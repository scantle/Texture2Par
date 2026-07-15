module m_grid
  use m_cellblocks
  use m_grid_kdtree, only: t_grid_kdtree
  implicit none

!-------------------------------------------------------------------------------------------------!
  type t_grid_search
    integer              :: n = 0
    logical, allocatable :: mask(:)
    type(t_grid_kdtree)  :: kdtree
  end type t_grid_search

!-------------------------------------------------------------------------------------------------!
  type t_Grid

    integer                   :: n
    integer                   :: dim
    real    , allocatable     :: coords(:,:)         ! xy[z] by n points
    type(t_grid_search),pointer :: tree => null()    ! Built search mask/rotation state

    ! Block Kriging Support
    type(t_cellblocks), pointer :: active_blocks => null()

    contains
      procedure,private        :: initialize
      procedure,private        :: finalize
      procedure,public         :: build_tree
      procedure,public         :: get_nnear
      procedure,public         :: burn_tree
  end type t_Grid
!-------------------------------------------------------------------------------------------------!
  contains
!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!

  function create_grid_object(dim, n) result(object)
    implicit none

    integer,intent(in)            :: dim, n
    class(t_Grid),pointer         :: object

    allocate(object)
    call object%initialize(dim, n)

  end function

!-------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------!
! GRID CLASS FILE TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

  subroutine initialize(this, dim, n)
    implicit none

    class(t_Grid)        :: this
    integer, intent(in)  :: dim, n

    this%n   = n
    this%dim = dim

    allocate(this%coords(dim,n))

    this%coords = 0.0

  end subroutine initialize

!-------------------------------------------------------------------------------------------------!

  subroutine build_tree(this, mask)
    use m_error_handler, only: error_handler
    implicit none
    class(t_Grid), intent(inout) :: this
    logical,intent(in),optional  :: mask(this%n)

    call this%burn_tree()
    allocate(this%tree)
    allocate(this%tree%mask(this%n))

    if (present(mask)) then
      this%tree%mask = mask
    else
      this%tree%mask = .true.
    end if
    this%tree%n = count(this%tree%mask)

    if (this%tree%n == 0) then
      call error_handler(9, opt_msg='Could not initialize search - no eligible data')
    end if

    call this%tree%kdtree%build(this%coords, this%tree%mask)

  end subroutine build_tree

!-------------------------------------------------------------------------------------------------!

  subroutine get_nnear(this, qpoint, nnear, idx, dis)
    use m_error_handler, only: error_handler
    implicit none
    ! Returns source indices and SQUARED distances, ordered by increasing
    ! distance with lower source index breaking exact ties.
    class(t_Grid), intent(inout) :: this
    real    ,intent(in)          :: qpoint(this%dim)
    integer, intent(in)          :: nnear
    integer, intent(out)         :: idx(nnear)
    real    ,intent(out)         :: dis(nnear)
    integer, allocatable         :: tied_idx(:)
    real, allocatable            :: tied_dis(:)
    real                         :: query(this%dim), cutoff
    integer                      :: neligible

    if (.not. associated(this%tree)) call this%build_tree()

    neligible = count(this%tree%mask)
    if (nnear < 1 .or. nnear > neligible) then
      call error_handler(9, opt_msg='Requested nearest-neighbor count exceeds eligible data')
    end if

    if (.not. this%tree%kdtree%is_ready()) then
      call get_nnear_brute_force(this, qpoint, nnear, idx, dis)
      return
    end if

    query = qpoint
    call this%tree%kdtree%n_nearest(query, nnear, idx, dis)
    call sort_neighbors(idx, dis)
    cutoff = dis(nnear)

    call this%tree%kdtree%radius(query, cutoff, tied_idx, tied_dis)
    call sort_neighbors(tied_idx, tied_dis)
    idx = tied_idx(:nnear)
    dis = tied_dis(:nnear)

  end subroutine get_nnear

!-------------------------------------------------------------------------------------------------!

  subroutine burn_tree(this)
    implicit none
    class(t_Grid), intent(inout)     :: this

    if (associated(this%tree)) then
      call this%tree%kdtree%clear()
      deallocate(this%tree)
    end if
  end subroutine burn_tree

!-------------------------------------------------------------------------------------------------!

  subroutine finalize(this)
    implicit none
    class(t_Grid), intent(inout) :: this

    if (allocated(this%coords))         deallocate(this%coords)
    if (associated(this%tree))          call this%burn_tree()
    if (associated(this%active_blocks)) deallocate(this%active_blocks)

    ! Include any other cleanup tasks here
  end subroutine finalize

!-------------------------------------------------------------------------------------------------!

  function squared_distance(point, qpoint) result(distance)
    real, intent(in)           :: point(:), qpoint(:)
    real                       :: distance
    real                       :: delta(size(point))

    delta = point - qpoint
    distance = sum(delta**2)
  end function squared_distance

!-------------------------------------------------------------------------------------------------!

  subroutine get_nnear_brute_force(this, qpoint, nnear, idx, dis)
    class(t_Grid), intent(in)       :: this
    real, intent(in)                :: qpoint(this%dim)
    integer, intent(in)             :: nnear
    integer, intent(out)            :: idx(nnear)
    real, intent(out)               :: dis(nnear)
    logical, allocatable            :: selected(:)
    integer                         :: i, j, best_idx
    real                            :: candidate_dis, best_dis

    allocate(selected(this%n))
    selected = .false.
    do j = 1, nnear
      best_idx = 0
      best_dis = huge(0.0)
      do i = 1, this%n
        if (selected(i) .or. .not. this%tree%mask(i)) cycle
        candidate_dis = squared_distance(this%coords(:,i), qpoint)
        if (candidate_dis < best_dis .or. &
            (candidate_dis == best_dis .and. (best_idx == 0 .or. i < best_idx))) then
          best_idx = i
          best_dis = candidate_dis
        end if
      end do
      selected(best_idx) = .true.
      idx(j) = best_idx
      dis(j) = best_dis
    end do
  end subroutine get_nnear_brute_force

!-------------------------------------------------------------------------------------------------!

  subroutine sort_neighbors(idx, dis)
    integer, intent(inout) :: idx(:)
    real, intent(inout)    :: dis(:)
    integer                :: i, j, idx_value
    real                   :: dis_value

    do i = 2, size(idx)
      idx_value = idx(i)
      dis_value = dis(i)
      j = i - 1
      do while (j >= 1)
        if (dis(j) < dis_value) exit
        if (dis(j) == dis_value .and. idx(j) < idx_value) exit
        idx(j+1) = idx(j)
        dis(j+1) = dis(j)
        j = j - 1
      end do
      idx(j+1) = idx_value
      dis(j+1) = dis_value
    end do
  end subroutine sort_neighbors

!-------------------------------------------------------------------------------------------------!
end module m_grid
