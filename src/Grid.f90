module m_grid
  use m_cellblocks
  implicit none

!-------------------------------------------------------------------------------------------------!
  type t_grid_search
    ! Temporary brute-force search state. The name `tree` is retained on
    ! t_grid so existing callers can continue to test association and read n.
    integer              :: n = 0
    logical, allocatable :: mask(:)
    real, allocatable    :: rotmat(:,:)
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
      procedure,public         :: get_r2near
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

  subroutine build_tree(this, mask, rotmat)
    use m_error_handler, only: error_handler
    implicit none
    class(t_Grid), intent(inout) :: this
    logical,intent(in),optional  :: mask(this%n)
    real    ,intent(in),optional :: rotmat(:,:) ! Legacy Rotation Matrix

    call this%burn_tree()
    allocate(this%tree)
    allocate(this%tree%mask(this%n))

    if (present(mask)) then
      this%tree%mask = mask
    else
      this%tree%mask = .true.
    end if
    this%tree%n = count(this%tree%mask)

    if (present(rotmat)) then
      if (size(rotmat,1) /= this%dim .or. size(rotmat,2) /= this%dim) then
        call error_handler(9, opt_msg='Search rotation matrix does not match grid dimension')
      end if
      allocate(this%tree%rotmat(this%dim,this%dim))
      this%tree%rotmat = rotmat
    end if

    if (this%tree%n == 0) then
      call error_handler(9, opt_msg='Could not initialize search - no eligible data')
    end if

  end subroutine build_tree

!-------------------------------------------------------------------------------------------------!

  subroutine get_nnear(this, qpoint, nnear, idx, dis, mask, rotmat)
    use m_error_handler, only: error_handler
    implicit none
    ! Returns source indices and SQUARED distances, ordered by increasing
    ! distance with lower source index breaking exact ties.
    class(t_Grid), intent(inout) :: this
    real    ,intent(in)          :: qpoint(this%dim)
    integer, intent(in)          :: nnear
    integer, intent(out)         :: idx(nnear)
    real    ,intent(out)         :: dis(nnear)
    logical,intent(in),optional  :: mask(this%n)
    real    ,intent(in),optional :: rotmat(:,:) ! Rotation Matrix
    logical, allocatable         :: selected(:)
    integer                      :: i, j, best_idx, neligible
    real                         :: candidate_dis, best_dis

    if (.not. associated(this%tree)) call this%build_tree(mask=mask, rotmat=rotmat)

    neligible = count(this%tree%mask)
    if (present(mask)) neligible = count(this%tree%mask .and. mask)
    if (nnear < 1 .or. nnear > neligible) then
      call error_handler(9, opt_msg='Requested nearest-neighbor count exceeds eligible data')
    end if

    allocate(selected(this%n))
    selected = .false.
    do j = 1, nnear
      best_idx = 0
      best_dis = huge(0.0)
      do i = 1, this%n
        if (selected(i) .or. .not. this%tree%mask(i)) cycle
        if (present(mask)) then
          if (.not. mask(i)) cycle
        end if

        if (present(rotmat)) then
          candidate_dis = squared_distance(this%coords(:,i), qpoint, rotmat)
        else if (allocated(this%tree%rotmat)) then
          candidate_dis = squared_distance(this%coords(:,i), qpoint, this%tree%rotmat)
        else
          candidate_dis = squared_distance(this%coords(:,i), qpoint)
        end if

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

  end subroutine get_nnear

!-------------------------------------------------------------------------------------------------!

  subroutine get_r2near(this, qpoint, r2, idx, dis, mask)
    implicit none
    ! Returns all eligible points within the SQUARED distance r2.
    class(t_Grid), intent(inout)     :: this
    real    ,intent(in)              :: qpoint(this%dim)
    real    ,intent(in)              :: r2
    integer, pointer                 :: idx(:)
    real    ,pointer                 :: dis(:)
    logical,intent(in),optional      :: mask(this%n)
    integer                          :: i, nfound
    real                             :: candidate_dis

    if (.not. associated(this%tree)) call this%build_tree(mask=mask)

    nfound = 0
    do i = 1, this%n
      if (.not. this%tree%mask(i)) cycle
      if (present(mask)) then
        if (.not. mask(i)) cycle
      end if
      if (allocated(this%tree%rotmat)) then
        candidate_dis = squared_distance(this%coords(:,i), qpoint, this%tree%rotmat)
      else
        candidate_dis = squared_distance(this%coords(:,i), qpoint)
      end if
      if (candidate_dis <= r2) nfound = nfound + 1
    end do

    if (associated(idx)) deallocate(idx)
    if (associated(dis)) deallocate(dis)
    allocate(idx(nfound), dis(nfound))

    nfound = 0
    do i = 1, this%n
      if (.not. this%tree%mask(i)) cycle
      if (present(mask)) then
        if (.not. mask(i)) cycle
      end if
      if (allocated(this%tree%rotmat)) then
        candidate_dis = squared_distance(this%coords(:,i), qpoint, this%tree%rotmat)
      else
        candidate_dis = squared_distance(this%coords(:,i), qpoint)
      end if
      if (candidate_dis <= r2) then
        nfound = nfound + 1
        idx(nfound) = i
        dis(nfound) = candidate_dis
      end if
    end do
    call sort_neighbors(idx, dis)

  end subroutine get_r2near

!-------------------------------------------------------------------------------------------------!

  subroutine burn_tree(this)
    implicit none
    class(t_Grid), intent(inout)     :: this

    if (associated(this%tree)) deallocate(this%tree)
  end subroutine burn_tree

!-------------------------------------------------------------------------------------------------!

  subroutine finalize(this)
    implicit none
    class(t_Grid), intent(inout) :: this

    if (allocated(this%coords))         deallocate(this%coords)
    if (associated(this%tree))          deallocate(this%tree)
    if (associated(this%active_blocks)) deallocate(this%active_blocks)

    ! Include any other cleanup tasks here
  end subroutine finalize

!-------------------------------------------------------------------------------------------------!

  function squared_distance(point, qpoint, rotmat) result(distance)
    real, intent(in)           :: point(:), qpoint(:)
    real, intent(in), optional :: rotmat(:,:)
    real                       :: distance
    real                       :: delta(size(point))

    delta = point - qpoint
    if (present(rotmat)) delta = matmul(rotmat, delta)
    distance = sum(delta**2)
  end function squared_distance

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
