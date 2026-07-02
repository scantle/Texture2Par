module m_grid
  use m_cellblocks
  use m_grid_kdtree, only: t_grid_kdtree
  implicit none

!-------------------------------------------------------------------------------------------------!
  type t_grid_search
    integer              :: n = 0
    logical, allocatable :: mask(:)
    real, allocatable    :: rotmat(:,:)
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

    if (allocated(this%tree%rotmat)) then
      call this%tree%kdtree%build(this%coords, this%tree%mask, this%tree%rotmat)
    else
      call this%tree%kdtree%build(this%coords, this%tree%mask)
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
    integer, allocatable         :: candidate_idx(:), tied_idx(:), valid_idx(:)
    real, allocatable            :: candidate_dis(:), tied_dis(:), valid_dis(:)
    real                         :: query(this%dim), cutoff
    integer                      :: ncandidate, neligible, nvalid

    if (.not. associated(this%tree)) call this%build_tree(mask=mask, rotmat=rotmat)

    neligible = count(this%tree%mask)
    if (present(mask)) neligible = count(this%tree%mask .and. mask)
    if (nnear < 1 .or. nnear > neligible) then
      call error_handler(9, opt_msg='Requested nearest-neighbor count exceeds eligible data')
    end if

    ! A query-time rotation can differ from the metric used to build the tree.
    ! Preserve that legacy override with the exact brute-force implementation.
    if (present(rotmat) .or. .not. this%tree%kdtree%is_ready()) then
      call get_nnear_brute_force(this, qpoint, nnear, idx, dis, mask, rotmat)
      return
    end if

    query = qpoint
    if (allocated(this%tree%rotmat)) query = matmul(this%tree%rotmat, qpoint)

    ! Dynamic masks (notably SGSIM and FlowModel active cells) are handled by
    ! progressively expanding the unmasked candidate set. Once enough eligible
    ! points are present, a radius query captures every exact tie at the cutoff.
    ncandidate = min(this%tree%n, max(nnear, 2*nnear))
    do
      allocate(candidate_idx(ncandidate), candidate_dis(ncandidate))
      call this%tree%kdtree%n_nearest(query, ncandidate, candidate_idx, candidate_dis)
      if (present(mask)) then
        nvalid = count(mask(candidate_idx))
      else
        nvalid = ncandidate
      end if
      if (nvalid >= nnear) exit
      deallocate(candidate_idx, candidate_dis)
      ncandidate = min(this%tree%n, 2*ncandidate)
    end do

    allocate(valid_idx(nvalid), valid_dis(nvalid))
    if (present(mask)) then
      valid_idx = pack(candidate_idx, mask(candidate_idx))
      valid_dis = pack(candidate_dis, mask(candidate_idx))
    else
      valid_idx = candidate_idx
      valid_dis = candidate_dis
    end if
    call sort_neighbors(valid_idx, valid_dis)
    cutoff = valid_dis(nnear)

    call this%tree%kdtree%radius(query, cutoff, tied_idx, tied_dis)
    if (present(mask)) then
      nvalid = count(mask(tied_idx))
      valid_idx = pack(tied_idx, mask(tied_idx))
      valid_dis = pack(tied_dis, mask(tied_idx))
    else
      nvalid = size(tied_idx)
      valid_idx = tied_idx
      valid_dis = tied_dis
    end if
    call sort_neighbors(valid_idx, valid_dis)
    idx = valid_idx(:nnear)
    dis = valid_dis(:nnear)

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
    integer, allocatable             :: found_idx(:)
    real, allocatable                :: found_dis(:)
    real                             :: query(this%dim)
    integer                          :: nfound

    if (.not. associated(this%tree)) call this%build_tree(mask=mask)

    if (associated(idx)) deallocate(idx)
    if (associated(dis)) deallocate(dis)

    if (.not. this%tree%kdtree%is_ready()) then
      call get_r2near_brute_force(this, qpoint, r2, idx, dis, mask)
      return
    end if

    query = qpoint
    if (allocated(this%tree%rotmat)) query = matmul(this%tree%rotmat, qpoint)
    call this%tree%kdtree%radius(query, r2, found_idx, found_dis)
    if (present(mask)) then
      nfound = count(mask(found_idx))
      allocate(idx(nfound), dis(nfound))
      idx = pack(found_idx, mask(found_idx))
      dis = pack(found_dis, mask(found_idx))
    else
      allocate(idx(size(found_idx)), dis(size(found_dis)))
      idx = found_idx
      dis = found_dis
    end if
    call sort_neighbors(idx, dis)

  end subroutine get_r2near

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

  subroutine get_nnear_brute_force(this, qpoint, nnear, idx, dis, mask, rotmat)
    class(t_Grid), intent(in)       :: this
    real, intent(in)                :: qpoint(this%dim)
    integer, intent(in)             :: nnear
    integer, intent(out)            :: idx(nnear)
    real, intent(out)               :: dis(nnear)
    logical, intent(in), optional   :: mask(this%n)
    real, intent(in), optional      :: rotmat(:,:)
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
  end subroutine get_nnear_brute_force

!-------------------------------------------------------------------------------------------------!

  subroutine get_r2near_brute_force(this, qpoint, r2, idx, dis, mask)
    class(t_Grid), intent(in)       :: this
    real, intent(in)                :: qpoint(this%dim), r2
    integer, pointer                :: idx(:)
    real, pointer                   :: dis(:)
    logical, intent(in), optional   :: mask(this%n)
    integer                         :: i, nfound
    real                            :: candidate_dis

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
  end subroutine get_r2near_brute_force

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
