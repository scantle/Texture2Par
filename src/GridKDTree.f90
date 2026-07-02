module m_grid_kdtree
  use kdtree2_module, only: kdkind, kdtree2, kdtree2_result, kdtree2_create, kdtree2_destroy, &
                            kdtree2_n_nearest, kdtree2_r_count, kdtree2_r_nearest
  implicit none
  private

  type, public :: t_grid_kdtree
    private
    type(kdtree2), pointer       :: handle => null()
    real(kdkind), pointer        :: data(:,:) => null()
    integer, allocatable         :: source_idx(:)
  contains
    procedure, public :: build => build_grid_kdtree
    procedure, public :: clear => clear_grid_kdtree
    procedure, public :: is_ready => grid_kdtree_is_ready
    procedure, public :: n_nearest => grid_kdtree_n_nearest
    procedure, public :: radius => grid_kdtree_radius
  end type t_grid_kdtree

contains

  subroutine build_grid_kdtree(this, coords, mask, rotmat)
    class(t_grid_kdtree), intent(inout) :: this
    real, intent(in)                    :: coords(:,:)
    logical, intent(in)                 :: mask(:)
    real, intent(in), optional          :: rotmat(:,:)
    integer                             :: i, j, n

    call this%clear()
    n = count(mask)

    ! KDTREE2 rejects data sets with fewer points than coordinate dimensions.
    ! Grid keeps its brute-force implementation for these very small searches.
    if (n < size(coords, 1)) return

    allocate(this%data(size(coords, 1), n), this%source_idx(n))
    j = 0
    do i = 1, size(mask)
      if (.not. mask(i)) cycle
      j = j + 1
      this%source_idx(j) = i
      if (present(rotmat)) then
        this%data(:,j) = matmul(rotmat, coords(:,i))
      else
        this%data(:,j) = coords(:,i)
      end if
    end do

    this%handle => kdtree2_create(this%data, sort=.false., rearrange=.true.)
  end subroutine build_grid_kdtree

  subroutine clear_grid_kdtree(this)
    class(t_grid_kdtree), intent(inout) :: this

    if (associated(this%handle)) call kdtree2_destroy(this%handle)
    nullify(this%handle)
    if (associated(this%data)) deallocate(this%data)
    nullify(this%data)
    if (allocated(this%source_idx)) deallocate(this%source_idx)
  end subroutine clear_grid_kdtree

  logical function grid_kdtree_is_ready(this)
    class(t_grid_kdtree), intent(in) :: this

    grid_kdtree_is_ready = associated(this%handle)
  end function grid_kdtree_is_ready

  subroutine grid_kdtree_n_nearest(this, qpoint, nnear, idx, dis)
    class(t_grid_kdtree), intent(inout) :: this
    real, intent(in)                    :: qpoint(:)
    integer, intent(in)                 :: nnear
    integer, intent(out)                :: idx(nnear)
    real, intent(out)                   :: dis(nnear)
    type(kdtree2_result), allocatable   :: results(:)
    real(kdkind), target                :: query(size(qpoint))
    integer                             :: i

    query = qpoint
    allocate(results(nnear))
    call kdtree2_n_nearest(this%handle, query, nnear, results)
    do i = 1, nnear
      idx(i) = this%source_idx(results(i)%idx)
      dis(i) = results(i)%dis
    end do
  end subroutine grid_kdtree_n_nearest

  subroutine grid_kdtree_radius(this, qpoint, r2, idx, dis)
    class(t_grid_kdtree), intent(inout) :: this
    real, intent(in)                    :: qpoint(:), r2
    integer, allocatable, intent(out)   :: idx(:)
    real, allocatable, intent(out)      :: dis(:)
    type(kdtree2_result), allocatable   :: results(:)
    real(kdkind), target                :: query(size(qpoint))
    integer                             :: i, nfound, nreported

    query = qpoint
    nfound = kdtree2_r_count(this%handle, query, real(r2, kdkind))
    allocate(idx(nfound), dis(nfound))
    if (nfound == 0) return

    allocate(results(nfound))
    call kdtree2_r_nearest(this%handle, query, real(r2, kdkind), nreported, nfound, results)
    do i = 1, nreported
      idx(i) = this%source_idx(results(i)%idx)
      dis(i) = results(i)%dis
    end do
  end subroutine grid_kdtree_radius

end module m_grid_kdtree
