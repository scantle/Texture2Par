module m_grid
  use m_global, only: dp
  use kdtree2_module
  implicit none
  
!-------------------------------------------------------------------------------------------------!
  type t_Grid
  
    integer                   :: n
    integer                   :: dim
    real(dp), allocatable     :: coords(:,:)    ! xy by n points
    type(kdtree2),pointer     :: tree           ! Tree for locating model nodes
    integer, allocatable      :: true_idx(:)    ! Mapping of indices from masked to full dataset
    
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
    implicit none
    class(t_Grid), intent(inout) :: this
    logical,intent(in),optional  :: mask(this%n)
    real(dp),optional            :: rotmat(:,:) ! Rotation Matrix
    integer                      :: j
    
    allocate(this%tree)
    
    if (present(mask)) then
      allocate(this%true_idx(count(mask)))
      this%true_idx = pack([(j, j = 1, this%n)], mask)
      this%tree => kdtree2_create(this%coords(:,this%true_idx), dim=this%dim, rotmat=rotmat)
    else
      this%tree => kdtree2_create(this%coords, dim=this%dim, rotmat=rotmat)
      allocate(this%true_idx(this%n))
      this%true_idx = [(j, j = 1, this%n)]
    end if
    
  end subroutine build_tree
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine get_nnear(this, qpoint, nnear, idx, dis, mask, rotmat)
    implicit none
    ! Uses tree to return n nearest grid indeces and distances
    class(t_Grid), intent(inout) :: this
    real(dp),intent(in)          :: qpoint(this%dim)
    integer, intent(in)          :: nnear
    integer, intent(out)         :: idx(nnear)
    real(dp),intent(out)         :: dis(nnear)
    type(kdtree2_result)         :: res(nnear)
    logical,intent(in),optional  :: mask(this%n)
    real(kdkind), target, optional  :: rotmat(:,:) ! Rotation Matrix

    if (associated(this%tree)==.false.) call this%build_tree(mask=mask)
    call kdtree2_n_nearest(this%tree, qpoint, nnear, res, rotmat)
    
    ! Worth getting fancier than copying? Pointers?
    if (allocated(this%true_idx)) then
      ! Map indices back to original dataset
      idx = this%true_idx(res(1:)%idx)
    else
      idx = res(:)%idx
    end if
    dis = res(:)%dis
    
  end subroutine get_nnear
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine get_r2near(this, qpoint, r2, idx, dis, mask)
    implicit none
    class(t_Grid), intent(inout)     :: this
    real(dp),intent(in)              :: qpoint(this%dim)
    real(dp),intent(in)              :: r2
    integer, pointer                 :: idx(:)
    real(dp),pointer                 :: dis(:)
    integer                          :: nfound
    type(kdtree2_result),allocatable :: res(:)
    logical,intent(in),optional      :: mask(this%n)
    
    if (associated(this%tree)==.false.) call this%build_tree(mask=mask)
    nfound = kdtree2_r_count(this%tree,qpoint,r2)
    if (associated(idx)) deallocate(idx)
    if (associated(dis)) deallocate(dis)
    allocate(idx(nfound), dis(nfound), res(nfound))
    call kdtree2_r_nearest(this%tree, qpoint, r2, nfound, nfound, res)
    
    ! Worth getting fancier than copying? Pointers?
    if (allocated(this%true_idx)) then
      ! Map indices back to original dataset
      idx = this%true_idx(res(:)%idx)
    else
      idx = res(:)%idx
    end if
    dis = res(:)%dis
    
  end subroutine get_r2near
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine burn_tree(this)
    ! sometimes you have to have a little fun with names
    implicit none
    class(t_Grid), intent(inout)     :: this
    if (associated(this%tree)) deallocate(this%tree)
  end subroutine burn_tree
  
!-------------------------------------------------------------------------------------------------!

  subroutine finalize(this)
    implicit none
    class(t_Grid), intent(inout) :: this
  
    if (allocated(this%coords)) deallocate(this%coords)
    if (associated(this%tree)) deallocate(this%tree)
  
    ! Include any other cleanup tasks here
  end subroutine finalize
  
!-------------------------------------------------------------------------------------------------!
end module m_grid