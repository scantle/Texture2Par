module m_grid
  use m_global, only: dp
  use kdtree2_module
  implicit none
  
!-------------------------------------------------------------------------------------------------!
  type t_Grid
  
    integer                   :: n
    integer                   :: dim
    real(dp), allocatable     :: coords(:,:)  ! xy by n points
    type(kdtree2),pointer     :: tree           ! Tree for locating model nodes
    
    contains
      procedure,private        :: initialize
      procedure,private        :: finalize
      procedure,private        :: build_tree
      procedure,public         :: get_nnear
      procedure,public         :: get_r2near
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
  
  subroutine build_tree(this)
    class(t_Grid), intent(inout) :: this
    
    allocate(this%tree)
    this%tree => kdtree2_create(this%coords, dim=this%dim)
    
  end subroutine build_tree
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine get_nnear(this, qpoint, nnear, idx, dis)
    ! Uses tree to return n nearest grid indeces and distances
    class(t_Grid), intent(inout) :: this
    real(dp),intent(in)          :: qpoint(this%dim)
    integer, intent(in)          :: nnear
    integer, intent(out)         :: idx(nnear)
    real(dp),intent(out)         :: dis(nnear)
    type(kdtree2_result)         :: res(nnear)
    
    if (associated(this%tree)==.false.) call this%build_tree()
    
    call kdtree2_n_nearest(this%tree, qpoint, nnear, res)
    
    ! Worth getting fancier than copying? Pointers?
    idx = res(:)%idx
    dis = res(:)%dis
    
  end subroutine get_nnear
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine get_r2near(this, qpoint, r2, idx, dis)
    class(t_Grid), intent(inout)     :: this
    real(dp),intent(in)              :: qpoint(this%dim)
    real(dp),intent(in)              :: r2
    integer, pointer                 :: idx(:)
    real(dp),pointer                 :: dis(:)
    integer                          :: nfound
    type(kdtree2_result),allocatable :: res(:)
    
    if (associated(this%tree)==.false.) call this%build_tree()
    nfound = kdtree2_r_count(this%tree,qpoint,r2)
    if (associated(idx)) deallocate(idx)
    if (associated(dis)) deallocate(dis)
    allocate(idx(nfound), dis(nfound), res(nfound))
    call kdtree2_r_nearest(this%tree, qpoint, r2, nfound, nfound, res)
    
    idx = res(:)%idx
    dis = res(:)%dis
    
  end subroutine get_r2near
  
!-------------------------------------------------------------------------------------------------!

  subroutine finalize(this)
    class(t_Grid), intent(inout) :: this
  
    if (allocated(this%coords)) deallocate(this%coords)
    if (associated(this%tree)) deallocate(this%tree)
  
    ! Include any other cleanup tasks here
  end subroutine finalize
  
!-------------------------------------------------------------------------------------------------!
end module m_grid