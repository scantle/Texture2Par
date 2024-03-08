module m_grid
  use m_global, only: dp
  implicit none
  
!-------------------------------------------------------------------------------------------------!
  type t_Grid
  
    integer                   :: n
    integer                   :: dim
    real(dp), allocatable     :: coords(:,:)  ! xy by n points
    
    contains
      procedure,private        :: initialize
      procedure,private        :: finalize
      !procedure,public        :: build_tree
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

  subroutine finalize(this)
    class(t_Grid), intent(inout) :: this
  
    ! Check if 'coords' is allocated and deallocate it
    if (allocated(this%coords)) then
      deallocate(this%coords)
    end if
  
    ! Include any other cleanup tasks here
  end subroutine finalize
  
!-------------------------------------------------------------------------------------------------!
end module m_grid