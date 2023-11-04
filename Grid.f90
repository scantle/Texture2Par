module m_grid
  implicit none
  
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,307)  ! 64-bit
  
!-------------------------------------------------------------------------------------------------!
  type t_Grid
  
    integer                   :: n
    integer                   :: dim
    real(dp), allocatable     :: xyz(:,:)  ! xyz by n points
    
    contains
      procedure,public        :: initialize
      !procedure,public        :: build_tree
  end type t_Grid
!-------------------------------------------------------------------------------------------------!
  
  contains
  
!-------------------------------------------------------------------------------------------------!
! GRID CLASS FILE TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine initialize(this, dim, n)
    implicit none
    
    class(t_Grid)        :: this
    integer, intent(in)  :: dim, n
  
    this%n   = n
    this%dim = dim
    
    allocate(this%xyz(dim,n))
    
  end subroutine initialize
  
!-------------------------------------------------------------------------------------------------!
end module m_grid