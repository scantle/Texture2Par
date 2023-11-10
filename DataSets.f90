module m_datasets
  use m_global, only: dp
  use m_grid, only: t_grid
  implicit none
!-------------------------------------------------------------------------------------------------!
  type :: t_dataset
    
    character(30)                  :: name             ! Texture class
    integer,allocatable            :: hsus(:)          ! Hydrostratigraphic Units for each data point
    integer,allocatable            :: nearest_node(:)  ! Nearest model node (cell) to data point
    real(dp),allocatable           :: values(:)        ! Data values associated with each point
    type(t_grid),pointer           :: grid
    
    contains
    
      procedure,public             :: initialize
  
  end type t_dataset
!-------------------------------------------------------------------------------------------------!

  contains
  
!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  function create_dataset_object(class_name, npoints) result(object)
    implicit none
    
    character(*),intent(in)          :: class_name
    integer,intent(inout)            :: npoints
    class(t_dataset),pointer         :: object
    
    call object%initialize(npoints)
    object%name = class_name
  
  end function
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine initialize(this, npoints)
    implicit none
    class(t_dataset),intent(inout)  :: this
    integer,intent(in)              :: npoints
    
    allocate(this%hsus         (npoints), &
             this%nearest_node (npoints), &
             this%values       (npoints)  )

    !this%grid => create_grid_object(3, npoints)  ! X,Y,Z coords by n 
  
  end subroutine initialize
  
!-------------------------------------------------------------------------------------------------!
end module m_datasets