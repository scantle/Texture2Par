module m_pilotpoints
  use m_datasets, only: t_dataset
  use m_global, only: dp
  use m_grid, only: t_grid, create_grid_object
  implicit none
  
  ! Module Variables
  integer,parameter           :: nproperties = 9  ! Kmin, Kmax, Ss, Sy, Aniso, Kd, KHp, KVp, STp
  
!-------------------------------------------------------------------------------------------------!
  type, extends(t_dataset) :: t_pilotpoints
    ! For reference: inherited from t_dataset
    !integer,allocatable           :: category(:)      ! Zones, Hydrostratigraphic Units (HSUs)
    !type(t_grid),pointer          :: grid             ! Grid
  
    integer                        :: class_id         ! Class ID
    integer,allocatable            :: lay_type(:)      ! Aquifer (1), Aquitard (-1), Unused (0)
    real(dp),allocatable           :: values(:,:)      ! Data values associated with each point
    
    contains
    
      procedure,public             :: initialize => init_pilotpoint
  
  end type t_pilotpoints
!-------------------------------------------------------------------------------------------------!

  contains
  
!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!

  function create_pilotpoint_object(class_id, npoints) result(object)
    implicit none
    
    integer,intent(in)               :: class_id
    integer,intent(inout)            :: npoints
    class(t_pilotpoints),pointer     :: object
    
    allocate(object)
    
    call object%initialize(npoints)
  
  end function
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! LAYER CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine init_pilotpoint(this, npoints)
    implicit none
    class(t_pilotpoints),intent(inout) :: this
    integer,intent(in)                 :: npoints
    
    allocate(this%category  (npoints), &
             this%lay_type  (npoints), &
             this%values    (nproperties, npoints))
    allocate(this%grid)

    ! Intialize with default values
    this%category = 1
    this%lay_type = 0
    
    this%grid => create_grid_object(2, npoints)  ! X,Y coords by n 
  
  end subroutine init_pilotpoint
  
!-------------------------------------------------------------------------------------------------!  
end module m_pilotpoints