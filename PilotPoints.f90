module m_pilotpoints
  use m_datasets, only: t_dataset
  use m_global, only: dp, NODATA
  use m_grid, only: t_grid, create_grid_object
  implicit none
  
  ! Module Variables
  integer,parameter           :: nproperties = 9  ! Kmin, Kmax, Ss, Sy, Aniso, Kd, KHp, KVp, STp
  
!-------------------------------------------------------------------------------------------------!
  type, extends(t_dataset) :: t_pilotpoints
    ! For reference: inherited from t_dataset
    !integer,allocatable           :: category(:)      ! Zones, Hydrostratigraphic Units (HSUs)
    !type(t_grid),pointer          :: grid             ! Grid
  
    integer,allocatable            :: lay_type(:)      ! Aquifer (1), Aquitard (2), Unused (0)
    integer,allocatable            :: zonelist(:)      ! List of all "categories" assigned to a pilot point
    real(dp),allocatable           :: values(:,:,:)    ! Data values by parameter (point, parameter, class_id)
    
    ! Counters
    integer                        :: max_a_cat        ! Max aquifer category id
    integer                        :: pp_split         ! First aquitard point index
    integer                        :: zonelist_split   ! The first aquitard zone in zonelist
    
    contains
    
      procedure,public             :: initialize => init_pilotpoint
  
  end type t_pilotpoints
!-------------------------------------------------------------------------------------------------!

  contains
  
!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!

  function create_pilotpoint_object(npoints, nclasses) result(object)
    implicit none
    
    integer,intent(inout)            :: npoints,nclasses
    class(t_pilotpoints),pointer     :: object
    
    allocate(object)
    
    call object%initialize(npoints)
    
    ! Have to allocate this seperately, since init cannot take more than npoints (inheritance)
    allocate(object%values(npoints, nproperties, nclasses))
    object%values   = NODATA
  
  end function
  
!-------------------------------------------------------------------------------------------------!
! LAYER CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine init_pilotpoint(this, npoints)
    implicit none
    class(t_pilotpoints),intent(inout) :: this
    integer,intent(in)                 :: npoints
    
    allocate(this%category  (npoints), &
             this%lay_type  (npoints)  )
    allocate(this%grid)

    ! Intialize with default values
    this%category       = 1
    this%lay_type       = 0
    this%max_a_cat      = 0
    this%pp_split       = 0
    this%zonelist_split = 0
    
    this%grid => create_grid_object(2, npoints)  ! X,Y coords by n 
  
  end subroutine init_pilotpoint
  
!-------------------------------------------------------------------------------------------------!  
end module m_pilotpoints