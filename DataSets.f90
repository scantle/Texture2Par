module m_datasets
  use m_global, only: dp
  use m_grid, only: t_grid
  implicit none
!-------------------------------------------------------------------------------------------------!
  type :: t_dataset  ! Abstract Base Class
    
    integer,allocatable            :: category(:)      ! Zones, Hydrostratigraphic Units (HSUs)
    type(t_grid),pointer           :: grid             ! Grid
    
    contains
    
      procedure,public             :: initialize
  
  end type t_dataset
!-------------------------------------------------------------------------------------------------!
!  type, extends(t_dataset) :: t_intervaldata   ! UNSURE IF NEEDED
!    
!    character(30)                  :: name             ! Texture class
!    integer                        :: class_id         ! Class ID
!    integer,allocatable            :: nearest_node(:)  ! Nearest model node (cell) to data point
!    real(dp),allocatable           :: values(:,:)      ! Data values associated with intervals at each point (interval, point)
!    
!    contains
!    
!      procedure,public             :: initialize => init_interval
!  
!  end type t_dataset
!-------------------------------------------------------------------------------------------------!
  type, extends(t_dataset) :: t_layerdata
    
    character(30)                  :: name             ! Class Name
    integer                        :: class_id         ! Class ID
    real(dp),allocatable           :: values(:)        ! Data values associated with each grid point
    !integer,allocatable            :: nearest_node(:)  ! Nearest model node (cell) to data point
    
    contains
    
      procedure,public             :: initialize => init_layerdata
  
  end type t_layerdata
!-------------------------------------------------------------------------------------------------!
! See also: t_pilotpoint derived class in PilotPoints.f90
!-------------------------------------------------------------------------------------------------!

  contains
  
!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  function create_dataset_object(npoints) result(object)
    implicit none
    
    integer,intent(inout)            :: npoints
    class(t_dataset),pointer         :: object
    
    call object%initialize(npoints)
  
  end function
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! BASE CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine initialize(this, npoints)
    implicit none
    class(t_dataset),intent(inout)  :: this
    integer,intent(in)              :: npoints
    
    allocate(this%category   (npoints))
    allocate(this%grid)

    !this%grid => create_grid_object(3, npoints)  ! X,Y,Z coords by n 
  
  end subroutine initialize
  
!-------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------!
! LAYER CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine init_layerdata(this, npoints)
    implicit none
    class(t_layerdata),intent(inout)  :: this
    integer,intent(in)                :: npoints
    
    allocate(this%category  (npoints), &
             this%values    (npoints)  )
    allocate(this%grid)

    !this%grid => create_grid_object(3, npoints)  ! X,Y,Z coords by n 
  
  end subroutine init_layerdata
  
!-------------------------------------------------------------------------------------------------!
  
end module m_datasets