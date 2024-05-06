module m_datasets
  use m_global, only: dp
  use m_grid, only: t_grid, create_grid_object
  implicit none
!-------------------------------------------------------------------------------------------------!
  type, abstract  :: t_dataset  ! Abstract Base Class
    
    integer,allocatable            :: category(:)      ! Zones, Hydrostratigraphic Units (HSUs)
    type(t_grid),pointer           :: grid             ! Grid
    ! Doesn't even have any data! Because array dimensions differ, subclasses are controversially
    ! allowed to have their own data storage
    
    contains
      procedure,nopass             :: initialize_abc
      generic,public               :: initialize => initialize_abc
  
  end type t_dataset
!-------------------------------------------------------------------------------------------------!
!  type, extends(t_dataset) :: t_intervaldata
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
!  end type t_intervaldata
!-------------------------------------------------------------------------------------------------!
  type, extends(t_dataset) :: t_layerdata
    ! Intended for input data
    character(30)                  :: name             ! Class Name
    integer                        :: class_id         ! Class ID
    real(dp),allocatable           :: values(:,:)        ! Data values associated with each grid point for each class (nclasses, npoints)
    integer,allocatable            :: fmodel_loc(:)    ! Containing model element/cell
    
    contains
      procedure                    :: init_layerdata
      generic                      :: initialize => init_layerdata
  
  end type t_layerdata
!-------------------------------------------------------------------------------------------------!
!  type, extends(t_dataset) :: t_layermultidata
!    ! Intended for intermediate data
!    character(30)                  :: name             ! Class Name
!    integer                        :: class_id         ! Class ID
!    real(dp),allocatable           :: values(:,:)      ! Data values associated with each grid point for each class (nclasses, npoints)
!    integer,allocatable            :: fmodel_loc(:)    ! Containing model element/cell
!    
!    contains
!      procedure                    :: init_layermultidata
!      generic                      :: initialize => init_layermultidata
!  
!  end type t_layermultidata
!-------------------------------------------------------------------------------------------------!
! See also: t_pilotpoint derived class in PilotPoints.f90
!-------------------------------------------------------------------------------------------------!

  contains
  
!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------!
! BASE CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

  subroutine initialize_abc(this)
    implicit none
    class(t_dataset),intent(inout)  :: this
  
  end subroutine initialize_abc

!-------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------!
! LAYER CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine init_layerdata(this, npoints, nclasses)
    implicit none
    class(t_layerdata),intent(inout)  :: this
    integer,intent(in)                :: npoints, nclasses
    
    allocate(this%category  (npoints), &
             this%values    (nclasses, npoints))
    allocate(this%grid)

    !this%grid = create_grid_object(2, npoints)  ! X,Y,Z coords by n 
  
  end subroutine init_layerdata
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! LAYER MULTI-DATA CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
!  
!  subroutine init_layermultidata(this, npoints, nclasses, grid)
!    implicit none
!    class(t_layermultidata),intent(inout)  :: this
!    integer,intent(in)                     :: npoints, nclasses
!    class(t_grid),pointer                  :: grid
!    
!    allocate(this%category  (npoints), &
!             this%values    (nclasses, npoints))
!    !allocate(this%grid)  ! Already should exist, avoid copying
!    this%grid => grid
!
!    !this%grid = create_grid_object(2, npoints)  ! X,Y,Z coords by n 
!  
!  end subroutine init_layermultidata
!  
!-------------------------------------------------------------------------------------------------!
  
end module m_datasets