module m_pilotpoints
  use m_datasets, only: t_dataset
  use m_global, only: dp, NODATA
  use m_grid, only: t_grid, create_grid_object
  implicit none
  
  ! Module Variables                              !    1     2   3   4      5   6    7    8    9
  integer,parameter           :: nproperties = 9  ! Kmin, Kmax, Ss, Sy, Aniso, Kd, KHp, KVp, STp
  
!-------------------------------------------------------------------------------------------------!
  type, extends(t_dataset) :: t_pilotpoints
    ! For reference: inherited from t_dataset
    !integer,allocatable           :: category(:)      ! Zones, Hydrostratigraphic Units (HSUs)
    !type(t_grid),pointer          :: grid             ! Grid
  
    integer,allocatable            :: lay_type(:)      ! Aquifer (1), Aquitard (2), Unused (0)
    integer,pointer                :: pp_id(:)         ! ID value of each pilot point (index matches pilot_points)
    integer,allocatable            :: zonelist(:)      ! List of all "categories" assigned to a pilot point
    real(dp),allocatable           :: values(:,:,:)    ! Data values by parameter (point, parameter, class_id)
    
    ! Counters
    integer                        :: npp              ! Number of pilot points
    integer                        :: nppzones         ! Number of pilot point zones
    integer                        :: max_a_cat        ! Max aquifer category id
    integer                        :: pp_split         ! First aquitard point index
    integer                        :: zonelist_split   ! The first aquitard zone in zonelist
    
    contains
    
      procedure,private            :: get_pp_index
      procedure,public             :: init_pilotpoint
      generic                      :: initialize => init_pilotpoint
      procedure,public             :: read_pp_line
  
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
    
    this%npp            = npoints

    ! Intialize with default values
    this%nppzones       = 1
    this%category       = 1
    this%lay_type       = 0
    this%max_a_cat      = 0
    this%pp_split       = npoints
    this%zonelist_split = 0
    
    this%grid => create_grid_object(2, npoints)  ! X,Y coords by n 
  
  end subroutine init_pilotpoint
  
!-------------------------------------------------------------------------------------------------!
  
  function get_pp_index(this, pp_id) result(index)
    implicit none
    class(t_pilotpoints), intent(inout)   :: this
    integer,intent(in)                    :: pp_id
    integer                               :: index

    do index=1, this%npp
      if (this%pp_id(index)==pp_id) exit
    end do
  
  end function get_pp_index
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine read_pp_line(this, id, pp_type, strings, length, class_names)
    use m_vstringlist, only: t_vstringlist, vstrlist_search, vstrlist_index
    use m_file_io, only: item2int, item2dp
    use m_error_handler, only: error_handler
    
    implicit none
    class(t_pilotpoints),intent(inout) :: this
    character(*),intent(in)            :: id, pp_type
    type(t_vstringlist),intent(in)     :: strings, class_names
    integer,intent(in)                 :: length
    integer                            :: pp_index, class_id
  
    ! Read line for pp_id and get internal pilot point index
    pp_index = this%get_pp_index(item2int(strings, 1))
    ! Now we react to type
    select case(pp_type)
      case("CLASS")
        class_id = vstrlist_search(class_names, vstrlist_index(strings, 2))
        if (pp_index < this%pp_split) then
        this%values(pp_index,1,class_id) = item2dp(strings, 3)  ! Kmin
        this%values(pp_index,2,class_id) = item2dp(strings, 4)  ! Kmax
        this%values(pp_index,3,class_id) = item2dp(strings, 5)  ! Ss
        this%values(pp_index,4,class_id) = item2dp(strings, 6)  ! Sy
        this%values(pp_index,5,class_id) = item2dp(strings, 7)  ! Aniso
        this%values(pp_index,6,class_id) = item2dp(strings, 8)  ! Kd
        else   ! Aquitard
        this%values(pp_index,1,class_id) = item2dp(strings, 3)  ! Kmin
        this%values(pp_index,2,class_id) = item2dp(strings, 4)  ! Kmax
        this%values(pp_index,5,class_id) = item2dp(strings, 5)  ! Aniso
        this%values(pp_index,6,class_id) = item2dp(strings, 6)  ! Kd
        end if
        ! K Processing & Checking
        if (this%values(pp_index,2,class_id) <= this%values(pp_index,1,class_id)) then
          call error_handler(1,"KMax <= Kmin for pilot point " // trim(id))
        else
          ! We krige the *difference* to ensure Kmax > Kmin everywhere
          ! Store this value as a difference
          this%values(pp_index,2,class_id) = this%values(pp_index,2,class_id) - this%values(pp_index,1,class_id)
        end if
      case("GLOBAL")
        ! Set for ALL classes
        this%values(pp_index,7, :) = item2dp(strings, 2)  ! KHp
        this%values(pp_index,8, :) = item2dp(strings, 3)  ! KVp
        this%values(pp_index,9, :) = item2dp(strings, 4)  ! STp
      case DEFAULT
        call error_handler(1,"Unknown Pilot Point Type: " // trim(pp_type))
    end select
  
  end subroutine read_pp_line
  
!-------------------------------------------------------------------------------------------------!  
end module m_pilotpoints