module m_T2P
  use m_vstringlist, only: t_vstringlist
  use m_options, only: t_options
  use m_flow_model, only: t_flow_model
  use m_pilotpoints, only: t_pilotpoints
  use m_vario, only: t_vgm
!-------------------------------------------------------------------------------------------------!
  type t_T2P
    type(t_options)                  :: opt  ! avoid conflict with Fortran intrinsic 'options'
    class(t_flow_model),pointer      :: fmodel => null()
    !type(t_Interpolator),pointer     :: interp_data
    !type(t_Interpolator),pointer     :: interp_pp
    !type(t_Dataset),pointer          :: datasets(:)
    type(t_vgm),pointer              :: variograms(:,:)  ! (stuctures, variogram)
    type(t_pilotpoints),pointer      :: pilot_points
    
    ! Counters
    integer                          :: nclasses   ! Number of data classes
    integer                          :: ndatafiles
    integer                          :: nprimary   ! Number of *primary* datasets (for which hydraulic properties exist)
    integer                          :: nhsus      ! Number of Hydrostratigraphic Units (HSUs)
                                                   ! (Formerly "geologic zones")
    integer                          :: nvario     ! Number of variograms
    integer,pointer                  :: nstruct(:) ! Number of structures for each variogram
    
    ! Trackers
    type(t_vstringlist),pointer      :: class_names
    integer,pointer                  :: class_connections(:,:)
    integer,pointer                  :: class_type(:)    ! Primary (1) or Secondary (2) for each class
    type(t_vstringlist)              :: dataset_files    ! dataset filenames
    type(t_vstringlist)              :: dataset_trans    ! dataset transformations
    integer,pointer                  :: nnear(:)         ! Number of nearest points used for interpolation
    
    ! Read trackers - used to check if a block has been read (and make sure it exists)
    logical                          :: has_opt
    logical                          :: has_classes
    logical                          :: has_interp
    logical                          :: has_fmodel
    logical                          :: has_pplocs
    logical                          :: has_pp
    logical                          :: has_datafiles
  
  contains
  
    procedure,public :: initialize
    procedure,public :: setup_data_classes
    procedure,public :: setup_variograms
    procedure,public :: setup_pplocs
    !procedure :: allocate_data_classes
    !procedure :: allocate
  
  end type t_T2P
!-------------------------------------------------------------------------------------------------!
  
  contains
  
  subroutine initialize(this)
    use m_file_io, only: open_file_writer, log_unit
    use m_vstringlist, only: vstrlist_new
    implicit none
    class(t_T2P), intent(inout)      :: this
    
    ! Initialize Read Trackers
    this%has_opt       = .false.
    this%has_classes   = .false.
    this%has_interp    = .false.
    this%has_fmodel    = .false.
    this%has_pplocs    = .false.
    this%has_pp        = .false.
    this%has_datafiles = .false.
    
    ! Initialize Elements
    call this%opt%initialize()
    call vstrlist_new (this%dataset_files)
    call vstrlist_new (this%dataset_trans)
    
    ! Intialize counters (useful to determine when a block is missing or has not been read yet)
    this%nclasses = 0
    this%nvario   = 0
    
    ! Set class defaults
    this%nhsus    = 1
    
  end subroutine initialize
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine setup_data_classes(this, num_primary, class_names, class_connections, class_type)
    implicit none
    class(t_T2P), intent(inout)   :: this
    integer,intent(in)            :: num_primary
    type(t_vstringlist),pointer   :: class_names
    integer,pointer               :: class_connections(:,:), class_type(:)
    integer                       :: i,j
    
    ! TODO needs to store different when INTERP_ORDER: AFTER is passed
    
    this%nprimary = num_primary
    this%nclasses = this%opt%DATA_CLASSES
    this%class_names => class_names
    this%class_connections => class_connections
    this%class_type => class_type
    
    ! TODO Maybe this belongs somewhere else, but this is when find out
    ! TODO handle INTERP_ORDER: AFTER
    this%nvario = maxval(this%class_connections) + 1  ! +1 for pilot points
    
    ! Subtract one if INFER because we don't krige - it's the remainder of the other class(es)
    if (this%opt%INFER_LAST_CLASS) this%nvario = this%nvario - 1  
    
    this%has_classes = .true.
  
  end subroutine setup_data_classes

!-------------------------------------------------------------------------------------------------!
  
  subroutine setup_variograms(this, variograms, nnear, nstruct)
    implicit none
    class(t_T2P), intent(inout)   :: this
    integer, pointer              :: nnear(:), nstruct(:)
    type(t_vgm),pointer           :: variograms(:,:)
    
    this%variograms => variograms
    this%nnear      => nnear
    this%nstruct    => nstruct
    
    this%has_interp = .true.
  
  end subroutine setup_variograms
!-------------------------------------------------------------------------------------------------!
  
  subroutine setup_pplocs(this, pp_id, nzones, zone_tracker)
    implicit none
    class(t_T2P), intent(inout)   :: this
    integer, pointer              :: pp_id(:)
    integer, intent(in)           :: nzones, zone_tracker(:)
    integer                       :: i
    
    this%pilot_points%pp_id => pp_id
    allocate(this%pilot_points%zonelist(nzones))
    
    do i=1, nzones
      this%pilot_points%zonelist(i) = zone_tracker(i)
    end do
    
    this%has_pplocs = .true.
  
  end subroutine setup_pplocs
!-------------------------------------------------------------------------------------------------!

end module m_T2P