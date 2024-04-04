module m_T2P
  use m_vstringlist, only: t_vstringlist, vstrlist_search
  use m_vstring, only: t_vstring
  use m_error_handler, only: error_handler,warnings
  use m_options, only: t_options
  use m_datasets, only: t_layerdata, t_dataset
  use m_flow_model, only: t_flow_model
  use m_pilotpoints, only: t_pilotpoints
  use m_interpolator, only: t_interpolator, t_idw
  use m_kriging, only: t_simpkrige, t_ordkrige
  use m_vario, only: t_vgm
!-------------------------------------------------------------------------------------------------!
  type t_T2P
    type(t_options)                  :: opt  ! avoid conflict with Fortran intrinsic 'options'
    class(t_flow_model),pointer      :: fmodel      => null()
    class(t_interpolator),pointer    :: interp_data => null()
    class(t_interpolator),pointer    :: interp_pp   => null()
    type(t_layerdata),pointer        :: datasets(:,:)     ! Datasets (dataset, layer)
    type(t_layerdata)                :: hydraulic_par(10) ! Hydraulic Parameters - KCMin, DeltaKC, KFMin, DeltaKF, SsC, SsF, SyC, Aniso
    type(t_vgm),pointer              :: variograms(:,:)   ! (stuctures, variogram)
    type(t_pilotpoints),pointer      :: pilot_points
    
    ! Counters
    integer                          :: nclasses   ! Number of data classes
    integer                          :: ndatasets  ! Number of data sets (one may be inferred)
    integer                          :: ndatafiles ! Number of data files to read in
    integer                          :: nprimary   ! Number of *primary* datasets (for which hydraulic properties exist)
    integer                          :: nhsus      ! Number of Hydrostratigraphic Units (HSUs)
                                                   ! (Formerly "geologic zones")
    integer                          :: nvario     ! Number of variograms
    integer,pointer                  :: nstruct(:) ! Number of structures for each variogram
    
    ! Trackers
    type(t_vstringlist),pointer      :: class_names
    integer,pointer                  :: class_connections(:,:)
    integer,pointer                  :: class_type(:)    ! Primary (1) or Secondary (2) for each class
    integer,allocatable              :: class_dataset(:) ! Which dataset object each class is stored in
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
    
    ! Groups (HSUs, PP Zones)
    logical                          :: has_hsus
    logical                          :: has_ppzones
  
  contains
    ! Setup procedures
    procedure,public   :: initialize
    procedure,public   :: process_options
    procedure,public   :: setup_flow_model
    procedure,public   :: setup_data_classes
    procedure,public   :: setup_variograms
    procedure,public   :: setup_pplocs
    !procedure :: allocate_data_classes
    !procedure :: allocate
    
    ! Main Procedures
    procedure,public   :: predict
    
    ! Access procedures
    procedure, public  :: get_class_id
  
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
    ! And group trackers...
    this%has_hsus      = .false.
    this%has_ppzones   = .false.
    
    ! Initialize Elements
    call this%opt%initialize()
    call vstrlist_new (this%dataset_files)
    call vstrlist_new (this%dataset_trans)
    
    ! Intialize counters (useful to determine when a block is missing or has not been read yet)
    this%nclasses = 0
    this%nvario   = 0
    
    ! Set class defaults
    this%nhsus    = 1
    
    ! Setup warnings
    call vstrlist_new (warnings)  ! belongs to m_error_handler
    
  end subroutine initialize

!-------------------------------------------------------------------------------------------------!
  
  subroutine process_options(this)
    ! Hook for after options are read
    implicit none
    class(t_T2P), intent(inout)   :: this
    
    ! Handle Interpolation Method
    select case (this%opt%INTERP_METHOD)
      case("OK")
        allocate(t_ordkrige::this%interp_data)
        allocate(t_ordkrige::this%interp_pp)
      case("SK")
        allocate(t_simpkrige::this%interp_data)
        allocate(t_simpkrige::this%interp_pp)
      case("IDW")
        allocate(t_idw::this%interp_data)
        allocate(t_idw::this%interp_pp)
      case DEFAULT
        call error_handler(1,opt_msg="Invalid Interpolation Method: " // this%opt%INTERP_METHOD)
    end select
    call this%interp_data%initialize()  ! Todo - custom parameters?
    call this%interp_pp%initialize()    ! Todo - custom parameters?
  
  end subroutine process_options
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine setup_flow_model(this)
    implicit none
    class(t_T2P), intent(inout)   :: this
    
    if (len_trim(this%fmodel%sim_file).or.len_trim(this%fmodel%grd_file)) then
      this%has_fmodel = .true.
    end if
    if (len_trim(this%fmodel%hsu_file))     this%has_hsus = .true.
    if (len_trim(this%fmodel%ppzone_file))  this%has_ppzones = .true.
    
  end subroutine setup_flow_model
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine setup_data_classes(this, num_primary, num_total, class_names, class_connections, class_type)
    implicit none
    class(t_T2P), intent(inout)   :: this
    integer,intent(in)            :: num_primary, num_total
    type(t_vstringlist),pointer   :: class_names
    integer,pointer               :: class_connections(:,:), class_type(:)
    integer                       :: i
    
    ! TODO needs to store different when INTERP_ORDER: AFTER is passed
    
    this%nprimary = num_primary
    this%nclasses = num_total
    this%ndatasets = num_total
    this%class_names => class_names
    this%class_connections => class_connections
    this%class_type => class_type
    
    ! TODO Maybe this belongs somewhere else, but this is when find out
    ! TODO handle INTERP_ORDER: AFTER
    this%nvario = maxval(this%class_connections) + 1  ! +1 for pilot points
    
    ! Subtract one if INFER because we don't krige - last class is the remainder of the other class(es)
    ! We also don't store it as a dataset
    if (this%opt%INFER_LAST_CLASS) then
      this%nvario = this%nvario - 1
      this%ndatasets = this%ndatasets - 1
    end if
    
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
  
  function get_class_id(this, string) result (id)
    ! Returns zero if not found
    implicit none
    class(t_T2P), intent(inout)   :: this
    type(t_vstring),intent(in)    :: string
    integer                       :: id
    id = vstrlist_search(this%class_names, string)
  end function get_class_id
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine predict(this)
    implicit none
    class(t_T2P), intent(inout)   :: this
    ! Assign PPTO_grid (data or fmodel)
    ! Krige PP to PPTO_grid
    ! Write PP weights
  
    ! INTERP_ORDER: BEFORE
    ! Krige textures to fmodel grid
    ! Write weights 
    ! Calc parameters at fmodel grid locations
  
    ! INTERP_ORDER: AFTER
    ! Calc parameters at data locations
    ! Krige parameters to fmodel grid
    ! Write weights
    
  end subroutine predict
!-------------------------------------------------------------------------------------------------!

end module m_T2P