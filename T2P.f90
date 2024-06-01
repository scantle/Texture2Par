module m_T2P
  use m_global, only: dp
  use m_vstringlist, only: t_vstringlist, vstrlist_search, vstrlist_append
  use m_vstring, only: t_vstring, vstring_toupper
  use m_error_handler, only: error_handler,warnings
  use m_options, only: t_options
  use m_datasets, only: t_layerdata, t_dataset
  use m_flow_model, only: t_flow_model, t_iwfm_model
  use m_pilotpoints, only: t_pilotpoints
  use m_interpolator, only: t_interpolator, t_idw
  use m_kriging, only: t_krige, t_simpkrige, t_ordkrige
  use m_categories, only: t_category
  use m_vario, only: t_vgm
!-------------------------------------------------------------------------------------------------!
  type t_T2P
    type(t_options)                  :: opt  ! avoid conflict with Fortran intrinsic 'options'
    class(t_flow_model),pointer      :: fmodel      => null()
    class(t_interpolator),pointer    :: interp_data => null()
    class(t_interpolator),pointer    :: interp_pp   => null()
    type(t_layerdata),pointer        :: datasets(:,:)           ! Datasets (dataset, layer) (1+nsecondary, layer)
    type(t_layerdata),pointer        :: dset_intermediate(:,:)  ! texture or parameters data interpolated to model grid (nsets, nlayers)
    type(t_layerdata),pointer        :: par(:)            ! Final model hydraulic parameters (nlayers)&values(npar, npoints)
    !type(t_vgm),pointer              :: variograms(:,:)   ! (stuctures, variogram)  ! Moved to interpolator
    type(t_pilotpoints),pointer      :: pilot_points
    type(t_pilotpoints),pointer      :: pilot_points_aquitard
    
    ! Counters
    integer                          :: npar        ! Kh, Kv, Ss, Sy
    integer                          :: nclasses    ! Number of data classes
    integer                          :: ntexclasses ! Texture datasets (for which hydraulic properties exist) (one may be inferred)
    integer                          :: nprimary    ! Number of active primary datasets for kriging (texture or parameter, see opt%INTERP_ORDER)
    integer                          :: nsecondary  ! Number of covariate data classes
    integer                          :: ndatasets   ! Number of dataset objects (one for primary datasets, n for secondary)
    integer                          :: nvario      ! Number of variograms, eventually passed to interpolator (if kriging)
    integer                          :: ndatafiles  ! Number of data files to read in
    integer                          :: nhsus       ! Number of Hydrostratigraphic Units (HSUs)
                                                    ! (Formerly "geologic zones")
    
    ! Trackers
    type(t_vstringlist),pointer      :: class_names           ! Data class names
    type(t_vstringlist)              :: parnames              ! Names of output parameters calculated
    integer,pointer                  :: nsecondary_by_id(:)   ! Connections (covariates) per [class, hyd property]
    integer,pointer                  :: connections(:,:)      ! (connections, [classes or hyd properties])
    integer,allocatable              :: class_dataset(:)      ! Translation from class id to dataset_id, ie datasets(dataset_id, layer)
    type(t_vstringlist)              :: dataset_files         ! dataset filenames
    type(t_vstringlist)              :: dataset_trans         ! dataset transformations (same order as filenames)
    type(t_vstringlist)              :: dataset_trans_params  ! dataset transformation parameters (same order as filenames) (unimplemented)
    type(t_category),pointer         :: dataset_ppzones       ! Pilot point zone dataset grid arrays & tracking (when INTERP_ORDER=AFTER)
    
    ! Read trackers - used to check if a block has been read (and make sure it exists)
    logical                          :: has_opt
    logical                          :: has_classes
    logical                          :: has_covar
    logical                          :: has_interp
    logical                          :: has_vario
    logical                          :: has_fmodel
    logical                          :: has_pplocs
    logical                          :: has_pplocs_aquitards
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
    procedure,public   :: setup_covariates
    procedure,public   :: setup_variograms
    procedure,public   :: setup_pplocs
    procedure,public   :: setup_dataset_ppzones
    procedure,public   :: check_major_inputs
    !procedure,private  :: check_pp_zone_layertype_assignments
    procedure,private  :: nvario_estimator
    procedure,private  :: no_covariates
    
    ! Main Procedures
    procedure,public   :: predict
    procedure,public   :: write_output
    procedure,private  :: interpolate_pp
    procedure,private  :: interpolate_model
    procedure,private  :: calc_block_parameters
    
    ! Access procedures
    procedure, public  :: get_class_id
    !procedure, public  :: get_covariate_id
    procedure, public  :: get_parameter_id
    procedure, public  :: get_primary_id
    procedure, public  :: get_active_texture_classes
    procedure, public  :: get_connection_list
    procedure, public  :: get_dataset_idx_list
  
  end type t_T2P
!-------------------------------------------------------------------------------------------------!
  
  contains

!-------------------------------------------------------------------------------------------------!
! SETUP PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine initialize(this)
    use m_vstringlist, only: vstrlist_new, vstrlist_length
    implicit none
    class(t_T2P), intent(inout)      :: this
    
    ! Initialize Read Trackers
    this%has_opt          = .false.
    this%has_classes      = .false.
    this%has_covar        = .false.
    this%has_interp       = .false.
    this%has_vario        = .false.
    this%has_fmodel       = .false.
    this%has_pplocs       = .false.
    this%has_pp           = .false.
    this%has_datafiles    = .false.
    this%has_pplocs_aquitards = .false.
    ! And group trackers...
    this%has_hsus         = .false.
    this%has_ppzones      = .false.
    
    ! Initialize Elements
    call this%opt%initialize()
    call vstrlist_new (this%dataset_files)
    call vstrlist_new (this%dataset_trans)
    
    ! Initialize Aquifer Parameter list - MUST BE IN CAPS
    call vstrlist_new(this%parnames)
    call vstrlist_append (this%parnames, "KH" )
    call vstrlist_append (this%parnames, "KV" )
    call vstrlist_append (this%parnames, "SY" )
    call vstrlist_append (this%parnames, "SS" )
    this%npar = vstrlist_length(this%parnames)
    
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
        this%has_interp = .true.
      case DEFAULT
        call error_handler(1,opt_msg="Invalid Interpolation Method: " // this%opt%INTERP_METHOD)
    end select
    call this%interp_data%initialize()  ! Todo - custom parameters?
    call this%interp_pp%initialize()    ! Todo - custom parameters?
    
    ! Invalid combinations...
    if (this%opt%INTERP_ORDER=="AFTER") then
      if (this%opt%ONLY_TEXTURE) call error_handler(1,opt_msg='ONLY_TEXTURE option is only valid for INTERP_ORDER=="BEFORE"')
      if (this%opt%INFER_LAST_CLASS) call error_handler(1,opt_msg='Cannot use INFER_LAST_CLASS with INTERP_ORDER=="AFTER"')
    end if
  
  end subroutine process_options
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine setup_flow_model(this)
    implicit none
    class(t_T2P), intent(inout)   :: this
    
    if (len_trim(this%fmodel%sim_file)>0.or.len_trim(this%fmodel%grd_file)>0) then
      this%has_fmodel = .true.
    end if
    if (len_trim(this%fmodel%hsu_file))     this%has_hsus = .true.
    if (len_trim(this%fmodel%ppzone_file))  this%has_ppzones = .true.
    
  end subroutine setup_flow_model
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine setup_data_classes(this, nclasses, class_names)
    implicit none
    class(t_T2P), intent(inout)   :: this
    integer,intent(in)            :: nclasses
    type(t_vstringlist),pointer   :: class_names   
    integer                       :: i,n
    
    this%nclasses = nclasses
    this%ntexclasses = nclasses
    this%class_names => class_names
    this%ndatasets = 1
    
    ! Allocate
    allocate(this%nsecondary_by_id(nclasses))
    this%nsecondary_by_id(1:) = 0
    
    call this%nvario_estimator()
    
    this%has_classes = .true.
  end subroutine setup_data_classes
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine setup_covariates(this, nsecondary, connections)
    implicit none
    class(t_T2P), intent(inout)   :: this
    integer,intent(in)            :: nsecondary
    type(t_vstringlist),pointer   :: covnames
    integer,pointer               :: connections(:,:)
    
    this%nsecondary   =  nsecondary
    this%nclasses     = this%nclasses + nsecondary
    this%ndatasets    = 1 + this%nsecondary
    this%connections => connections

    ! Now we can setup the dataset-class connection
    allocate(this%class_dataset(this%nclasses))
    this%class_dataset(1:) = 0
    
    call this%nvario_estimator()
    
    if (maxval(this%nsecondary_by_id) > 0) then
      this%has_covar = .true.
    end if
  end subroutine setup_covariates

!-------------------------------------------------------------------------------------------------!
  
  subroutine setup_dataset_ppzones(this)
    ! Some copying b/c we need PPZONE data at dataset points when opt%INTERP_ORDER=="AFTER"
    implicit none
    class(t_T2P), intent(inout)   :: this
    integer :: nlayers, nnodes, i, j, idx
    
    ! Initialize...
    call this%dataset_ppzones%initialize(this%fmodel%ppzones%narr, this%fmodel%nlayers, this%fmodel%nnodes(1), 0)
    
    ! Copy over all categories
    this%dataset_ppzones%lay2arr = this%fmodel%ppzones%lay2arr
    this%dataset_ppzones%from_file = this%fmodel%ppzones%from_file
    ! Loop over each array and add each category stored in arrcat
    do i = 1, this%fmodel%ppzones%narr
      do j = 1, this%fmodel%ppzones%arrncat(i)
        call this%dataset_ppzones%add_category_int(this%fmodel%ppzones%arrcat(j, i), i, idx)
      end do
    end do
  
  end subroutine setup_dataset_ppzones

!-------------------------------------------------------------------------------------------------!
  
  subroutine nvario_estimator(this)
    ! Need x(x-1)/2 variograms per class where x is classes involved (prime + secondary)
    ! Called during setup_data_classes AND setup_covariates as more info comes in
    ! LATER nvario is passed to the interpolator (if kriging)
    implicit none
    class(t_T2P), intent(inout)   :: this
    integer                       :: n, m, i, j, x, vidx
    
    ! Get problem size
    if (this%opt%INTERP_ORDER=="BEFORE") then
      n = this%ntexclasses
      m = this%nclasses 
    else ! this%opt%INTERP_ORDER=="AFTER")
      n = this%npar
      m = this%npar + this%nsecondary
    end if
    
    ! Defaults
    this%nprimary = n
    this%nvario = m
    
    ! Allocate variogram matrix for data
    ! TODO move to kriging class, call directly using p%
    ! (Requires lots of revising, t2p is setup to manage all of the connections/nvario stuff and only later pass it)
    select type(p => this%interp_data)
      class is(t_krige)
        if (allocated(p%global_vario_idx)) deallocate(p%global_vario_idx)
        allocate(p%global_vario_idx(m,m))
        p%global_vario_idx = 0
    
        do i=1, this%nvario
          p%global_vario_idx(i,i) = i
        end do
    
        ! Subtract one if INFER because we don't krige - last class is the remainder of the other class(es)
        if (this%opt%INFER_LAST_CLASS.and.this%opt%INTERP_ORDER=="BEFORE") then
          p%global_vario_idx(n,n) = 0
          p%global_vario_idx(n+1:m,n+1:m) = p%global_vario_idx(n+1:m,n+1:m) - 1
          this%nprimary = this%nprimary - 1
          this%nvario = this%nvario - 1
        end if
      
        if (ASSOCIATED(this%connections)) then
          vidx = this%nvario
          ! Loop over connections
          do i=1, this%nprimary
            ! Fill matrix
            do j=1, this%nsecondary_by_id(i)
              vidx = vidx + 1
              p%global_vario_idx(i,this%connections(j+1,i)) = vidx
              p%global_vario_idx(this%connections(j+1,i),i) = vidx ! symmetry          
            end do
          end do
          this%nvario = vidx
        end if
      class default
        if (m /= n) call error_handler(4,opt_msg="Covariates are not supported for non-kriging options.")
    end select
    
    ! Add Pilot Point variogram
    this%nvario = this%nvario + 1  
    
    ! Allocate variogram matrix for pp
    select type(p => this%interp_pp)
      class is(t_krige)
        if (allocated(p%global_vario_idx)==.false.) allocate(p%global_vario_idx(1,1))
        p%global_vario_idx(1,1) = 1  ! Is passed only this variogram
    end select
  
  end subroutine nvario_estimator
  
!-------------------------------------------------------------------------------------------------!
  subroutine no_covariates(this)
    implicit none
    class(t_T2P), intent(inout)   :: this
    integer                       :: i
    
    allocate(this%connections(1,this%nprimary))
    this%nsecondary_by_id(:) = 0
    
    ! Setup matrix, values indicate variogram ID. Can start off assuming each class has a variogram (fill diagonal)
    this%connections = 0
    do i = 1, this%nprimary
      this%connections(1, i) = i
    end do
  
  end subroutine no_covariates
!-------------------------------------------------------------------------------------------------!
  
  subroutine setup_variograms(this, variograms, nnear, nstruct)
    implicit none
    class(t_T2P), intent(inout)   :: this
    integer, pointer              :: nnear(:), nstruct(:)
    type(t_vgm),pointer           :: variograms(:,:)
    class(*),pointer              :: p
    integer                       :: nsets
    
    nsets = this%nprimary+this%nsecondary
    
    ! Place directly in the kriging interpolator objects
    select type(p => this%interp_data)
      class is(t_krige)
        ! First n are for classes/properties
        p%nvario     = this%nvario -1
        p%nstruct    => nstruct(1:this%nvario-1)
        p%variograms => variograms(:,1:p%nvario)
        p%nnear      => nnear(1:nsets)
    end select
    ! PP is the last one
    select type(p => this%interp_pp)
      class is(t_krige)
        p%nvario     = 1
        p%nstruct    => nstruct(this%nvario:this%nvario)
        p%variograms => variograms(:,this%nvario:this%nvario)
        p%nnear      => nnear(nsets+1:nsets+1)
    end select
    
    this%has_interp = .true.
    this%has_vario  = .true.
  end subroutine setup_variograms

!-------------------------------------------------------------------------------------------------!
  
  subroutine setup_pplocs(this, pp_id, nzones, zone_tracker)
    implicit none
    class(t_T2P), intent(inout)   :: this
    integer, pointer              :: pp_id(:)
    integer, intent(in)           :: nzones, zone_tracker(:)
    integer                       :: i
    
    ! Called for both aquifer & aquitard pilot points
    
    ! AQUIFER
    if (associated(this%pilot_points)) then
      this%pilot_points%pp_id => pp_id
      allocate(this%pilot_points%zonelist(nzones))
    
      this%pilot_points%nppzones = nzones
      do i=1, nzones
        this%pilot_points%zonelist(i) = zone_tracker(i)
      end do
    
      this%has_pplocs = .true.
      
    ! AQUITARD
    elseif (associated(this%pilot_points_aquitard)) then
      this%pilot_points_aquitard%pp_id => pp_id
      allocate(this%pilot_points_aquitard%zonelist(nzones))
    
      this%pilot_points_aquitard%nppzones = nzones
      do i=1, nzones
        this%pilot_points_aquitard%zonelist(i) = zone_tracker(i)
      end do
    
      this%has_pplocs_aquitards = .true.
    end if
  
  end subroutine setup_pplocs
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine check_major_inputs(this)
    implicit none
    class(t_T2P), intent(inout)   :: this
    
    ! Do we have everything we need to run the model?
    if(this%has_opt      ==.false.) call error_handler(1,opt_msg="Missing Required OPTIONS Block")
    if(this%has_classes  ==.false.) call error_handler(1,opt_msg="Missing Required CLASSES Block")
    if(this%has_fmodel   ==.false.) call error_handler(1,opt_msg="Missing Required FLOW_MODEL Block")
    if(this%has_pplocs   ==.false.) call error_handler(1,opt_msg="Missing Required PP_LOCS Block")
    if(this%has_pp       ==.false.) call error_handler(1,opt_msg="Missing Required PP_PARAMETERS Block")
    if(this%has_datafiles==.false.) call error_handler(1,opt_msg="Missing Required DATASET Block(s)")
    if(this%has_covar    ==.false.) call this%no_covariates()  ! Allowable, just needs to be handled
    if(this%has_interp   ==.false.) then
      if (this%opt%INTERP_METHOD=="OK".or.this%opt%INTERP_METHOD=="SK") then
        call error_handler(1,opt_msg="Missing Required VARIOGRAMS Block")
      else  ! Unknown interpolator parameters required??
        call error_handler(1,opt_msg="Missing Required Interpolation Method Parameters")
      end if
    end if
    
    ! Check pilot point zones
    ! Doesn't ensure they are the same ID values...
    if (this%pilot_points%nppzones /= this%fmodel%get_active_ppzones()) then
      call error_handler(1,opt_msg="Mismatch between Pilot Point Zones assigned to Pilot Points vs Model Grid")
    end if
    
    ! Make sure pilot points and flow model agree on existance of aquitards
    ! TODO - more robust? Nothing checks to ensure values at this points
    if (this%has_pplocs_aquitards /= this%fmodel%has_aquitards) then
      call error_handler(1,opt_msg="Mismatch in aquitard existance between Pilot Points and Model")
    end if
    
    !call this%check_pp_zone_layertype_assignments()
    
  end subroutine check_major_inputs
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! MAIN PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine predict(this)
    use m_grid, only: t_grid
    implicit none
    class(t_T2P), intent(inout)   :: this
    class(t_grid),pointer         :: pp_to                                      ! Where we interpolate the pilot points to
    type(t_category),pointer      :: cat_to                                     ! Categories the pilot points are interpolated using
    type(t_pilotpoints)           :: pp_intermediate(this%fmodel%ppzones%narr)  ! The pilot point data interpolated to pp_to (by ppzone array)
    integer                       :: i,j
    
    ! Allocate the final parameter objects
    allocate(this%par(this%fmodel%nlayers))
    
    if (this%opt%ONLY_TEXTURE) then
      ! If just texture, our jobs is very easy
      allocate(this%dset_intermediate(1,this%fmodel%nlayers))
      call this%interpolate_model(this%datasets, this%dset_intermediate(1,1:))
      RETURN
    end if
    
    ! Assign PPTO_grid (data or fmodel)
    if (this%opt%INTERP_ORDER=="BEFORE") then
      ! T2P 1.0 Method
      allocate(this%dset_intermediate(1,this%fmodel%nlayers))
      ! Uses layer 1 grid
      pp_to => this%fmodel%grid(1)
      cat_to => this%fmodel%ppzones
    else if (this%opt%INTERP_ORDER=="AFTER") then
      allocate(this%dset_intermediate(size(this%datasets,1),this%fmodel%nlayers))
      this%dset_intermediate(2:,1:) => this%datasets(2:,1:)  ! All secondary data gets pointed over
      ! Uses primary dataset layer 1 grid
      ! Again, assumes uniformity of grid by layer
      pp_to => this%datasets(1,1)%grid
      cat_to => this%dataset_ppzones
    end if
  
    ! Interpolate PP to PPTO_grid
    call this%interpolate_pp(pp_to, cat_to, pp_intermediate)
  
    if (this%opt%INTERP_ORDER=="BEFORE") then
      ! Krige textures to fmodel grid
      call this%interpolate_model(this%datasets, this%dset_intermediate(1,1:))
      ! Write weights 
      ! Calc parameters at fmodel grid locations (pp,texture,result)  Kh, Kv, Ss, Sy
      call this%calc_block_parameters(pp_intermediate, this%dset_intermediate(1,1:), this%par)
      
    else if (this%opt%INTERP_ORDER=="AFTER") then
      ! Calc parameters at data locations
      call this%calc_block_parameters(pp_intermediate, this%datasets(1,1:), this%dset_intermediate(1,1:))
      ! Krige parameters to fmodel grid
      call this%interpolate_model(this%dset_intermediate, this%par)
      ! Write weights
    end if
    
  end subroutine predict
!-------------------------------------------------------------------------------------------------!
  
  subroutine interpolate_pp(this, to_grid, to_categories, results)
    use m_grid, only: t_grid
    use m_weights, only: t_weights
    use m_pilotpoints, only: nproperties, naquitard_props, aquitard_props
    ! Pilot Point type in, Pilot Point type out
    implicit none
    class(t_T2P), intent(inout)     :: this
    type(t_pilotpoints)             :: results(this%fmodel%ppzones%narr)
    class(t_grid),pointer           :: to_grid
    type(t_category)                :: to_categories
    type(t_weights)                 :: weights
    integer                         :: i,j,k
    integer                         :: class_idx(1)
    type(t_pilotpoints)             :: obs_data(1)  ! a little trick to appease calc_weights
    
    obs_data(1) = this%pilot_points
    class_idx(1) = 1
    call weights%initialize(this%pilot_points%grid%n, to_grid%n)  !input_n, output_n
    
    ! Go by pilot point zone arrays (up to nlayers)
    do k=1, this%fmodel%ppzones%narr
    
      ! Initialize results for this array
      call results(k)%init_pilotpoint(to_grid%n, to_grid) ! Copies grid over
      allocate(results(k)%values(to_grid%n, size(this%pilot_points%values,2), this%ntexclasses))
    
      ! Check if we have precalculated weights for the pilot point interpolation
      if (this%opt%READ_PP_WEIGHTS) then
        ! TODO READ
      else
        call this%interp_pp%calc_weights(obs_data, &
                                         class_idx,  &
                                         to_grid,  &
                                         weights%values, &
                                         to_categories%arrcat(1:to_categories%arrncat(k),k), &
                                         to_categories%array(:,k),&
                                         this%opt)
      end if
    
      ! Write PP weights
      if (this%opt%WRITE_PP_WEIGHTS) then
        ! TODO WRITE
      end if
    
      ! Apply weights and create results
      do i=1, this%ntexclasses
        do j=1, nproperties
          results(k)%values(:,j,i) = matmul(weights%values,this%pilot_points%values(:,j,i))
        end do
      end do
    end do
    
    ! Aquitards
    if (this%fmodel%has_aquitards) then
      obs_data(1) = this%pilot_points_aquitard
      call weights%initialize(this%pilot_points%grid%n, to_grid%n)  !input_n, output_n
    
      ! Go by pilot point zone arrays (up to nlayers)
      do k=1, this%fmodel%ppzones%narr
    
        ! Check if we have precalculated weights for the pilot point interpolation
        if (this%opt%READ_PP_WEIGHTS) then
          ! TODO READ
        else
          call this%interp_pp%calc_weights(obs_data, &
                                           class_idx,  &
                                           to_grid,  &
                                           weights%values, &
                                           to_categories%arrcat(1:to_categories%arrncat(k),k), &
                                           to_categories%array(:,k),&
                                           this%opt)
        end if
    
        ! Write PP weights
        if (this%opt%WRITE_PP_WEIGHTS) then
          ! TODO WRITE
        end if
    
        ! Apply weights and create results - subset of properties
        do i=1, this%ntexclasses
          do j=1, naquitard_props
            results(k)%values(:,aquitard_props(j),i) = matmul(weights%values,this%pilot_points%values(:,j,i))
          end do
        end do
      end do
    end if
    
  end subroutine interpolate_pp

!-------------------------------------------------------------------------------------------------!
  
  subroutine interpolate_model(this, obs, dout)
    use m_grid, only: t_grid
    use m_weights, only: t_weights
    
    implicit none
    class(t_T2P), intent(inout)   :: this
    class(t_layerdata)            :: obs(:,:)  ! (nprimary, nlayers)
    type(t_layerdata)             :: dout(this%fmodel%nlayers)
    integer                       :: k, i, j, x, wsize
    integer                       :: wstart, wend
    integer,allocatable           :: class_idx(:), dataset_idx(:)
    real(dp)                      :: v_adjust
    type(t_weights)               :: weights
    
    ! Loop over model layers
    do k=1, this%fmodel%nlayers
    
      ! Check if layer is active (generally IWFM: means it has thickness)
      if (this%fmodel%is_active_layer(k)==.false.) cycle
      
      ! Initialize results
      if (this%opt%INFER_LAST_CLASS) then
        call dout(k)%initialize(this%fmodel%grid(1)%n, this%nprimary+1)
      else
        call dout(k)%initialize(this%fmodel%grid(1)%n, this%nprimary)
      end if
      dout(k)%values = 0.0
    
      ! Loop over primary classes
      do i=1, this%nprimary
    
        ! Check if we have precalculated weights for the model grid interpolation
        if (this%opt%READ_GRID_WEIGHTS) then
          ! TODO READ
        else
          ! Interpolate
          class_idx   = this%get_connection_list(i)           ! Classes involved
          dataset_idx = this%get_dataset_idx_list(class_idx)  ! What datasets those classes are in
          wsize = 0
          do j=1, size(dataset_idx)
            wsize = wsize + obs(dataset_idx(j),k)%grid%n
          end do
          call weights%initialize(wsize, this%fmodel%grid(1)%n)  !input_n, output_n
          ! Multiple weights objects?
          !call weights%initialize(obs(1,1)%grid%n, this%fmodel%grid(1)%n)  !input_n, output_n
          !call weights2%initialize(obs(2,1)%grid%n, this%fmodel%grid(1)%n)
          
          x = this%fmodel%hsus%lay2arr(k)  ! Array id for this layer
          call this%interp_data%calc_weights(obs(dataset_idx,k), &
                                            class_idx,  &
                                            this%fmodel%grid(1),  &
                                            weights%values, &
                                            this%fmodel%hsus%arrcat(1:this%fmodel%hsus%arrncat(x),x), &
                                            this%fmodel%hsus%array(1:,x),&
                                            this%opt)
        end if
    
        ! Write weights
        if (this%opt%WRITE_GRID_WEIGHTS) then
          ! TODO WRITE
        end if
        ! Apply weights and create results
        wsize = 0
        do j=1, size(dataset_idx)
          wstart = wsize+1
          wend = wsize+obs(dataset_idx(j),k)%grid%n
          v_adjust = 0.0_dp
          if (j>1) v_adjust = obs(1,k)%mean(i) - obs(dataset_idx(j),k)%mean(1)  ! secondary data standardization factor
          ! Weights%values are (output_points, input_points)
          !                     ngrid          nobs (all datasets)
          dout(k)%values(i,1:) =  dout(k)%values(i,1:) + matmul(weights%values(:,wstart:wend),obs(dataset_idx(j),k)%values(i,1:)+v_adjust)
          wsize = wsize + obs(dataset_idx(j),k)%grid%n
        end do
        ! Enforce values must be greater than zero
        where(dout(k)%values(i,1:) < 0.0) dout(k)%values(i,1:) = 0.0_dp
        ! If we're kriging texture - must be below 1.0
        if (this%opt%INTERP_ORDER=="BEFORE") where(dout(k)%values(i,1:) > 1.0) dout(k)%values(i,1:) = 1.0_dp
      end do
      ! If we need to, infer the last class
      if (this%opt%INFER_LAST_CLASS) dout(k)%values(this%nprimary+1,1:) = 1.0_dp - sum(dout(k)%values(1:this%nprimary,1:), dim=1)
    end do
    
  end subroutine interpolate_model
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine calc_block_parameters(this, pp_obj, tex_obs, par_out)
    implicit none
    class(t_T2P), intent(inout)    :: this
    type(t_pilotpoints),intent(in) :: pp_obj(this%fmodel%ppzones%narr) ! Kmin, Kmax, Ss, Sy, Aniso, Kd, KHp, KVp, STp
    type(t_layerdata),intent(in)   :: tex_obs(this%fmodel%nlayers)  !  values(ntexclasses, npoints)
    type(t_layerdata),intent(out)  :: par_out(this%fmodel%nlayers)  !values([Kh, Kv, Ss, Sy], npoints)
    integer                        :: i, k, r
    real(dp),allocatable           :: depths(:),KhD(:,:),KvD(:,:)
    
    write(*,*) 'here'
    
    ! Loop over layers
    do k=1, this%fmodel%nlayers
    
      ! Check if layer is active (generally IWFM: means it has thickness)
      if (this%fmodel%is_active_layer(k)==.false.) cycle
    
      ! Initialize Par
      call par_out(k)%initialize(this%fmodel%nnodes(k), this%npar)
      par_out(k)%grid => this%fmodel%grid(k)
      
      ! Allocate arrays
      allocate(KhD(this%fmodel%nnodes(k), this%ntexclasses))
      allocate(KvD(this%fmodel%nnodes(k), this%ntexclasses))
    
      ! Get node midpoint depths
      depths = this%fmodel%get_midpoint_depths(k)
      
      ! Get pilot point array  (point, parameter, class_id)
      r = this%fmodel%ppzones%lay2arr(k)
      
      do i=1, this%fmodel%nnodes(k)
        ! Kh - Depth adjusted
        KhD(i,:) = pp_obj(r)%values(i, 1, :) + pp_obj(r)%values(i, 2, :) * exp(-1_dp *pp_obj(r)%values(1, 6, :) * depths(i))
        ! Kv - Depth adjusted, calculated from anisotropy
        KvD(i,:) = KhD(i,:) / pp_obj(r)%values(i, 5, :)
        ! Kh
        par_out(k)%values(1,i) = DOT_PRODUCT(KhD(i,:)**pp_obj(r)%values(i, 7, :), tex_obs(k)%values(:,i)) ** (1.0_dp / pp_obj(r)%values(i, 7, 1))
        ! Kv
        par_out(k)%values(2,i) = DOT_PRODUCT(KvD(i,:)**pp_obj(r)%values(i, 8, :), tex_obs(k)%values(:,i)) ** (1.0_dp / pp_obj(r)%values(i, 8, 1))
        ! Ss
        par_out(k)%values(3,i) = DOT_PRODUCT(pp_obj(r)%values(i,3,:)**pp_obj(r)%values(i, 9, :), tex_obs(k)%values(:,i)) ** (1.0_dp / pp_obj(r)%values(i, 9, 1))
        ! Sy
        par_out(k)%values(4,i) = DOT_PRODUCT(pp_obj(r)%values(i,4,:)**pp_obj(r)%values(i, 9, :), tex_obs(k)%values(:,i)) ** (1.0_dp / pp_obj(r)%values(i, 9, 1))
      end do
      
      deallocate(KhD, KvD)
    
    end do
  
  end subroutine calc_block_parameters
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine write_output(this)
    implicit none
    class(t_T2P), intent(inout)   :: this
    type(t_layerdata),pointer     :: p(:)
    
    ! Write parameters by node/cell, and if texture was kriged, textures by node/cell
    if (this%opt%WRITE_NODE_FILES) call this%fmodel%write_node_files(this%par, this%parnames)
    if (this%opt%WRITE_NODE_FILES.and.this%opt%INTERP_ORDER=="BEFORE") then
      p=>this%dset_intermediate(1,:)
      call this%fmodel%write_node_files(p, this%class_names)
    end if
    
    if (this%opt%WRITE_DATASET_FILES) write(*,*) 'todo'
    
    if (this%opt%WRITE_NONE) return  ! Skips model output files
    
    ! Write out flow model input files
    call this%fmodel%write_fmodel_input(this%par)
  
  end subroutine write_output
  
!-------------------------------------------------------------------------------------------------!
! ACCESS PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  function get_class_id(this, string, adjust_for_infer) result (id)
    ! Returns zero if not found
    implicit none
    class(t_T2P), intent(inout)   :: this
    type(t_vstring),intent(in)    :: string
    integer                       :: id
    logical,optional,intent(in)   :: adjust_for_infer
    id = vstrlist_search(this%class_names, vstring_toupper(string), exact=.true.)
    if (present(adjust_for_infer)) then
      ! If you're inferring the last class, it gets off past one after ntexclasses
      if (adjust_for_infer.and.this%opt%INFER_LAST_CLASS.and.(id>this%ntexclasses)) id = id -1
    end if
    return
  end function get_class_id

!-------------------------------------------------------------------------------------------------!
  
  function get_parameter_id(this, string) result (id)
    ! Returns zero if not found
    implicit none
    class(t_T2P), intent(inout)   :: this
    type(t_vstring),intent(in)    :: string
    integer                       :: id
    id = vstrlist_search(this%parnames, vstring_toupper(string), exact=.true.)
    return
  end function get_parameter_id

!-------------------------------------------------------------------------------------------------!
  
  function get_primary_id(this, string) result(prime_id)
    ! Returns zero if not found
    implicit none
    class(t_T2P), intent(inout)   :: this
    type(t_vstring),intent(in)    :: string
    integer                       :: prime_id
    
    if (this%opt%INTERP_ORDER=="BEFORE") then
      prime_id = this%get_class_id(string, adjust_for_infer=.true.)
    else ! this%opt%INTERP_ORDER=="AFTER")
      prime_id = this%get_parameter_id(string)
    end if
  
  end function get_primary_id
  
!-------------------------------------------------------------------------------------------------!
  
  function get_connection_list(this,prime) result(connlist)
    implicit none
    class(t_T2P)            :: this
    integer,intent(in)      :: prime
    integer,allocatable     :: connlist(:)
    
    if (allocated(connlist)) deallocate(connlist)
    allocate(connlist(this%nsecondary_by_id(prime)))
    connlist = this%connections(1:1+this%nsecondary_by_id(prime),prime) 
    return
  end function

!-------------------------------------------------------------------------------------------------!
  
  function get_dataset_idx_list(this,class_ids) result(dataset_idx)
    implicit none
    class(t_T2P)            :: this
    integer,intent(in)      :: class_ids(:)
    integer,allocatable     :: dataset_idx(:)
    
    if (allocated(dataset_idx)) deallocate(dataset_idx)
    allocate(dataset_idx(size(class_ids)))
    dataset_idx(:) = this%class_dataset(class_ids)
    return
  end function
  
!-------------------------------------------------------------------------------------------------!
  
  function get_active_texture_classes(this) result(n)
    implicit none
    class(t_T2P)            :: this
    integer                 :: n
    n = this%ntexclasses
    if (this%opt%INFER_LAST_CLASS) n = n-1
    return
  end function get_active_texture_classes

  !-------------------------------------------------------------------------------------------------!
  
end module m_T2P