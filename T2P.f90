module m_T2P
  use m_vstringlist, only: t_vstringlist
  use m_options, only: t_options
  use m_flow_model, only: t_flow_model
  use m_vario, only: t_vgm
!-------------------------------------------------------------------------------------------------!
  type t_T2P
    type(t_options)                   :: opt  ! avoid conflict with Fortran intrinsic 'options'
    class(t_flow_model),pointer       :: fmodel => null()
    !type(t_Interpolator),pointer      :: interp_data
    !type(t_Interpolator),pointer      :: interp_pp
    !type(t_Dataset),pointer           :: datasets(:)
    type(t_vgm),pointer               :: variograms(:,:)  ! (stuctures, variogram)
    !type(t_PilotPoint),pointer        :: pilot_points
    
    ! Counters
    integer                          :: nclasses   ! Number of data classes
    integer                          :: nprimary   ! Number of *primary* datasets (for which hydraulic properties exist)
    integer                          :: nppzones   ! Number of Pilot Point Zones
    integer                          :: nhsus      ! Number of Hydrostratigraphic Units (HSUs)
                                                   ! (Formerly "geologic zones")
    integer                          :: nvario     ! Number of variograms
    integer,allocatable              :: nstruct(:) ! Number of structures for each variogram
    
    ! Trackers
    type(t_vstringlist),pointer      :: class_names
    integer,pointer                  :: class_connections(:,:)
    integer,pointer                  :: class_type(:)    ! Primary (1) or Secondary (2) for each class
    type(t_vstringlist)              :: dataset_files    ! dataset filenames
    type(t_vstringlist)              :: dataset_trans    ! dataset transformations
  
  contains
  
    procedure :: initialize
    procedure :: setup_data_classes
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
    
    ! Initialize Elements
    call this%opt%initialize()
    call vstrlist_new (this%dataset_files)
    call vstrlist_new (this%dataset_trans)
    
    ! Intialize counters (useful to determine when a block is missing or has not been read yet)
    this%nclasses = 0
    this%nvario   = 0
    
    ! Set class defaults
    this%nppzones = 1
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
    
    this%nprimary = num_primary
    this%nclasses = this%opt%DATA_CLASSES
    this%class_names => class_names
    this%class_connections => class_connections
    this%class_type => class_type
  
  end subroutine setup_data_classes

!-------------------------------------------------------------------------------------------------!
  
  subroutine setup_variograms(this)
    implicit none
    class(t_T2P), intent(inout)   :: this
  
  end subroutine setup_variograms
!-------------------------------------------------------------------------------------------------!
  
end module m_T2P