module m_flow_model
  use m_file_io
  use m_grid, only: t_grid, create_grid_object
  use m_error_handler, only: error_handler
  implicit none
!-------------------------------------------------------------------------------------------------!
  ! Module variables
  type(t_file_reader), pointer :: reader          ! For reading input files
!-------------------------------------------------------------------------------------------------!
  
  type :: t_flow_model  ! Abstract Base Class
  ! Base class contains some derived type variables (e.g., the IWFM preproc file) since otherwise
  ! the compiler gets confused when you try to set them on the polymorphic t2p%fmodel type
  ! The alternative was to write a lot of get/set-type methods that would have to be overwritten
  ! by each derived class.
  
    integer                       :: model_type    ! 1 - grid, 2 - IWFM, 3 - MODFLOW 2000
    integer                       :: nlayers       ! Number of model layers
    integer,allocatable           :: nnodes(:)     ! Number of nodes (points at which flow is solved) by layer
    character(100)                :: preproc_file  ! IWFM main discretization file (empty for others)
    character(100)                :: sim_file      ! Simulation File or MODFLOW name file
    character(100)                :: ppzone_file   ! File of pilot point zones defined for each node
    character(100)                :: hsu_file      ! File of Hydrostratigraphic Units (HSUs) defined for each node
    character(100)                :: template_file ! File copied for output_file
    character(100)                :: output_file   ! File where aquifer parameters are written
    type(t_Grid), pointer         :: layer(:)      ! Grid for ground surface (0) + each layer bottom
    logical                       :: laysamexy     ! Flag for when X,Y coords can be assumed same for each layer
  
    contains
      procedure,public            :: initialize
      procedure,public            :: read_model
      procedure,private           :: alloc_layers
      procedure,private           :: init_layers
      !procedure,private     :: write_model_input
  
  end type t_flow_model
!-------------------------------------------------------------------------------------------------!
  type, extends(t_flow_model) :: t_mgrid
  
  contains
    procedure,public          :: initialize => init_mgrid
  
  end type 
!-------------------------------------------------------------------------------------------------!
  type, extends(t_flow_model) :: t_iwfm_model   ! I DO realize the "m" in IWFM is for model
  
    character(100)            :: elem_file
    character(100)            :: node_file
    character(100)            :: strt_file
    logical,allocatable       :: lay_exist(:)   ! Tracks whether layer has thickness > 0.00
  
  contains
    procedure,public          :: initialize => init_iwfm
    procedure,public          :: read_model => read_iwfm
    procedure,private         :: alloc_layers => alloc_layers_iwfm
  
  end type 
!-------------------------------------------------------------------------------------------------!
  type, extends(t_flow_model) :: t_modflow_model
  
  logical                     :: structured    ! T/F Structured (row/col) model
  
  contains
    procedure,public          :: initialize => init_modflow2000
  
  end type 
!-------------------------------------------------------------------------------------------------!  
  
  contains
  
!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  function create_flow_model_object(model_type_char) result(object)
    implicit none
    
    character(*),intent(in)        :: model_type_char
    class(t_flow_model),pointer    :: object
    
    select case(model_type_char)
      case ("GRID")
        allocate(t_mgrid::object)
      case ("IWFM")
        allocate(t_iwfm_model::object)
      case ("MODFLOW")
        allocate(t_modflow_model::object)
      case DEFAULT
        call error_handler(1,opt_msg="Unknown Flow Model Type Name: " // trim(model_type_char))
    end select
    
    ! Initialize new object
    call object%initialize()
  
  end function create_flow_model_object
  
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine check_id_validity(id, fmodel, valid_type)
    implicit none
    ! Checks (A) if flow model type has been set and (B) optionally if it is the correct type of
    ! input for the model type (valid_type)
    ! Used while reading in options
    class(t_flow_model),pointer  :: fmodel
    character(30)                :: id
    integer,optional             :: valid_type
    
    !select type (fmodel)
    !  type is (t_iwfm_model)
    !    write(*,*) 'IWFM MODEL'
    !  type is (t_flow_model)
    !    write(*,*) 'GENERIC!'
    !end select
    
    ! Is type set?
    if (ASSOCIATED(fmodel)==.false.) then
      call error_handler(2,opt_msg="Must specify flow model type before other settings")
    else if (present(valid_type)) then
      ! Is it the right type for this setting?
      if (valid_type /= fmodel%model_type) then
        call error_handler(2,opt_msg="Wrong flow model type for option: " // id)
      end if
    end if
    
  end subroutine check_id_validity
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! BASE CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine initialize(this)
    ! Should not be run
    implicit none
    class(t_flow_model)       :: this
    
    this%model_type    = 0
    this%laysamexy     = .false.
    this%preproc_file  = ''
    this%sim_file      = ''
    this%ppzone_file   = ''
    this%hsu_file      = ''
    this%template_file = ''
    this%output_file   = ''
    
  end subroutine initialize
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine read_model(this)
    ! Base class method for reading in [nonexistant] model input file(s)
    ! Should not be run
    implicit none
    class(t_flow_model)       :: this
    
  end subroutine read_model
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine alloc_layers(this)
    implicit none
    class(t_flow_model)       :: this
    
    allocate(this%nnodes(this%nlayers), &
             this%layer (0:this%nlayers))  ! +1 for GSE (layer(0))
    
  end subroutine alloc_layers
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine init_layers(this)
    implicit none
    class(t_flow_model)       :: this
    integer                   :: i
    
    do i=1, this%nlayers
      this%layer(i) = create_grid_object(3, this%nnodes(i))
    end do
    
  end subroutine init_layers
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! GRID CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

  subroutine init_mgrid(this)
    implicit none
    class(t_mgrid)       :: this
    
    call initialize(this)
    this%laysamexy   = .false.  ! Todo
    this%model_type  = 1
    
  end subroutine init_mgrid  
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! IWFM CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

  subroutine init_iwfm(this)
    implicit none
    class(t_iwfm_model)       :: this
    
    call initialize(this)
    this%laysamexy   = .true.
    this%model_type  = 2
    
  end subroutine init_iwfm  

!-------------------------------------------------------------------------------------------------!
  
  subroutine alloc_layers_iwfm(this)
    implicit none
    class(t_iwfm_model)       :: this
    
    allocate(this%lay_exist(this%nlayers))
    
    ! Call super
    call alloc_layers(this)
    
  end subroutine alloc_layers_iwfm
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine read_iwfm(this)
    use m_global, only: dp
    use fpath
    use m_vstringlist, only: t_vstringlist
    implicit none
    type(t_file_reader), pointer :: reader2
    class(t_iwfm_model)          :: this
    character(200)               :: path
    real(dp)                     :: factor
    integer                      :: itemp
    
    ! Standard file reader variables
    integer                    :: status, length, eof
    character(200)             :: value
    type(t_vstringlist)        :: strings
    
    ! Get relative path to later be used
    call fpath_strip(this%preproc_file, path, value)
    
    ! Preprocessor file
    reader => open_file_reader(this%preproc_file)
    call reader%read_to_next_iwfm_line(eof)
    call reader%skip(3, eof)                  ! Title
    call reader%read_to_next_iwfm_line(eof)
    call reader%skip(1, eof)                  ! Binary output
    call reader%next_item(eof, strings)       ! Element File
    call item2char(strings, 1, value)
    call fpath_join(path, value, this%elem_file)
    call reader%next_item(eof, strings)       ! Node File
    call item2char(strings, 1, value)
    call fpath_join(path, value, this%node_file)
    call reader%next_item(eof, strings)       ! Stratigraphy File
    call item2char(strings, 1, value)
    call fpath_join(path, value, this%strt_file)
    call reader%close_file()
    
    ! Stratigraphy file (only to nlayers)
    reader2 => open_file_reader(this%strt_file)
    call reader%read_to_next_iwfm_line(eof)
    call reader%next_item(eof, strings)       ! nlayers
    this%nlayers = item2int(strings, 1) * 2   ! Each layer has overlying aquitard (to be treated as layer)
    
    ! Allocate layers
    call this%alloc_layers()
    
    ! Node file
    reader => open_file_reader(this%node_file)
    call reader%read_to_next_iwfm_line(eof)
    call reader%next_item(eof, strings)       ! nnodes
    this%nnodes(:) = item2int(strings, 1)     ! Node count does not vary by layer in IWFM
    call reader%next_item(eof, strings)       ! xy mult factor
    factor = item2dp(strings, 1)
    
    ! Allocate nodes
    call this%init_layers()
    
    ! Continue Node file read
    
    
  end subroutine read_iwfm
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! MODFLOW 2000 CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine init_modflow2000(this)
    implicit none
    class(t_modflow_model)       :: this
    
    call initialize(this)
    this%model_type  = 3
    this%structured  = .true.
    this%laysamexy   = .true.
    
  end subroutine init_modflow2000
  
!-------------------------------------------------------------------------------------------------!
end module m_flow_model