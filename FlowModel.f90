module m_flow_model
  use m_file_io, only: t_file_reader, t_file_writer
  use m_grid, only: t_grid
  use m_error_handler, only: error_handler
  implicit none
!-------------------------------------------------------------------------------------------------!
  type :: t_flow_model  ! Abstract Base Class
  ! Base class contains some derived type variables (e.g., the IWFM preproc file) since otherwise
  ! the compiler gets confused when you try to set them on the polymorphic t2p%fmodel type
  ! The alternative was to write a lot of get/set-type methods that would have to be overwritten
  ! by each derived class.
  
    integer                       :: model_type    ! 1 - grid, 2 - IWFM, 3 - MODFLOW 2000
    integer                       :: nnodes        ! Number of nodes (points at which flow is solved)
    integer                       :: nlayers       ! Number of model layers
    character(100)                :: preproc_file  ! IWFM main discretization file (empty for others)
    character(100)                :: sim_file      ! Simulation File or MODFLOW name file
    character(100)                :: ppzone_file   ! File of pilot point zones defined for each node
    character(100)                :: hsu_file      ! File of Hydrostratigraphic Units (HSUs) defined for each node
    character(100)                :: template_file ! File copied for output_file
    character(100)                :: output_file   ! File where aquifer parameters are written
    type(t_Grid), pointer         :: layer(:)      ! Grid for ground surface + each layer bottom
  
    contains
      procedure,public            :: initialize
      !procedure,private     :: write_model_input
      !procedure,private     ::
  
  end type t_flow_model
!-------------------------------------------------------------------------------------------------!
  type, extends(t_flow_model) :: t_mgrid
  
  contains
    procedure,public          :: initialize => init_mgrid
  
  end type 
!-------------------------------------------------------------------------------------------------!
  type, extends(t_flow_model) :: t_iwfm_model   ! I DO realize the "m" in IWFM is for model
  
  contains
    procedure,public          :: initialize => init_iwfm
  
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
    this%preproc_file  = ''
    this%sim_file      = ''
    this%ppzone_file   = ''
    this%hsu_file      = ''
    this%template_file = ''
    this%output_file   = ''
    
  end subroutine initialize
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! GRID CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

  subroutine init_mgrid(this)
    implicit none
    class(t_mgrid)       :: this
    
    call initialize(this)
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
    this%model_type  = 2
    
  end subroutine init_iwfm  
  
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
    
  end subroutine init_modflow2000
  
!-------------------------------------------------------------------------------------------------!
end module m_flow_model