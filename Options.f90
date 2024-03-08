module m_options
!-------------------------------------------------------------------------------------------------!  
! HOW TO ADD AN OPTION
!  1. Add it as a variable in the type t_options
!  2. Add it to the initialize routine and give it a default value
!  3. Create a new case for it in the read_options routine. Write code to handle assignment.
!-------------------------------------------------------------------------------------------------!
  type t_options
    ! Required Options:
    integer                  :: MAX_VSTRUCT           ! Maximum variogram structures
    integer                  :: DATA_CLASSES          ! Number of texture & covariate classes
    ! "Optional" Options - Default Value and Settings
    real                     :: NO_DATA_VALUE         !
    real                     :: IDW_POWER             ! Inverse Distance Weighting Power Parameter
    character(10)            :: INTERP_DIM            ! 2D, 3D
    character(10)            :: INTERP_METHOD         ! KRIGING, IDW
    character(10)            :: INTERP_ORDER          ! BEFORE, AFTER
    character(60)            :: MAIN_OUTPUT           ! 
    !character(100)           :: OUTPUT_DIR            ! Directory where output gets written
    !character(100)           :: INPUT_DIR             ! Direction where input gets read
    logical                  :: SILENT_MODE           !
    logical                  :: DEBUG_MODE            !
    logical                  :: WRITE_MATRIX          !
    logical                  :: WRITE_PP_WEIGHTS      !
    logical                  :: WRITE_GRID_WEIGHTS    !
    logical                  :: READ_PP_WEIGHTS       !
    logical                  :: READ_GRID_WEIGHTS     !
    logical                  :: WRITE_NONE            !
    logical                  :: WRITE_NODE_FILES      !
    logical                  :: WRITE_DATASET_FILES   !
    logical                  :: OUTPUT_GRID           !
    logical                  :: OUTPUT_ARRAY          !
    logical                  :: LOG_KRIGE_K           !
    logical                  :: INFER_LAST_CLASS      ! For classic T2P, where the data is only percent coarse
    
    contains
      procedure, public :: initialize
      procedure, public :: read_options
    
  end type t_options

!-------------------------------------------------------------------------------------------------!
  
  contains

!-------------------------------------------------------------------------------------------------!

  subroutine initialize(this)
    implicit none
    class(t_options),intent(inout)         :: this

    this%MAX_VSTRUCT          = -1
    this%DATA_CLASSES         = -1
    this%NO_DATA_VALUE        = -999.0
    this%IDW_POWER            =  2.0
    this%INTERP_DIM           = '2D'
    this%INTERP_METHOD        = 'KRIGING'
    this%INTERP_ORDER         = 'BEFORE'
    this%MAIN_OUTPUT          = 't2p_log.txt'
    !this%OUTPUT_DIR           = './'
    !this%INPUT_DIR            = './'
    this%SILENT_MODE          = .false.
    this%DEBUG_MODE           = .false.
    this%WRITE_MATRIX         = .false.
    this%WRITE_PP_WEIGHTS     = .false.
    this%WRITE_GRID_WEIGHTS   = .false.
    this%READ_PP_WEIGHTS      = .false.
    this%READ_GRID_WEIGHTS    = .false.
    this%WRITE_NONE           = .false.
    this%WRITE_NODE_FILES     = .false.
    this%WRITE_DATASET_FILES  = .false.
    this%OUTPUT_GRID          = .false.
    this%OUTPUT_ARRAY         = .false.
    this%LOG_KRIGE_K          = .false.
    this%INFER_LAST_CLASS     = .false.
  end subroutine initialize

!-------------------------------------------------------------------------------------------------!
  
subroutine read_options(this, reader)
  use m_file_io, only: t_file_reader, item2int, item2char, item2real
  use m_vstringlist, only: t_vstringlist
  use m_error_handler, only: error_handler
  implicit none
  class(t_options)           :: this
  type(t_file_reader), pointer :: reader
  integer                    :: status, length
  character(30)              :: id
  type(t_vstringlist)        :: strings
  
  ! To hold anything read
  integer                    :: in_int
  real                       :: in_real
  character(30)              :: in_char

  do
    call reader%next_block_item(status, id, strings, length)
    if (status /= 0) exit  ! exit if end of block or end of file
    select case(trim(id))
      case("MAX_VSTRUCT")
        ! TODO LENGTH CHECK
        this%MAX_VSTRUCT          = item2int(strings, 2)
        ! WHAT HAPPENS IF THIS FAILS?
      case("DATA_CLASSES")
        this%DATA_CLASSES         = item2int(strings, 2)
      case("INTERP_DIM")
        call item2char(strings, 2, this%INTERP_DIM)
      case("INTERP_METHOD")
        call item2char(strings, 2, this%INTERP_METHOD)
      case("INTERP_ORDER")
        call item2char(strings, 2, this%INTERP_ORDER)
      case("MAIN_OUTPUT")
        call item2char(strings, 2, this%MAIN_OUTPUT)
      !case("OUTPUT_DIR")
      !  call item2char(strings, 2, this%OUTPUT_DIR)
      !case("INPUT_DIR")
      !  call item2char(strings, 2, this%INPUT_DIR)
      case("NO_DATA_VALUE")
        this%NO_DATA_VALUE        = item2real(strings, 2)
      case("IDW_POWER")
        this%IDW_POWER            = item2real(strings, 2)
      case("SILENT_MODE")
        this%SILENT_MODE          = .true.
      case("DEBUG_MODE")
        this%DEBUG_MODE           = .true.
      case("WRITE_MATRIX")
        this%WRITE_MATRIX         = .true.
      case("WRITE_PP_WEIGHTS")
        this%WRITE_PP_WEIGHTS     = .true.
      case("WRITE_GRID_WEIGHTS")
        this%WRITE_GRID_WEIGHTS   = .true.
      case("READ_PP_WEIGHTS")
        this%READ_PP_WEIGHTS      = .true.
      case("READ_GRID_WEIGHTS")
        this%READ_GRID_WEIGHTS    = .true.
      case("WRITE_NONE")
        this%WRITE_NONE           = .true.
      case("WRITE_NODE_FILES")
        this%WRITE_NODE_FILES     = .true.
      case("WRITE_DATASET_FILES")
        this%WRITE_DATASET_FILES  = .true.
      case("OUTPUT_GRID")
        this%OUTPUT_GRID          = .true.
      case("OUTPUT_ARRAY")
        this%OUTPUT_ARRAY         = .true.
      case("LOG_KRIGE_K")
        this%LOG_KRIGE_K          = .true.
      case("INFER_LAST_CLASS")
        this%INFER_LAST_CLASS     = .true.
      case DEFAULT
        call error_handler(1,reader%file,"Unknown Option Name: " // trim(id))
    end select
  end do
  
  ! Handle EOF
    
end subroutine read_options

end module m_options
