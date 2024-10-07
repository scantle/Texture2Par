module m_options
!-------------------------------------------------------------------------------------------------!  
! HOW TO ADD AN OPTION
!  1. Add it as a variable in the type t_options
!  2. Add it to the initialize routine and give it a default value
!  3. Create a new case for it in the read_options routine. Write code to handle assignment.
!  4. Add it to the write_input_summary routine
!  5. If necessary, add a handler for the option to the t2p%process_options() hook that is called
!     after the options are read.
!-------------------------------------------------------------------------------------------------!
  
  type t_options
    ! Required Options:
    integer                  :: MAX_VSTRUCT           ! Maximum variogram structures
    !integer                  :: DATA_CLASSES          ! Number of texture & covariate classes
    ! "Optional" Options - Default Value and Settings
    real                     :: NO_DATA_VALUE         !
    real                     :: IDW_POWER             ! Inverse Distance Weighting Power Parameter
    real                     :: MAX_SECONDARY_DIST
    integer                  :: INTERP_DIM            ! 2D, 3D, stored as integer, read as string
    integer                  :: NSIM                  ! number of sequential gaussian simulation; default is zero which is Kriging
    integer                  :: SEED                  ! seed number for generating random path and samples for SGSIM
    character(10)            :: INTERP_METHOD         ! KRIGING, IDW
    character(10)            :: INTERP_ORDER          ! Interpolate to grid: [BEFORE, AFTER] upscaling/aggegating?
    character(60)            :: MAIN_OUTPUT           ! 
    !character(100)           :: OUTPUT_DIR            ! Directory where output gets written
    !character(100)           :: INPUT_DIR             ! Direction where input gets read
    logical                  :: SILENT_MODE           !
    logical                  :: DEBUG_MODE            !
    logical                  :: ONLY_TEXTURE          ! Stops after calculating texture
    logical                  :: WRITE_MATRIX          !
    logical                  :: WRITE_PP_WEIGHTS      !
    logical                  :: WRITE_GRID_WEIGHTS    !
    logical                  :: READ_PP_WEIGHTS       !
    logical                  :: READ_GRID_WEIGHTS     !
    logical                  :: WRITE_NONE            ! Skips model output files (as implemented...)
    logical                  :: WRITE_NODE_FILES      ! Writes parameters by node/cell
    logical                  :: WRITE_DATASET_FILES   ! Writes out datasets averaged by layer
    logical                  :: OUTPUT_GRID           !
    logical                  :: OUTPUT_ARRAY          !
    logical                  :: LOG_KRIGE_K           !
    logical                  :: INFER_LAST_CLASS      ! For classic T2P, where the data is only percent coarse
    logical                  :: CORRECT_WEIGHTS       ! Zero negative (kriging) weights and rebalance if true, see Deutsch (1996) (only for OK, COK)
    
    contains
      procedure, public :: initialize
      procedure, public :: read_options
      procedure, public :: write_input_summary
    
  end type t_options

!-------------------------------------------------------------------------------------------------!
  
  contains

!-------------------------------------------------------------------------------------------------!

  subroutine initialize(this)
    implicit none
    class(t_options),intent(inout)         :: this

    this%MAX_VSTRUCT          = -1
    !this%DATA_CLASSES         = -1
    this%NO_DATA_VALUE        = -999.0
    this%IDW_POWER            = 2.0
    this%MAX_SECONDARY_DIST   = 0.0 
    this%INTERP_DIM           = 2
    this%NSIM                 = 0
    this%SEED                 = 0
    this%INTERP_METHOD        = 'OK'
    this%INTERP_ORDER         = 'BEFORE'
    this%MAIN_OUTPUT          = 't2p_log.txt'
    !this%OUTPUT_DIR           = './'
    !this%INPUT_DIR            = './'
    this%SILENT_MODE          = .false.
    this%DEBUG_MODE           = .false.
    this%ONLY_TEXTURE         = .false.
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
    this%CORRECT_WEIGHTS      = .false.
  end subroutine initialize

!-------------------------------------------------------------------------------------------------!
  
subroutine read_options(this, reader)
  use m_file_io, only: t_file_reader, item2int, item2char, item2real
  use m_vstringlist, only: t_vstringlist
  use m_error_handler, only: error_handler
  use tools, only: random_seed_initialize
  implicit none
  class(t_options)           :: this
  type(t_file_reader), pointer :: reader
  integer                    :: status, length
  character(30)              :: id, temp
  type(t_vstringlist)        :: strings

  do
    call reader%next_block_item(status, id, strings, length)
    if (status /= 0) exit  ! exit if end of block or end of file
    select case(trim(id))
      case("MAX_VSTRUCT")
        this%MAX_VSTRUCT          = item2int(strings, 2)
      !case("DATA_CLASSES")
      !  this%DATA_CLASSES         = item2int(strings, 2)
      case("INTERP_DIM")
        call item2char(strings, 2, temp,toupper=.true.)
        select case (trim(temp))
          case("2","2D")
            this%INTERP_DIM = 2
          case("3","3D")
            this%INTERP_DIM = 3
          case DEFAULT
            call error_handler(1,reader%file,"Invalid INTERP_DIM: " // trim(id))
        end select
      case("SGSIM")
        this%NSIM              = 1
      case("SEED")
        this%SEED              = item2int(strings, 2)
      case("INTERP_METHOD")
        call item2char(strings, 2, this%INTERP_METHOD,toupper=.true.)
      case("INTERP_ORDER")
        call item2char(strings, 2, this%INTERP_ORDER,toupper=.true.)
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
      case("MAX_SECONDARY_DIST")
        this%MAX_SECONDARY_DIST   = item2real(strings, 2)
      case("SILENT_MODE")
        this%SILENT_MODE          = .true.
      case("DEBUG_MODE")
        this%DEBUG_MODE           = .true.
      case("ONLY_TEXTURE")
        this%ONLY_TEXTURE         = .true.
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
      case("WRITE_NODE_FILES","WRITE_CELL_FILES")
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
      case("CORRECT_WEIGHTS")
        this%CORRECT_WEIGHTS      = .true.
      case DEFAULT
        call error_handler(1,reader%file,"Unknown Option Name: " // trim(id))
    end select
  end do
  !call random_seed_initialize(this%SEED)
  ! Handle EOF
    
end subroutine read_options

!-------------------------------------------------------------------------------------------------!

subroutine write_input_summary(this)
  use m_global, only: log
  implicit none
  class(t_options)           :: this
  character(200)             :: logic_true, interp_text
  integer                    :: no_true_len

    call log%write_valueline('Max Variogram Structures:', this%MAX_VSTRUCT)
    call log%write_valueline('No Data Value:',this%NO_DATA_VALUE)
    call log%write_valueline('IDW Power Value:',this%IDW_POWER)
    call log%write_valueline('Max Secondary Data Dist:', this%MAX_SECONDARY_DIST)
    !call log%write_valueline('',this%INTERP_DIM)  ! unimplemented
    ! Going to do something a little nicer for the Interpolation Method
    select case(this%INTERP_METHOD)
      case("OK")
        interp_text = "Ordinary Kriging (OK)"
      case("SK")
        interp_text = "Simple Kriging (SK)"
      case("IDW")
        interp_text = "Inverse Distance Weighting (IDW)"
      case DEFAULT  ! Shouldn't be able to get here
        interp_text = this%INTERP_METHOD
    end select
    call log%write_valueline('Interpolation Method:',trim(interp_text))
    call log%write_valueline('Interpolation Order:',this%INTERP_ORDER)
    write(interp_text, *) this%NSIM; call log%write_valueline('Realization Number:', trim(interp_text))
    write(interp_text, *) this%SEED; call log%write_valueline('Random Seed Number:', trim(interp_text))
    !call log%write_valueline('',this%MAIN_OUTPUT) ! unneeded (we're writing to that file!)
    !call log%write_valueline('',this%OUTPUT_DIR)  ! unimplemented
    !call log%write_valueline('',this%INPUT_DIR)   ! unimplemented
    
    ! Logical variables
    logic_true = ''
    no_true_len = 0 !len(trim(logic_true))
    
    if (this%SILENT_MODE        ) logic_true = trim(logic_true) // "SILENT_MODE"         // ","
    if (this%DEBUG_MODE         ) logic_true = trim(logic_true) // "DEBUG_MODE"          // ","
    if (this%ONLY_TEXTURE       ) logic_true = trim(logic_true) // "ONLY_TEXTURE"        // ","
    if (this%WRITE_MATRIX       ) logic_true = trim(logic_true) // "WRITE_MATRIX"        // ","
    if (this%WRITE_PP_WEIGHTS   ) logic_true = trim(logic_true) // "WRITE_PP_WEIGHTS"    // ","
    if (this%WRITE_GRID_WEIGHTS ) logic_true = trim(logic_true) // "WRITE_GRID_WEIGHTS"  // ","
    if (this%READ_PP_WEIGHTS    ) logic_true = trim(logic_true) // "READ_PP_WEIGHTS"     // ","
    if (this%READ_GRID_WEIGHTS  ) logic_true = trim(logic_true) // "READ_GRID_WEIGHTS"   // ","
    if (this%WRITE_NONE         ) logic_true = trim(logic_true) // "WRITE_NONE"          // ","
    if (this%WRITE_NODE_FILES   ) logic_true = trim(logic_true) // "WRITE_NODE_FILES"    // ","
    if (this%WRITE_DATASET_FILES) logic_true = trim(logic_true) // "WRITE_DATASET_FILES" // ","
    if (this%OUTPUT_GRID        ) logic_true = trim(logic_true) // "OUTPUT_GRID"         // ","
    if (this%OUTPUT_ARRAY       ) logic_true = trim(logic_true) // "OUTPUT_ARRAY"        // ","
    if (this%LOG_KRIGE_K        ) logic_true = trim(logic_true) // "LOG_KRIGE_K"         // ","
    if (this%INFER_LAST_CLASS   ) logic_true = trim(logic_true) // "INFER_LAST_CLASS"    // ","
    if (this%CORRECT_WEIGHTS    ) logic_true = trim(logic_true) // "CORRECT_WEIGHTS"     // ","
    
    if (len(trim(logic_true)) > no_true_len) then
      !call log%write_line(logic_true)
      call log%write_valueline('Logical Options Set to TRUE:',logic_true)
    else
      !call log%write_line("All Logical Options Set to FALSE")
      call log%write_valueline('All Logical Options Set to','FALSE')
    end if
  
end subroutine write_input_summary
!-------------------------------------------------------------------------------------------------!

end module m_options
