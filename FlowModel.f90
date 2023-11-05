module m_flow_model
  use m_file_io, only: t_file_reader, t_file_writer
  use m_grid, only: t_grid
  implicit none
!-------------------------------------------------------------------------------------------------!
  type t_flow_model  ! Abstract Base Class, also used for generic grid output
  
    integer                       :: nnodes        ! Number of nodes (points at which flow is solved)
    integer                       :: nlayers       ! Number of model layers
    logical                       :: structured    ! T/F Structured (row/col) model
    type(t_file_reader), pointer  :: template_file ! File copied for output_file
    type(t_file_writer), pointer  :: output_file   ! File where aquifer parameters are written
    type(t_Grid), allocatable     :: layer(:)      ! Grid for ground surface + each layer bottom
  
    contains
      procedure,public              :: initialize
      !procedure,private     :: write_model_input
      !procedure,private     ::
  
  end type t_flow_model
!-------------------------------------------------------------------------------------------------!
  type, extends(t_flow_model) :: t_iwfm_model
  
  contains
    procedure,public              :: initialize => init_iwfm
  
  end type 
!-------------------------------------------------------------------------------------------------!
  type, extends(t_flow_model) :: t_modflow_struct_model
  
  contains
    procedure,public              :: initialize => init_modflow_struct
  
  end type 
!-------------------------------------------------------------------------------------------------!  
  
  contains
  
!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine read_FLOW_MODEL_block(reader)
    use m_vstringlist, only: t_vstringlist
    implicit none
  ! Called during main input file read
    type(t_file_reader), pointer :: reader
    ! Standard file reader variables
    integer                    :: status, length
    character(30)              :: id
    type(t_vstringlist)        :: strings
    
    
    call reader%next_block_item(status, id, strings, length)
    do
      if (status /= 0) exit  ! exit if end of block or end of file
      select case(trim(id))
        case("LOG_KRIGE_K")
          write(*,*) 'TESTING'
        case DEFAULT
          ! TODO! Raise error!
          write(*,*) "BAD BLOCK!"
          stop
      end select
    end do
  
  end subroutine read_FLOW_MODEL_block
  
  ! I think we'll need a subroutine to create the pointer
  
!-------------------------------------------------------------------------------------------------!
! BASE CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine initialize(this)
    implicit none
    class(t_flow_model)       :: this
    
  end subroutine initialize

!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! IWFM CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

  subroutine init_iwfm(this)
    implicit none
    class(t_iwfm_model)       :: this
    
  end subroutine init_iwfm  
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! MODFLOW STRUCTURED CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine init_modflow_struct(this)
    implicit none
    class(t_modflow_struct_model)       :: this
    
  end subroutine init_modflow_struct
  
!-------------------------------------------------------------------------------------------------!
end module m_flow_model