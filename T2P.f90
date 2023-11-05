module m_T2P
  use m_vstringlist, only: t_vstringlist
  use m_options, only: t_options
  use m_flow_model, only: t_flow_model
!-------------------------------------------------------------------------------------------------!
  type t_T2P
!    type(t_file_writer), pointer     :: log   ! Runtime Log
    type(t_options)                  :: opt  ! avoid conflict with Fortran intrinsic 'options'
    type(t_flow_model), pointer      :: fmodel
    !type(t_Interpolator)             :: interpolator
    !type(t_Dataset),allocatable()    :: datasets
    !type(t_Variogram),allocatable()  :: variograms
    !type(t_PilotPoint),allocatable() :: pilot_points
    
    ! Counters
    integer                          :: nppzones  ! Number of Pilot Point Zones
    integer                          :: nhsus     ! Number of Hydrostratigraphic Units (HSUs)
                                                  ! (Formerly "geologic zones")
  
  contains
  
    procedure :: initialize
    !procedure :: allocate
  
  end type t_T2P
!-------------------------------------------------------------------------------------------------!
  
  contains
  
  subroutine initialize(this)
    use m_file_io, only: open_file_writer, log_unit
    implicit none
    class(t_T2P), intent(inout)      :: this
    
    ! Initialize Elements
    call this%opt%initialize()
    
    ! Set class defaults
    this%nppzones = 1
    this%nhsus    = 1
    
  end subroutine initialize
  
!-------------------------------------------------------------------------------------------------!
  

!-------------------------------------------------------------------------------------------------!
  
end module m_T2P