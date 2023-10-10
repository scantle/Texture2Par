module m_T2P
  use m_vstringlist, only: t_vstringlist
  use m_options, only: t_options
  use m_file_io, only: t_file_writer
!-------------------------------------------------------------------------------------------------!
  type t_T2P
    type(t_file_writer), pointer     :: log   ! Runtime Log
    type(t_options)                  :: opt  ! avoid conflict with Fortran intrinsic 'options'
    type(t_vstringlist)              :: warnings
    !type(t_FlowModel)                :: fmodel
    !type(t_Interpolator)             :: interpolator
    !type(t_Dataset),allocatable()    :: datasets
    !type(t_Variogram),allocatable()  :: variograms
    !type(t_PilotPoint),allocatable() :: pilot_points
  
  contains
  
    procedure :: initialize
    procedure :: log_init
    !procedure :: allocate
    
    ! Exit Routines (Hooks)
    procedure :: early_exit
    procedure :: normal_exit
    
    ! Error Routines
    procedure :: error_handler
    procedure :: warning_handler
  
  end type t_T2P
!-------------------------------------------------------------------------------------------------!
  
  contains
  
  subroutine initialize(this)
    use m_file_io, only: open_file_writer, log_unit
    implicit none
    class(t_T2P), intent(inout)      :: this
    
    ! Initialize Log
    this%log => open_file_writer("log.t2p")
    log_unit = this%log%unit
    
    ! Initialize Elements
    call this%opt%initialize()
    
  end subroutine initialize
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine early_exit(this)
    use m_file_io, only: close_file
    implicit none
    class(t_T2P), intent(inout)      :: this
    
    call close_file(this%log)
    
    write(*,'(a)') 'Texture2Par Early Termination.'
    error stop
    
  end subroutine early_exit
    
!-------------------------------------------------------------------------------------------------!
  
  subroutine normal_exit(this)
    use m_file_io, only: close_file
    implicit none
    class(t_T2P), intent(inout)      :: this
    
    call close_file(this%log)
    
    write(*,'(a)') 'Texture2Par Normal Termination.'
    
  end subroutine normal_exit
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine error_handler(this,err,filename,opt_msg)
    implicit none
    class(t_T2P), intent(inout)      :: this
    integer, intent(in)              :: err
    character(*),intent(in),optional :: filename
    character(*),intent(in),optional :: opt_msg
    character(200)                   :: err_msg
    
    ! Start message
    err_msg = "Error -"
    
    ! All available errors
    select case(err)
      case(0) ! Not an error
        ! Important to cover all your bases :)
        return ! Returns with no error
      case(1) ! Generic Input Error
        err_msg = trim(err_msg) // " Invalid Input -"
      case(2) ! Missing Data Error
        err_msg = trim(err_msg) // " Missing Required Data -"
      case(3) ! Invalid data, unable to convert
        err_msg = trim(err_msg) // " Data Format/Conversion Failure -"
      case DEFAULT
        write(*,*) "An unknown error has occured. Please contact the developer: lelands@sspa.com"
    end select
    
    ! Add optional stuff
    if (present(opt_msg)) then
      err_msg = trim(err_msg) //" "// trim(opt_msg) // "."
    end if
    if (present(filename)) then
      err_msg = trim(err_msg) // " File: " // trim(filename)
    end if
    
    ! Write the error
    write(*,'(a)') err_msg
    
    ! Call early_exit handler
    call this%early_exit()
  
  end subroutine error_handler
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine warning_handler(this, warning)
    use m_vstringlist, only: vstrlist_append
    implicit none
    class(t_T2P), intent(inout)      :: this
    character(*), intent(in)         :: warning
    
    call vstrlist_append(this%warnings, warning)
  
  end subroutine warning_handler

!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
  subroutine log_init(this)
    implicit none
    class(t_T2P), intent(inout)      :: this
  
    call this%log%write_line('    _________   _______  ________   ')
    call this%log%write_line('--|\___   ___\/  ___  \|\   __  \ ----')
    call this%log%write_line('   \|___ \  \_/__/|_/  /\ \  \|\  \ ')
    call this%log%write_line('------- \ \  \|__|//  / /\ \   ____\----')
    call this%log%write_line('         \ \  \   /  /_/__\ \  \___|')
    call this%log%write_line('--------- \ \__\ |\________\ \__\ --------')
    call this%log%write_line('           \|__|  \|_______|\|__|   ', blank_end=1)
    call this%log%write_line('             Texture2Par v2.0.0             ')
    call this%log%write_line( '----------------------------------------------')
    call this%log%write_line('Copyright S.S. Papadopulos & Associates Inc.', blank_end=1)
    
    ! TODO write time of start
  
  end subroutine log_init
!-------------------------------------------------------------------------------------------------!  
  
end module m_T2P