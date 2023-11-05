module m_log
  use m_file_io, only: open_file_writer, t_file_writer
  
  contains
  
!-------------------------------------------------------------------------------------------------!
  subroutine log_init(log)
    use m_file_io, only: log_unit
    implicit none
  
    type(t_file_writer),pointer  :: log
  
    ! Initialize
    log => open_file_writer("log.t2p")
    ! To avoid a circular dependency with the FileIO class
    log_unit = log%unit
  
    call log%write_line('    _________   _______  ________   ')
    call log%write_line('--|\___   ___\/  ___  \|\   __  \ ----')
    call log%write_line('   \|___ \  \_/__/|_/  /\ \  \|\  \ ')
    call log%write_line('------- \ \  \|__|//  / /\ \   ____\----')
    call log%write_line('         \ \  \   /  /_/__\ \  \___|')
    call log%write_line('--------- \ \__\ |\________\ \__\ --------')
    call log%write_line('           \|__|  \|_______|\|__|   ', blank_end=1)
    call log%write_line('             Texture2Par v2.0.0             ')
    call log%write_line( '----------------------------------------------')
    call log%write_line('Copyright S.S. Papadopulos & Associates Inc.', blank_end=1)
    
    ! TODO write time of start
  
  end subroutine log_init
!-------------------------------------------------------------------------------------------------!  
  
  
end module m_log