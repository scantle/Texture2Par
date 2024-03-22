module m_log
  use m_file_io, only: open_file_writer, t_file_writer
  use t_time
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
    call write_time_start(log_unit,'Started at: ')
    call start_timer()
  
  end subroutine log_init

!-------------------------------------------------------------------------------------------------!  
  
  subroutine write_elapsed_time(log)
    implicit none
    type(t_file_writer),pointer  :: log
    integer                      :: elapsed_time
    integer                      :: hours, minutes, seconds
    character(100)               :: line
    
    elapsed_time = get_elapsed_time()
    
    hours = int(elapsed_time / 3600)
    minutes = int(mod(elapsed_time, 3600) / 60)
    seconds = int(mod(elapsed_time, 60))

    if (hours > 0) then
      write(line,'(a,i3,a,(i2,a))') "Elapsed time: ", hours, " hours, ", minutes, " minutes, ", seconds, " seconds."
    else if (minutes > 0) then
      write(line,'(a,2(i2,a))') "Elapsed time: ", minutes, " minutes, ", seconds, " seconds."
    else
      write(line,'(a,i2,a)') "Elapsed time: ", seconds, " seconds."
    endif
    
    call log%write_line(line, blank_start=1)
  
  end subroutine write_elapsed_time

!-------------------------------------------------------------------------------------------------!
  
  subroutine write_warnings(log)
    use m_error_handler, only: warnings
    use m_vstringlist, only: vstrlist_length, vstrlist_index
    use m_vstring
    implicit none
    type(t_file_writer),pointer  :: log
    integer           :: i
  
    if (vstrlist_length(warnings) < 1) return   ! Early return if no warnings
    
    call log%write_line('WARNINGS:', blank_start=1)
    do i=1, vstrlist_length(warnings)
      call log%write_line(vstrlist_index(warnings,i))
    end do
  
  end subroutine write_warnings
  
!-------------------------------------------------------------------------------------------------!
  
end module m_log