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
  
  subroutine write_input_summary(log, t2p)
    use m_T2P, only: t_t2p
    use m_flow_model
    implicit none
    type(t_file_writer),pointer  :: log
    type(t_t2p)                  :: t2p
    
    call log%write_line('--- INPUT SUMMARY ---', blank_start=1)

    ! T2P
    call log%write_valueline('Total Data Classes:',   t2p%nclasses)
    call log%write_valueline('Primary Data Classes:',   t2p%nprimary)
    call log%write_valueline('Covariate Data Classes:', t2p%nsecondary)
    call log%write_valueline('Number of Dataset Files:',t2p%ndatafiles)
    call log%write_valueline('Number of HSUs:',         t2p%nhsus)
    call log%write_valueline('Number of Variograms:',   t2p%nvario)
    
    ! Flow Model
    call t2p%fmodel%write_input_summary()
    
    ! Pilot Points
    call log%write_valueline('Number of Pilot Points',      t2p%pilot_points%npp)
    if (t2p%pilot_points%has_aquitard_pp) then
      call log%write_valueline('Number of Aquifer Pilot Points',   count(t2p%pilot_points%is_aquifer))
      call log%write_valueline('Number of Aquitard Pilot Points',  count(t2p%pilot_points%is_aquitard))
    end if
    !call log%write_valueline('Number of Pilot Point Zones', t2p%pilot_points%nppzones)
    
    ! Options
    call t2p%opt%write_input_summary() ! Move here?
    
  end subroutine write_input_summary

!-------------------------------------------------------------------------------------------------!
  
  subroutine write_warnings(log)
    use m_error_handler, only: warnings
    use m_vstringlist, only: vstrlist_length, vstrlist_index
    use m_vstring
    implicit none
    type(t_file_writer),pointer  :: log
    integer           :: i
  
    if (vstrlist_length(warnings) < 1) return   ! Early return if no warnings
    
    call log%write_line('--- WARNINGS ---', blank_start=1)
    do i=1, vstrlist_length(warnings)
      call log%write_line(vstrlist_index(warnings,i))
    end do
  
  end subroutine write_warnings
  
!-------------------------------------------------------------------------------------------------!

  
end module m_log