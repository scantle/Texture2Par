program Texture2Par
  use m_global, only: log
  use m_T2P, only: t_t2p
  use m_file_io, only: io_initialize, get_command_args
  use m_read_main_input, only: read_main_input
  use m_log, only: log_init, write_warnings, write_elapsed_time
  use m_error_handler, only: normal_exit
  use m_read_datasets, only: read_datasets
  implicit none

!-----------------------------------------------------------------------------!
!                             Texture2Par
!                         ********************
!                               V 2.0.0
!
! Converts well percent course with depth data to aquifer parameters
! at specified nodes. Intended to work with input/output files of the
! Integrated Water Flow Model (IWFM) with added support for MODFLOW 2000/2005.
! Written by Leland Scantlebury and Marinko Karanovic of S.S. Papadopulos &
! Associates, based on original method and VBA program created by Timothy J.
! Durbin.
!
! Requires Input File: Texture2Par.in
!
! Copyright 2024 S.S. Papadopulos & Associates. All rights reserved.
!-----------------------------------------------------------------------------!

  character(100)      :: main_input_file = 'Texture2Par.in'
  type(t_T2P)         :: t2p   ! main Texture2Par storage
  
  ! Write Flag
  write(*,'(/,4x,a)') '             Texture2Par v2.0.0             '
  write(*,'(3x,a)') '----------------------------------------------'
  write(*,'(4x,a)')   'Copyright S.S. Papadopulos & Associates Inc.'
    
  ! Pre-read initializations
  call get_command_args(main_input_file)
  call io_initialize()
  call log_init(log)
  call t2p%initialize()
  
  ! Read Main Input
  call read_main_input(t2p, main_input_file)
  call t2p%fmodel%read_model()
  ! Read datasets
  call read_datasets(t2p)
  call write_elapsed_time(log)
  ! Check Major Inputs/write warnings
  
  call write_warnings(log)
    
  ! Main
  call t2p%predict()
    
  ! Output
  !call t2p%write_output()
    
  ! Exit
  call write_elapsed_time(log)
  call normal_exit()

  end program Texture2Par
!-----------------------------------------------------------------------------!