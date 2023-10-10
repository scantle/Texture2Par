program Texture2Par
  use m_file_io, only: io_initialize
  use m_read_main_input, only: read_main_input
  use m_global, only: t2p
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
! Copyright 2022 S.S. Papadopulos & Associates. All rights reserved.
!-----------------------------------------------------------------------------!

  character(100)   :: main_input_file = 'Texture2Par.in'
  
  ! Write Flag
  write(*,'(/,4x,a)') '             Texture2Par v2.0.0             '
  write(*,'(3x,a)') '----------------------------------------------'
  write(*,'(4x,a)')   'Copyright S.S. Papadopulos & Associates Inc.'
    
  ! Pre-read initializations
  call io_initialize()
  call t2p%initialize()
  
  ! Read Main Input
  call read_main_input(main_input_file)
    
  ! Check Major Inputs
    
    
  ! Allocate
    
    
  ! Read Sub-Input Files
    
    
  ! Check Sub Inputs
    
    
  ! Main
    
    
  ! Output
    
    
  ! Exit
  call t2p%normal_exit()

  end program Texture2Par
!-----------------------------------------------------------------------------!