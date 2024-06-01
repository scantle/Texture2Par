module m_write_fmodel
  use m_global, only: dp

  contains
!-----------------------------------------------------------------------------!
! IWFM File copying subroutines
!-----------------------------------------------------------------------------!  
subroutine writeIWFMgwfile(gwtemp, gwfile, par, nlayers, nnodes, nodeid, is_aquitard, is_active_layer)
  use m_file_io, only: t_file_reader, t_file_writer, open_file_reader, open_file_writer
  use m_datasets, only: t_layerdata
  implicit none

  !integer, intent(in)       :: nskip
  character(100),intent(in)  :: gwtemp, gwfile
  integer,intent(in)         :: nlayers, nnodes
  integer,intent(in)         :: nodeid(nnodes)
  type(t_layerdata),pointer  :: par(:)            ! Final model hydraulic parameters (nlayers)&values([Kh, Kv, Ss, Sy], npoints)
  logical,intent(in)         :: is_aquitard(nlayers), is_active_layer(nlayers)
  integer                    :: i, k, ierr, nouth, noutf, ngroup, ru, wu
  character(256)             :: cerr
  character(1000)            :: line
  real(dp)                   :: KvB_aquitard
  type(t_file_reader),pointer:: fread
  type(t_file_writer),pointer:: fwrite
  

  write(*,'(2x,a)') 'Reading gw temp file, writing new gw file'

  !open(ru, file=trim(gwtemp), status='old', iostat = ierr, iomsg=cerr)   ! Temp File
  !call iomsghandler(ierr, cerr)
  !open(12, file=trim(GWFile), status='replace', iostat = ierr, iomsg=cerr)   ! New File
  !call iomsghandler(ierr, cerr)
  fread  => open_file_reader(gwtemp)
  ru = fread%unit
  fwrite => open_file_writer(gwfile)
  wu = fwrite%unit

  ! Copy Section Number
  call CopySettingsSection(ru, wu, 1)
  ! Copy over the file header
  call CopyCommentSection(ru, wu)
  ! Copy over main groundwater file settings (19 lines in V4.0)
  call CopySettingsSection(ru, wu)
  ! Copy Debugging
  call CopyCommentSection(ru, wu)
  call CopySettingsSection(ru, wu)
  ! Copy Groundwater Hydrograph Output Settings
  call CopyCommentSection(ru, wu)
  read(ru, '(a1000)', iostat = ierr) line
  read(line, *) nouth
  write(wu,'(a)') trim(line)
  call CopySettingsSection(ru, wu)
  ! Copy Hydrogrph data
  call CopyCommentSection(ru, wu)
  call CopySettingsSection(ru, wu, nouth)
  ! Copy Element Face Flow Output Settings
  call CopyCommentSection(ru, wu)
  read(ru, '(a1000)', iostat = ierr) line
  read(line, *) noutf
  write(wu,'(a)') trim(line)
  call CopySettingsSection(ru, wu)
  ! Copy Element Face Data
  call CopyCommentSection(ru, wu)
  call CopySettingsSection(ru, wu, noutf)
  ! Copy Aquifer Parameters, Handle NGROUP setting
  call CopyCommentSection(ru, wu)
  read(ru, '(a1000)', iostat = ierr) line
  write(wu,'(a)') trim(line)
  read(line, *) ngroup
  if (ngroup > 0) then
    write(*,'(a)') 'ERROR - NGROUP > 0 in Groundwater file (parametric grid)'
    write(*,'(a)') 'Only NGROUP = 0 is supported (node parameters)'
    stop
  end if
  ! Copy all the way down to aquifer parameter section
  call CopyCommentSection(ru, wu)
  call CopySettingsSection(ru, wu) ! Conversion Factors
  call CopyCommentSection(ru, wu)
  call CopySettingsSection(ru, wu) ! Units
  call CopyCommentSection(ru, wu)

  ! Finally ready to write the aquifer parameters!
  do i=1, nnodes
    do k = 2, nlayers, 2  ! Aquitards are setup as layers on top of each aquifer layer
      KvB_aquitard = 0.0
      if (is_aquitard(k-1).and.is_active_layer(k-1))  KvB_aquitard = par(k-1)%values(2,i)
      if (k.eq.2) then  ! nodeid, Kh, Ss, Sy, KvB_aquitard KvB
        write(wu,'(i13,f14.8,e12.3,f14.8,f14.8,f14.8)') nodeid(i), par(k)%values(1,i), par(k)%values(3:4,i), KvB_aquitard, par(k)%values(2,i)
      else
        write(wu,'(13x,f14.8,e12.3,f14.8,f14.8,f14.8)') par(k)%values(1,i), par(k)%values(3:4,i), KvB_aquitard, par(k)%values(2,i)
      end if
      read(ru,*)
    end do
  end do

  ! Copy the rest of the file
  do
    read(ru, '(a1000)', iostat = ierr) line
    if (ierr /= 0) exit
    write(wu,'(a)') adjustl(trim(line))
  end do

  call fread%close_file()
  call fwrite%close_file()

end subroutine writeIWFMgwfile

!-----------------------------------------------------------------------------!
! IWFM File copying subroutines
!-----------------------------------------------------------------------------!

subroutine CopyCommentSection(copyunit, newunit)
  implicit none

  integer,intent(in)       :: copyunit, newunit
  character(1)             :: compare
  character(1000)          :: line

  ! Copy lines until there are no more comment lines
  do
    read(copyunit, '(a1000)') line
    compare = ADJUSTL(line)
    if ((compare == 'C').or.(compare == 'c').or.&
        (compare == '*').or.(compare == '#')) then
      write(newunit,'(a)') trim(line)
    else
      ! Non comment found, back one line and exit
      backspace(copyunit)
      exit
    end if
  end do

end subroutine
!-----------------------------------------------------------------------------!
subroutine CopySettingsSection(copyunit, newunit, nlines)
  implicit none

  integer,intent(in)       :: copyunit, newunit
  integer                  :: i
  character(1)             :: compare
  character(1000)          :: line
  integer,intent(in),optional  :: nlines

  ! If nlines is passed, just copy over that many lines
  if (present(nlines)) then
    do i=1, nlines
      read(copyunit, '(a1000)') line
      write(newunit,'(a)') trim(line)
    end do
  else
    ! Copy lines until a comment line is found
    do
      read(copyunit, '(a1000)') line
      compare = adjustl(line)
      if ((compare == 'C').or.(compare == 'c').or.&
          (compare == '*').or.(compare == '#')) then
        ! Non comment found, back one line and exit
        backspace(copyunit)
        exit
      else
        write(newunit,'(a)') trim(line)
      end if
    end do
  end if

end subroutine

!-----------------------------------------------------------------------------!

end module m_write_fmodel