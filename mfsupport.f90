module MFSupport
  implicit none
  
! Added in v1.3.1
  
  SAVE
  real                      :: xoff, yoff, mfrot
  real,allocatable          :: mftop(:,:,:), mfbot(:,:,:), locx(:), locy(:)
  integer                   :: mfmaxunit, layfile, ncbd
  character(3)              :: gwftyp
  character(10)             :: mfzonename
  character(200)            :: layfiletemp
  logical                   :: lop
  
  CONTAINS

!-----------------------------------------------------------------------------!
  subroutine read_modflow(namfile)
    
! MODFLOW globals to import
    use GLOBAL, ONLY : IUNSTR,NLAY,NROW,NCOL,DELR,DELC,TOP,BOT,LAYCBD,IUNIT
    use PARAMMODULE, ONLY: NMLTAR, NZONAR
    use GWFBCFMODULE, ONLY: VKA, SC1, SC2
    
! Texture2Par modules
    use MakePar
    use errorhandle
    
    implicit none
    
! Local variables
    character(200)         :: namfile
    integer                :: ierr, i, j, k
    
! Open modflow files
    call init_modflow_usgs(namfile,ierr,mfmaxunit)
    call iostathandler(ierr)
    if (IUNSTR==0) then
!      write(*,*) 'Structured Model'
    else if (IUNSTR==1) then
      call USGunimplemented()
    end if

! Set MakePar module values
    nnodes = nrow * ncol
    nlayers = nlay
    nzones = NZONAR
    
! Handle MODFLOW settings
    ! UPW or LPF (or none)
    layfile = 0
    if (IUNIT(23) > 0) then
      layfile = 23
      gwftyp = 'LPF'
    end if
    if (IUNIT(45) > 0) then
      layfile = 45
      gwftyp = 'UPW'
    end if
    if (layfile==0) call NotInNam('LPF or UPW')
    
! Allocate
    call allocate_modflow()
    
! Get Top & Bot arrays into node format
    ncbd = 0
    Elevation(1:nnodes) = top(1:nnodes)
    do k=1, nlay
      topelev(1:nnodes,k) = top((k-1)*nnodes+1 : k*nnodes)
      botelev(1:nnodes,k) = bot((k-1)*nnodes+1 : k*nnodes)
      ! Quasi-3d confining beds == Aquitards
      if (laycbd(k) > 0) then
        ncbd = ncbd + 1
        aqtardtopelev(1:nnodes,k) = bot((k-1)*nnodes+1 : k*nnodes)
        aqtardbotelev(1:nnodes,k) = top(k*nnodes+1 : (k+1)*nnodes)
      else
        ! Fill with junk data
        aqtardtopelev(1:nnodes,k) = NODATA
        aqtardbotelev(1:nnodes,k) = NODATA
      end if
    end do
    
  end subroutine read_modflow
!-----------------------------------------------------------------------------!
  
!-----------------------------------------------------------------------------!
  subroutine calcCellCenters_global()
    use GLOBAL, ONLY : NROW, NCOL, DELR, DELC
    use MakePar
    ! Calculates cell centers in global coordinates
    
    integer                :: i, j, n
    real                   :: xbase, ybase
    
    ybase = sum(delc)
    do i=1, nrow
      xbase = 0.0
      do j=1, ncol
        n = (i-1)*ncol+j
        locx(n) = xbase + 0.5 * delr(j)
        locy(n) = ybase - 0.5 * delc(i)
        call loc2glo(locx(n), locy(n), xoff, yoff, mfrot, nodex(n), nodey(n))
        xbase = xbase + delr(j)
      end do
      ybase = ybase - delc(i)
    end do
    
  end subroutine calcCellCenters_global  
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
SUBROUTINE inputwells_modflow()
    use MakePar
    use kd_tree
    implicit none

!-----------------------------------------------------------------------------!
! Inputs well information.
! Uses NodeTree to locate closest nodes and pass this information to
! subroutine that locates cell/element of well.
!-----------------------------------------------------------------------------!
    
    integer, parameter       :: nnear=4
    integer                  :: i, iwell, izone, tpc, currwellelem, node
    integer                  :: closenodes(nnear)
    real                     :: depth, tx, ty, tz, td
    real                     :: nodedists(nnear), queryxy(2)
    character(30)            :: wellname

! Read well file
    write(*,'(2x,2a)') 'Reading Well File ', trim(wellfile)
    open(10, file =trim(wellfile), status='old', action='read')
    read(10,*)                        ! Header
    do i=1, wfilelen
      read(10,*) wellname, iwell, izone, tpc, tx, ty, tz, td
      ! Well "zone" is discrete depth interval and should not be confused with
      ! other uses of zone (e.g. pilot points)
      PcWellZone(iwell, izone) = tpc
      numZone(iwell) = izone
      if (izone == 1) then
        Xwell(iwell) = tx
        Ywell(iwell) = ty
        Zland(iwell) = tz
        ! Use node tree to find closest node
        queryxy(1) = tx
        queryxy(2) = ty
        call n_nearest_to(nodetree, queryxy, nnear, closenodes, nodedists)
        ! Handle Model Differences
        
        ! MODFLOW
        node = get_well_cell(tx, ty, closenodes, nnear)
        wellelem(iwell) = node
        ! While looping over wells, might as well assign modflow
        ! elevations to wells, as no interpolation is required
        wellelevtop(iwell,1:nlayers) = topelev(node,1:nlayers)
        wellelevbot(iwell,1:nlayers) = botelev(node,1:nlayers)
        if (ncbd > 0) then
          wellaqtardelevtop(iwell,1:nlayers) = aqtardtopelev(node,1:nlayers)
          wellaqtardelevbot(iwell,1:nlayers) = aqtardbotelev(node,1:nlayers)
        end if
        
      end if
      Zwell(iwell, izone) = tz - td
    end do
    close(10)
    
    !write(*,'(4x,a)') 'Well file read successfully.'

end subroutine inputwells_modflow
!-----------------------------------------------------------------------------!
  
!-----------------------------------------------------------------------------!
  subroutine inputZones_modflow()
    use GLOBAL, ONLY : IUNSTR, NLAY, NROW, NCOL, DELR, DELC, TOP, BOT, NODES
    use PARAMMODULE, ONLY: NZONAR, IZON, ZONNAM
    use MakePar
    !
    
    integer                :: i, j, k, zone, namindex
    logical                :: newzone
    
    ! Use MODFLOW Zone array specified in input file to fill zones and
    ! nodezones arrays defined in MakePar
    nzones = 0
    
    do i=1, nzonar
      if (trim(ZONNAM(i))==trim(mfzonename)) then
        namindex = i
      end if
    end do
    
    ! IZON is (NCOL,NROW,NZONAR))
    do i=1, nrow
      do j=1, ncol
        zone = IZON(j,i,namindex)
        ! Is this a new zone?
        newzone = .true.
        do k=1, nzones
          if (zone == zones(k)) newzone = .false.
        end do
        ! Add to Array
        nodezone((i-1)*ncol+j) = zone
        ! Add if new
        if (newzone) then
          nzones = nzones + 1
          zones(nzones) = zone
        end if
      end do
    end do

    
  end subroutine inputZones_modflow  
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
integer function get_well_cell(x, y, closenodes, nnear) result(cell)
  implicit none
  
  real, intent(in)         :: x, y
  integer, intent(in)      :: nnear, closenodes(nnear)
  integer                  :: i
  real                     :: xl, yl, xmin, xmax, ymin, ymax
  
  cell = 0
  ! Check each cell center to see if well belongs
  do i=1, nnear
    call getCellBounds_local(closenodes(i), xmin, xmax, ymin, ymax)
    call glo2loc(x, y, xoff, yoff, mfrot, xl, yl)
    if ((xl >= xmin).and.(xl < xmax).and.(yl >= ymin).and.(yl < ymax)) then
      cell = closenodes(i)
      exit
    end if
  end do
  
  if (cell == 0) then
    ! Maybe this should just be a warning?
    write(*,*) 'ERROR well not found in a MODFLOW cell'
    write(*,*) x, y
    stop
  end if

end function get_well_cell
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
  subroutine getCellBounds_local(node, xmin, xmax, ymin, ymax)
    use MakePar
    use GLOBAL, ONLY : DELR, DELC
    
    integer,intent(in)     :: node
    real,intent(out)       :: xmin,xmax,ymin,ymax
    integer                :: row, col
    
    call node2rowcol(node, row, col)
    
    xmin = locx(node) - 0.5 * delr(col)
    xmax = locx(node) + 0.5 * delr(col)
    ymin = locy(node) - 0.5 * delc(row)
    ymax = locy(node) + 0.5 * delc(row)
    
  end subroutine getCellBounds_local  
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
  subroutine getCellBounds_global(node, xmin, xmax, ymin, ymax)
    use MakePar
    use GLOBAL, ONLY : DELR, DELC
    
    ! Not compatible with rotated grids - not currently used in Texture2Par
    
    integer,intent(in)     :: node
    real,intent(out)       :: xmin,xmax,ymin,ymax
    integer                :: row, col
    
    call node2rowcol(node, row, col)
    
    xmin = nodex(node) - 0.5 * delr(col)
    xmax = nodex(node) + 0.5 * delr(col)
    ymin = nodey(node) - 0.5 * delc(row)
    ymax = nodey(node) + 0.5 * delc(row)
    
  end subroutine getCellBounds_global  
!-----------------------------------------------------------------------------!
  
!-----------------------------------------------------------------------------!
  subroutine node2rowcol(node, row, col)
    use GLOBAL, ONLY : NROW, NCOL
    
    integer,intent(in)     :: node
    integer,intent(out)    :: row, col
    
    row = CEILING(real(node)/real(ncol))
    col = node - (row-1)*ncol
    
  end subroutine node2rowcol  
!-----------------------------------------------------------------------------!
  
!-----------------------------------------------------------------------------!
  subroutine allocate_modflow()
    use GLOBAL, ONLY : IUNSTR, NLAY, NROW, NCOL, DELR, DELC, TOP, BOT, NODES
    use MakePar
    implicit none
    
    allocate(locx(nnodes), &
             locy(nnodes)  )
    
    ! May be better to put into a MakePar routine for allocation
    ! Since these are all use for both models
    allocate(nodex(nnodes),              &
             nodey(nnodes),              &
             nodexy(nnodes,2),           &
             Elevation(nnodes),          &
             topelev(nnodes,nlay),       &
             botelev(nnodes,nlay),       &
             aqtardtopelev(nnodes,nlay), &
             aqtardbotelev(nnodes,nlay)  )

  end subroutine allocate_modflow
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine writeWellFile_modflow(v, filename)
  use MakePar
  implicit none

  real,intent(in)          :: v(nwell,nlayers)
  character(*),intent(in)  :: filename
  integer                  :: iwell, ilay, row, col
  character(60)            :: fmt(2)
  
  write(6,'("+",4x,a,a35)') 'Writing:', filename

  write(fmt(1),'("(3a10,",i5,"i14)")') nlayers
  write(fmt(2),'("(3i10,",i5,"es14.5)")') nlayers

  open(20, file=trim(filename), status='replace', action='write', buffered='YES')
  ! Header
  write(20,fmt(1)) 'Well', 'Row', 'Col', (ilay, ilay=1, nlayers)
  ! Values
  do iwell = 1, nwell
    call node2rowcol(wellelem(iwell), row, col)
    write(20,fmt(2)) iwell, row, col, (v(iwell,ilay), ilay=1, nlayers)
  end do
  
  close(20)

end subroutine writeWellFile_modflow
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine writeCellFile(v, filename)
  use MakePar
  implicit none

  real,intent(in)          :: v(nlayers,nnodes)
  character(*),intent(in)  :: filename
  integer                  :: inode, ilay, row, col
  character(60)            :: fmt(2)
  
  write(6,'("+",4x,a,a35)') 'Writing :', filename

  write(fmt(1),'("(3a8,2a18,",i5,"i14)")') nlayers
  write(fmt(2),'("(3i8,2f18.5,",i5,"e14.5)")') nlayers

  open(20, file=trim(filename), status='replace', action='write', buffered='YES')
  ! Header
  write(20,fmt(1)) 'Row', 'Col', 'Zone', 'X', 'Y', (ilay, ilay=1, nlayers)
  ! Values
  do inode = 1, nnodes
    call node2rowcol(inode, row, col)
    write(20,fmt(2)) row, col, nodezone(inode), nodex(inode), nodey(inode), &
                     (v(ilay, inode), ilay=1,nlayers)
  end do
  
  close(20)

end subroutine writeCellFile
!-----------------------------------------------------------------------------!
  
!-----------------------------------------------------------------------------!
  subroutine readwriteFlowPackage()
    use GLOBAL, ONLY:IUNIT
    use errorhandle
    implicit none
    
    integer                :: ierr, upwflag
    character(300)         :: fname
    logical                :: lop
    
    write(*,'(2x,2(a,3a))') 'Reading temp ', gwftyp, ', writing new ', gwftyp
    
    ! If LPF/UPW is open, close it. We need to do a bait-and-switch on MODFLOW
    inquire(unit=IUNIT(layfile), opened=lop, name=fname)
    if (lop) close(IUNIT(layfile))
    
    ! Initialize MODFLOW LPF values with template LPF file (user provided)
    open(IUNIT(layfile),file=trim(layfiletemp),iostat=ierr,status='old')
    call iostathandler(ierr,layfiletemp)
    upwflag = 0
    if (gwftyp=='UPW') upwflag = 1
    CALL GWF2BCFU1AR(IUNIT(1),IUNIT(22),IUNIT(layfile),IUNIT(64), upwflag)
    close(IUNIT(layfile))
    
    ! Re-open actual LPF file and OVERWRITE
    open(IUNIT(layfile),file=trim(fname),iostat=ierr,status='replace')
    call iostathandler(ierr,fname)
    call writeGWFP(IUNIT(layfile))
    close(IUNIT(layfile))
    
  end subroutine readwriteFlowPackage
!-----------------------------------------------------------------------------!

end module MFSupport