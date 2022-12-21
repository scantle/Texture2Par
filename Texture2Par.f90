program Texture2Par
  use MakePar
  use fpath
  use IWFM
  use Krige
  use PilotPoints
  use errorhandle
  use pddata
  use str_utls
  use MFSupport
  implicit none

!-----------------------------------------------------------------------------!
!                             Texture2Par
!                         ********************
!                               V 1.1.0
!
! Converts well percent course with depth data to aquifer parameters
! at specified nodes. Intended to work with input/output files of the
! Integrated Water Flow Model (IWFM) with added support for MODFLOW 2000/2005.
! Written by Leland Scantlebury and Marinko Karanovic of S.S. Papadopulos &
! Associates, based on original method and VBA program created by Timothy J.
! Durbin.
!
!
! Requires Input File: Texture2Par.in
!
!-----------------------------------------------------------------------------!

    integer                    :: ioerr, ilay, ifail, ipp, &
                                  i, j, ierr, iarr, twell, wpoint, zone, &
                                  modflowflag, zoneflag, geoflag, &
                                  vi, vj, ind, maxint, flen
    real                       :: a_hmin
    character(150)             :: jnk,fpath_line
    character(256)             :: cerr
    character(1000)            :: line
    logical                    :: writenodefiles

! Intel-Fortran specific way of writing to same command line line
!    open(6, carriagecontrol='fortran')  ! Uses removed by VB

    ! Write Flag
    write(*,'(/,8x,a)') '   Texture2Par v1.1.0'
    write(*,'(8x,a)')   '------------------------'

    ! Read input file
    ierr = 0
    writeNodeFiles = .false.
    write(*,'(/,2x,a)') 'Reading input file'
    open(10, file='Texture2Par.in', status='old', action='read', &
         iostat = ierr,  iomsg=cerr)
    call iomsghandler(ierr, cerr)

    ! Read in model type and well file
    call ReadToSymbol(10, '*')
    read(10,*) line
    read(10,*) wellfile
    read(10,*) geozonefile
    read(10,*) ndata
    allocate(cokrigefiles(ndata-1), codata(ndata-1), knnear(ndata))
    do i=1, (ndata-1)
      read(10,*) cokrigefiles(i)
    end do
    read(10,*) ntclasses
    
    ! Cokriging with classes check (unimplemented)
    if (ndata>1.and.ntclasses>1) then
      write(*,*) 'Error - Currently only cokriging OR multiple datasets supported.'
      stop
    end if
    call ReadToSymbol(10, '*')

    ! Check if MODFLOW or IWFM
    if (index(to_upper(trim(line)),'IWFM')>0) then
      modflowflag = 0
    else if (index(to_upper(trim(line)), 'MODFLOW')>0) then
      modflowflag = 1
    else
      write(*,'(a)') 'Error - Invalid Model Type. Stopping.'
      write(*,'(a)') 'Valid Model Types: IWFM, MODFLOW'
      stop
    end if

    ! Different order/different files for each model
    if (modflowflag==0) then
      ! IWFM
      read(10,*) simFile
      read(10,*) fpath_line    !preSimFile
      call fpath_strip(fpath_line, prefolder, preSimFile)
      read(10,*) gwtemp
      read(10,*) zonefile
      ! MODFLOW Confining flag, set to 0 when IWFM
      ncbd = 0
    else if (modflowflag==1) then
      ! MODFLOW
      read(10,*) simFile
      read(10,*) layfiletemp
      read(10,*) zonefile
      read(10,*) xoff
      read(10,*) yoff
      read(10,*) mfrot
    end if

    ! Node output setting
    call ReadToSymbol(10, '*')
    read(10,*) line
    if (index(to_upper(trim(line)), 'TRUE')>0) writeNodeFiles = .true.

    ! Variograms
    call ReadToSymbol(10, '*')
    read(10,*) line
		all_data_flag = 0  ! Assume nearest neighbors (NN)
		! But check for ALL
		if (index(to_upper(trim(line)), 'ALL')>0) all_data_flag = 1
		! TODO: What about improper entries??
		
		if (all_data_flag==0) then
		  ! For NN we need the NN values all data sets, as well as the max secondary data distance
      do i=1, ndata
        read(10,*) knnear(i)
      end do
			! Only need max secondary if secondary data exists
			if (ndata > 1) then
        read(10,*) max_secondary_dist
			end if
		end if
		
    read(10,*) maxstruct
    
    ! Variogram-related calcs and allocation
    if (ntclasses==1) then             ! NOT using multiple data classes
      nvario = (ndata*(ndata+1)/2) + 1 ! Calculate number of variograms (+1 for pps)
    else
      nvario = ntclasses +1            ! Multiple tex classes = 1 vario per class +1 for pps
    end if
    call allocate_krige()              ! requires nvario & maxstruct
    
    ! Read variograms
    do i=1, nvario
      read(10,*) vi, vj
      if (vi==0) then
        ind=nvario
      else if (vi==vj) then
        ind=vi
      else
        ind=ndata + (vi-1) + (vj-1)
      end if
      read(10,*) nstruct(ind), vario(1,ind)%nugget
      do j=1, nstruct(ind)
        read(10, *) vario(j,ind)%itype, vario(j,ind)%sill, vario(j,ind)%ang1
        read(10, *) vario(j,ind)%a_hmax, vario(j,ind)%a_hmin
      end do
    end do
    ! Process/setup
    call init_krige()

    ! Sediment Settings
    call ReadToSymbol(10, '*')
    if (ntclasses==1) then
      read(10,*) KCk
      read(10,*) KFk
    else
      allocate(KkCvals(ntclasses))
      do i=1, ntclasses
        read(10,*) KkCvals(i)
      end do
    end if
    read(10,*) KHp
    read(10,*) KVp
    read(10,*) Stp
    ! Handle p==0 - power-law function is discountinuous at 0
    if (KHp == 0.0) KHp = tiny(KHp)
    if (KVp == 0.0) KVp = tiny(KVp)
    
    ! Aquifer Pilot Point section
    ! TODO: Add descriptive error when not all pilot point parameters are present
    call ReadToSymbol(10, '*')
    call FindSectionLength(10, ppn, '*')
    if (ntclasses==1) then
      call AquiferPPAllocate()
      do ipp = 1, ppn
        read(10,*) ppx(ipp), ppy(ipp), ppKCMin(ipp), ppdKC(ipp), ppKFMin(ipp), &
                   ppdKF(ipp), ppSsC(ipp), ppSsF(ipp), ppSyC(ipp),ppSyF(ipp), &
                   ppAnisoC(ipp), ppAnisoF(ipp), ppZone(ipp)
      end do

      ! Aquitard Pilot Point section
      !read(10,*)    !call ReadToSymbol(10, '*')
      call ReadToSymbol(10, '*')
      call FindSectionLength(10, atppn, '*')
      call AquitardPPAllocate()
      do ipp = 1, atppn
        read(10,*) ppatx(ipp), ppaty(ipp), ppatKCMin(ipp), ppatdKC(ipp), &
                   ppatKFMin(ipp), ppatdKF(ipp), ppatAnisoC(ipp), &
                   ppatAnisoF(ipp), ppatZone(ipp)
      end do
      close(10)
    else
      ! Multiple texture classes on each pilot point property
      call AquiferPPAllocate_multiclass(ntclasses)
      do ipp=1, ppn
        read(10,*) ppx(ipp), ppy(ipp), (ppKMinCvals(ipp,i), ppdKCvals(ipp,i), i=1, ntclasses), ppSsCvals(ipp,:), ppSyCvals(ipp,:), ppAnisoCvals(ipp,:), ppZone(ipp)
      end do
    end if

    ! MODEL SPECIFIC
    if (modflowflag == 0) then
      write(*,'(2x,a)') 'Reading IWFM Input Files'
      call readIWFM()
      call create_nodetree()
    else if (modflowflag == 1) then
      write(*,'(2x,a)') 'Reading MODFLOW Input Files'
      call read_modflow(simfile)
      call calcCellCenters_global()
      call create_nodetree()
    end if

    ! Pre-read Well file for file length, nwells, maxWPoints
    open(10, file = trim(wellfile), status='old', iostat = ierr,  iomsg=cerr)
    call iomsghandler(ierr, cerr)
    read(10,*)                                    ! Header
    ioerr = 0
    wfilelen = 0
    nwell = 0
    maxWPoints = 0
    do while (ioerr == 0)
      read(10,*, iostat = ioerr) jnk, twell, wpoint
      wfilelen = wfilelen + 1
      if (twell > nwell) nwell = twell
      if (wpoint > maxWPoints) maxWPoints = wpoint
    end do
    wfilelen = wfilelen - 1
    close(10)

    ! Read Node Zone file, unless zonefile was set to NONE
    allocate(pp_lay2arr(nlayers+natlayers))
    if (to_upper(trim(zonefile))=='NONE') then
      zoneflag = 0
    else
      ! External file
      zoneflag = 1
      if (modflowflag == 0) then
        call inputZones_IWFM()
      else
        call inputZones_modflow()
      end if
    end if

    ! Handle No Pilot Point Zones (same as one zone)
    if (zoneflag == 0) then
      call no_ppzone_setup(nnodes, nlayers, natlayers)
    end if

    ! Allocate now that many dimensions are known
    call WellAllocate()
    call NodeAllocate()

    ! Read geologic units file
    if (to_upper(trim(geozonefile))=='NONE') then
      ! No geologic zones
      ! perhaps a different subroutine should be used instead of dummy values?
      geoflag = 0
      ngeozones = 1
      ngeozones_lay(1:nlayers) = 1
      geozones(1:nnodes,1:nlayers) = '[NONE]'
      geonames_lay(1,1:nlayers) = '[NONE]'
    else
      geoflag = 1
      if (modflowflag == 0) then
        call readNodeGeoZones_IWFM()
      else if (modflowflag == 1) then
        call readNodeGeoZones_modflow()
      end if
    end if

    ! Write Read Summary
    write(*,'(2x,a)') 'Setup Summary:                          '
    write(*,'(4x,1a36,i8)') 'Number of aquifer pilot points: ', ppn
    write(*,'(4x,1a36,i8)') 'Number of aquitard pilot points: ', atppn
    write(*, '(4x,1a36,i8)') 'Number of Layers: ', nlayers
    if (modflowflag == 0) then
      write(*, '(4x,1a36,i8)') 'Number of Nodes: ', nnodes
      write(*, '(4x,1a36,i8)') 'Number of Elements: ', nelements
    else if (modflowflag == 1) then
      ! Perhaps "cells per layer"
      write(*, '(4x,1a36,i8)') 'Number of Cells: ', nnodes
      write(*, '(4x,1a36,i8)') 'Quasi-3D Confining Beds: ', ncbd
    end if
    !write(*,'(4x,1a36,i)') 'Well File Length: ', wfilelen
    write(*,'(4x,1a36,i8)') 'Number of Texture Classes: ', ntclasses
    write(*,'(4x,1a36,i8)') 'Number of Well Logs: ', nwell
    write(*,'(4x,1a36,i8)') 'Number of Pilot Point Zones: ', maxval(pp_arr_nzones)
    write(*,'(4x,1a36,i8,/)') 'Number of Geologic Units: ', ngeozones

! MODEL SPECIFIC
! Call input wells subroutine
    if (modflowflag==0) then
      call inputwells_IWFM(geoflag)
    else if (modflowflag==1) then
      call inputwells_modflow(geoflag)
    end if
    
! Handle additional datasets
    do i=1, (ndata-1)
      call pre_read_pddata(11, cokrigefiles(i), flen, codata(i)%npdat, maxint)
      call allocate_pddataset(codata(i), maxint, nlayers)
      call read_pddata_IWFM(codata(i), 11, cokrigefiles(i), flen, geoflag)
    end do

! Call Layer Elevation subroutine
    if (modflowflag == 0) then
      write(*,'(2x,a)') 'Calculating Layer Elevations at Data Points'
      call layerelevation(TopElev, BotElev, wellelevtop, wellelevbot)
      ! For aquitard
      call layerelevation(AqTardTopElev, AqTardBotElev, &
                          wellaqtardelevtop, wellaqtardelevbot)
      call pd_elev_at_points(ndata-1, TopElev, BotElev, AqTardTopElev, AqTardBotElev)
    end if
    ! Modflow Well elevations are assigned during well read

! Move co-located wells (if kriging still crashes, well locations need to be revised)
    call colocate_corrector(xwell, ywell, nwell)

! Call Layer Texture subroutine
    write(*,'(2x,a)') 'Calculating Percent Coarse at Wells'
    do i=1, ntclasses
      call layertexture(wellelevtop, wellelevbot, i, PcWellElem)
      ! TODO find a more reasonable place for this
      wellmean = sum(PcWellElem(:,:,i),PcWellElem(:,:,i)>0.0)/count(PcWellElem(:,:,i)>0.0)
      ! If IWFM or if there are MODFLOW confining beds
      if ((modflowflag == 0).or.(ncbd > 0 )) then
        call layertexture(wellaqtardelevtop, wellaqtardelevbot, 1, PcWellAqtardElem)
      end if
    end do
      
    ! Average pddatasets by layer
    do i=1, (ndata-1)
      call pd_avg2lay(codata(i), codata(i)%npdat, codata(i)%laytop, codata(i)%laybot, codata(i)%layv)
      ! TODO find a more reasonable place for this (all-layer average)
      codata(i)%v_mean = sum(codata(i)%layv,codata(i)%layv>NODATA)/count(codata(i)%layv>NODATA)
      if ((modflowflag == 0).or.(ncbd > 0 )) then
        call pd_avg2lay(codata(i), codata(i)%npdat, codata(i)%aqttop, codata(i)%aqtbot, codata(i)%aqtv)
      end if
    end do

    ! Krige from pilot points to nodes
    ! Aquifer Parameters - KCMin, DeltaKC, KFMin, DeltaKF, SsC, SsF, SyC, Aniso
    if (ntclasses==1) then
      write(*,'(2x,a,i3)') 'Kriging Aquifer Pilot Point Parameters'
      call ppaq_parkrige_arr(nnodes, nodex, nodey, ppKCMin,  nvario, nKCMin )
      call ppaq_parkrige_arr(nnodes, nodex, nodey, ppdKC,    nvario, dKC    )
      call ppaq_parkrige_arr(nnodes, nodex, nodey, ppKFMin,  nvario, nKFMin )
      call ppaq_parkrige_arr(nnodes, nodex, nodey, ppdKF,    nvario, dKF    )
      call ppaq_parkrige_arr(nnodes, nodex, nodey, ppSsC,    nvario, nSsC   )
      call ppaq_parkrige_arr(nnodes, nodex, nodey, ppSsF,    nvario, nSsF   )
      call ppaq_parkrige_arr(nnodes, nodex, nodey, ppSyC,    nvario, nSyC   )
      call ppaq_parkrige_arr(nnodes, nodex, nodey, ppSyF,    nvario, nSyF   )
      call ppaq_parkrige_arr(nnodes, nodex, nodey, ppAnisoC, nvario, nAnisoC)
      call ppaq_parkrige_arr(nnodes, nodex, nodey, ppAnisoF, nvario, nAnisoF)

      ! Krige Aquitard parameters from pilot points to nodes
      if ((modflowflag == 0).or.(ncbd > 0 )) then
        if (atppn > 0) then
          ! Aquitard Parameters - KCMin, DeltaKC, KFMin, DeltaKF, Aniso
          write(*,'(2x,a,i3)') 'Kriging Aquitard Pilot Point Parameters'
          call ppat_parkrige_arr(nnodes, nodex, nodey, ppatKCMin,  nvario, nKCMin )
          call ppat_parkrige_arr(nnodes, nodex, nodey, ppatdKC,    nvario, dKC    )
          call ppat_parkrige_arr(nnodes, nodex, nodey, ppatKFMin,  nvario, nKFMin )
          call ppat_parkrige_arr(nnodes, nodex, nodey, ppatdKF,    nvario, dKF    )
          call ppat_parkrige_arr(nnodes, nodex, nodey, ppatAnisoC, nvario, nAnisoC)
          call ppat_parkrige_arr(nnodes, nodex, nodey, ppatAnisoF, nvario, nAnisoF)
        else
          write(*,'(2x,a,/)') 'No Aquitard Pilot Points - Skipping       '
          nKCMin(:,ppaq_narr+1:pp_narr) = -999
          dKC   (:,ppaq_narr+1:pp_narr) = -999
          nKFMin(:,ppaq_narr+1:pp_narr) = -999
          dKF   (:,ppaq_narr+1:pp_narr) = -999
        end if
      end if
    else
      ! Multiple Classes
      do i=1, ntclasses
        call ppaq_parkrige_arr(nnodes, nodex, nodey, ppKMinCvals(:,i),   nvario,   KMinCvals(:,:,i))
        call ppaq_parkrige_arr(nnodes, nodex, nodey, ppdKCvals(:,i),     nvario,     dKCvals(:,:,i))
        call ppaq_parkrige_arr(nnodes, nodex, nodey, ppSsCvals(:,i),     nvario,     SsCvals(:,:,i))
        call ppaq_parkrige_arr(nnodes, nodex, nodey, ppSyCVals(:,i),     nvario,     SyCVals(:,:,i))
        call ppaq_parkrige_arr(nnodes, nodex, nodey, ppAnisoCvals(:,i),  nvario, nAnisoCvals(:,:,i))
      end do
    end if
      
! Call interpolate subroutines (Well PC to Node Aquifer parameter)
    write(*,'(2x,a)') 'Calculating Node Aquifer Parameters           '
    if (ntclasses==1) then
      call interpolate()
      if ((modflowflag == 0).or.(ncbd > 0 )) then
        call interpolate_aqtard()
      end if
    else
      call interpolate_multiclass()
      ! AQUITARD INTERPOLATE UNIMPLEMENTED FOR MULTICLASS
      KvB_aqtard = 0.0
    end if

! Write Output files, routines called depend on model type
    if (modflowflag == 0) then
      if (writeNodeFiles) then
        ! Only calculate if writing the file
        call calc_interbed_thickness()
        write(*,'(2x,a)') 'Writing Node Parameter Files'
        do i=1, ntclasses
          write(line,'(a,i1,a)') 't2p_WellPC_class', i, '.out'
          call writeWellFile_IWFM(PcWellElem(:,:,i), line)
        end do
        call writeWellFile_IWFM(PcWellAqtardElem, 't2p_WellAqTardPC.out')
        ! Future output selection
        !call writeWellFile_IWFM(TopElev, 't2p_WellTop.out')
        !call writeWellFile_IWFM(BotElev, 't2p_WellBot.out')
        if (ntclasses==1) then
          call writeNodeFile(PcNode, 't2p_NodePC.out')
        else if (ntclasses > 1) then
          do i=1, ntclasses
            write(line,'(a,i1,a)') 't2p_NodePc_class', i, '.out'
            call writeNodeFile(nodecval(:,:,i), line)
          end do
        end if
        call writeNodeFile(KhB,    't2p_KhB.out'   )
        call writeNodeFile(KvB,    't2p_KvB.out'   )
        call writeNodeFile(SsB,    't2p_SsB.out'   )
        call writeNodeFile(SyB,    't2p_SyB.out'   )
        call writeNodeFile(interbed_thick,'t2p_InterbedThickness.out')
        call writeNodeFile(PcNode_aqtard, 't2p_NodeAqTardPC.out')
        call writeNodeFile(KvB_aqtard,    't2p_AqTardKvB.out'   )
      else
        write(*,'(2x,a)') 'WriteNodeFiles set to FALSE - no parameter files written'
      end if
        call writeIWFMgwfile()

    else if (modflowflag == 1) then
      if (writeNodeFiles) then
        ! Only calculate if writing the file
        call calc_interbed_thickness()
        write(*,'(2x,a)') 'Writing Cell Parameter Files'
        call writeWellFile_modflow(PcWellElem,       't2p_WellPC.out')
        call writeCellFile(PcNode, 't2p_NodePC.out')
        call writeCellFile(KhB,    't2p_KhB.out'   )
        call writeCellFile(KvB,    't2p_KvB.out'   )
        call writeCellFile(SsB,    't2p_SSB.out'   )
        call writeCellFile(SyB,    't2p_SyB.out'   )
        call writeCellFile(interbed_thick,'t2p_InterbedThickness.out')
        if (ncbd > 0 ) then
          call writeWellFile_modflow(PcWellAqtardElem, 't2p_WellAqTardPC.out')
          call writeCellFile(PcNode_aqtard, 't2p_AqTardNodePC.out')
          call writeCellFile(KvB_aqtard,    't2p_AqTardKvB.out'   )
        end if
      else
        write(*,'(2x,a)') 'WriteCellFiles set to FALSE - no parameter files written'
      end if
      call readwriteFlowPackage()
      call clear_modflow_memory()
    end if
    
    !TODO Poorly Implemented
    do i=1, (ndata-1)
      write(jnk,'(a,i1,a)') 't2p_cokrige_',i,'.out'
      call write_pd_file(codata(i), jnk, nlayers)
    end do

    write(*,'(/,2x,a)') 'Program finished!'

  end program Texture2Par
!-----------------------------------------------------------------------------!