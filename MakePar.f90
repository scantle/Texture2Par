MODULE MakePar
  use kdtree2_module
  use pddata
  implicit none

!-----------------------------------------------------------------------------!
! MakePar Module for Texture2Par
! Some subroutines were originally VBA subroutines written by Tim Durbin
! Author: Leland Scantlebury of S.S. Papadopulos & Associates
!-----------------------------------------------------------------------------!

    SAVE
    ! Parameters
    integer, parameter         :: max_node_elem_members = 15, &                 ! Maximum number of elements a node is a part of (will stop w/ error if too low)
                                  max_zones = 100, &                            ! Maximum number of zones (used for both pilot point and geologic)
                                  USE_IDW = 0, &
                                  max_nodeid=5E6                                ! Maximum nodeid value (to prevent enormous array allocation)
    real, parameter            :: NODATA=-999.0
                               
    real,allocatable           :: wellelevbot(:,:), ssb(:,:), pcwellelem(:,:,:), &
                                  xwell(:), wellaqtardelevtop(:,:), &
                                  wellelevtop(:,:),  &
                                  ywell(:), &
                                  zland(:), aqtardbotelev(:,:), &
                                  nodex(:), nodey(:), &
                                  zwell(:,:), khb(:,:), &
                                  pcwellaqtardelem(:,:,:), &
                                  aqtardtopelev(:,:),  kvb(:,:), &
                                  syb(:,:), elevation(:), botelev(:,:), &
                                  wellaqtardelevbot(:,:), topelev(:,:), &
                                  pcnode(:,:), pcnode_aqtard(:,:), &
                                  kvb_aqtard(:,:), pcwellzone(:,:,:), &
                                  Thick(:,:), AqTardThick(:,:), &
                                  nodexy(:,:), &
                                  interbed_thick(:,:), &
                                  nodecval       (:,:,:), &
                                  nodecval_aqtard(:,:,:)
    integer,allocatable        :: wellelem(:), numzone(:), &
                                  elements(:,:), nodeelem(:,:), nnodeelem(:), &
                                  ngeozones_lay(:), nodeid(:), nodebyid(:), &
                                  knnear(:)
    integer                    :: nlayers, nnodes, nwell, wfilelen, &
                                  nelements, maxWPoints, &
                                  ngeozones, natlayers, all_data_flag, ntclasses
    real                       :: kfk, stp, khp, kvp, kck, wellmean
    character(30),allocatable  :: geozones(:,:), wellgeozones(:,:), &
                                  geonames_lay(:,:)
    character(30)              :: geonames(max_zones)
    character(60)              :: outfile
    character(150)             :: wellfile, geozonefile
    character(150),allocatable :: cokrigefiles(:)
    logical                    :: colocated=.false.
                               
    ! Pilot Point Parameters   
    real,allocatable           :: nKCMin(:,:), dKC(:,:), nKFMin(:,:),&
                                  dKF(:,:), nSsC(:,:), nSsF(:,:),    &
                                  nSyC(:,:), nSyF(:,:), nAnisoC(:,:),&
                                  nAnisoF(:,:)

    ! Cokriging data
    type(pddataset),allocatable :: codata(:)
    
    ! Alternatively, multi-class (texture class) data
    type(pddataset)             :: tcdata
    real,allocatable            :: KkCvals(:),KMinCvals(:,:,:), &
                                   dKCvals(:,:,:),SsCvals(:,:,:),       &
                                   SyCVals(:,:,:),nAnisoCvals(:,:,:)    
    
    ! Global Trees
    type(kdtree2), pointer :: nodetree

    CONTAINS

!-----------------------------------------------------------------------------!
SUBROUTINE layerelevation(top_elev, bot_elev, result_top, result_bot)
    implicit none
!-----------------------------------------------------------------------------!
! Find layer elevations at each well using node elevation data.
! Inverse Distance Weighting is used to determine the elevations.
!
! Arguments:
!    - topelev    - Top elevation at nodes
!    - botelev    - Bottom elevation at nodes
!    - result_top - Layer top elevations at each well
!    - result_bot - Layer bottom elevations at each well
!-----------------------------------------------------------------------------!
    real,intent(in)          :: top_elev(nnodes, nlayers), &
                                bot_elev(nnodes, nlayers)
    real,intent(inout)       :: result_top(nwell, nlayers), &
                                result_bot(nwell, nlayers)

    result_top = 0.0
    result_bot = 0.0

    call IDW(nwell, xwell, ywell, wellelem, result_top, nnodes, nlayers, &
             nodex, nodey, top_elev, nelements, elements)
    call IDW(nwell, xwell, ywell, wellelem, result_bot, nnodes, nlayers, &
             nodex, nodey, bot_elev, nelements, elements)

end subroutine layerelevation
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
SUBROUTINE pd_elev_at_points(nsets, top_elev, bot_elev, top_tard, bot_tard)
    implicit none
!-----------------------------------------------------------------------------!
! Find layer elevations at each well using node elevation data.
! Inverse Distance Weighting is used to determine the elevations.
!
! Arguments:
!    - top_elev    - Top elevation at nodes
!    - bot_elev    - Bottom elevation at nodes
!    - top_tard    - Top elevation at aquitard layers
!    - bot_tard    - Bottom elevation at aquitard layers
!-----------------------------------------------------------------------------!
    real,intent(in)               :: top_elev(nnodes, nlayers), &
                                     bot_elev(nnodes, nlayers), &
                                     top_tard(nnodes, nlayers), &
                                     bot_tard(nnodes, nlayers)
    integer, intent(in)           :: nsets
    integer                       :: i
    
    do i=1, nsets
      codata(i)%laytop = 0.0
      codata(i)%laybot = 0.0
      call IDW(codata(i)%npdat, codata(i)%xy(1,:), codata(i)%xy(2,:), &
               codata(i)%node, codata(i)%laytop, nnodes, nlayers, nodex, nodey, &
               top_elev, nelements, elements)
      call IDW(codata(i)%npdat, codata(i)%xy(1,:), codata(i)%xy(2,:), &
               codata(i)%node, codata(i)%laybot, nnodes, nlayers, nodex, nodey, &
               bot_elev, nelements, elements)
      call IDW(codata(i)%npdat, codata(i)%xy(1,:), codata(i)%xy(2,:), & 
               codata(i)%node, codata(i)%aqttop, nnodes, nlayers, nodex, nodey, &
               top_tard, nelements, elements)
      call IDW(codata(i)%npdat, codata(i)%xy(1,:), codata(i)%xy(2,:), &
               codata(i)%node, codata(i)%aqtbot, nnodes, nlayers, nodex, nodey, &
               bot_tard, nelements, elements)
    end do

end subroutine pd_elev_at_points
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
SUBROUTINE layertexture(well_topelev, well_botelev, pcwellout)
    implicit none
!-----------------------------------------------------------------------------!
! This subprogram calculates the proportion of coarse-grained material at each
! well within an element layer, where the well at least partially penetrates
! the element layer.
!
! Arguments:
!    - well_topelev - Top elevation at every well
!    - well_botelev - Bottom elevation at every well
!    - pcwellout    - Percent coarse at each well element
!-----------------------------------------------------------------------------!

    real,intent(in)          :: well_topelev(nwell, nlayers), &
                                well_botelev(nwell, nlayers)
    real,intent(inout)       :: pcwellout(nwell, nlayers,1)
    real                     :: elemtop, elembot, zonetop, zonebot, &
                                ztop, zbot, lengthin, sumlengthin, &
                                sumpcwellzone
    integer                  :: i, iwell, ilay, izone, irow
    character(60)            :: fmt(2)

    !write(*,'(4x, a)') 'Calculating coarse-grained % for each well in each layer'
    do ilay = 1, nlayers
      do iwell = 1, nwell
        elemTop = well_topelev(iwell, ilay)
        elemBot = well_botelev(iwell, ilay)
!     Interval within layer
        SumPcWellZone = 0
        SumLengthIn = 0
        do izone = 1, numZone(iwell)
          if (izone == 1) then
            zoneBot = Zwell(iwell, izone)
            zoneTop = Zland(iwell)
          else
            zoneBot = Zwell(iwell, izone)
            zoneTop = Zwell(iwell, izone - 1)
          end if
!       Calculate vertical weight for screen interval within element
          if ((zoneTop >= elemBot).and.(zoneBot <= elemTop).and.&
              (PcWellZone(iwell, izone, 1) >= 0)) then
            ! Find min top/bot values
            Ztop = min(elemTop, zoneTop)
            Zbot = max(elemBot, zoneBot)
            LengthIn = Ztop - Zbot
            SumLengthIn = SumLengthIn + LengthIn
            SumPcWellZone = SumPcWellZone + LengthIn * PcWellZone(iwell, izone, 1)
          end if
        end do
        pcwellout(iwell, ilay, 1) = NODATA
        if (SumLengthIn > 0.0) then
          pcwellout(iwell, ilay, 1) = SumPcWellZone / SumLengthIn
        end if
      end do
    end do

end subroutine layertexture
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
SUBROUTINE pd_avg2lay(pdset, n, topelevs, botelevs, layvarray)
    
  implicit none
!-----------------------------------------------------------------------------!
! This subprogram calculates the proportion of coarse-grained material at each
! well within an element layer, where the well at least partially penetrates
! the element layer.
!
! Arguments:
!-----------------------------------------------------------------------------!

  type(pddataset),intent(in) :: pdset
  integer,intent(in)         :: n
  real,intent(in)            :: topelevs(n, nlayers), &
                                botelevs(n, nlayers)
  real,intent(inout)         :: layvarray(n, nlayers)
  real                       :: intv_top, intv_bot, lengthin, sumlengthin, intv_sum
  integer                    :: i, p, ilay, intv
    
  ! TODO revise with WHERE/FORALL statement(s)
  do ilay = 1, nlayers
  ! Loop over data points
    do i = 1, n
      ! Loop over depth intervals, get top&bot
      intv_sum = 0.0
      SumLengthIn = 0.0
      do intv = 1, pdset%ndint(i)
        if (intv == 1) then 
          intv_bot = pdset%z(intv, i)
          intv_top = pdset%zland(intv)
        else
          intv_bot = pdset%z(intv, i)
          intv_top = pdset%z(intv - 1,i)
        end if
  !       Calculate vertical weight for screen interval within element
        if ((intv_top >= botelevs(i,ilay)).and.&
            (intv_bot <= topelevs(i,ilay)).and.&
            (pdset%v(intv,i) >= 0.0)) then
          ! Find min top/bot values
          LengthIn = min(topelevs(i,ilay), intv_top)-max(botelevs(i,ilay), intv_bot)
          intv_sum = intv_sum + LengthIn * pdset%v(intv,i)
          SumLengthIn = SumLengthIn + LengthIn
        end if
      end do
      layvarray(i, ilay) = NODATA
      if (SumLengthIn > 0.0) then
        layvarray(i, ilay) = intv_sum / SumLengthIn
      end if
    end do
  end do

end subroutine pd_avg2lay
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
SUBROUTINE interpolate()
    use Krige
    use PilotPoints
    use kdtree2_module
    implicit none
!-----------------------------------------------------------------------------!
! Interpolates from wells to nodes 
!
! Well GeoZones Loop added 7/1/2019
! Utilizes arrays geonames_lay and ngeozones_lay which list geologic zones
! ("Units") present in each model layer, based upon geologic node zones
! -LS
!-----------------------------------------------------------------------------!
    real                      :: sum1, sum2, weight, dcol
    real,allocatable          :: khc(:,:), khf(:,:), kvc(:,:), kvf(:,:)
    integer                   :: iwell, ilay, inode, iarr, igeo, i, wcount, count
    real,allocatable          :: wellxy(:,:), wellv(:),temppc(:)
    character(30),allocatable :: tempzones(:)
    character(60),allocatable :: fmt(:)
    type(pddataset)           :: codata_nona_lay(ndata-1)
		type(kdata)								:: cokdata(ndata)
    type(tree_array),allocatable :: dtrees(:)

!  Avoiding stack overflow by allocating after declaration
    allocate(KhC(nnodes, nlayers), &
             KhF(nnodes, nlayers), &
             KvC(nnodes, nlayers), &
             KvF(nnodes, nlayers)  )
    allocate(wellxy(2,nwell), wellv(nwell),temppc(nnodes))
    allocate(tempzones(nnodes),fmt(2))
    allocate(dtrees(ndata))

    if (USE_IDW == 0) then
      ! Calculate percent coarse at each node by kriging well data by layer
      do ilay=1, nlayers
        tempzones = geozones(1:nnodes,ilay)
        ! Compile wells (x,y,val) that have PC for each layer/geo zone
        ! Use to build kd-tree for kriging subroutine
        do igeo=1, ngeozones_lay(ilay)
          wcount = 0
          do iwell=1, nwell 
            if ((pcwellelem(iwell, ilay,1) >= 0.0d0).and.&
                (wellgeozones(iwell,ilay) == geonames_lay(igeo,ilay))) then
              wcount = wcount + 1
              wellv(wcount) = pcwellelem(iwell, ilay,1)
              wellxy(1,wcount) = xwell(iwell)
              wellxy(2,wcount) = ywell(iwell)
            end if
          end do
          if (wcount > 0) then
					  ! All Data
					  if (all_data_flag==1) then
              ! Write some output so the user knows T2P is working (can be slow!)
              write(*,'(a,i2)') 'LAYER', ilay
						  ! Sure is great to copy this around a million times
						  cokdata(1)%n = wcount
						  call kdata_allocator(cokdata(1))
							cokdata(1)%xy(:,:) = wellxy(:,:)
							cokdata(1)%v(:) = wellv(:)
							cokdata(1)%v_mean = wellmean
              do i=2, ndata
                !TODO zone subset
						    call pd2kdata(codata(i-1), cokdata(i), ilay, 0)
							end do
						  call spcokrige_group(cokdata, ndata, nnodes, nodex, nodey, temppc, 1, &
													 0, tempzones, geonames_lay(igeo,ilay))
						else
						! Nearest Neightbors
              if (wcount < knnear(1)) then
                ! Warn if less wells than expected
                write(*,'(4x,a,i8,3a,i3,a)') 'Warning - Layer', ilay, ' - ', &
                         trim(geonames_lay(igeo,ilay)), ', only ', wcount, ' wells'
              end if
              dtrees(1)%tp => kdtree2_create(wellxy(:,1:wcount))
              ! Pass well tree to kriging routine
              ! Will ignore nodes not in the target group
              if (ndata>1) then
                ! Create cokriging trees
                ! ZONES UNIMPLEMENTED
                do i=2, ndata
                  !TODO zone subset
                  call pd4tree(codata(i-1), codata_nona_lay(i-1), ilay, 0)
                  dtrees(i)%tp => kdtree2_create(codata_nona_lay(i-1)%xy(:,:))
                end do
                !cokrige
                call spcokrige_tree_group(codata_nona_lay, ndata, nnodes, nodex, nodey, temppc, 1, &
                                          dtrees, knnear, 0, &
                                          wcount, wellv, wellmean, tempzones, geonames_lay(igeo,ilay))
              else
                ! One data -> ordinary kriging (ind=1, first variogram)
                call spkrige_tree_group(nnodes, nodex, nodey, dtrees(1)%tp, wcount, &
                                  wellv, 1, temppc, knnear(1), &
                                  ggroup=tempzones, targetgroup=geonames_lay(igeo,ilay))
              end if
            end if
          else
            write(*,'(4x,a,i8,2a)') 'Warning - No well data in layer ', ilay, &
                                   'for geologic zone', geonames_lay(igeo,ilay)
            ! Only affect one zone
            where (tempzones == geonames_lay(igeo,ilay)) temppc = NODATA
          end if
        end do
        ! Enforce 0-1 boundaries of percent coarse
        where(temppc < 0.0)
          temppc = 0
        elsewhere(temppc > 1.0)
          temppc = 1
        end where
        ! Copy over
        PcNode(1:nnodes,ilay) = temppc
        ! Destroy trees & temp pddatasets
        do i=1, ndata
          if (all_data_flag==0) then
            call kdtree2_destroy(dtrees(i)%tp)
            if (i>1) then
              call deallocate_pddataset(codata_nona_lay(i-1))
            end if
          else !== alldata
            call kdata_deallocator(cokdata(i))
          end if
        end do
      end do
    else  ! i.e. USE_IDW == 1
      ! Calculate percent coarse for each node based on normalized inverse-distance
      write(*,'(4x,a)') 'Calculating inverse-distance well-node weights'
      do ilay = 1, nlayers
        write(*,'(4x,a,i2)') 'Processing layer ', ilay
        do inode = 1, nnodes
          sum1 = 0.0
          sum2 = 0.0
          count = 0            ! LS Added to avoid "no well data" when PC is 0.0 from data
          do iwell = 1, nwell
            if ((PcWellElem(iwell, ilay,1) >= 0).and.&
                       (geozones(inode,ilay) == wellgeozones(iwell,ilay))) then
              count = count + 1
              weight = (((Xwell(Iwell) - nodex(inode)) ** 2 &
                       + (Ywell(Iwell) - nodey(inode)) ** 2) ** 0.5) ** 2
              weight = 1 / weight
              sum1 = sum1 + weight * PcWellElem(iwell, ilay,1)
              sum2 = sum2 + weight
            end if
          end do
          if (count > 0) then                               ! LS sum2 > 0
            PcNode(inode, ilay) = sum1 / sum2
          else
            PcNode(inode, ilay) = NODATA
            write(*,'(4x,a,i8)') 'No well data for node ', inode
          end if
        end do
      end do
    end if

! Use input values to calculate aquifer parameters for each node
    !write(*,'(4x,a)') 'Calculating element aquifer parameters'
    do ilay = 1, nlayers
      iarr = pp_lay2arr(ilay)  ! Get array index
      do inode = 1, nnodes
        Dcol = elevation(inode) - ((topelev(inode,ilay) + botelev(inode, ilay)) / 2.0)
        if (PcNode(inode, ilay) >= 0.0d0) then
            KhC(inode, ilay) = nKCMin(inode, iarr) + (dKC(inode,iarr) * &
                               exp(-1.0 * KCk * Dcol))
            KhF(inode, ilay) = nKFMin(inode, iarr) + (dKF(inode,iarr) * &
                               exp(-1.0 * KFk * Dcol))
            KvC(inode, ilay) = KhC(inode, ilay) / nAnisoC(inode,iarr)
            KvF(inode, ilay) = KhF(inode, ilay) / nAnisoF(inode,iarr)
            KhB(inode, ilay) = ((PcNode(inode, ilay) * &
                                (KhC(inode, ilay) ** KHp)) + &
                                ((1 - PcNode(inode, ilay)) * &
                                 (KhF(inode, ilay) ** KHp))) ** (1.0 / KHp)
            KvB(inode, ilay) = ((PcNode(inode, ilay) * &
                                (KvC(inode, ilay) ** KVp)) + &
                                ((1 - PcNode(inode, ilay)) * &
                                 (KvF(inode, ilay) ** KVp))) ** (1.0 / KVp)

!            SsB(inode, ilay) = (PcNode(inode, ilay) * (nSsC(inode))) + &
!                               ((1 - PcNode(inode, ilay)) * (nSsF(inode)))
! LS 1/6/2022 Storage parameters now share a common power law term, stp
            SsB(inode, ilay) = ((PcNode(inode, ilay) * (nSsC(inode, iarr) ** Stp)) + &
                                ((1 - PcNode(inode, ilay)) * (nSsF(inode, iarr) ** &
                                  Stp))) ** (1.0 / Stp)
            SyB(inode, ilay) = ((PcNode(inode, ilay) * (nSyC(inode, iarr) ** Stp)) + &
                                ((1 - PcNode(inode, ilay)) * (nSyF(inode, iarr) ** &
                                  Stp))) ** (1.0 / Stp)

! Removing these, not used in model - LS 9/7/2017
!            SskeB(inode, ilay) = (PcNode(inode, ilay) * SskeC) + &
!                                 ((1 - PcNode(inode, ilay)) * SskeF)
!            SskvB(inode, ilay) = (PcNode(inode, ilay) * SskvC) + &
!                                 ((1 - PcNode(inode, ilay)) * SskvF)
        else
            KhB(inode, ilay) = NODATA
            KvB(inode, ilay) = NODATA
            SsB(inode, ilay) = NODATA
            SyB(inode, ilay) = NODATA
!            SskeB(inode, ilay) = NODATA
!            SskvB(inode, ilay) = NODATA
        end if
      end do
    end do

    ! Deallocate
    deallocate(KhC, KhF, KvC, KvF)
    deallocate(wellxy, wellv, temppc)
    deallocate(tempzones,fmt)

end subroutine interpolate
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
SUBROUTINE interpolate_aqtard()
  use Krige
  use PilotPoints
  use kdtree2_module
    implicit none
!-----------------------------------------------------------------------------!
! Interpolates from wells to nodes using inverse-distance method
!-----------------------------------------------------------------------------!

    real                     :: sum1, sum2, weight, dcol
    real,allocatable         :: khc(:,:), khf(:,:), kvc(:,:), kvf(:,:), pc(:)
    integer                  :: i,iwell, ilay, inode, iarr, wcount, count
    real,allocatable         :: wellxy(:,:), wellv(:)
    character(60),allocatable:: fmt(:)
    type(pddataset)           :: codata_nona_lay(ndata-1)
		type(kdata)								:: cokdata(ndata-1)
    type(tree_array),allocatable :: dtrees(:)
    allocate(dtrees(ndata))

    allocate(KhC(nnodes, nlayers), &
             KhF(nnodes, nlayers), &
             KvC(nnodes, nlayers), &
             KvF(nnodes, nlayers), &
             pc(nnodes))
    allocate(wellxy(2,nwell), wellv(nwell))
    allocate(fmt(2))
    if (USE_IDW == 0) then
      ! Calculate percent coarse at each node by kriging well data by layer
      do ilay=1, nlayers
        ! Get wells (x,y,val) that have PC for this layer to build kd-tree
        wcount = 0
        do iwell=1, nwell
          if (PcWellAqtardElem(iwell, ilay,1) >= 0.0d0) then
            wcount = wcount + 1
            wellv(wcount) = PcWellAqtardElem(iwell, ilay,1)
            wellxy(1,wcount) = xwell(iwell)
            wellxy(2,wcount) = ywell(iwell)
          end if
        end do

        ! Only build tree & krige if there's data for the layer
        if (wcount > 0) then
					! All Data
					if (all_data_flag==1) then
            do i=2, ndata
              !TODO zone subset
						  call pd2kdata(codata(i-1), cokdata(i-1), ilay, 1)
						end do
						call spcokrige_group(cokdata, ndata, nnodes, nodex, nodey, pc, 1, 1)
					else
					  ! Nearest Neightbors
            dtrees(1)%tp => kdtree2_create(wellxy(:,1:wcount))
            if (ndata>1) then
              ! Create cokriging trees
              do i=2, ndata
                call pd4tree(codata(i-1), codata_nona_lay(i-1), ilay, 1)
                dtrees(i)%tp => kdtree2_create(codata(i)%xy(:,:))
              end do
              !cokrige
              call spcokrige_tree_group(codata_nona_lay, ndata, nnodes, nodex, nodey, pc, 1, &
                                        dtrees, knnear, 1, wcount, wellv, wellmean)
            else
              call spkrige_tree(nnodes, nodex, nodey, dtrees(i)%tp, wcount, wellv, &
                                1, pc, knnear(1))
            end if
            ! Enforce 0-1 boundaries of percent coarse
            where(pc < 0.0)
              pc = 0
            elsewhere(pc > 1.0)
              pc = 1
            end where
            PcNode_aqtard(1:nnodes,ilay) = pc(1:nnodes)
            ! Destroy trees & temp pddatasets
            do i=1, ndata
              if (all_data_flag==0) then
                call kdtree2_destroy(dtrees(i)%tp)
                if (i>1) then
                  call deallocate_pddataset(codata_nona_lay(i-1))
                end if
              else !== alldata
                call kdata_deallocator(cokdata(i))
              end if
            end do
					end if
        else
          ! Check for layer thickness
          do inode = 1, nnodes
            if ((aqtardtopelev(inode, ilay) - aqtardbotelev(inode, ilay)) > 0.0d0) then
              ! No wells in layer, despite having thickness
              pcnode_aqtard(inode, ilay) = NODATA
            else
              ! No thickness, no problem. PC is 0.0
              pcnode_aqtard(inode, ilay) = 0.0d0
            end if
          end do
          if (minval(pcnode_aqtard(:,ilay))<0.0d0) then
            ! Need to warn if NODATA values were written
            write(*,'(4x,a,i8)') 'No well data for aquitard layer ', ilay
          end if
        end if
      end do
    else
      ! Calculate percent course for each node based on normalized inverse-distance
      write(*,'(4x,a)') 'Calculating inverse-distance well-node weights'
      do ilay = 1, nlayers
  !      write(*,'(4x,a,i2)') 'Processing layer ', ilay
        do inode = 1, nnodes
          ! Check to make sure there is aquifer thickness
          if ((aqtardtopelev(inode, ilay) - aqtardbotelev(inode, ilay)) > 0.0d0) then
            sum1 = 0.0
            sum2 = 0.0
            count = 0
            do iwell = 1, nwell
              if (PcWellAqtardElem(iwell, ilay,1) > -1e-5) then                   ! PcWellAqtardElem is -999 if blank
                count = count + 1
                weight = (((Xwell(Iwell) - nodex(inode)) ** 2 &
                         + (Ywell(Iwell) - nodey(inode)) ** 2) ** 0.5) ** 2
                weight = 1 / weight
                sum1 = sum1 + weight * PcWellAqtardElem(iwell, ilay,1)
                sum2 = sum2 + weight
              end if
            end do
            if (count > 0) then
              pcnode_aqtard(inode, ilay) = sum1 / sum2
            else
              pcnode_aqtard(inode, ilay) = NODATA
              write(*,'(4x,a,i8)') 'No well data for node ', inode
            end if
          else
            ! No thickness - no error, just zero
            pcnode_aqtard(inode, ilay) = 0.0
          end if
        end do
      end do
    end if

! Use input values to calculate aquifer parameters for each node
    !write(*,'(4x,a)') 'Calculating element aquitard parameters'
    do ilay = 1, natlayers
      iarr = pp_lay2arr(nlayers + ilay) ! Get array index
      do inode = 1, nnodes
        Dcol = elevation(inode) - ((aqtardtopelev(inode,ilay) + aqtardbotelev(inode, ilay)) / 2.0)
        if (pcnode_aqtard(inode, ilay) >= 0.0d0) then
            KhC(inode, ilay) = nKCMin(inode, iarr) + (dKC(inode,iarr) * &
                               exp(-1.0 * KCk * Dcol))
            KhF(inode, ilay) = nKFMin(inode, iarr) + (dKF(inode,iarr) * &
                               exp(-1.0 * KFk * Dcol))
            KvC(inode, ilay) = KhC(inode, ilay) / nAnisoC(inode,iarr)
            KvF(inode, ilay) = KhF(inode, ilay) / nAnisoF(inode,iarr)
            if ((aqtardtopelev(inode, ilay) - aqtardbotelev(inode, ilay)) > 1e-5) then
              ! Has thickness
              KvB_aqtard(inode, ilay) = ((pcnode_aqtard(inode, ilay) * &
                                  (KvC(inode, ilay) ** KVp)) + &
                                  ((1 - pcnode_aqtard(inode, ilay)) * &
                                   (KvF(inode, ilay) ** KVp))) ** (1.0 / KVp)
            else
              ! No Thickness
              KvB_aqtard(inode, ilay) = 0.0d0
            end if
        else
            KvB_aqtard(inode, ilay) = NODATA
        end if
      end do
    end do

    ! Deallocate
    deallocate(KhC, KhF, KvC, KvF)
    deallocate(wellxy, wellv)
    deallocate(fmt)

end subroutine interpolate_aqtard
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
SUBROUTINE interpolate_multiclass()
    use Krige
    use PilotPoints
    use kdtree2_module
    implicit none
!-----------------------------------------------------------------------------!
! Interpolates from wells to nodes
!
! Well GeoZones Loop added 7/1/2019
! Utilizes arrays geonames_lay and ngeozones_lay which list geologic zones
! ("Units") present in each model layer, based upon geologic node zones
! -LS
!-----------------------------------------------------------------------------!
    real                      :: sum1, sum2, weight, dcol
    real,allocatable          :: khD(:,:,:), kvd(:,:,:), tcvals(:,:)
    integer                   :: iwell, ilay, iclass, inode, iarr, igeo, i, wcount, count
    real,allocatable          :: wellxy(:,:), wellv(:),tempv(:)
    character(30),allocatable :: tempzones(:)
    character(60),allocatable :: fmt(:)
    type(pddataset)           :: codata_nona_lay(ndata-1)
		type(kdata)								:: cokdata(ndata)
    type(tree_array),allocatable :: dtrees(:)

!  Avoiding stack overflow by allocating after declaration
    allocate(KhD(nnodes, nlayers, ntclasses), &
             KvD(nnodes, nlayers, ntclasses), &
             tcvals(nnodes, ntclasses))
    allocate(wellxy(2,nwell), wellv(nwell),tempv(nnodes))
    allocate(tempzones(nnodes),fmt(2))
    allocate(dtrees(ndata))

    ! Loop over classes
    do iclass=1, ntclasses
      ! Calculate class value at each node by kriging well data by layer
      do ilay=1, nlayers
        tempzones = geozones(1:nnodes,ilay)
        ! Compile wells (x,y,val) that have PC for each layer/geo zone
        ! Use to build kd-tree for kriging subroutine
        do igeo=1, ngeozones_lay(ilay)
          wcount = 0
          do iwell=1, nwell 
            if ((pcwellelem(iwell, ilay,1) >= 0.0d0).and.&
                (wellgeozones(iwell,ilay) == geonames_lay(igeo,ilay))) then
              wcount = wcount + 1
              wellv(wcount) = pcwellelem(iwell, ilay, iclass)
              wellxy(1,wcount) = xwell(iwell)
              wellxy(2,wcount) = ywell(iwell)
            end if
          end do
          if (wcount > 0) then
            if (wcount < knnear(1)) then
              ! Warn if less wells than expected
              write(*,'(4x,a,i8,3a,i3,a)') 'Warning - Layer', ilay, ' - ', &
                        trim(geonames_lay(igeo,ilay)), ', only ', wcount, ' wells'
            end if
            dtrees(1)%tp => kdtree2_create(wellxy(:,1:wcount))
            ! Pass well tree to kriging routine
            ! Will ignore nodes not in the target group
            ! One data -> ordinary kriging (ind=1, first variogram)
            call spkrige_tree_group(nnodes, nodex, nodey, dtrees(1)%tp, wcount, &
                              wellv, 1, tempv, knnear(1), &
                              ggroup=tempzones, targetgroup=geonames_lay(igeo,ilay))
          else
            write(*,'(4x,a,i8,2a)') 'Warning - No well data in layer ', ilay, &
                                    'for geologic zone', geonames_lay(igeo,ilay)
            ! Only affect one zone
            where (tempzones == geonames_lay(igeo,ilay)) tempv = NODATA
          end if
        end do
        ! Copy over
        nodecval(:,ilay,iclass) = tempv
        ! Destroy trees & temp pddatasets
        do i=1, ndata
          if (all_data_flag==0) then
            call kdtree2_destroy(dtrees(i)%tp)
            if (i>1) then
              call deallocate_pddataset(codata_nona_lay(i-1))
            end if
          else !== alldata
            call kdata_deallocator(cokdata(i))
          end if
        end do
      end do
    end do

! Use input values to calculate aquifer parameters for each node
    !write(*,'(4x,a)') 'Calculating element aquifer parameters'
    do ilay = 1, nlayers
      iarr = pp_lay2arr(ilay)  ! Get array index
      do inode = 1, nnodes
        Dcol = elevation(inode) - ((topelev(inode,ilay) + botelev(inode, ilay)) / 2.0)
        
        ! Ensure classes add to one
        ! TODO add dominant classes option
        nodecval(inode,ilay,:) = nodecval(inode,ilay,:)/sum(nodecval(inode,ilay,:))
        
        ! Calculate Depth-decay hydraulic conductivity values...
        KhD(inode,ilay,:) = KMinCvals(inode,iarr,:) + dKCvals(inode,iarr,:) * exp(-1.0 * KkCvals(:) * Dcol)
        KvD(inode,ilay,:) = KhD(inode,iarr,:) / nAnisoCvals(inode,iarr,:)
        
        ! Calculate bulk conductivity values
        KhB(inode,ilay) = DOT_PRODUCT(KhD(inode,ilay,:)**KHp, nodecval(inode,ilay,:)) ** (1.0 / KHp)
        KvB(inode,ilay) = DOT_PRODUCT(KvD(inode,ilay,:)**Kvp, nodecval(inode,ilay,:)) ** (1.0 / Kvp)
        
        ! Calculate bulk storage values
        SsB(inode,ilay) = DOT_PRODUCT(SsCvals(inode,ilay,:)**Stp, nodecval(inode,ilay,:)) ** (1.0 / Stp)
        SyB(inode,ilay) = DOT_PRODUCT(SyCvals(inode,ilay,:)**Stp, nodecval(inode,ilay,:)) ** (1.0 / Stp)
      end do
    end do

    ! Deallocate
    deallocate(KhD, KvD, tcvals)
    deallocate(wellxy, wellv, tempv)
    deallocate(tempzones,fmt)

end subroutine interpolate_multiclass
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
! - IWFM KD-tree related subroutines
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!

subroutine create_nodetree()
  use kdtree2_module
  implicit none

  integer                  :: i

  ! Reassignment into a single xy array
  ! memory usage could be improved by modifying KDTree
  do i=1, nnodes
    nodexy(1,i) = nodex(i)
    nodexy(2,i) = nodey(i)
  end do

  nodetree => kdtree2_create(nodexy, dim=2)

end subroutine create_nodetree

!-----------------------------------------------------------------------------!

integer function get_well_element(x, y, iwell, closenodes, nnear) result(element)
  implicit none

  integer                  :: i, j, inode, elem, points, count, node
  integer, intent(in)      :: nnear, iwell, closenodes(nnear)
  real, intent(in)         :: x, y
  real                     :: elemx(4), elemy(4), inelem, pinpol

  element = 0

  ! Loop over elements of closest node
outerloop:  do inode=1, nnear
    node = closenodes(inode)
    do i=1, nnodeelem(node)
      elem = nodeelem(node,i)
      count = 0
      ! Get Element node x,y coordinates
      do j=1, 4
        if (elements(elem,j) > 0) then
          elemx(j) = nodex(elements(elem,j))
          elemy(j) = nodey(elements(elem,j))
          count = count + 1
        end if
      end do
      ! Check if the point is in this polygon
      inelem = pinpol(x, y, elemx, elemy, count)
      if (inelem >= 0.0d0) then             ! Greedy, treats on line as in elem
        element = elem
        exit outerloop
      end if
    end do
  end do outerloop

  if (element == 0) then
    write(*,'(4x,a)') 'ERROR - Well not inside element. Likely located out of model.'
    write(*,'(4x,a8,2a18)') 'WellNum', 'X', 'Y'
    write(*,'(4x,i8,2f18.5)') iwell, x, y
    stop
  end if

  return

end function get_well_element

!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
! - Extra Subroutines
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine calc_interbed_thickness()
  implicit none

  ! Estimates interbed thickness as percent fine * layer thickness
  ! Added by Leland Scantlebury July 2019
  !
  ! For IWFM, layer thickness is directly read in in variable thick
  ! However, for MODFLOW support we're just going to use the layer top/bot

  integer                  :: ilay, inode
  real, allocatable        :: laythick(:,:)

  allocate(interbed_thick(nnodes,nlayers),laythick(nnodes,nlayers))

  laythick = topelev - botelev
  interbed_thick = (1 - pcnode) * laythick

end subroutine calc_interbed_thickness
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine colocate_corrector(xcoords, ycoords, ncoords)
implicit none
! Checks for duplicate X,Y values and moves them away from each other.
!
! Added by Leland Scantlebury July 2019
  integer, intent(in)     :: ncoords
  real, intent(inout)     :: xcoords(ncoords), ycoords(ncoords)
  integer                 :: i, j, k, ncoloc, max
  real, parameter         :: move=0.01

  ! Setup
  max = 1

  ! Outer Loop over coordinates
  do i=1, ncoords
    ncoloc = 0
    ! Inner Loop over coordinates to search for pairs
    do j=1, ncoords
      ! Skip the outer loop coordinate
      if (i==j) cycle
      ! Check for x-y match
      if (xcoords(i) == xcoords(j).and.ycoords(i) == ycoords(j)) then
        ! Warn if this is our first match
        if (colocated == .false.) then
          colocated = .true.
          write(*,'(4x,a)') 'Warning - co-located wells. Attempting to move apart for kriging'
        end if
        ! Add to count
        ncoloc = ncoloc + 1
        ! See if greater than max
        if (ncoloc > max) max = ncoloc
        ! Move
        xcoords(j) = xcoords(j) + move * ncoloc
        ycoords(j) = ycoords(j) + move * ncoloc
      end if
    end do
  end do

  if (colocated) then
    write(*,'(4x,a,i3,a,f6.3)') 'Max co-located wells:',max,' | Max dist moved:', sqrt(2*(max*move)**2)
  end if

end subroutine colocate_corrector
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
! Allocation Subroutine(s)
!-----------------------------------------------------------------------------!

subroutine WellAllocate()
  implicit none

    allocate(Xwell(nwell), &
             Ywell(nwell), &
             Zland(nwell), &
             Zwell(nwell, maxWPoints), &
             wellelem(nwell), &
             numZone(nwell), &
             PcWellZone(nwell, maxWPoints, ntclasses), &
             PcWellElem(nwell, nlayers,ntclasses), &
             PcWellAqtardElem(nwell,nlayers,ntclasses), &
             wellelevtop(nwell, nlayers), &
             wellelevbot(nwell, nlayers), &
             wellaqtardelevtop(nwell,nlayers), &
             wellaqtardelevbot(nwell,nlayers), &
             wellgeozones(nwell, nlayers))

end subroutine
!-----------------------------------------------------------------------------!
subroutine NodeAllocate()
  use PilotPoints
  implicit none

  allocate(PcNode       (nnodes, nlayers), &
           pcnode_aqtard(nnodes, nlayers), &
           nKCMin       (nnodes, pp_narr), &
           dKC          (nnodes, pp_narr), &
           nKFMin       (nnodes, pp_narr), &
           dKF          (nnodes, pp_narr), &
           nAnisoC      (nnodes, pp_narr), &
           nAnisoF      (nnodes, pp_narr), &
           nSsC         (nnodes, ppaq_narr), &
           nSsF         (nnodes, ppaq_narr), &
           nSyC         (nnodes, ppaq_narr), &
           nSyF         (nnodes, ppaq_narr), &
           KhB          (nnodes, nlayers), &
           KvB          (nnodes, nlayers), &
           SsB          (nnodes, nlayers), &
           SyB          (nnodes, nlayers), &
           KvB_aqtard   (nnodes, nlayers), &
           geozones     (nnodes, nlayers), &
           ngeozones_lay(nlayers), &
           geonames_lay (max_zones,nlayers))
           
  if (ntclasses > 1) then
    allocate(nodecval       (nnodes,nlayers,ntclasses), &
             nodecval_aqtard(nnodes,nlayers,ntclasses), &
               KMinCvals    (nnodes,pp_narr,ntclasses), &
                 dKCvals    (nnodes,pp_narr,ntclasses), &
                 SsCvals    (nnodes,pp_narr,ntclasses), &
                 SyCVals    (nnodes,pp_narr,ntclasses), &
             nAnisoCvals    (nnodes,pp_narr,ntclasses)  )
  end if

end subroutine
!-----------------------------------------------------------------------------!

END MODULE MakePar