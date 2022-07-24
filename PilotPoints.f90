module PilotPoints
  implicit none
!-----------------------------------------------------------------------------!
! Pilot Point Module
! 
! Author: Leland Scantlebury of S.S. Papadopulos & Assoc.
! Improved/Moved by LS July 2022
!-----------------------------------------------------------------------------!

! Variables related to pilot points and their interpolated (kriged) nodal values 
  SAVE
  integer                 :: ppn, atppn             ! Counts for aquifer & aquitard pilot points 
  integer,allocatable     :: ppZone(:),     &       ! Zone of each aquifer pilot point
                             ppatZone(:)            ! Zone of each aquitard pilot point
  real, allocatable       :: ppx(:),        &       ! X-coordinate of aquifer pilot points
                             ppy(:),        &       ! Y-coordinate of aquifer pilot points
                             ppKCMin(:),    &       ! Coarse hydraulic conductivity minimum
                             ppdKC(:),      &       ! Coarse hydraulic conductivity maximum change (delta)
                             ppKFMin(:),    &       ! Fine hydraulic conductivity minimum
                             ppdKF(:),      &       ! Fine hydraulic conductivity maximum change (delta)
                             ppSsC(:),      &       ! Coarse specific storage
                             ppSsF(:),      &       ! Fine specific storage
                             ppSyC(:),      &       ! Coarse specific yield
                             ppSyF(:),      &       ! Fine specific yield
                             ppAnisoC(:),   &       ! Coarse Anisotropy (Kh/Kv)
                             ppAnisoF(:),   &       ! Fine Anisotropy (Kh/Kv)
                             ppatx(:),      &       ! X-coordinate of aquitard pilot points
                             ppaty(:),      &       ! Y-coordinate of aquitard pilot points
                             ppatKCMin(:),  &       ! Aquitard coarse hydraulic conductivity minimum
                             ppatdKC(:),    &       ! Aquitard coarse hydraulic conductivity maximum change (delta)
                             ppatKFMin(:),  &       ! Aquitard fine hydraulic conductivity minimum
                             ppatdKF(:),    &       ! Aquitard fine hydraulic conductivity maximum change (delta)
                             ppatAnisoC(:), &       ! Aquitard coarse Anisotropy (Kh/Kv)
                             ppatAnisoF(:)          ! Aquitard fine Anisotropy (Kh/Kv)

! Variables related to pilot point zone arrays (nodes/cells)
  character(150)            :: zonefile              ! File with zones, read differently for IWFM & MF
  integer                   :: pp_narr,           &  ! Number of pilot point zone arrays
                               ppaq_narr,         &  ! Number of aquifer pilot point zone arrays
                               ppat_narr             ! Number of aquitard pilot point zone arrays
  integer,allocatable       :: pp_lay2arr(:),     &  ! Map from layers to arrays, after nlayers is aquitards
                               pp_node_zone(:,:), &  ! Zone membership by nodes (cells) (nnodes, pp_narr)
                               pp_arr_zones(:,:), &  ! Zones present in each array (nzones, pp_narr)
                               pp_arr_nzones(:)          ! Count of zones present in each array

CONTAINS  !-- Module Subroutines

!-----------------------------------------------------------------------------!
subroutine ppaq_parkrige_arr(nnodes, nodex, nodey, pp_par_vals, kriged)
  use krige, only: spkrige_group
  implicit none
!-----------------------------------------------------------------------------!
! Loops over zones kriging pilot point parameter values to nodes
!-----------------------------------------------------------------------------!

  integer,intent(in)       :: nnodes
  real,intent(in)          :: nodex(nnodes), nodey(nnodes), pp_par_vals(nnodes)
  real,intent(out)         :: kriged(nnodes, pp_narr)

  integer                  :: i, iarr, ipp, ppcount, zone
  real                     :: pzx(ppn), pzy(ppn), pzv(ppn)

  ! Loop over aquifer arrays
  do iarr=1, ppaq_narr
    ! Loop over zones present in array
    do i=1, pp_arr_nzones(iarr)
      zone = pp_arr_zones(i,iarr)    ! Zones do not necessarily have to be in any order
      ! Collect pilot points in zone
      ppcount = 0
      do ipp=1, ppn
        if (ppZone(ipp) == zone) then
          ppcount = ppcount + 1
          pzx(ppcount) = ppx(ipp)
          pzy(ppcount) = ppy(ipp)
          pzv(ppcount) = pp_par_vals(ipp)
        end if
      end do ! pilot points
      if (ppcount > 0) then
        ! Krige
        call spkrige_group(ppcount, nnodes, nodex, nodey, pzx, pzy, pzv, &
                           kriged(:, iarr), pp_node_zone, zone)
      else
        write(*,'(2(a,i3))') 'ERROR - no pilot points for pilot point zone ', &
                   zone, ' in array ', iarr
        stop
      end if
    end do ! Zones
  end do ! Arrays

end subroutine ppaq_parkrige_arr
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine ppat_parkrige_arr(nnodes, nodex, nodey, pp_par_vals, kriged)
  use krige, only: spkrige_group
  implicit none
!-----------------------------------------------------------------------------!
! Loops over zones kriging pilot point parameter values to nodes
!-----------------------------------------------------------------------------!

  integer,intent(in)       :: nnodes
  real,intent(in)          :: nodex(nnodes), nodey(nnodes), pp_par_vals(nnodes)
  real,intent(out)         :: kriged(nnodes, pp_narr)

  integer                  :: i, iarr, ipp, ppcount, zone
  real                     :: pzx(ppn), pzy(ppn), pzv(ppn)

  ! Loop over aquifer arrays
  do iarr=ppaq_narr, pp_narr
    ! Loop over zones present in array
    do i=1, pp_arr_nzones(iarr)
      zone = pp_arr_zones(i,iarr)    ! Zones do not necessarily have to be in any order
      ! Collect pilot points in zone
      ppcount = 0
      do ipp=1, atppn
        if (ppatZone(ipp) == zone) then
          ppcount = ppcount + 1
          pzx(ppcount) = ppatx(ipp)
          pzy(ppcount) = ppaty(ipp)
          pzv(ppcount) = pp_par_vals(ipp)
        end if
      end do ! pilot points
      if (ppcount > 0) then
        ! Krige
        call spkrige_group(ppcount, nnodes, nodex, nodey, pzx, pzy, pzv, &
                           kriged(:, iarr), pp_node_zone, zone)
      else
        write(*,'(2(a,i3))') 'ERROR - no pilot points for pilot point zone ', &
                   zone, ' in array ', iarr
        stop
      end if
    end do ! Zones
  end do ! Arrays

end subroutine ppat_parkrige_arr
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine PPZoneAllocate(nzones, nnodes)
  implicit none
  ! Called DURING zone reading subroutines
  integer,intent(in) :: nzones, nnodes
  
  ! Must be set: nlayers, pp_narr
  allocate(pp_arr_nzones   (pp_narr),    &
           pp_node_zone(nnodes,pp_narr), &
           pp_arr_zones(nzones,pp_narr)  )
  
end subroutine PPZoneAllocate

!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine AquiferPPAllocate()
  implicit none
  allocate(ppx(ppn), ppy(ppn), ppKCMin(ppn), ppdKC(ppn), ppKFMin(ppn), &
          ppdKF(ppn), ppSsC(ppn), ppSsF(ppn), ppSyC(ppn), ppSyF(ppn), &
          ppAnisoC(ppn), ppAnisoF(ppn), ppZone(ppn))
end subroutine AquiferPPAllocate
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine AquitardPPAllocate()
  implicit none
  allocate(ppatx(atppn), ppaty(atppn), ppatKCMin(atppn), ppatdKC(atppn), &
            ppatKFMin(atppn), ppatdKF(atppn), ppatAnisoC(atppn), &
            ppatAnisoF(atppn), ppatZone(atppn))
end subroutine AquitardPPAllocate
!-----------------------------------------------------------------------------!

end module