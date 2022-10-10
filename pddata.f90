module pddata
  implicit none
!-----------------------------------------------------------------------------!
! Point depth data
! Data, at a point location, containing data with depth
! e.g., a well log
!-----------------------------------------------------------------------------!
  SAVE
!-----------------------------------------------------------------------------!
! Module Data Structures
!-----------------------------------------------------------------------------!
  type :: pddataset

    integer                  :: npdat              ! Number of datapoints
    integer, allocatable     :: ndint(:),    &      ! Number of depth intervals for each point
                                zone(:,:),   &      ! GeoZone, ie HSU of (points, layers)
                                node(:)             ! Containing model node/element
    real                     :: v_mean
    real, allocatable        :: xy(:,:),        &   ! Coordinates (x,y) of points
                                zland(:),    &      ! Land surface elevation 
                                z(:,:),      &      ! Interval Bottom Elevation (interval, point)
                                v(:,:),      &      ! Values (interval, point)
                                laytop(:,:), &      ! Model layer top
                                laybot(:,:), &      ! Model layer bottom
                                aqttop(:,:), &      ! Model aquitard top
                                aqtbot(:,:), &      ! Model aquitard bottom
                                layv(:,:),   &      ! Values averaged by layer
                                aqtv(:,:)           ! Values averaged by aquitard layer
  end type pddataset
!-----------------------------------------------------------------------------!
  
!-----------------------------------------------------------------------------!
! Module Variables
!-----------------------------------------------------------------------------!
  
!-----------------------------------------------------------------------------!
! Module Subroutines & Functions
!-----------------------------------------------------------------------------!
  CONTAINS

subroutine allocate_pddataset(pdset, maxint, nlay, modflowflag)
  implicit none
  
  type(pddataset),intent(inout)  :: pdset
  integer,intent(in)          :: maxint, nlay, modflowflag
  
  allocate(pdset%ndint   (pdset%npdat),    &
           pdset%xy      (pdset%npdat,2),  &
           pdset%zland   (pdset%npdat),    &
           pdset%node    (pdset%npdat),    &
           pdset%z(maxint,pdset%npdat),    &
           pdset%v(maxint,pdset%npdat),    &
           pdset%zone  (pdset%npdat,nlay), &
           pdset%laytop(pdset%npdat,nlay), &
           pdset%laybot(pdset%npdat,nlay), &
           pdset%aqttop(pdset%npdat,nlay), &
           pdset%aqtbot(pdset%npdat,nlay), &
           pdset%layv  (pdset%npdat,nlay), &
           pdset%aqtv  (pdset%npdat,nlay)  )
           !kd tree poor data setup   pdset%xy    (2,pdset%npdat),    &

end subroutine allocate_pddataset
  
!-----------------------------------------------------------------------------!

subroutine pre_read_pddata(unit, file, filelen, npoints, maxint)
  use errorhandle
  implicit none
!-----------------------------------------------------------------------------!
!
!-----------------------------------------------------------------------------!
  integer,intent(in)             :: unit
  character(150),intent(in)      :: file
  character(256)                 :: cerr
  character(30)                  :: jnk
  integer,intent(out)            :: filelen, npoints, maxint
  integer                        :: ioerr, pid, interval
  
  open(unit, file = trim(file), status='old', iostat = ioerr, iomsg=cerr)
  call iomsghandler(ioerr, cerr)
  read(unit,*)                                    ! Header
  ioerr = 0
  filelen = 0
  npoints = 0
  maxint = 0
  do while (ioerr == 0)
    read(unit,*, iostat = ioerr) jnk, pid, interval  ! name, point id, point depth interval
    filelen = filelen + 1
    if (pid > npoints) npoints = pid
    if (interval > maxint) maxint = interval
  end do
  filelen = filelen - 1
  close(unit)
end subroutine pre_read_pddata

!-----------------------------------------------------------------------------!

subroutine write_pd_file(pdset, filename)
  type(pddataset) :: pdset
  character(150), intent(in) :: filename
  
  open(20, file=trim(filename), status='replace', action='write', buffered='YES')
  write(*,*) pdset%layv
  close(20)
  
end subroutine write_pd_file

!-----------------------------------------------------------------------------!

end module pddata