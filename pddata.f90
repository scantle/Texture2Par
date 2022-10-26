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

subroutine allocate_pddataset(pdset, maxint, nlay)
  implicit none
  
  type(pddataset),intent(inout)  :: pdset
  integer,intent(in)          :: maxint, nlay
  
  allocate(pdset%ndint   (pdset%npdat),    &
           pdset%xy    (2,pdset%npdat),    &
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

end subroutine allocate_pddataset
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine deallocate_pddataset(pdset)
  implicit none
  
  type(pddataset),intent(inout)  :: pdset
  
  deallocate(pdset%ndint,   &
             pdset%xy,      &
             pdset%zland,   &
             pdset%node,    &
             pdset%z,   &
             pdset%v,   &
             pdset%zone  ,  &
             pdset%laytop,  &
             pdset%laybot,  &
             pdset%aqttop,  &
             pdset%aqtbot,  &
             pdset%layv  ,  &
             pdset%aqtv     )

end subroutine deallocate_pddataset
!-----------------------------------------------------------------------------!

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

subroutine pd4tree(pdset, pdsetout, lay, aquitard)
  implicit none
  
  type(pddataset),intent(inout) :: pdset, pdsetout
  integer,intent(in)            :: lay, aquitard
  integer                       :: i, p
  real,parameter                :: cutoff=-1.0
  
  if (aquitard==1) then
    pdsetout%npdat = count(pdset%aqtv(:,lay) > cutoff)
  else
    pdsetout%npdat = count(pdset%layv(:,lay) > cutoff)
  end if
  
  call allocate_pddataset(pdsetout, 1, 1)
  
  p = 1
  if (aquitard==1) then
    do i=1, pdset%npdat
      if (pdset%aqtv(i,lay) > cutoff) then
        pdsetout%aqtv(p,1) = pdset%aqtv(i,lay)
        pdsetout%xy(:,p) = pdset%xy(:,i)
        pdsetout%zone(p,1) = pdset%zone(i,lay)
        p = p+1
      end if
    end do
    !pdsetout%v_mean = sum(pdsetout%aqtv(:,1)) / pdsetout%npdat
  else
    do i=1, pdset%npdat
      if (pdset%layv(i,lay) > cutoff) then
        pdsetout%layv(p,1) = pdset%layv(i,lay)
        pdsetout%xy(:,p) = pdset%xy(:,i)
        pdsetout%zone(p,1) = pdset%zone(i,lay)
        p = p+1
      end if
    end do
    !pdsetout%v_mean = sum(pdsetout%layv(:,1)) / pdsetout%npdat
  end if
  
  ! maintain mean
  pdsetout%v_mean = pdset%v_mean    

end subroutine pd4tree

!-----------------------------------------------------------------------------!

subroutine write_pd_file(pdset, filename, nlayers)
  type(pddataset) :: pdset
  character(150), intent(in) :: filename
  integer, intent(in)        :: nlayers
  integer :: i, ilay
  character(60)            :: fmt(2)

  write(*,'(4x,a,a35)') 'Writing:', filename

  write(fmt(1),'("(2a8,2a18,",i5,"i14)")') nlayers
  write(fmt(2),'("(2i8,2f18.5,",i5,"e14.5)")') nlayers
  
  open(20, file=trim(filename), status='replace', action='write', buffered='YES')
  write(20,fmt(1)) 'Point', 'Element', 'X', 'Y', (ilay, ilay=1, nlayers)
  ! Values
  do i = 1, pdset%npdat
    write(20,fmt(2)) i, pdset%node(i), pdset%xy(1:2,i), (pdset%layv(i,ilay), ilay=1, nlayers)
  end do
  
  close(20)
  
end subroutine write_pd_file

!-----------------------------------------------------------------------------!

end module pddata