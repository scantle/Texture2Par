module pddata
  implicit none
!-----------------------------------------------------------------------------!
! Point depth data
! Data, at a point, containing data with depth
! e.g., a well log
!-----------------------------------------------------------------------------!
  SAVE
  
!-----------------------------------------------------------------------------!
! Module Data Structures
!-----------------------------------------------------------------------------!
  type :: pddataset

    integer                  :: npdat            ! Number of datapoints
    integer, allocatable     :: ndint(:), &      ! Number of depth intervals for each point
                                zone(:,:)        ! GeoZone, ie HSU of (points, layers)
    real, allocatable        :: x(:),     &      ! X-coordinates of points
                                y(:),     &      ! Y-coordinates of points
                                z(:),     &      ! Elevation above datum (zland - depth)
                                zland(:), &      ! Land surface elevation 
                                v(:,:)           ! Values (interval, point)
  end type pddataset
!-----------------------------------------------------------------------------!

  
!-----------------------------------------------------------------------------!
! Module Variables
!-----------------------------------------------------------------------------!
  
!-----------------------------------------------------------------------------!
! Module Subroutines & Functions
!-----------------------------------------------------------------------------!
  CONTAINS

subroutine pre_read_pddata(unit, file, filelen, npoints, maxint)
  implicit none
!-----------------------------------------------------------------------------!
!
!-----------------------------------------------------------------------------!
  integer,intent(in)             :: unit
  character(150),intent(in)      :: file, jnk
  integer,intent(out)            :: filelen, npoints, maxint
  integer                        :: ioerr, pid, interal
  
  
  open(unit, file = trim(wellfile), status='old', iostat = ierr,  iomsg=cerr)
  call iomsghandler(ierr, cerr)
  read(unit,*)                                    ! Header
  ioerr = 0
  filelen = 0
  npoints = 0
  maxint = 0
  do while (ioerr == 0)
    read(unit,*, iostat = ioerr) jnk, pid, interval  ! name, point id, point depth interval
    filelen = filelen + 1
    if (id > npoints) npoints = pid
    if (interval > maxint) maxint = interval
  end do
  filelen = filelen - 1
  close(unit)
end subroutine pre_read_pddata

!-----------------------------------------------------------------------------!
    
end module pddata