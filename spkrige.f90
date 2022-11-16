MODULE Krige
  implicit none
!-----------------------------------------------------------------------------!
! Kriging/Cokriging Module
!
! 
! Author: Leland Scantlebury of S.S. Papadopulos & Assoc.
! (with many subroutines/functions by others, credited below, where known)
!
! Note: Z-values are all set to zero. There is no interpolation done in the
!       z-direction since well pc values and pilot points effectively
!       have no z-value
!
! Note: number of rotation matricies is currently hardcoded to 1 (maxrot)
!    (more than that number would need to be changed to implement multiple)
!-----------------------------------------------------------------------------!
  SAVE

!-----------------------------------------------------------------------------!
! Module Data Structures
!-----------------------------------------------------------------------------!
  type :: vgm_struct
    ! Structure for storing variograms or cross variograms
    integer                  :: itype        ! Variogram model (0-linear, 1-spherical, 2-exponential)
    real                     :: sill,     &  ! Contribution of each structure, when summed (with nugget) the sill
                                a_hmax,   &  ! Range: distance where data are no longer correlated
                                a_hmin,   &  ! Range min (equal to amax when anis1/anis2 ==1)
                                nugget,   &  ! Nugget: sill immediately at distances > 0
                                ang1,     &  ! Angle 1 (azimuth correction)
                                ang2,     &  ! Angle 2 (dip correction)
                                ang3,     &  ! Angle 3 (plunge correction)
                                anis1,    &  ! Anisotropy 1 (horizontal)
                                anis2        ! Anisotropy 2 (vertical)
  end type vgm_struct
!-----------------------------------------------------------------------------!
  type :: kdata
    ! Structure for storing primary/secondary data
    integer                  :: n             ! Number of points
    real                     :: v_mean        ! Dataset arithmetic mean
    real, allocatable        :: xy(:,:), &    ! X-coordinates
                                v(:)          ! Values
  end type kdata
!-----------------------------------------------------------------------------!
  type :: grid2D
    ! Structure for storing grid points
    integer                  :: n             ! Number of points
    real, allocatable        :: x(:), &       ! X-coordinates
                                y(:)          ! Y-coordinates
  end type grid2D
!-----------------------------------------------------------------------------!  
    
!-----------------------------------------------------------------------------!
! Module Variables
!-----------------------------------------------------------------------------!
  integer                      :: ndata,         &
                                  nvario,        &
                                  maxstruct
  real                         :: max_secondary_dist
  integer, allocatable         :: nstruct(:)
  real, allocatable            :: rotmat(:,:,:,:)
  type(grid2D)                 :: kgrid
  type(kdata), allocatable     :: datasets(:)
  type(vgm_struct), allocatable :: vario(:,:)
  
!-----------------------------------------------------------------------------!

CONTAINS  !-- Module Subroutines

!-----------------------------------------------------------------------------!
! Allocation routines
!-----------------------------------------------------------------------------!

subroutine allocate_krige()
  implicit none

  allocate(vario    (maxstruct, nvario), &
           nstruct             (nvario), &
           rotmat(maxstruct,nvario,3,3)  )

end subroutine allocate_krige
!-----------------------------------------------------------------------------!

  subroutine grid_allocator(obj)
    implicit none
    type(grid2D), intent(inout)  :: obj
    allocate(obj%x(obj%n), obj%y(obj%n))
  end subroutine grid_allocator
  
!-----------------------------------------------------------------------------!
  
  subroutine kdata_allocator(obj)
    implicit none
    type(kdata), intent(inout)  :: obj
    allocate(obj%xy(2,obj%n), obj%v(obj%n))
  end subroutine kdata_allocator

!-----------------------------------------------------------------------------!
	
subroutine pd2kdata(pdset, kset, lay, aquitard)
  use pddata
  implicit none
  
  type(pddataset),intent(in)    :: pdset
	type(kdata),intent(inout)			:: kset
  integer,intent(in)            :: lay, aquitard
  integer                       :: i, p
  real,parameter                :: cutoff=-1.0
  
  if (aquitard==1) then
    kset%n = count(pdset%aqtv(:,lay) > cutoff)
  else
    kset%n = count(pdset%layv(:,lay) > cutoff)
  end if
  
  call kdata_allocator(kset)
  
  p = 1
  if (aquitard==1) then
    do i=1, pdset%npdat
      if (pdset%aqtv(i,lay) > cutoff) then
        kset%v(p) = pdset%aqtv(i,lay)
        kset%xy(:,p) = pdset%xy(:,i)
        p = p + 1
      end if
    end do
  else
    do i=1, pdset%npdat
      if (pdset%layv(i,lay) > cutoff) then
        kset%v(p) = pdset%layv(i,lay)
        kset%xy(:,p) = pdset%xy(:,i)
        p = p + 1
      end if
    end do
    !pdsetout%v_mean = sum(pdsetout%layv(:,1)) / pdsetout%npdat
  end if
  
  ! maintain mean
  kset%v_mean = pdset%v_mean    

end subroutine pd2kdata

!-----------------------------------------------------------------------------!
subroutine init_krige()
  implicit none
  ! To be called AFTER allocate_krige()
  
  integer       :: i,j

  ! Loop over variograms
  do i=1, nvario
    do j=1, nstruct(i)
      ! Sets anis1, anis2 ratios based on ranges aa1, aa2
      vario(j,i)%anis1 = vario(j,i)%a_hmin/vario(j,i)%a_hmax
      vario(j,i)%anis2 = 1 ! LS Hardcoded OUT a_vert/vario(j,i)%a_hmax
    
      ! Exist only in case we ever move to true 3d kriging
      vario(j,i)%ang2 = 0.0
      vario(j,i)%ang3 = 0.0
      
      ! Zero out nuggets on subsequent structures
      if (j > 1) vario(j,i)%nugget = 0.0
    
      ! Now the rotation matrices get initialized
      call setrot(i, j)
    end do
  end do

end subroutine init_krige
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
! Kriging routines
!-----------------------------------------------------------------------------!
subroutine spkrige(nobs, ngrid, gx, gy, obsx, obsy, obsv, ind, kriged_val)
    implicit none
!-----------------------------------------------------------------------------!
!                       Simple Point Kriging Subroutine
!                     ***********************************
!
! Requires: PSET, decomp, variogram, doolittle
!
! Arguments:
!   nobs (int) number of observations
!   ngrid (int) number of grid points
!   gx (real array) grid X values
!   gy (real array) grid Y values
!   itype (int) variogram shape (0-linear, 1-spherical, 2-exponential)
!   ind (int) matrix indicator to initialize
!   sill
!   ang1, ang2, ang3 (real) azimuth, dip plunge
!   anis1, anis2 (real) First & second anisotropy ratios
!   maxrot (int) maximum number of rotation matrices dimensioned
!   rotmat (int) rotation matrices
!   a
!   nugget
!   kriged_val (real array) kriging results
!-----------------------------------------------------------------------------!

    integer                  :: maxpt, ntot, i, g, n, ier, j, ind
    integer, allocatable     :: id(:)
    integer, intent(in)      :: nobs, ngrid
    real                     :: value, zkrig
    real, allocatable        :: p(:), r(:), w(:)
    real, intent(in)         :: gx(ngrid), gy(ngrid), obsx(nobs), obsy(nobs), &
                                obsv(nobs)
    real, intent(inout)      :: kriged_val(ngrid)
    real, parameter          :: nodataval = -999.0

    ntot = nobs + 1              ! +1 for ordinary
    maxpt = ntot * (ntot + 1) / 2
    allocate(p(maxpt),r(ntot),w(ntot),id(ntot))
    p = 0.d0
    id = 0
    call PSET(ind,p,obsx,obsy,id,nobs,ntot)
    call decomp(p,id,ntot, ier)

    if (ier == 1) then
      write(*,'(2x,a)') 'Error - Kriging Matrix Decomposition Failure.'
      write(*,'(2x,a)') 'This error is usually caused by co-located values.'
      stop
    end if

    ! Loop over each grid point
    do i = 1, ngrid
      r = 0.d0
      w = 0.d0
      r(ntot) = 1.0d0              ! commented = simple, uncommented = ordinary
      ! Loop over observations
      do n = 1, nobs
        do j=1, nstruct(ind)
          call variogram(ind,j,gx(i),gy(i),0.0d0,obsx(n),obsy(n),0.0d0,value)
          r(n) = r(n) + value
        end do
      end do
      call doolittle(r,p,id,w,ntot)
      value = 0.0d0
      zkrig = 0.0d0
      do n=1, nobs
        zkrig = zkrig + w(n) * obsv(n)
      end do
      zkrig = zkrig + value
      kriged_val(i) = zkrig
    end do

    deallocate(p, r, w, id)

end subroutine spkrige
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine spkrige_tree(ngrid, gx, gy, obstree, nobs, obsval, ind, kriged_val, nnear)
    use kdtree2_module
    implicit none
!-----------------------------------------------------------------------------!
!                       Simple Point Kriging Subroutine
!                     ***********************************
! Requires: PSET, decomp, variogram, doolittle
!-----------------------------------------------------------------------------!

    integer                  :: maxpt, ntot, i, j, g, n, &
                                usenear, ier
    integer, allocatable     :: id(:)
    integer, intent(in)      :: ngrid, nnear, nobs, ind
    real                     :: value, zkrig, queryxy(2), &
                                obsx(nnear), obsy(nnear)
    real, allocatable        :: p(:), r(:), w(:)
    real, intent(in)         :: gx(ngrid), gy(ngrid), obsval(nobs)
    real, intent(inout)      :: kriged_val(ngrid)

    type(kdtree2_result)     :: closeobs(nnear)
    type(kdtree2), intent(in), pointer :: obstree

    ! If there's not many datapoints
    if (obstree%N < nnear) then
      write(*,'(4x,a)') 'Well Data < Kriging Tree Search. Adjusting...'
      usenear = obstree%N
    else
      usenear = nnear
    end if

    ntot = usenear + 1              ! +1 for ordinary
    maxpt = ntot * (ntot + 1) / 2
    allocate(p(maxpt),r(ntot),w(ntot),id(ntot))

    ! Loop over each grid point
    do i = 1, ngrid
      ! Get nnear closest observation points
      queryxy = (/ gx(i), gy(i) /)
      call kdtree2_n_nearest(obstree, queryxy, usenear, closeobs)
      do j=1, usenear
        obsx(j) = obstree%THE_DATA(1,closeobs(j)%idx)
        obsy(j) = obstree%THE_DATA(2,closeobs(j)%idx)
      end do

      ! Reset values, matrices for new grid point
      p = 0.d0
      id = 0
      call PSET(ind,p,obsx,obsy,id,usenear,ntot)
      call decomp(p,id,ntot, ier)

      if (ier == 1) then
        write(*,'(2x,a)') 'Error - Kriging Matrix Decomposition Failure.'
        write(*,'(2x,a)') 'This error is usually caused by co-located wells.'
        stop
      end if
      r = 0.d0
      w = 0.d0
      r(ntot) = 1.0d0              ! commented = simple, uncommented = ordinary

      ! Loop over observations
      do n = 1, usenear
        do j=1, nstruct(ind)
          call variogram(ind,j,gx(i),gy(i),0.0d0,obsx(n),obsy(n),0.0d0,value)
          r(n) = r(n) + value
        end do
      end do
      call doolittle(r,p,id,w,ntot)

      ! Correct for negative weights, if necessary
      if (minval(w(1:usenear)) < 0.0d0) then
        call ok_negative_correction(w, usenear)
      end if

      value = 0.0d0
      zkrig = 0.0d0
      do n=1, usenear
        zkrig = zkrig + w(n) * obsval(closeobs(n)%idx)
      end do
      zkrig = zkrig + value
      kriged_val(i) = zkrig
    end do

    deallocate(p, r, w, id)

end subroutine spkrige_tree
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine spkrige_tree_group(ngrid, gx, gy, obstree, nobs, obsval, ind, &
                              kriged_val, nnear, ggroup, targetgroup)
    use kdtree2_module
    implicit none
!-----------------------------------------------------------------------------!
!                       Simple Point Kriging Subroutine
!                     ***********************************
!
! Groups are set as CHARACTERS(30)
!-----------------------------------------------------------------------------!

    integer                  :: maxpt, ntot, i, j, g, n, &
                                usenear, ier
    integer, allocatable     :: id(:)
    integer, intent(in)      :: ngrid, nnear, nobs, ind
    real                     :: value, zkrig, queryxy(2), &
                                obsx(nnear), obsy(nnear)
    real, allocatable        :: p(:), r(:), w(:)
    real, intent(in)         :: gx(ngrid), gy(ngrid), obsval(nobs)
    real, intent(inout)      :: kriged_val(ngrid)
    character(30),intent(in) :: ggroup(ngrid), targetgroup

    type(kdtree2_result)     :: closeobs(nnear)
    type(kdtree2), intent(in), pointer :: obstree

    ! If there's not many datapoints
    if (obstree%N < nnear) then
      ! Warning moved outside routine
      !write(*,'(4x,2a)') 'Well Data < Kriging Tree Search in group: ', targetgroup
      usenear = obstree%N
    else
      usenear = nnear
    end if

    ntot = usenear + 1              ! +1 for ordinary
    maxpt = ntot * (ntot + 1) / 2
    allocate(p(maxpt),r(ntot),w(ntot),id(ntot))

    ! Loop over each grid point
    do i = 1, ngrid

      ! Skip if grid point not in group
      if (ggroup(i) /= targetgroup) then
        !kriged_val(i) = nodataval
        cycle
      end if

      ! Get nnear closest observation points
      queryxy = (/ gx(i), gy(i) /)
      call kdtree2_n_nearest(obstree, queryxy, usenear, closeobs)
      do j=1, usenear
        obsx(j) = obstree%THE_DATA(1,closeobs(j)%idx)
        obsy(j) = obstree%THE_DATA(2,closeobs(j)%idx)
      end do

      ! Reset values, matricies for new grid point
      p = 0.d0
      id = 0
      call PSET(ind,p,obsx,obsy,id,usenear,ntot)
      call decomp(p,id,ntot, ier)

      if (ier == 1) then
        write(*,'(2x,a)') 'Error - Kriging Matrix Decomposition Failure.'
        write(*,'(2x,a)') 'This error is usually caused by co-located wells.'
        stop
      end if
      r = 0.d0
      w = 0.d0
      r(ntot) = 1.0d0              ! commented = simple, uncommented = ordinary

      ! Loop over observations
      do n = 1, usenear
        do j=1, nstruct(ind)
          call variogram(ind,j,gx(i),gy(i),0.0d0,obsx(n),obsy(n),0.0d0,value)
          r(n) = r(n) + value
        end do
      end do
      call doolittle(r,p,id,w,ntot)

      ! Correct for negative weights, if necessary
      if (minval(w(1:usenear)) < 0.0d0) then
        call ok_negative_correction(w, usenear)
      end if

      value = 0.0d0
      zkrig = 0.0d0
      do n=1, usenear
        zkrig = zkrig + w(n) * obsval(closeobs(n)%idx)
      end do
      zkrig = zkrig + value
      kriged_val(i) = zkrig
    end do

    deallocate(p, r, w, id)

end subroutine spkrige_tree_group
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine spkrige_group(nobs, ngrid, gx, gy, obsx, obsy, obsv, ind, kriged_val, &
                   ggroup, targetgroup)
    implicit none
!-----------------------------------------------------------------------------!
!                       Simple Point Kriging Subroutine
!                     ***********************************
!
! Groups are set as INTEGERS
!-----------------------------------------------------------------------------!

    integer                  :: maxpt, ntot, i, j, g, n, ier
    integer, allocatable     :: id(:)
    integer, intent(in)      :: nobs, ngrid, ind
    real                     :: value, zkrig
    real, allocatable        :: p(:), r(:), w(:)
    real, intent(in)         :: gx(ngrid), gy(ngrid), obsx(nobs), obsy(nobs), &
                                obsv(nobs)
    real, intent(inout)      :: kriged_val(ngrid)
    real, parameter          :: nodataval = -999.0
    ! Arguments for grouped kriging
    integer,intent(in)       :: ggroup(ngrid), targetgroup

    ntot = nobs + 1              ! +1 for ordinary
    maxpt = ntot * (ntot + 1) / 2
    allocate(p(maxpt),r(ntot),w(ntot),id(ntot))
    id = 0
    call PSET(ind,p,obsx,obsy,id,nobs,ntot)
    call decomp(p,id,ntot, ier)

    if (ier == 1) then
      write(*,'(2x,a)') 'Error - Kriging Matrix Decomposition Failure.'
      write(*,'(2x,a)') 'This error is usually caused by co-located pilot points.'
      stop
    end if

    ! Loop over each grid point
    do i = 1, ngrid
      ! Skip if grid point not in group
      if (ggroup(i) /= targetgroup) then
        !kriged_val(i) = nodataval
        cycle
      end if

      r = 0.d0
      w = 0.d0
      r(ntot) = 1.0d0              ! commented = simple, uncommented = ordinary
      ! Loop over observations
      do n = 1, nobs
        do j=1, nstruct(ind)
          call variogram(ind,j,gx(i),gy(i),0.0d0,obsx(n),obsy(n),0.0d0,value)
          r(n) = r(n) + value
        end do
      end do
      call doolittle(r,p,id,w,ntot)

      ! Correct for negative weights, if necessary
      if (minval(w(1:nobs)) < 0.0d0) then
        call ok_negative_correction(w, nobs)
      end if

      value = 0.0d0
      zkrig = 0.0d0
      do n=1, nobs
        zkrig = zkrig + w(n) * obsv(n)
      end do
      zkrig = zkrig + value
      kriged_val(i) = zkrig
    end do

    deallocate(p, r, w, id)

end subroutine
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine spcokrige(dsets, ndsets, grid, kriged_val)
    use kdtree2_module
    implicit none
!-----------------------------------------------------------------------------!
!                       Standardized Ordinary Cokriging
!                     ***********************************
!
! Requires: PSET, decomp, variogram, doolittle
!
!
!-----------------------------------------------------------------------------!
    integer, intent(in)      :: ndsets
    type(grid2D), intent(in) :: grid
    type(kdata),  intent(in) :: dsets(ndsets)
    real, intent(inout)      :: kriged_val(grid%n)
    integer                  :: maxpt, npoints, ntot, i, j, g, n,nn, ier, ind, rcount, d
    integer, allocatable     :: id(:)
    real                     :: value, zkrig
    real, allocatable        :: p(:), r(:), w(:), obsv(:)
    
    npoints = 0 + sum(dsets(:)%n) !(dsets(i)%n, i=1,ndsets)
    ntot = npoints + 1                         ! +1 for unbias
    maxpt = ntot * (ntot + 1) / 2
    allocate(p(maxpt),r(ntot),w(ntot),id(ntot),obsv(ntot))
    id = 0
    call pset_extreme(dsets, ndsets, ntot, p, id, maxpt)
    call decomp(p,id,ntot, ier)

    if (ier == 1) then
      write(*,'(2x,a)') 'Error - Kriging Matrix Decomposition Failure.'
      write(*,'(2x,a)') 'This error is usually caused by co-located values.'
      stop
    end if
    
    ! Pre-build observation dataset (constant over points)
    obsv = 0.0d0
    obsv(1:dsets(1)%n) = dsets(1)%v
    do d = 2, ndsets
      rcount = sum(dsets(1:(d-1))%n)
      ! Standardize if secondary data
      obsv(rcount+1:rcount+dsets(d)%n) = dsets(d)%v - dsets(d)%v_mean + dsets(1)%v_mean
    end do

    ! Loop over each grid point
    do i = 1, grid%n
      r = 0.d0
      w = 0.d0
      r(ntot) = 1.0d0              ! commented = simple, uncommented = ordinary
      rcount = 1
      ind = 1
      ! Loop over datasets
      do d = 1, ndsets
        if (d>1) ind = ind + ndsets
        ! loop over data points
        do n=1, dsets(d)%n
          do j=1, nstruct(ind)
            call variogram(ind,j,grid%x(i),grid%y(i),0.0d0,dsets(d)%xy(1,n),dsets(d)%xy(2,n),0.0d0,value)
            r(rcount) = r(rcount) + value
          end do
          rcount = rcount + 1
        end do
      end do
      call doolittle(r,p,id,w,ntot)
      kriged_val(i) = DOT_PRODUCT(w,obsv)
    end do

    deallocate(p, r, w, id, obsv)

end subroutine spcokrige
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine spcokrige_group(dsets, ndsets, ngrid, gx, gy, kriged_val, lay, &
													 aquitard, ggroup, targetgroup)
		use pddata
    use kdtree2_module
    implicit none
!-----------------------------------------------------------------------------!
!                       Standardized Ordinary Cokriging
!                     ***********************************
!
! Requires: PSET, decomp, variogram, doolittle
!
!
!-----------------------------------------------------------------------------!
    integer, intent(in)      :: ndsets, ngrid, lay, aquitard
    type(kdata), intent(in)  :: dsets(ndsets)
		real,intent(in)          :: gx(ngrid), gy(ngrid)
    real, intent(inout)      :: kriged_val(ngrid)
    integer                  :: maxpt, npoints, ntot, i, j, g, n,nn, ier, ind, rcount, d
    integer, allocatable     :: id(:)
    real                     :: value, zkrig
    real, allocatable        :: p(:), r(:), w(:), obsv(:)
		
		! Group arguments
		character(30),intent(in),optional :: ggroup(ngrid), targetgroup

		npoints = sum(dsets(:)%n)
    ntot = npoints + 1                         ! +1 for unbias
    maxpt = ntot * (ntot + 1) / 2
    allocate(p(maxpt),r(ntot),w(ntot),id(ntot),obsv(ntot))
    id = 0
		
		call pset_extreme(dsets, ndsets, ntot, p, id, maxpt)
    call decomp(p,id,ntot, ier)
    
    ! Pre-build observation dataset (constant over points)
    obsv = 0.0d0
		
		! Build observation array (obsv)
		rcount = 0
    do d = 1, ndsets
		  if (d>1) then
        ! Standardize if secondary data
        obsv(rcount+1:rcount+dsets(d)%n) = dsets(d)%v - dsets(d)%v_mean + dsets(1)%v_mean
			else
			  obsv(rcount+1:rcount+dsets(d)%n) = dsets(d)%v
			end if
      rcount = rcount + dsets(d)%n
    end do

    if (ier == 1) then
      write(*,'(2x,a)') 'Error - Kriging Matrix Decomposition Failure.'
      write(*,'(2x,a)') 'This error is usually caused by co-located values.'
      stop
    end if

    ! Loop over each grid point
    do i = 1, ngrid
			if (present(targetgroup)) then
        ! Skip if grid point not in group
        if (ggroup(i) /= targetgroup) then
          !kriged_val(i) = nodataval
          cycle
        end if
      end if

      r = 0.d0
      w = 0.d0
      r(ntot) = 1.0d0              ! commented = simple, uncommented = ordinary
      rcount = 1
      ind = 1
      ! Loop over datasets
      do d = 1, ndsets
      ! loop over data points
				ind = ind + ndsets
        do n=1, dsets(d)%n
          do j=1, nstruct(ind)
            call variogram(ind,j,gx(i),gy(i),0.0d0,dsets(d)%xy(1,n),dsets(d)%xy(2,n),0.0d0,value)
            r(rcount) = r(rcount) + value
          end do
          rcount = rcount + 1
        end do
      end do
      call doolittle(r,p,id,w,ntot)
      kriged_val(i) = DOT_PRODUCT(w,obsv)
    end do

    deallocate(p, r, w, id, obsv)

end subroutine spcokrige_group
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine spcokrige_tree_group(dsets, ndsets, ngrid, gx, gy, kriged_val, lay, &
                                dtrees, knnear, aquitard, &
                                nobs, obsval, obsmean, ggroup, targetgroup)
    use kdtree2_module
    use pddata
    implicit none
!-----------------------------------------------------------------------------!
!                       Standardized Ordinary Cokriging
!                     ***********************************
!
! Requires: PSET, decomp, variogram, doolittle
!
!
!-----------------------------------------------------------------------------!
    integer, intent(in)      :: ndsets, ngrid, nobs, aquitard, lay
    type(pddataset),intent(in) :: dsets(ndsets-1)
    type(kdata)              :: nearsets(ndsets)
    real,intent(in)          :: gx(ngrid), gy(ngrid), obsval(nobs),obsmean
    real, intent(inout)      :: kriged_val(ngrid)
    integer                  :: maxpt, npoints, ntot, i, j, k, g, n,nn, ier, ind, rcount, nearcount, d
    integer                  :: maxnear(ndsets)
    integer, allocatable     :: id(:)
    real                     :: value, zkrig, queryxy(2)
    real, allocatable        :: p(:), r(:), w(:), obsv(:), obsdists(:)
    
    ! Tree/Group arguments
    integer,intent(in)                :: knnear(ndsets)
    type(kdtree2_result),allocatable  :: closeobs(:)
    type(tree_array), intent(in)      :: dtrees(ndsets)
    character(30),intent(in),optional :: ggroup(ngrid), targetgroup
    
    ! Warn if less data than expected
    maxnear = knnear
    do i=1, ndsets
      if (dtrees(i)%tp%n < knnear(i)) then
        write(*,'(4x,a,i3,a)') 'Dataset' ,i,' - Well Data < Kriging Tree Search. Adjusting...'
        maxnear(i) = dtrees(i)%tp%n
      end if
			!TODO use kdata_allocator(obj)
      allocate(nearsets(i)%xy(2,maxnear(i)), & 
               nearsets(i)%v(maxnear(i))        )
    end do
    !npoints = sum(maxnear)

    do i=1, ngrid
      if (present(targetgroup)) then
        ! Skip if grid point not in group
        if (ggroup(i) /= targetgroup) then
          !kriged_val(i) = nodataval
          cycle
        end if
      end if
      
      ! Get nnear closest observation points
      !obsv = 0.0d0
      queryxy = (/ gx(i), gy(i) /)
      do j=1, ndsets
        ! Gather data
        allocate(closeobs(maxnear(j)), obsdists(maxnear(j)))
        call kdtree2_n_nearest(dtrees(j)%tp, queryxy, maxnear(j), closeobs)
        !nearsets(j)%xy(:,1:maxnear(j)) = dtrees(j)%tp%THE_DATA(:,closeobs(1:maxnear(j))%idx)
        !nearsets(j)%n = maxnear(j)
        nearcount     = 0
        do k=1, maxnear(j)
          if (j==1) then
            nearcount = nearcount + 1
            nearsets(j)%v(nearcount) = obsval(closeobs(k)%idx)
            !obsv         (nearcount) = obsval(closeobs(k)%idx)
          else if (sqrt(closeobs(k)%dis) > max_secondary_dist) then
            ! Skip
            cycle
          else if (aquitard==1) then
            nearcount = nearcount + 1
            nearsets(j)%v(nearcount) = dsets(j-1)%aqtv(closeobs(k)%idx,lay)
          else
            nearcount = nearcount + 1
            nearsets(j)%v(nearcount) = dsets(j-1)%layv(closeobs(k)%idx,lay)
          end if
          nearsets(j)%xy(:,nearcount) = dtrees(j)%tp%THE_DATA(:,closeobs(nearcount)%idx)
        end do
        nearsets(j)%n = nearcount
        ! Move secondary data to observation array
        !if (j>1) then
        !  !rcount = sum(maxnear(1:(j-1)))
        !  rcount = sum(nearsets(1:(j-1))%n)
        !  ! Standardize if secondary data
        !  !obsv(rcount+1:rcount+maxnear(j)) = nearsets(j)%v - dsets(j)%v_mean + dsets(1)%v_mean
        !  !obsv(rcount+1:rcount+maxnear(j)) = nearsets(j)%v - dsets(j-1)%v_mean + obsmean
        !  obsv(rcount+1:rcount+nearcount) = nearsets(j)%v(1:nearcount) - dsets(j-1)%v_mean + obsmean
        !end if
        deallocate(closeobs, obsdists)
      end do
      
      npoints = sum(nearsets(:)%n)
      ntot = npoints + 1                         ! +1 for unbias
      maxpt = ntot * (ntot + 1) / 2
      allocate(p(maxpt),r(ntot),w(ntot),id(ntot),obsv(ntot))
      
      ! Fill obsv
      obsv = 0.0d0
      rcount = 0
      do j=1, ndsets
        if (j==1) then
          obsv(rcount+1:rcount+nearsets(j)%n) = nearsets(j)%v(1:nearsets(j)%n)
        else  !(j>1)
          ! Standardize secondary data
          obsv(rcount+1:rcount+nearsets(j)%n) = nearsets(j)%v(1:nearsets(j)%n) - dsets(j-1)%v_mean + obsmean
        end if
        rcount = rcount + nearsets(j)%n
      end do
      
      ! Build correlation matrix
      id = 0
      call pset_extreme(nearsets, ndsets, ntot, p, id, maxpt)
      call decomp(p,id,ntot, ier)
      
      if (ier == 1) then
        write(*,'(2x,a)') 'Error - Kriging Matrix Decomposition Failure.'
        write(*,'(2x,a)') 'This error is usually caused by co-located wells.'
        stop
      end if

      r = 0.d0
      w = 0.d0
      r(ntot) = 1.0d0
      rcount = 1
      ind = 1
      ! Loop over datasets
      do d = 1, ndsets
        if (d>1) ind = ind + ndsets
        ! loop over data points
        do n=1, nearsets(d)%n
          do j=1, nstruct(ind)
            call variogram(ind,j,gx(i),gy(i),0.0d0,nearsets(d)%xy(1,n),nearsets(d)%xy(2,n),0.0d0,value)
            r(rcount) = r(rcount) + value
          end do
          rcount = rcount + 1
        end do
      end do
      call doolittle(r,p,id,w,ntot)
      kriged_val(i) = DOT_PRODUCT(w,obsv)
      deallocate(p, r, w, id, obsv)
    end do


end subroutine spcokrige_tree_group
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
  subroutine ok_negative_correction(wts, nwts)
    implicit none
!-----------------------------------------------------------------------------!
! Use for correcting issues with negative weights, presumably due to the
! screen effect. Follows the advice of Deutsch (1995) "Correcting for Negative
! Weights in Ordinary Kriging" and assigns negative weights to 0 and normalizes
! the weight sum back to one
!-----------------------------------------------------------------------------!

    real, intent(inout)    :: wts(nwts)
    integer, intent(in)    :: nwts
    integer                :: i
    real                   :: wtsum

    ! Loop over weights, correct negative to positive, add up sum
    wtsum = 0.0
    do i=1, nwts
      if (wts(i) < 0.0d0) wts(i) = 0.0d0
      wtsum = wtsum + wts(i)
    end do

    ! Normalize weight sum back to unity
    do i=1, nwts
      wts(i) = wts(i) / wtsum
    end do

  end subroutine ok_negative_correction
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
      subroutine setrot(i,j)
      implicit none
!
!              Sets up an Anisotropic Rotation Matrix
!              **************************************
!
! Sets up the matrix to transform cartesian coordinates to coordinates
! accounting for angles and anisotropy (see manual for a detailed
! definition):
!
!
! INPUT PARAMETERS: (moved to module definition)
!
!   ang1             Azimuth angle for principal direction
!   ang2             Dip angle for principal direction
!   ang3             Third rotation angle
!   anis1            First anisotropy ratio
!   anis2            Second anisotropy ratio
!   ind              The matrix indicator to initialize
!   MAXROT           The maximum number of rotation matrices dimensioned
!   rotmat           The rotation matrices
!
! Author: C. Deutsch                                Date: September 1989
! Modified for T2P by Leland Scantlebury
!-----------------------------------------------------------------------
      integer, intent(in)       :: i, j
      real*8, parameter         :: DEG2RAD = 3.141592654/180.0, &
                                   EPSLON = 1.e-10
      real*8                    :: alpha, beta, theta,&
                                   sina, sinb, sint, &
                                   cosa, cosb, cost, &
                                   afac1, afac2
!
! Converts the input angles to three angles which make more
!  mathematical sense:
!
!         alpha   angle between the major axis of anisotropy and the
!                 E-W axis. Note: Counter clockwise is positive.
!         beta    angle between major axis and the horizontal plane.
!                 (The dip of the ellipsoid measured positive down)
!         theta   Angle of rotation of minor axis about the major axis
!                 of the ellipsoid.
!

      if(vario(j,i)%ang1.ge.0.0d0.and.vario(j,i)%ang1.lt.270.0d0) then
          alpha = (90.0d0   - vario(j,i)%ang1) * DEG2RAD
      else
          alpha = (450.0d0  - vario(j,i)%ang1) * DEG2RAD
      endif
      beta  = -1.0d0 * vario(j,i)%ang2 * DEG2RAD
      theta =          vario(j,i)%ang3 * DEG2RAD
!
! Get the required sines and cosines:
!
      sina = sin(alpha)
      sinb = sin(beta)
      sint = sin(theta)
      cosa = cos(alpha)
      cosb = cos(beta)
      cost = cos(theta)
!
! Construct the rotation matrix in the required memory:
!
      afac1 = 1.0d0 / max(vario(j,i)%anis1,EPSLON)
      afac2 = 1.0d0 / max(vario(j,i)%anis2,EPSLON)
      rotmat(j,i,1,1) =       (cosb * cosa)
      rotmat(j,i,1,2) =       (cosb * sina)
      rotmat(j,i,1,3) =       (-sinb)
      rotmat(j,i,2,1) = afac1*(-cost*sina + sint*sinb*cosa)
      rotmat(j,i,2,2) = afac1*(cost*cosa + sint*sinb*sina)
      rotmat(j,i,2,3) = afac1*( sint * cosb)
      rotmat(j,i,3,1) = afac2*(sint*sina + cost*sinb*cosa)
      rotmat(j,i,3,2) = afac2*(-sint*cosa + cost*sinb*sina)
      rotmat(j,i,3,3) = afac2*(cost * cosb)
!
! Return to calling program:
!
      return
    end subroutine setrot
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
    SUBROUTINE PSET (ind,P,X,Y,ID,N,ntot)
      implicit none
!***  FROM SUBROUTINE WRITTEN BY SKRIVAN AND KARLINGER AND SLIGHTLY
!***  MODIFIED BY W.D. GRUNDY
!***  CALCULATE COEFFICIENTS IN THE SMETRIC P MATRIX WHICH WILL
!***  CONSIST OF THE LOWER TRIANGULAR PORTION STORED BY COLUMNS.
!***  ID IS THE POINTER VECTOR GIVING THE LOCATION OF THE DIAGONAL
!***  ELEMENTS IN P.
! Multiple Structures implemented by LS
      integer, intent(in) :: ind
      integer    :: id,ntot,np1,N,MAXPT,i,j,js,jj
      real       :: p,x,y,x1,x2,y1,y2,value,value_sum
      DIMENSION  :: P(1),X(1),Y(1),ID(1)
!      SAVE
!      NTOT=N+1
      NP1=N+1
!***  ZERO OUT ELEMENTS OF P-MATRIX NEEDED FOR KRIGING
      MAXPT=NTOT*(NTOT+1)/2
      DO 10 I=1,MAXPT
      P(I)=0.
   10 CONTINUE
!***  CALCULATE THE POINTERS FOR ID
      ID(1)=1
      DO 20 I=2,NTOT
   20 ID(I)=ID(I-1)+NTOT+2-I
!***  CALCULATE THE COEFFICIENTS OF THE COVARIANCE OF POINTS I AND J
      DO 30 J=1,N
        JS=ID(J)-J
      DO 30 I=J,N
        x1=x(i)
        y1=y(i)
        x2=x(j)
        y2=y(j)
        value_sum=0.0d0
        do jj=1, nstruct(ind)
          call variogram(ind,jj,x1,y1,0.0d0,x2,y2,0.0d0,value)
          value_sum = value_sum + value
        end do
        P(JS+I)=value_sum
   30 CONTINUE
      DO 40 I=1,N
   40 P(ID(I)+NP1-I)=1
      RETURN
      END SUBROUTINE PSET
!-----------------------------------------------------------------------------!
      
  subroutine pset_extreme(dsets, ndsets, ntot, p, id, maxpt)
    implicit none
    ! Rewrite of PSET routine above for cokriging
    ! Author Leland Scantlebury
    
    ! Inputs
    integer,intent(in)     :: ndsets, ntot, maxpt
    integer,intent(inout)  :: id(ntot)
    type(kdata),intent(in) :: dsets(ndsets)
    real,intent(inout)     :: p(maxpt)
    
    ! Subroutine variables
    integer             :: i,ii,j,jj,n,ind,JS, rowcount
    integer,allocatable :: colcount(:)
    real                :: x1,x2,y1,y2,value,value_sum
    
    allocate(colcount(ndsets))
    
    colcount = 0
    do i=1, (ndsets-1)
      colcount(i+1) = colcount(i) + dsets(i)%n
    end do
    
!***  ZERO OUT ELEMENTS OF P-MATRIX NEEDED FOR KRIGING
    P=0.0
!***  CALCULATE THE POINTERS FOR ID
    ID(1)=1
    DO I=2,NTOT
      ID(I)=ID(I-1)+NTOT+2-I
    END DO
    ! Outer Loop over datasets (columns)
    do j=1, ndsets
      ! Inner loop over datasets (rows)
      rowcount = 0
      do i=j, ndsets
        ! Variogram index for section
        ind = j
        if (i/=j) then
          ind = ndsets + j*(i-1)  ! ndsets + j-1 + i-1 ??
        end if
        ! Outer loop over column points
        do jj=1, dsets(j)%n
          JS=ID(colcount(j)+jj)-(colcount(j)+jj)
          x2=dsets(j)%xy(1,jj)
          y2=dsets(j)%xy(2,jj)
          ! Inner loop over row points
          ii = jj
          if (i/=j) ii=1
          do while (ii <= dsets(i)%n)
            x1=dsets(i)%xy(1,ii)
            y1=dsets(i)%xy(2,ii)
            value_sum=0.0d0
            do n=1, nstruct(ind)
              call variogram(ind,n,x1,y1,0.0d0,x2,y2,0.0d0,value)
              value_sum = value_sum + value
            end do
            P(JS+rowcount+colcount(j)+ii)=value_sum
            ii = ii + 1
          end do
        end do
        rowcount = rowcount + (ii-1)
      end do
    end do
    ! Fill last row with unbias
    DO I=1,ntot-1
      P(ID(I)+ntot-I)=1
    END DO
      
    deallocate(colcount)
    RETURN

  end subroutine pset_extreme

!-----------------------------------------------------------------------------!
subroutine decomp(p,id,ntot,ier)
    implicit none
!*** DECOMPOSITION TO GET LOWER TRIANGULAR S MATRIX
    integer                :: id,ier,ks,jm1,k,js,is,is1,j,i,i1,ntot
    real                   :: p,sum,p1,p2
    dimension              :: p(1),id(1)

          IER=0
          P(1)=1./P(1)
          DO 65 J=2,NTOT
          DO 65 I=J,NTOT
            KS=ID(I)
            JM1=J-1
            SUM=0.
            DO 55 K=1,JM1
              JS=ID(K)
              IS=JS+I-K
              IS1=JS+J-K
              P1=P(IS)
              P2=P(IS1)
              IF (P1) 39,55,39
   39         IF (P2) 50,55,50
   50         SUM=SUM+P1*P2*P(JS)
   55         CONTINUE
            I1=ID(J)+I-J
            P(I1)=P(I1)-SUM
            IF (I-J) 65,70,65
   70       IF(P(KS)) 90,91,90
   90       P(KS)=1./P(KS)
   65       CONTINUE
          GOTO 92
   91     IER=1
   92     CONTINUE
      return
end subroutine decomp
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------
subroutine variogram(ind,str,x1,y1,z1,x2,y2,z2,value)
      implicit none
!      common block/vario/slope,sill,a
      integer,intent(in) :: ind, str
      real           :: x1,y1,z1,x2,y2,z2,value,h,tmp
      real,parameter :: epsilon=1.e-10

      h = sqdist(x1,y1,z1,x2,y2,z2,ind,str)

  ! Check for "zero" distance, return with maximum value
      if (h < epsilon) then
        value = vario(str,ind)%nugget + vario(str,ind)%sill
        return
      end if
  ! linear variogram
      if(vario(str,ind)%itype.eq.0) then
        h=sqrt(h)*vario(str,ind)%sill
        !write(11,*) 'Linear', h
      end if
  ! spherical variogram
      if(vario(str,ind)%itype.eq.1) then
        H=SQRT(h)
        tmp=h/vario(str,ind)%a_hmax
        if(tmp.ge.1.) then
          h=vario(str,ind)%sill
        else
          h=vario(str,ind)%sill*(1.5*tmp-0.5*tmp*tmp*tmp)
        endif
      endif
  ! exponential variogram
      if (vario(str,ind)%itype.eq.2) then
        H=SQRT(h)
        h=vario(str,ind)%sill*(1.d0-exp(-3.d0*h/vario(str,ind)%a_hmax))
      endif
      value=vario(str,ind)%sill-h

end subroutine variogram
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine doolittle(r,p,id,w,ntot)

    implicit none
          integer   :: id,ntot,i,is,js,im1,j,ii,iip1,iis
          real      :: r,p,w,sum
      dimension     :: r(1),p(1),w(1),id(1)
!***  START OF DOOLITTLE ***
!***  FORWARD SUBSTITUTION FOR SYSTEM S*W = R
          W(1)=R(1)*P(1)
          DO 130 I=2,NTOT
            IS=ID(I)
            SUM=0.
            IM1=I-1
            DO 120 J=1,IM1
              JS=ID(J)+I-J
  120       SUM=SUM+W(J)*P(JS)
  130     W(I)=(R(I)-SUM)*P(IS)
!***  BACKWARD SUBSTITUTION FOR SYSTEM T*V = W
!***  WHERE V IS STORED IN W
          DO 150 I=2,NTOT
            SUM=0.
            II=NTOT-I+1
            JS=ID(II)
            IIP1=II+1
            DO 140 J=IIP1,NTOT
              IIS=JS+J-II
  140       SUM=SUM+W(J)*P(IIS)
  150     W(II)=W(II)-SUM*P(JS)
!*** END OF DOOLITTLE ****
      return

end subroutine doolittle
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
real function sqdist(x1,y1,z1,x2,y2,z2,ind,str)
      implicit none
!
!    Squared Anisotropic Distance Calculation Given Matrix Indicator
!    ***************************************************************
!
! This routine calculates the anisotropic distance between two points
!  given the coordinates of each point and a definition of the
!  anisotropy.
!
! INPUT VARIABLES:
!   x1,y1,z1         Coordinates of first point
!   x2,y2,z2         Coordinates of second point
!   ind              The matrix indicator to initialize
!   MAXROT           The maximum number of rotation matrices dimensioned
!   rotmat           The rotation matrices
!
!
! OUTPUT VARIABLES:
!   sqdist           The squared distance accounting for the anisotropy
!                      and the rotation of coordinates (if any).
!
!
! Author: C. Deutsch                                Date: September 1989
      integer   :: i,ind, str
      real      :: dx,dy,dz,x1,y1,x2,y2,z1,z2,cont
!
! Compute component distance vectors and the squared distance:
!
      dx = x1 - x2
      dy = y1 - y2
      dz = z1 - z2
          sqdist = 0.0d0
      do i=1,3
        cont   = rotmat(str,ind,i,1) * dx &
               + rotmat(str,ind,i,2) * dy &
               + rotmat(str,ind,i,3) * dz
        sqdist = sqdist + cont * cont
      end do
      return

end function sqdist
!-----------------------------------------------------------------------------!

  END MODULE Krige