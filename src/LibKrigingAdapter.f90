! Texture2Par-side boundary for the native Fortran libkriging API.
!
! Texture2Par owns input parsing, categories, target/observation packing, and
! scattering packed results back to model arrays.  libkriging is intended to
! own covariance/variogram evaluation, neighbour search, the kriging solve,
! supported caching/reuse, and (in a later stage) simulation.
!
! The legacy Texture2Par dip/rotation convention may intentionally differ
! from libkriging's GSLIB-style convention.  This adapter does not attempt to
! reconcile those conventions; that mapping belongs to a later cutover stage.
module m_libkriging_adapter
  use m_interpolation_types, only: t_interp_request, t_interp_result
  use m_vario,              only: t_vgm
  use m_error_handler,      only: error_handler
  use kriging,              only: t_kriging
  use kriging_err,          only: kriging_clear_error, kriging_failed, kriging_last_error
  implicit none
  private

  public :: libkriging_adapter_smoke_test
  public :: libkriging_solve
  public :: libkriging_interpolate_many

contains

!-------------------------------------------------------------------------------------------------!
  ! Translate a backend-neutral batch request plus Texture2Par variograms into one native
  ! libkriging solve. OK and SK share this path; unbias/sk_mean are configuration, not classes.
  subroutine libkriging_solve(request, variograms, nstruct, unbias, sk_mean, result)
    type(t_interp_request), intent(in)  :: request
    type(t_vgm), intent(in)             :: variograms(:)
    integer, intent(in)                 :: nstruct, unbias
    real, intent(in)                    :: sk_mean
    type(t_interp_result), intent(out)  :: result

    type(t_kriging)        :: solver
    character(3)           :: vtype
    character(:), allocatable :: failure
    integer                :: i, ntarget
    logical                :: initialized, succeeded

    initialized = .false.
    succeeded = .false.
    ntarget = size(request%target_coordinates, 2)
    allocate(result%estimates(ntarget), result%variance(ntarget), result%target_index(ntarget))
    result%estimates = 0.0
    result%variance = 0.0
    result%target_index = request%target_index
    if (ntarget == 0) return

    if (request%dimension < 2 .or. request%dimension > 3) then
      call error_handler(4, opt_msg='OK_LIB supports only 2D and 3D interpolation')
    end if
    if (nstruct < 1 .or. nstruct > size(variograms)) then
      call error_handler(4, opt_msg='OK_LIB received an invalid variogram structure count')
    end if

    call kriging_clear_error()
    call solver%initialize(ndim=request%dimension, nvar=1, unbias=unbias, &
                           anisotropic_search=request%search%anisotropic_search, verbose=.false.)
    initialized = .true.
    if (kriging_failed()) goto 900

    if (unbias == 0) then
      call solver%set_obs(1, request%observation_coordinates, request%observation_values, &
                          request%observation_variance, sk_mean=sk_mean)
    else
      call solver%set_obs(1, request%observation_coordinates, request%observation_values, &
                          request%observation_variance)
    end if
    if (kriging_failed()) goto 900

    call solver%set_grid(request%target_coordinates)
    if (kriging_failed()) goto 900

    do i = 1, nstruct
      vtype = libkriging_variogram_type(variograms(i)%itype)
      call solver%set_vgm(1, 1, vtype=vtype, nugget=variograms(i)%nugget, &
                          sill=variograms(i)%sill, a_major=variograms(i)%a_hmax, &
                          a_minor1=variograms(i)%a_hmin, a_minor2=variograms(i)%a_vert, &
                          azimuth=variograms(i)%ang1, dip=variograms(i)%ang2, &
                          plunge=variograms(i)%ang3)
      if (kriging_failed()) goto 900
    end do

    ! libkriging uses its GSLIB-style angle convention directly. Texture2Par's legacy
    ! rotation matrix is intentionally not passed into this backend.
    if (request%search%anisotropic_search) then
      call set_search_anisotropic(solver, request, variograms(1))
    else
      call set_search_isotropic(solver, request)
    end if
    if (kriging_failed()) goto 900

    call solver%solve(nthread=1)
    if (kriging_failed() .or. .not. solver%solved) goto 900

    result%estimates = solver%block%value(1,1,:)
    result%variance  = solver%block%variance(1,1,:)
    result%status = 0
    succeeded = .true.

900 continue
    if (.not. succeeded) failure = kriging_last_error()
    if (initialized) call solver%finalize()
    if (.not. succeeded) then
      if (len_trim(failure) == 0) failure = 'libkriging solve did not complete'
      call error_handler(4, opt_msg='OK_LIB: '//trim(failure))
    end if
  end subroutine libkriging_solve
!-------------------------------------------------------------------------------------------------!

  ! Interpolate many value columns over one packed observation/target geometry.
  ! value_matrix is (observation, value-column); estimates and variance are
  ! (target, value-column).
  subroutine libkriging_interpolate_many(request, variograms, nstruct, unbias, sk_mean, &
                                         value_matrix, estimates, variance)
    type(t_interp_request), intent(in)   :: request
    type(t_vgm), intent(in)              :: variograms(:)
    integer, intent(in)                  :: nstruct, unbias
    real, intent(in)                     :: sk_mean
    real, intent(in)                     :: value_matrix(:,:)
    real, allocatable, intent(out)       :: estimates(:,:), variance(:,:)

    type(t_kriging)        :: solver
    character(3)           :: vtype
    character(:), allocatable :: failure
    integer                :: column, i, nobs, ntarget, nvalues
    logical                :: initialized, succeeded

    initialized = .false.
    succeeded = .false.
    nobs = size(request%observation_coordinates, 2)
    ntarget = size(request%target_coordinates, 2)
    nvalues = size(value_matrix, 2)
    if (size(value_matrix, 1) /= nobs) then
      call error_handler(4, opt_msg='OK_LIB many-value observation geometry/value shape mismatch')
    end if

    allocate(estimates(ntarget, nvalues), variance(ntarget, nvalues))
    estimates = 0.0
    variance = 0.0
    if (ntarget == 0 .or. nvalues == 0) return

    if (request%dimension < 2 .or. request%dimension > 3) then
      call error_handler(4, opt_msg='OK_LIB supports only 2D and 3D interpolation')
    end if
    if (nstruct < 1 .or. nstruct > size(variograms)) then
      call error_handler(4, opt_msg='OK_LIB received an invalid variogram structure count')
    end if

    call kriging_clear_error()
    ! store_weight retains the first solve's neighbours, weights, and variances
    ! for value columns sharing the same observation and target geometry.
    call solver%initialize(ndim=request%dimension, nvar=1, unbias=unbias, &
                           anisotropic_search=request%search%anisotropic_search, &
                           store_weight=nvalues > 1, verbose=.false.)
    initialized = .true.
    if (kriging_failed()) goto 900

    if (unbias == 0) then
      call solver%set_obs(1, request%observation_coordinates, value_matrix(:,1), &
                          request%observation_variance, sk_mean=sk_mean)
    else
      call solver%set_obs(1, request%observation_coordinates, value_matrix(:,1), &
                          request%observation_variance)
    end if
    if (kriging_failed()) goto 900

    call solver%set_grid(request%target_coordinates)
    if (kriging_failed()) goto 900

    do i = 1, nstruct
      vtype = libkriging_variogram_type(variograms(i)%itype)
      call solver%set_vgm(1, 1, vtype=vtype, nugget=variograms(i)%nugget, &
                          sill=variograms(i)%sill, a_major=variograms(i)%a_hmax, &
                          a_minor1=variograms(i)%a_hmin, a_minor2=variograms(i)%a_vert, &
                          azimuth=variograms(i)%ang1, dip=variograms(i)%ang2, &
                          plunge=variograms(i)%ang3)
      if (kriging_failed()) goto 900
    end do

    if (request%search%anisotropic_search) then
      call set_search_anisotropic(solver, request, variograms(1))
    else
      call set_search_isotropic(solver, request)
    end if
    if (kriging_failed()) goto 900

    call solver%solve(nthread=1)
    if (kriging_failed() .or. .not. solver%solved) goto 900
    estimates(:,1) = solver%block%value(1,1,:)
    variance(:,1) = solver%block%variance(1,1,:)

    do column = 2, nvalues
      ! update_obs_value preserves observation variance and activates reuse of
      ! the populated in-memory weights. External Texture2Par READ/WRITE weight
      ! files remain a separate compatibility task.
      call solver%update_obs_value(1, value_matrix(:,column))
      if (kriging_failed()) goto 900
      call solver%solve(nthread=1)
      if (kriging_failed() .or. .not. solver%solved) goto 900
      estimates(:,column) = solver%block%value(1,1,:)
      variance(:,column) = solver%block%variance(1,1,:)
    end do
    succeeded = .true.

900 continue
    if (.not. succeeded) failure = kriging_last_error()
    if (initialized) call solver%finalize()
    if (.not. succeeded) then
      if (len_trim(failure) == 0) failure = 'libkriging many-value solve did not complete'
      call error_handler(4, opt_msg='OK_LIB: '//trim(failure))
    end if
  end subroutine libkriging_interpolate_many
!-------------------------------------------------------------------------------------------------!

  subroutine set_search_isotropic(solver, request)
    type(t_kriging), intent(inout)       :: solver
    type(t_interp_request), intent(in)   :: request

    ! TODO: expose libkriging sector search after Texture2Par has sector-count controls.
    if (valid_maxdist(request%search%max_primary_distance)) then
      call solver%set_search(1, nmax=request%search%n_nearest_primary, &
                             maxdist=request%search%max_primary_distance, sector_search=.false.)
    else
      call solver%set_search(1, nmax=request%search%n_nearest_primary, sector_search=.false.)
    end if
  end subroutine set_search_isotropic
!-------------------------------------------------------------------------------------------------!

  subroutine set_search_anisotropic(solver, request, variogram)
    type(t_kriging), intent(inout)       :: solver
    type(t_interp_request), intent(in)   :: request
    type(t_vgm), intent(in)              :: variogram

    if (valid_maxdist(request%search%max_primary_distance)) then
      call solver%set_search(1, anis1=variogram%anis1, anis2=variogram%anis2, &
                             azimuth=variogram%ang1, dip=variogram%ang2, &
                             plunge=variogram%ang3, nmax=request%search%n_nearest_primary, &
                             maxdist=request%search%max_primary_distance, sector_search=.false.)
    else
      call solver%set_search(1, anis1=variogram%anis1, anis2=variogram%anis2, &
                             azimuth=variogram%ang1, dip=variogram%ang2, &
                             plunge=variogram%ang3, nmax=request%search%n_nearest_primary, &
                             sector_search=.false.)
    end if
  end subroutine set_search_anisotropic
!-------------------------------------------------------------------------------------------------!

  logical function valid_maxdist(maxdist)
    real, intent(in) :: maxdist
    valid_maxdist = maxdist > 0.0 .and. maxdist < sqrt(huge(maxdist))
  end function valid_maxdist
!-------------------------------------------------------------------------------------------------!

  function libkriging_variogram_type(itype) result(vtype)
    integer, intent(in) :: itype
    character(3)        :: vtype

    select case (itype)
      case (1); vtype = 'sph'
      case (2); vtype = 'exp'
      case (3); vtype = 'hol'
      case (4); vtype = 'gau'
      case (5); vtype = 'pow'
      case (6); vtype = 'cir'
      case (7); vtype = 'lin'
      case default
        call error_handler(4, opt_msg='OK_LIB received an unsupported variogram type')
    end select
  end function libkriging_variogram_type
!-------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------!
  ! Exercise a real native ordinary-kriging solve without entering Texture2Par's runtime path.
  ! Keeping this routine in the executable build forces the linker to resolve the libkriging API.
  subroutine libkriging_adapter_smoke_test(estimate, kriging_variance, succeeded)
    real,    intent(out), optional :: estimate, kriging_variance
    logical, intent(out), optional :: succeeded

    type(t_interp_request) :: request
    type(t_interp_result)  :: result
    type(t_kriging)        :: solver
    logical                :: initialized

    if (present(estimate))         estimate = 0.0
    if (present(kriging_variance)) kriging_variance = 0.0
    if (present(succeeded))        succeeded = .false.
    initialized = .false.

    ! Neutral Texture2Par storage and libkriging both use (dimension, point).
    request%dimension = 2
    allocate(request%observation_coordinates(2,4))
    request%observation_coordinates(:,1) = [0.0, 0.0]
    request%observation_coordinates(:,2) = [1.0, 0.0]
    request%observation_coordinates(:,3) = [0.0, 1.0]
    request%observation_coordinates(:,4) = [1.0, 1.0]
    request%observation_values   = [1.0, 2.0, 3.0, 4.0]
    request%observation_variance = [0.0, 0.0, 0.0, 0.0]
    allocate(request%target_coordinates(2,1))
    request%target_coordinates(:,1) = [0.5, 0.5]

    call kriging_clear_error()
    call solver%initialize(ndim=request%dimension, nvar=1, unbias=1, verbose=.false.)
    initialized = .true.
    if (kriging_failed()) goto 900

    call solver%set_obs(1, request%observation_coordinates, &
                        request%observation_values, request%observation_variance)
    if (kriging_failed()) goto 900

    call solver%set_grid(request%target_coordinates)
    if (kriging_failed()) goto 900

    call solver%set_vgm(1, 1, vtype='sph', nugget=0.0, sill=1.0, a_major=2.0)
    if (kriging_failed()) goto 900

    call solver%set_search(1, nmax=size(request%observation_values))
    if (kriging_failed()) goto 900

    call solver%solve(nthread=1)
    if (kriging_failed() .or. .not. solver%solved) goto 900

    allocate(result%estimates(1), result%variance(1))
    result%estimates(1) = solver%block%value(1,1,1)
    result%variance(1)  = solver%block%variance(1,1,1)
    if (present(estimate))         estimate = result%estimates(1)
    if (present(kriging_variance)) kriging_variance = result%variance(1)
    if (present(succeeded))        succeeded = .true.

900 continue
    if (initialized) call solver%finalize()
  end subroutine libkriging_adapter_smoke_test
!-------------------------------------------------------------------------------------------------!

end module m_libkriging_adapter
