! Texture2Par-side boundary for the native Fortran libkriging API.
!
! Texture2Par owns input parsing, categories, target/observation packing, and
! scattering packed results back to model arrays.  libkriging is intended to
! own covariance/variogram evaluation, neighbour search, the kriging solve,
! supported caching/reuse, and simulation.
!
! The legacy Texture2Par dip/rotation convention may intentionally differ
! from libkriging's GSLIB-style convention.  This adapter does not attempt to
! reconcile those conventions; that mapping belongs to a later cutover stage.
module m_libkriging_adapter
  use m_interpolation_types, only: t_interp_request, t_interp_result, &
                                   PREDICT_ESTIMATE, PREDICT_SGSIM
  use m_vario,              only: t_vgm
  use m_error_handler,      only: error_handler
  use kriging,              only: t_kriging
  use kriging_err,          only: kriging_clear_error, kriging_failed, kriging_last_error
  use utils,                only: set_seq, r8vec_normal_01
  implicit none
  private

  public :: libkriging_adapter_smoke_test
  public :: libkriging_solve
  public :: libkriging_interpolate_many

  type, public :: t_libkriging_variable
    real, allocatable :: coordinates(:,:)
    real, allocatable :: values(:)
    real, allocatable :: variance(:)
    integer           :: nmax = 0
    real              :: max_distance = -1.0
    real              :: sk_mean = 0.0
  end type t_libkriging_variable

  interface libkriging_solve
    module procedure libkriging_solve_single
    module procedure libkriging_solve_variables
  end interface libkriging_solve

contains

!-------------------------------------------------------------------------------------------------!
  ! Translate a backend-neutral batch request plus Texture2Par variograms into one native
  ! libkriging solve. OK and SK share this path; unbias/sk_mean are configuration, not classes.
  subroutine libkriging_solve_single(request, variograms, nstruct, unbias, sk_mean, result)
    type(t_interp_request), intent(in)  :: request
    type(t_vgm), intent(in)             :: variograms(:)
    integer, intent(in)                 :: nstruct, unbias
    real, intent(in)                    :: sk_mean
    type(t_interp_result), intent(out)  :: result

    type(t_libkriging_variable) :: variables(1)
    type(t_vgm), allocatable    :: variogram_matrix(:,:)
    integer                    :: nstruct_by_vgm(1), vgm_index(1,1)

    variables(1)%coordinates = request%observation_coordinates
    variables(1)%values = request%observation_values
    variables(1)%variance = request%observation_variance
    variables(1)%nmax = request%search%n_nearest_primary
    variables(1)%max_distance = request%search%max_primary_distance
    variables(1)%sk_mean = sk_mean
    allocate(variogram_matrix(size(variograms),1))
    variogram_matrix(:,1) = variograms
    nstruct_by_vgm(1) = nstruct
    vgm_index(1,1) = 1

    call libkriging_solve_variables(request, variables, variogram_matrix, &
                                    nstruct_by_vgm, vgm_index, unbias, .true., result)
  end subroutine libkriging_solve_single
!-------------------------------------------------------------------------------------------------!

  subroutine libkriging_solve_variables(request, variables, variograms, nstruct, &
                                        vgm_index, unbias, nonstd_ck_mean_correction, result)
    type(t_interp_request), intent(in)      :: request
    type(t_libkriging_variable), intent(in) :: variables(:)
    type(t_vgm), intent(in)                 :: variograms(:,:)
    integer, intent(in)                     :: nstruct(:), vgm_index(:,:), unbias
    logical, intent(in)                     :: nonstd_ck_mean_correction
    type(t_interp_result), intent(out)      :: result

    type(t_kriging)           :: solver
    character(3)              :: vtype
    character(:), allocatable :: failure
    integer                   :: i, isim, ivar, jvar, nvar, ntarget, vidx
    integer, allocatable      :: simulation_path(:)
    real, allocatable         :: simulation_sample(:,:,:), draw(:)
    logical                   :: initialized, succeeded

    initialized = .false.
    succeeded = .false.
    nvar = size(variables)
    ntarget = size(request%target_coordinates, 2)
    allocate(result%estimates(ntarget), result%variance(ntarget), result%target_index(ntarget))
    result%estimates = 0.0
    result%variance = 0.0
    result%target_index = request%target_index
    if (ntarget == 0) return

    if (request%dimension < 2 .or. request%dimension > 3) &
      call error_handler(4, opt_msg='libkriging supports only 2D and 3D interpolation')
    if (nvar < 1 .or. size(vgm_index,1) /= nvar .or. size(vgm_index,2) /= nvar) &
      call error_handler(4, opt_msg='libkriging received invalid cokriging variable metadata')

    do ivar = 1, nvar
      if (.not. allocated(variables(ivar)%coordinates) .or. &
          .not. allocated(variables(ivar)%values) .or. &
          .not. allocated(variables(ivar)%variance)) &
        call error_handler(4, opt_msg='libkriging received incomplete observation data')
      if (size(variables(ivar)%coordinates,1) /= request%dimension .or. &
          size(variables(ivar)%coordinates,2) /= size(variables(ivar)%values) .or. &
          size(variables(ivar)%variance) /= size(variables(ivar)%values)) &
        call error_handler(4, opt_msg='libkriging observation geometry/value shape mismatch')
      if (variables(ivar)%nmax < 1) &
        call error_handler(4, opt_msg='libkriging requires NNEAR to be at least one for every variable')
    end do

    call kriging_clear_error()
    select case (request%prediction_mode)
      case (PREDICT_ESTIMATE)
        call solver%initialize(ndim=request%dimension, nvar=nvar, unbias=unbias, &
                               std_ck=nvar == 1, &
                               nonstd_ck_mean_correction=nonstd_ck_mean_correction, &
                               anisotropic_search=request%search%anisotropic_search, verbose=.false.)
      case (PREDICT_SGSIM)
        if (request%simulation_count < 1) &
          call error_handler(4, opt_msg='libkriging SGSIM requires at least one realization')
        if (request%request_block_support) &
          call error_handler(1, opt_msg='BLOCK_MODE with SGSIM is not supported in this migration stage')
        call solver%initialize(ndim=request%dimension, nvar=nvar, unbias=unbias, &
                               nsim=request%simulation_count, seed=request%random_seed, &
                               std_ck=nvar == 1, &
                               nonstd_ck_mean_correction=nonstd_ck_mean_correction, &
                               anisotropic_search=request%search%anisotropic_search, verbose=.false.)
      case default
        call error_handler(4, opt_msg='libkriging received an invalid prediction mode')
    end select
    initialized = .true.
    if (kriging_failed()) goto 900

    ! Texture2Par selects libkriging's local CK mean correction or its historical
    ! global secondary standardization before entering this shared solve path.

    do ivar = 1, nvar
      if (unbias == 0) then
        call solver%set_obs(ivar, variables(ivar)%coordinates, variables(ivar)%values, &
                            variables(ivar)%variance, sk_mean=variables(ivar)%sk_mean)
      else
        call solver%set_obs(ivar, variables(ivar)%coordinates, variables(ivar)%values, &
                            variables(ivar)%variance)
      end if
      if (kriging_failed()) goto 900
    end do

    call set_solver_grid(solver, request)
    if (kriging_failed()) goto 900

    do ivar = 1, nvar
      do jvar = ivar, nvar
        vidx = vgm_index(ivar,jvar)
        if (vidx < 1 .or. vidx > size(nstruct) .or. vidx > size(variograms,2)) &
          call error_handler(4, opt_msg='libkriging could not map an auto/cross variogram')
        if (nstruct(vidx) < 1 .or. nstruct(vidx) > size(variograms,1)) &
          call error_handler(4, opt_msg='libkriging received an invalid variogram structure count')
        do i = 1, nstruct(vidx)
          vtype = libkriging_variogram_type(variograms(i,vidx)%itype)
          call solver%set_vgm(ivar, jvar, vtype=vtype, nugget=variograms(i,vidx)%nugget, &
                              sill=variograms(i,vidx)%sill, a_major=variograms(i,vidx)%a_hmax, &
                              a_minor1=variograms(i,vidx)%a_hmin, a_minor2=variograms(i,vidx)%a_vert, &
                              azimuth=variograms(i,vidx)%ang1, dip=variograms(i,vidx)%ang2, &
                              plunge=variograms(i,vidx)%ang3)
          if (kriging_failed()) goto 900
        end do
      end do
    end do

    if (request%prediction_mode == PREDICT_SGSIM) then
      ! Supplying the path and samples keeps simulation state in memory; the
      ! libkriging default otherwise writes diagnostic sgs_path/sample files.
      allocate(simulation_path(ntarget), &
               simulation_sample(request%simulation_count,nvar,ntarget), draw(ntarget))
      call set_seq(simulation_path, ntarget, .true.)
      do isim = 1, request%simulation_count
        do ivar = 1, nvar
          call r8vec_normal_01(ntarget, draw)
          simulation_sample(isim,ivar,:) = draw
        end do
      end do
      call solver%set_sim(randpath=simulation_path, sample=simulation_sample)
      if (kriging_failed()) goto 900
    end if

    ! Search anisotropy is configured per variable from its first auto-variogram structure.
    do ivar = 1, nvar
      vidx = vgm_index(ivar,ivar)
      call set_variable_search(solver, ivar, variables(ivar), &
                               request%search%anisotropic_search, variograms(1,vidx))
      if (kriging_failed()) goto 900
    end do

    call solver%solve(nthread=1)
    if (kriging_failed() .or. .not. solver%solved) goto 900

    ! Texture2Par's current output contract stores one primary realization.
    ! In estimate mode the first value plane is the deterministic estimate.
    result%estimates = solver%block%value(1,1,:)
    result%variance  = solver%block%variance(1,1,:)
    result%status = 0
    succeeded = .true.

900 continue
    if (.not. succeeded) failure = kriging_last_error()
    if (initialized) call solver%finalize()
    if (.not. succeeded) then
      if (len_trim(failure) == 0) failure = 'libkriging solve did not complete'
      call error_handler(4, opt_msg='libkriging: '//trim(failure))
    end if
  end subroutine libkriging_solve_variables
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
      call error_handler(4, opt_msg='libkriging many-value observation geometry/value shape mismatch')
    end if

    allocate(estimates(ntarget, nvalues), variance(ntarget, nvalues))
    estimates = 0.0
    variance = 0.0
    if (ntarget == 0 .or. nvalues == 0) return

    if (request%dimension < 2 .or. request%dimension > 3) then
      call error_handler(4, opt_msg='libkriging supports only 2D and 3D interpolation')
    end if
    if (nstruct < 1 .or. nstruct > size(variograms)) then
      call error_handler(4, opt_msg='libkriging received an invalid variogram structure count')
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

    call set_solver_grid(solver, request)
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
      call error_handler(4, opt_msg='libkriging: '//trim(failure))
    end if
  end subroutine libkriging_interpolate_many
!-------------------------------------------------------------------------------------------------!

  subroutine set_solver_grid(solver, request)
    type(t_kriging), intent(inout)       :: solver
    type(t_interp_request), intent(in)   :: request

    if (request%request_block_support) then
      if (.not. allocated(request%target_block_size)) &
        call error_handler(4, opt_msg='libkriging Gaussian block dimensions are incomplete')
      if (size(request%target_block_size,1) /= request%dimension .or. &
          size(request%target_block_size,2) /= size(request%target_coordinates,2)) &
        call error_handler(4, opt_msg='libkriging Gaussian block dimension shape mismatch')
      call solver%set_grid(request%target_coordinates, block_type=-4, &
                           blocksize=request%target_block_size)
    else
      call solver%set_grid(request%target_coordinates)
    end if
  end subroutine set_solver_grid
!-------------------------------------------------------------------------------------------------!

  subroutine set_variable_search(solver, ivar, variable, anisotropic_search, variogram)
    type(t_kriging), intent(inout)          :: solver
    integer, intent(in)                     :: ivar
    type(t_libkriging_variable), intent(in) :: variable
    logical, intent(in)                     :: anisotropic_search
    type(t_vgm), intent(in)                 :: variogram

    if (anisotropic_search) then
      if (valid_maxdist(variable%max_distance)) then
        call solver%set_search(ivar, anis1=variogram%anis1, anis2=variogram%anis2, &
                               azimuth=variogram%ang1, dip=variogram%ang2, &
                               plunge=variogram%ang3, nmax=variable%nmax, &
                               maxdist=variable%max_distance, sector_search=.false.)
      else
        call solver%set_search(ivar, anis1=variogram%anis1, anis2=variogram%anis2, &
                               azimuth=variogram%ang1, dip=variogram%ang2, &
                               plunge=variogram%ang3, nmax=variable%nmax, &
                               sector_search=.false.)
      end if
    else
      if (valid_maxdist(variable%max_distance)) then
        call solver%set_search(ivar, nmax=variable%nmax, maxdist=variable%max_distance, &
                               sector_search=.false.)
      else
        call solver%set_search(ivar, nmax=variable%nmax, sector_search=.false.)
      end if
    end if
  end subroutine set_variable_search
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
        call error_handler(4, opt_msg='libkriging received an unsupported variogram type')
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
