! Batch-oriented Texture2Par backend for native libkriging OK and SK.
module m_libkriging_interpolator
  use m_kriging_base, only: t_kriging_base
  use m_options, only: t_options, CK_MEAN_LIBKRIGING, CK_MEAN_LEGACY_STANDARDIZED
  use m_interpolation_types, only: t_interp_request, t_interp_result, INTERP_PURPOSE_MODEL, &
                                   PREDICT_ESTIMATE, PREDICT_SGSIM
  use m_libkriging_adapter, only: libkriging_solve, libkriging_interpolate_many, &
                                  t_libkriging_variable
  use m_error_handler, only: error_handler
  implicit none
  private

  integer, parameter, public :: KRIGING_MODE_ORDINARY = 1
  integer, parameter, public :: KRIGING_MODE_SIMPLE   = 2

  integer, parameter, public :: KRIGING_PREDICT_ESTIMATE = PREDICT_ESTIMATE
  integer, parameter, public :: KRIGING_PREDICT_SGSIM    = PREDICT_SGSIM

  type, extends(t_kriging_base), public :: t_kriging_interpolator
    ! OK and SK are configurations of this one backend, matching libkriging's t_kriging API.
    integer :: kriging_mode = KRIGING_MODE_ORDINARY
    integer :: prediction_mode = KRIGING_PREDICT_ESTIMATE
    integer :: simulation_sequence = 0
    integer :: nthreads = 1
  contains
    procedure :: interpolate => interpolate_kriging
    procedure, private :: interpolate_many => interpolate_many_libkrige
  end type t_kriging_interpolator

contains

  subroutine configure_ck_mean_mode(opt, shift_secondary_mean, local_mean_correction)
    type(t_options), intent(in) :: opt
    logical, intent(out)        :: shift_secondary_mean, local_mean_correction

    select case (opt%CK_MEAN_MODE)
      case (CK_MEAN_LIBKRIGING)
        shift_secondary_mean = .false.
        local_mean_correction = .true.
      case (CK_MEAN_LEGACY_STANDARDIZED)
        shift_secondary_mean = .true.
        local_mean_correction = .false.
      case DEFAULT
        call error_handler(4, opt_msg='Invalid CK mean mode in libkriging interpolator')
    end select
  end subroutine configure_ck_mean_mode
!-------------------------------------------------------------------------------------------------!

  subroutine interpolate_kriging(this, request, values, opt, result)
    class(t_kriging_interpolator), intent(inout) :: this
    type(t_interp_request), intent(inout)         :: request
    real, intent(in)                              :: values(:,:)
    type(t_options), intent(in)                   :: opt
    type(t_interp_result), intent(inout)          :: result
    real, allocatable                             :: estimates(:,:), variance(:,:)

    if (request%purpose == INTERP_PURPOSE_MODEL) then
      call interpolate_model_kriging(this, request, values, opt, result)
      return
    end if

    call this%interpolate_many(request, values, opt, estimates, variance)
    if (allocated(result%estimates_many)) deallocate(result%estimates_many)
    if (allocated(result%variance_many)) deallocate(result%variance_many)
    if (allocated(result%target_index)) deallocate(result%target_index)
    call move_alloc(estimates, result%estimates_many)
    call move_alloc(variance, result%variance_many)
    allocate(result%target_index(size(request%target_index)))
    result%target_index = request%target_index
    result%status = 0
  end subroutine interpolate_kriging
!-------------------------------------------------------------------------------------------------!

  subroutine interpolate_model_kriging(this, request, values, opt, result)
    class(t_kriging_interpolator), intent(inout) :: this
    type(t_interp_request), intent(inout)         :: request
    real, intent(in)                              :: values(:,:)
    type(t_options), intent(in)                   :: opt
    type(t_interp_result), intent(inout)          :: result
    type(t_libkriging_variable), allocatable      :: variables(:)
    type(t_interp_result)                         :: solve_result
    integer, allocatable                          :: indices(:), vgm_index(:,:)
    logical, allocatable                          :: variable_mask(:)
    integer                                       :: d, i, nvar, unbias, zone, vgm_idx
    logical                                       :: shift_secondary_mean, local_mean_correction

    if (.not. allocated(request%variable_class_id) .or. &
        .not. allocated(request%variable_mean) .or. &
        .not. allocated(request%observation_variable_id)) &
      call error_handler(4, opt_msg='libkriging model interpolation requires packed variable metadata')
    if (size(values,1) /= size(request%observation_coordinates,2) .or. size(values,2) /= 1) &
      call error_handler(4, opt_msg='libkriging model observation geometry/value shape mismatch')

    nvar = size(request%variable_class_id)
    if (size(request%variable_mean) /= nvar) &
      call error_handler(4, opt_msg='libkriging model variable metadata shape mismatch')
    if (nvar > 1 .and. opt%ASYMMETRIC_COV) &
      call error_handler(4, opt_msg='libkriging CK requires symmetric cross-variograms; ASYMMETRIC_COV is unsupported')
    if (opt%NSIM > 1) &
      call error_handler(1, opt_msg='Texture2Par SGSIM output currently supports exactly one realization')

    this%nthreads = opt%NTHREADS
    do d = 1, nvar
      call validate_metadata(this, request%variable_class_id(d))
    end do
    select case (this%kriging_mode)
      case (KRIGING_MODE_ORDINARY)
        unbias = 1
        call configure_ck_mean_mode(opt, shift_secondary_mean, local_mean_correction)
      case (KRIGING_MODE_SIMPLE)
        unbias = 0
        shift_secondary_mean = .false.
        local_mean_correction = .false.
      case default
        call error_handler(4, opt_msg='libkriging has an invalid kriging mode')
    end select

    allocate(variables(nvar), vgm_index(nvar,nvar))
    allocate(variable_mask(size(values,1)))
    do d = 1, nvar
      variable_mask = request%observation_variable_id == d
      indices = pack([(i,i=1,size(values,1))], variable_mask)
      if (size(indices) == 0) &
        call error_handler(4, opt_msg='libkriging model request contains an empty variable')
      variables(d)%coordinates = request%observation_coordinates(:,indices)
      variables(d)%values = values(indices,1)
      allocate(variables(d)%variance(size(indices)))
      if (allocated(request%observation_variance)) then
        variables(d)%variance = request%observation_variance(indices)
      else
        variables(d)%variance = 0.0
      end if
      if (d > 1 .and. shift_secondary_mean) variables(d)%values = variables(d)%values + &
        request%variable_mean(1) - request%variable_mean(d)
      if (this%kriging_mode == KRIGING_MODE_SIMPLE) &
        variables(d)%sk_mean = request%variable_mean(d)
      variables(d)%nmax = min(this%nnear(request%variable_class_id(d)), size(indices))
    end do

    zone = 0
    if (allocated(this%cat2zone)) then
      if (request%category_id < 1 .or. request%category_id > size(this%cat2zone)) &
        call error_handler(4, opt_msg='libkriging OK category is outside the category-to-zone map')
      zone = this%cat2zone(request%category_id)
    end if
    vgm_index = this%global_vario_idx(zone, request%variable_class_id, request%variable_class_id)
    if (any(vgm_index < 1) .or. any(vgm_index > this%nvario)) &
      call error_handler(4, opt_msg='libkriging OK could not map all auto/cross variograms')

    variables(1)%max_distance = opt%MAX_PRIMARY_DIST
    do d = 2, nvar
      if (opt%MAX_SECONDARY_DIST > 0.1) then
        variables(d)%max_distance = opt%MAX_SECONDARY_DIST
      else
        vgm_idx = vgm_index(d,d)
        variables(d)%max_distance = 3.0 * &
          maxval(this%variograms(1:this%nstruct(vgm_idx),vgm_idx)%a_hmax)
      end if
    end do
    if (any(variables%nmax < 1)) &
      call error_handler(4, opt_msg='libkriging OK requires NNEAR to be at least one for every variable')

    request%zone_id = zone
    request%search%n_nearest_primary = variables(1)%nmax
    request%search%max_primary_distance = variables(1)%max_distance
    if (nvar > 1) then
      request%search%n_nearest_secondary = variables(2)%nmax
      request%search%max_secondary_distance = variables(2)%max_distance
    end if
    request%search%anisotropic_search = opt%ANISOTROPIC_SEARCH
    request%prediction_mode = this%prediction_mode
    if (request%prediction_mode == PREDICT_SGSIM) then
      this%simulation_sequence = this%simulation_sequence + 1
      request%simulation_count = opt%NSIM
      request%random_seed = opt%SEED
      if (request%random_seed > 0) &
        request%random_seed = request%random_seed + this%simulation_sequence - 1
    end if

    call libkriging_solve(request, variables, this%variograms, this%nstruct, vgm_index, &
                          unbias, local_mean_correction, opt%CORRECT_WEIGHTS, &
                          opt%WRITE_MATRIX, this%nthreads, solve_result)
    if (allocated(result%estimates_many)) deallocate(result%estimates_many)
    if (allocated(result%variance_many)) deallocate(result%variance_many)
    if (allocated(result%target_index)) deallocate(result%target_index)
    allocate(result%estimates_many(size(solve_result%estimates),1))
    allocate(result%variance_many(size(solve_result%variance),1))
    allocate(result%target_index(size(solve_result%target_index)))
    result%estimates_many(:,1) = solve_result%estimates
    result%variance_many(:,1) = solve_result%variance
    result%target_index = solve_result%target_index
    result%status = solve_result%status
  end subroutine interpolate_model_kriging
!-------------------------------------------------------------------------------------------------!

  ! Batch values share one packed geometry. value_matrix is (observation,
  ! value-column); estimates and variance are returned as (target,value-column).
  subroutine interpolate_many_libkrige(this, request, value_matrix, opt, estimates, variance)
    class(t_kriging_interpolator), intent(inout)  :: this
    type(t_interp_request), intent(inout) :: request
    real, intent(in)                  :: value_matrix(:,:)
    type(t_options), intent(in)       :: opt
    real, allocatable, intent(out)    :: estimates(:,:), variance(:,:)

    integer :: zone, vgm_idx, nmax, unbias
    real    :: sk_mean

    call validate_metadata(this, request%class_id)
    this%nthreads = opt%NTHREADS
    if (size(value_matrix,1) /= size(request%observation_coordinates,2)) &
      call error_handler(4, opt_msg='libkriging OK many-value observation geometry/value shape mismatch')
    if (.not. allocated(request%observation_variance)) then
      allocate(request%observation_variance(size(value_matrix,1)))
      request%observation_variance = 0.0
    end if

    select case (this%kriging_mode)
      case (KRIGING_MODE_ORDINARY)
        unbias = 1
        sk_mean = 0.0
      case (KRIGING_MODE_SIMPLE)
        call error_handler(4, opt_msg='libkriging SK pilot-point/interbed interpolation is not enabled; '// &
                                     'those workflows do not define an SK mean convention')
      case default
        call error_handler(4, opt_msg='libkriging has an invalid kriging mode')
    end select

    zone = 0
    if (allocated(this%cat2zone)) then
      if (request%category_id < 1 .or. request%category_id > size(this%cat2zone)) &
        call error_handler(4, opt_msg='libkriging OK category is outside the category-to-zone map')
      zone = this%cat2zone(request%category_id)
    end if
    if (zone < lbound(this%global_vario_idx,1) .or. &
        zone > ubound(this%global_vario_idx,1)) &
      call error_handler(4, opt_msg='libkriging OK category has no variogram zone mapping')
    vgm_idx = this%global_vario_idx(zone, request%class_id, request%class_id)
    if (vgm_idx < 1 .or. vgm_idx > this%nvario) &
      call error_handler(4, opt_msg='libkriging OK could not map the category and class to a variogram')

    nmax = min(this%nnear(request%class_id), size(value_matrix,1))
    if (nmax < 1) call error_handler(4, opt_msg='libkriging OK requires NNEAR to be at least one')
    request%zone_id = zone
    request%search%n_nearest_primary = nmax
    request%search%max_primary_distance = opt%MAX_PRIMARY_DIST
    request%search%anisotropic_search = opt%ANISOTROPIC_SEARCH

    call libkriging_interpolate_many(request, this%variograms(:,vgm_idx), &
                                      this%nstruct(vgm_idx), unbias, sk_mean, &
                                      opt%CORRECT_WEIGHTS, opt%WRITE_MATRIX, this%nthreads, &
                                      value_matrix, estimates, variance)
  end subroutine interpolate_many_libkrige
!-------------------------------------------------------------------------------------------------!

  subroutine validate_metadata(this, class_id)
    class(t_kriging_interpolator), intent(in) :: this
    integer, intent(in)            :: class_id

    if (.not. associated(this%variograms) .or. .not. associated(this%nstruct) .or. &
        .not. associated(this%nnear)) then
      call error_handler(4, opt_msg='libkriging OK variogram/search metadata was not initialized')
    end if
    if (.not. allocated(this%global_vario_idx)) &
      call error_handler(4, opt_msg='libkriging OK variogram mapping was not initialized')
    if (class_id < 1 .or. class_id > size(this%global_vario_idx,2)) &
      call error_handler(4, opt_msg='libkriging OK class is outside the variogram map')
  end subroutine validate_metadata
!-------------------------------------------------------------------------------------------------!

end module m_libkriging_interpolator
