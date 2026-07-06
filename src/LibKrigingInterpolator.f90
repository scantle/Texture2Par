! Batch-oriented Texture2Par backend for native libkriging OK and SK.
module m_libkriging_interpolator
  use m_kriging_base, only: t_kriging_base
  use m_datasets, only: t_layerdata
  use m_grid, only: t_grid
  use m_options, only: t_options, CK_MEAN_LIBKRIGING, CK_MEAN_LEGACY_STANDARDIZED
  use m_interpolation_types, only: t_interp_request, t_interp_result, INTERP_PURPOSE_MODEL, &
                                   PREDICT_ESTIMATE, PREDICT_SGSIM
  use m_interpolation_packers, only: pack_grid_targets, pack_layerdata_observations, &
                                     validate_unique_observations
  use m_libkriging_adapter, only: libkriging_solve, libkriging_interpolate_many, &
                                  t_libkriging_variable
  use m_error_handler, only: error_handler
  use m_global, only: log
  implicit none
  private

  integer, parameter, public :: LIBKRIGE_MODE_ORDINARY = 1
  integer, parameter, public :: LIBKRIGE_MODE_SIMPLE   = 2

  integer, parameter, public :: LIBKRIGE_PREDICT_ESTIMATE = PREDICT_ESTIMATE
  integer, parameter, public :: LIBKRIGE_PREDICT_SGSIM    = PREDICT_SGSIM

  type, extends(t_kriging_base), public :: t_libkrige
    ! OK and SK are configurations of this one backend, matching libkriging's t_kriging API.
    integer :: kriging_mode = LIBKRIGE_MODE_ORDINARY
    integer :: prediction_mode = LIBKRIGE_PREDICT_ESTIMATE
    integer :: simulation_sequence = 0
    integer :: nthreads = 1
  contains
    procedure :: calc => calc_libkrige
    procedure :: interpolate_many => interpolate_many_libkrige
  end type t_libkrige

contains

!-------------------------------------------------------------------------------------------------!
  subroutine calc_libkrige(this, obs_data, obs_idx, val_idx, to_grid, active, results, variance, &
                           catlist, to_cats, opt, layer)
    class(t_libkrige), intent(inout) :: this
    class(t_layerdata)               :: obs_data(:)
    integer, intent(in)              :: obs_idx(:), val_idx(:)
    type(t_grid), intent(in)         :: to_grid
    real                            :: results(:), variance(:)
    logical, intent(in)             :: active(:)
    integer, intent(in)             :: catlist(:), to_cats(:)
    type(t_options), intent(in)      :: opt
    integer, intent(in)              :: layer

    type(t_interp_request) :: request, variable_request
    type(t_interp_result)  :: batch
    type(t_libkriging_variable), allocatable :: variables(:)
    integer, allocatable   :: vgm_index(:,:)
    integer                :: i, d, cat, zone, vgm_idx, nmax, unbias
    logical                :: missing_variable, shift_secondary_mean, local_mean_correction

    call validate_libkriging_scope(obs_data, obs_idx, val_idx, opt)
    this%nthreads = opt%NTHREADS
    do d = 1, size(obs_idx)
      call validate_metadata(this, obs_idx(d))
    end do

    results = 0.0
    variance = 0.0

    select case (this%kriging_mode)
      case (LIBKRIGE_MODE_ORDINARY)
        unbias = 1
      case (LIBKRIGE_MODE_SIMPLE)
        unbias = 0
      case default
        call error_handler(4, opt_msg='libkriging has an invalid kriging mode')
    end select
    if (this%kriging_mode == LIBKRIGE_MODE_ORDINARY) then
      call configure_ck_mean_mode(opt, shift_secondary_mean, local_mean_correction)
    else
      shift_secondary_mean = .false.
      local_mean_correction = .false.
    end if

    allocate(variables(size(obs_data)), vgm_index(size(obs_data),size(obs_data)))

    do i = 1, size(catlist)
      cat = catlist(i)
      missing_variable = .false.
      do d = 1, size(obs_data)
        call pack_layerdata_observations(variable_request, obs_data(d:d), obs_idx(d:d), &
                                         val_idx(d:d), opt%INTERP_DIM, cat)
        if (size(variable_request%observation_values) == 0) then
          missing_variable = .true.
          cycle
        end if
        call validate_unique_observations(variable_request)
        variables(d)%coordinates = variable_request%observation_coordinates
        variables(d)%values = variable_request%observation_values
        if (allocated(variable_request%observation_variance)) then
          variables(d)%variance = variable_request%observation_variance
        else
          allocate(variables(d)%variance(size(variables(d)%values)))
          variables(d)%variance = 0.0
        end if
        if (d > 1 .and. shift_secondary_mean) variables(d)%values = variables(d)%values + &
          obs_data(1)%mean(val_idx(1)) - obs_data(d)%mean(val_idx(d))
        if (this%kriging_mode == LIBKRIGE_MODE_SIMPLE) &
          variables(d)%sk_mean = obs_data(d)%mean(val_idx(d))
        variables(d)%nmax = min(this%nnear(obs_idx(d)), size(variables(d)%values))
      end do
      if (missing_variable) then
        write(*,'(2(A,I0))') 'Warning - no data for layer ', layer, ' in category ', cat
        write(log%unit,'(2(A,I0))') 'Warning - no data for layer ', layer, ' in category ', cat
        cycle
      end if

      call pack_grid_targets(request, to_grid, active, opt%INTERP_DIM, to_cats, cat, &
                             include_block_support=opt%BLOCK_KRIGE > 0)
      if (size(request%target_index) == 0) cycle

      zone = 0
      if (allocated(this%cat2zone)) then
        if (cat < 1 .or. cat > size(this%cat2zone)) &
          call error_handler(4, opt_msg='libkriging OK category is outside the category-to-zone map')
        zone = this%cat2zone(cat)
      end if
      vgm_index = this%global_vario_idx(zone, obs_idx, obs_idx)
      if (any(vgm_index < 1) .or. any(vgm_index > this%nvario)) &
        call error_handler(4, opt_msg='libkriging OK could not map all auto/cross variograms')

      variables(1)%max_distance = opt%MAX_PRIMARY_DIST
      do d = 2, size(variables)
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
      request%purpose = INTERP_PURPOSE_MODEL
      request%model_layer = layer
      request%class_id = obs_idx(1)
      request%zone_id = zone
      request%store_weight = opt%WRITE_GRID_WEIGHTS
      request%use_old_weight = opt%READ_GRID_WEIGHTS
      if (request%store_weight .or. request%use_old_weight) &
        write(request%weight_file, '(A,I2.2,A,I0,A,I0,A)') &
          trim(opt%PROJECT_NAME)//'_grdweights_arr', layer, '_obs', obs_idx(1), &
          '_cat', cat, '.wts'
      request%search%n_nearest_primary = variables(1)%nmax
      request%search%max_primary_distance = opt%MAX_PRIMARY_DIST
      if (size(variables) > 1) then
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

      call libkriging_solve(request, variables, this%variograms, this%nstruct, &
                            vgm_index, unbias, local_mean_correction, opt%CORRECT_WEIGHTS, &
                            opt%WRITE_MATRIX, this%nthreads, batch)
      results(batch%target_index) = batch%estimates
      variance(batch%target_index) = batch%variance
    end do

    ! TODO: when nnear >= nobs, use a cached/global solve if libkriging exposes
    ! that optimization while preserving kriging_mode (never silently convert OK to SK).
  end subroutine calc_libkrige
!-------------------------------------------------------------------------------------------------!

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

  ! Batch values share one packed geometry. value_matrix is (observation,
  ! value-column); estimates and variance are returned as (target,value-column).
  subroutine interpolate_many_libkrige(this, request, value_matrix, opt, estimates, variance)
    class(t_libkrige), intent(inout)  :: this
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
    call validate_unique_observations(request)
    if (.not. allocated(request%observation_variance)) then
      allocate(request%observation_variance(size(value_matrix,1)))
      request%observation_variance = 0.0
    end if

    select case (this%kriging_mode)
      case (LIBKRIGE_MODE_ORDINARY)
        unbias = 1
        sk_mean = 0.0
      case (LIBKRIGE_MODE_SIMPLE)
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

  subroutine validate_libkriging_scope(obs_data, obs_idx, val_idx, opt)
    class(t_layerdata)          :: obs_data(:)
    integer, intent(in)         :: obs_idx(:), val_idx(:)
    type(t_options), intent(in) :: opt

    if (size(obs_data) /= size(obs_idx) .or. size(obs_data) /= size(val_idx)) &
      call error_handler(4, opt_msg='libkriging received inconsistent cokriging variable metadata')
    if (size(obs_data) > 1 .and. opt%ASYMMETRIC_COV) &
      call error_handler(4, opt_msg='libkriging CK requires symmetric cross-variograms; ASYMMETRIC_COV is unsupported')
    if (opt%NSIM > 1) &
      call error_handler(1, opt_msg='Texture2Par SGSIM output currently supports exactly one realization')
  end subroutine validate_libkriging_scope
!-------------------------------------------------------------------------------------------------!

  subroutine validate_metadata(this, class_id)
    class(t_libkrige), intent(in) :: this
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
