! Experimental batch-oriented Texture2Par backend for native libkriging.
module m_libkriging_interpolator
  use m_kriging_base, only: t_kriging_base
  use m_datasets, only: t_layerdata
  use m_grid, only: t_grid
  use m_options, only: t_options
  use m_interpolation_types, only: t_interp_request, t_interp_result, INTERP_PURPOSE_MODEL
  use m_interpolation_packers, only: pack_grid_targets, pack_layerdata_observations, &
                                     validate_unique_observations
  use m_libkriging_adapter, only: libkriging_solve, libkriging_interpolate_many
  use m_error_handler, only: error_handler
  use m_global, only: log
  implicit none
  private

  integer, parameter, public :: LIBKRIGE_MODE_ORDINARY = 1
  integer, parameter, public :: LIBKRIGE_MODE_SIMPLE   = 2

  type, extends(t_kriging_base), public :: t_libkrige
    ! OK and SK are configurations of this one backend, matching libkriging's t_kriging API.
    integer :: kriging_mode = LIBKRIGE_MODE_ORDINARY
    real    :: sk_mean = 0.0
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

    type(t_interp_request) :: request
    type(t_interp_result)  :: batch
    integer                :: i, cat, zone, vgm_idx, nmax, unbias

    call validate_ok_lib_scope(obs_data, obs_idx, val_idx, opt)
    call validate_metadata(this, obs_idx(1))

    results = 0.0
    variance = 0.0

    select case (this%kriging_mode)
      case (LIBKRIGE_MODE_ORDINARY)
        unbias = 1
      case (LIBKRIGE_MODE_SIMPLE)
        unbias = 0
      case default
        call error_handler(4, opt_msg='OK_LIB has an invalid kriging mode')
    end select

    do i = 1, size(catlist)
      cat = catlist(i)
      call pack_layerdata_observations(request, obs_data, obs_idx, val_idx, &
                                       opt%INTERP_DIM, cat)
      if (size(request%observation_values) == 0) then
        write(*,'(2(A,I0))') 'Warning - no data for layer ', layer, ' in category ', cat
        write(log%unit,'(2(A,I0))') 'Warning - no data for layer ', layer, ' in category ', cat
        cycle
      end if
      call validate_unique_observations(request)
      if (.not. allocated(request%observation_variance)) then
        allocate(request%observation_variance(size(request%observation_values)))
        request%observation_variance = 0.0
      end if

      call pack_grid_targets(request, to_grid, active, opt%INTERP_DIM, to_cats, cat)
      if (size(request%target_index) == 0) cycle

      zone = 0
      if (allocated(this%cat2zone)) then
        if (cat < 1 .or. cat > size(this%cat2zone)) &
          call error_handler(4, opt_msg='OK_LIB category is outside the category-to-zone map')
        zone = this%cat2zone(cat)
      end if
      vgm_idx = this%global_vario_idx(zone, obs_idx(1), obs_idx(1))
      if (vgm_idx < 1 .or. vgm_idx > this%nvario) &
        call error_handler(4, opt_msg='OK_LIB could not map the category and class to a variogram')

      nmax = min(this%nnear(obs_idx(1)), size(request%observation_values))
      if (nmax < 1) call error_handler(4, opt_msg='OK_LIB requires NNEAR to be at least one')
      request%purpose = INTERP_PURPOSE_MODEL
      request%model_layer = layer
      request%class_id = obs_idx(1)
      request%zone_id = zone
      request%search%n_nearest_primary = nmax
      request%search%max_primary_distance = opt%MAX_PRIMARY_DIST
      request%search%anisotropic_search = opt%ANISOTROPIC_SEARCH

      call libkriging_solve(request, this%variograms(:,vgm_idx), this%nstruct(vgm_idx), &
                            unbias, this%sk_mean, batch)
      results(batch%target_index) = batch%estimates
      variance(batch%target_index) = batch%variance
    end do

    ! TODO(Stage 6C+): when nnear >= nobs, use a cached/global solve if libkriging exposes
    ! that optimization while preserving kriging_mode (never silently convert OK to SK).
  end subroutine calc_libkrige
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

    call validate_metadata(this, request%class_id)
    if (size(value_matrix,1) /= size(request%observation_coordinates,2)) &
      call error_handler(4, opt_msg='OK_LIB many-value observation geometry/value shape mismatch')
    call validate_unique_observations(request)
    if (.not. allocated(request%observation_variance)) then
      allocate(request%observation_variance(size(value_matrix,1)))
      request%observation_variance = 0.0
    end if

    select case (this%kriging_mode)
      case (LIBKRIGE_MODE_ORDINARY); unbias = 1
      case (LIBKRIGE_MODE_SIMPLE);   unbias = 0
      case default
        call error_handler(4, opt_msg='OK_LIB has an invalid kriging mode')
    end select

    zone = 0
    if (allocated(this%cat2zone)) then
      if (request%category_id < 1 .or. request%category_id > size(this%cat2zone)) &
        call error_handler(4, opt_msg='OK_LIB category is outside the category-to-zone map')
      zone = this%cat2zone(request%category_id)
    end if
    if (zone < lbound(this%global_vario_idx,1) .or. &
        zone > ubound(this%global_vario_idx,1)) &
      call error_handler(4, opt_msg='OK_LIB category has no variogram zone mapping')
    vgm_idx = this%global_vario_idx(zone, request%class_id, request%class_id)
    if (vgm_idx < 1 .or. vgm_idx > this%nvario) &
      call error_handler(4, opt_msg='OK_LIB could not map the category and class to a variogram')

    nmax = min(this%nnear(request%class_id), size(value_matrix,1))
    if (nmax < 1) call error_handler(4, opt_msg='OK_LIB requires NNEAR to be at least one')
    request%zone_id = zone
    request%search%n_nearest_primary = nmax
    request%search%max_primary_distance = opt%MAX_PRIMARY_DIST
    request%search%anisotropic_search = opt%ANISOTROPIC_SEARCH

    call libkriging_interpolate_many(request, this%variograms(:,vgm_idx), &
                                      this%nstruct(vgm_idx), unbias, this%sk_mean, &
                                      value_matrix, estimates, variance)
  end subroutine interpolate_many_libkrige
!-------------------------------------------------------------------------------------------------!

  subroutine validate_ok_lib_scope(obs_data, obs_idx, val_idx, opt)
    class(t_layerdata)          :: obs_data(:)
    integer, intent(in)         :: obs_idx(:), val_idx(:)
    type(t_options), intent(in) :: opt

    if (size(obs_data) /= 1 .or. size(obs_idx) /= 1 .or. size(val_idx) /= 1) then
      call error_handler(4, opt_msg='OK_LIB currently supports one primary variable; cokriging is not supported')
    end if
    if (opt%NSIM > 0) &
      call error_handler(1, opt_msg='OK_LIB does not currently support SGSIM')
    if (opt%BLOCK_KRIGE > 0) &
      call error_handler(1, opt_msg='OK_LIB does not currently support block kriging')
    if (opt%READ_GRID_WEIGHTS .or. opt%WRITE_GRID_WEIGHTS) &
      call error_handler(1, opt_msg='OK_LIB does not support READ_GRID_WEIGHTS or WRITE_GRID_WEIGHTS currently')
  end subroutine validate_ok_lib_scope
!-------------------------------------------------------------------------------------------------!

  subroutine validate_metadata(this, class_id)
    class(t_libkrige), intent(in) :: this
    integer, intent(in)            :: class_id

    if (.not. associated(this%variograms) .or. .not. associated(this%nstruct) .or. &
        .not. associated(this%nnear)) then
      call error_handler(4, opt_msg='OK_LIB variogram/search metadata was not initialized')
    end if
    if (.not. allocated(this%global_vario_idx)) &
      call error_handler(4, opt_msg='OK_LIB variogram mapping was not initialized')
    if (class_id < 1 .or. class_id > size(this%global_vario_idx,2)) &
      call error_handler(4, opt_msg='OK_LIB class is outside the variogram map')
  end subroutine validate_metadata
!-------------------------------------------------------------------------------------------------!

end module m_libkriging_interpolator
