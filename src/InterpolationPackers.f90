! Helpers for translating existing Texture2Par storage into batch requests.
!
! These routines are the one-way adapter from Texture2Par grid/dataset types
! to the backend-neutral batch interpolation request.
module m_interpolation_packers
  use m_interpolation_types, only: t_interp_request
  use m_datasets, only: t_layerdata
  use m_grid, only: t_grid
  use m_pilotpoints, only: t_pilotpoints
  use m_categories, only: CAT_ALL
  use m_global, only: NODATA
  use m_error_handler, only: error_handler
  implicit none
  private

  public :: pack_grid_targets
  public :: pack_layerdata_observations
  public :: pack_pilotpoint_observations

contains

!-------------------------------------------------------------------------------------------------!

  subroutine pack_grid_targets(request, grid, active, dimension, target_categories, category_id, &
                               include_block_support)
    type(t_interp_request), intent(inout) :: request
    type(t_grid), intent(in)              :: grid
    logical, intent(in)                   :: active(:)
    integer, intent(in)                   :: dimension
    integer, intent(in), optional         :: target_categories(:)
    integer, intent(in), optional         :: category_id
    logical, intent(in), optional         :: include_block_support
    logical, allocatable                  :: pack_mask(:)
    logical                               :: pack_blocks
    integer                               :: i, j, ntarget

    if (dimension < 2 .or. dimension > 3 .or. dimension > grid%dim) then
      call packer_error('target dimension is incompatible with the grid')
    end if
    if (size(active) /= grid%n) then
      call packer_error('target active mask size does not match the grid')
    end if
    if (present(target_categories)) then
      if (.not. present(category_id)) then
        call packer_error('target categories require a category id')
      end if
      if (size(target_categories) /= grid%n) then
        call packer_error('target category array size does not match the grid')
      end if
    end if

    allocate(pack_mask(grid%n))
    pack_mask = active
    if (present(target_categories)) then
      ! CAT_ALL is an observation wildcard only; target categories match exactly.
      pack_mask = pack_mask .and. (target_categories == category_id)
    end if

    ntarget = count(pack_mask)
    call clear_targets(request)
    allocate(request%target_coordinates(dimension, ntarget), request%target_index(ntarget))

    j = 0
    do i = 1, grid%n
      if (.not. pack_mask(i)) cycle
      j = j + 1
      request%target_coordinates(:,j) = grid%coords(1:dimension,i)
      request%target_index(j) = i
    end do

    pack_blocks = .false.
    if (present(include_block_support)) pack_blocks = include_block_support
    request%request_block_support = pack_blocks
    if (request%request_block_support) then
      if (.not. associated(grid%active_blocks)) &
        call packer_error('block dimensions were requested but are not attached to the target grid')
      if (.not. allocated(grid%active_blocks%block_size)) &
        call packer_error('block target dimensions are incomplete')
      if (size(grid%active_blocks%block_size,1) < dimension .or. &
          size(grid%active_blocks%block_size,2) /= grid%n) &
        call packer_error('block target dimensions do not match the target grid')

      allocate(request%target_block_size(dimension,ntarget))
      do j = 1, ntarget
        request%target_block_size(:,j) = &
          grid%active_blocks%block_size(1:dimension,request%target_index(j))
      end do
      if (any(request%target_block_size <= 0.0)) &
        call packer_error('Gaussian block dimensions must be positive')
    end if

    request%dimension = dimension
    request%target_total_count = grid%n
    if (present(category_id)) request%category_id = category_id
  end subroutine pack_grid_targets

!-------------------------------------------------------------------------------------------------!

  subroutine pack_layerdata_observations( &
    request, observations, class_idx, dataval_idx, dimension, category_id)
    type(t_interp_request), intent(inout) :: request
    type(t_layerdata), intent(in)         :: observations(:)
    integer, intent(in)                   :: class_idx(:)
    integer, intent(in)                   :: dataval_idx(:)
    integer, intent(in)                   :: dimension, category_id
    logical, allocatable                  :: pack_mask(:)
    integer                               :: d, i, j, nobs, total_count

    if (dimension < 2 .or. dimension > 3) then
      call packer_error('observation dimension must be 2 or 3')
    end if
    if (size(class_idx) /= size(observations) .or. &
        size(dataval_idx) /= size(observations)) then
      call packer_error('class and data-value index arrays must match the datasets')
    end if

    nobs = 0
    total_count = 0
    do d = 1, size(observations)
      call check_layerdata(observations(d), dataval_idx(d), dimension)
      allocate(pack_mask(observations(d)%grid%n))
      call fill_category_value_mask(observations(d), dataval_idx(d), category_id, pack_mask)
      nobs = nobs + count(pack_mask)
      total_count = total_count + observations(d)%grid%n
      deallocate(pack_mask)
    end do

    call clear_observations(request)
    allocate(request%observation_coordinates(dimension, nobs))
    allocate(request%observation_values(nobs))
    allocate(request%observation_dataset_id(nobs))
    allocate(request%observation_class_id(nobs))
    allocate(request%observation_variable_id(nobs))
    allocate(request%observation_source_index(nobs))
    allocate(request%observation_variance(nobs))
    allocate(request%variable_class_id(size(observations)))
    allocate(request%variable_mean(size(observations)))
    request%observation_variance = 0.0
    j = 0
    do d = 1, size(observations)
      request%variable_class_id(d) = class_idx(d)
      request%variable_mean(d) = observations(d)%mean(dataval_idx(d))
      allocate(pack_mask(observations(d)%grid%n))
      call fill_category_value_mask(observations(d), dataval_idx(d), category_id, pack_mask)
      do i = 1, observations(d)%grid%n
        if (.not. pack_mask(i)) cycle
        j = j + 1
        request%observation_coordinates(:,j) = observations(d)%grid%coords(1:dimension,i)
        request%observation_values(j) = observations(d)%values(dataval_idx(d),i)
        if (allocated(observations(d)%variance)) then
          request%observation_variance(j) = observations(d)%variance(dataval_idx(d),i)
        end if
        request%observation_dataset_id(j) = observations(d)%id
        request%observation_class_id(j) = class_idx(d)
        request%observation_variable_id(j) = d
        request%observation_source_index(j) = i
      end do
      deallocate(pack_mask)
    end do

    request%dimension = dimension
    request%observation_total_count = total_count
    request%category_id = category_id
    if (size(class_idx) == 1) request%class_id = class_idx(1)
  end subroutine pack_layerdata_observations

!-------------------------------------------------------------------------------------------------!

  subroutine pack_pilotpoint_observations( &
    request, observations, parameter_idx, class_id, category_id, dataset_id)
    type(t_interp_request), intent(inout) :: request
    type(t_pilotpoints), intent(in)       :: observations
    integer, intent(in)                   :: parameter_idx, class_id, category_id
    integer, intent(in), optional         :: dataset_id
    logical, allocatable                  :: pack_mask(:)
    integer                               :: i, j, nobs, source_dataset_id

    if (.not. associated(observations%grid)) then
      call packer_error('pilot-point grid is not associated')
    end if
    if (.not. allocated(observations%values) .or. &
        .not. allocated(observations%category)) then
      call packer_error('pilot-point values or categories are not allocated')
    end if
    if (size(observations%category) /= observations%grid%n .or. &
        size(observations%values,1) /= observations%grid%n) then
      call packer_error('pilot-point arrays do not match the source grid')
    end if
    if (parameter_idx < 1 .or. parameter_idx > size(observations%values,2) .or. &
        class_id < 1 .or. class_id > size(observations%values,3)) then
      call packer_error('pilot-point parameter or class index is out of range')
    end if

    allocate(pack_mask(observations%grid%n))
    pack_mask = ((observations%category == category_id) .or. &
                 (observations%category == CAT_ALL)) .and. &
                (observations%values(:,parameter_idx,class_id) /= NODATA)
    nobs = count(pack_mask)

    source_dataset_id = 1
    if (present(dataset_id)) source_dataset_id = dataset_id

    call clear_observations(request)
    allocate(request%observation_coordinates(observations%grid%dim, nobs))
    allocate(request%observation_values(nobs))
    allocate(request%observation_dataset_id(nobs))
    allocate(request%observation_class_id(nobs))
    allocate(request%observation_variable_id(nobs))
    allocate(request%observation_source_index(nobs))

    j = 0
    do i = 1, observations%grid%n
      if (.not. pack_mask(i)) cycle
      j = j + 1
      request%observation_coordinates(:,j) = observations%grid%coords(:,i)
      request%observation_values(j) = observations%values(i,parameter_idx,class_id)
      request%observation_dataset_id(j) = source_dataset_id
      request%observation_class_id(j) = class_id
      request%observation_variable_id(j) = 1
      request%observation_source_index(j) = i
    end do

    request%dimension = observations%grid%dim
    request%observation_total_count = observations%grid%n
    request%category_id = category_id
    request%class_id = class_id
  end subroutine pack_pilotpoint_observations

!-------------------------------------------------------------------------------------------------!

  subroutine fill_category_value_mask(observations, dataval_idx, category_id, mask)
    type(t_layerdata), intent(in) :: observations
    integer, intent(in)           :: dataval_idx, category_id
    logical, intent(out)          :: mask(:)
    integer                       :: i

    do i = 1, size(mask)
      mask(i) = ((observations%category(i) == category_id) .or. &
                 (observations%category(i) == CAT_ALL)) .and. &
                (observations%values(dataval_idx,i) /= NODATA)
    end do
  end subroutine fill_category_value_mask

!-------------------------------------------------------------------------------------------------!

  subroutine check_layerdata(observations, dataval_idx, dimension)
    type(t_layerdata), intent(in) :: observations
    integer, intent(in)           :: dataval_idx, dimension

    if (.not. associated(observations%grid)) then
      call packer_error('layerdata grid is not associated')
    end if
    if (.not. allocated(observations%values) .or. &
        .not. allocated(observations%category)) then
      call packer_error('layerdata values or categories are not allocated')
    end if
    if (dimension > observations%grid%dim) then
      call packer_error('observation dimension is incompatible with a source grid')
    end if
    if (dataval_idx < 1 .or. dataval_idx > size(observations%values,1)) then
      call packer_error('data-value index is out of range')
    end if
    if (size(observations%category) /= observations%grid%n .or. &
        size(observations%values,2) /= observations%grid%n) then
      call packer_error('layerdata arrays do not match the source grid')
    end if
    if (allocated(observations%variance)) then
      if (size(observations%variance,1) < dataval_idx .or. &
          size(observations%variance,2) /= observations%grid%n) then
        call packer_error('layerdata variance does not match the values array')
      end if
    end if
  end subroutine check_layerdata

!-------------------------------------------------------------------------------------------------!

  subroutine clear_observations(request)
    type(t_interp_request), intent(inout) :: request

    if (allocated(request%observation_coordinates)) deallocate(request%observation_coordinates)
    if (allocated(request%observation_values)) deallocate(request%observation_values)
    if (allocated(request%observation_variance)) deallocate(request%observation_variance)
    if (allocated(request%observation_dataset_id)) deallocate(request%observation_dataset_id)
    if (allocated(request%observation_class_id)) deallocate(request%observation_class_id)
    if (allocated(request%observation_variable_id)) deallocate(request%observation_variable_id)
    if (allocated(request%observation_source_index)) deallocate(request%observation_source_index)
    if (allocated(request%variable_class_id)) deallocate(request%variable_class_id)
    if (allocated(request%variable_mean)) deallocate(request%variable_mean)
  end subroutine clear_observations

!-------------------------------------------------------------------------------------------------!

  subroutine clear_targets(request)
    type(t_interp_request), intent(inout) :: request

    if (allocated(request%target_coordinates)) deallocate(request%target_coordinates)
    if (allocated(request%target_index)) deallocate(request%target_index)
    if (allocated(request%target_block_size)) deallocate(request%target_block_size)
    if (allocated(request%target_active)) deallocate(request%target_active)
    request%request_block_support = .false.
  end subroutine clear_targets

!-------------------------------------------------------------------------------------------------!

  subroutine packer_error(message)
    character(*), intent(in) :: message

    call error_handler(9, opt_msg='Interpolation packer: '//trim(message))
  end subroutine packer_error

end module m_interpolation_packers
