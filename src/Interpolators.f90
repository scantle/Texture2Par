module m_interpolator
  use m_weights, only: t_weights
  use m_error_handler, only: error_handler
  use m_options, only: t_options
  use m_interpolation_types, only: t_interp_request, t_interp_result, INTERP_PURPOSE_MODEL
  implicit none

!-------------------------------------------------------------------------------------------------!
! Every interpolation backend accepts packed geometry and one or more value columns. Model/data,
! pilot-point, and interbed workflows all use this interface.
!-------------------------------------------------------------------------------------------------!
  type, abstract :: t_interpolator
  contains
    procedure, private :: initialize_abc
    generic :: initialize => initialize_abc
    procedure(interpolate_batch_ifc), deferred, public :: interpolate
    procedure :: write_summary
    procedure :: finalize
  end type t_interpolator

  type, extends(t_interpolator) :: t_idw
    real, allocatable :: cached_weights(:,:)
    character(len=1024) :: cached_weight_file = ""
    logical :: cached_weights_loaded = .false.
  contains
    procedure :: interpolate => interpolate_idw
    procedure :: finalize => finalize_idw
  end type t_idw

  abstract interface
    subroutine interpolate_batch_ifc(this, request, values, opt, result)
      import :: t_interpolator, t_interp_request, t_interp_result, t_options
      class(t_interpolator), intent(inout) :: this
      type(t_interp_request), intent(inout) :: request
      real, intent(in) :: values(:,:)
      type(t_options), intent(in) :: opt
      type(t_interp_result), intent(inout) :: result
    end subroutine interpolate_batch_ifc
  end interface

contains

!-------------------------------------------------------------------------------------------------!

  function dist(point1, point2) result(distance)
    real, intent(in) :: point1(:), point2(:)
    real :: distance

    distance = sqrt(sum((point2 - point1)**2))
  end function dist

!-------------------------------------------------------------------------------------------------!

  function dist_vector(p, point_array) result(distances)
    real, intent(in) :: p(:)
    real, intent(in) :: point_array(:,:)
    real :: distances(size(point_array,2))
    integer :: i

    distances = 0.0
    do i = 1, size(p)
      distances = distances + (p(i) - point_array(i,:))**2
    end do
    distances = sqrt(distances)
  end function dist_vector

!-------------------------------------------------------------------------------------------------!

  subroutine initialize_abc(this)
    class(t_interpolator), intent(in) :: this
  end subroutine initialize_abc

!-------------------------------------------------------------------------------------------------!

  subroutine write_summary(this)
    class(t_interpolator), intent(inout) :: this
  end subroutine write_summary

!-------------------------------------------------------------------------------------------------!

  subroutine finalize(this)
    class(t_interpolator), intent(inout) :: this
  end subroutine finalize

!-------------------------------------------------------------------------------------------------!

  subroutine interpolate_idw(this, request, values, opt, result)
    class(t_idw), intent(inout) :: this
    type(t_interp_request), intent(inout) :: request
    real, intent(in) :: values(:,:)
    type(t_options), intent(in) :: opt
    type(t_interp_result), intent(inout) :: result
    type(t_weights) :: weights_file
    integer, allocatable :: observation_index(:), cache_index(:)
    logical, allocatable :: observation_mask(:)
    real, allocatable :: weights(:), distances(:)
    real :: total_weight
    integer :: i, j, zero_count
    logical :: reset_cache

    if (.not. allocated(request%observation_coordinates) .or. &
        .not. allocated(request%target_coordinates) .or. &
        .not. allocated(request%target_index)) &
      call error_handler(4, opt_msg='IDW batch interpolation requires observation and target coordinates')
    if (size(values,1) /= size(request%observation_coordinates,2)) &
      call error_handler(4, opt_msg='IDW batch observation geometry/value shape mismatch')

    allocate(observation_mask(size(values,1)))
    observation_mask = .true.
    if (request%purpose == INTERP_PURPOSE_MODEL .and. &
        allocated(request%observation_variable_id)) then
      ! IDW is univariate. Connected secondary variables are packed for the
      ! common model workflow but do not participate in an IDW estimate.
      observation_mask = request%observation_variable_id == 1
    end if
    observation_index = pack([(j,j=1,size(values,1))], observation_mask)
    if (size(observation_index) == 0) &
      call error_handler(4, opt_msg='IDW batch interpolation has no observations')

    if (allocated(result%estimates_many)) deallocate(result%estimates_many)
    if (allocated(result%variance_many)) deallocate(result%variance_many)
    if (allocated(result%target_index)) deallocate(result%target_index)
    allocate(result%estimates_many(size(request%target_coordinates,2),size(values,2)))
    allocate(result%variance_many(size(request%target_coordinates,2),size(values,2)))
    allocate(result%target_index(size(request%target_index)))
    result%variance_many = 0.0
    result%target_index = request%target_index
    result%status = 0

    if (request%store_weight .or. request%use_old_weight) then
      if (request%observation_total_count < 1 .or. request%target_total_count < 1 .or. &
          .not. allocated(request%observation_source_index)) &
        call error_handler(4, opt_msg='IDW saved weights require full source and target index metadata')
      cache_index = request%observation_source_index(observation_index)
      reset_cache = trim(this%cached_weight_file) /= trim(request%weight_file)
      if (allocated(this%cached_weights)) then
        reset_cache = reset_cache .or. &
          size(this%cached_weights,1) /= request%target_total_count .or. &
          size(this%cached_weights,2) /= request%observation_total_count
      end if
      if (reset_cache .or. .not. allocated(this%cached_weights)) then
        if (allocated(this%cached_weights)) deallocate(this%cached_weights)
        allocate(this%cached_weights(request%target_total_count,request%observation_total_count))
        this%cached_weights = 0.0
        this%cached_weight_file = request%weight_file
        this%cached_weights_loaded = .false.
      end if
      if (request%use_old_weight .and. .not. this%cached_weights_loaded) then
        call weights_file%initialize(request%observation_total_count,request%target_total_count)
        call weights_file%load_from_file(trim(request%weight_file),request%weight_file_text)
        this%cached_weights = weights_file%values
        this%cached_weights_loaded = .true.
      end if
    end if

    allocate(weights(size(observation_index)), distances(size(observation_index)))
    do i = 1, size(request%target_coordinates,2)
      if (request%use_old_weight) then
        weights = this%cached_weights(request%target_index(i),cache_index)
      else
        distances = dist_vector(request%target_coordinates(:,i), &
                                request%observation_coordinates(:,observation_index))
        zero_count = count(distances == 0.0)
        if (zero_count > 0) then
          weights = 0.0
          where (distances == 0.0) weights = 1.0 / zero_count
        else
          weights = 1.0 / distances**opt%IDW_POWER
          total_weight = sum(weights)
          weights = weights / total_weight
        end if
        if (request%store_weight) &
          this%cached_weights(request%target_index(i),cache_index) = weights
      end if
      result%estimates_many(i,:) = matmul(weights,values(observation_index,:))
    end do

    if (request%store_weight) then
      call weights_file%initialize(request%observation_total_count,request%target_total_count)
      weights_file%values = this%cached_weights
      call weights_file%save_to_file(trim(request%weight_file),request%weight_file_text)
    end if
  end subroutine interpolate_idw

!-------------------------------------------------------------------------------------------------!

  subroutine finalize_idw(this)
    class(t_idw), intent(inout) :: this

    if (allocated(this%cached_weights)) deallocate(this%cached_weights)
    this%cached_weight_file = ""
    this%cached_weights_loaded = .false.
  end subroutine finalize_idw

!-------------------------------------------------------------------------------------------------!
end module m_interpolator
