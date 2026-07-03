! Backend-neutral batch interpolation data types.
!
! These types establish the seam used by interpolation backends without
! depending on Texture2Par's dataset, grid, or point-wise interpolator classes.
! They remain data-only while carrying the current libkriging runtime requests.
! Real components intentionally use Texture2Par's prevailing default real
! kind. Any explicit kind conversion belongs at the backend adapter.
module m_interpolation_types
  implicit none
  private

  integer, parameter, public :: INTERP_PURPOSE_UNSPECIFIED = 0
  integer, parameter, public :: INTERP_PURPOSE_MODEL       = 1
  integer, parameter, public :: INTERP_PURPOSE_PILOT_POINT = 2
  integer, parameter, public :: INTERP_PURPOSE_INTERBED    = 3

  integer, parameter, public :: PREDICT_ESTIMATE = 1
  integer, parameter, public :: PREDICT_SGSIM    = 2

  type, public :: t_search_spec
    ! Current Texture2Par search controls.
    integer :: n_nearest_primary   = 0
    integer :: n_nearest_secondary = 0
    real    :: max_primary_distance   = -1.0
    real    :: max_secondary_distance = -1.0
    logical :: anisotropic_search = .false.

    ! Reserved for a later libkriging sector-search input extension.
    logical :: sector_search    = .false.
    integer :: n_sectors        = 0
    integer :: min_per_sector   = 0
    integer :: max_per_sector   = 0
  end type t_search_spec

  type, public :: t_interp_request
    ! Problem metadata. Zero denotes not applicable or not yet assigned.
    integer :: dimension   = 0
    integer :: model_layer = 0
    integer :: category_id = 0
    integer :: zone_id     = 0
    integer :: class_id    = 0
    integer :: purpose     = INTERP_PURPOSE_UNSPECIFIED

    ! Prediction and simulation configuration.  SGSIM uses the same packed
    ! observations, targets, categories, variograms, and search specification
    ! as an estimate; only the backend solve mode and result extraction differ.
    integer :: prediction_mode = PREDICT_ESTIMATE
    integer :: simulation_count = 0
    integer :: random_seed = 0

    type(t_search_spec) :: search

    ! Coordinates use (dimension, point) storage to match existing T2P grids.
    real, allocatable :: observation_coordinates(:,:)
    real, allocatable :: observation_values(:)
    real, allocatable :: observation_variance(:)
    integer, allocatable :: observation_dataset_id(:)
    integer, allocatable :: observation_class_id(:)
    integer, allocatable :: observation_source_index(:)

    ! The preferred path packs only targets that should be interpolated.
    real, allocatable :: target_coordinates(:,:)
    ! Maps each packed target back to its full model/grid/node index.
    integer, allocatable :: target_index(:)
    ! Optional half-dimensions for native libkriging Gaussian block targets.
    real, allocatable :: target_block_size(:,:)
    ! Reserved for a future full-grid request; not required for packed targets.
    logical, allocatable :: target_active(:)

    ! Backend capability requests.
    logical :: request_block_support   = .false.
    logical :: request_reusable_weights = .false.
    integer :: weight_cache_id = 0
    logical :: store_weight = .false.
    logical :: use_old_weight = .false.
    character(len=1024) :: weight_file = ""
  end type t_interp_request

  type, public :: t_interp_result
    real, allocatable :: estimates(:)
    real, allocatable :: variance(:)

    ! Optional per-target status and mapping back to model nodes/cells.
    integer, allocatable :: target_status(:)
    integer, allocatable :: target_index(:)
    integer :: status = 0
  end type t_interp_result

end module m_interpolation_types
