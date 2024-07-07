module m_interpolator
  use m_global, only: dp
  use m_weights, only: t_weights
  use m_grid, only: t_grid, create_grid_object
  use m_datasets, only: t_dataset, t_layerdata
  use m_error_handler, only: error_handler
  use m_options, only: t_options
  use m_vario, only: t_vgm
  use m_categories, only: t_category
  implicit none
!-------------------------------------------------------------------------------------------------!
! ADDING A NEW INTERPOLATOR
! 
! Interpolators must override these methods:
!  - setup(this, obs_data, obs_idx, to_grid, cat, to_cats, weights, opt)
!  - predict_point(this, idx, p, obs_data, cat, weights, widx, ppd)
! 
! Optionally, they can also override these:
!  - wsize_estimate(this, obs_data, obs_idx)
!  - finalize(this)
! 
! `setup` is called before each loop over a set of points in a given category (pp zone, HSU).
! The intent is for it to be used to allocate arrays, etc.
! 
! `predict_point` is called for each point on the output grid (MF, IWFM, GSF...). The intent
! is to return a series of a weights for each input observation (pp or class data). In the
! interpolate subroutine `calc_weights` it returns weights for the entire observation set.
! In `calc` however (currently used for class interpolation) a sparse weight representation
! is used; the predict_point() optional arguments ppd and widx are used to return the number
! of observations and the observation indices, respectively.
!
! `wsize_estimator` returns the maximum size of the weight array (max obs points used)
! `finalize` deallocates arrays and nullifies pointers
! 
! The arguments of these methods are documented in their abstract base class placeholders
! below.
! 
! To add the new interpolator as a valid option in the input file, add a new case
! statement to T2P.f90 in the `process_options` method. Their the interpolator can be setup
! to be allocated as the correct subclass. 
! 
! If your interpolator has additional options that need to be included via the input file,
! their are a couple possible routes. New options can be added to the options class in
! options.f90, and then when the interpolator `setup()` method is called it will have access
! to the new options. If the interpolator requires it's own block in the input file (as
! kriging did with variograms) you'll need to:
!  1. Add a new block in the main input file with a new reader in the Read_Main_Input.f90 
!     file
!  2. Setup the OPTIONS block to be required before your new block in `pre_block_read_check()`
!  3. Write a subroutine in T2P to be called after your block is read. Pass the necessary
!     data or options to t2p%interp_pp and t2p%interp_data - see setup_variograms in T2P.f90
! 	for an example
! Other less drastic options include:
!  1. Adding an input file to the options block your `setup()` routine can read
!     Downside: `setup()` can't distinguish between pp and data interpolation easily
!  2. Modify the interpolator `initialize()` routine to accept your parameter as an input.
!     Override `initialize()` for your interpolator so it properly handles the parameter
! 	The `initialize()` calls are currently in `process_options()` but could be moved
! 	to run after everything is read - they currently do nothing!
!-------------------------------------------------------------------------------------------------!
  type, abstract:: t_interpolator
  
  ! Abstract Base Class
  contains
    procedure,private      :: initialize_abc
    generic                :: initialize => initialize_abc       ! placeholder
    procedure,public       :: calc_weights
    procedure,public       :: calc
    procedure,private      :: single_point_full_weight
    procedure              :: setup
    procedure              :: predict_point
    procedure              :: wsize_estimate
    procedure              :: finalize
  end type t_interpolator
!-------------------------------------------------------------------------------------------------!
  
  type, extends(t_interpolator) :: t_idw
    real(dp)               :: power
    integer                :: nobs
    logical,allocatable    :: catmask(:)
    
  contains
    procedure              :: setup => setup_idw
    procedure              :: predict_point => predict_point_idw
    procedure              :: finalize => finalize_idw
  end type t_idw
  
!-------------------------------------------------------------------------------------------------!
! Kriging class/routines are located in kriging.f90. They inheret from Interpolator
!-------------------------------------------------------------------------------------------------!

  contains

!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  function dist_2D(point1, point2) result(distance)
    real(dp), intent(in) :: point1(2), point2(2)  ! Input: Coordinates of the two points
    real(dp)             :: distance              ! Output: Euclidean distance

    distance = sqrt((point2(1) - point1(1))**2 + (point2(2) - point1(2))**2)
  end function dist_2D
  
!-------------------------------------------------------------------------------------------------!
  
  function dist_2D_vector(p, point_array) result(distances)
    real(dp), intent(in)                      :: p(2)          
    real(dp), dimension(:, :), intent(in)     :: point_array
    real(dp), dimension(size(point_array, 2)) :: distances
    distances = sqrt((point_array(1, :) - p(1)) ** 2 + (point_array(2, :) - p(2)) ** 2)
  end function dist_2D_vector
  
!-------------------------------------------------------------------------------------------------!
! BASE CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine initialize_abc(this)
    implicit none
    class(t_interpolator), intent(in)     :: this
        
  end subroutine initialize_abc
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine calc_weights(this, obs_data, obs_idx, to_grid, weights, catlist, to_cats, opt)
    implicit none
    class(t_interpolator),intent(inout)   :: this
    class(t_dataset)                      :: obs_data(:)
    integer,intent(in)                    :: obs_idx(:)
    type(t_grid),intent(in)               :: to_grid
    real(dp)                              :: weights(:,:) !(output_points, input_points)
    integer,intent(in)                    :: catlist(:), to_cats(:)
    type(t_options),intent(in)            :: opt
    integer                               :: i, j, cat, count_cat
    
    do i=1, size(catlist)
      ! Setup for each category
      cat = catlist(i)
      count_cat = 0
      ! Count occurrences across datasets
      do j = 1, size(obs_data)
        count_cat = count_cat + count(obs_data(j)%category == cat)
        if (count_cat > 1) exit
      end do
      
      if (count_cat>1) then ! more than one obs point
        call this%setup(obs_data, obs_idx, to_grid, cat, to_cats, weights, opt)
        do j=1, to_grid%n
          ! Interpolate for each point
          if (to_cats(j)==cat) call this%predict_point(j, to_grid%coords(1:2,j), obs_data, cat, weights(j,:))
        end do
      else  ! only one point!
        call this%single_point_full_weight(obs_data, cat, to_cats, weights)
      end if
    end do
    
  end subroutine calc_weights
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine calc(this, obs_data, obs_idx, val_idx, to_grid, results, catlist, to_cats, opt)
    use m_weights, only: t_sosa_weights
    implicit none
    class(t_interpolator),intent(inout)   :: this
    class(t_layerdata)                    :: obs_data(:)
    integer,intent(in)                    :: obs_idx(:), val_idx(:)
    type(t_grid),intent(in)               :: to_grid
    real(dp)                              :: results(:)
    integer,intent(in)                    :: catlist(:), to_cats(:)
    type(t_options),intent(in)            :: opt
    type(t_sosa_weights)                  :: weights
    integer                               :: i, j, k, woffset, cat, wlen
    real(dp)                              :: fake_weights(1,1), v_adjust
    integer                               :: ppd(size(obs_data))
    
    results = 0.0_dp  ! necessary?
    call weights%initialize(this%wsize_estimate(obs_data, obs_idx))
    
    do i=1, size(catlist)
      ! Setup for each category
      cat = catlist(i)
      call this%setup(obs_data, obs_idx, to_grid, cat, to_cats, fake_weights, opt)
      
      do j=1, to_grid%n
        ! Interpolate for each point
        if (to_cats(j)==cat) then
          call this%predict_point(j, to_grid%coords(1:2,j), obs_data, cat, weights%values, weights%idx, ppd)
          ! Calculate for each point
          woffset = 0
          do k=1, size(obs_data)
            v_adjust = 0.0_dp
            if (k>1) v_adjust = obs_data(1)%mean(i) - obs_data(k)%mean(val_idx(1))
            results(j) = results(j) + dot_product(weights%values(woffset+1:woffset+ppd(k)), obs_data(k)%values(val_idx(k), weights%idx(woffset+1:woffset+ppd(k)))+v_adjust)
            woffset = woffset + ppd(k)
          end do
        end if
      end do
    end do
    
  end subroutine calc

!-------------------------------------------------------------------------------------------------!
subroutine single_point_full_weight(this, obs_data, cat, to_cats, weights)
  implicit none
  class(t_interpolator), intent(inout) :: this
  class(t_dataset)                     :: obs_data(:)
  real(dp)                             :: weights(:,:)
  integer, intent(in)                  :: cat, to_cats(:)
  integer                              :: i, j, k, obs_idx, w_offset

  w_offset = 0

  ! Loop over datasets to find the index of the ASSUMED SINGLE observation that matches cat
  do i = 1, size(obs_data)
    obs_idx = findloc(obs_data(i)%category, cat, dim=1)
    if (obs_idx > 0) then  ! If found in dataset 'i'
      obs_idx = obs_idx + w_offset
      ! set weights for the found observation
      do k = 1, size(to_cats)
        if (to_cats(k) == cat) weights(k, obs_idx) = 1.0_dp
      end do
      exit
    endif
    w_offset = w_offset + obs_data(i)%grid%n
  end do
end subroutine single_point_full_weight
  
!-------------------------------------------------------------------------------------------------!
!> Initializes the interpolator before each loop over a set of points in a given category.
!! This is typically used to allocate necessary arrays and set up the environment.
!! 
!! @param this The interpolator object.
!! @param obs_data Observed datasets to be used for interpolation.
!! @param obs_idx Indices of the observed data points.
!! @param to_grid Grid object where the interpolated values will be placed.
!! @param cat Category data.
!! @param to_cats Categories for the target grid points.
!! @param weights Weights to be calculated during interpolation.
!! @param opt T2P Options object.
  subroutine setup(this, obs_data, obs_idx, to_grid, cat, to_cats, weights, opt)
    ! Abstract Base Class method - should not be called
    implicit none
    class(t_interpolator),intent(inout)   :: this
    class(t_dataset)                      :: obs_data(:)
    integer,intent(in)                    :: obs_idx(:)
    type(t_grid),intent(in)               :: to_grid
    real(dp)                              :: weights(:,:)
    integer,intent(in)                    :: cat
    integer,intent(in)                    :: to_cats(:)
    type(t_options),intent(in)            :: opt
    
  end subroutine setup
  
!-------------------------------------------------------------------------------------------------!
!> Predicts the value at a specific point on the output grid based on the observed data and weights.
!! This is called for each point on the output grid to compute the interpolation weights.
!! 
!! @param this The interpolator object.
!! @param idx Index of the point to predict.
!! @param p Coordinates of the point to predict.
!! @param obs_data Observed data used for prediction.
!! @param cat Category data.
!! @param weights Weights for the prediction.
!! @param widx Indices of the weights (optional).
!! @param ppd Number of observations (optional).
  subroutine predict_point(this, idx, p, obs_data, cat, weights, widx, ppd)
    ! Abstract Base Class method - should not be called
    implicit none
    class(t_interpolator),intent(inout)  :: this
    integer,intent(in)                   :: idx, cat
    real(dp),intent(in)                  :: p(2)
    class(t_dataset)                     :: obs_data(:)
    real(dp)                             :: weights(:)
    integer, optional                    :: widx(:)              ! idx associated with weights
    integer, optional                    :: ppd(size(obs_data))  ! points per dataset
  
  end subroutine predict_point
  
!-------------------------------------------------------------------------------------------------!
!> Estimates the maximum size of the weight array for the interpolation.
!! This subroutine returns the maximum number of observation points used.
!! 
!! @param this The interpolator object.
!! @param obs_data Observed data used for interpolation.
!! @param obs_idx Indices of the observed data points.
!! @return Maximum size of the weight array.
  function wsize_estimate(this, obs_data, obs_idx)
    ! To be overwritten where necessary
    implicit none
    class(t_interpolator),intent(inout)  :: this
    class(t_layerdata)                   :: obs_data(:)
    integer,intent(in)                   :: obs_idx(:)
    integer                              :: i,wsize_estimate
    wsize_estimate = 0
    do i=1, size(obs_data)
      wsize_estimate = wsize_estimate + obs_data(i)%grid%n
    end do
    return
  end function wsize_estimate
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine finalize(this)
    implicit none
    class(t_interpolator),intent(inout)  :: this
    ! Nothing for the base class
  end subroutine finalize
    
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! INVERSE DISTANCE WEIGHTING (IDW) CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine setup_idw(this, obs_data, obs_idx, to_grid, cat, to_cats, weights, opt)
    implicit none
    class(t_idw),intent(inout)    :: this
    class(t_dataset)              :: obs_data(:)
    integer,intent(in)            :: obs_idx(:)
    type(t_grid),intent(in)       :: to_grid
    real(dp)                      :: weights(:,:)
    integer,intent(in)            :: cat
    integer,intent(in)            :: to_cats(:)
    type(t_options),intent(in)    :: opt
    type(t_layerdata),pointer     :: p
    logical                       :: na_mask(obs_data(1)%grid%n)
    
    this%nobs = obs_data(1)%grid%n
    this%power = opt%IDW_POWER
    na_mask = .true.
    
    allocate(this%catmask(this%nobs))
    ! For layer dataset, ignore all NA values
    select type(p=>obs_data(1))
      class is (t_layerdata)
        na_mask = .not. any(p%values(1:, :) == -999.0, dim=1)
    end select
    this%catmask = (obs_data(1)%category == cat) .and. (na_mask)
  end subroutine setup_idw

!-------------------------------------------------------------------------------------------------!
  
  subroutine predict_point_idw(this, idx, p, obs_data, cat, weights, widx, ppd)
    implicit none
    class(t_idw),intent(inout)    :: this
    integer,intent(in)            :: idx, cat
    real(dp),intent(in)           :: p(2)
    class(t_dataset)              :: obs_data(:)
    real(dp)                      :: weights(:)
    real(dp)                      :: distances(this%nobs)
    real(dp)                      :: total_weight
    integer                       :: zero_count, i
    integer, optional             :: widx(:)              ! idx associated with weights
    integer, optional             :: ppd(size(obs_data))  ! points per dataset
    
    ! Calculate distances for points matching the mask
    distances = huge(dp)
    distances(this%catmask) = dist_2D_vector(p, obs_data(1)%grid%coords(:, this%catmask))
    where(this%catmask)
      weights = 1.0_dp / distances**this%power
    elsewhere
      weights = 0.0_dp
    end where
    
     ! Check if there are any zero distances
    zero_count = count(distances == 0.0_dp)
    if (zero_count > 0) then
      ! Assign equal weight to all zero-distance points (TODO should this produce a warning?)
      weights(:) = 0.0_dp
      where (distances == 0.0_dp)
        weights = 1.0_dp / zero_count
      end where
    else
      ! Normalize weights to sum up to 1
      total_weight = sum(weights)
      weights = weights / total_weight
    end if
    
    ! Handle widx & ppd presence
    if (present(widx)) then
      ! todo oneliner?
      do i=1,obs_data(1)%grid%n
        widx = i
      end do
    end if
    if (present(ppd)) then
      do i=1,size(obs_data)
        ppd(i) = obs_data(i)%grid%n
      end do
    end if
    
  end subroutine predict_point_idw
!-------------------------------------------------------------------------------------------------!
  
  subroutine finalize_idw(this)
    implicit none
    class(t_idw),intent(inout)  :: this
    
    if (allocated(this%catmask)) deallocate(this%catmask)
  end subroutine finalize_idw
  
!-------------------------------------------------------------------------------------------------!
end module m_interpolator