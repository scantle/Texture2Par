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
  type, abstract:: t_interpolator
  
  ! Abstract Base Class
  contains
    procedure,private      :: initialize_abc
    generic                :: initialize => initialize_abc       ! placeholder
    procedure,public       :: calc_weights
    procedure,public       :: calc
    procedure              :: setup
    procedure              :: predict_point
    procedure              :: wsize_estimate
  end type t_interpolator
!-------------------------------------------------------------------------------------------------!
  
  type, extends(t_interpolator) :: t_idw
    real(dp)               :: power
    integer                :: nobs
    logical,allocatable    :: catmask(:)
    
  contains
    procedure               :: setup => setup_idw
    procedure               :: predict_point => predict_point_idw
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
    real(dp)                              :: weights(:,:)
    integer,intent(in)                    :: catlist(:), to_cats(:)
    type(t_options),intent(in)            :: opt
    integer                               :: i, j, cat
    
    do i=1, size(catlist)
      ! Setup for each category
      cat = catlist(i)
      call this%setup(obs_data, obs_idx, to_grid, cat, to_cats, weights, opt)
      do j=1, to_grid%n
        ! Interpolate for each point
        if (to_cats(j)==cat) call this%predict_point(j, to_grid%coords(1:2,j), obs_data, cat, weights(j,:))
      end do
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
    call weights%initialize(this%wsize_estimate(obs_data))
    
    do i=1, size(catlist)
      ! Setup for each category
      cat = catlist(i)
      call this%setup(obs_data, obs_idx, to_grid, cat, to_cats, fake_weights, opt)
      
      do j=1, to_grid%n
        ! Interpolate for each point
        if (to_cats(j)==cat) call this%predict_point(j, to_grid%coords(1:2,j), obs_data, cat, weights%values, weights%idx, ppd)
        ! Calculate for each point
        woffset = 0
        do k=1, size(obs_data)
          v_adjust = 0.0_dp
          if (k>1) v_adjust = obs_data(1)%mean(i) - obs_data(k)%mean(val_idx(1))
          results(j) = results(j) + dot_product(weights%values(woffset+1:woffset+ppd(k)), obs_data(k)%values(val_idx(k), weights%idx(woffset+1:woffset+ppd(k)))+v_adjust)
          woffset = woffset + ppd(k)
        end do
      end do
    end do
    
  end subroutine calc
  
!-------------------------------------------------------------------------------------------------!
  
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
  
  function wsize_estimate(this, obs_data)
    implicit none
    class(t_interpolator),intent(inout)  :: this
    class(t_layerdata)                   :: obs_data(:)
    integer                              :: i,wsize_estimate
    wsize_estimate = 0
    do i=1, size(obs_data)
      wsize_estimate = wsize_estimate + obs_data(i)%grid%n
    end do
    return
  end function wsize_estimate
  
!-------------------------------------------------------------------------------------------------!
! INVERSE DISTANCE WEIGHTING (IDW) CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  !subroutine initialize_idw(this)
  !  implicit none
  !  class(t_idw)                  :: this
  !  
  !end subroutine initialize_idw
  
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
end module m_interpolator