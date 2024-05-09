module m_interpolator
  use m_global, only: dp
  use m_weights, only: t_weights, create_weight_object
  use m_grid, only: t_grid, create_grid_object
  use m_datasets, only: t_dataset, t_layerdata
  use m_error_handler, only: error_handler
  use m_options, only: t_options
  use m_vario, only: t_vgm
  use m_categories, only: t_category
  implicit none
!-------------------------------------------------------------------------------------------------!
  type, abstract :: t_interpolator
  ! Abstract Base Class
  contains
    procedure,private      :: initialize_abc
    generic                :: initialize => initialize_abc       ! placeholder
    procedure,public       :: calc_weights
    procedure,private      :: setup
    procedure,private      :: predict_point
  end type t_interpolator
!-------------------------------------------------------------------------------------------------!
  
  type, extends(t_interpolator) :: t_idw
    real(dp)               :: power
    integer                :: nobs
    logical,allocatable    :: catmask(:)
    
  contains
    procedure               :: setup => setup_idw
    procedure               :: predict_point => predict_point_idw
    !procedure              :: initialize_idw
    !generic                :: initialize => initialize_idw
    !procedure,public       :: calc_weights => calc_weights_idw
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

    ! Calculate the differences using appropriate indexing
    distances = sqrt((point_array(1, :) - p(1)) ** 2 + &
                      (point_array(2, :) - p(2)) ** 2)
  end function dist_2D_vector
  
!-------------------------------------------------------------------------------------------------!
! BASE CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine initialize_abc(this)
    implicit none
    class(t_interpolator), intent(in)     :: this
        
  end subroutine initialize_abc
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine calc_weights(this, from_grids, to_grid, weights, catlist, from_cats, to_cats, opt)
    implicit none
    class(t_interpolator), intent(in)     :: this
    type(t_grid),intent(in)               :: from_grids(:), to_grid
    real(dp)                              :: weights(:,:)
    integer,intent(in)                    :: catlist(:), from_cats(:), to_cats(:)
    type(t_options),intent(in)            :: opt
    integer                               :: i, j, cat
    
    write(*,*) 'meely meely'
    
    ! Need to do anything to setup weights?
    ! Make sure they're zero'd?
    
    do i=1, size(catlist)
      ! Setup for each category
      cat = catlist(i)
      call this%setup(from_grids, to_grid, cat, from_cats, to_cats, weights, opt)
      do j=1, to_grid%n
        ! Interpolate for each point
        if (to_cats(j)==cat) call this%predict_point(j, to_grid%coords(1:2,j), from_grids, cat, to_cats, weights(:,j))
      end do
    end do
    
    !call this%cleanup() needed?
    
  end subroutine calc_weights
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine setup(this, from_grids, to_grid,  cat, from_cats, to_cats, weights, opt)
    implicit none
    class(t_interpolator)                 :: this
    type(t_grid),intent(in)               :: from_grids(:), to_grid
    real(dp)                              :: weights(:,:)
    integer,intent(in)                    :: cat
    integer,intent(in)                    :: from_cats(:), to_cats(:)
    type(t_options),intent(in)            :: opt
    
  end subroutine setup
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine predict_point(this, idx, p, from_grids, cat, from_cats, weights)
    implicit none
    class(t_interpolator), intent(in) :: this
    integer,intent(in)                :: idx, cat
    real(dp),intent(in)               :: p(2)
    type(t_grid),intent(in)           :: from_grids(:)
    integer,intent(in)                :: from_cats(:)
    real(dp)                          :: weights(:)
  
  end subroutine predict_point
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! INVERSE DISTANCE WEIGHTING (IDW) CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  !subroutine initialize_idw(this)
  !  implicit none
  !  class(t_idw)                  :: this
  !  
  !end subroutine initialize_idw
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine setup_idw(this, from_grids, to_grid, cat, from_cats, to_cats, weights, opt)
    implicit none
    class(t_idw)                  :: this
    type(t_grid),intent(in)       :: from_grids(:), to_grid
    real(dp)                      :: weights(:,:)
    integer,intent(in)            :: cat
    integer,intent(in)            :: from_cats(:), to_cats(:)
    type(t_options),intent(in)    :: opt
    
    ! Do elsewhere?
    if (size(from_grids) > 1) call error_handler(4,opt_msg="IDW recieved multiple grids... covariates not supported")
    
    this%nobs = from_grids(1)%n
    this%power = opt%IDW_POWER
    
    allocate(this%catmask(this%nobs))
    this%catmask = (from_cats == cat)
    
  end subroutine setup_idw

!-------------------------------------------------------------------------------------------------!
  
  subroutine predict_point_idw(this, idx, p, from_grids, cat, from_cats, weights)
    implicit none
    class(t_idw), intent(in)      :: this
    integer,intent(in)            :: idx, cat
    real(dp),intent(in)           :: p(2)
    type(t_grid),intent(in)       :: from_grids(:)
    integer,intent(in)            :: from_cats(:)
    real(dp)                      :: weights(:)
    real(dp)                      :: distances(this%nobs)
    real(dp)                      :: total_weight
    integer                       :: zero_count
    
    ! Calculate distances for points matching the mask
    distances = huge(dp)
    distances(this%catmask) = dist_2D_vector(p, from_grids(1)%coords(:, this%catmask))
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
    
  end subroutine predict_point_idw
!-------------------------------------------------------------------------------------------------!
end module m_interpolator