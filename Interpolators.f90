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
    procedure              :: calc_weights_by_category_pp
    generic                :: calc_weights_by_category => calc_weights_by_category_pp

  end type t_interpolator
!-------------------------------------------------------------------------------------------------!
  
  type, extends(t_interpolator) :: t_idw
    integer                :: p
  contains
    procedure              :: initialize_idw
    generic                :: initialize => initialize_idw
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
  
!-------------------------------------------------------------------------------------------------!
! BASE CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine initialize_abc(this)
    implicit none
    class(t_interpolator), intent(in)     :: this
        
  end subroutine initialize_abc
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine calc_weights_by_category_pp(this, pilot_points, output_grid, weights, output_categories, arr_id, variograms, nstruct)
    use m_pilotpoints, only: t_pilotpoints
    implicit none
    class(t_interpolator)             :: this
    type(t_pilotpoints)               :: pilot_points
    type(t_grid)                      :: output_grid
    type(t_weights)                   :: weights
    type(t_category),optional         :: output_categories
    integer,intent(in),optional       :: arr_id, nstruct
    type(t_vgm),optional              :: variograms(:)
    
    integer                           :: i
    
    do i=1, output_categories%ncat
      where (output_categories%array(:,arr_id)== output_categories%catlist(i))
        ! Passes to the method specific, grid-only routine
        !call this%calc_weights(this%pilot_points%grid, this%output_grid, weights, variograms, nstruct)
      end where
    end do
    
  end subroutine calc_weights_by_category_pp
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! INVERSE DISTANCE WEIGHTING (IDW) CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine initialize_idw(this, opt)
    implicit none
    class(t_idw)                  :: this
    class(t_options)              :: opt
    
    this%p = opt%IDW_POWER
    
  end subroutine initialize_idw
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
end module m_interpolator