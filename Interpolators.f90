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
  end type t_interpolator
!-------------------------------------------------------------------------------------------------!
  
  type, extends(t_interpolator) :: t_idw
    integer                :: p
  contains
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
  
!-------------------------------------------------------------------------------------------------!
! BASE CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine initialize_abc(this)
    implicit none
    class(t_interpolator), intent(in)     :: this
        
  end subroutine initialize_abc
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine calc_weights(this, from_grid, to_grid, weights, catid)
    ! Abstract Base Class Routine - not meant to be used
    implicit none
    class(t_interpolator), intent(in)     :: this
    type(t_grid),intent(in)               :: from_grid, to_grid
    type(t_weights)                       :: weights
    
    integer,intent(in)                    :: catid
    
  end subroutine calc_weights
  
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

    subroutine calc_weights_idw(this, from_grid, to_grid, weights)
    implicit none
    class(t_idw), intent(in)        :: this
    type(t_grid),intent(in)         :: from_grid, to_grid
    type(t_weights)                 :: weights
    
  end subroutine calc_weights_idw
  
!-------------------------------------------------------------------------------------------------!
end module m_interpolator