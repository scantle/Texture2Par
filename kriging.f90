module m_kriging
  use m_global, only: dp
  use m_interpolator
  use m_weights, only: t_weights, create_weight_object
  use m_grid, only: t_grid, create_grid_object
  use m_vario, only: t_vgm
  use m_datasets, only: t_dataset, t_layerdata
  use m_error_handler, only: error_handler
  use m_options, only: t_options
  implicit none
!-------------------------------------------------------------------------------------------------!

  type, extends(t_interpolator) :: t_krige
    integer                          :: nvario            ! Number of variograms stored for interpolator
    integer,pointer                  :: nnear(:)          ! Number of nearest points used for interpolation
    integer,pointer                  :: nstruct(:)        ! Number of structures for each variogram
    type(t_vgm),pointer              :: variograms(:,:)   ! (stuctures, variogram)
    ! Attributes for KD-tree and search parameters could be added here
  contains
    procedure                 :: intialize => initialize_kriging

  end type t_krige
!-------------------------------------------------------------------------------------------------!
  
  type, extends(t_krige) :: t_simpkrige
  contains
    !procedure              :: calc_weights => calc_weights_sk
    !procedure              :: apply_weights => apply_weights_sk
  end type t_simpkrige
!-------------------------------------------------------------------------------------------------!
  
  type, extends(t_krige) :: t_ordkrige
  contains
    !procedure              :: calc_weights => calc_weights_ok
    !procedure              :: apply_weights => apply_weights_ok
  end type t_ordkrige
!-------------------------------------------------------------------------------------------------!

  contains

!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------!
! KRIGING CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine initialize_kriging(this)
    implicit none
    class(t_krige)                :: this
    
    
  end subroutine initialize_kriging
  
  
!-------------------------------------------------------------------------------------------------!
    
!-------------------------------------------------------------------------------------------------!
! SIMPLE KRIGING CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

    subroutine calc_weights_sk(this, from_grid, to_grid, weights)
    implicit none
    class(t_simpkrige), intent(in)  :: this
    type(t_grid),intent(in)         :: from_grid, to_grid
    type(t_weights)                 :: weights
    
  end subroutine calc_weights_sk
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! ORDINARY KRIGING CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

    subroutine calc_weights_ok(this, from_grid, to_grid, weights)
    implicit none
    class(t_ordkrige), intent(in)   :: this
    type(t_grid),intent(in)         :: from_grid, to_grid
    type(t_weights)                 :: weights
    
  end subroutine calc_weights_ok
  
!-------------------------------------------------------------------------------------------------!
  
end module m_kriging