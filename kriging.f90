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
    
    class(t_vgm), allocatable :: vgm
    ! Attributes for KD-tree and search parameters could be added here
  contains
    procedure                 :: intialize => initialize_kriging

  end type t_krige
!-------------------------------------------------------------------------------------------------!
  
  type, extends(t_krige) :: t_simpkrige
  contains
    !procedure              :: calculate_weights => calculate_weights_sk
    !procedure              :: apply_weights => apply_weights_sk
  end type t_simpkrige
!-------------------------------------------------------------------------------------------------!
  
  type, extends(t_krige) :: t_ordkrige
  contains
    !procedure              :: calculate_weights => calculate_weights_ok
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
  
  subroutine initialize_kriging(this, vgm)
    implicit none
    class(t_krige)                :: this
    class(t_vgm)                  :: vgm
    
    this%vgm = vgm
    
  end subroutine initialize_kriging
  
  
!-------------------------------------------------------------------------------------------------!
    
!-------------------------------------------------------------------------------------------------!
! SIMPLE KRIGING CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! ORDINARY KRIGING CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------!
  
end module m_kriging