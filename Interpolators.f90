module m_interpolator
  use m_weights
  use m_grid
  use m_vario
  implicit none
!-------------------------------------------------------------------------------------------------!
  type, abstract :: Interpolator
  ! Abstract Base Class
  contains
    procedure              :: calculate_weights
    procedure              :: apply_weights
  end type Interpolator
!-------------------------------------------------------------------------------------------------!
  type, extends(Interpolator) :: Kriging
    class(t_vgm), allocatable :: variogramModel
    ! Attributes for KD-tree and search parameters could be added here
 ! contains

  end type Kriging
!-------------------------------------------------------------------------------------------------!
  type, extends(Kriging) :: SK
  contains
    procedure              :: calculate_weights => calculate_weights_sk
    procedure              :: apply_weights => apply_weights_sk
  end type SK
!-------------------------------------------------------------------------------------------------!
  type, extends(Kriging) :: OK
  contains
    procedure              :: calculate_weights => calculate_weights_ok
    procedure              :: apply_weights => apply_weights_ok
  end type OK
!-------------------------------------------------------------------------------------------------!
  type, extends(Interpolator) :: IDW
  contains
    procedure              :: calculate_weights => calculate_weights_idw
    procedure              :: apply_weights => apply_weights_idw
  end type IDW
!-------------------------------------------------------------------------------------------------!

  contains

!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! BASE CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  subroutine calculate_weights(this, inputGrids, outputGrid, dataArrays)
    class(Interpolator), intent(in) :: this
    type(t_grid), intent(in) :: inputGrids(:)
    type(t_grid), intent(in) :: outputGrid
    real(dp), intent(in) :: dataArrays(:)
    ! Implementation...
  end subroutine calculate_weights
!-------------------------------------------------------------------------------------------------!
  subroutine apply_weights(this, dataArrays, weightsObject, outputValues)
    class(Interpolator), intent(in) :: this
    real(dp), intent(in) :: dataArrays(:)
    type(t_weights), intent(in) :: weightsObject
    real, intent(out) :: outputValues(:)
    ! Implementation...
  end subroutine apply_weights
  
!-------------------------------------------------------------------------------------------------!
! KRIGING CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! SIMPLE KRIGING CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine calculate_weights_sk(this, inputGrids, outputGrid, dataArrays)
    class(SK), intent(in) :: this
    type(t_grid), intent(in) :: inputGrids(:)
    type(t_grid), intent(in) :: outputGrid
    real(dp), intent(in) :: dataArrays(:)
    ! Implementation...
  end subroutine calculate_weights_sk
!-------------------------------------------------------------------------------------------------!
  subroutine apply_weights_sk(this, dataArrays, weightsObject, outputValues)
    class(SK), intent(in) :: this
    real(dp), intent(in) :: dataArrays(:)
    type(t_weights), intent(in) :: weightsObject
    real, intent(out) :: outputValues(:)
    ! Implementation...
  end subroutine apply_weights_sk

!-------------------------------------------------------------------------------------------------!
! ORDINARY KRIGING CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine calculate_weights_ok(this, inputGrids, outputGrid, dataArrays)
    class(OK), intent(in) :: this
    type(t_grid), intent(in) :: inputGrids(:)
    type(t_grid), intent(in) :: outputGrid
    real(dp), intent(in) :: dataArrays(:)
    ! Implementation...
  end subroutine calculate_weights_ok
!-------------------------------------------------------------------------------------------------!
  subroutine apply_weights_ok(this, dataArrays, weightsObject, outputValues)
    class(OK), intent(in) :: this
    real(dp), intent(in) :: dataArrays(:)
    type(t_weights), intent(in) :: weightsObject
    real, intent(out) :: outputValues(:)
    ! Implementation...
  end subroutine apply_weights_ok

!-------------------------------------------------------------------------------------------------!
! INVERSE DISTANCE WEIGHTING (IDW) CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine calculate_weights_idw(this, inputGrids, outputGrid, dataArrays)
    class(IDW), intent(in) :: this
    type(t_grid), intent(in) :: inputGrids(:)
    type(t_grid), intent(in) :: outputGrid
    real(dp), intent(in) :: dataArrays(:)
    ! Implementation...
  end subroutine calculate_weights_idw
!-------------------------------------------------------------------------------------------------!
  subroutine apply_weights_idw(this, dataArrays, weightsObject, outputValues)
    class(IDW), intent(in) :: this
    real(dp), intent(in) :: dataArrays(:)
    type(t_weights), intent(in) :: weightsObject
    real, intent(out) :: outputValues(:)
    ! Implementation...
  end subroutine apply_weights_idw

!-------------------------------------------------------------------------------------------------!
end module m_interpolator