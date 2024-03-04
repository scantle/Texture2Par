module m_interpolator
  use m_weights
  use m_grid
  use m_vario
  implicit none
!-------------------------------------------------------------------------------------------------!
  type, abstract :: Interpolator
  ! Abstract Base Class
  contains
    procedure,public       :: predict
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
    !procedure              :: calculate_weights => calculate_weights_sk
    !procedure              :: apply_weights => apply_weights_sk
  end type SK
!-------------------------------------------------------------------------------------------------!
  
  type, extends(Kriging) :: OK
  contains
    !procedure              :: calculate_weights => calculate_weights_ok
    !procedure              :: apply_weights => apply_weights_ok
  end type OK
!-------------------------------------------------------------------------------------------------!
  
  type, extends(Interpolator) :: IDW
  contains
    !procedure              :: calculate_weights => calculate_weights_idw
    !procedure              :: apply_weights => apply_weights_idw
  end type IDW
!-------------------------------------------------------------------------------------------------!

  contains

!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! BASE CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine predict(this, input_grids, output_grid, in_values, out_values, ngroups, groups, weights)
    implicit none
    class(Interpolator), intent(in)     :: this
    type(t_grid), intent(inout)         :: input_grids(:)
    type(t_grid), intent(inout)         :: output_grid
    real(dp),intent(inout)              :: in_values(:)       ! placeholder
    real(dp),intent(out)                :: out_values(:)
    integer ,intent(in)                 :: ngroups
    integer ,intent(in)                 :: groups(ngroups)
    type(t_weights),optional            :: weights
    
    ! Get weights
    if (present(weights)==.false.) then
      weights = this%calculate_weights(input_grids, output_grid)
    end if
    
    ! Use weights to make prediction
    out_values = this%apply_weights(in_values, weights)
  
  end subroutine predict
  
!-------------------------------------------------------------------------------------------------!
  
  function calculate_weights(this, inputGrids, outputGrid) result(weights)
    implicit none
    class(Interpolator), intent(in) :: this
    type(t_grid), intent(in)        :: inputGrids(:)
    type(t_grid), intent(in)        :: outputGrid
    type(t_weights)                 :: weights
    
    ! Implementation...
    
  end function calculate_weights
  
!-------------------------------------------------------------------------------------------------!
  
  function apply_weights(this, in_values, weights) result (out_values)
    implicit none
    class(Interpolator), intent(in) :: this
    real(dp), intent(in)            :: in_values(:)
    type(t_weights), intent(in)     :: weights
    real(dp)                        :: out_values(weights%nout)
    
    ! Implementation...
    
  end function apply_weights
  
!-------------------------------------------------------------------------------------------------!
! KRIGING CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! SIMPLE KRIGING CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  

  
!-------------------------------------------------------------------------------------------------!
  


!-------------------------------------------------------------------------------------------------!
! ORDINARY KRIGING CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  

  
!-------------------------------------------------------------------------------------------------!
  


!-------------------------------------------------------------------------------------------------!
! INVERSE DISTANCE WEIGHTING (IDW) CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  

  
!-------------------------------------------------------------------------------------------------!
  


!-------------------------------------------------------------------------------------------------!
end module m_interpolator