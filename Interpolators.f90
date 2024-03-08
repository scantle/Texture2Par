module m_interpolator
  use m_global, only: dp
  use m_weights, only: t_weights, create_weight_object
  use m_grid, only: t_grid, create_grid_object
  use m_datasets, only: t_dataset, t_layerdata
  use m_error_handler, only: error_handler
  use m_options, only: t_options
  implicit none
!-------------------------------------------------------------------------------------------------!
  type, abstract :: Interpolator
  ! Abstract Base Class
  contains
    procedure,pass         :: initialize       ! placeholder
    procedure              :: predict
    procedure              :: predict_coords
    procedure              :: predict_grouped
    procedure              :: calculate_weights
    procedure              :: apply_weights => apply_weights_t_layerdata, apply_weights_real
  end type Interpolator
!-------------------------------------------------------------------------------------------------!
  
  type, extends(Interpolator) :: IDW
    integer                 :: p
  contains
    procedure              :: initalize => initialize_idw
    procedure              :: calculate_weights => calculate_weights_idw
    !procedure              :: apply_weights => apply_weights_idw
  end type IDW
  
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
! BASE CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine initialize(this)
    implicit none
    class(Interpolator), intent(in)     :: this
        
  end subroutine initialize
  
!-------------------------------------------------------------------------------------------------!
    
!-------------------------------------------------------------------------------------------------!
  
  subroutine predict(this, input_grids, output_grid, in_values, out_values, weights)
    implicit none
    class(Interpolator), intent(in)     :: this
    type(t_grid), intent(inout)         :: input_grids(:)
    type(t_grid), intent(inout)         :: output_grid
    type(t_layerdata),intent(inout)     :: in_values(:)
    real(dp),intent(out)                :: out_values(:)
    type(t_weights),pointer             :: weights
    
    ! Get weights
    if (allocated(weights%values)==.false.) then  ! if this doesn't work, try checking if weights%values is allocated
      weights = this%calculate_weights(input_grids, output_grid)
    end if
    
    ! Use weights to make prediction
    out_values = this%apply_weights(in_values, weights)
  
  end subroutine predict

!-------------------------------------------------------------------------------------------------!
  
  subroutine predict_coords(this, input_grids, coords, in_values, out_values)
    implicit none
    class(Interpolator), intent(in)     :: this
    type(t_grid), intent(inout)         :: input_grids(:)
    real(dp)                            :: coords(:,:)
    real(dp),intent(inout)              :: in_values(:)
    real(dp),intent(out)                :: out_values(:)
    type(t_weights)                     :: weights
    type(t_grid),pointer                :: temp
    
    ! A little awkward - overloading could be used if this turns out to be SLOW
    temp = create_grid_object(size(coords,1), size(coords,2))
    temp%coords = coords
    
    ! Get weights
    weights = this%calculate_weights(input_grids, temp)
    
    ! Use weights to make prediction
    out_values = this%apply_weights_real(in_values, weights)
  
  end subroutine predict_coords
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine predict_grouped(this, input_grids, output_grid, in_values, out_values, group_ids, groups, weights)
    implicit none
    class(Interpolator), intent(in)     :: this
    type(t_grid), intent(inout)         :: input_grids(:)
    type(t_grid), intent(inout)         :: output_grid
    type(t_layerdata),intent(inout)     :: in_values(:)
    real(dp),intent(out)                :: out_values(:)
    integer,intent(in)                  :: group_ids(:), groups(:)
    type(t_weights),optional            :: weights
    
    ! Get weights
    !if (present(weights)==.false.) then
    !  weights = this%calculate_weights(input_grids%coords, output_grid%coords)
    !end if
    
    ! Use weights to make prediction
    out_values = this%apply_weights(in_values, weights)
  
  end subroutine predict_grouped
  
!-------------------------------------------------------------------------------------------------!
  
  function calculate_weights(this, input_grids, output_grid) result(weights)
    ! Abstract base class method - should not be used
    implicit none
    class(Interpolator)             :: this
    type(t_grid), intent(in)        :: input_grids(:)
    type(t_grid), intent(in)        :: output_grid
    type(t_weights)                 :: weights
    
    ! Initialize, Default values
    weights = create_weight_object(sum(input_grids(:)%n), output_grid%n)
    weights%values = 0
    
  end function calculate_weights
  
!-------------------------------------------------------------------------------------------------!

  function apply_weights_t_layerdata(this, in_values, weights) result (out_values)
    implicit none
    class(Interpolator)             :: this
    type(t_layerdata),intent(inout) :: in_values(:)
    type(t_weights), intent(in)     :: weights
    real(dp)                        :: out_values(weights%nout)
    integer                         :: i, start_col, end_col
    
    ! Initialize
    out_values = 0.0_dp
    start_col = 1
    do i = 1, size(in_values)
      end_col = start_col + size(in_values(i)%values) - 1

      ! Aggregate the result into out_values
      out_values = out_values + matmul(weights%values(:, start_col:end_col), in_values(i)%values)

      start_col = end_col + 1
    end do
    
  end function apply_weights_t_layerdata
  
!-------------------------------------------------------------------------------------------------!

  function apply_weights_real(this, in_values, weights) result (out_values)
    implicit none
    class(Interpolator)             :: this
    real(dp),intent(inout)          :: in_values(:)
    type(t_weights), intent(in)     :: weights
    real(dp)                        :: out_values(weights%nout)

    ! Do
    out_values = matmul(weights%values, in_values)
    
  end function apply_weights_real
  
!-------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------!
! INVERSE DISTANCE WEIGHTING (IDW) CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine initialize_idw(this, opt)
    implicit none
    class(IDW)                    :: this
    class(t_options)              :: opt
    
    this%p = opt%IDW_POWER
    
  end subroutine initialize_idw
  
!-------------------------------------------------------------------------------------------------!
  
  function calculate_weights_idw(this, input_grids, output_grid) result(weights)
    implicit none
    class(IDW)                      :: this
    type(t_grid), intent(in)        :: input_grids(:)
    type(t_grid), intent(in)        :: output_grid
    type(t_weights)                 :: weights
    real(dp)                        :: dist
    integer :: i, j
    
    ! This method can only handle one grid
    if (size(input_grids, 1)) call error_handler(9,opt_msg="Too many grids for IDW - can only handle one input grid")
    
    ! Create weights object
    weights = create_weight_object(sum(input_grids(:)%n), output_grid%n)
    
    ! Loop over each output point
    do i = 1, size(output_grid%coords, 2)
      ! Loop over each input point
      do j = 1, size(input_grids(1)%coords, 2)
        ! Calculate distance, max to avoid divide by zero
        dist = max(dist_2D(output_grid%coords(:,i), input_grids(1)%coords(:,j)), 1.0e-8)
        ! Compute weight
        weights%values(i, j) = 1.0 / (dist**this%p)
      end do
    end do

    ! Normalize weights for each output point so they sum to 1
    do i = 1, size(weights%values, 1)
      weights%values(i, :) = weights%values(i, :) / sum(weights%values(i, :))
    end do
  end function calculate_weights_idw
  
!-------------------------------------------------------------------------------------------------!
end module m_interpolator