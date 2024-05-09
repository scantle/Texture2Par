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
    integer,allocatable              :: vario_idx(:,:)        ! Variogram index given (class_id, class_id)
    ! Attributes for KD-tree and search parameters could be added here
  contains
    procedure                 :: intialize => initialize_kriging

  end type t_krige
!-------------------------------------------------------------------------------------------------!
  
  type, extends(t_krige) :: t_simpkrige
  contains

  end type t_simpkrige
!-------------------------------------------------------------------------------------------------!
  
  type, extends(t_krige) :: t_ordkrige
  contains

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

  subroutine setup_sk(this, from_grids, to_grid,  cat, from_cats, to_cats, weights, opt)
    implicit none
    class(t_simpkrige)                    :: this
    type(t_grid),intent(in)               :: from_grids(:), to_grid
    real(dp)                              :: weights(:,:)
    integer,intent(in)                    :: cat
    integer,intent(in)                    :: from_cats(:), to_cats(:)
    type(t_options),intent(in)            :: opt
    
  end subroutine setup_sk
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine predict_point_sk(this, idx, p, from_grids, cat, from_cats, weights)
    implicit none
    class(t_simpkrige), intent(in)    :: this
    integer,intent(in)                :: idx, cat
    real(dp),intent(in)               :: p(2)
    type(t_grid),intent(in)           :: from_grids(:)
    integer,intent(in)                :: from_cats(:)
    real(dp)                          :: weights(:)
  
  end subroutine predict_point_sk
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! ORDINARY KRIGING CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

  subroutine setup_ok(this, from_grids, to_grid,  cat, from_cats, to_cats, weights, opt)
    implicit none
    class(t_simpkrige)                    :: this
    type(t_grid),intent(in)               :: from_grids(:), to_grid
    real(dp)                              :: weights(:,:)
    integer,intent(in)                    :: cat
    integer,intent(in)                    :: from_cats(:), to_cats(:)
    type(t_options),intent(in)            :: opt
    
  end subroutine setup_ok
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine predict_point_ok(this, idx, p, from_grids, cat, from_cats, weights)
    implicit none
    class(t_simpkrige), intent(in)    :: this
    integer,intent(in)                :: idx, cat
    real(dp),intent(in)               :: p(2)
    type(t_grid),intent(in)           :: from_grids(:)
    integer,intent(in)                :: from_cats(:)
    real(dp)                          :: weights(:)
  
  end subroutine predict_point_ok
  
!-------------------------------------------------------------------------------------------------!
  
end module m_kriging