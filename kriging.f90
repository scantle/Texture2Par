module m_kriging
  use m_global, only: dp
  use m_interpolator
  use m_weights, only: t_weights, create_weight_object
  use m_grid, only: t_grid, create_grid_object
  use m_vario, only: t_vgm
  use m_datasets, only: t_dataset, t_layerdata
  use m_error_handler, only: error_handler
  use m_options, only: t_options
  use m_sparse_symmetric, only: t_sparse_sym
  implicit none
!-------------------------------------------------------------------------------------------------!

  type, extends(t_interpolator) :: t_krige
    integer                          :: nvario            ! Number of variograms stored for interpolator
    integer,pointer                  :: nnear(:)          ! Number of nearest points used for interpolation
    integer,pointer                  :: nstruct(:)        ! Number of structures for each variogram
    type(t_vgm),pointer              :: variograms(:,:)   ! (stuctures, variogram)
    integer,allocatable              :: vario_idx(:,:)    ! Variogram index given (class_id, class_id)
    integer                          :: nsets             ! Number of datasets involved
    real(dp),allocatable             :: w(:)              ! Internally calculated weights
    ! Variables from PPSGS (slightly adjusted)
    integer                          :: ndim, ngrid, ndrift, seed, unbias, nsim
    integer,allocatable              :: nobs(:), nmax(:), near_idx(:)
    real(dp),allocatable             :: matA(:,:), rhsB(:)
    ! Attributes for KD-tree and search parameters could be added here
  contains
    procedure               :: intialize => initialize_kriging
    procedure               :: setup => setup_krige

  end type t_krige
!-------------------------------------------------------------------------------------------------!
  
  type, extends(t_krige) :: t_simpkrige
  contains
    procedure               :: setup => setup_sk
    procedure               :: predict_point => predict_point_sk
  end type t_simpkrige
!-------------------------------------------------------------------------------------------------!
  
  type, extends(t_krige) :: t_ordkrige
    integer                 :: maxnear
    integer,allocatable     :: usenear(:)
  contains
    procedure               :: setup => setup_ok
    procedure               :: predict_point => predict_point_ok
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
  
  subroutine setup_krige(this, obs_data, to_grid, cat, to_cats, weights, opt)
    ! Does everything common to all kriging
    ! Intended to be called by subclasses during their setup
    implicit none
    class(t_krige)              :: this
    class(t_dataset)            :: obs_data(:)
    type(t_grid),intent(in)     :: to_grid
    real(dp)                    :: weights(:,:)
    integer,intent(in)          :: cat
    integer,intent(in)          :: to_cats(:)
    type(t_options),intent(in)  :: opt
    integer                     :: i
    
    ! Defaults
    this%ndrift = 0  ! Hardcoded out for now
    this%unbias = 0 
    this%nsets  = size(obs_data)
    this%ngrid  = to_grid%n
    
    !call random_seed_initialize(seed)
    !call set_samples()
    
    !if (nsim>0) then
    !  call set_randpath()
    !  grid = grid(:, irandpath)
    !else
    !  irandpath=[(ig, ig=1, ngrid)]
    !end if
    
    
  end subroutine setup_krige
  
!-------------------------------------------------------------------------------------------------!
    
!-------------------------------------------------------------------------------------------------!
! SIMPLE KRIGING CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

  subroutine setup_sk(this, obs_data, to_grid, cat, to_cats, weights, opt)
    implicit none
    class(t_simpkrige)          :: this
    class(t_dataset)            :: obs_data(:)
    type(t_grid),intent(in)     :: to_grid
    real(dp)                    :: weights(:,:)
    integer,intent(in)          :: cat
    integer,intent(in)          :: to_cats(:)
    type(t_options),intent(in)  :: opt
    ! Local
    integer                     :: i,matsize
    
    ! Call super setup
    call setup_krige(this, obs_data, to_grid, cat, to_cats, weights, opt)
    
    matsize = 0
    do i=1, this%nsets
      matsize = obs_data(i)%grid%n
    end do
    
    matsize = to_grid%n+matsize+this%unbias+this%ndrift
    allocate(this%w(matsize))
    allocate(this%matA(matsize, matsize), this%rhsB(matsize))
    
    !pre-calculate distances
    !call setdist()
    
  end subroutine setup_sk
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine predict_point_sk(this, idx, p, obs_data, cat, weights)
    implicit none
    class(t_simpkrige), intent(in)    :: this
    integer,intent(in)                :: idx, cat
    real(dp),intent(in)               :: p(2)
    type(t_dataset)                   :: obs_data(:)
    real(dp)                          :: weights(:)
  
  end subroutine predict_point_sk
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! ORDINARY KRIGING CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

  subroutine setup_ok(this, obs_data, to_grid, cat, to_cats, weights, opt)
    implicit none
    class(t_ordkrige)          :: this
    class(t_dataset)            :: obs_data(:)
    type(t_grid),intent(in)     :: to_grid
    real(dp)                    :: weights(:,:)
    integer,intent(in)          :: cat
    integer,intent(in)          :: to_cats(:)
    type(t_options),intent(in)  :: opt
    ! Local
    integer                     :: i,matsize,vidx
    
    ! Call super setup
    call setup_krige(this, obs_data, to_grid, cat, to_cats, weights, opt)
    
    ! Setup nnearest tracker (usenear) (sometimes there is not this%nnear(i) values in the dataset)
    allocate(this%usenear(this%nsets))
    
    ! Build search trees in observed grid objects
    ! Also, figure out number of datasets in each category
    do i=1, this%nsets
      vidx = this%vario_idx(i,i)              ! Should always equal i, but to be safe
      call obs_data(i)%build_tree_by_category(cat)   ! rotmat passed during n_nearest call
      this%usenear(i) = min(this%nnear(i), obs_data(i)%grid%tree%n)
      ! TODO Warn that there is so little data? Do we even have enough information to warn?
      ! Warnings only get printed before the calculations
      ! In this scope we don't if it's pp or data...
    end do
    
    this%maxnear = maxval(this%nnear)
    
    matsize = sum(this%nnear)+this%unbias+this%ndrift
    allocate(this%w(matsize))
    allocate(this%matA(matsize, matsize), this%rhsB(matsize))
    
  end subroutine setup_ok
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine predict_point_ok(this, idx, p, obs_data, cat, weights)
    implicit none
    class(t_ordkrige), intent(in)     :: this
    integer,intent(in)                :: idx, cat
    real(dp),intent(in)               :: p(2)
    type(t_dataset)                   :: obs_data(:)
    real(dp)                          :: weights(:)
    ! local
    integer                           :: i
    integer                           :: inear(this%maxnear,this%nsets)
    real(dp)                          :: dist(this%maxnear,this%nsets)
    
    ! Get nearest for each dataset
    do i=1, this%nsets
      call obs_data(i)%grid%get_nnear(p, this%usenear(i), inear(:,i), dist(:,i))
    end do
    
    
  
  end subroutine predict_point_ok
  
!-------------------------------------------------------------------------------------------------!
  
end module m_kriging