module m_weights
!-------------------------------------------------------------------------------------------------!  
  implicit none

  ! Define a custom type for managing weights
  type :: t_weights
      integer                :: nin, nout     ! grid input point count, grid output point count
      real, allocatable      :: values(:,:)  ! The weights (output_points, input_points)
      
  contains
      procedure,public  :: initialize
      procedure,public  :: save_to_file
      procedure,public  :: load_from_file
      procedure,private :: finalize
  end type t_weights
!-------------------------------------------------------------------------------------------------!  
! Alternatively: single-output sparse (SOSA)
  
  type :: t_sosa_weights
    integer             :: nin
    integer             :: unit
    integer,allocatable :: idx(:)
    real,allocatable    :: values(:)
  contains
      procedure,public  :: initialize => initialize_sosa
      procedure,public  :: open_file => open_file_sosa
      procedure,public  :: save_to_file => save_to_file_sosa
      procedure,public  :: load_from_file => load_from_file_sosa
      procedure,private :: finalize => finalize_sosa
  end type t_sosa_weights
  
contains

!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------!
! WEIGHT CLASS PROCEDURES
!-------------------------------------------------------------------------------------------------!  

  ! Initialize the weights object
  subroutine initialize(this, input_n, output_n)
    implicit none
    class(t_weights), intent(inout) :: this
    integer,intent(in)              :: input_n, output_n
    
    this%nin  = input_n
    this%nout = output_n
    if (allocated(this%values)) deallocate(this%values)
    allocate(this%values (output_n, input_n))
      
  end subroutine initialize
  
!-------------------------------------------------------------------------------------------------!
  
  ! Save weights to a file
  subroutine save_to_file(this, filename)
    implicit none
    class(t_weights), intent(in) :: this
    character(len=*), intent(in) :: filename

      ! Implementation to serialize and save weights to a file
  end subroutine save_to_file

!-------------------------------------------------------------------------------------------------!
  
  ! Load weights from a file
  subroutine load_from_file(this, filename)
    implicit none
    class(t_weights), intent(inout) :: this
    character(len=*), intent(in) :: filename

      ! Implementation to deserialize and load weights from a file
  end subroutine load_from_file

!-------------------------------------------------------------------------------------------------!
  
  subroutine finalize(this)
    class(t_weights), intent(inout) :: this
    if (allocated(this%values)) deallocate(this%values)
  end subroutine finalize
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! SOSA CLASS PROCEDURES
!-------------------------------------------------------------------------------------------------!  
  
  ! Initialize the weights object
  subroutine initialize_sosa(this, input_n, output_n)
    implicit none
    class(t_sosa_weights), intent(inout) :: this
    integer,intent(in)              :: input_n, output_n
    
    this%nin  = input_n
    if (allocated(this%idx)) deallocate(this%idx)
    if (allocated(this%values)) deallocate(this%values)
    allocate(this%idx(input_n), this%values(input_n))
      
  end subroutine initialize_sosa
  
!-------------------------------------------------------------------------------------------------!

    ! Save weights to a file
  subroutine open_file_sosa(this, filename)
    implicit none
    class(t_sosa_weights), intent(in) :: this
    character(len=*), intent(in) :: filename

      ! Implementation to serialize and save weights to a file
  end subroutine open_file_sosa

  !-------------------------------------------------------------------------------------------------!
  
  ! Save weights to a file
  subroutine save_to_file_sosa(this, filename)
    implicit none
    class(t_sosa_weights), intent(in) :: this
    character(len=*), intent(in) :: filename

      ! Implementation to serialize and save weights to a file
  end subroutine save_to_file_sosa

!-------------------------------------------------------------------------------------------------!
  
  ! Load weights from a file
  subroutine load_from_file_sosa(this, filename)
    implicit none
    class(t_sosa_weights), intent(inout) :: this
    character(len=*), intent(in) :: filename

      ! Implementation to deserialize and load weights from a file
  end subroutine load_from_file_sosa

!-------------------------------------------------------------------------------------------------!
  
  subroutine finalize_sosa(this)
    class(t_sosa_weights), intent(inout) :: this
    if (allocated(this%values)) deallocate(this%values)
    if (allocated(this%values)) deallocate(this%values)
  end subroutine finalize_sosa

!-------------------------------------------------------------------------------------------------!
  
end module m_weights
