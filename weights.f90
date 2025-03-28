module m_weights
  use m_file_io !, only: open_file_writer, open_file_writer, t_file_reader, t_file_writer
  use m_error_handler, only: error_handler
  
!-------------------------------------------------------------------------------------------------!  
  implicit none

  ! Define a custom type for managing weights
  type :: t_weights
      integer                :: nin, nout     ! grid input point count, grid output point count
      real    , allocatable  :: values(:,:)  ! The weights (output_points, input_points)
      
  contains
      procedure,public   :: initialize
      procedure,public   :: save_to_file
      procedure,public   :: load_from_file
      procedure,private  :: finalize
  end type t_weights
!-------------------------------------------------------------------------------------------------!  
! Alternatively: single-output sparse (SOSA)
  
  type :: t_sosa_weights
    integer                        :: nin
    class(t_file_handler),pointer  :: f => null()
    integer,allocatable            :: ppd(:)
    integer,allocatable            :: idx(:)
    real    ,allocatable           :: values(:)
  contains
      procedure,public   :: initialize => initialize_sosa
      procedure,public   :: open_file => open_file_sosa
      procedure,public   :: close_file => close_file_sosa
      procedure,public   :: save_to_file => save_to_file_sosa
      procedure,public   :: load_from_file => load_from_file_sosa
      procedure,private  :: finalize => finalize_sosa
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
    this%values = 0.0
      
  end subroutine initialize
  
!-------------------------------------------------------------------------------------------------!
  
  ! Save weights to a file
  subroutine save_to_file(this, filename)
    implicit none
    class(t_weights), intent(in) :: this
    character(len=*), intent(in) :: filename
    type(t_file_writer),pointer   :: f

    f => open_file_writer(filename, binary=.true.)
    write(f%unit) this%nin, this%nout
    write(f%unit) this%values(:,:)
    call f%close_file()
    
  end subroutine save_to_file

!-------------------------------------------------------------------------------------------------!
  
  ! Load weights from a file
  subroutine load_from_file(this, filename)
    implicit none
    class(t_weights), intent(inout) :: this
    character(len=*), intent(in)    :: filename
    type(t_file_reader),pointer     :: f
    integer                         :: input_n, output_n
    character(256)                  :: iomsg

    f => open_file_reader(filename, binary=.true.)
    read(f%unit) input_n, output_n
    ! Check input
    if ((input_n /= this%nin).or.(output_n /= this%nout)) then
      call error_handler(1,filename,opt_msg="Weight file dimension mismatch")
    end if
    read(f%unit, iomsg=iomsg) this%values(:,:)
    call f%iomsg_handler(iomsg)  ! in case something goes wrong with reading the values
    call f%close_file()
      
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
  subroutine initialize_sosa(this, input_n, classes_n)
    implicit none
    class(t_sosa_weights), intent(inout) :: this
    integer,intent(in)              :: input_n, classes_n
    
    this%nin  = input_n
    if (allocated(this%idx)) deallocate(this%idx)
    if (allocated(this%values)) deallocate(this%values)
    if (allocated(this%ppd)) deallocate(this%ppd)
    allocate(this%idx(input_n), this%values(input_n), this%ppd(classes_n))
      
  end subroutine initialize_sosa
  
!-------------------------------------------------------------------------------------------------!

  subroutine open_file_sosa(this, filename, writing, reading)
    implicit none
    class(t_sosa_weights), intent(inout) :: this
    character(len=*), intent(in)         :: filename
    logical,intent(in)                   :: writing, reading
    integer                              :: arrdim(2)
    
    if (reading) then
      allocate(t_file_reader::this%f)
      this%f => open_file_reader(filename, .true.)
      ! Check dimensions
      read(this%f%unit) arrdim
      if ((arrdim(1) /= size(this%ppd)).or.(arrdim(2) /= this%nin)) then
        call error_handler(1,filename,opt_msg="Weight file dimension mismatch")
      end if
    else if (writing) then
      allocate(t_file_writer::this%f)
      this%f => open_file_writer(filename, .true.)
      ! Write some dimensions to be checked by read
      write(this%f%unit) size(this%ppd), this%nin
    end if
    
  end subroutine open_file_sosa

  !-------------------------------------------------------------------------------------------------!
  
  subroutine close_file_sosa(this)
    implicit none
    class(t_sosa_weights), intent(inout) :: this

    call this%f%close_file()
    
  end subroutine close_file_sosa

  !-------------------------------------------------------------------------------------------------!
  
  ! Save weights to a file - one grid point at a time
  subroutine save_to_file_sosa(this, grid_idx)
    implicit none
    class(t_sosa_weights), intent(in) :: this
    integer,intent(in)                :: grid_idx

    write(this%f%unit) grid_idx, this%ppd(1:), this%idx(1:), this%values(1:)
    
  end subroutine save_to_file_sosa

!-------------------------------------------------------------------------------------------------!
  
  ! Load weights from a file
  subroutine load_from_file_sosa(this, grid_idx)
    implicit none
    class(t_sosa_weights), intent(inout) :: this
    integer,intent(inout)                :: grid_idx

    ! Read it in...
    read(this%f%unit) grid_idx, this%ppd(1:), this%idx(1:), this%values(1:)
    
  end subroutine load_from_file_sosa

!-------------------------------------------------------------------------------------------------!
  
  subroutine finalize_sosa(this)
    class(t_sosa_weights), intent(inout) :: this
    if (allocated(this%idx))    deallocate(this%idx)
    if (allocated(this%values)) deallocate(this%values)
    if (allocated(this%ppd))    deallocate(this%ppd)
  end subroutine finalize_sosa

!-------------------------------------------------------------------------------------------------!
  
end module m_weights
