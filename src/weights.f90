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
  subroutine save_to_file(this, filename, text_file)
    implicit none
    class(t_weights), intent(in) :: this
    character(len=*), intent(in) :: filename
    logical,intent(in)           :: text_file
    type(t_file_writer),pointer  :: f

    f => open_file_writer(filename, binary=.not.text_file)
    if (f%binary) then
      write(f%unit) this%nin, this%nout
      write(f%unit) this%values(:,:)
    else
      write(f%unit,*) this%nin, this%nout
      write(f%unit,*) this%values(:,:)
    end if
    call f%close_file()
    
  end subroutine save_to_file

!-------------------------------------------------------------------------------------------------!
  
  ! Load weights from a file
  subroutine load_from_file(this, filename, text_file)
    implicit none
    class(t_weights), intent(inout) :: this
    character(len=*), intent(in)    :: filename
    logical,intent(in)              :: text_file
    type(t_file_reader),pointer     :: f
    integer                         :: input_n, output_n
    character(256)                  :: iomsg

    f => open_file_reader(filename, binary=.not.text_file)
    if (f%binary) then
      read(f%unit) input_n, output_n
    else
      read(f%unit,*) input_n, output_n
    end if
    ! Check input
    if ((input_n /= this%nin).or.(output_n /= this%nout)) then
      call error_handler(1,filename,opt_msg="Weight file dimension mismatch")
    end if
    if (f%binary) then
      read(f%unit, iomsg=iomsg) this%values(:,:)
    else
      read(f%unit,*, iomsg=iomsg) this%values(:,:)
    end if
    call f%iomsg_handler(iomsg)  ! in case something goes wrong with reading the values
    call f%close_file()
      
  end subroutine load_from_file

!-------------------------------------------------------------------------------------------------!
  
  subroutine finalize(this)
    class(t_weights), intent(inout) :: this
    if (allocated(this%values)) deallocate(this%values)
  end subroutine finalize
  
!-------------------------------------------------------------------------------------------------!
  
end module m_weights
