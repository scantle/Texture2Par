module m_weights
!-------------------------------------------------------------------------------------------------!  
  implicit none

  ! Define a custom type for managing weights
  type :: t_weights
      integer                :: nin, nout     ! grid input point count, grid output point count
      real, allocatable      :: weights(:,:)  ! The weights (input points, output points)
      integer, allocatable   :: ids(:,:)      ! index of weights within their dataset
      integer, allocatable   :: group_id(:)   ! Group ID for each weight
      ! Additional metadata can be added here as needed
  contains
      procedure,private :: initialize
      procedure,public :: save_to_file
      procedure,public :: load_from_file
      procedure,private :: finalize
  end type t_weights
!-------------------------------------------------------------------------------------------------!  
contains

!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!

function create_weight_object(input_n, output_n) result(object)
  implicit none
  integer,intent(inout)    :: input_n, output_n
  type(t_weights), pointer :: object
  
  allocate(object)
  call object%initialize(input_n,output_n)

end function create_weight_object

!-------------------------------------------------------------------------------------------------!
! WEIGHT CLASS PROCEDURES
!-------------------------------------------------------------------------------------------------!  

  ! Initialize the weights object
  subroutine initialize(self, input_n, output_n)
    implicit none
    integer,intent(in)              :: input_n, output_n
    class(t_weights), intent(inout) :: self
    
    self%nin  = input_n
    self%nout = output_n
    allocate(self%weights (input_n, output_n), &
             self%ids     (input_n, output_n), &
             self%group_id(output_n)           )
      
  end subroutine initialize
  
!-------------------------------------------------------------------------------------------------!
  
  ! Save weights to a file
  subroutine save_to_file(self, filename)
    implicit none
    class(t_weights), intent(in) :: self
    character(len=*), intent(in) :: filename

      ! Implementation to serialize and save weights to a file
  end subroutine save_to_file

!-------------------------------------------------------------------------------------------------!
  
  ! Load weights from a file
  subroutine load_from_file(self, filename)
    implicit none
    class(t_weights), intent(inout) :: self
    character(len=*), intent(in) :: filename

      ! Implementation to deserialize and load weights from a file
  end subroutine load_from_file

!-------------------------------------------------------------------------------------------------!
  subroutine finalize(self)
    class(t_weights), intent(inout) :: self
    if (allocated(self%weights)) deallocate(self%weights)
    if (allocated(self%ids)) deallocate(self%ids)
    if (allocated(self%group_id)) deallocate(self%group_id)
  end subroutine finalize
!-------------------------------------------------------------------------------------------------!
  
end module m_weights
