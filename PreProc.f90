module m_preproc
  use m_global, only: NODATA
  use m_vstringlist
  use m_error_handler, only: error_handler,warning_handler,warnings
!-------------------------------------------------------------------------------------------------!
! Preprocessor module: convert non-texture data to a texture *after* interpolation to model grid
! Available preprocessors:
!  - RES2TEXTURE
!-------------------------------------------------------------------------------------------------!
  implicit none
  
  type :: t_preprocessor  ! abstract base class
    character(100)                   :: filename
    integer                          :: nclasses     ! Number of *preprocessor* classes (usually 1)
    integer                          :: ntexclasses  ! Texture datasets (in T2P)
    type(t_vstringlist)              :: class_names  ! *preprocessor* data class names
    
    contains
      procedure,private              :: initialize
    
  end type
  
  type,  extends(t_preprocessor) :: t_res2texture

  end type
!-------------------------------------------------------------------------------------------------!
  contains
!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  function create_preproc_object(preproc_mode_char, ntexclasses) result(object)
    implicit none

    character(*),intent(in)        :: preproc_mode_char
    integer, intent(in)            :: ntexclasses
    class(t_preprocessor),pointer  :: object

    select case(preproc_mode_char)
      case ("RES2TEXTURE")
        allocate(t_res2texture::object)
      case DEFAULT
        call error_handler(1,opt_msg="Unknown Preprocessor Mode Name: " // trim(preproc_mode_char))
    end select

    ! Initialize new object
    call object%initialize(ntexclasses)

  end function create_preproc_object
  
!-------------------------------------------------------------------------------------------------!
! BASE CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine initialize(this, ntexclasses)
    implicit none
    class(t_preprocessor), intent(inout) :: this
    integer, intent(in)                  :: ntexclasses
    
    call vstrlist_new(this%class_names)
    this%nclasses    = 0
    this%ntexclasses = ntexclasses
  
  end subroutine initialize
  
!-------------------------------------------------------------------------------------------------!
! RES2TEXTURE TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
end module m_preproc