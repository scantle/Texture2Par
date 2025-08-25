module m_preproc
  use m_global, only: NODATA
  use m_vstringlist
  use m_vstring, only: t_vstring, vstring_toupper
  use m_error_handler, only: error_handler,warning_handler,warnings
!-------------------------------------------------------------------------------------------------!
! Preprocessor module: convert non-texture data to a texture *after* interpolation to model grid
! Available preprocessors:
!  - RES2TEXTURE
!-------------------------------------------------------------------------------------------------!
  implicit none
  
  type :: t_preprocessor  ! abstract base class
    character(100)                   :: filename, mult_filename
    integer                          :: nclasses     ! Number of *preprocessor* classes (usually 1)
    integer                          :: ntexclasses  ! Texture datasets (in T2P)
    type(t_vstringlist)              :: class_names  ! *preprocessor* data class names
    
    contains
      procedure,private              :: initialize
      procedure,public               :: read_input_files
      procedure,public               :: to_texture
      procedure,public               :: get_class_id
      procedure,public               :: get_dataset_idx_list
      procedure,public               :: get_dataval_idx_list
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
  
  subroutine read_input_files(this)
    implicit none
    class(t_preprocessor), intent(inout) :: this
    
  end subroutine read_input_files
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine to_texture(this)
    implicit none
    class(t_preprocessor), intent(inout) :: this
    
  end subroutine to_texture
  
!-------------------------------------------------------------------------------------------------!
! T2P class replicating functions  
!-------------------------------------------------------------------------------------------------!
  
  function get_class_id(this, string, adjust_for_infer) result (id)
    ! Returns zero if not found
    implicit none
    class(t_preprocessor), intent(inout)   :: this
    type(t_vstring),intent(in)    :: string
    integer                       :: id
    logical,optional,intent(in)   :: adjust_for_infer  ! not used, just to keep compatability with t2p%get_class_id
    id = vstrlist_search(this%class_names, vstring_toupper(string), exact=.true.)
    return
  end function get_class_id

!-------------------------------------------------------------------------------------------------!
  
  function get_dataset_idx_list(this,class_ids) result(dataset_idx)
    implicit none
    class(t_preprocessor)   :: this
    integer,intent(in)      :: class_ids(:)
    integer,allocatable     :: dataset_idx(:)

    if (allocated(dataset_idx)) deallocate(dataset_idx)
    allocate(dataset_idx(size(class_ids)))
    dataset_idx(:) = 1    ! Will need to revise if we allow multiple preproc classes
    return
  end function

!-------------------------------------------------------------------------------------------------!

  function get_dataval_idx_list(this, class_ids) result (dataval_idx)
    implicit none
    class(t_preprocessor)   :: this
    integer,intent(in)      :: class_ids(:)
    integer,allocatable     :: dataval_idx(:)
    integer                 :: nactive

    if (allocated(dataval_idx)) deallocate(dataval_idx)
    nactive = 1    ! Will need to revise if we allow multiple preproc classes/cokriging
    dataval_idx = class_ids
    where (dataval_idx>nactive) dataval_idx = 1  ! All covariates use idx=1 for values

  end function get_dataval_idx_list
  
!-------------------------------------------------------------------------------------------------!
! RES2TEXTURE TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
end module m_preproc