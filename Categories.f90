module m_categories
  use m_vstring, only: t_vstring
  use m_vstringlist, only: t_vstringlist, vstrlist_new, vstrlist_append, vstrlist_search
  implicit none
  
  integer,parameter   :: DEFAULT_CATLIST_SIZE=10
!-------------------------------------------------------------------------------------------------!
  type :: t_category
    integer                    :: dtype          ! Data type, 0-int, 1-str
    integer                    :: narr           ! Number of unique arrays
    integer                    :: ncat           ! Number of categories
    integer,allocatable        :: catlist(:)     ! List of category values
    type(t_vstringlist)        :: strlist        ! Category values, but when string
    integer,allocatable        :: lay2arr(:)     ! Translation from layer to array
    integer,allocatable        :: array(:,:)     ! Arrays of category values (nvalues, narr)
    
    contains
      procedure, public        :: initialize
      procedure, public        :: add_category => add_category_int, add_category_str
  end type
!-------------------------------------------------------------------------------------------------!
  contains
!-------------------------------------------------------------------------------------------------!
! CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  subroutine initialize(this, narr, fmodel, datatype)
    use m_flow_model, only: t_flow_model
    implicit none
    class(t_category)     :: this
    integer,intent(in)    :: narr, datatype
    type(t_flow_model)    :: fmodel
    
    this%dtype = datatype
    this%narr = narr
    allocate(this%lay2arr(fmodel%nlayers      ), &
             this%array(fmodel%nnodes(1), narr), &  ! Assumes uniform size layers
             this%catlist(DEFAULT_CATLIST_SIZE)  )
    if (this%dtype==1) call vstrlist_new(this%strlist)
  
  end subroutine initialize

!-------------------------------------------------------------------------------------------------!

  subroutine add_category_int(this, value, index)
    class(t_category), intent(inout) :: this
    integer, intent(in)              :: value
    logical                          :: is_new = .true.
    integer, intent(out)             :: index

    ! Check if new_cat is a new unique category
    do index = 1, this%ncat
      if (this%catlist(index) == value) then
        is_new = .false.
        exit             ! Early exit
      end if
    end do

    if (is_new) then
      if (this%ncat >= size(this%catlist)) then
        ! Need to reallocate catlist to make room for new category
        call reallocate_catlist(this)
      end if
      this%ncat = this%ncat + 1
      index = this%ncat
      this%catlist(this%ncat) = value
    end if
  end subroutine add_category_int
  
!-------------------------------------------------------------------------------------------------!

  subroutine add_category_str(this, value, index)
    class(t_category), intent(inout) :: this
    type(t_vstring), intent(in)      :: value
    integer, intent(out)             :: index

    ! Check if new_cat is a new unique category
    index = vstrlist_search(this%strlist, value)
    if (index == 0) then
      call vstrlist_append(this%strlist, value)
      this%ncat = this%ncat + 1
      index = this%ncat
    end if
  end subroutine add_category_str
  
!-------------------------------------------------------------------------------------------------!

  subroutine reallocate_catlist(this)
    class(t_category), intent(inout) :: this
    integer, allocatable :: temp_catlist(:)
    integer :: new_size

    new_size = size(this%catlist) + 5  ! Increase size by 5, adjust as necessary
    allocate(temp_catlist(new_size))
    temp_catlist(:size(this%catlist)) = this%catlist
    call move_alloc(temp_catlist, this%catlist)
  end subroutine reallocate_catlist
  
!-------------------------------------------------------------------------------------------------!
end module m_categories