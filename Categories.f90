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
      procedure, public        :: read_file
      procedure, private       :: initialize
      procedure, private       :: add_category => add_category_int, add_category_str
  end type
!-------------------------------------------------------------------------------------------------!
  contains
  
!-------------------------------------------------------------------------------------------------!
! CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  subroutine initialize(this, narr, nlayers, nnodes, datatype)
    ! Called by this%read_file()
    implicit none
    class(t_category)     :: this
    integer,intent(in)    :: narr, datatype, nlayers, nnodes
    
    this%dtype = datatype
    this%narr = narr
    allocate(this%lay2arr(nlayers      ), &
             this%array   (nnodes, narr), &  ! Assumes uniform size layers
             this%catlist(DEFAULT_CATLIST_SIZE)  )
    if (this%dtype==1) call vstrlist_new(this%strlist)
    
    ! By Default
    this%lay2arr =  1
    this%array   =  1
    this%catlist = -1
  
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
  
  subroutine read_file(this, file, nlayers, nnodes, datatype)
    use m_file_io, only: t_file_reader, open_file_reader, item2int, item2char
    use m_error_handler, only: error_handler
    implicit none
    class(t_category), intent(inout) :: this
    character(*), intent(in)         :: file
    integer,intent(in)               :: datatype, nlayers, nnodes
    type(t_file_reader), pointer     :: reader
    type(t_vstringlist)              :: strings
    integer                          :: lay2arr(nlayers)
    integer                          :: eof, narr
    
    ! Check if a file was actually passed. If not, we'll just initialize and be done
    if (len(trim(file))>0) then
      reader => open_file_reader(file)
      call reader%next_item(eof, strings)
      
      narr = maxval(lay2arr)
      call this%initialize(narr, nlayers, nnodes, datatype)
    else
      narr = 1
      call this%initialize(narr, nlayers, nnodes, datatype)
      return                   ! Early return if no file
    end if
    
    ! Continue if file passed
    
  end subroutine read_file
  
!-------------------------------------------------------------------------------------------------!
end module m_categories