module m_categories
  use m_vstring, only: t_vstring
  use m_vstringlist
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
    logical                    :: from_file      ! Categories were read from file
    
    contains
      procedure, public        :: read_file
      procedure, private       :: initialize
      procedure, private       :: no_category_setup
      procedure, pass          :: add_category_int
      procedure, pass          :: add_category_str
      procedure, pass          :: add_category_char
      generic, private         :: add_category => add_category_int,add_category_str,add_category_char
      procedure, pass          :: get_id_char
      generic, public          :: get_id => get_id_char
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
    this%ncat       =  0
    this%lay2arr    =  1
    this%array      =  1
    this%catlist    = -1
  
  end subroutine initialize

!-------------------------------------------------------------------------------------------------!
  
  subroutine no_category_setup(this)
    class(t_category), intent(inout) :: this
    
    this%ncat = 1
    this%catlist(1) = 1
    if (this%dtype==1) then
      call vstrlist_append(this%strlist, '')  ! Necessary?
    end if
  
  end subroutine no_category_setup
  
!-------------------------------------------------------------------------------------------------!

  subroutine add_category_int(this, value, index)
    class(t_category), intent(inout) :: this
    integer, intent(in)              :: value
    logical                          :: is_new
    integer, intent(out)             :: index

    ! Check if new_cat is a new unique category
    is_new = .true.
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
      this%catlist(this%ncat) = this%ncat
      index = this%ncat
    end if
  end subroutine add_category_str
  
!-------------------------------------------------------------------------------------------------!

  subroutine add_category_char(this, value, index)
    class(t_category), intent(inout) :: this
    character(*), intent(in)         :: value
    integer, intent(out)             :: index

    ! Check if new_cat is a new unique category
    index = vstrlist_search(this%strlist, trim(value))
    if (index == 0) then
      call vstrlist_append(this%strlist, trim(value))
      this%ncat = this%ncat + 1
      this%catlist(this%ncat) = this%ncat
      index = this%ncat
    end if
  end subroutine add_category_char
  
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
    !type(t_vstring)                  :: temp
    integer                          :: lay2arr(nlayers)
    integer                          :: eof, narr, cidx, i, j, inode
    character(30)                    :: chartemp(nnodes)
    
    lay2arr = 0
    
    ! Check if a file was actually passed. If not, we'll just initialize and be done
    if (len(trim(file))>0) then
      this%from_file = .true.
      reader => open_file_reader(file)
      call reader%next_item(eof, strings)
      
      if (vstrlist_length(strings) < nlayers) then
        call error_handler(2,filename=file,opt_msg="Not enough array indices assigned: need one index for each layer")
      end if
      
      ! Move to temp lay2arr
      do i=1, nlayers
        lay2arr(i) = item2int(strings, i)
      end do
      
      narr = maxval(lay2arr)
      call this%initialize(narr, nlayers, nnodes, datatype)
      ! Move to new lay2array
      this%lay2arr = lay2arr
    else
      this%from_file = .false.
      narr = 1
      call this%initialize(narr, nlayers, nnodes, datatype)
      call this%no_category_setup()
      return                   ! Early return if no file
    end if
    
    ! Continue if file passed - Loop over arrays, skipping header
    do i=1, narr
      call reader%read_to_next_line(eof)  ! Skip header
      if (datatype==0) then ! INT
        read(reader%unit,*) this%array(:,i)
        do j=1, nnodes
          call this%add_category(this%array(j,i), cidx)
          !this%array(j,i) = cidx  ! Excessively confusing, e.g. cat 4 is 3...
        end do
      else if (datatype==1) then  ! STR/CHAR
        read(reader%unit,*) chartemp
        do j=1, nnodes
          call this%add_category(chartemp(j), cidx)
          this%array(j,i) = cidx
        end do
        ! Too slow with vstrlists
        !inode = 1
        !do while (inode < nnodes)
        !  call reader%next_item(eof, strings)
        !  do j=1, vstrlist_length(strings)
        !    temp = vstrlist_index(strings, j)
        !    call this%add_category(temp, cidx)
        !    this%array(inode,i) = cidx
        !    inode = inode + 1
        !  end do
        !end do
      end if
    end do
    
  end subroutine read_file
  
!-------------------------------------------------------------------------------------------------!
  
  function get_id_char(this, value) result(index)
    ! Gets strlist index of pass character category
    implicit none
    class(t_category), intent(inout) :: this
    character(*),intent(in)          :: value
    integer                          :: index
    
    index = vstrlist_search(this%strlist, trim(value))
  
  end function get_id_char
  
!-------------------------------------------------------------------------------------------------!
end module m_categories