module m_file_io
  use str_module
  use m_vstring
  use m_vstringlist
  use list_module, only: stack, list_node
  implicit none

!-------------------------------------------------------------------------------------------------!
  type t_file_handler  ! Abstract Base Class
    
    character(100)           :: file
    integer                  :: unit
    integer                  :: ioerr
    
    contains
      procedure,private      :: iomsg_handler
      procedure,public       :: close_file
    
  end type t_file_handler
!-------------------------------------------------------------------------------------------------!
  type, extends(t_file_handler) :: t_file_reader
  
    type(stack)              :: backstack
  
    contains
      procedure,private      :: jump_to_external
      procedure,private      :: handle_eof
      procedure,private      :: read_to_next_line
      procedure,public       :: next_block_id
      procedure,public       :: next_block_item
      procedure,public       :: get_block_len
      procedure,public       :: back
    
  end type t_file_reader
!-------------------------------------------------------------------------------------------------!
    type, extends(t_file_handler) :: t_file_writer
  
    contains
      procedure,public      :: write_line
    
  end type t_file_writer

!-------------------------------------------------------------------------------------------------!

  integer                     :: max_unit   ! Maximum unit currently open (starts at 9)
                                            ! (notably, this is not the maximum units allowed open)
  integer                     :: log_unit   ! Unit the log file is open on
                                            ! (so IO errors can be written to it)
  character(500)              :: line       ! For reading in lines from files
  
  contains

!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine io_initialize()
    implicit none
    
    ! Sets the initial unit number used for opening a file
    
    max_unit = 9
    log_unit = 0
    
  end subroutine io_initialize

!-------------------------------------------------------------------------------------------------!
  
  integer function get_next_iunit() result(unit_no)
    implicit none
    logical      :: isopen
    isopen = .true.
    unit_no = max_unit
    do while (isopen==.true.)
      unit_no = unit_no + 1
      inquire(unit=unit_no, opened=isopen)
    end do
    return
    max_unit = unit_no
  end function get_next_iunit

!-------------------------------------------------------------------------------------------------!
  
  subroutine return_iunit(unit_no)  ! should re-name, maybe "de-register" ?
    implicit none
    integer, intent(in) :: unit_no
    
    ! TODO - more robust handling
    max_unit = unit_no - 1
    
  end subroutine return_iunit
  
!-------------------------------------------------------------------------------------------------!
  
  function open_file_reader(file) result(f)
    implicit none
    
    type(t_file_reader),pointer  :: f
    character(*), intent(in)   :: file
    character(256)             :: iomsg
    
    allocate(f)
    
    f%unit = get_next_iunit()
    f%file = trim(file)
    
    open(f%unit, file=f%file, action='READ', STATUS='OLD', iostat=f%ioerr, iomsg=iomsg)
    write(log_unit,'(3a)') 'Opening file', trim(f%file), ' for reading.'
    call f%iomsg_handler(iomsg)
    
    
  end function open_file_reader

!-------------------------------------------------------------------------------------------------!
  
  function open_file_writer(file) result(f)
    implicit none
    
    type(t_file_writer),pointer  :: f
    character(*), intent(in)   :: file
    character(256)             :: iomsg
    
    allocate(f)
    
    f%unit = get_next_iunit()
    f%file = trim(file)
    
    open(f%unit, file=f%file, action='WRITE', STATUS='REPLACE', iostat=f%ioerr, iomsg=iomsg)
    if (log_unit > 0) write(log_unit,'(3a)') 'Opening file', trim(f%file), ' for writing.'
    call f%iomsg_handler(iomsg)
    
  end function open_file_writer

!-------------------------------------------------------------------------------------------------!
  
    ! Should these be here or in a different module??
  function item2int(strings, i) result(out_int)
    implicit none
    class(t_vstringlist),intent(in)   :: strings
    integer, intent(in)               :: i
    integer                           :: out_int
    
    call vstring_cast(vstrlist_index(strings, i), out_int)
  end function item2int
  
!-------------------------------------------------------------------------------------------------!
  
    ! Should these be here or in a different module??
  function item2real(strings, i) result(out_real)
    implicit none
    class(t_vstringlist),intent(in)   :: strings
    integer, intent(in)               :: i
    real                              :: out_real
    
    call vstring_cast(vstrlist_index(strings, i), out_real)
  end function item2real
  
!-------------------------------------------------------------------------------------------------!
  
    ! Should these be here or in a different module??
  subroutine item2char(strings, i, out_char, toupper)
    implicit none
    class(t_vstringlist),intent(in)   :: strings
    integer, intent(in)               :: i
    logical, optional                 :: toupper
    character(*),intent(out)          :: out_char
    
    if (present(toupper)) then
      if (toupper) then
        call vstring_cast(vstring_toupper(vstrlist_index(strings, i)), out_char)
      end if
    else
      call vstring_cast(vstrlist_index(strings, i), out_char)
    end if
  end subroutine item2char
  
!-------------------------------------------------------------------------------------------------!
  
    ! Should these be here or in a different module??
  function find_string_index_in_list(strings, i, list, toupper) result(index)
    implicit none
    ! Given input (strings(i)), find the index of it in the comma delimited character list "list"
    ! Useful when you have an input (e.g., TYPE IWFM) and you want to know which valid entry in
    ! a character list ("GRID,IWFM,MODFLOW") it is (in this case, IWFM is 2). Then a variable
    ! can be set using this value. If the input is not in the character list, it returns zero.
    class(t_vstringlist),intent(in)   :: strings
    integer, intent(in)               :: i
    character(*), intent(in)          :: list
    logical, optional                 :: toupper
    type(t_vstring)                   :: temp_str
    integer                           :: index
    
    call vstring_new(temp_str, list)
    if (present(toupper)) then
      if (toupper) then
        index = vstrlist_search(vstrlist_split(temp_str, ","), vstring_toupper(vstrlist_index(strings, i)) )
      end if
    else
      index = vstrlist_search(vstrlist_split(temp_str, ","), vstrlist_index(strings, i) )
    end if
    
  end function find_string_index_in_list
  
!-------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------!
! BASE CLASS FILE TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine iomsg_handler(this,iomsg)
    implicit none
    class(t_file_handler)        :: this
    character(256),intent(in)    :: iomsg
    
    if (this%ioerr /= 0) then
      write(*,'(a)') trim(iomsg)
      write(log_unit,'(3a)') 'IO Error -', trim(iomsg)
      close(this%unit)
      error stop
    end if
    
  end subroutine iomsg_handler

!-------------------------------------------------------------------------------------------------!

  subroutine close_file(this)
    implicit none
    class(t_file_handler)        :: this
    close(this%unit)
    
    call return_iunit(this%unit)
    
    !deallocate(f)
    !nullify(f)
  end subroutine close_file
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! READER TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine read_to_next_line(this, eof, skipped)
    implicit none
    class(t_file_reader)              :: this
    type(t_file_reader),pointer       :: next_file
    integer                           :: i, ierr
    integer,intent(out)               :: eof
    integer,intent(inout),optional    :: skipped
    
    ! For line splitting
    integer             :: left(2), right(2)

    eof = 0
    i = 0
    do while (eof == 0)
      read(this%unit, '(a500)', iostat=ierr) line
      line = adjustl(line)
      ! Handle EOF
      if (IS_IOSTAT_END(ierr)) then
        call this%handle_eof(eof)
      else if ((line(1:1) == '#').or.(line(1:1) == '*').or.(line(1:1) == '')) then
        ! Comment line, blank line
        i = i + 1
        continue
      else if (index(line, 'EXTERNAL')>0) then
        ! External command, continue reading in new file
        call multisplit(ierr,2,left,right,line)
        call this%jump_to_external(line(left(2):right(2)))
        i = 0  ! No lines skipped if it's a new file
        !call this%read_to_next_line(eof)
      else
        ! Data
        backspace(this%unit)
        exit
      end if
    end do
    
    if (present(skipped)) skipped = i
  
  end subroutine
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine jump_to_external(this, file)
    implicit none
    class(t_file_reader)  :: this
    character(*), intent(in) :: file
    character(256)           :: iomsg
    
    ! Store old iunit
    call this%backstack%push(this%unit)
    
    ! Update where we're reading from
    this%unit = get_next_iunit()
  
    open(this%unit, file=file, action='READ', STATUS='OLD', iostat=this%ioerr)
    write(log_unit,'(3a)') 'Opening external file', trim(file), ' for reading.'
    call this%iomsg_handler(iomsg)
    
  end subroutine jump_to_external
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine handle_eof(this, eof)
    implicit none
    class(t_file_reader)     :: this
    type(list_node)        :: prev
    integer, intent(inout) :: eof
    
    if (this%backstack%size > 0) then
      close(this%unit)
      prev = this%backstack%pop()
      this%unit  = prev%value
      eof = 0
    else
      eof = -1
    end if
    
  end subroutine handle_eof
  
!-------------------------------------------------------------------------------------------------!
    subroutine next_item(this, eof, strings)
    implicit none
    class(t_file_reader)              :: this
    integer,intent(out)               :: eof
    type(t_vstringlist), intent(out)  :: strings
    type(t_vstring)                   :: line_string
    integer                           :: ierr

    ! Get to the next line
    call this%read_to_next_line(eof)
    if (eof == 0) then
    
      ! Read string, ensure block -> pass to block reader
      read(this%unit, '(a500)', iostat=ierr) line
      ! TODO Error handling
      call vstring_new(line_string, trim(adjustl(line)))
      strings = vstrlist_split(line_string)
    
    else
      ! End of file - just return
    end if
  
  end subroutine next_item
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine next_block_id(this, eof, id)  
    implicit none
    class(t_file_reader)        :: this
    integer,intent(out)       :: eof
    character(*),intent(out)  :: id
    type(t_vstringlist)       :: slist

    ! Get next item in file
    call next_item(this, eof, slist)
    if (eof == 0) then
      ! Make sure there's enough there...
      if (vstrlist_length(slist) < 2) then
        ! ERROR
        stop
      end if
      !! Get block name
      if (vstring_equals(vstring_toupper(vstrlist_index(slist, 1)), "BEGIN")) then
        call vstring_cast(vstring_toupper(vstrlist_index(slist, 2)), id)
      end if
    else
      ! End of file - just return
    end if
  
  end subroutine next_block_id

!-------------------------------------------------------------------------------------------------!
  
  subroutine next_block_item(this, status, id, strings, length)
    implicit none
    class(t_file_reader)      :: this
    integer,intent(out)       :: status   ! 0 = Fine. -1 = EOF. +1 = End of Block
    integer,intent(out)       :: length
    character(*),intent(out)  :: id
    type(t_vstringlist)       :: strings

    ! Get next item in file
    call next_item(this, status, strings)
    if (status == 0) then
      
      !! Check for end
      call vstring_cast(vstring_toupper(vstrlist_index(strings, 1)), id)
      if (trim(id)=="END") then
        status = 1
        length = 1
      else
        length = vstrlist_length(strings)
      end if
    else
      status = -1
    end if
  
  end subroutine next_block_item

!-------------------------------------------------------------------------------------------------!
  
  function get_block_len(this) result(block_len)
    implicit none
    class(t_file_reader), intent(in) :: this
    integer                         :: i, block_len, eof, ierr
    integer                         :: status, skipped, total_skipped, length
    character(30)                   :: id
    type(t_vstringlist)             :: strings
    
    ! TODO what if file changes in block? (backspace may fail)

    status = 0
    total_skipped = 0
    block_len = -1  ! To account for reading END line
          
    ! Loop through the block items and count non-comment lines
    call this%read_to_next_line(eof, skipped=total_skipped)
    !read(this%unit, '(a500)', iostat=ierr) line
    if (index(line, "END") > 0.or.eof < 0) status = 1
    do while (status == 0)
      block_len = block_len + 1
      call this%read_to_next_line(eof, skipped=skipped)
      read(this%unit, '(a500)', iostat=ierr) line
      total_skipped = total_skipped + skipped + 1
      if (index(line, "END") > 0.or.eof < 0) status = 1
    end do
    
    ! Rewind!
    do i=1, total_skipped
      backspace(this%unit)
    end do
          
  end function get_block_len

!-------------------------------------------------------------------------------------------------!
  
  subroutine back(this, nlines)
    implicit none
    class(t_file_reader), intent(in) :: this
    integer,intent(in)               :: nlines
    integer                          :: i
    
    do i=1, nlines
      backspace(this%unit)
    end do
    
  end subroutine back
  
!-------------------------------------------------------------------------------------------------!
! WRITER TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine write_line(this, outline, blank_start, blank_end)
    implicit none
    class(t_file_writer), intent(inout) :: this
    character(len=*), intent(in)        :: outline
    integer,intent(in),optional         :: blank_start, blank_end
    integer                             :: i
    
    if (present(blank_start)) then
      do i=1, blank_start
        write(this%unit,*)
      end do
    end if
    
    write(this%unit, '(A)') outline
    
    if (present(blank_end)) then
      do i=1, blank_end
        write(this%unit,*)
      end do
    end if
    
  end subroutine write_line

!-------------------------------------------------------------------------------------------------!

end module m_file_io