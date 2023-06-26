module m_file_reader
  use str_module
  use m_vstring
  use m_vstringlist
  use list_module, only: stack, list_node
  implicit none

!-------------------------------------------------------------------------------------------------!
  type file_reader
    
    character(100)           :: file
    integer                  :: unit
    integer                  :: ioerr
    type(stack)              :: backstack
  
    contains
      procedure,private      :: jump_to_external
      procedure,private      :: handle_eof
      procedure,private      :: read_to_next_line
      procedure,public       :: next_block_id
      procedure,public       :: next_block_item
      !procedure,public :: get_block_length
    
  end type file_reader
!-------------------------------------------------------------------------------------------------!

  integer                     :: max_unit
  character(500)              :: line
  
  contains

!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine file_reader_module_initialize()
    implicit none
    
    max_unit = 8

  end subroutine file_reader_module_initialize

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
  
  subroutine return_iunit(unit_no)
    implicit none
    integer, intent(in) :: unit_no
    
    ! TODO - more robust handling
    max_unit = unit_no - 1
    
  end subroutine return_iunit
  
!-------------------------------------------------------------------------------------------------!
  
  function open_file(file) result(f)
    implicit none
    
    type(file_reader),pointer  :: f
    character(*), intent(in) :: file
    
    allocate(f)
    
    f%unit = get_next_iunit()
    f%file = trim(file)
    
    open(f%unit, file=f%file, action='READ', STATUS='OLD', iostat=f%ioerr)
    ! TODO ADD IOMSG & ISOMSG HANLDR
    ! TODO RECORD
    
  end function open_file

!-------------------------------------------------------------------------------------------------!

  subroutine close_file(f)
    implicit none
    type(file_reader),pointer :: f
    close(f%unit)
    
    call return_iunit(f%unit)
    
    !deallocate(f)
    !nullify(f)
  end subroutine close_file
  
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
  subroutine item2char(strings, i, out_char)
    implicit none
    class(t_vstringlist),intent(in)   :: strings
    integer, intent(in)               :: i
    character(*),intent(out)             :: out_char
    
    call vstring_cast(vstrlist_index(strings, i), out_char)
  end subroutine item2char
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine read_to_next_line(this, eof)
    implicit none
    class(file_reader)  :: this
    type(file_reader),pointer   :: next_file
    integer             :: ierr
    integer,intent(out) :: eof
    
    ! For line splitting
    integer             :: left(2), right(2)

    eof = 0
    do while (eof == 0)
      read(this%unit, '(a500)', iostat=ierr) line
      line = adjustl(line)
      ! Handle EOF
      if (IS_IOSTAT_END(ierr)) then
        call this%handle_eof(eof)
      else if ((line(1:1) == '#').or.(line(1:1) == '*').or.(line(1:1) == '')) then
        ! Comment line, blank line
        continue
      else if (index(line, 'EXTERNAL')>0) then
        ! External command, continue reading in new file
        call multisplit(ierr,2,left,right,line)
        call this%jump_to_external(line(left(2):right(2)))
        !call this%read_to_next_line(eof)
      else
        ! Data
        backspace(this%unit)
        exit
      end if
    end do
  
  end subroutine
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine jump_to_external(this, file)
    implicit none
    class(file_reader)  :: this
    character(*), intent(in) :: file
    
    ! Store old iunit
    call this%backstack%push(this%unit)
    
    ! Update where we're reading from
    this%unit = get_next_iunit()
  
    open(this%unit, file=file, action='READ', STATUS='OLD', iostat=this%ioerr)
    ! TODO ADD IOMSG & ISOMSG HANLDR
    ! TODO RECORD
    
  end subroutine jump_to_external
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine handle_eof(this, eof)
    implicit none
    class(file_reader)     :: this
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
    class(file_reader)                :: this
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
      ! TODO Log
      call vstring_new(line_string, trim(adjustl(line)))
      strings = vstrlist_split(line_string)
    
    else
      ! End of file - just return
    end if
  
  end subroutine next_item
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine next_block_id(this, eof, id)  
    implicit none
    class(file_reader)        :: this
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
    class(file_reader)        :: this
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
      ! End of file - just return
    end if
  
  end subroutine next_block_item

!-------------------------------------------------------------------------------------------------!
  
end module m_file_reader