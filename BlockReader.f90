module block_reader_module
  use str_module, only: to_upper,count_fragments,multisplit
  implicit none

  contains
  
!-------------------------------------------------------------------------------------------------!



!-------------------------------------------------------------------------------------------------!
  
  subroutine get_block_line(block_name, unit, line, n)
    implicit none
    
    character(*), intent(in)       :: block_name
    integer, intent(in)            :: unit
    integer                        :: ierr
    character(300), intent(out)    :: line
    integer, intent(out)           :: n
    
    ! Initialize
    n = 0
    
    do
      read(unit,'(a300)', iostat=ierr) line
      line = trim(line)
      if (ierr /= 0) then                          ! End of file, etc TODO: Handle error 
        stop 
      else if (line== "".or.line(1:1)=="#") then   ! Blank or Comment Line
        cycle
      else if (line(1:3)=="END") then              ! End of block (return negative value for n)
        n = -1
      else                                         ! Hopefully data
        n = count_fragments(line)                  ! Get number of items
      end if
    end do
  
  end subroutine get_block_line
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine generic_read(unit)
    implicit none
    
    integer, intent(in)      :: unit
    integer                  :: n, arrlen, ierr
    integer,allocatable      :: lw(:), rw(:)      ! TODO replace with pointers
    character(300)           :: line
    
    ! Setup left and right character fragment location arrays (see multsplit)
    arrlen=15
    allocate(lw(arrlen), rw(arrlen))
    
    do
      call get_block_line("BLOCK", unit, line, n)
      if (n > arrlen) then
        ! Unwarranted robustness
        arrlen = n*2
        deallocate(lw, rw)
        allocate(lw(arrlen), rw(arrlen))
      end if
      call multisplit(ierr,n,lw,rw,line)
      ! Begin actual data conditions
      
    end do
  
  end subroutine generic_read
  
!-------------------------------------------------------------------------------------------------!

end module block_reader_module