module read_main_input_module
  use file_reader_module, only: file_reader, open_file, close_file
  implicit none

  type(file_reader), pointer :: reader
  
  contains

!-------------------------------------------------------------------------------------------------!

  subroutine read_main_input(file)
    implicit none
    character(100) :: file
    integer        :: eof_check
    character(500) :: line
  
    reader => open_file(file)
    
    eof_check = 0
    do
      call reader%read_to_next_line(eof_check)
      ! TODO write function to identify/grab blocks
      if (eof_check == 0) then
        read(9, '(a500)') line
        write(*,*) trim(line)
      else
        exit
      end if
    end do
    
    call close_file(reader)
  
  end subroutine read_main_input

!-------------------------------------------------------------------------------------------------!
end module read_main_input_module