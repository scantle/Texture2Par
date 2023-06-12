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
    character(100) :: id
  
    reader => open_file(file)
    
    eof_check = 0
    do
      call reader%read_next_block_id(eof_check, id)
      if (eof_check == 1) exit
      call block_redirect(id)
    end do
    write(*,*) 'DONE READING!'
    
    call close_file(reader)
  
  end subroutine read_main_input

!-------------------------------------------------------------------------------------------------!
  
    subroutine block_redirect(id)
    implicit none

    character(*), intent(in)        :: id

    select case(trim(id))
      case("OPTIONS")
        write(*,*) 'yeahhh'
        !CALL read_OPTIONS_block(unit)
      case("FLOW_MODEL")
        !CALL read_FLOW_MODEL_block(unit)
      case("CLASSES")
        !CALL read_class_block(unit)
      case("DATASET")
        !CALL read_dataset_block(unit)
      case("VARIOGRAMS")
        !CALL read_variogram_block(unit)
      case("PP_LOCS")
        !CALL read_pploc_block(unit)
      case("PP_PARAMETERS")
        !CALL read_ppparm_block(unit)
      case DEFAULT
        ! TODO! Raise error!
        write(*,*) "BAD BLOCK!"
        stop
    end select

  end subroutine block_redirect
  
!-------------------------------------------------------------------------------------------------!
  
end module read_main_input_module