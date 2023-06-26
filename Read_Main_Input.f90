module m_read_main_input
  use m_file_reader, only: file_reader, open_file, close_file
  implicit none

  type(file_reader), pointer :: reader
  
  contains

!-------------------------------------------------------------------------------------------------!

  subroutine read_main_input(file)
    implicit none
    character(100), intent(in)  :: file
    integer                     :: eof_check
    character(100)              :: id
  
    reader => open_file(file)
    
    eof_check = 0
    do
      call reader%next_block_id(eof_check, id)
      if (eof_check == -1) exit
      call block_redirect(id)
    end do
    write(*,*) 'DONE READING!'
    
    call close_file(reader)
  
  end subroutine read_main_input

!-------------------------------------------------------------------------------------------------!
  
    subroutine block_redirect(id)
    use m_global, only: t2p
    implicit none

    character(*), intent(in)        :: id

    select case(trim(id))
      case("OPTIONS")
        write(*,*) 'yeahhh'
        call t2p%opt%read_options(reader)
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
  
end module m_read_main_input