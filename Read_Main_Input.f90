module m_read_main_input
  use m_file_io, only: t_file_reader, open_file_reader, close_file
  use m_flow_model, only: read_FLOW_MODEL_block
  use m_error_handler, only: error_handler
  implicit none

  type(t_file_reader), pointer :: reader
  
  contains

!-------------------------------------------------------------------------------------------------!

  subroutine read_main_input(t2p, file)
    use m_T2P, only: t_t2p
    implicit none
    
    type(t_t2p),intent(inout)       :: t2p
    character(100), intent(in)  :: file
    integer                     :: eof_check
    character(100)              :: id
  
    reader => open_file_reader(file)
    
    eof_check = 0
    do
      call reader%next_block_id(eof_check, id)
      if (eof_check == -1) exit
      call block_redirect(t2p, id)
    end do
    write(*,*) 'DONE READING!'
    
    call close_file(reader)
  
  end subroutine read_main_input

!-------------------------------------------------------------------------------------------------!
  
  subroutine block_redirect(t2p, id)
    use m_T2P, only: t_t2p
    implicit none

    type(t_t2p),intent(inout)       :: t2p
    character(*), intent(in)        :: id

    select case(trim(id))
      case("OPTIONS")
        call t2p%opt%read_options(reader)
      case("FLOW_MODEL")
        call read_FLOW_MODEL_block(reader)
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
        call error_handler(1,reader%file,"Bad Block Name: " // trim(id))
    end select

  end subroutine block_redirect
  
!-------------------------------------------------------------------------------------------------!
  
end module m_read_main_input