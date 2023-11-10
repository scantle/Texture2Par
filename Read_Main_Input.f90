module m_read_main_input
  use m_T2P, only: t_t2p
  use m_vstring, only: t_vstring, vstring_new
  use m_vstringlist, only: t_vstringlist
  use m_file_io
  use m_error_handler, only: error_handler
  implicit none

  type(t_file_reader), pointer :: reader
  
  contains

!-------------------------------------------------------------------------------------------------!

  subroutine read_main_input(t2p, file)
    implicit none
    
    type(t_t2p),intent(inout)   :: t2p
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
    implicit none

    type(t_t2p),intent(inout)       :: t2p
    character(*), intent(in)        :: id

    select case(trim(id))
      case("OPTIONS")
        call t2p%opt%read_options(reader)
      case("FLOW_MODEL")
        call read_FLOW_MODEL_block(t2p, reader)
      case("CLASSES")
        CALL read_class_block(t2p, reader)
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
        call error_handler(1,reader%file,"Unknown Block Name: " // trim(id))
    end select

  end subroutine block_redirect
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine read_FLOW_MODEL_block(t2p, reader)
    use m_flow_model, only: t_flow_model, create_flow_model_object, check_id_validity
    implicit none
    type(t_t2p),intent(inout)       :: t2p
    type(t_file_reader), pointer    :: reader
    
    ! Standard file reader variables
    integer                    :: status, length
    character(30)              :: id, value
    type(t_vstringlist)        :: strings
    
    ! Loop until end of block or end of file
    do
      call reader%next_block_item(status, id, strings, length)
      if (status /= 0) exit  ! exit if end of block or end of file
      select case(trim(id))
        case("TYPE")
          call item2char(strings, 2, value, toupper=.true.)
          t2p%fmodel => create_flow_model_object(value)
          !call create_flow_model_object(value, t2p%fmodel)
        case("SIM_FILE")
          call check_id_validity(id, t2p%fmodel, 2) ! Only IWFM
          call item2char(strings, 2, t2p%fmodel%sim_file)
        case("PREPROC_FILE")
          call check_id_validity(id, t2p%fmodel, 2) ! Only IWFM
          call item2char(strings, 2, t2p%fmodel%preproc_file)
        case("NAME_FILE")
          call check_id_validity(id, t2p%fmodel, 3) ! only MODFLOW2000
          call item2char(strings, 2, t2p%fmodel%sim_file)
        case("TEMPLATE_FILE")
          call check_id_validity(id, t2p%fmodel)
          call item2char(strings, 2, t2p%fmodel%template_file)
        case("PP_ZONE_FILE")
          call check_id_validity(id, t2p%fmodel)
          call item2char(strings, 2, t2p%fmodel%ppzone_file)
        case("HSU_FILE")
          call check_id_validity(id, t2p%fmodel)
          call item2char(strings, 2, t2p%fmodel%hsu_file)
        case DEFAULT
          call error_handler(1,reader%file,"Unknown Flow Model option: " // trim(id))
      end select
    end do
  
  end subroutine read_FLOW_MODEL_block
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine read_class_block(t2p, reader)
    
    implicit none
    type(t_t2p),intent(inout)       :: t2p
    type(t_file_reader), pointer    :: reader
    
    ! Standard file reader variables
    integer                    :: status, length
    character(30)              :: id, value
    type(t_vstringlist)        :: strings
    ! sub variables
    integer                       :: i, j
    type(t_vstringlist),pointer   :: class_names
    integer,pointer               :: class_connections(:)
    
    allocate(class_names)
    call vstrlist_new ( class_names )
    allocate(class_connections(t2p%opt%DATA_CLASSES))
    i = 0  ! counts primary datasets
    j = 1  ! tracks place in class_connections
    
    ! Loop until end of block or end of file
    do
      i = i + 1
      call reader%next_block_item(status, id, strings, length)
      if (status /= 0) exit  ! exit if end of block or end of file
       if (j > t2p%opt%DATA_CLASSES) then
        call error_handler(1,reader%file,"Too many data classes listed in CLASSES block. Increase DATA_CLASSES")
      end if
      ! Add items to arrays
      call vstrlist_append(class_names, strings)
      class_connections(j) = i
      if (length > 1) class_connections(j+1:j+length-1) = j
      j = j + length
    end do
    
    ! Some early validation
    if (j < 1) then
      call error_handler(2,reader%file,"No Data Classes Listed")
    else if (j < t2p%opt%DATA_CLASSES) then
      call error_handler(2,reader%file,"Too Few Data Classes Listed")
    end if
    
    ! Pass to T2P
    call t2p%setup_data_classes(class_names, class_connections)
    
  end subroutine read_class_block
  
!-------------------------------------------------------------------------------------------------!
  
end module m_read_main_input