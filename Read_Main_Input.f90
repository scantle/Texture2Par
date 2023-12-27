module m_read_main_input
  use m_T2P, only: t_t2p
  use m_vstring, only: t_vstring, vstring_new
  use m_vstringlist, only: t_vstringlist, vstrlist_search, vstrlist_split
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

    ! Regardless of block
    call pre_block_read_check(t2p, id)
    
    select case(trim(id))
      case("OPTIONS")
        ! This is the only block we let read it's own entries and set it's own read tracker (has_opt)
        ! In the future we may need a post-read hook to process this info, but for now we keep it simple
        call t2p%opt%read_options(reader)
        t2p%has_opt = .true.
      case("FLOW_MODEL")
        call read_FLOW_MODEL_block(t2p, reader)
      case("CLASSES")
        CALL read_class_block(t2p, reader)
      case("DATASET")
        CALL read_dataset_block(t2p, reader)
      case("VARIOGRAMS")
        CALL read_variogram_block(t2p, reader)
      case("PP_LOCS")
        call pre_block_read_check(t2p, id)
        !CALL read_pploc_block(t2p, reader)
      case("PP_PARAMETERS")
        call pre_block_read_check(t2p, id)
        !CALL read_ppparm_block(t2p, reader)
      case DEFAULT
        ! TODO! Raise error!
        call error_handler(1,reader%file,"Unknown Block Name: " // trim(id))
    end select

  end subroutine block_redirect
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine pre_block_read_check(t2p, block_id)
    ! Generic hook for making sure a block has its prerequisite data read in.
    implicit none
    type(t_t2p),intent(inout)       :: t2p
    character(*), intent(in)        :: block_id
    
    select case(trim(block_id))
      case("CLASSES")
        if (t2p%has_opt==.false.) call error_handler(2,reader%file,"OPTIONS block must be declared before CLASSES")
      case("VARIOGRAMS")
        if (t2p%has_classes==.false.) call error_handler(2,reader%file,"CLASSES block must be declared before VARIOGRAMS")
      case("PP_PARAMETERS")
        if (t2p%has_pplocs==.false.) call error_handler(2,reader%file,"PP_LOCS block must be declared before PP_PARAMETERS")
      case DEFAULT
      ! Do nothing
    end select
  
  end subroutine pre_block_read_check
  
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
    integer                       :: i, j, k, n, prime_id, class_id
    type(t_vstringlist),pointer   :: class_names
    integer,pointer               :: class_connections(:,:), class_type(:)
    
    allocate(class_names)
    call vstrlist_new (class_names)
    n = t2p%opt%DATA_CLASSES
    allocate(class_connections(n, n)) ! symmetric n x n
    allocate(class_type(n))
    
    ! Setup matrix, values indicate variogram ID. Can start off assuming each class has a variogram (fill diagonal)
    class_connections = 0
    do i = 1, n
      class_connections(i, i) = i
    end do
    
    i = 0   ! counts primary datasets
    j = n+1 ! tracks relationship id
    
    ! Loop until end of block or end of file
    do
      i = i + 1
      call reader%next_block_item(status, id, strings, length)
      if (status /= 0) exit  ! exit if end of block or end of file
      if (i+(j-n-1) > t2p%opt%DATA_CLASSES) then
        call error_handler(1,reader%file,"Too many data classes listed in CLASSES block. Increase DATA_CLASSES")
      end if
      ! Add items to arrays
      call vstrlist_append(class_names, strings)
      prime_id = vstrlist_search(class_names, vstrlist_index(strings, 1))
      class_type(prime_id) = 1
      ! Loop over specifying the cross relationships
      do k=2, length
        class_id = vstrlist_search(class_names, vstrlist_index(strings, k))
        class_connections(prime_id,class_id) = j
        class_connections(class_id,prime_id) = j ! symmetric
        class_type(class_id) = 2
        j = j + 1
      end do
    end do
    
    ! Some early validation
    if (j < 1) then
      call error_handler(2,reader%file,"No Data Classes Listed")
    else if (j < t2p%opt%DATA_CLASSES) then
      call error_handler(2,reader%file,"Too Few Data Classes Listed")
    end if
    
    ! Pass to T2P
    call t2p%setup_data_classes(i - 1, class_names, class_connections, class_type)
    
  end subroutine read_class_block
  
!-------------------------------------------------------------------------------------------------!
  
subroutine read_dataset_block(t2p, reader)
  implicit none
  type(t_t2p),intent(inout)       :: t2p
  type(t_file_reader), pointer    :: reader
  
  ! Standard file reader variables
  integer                    :: status, length
  character(30)              :: id, value
  type(t_vstringlist)        :: strings
  
  ! Sub variables
  logical                    :: fname_set
  type(t_vstring)            :: fname, ftrans
  
  ! Default values
  call vstring_new (ftrans , "NONE" )
    
  ! Loop until end of block or end of file
  do
    call reader%next_block_item(status, id, strings, length)
    if (status /= 0) exit  ! exit if end of block or end of file
    select case(trim(id))
      case("FILE")
        fname = vstrlist_index(strings , 2)
        fname_set = .true.
      case("TRANSFORM")
        ftrans = vstring_toupper(vstrlist_index(strings , 2))
      case DEFAULT
        call error_handler(1,reader%file,"Unknown Dataset option: " // trim(id))
    end select
  end do
  
  ! Did at least filename get set??
  if (fname_set == .false.) then
    call error_handler(2,reader%file,"Dataset Block must include a FILE")
  end if
  
  ! Add to T2P
  call vstrlist_append(t2p%dataset_files, fname)
  call vstrlist_append(t2p%dataset_trans, ftrans)
      
end subroutine read_dataset_block

!-------------------------------------------------------------------------------------------------!
  
subroutine read_variogram_block(t2p, reader)
  use m_vario, only: t_vgm, get_vgm_from_line
  implicit none
  type(t_t2p),intent(inout)       :: t2p
  type(t_file_reader), pointer    :: reader
  
  ! Standard file reader variables
  integer                    :: status, length
  character(30)              :: id, value
  type(t_vstringlist)        :: strings
  
  ! Sub variables
  type(t_vgm),pointer        :: vgms(:,:)
  integer,pointer            :: nnear(:), nstruct(:)
  integer                    :: class_id, cross_id, vario_id, struct_id, vcount
  logical                    :: iscross, ispp
  type(t_vstring)            :: class_name, cross_name
  type(t_vstringlist)        :: class_temp
  character(100)             :: msg
  
  ! Allocate variograms
  ! nvario is currently set during the class block read, max structures in settings
  allocate(vgms(t2p%opt%max_vstruct,t2p%nvario), nnear(t2p%nvario), nstruct(t2p%nvario))
  vcount = 0
  nnear = 0
  nstruct = 0
  
  ! Loop until end of block or end of file
  do
    call reader%next_block_item(status, id, strings, length)
    if (status /= 0) exit  ! exit if end of block or end of file
    select case(trim(id))
      case("CLASS")
        ! Now we need to figure out what variogram we have, and read all the structures for it
        iscross = .false.
        ispp    = .false.
        class_name = vstrlist_index(strings , 2)
        
        ! Is this a single class covariogram or is this a cross covariogram?
        if (length > 3 .and. vstrlist_search(strings,":")) then
          !Covar case 1 class_name : cross_name
          cross_name = vstrlist_index(strings,4)
          ! TODO if the cross name doesn't exist
          iscross = .true.
        else if (vstring_match(class_name,"*:*")) then
          !Covar case 2 No space class_name:cross_name
          class_temp = vstrlist_split(class_name, ":")
          class_name = vstrlist_index(class_temp, 1)
          cross_name = vstrlist_index(class_temp, 2)
          ! TODO if the cross name doesn't exist
          iscross = .true.
        !else if (vstring_match(class_name  ! WHAT ABOUT PPs?
        end if
        
        ! Get variogram ids
        ! Check for pilot point class
        if (vstring_match(class_name,"PilotPoints",nocase=.true.)) then
          ispp = .true.
          vario_id = t2p%nvario  ! PP vario is max
        else
        ! Not PP
          class_id = vstrlist_search(t2p%class_names, class_name)
          ! TODO if the class name doesn't exist
          if (iscross) then
            cross_id = vstrlist_search(t2p%class_names, cross_name)
          else
            cross_id = class_id
          end if
          vario_id = t2p%class_connections(class_id, cross_id)
          !TODO Validate vario_id
        end if
        
        ! Now we know what variogram, loop and read in structures
        do
          call reader%next_block_item(status, id, strings, length)
          if (status /= 0) exit       ! exit if end of block or end of file
          if (trim(id)=="CLASS") then ! move back one line & exit if it's the next class
            call reader%back(1)
            exit 
          end if
          ! Read structure
          struct_id = item2int(strings, 1)
          nstruct(vario_id) = struct_id
          !TODO validate struct_id
          ! Get full variogram item from the rest of the line
          vgms(struct_id,vario_id) = get_vgm_from_line(strings, struct_id)
          ! NNear on first struct read, advance vcount
          if (struct_id == 1) then
            nnear(vario_id) = item2int(strings, 11)
            vcount = vcount + 1
          end if
        end do
        ! End of variogram structure reading, moves to next class variogram or end of block
      case DEFAULT
        call error_handler(1,reader%file,"Unknown Variogram option: " // trim(id))
    end select
  end do ! end of variogram reading loop
  
  ! Check we have all variograms, send to t2p object
  if (vcount < t2p%nvario) then
    write(msg, '(2(a,i3))') "Too Few Variograms, expected ", t2p%nvario, ", found: ", vcount
    call error_handler(2,reader%file,msg)
  end if
  call t2p%setup_variograms(vgms, nnear, nstruct)

end subroutine read_variogram_block



subroutine read_pploc_block(t2p, reader)
  implicit none
  type(t_t2p),intent(inout)       :: t2p
  type(t_file_reader), pointer    :: reader
  
  ! Standard file reader variables
  integer                    :: status, length
  character(30)              :: id, value
  type(t_vstringlist)        :: strings
  
  

end subroutine read_pploc_block

end module m_read_main_input