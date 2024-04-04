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
        call t2p%opt%read_options(reader)
        call t2p%process_options()
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
        CALL read_pploc_block(t2p, reader)
      case("PP_PARAMETERS")
        call pre_block_read_check(t2p, id)
        CALL read_ppparm_block(t2p, reader)
      case ("END")
        ! Do nothing - errant line? Block reader ended early?
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
      case("PP_LOCS")
        if (t2p%has_classes==.false.) call error_handler(2,reader%file,"CLASSES block must be declared before PP_LOCS")
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
        case("SIM_FILE")
          call check_id_validity(id, t2p%fmodel, 2) ! Only IWFM
          call item2char(strings, 2, t2p%fmodel%sim_file)
        case("PREPROC_FILE")
          call check_id_validity(id, t2p%fmodel, 2) ! Only IWFM
          call item2char(strings, 2, t2p%fmodel%preproc_file)
        case("NAME_FILE", "NAM_FILE")
          call check_id_validity(id, t2p%fmodel, 3) ! only MODFLOW2000
          call item2char(strings, 2, t2p%fmodel%sim_file)
        case("GSF_FILE")
          call check_id_validity(id, t2p%fmodel)    ! Grid or MFUSG
          call item2char(strings, 2, t2p%fmodel%grd_file)
          t2p%fmodel%grd_file_type = 'gsf'
        case("TEMPLATE_FILE")
          call check_id_validity(id, t2p%fmodel)
          call item2char(strings, 2, t2p%fmodel%template_file)
        case("PP_ZONE_FILE")
          call check_id_validity(id, t2p%fmodel)
          call item2char(strings, 2, t2p%fmodel%ppzone_file)
        case("HSU_FILE")
          call check_id_validity(id, t2p%fmodel)
          call item2char(strings, 2, t2p%fmodel%hsu_file)
        case("XOFF","XOFFSET")
          call check_id_validity(id, t2p%fmodel, 3) ! only MODFLOW2000
          t2p%fmodel%xoff = item2dp(strings, 2)
        case("YOFF","YOFFSET")
          call check_id_validity(id, t2p%fmodel, 3) ! only MODFLOW2000
          t2p%fmodel%yoff = item2dp(strings, 2)
        case("ROT","ROTATION")
          call check_id_validity(id, t2p%fmodel, 3) ! only MODFLOW2000
          t2p%fmodel%rot= item2dp(strings, 2)
        case DEFAULT
          call error_handler(1,reader%file,"Unknown Flow Model option: " // trim(id))
      end select
    end do
    
    ! Hook for T2P to process flowmodel entries
    call t2p%setup_flow_model()
  
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
    integer                       :: i, j, k, n, nclasses, prime_id, class_id
    type(t_vstringlist),pointer   :: class_names
    integer,pointer               :: class_connections(:,:), class_type(:)
    
    allocate(class_names)
    call vstrlist_new (class_names)
    
    ! Pre-read for number of classes
    call reader%get_block_dim(n, nclasses)  !t2p%opt%DATA_CLASSES
    if (n < 1) then
      call error_handler(2,reader%file,"No Data Classes Listed")
    end if
    allocate(class_connections(nclasses, nclasses)) ! symmetric n x n
    allocate(class_type(nclasses))
    
    ! Setup matrix, values indicate variogram ID. Can start off assuming each class has a variogram (fill diagonal)
    class_connections = 0
    do i = 1, n
      class_connections(i, i) = i
    end do
    
    j = n+1 ! tracks relationship id
    
    ! Loop over data class lines
    do i=1, n
      call reader%next_block_item(status, id, strings, length)
      if (status /= 0) exit  ! exit if end of block or end of file
      !if (i+(j-n-1) > t2p%opt%DATA_CLASSES) then
      !  call error_handler(1,reader%file,"Too many data classes listed in CLASSES block. Increase DATA_CLASSES")
      !end if
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
    read(reader%unit,*)  ! Read past end block
    
    ! Pass to T2P
    call t2p%setup_data_classes(n,nclasses, class_names, class_connections, class_type)
    
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
  t2p%ndatafiles = t2p%ndatafiles + 1
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
  status = 0
  do
    if (status /= 0) exit  ! Exit if inner loop variogram read loop found the end
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
          ! NNear on first struct read, if not cross, advance vcount
          if (struct_id == 1) then
            if (iscross) then
              nnear(vario_id) = 0 ! Can't get n cross points...
            else
              nnear(vario_id) = item2int(strings, 11)
            end if
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

!-------------------------------------------------------------------------------------------------!

subroutine read_pploc_block(t2p, reader)
  use m_pilotpoints, only: t_pilotpoints, create_pilotpoint_object
  implicit none
  type(t_t2p),intent(inout)       :: t2p
  type(t_file_reader), pointer    :: reader
  
  ! Standard file reader variables
  integer                    :: status, length
  character(30)              :: id, value
  type(t_vstringlist)        :: strings
  
  ! Sub variables
  integer,pointer            :: pp_id(:)
  integer,allocatable        :: zone_tracker(:)
  integer                    :: i, j, npp, pid, nzones
  logical                    :: zone_counted
  
  ! Assume block rows = # of pilot points
  npp = reader%get_block_len()
  
  ! Allocate & setup PP object
  allocate(pp_id       (npp), &
           zone_tracker(npp)  )
  t2p%pilot_points => create_pilotpoint_object(npp, t2p%nclasses)
  nzones = 0
  zone_tracker = 0
  
  ! Now read the pilot point locations
  do i=1, npp
    call reader%next_block_item(status, id, strings, length)
    if (status /= 0) exit  ! exit if end of block or end of file
    ! Should this be read with it's own function?
    pp_id(i)                           = item2int(strings, 1)
    t2p%pilot_points%grid%coords(1,i)  = item2dp (strings, 2)
    t2p%pilot_points%grid%coords(2,i)  = item2dp (strings, 3)
    t2p%pilot_points%category(i)       = item2int(strings, 4)
    t2p%pilot_points%lay_type(i)       = find_string_index_in_list(strings, 5, "A,T", .true.)
    
    ! Early input checking & handling
    do j=1, i-1
      if (pp_id(j) == pp_id(i)) call error_handler(1,reader%file,"All pilot point IDs must be unique")
    end do
    if (i == 1) then  ! Make sure first pilot point wasn't aquitard
      if(t2p%pilot_points%lay_type(i)==2) then
        call error_handler(1,reader%file,"Must list all aquifer pilot points (A) before aquitard pilot points (T)")
      end if
    else  ! Check to make sure all aquifer pps listed before aquitard & store max aquifer pp
      if (t2p%pilot_points%lay_type(i)==1 .and. t2p%pilot_points%lay_type(i-1)==2) then
        call error_handler(1,reader%file,"Must list all aquifer pilot points (A) before aquitard pilot points (T)")
      else if (t2p%pilot_points%lay_type(i)==2 .and. t2p%pilot_points%lay_type(i-1)==1) then
        ! First aquitard pilot point - record (TODO: All Necessary?)
        t2p%pilot_points%max_a_cat = maxval(t2p%pilot_points%category(:))
        t2p%pilot_points%pp_split = i
        t2p%pilot_points%zonelist_split = nzones + 1
      end if
      
      ! Adjust aquitard zone id for internal use - offset by maximum aquifer category
      if (t2p%pilot_points%lay_type(i)==2) then
        t2p%pilot_points%category(i) = t2p%pilot_points%max_a_cat + t2p%pilot_points%category(i)
      end if

    end if
    ! Zone tracking
    zone_counted = .false.
    do j=1, nzones
      if (zone_tracker(j) == t2p%pilot_points%category(i)) zone_counted = .true.
    end do
    if (zone_counted == .false.) then
      nzones = nzones + 1
      zone_tracker(nzones) = t2p%pilot_points%category(i)
    end if
  end do
  
  ! Should be at block end, check to be sure!
  call reader%next_block_item(status, id, strings, length)
  if (trim(id)/="END") then
    ! Assumption is the input file is formatted incorrectly... hopefully that's why
    call error_handler(1,reader%file,"PP_LOCS section improperly formatted.")
  end if
  
  ! Run pplocs setup hook
  call t2p%setup_pplocs(pp_id, nzones, zone_tracker)
  
end subroutine read_pploc_block

!-------------------------------------------------------------------------------------------------!
  
subroutine read_ppparm_block(t2p, reader)
  implicit none
  type(t_t2p),intent(inout)       :: t2p
  type(t_file_reader), pointer    :: reader
  
  ! Standard file reader variables
  integer                    :: status, length
  character(30)              :: id, value
  type(t_vstringlist)        :: strings
  
  ! Subroutine variables
  integer                    :: pp_index
  
  ! Loop until end of block or end of file
  status = 0
  do
    if (status /= 0) exit  ! Exit if inner loop pp read loop found the end
    call reader%next_block_item(status, id, strings, length)
    if (status /= 0) exit  ! exit if end of block or end of file
    select case(trim(id))
      case("TYPE")  ! Only valid id (for now)
        call item2char(strings, 2, value, toupper=.true.)  ! Get type
        ! Internal Loop over pilot points, whatever type they are (global or class)
        do
          call reader%next_block_item(status, id, strings, length)
          if (status /= 0) exit  ! exit if end of block or end of file
          if (trim(id)=="TYPE") then ! move back one line & exit if it's the next class
            call reader%back(1)
            exit
          end if
          ! Read the line, save
          call t2p%pilot_points%read_pp_line(id, value, strings, length, t2p%class_names)
        end do ! End of pp reading loop
      case DEFAULT
        call error_handler(1,reader%file,"Unknown Pilot Point Parameter option: " // trim(id))
    end select
  end do
  
  ! May need to be replaced with setup hook
  t2p%has_pp = .true.
  
end subroutine read_ppparm_block

!-------------------------------------------------------------------------------------------------!
  
end module m_read_main_input