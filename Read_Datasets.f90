module m_read_datasets
  use m_datasets
  use m_T2P, only: t_t2p
  use m_vstring, only: t_vstring, vstring_is
  use m_vstringlist, only: t_vstringlist, vstrlist_search, vstrlist_split, vstrlist_length,vstrlist_index
  use m_file_io
  use m_error_handler, only: error_handler, warning_handler
  implicit none

  type(t_file_reader), pointer :: reader
  
  contains

!-------------------------------------------------------------------------------------------------!
  subroutine read_datasets(t2p)
    implicit none
    class(t_T2P), intent(inout)   :: t2p
    type(t_vstringlist)           :: strings
    character(200)                :: file
    integer                       :: i, j, ndat, nlines, npoints, nsets, eof, idx
    integer,pointer               :: ids(:)
    
    ! Loop over files
    do i=1, t2p%ndatafiles
      call item2char(t2p%dataset_files, i, file)
      reader => open_file_reader(file)
      
      ! Read in header
      call reader%next_item(eof, strings)
      call process_header(t2p, file, strings, nsets, ids)
      
      ! Get data length
      nlines = reader%get_data_len()
      ! Get npoints
      backspace(reader%unit)  ! EOF
      backspace(reader%unit)  ! Last line of data
      call reader%next_item(eof, strings)
      npoints = item2int(strings,2)  ! Reasonable assumption? 
      ! Allocate dataset(s)
      do j=1, nsets
        call t2p%datasets(ids(j))%initialize(npoints)
      end do
      
      ! Rewind
      rewind(reader%unit)
      ! Skip header this time...
      call reader%skip(1, eof)
      
      ! The real work begins
      do j=1, nlines
        
      end do
      
    end do
  
  end subroutine read_datasets
!-------------------------------------------------------------------------------------------------!

  subroutine process_header(t2p, file, header, nsets, ids)
    implicit none
    class(t_T2P), intent(inout)        :: t2p
    character(*),intent(in)            :: file
    type(t_vstringlist),intent(in)     :: header
    integer,intent(out)                :: nsets
    integer,pointer                    :: ids(:)
    integer,pointer                    :: temp(:)
    integer                            :: i, count
    character(30)                      :: badclass
    character(5)                       :: cint
    
    ! Columns of header should be Location ID n X Y Zland Depth follwed by class data, then layers
    !                             1        2  3 4 5 6     7
    nsets = vstrlist_length(header) - 7
    if (t2p%has_hsus) nsets = nsets - t2p%fmodel%nlayers
    ! TODO check to make zure nsets > 0
    allocate(temp(nsets))
    count = 0
    
    do i=1, nsets
      temp(i) = t2p%get_class_id(vstrlist_index(header,7+i))
      if (temp(i)==0) then
        if (vstring_is(vstrlist_index(header,7+i),'integer')) then
          ! Assumption - an integer indicates this is an (unused) HSU zone == done with classes
          write(cint,'(i5)') count
          call warning_handler("Encountered integer in header for "//trim(file)//" - assuming "//trim(adjustl(cint))//" classes in file")
          exit
        end if
        call item2char(header,7+i,badclass)
        call error_handler(1,filename=file,opt_msg="Bad class name = " // badclass)
      end if
      count = count + 1
    end do
    
    nsets = count
    allocate(ids(nsets))
    ids(:) = temp(:nsets)
    deallocate(temp)
  
  end subroutine process_header
  
!-------------------------------------------------------------------------------------------------!
  
  
!-------------------------------------------------------------------------------------------------!
  
end module m_read_datasets