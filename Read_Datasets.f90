module m_read_datasets
  use m_datasets
  use m_T2P, only: t_t2p
  use m_vstring, only: t_vstring, vstring_new
  use m_vstringlist, only: t_vstringlist, vstrlist_search, vstrlist_split, vstrlist_length,vstrlist_index
  use m_file_io
  use m_error_handler, only: error_handler
  implicit none

  type(t_file_reader), pointer :: reader
  
  contains

!-------------------------------------------------------------------------------------------------!
  subroutine read_datasets(t2p)
    implicit none
    class(t_T2P), intent(inout)   :: t2p
    type(t_vstringlist)           :: strings
    character(200)                :: file
    integer                       :: i, ndat, nlines, npoints, nsets, eof
    integer,pointer               :: ids(:)
    type(t_layerdata),pointer     :: dataset
    
    ! Loop over files
    do i=1, ndat
      call item2char(t2p%dataset_files, i, file)
      reader => open_file_reader(file)
      
      ! Read in header
      call reader%next_item(eof, strings)
      
      ! Process header
      call process_header(t2p, file, strings, nsets, ids)
      
      ! Get data length
      nlines = reader%get_data_len()
      ! Get npoints
      backspace(reader%unit)
      call reader%next_item(eof, strings)
      npoints = item2int(strings,2)
      
      ! Setup object(s)
      dataset = create_layerdata_object(npoints)
      
      ! Rewind
      rewind(reader%unit)
      
      ! Check
      read(reader%unit) file
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
    integer                            :: i
    character(30)                      :: badclass
    
    ! Columns of header should be Location ID n X Y Zland Depth follwed by class data
    !                             1        2  3 4 5 6     7
    nsets = vstrlist_length(header) - 7
    allocate(ids(nsets))
    
    do i=1, nsets
      ids(i) = t2p%get_class_id(vstrlist_index(header,7+i))
      if (ids(i)==0) then
        call item2char(header,7+i,badclass)
        call error_handler(1,filename=file,opt_msg="Bad class name = " // badclass)
      end if
    end do
  
  end subroutine process_header
  
!-------------------------------------------------------------------------------------------------!
  
  
!-------------------------------------------------------------------------------------------------!
  
end module m_read_datasets