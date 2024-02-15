module m_read_datasets
  use m_T2P, only: t_t2p
  use m_vstring, only: t_vstring, vstring_new
  use m_vstringlist, only: t_vstringlist, vstrlist_search, vstrlist_split, vstrlist_length
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
    integer                       :: i, ndat, eof
    
    ! Loop over files
    do i=1, ndat
      call item2char(t2p%dataset_files, i, file)
      reader => open_file_reader(file)
      ! Read in header
      call reader%next_item(eof, strings)
      
      ! Get data length
      ndat = reader%get_data_len()
      
      ! Allocate
      
      ! Rewind
      rewind(reader%unit)
      
      ! Check
      read(reader%unit) file
    end do
  
  end subroutine read_datasets
!-------------------------------------------------------------------------------------------------!
  
end module m_read_datasets