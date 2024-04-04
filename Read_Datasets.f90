module m_read_datasets
  use m_global, only: dp, NODATA
!  use m_datasets
  use m_T2P, only: t_t2p
  use m_vstring, only: t_vstring, vstring_is
  use m_vstringlist, only: t_vstringlist, vstrlist_search, vstrlist_split, vstrlist_length,vstrlist_index
  use m_file_io, only: t_file_reader, open_file_reader, item2int, item2char
  use m_error_handler, only: error_handler, warning_handler
  use m_grid, only: t_grid, create_grid_object
  implicit none

  type(t_file_reader), pointer :: reader
  
  contains

!-------------------------------------------------------------------------------------------------!
  subroutine read_datasets(t2p)
    implicit none
    class(t_T2P), intent(inout)   :: t2p
    type(t_vstringlist)           :: strings
    character(200)                :: file, name
    character(30),allocatable     :: hsus(:)
    integer                       :: i, j, k, l, ndat, nlines, npoints, nsets, eof
    integer                       :: prev_idx, idx, iint, fmodel_loc
    real(dp)                      :: x, y, z, prev_depth, depth, intv_top, intv_bot, thick, frac_thick
    real(dp),allocatable          :: wvalues(:), wv_accumulator(:,:), local_elevs(:), total_thick(:,:)
    integer,pointer               :: ids(:)
    type(t_grid),pointer          :: grid
    
    ! Allocate datasets (is this the right place for this to happen?)
    allocate(t2p%datasets(t2p%ndatasets, t2p%fmodel%nlayers))
    
    ! Loop over files
    do i=1, t2p%ndatafiles
      call item2char(t2p%dataset_files, i, file)
      reader => open_file_reader(file)
      
      ! Read in header
      call reader%next_item(eof, strings)
      call process_header(t2p, file, strings, nsets, ids)
      
      ! Allocate for reading
      allocate(wvalues      (nsets), &
              total_thick   (nsets,t2p%fmodel%nlayers), &
              wv_accumulator(nsets,t2p%fmodel%nlayers), &
              hsus(0))  ! TODO handle nlayers of hsus when active
      if (allocated(local_elevs)==.false.) allocate(local_elevs(0:t2p%fmodel%nlayers))
      
      ! Get data length (not file length, since we've gone past the header)
      nlines = reader%get_data_len(rewind_to_last_line=.true.)
      ! Get npoints
      !backspace(reader%unit)  ! EOF
      !backspace(reader%unit)  ! Last line of data
      call reader%next_item(eof, strings)
      npoints = item2int(strings,2)  ! Reasonable assumption? 
      ! Allocate shared grid
      grid => create_grid_object(2, npoints)
      ! Allocate dataset(s) for all layers, set default value
      do k=1, t2p%fmodel%nlayers
        do j=1, nsets
          call t2p%datasets(ids(j),k)%initialize(npoints)
          t2p%datasets(ids(j),k)%values = NODATA
          ! Point all datasets to the same grid object
          t2p%datasets(ids(j),k)%grid => grid
        end do
      end do
      
      ! Rewind
      rewind(reader%unit)
      ! Skip header this time...
      call reader%skip(1, eof)
      
      ! The real work begins
      do j=1, nlines
        read(reader%unit, *) name, idx, iint, x, y, z, depth, wvalues(1:nsets), hsus(1:0)
        if (iint==1) then
          ! Store previous log values...
          if (j>1) then
            do k=1, t2p%fmodel%nlayers
              do l=1, nsets
                if (total_thick(l,k)>0.0) then
                  t2p%datasets(ids(l),k)%values(prev_idx) = wv_accumulator(l,k) / total_thick(l,k)
                end if
              end do
            end do
          end if
          ! And reset for new log
          wv_accumulator = 0.0_dp
          total_thick = 0.0_dp
          intv_top = z
          thick = depth
          grid%coords(1,idx) = x
          grid%coords(2,idx) = y
          ! TODO HSUS
          ! Get location in model, get elevation at model layers in model
          fmodel_loc = t2p%fmodel%find_in_grid(grid%coords(:,idx))
          call t2p%fmodel%get_elev_at_point(grid%coords(:,idx), local_elevs(:), fmodel_loc)
          prev_idx = idx
        else
          thick = depth - prev_depth
          intv_top = z - prev_depth
        end if
        intv_bot = z - depth
        ! Accumulate wvalues * interval values for averaging... by layer. It's fun.
        do k=1, t2p%fmodel%nlayers
          if ((intv_top >= local_elevs(k)).and.(intv_bot <= local_elevs(k-1))) then
            ! Find min top/bot values
            frac_thick = min(local_elevs(k-1), intv_top)-max(local_elevs(k), intv_bot)
            do l=1, nsets
              if (wvalues(l)>=0) then
                wv_accumulator(l,k) = wv_accumulator(l,k) + wvalues(l) * frac_thick
                total_thick(l,k) = total_thick(l,k) + frac_thick
              end if
            end do
          end if
        end do
        prev_depth = depth
      end do
      
      ! Store final log values
      do k=1, t2p%fmodel%nlayers
        do l=1, nsets
          if (total_thick(l,k)>0.0) then
            t2p%datasets(ids(l),k)%values(prev_idx) = wv_accumulator(l,k) / total_thick(l,k)
          end if
        end do
      end do
      
      ! Cleanup loop variables
      deallocate(wvalues, total_thick, wv_accumulator, hsus)
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
    
    ! Columns of header should be Location ID n X Y Zland Depth follwed by class data, then layer groups
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