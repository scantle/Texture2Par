module m_read_datasets
  use m_global, only: dp, NODATA
!  use m_datasets
  use m_T2P, only: t_t2p
  use m_vstring, only: t_vstring, vstring_is
  use m_vstringlist, only: t_vstringlist, vstrlist_search, vstrlist_split, vstrlist_length,vstrlist_index
  use m_file_io, only: t_file_reader, open_file_reader, item2int, item2char
  use m_error_handler, only: error_handler, warning_handler
  use m_grid, only: t_grid, create_grid_object
  use m_data_transformer, only: t_data_transformer
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
    integer                       :: i, j, k, l, ndat, nlines, npoints, nsets, ncats, eof
    integer                       :: prev_idx, idx, iint, fmodel_loc, didx
    integer,allocatable           :: cids(:)
    real(dp)                      :: x, y, z, prev_depth, depth, intv_top, intv_bot, thick, frac_thick
    real(dp),allocatable          :: wvalues(:), wv_accumulator(:,:), local_elevs(:), total_thick(:,:)
    type(t_grid),pointer          :: grid
    type(t_data_transformer)      :: data_transformer
    
    ! Allocate datasets (is this the right place for this to happen?)
    allocate(t2p%datasets(t2p%ndatasets, t2p%fmodel%nlayers))
    
    ! Setup for HSUs - if a HSU file was given/read, then they are expected in the dataset files
    ! Should we skip if ncat is 1???
    if (t2p%fmodel%hsus%from_file) then
      ncats = t2p%fmodel%nlayers
    else
      ncats = 0
    end if
    
    ! Loop over files
    do i=1, t2p%ndatafiles
      call item2char(t2p%dataset_files, i, file)
      reader => open_file_reader(file)
      
      ! Setup data transformer based on information passed to main input file
      ! TODO implement additional parameter(s)
      call data_transformer%set_transform(vstrlist_index(t2p%dataset_trans,i))
      
      ! Read in header
      call reader%next_item(eof, strings)
      call process_header(t2p, file, strings, nsets, didx, cids)
      
      ! Allocate for reading
      allocate(wvalues      (nsets), &
              total_thick   (nsets,t2p%fmodel%nlayers), &
              wv_accumulator(nsets,t2p%fmodel%nlayers), &
              hsus(ncats))
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
        call t2p%datasets(didx,k)%initialize(npoints, nsets)
        t2p%datasets(didx,k)%values = NODATA
        t2p%datasets(didx,k)%grid => grid
      end do
      
      ! For texture classes, if opt%INTERP_ORDER=AFTER we need to store ppzones for these points, too!
      if (didx==1.and.t2p%opt%INTERP_ORDER=="AFTER") then
        call t2p%setup_dataset_ppzones()
      end if
      
      ! Rewind
      rewind(reader%unit)
      ! Skip header this time...
      call reader%skip(1, eof)
      
      ! The real work begins
      do j=1, nlines
        read(reader%unit, *) name, idx, iint, x, y, z, depth, wvalues(1:nsets), hsus(1:ncats)
        
        if (iint==1) then
          ! Store previous log values...
          if (j>1) then
            do k=1, t2p%fmodel%nlayers
              do l=1, nsets
                if (total_thick(l,k)>0.0) then
                  t2p%datasets(didx,k)%values(cids(l), prev_idx) = wv_accumulator(l,k) / total_thick(l,k)
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
          ! Get location in model, get elevation at model layers in model
          fmodel_loc = t2p%fmodel%find_in_grid(grid%coords(:,idx))
          ! If using, get ppzone in model
          if (t2p%opt%INTERP_ORDER=="AFTER") then
            t2p%dataset_ppzones%array(idx,1:) = t2p%dataset_ppzones%array(fmodel_loc,1:)
          end if
          ! Store HSUs for each layer
          do k=1, ncats
            do l=1, nsets
              t2p%datasets(didx,k)%category(idx) = t2p%fmodel%hsus%get_id(hsus(k))
              if (t2p%datasets(didx,k)%category(idx)==0) then
                call error_handler(1,file,"Invalid HSU: "//trim(hsus(k))//" (No matching category in HSU file)")
              end if
            end do
          end do
          call t2p%fmodel%get_elev_at_point(grid%coords(:,idx), local_elevs(:), fmodel_loc)
          prev_idx = idx
        else
          thick = depth - prev_depth
          intv_top = z - prev_depth
        end if  ! End of new log hanlding
        
        ! Start of accumulator
        ! Transform values
        call data_transformer%apply(wvalues)
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
            t2p%datasets(didx,k)%values(cids(l),prev_idx) = wv_accumulator(l,k) / total_thick(l,k)
          end if
        end do
      end do
      
      ! Cleanup loop variables
      deallocate(wvalues, total_thick, wv_accumulator, hsus)
    end do
  
  end subroutine read_datasets
!-------------------------------------------------------------------------------------------------!

  subroutine process_header(t2p, file, header, nsets, dataset_idx, class_ids)
    implicit none
    class(t_T2P), intent(inout)        :: t2p
    character(*),intent(in)            :: file
    type(t_vstringlist),intent(in)     :: header
    integer,intent(out)                :: nsets, dataset_idx
    integer,allocatable                :: temp(:)
    integer,allocatable,intent(out)    :: class_ids(:)
    integer                            :: i, count
    character(30)                      :: badclass
    character(5)                       :: cint
    
    ! Columns of header should be Location ID n X Y Zland Depth follwed by class data, then layer groups
    !                             1        2  3 4 5 6     7
    nsets = vstrlist_length(header) - 7
    if (t2p%has_hsus) nsets = nsets - t2p%fmodel%nlayers
    if (nsets < 1) then
       call error_handler(1,file,"Invalid number of entries in dataset header")
    end if
    
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
        call error_handler(1,file,"Bad class name = " // badclass)
      end if
      count = count + 1
    end do
    
    ! Move to output variables
    nsets = count
    allocate(class_ids(nsets))
    class_ids(:) = temp(1:nsets)
    deallocate(temp)
    
    ! Validate datasets either (1) are all the primary datasets or (2) are a secondary dataset
    ! Return index of t2p dataset object to store them in
    if (maxval(class_ids) <= t2p%get_active_texture_classes()) then
      ! is primary...
      if (nsets==t2p%get_active_texture_classes()) then
        dataset_idx = 1
        t2p%class_dataset(class_ids) = dataset_idx
      else
        call error_handler(1,file,"All primary datasets must be in same dataset file")
      end if
    else
      ! Secondary: we adjust all indexes relative to the secondary datasets
      ! (Each in their own dataset)
      if (nsets==1) then
        dataset_idx = class_ids(1) - t2p%get_active_texture_classes()
        t2p%class_dataset(class_ids(1)) = dataset_idx
        class_ids(:) = 1
      else
        call error_handler(1,file,"Each secondary dataset must have their own dataset file")
      end if
    end if
  end subroutine process_header
  
!-------------------------------------------------------------------------------------------------!
  
end module m_read_datasets