module m_read_datasets
  use m_global ,only: NODATA, log
!  use m_datasets
  use m_T2P, only: t_t2p
  use str_module, only: multisplit, span2int, span2real, span2char
  use m_vstring, only: t_vstring, vstring_is, vstring_match
  use m_vstringlist, only: t_vstringlist, vstrlist_search, vstrlist_length,vstrlist_index
  use m_file_io, only: t_file_reader, open_file_reader, item2int, item2char
  use m_error_handler, only: error_handler, warning_handler
  use m_grid, only: t_grid, create_grid_object
  use m_data_transformer, only: t_data_transformer
  use m_subsidence, only: accumulate_interbeds, reset_interbeds, build_interbed_mask, is_interbed_interval
  use m_log, only: write_colmap
  use m_colmap, only: t_dataset_colmap
  implicit none

  ! Module variables
  type(t_file_reader), pointer :: reader
  type(t_dataset_colmap)       :: colmap

  contains

!-------------------------------------------------------------------------------------------------!
  subroutine read_datasets(t2p)
    implicit none
    class(t_T2P), intent(inout)   :: t2p

    if (t2p%opt%READ_GRID) then ! Data already interpolated
      call read_gridded_data(t2p)
    else
      ! Business as usual
      if (t2p%opt%INTERP_DIM==2) then
        call read_datasets_2D(t2p)
      else if (t2p%opt%INTERP_DIM==3) then
        call read_datasets_3D(t2p)
      else
        ! Shouldn't ever reach here, input checked during m_options::read_options() subroutine
      end if
    end if

  end subroutine read_datasets

!-------------------------------------------------------------------------------------------------!

  subroutine read_datasets_2D(t2p)
    implicit none
    class(t_T2P), intent(inout)   :: t2p
    type(t_vstringlist)           :: strings
    character(200)                :: file, name
    character(30),allocatable     :: hsus(:)
    character(256)                :: iomsg
    character(2048)               :: line
    character(20)                 :: tmp,tmp2
    integer                       :: i, j, k, l, ndat, nlines, npoints, nsets, ncats, eof, hsu_id
    integer                       :: idx, idx_offset, iint, prev_idx, prev_iint, fmodel_loc, didx, ios
    integer,allocatable           :: cids(:), all_count(:), lw(:), rw(:)
    real                          :: x, y, z, fmodel_z, depth0, depth, intv_top, intv_bot, frac_thick
    real    ,allocatable          :: wvalues(:), wv_accumulator(:,:), local_elevs(:), total_thick(:,:), wvariance(:), var_accumulator(:,:)
    real    ,allocatable          :: all_sum(:), all_depth(:)
    type(t_grid),pointer          :: grid
    type(t_data_transformer)      :: data_transformer
    ! subsidence variables
    real,allocatable              :: kbed_run(:)
    logical,allocatable           :: kbed_open(:), interbed_mask(:)
    logical                       :: is_interbed, inferred_fine

    allocate(t2p%datasets(t2p%ndatasets, t2p%fmodel%nlayers))
    if (t2p%opt%SUBSIDENCE_ON) then
      allocate(t2p%interbed_count(t2p%fmodel%nlayers),&
               t2p%interbed_thick(t2p%fmodel%nlayers),&
               kbed_run(t2p%fmodel%nlayers),&
               kbed_open(t2p%fmodel%nlayers))
    end if

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
      call process_header(t2p, file, strings, nsets, didx, cids, colmap, t2p%dataset_colmap(i))
      call write_colmap(log, colmap, t2p, t2p%dataset_colmap(i)%nsets>0)

      ! Allocate for reading
      allocate(wvalues      (nsets))
      allocate(wvariance    (nsets))
      allocate(all_count    (nsets))
      allocate(all_sum      (nsets))
      allocate(all_depth    (nsets))
      allocate(total_thick    (nsets,t2p%fmodel%nlayers))
      allocate(wv_accumulator (nsets,t2p%fmodel%nlayers))
      allocate(var_accumulator(nsets,t2p%fmodel%nlayers))
      allocate(hsus(ncats))
      allocate(rw(colmap%maxcol), lw(colmap%maxcol))

      var_accumulator = 0.0
      wvariance = 0.0
      if (.not. allocated(local_elevs)) allocate(local_elevs(0:t2p%fmodel%nlayers))
      all_count = 0
      all_sum = 0
      all_depth = 0
      prev_idx = 0
      prev_iint = 0
      ! interbeds (allocated only when subsidence is on)
      if (t2p%opt%SUBSIDENCE_ON) then
        kbed_run = 0.0
        kbed_open = .false.
      end if
      inferred_fine = .false.

      if (t2p%opt%SUBSIDENCE_ON.and.didx==1) then
        ! Handle interbed mask
        allocate(interbed_mask(nsets))
        call build_interbed_mask(nsets, t2p%fine_classes, cids, interbed_mask)
        inferred_fine = t2p%is_inferred_class_fine()
      end if

      ! Get data length (not file length, since we've gone past the header)
      call reader%next_item(eof, strings)
      idx_offset = item2int(strings,2) - 1          ! First well id, must be sequential after
      nlines = reader%get_data_len(rewind_to_last_line=.true.) + 1
      ! Get npoints
      !backspace(reader%unit)  ! EOF
      !backspace(reader%unit)  ! Last line of data
      call reader%next_item(eof, strings)
      npoints = item2int(strings,2) - idx_offset    ! Last - (First - 1) = npoints
      ! Allocate shared grid
      grid => create_grid_object(2, npoints)
      ! Allocate dataset(s) for all layers, set default value
      do k=1, t2p%fmodel%nlayers
        call t2p%datasets(didx,k)%initialize(npoints, nsets)
        t2p%datasets(didx,k)%values = NODATA
        t2p%datasets(didx,k)%grid => grid
        t2p%datasets(didx,k)%id = didx
        if (t2p%opt%SUBSIDENCE_ON.and.didx==1) then
          call t2p%interbed_count(k)%initialize(npoints, 1)
          call t2p%interbed_thick(k)%initialize(npoints, 2)
          t2p%interbed_count(k)%values = NODATA
          t2p%interbed_count(k)%grid => grid
          t2p%interbed_count(k)%id = didx
          t2p%interbed_thick(k)%values = NODATA
          t2p%interbed_thick(k)%grid => grid
          t2p%interbed_thick(k)%id = didx
        end if
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
        write(tmp, '(i)') j+1
        read(reader%unit, "(A)", iostat=reader%ioerr, iomsg=iomsg) line
        call reader%iomsg_handler(iomsg)
        call multisplit(ios, colmap%maxcol, lw, rw, line)
        if (ios /= 0) call error_handler(1,file,"Too few columns on line "//adjustl(trim(tmp)))
        !read(line, *, iostat=reader%ioerr, iomsg=iomsg) name, idx, iint, x, y, z, depth0, depth
        ! --- First 8 fixed columns
        call span2char(line, lw(1), rw(1), name)
        idx   = span2int (line, lw(2), rw(2), ios); if (ios/=0) call error_handler(1,file,"Bad location id on line "//adjustl(trim(tmp)))
        iint  = span2int (line, lw(3), rw(3), ios); if (ios/=0) call error_handler(1,file,"Bad interval id on line "//adjustl(trim(tmp)))
        x     = span2real(line, lw(4), rw(4), ios); if (ios/=0) call error_handler(1,file,"Bad X on line "//adjustl(trim(tmp)))
        y     = span2real(line, lw(5), rw(5), ios); if (ios/=0) call error_handler(1,file,"Bad Y on line "//adjustl(trim(tmp)))
        z     = span2real(line, lw(6), rw(6), ios); if (ios/=0) call error_handler(1,file,"Bad Z on line "//adjustl(trim(tmp)))
        depth0= span2real(line, lw(7), rw(7), ios); if (ios/=0) call error_handler(1,file,"Bad top depth on line "//adjustl(trim(tmp)))
        depth = span2real(line, lw(8), rw(8), ios); if (ios/=0) call error_handler(1,file,"Bad bot depth on line "//adjustl(trim(tmp)))
        ! Then the data columns
        do l = 1, colmap%nsets
          wvalues(l) = span2real(line, lw(colmap%col_val(l)), rw(colmap%col_val(l)), ios)
          if (ios /= 0) call error_handler(1,file,"Bad class value on line "//adjustl(trim(tmp)))
        end do
        if (t2p%opt%KVME) then
          do l = 1, colmap%nsets
            wvariance(l) = span2real(line, lw(colmap%col_var(l)), rw(colmap%col_var(l)), ios)
            if (ios /= 0) call error_handler(1,file,"Bad class variance on line "//adjustl(trim(tmp)))
          end do
        end if
        ! Then the HSUS
        if (t2p%has_hsus) then
          do k = 1, ncats
            call span2char(line, lw(colmap%col_hsu(k)), rw(colmap%col_hsu(k)), hsus(k))
          end do
        end if

        if ((.not. t2p%opt%NO_TEX_BOUNDS).and.didx==1.and.any(wvalues > 1.0)) call error_handler(1,file,"Texture observations cannot be > 1.0, see line: "//adjustl(trim(tmp)))
        all_depth = all_depth + depth

        if (iint==1) then
          idx = idx - idx_offset
          if (idx <= prev_idx) call error_handler(1,file,"Loc IDs must be sequential - Check Loc "//trim(name)//" on line "//adjustl(trim(tmp)))
          ! Check for duplicate well locations
          do k = 1, idx-1
            if (x==grid%coords(1,k).and.y==grid%coords(2,k)) then
              write(tmp , '(i)') idx+idx_offset
              write(tmp2, '(i)') k+idx_offset
              call error_handler(1,file,"Co-located data point at "//trim(name)//". Loc ID "//trim(adjustl(tmp))//" has same XY as ID "//trim(adjustl(tmp2)))
            end if
          end do
          ! Store previous log values...
          if (j>1) then
            do k=1, t2p%fmodel%nlayers
              do l=1, nsets
                if (total_thick(l,k)>0.0) then
                  t2p%datasets(didx,k)%values(cids(l), prev_idx) = wv_accumulator(l,k) / total_thick(l,k)
                  if (t2p%opt%KVME) t2p%datasets(didx,k)%variance(cids(l), prev_idx) = var_accumulator(l,k) / total_thick(l,k)**2
                  all_sum(l) = all_sum(l) + wv_accumulator(l,k) / total_thick(l,k)
                  all_count(l) = all_count(l) + 1
                end if
              end do
            end do
          end if
          ! Store/close open interbeds...
          if (t2p%opt%SUBSIDENCE_ON.and.didx==1) then
            if (j>1) call reset_interbeds(prev_idx, kbed_open, kbed_run, t2p%interbed_thick)
          end if
          ! And reset for new log
          wv_accumulator = 0.0
          var_accumulator = 0.0
          total_thick = 0.0
          intv_top = z
          grid%coords(1,idx) = x
          grid%coords(2,idx) = y
          prev_iint = 0
          ! Get location in model, get elevation at model layers in model
          fmodel_loc = t2p%fmodel%find_in_grid([x,y], t2p%opt%MAX_OUTSIDE_DIST)
          if (fmodel_loc==0) then
            if (t2p%opt%MAX_OUTSIDE_DIST > 0.0) then
              call error_handler(1,file,"Well distance from model greater than MAX_OUTSIDE_DIST: "//name)
            else
              call error_handler(1,file,"Well not found in model domain: "//name)
            end if
          end if
          ! Store location
          do k=1, t2p%fmodel%nlayers
            t2p%datasets(didx,k)%fmodel_loc(idx) = fmodel_loc
            if (t2p%opt%SUBSIDENCE_ON.and.didx==1) then
              t2p%interbed_count(k)%fmodel_loc(idx) = fmodel_loc
              t2p%interbed_thick(k)%fmodel_loc(idx) = fmodel_loc
            end if
          end do
          ! If using, get ppzone in model
          if (t2p%opt%INTERP_ORDER=="AFTER") then
            ! TODO for IWFM this is attempting to match element (fmodel_Loc) to node (ppzones)
            ! Should replace with a fmodel routine get_ppzone_by_loc. For modflow, we'll use the cell, but from iwfm we'll use the nearest node?
            t2p%dataset_ppzones%array(idx,1:) = t2p%dataset_ppzones%array(fmodel_loc,1:)
          end if
          ! Store HSUs for each layer
          do k=1, ncats
            hsu_id = t2p%fmodel%hsus%get_id(hsus(k))
            t2p%datasets(didx,k)%category(idx) = hsu_id
            if (t2p%opt%SUBSIDENCE_ON.and.didx==1) then
              t2p%interbed_count(k)%category(idx) = hsu_id
              t2p%interbed_thick(k)%category(idx) = hsu_id
            end if
            if (t2p%datasets(didx,k)%category(idx)==0) then
              call error_handler(1,file,"Invalid HSU: "//trim(hsus(k))//" (No matching category in HSU file)")
            end if
          end do
          call t2p%fmodel%get_elev_at_point(grid%coords(:,idx), local_elevs(:), fmodel_loc)
          fmodel_z = local_elevs(0)
          if (t2p%opt%USE_MODEL_GSE) intv_top = fmodel_z
          prev_idx = idx
        else if (iint <= prev_iint) then
          call error_handler(1,file,"Well intervals must be sequential - Check Well "//trim(name)//" on line "//adjustl(trim(tmp)))
        end if  ! End of new log handling
        if (depth0 >= depth)  call error_handler(1,file,"Bad well depth: Well "//trim(name)//" on line "//adjustl(trim(tmp))//new_line("")//trim(line)//new_line(""))

        if (t2p%opt%USE_MODEL_GSE) then
          intv_top = fmodel_z - depth0
          intv_bot = fmodel_z - depth
        else
          intv_top = z - depth0
          intv_bot = z - depth
        end if

        ! Interbed thickness accumulator...
        if (t2p%opt%SUBSIDENCE_ON.and.didx==1) then
          is_interbed = is_interbed_interval(wvalues, nsets, interbed_mask, inferred_fine)
          call accumulate_interbeds(local_elevs, intv_top, intv_bot, is_interbed, kbed_open, kbed_run, prev_idx, t2p%interbed_count, t2p%interbed_thick)
        end if

        ! Start of accumulator
        ! Transform values
        if (minval(wvalues)>NODATA) call data_transformer%apply(wvalues)

        ! Accumulate wvalues * interval values for averaging... by layer. It's fun.
        do k=1, t2p%fmodel%nlayers
          if ((intv_top >= local_elevs(k)).and.(intv_bot <= local_elevs(k-1))) then
            ! Find min top/bot values
            frac_thick = min(local_elevs(k-1), intv_top)-max(local_elevs(k), intv_bot)
            do l=1, nsets
              if (wvalues(l)>NODATA) then
                wv_accumulator(l,k) = wv_accumulator(l,k) + wvalues(l) * frac_thick
                if (t2p%opt%KVME) var_accumulator(l,k) = var_accumulator(l,k) + wvariance(l) * frac_thick**2
                total_thick(l,k) = total_thick(l,k) + frac_thick
              end if
            end do
          end if
        end do
        ! NEXT
        prev_iint = iint
      end do  ! end of loop over file lines

      ! Store final log values
      do k=1, t2p%fmodel%nlayers
        do l=1, nsets
          if (total_thick(l,k)>0.0) then
            t2p%datasets(didx,k)%values(cids(l),prev_idx) = wv_accumulator(l,k) / total_thick(l,k)
            if (t2p%opt%KVME) t2p%datasets(didx,k)%variance(cids(l), prev_idx) = var_accumulator(l,k) / total_thick(l,k)**2
            all_sum(l) = all_sum(l) + wv_accumulator(l,k) / total_thick(l,k)
            all_count(l) = all_count(l) + 1
          end if
        end do
      end do

      ! Store/close open interbeds...
      if (t2p%opt%SUBSIDENCE_ON.and.didx==1) then
        call reset_interbeds(prev_idx, kbed_open, kbed_run, t2p%interbed_thick)
      end if

      ! Store means
      do l=1, nsets
        do k=1, t2p%fmodel%nlayers
          t2p%datasets(didx,k)%mean(l) = (all_sum(l) / all_count(l))
        end do
      end do

      ! Cleanup loop variables
      deallocate(wvalues, wvariance, all_count, all_sum, all_depth, total_thick, wv_accumulator, var_accumulator, hsus, rw, lw)
      if (allocated(interbed_mask)) deallocate(interbed_mask)
    end do
    if (allocated(kbed_run)) deallocate(kbed_run)
    if (allocated(kbed_open)) deallocate(kbed_open)
  end subroutine read_datasets_2D

!-------------------------------------------------------------------------------------------------!

  subroutine read_datasets_3D(t2p)
    implicit none
    class(t_T2P), intent(inout)   :: t2p
    type(t_vstringlist)           :: strings
    character(2048)               :: line
    character(200)                :: file, name
    character(30),allocatable     :: hsus(:)
    character(256)                :: iomsg
    character(20)                 :: tmp
    integer                       :: i, j, k, l, nlines, npoints, nsets, ncats, eof, idx_offset, ios, hsu_id
    integer                       :: idx, iint, ipoint, prev_idx, prev_iint, fmodel_loc, didx, resized_intervals
    integer                       :: all_count
    integer,allocatable           :: cids(:), lw(:), rw(:)
    real                          :: x, y, z, depth0, depth, re_depth, thick, fmodel_z, intv_top, intv_bot
    real    ,allocatable          :: wvalues(:), local_elevs(:), all_thick(:), all_sum(:), wvariance(:)
    type(t_grid),pointer          :: grid, interbed_grid
    type(t_data_transformer)      :: data_transformer
    ! subsidence variables
    real,allocatable              :: kbed_run(:)
    logical,allocatable           :: kbed_open(:), interbed_mask(:)
    logical                       :: is_interbed, inferred_fine

    if (t2p%opt%SUBSIDENCE_ON) then
      allocate(t2p%interbed_count(t2p%fmodel%nlayers),&
               t2p%interbed_thick(t2p%fmodel%nlayers),&
               kbed_run(t2p%fmodel%nlayers),&
               kbed_open(t2p%fmodel%nlayers))
    end if

    allocate(t2p%datasets(t2p%ndatasets, 1)) ! For 3D, we store it all in one layer

    ! Setup for HSUs - if a HSU file was given/read, then they are expected in the dataset files
    ! Should we skip if ncat is 1???
    if (t2p%fmodel%hsus%from_file) then
      ncats = 1
    else
      ncats = 0
    end if

    ! Loop over files
    fileloop: do i=1, t2p%ndatafiles
      call item2char(t2p%dataset_files, i, file)
      reader => open_file_reader(file)

      ! Setup data transformer based on information passed to main input file
      ! TODO implement additional parameter(s)
      call data_transformer%set_transform(vstrlist_index(t2p%dataset_trans,i))

      ! Read in header
      call reader%next_item(eof, strings)
      call process_header(t2p, file, strings, nsets, didx, cids, colmap, t2p%dataset_colmap(i))
      call write_colmap(log, colmap, t2p, t2p%dataset_colmap(i)%nsets>0)

      ! Allocate for reading
      allocate(wvalues      (nsets))
      allocate(wvariance    (nsets))
      allocate(all_sum      (nsets))
      allocate(all_thick    (nsets))
      allocate(hsus(ncats))
      allocate(rw(colmap%maxcol), lw(colmap%maxcol))

      wvariance = 0.0
      if (.not. allocated(local_elevs)) allocate(local_elevs(0:t2p%fmodel%nlayers))
      nlines = 0
      all_count = 0
      all_sum = 0
      all_thick = 0
      prev_idx = 0
      prev_iint = 0
      npoints = 0    ! in 3D, npoints =/ nwells. Instead, its the total number of intervals.
      ! interbeds (allocated only when subsidence is on)
      if (t2p%opt%SUBSIDENCE_ON) then
        kbed_run = 0.0
        kbed_open = .false.
      end if
      inferred_fine = .false.

      if (t2p%opt%SUBSIDENCE_ON.and.didx==1) then
        ! Handle interbed mask
        allocate(interbed_mask(nsets))
        call build_interbed_mask(nsets, t2p%fine_classes, cids, interbed_mask)
        inferred_fine = t2p%is_inferred_class_fine()
      end if

      ! Pre-read to get npoints, given opt%MAX_LOG_LENGTH. Logs are to be split evenly into length/MAX_LOG_LENGTH
      do
        read(reader%unit, *, iostat=reader%ioerr, iomsg=iomsg) name, idx, iint, x, y, z, depth0, depth
        if (reader%ioerr /= 0) exit   ! done
        nlines = nlines + 1
        if (nlines==1) idx_offset = idx - 1
        ! How many points will this interval become?
        npoints = npoints + ceiling((depth - depth0) / t2p%opt%MAX_LOG_LENGTH(i))
        ! write(*,*) nlines, npoints
      end do

      ! Initialize, setup
      grid => create_grid_object(3, npoints)
      call t2p%datasets(didx,1)%initialize(npoints, nsets)
      t2p%datasets(didx,1)%values = NODATA
      t2p%datasets(didx,1)%grid => grid
      t2p%datasets(didx,1)%id = didx
      ! Interbeds are still by layer
      if (t2p%opt%SUBSIDENCE_ON.and.didx==1) then
        ! Convert npoints back to it's 2D meaning...
        npoints = idx - idx_offset
        interbed_grid => create_grid_object(2, npoints)
        do k=1, t2p%fmodel%nlayers
          call t2p%interbed_count(k)%initialize(npoints, 1)
          call t2p%interbed_thick(k)%initialize(npoints, 2)
          t2p%interbed_count(k)%values = NODATA
          t2p%interbed_count(k)%grid => interbed_grid
          t2p%interbed_count(k)%id = didx
          t2p%interbed_thick(k)%values = NODATA
          t2p%interbed_thick(k)%grid => interbed_grid
          t2p%interbed_thick(k)%id = didx
        end do
      end if

      ! For texture classes, if opt%INTERP_ORDER=AFTER we need to store ppzones for these points, too!
      if (didx==1.and.t2p%opt%INTERP_ORDER=="AFTER") then
        call t2p%setup_dataset_ppzones()
      end if

      ! Rewind
      rewind(reader%unit)
      ! Skip header this time...
      call reader%skip(1, eof)

      ! The real work begins
      ipoint = 1
      hsu_id = 1
      lineloop: do j=1, nlines
        write(tmp, '(i)') j+1
        read(reader%unit, "(A)", iostat=reader%ioerr, iomsg=iomsg) line
        call reader%iomsg_handler(iomsg)
        call multisplit(ios, colmap%maxcol, lw, rw, line)
        if (ios /= 0) call error_handler(1,file,"Too few columns on line "//adjustl(trim(tmp)))
        !read(line, *, iostat=reader%ioerr, iomsg=iomsg) name, idx, iint, x, y, z, depth0, depth
        ! --- First 8 fixed columns
        call span2char(line, lw(1), rw(1), name)
        idx   = span2int (line, lw(2), rw(2), ios); if (ios/=0) call error_handler(1,file,"Bad location id on line "//adjustl(trim(tmp)))
        iint  = span2int (line, lw(3), rw(3), ios); if (ios/=0) call error_handler(1,file,"Bad interval id on line "//adjustl(trim(tmp)))
        x     = span2real(line, lw(4), rw(4), ios); if (ios/=0) call error_handler(1,file,"Bad X on line "//adjustl(trim(tmp)))
        y     = span2real(line, lw(5), rw(5), ios); if (ios/=0) call error_handler(1,file,"Bad Y on line "//adjustl(trim(tmp)))
        z     = span2real(line, lw(6), rw(6), ios); if (ios/=0) call error_handler(1,file,"Bad Z on line "//adjustl(trim(tmp)))
        depth0= span2real(line, lw(7), rw(7), ios); if (ios/=0) call error_handler(1,file,"Bad top depth on line "//adjustl(trim(tmp)))
        depth = span2real(line, lw(8), rw(8), ios); if (ios/=0) call error_handler(1,file,"Bad bot depth on line "//adjustl(trim(tmp)))
        ! Then the data columns
        do l = 1, colmap%nsets
          wvalues(l) = span2real(line, lw(colmap%col_val(l)), rw(colmap%col_val(l)), ios)
          if (ios /= 0) call error_handler(1,file,"Bad class value on line "//adjustl(trim(tmp)))
        end do
        if (t2p%opt%KVME) then
          do l = 1, colmap%nsets
            wvariance(l) = span2real(line, lw(colmap%col_var(l)), rw(colmap%col_var(l)), ios)
            if (ios /= 0) call error_handler(1,file,"Bad class variance on line "//adjustl(trim(tmp)))
          end do
        end if
        ! Then the HSUS
        if (t2p%has_hsus) then
          do k = 1, ncats
            call span2char(line, lw(colmap%col_hsu(k)), rw(colmap%col_hsu(k)), hsus(k))
          end do
        end if

        if ((.not. t2p%opt%NO_TEX_BOUNDS).and.didx==1.and.any(wvalues > 1.0)) call error_handler(1,file,"Texture observations cannot be > 1.0, see line: "//adjustl(trim(tmp)))
        if (iint==1) then
          idx = idx - idx_offset
          if (idx <= prev_idx) call error_handler(1,file,"Well IDs must be sequential - Check Well "//trim(name)//" on line "//adjustl(trim(tmp)))
          ! Get location in model only if need to align with grid's top elevation OR if subsidence is on (needed for interbeds)
          if (t2p%opt%USE_MODEL_GSE.or.t2p%opt%SUBSIDENCE_ON) then
            fmodel_loc = t2p%fmodel%find_in_grid([x,y], t2p%opt%MAX_OUTSIDE_DIST)
            if (fmodel_loc==0) then
              if (t2p%opt%MAX_OUTSIDE_DIST > 0.0) then
                call error_handler(1,file,"Well distance from model greater than MAX_OUTSIDE_DIST: "//name)
              else
                call error_handler(1,file,"Well not found in model domain: "//name)
              end if
            end if
            call t2p%fmodel%get_elev_at_point([x,y], local_elevs(:), fmodel_loc)
            if (t2p%opt%SUBSIDENCE_ON.and.didx==1) then
              do k=1, t2p%fmodel%nlayers
                t2p%interbed_count(k)%fmodel_loc(idx) = fmodel_loc
                t2p%interbed_thick(k)%fmodel_loc(idx) = fmodel_loc
              end do
            end if
            if (t2p%opt%USE_MODEL_GSE) fmodel_z = local_elevs(0)
          end if
          ! Store/close open interbeds, setup grid
          if (t2p%opt%SUBSIDENCE_ON.and.didx==1) then
            interbed_grid%coords(1, idx) = x
            interbed_grid%coords(2, idx) = y
            if (j>1) call reset_interbeds(prev_idx, kbed_open, kbed_run, t2p%interbed_thick)
          end if
          prev_idx = idx
        else if (iint <= prev_iint) then
          call error_handler(1,file,"Well intervals must be sequential - Check Well "//trim(name)//" on line "//adjustl(trim(tmp)))
        end if  ! End new well if
        if (depth0 >= depth) call error_handler(1,file,"Bad well depth: Well "//trim(name)//" on line "//adjustl(trim(tmp)))

        thick = (depth - depth0)
        resized_intervals = ceiling(thick / t2p%opt%MAX_LOG_LENGTH(i))

        ! Elevations
        if (t2p%opt%USE_MODEL_GSE) then
          intv_top = fmodel_z - depth0
          intv_bot = fmodel_z - depth
        else
          intv_top = z - depth0
          intv_bot = z - depth
        end if

        ! Interbed handling
        if (t2p%opt%SUBSIDENCE_ON.and.didx==1) then
          is_interbed = is_interbed_interval(wvalues, nsets, interbed_mask, inferred_fine)
          call accumulate_interbeds(local_elevs, intv_top, intv_bot, is_interbed, kbed_open, kbed_run, prev_idx, t2p%interbed_count, t2p%interbed_thick)
        end if

        ! Transform values
        if (minval(wvalues)>NODATA) call data_transformer%apply(wvalues)

        do k=1, nsets
          if (wvalues(k) > NODATA) then
            all_sum(k) = all_sum(k) + wvalues(k) * thick
            all_thick(k) = all_thick(k) + thick
          end if
        end do

        ! Store HSU category for this interval - if multiple HSUs, we just take the first one, but this can be expanded in the future if needed
        if (t2p%has_hsus) then
          hsu_id = t2p%fmodel%hsus%get_id(hsus(1))
          if (hsu_id==0) then
            call error_handler(1,file,"Invalid HSU: "//trim(hsus(1))//" (No matching category in HSU file)")
          end if
        end if
        do k=1, resized_intervals
          grid%coords(1, ipoint) = x
          grid%coords(2, ipoint) = y
          re_depth = depth0+(thick/resized_intervals)*(k-0.5)
          if (t2p%opt%USE_MODEL_GSE) then
            grid%coords(3, ipoint) = fmodel_z - re_depth
          else
            grid%coords(3, ipoint) = z - re_depth
          end if
          t2p%datasets(didx,1)%values(:, ipoint) = wvalues(1:nsets)
          t2p%datasets(didx,1)%category(ipoint) = hsu_id
          if (t2p%opt%KVME) t2p%datasets(didx,1)%variance(:, ipoint) = wvariance(1:nsets) * resized_intervals
          ! Store location
          t2p%datasets(didx,1)%fmodel_loc(ipoint) = fmodel_loc
          ipoint = ipoint + 1
        end do

        if (depth>0) all_count = all_count + 1
        prev_iint = iint
      end do lineloop

      ! Store/close open interbeds...
      if (t2p%opt%SUBSIDENCE_ON.and.didx==1) then
        call reset_interbeds(prev_idx, kbed_open, kbed_run, t2p%interbed_thick)
      end if

      ! Store means
      do j=1, nsets
        if (all_thick(j)>0) then
          t2p%datasets(didx,1)%mean(j) = (all_sum(j) / all_thick(j))
        end if
      end do

      ! Cleanup loop variables
      deallocate(wvalues, all_sum, all_thick, hsus, wvariance, rw, lw)
      if (allocated(interbed_mask)) deallocate(interbed_mask)
    end do fileloop
    if (allocated(kbed_run)) deallocate(kbed_run)
    if (allocated(kbed_open)) deallocate(kbed_open)
  end subroutine read_datasets_3D

!-------------------------------------------------------------------------------------------------!

  subroutine process_header(t2p, file, header, nsets, dataset_idx, class_ids, colmap, colmap_user)
    implicit none
    class(t_T2P), intent(inout)           :: t2p
    character(*),intent(in)               :: file
    type(t_vstringlist),intent(in)        :: header
    integer,intent(out)                   :: nsets, dataset_idx
    integer,allocatable,intent(out)       :: class_ids(:)
    type(t_dataset_colmap), intent(inout) :: colmap
    type(t_dataset_colmap), intent(in)    :: colmap_user
    integer,allocatable                   :: tmp_ids(:), tmp_cols(:)
    integer                               :: i, count, nhead, ncol, cid, nactive, idx_var, idx_hsu
    character(64) :: cname
    character(30)                         :: badclass
    character(5)                          :: cint

    ! Columns of header should be Location ID n X Y Zland DepthTop DepthBot followed by data
    !                             1        2  3 4 5 6     7        8
    nhead = 8
    ncol = vstrlist_length(header)
    if ((ncol - nhead) < 1) then
       call error_handler(1,file,"Invalid number of entries in dataset header")
    end if

    ! If the user colmap is set, we can just do a little validation and then get out
    if (colmap_user%nsets>0) then

      ! Move it over!
      colmap = colmap_user

      ! Set some varibles we use outside
      nsets = colmap_user%nsets
      class_ids = colmap_user%class_ids

      if (colmap%maxcol > ncol) then
        call error_handler(1,file,"Explicit column mapping refers to column beyond header length")
      end if

      ! The only thing we couldn't verify during specified colmap read was HSUs
      ! Also, if they only passed one and there are multiple layers, that means they passed the first column
      if (t2p%has_hsus) then
        if (size(colmap%col_hsu)>1) then
          if (size(colmap%col_hsu)/=t2p%fmodel%nlayers) then
            call error_handler(2,opt_msg="DATASET explicit column mapping does not list all HSU columns (needs either first, or one for each layer)")
          end if
        else ! It's only one long
          if (t2p%fmodel%nlayers>1) then ! Assume sequential
            allocate(tmp_cols(t2p%fmodel%nlayers))
            do i=1, t2p%fmodel%nlayers
              tmp_cols(i) = colmap%col_hsu(1) + (i-1)
            end do
            ! and replace
            colmap%col_hsu = tmp_cols
            colmap%maxcol = max(colmap%maxcol,maxval(colmap%col_hsu))
          end if
        end if
      end if

    else ! Autoread

      allocate(tmp_ids(ncol - nhead))
      allocate(tmp_cols(ncol - nhead))
      tmp_ids  = 0
      tmp_cols = 0
      count = 0

      ! Scan ALL columns after nhead and pick those that match known classes.
      ! This makes class columns order-independent and tolerant to extra columns.
      do i = nhead+1, ncol
        cid = t2p%get_class_id(vstrlist_index(header, i))
        if (cid > 0) then
          if (any(tmp_ids(1:count) == cid)) then
            call item2char(header, i, cname)
            call error_handler(1,file,"Duplicate class column in header: "//trim(cname))
          end if
          count = count + 1
          tmp_ids(count)  = cid
          tmp_cols(count) = i
        end if
      end do

      if (count < 1) then
        call error_handler(1,file,"No recognizable data class columns found after location columns")
      end if

      nsets = count
      allocate(class_ids(nsets))
      class_ids(:) = tmp_ids(1:nsets)

      ! Allocate/Init colmap
      colmap%nsets = nsets
      if (allocated(colmap%class_ids)) deallocate(colmap%class_ids)
      if (allocated(colmap%col_val)) deallocate(colmap%col_val)
      allocate(colmap%class_ids(nsets))
      allocate(colmap%col_val(nsets))
      colmap%class_ids(:)  = class_ids(:)
      colmap%col_val(:)    = tmp_cols(1:nsets)
      colmap%maxcol = maxval(colmap%col_val)

      deallocate(tmp_ids, tmp_cols)

      ! KVME variance columns: locate var_<classname> for each class column
      if (t2p%opt%KVME) then
        if (allocated(colmap%col_var)) deallocate(colmap%col_var)
        allocate(colmap%col_var(nsets))
        colmap%col_var = 0

        do i = 1, nsets
          call item2char(header, colmap%col_val(i), cname)
          idx_var = vstrlist_search(header, "var_"//trim(cname), first=nhead+1)
          if (idx_var == 0) then
            call error_handler(1,file,"KVME on: missing variance column 'var_"//trim(cname)//"'")
          end if
          colmap%col_var(i) = idx_var
        end do
        colmap%maxcol = max(colmap%maxcol, maxval(colmap%col_var))
      else
        if (allocated(colmap%col_var)) deallocate(colmap%col_var)
      end if

      ! HSU columns (if enabled)
      if (t2p%has_hsus) then
        if (allocated(colmap%col_hsu)) deallocate(colmap%col_hsu)
        allocate(colmap%col_hsu(t2p%fmodel%nlayers))
        colmap%col_hsu = 0

        ! Prefer explicit header names: hsu1..hsuN
        do i = 1, t2p%fmodel%nlayers
          write(cname,'("hsu",i0)') i
          idx_hsu = vstrlist_search(header, trim(cname), first=nhead+1)
          if (idx_hsu /= 0) colmap%col_hsu(i) = idx_hsu
        end do

        if (any(colmap%col_hsu /= 0)) then
          ! If user uses hsu# naming, require ALL layers to be present
          if (any(colmap%col_hsu == 0)) then
            call error_handler(1,file,"HSU columns: found some 'hsu#' headers but not all 1..N. Use hsu1..hsuN for all layers or none.")
          end if

        else
          ! Fallback: assume HSUs are the last nlayers columns
          if (ncol < nhead + nsets + 1) then
            write(cname,'(i0)') ncol
            call error_handler(1,file,"HSUs expected but header does not have enough columns to contain all layers. "//trim(cname)//" columns found")
          end if
          do i = 1, t2p%fmodel%nlayers
            colmap%col_hsu(i) = ncol - t2p%fmodel%nlayers + i
          end do
        end if
        colmap%maxcol = max(colmap%maxcol, maxval(colmap%col_hsu))
      else
        if (allocated(colmap%col_hsu)) deallocate(colmap%col_hsu)
      end if
    end if ! end of user passed or not if

    ! Validate primary/secondary dataset classes, regardless of user/auto colmap
    nactive = t2p%get_active_texture_classes()

    if (all(class_ids <= t2p%ntexclasses)) then
      ! Primary dataset: must include all active texture classes in one file.
      ! When INFER_LAST_CLASS is active, the inferred final class may still be
      ! present in the source file, but it should not make the file secondary.
      do i = 1, nactive
        if (.not. any(class_ids == i)) then
          call error_handler(1,file,"All active primary datasets must be in same dataset file")
        end if
      end do
      dataset_idx = 1
      t2p%class_dataset(class_ids) = dataset_idx
    else
      ! Secondary dataset: only one class per file
      if (nsets == 1) then
        dataset_idx = class_ids(1) - t2p%ntexclasses + 1
        t2p%class_dataset(class_ids(1)) = dataset_idx
        class_ids(1) = 1
        !colmap%class_ids(1) = 1  colmap needs global class id!
      else
        call error_handler(1,file,"Each secondary dataset must have their own dataset file")
      end if
    end if

  end subroutine process_header

!-------------------------------------------------------------------------------------------------!

  subroutine read_gridded_data(t2p)
  use m_global ,only: NODATA, log
  use str_module, only: multisplit, span2real
  use m_vstringlist, only: t_vstringlist, vstrlist_length, vstrlist_index
  use m_file_io, only: t_file_reader, open_file_reader, item2char
  use m_error_handler, only: error_handler, warning_handler
  implicit none
  class(t_T2P), intent(inout) :: t2p

  type(t_vstringlist)        :: header
  character(200)             :: file
  character(2048)            :: line
  character(256)             :: iomsg
  character(64)              :: cname
  character(20)              :: tmp
  integer                    :: eof, ios
  integer                    :: ncol, nhead, nlines, nactive
  integer                    :: i, j, k
  integer                    :: cid
  integer                    :: lay, node, nval
  integer, allocatable       :: lw(:), rw(:)
  integer, allocatable       :: col(:), cidx(:)   ! mapped class columns -> primary index
  integer, allocatable       :: seen(:)
  integer                    :: nfound, nn, total_expected, total_read
  real                       :: v

  ! Enforce: only one dataset file (no secondary data supported)
  if (t2p%ndatafiles /= 1) then
    write(tmp,'(i0)') t2p%ndatafiles
    call error_handler(1, opt_msg="Only one dataset file can be present when using READ_GRID "//trim(tmp))
  end if

  call item2char(t2p%dataset_files, 1, file)
  reader => open_file_reader(file)

  ! Header + column mapping
  call reader%next_item(eof, header)
  if (eof /= 0) call error_handler(1, file, "Empty gridded dataset file")

  ncol  = vstrlist_length(header)
  nhead = 1 + t2p%opt%INTERP_DIM   ! node id + (x,y[,z])

  if (ncol <= nhead) then
    call error_handler(1, file, "Gridded dataset header has no data columns after node/coords")
  end if

  allocate(lw(ncol), rw(ncol))
  allocate(col(ncol), cidx(ncol))
  col  = 0
  cidx = 0
  nfound = 0

  ! Track presence of each primary class in header
  nactive = t2p%get_active_texture_classes()
  allocate(seen(t2p%nprimary))
  seen = 0

  do i = nhead+1, ncol
    call item2char(header, i, cname)
    cid = t2p%get_class_id(vstrlist_index(header, i), adjust_for_infer=.false.)
    if (cid > 0) then
      if (cid > nactive) then
        call error_handler(1, file, "Secondary/inferred class column not allowed: "//trim(cname))
      end if
      if (seen(cid) == 1) then
        call error_handler(1, file, "Duplicate primary class column in header: "//trim(cname))
      end if
      nfound = nfound + 1
      col(nfound)  = i
      cidx(nfound) = cid
      seen(cid)    = 1
    end if
  end do

  if (nfound < 1) then
    call error_handler(1, file, "No recognizable primary texture class columns found in gridded dataset header")
  end if

  ! If not inferring last class, require all primaries present
  if (.not. t2p%opt%INFER_LAST_CLASS) then
    if (any(seen == 0)) then
      call error_handler(1, file, "READ_GRID: missing one or more primary texture class columns in header")
    end if
  end if

  deallocate(seen)

  ! Allocate output: dset_intermediate(1, nlayers)
  allocate(t2p%dset_intermediate(1, t2p%fmodel%nlayers))

  do k = 1, t2p%fmodel%nlayers
    nn = t2p%fmodel%nnodes(k)
    if (t2p%opt%INFER_LAST_CLASS) then
      call t2p%dset_intermediate(1,k)%initialize(nn, t2p%nprimary+1, default_value=NODATA)
    else
      call t2p%dset_intermediate(1,k)%initialize(nn, t2p%nprimary, default_value=NODATA)
    end if
    t2p%dset_intermediate(1,k)%grid => t2p%fmodel%grid(t2p%fmodel%lay2grid(k))%ptrgrid
    t2p%dset_intermediate(1,k)%id   = 1
  end do

  ! Read body: sequentially fill (layer,node) assuming perfect order
  nlines = reader%get_data_len(rewind_to_last_line=.true.) + 1

  rewind(reader%unit)
  call reader%skip(1, eof)  ! skip header

  lay = 1
  node = 1
  total_read = 0
  total_expected = 0
  do k=1, t2p%fmodel%nlayers
    if (t2p%fmodel%is_active_layer(k)) then
      total_expected = total_expected + t2p%fmodel%nnodes(k)
    else
      if (lay==k) lay = k+1  ! find minimum active layer
    end if
  end do

  do j = 1, nlines
    write(tmp, '(i)') j+1
    read(reader%unit, "(A)", iostat=reader%ioerr, iomsg=iomsg) line
    if (reader%ioerr /= 0) exit
    call reader%iomsg_handler(iomsg)

    call multisplit(ios, ncol, lw, rw, line)
    if (ios /= 0) call error_handler(1, file, "Too few columns on line "//adjustl(trim(tmp)))
    ! Stop if we've already filled the full model
    if (lay > t2p%fmodel%nlayers) exit

    ! Fill only recognized class columns; ignore everything else
    do i = 1, nfound
      v = span2real(line, lw(col(i)), rw(col(i)), ios)
      if (ios /= 0) call error_handler(1, file, "Bad class value on line "//adjustl(trim(tmp)))
      t2p%dset_intermediate(1,lay)%values(cidx(i), node) = v
    end do

    total_read = total_read + 1

    ! Advance node/layer
    node = node + 1
    if (node > t2p%fmodel%nnodes(lay)) then
      lay = lay + 1
      ! Make sure it's an active layer
      do
        if (lay > t2p%fmodel%nlayers) exit
        if (.not.t2p%fmodel%is_active_layer(lay)) then
          lay = lay + 1
        else
          exit
        end if
      end do
      node = 1
    end if
  end do

  ! Post-processing: bounds, infer last class, inactive -> NODATA, means
  do k = 1, t2p%fmodel%nlayers
    if (.not. t2p%fmodel%is_active_layer(k)) cycle
    nval = size(t2p%dset_intermediate(1,k)%values, dim=1)

    if (.not. t2p%opt%NO_TEX_BOUNDS) then
      where(t2p%dset_intermediate(1,k)%values < 0.0) t2p%dset_intermediate(1,k)%values = 0.0
      where(t2p%dset_intermediate(1,k)%values > 1.0) t2p%dset_intermediate(1,k)%values = 1.0
    end if

    if (t2p%opt%INFER_LAST_CLASS) then
      do node = 1, size(t2p%dset_intermediate(1,k)%values, dim=2)
        if (t2p%fmodel%active_nodes_bylay(node,k)) then
          t2p%dset_intermediate(1,k)%values(t2p%nprimary+1, node) = 1.0 - sum(t2p%dset_intermediate(1,k)%values(1:t2p%nprimary, node), dim=1)
        end if
      end do
    end if

    ! Inactive to NODATA
    do node = 1, size(t2p%dset_intermediate(1,k)%values, dim=2)
      if (.not. t2p%fmodel%active_nodes_bylay(node,k)) then
        t2p%dset_intermediate(1,k)%values(:, node) = NODATA
      end if
    end do
  end do

  deallocate(lw, rw, col, cidx)

end subroutine read_gridded_data

!-------------------------------------------------------------------------------------------------!

end module m_read_datasets
