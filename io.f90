
module io
  character(1024)               :: prefix='', ipcfile='', extlpcfile(100)='', mfdis='', mf6dis='', mf6disv='', mf6disu=''
  character(24)                 :: fmtHK,fmtVK,fmtSS,fmtSY
  integer                       :: count_extlpc = 0
  logical                       :: mf6=.false., kriging=.false., nonexistok=.true.
  namelist /options/ prefix,fmtHK,fmtVK,fmtSS,fmtSY,extlpcfile,ipcfile,mfdis,mf6dis,mf6disv,mf6disu
contains


  subroutine read_options(nml)
    character(*)                  :: nml
    integer                       :: ifile, ierr
    character(1024)               :: ioerrmsg

    fmtHK = '10es10.3'
    fmtVK = '10es10.3'
    fmtSS = '10es10.3'
    fmtSY = '10f10.6'

    open (status='old', file=trim(nml), newunit=ifile)
    read (ifile, nml=options, iostat=ierr, iomsg=ioerrmsg)
    if (ierr /= 0) then
      print*, 'Namelist file error: '//ioerrmsg
      stop
    end if
    count_extlpc = count(extlpcfile /= '')
    if (count_extlpc>0 .and. ipcfile == '') then
      print *, 'ipcfile is not defined in namelist.'
      stop
    end if
    mf6 = (mf6dis/='' .or. mf6disu/=''.or. mf6disv/='')
  end subroutine

!-----------------------------------------------------------------------------!
  SUBROUTINE read_external_pc(modflowflag)
    use MakePar
    use MFSupport, only : get_well_cell
    use utl
    use InputOutputModule, only : UPCASE, getunit
    implicit none
    integer, intent(in)         :: modflowflag              ! if this is a MODFLOW model
    integer, parameter          :: nnear=4
    real                        :: x, y, z, pc
    real                        :: nodedists(nnear), queryxy(2)
    integer                     :: closenodes(nnear), iextl, ifile, inode, iline, idummy, ilayer, ierr
    integer, allocatable        :: ipc(:,:), nextlpc(:,:)
    character(1024)             :: line, tmpstr

    pcnode = NODATA
    allocate(ipc(nnodes,nlayers), nextlpc(nnodes,nlayers))
    ! read ibound like ipc
    tmpstr = adjustl(trim(ipcfile))
    call UPCASE(tmpstr)
    if (tmpstr(:8) == 'CONSTANT') then
      read(tmpstr(9:), *) idummy
      ipc = idummy
    else
      ifile = getunit() ! Assign a free unopened unit number to the iu dummy argument.
      open(unit=ifile, file=trim(ipcfile), status='old')
      read(ifile, *) ipc
      close(ifile)
    end if
    kriging = any(ipc==0)
    nextlpc = 0
    print *, 'Reading External Percent Coarse Files:'
    do iextl = 1, count_extlpc
      print *, ' - '//trim(extlpcfile(iextl))
      open(newunit=ifile, file=trim(extlpcfile(iextl)), status='old')
      read(ifile,*,iostat=ierr) line ! Header
      where(ipc==iextl)
        nextlpc = 0
        pcnode = 0
      end where
      iline = 1
      do
        read(ifile,*,iostat=ierr) x, y, z, pc ! coordinates and
        iline = iline + 1
        if (ierr < 0) then
          exit
        else if (pc<0 .or. pc>1. .or. ierr>0) then
          write(*,'(A,I9)'), 'Incorrect percent coarse identified in file '//trim(extlpcfile(iextl))//' Line ', iline
          stop
        end if
        ! find the point location (from inputwells_IWFM and inputwells_modflow, which these two should be combined)
        queryxy(1) = x
        queryxy(2) = y
        call n_nearest_to(nodetree, queryxy, nnear, closenodes, nodedists)
        if (modflowflag==0) then
          inode = get_well_element(x, y, idummy, closenodes, nnear)
        else if (modflowflag==1) then
          inode = get_well_cell(x, y, closenodes, nnear, .true.)
        end if
        ! find layer the texture is located in
        ilayer = binarySearch(topelev(inode,:), z, outside_ok=.true.)
        ! add the percent value and count
        pcnode(inode,ilayer)  = pcnode(inode,ilayer) + pc
        nextlpc(inode,ilayer) = nextlpc(inode,ilayer) + 1
      end do
      ! check if all the nodes are defined externally
      write(tmpstr, *) count(nextlpc>0)
      write(*,'(A)'), '  - The percent coarse in '//trim(adjustl(tmpstr))//' model cells have been defined.'
      if (any(pack(nextlpc, ipc==iextl)==0)) then
        write(*,'(A)'), '    - Some model cells do not contain percent coarse in the file.'
        if (nonexistok) then
          write(*,'(A)'), '    - Percent coarse in these cells will be kriged from well logs.'
          !where (ipc==iextl .and. nextlpc==0) ipc = 0 ! convert to kriging of no percent coarse is found in a cell
        else
          stop
        end if
      end if
      close(ifile)
      write(*,'(A)'), 'Some model cells do not contain percent coarse in file '//trim(extlpcfile(iextl))
    end do
    where(ipc > 0 .and. nextlpc>0)
      ! calculate the average
      pcnode = pcnode / nextlpc
      ! create fake geozones for external
      geozones = geoname_extlpc
    end where
  end subroutine read_external_pc

  subroutine read_mf6dis()
    use BaseDisModule,          only:  DisBaseType
    use GwfDisModule,           only:  dis_cr, GwfDisType
    use GwfDisvModule,          only:  disv_cr, GwfDisvType
    use GwfDisuModule,          only:  disu_cr, GwfDisuType
    use InputOutputModule,      only:  openfile, getunit
    use MakePar,                only:  nnodes,nlayers,nodex,nodey,nzones,botelev,topelev,elevation,aqtardtopelev,aqtardbotelev
    use MFSupport,              only:  allocate_modflow, xoff, yoff, mfrot
    use GLOBAL,                 ONLY:  NLAY, NROW, NCOL, DELR, DELC,LAYCBD
    ! dummy
    class(DisBaseType), pointer     :: dis  => null()
    integer                         :: indis
    integer                         :: iout = 6
    integer                         :: k
    character(:),allocatable        :: name_model

    name_model = 'T2P_MF6'

    indis = getunit() ! Assign a free unopened unit number to the iu dummy argument.
    if(mf6dis /= '') then
      open(unit=indis, file=mf6dis, status='old')
      call dis_cr(dis, name_model, indis, iout)
    elseif(mf6disv /= '') then
      open(unit=indis, file=mf6disv, status='old')
      call disv_cr(dis, name_model, indis, iout)
    elseif(mf6disu /= '') then
      open(unit=indis, file=mf6disu, status='old')
      call disu_cr(dis, name_model, indis, iout)
    else
      return
    endif
    ! read file
    call dis%dis_df()
    close(indis)
    ! prepare data for T2P
    xoff = dis%xorigin
    yoff = dis%yorigin
    mfrot = dis%angrot

    select type(dis)
    type is (GwfDisType)
      print*, 'Regular MF6 Grid with Layer, Row, Column: ', dis%mshape
      ! todo: prepare data for T2P
      nlayers = dis%nlay
      nnodes = dis%nrow * dis%ncol

      if (.not. associated(NLAY))  allocate(NLAY)
      if (.not. associated(NROW))  allocate(NROW)
      if (.not. associated(NCOL))  allocate(NCOL)
      if (.not. allocated(DELR))   allocate(DELR(NCOL))
      if (.not. allocated(DELC))   allocate(DELC(NROW))
      if (.not. allocated(LAYCBD)) allocate(LAYCBD(NLAY))

      NLAY = dis%nlay
      NROW = dis%nrow
      NCOL = dis%ncol

      DELR = dis%delr
      DELC = dis%delc

      call allocate_modflow()
      Elevation            = reshape(dis%top2d,(/nnodes/))
      do k = 1, NLAY
        botelev(:,k)       = reshape(dis%bot3d(:,:,k), (/nnodes/))
      end do
      topelev(:,1)         = Elevation
      topelev(:,2:nlayers) = botelev(:,:(nlayers-1))
      aqtardtopelev = NODATA
      aqtardbotelev = NODATA
    type is (GwfDisvType)
      print*, 'DISv Grid Layer, Not supported yet'
      stop
    type is (GwfDisuType)
      print*, 'Unstructured Grid, Not supported yet'
    end select
  end subroutine

  subroutine write_data()
    use MakePar, only: KhB,KvB,SsB,SyB,pcnode
    if (prefix/='') then
      open(newunit=ifile, file=trim(prefix)//'.hk', status='replace')
      write(ifile, '('//trim(fmtHK)//')') KhB
      close(ifile)
      open(newunit=ifile, file=trim(prefix)//'.vk', status='replace')
      write(ifile, '('//trim(fmtVK)//')') KvB
      close(ifile)
      open(newunit=ifile, file=trim(prefix)//'.ss', status='replace')
      write(ifile, '('//trim(fmtSS)//')') SsB
      close(ifile)
      open(newunit=ifile, file=trim(prefix)//'.sy', status='replace')
      write(ifile, '('//trim(fmtSy)//')') SyB
      close(ifile)
      open(newunit=ifile, file=trim(prefix)//'.pc', status='replace')
      write(ifile, '(20f10.6)') pcnode
      close(ifile)
    end if
  end subroutine
end module
