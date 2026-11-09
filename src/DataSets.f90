module m_datasets
  use m_grid, only: t_grid, create_grid_object
  implicit none
!-------------------------------------------------------------------------------------------------!
  type  :: t_dataset  ! Abstract Base Class

    integer,allocatable            :: category(:)      ! Zones, Hydrostratigraphic Units (HSUs)
    type(t_grid),pointer           :: grid             ! Grid
    ! Doesn't even have any data! Because array dimensions differ, subclasses are controversially
    ! allowed to have their own data storage

    contains
      procedure,nopass             :: initialize_abc
      generic,public               :: initialize => initialize_abc
      procedure,public             :: finalize
  end type t_dataset
!-------------------------------------------------------------------------------------------------!
!  type, extends(t_dataset) :: t_intervaldata
!
!    character(30)                  :: name             ! Texture class
!    integer                        :: class_id         ! Class ID
!    integer,allocatable            :: nearest_node(:)  ! Nearest model node (cell) to data point
!    real    ,allocatable           :: values(:,:)      ! Data values associated with intervals at each point (interval, point)
!
!    contains
!
!      procedure,public             :: initialize => init_interval
!
!  end type t_intervaldata
!-------------------------------------------------------------------------------------------------!
  type, extends(t_dataset) :: t_layerdata
    ! Intended for input data
    character(30)                  :: name             ! Class Name  (UNUSED? TODO)
    integer                        :: id               ! Index within t2p%datasets
    real    ,allocatable           :: values(:,:)      ! Data values associated with each grid point for each class (nclasses, npoints)
    real    ,allocatable           :: variance(:,:)    ! Variance of each class dataset (nclasses, npoints). Can be a user input or interpolation result.
    integer ,allocatable           :: fmodel_loc(:)    ! Containing model element/cell
    real    ,allocatable           :: mean(:)          ! Mean of values, by class, for cokriging standardization (globally, not by layer)

    contains
      procedure                    :: init_layerdata
      generic                      :: initialize => init_layerdata
      procedure,public             :: finalize => finalize_layerdata
  end type t_layerdata
!-------------------------------------------------------------------------------------------------!
!  type, extends(t_dataset) :: t_layermultidata
!    ! Intended for intermediate data
!    character(30)                  :: name             ! Class Name
!    integer                        :: class_id         ! Class ID
!    real    ,allocatable           :: values(:,:)      ! Data values associated with each grid point for each class (nclasses, npoints)
!    integer,allocatable            :: fmodel_loc(:)    ! Containing model element/cell
!
!    contains
!      procedure                    :: init_layermultidata
!      generic                      :: initialize => init_layermultidata
!
!  end type t_layermultidata
!-------------------------------------------------------------------------------------------------!
! See also: t_pilotpoint derived class in PilotPoints.f90
!-------------------------------------------------------------------------------------------------!

  contains

!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!

  subroutine write_layerdata(&
    fname,         & ! output filename
    layerdata,     & ! array of layerdata
    coords,        & ! coordinates
    kstart,        & ! starting layer for ourput
    kintrvl,       & ! layer step for output (used for IWFM)
    classid,       & ! index of values
    write_variance,&
    gridncol,      & ! optional; number of columns for structured grid; default is 0 for unstructured grid
    lname,         & ! optional; prefix of layer name
    rownames,      & ! optional; row names written in the first column in the output file
    modelnodes,    & ! optional; index of the node numbers of the output, used when output corresponding cell of the wells
    mdim,          & ! optional; number of dimensions; default is 2
    nodecol        & ! optional; column name of the node number
    
  )
    implicit none
    character(*)               :: fname
    type(t_layerdata),pointer  :: layerdata(:)  ! (nlayers)
    integer,intent(in)         :: kstart, kintrvl, classid
    real   ,intent(in)         :: coords(:,:)
    integer,intent(in),optional:: gridncol, mdim, modelnodes(:)
    logical,intent(in)         :: write_variance
    character(*),intent(in),optional :: lname, rownames(0:), nodecol
    ! local
    integer                    :: inode, igrid, i, k, row, col, unit, ncol, ndim, nlayers
    integer, allocatable       :: inodes(:)
    integer                    :: nnodes
    character(:),allocatable   :: header, fmt, nodename
    character(20)              :: tmpstr

    nlayers = size(layerdata)
    nodename = "Node"
    if (present(nodecol)) nodename = trim(nodecol)
    if (present(modelnodes)) then
      inodes = modelnodes
    else
      nnodes = size(layerdata(1)%values, dim=2)
      allocate(inodes(nnodes))
      !inodes = [(inode, inode=1,size(layerdata(1)%values, dim=2))]
      do i = 1, nnodes
        inodes(i) = i
      end do
    end if
    ndim = 2
    if (present(mdim)) ndim = mdim
    nnodes = size(inodes)
    ncol = 0
    if (present(gridncol)) then
      ncol = gridncol
    end if
    if (ncol>0) then
      header="Row,Column,X,Y"
      fmt="I0,',',I0,2(',',G0.10)"
    else
      header=nodename//",X,Y"
      fmt="I0,2(',',G0.10)"
    end if
    if (present(rownames)) then
      header = trim(rownames(0))//","//header
      fmt = "A,',',"//fmt
    end if
    if (ndim==3) then
      header = header//",Z"
      fmt = fmt//",',',G0.6"
    end if

    do k=1, int(nlayers/kintrvl)
      if (present(lname)) then
        write(tmpstr, "(2A,I0)") ",", trim(lname), k
      else
        write(tmpstr, "(2A,I0)") ",","Layer", k
      end if
      header = header//trim(tmpstr)
    end do
    fmt = "("//fmt//",*(:',',G0.6))"
    
    if (write_variance) then
      do k = 1, int(nlayers/kintrvl)
        write(tmpstr,"(A,I0)") ",Var_Layer",k
        header = header//trim(tmpstr)
      end do
      !fmt = "("//fmt//",*(:',',G0.6))"
    end if

    open(newunit=unit, file=trim(fname), status='replace')

    ! Header
    write(unit,"(A)") header

    ! Values
    do inode = 1, nnodes
      igrid = inodes(inode)
      if (ncol>0) then
        row = int((igrid-1)/ncol) + 1
        col = igrid - (row-1)*ncol
        if (present(rownames)) then
          if (write_variance) then
            write(unit,fmt) trim(rownames(inode)), row, col, coords(1:ndim,inode), (layerdata(k)%values(classid,inode),  k=kstart,nlayers,kintrvl), (layerdata(k)%variance(classid,inode), k=kstart,nlayers,kintrvl)
          else
            write(unit,fmt) trim(rownames(inode)), row, col, coords(1:ndim,inode), (layerdata(k)%values(classid,inode),  k=kstart,nlayers,kintrvl)
          end if
        else
          if (write_variance) then
            write(unit,fmt) row, col, coords(1:ndim,inode), (layerdata(k)%values(classid,inode),  k=kstart,nlayers,kintrvl), (layerdata(k)%variance(classid,inode), k=kstart,nlayers,kintrvl)
          else
            write(unit,fmt) row, col, coords(1:ndim,inode), (layerdata(k)%values(classid,inode),  k=kstart,nlayers,kintrvl)
          end if
        end if
      else
        if (present(rownames)) then
          if (write_variance) then
            write(unit,fmt) trim(rownames(inode)), igrid, coords(1:ndim,inode), (layerdata(k)%values(classid,inode),  k=kstart,nlayers,kintrvl), (layerdata(k)%variance(classid,inode), k=kstart,nlayers,kintrvl)
          else
            write(unit,fmt) trim(rownames(inode)), igrid, coords(1:ndim,inode), (layerdata(k)%values(classid,inode),  k=kstart,nlayers,kintrvl)
          end if
        else
          if (write_variance) then
            write(unit,fmt) igrid, coords(1:ndim,inode), (layerdata(k)%values(classid,inode),  k=kstart,nlayers,kintrvl), (layerdata(k)%variance(classid,inode), k=kstart,nlayers,kintrvl)
          else
            write(unit,fmt) igrid, coords(1:ndim,inode), (layerdata(k)%values(classid,inode),  k=kstart,nlayers,kintrvl)
          end if
        end if
      end if
    end do
    
    close(unit)
  end subroutine write_layerdata

!-------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------!
! BASE CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

  subroutine initialize_abc(this)
    implicit none
    class(t_dataset),intent(inout)  :: this

  end subroutine initialize_abc

!-------------------------------------------------------------------------------------------------!

  subroutine finalize(this)
    implicit none
    class(t_dataset),intent(inout) :: this

    if (allocated(this%category)) deallocate(this%category)
    if (associated(this%grid)) nullify(this%grid)

  end subroutine finalize

!-------------------------------------------------------------------------------------------------!
! LAYER CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

  subroutine init_layerdata(this, npoints, nclasses, default_value)
    implicit none
    class(t_layerdata),intent(inout)  :: this
    integer,intent(in)                :: npoints, nclasses
    real    ,optional                 :: default_value

    allocate(this%category  (npoints), &
             this%fmodel_loc(npoints), &
             this%values    (nclasses, npoints), &
             this%variance  (nclasses, npoints), &
             this%mean      (nclasses))
    allocate(this%grid)

    ! By default
    this%category = 1
    this%fmodel_loc = 0

    if (present(default_value)) this%values = default_value
    this%variance = 0.0

    !this%grid = create_grid_object(2, npoints)  ! X,Y,Z coords by n

  end subroutine init_layerdata

!-------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------!
! LAYER MULTI-DATA CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
!
!  subroutine init_layermultidata(this, npoints, nclasses, grid)
!    implicit none
!    class(t_layermultidata),intent(inout)  :: this
!    integer,intent(in)                     :: npoints, nclasses
!    class(t_grid),pointer                  :: grid
!
!    allocate(this%category  (npoints), &
!             this%values    (nclasses, npoints))
!    !allocate(this%grid)  ! Already should exist, avoid copying
!    this%grid => grid
!
!    !this%grid = create_grid_object(2, npoints)  ! X,Y,Z coords by n
!
!  end subroutine init_layermultidata
!
!-------------------------------------------------------------------------------------------------!

  subroutine finalize_layerdata(this)
    implicit none
    class(t_layerdata),intent(inout) :: this

    if (allocated(this%category  )) deallocate(this%category  )
    if (allocated(this%values    )) deallocate(this%values    )
    if (allocated(this%variance  )) deallocate(this%variance  )
    if (allocated(this%fmodel_loc)) deallocate(this%fmodel_loc)
    if (allocated(this%mean      )) deallocate(this%mean      )
    if (associated(this%grid     )) nullify(this%grid)

  end subroutine finalize_layerdata

!-------------------------------------------------------------------------------------------------!

end module m_datasets
