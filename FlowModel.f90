module m_flow_model
  use m_file_io
  use m_grid, only: t_grid, create_grid_object
  use m_error_handler, only: error_handler
  use m_categories, only: t_category
  !use kdtree2_module, only: kdtree2, kdtree2_create
  implicit none
!-------------------------------------------------------------------------------------------------!
  ! Module variables
  type(t_file_reader), pointer :: reader          ! For reading input files
  integer,parameter            :: max_node_elem_members = 30
!-------------------------------------------------------------------------------------------------!
  
  type :: t_flow_model  ! Abstract Base Class
  ! Base class contains some derived type variables (e.g., the IWFM preproc file) since otherwise
  ! the compiler gets upset when you try to set them on the polymorphic t2p%fmodel type
  ! The alternative was to write a lot of get/set-type methods that would have to be overwritten
  ! by each derived class.
    
    integer                       :: model_type     ! 1 - grid, 2 - IWFM, 3 - MODFLOW 2000
    real(dp)                      :: xoff,yoff,rot  ! Conversion from/to global coordinates
    integer                       :: nlayers        ! Number of model layers
    integer                       :: ngrids         ! Number of unique model grid (xy) layers
    integer,allocatable           :: nnodes(:)      ! Number of nodes (points at which flow is solved) by layer
    integer,pointer               :: lay2grid(:)    ! Which grid each layer uses (e.g., same xy between layers)
    character(30)                 :: name           ! For printing nicely to log/screen
    character(100)                :: preproc_file   ! IWFM main discretization file (empty for others)
    character(100)                :: sim_file       ! Simulation File or MODFLOW name file
    character(100)                :: ppzone_file    ! File of pilot point zones defined for each node
    character(100)                :: hsu_file       ! File of Hydrostratigraphic Units (HSUs) defined for each node
    character(100)                :: template_file  ! File copied for output_file
    character(3)                  :: grd_file_type  ! Character string of grid file type (GSF, ...)
    character(100)                :: grd_file       ! Grid specification file for grid type or MODFLOW-USG
    character(100)                :: output_file    ! File where aquifer parameters are written
    real(dp), allocatable         :: elev(:,:)      ! Bottom Elevations for each layer (node, layer) (indexed from zero = model top)
    type(t_grid),pointer          :: grid(:)        ! xy grids for layers
    type(t_category)              :: ppzones        ! Pilot point zone fmodel grid arrays & tracking
    type(t_category)              :: hsus           ! HSU fmodel grid arrays & tracking
    logical                       :: has_aquitards  ! For now, only in IWFM - but helps T2P w/ aquitard pilot points
  
    contains
      procedure,public            :: initialize
      procedure,public            :: read_model
      procedure,public            :: read_category_files
      procedure,public            :: write_input_summary
      procedure,public            :: find_in_grid
      procedure,public            :: get_elev_at_point
      procedure,public            :: get_midpoint_depths
      procedure,public            :: get_active_ppzones
      procedure,public            :: get_layer_id
      procedure,public            :: get_layer_no   ! For when layer =/ layer_id (e.g., IWFM aquitards)
      procedure,public            :: is_active_layer
      procedure,public            :: is_aquitard_layer
      procedure,public            :: write_fmodel_input
      procedure,public            :: write_node_files
      procedure,public            :: write_node_file
      procedure,private           :: alloc_layers
      procedure,private           :: alloc_nodes
  
  end type t_flow_model
!-------------------------------------------------------------------------------------------------!
  type,extends(t_flow_model) :: t_flow_fdmodel    ! Finite Difference Model
    ! Ineligant solution to both MF-USG and Grid using the GSF
    logical                   :: structured    ! T/F Structured (row/col) model
    integer                   :: nvert
    ! I am not sure we actually need either of these vert arrays since we have node coordinates
    !integer,allocatable       :: vertnnode(:)   ! Number of nodes that include each vertex
    !integer,allocatable       :: vert2node(:,:) ! Nodes that include each vertex
    ! These two are used to confirm a point is within a node (cell) in find_in_grid_gsf
    integer,allocatable       :: nodenvert(:)   ! Number of verticies attached to each node (vert)
    integer,allocatable       :: node2vert(:,:) ! ID of verticies attached to each node (vert, node)
    type(t_grid),pointer      :: vertgrid       ! xyz grid for vertices
    
    contains
      procedure,private       :: read_gsf
      !procedure,private       :: add_node2vert2node
      procedure,public        :: get_elev_at_point => get_elev_at_point_cell
  end type t_flow_fdmodel
!-------------------------------------------------------------------------------------------------!
  type, extends(t_flow_fdmodel) :: t_mgrid
  
  contains
    procedure,public          :: initialize => init_mgrid
    procedure,public          :: read_model => read_mgrid
    procedure,public          :: find_in_grid => find_in_grid_gsf
  
  end type 
!-------------------------------------------------------------------------------------------------!
  type, extends(t_flow_model) :: t_iwfm_model   ! I DO realize the "m" in IWFM is for model
  
    character(100)            :: elem_file
    character(100)            :: node_file
    character(100)            :: strt_file
    integer                   :: nelements
    integer,allocatable       :: nodeid(:)      ! IWFM Node ID for each node
    integer,allocatable       :: nodebyid(:)    ! Given an IWFM node ID, contains internal t2p node index
    integer,allocatable       :: elements(:,:)  ! Comprised of 3-4 nodes (inode, nelem)
    integer,allocatable       :: elemid(:)      ! IWFM Elem ID for each element
    integer,allocatable       :: elemnnodes(:)  ! Number of nodes forming each element
    integer,allocatable       :: nodenelem(:)   ! Number of elements that include each node
    integer,allocatable       :: node2elem(:,:) ! Elements that include each node (elemid,nodes)
    logical,allocatable       :: lay_exist(:)   ! Tracks whether layer has thickness > 0.00
    logical,allocatable       :: lay_aqtard(:)  ! Tracks whether layer is an aquitard
  
  contains
    procedure,public          :: initialize => init_iwfm
    procedure,public          :: read_model => read_iwfm
    procedure,public          :: find_in_grid => find_in_grid_elem
    procedure,public          :: get_elev_at_point =>  get_elev_at_point_elem
    procedure,public          :: get_layer_id => get_layer_id_iwfm
    procedure,public          :: get_layer_no => get_layer_no_iwfm
    procedure,public          :: is_active_layer => is_active_layer_iwfm
    procedure,public          :: is_aquitard_layer => is_aquitard_layer_iwfm
    !procedure,public          :: get_active_ppzones => get_active_ppzones_iwfm
    procedure,public          :: write_fmodel_input => write_fmodel_input_iwfm
    procedure,public          :: write_node_file => write_node_file_iwfm
    procedure,private         :: iwfm_elem_elev_idw
    procedure,private         :: alloc_layers => alloc_layers_iwfm
    procedure,private         :: alloc_nodes  => alloc_nodes_iwfm
    procedure,private         :: alloc_elems
    procedure,private         :: init_nodebyid
    procedure,private         :: process_iwfm_layers
    procedure,private         :: add_elem2nodeelem
  
  end type 
!-------------------------------------------------------------------------------------------------!
  type, extends(t_flow_fdmodel) :: t_modflow_model
  
  integer                     :: layfile       ! LPF/UPW IUNIT number
  character(3)                :: gwftyp        ! LPF or UPW, for printing
  
  contains
    procedure,public          :: initialize => init_modflow2000
    procedure,public          :: read_model => read_modflow2000
    procedure,public          :: find_in_grid => find_in_grid_cell
    procedure,public          :: write_node_file => write_node_file_mf
    procedure,public          :: cell2rowcol
    procedure,private         :: get_cell_local
    procedure,private         :: calc_mfcell_centers_global
  
  end type 
!-------------------------------------------------------------------------------------------------!  
  
  contains
  
!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  function create_flow_model_object(model_type_char) result(object)
    implicit none
    
    character(*),intent(in)        :: model_type_char
    class(t_flow_model),pointer    :: object
    
    select case(model_type_char)
      case ("GRID")
        allocate(t_mgrid::object)
      case ("IWFM")
        allocate(t_iwfm_model::object)
      case ("MODFLOW")
        allocate(t_modflow_model::object)
      case DEFAULT
        call error_handler(1,opt_msg="Unknown Flow Model Type Name: " // trim(model_type_char))
    end select
    
    ! Initialize new object
    call object%initialize()
  
  end function create_flow_model_object
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine check_id_validity(id, fmodel, valid_type)
    implicit none
    ! Checks (A) if flow model type has been set and (B) optionally if it is the correct type of
    ! input for the model type (valid_type)
    ! Used while reading in options
    class(t_flow_model),pointer  :: fmodel
    character(30)                :: id
    integer,optional             :: valid_type
    
    ! Is type set?
    if (ASSOCIATED(fmodel)==.false.) then
      call error_handler(2,opt_msg="Must specify flow model type before other settings")
    else if (present(valid_type)) then
      ! Is it the right type for this setting?
      if (valid_type /= fmodel%model_type) then
        call error_handler(2,opt_msg="Wrong flow model type for option: " // id)
      end if
    end if
    
  end subroutine check_id_validity
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! BASE CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine initialize(this)
    ! Should not be run
    implicit none
    class(t_flow_model)       :: this
    
    this%model_type    = 0
    this%nlayers       = 0
    this%preproc_file  = ''
    this%sim_file      = ''
    this%ppzone_file   = ''
    this%hsu_file      = ''
    this%template_file = ''
    this%output_file   = ''
    this%grd_file      = ''
    this%has_aquitards = .false.
    
  end subroutine initialize
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine read_model(this)
    ! Base class method for reading in [nonexistant] model input file(s)
    ! Should not be run
    implicit none
    class(t_flow_model)       :: this
    
  end subroutine read_model
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine alloc_layers(this, lay2grid)
    implicit none
    class(t_flow_model)       :: this
    integer                   :: lay2grid(this%nlayers)
    integer                   :: i

    ! Allocate grid based on unique grids
    this%ngrids   = maxval(lay2grid)
    
    allocate(this%nnodes  (0:this%nlayers), &
             this%lay2grid(  this%nlayers), &  
             this%grid    (  this%ngrids )  )
             
    this%lay2grid = lay2grid
    
  end subroutine alloc_layers
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine alloc_nodes(this)
    implicit none
    class(t_flow_model)       :: this
    integer                   :: i
    
    allocate(this%elev(maxval(this%nnodes(:)), 0:this%nlayers))
    
    ! Loop over layers creating grid object
    do i=1, this%ngrids
      this%grid(i) = create_grid_object(2, this%nnodes(this%lay2grid(i)))
    end do
    
  end subroutine alloc_nodes
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine write_input_summary(this)
    use m_global, only: log
    implicit none
    class(t_flow_model)     :: this
    
    select type (this)
      class is (t_iwfm_model)
        call log%write_valueline('Number of Model Nodes',    this%nnodes(1))
        call log%write_valueline('Number of Model Elements', this%nelements)
        call log%write_valueline('Number of Aquifer Layers', count(this%lay_exist .and. .not. this%lay_aqtard))
        call log%write_valueline('Number of Active Aquitard Layers', count(this%lay_exist .and. this%lay_aqtard))
      class is (t_flow_fdmodel)
        call log%write_valueline('Number of Model Cells (All Layers)', sum(this%nnodes))
        call log%write_valueline('Number of Model Layers',      this%nlayers)
    end select
  
  end subroutine write_input_summary
  
!-------------------------------------------------------------------------------------------------!
  
  function find_in_grid(this, coords) result(res)
    ! Base class method, should not be run
    implicit none
    class(t_flow_model)     :: this
    real(dp)                :: coords(:)
    integer                 :: res
    
    res = 0
    
  end function find_in_grid

!-------------------------------------------------------------------------------------------------!
  
  subroutine get_elev_at_point(this, coords, point_elevs, locid)
    ! Base class method, should not be run
    implicit none
    class(t_flow_model)     :: this
    real(dp), intent(in)    :: coords(:)
    real(dp), intent(out)   :: point_elevs(this%nlayers)
    integer, optional       :: locid    ! Cell id, element id
    
    point_elevs = 0
  
  end subroutine get_elev_at_point

!-------------------------------------------------------------------------------------------------!
  
  function get_midpoint_depths(this, k) result(depths)
    ! To be overloaded where necessary (should be the same for all models...
    implicit none
    class(t_flow_model)                :: this
    integer,intent(in)                 :: k
    real(dp)                           :: depths(this%nnodes(k))
    
    depths = this%elev(:,0) - ((this%elev(:,k-1) + this%elev(:,k)) / 2.0)
    
  end function get_midpoint_depths
  
!-------------------------------------------------------------------------------------------------!
  
  function get_active_ppzones(this) result(nactive)
    ! To be overloaded where necessary
    implicit none
    class(t_flow_model)     :: this
    integer                 :: nactive
    nactive = this%ppzones%ncat
  end function get_active_ppzones

!-------------------------------------------------------------------------------------------------!
  
  function get_layer_id(this, k) result (id)
    ! To be overloaded where necessary
    implicit none
    class(t_flow_model)     :: this
    integer,intent(in)      :: k
    integer                 :: id
    id = k
    return
  end function get_layer_id

!-------------------------------------------------------------------------------------------------!

  function get_layer_no(this, k) result(no)
    ! To be overloaded where necessary
    implicit none
    class(t_flow_model)     :: this
    integer,intent(in)      :: k
    integer                 :: no
    no = k
    return
  end function get_layer_no
  
!-------------------------------------------------------------------------------------------------!

  subroutine read_category_files(this)
    ! Called after model read
    ! Assumes all layers have same number of nodes
    ! (Could be modified in the future by passing all this%nnodes)
    implicit none
    class(t_flow_model)     :: this
    
    ! Read ppzones file
    call this%ppzones%read_file(this%ppzone_file, this%nlayers, this%nnodes(1), 0)
    
    ! Read hsus file
    call this%hsus%read_file(this%hsu_file, this%nlayers, this%nnodes(1), 1)
    
  end subroutine read_category_files
  
!-------------------------------------------------------------------------------------------------!
  
  function is_active_layer(this, k)
    ! To be overloaded where necessary (IWFM has zero-thickness layers, typically aquitards)
    implicit none
    class(t_flow_model)     :: this
    integer,intent(in)      :: k
    logical                 :: is_active_layer
    
    is_active_layer = .true.
    
  end function is_active_layer

!-------------------------------------------------------------------------------------------------!
  
  function is_aquitard_layer(this, k)
    ! To be overloaded where necessary (IWFM has aquitards)
    implicit none
    class(t_flow_model)     :: this
    integer,intent(in)      :: k
    logical                 :: is_aquitard_layer
    
    is_aquitard_layer = .false.
    
  end function is_aquitard_layer

!-------------------------------------------------------------------------------------------------!
  
  subroutine write_fmodel_input(this, par)
    ! Base class method, should not be run
    use m_datasets, only: t_layerdata
    implicit none
    class(t_flow_model)        :: this
    type(t_layerdata),pointer  :: par(:)
  
  end subroutine write_fmodel_input
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine write_node_files(this, datasets, dat_names)
    use m_vstringlist, only: t_vstringlist, vstrlist_length
    ! Calls methods overloaded by subclasses
    use m_datasets, only: t_layerdata
    implicit none
    class(t_flow_model)        :: this
    type(t_layerdata),pointer  :: datasets(:)  ! (nlayers)
    type(t_vstringlist)        :: dat_names
    character(100)             :: fname, temp
    integer                    :: i
    
    do i=1, vstrlist_length(dat_names)  ! ndat, hopefully
      call item2char(dat_names, i, temp)
      write(fname,'(3a)') 't2p_', trim(temp), '.out'
      call this%write_node_file(fname, datasets, i)
    end do
  
  end subroutine write_node_files
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine write_node_file(this, fname, par, par_id)
    ! Base class method, should not be run
    use m_datasets, only: t_layerdata
    implicit none
    class(t_flow_model)        :: this
    character(*)               :: fname
    type(t_layerdata),pointer  :: par(:)
    integer,intent(in)         :: par_id
  
  end subroutine write_node_file
  
!-------------------------------------------------------------------------------------------------!
  
  !subroutine add_node2vert2node(this, node, nvert)
  !  ! Not currently used
  !  implicit none
  !  class(t_flow_fdmodel)     :: this
  !  integer                   :: i, j, k
  !  integer, intent(in)       :: node, nvert
  !  
  !  do j=1, nvert
  !    ! Increment node-element membership array
  !    this%vertnnode(this%node2vert(j)) = this%vertnnode(this%node2vert(j)) + 1
  !    ! Ensure we're not going over our max
  !    if (this%vertnnode(this%node2vert(j)) > max_node_elem_members) then
  !      write(*,*) 'ERROR - max_vert_node_members exceeded!'     ! a slight lie in variable name
  !      write(*,'(a,i7,a,i3)') 'vert ', this%node2vert(j), 'is in >', max_node_elem_members
  !      stop
  !    end if
  !    ! Otherwise, add to the array
  !    this%vert2node(this%vertnnode(this%node2vert(j)),this%node2vert(j)) = node
  !  end do
  !
  !end subroutine add_node2vert2node
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine read_gsf(this)
    implicit none
    class(t_flow_fdmodel)     :: this
    type(t_vstring)           :: strtemp
    integer                   :: i, j, total_nodes, ngrids, nvert, node, lay, grid, itemp, &
                                 gridcount, nodevert
    integer,allocatable       :: nnodes(:), grid_count(:), lay2grid(:),temp_node_vert(:)
    real(dp)                  :: rtemp(3)
    character(200)            :: ctemp
    logical                   :: found
    
    ! Standard file reader variables
    integer                    :: eof
    type(t_vstringlist)        :: strings
  
    reader => open_file_reader(this%grd_file)
    
    call reader%next_item(eof, strings) ! Line 1 structured or unstructured (unused)
    ! Next couple lines are read if t_mgrid, but not if MODFLOW
    ! Additionally, arrays SHOULD already be allocated for MODFLOW
    select type (this)
      type is (t_mgrid)
        ! Line 2 Nodes for all layers, nlay, iz, ic (used if mgrid, not MF-USG)
        read(reader%unit, *) total_nodes, this%nlayers
        allocate(nnodes(this%nlayers))
        ! Line 3 (Custom T2P) nodes per layer (NODELAY in USG)
        call reader%next_item(eof, strings)
        ngrids = 1
        allocate(grid_count(this%nlayers), lay2grid(this%nlayers))
        if (vstring_match(vstrlist_index(strings,1),'CONSTANT')) then
          nnodes(:) = item2int(strings, 2)              ! same nnodes each layer
          lay2grid = 1
        else                                                 ! The hard way
          do i=1, this%nlayers
            nnodes(i) = item2int(strings, i)
          end do
          grid_count(1) = nnodes(1)
          do i=2, this%nlayers
            found = .false.
            do j=1, ngrids
              if (nnodes(i) == grid_count(j)) then
                found = .true.
                lay2grid(i) = j
                exit
              end if
            end do
            if (found == .false.) then
              ngrids = ngrids + 1
              grid_count(ngrids) = nnodes(i)
              lay2grid(i) = ngrids
            end if
          end do
        end if
      type is (t_modflow_model)
        call reader%skip(2, eof)
    end select  ! Back to the easy stuff
    
    ! Allocate layers, nodes
    this%ngrids  = ngrids
    call this%alloc_layers(lay2grid)
    this%nnodes = nnodes
    call this%alloc_nodes()
    
    ! Line 4 Number of vertices
    read(reader%unit, *) nvert
    this%nvert = nvert
    allocate(this%vertgrid)
    this%vertgrid = create_grid_object(3, this%nvert)
    write(*,*) this%nnodes
        
    ! Read vertices
    do i=1, nvert
      read(reader%unit, *) this%vertgrid%coords(1:3,i)  ! xyz
    end do
    
    ! Read to get number of verticies associated with each node
    allocate(this%nodenvert(total_nodes))
    this%nodenvert = 0
    grid = 1
    gridcount = 0
    do i=1, total_nodes
      read(reader%unit, *) node, rtemp(1:3), lay, nodevert
      this%nodenvert(i) = nodevert
    end do
    
    ! Allocate vertex lists...
    !allocate(this%vertnnode(nvert), &
    !         this%vert2node(max_node_elem_members,nvert), &
    !         this%node2vert(maxval(this%nodenvert(:))))
    !this%vertnnode = 0
    !this%vert2node = 0
    allocate(temp_node_vert(maxval(this%nodenvert(:))), &
             this%node2vert(maxval(this%nodenvert(:)), total_nodes))
    this%node2vert = 0
    
    ! Now we backup... and read for vert2node, node x, y, [and a calculated z if grid, not MF-USG]
    call reader%back(total_nodes)
    grid = 1
    gridcount = 0
    do i=1, total_nodes
      temp_node_vert = 0
      read(reader%unit, *) itemp, rtemp(1:3), lay, nodevert, temp_node_vert(1:this%nodenvert(i))
      node = merge(itemp - sum(this%nnodes(2:lay)), itemp, lay > 1)
      this%node2vert(1:this%nodenvert(i), node) = temp_node_vert(1:this%nodenvert(i))
      !call this%add_node2vert2node(i, nodevert/2, this%node2vert)  ! Only connecting to top vertices
      if (lay2grid(lay) == grid) then
        ! Process node
        this%grid(grid)%coords(1:2,i) = rtemp(1:2)
      end if
      ! Get mean top/bot if grid
      select type (this)
        type is (t_mgrid)
          if (node==1) write(*,*) 'REAL NODE FOUND LAYER', lay, itemp, node
          if (lay==1) then  ! top
            this%elev(node, lay-1) = sum(this%vertgrid%coords(3,temp_node_vert(1:nodevert/2))) / (nodevert/2)
          end if            ! bot
            this%elev(node, lay)   = sum(this%vertgrid%coords(3,temp_node_vert(nodevert/2+1:nodevert  ))) / (nodevert/2)
          ! Make sure top/bot make sense
          if (this%elev(node,lay) > this%elev(node,lay-1)) then
            write(ctemp,*) itemp
            call error_handler(2,opt_msg='Bad vertex elevation or sorting: Layer elevation calc error for node' // ctemp)
          end if
      end select
      if (node==this%nnodes(lay))then
      grid = grid+1                       ! just advance to next grid
      write(*,*) 'GRID FOUND LAYER', lay, grid
      end if
    end do
    
    call reader%close_file()
    
  end subroutine
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! GRID CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

  subroutine init_mgrid(this)
    implicit none
    class(t_mgrid)       :: this
    
    call initialize(this)
    this%model_type  = 1
    this%name = "GSF"
    
  end subroutine init_mgrid  
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine read_mgrid(this)
    implicit none
    class(t_mgrid)       :: this
    
    ! Potentially multiple types of grid files (although unimplemented so far)
    select case (this%grd_file_type)
      case('gsf')
        call this%read_gsf()
      case DEFAULT
        call error_handler(2,opt_msg='Unknown Grid Type for file: ' // this%grd_file)
    end select
  
  end subroutine read_mgrid
  
!-------------------------------------------------------------------------------------------------!
  
  function find_in_grid_gsf(this, coords) result(node)
    use tools, only: pinpol
    ! Assumes coords is in GLOBAL coordinates
    ! Assumes uniform xy layers
    ! Returns 0 if not found
    ! in T2P 1.0, it gave an error within here if the cell was not found. However, seems best to
    ! let the calling routine handle the error - it has more context.
    implicit none
    class(t_mgrid)   :: this
    real(dp)                :: coords(:)
    integer                 :: i,j,node
    integer                 :: nodeids(5)
    real(dp)                :: in_node
    real(dp)                :: nodedist(5), elemx(max_node_elem_members), elemy(max_node_elem_members)
    
    node = 0
    call this%grid(this%lay2grid(1))%get_nnear(coords, 5, nodeids, nodedist)
    
    ! Loop over elements of closest node
    do i=1, 5
      do j=1, this%nodenvert(nodeids(i))/2  ! Only the top vertices
        elemx(j) = this%vertgrid%coords(1, this%node2vert(j, nodeids(i)))
        elemy(j) = this%vertgrid%coords(2, this%node2vert(j, nodeids(i)))
      end do
      ! Check if the point is in this polygon
      in_node = pinpol(coords(1), coords(2), elemx, elemy, this%nodenvert(nodeids(i))/2)
      if (in_node >= 0.0d0) then             ! Greedy, treats on line as in elem
        node = nodeids(i)
        exit
      end if
    end do
    
  end function find_in_grid_gsf
  
!-------------------------------------------------------------------------------------------------!

  
!-------------------------------------------------------------------------------------------------!
! IWFM CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

  subroutine init_iwfm(this)
    implicit none
    class(t_iwfm_model)       :: this
    
    call initialize(this)
    this%model_type  = 2
    this%name = "IWFM"
    
  end subroutine init_iwfm  

!-------------------------------------------------------------------------------------------------!
  
  subroutine alloc_layers_iwfm(this, lay2grid)
    implicit none
    class(t_iwfm_model)       :: this
    integer                   :: i, lay2grid(this%nlayers)
    
    allocate(this%lay_exist (this%nlayers),&
             this%lay_aqtard(this%nlayers) )
    this%lay_exist = .true.
    this%lay_aqtard = .false.
    
    ! Set lay_aqtard to true for every odd layer
    this%lay_aqtard = mod((/(i, i=1, this%nlayers)/), 2) == 1
    
    ! Call super
    call alloc_layers(this, lay2grid)
    
  end subroutine alloc_layers_iwfm
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine alloc_nodes_iwfm(this)
    ! Additional arrays required for IWFM due to Finite-element
    implicit none
    class(t_iwfm_model)       :: this
  
    allocate(this%node2elem(max_node_elem_members,this%nnodes(1)), &
             this%nodenelem(this%nnodes(1)),&
             this%nodeid   (this%nnodes(1)) )
    
    ! Call super
    call alloc_nodes(this)
    
    ! Initialize
    this%node2elem = 0
    this%nodenelem = 0
    
  end subroutine alloc_nodes_iwfm
  
!-------------------------------------------------------------------------------------------------!
  
  function find_in_grid_elem(this, coords) result(element)
    use tools, only: pinpol
    ! Assumes uniform xy layers
    ! Returns 0 if not found
    ! in T2P 1.0, it gave an error within here if the element was not found. However, seems best to
    ! let the calling routine handle the error - it has more context.
    implicit none
    class(t_iwfm_model)     :: this
    real(dp)                :: coords(:)
    integer                 :: element
    integer                 :: nodeids(5)
    real(dp)                :: nodedist(5), elemx(4), elemy(4), inelem
    integer                 :: inode, i, j, count, elem
    
    element = 0
    call this%grid(this%lay2grid(1))%get_nnear(coords, 5, nodeids, nodedist)
    
    ! Loop over elements of closest node
    outerloop:  do inode=1, 5
      do i=1, this%nodenelem(nodeids(inode))
        elem = this%node2elem(i,nodeids(inode))
        count = 0
        ! Get Element node x,y coordinates
        do j=1, 4
          if (this%elements(j, elem) > 0) then
            elemx(j) = this%grid(this%lay2grid(1))%coords(1, this%elements(j, elem))
            elemy(j) = this%grid(this%lay2grid(1))%coords(2, this%elements(j, elem))
            count = count + 1
          end if
        end do
        ! Check if the point is in this polygon
        inelem = pinpol(coords(1), coords(2), elemx, elemy, count)
        if (inelem >= 0.0d0) then             ! Greedy, treats on line as in elem
          element = elem
          exit outerloop
        end if
      end do
    end do outerloop
    
  end function find_in_grid_elem
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine get_elev_at_point_elem(this, coords, point_elevs, locid)
    ! Returns elevations for each layer given coordinates.
    implicit none
    class(t_iwfm_model)   :: this
    real(dp), intent(in)    :: coords(:)
    real(dp), intent(out)   :: point_elevs(0:this%nlayers)
    integer, optional       :: locid      ! Containing cell id if known
    integer                 :: elemid
    
    if (present(locid)) then
      elemid = locid
    else
      elemid = this%find_in_grid(coords)
    end if

    call this%iwfm_elem_elev_idw(coords,elemid,point_elevs)
    
  end subroutine get_elev_at_point_elem

!-------------------------------------------------------------------------------------------------!
  
  subroutine alloc_elems(this)
    implicit none
    class(t_iwfm_model)       :: this
    integer                   :: i
    
    allocate(this%elements(4,this%nelements), & 
             this%elemnnodes(this%nelements), &
             this%elemid    (this%nelements)  )

    this%elements(4,this%nelements) = -1  ! initialize to obvious non-node value
  
  end subroutine alloc_elems

!-------------------------------------------------------------------------------------------------!
  
  subroutine init_nodebyid(this, maxnodeid)
    implicit none
    class(t_iwfm_model)       :: this
    integer,intent(in)        :: maxnodeid
    integer                   :: i
    
    allocate(this%nodebyid(maxnodeid))

    ! Initialize all values to a negative number so we can easily check whether
    ! a nodeid is valid
    this%nodebyid = -1

    ! Fill array
    do i=1, this%nnodes(i)
      this%nodebyid(this%nodeid(i)) = i
    end do
  
  end subroutine init_nodebyid

!-------------------------------------------------------------------------------------------------!
  
  subroutine process_iwfm_layers(this)
    implicit none
    class(t_iwfm_model)       :: this
    integer :: i
    
    ! Calculate maximum thickness (differences) between consecutive layers
    ! Determine if layer exists based on max thickness > 0.0
    this%lay_exist = maxval(abs(this%elev(:, 1:this%nlayers) - this%elev(:, 0:this%nlayers-1)), dim=1) > 0.0_dp
    
    ! See if any of these are aquitard layers
    do i=1, this%nlayers
      if (this%lay_exist(i).and.this%lay_aqtard(i)) this%has_aquitards = .true.
    end do
    
  end subroutine process_iwfm_layers
  
!-------------------------------------------------------------------------------------------------!
  
  function is_active_layer_iwfm(this, k)
    ! To be overloaded where necessary (IWFM has zero-thickness layers, typically aquitards)
    implicit none
    class(t_iwfm_model)     :: this
    integer,intent(in)      :: k
    logical                 :: is_active_layer_iwfm
    
    is_active_layer_iwfm = this%lay_exist(k)
    
  end function is_active_layer_iwfm

!-------------------------------------------------------------------------------------------------!
  
  function is_aquitard_layer_iwfm(this, k)
    ! To be overloaded where necessary (IWFM has aquitards)
    implicit none
    class(t_iwfm_model)     :: this
    integer,intent(in)      :: k
    logical                 :: is_aquitard_layer_iwfm
    
    is_aquitard_layer_iwfm = this%lay_aqtard(k)
    
  end function is_aquitard_layer_iwfm
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine add_elem2nodeelem(this, ielem)
    implicit none
    class(t_iwfm_model)       :: this
    integer                   :: i, j, k, node
    integer, intent(in)       :: ielem

    do j=1, 4
      node = this%elements(j,ielem)
      if (node > 0) then
        ! Increment node-element membership array
        this%nodenelem(node) = this%nodenelem(node) + 1
        ! Ensure we're not going over our max
        if (this%nodenelem(node) > max_node_elem_members) then
          write(*,*) 'ERROR - max_node_elem_members exceeded!'
          write(*,'(a,i7,a,i3)') 'node ', node, 'is in >', max_node_elem_members
          stop
        end if
        ! Otherwise, add to the array
        this%node2elem(this%nodenelem(node),node) = ielem
      end if
    end do

  end subroutine add_elem2nodeelem

!-------------------------------------------------------------------------------------------------!
  
  subroutine read_iwfm(this)
    use m_global, only: dp
    use fpath
    use m_vstringlist, only: t_vstringlist
    implicit none
    type(t_file_reader), pointer :: reader2
    class(t_iwfm_model)          :: this
    character(200)               :: path, ctemp
    real(dp)                     :: factor, rtemp1, rtemp2
    real(dp),allocatable         :: rarray(:)
    integer                      :: i, j, itemp, inode, maxnodeid
    integer, allocatable         :: lay2grid(:)
    
    ! Standard file reader variables
    integer                    :: status, length, eof
    character(200)             :: value
    type(t_vstringlist)        :: strings
    
    ! Get relative path to later be used
    call fpath_strip(this%preproc_file, path, value)
    
    ! Preprocessor file
    reader => open_file_reader(this%preproc_file)
    call reader%read_to_next_iwfm_line(eof)
    call reader%skip(3, eof)                  ! Title
    call reader%read_to_next_iwfm_line(eof)
    call reader%skip(1, eof)                  ! Binary output
    call reader%next_item(eof, strings)       ! Element File
    call item2char(strings, 1, value)
    call fpath_join(path, value, this%elem_file)
    call reader%next_item(eof, strings)       ! Node File
    call item2char(strings, 1, value)
    call fpath_join(path, value, this%node_file)
    call reader%next_item(eof, strings)       ! Stratigraphy File
    call item2char(strings, 1, value)
    call fpath_join(path, value, this%strt_file)
    call reader%close_file()
    
    ! Stratigraphy file (only to nlayers)
    reader2 => open_file_reader(this%strt_file)
    call reader%read_to_next_iwfm_line(eof)
    call reader%next_item(eof, strings)       ! nlayers
    this%nlayers = item2int(strings, 1) * 2   ! Each layer has overlying aquitard (to be treated as layer)
    
    ! Allocate layers
    allocate(lay2grid(this%nlayers))
    lay2grid(:) = 1                           ! all layers use the same xy grid in IWFM
    call this%alloc_layers(lay2grid)
    
    ! Node file
    reader => open_file_reader(this%node_file)
    call reader%read_to_next_iwfm_line(eof)
    call reader%next_item(eof, strings)       ! nnodes
    this%nnodes(:) = item2int(strings, 1)     ! Node count does not vary by layer in IWFM
    call reader%next_item(eof, strings)       ! xy mult factor
    factor = item2dp(strings, 1)
    
    ! Allocate nodes
    call this%alloc_nodes()
    
    ! Continue Node file read
    call reader%read_to_next_iwfm_line(eof)
    maxnodeid = 0
    do i=1,this%nnodes(1)
      ! This marks the spot I realized that my reader class is far too slow for reading long datasets
      read(reader%unit,*) this%nodeid(i), rtemp1, rtemp2
      this%grid(1)%coords(1,i) = rtemp1 * factor
      this%grid(1)%coords(2,i) = rtemp2 * factor
    end do
    call this%init_nodebyid(maxval(this%nodeid)) ! Setup mapping from nodeid->array index
    call reader%close_file()                     ! Close Node file 
    
    ! Continue stratigraphy file read
    call reader2%next_item(eof, strings)         ! elev, thickness mult factor
    factor = item2dp(strings, 1)
    call reader2%read_to_next_iwfm_line(eof)
    allocate(rarray(this%nlayers+1))
    do i=1,this%nnodes(1)
      read(reader2%unit,*) itemp, rarray
      inode = this%nodebyid(itemp)
      this%elev(inode,:) = rarray(1) * factor
      do j=1, this%nlayers                      ! Loop over layers subtracting thicknesses
        this%elev(inode,j:this%nlayers) = this%elev(inode,j:this%nlayers) - (rarray(1+j) * factor)
      end do
    end do
    call reader2%close_file()                   ! Close Stratigraphy file
    call this%process_iwfm_layers()             ! Find layers of zero thickness
    
    ! Elements file
    reader => open_file_reader(this%elem_file)
    call reader%read_to_next_iwfm_line(eof)
    call reader2%next_item(eof, strings)       ! Number of elements
    this%nelements = item2int(strings, 1)
    call reader%read_to_next_iwfm_line(eof)
    call reader2%next_item(eof, strings)       ! Number of regions
    itemp = item2int(strings, 1)
    call reader%read_to_next_iwfm_line(eof)
    call reader%skip(itemp, eof)               ! Skip region descriptions
    
    ! Allocate element arrays
    call this%alloc_elems()
    
    ! Read elements
    call reader%read_to_next_iwfm_line(eof)
    do i=1, this%nelements
      read(reader2%unit,*) this%elemid(i), this%elements(1:4,i)
      this%elemnnodes(i) = count(this%elements(1:4,i) > 0)
  
      ! Translate/Validate nodeid to node index for each connection
      do j=1,4
        if (this%elements(j,i) > 0) then
          if (this%nodebyid(this%elements(j,i)) < 1) then
            write(ctemp,'(i)') this%elemid(i)
            call error_handler(2,opt_msg='Invalid node ID in Element File. Element Number = ' // trim(ctemp))
          end if
          this%elements(j,i) = this%nodebyid(this%elements(j,i))
        end if
      end do
      call this%add_elem2nodeelem(i)
    end do
    call reader%close_file()                   ! Close Element file
    
    ! Simulation file
    reader => open_file_reader(this%sim_file)
    call fpath_strip(this%sim_file, path, value)
    
    call reader%read_to_next_iwfm_line(eof)
    call reader%skip(3, eof)                   ! Title
    call reader%read_to_next_iwfm_line(eof)
    call reader%skip(1, eof)                   ! Binary file
    call reader%next_item(eof, strings)        ! GW File
    call item2char(strings, 1, value)
    call fpath_join(path, value, this%output_file)
    call reader%close_file()                   ! Close Simulation file
    
  end subroutine read_iwfm
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine iwfm_elem_elev_idw(this,coords,ObsElem,InterpValues)
    use m_interpolator, only: dist_2D
  ! Author: Marinko Karanovic, SSP&A
  ! Reworked for T2P v2 by Leland Scantlebury
  ! p hardcoded to 1, hardcoded to layer elevations
    implicit none
    class(t_iwfm_model)   :: this
    integer               :: j,k,nodeID
    integer               :: ObsElem
    real(dp)              :: dist
    real(dp)              :: coords(2),InterpValues(:)
    real(dp)              :: wgt_tmp,wgt
    InterpValues=0
    wgt=0
    do j=1,4
      if (ObsElem.gt.0) then
        nodeID = this%elements(j,ObsElem)
        if (nodeID.gt.0) then
          dist = dist_2D(coords, this%grid(1)%coords(:,nodeID)) ! Assumes all layers share xy
          wgt_tmp = 1.0 / dist
	        wgt = wgt + wgt_tmp
          InterpValues(:) =InterpValues(:) + wgt_tmp * this%elev(nodeID,:)
        end if
      end if
    end do
    InterpValues(:) =InterpValues(:) / wgt
  end subroutine iwfm_elem_elev_idw
  
!-------------------------------------------------------------------------------------------------!
  
  function get_layer_id_iwfm(this, k) result (id)
    implicit none
    class(t_iwfm_model)     :: this
    integer,intent(in)      :: k
    integer                 :: id
    id = ceiling(real(k)/real(2))
    return
  end function get_layer_id_iwfm

!-------------------------------------------------------------------------------------------------!

  function get_layer_no_iwfm(this, k) result(no)
    implicit none
    class(t_iwfm_model)     :: this
    integer,intent(in)      :: k
    integer                 :: no
    no = k*2
    return
  end function get_layer_no_iwfm
    
!-------------------------------------------------------------------------------------------------!

  subroutine write_fmodel_input_iwfm(this, par)
    use m_write_fmodel, only: writeIWFMgwfile
    use m_datasets, only: t_layerdata
    implicit none
    class(t_iwfm_model)        :: this
    type(t_layerdata),pointer  :: par(:)
    
    ! External for now
    call writeIWFMgwfile(this%template_file, this%output_file, par, this%nlayers, this%nnodes(1), this%nodeid, this%lay_aqtard, this%lay_exist)
  
  end subroutine write_fmodel_input_iwfm
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine write_node_file_iwfm(this, fname, par, par_id)
    ! Base class method, should not be run
    use m_datasets, only: t_layerdata
    implicit none
    class(t_iwfm_model)        :: this
    character(*)               :: fname
    type(t_layerdata),pointer  :: par(:)
    integer,intent(in)         :: par_id
    integer                    :: i, k, unit
    character(60)              :: fmt(2)
    type(t_file_writer),pointer:: fwrite

    !write(*,'(4x,a,a35)') 'Writing:', trim(filename)

    write(fmt(1),'("(1a8,2a18,",i5,"i14)")') this%nlayers
    write(fmt(2),'("(1i8,2f18.5,",i5,"e14.5)")') this%nlayers

    !open(20, file=trim(filename), status='replace', action='write', buffered='YES')
    fwrite => open_file_writer(fname)
    unit = fwrite%unit
    
    ! Header
    write(unit,fmt(1)) 'Node', 'X', 'Y', (k, k=1, (this%nlayers/2))
    ! Values
    do i = 1, this%nnodes(1)
      write(unit,fmt(2)) this%nodeid(i), this%grid(1)%coords(1:2,i), (par(k)%values(par_id,i), k=2, this%nlayers, 2)  ! Only writes aquifer params
    end do

    call fwrite%close_file()
  
  end subroutine write_node_file_iwfm
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! FINITE-DIFFERENCE CLASS TYPE-BOUND PROCEDURES (MODFLOW, GSF)
!-------------------------------------------------------------------------------------------------!
  
  subroutine get_elev_at_point_cell(this, coords, point_elevs, locid)
    ! Returns elevations for each layer given coordinates.
    implicit none
    class(t_flow_fdmodel)   :: this
    real(dp), intent(in)    :: coords(:)
    real(dp), intent(out)   :: point_elevs(0:this%nlayers)
    integer, optional       :: locid      ! Containing cell id if known
    integer                 :: cellid
    
    if (present(locid)) then
      cellid = locid
    else
      cellid = this%find_in_grid(coords)
    end if
    
    point_elevs = this%elev(cellid,:)
  
  end subroutine get_elev_at_point_cell
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! MODFLOW 2000 CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine init_modflow2000(this)
    implicit none
    class(t_modflow_model)       :: this
    
    call initialize(this)
    this%model_type  = 3
    this%structured  = .true.
    this%name = "MODFLOW"
    
  end subroutine init_modflow2000
  
!-------------------------------------------------------------------------------------------------!
  subroutine calc_mfcell_centers_global(this)
    use tools, only: loc2glo
    use GLOBAL, ONLY : NROW, NCOL, DELR, DELC
    implicit none
    ! Calculates cell centers in global coordinates  
    class(t_modflow_model)       :: this
    integer                :: i, j, n
    real(dp)               :: xbase, ybase, locx(this%nnodes(1)), locy(this%nnodes(1))
  
    ybase = sum(delc)
    do i=1, nrow
      xbase = 0.0
      do j=1, ncol
        n = (i-1)*ncol+j
        locx(n) = xbase + 0.5 * delr(j)
        locy(n) = ybase - 0.5 * delc(i)
        call loc2glo(locx(n), locy(n), this%xoff, this%yoff, this%rot, this%grid(1)%coords(1,n), this%grid(1)%coords(2,n))
        xbase = xbase + delr(j)
      end do
      ybase = ybase - delc(i)
    end do
  
  end subroutine calc_mfcell_centers_global
!-------------------------------------------------------------------------------------------------!
  
  function find_in_grid_cell(this, coords) result(cell)
    use tools, only: glo2loc
    ! Assumes coords is in GLOBAL coordinates
    ! Assumes uniform xy layers
    ! Returns 0 if not found
    ! in T2P 1.0, it gave an error within here if the cell was not found. However, seems best to
    ! let the calling routine handle the error - it has more context.
    implicit none
    class(t_modflow_model)  :: this
    real(dp)                :: coords(:)
    integer                 :: cellids(4)
    real(dp)                :: celldist(4), xmin,xmax,ymin,ymax,xl,yl
    integer                 :: cell, i
    
    cell = 0
    call this%grid(this%lay2grid(1))%get_nnear(coords, 4, cellids, celldist)  ! do we need more than 4 for USG?
    
    if (this%structured) then
      ! Check each cell center to see if well belongs
      do i=1, 4
        call this%get_cell_local(cellids(i), xmin, xmax, ymin, ymax)
        call glo2loc(coords(1), coords(2), this%xoff, this%yoff, this%rot, xl, yl)
        if ((xl >= xmin).and.(xl < xmax).and.(yl >= ymin).and.(yl < ymax)) then
          cell = cellids(i)
          exit
        end if
      end do
    else  ! USG
      !cell = find_in_grid_gsf(coord) 
    end if
    
  end function find_in_grid_cell
  
  !-------------------------------------------------------------------------------------------------!

  subroutine cell2rowcol(this, cell, row, col)
    use GLOBAL, ONLY : NCOL
    implicit none
    class(t_modflow_model)  :: this
    integer,intent(in)     :: cell
    integer,intent(out)    :: row, col

    row = CEILING(real(cell) / real(NCOL))
    col = cell - (row-1) * NCOL

  end subroutine cell2rowcol
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine get_cell_local(this, cell, xmin, xmax, ymin, ymax)
    ! Assumes uniform xy layers
    use GLOBAL, ONLY : DELR, DELC
    implicit none
    class(t_modflow_model)  :: this
    integer,intent(in)     :: cell
    real(dp),intent(out)       :: xmin,xmax,ymin,ymax
    integer                :: row, col
    
    if (this%structured) then
      call this%cell2rowcol(cell, row, col)
      xmin = this%grid(this%lay2grid(1))%coords(1,cell) - 0.5 * delr(col)
      xmax = this%grid(this%lay2grid(1))%coords(1,cell) + 0.5 * delr(col)
      ymin = this%grid(this%lay2grid(1))%coords(2,cell) - 0.5 * delc(row)
      ymax = this%grid(this%lay2grid(1))%coords(2,cell) + 0.5 * delc(row)
    else
      ! TODO
    end if      

  end subroutine get_cell_local

!-------------------------------------------------------------------------------------------------!
  
  subroutine read_modflow2000(this)
    use m_file_io, only: log_unit
    ! MODFLOW globals to import
    use GLOBAL, ONLY : IUNSTR,NLAY,NROW,NCOL,DELR,DELC,TOP,BOT,LAYCBD,IUNIT,NODLAY
    use PARAMMODULE, ONLY: NMLTAR, NZONAR
    use GWFBCFMODULE, ONLY: VKA, SC1, SC2
    implicit none
    class(t_modflow_model)       :: this
    integer                      :: ierr, ngrids, i, j, k, n, mfmaxunit
    integer,allocatable          :: grid_count(:), lay2grid(:)
    character(256)               :: iomsg
    logical                      :: found
    
    ! MODFLOW-USG does the heavy lifting
    call init_modflow_usg(this%sim_file,ierr,iomsg,mfmaxunit)
    
    ! Error handling for external reader...
    if (ierr /= 0) then
      ! Todo replace with actual log writing function
      write(*,'(a)') trim(iomsg)
      write(log_unit,'(3a)') 'IO Error - ', trim(iomsg)
      error stop
    end if
    
    ! Handle MODFLOW settings
    ! UPW or LPF (or none)
    this%layfile = 0
    if (IUNIT(23) > 0) then
      this%layfile = 23
      this%gwftyp = 'LPF'
    end if
    if (IUNIT(45) > 0) then
      this%layfile = 45
      this%gwftyp = 'UPW'
    end if
    if (this%layfile==0) call error_handler(2,'MODFLOW NAM file must have LPF or UPW package present',this%sim_file)
    
    allocate(lay2grid(nlay))
    lay2grid = 1
    
    ! Check out if/how unstructured the model is
    ngrids = 1
    if (IUNSTR) then
      this%structured = .false.
      ! Equal nodes each layer?
      allocate(grid_count(nlay))
      grid_count(1) = nodlay(1)
      do i=2, nlay
        found = .false.
        n = nodlay(i)-nodlay(i-1)
        do j=1, ngrids
          if (n == grid_count(j)) then
            found = .true.
            lay2grid(i) = j
            exit
          end if
        end do
        if (found == .false.) then
          ngrids = ngrids + 1
          grid_count(ngrids) = n
          lay2grid(i) = ngrids
        end if
      end do
    end if
    
    ! Set Layers
    this%nlayers = nlay
    this%ngrids  = ngrids
    call this%alloc_layers(lay2grid)
    
    ! Count nodes
    this%nnodes = (/ (nodlay(i) - nodlay(i-1), i=1, nlay) /)
    call this%alloc_nodes()
    
    ! Move elevation arrays
    this%elev(:,0) = top(1:this%nnodes(1))
    do k=1, nlay
      this%elev(1:this%nnodes(k),k) = bot(nodlay(k-1)+1: nodlay(k))
    end do
    
    ! Get cell/node centers
    if (IUNSTR) then
      ! TODO
    else
      call this%calc_mfcell_centers_global()
    end if
    
    ! Clean up
    if (allocated(grid_count)) deallocate(grid_count)
    
  end subroutine read_modflow2000
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine write_node_file_mf(this, fname, par, par_id)
    ! Base class method, should not be run
    use m_datasets, only: t_layerdata
    implicit none
    class(t_modflow_model)     :: this
    character(*)               :: fname
    type(t_layerdata),pointer  :: par(:)
    integer,intent(in)         :: par_id
    integer                    :: i, k, row, col, unit
    character(60)              :: fmt(2)
    type(t_file_writer),pointer:: fwrite

    !write(*,'(4x,a,a35)') 'Writing:', trim(filename)
    i = 1
    if (this%structured) i = 2
    write(fmt(1),'("(",i1,"a8,2a18,",i5,"i14)")') i, this%nlayers
    write(fmt(2),'("(",i1,"i8,2f18.5,",i5,"e14.5)")') i, this%nlayers

    !open(20, file=trim(filename), status='replace', action='write', buffered='YES')
    fwrite => open_file_writer(fname)
    unit = fwrite%unit
    ! Header
    if (this%structured) then
      write(unit,fmt(1)) 'Row','Col', 'X', 'Y', (k, k=1, this%nlayers)
    else
      write(unit,fmt(1)) 'Node', 'X', 'Y', (k, k=1, this%nlayers)
    end if
    ! Values
    do i = 1, this%nnodes(1)
      if (this%structured) then
        call this%cell2rowcol(i, row, col)
        write(unit,fmt(2)) row, col, this%grid(1)%coords(1:2,i), (par(k)%values(par_id,i), k=1,this%nlayers)
      else
        write(unit,fmt(2)) i, this%grid(1)%coords(1:2,i), (par(k)%values(par_id,i), k=1,this%nlayers)
      end if
    end do

    call fwrite%close_file()
  
  end subroutine write_node_file_mf
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
  
end module m_flow_model