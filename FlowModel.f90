module m_flow_model
  use m_file_io
  use m_grid, only: t_grid, create_grid_object
  use m_error_handler, only: error_handler
  implicit none
!-------------------------------------------------------------------------------------------------!
  ! Module variables
  type(t_file_reader), pointer :: reader          ! For reading input files
  integer,parameter         :: max_node_elem_members = 30
!-------------------------------------------------------------------------------------------------!
  
  type :: t_flow_model  ! Abstract Base Class
  ! Base class contains some derived type variables (e.g., the IWFM preproc file) since otherwise
  ! the compiler gets confused when you try to set them on the polymorphic t2p%fmodel type
  ! The alternative was to write a lot of get/set-type methods that would have to be overwritten
  ! by each derived class.
  
    integer                       :: model_type     ! 1 - grid, 2 - IWFM, 3 - MODFLOW 2000
    integer                       :: nlayers        ! Number of model layers
    integer                       :: ngrids         ! Number of unique model grid (xy) layers
    integer,allocatable           :: nnodes(:)      ! Number of nodes (points at which flow is solved) by layer
    integer,pointer               :: lay2grid(:)    ! Which grid each layer uses (e.g., same xy between layers)
    character(100)                :: preproc_file   ! IWFM main discretization file (empty for others)
    character(100)                :: sim_file       ! Simulation File or MODFLOW name file
    character(100)                :: ppzone_file    ! File of pilot point zones defined for each node
    character(100)                :: hsu_file       ! File of Hydrostratigraphic Units (HSUs) defined for each node
    character(100)                :: template_file  ! File copied for output_file
    character(100)                :: output_file    ! File where aquifer parameters are written
    real(dp), allocatable         :: elev(:,:)      ! Elevations for each layer (node, layer)
    type(t_grid),allocatable :: grid(:)        ! xy grids for layers
  
    contains
      procedure,public            :: initialize
      procedure,public            :: read_model
      procedure,private           :: alloc_layers
      procedure,private           :: alloc_nodes
      !procedure,private     :: write_model_input
  
  end type t_flow_model
!-------------------------------------------------------------------------------------------------!
  type, extends(t_flow_model) :: t_mgrid
  
  contains
    procedure,public          :: initialize => init_mgrid
  
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
    integer,allocatable       :: node2elem(:,:) ! Elements that include each node
    logical,allocatable       :: lay_exist(:)   ! Tracks whether layer has thickness > 0.00
    logical,allocatable       :: lay_aqtard(:)  ! Tracks whether layer is an aquitard
  
  contains
    procedure,public          :: initialize => init_iwfm
    procedure,public          :: read_model => read_iwfm
    procedure,private         :: alloc_layers => alloc_layers_iwfm
    procedure,private         :: alloc_nodes  => alloc_nodes_iwfm
    procedure,private         :: alloc_elems
    procedure,private         :: init_nodebyid
    procedure,private         :: process_iwfm_layers
    procedure,private         :: add_elem2nodeelem
  
  end type 
!-------------------------------------------------------------------------------------------------!
  type, extends(t_flow_model) :: t_modflow_model
  
  logical                     :: structured    ! T/F Structured (row/col) model
  
  contains
    procedure,public          :: initialize => init_modflow2000
  
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
    
    !select type (fmodel)
    !  type is (t_iwfm_model)
    !    write(*,*) 'IWFM MODEL'
    !  type is (t_flow_model)
    !    write(*,*) 'GENERIC!'
    !end select
    
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
    this%preproc_file  = ''
    this%sim_file      = ''
    this%ppzone_file   = ''
    this%hsu_file      = ''
    this%template_file = ''
    this%output_file   = ''
    
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
  
!-------------------------------------------------------------------------------------------------!
! GRID CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

  subroutine init_mgrid(this)
    implicit none
    class(t_mgrid)       :: this
    
    call initialize(this)
    this%model_type  = 1
    
  end subroutine init_mgrid  
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! IWFM CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!

  subroutine init_iwfm(this)
    implicit none
    class(t_iwfm_model)       :: this
    
    call initialize(this)
    this%model_type  = 2
    
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
    
  end subroutine process_iwfm_layers
  
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
  
!-------------------------------------------------------------------------------------------------!
! MODFLOW 2000 CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine init_modflow2000(this)
    implicit none
    class(t_modflow_model)       :: this
    
    call initialize(this)
    this%model_type  = 3
    this%structured  = .true.
    
  end subroutine init_modflow2000
  
!-------------------------------------------------------------------------------------------------!
end module m_flow_model