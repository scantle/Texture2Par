module m_preproc
  use m_global, only: NODATA
  use m_vstringlist
  use m_vstring, only: t_vstring, vstring_toupper
  use m_file_io, only: t_file_reader, open_file_reader
  use m_error_handler, only: error_handler
  use m_flow_model, only: t_flow_model
  use m_datasets, only: t_layerdata
!-------------------------------------------------------------------------------------------------!
! Convert a single non-texture primary variable to texture after interpolation to the model grid.
! Model-sized arrays are allocatable to keep them off the stack.
!-------------------------------------------------------------------------------------------------!
  implicit none
  private

  public :: t_preprocessor, create_preproc_object

  type :: t_preprocessor
    character(500)                  :: filename = ''
    character(500)                  :: mult_filename = ''
    integer                         :: nclasses = 0
    integer                         :: ntexclasses = 0
    type(t_vstringlist)             :: class_names
    class(t_flow_model),pointer     :: fmodel => null()
    real,allocatable                :: mult(:,:)

    contains
      procedure,private             :: initialize
      procedure,private             :: read_multiplier_file
      procedure,private             :: apply_mult
      procedure,public              :: read_input_files
      procedure,public              :: to_texture
      procedure,public              :: get_class_id
  end type t_preprocessor

  type, extends(t_preprocessor) :: t_res2texture
    real,allocatable                :: shapes(:,:,:)
    real,allocatable                :: scales(:,:,:)

    contains
      procedure,public              :: read_input_files => read_input_files_r2t
      procedure,public              :: to_texture       => to_texture_r2t
      procedure,private             :: calc_normalized_likelihoods
  end type t_res2texture

  contains

!-------------------------------------------------------------------------------------------------!

  function create_preproc_object(preproc_mode_char, ntexclasses, fmodel) result(object)
    implicit none
    character(*),intent(in)                :: preproc_mode_char
    integer,intent(in)                     :: ntexclasses
    class(t_flow_model),pointer,intent(in) :: fmodel
    class(t_preprocessor),pointer          :: object

    select case(trim(preproc_mode_char))
      case('RES2TEXTURE')
        allocate(t_res2texture::object)
      case default
        call error_handler(1,opt_msg='Unknown Preprocessor Mode Name: '//trim(preproc_mode_char))
    end select
    call object%initialize(ntexclasses, fmodel)
  end function create_preproc_object

!-------------------------------------------------------------------------------------------------!

  subroutine initialize(this, ntexclasses, fmodel)
    implicit none
    class(t_preprocessor),intent(inout)     :: this
    integer,intent(in)                     :: ntexclasses
    class(t_flow_model),pointer,intent(in) :: fmodel

    call vstrlist_new(this%class_names)
    this%nclasses = 0
    this%ntexclasses = ntexclasses
    this%filename = ''
    this%mult_filename = ''
    this%fmodel => fmodel
  end subroutine initialize

!-------------------------------------------------------------------------------------------------!

  subroutine read_input_files(this)
    implicit none
    class(t_preprocessor),intent(inout) :: this

    call this%read_multiplier_file()
  end subroutine read_input_files

!-------------------------------------------------------------------------------------------------!

  subroutine read_multiplier_file(this)
    implicit none
    class(t_preprocessor),intent(inout) :: this
    type(t_file_reader),pointer         :: reader
    logical,allocatable                 :: seen(:,:)
    integer                             :: ierr, inode, lay, k
    real                                :: multiplier
    character(500)                      :: header

    if (len_trim(this%mult_filename) == 0) return

    allocate(this%mult(maxval(this%fmodel%nnodes),this%fmodel%nlayers), &
             seen(maxval(this%fmodel%nnodes),this%fmodel%nlayers))
    this%mult = 0.0
    seen = .false.

    reader => open_file_reader(this%mult_filename)
    read(reader%unit,'(A)',iostat=ierr) header
    if (ierr /= 0) call error_handler(3,this%mult_filename,'Unable to read multiplier header')
    do
      read(reader%unit,*,iostat=ierr) inode, lay, multiplier
      if (ierr < 0) exit
      if (ierr > 0) call error_handler(3,this%mult_filename,'Invalid multiplier row')
      if (lay < 1 .or. lay > this%fmodel%nlayers) &
        call error_handler(1,this%mult_filename,'Multiplier layer is outside the model')
      if (inode < 1 .or. inode > this%fmodel%nnodes(lay)) &
        call error_handler(1,this%mult_filename,'Multiplier node is outside the model layer')
      if (seen(inode,lay)) call error_handler(1,this%mult_filename,'Duplicate multiplier node/layer')
      this%mult(inode,lay) = multiplier
      seen(inode,lay) = .true.
    end do
    call reader%close_file()

    do k=1, this%fmodel%nlayers
      if (any(.not.seen(1:this%fmodel%nnodes(k),k))) &
        call error_handler(2,this%mult_filename,'Multiplier file does not cover every model cell')
    end do
    deallocate(seen)
  end subroutine read_multiplier_file

!-------------------------------------------------------------------------------------------------!

  subroutine apply_mult(this, layerdata, class_idx)
    implicit none
    class(t_preprocessor),intent(in) :: this
    type(t_layerdata),intent(inout)  :: layerdata(:)
    integer,intent(in)               :: class_idx
    integer                          :: i, k

    if (.not.allocated(this%mult)) return
    do k=1, this%fmodel%nlayers
      do i=1, this%fmodel%nnodes(k)
        if (.not.this%fmodel%active_nodes_bylay(i,k)) cycle
        if (layerdata(k)%values(class_idx,i) <= NODATA) cycle
        layerdata(k)%values(class_idx,i) = layerdata(k)%values(class_idx,i) + &
          sqrt(max(layerdata(k)%variance(class_idx,i),0.0)) * this%mult(i,k)
      end do
    end do
  end subroutine apply_mult

!-------------------------------------------------------------------------------------------------!

  subroutine to_texture(this, layerdata)
    implicit none
    class(t_preprocessor),intent(inout) :: this
    type(t_layerdata),intent(inout)     :: layerdata(:)

    call this%apply_mult(layerdata, 1)
  end subroutine to_texture

!-------------------------------------------------------------------------------------------------!

  function get_class_id(this, string, adjust_for_infer) result(id)
    implicit none
    class(t_preprocessor),intent(in) :: this
    type(t_vstring),intent(in)       :: string
    logical,optional,intent(in)      :: adjust_for_infer
    integer                          :: id

    id = vstrlist_search(this%class_names, vstring_toupper(string), exact=.true.)
  end function get_class_id

!-------------------------------------------------------------------------------------------------!

  subroutine read_input_files_r2t(this)
    implicit none
    class(t_res2texture),intent(inout) :: this
    type(t_file_reader),pointer        :: reader
    logical,allocatable                :: seen(:,:)
    real,allocatable                   :: scales(:), shapes(:)
    integer                            :: ierr, inode, lay, k
    character(1000)                    :: header

    call read_input_files(this)
    if (len_trim(this%filename) == 0) &
      call error_handler(2,opt_msg='RES2TEXTURE requires SCALE_FILE')

    allocate(this%scales(this%ntexclasses,maxval(this%fmodel%nnodes),this%fmodel%nlayers), &
             this%shapes(this%ntexclasses,maxval(this%fmodel%nnodes),this%fmodel%nlayers), &
             scales(this%ntexclasses), shapes(this%ntexclasses), &
             seen(maxval(this%fmodel%nnodes),this%fmodel%nlayers))
    this%scales = NODATA
    this%shapes = NODATA
    seen = .false.

    reader => open_file_reader(this%filename)
    read(reader%unit,'(A)',iostat=ierr) header
    if (ierr /= 0) call error_handler(3,this%filename,'Unable to read texture-distribution header')
    do
      read(reader%unit,*,iostat=ierr) inode, lay, scales, shapes
      if (ierr < 0) exit
      if (ierr > 0) call error_handler(3,this%filename,'Invalid texture-distribution row')
      if (lay < 1 .or. lay > this%fmodel%nlayers) &
        call error_handler(1,this%filename,'Texture-distribution layer is outside the model')
      if (inode < 1 .or. inode > this%fmodel%nnodes(lay)) &
        call error_handler(1,this%filename,'Texture-distribution node is outside the model layer')
      if (seen(inode,lay)) call error_handler(1,this%filename,'Duplicate texture-distribution node/layer')
      if (any(scales <= 0.0) .or. any(shapes <= 0.0)) &
        call error_handler(1,this%filename,'Lognormal scale and shape values must be positive')
      this%scales(:,inode,lay) = scales
      this%shapes(:,inode,lay) = shapes
      seen(inode,lay) = .true.
    end do
    call reader%close_file()

    do k=1, this%fmodel%nlayers
      if (any(.not.seen(1:this%fmodel%nnodes(k),k))) &
        call error_handler(2,this%filename,'Texture-distribution file does not cover every model cell')
    end do
    deallocate(scales,shapes,seen)
  end subroutine read_input_files_r2t

!-------------------------------------------------------------------------------------------------!

  subroutine to_texture_r2t(this, layerdata)
    implicit none
    class(t_res2texture),intent(inout) :: this
    type(t_layerdata),intent(inout)    :: layerdata(:)
    real,allocatable                   :: texture_lhoods(:)
    integer                            :: i, k

    if (size(layerdata) /= this%fmodel%nlayers) &
      call error_handler(9,opt_msg='Preprocessor layer count does not match flow model')
    allocate(texture_lhoods(this%ntexclasses))
    call this%apply_mult(layerdata, 1)

    do k=1, this%fmodel%nlayers
      if (size(layerdata(k)%values,1) < this%ntexclasses) &
        call error_handler(9,opt_msg='Preprocessor output does not have all texture classes')
      do i=1, this%fmodel%nnodes(k)
        if (.not.this%fmodel%active_nodes_bylay(i,k) .or. &
            layerdata(k)%values(1,i) <= NODATA) then
          texture_lhoods = NODATA
        else
          call this%calc_normalized_likelihoods(exp(layerdata(k)%values(1,i)), &
            this%shapes(:,i,k), this%scales(:,i,k), texture_lhoods)
        end if
        layerdata(k)%values(1:this%ntexclasses,i) = texture_lhoods
      end do
    end do
    deallocate(texture_lhoods)
  end subroutine to_texture_r2t

!-------------------------------------------------------------------------------------------------!

  subroutine calc_normalized_likelihoods(this, value, ln_shape, ln_scale, lhood)
    implicit none
    class(t_res2texture),intent(in) :: this
    real,intent(in)                 :: value
    real,intent(in)                 :: ln_shape(this%ntexclasses), ln_scale(this%ntexclasses)
    real,intent(out)                :: lhood(this%ntexclasses)
    real,parameter                  :: delta=1.0, cutoff=1.0e-10
    integer                         :: i

    lhood = 0.0
    if (value < minval(ln_scale)) then
      lhood(minloc(ln_scale,dim=1)) = 1.0
    else if (value > maxval(ln_scale)) then
      lhood(maxloc(ln_scale,dim=1)) = 1.0
    else
      do i=1, this%ntexclasses
        lhood(i) = lognorm_prob(value, delta, ln_scale(i), ln_shape(i))
        if (lhood(i) < cutoff) lhood(i) = 0.0
      end do
      if (sum(lhood) <= 0.0) &
        call error_handler(9,opt_msg='Probabilities are zero for every texture class')
      lhood = lhood / sum(lhood)
    end if
  end subroutine calc_normalized_likelihoods

!-------------------------------------------------------------------------------------------------!
! These scalar routines preserve the dev_res2par likelihood calculation without placing it in the
! general Tools module.
!-------------------------------------------------------------------------------------------------!

  function lognorm_mu(norm_mu, norm_sigma) result(mu)
    implicit none
    real,intent(in) :: norm_mu, norm_sigma
    real            :: mu

    mu = log(norm_mu**2 / sqrt(norm_sigma**2 + norm_mu**2))
  end function lognorm_mu

!-------------------------------------------------------------------------------------------------!

  subroutine normal_01_cdf(x, cdf)
    implicit none
    real,intent(in)  :: x
    real,intent(out) :: cdf
    real,parameter   :: a1=0.398942280444, a2=0.399903438504, a3=5.75885480458, &
                        a4=29.8213557808, a5=2.62433121679, a6=48.6959930692, &
                        a7=5.92885724438, b0=0.398942280385, b1=3.8052e-08, &
                        b2=1.00000615302, b3=3.98064794e-04, b4=1.98615381364, &
                        b5=0.151679116635, b6=5.29330324926, b7=4.8385912808, &
                        b8=15.1508972451, b9=0.742380924027, b10=30.789933034, &
                        b11=3.99019417011
    real             :: q, y

    if (abs(x) <= 1.28) then
      y = 0.5*x*x
      q = 0.5 - abs(x)*(a1-a2*y/(y+a3-a4/(y+a5+a6/(y+a7))))
    else if (abs(x) <= 12.7) then
      y = 0.5*x*x
      q = exp(-y)*b0/(abs(x)-b1+b2/(abs(x)+b3+b4/(abs(x)-b5+ &
        b6/(abs(x)+b7-b8/(abs(x)+b9+b10/(abs(x)+b11))))))
    else
      q = 0.0
    end if
    if (x < 0.0) then
      cdf = q
    else
      cdf = 1.0-q
    end if
  end subroutine normal_01_cdf

!-------------------------------------------------------------------------------------------------!

  subroutine log_normal_cdf(x, mu, sigma, cdf)
    implicit none
    real,intent(in)  :: x, mu, sigma
    real,intent(out) :: cdf

    if (x <= 0.0) then
      cdf = 0.0
    else
      call normal_01_cdf((log(x)-lognorm_mu(mu,sigma))/sigma,cdf)
    end if
  end subroutine log_normal_cdf

!-------------------------------------------------------------------------------------------------!

  function lognorm_prob(x, delta, mu, sigma) result(prob)
    implicit none
    real,intent(in) :: x, delta, mu, sigma
    real            :: right, left, prob

    call log_normal_cdf(x+delta,mu,sigma,right)
    call log_normal_cdf(x-delta,mu,sigma,left)
    prob = right-left
  end function lognorm_prob

!-------------------------------------------------------------------------------------------------!

end module m_preproc
