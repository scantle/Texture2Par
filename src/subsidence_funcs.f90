module m_subsidence
  ! Extra routines for subsidence - generally to hide away interbed calculations
  use m_global, only: log, NODATA
  use m_datasets, only: t_layerdata
  use m_interpolator, only: t_interpolator, t_idw
  use m_kriging_base, only: t_kriging_base
  use m_libkriging_interpolator, only: t_kriging_interpolator
  use m_vario, only: t_vgm
  use m_options, only: t_options
  use m_error_handler, only: error_handler
  use m_grid, only: t_grid

  implicit none

  real                             :: interbed_cutoff = 1.0
  class(t_interpolator),pointer    :: interp_ib => null()

  logical                          :: has_ib_variogram = .false.
  logical                          :: has_ib_classes   = .false.

  contains

!-------------------------------------------------------------------------------------------------!

  subroutine build_interbed_mask(nsets, interbed_class_ids, file_class_ids, interbed_mask)
    implicit none
    integer,intent(in)    :: nsets, interbed_class_ids(:), file_class_ids(nsets)
    logical,intent(out)   :: interbed_mask(nsets)
    integer               :: i, j

    interbed_mask = .false.

    do i=1, size(file_class_ids)
      do j=1, size(interbed_class_ids)
        if (file_class_ids(i)==interbed_class_ids(j)) interbed_mask(i) = .true.
      end do
    end do

  end subroutine

!-------------------------------------------------------------------------------------------------!

  subroutine setup_ib_variogram(variograms, nnear, nstruct, opt)
    implicit none
    integer, pointer              :: nnear(:), nstruct(:)
    type(t_vgm),pointer           :: variograms(:,:)
    class(*),pointer              :: p
    type(t_options),intent(in)    :: opt

    ! Handle Interpolation Method
    select case (opt%IB_INTERP_METHOD)
      case("OK")
        allocate(t_kriging_interpolator::interp_ib)
      case("IDW")
        allocate(t_idw::interp_ib)
      case DEFAULT
        call error_handler(1,opt_msg="Invalid IB_INTERP_METHOD: " // opt%IB_INTERP_METHOD)
    end select
    call interp_ib%initialize()

    select type(p => interp_ib)
      class is(t_kriging_base)
        p%nvario     = 1
        p%nstruct    => nstruct
        p%variograms => variograms
        p%nnear      => nnear
        if (allocated(p%global_vario_idx)) deallocate(p%global_vario_idx)
        allocate(p%global_vario_idx(1,1,1))
        p%global_vario_idx = 1
    end select

    has_ib_variogram = .true.
  end subroutine setup_ib_variogram

!-------------------------------------------------------------------------------------------------!

  subroutine check_subsidence_inputs(opt)
    ! Run in main Texture2Par file, currently before datasets are read but after fmodel/categories
    implicit none
    type(t_options),intent(in)    :: opt

    if(.not. has_ib_classes) call error_handler(1,opt_msg="Missing FINE_CLASSES in SUBSIDENCE block")

    if (opt%IB_INTERP_METHOD=="OK".and. &
        .not.has_ib_variogram) then
      call error_handler(1,opt_msg="Missing Required VARIOGRAM argument in SUBSIDENCE block")
    end if

  end subroutine check_subsidence_inputs

!-------------------------------------------------------------------------------------------------!

logical function is_interbed_interval(wvalues, nsets, interbed_mask, inferred_fine)
  implicit none
  real   , intent(in) :: wvalues(:)
  integer, intent(in) :: nsets
  logical, intent(in) :: interbed_mask(:)   ! size >= nsets, marks which texture columns count as interbed
  logical, intent(in) :: inferred_fine      ! when true, count (1 - sum(known)) as a candidate interbed class
  integer :: i
  real :: s_known, s_interbed, inferred

  s_known = 0.0
  s_interbed = 0.0
  do i=1,nsets
    if (wvalues(i) > NODATA) then
      s_known = s_known + wvalues(i)
      if (interbed_mask(i)) s_interbed = s_interbed + wvalues(i)
    end if
  end do
  if (inferred_fine) then
    inferred = max(0.0, 1.0 - s_known)
    s_interbed = s_interbed + inferred
  end if
  is_interbed_interval = (s_interbed >= interbed_cutoff)
end function is_interbed_interval

!-------------------------------------------------------------------------------------------------!

subroutine accumulate_interbeds(local_elevs, intv_top, intv_bot, is_interbed, kbed_open, kbed_run, iw, interbed_count, interbed_thick)
  implicit none
  real   , intent(in)               :: local_elevs(0:)    ! layer tops: 0..nlayers
  real   , intent(in)               :: intv_top, intv_bot ! elevation coords
  logical, intent(in)               :: is_interbed
  logical, intent(inout)            :: kbed_open(:)       ! (nlayers)
  real   , intent(inout)            :: kbed_run(:)        ! (nlayers)
  integer, intent(in)               :: iw                 ! well index (1..npoints)
  type(t_layerdata), intent(inout)  :: interbed_count(:)  ! interbed counts by (layer) with %values(:) of counts by location
  type(t_layerdata), intent(inout)  :: interbed_thick(:)  ! interbed thickness by (layer) with %values(2,:) (thicknesses, squared thicknesses)
  integer                           :: k
  real                              :: frac_thick

  do k=1, size(kbed_open)
    if ((intv_top >= local_elevs(k)) .and. (intv_bot <= local_elevs(k-1))) then
      frac_thick = min(local_elevs(k-1),intv_top) - max(local_elevs(k),intv_bot)
      if (frac_thick > 0.0) then
        ! If it's the first value, we need to initialize
        if (interbed_count(k)%values(1,iw) <= NODATA) then
          interbed_count(k)%values(1,iw)  = 0.0
          interbed_thick(k)%values(:,iw)  = 0.0
        end if
        ! Business as Usual
        if (is_interbed) then
          if (.not. kbed_open(k)) then
            kbed_open(k) = .true.
            interbed_count(k)%values(1,iw) = interbed_count(k)%values(1,iw) + 1
          end if
          kbed_run(k) = kbed_run(k) + frac_thick
        else
          if (kbed_open(k)) then
            interbed_thick(k)%values(1,iw) = interbed_thick(k)%values(1,iw) + kbed_run(k)
            interbed_thick(k)%values(2,iw) = interbed_thick(k)%values(2,iw) + kbed_run(k)**2
            kbed_run(k) = 0.0
            kbed_open(k) = .false.
          end if
        end if
      end if
    end if
  end do
end subroutine accumulate_interbeds

!-------------------------------------------------------------------------------------------------!

subroutine reset_interbeds(iw, kbed_open, kbed_run, interbed_thick)
  implicit none
  integer, intent(in)    :: iw
  logical, intent(inout) :: kbed_open(:)
  real   , intent(inout) :: kbed_run(:)
  type(t_layerdata), intent(inout)  :: interbed_thick(:)
  integer :: k
  do k=1, size(kbed_open)
    if (kbed_open(k)) then
      interbed_thick(k)%values(1,iw) = interbed_thick(k)%values(1,iw) + kbed_run(k)
      interbed_thick(k)%values(2,iw) = interbed_thick(k)%values(2,iw) + kbed_run(k)**2
    end if
    kbed_run(k) = 0.0
    kbed_open(k) = .false.
  end do
end subroutine reset_interbeds

!-------------------------------------------------------------------------------------------------!

end module m_subsidence
