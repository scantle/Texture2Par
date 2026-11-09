module m_kriging_base
  use m_global, only: log
  use m_interpolator, only: t_interpolator
  use m_vario, only: t_vgm
  use m_options, only: t_options
  use m_error_handler, only: error_handler
  use m_vstringlist
  implicit none

  type, abstract, extends(t_interpolator) :: t_kriging_base
    integer                   :: nvario
    integer                   :: nzones
    integer,pointer           :: nstruct(:)
    integer,pointer           :: nnear(:)
    type(t_vstringlist)       :: zones
    integer,allocatable       :: cat2zone(:)
    type(t_vgm),pointer       :: variograms(:,:)
    integer,allocatable       :: global_vario_idx(:,:,:)
  contains
    procedure :: add_zonal_variogram
    procedure :: init_zones
    procedure :: write_summary => write_summary_kriging_base
    procedure :: finalize => finalize_kriging_base
  end type t_kriging_base

contains

  function add_zonal_variogram(this, zone_name, class1, class2, opt, glo_nvario) result(zone_id)
    use m_vstring, only: t_vstring
    class(t_kriging_base), intent(inout) :: this
    type(t_vstring), intent(in)          :: zone_name
    integer, intent(in)                  :: class1, class2
    integer, intent(inout)               :: glo_nvario
    type(t_options), intent(in)          :: opt
    integer                              :: zone_id
    integer                              :: i, j
    integer                              :: vgm_idx, old_nzone
    integer, allocatable                 :: tmp(:,:,:)

    if (.not.vstrlist_exists(this%zones)) then
      call vstrlist_new(this%zones)
      this%nzones = 0
    end if

    i = min(class1, class2)
    j = max(class1, class2)
    zone_id = vstrlist_search(this%zones,zone_name)

    if (zone_id == VSTRINGLIST_INDEX_UNKNOWN) then
      call vstrlist_append(this%zones, zone_name)
      zone_id = vstrlist_length(this%zones)

      old_nzone = this%nzones
      this%nzones = zone_id

      allocate(tmp(0:this%nzones, size(this%global_vario_idx,2), size(this%global_vario_idx,3)))
      tmp = -1
      tmp(0:old_nzone,:,:) = this%global_vario_idx(0:old_nzone,:,:)
      tmp(zone_id,:,:) = tmp(0,:,:)
      call move_alloc(tmp, this%global_vario_idx)
    end if

    glo_nvario = glo_nvario + 1
    if (opt%ONLY_TEXTURE) then
      vgm_idx = glo_nvario
    else
      vgm_idx = glo_nvario - 1
    end if

    this%global_vario_idx(zone_id, i, j) = vgm_idx
    this%global_vario_idx(zone_id, j, i) = vgm_idx
  end function add_zonal_variogram

  subroutine init_zones(this, glo_zone_names)
    use m_file_io, only: item2char
    use m_vstring, only: vstring_toupper
    class(t_kriging_base), intent(inout) :: this
    type(t_vstringlist), intent(in)       :: glo_zone_names
    integer                               :: i, cat_id
    character(32)                         :: badzone

    if (.not. vstrlist_exists(this%zones)) return

    do i=1, vstrlist_length(this%zones)
      cat_id = vstrlist_search(glo_zone_names, vstring_toupper(vstrlist_index(this%zones,i)))
      if (cat_id==VSTRINGLIST_INDEX_UNKNOWN) then
        call item2char(this%zones, i, badzone)
        call error_handler(1,opt_msg="Invalid HSU Variogram Used: "//badzone)
      else
        this%cat2zone(cat_id) = i
      end if
    end do
  end subroutine init_zones

  subroutine write_summary_kriging_base(this)
    class(t_kriging_base), intent(inout) :: this
    integer                              :: i, j
    character(len=512)                   :: line, struct
    integer, dimension(3)                :: idx

    do i=1, this%nvario
      idx = findloc(this%global_vario_idx, i)
      write(line, '(A9,I2,A6,I2,A13,I2,A1,I2,A1)') &
        " * Vario ", i, ": HSU=", idx(1), " Class Pair=(", idx(2), ",", idx(3), ")"
      do j=1, this%nstruct(i)
        if (j==1) then
          write(struct, '(i5,2x,a)') j, trim(adjustl(this%variograms(j,i)%write_variogram()))
          call log%write_line(trim(line) // ' - ' // trim(struct))
        else
          write(struct, '(41x,i5,2x,a)') j, this%variograms(j,i)%write_variogram()
          call log%write_line(trim(struct))
        end if
      end do
    end do
  end subroutine write_summary_kriging_base

  subroutine finalize_kriging_base(this)
    class(t_kriging_base), intent(inout) :: this

    if (associated(this%nnear)) deallocate(this%nnear)
    if (associated(this%nstruct)) deallocate(this%nstruct)
    if (associated(this%variograms)) deallocate(this%variograms)
    if (allocated(this%global_vario_idx)) deallocate(this%global_vario_idx)
  end subroutine finalize_kriging_base

end module m_kriging_base
