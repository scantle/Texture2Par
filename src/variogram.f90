module m_vario

! Variogram configuration read from Texture2Par input. Covariance evaluation,
! anisotropy transforms, and kriging are implemented by libkriging.
type :: t_vgm
  integer :: itype       ! Variogram model identifier
  integer :: ndim        ! Number of dimensions
  real    :: sill        ! Partial sill for this structure
  real    :: a_hmax      ! Major horizontal range
  real    :: a_hmin      ! Minor horizontal range
  real    :: a_vert      ! Vertical range
  real    :: nugget      ! Nugget contribution
  real    :: ang1        ! Azimuth
  real    :: ang2        ! Dip
  real    :: ang3        ! Plunge
contains
  procedure, public :: write_variogram
end type t_vgm

contains

  function write_variogram(this) result(line)
    implicit none
    class(t_vgm), intent(in) :: this
    character(len=256)       :: line

    character(len=3), dimension(7) :: model_names
    character(len=3) :: model
    integer :: id

    model_names = ["SPH", "EXP", "HOL", "GAU", "POW", "CIR", "LIN"]
    id = this%itype

    if (id >= 1 .and. id <= size(model_names)) then
      model = model_names(id)
    else
      model = "UNK"
    end if

    if (this%ndim == 2) then
      write(line, '(A,1X,F7.4,1X,F7.4,1X,ES12.4,1X,ES12.4,1X,F8.2)') &
          trim(model), this%nugget, this%sill, this%a_hmin, this%a_hmax, this%ang1
    else if (this%ndim == 3) then
      write(line, '(A,1X,F7.4,1X,F7.4,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,F8.2,1X,F8.2,1X,F8.2)') &
          trim(model), this%nugget, this%sill, this%a_hmin, this%a_hmax, this%a_vert, &
          this%ang1, this%ang2, this%ang3
    end if
  end function write_variogram

end module m_vario