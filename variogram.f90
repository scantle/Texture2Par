module m_vario
!-------------------------------------------------------------------------------------------------!

type :: t_vgm
! Variogram or cross variogram structure
  integer                  :: itype       ! Variogram model (0-linear, 1-spherical, 2-exponential)
  real                     :: sill        ! Contribution of each structure, when summed (with nugget) the sill
  real                     :: a_hmax      ! Range: distance where data are no longer correlated
  real                     :: a_hmin      ! Range min (equal to amax when anis1/anis2 ==1)
  real                     :: a_vert      ! Vertical range
  real                     :: nugget      ! Nugget: sill immediately at distances > 0
  real                     :: ang1        ! Angle 1 (azimuth correction)
  real                     :: ang2        ! Angle 2 (dip correction)
  real                     :: ang3        ! Angle 3 (plunge correction)
  real                     :: anis1       ! Anisotropy 1 (horizontal)
  real                     :: anis2       ! Anisotropy 2 (vertical)

  contains
    procedure,public        :: initialize
  
end type t_vgm

!-------------------------------------------------------------------------------------------------!

  contains
  
!-------------------------------------------------------------------------------------------------!
! MODULE PROCEDURES
!-------------------------------------------------------------------------------------------------!

  function get_vgm_from_line(strings, struct_id) result(vgm)
    use m_file_io
    use m_vstringlist, only: t_vstringlist
    use m_error_handler, only: error_handler
    implicit none
    
    type(t_vstringlist),intent(in)   :: strings
    integer, intent(in)              :: struct_id
    class(t_vgm),pointer             :: vgm
    
    allocate(vgm)
    
    vgm%itype  = find_string_index_in_list(strings, 2, "LIN,SPH,EXP", toupper=.true.)
    if (vgm%itype==0) then
      call error_handler(1,"Invalid variogram type - must be one of Lin, Sph, Exp")
    else
      ! Convert to zero index (see itype definiton in t_vgm type)
      vgm%itype = vgm%itype - 1
    end if
    
    vgm%nugget = item2real(strings, 3)
    vgm%sill   = item2real(strings, 4)
    vgm%a_hmin = item2real(strings, 5)
    vgm%a_hmax = item2real(strings, 6)
    vgm%a_vert = item2real(strings, 7)
    vgm%ang1   = item2real(strings, 8)
    vgm%ang2   = item2real(strings, 9)
    vgm%ang3   = item2real(strings,10)
    !nnear read outside of routine
    
    ! initialize
    call vgm%initialize(struct_id)
    
  end function get_vgm_from_line
  
!-------------------------------------------------------------------------------------------------!
  
!-------------------------------------------------------------------------------------------------!
! BASE CLASS TYPE-BOUND PROCEDURES
!-------------------------------------------------------------------------------------------------!
  
  subroutine initialize(this, sid)
    implicit none
    class(t_vgm)         :: this
    integer, intent(in)  :: sid
    integer              :: i,j

    ! Sets anis1, anis2 ratios based on ranges aa1, aa2
    this%anis1 = this%a_hmin/this%a_hmax
    this%anis2 = this%a_vert/this%a_hmax
    
    ! Correct for zero a_vert (may not be using 3D kriging)
    if (this%anis2 < tiny(this%anis2)) this%anis2 = 1.0
      
    ! Zero out nuggets on structures above 1
    if (sid > 1) this%nugget = 0.0
    
    ! Now the rotation matrices get initialized
    !TODO call setrot(i, j)
    
  end subroutine initialize

!-------------------------------------------------------------------------------------------------!
  
end module m_vario