module m_T2P
  use m_options, only: t_options
!-------------------------------------------------------------------------------------------------!
  type t_T2P
    type(t_options)                  :: opt  ! avoid conflict with Fortran intrinsic 'options'
    !type(t_vstringlist)              :: errors
    !type(t_FlowModel)                :: fmodel
    !type(t_Interpolator)             :: interpolator
    !type(t_Dataset),allocatable()    :: datasets
    !type(t_Variogram),allocatable()  :: variograms
    !type(t_PilotPoint),allocatable() :: pilot_points
  
  contains
  
    procedure :: initialize
    !procedure :: allocate
  
  end type t_T2P
!-------------------------------------------------------------------------------------------------!
  
  contains
  
  subroutine initialize(this)
    implicit none
    class(t_T2P), intent(inout)      :: this
    
    ! Initialize Elements
    call this%opt%initialize()
    
  end subroutine initialize
  
end module m_T2P