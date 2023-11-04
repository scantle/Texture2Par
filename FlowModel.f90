module m_flow_model
  use m_file_io, only: t_file_reader, t_file_writer
  use m_grid, only: t_grid
  implicit none
!-------------------------------------------------------------------------------------------------!
  type t_flow_model  ! Abstract Base Class
  
    integer                       :: nnodes        ! Number of nodes (points at which flow is solved)
    integer                       :: nlayers       ! Number of model layers
    type(t_file_reader), pointer  :: template_file ! File copied for output_file
    type(t_file_writer), pointer  :: output_file   ! File where aquifer parameters are written
    type(t_Grid), allocatable     :: layer(:)      ! Grid for ground surface + each layer bottom
  
    contains
      !procedure,private     :: write_model_input
      !procedure,private     ::
  
  end type t_flow_model
!-------------------------------------------------------------------------------------------------!
end module m_flow_model