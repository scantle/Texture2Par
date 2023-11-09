module m_global
  use m_file_io, only: t_file_writer
  
  INTEGER, PARAMETER               :: dp = SELECTED_REAL_KIND(15,307)  ! 64-bit
  type(t_file_writer), pointer     :: log   ! Runtime Log
  
end module m_global