module m_global
  use m_T2P, only: t_T2P
  
  type(t_T2P),save         :: t2p   ! main Texture2Par storage
  
end module m_global