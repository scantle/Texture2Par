module m_libkriging_test
  use kriging, only: t_kriging
  implicit none
contains

  subroutine libkriging_available(ok)
    logical, intent(out) :: ok
    type(t_kriging) :: k

    ok = .true.
  end subroutine libkriging_available

end module m_libkriging_test