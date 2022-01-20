
!-----------------------------------------------------------------------------!
module utl

  contains

  function binarySearch(sorted, x, outside_ok, ascending) result(iloc)
    ! search the location to insert x in a sorted (descending) array
    integer               :: imid, ileft, iright, itmp
    real, intent(in)      :: sorted(:), x
    real, allocatable     :: array(:)
    real                  :: xx
    integer               :: iloc
    logical, optional     :: outside_ok, ascending
    logical               :: found, outsideok, descending
    found = .false.

    if (present(ascending)) then
      descending = .not. ascending
    else
      descending = .true.
    end if

    if (present(outside_ok)) then
      outsideok = outside_ok
    else
      outsideok = .false.
    end if

    if (descending) then
      xx = -x
      array = -sorted
    else
      xx = x
      array = sorted
    end if

    ileft =  1
    iright = size(array)

    iloc = -1
    if (xx < array(ileft) .or. xx > array(iright)) then
      if (outsideok) then
        if (xx < array(ileft)) then
          iloc = 1
        else
          iloc = iright
        end if
      else
        print *, 'binarySearch error: x outside of sorted array'
      end if
      return
    end if

    do while( iright - ileft > 1)
      imid = int((iright + ileft) /2)
      if (xx > array(imid)) then
        ileft = imid
      else
        iright = imid
      end if
    end do
    iloc = ileft

  end function


  subroutine ReadToSymbol(unitno, symchar, kerr)
    use ErrorHandle
    implicit none

    ! Reads file until a specified symbol is found at the start of a line,
    ! Checks that next line doesn't have symbol, then backspaces
    ! Effictively giving you the next "non-commented" line for reading
    ! Useful for jumping past commented lines
    ! LS 7/10/2018

    integer, intent(in)      :: unitno
    integer, optional        :: kerr
    integer                  :: ierr
    character(1), intent(in) :: symchar
    character(5)             :: line

    ! Find Symbol
    do
      read(unitno, *, iostat=ierr) line
      if (present(kerr)) then
          kerr = ierr
      else
        call iostathandler(ierr)
      end if
      if (index(line,symchar) > 0) exit
    end do

    ! Make sure following lines don't have symbol
    do
      read(unitno, *, iostat=ierr) line
      if (present(kerr)) kerr = ierr
      if (ierr /= 0) exit
      if (index(line,symchar) == 0) then
        ! Back one line
        backspace(unitno)
        exit
      end if
    end do
    ! End of File is NOT handled in this subroutine


  end subroutine ReadToSymbol
end module utl
