module tools
  use m_global, only: dp
  contains
!-----------------------------------------------------------------------------!
  subroutine loc2glo(locx, locy, xoff, yoff, rot, glox, gloy)
    implicit none

    real(dp),intent(in)        :: locx, locy, xoff, yoff, rot
    real(dp),intent(out)       :: glox, gloy
    real(dp), parameter        :: DEG2RAD = 3.141592654/180.0

    glox = xoff + locx * cos(rot*DEG2RAD) - locy * sin(rot*DEG2RAD)
    gloy = yoff + locy * cos(rot*DEG2RAD) + locx * sin(rot*DEG2RAD)

  end subroutine loc2glo
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
  subroutine glo2loc(glox, gloy, xoff, yoff, rot, locx, locy)
    implicit none

    real(dp),intent(in)        :: glox, gloy, xoff, yoff, rot
    real(dp),intent(out)       :: locx, locy
    real(dp), parameter        :: DEG2RAD = 3.141592654/180.0

    locx = (glox-xoff) * cos(rot*DEG2RAD) + (gloy-yoff) * sin(rot*DEG2RAD)
    locy = (gloy-yoff) * cos(rot*DEG2RAD) - (glox-xoff) * sin(rot*DEG2RAD)

  end subroutine glo2loc
!-----------------------------------------------------------------------------!
  
!-----------------------------------------------------------------------------!
function pinpol(px, py, polyx, polyy, npoly) result(inpoly)
  implicit none
  !---------------------------------------------------------------------------!
  !                                 PINPOL
  !                               *********
  !                 Checks if a point is inside a polygon.
  !      "An improved version of the algorithm of Nordbeck and Rydstedt"
  !
  ! Arguments:
  !   -- px is the x-coordinate of the point to be tested
  !   -- py is the y-coordinate of the point to be tested
  !   -- polyx is a vector of nodal x-coordinates in counter-clockwise order
  !   -- polyy is a vector of nodal y-coordinates in counter-clockwise order
  !   -- npoly is the number of polygon verticies
  !
  ! Requires:
  !   -- Function tridet to calculate triangle determinant
  !
  ! Returns:
  !   -- A real value indicating the point's position relative to the polygon:
  !      - Exact Zero (0.0): The point lies exactly on the boundary of the polygon.
  !      - Positive Value: The point is outside the polygon.
  !      - Negative Value: The point is inside the polygon.
  !
  ! Adapted from Sloan, S. W., 1985 A point-in-polygon program, Adv. Eng.
  !                  Software, Vol 7., No. 1
  !
  !               Implemented in f90 by Leland Scantlebury
  !
  ! Copyright 2018-2024 Leland Scantlebury
  !---------------------------------------------------------------------------!

  integer, intent(in)      :: npoly
  integer                  :: i, next, prev, j
  real(dp), intent(in)     :: px, py, polyx(npoly), polyy(npoly)
  real(dp)                 :: x1, y1, x21, y21, x1p, y1p, t, d, dx, dy, area
  real(dp), parameter      :: smalld = 1e-6
  logical                  :: snear
  real(dp)                 :: inpoly

  inpoly = HUGE(1.0)

  ! Loop over each side defining the polygon
  do i=1, npoly
    ! Side coordinates, length, distance from x,y
    next = i + 1
    if (i == npoly) next = 1
    x1 = polyx(i)
    y1 = polyy(i)
    x21 = polyx(next)-x1
    y21 = polyy(next)-y1
    x1p = x1-px
    y1p = y1-py

    ! Find where normal of px,py intersects infinite line
    t = -(x1p*x21 + y1p*y21)/(x21**2 + y21**2)
    if (t < 0.0d0) then
      ! Normal does not intersect side, point is closer to (x1, y1)
      ! Compute square distance to vertex
      d = x1p**2 + y1p**2
      if (d < inpoly) then
        ! Smallest distance yet
        snear = .false.
        inpoly = d
        j = i
      end if
    else if (t < 1.0d0) then
      ! Normal intersects the side
      dx = x1p + t * x21
      dy = y1p + t * y21
      d = dx**2 + dy**2
      if (d < inpoly) then
        ! Smallest distance yet
        snear = .true.
        inpoly = d
        j = i
      end if
    else
      ! Point is closer to the next vertex, continue on to next side
      cycle
    end if
  end do

  if (inpoly < smalld) then
    ! Point lies on the side of the polygon
    inpoly = 0.0d0
  else
    next = j + 1
    prev = j - 1
    if (j == 1) then
      prev = npoly
    else if (j == npoly) then
      next = 1
    else
      continue
    end if
    if (snear) then
      ! Point is closer to side, check if right or left of polygon
      ! (using determinant) if left, point is in polygon
      area = tridet(polyx(j),polyy(j),polyx(next),polyy(next),px,py)
      inpoly = sign(inpoly, area)
    else
      ! Point is closer to node. Check if nearest vertex is concave
      ! If concave, point is inside the polygon
      !if (j == 1) j = npoly + 1
      area = tridet(polyx(next),polyy(next),polyx(j),polyy(j),polyx(prev),polyy(prev))
      inpoly = sign(inpoly, area)
    end if
  end if

  return

  end function pinpol

!-----------------------------------------------------------------------------!

function tridet(x1, y1, x2, y2, x3, y3) result(det)
  implicit none
!-----------------------------------------------------------------------------!
! Computes twice the area of the triangle defined by coordinates
! (x1,y1) (x2,y2) (x3,y3) using determinate formula
!
! If the area is positive, the points are counter-clockwise
! If the area is negative, the points are clockwise
! If the area is zero, two or more of the points are co-located or all three
!    points are collinear
!
! Useful for determining what side of a line a point (x3, y3) is on.
!-----------------------------------------------------------------------------!
  real(dp), intent(in)     :: x1, y1, x2, y2, x3, y3
  real(dp)                 :: det

  det = (x1-x3)*(y2-y3)-(x2-x3)*(y1-y3)

  end function tridet
!-----------------------------------------------------------------------------!
  
end module tools 