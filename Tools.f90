module tools
  contains
!-----------------------------------------------------------------------------!
  subroutine loc2glo(locx, locy, xoff, yoff, rot, glox, gloy)
    implicit none

    real,intent(in)        :: locx, locy, xoff, yoff, rot
    real,intent(out)       :: glox, gloy
    real, parameter        :: DEG2RAD = 3.141592654/180.0

    glox = xoff + locx * cos(rot*DEG2RAD) - locy * sin(rot*DEG2RAD)
    gloy = yoff + locy * cos(rot*DEG2RAD) + locx * sin(rot*DEG2RAD)

  end subroutine loc2glo
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
  subroutine glo2loc(glox, gloy, xoff, yoff, rot, locx, locy)
    implicit none

    real,intent(in)        :: glox, gloy, xoff, yoff, rot
    real,intent(out)       :: locx, locy
    real, parameter        :: DEG2RAD = 3.141592654/180.0

    locx = (glox-xoff) * cos(rot*DEG2RAD) + (gloy-yoff) * sin(rot*DEG2RAD)
    locy = (gloy-yoff) * cos(rot*DEG2RAD) - (glox-xoff) * sin(rot*DEG2RAD)

  end subroutine glo2loc
!-----------------------------------------------------------------------------!
  
end module tools 