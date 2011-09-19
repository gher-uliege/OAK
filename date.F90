!
!  OAK, Ocean Assimilation Kit
!  Copyright(c) 2002-2011 Alexander Barth and Luc Vandenblucke
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!

! include the fortran preprocessor definitions
#include "ppdef.h"

module date
contains

function JulianDay(d,m,y)
implicit none
integer, intent(in) :: d,m,y
real(kind=8) :: JulianDay

JulianDay = ( 1461 * ( y + 4800 + ( m - 14 ) / 12 ) ) / 4 +      &
          ( 367 * ( m - 2 - 12 * ( ( m - 14 ) / 12 ) ) ) / 12 -  &
          ( 3 * ( ( y + 4900 + ( m - 14 ) / 12 ) / 100 ) ) / 4 + &
          d - 32075.5_8

end function

!_______________________________________________________
!

function ChronologicalJulianDay(d,m,y)
implicit none
integer, intent(in) :: d,m,y
integer :: ChronologicalJulianDay

! Mathematicians and programmers have naturally 
! interested themselves in mathematical and computational 
! algorithms to convert between Julian day numbers and 
! Gregorian dates. The following conversion algorithm is due 
! to Henry F. Fliegel and Thomas C. Van Flandern: 
! The Julian day (jd) is computed from Gregorian day, month and year (d, m, y) as follows:
! http://hermetic.magnet.ch/cal_stud/jdn.htm

ChronologicalJulianDay = ( 1461 * ( y + 4800 + ( m - 14 ) / 12 ) ) / 4 +      &
          ( 367 * ( m - 2 - 12 * ( ( m - 14 ) / 12 ) ) ) / 12 -  &
          ( 3 * ( ( y + 4900 + ( m - 14 ) / 12 ) / 100 ) ) / 4 + &
          d - 32075

end function

!_______________________________________________________
!

function ModifiedJulianDay(d,m,y)
implicit none
integer, intent(in) :: d,m,y
integer :: ModifiedJulianDay

! Mathematicians and programmers have naturally 
! interested themselves in mathematical and computational 
! algorithms to convert between Julian day numbers and 
! Gregorian dates. The following conversion algorithm is due 
! to Henry F. Fliegel and Thomas C. Van Flandern: 
! The Julian day (jd) is computed from Gregorian day, month and year (d, m, y) as follows:
! http://hermetic.magnet.ch/cal_stud/jdn.htm

! ModifiedJulianDay = 0 for 1858-11-17 CE.

ModifiedJulianDay = ( 1461 * ( y + 4800 + ( m - 14 ) / 12 ) ) / 4 +      &
          ( 367 * ( m - 2 - 12 * ( ( m - 14 ) / 12 ) ) ) / 12 -  &
          ( 3 * ( ( y + 4900 + ( m - 14 ) / 12 ) / 100 ) ) / 4 + &
          d - 32075 - 2400001
end function

!_______________________________________________________
!

subroutine GregorianDate(cjd,d,m,y)
implicit none
integer, intent(in)  :: cjd
integer, intent(out) :: d,m,y

integer              :: l,n,i,j

! Converting from the chronological Julian day number to the Gregorian 
! date is performed thus:


        l = cjd + 68569
        n = ( 4 * l ) / 146097
        l = l - ( 146097 * n + 3 ) / 4
        i = ( 4000 * ( l + 1 ) ) / 1461001
        l = l - ( 1461 * i ) / 4 + 31
        j = ( 80 * l ) / 2447
        d = l - ( 2447 * j ) / 80
        l = j / 11
        m = j + 2 - ( 12 * l )
        y = 100 * ( n - 49 ) + i + l

end subroutine


!_______________________________________________________
!

      function mjd(y,m,d,s)
      implicit none
      integer d,m,y
      real s
      real*8 mjd

! Mathematicians and programmers have naturally 
! interested themselves in mathematical and computational 
! algorithms to convert between Julian day numbers and 
! Gregorian dates. The following conversion algorithm is due 
! to Henry F. Fliegel and Thomas C. Van Flandern: 
! The Julian day (jd) is computed from Gregorian day, month and year (d, m, y) as follows:
! http://hermetic.magnet.ch/cal_stud/jdn.htm

! ModifiedJulianDay = 0 for 1858-11-17 CE.

      mjd = (( 1461 * ( y + 4800 + ( m - 14 ) / 12 ) ) / 4 +        &
            ( 367 * ( m - 2 - 12 * ( ( m - 14 ) / 12 ) ) ) / 12 -   &
            ( 3 * ( ( y + 4900 + ( m - 14 ) / 12 ) / 100 ) ) / 4 +  &
            d - 32075 - 2400001)*1d0 + s/(24*60*60d0)               
      end function mjd


end module
!_______________________________________________________
!
