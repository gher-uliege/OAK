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

program applyobsoper
 use matoper
 use rrsqrt
 use date
 use ufileformat
 use assimilation
 use initfile

  implicit none

  real, pointer      :: x(:)
  type(SparseMatrix) :: H
  type(MemLayout) :: ML
  character(256) :: str
  integer iargc

  if (iargc().ne.1) then
     write(stderr,*) 'Usage: applyoper <init file>'
     stop
  end if


  call getarg(1,str); initfname = str

  call MemoryLayout('Model.',ML)
  allocate(x(ML%effsize))

!  call loadVector('In.value',ML,x)
  call loadSparseMatrix('Operator',ML,ML,H)
!  call saveSparseMatrix('filter1.u',ML,ML,H)
   call loadVector('In.value',ML,x)

  call saveVector('Out.value',ML,H.x.x)
end program

