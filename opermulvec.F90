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


program opermul
  use matoper
  use rrsqrt
  use ufileformat
  use date
  use initfile
  use assimilation

  implicit none

  real, pointer      :: x(:)
  type(SparseMatrix) :: H
  type(MemLayout) :: ML1,ML2,ML3
  character(256) :: str
  integer iargc
  logical, pointer :: valid1(:),valid2(:)

  if (iargc().ne.1) then
     write(stderr,*) 'Usage: opermulvec <init file>'
     stop
  end if

  call getarg(1,str); initfname = str

  call MemoryLayout('Space1.',ML1)
  call MemoryLayout('Space2.',ML2)

  allocate(valid1(ML1%effsize),valid2(ML3%effsize),x(ML1%effsize))

  call loadSparseMatrix('Operator',ML2,ML1,H)
  write(6,*) 'H%n,H%m ',H%n,H%m
  call loadVector('In.value',ML1,x)

  call saveVector('Out.value',ML2,H.x.x)
end program
