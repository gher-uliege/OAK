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


program anamorph
  use matoper
  use rrsqrt
  use ufileformat
  use date
  use initfile
  use assimilation

  implicit none

  real, pointer      :: v1(:),v2(:),w(:)
  character(256) :: str
  integer iargc,TransformDir

  if (iargc().ne.1) then
     write(stderr,*) 'Usage: anamorph <init file>'
     stop
  end if

  call getarg(1,str); initfname = str

  call init(initfname)

!  call MemoryLayout('Model.',ModML)
  call getInitValue(initfname,'TransformDir',TransformDir)

  allocate(v1(ModML%effsize))

  call loadVector('Vec1.value',ModML,v1)

  if (TransformDir.eq.0) then
    call anamorphosisTransform(v1)
  else
    call invanamorphosisTransform(v1)
  end if
  call saveVector('Vec2.value',ModML,v1)

  deallocate(v1)
end program
