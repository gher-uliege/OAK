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
 use ufileformat
 use date
 use initfile
 use assimilation

 implicit none

 real, pointer      :: xf(:), invsqrtR(:), Hshift(:)
 integer            :: iargc,ntime,startntime,endntime
 character(len=124)   :: ntimeindex
 character(len=124) :: str,path
 character(len=256) :: prefix, Oprefix, Dprefix
 integer            :: m,n,k
 type(SparseMatrix) :: H
 type(MemLayout)    :: ObsML

 if (iargc().ne.2.and.iargc().ne.3) then
   write(stderr,*) 'Usage: applyobsoper <init file> <time index> '
   write(stderr,*) 'Usage: applyobsoper <init file> <start time index> <end time index> '
   stop
 end if

 call getarg(1,str); call init(str)
 call getarg(2,str); read(str,*) startntime  
 if (iargc().eq.3) then
   call getarg(3,str); read(str,*) endntime  
 else
   endntime = startntime
 end if

 do ntime=startntime,endntime

!   write(ntimeindex,'(I3.3,A)') ntime,'.'
!   write(Oprefix,'(A,I3.3,A)') 'Obs',ntime,'.'
!   write(Dprefix,'(A,I3.3,A)') 'Diag',ntime,'.'

   call fmtIndex('',ntime,'.',ntimeindex)
   Oprefix = 'Obs'//trim(ntimeindex)
   Dprefix = 'Diag'//trim(ntimeindex)


   allocate(xf(ModMLParallel%effsize))

   call loadVector('Forecast'//trim(ntimeindex)//'value',ModMLParallel,xf)
   call MemoryLayout(Oprefix,ObsML)

   allocate(invsqrtR(ObsML%effsize),Hshift(ObsML%effsize))
   invsqrtR = 1.

   !
   ! load the obervation matrix. All points out of grid will have a 
   ! weight (invsqrtR) = 0
   !

   call loadObservationOper(ntime,ObsML,H,Hshift,invsqrtR)

   call saveVector(trim(Dprefix)//'Hxf',ObsML,(H.x.xf)+Hshift,invsqrtR.ne.0.)

   if (presentInitValue(initfname,trim(Dprefix)//'invsqrtR')) then
     call saveVector(trim(Dprefix)//'invsqrtR',ObsML,invsqrtR)
   end if

   deallocate(xf,H%i,H%j,H%s,invsqrtR,Hshift)
   call MemoryLayoutDone(ObsML)
 end do
end program applyobsoper
