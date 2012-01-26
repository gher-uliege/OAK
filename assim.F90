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

program assimtest
 use matoper
 use rrsqrt
 use ufileformat
 use initfile
 use assimilation
 use parall

 implicit none

 real, allocatable  :: xf(:), xa(:), &
      Sf(:,:), Sa(:,:)

 integer            :: iargc, ntime, startntime, endntime
 character(len=124) :: str
 character(len=4) :: ntimeindex

 if (iargc().ne.2.and.iargc().ne.3) then
   write(stderr,*) 'Usage: assim <init file> <time index> '
   write(stderr,*) 'Usage: assim <init file> <start time index> <end time index> '
   ERROR_STOP
 end if

#ifdef ASSIM_PARALLEL
   call parallInit()
#endif

 call getarg(1,str); call init(str)
 call getarg(2,str); read(str,*) startntime  
 if (iargc().eq.3) then
   call getarg(3,str); read(str,*) endntime  
 else
   endntime = startntime
 end if

 allocate(xf(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel), &
      xa(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel), &
      Sa(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel,ErrorSpaceDim), &
      Sf(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel,ErrorSpaceDim))

 call loadErrorSpace('ErrorSpace.init',Sf,xf)

 do ntime=startntime,endntime
   write(ntimeindex,'(I3.3,A)') ntime,'.'

   if (presentInitValue(initfname,'Forecast'//ntimeindex//'value')) then
     write(stdlog,*) 'Use forecast xf and ignore ensemble mean'
     call loadStateVector('Forecast'//ntimeindex//'value',xf)
   else
     write(stdlog,*) 'Use ensemble mean for xf'
   end if

!$omp parallel

!$omp critical (writeStdout)
   !   write(6,*) 'sum(xf) ',i,sum(xf),omp_get_thread_num()
!$omp end critical (writeStdout)


   call assim(ntime,xf,Sf,xa,Sa)
!$omp end parallel


   if (presentInitValue(initfname,'Analysis'//ntimeindex//'value'))  &
        call saveStateVector('Analysis'//ntimeindex//'value',xa)

 end do

 call done()

#ifdef ASSIM_PARALLEL
 call parallDone()
#endif


 deallocate(xf,xa,Sa,Sf)
end program assimtest

