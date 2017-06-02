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
      Sf(:,:), Sa(:,:), mE(:)

 integer            :: iargc, ntime, enstype, k
 character(len=124) :: str
 character(len=124) :: ntimeindex
 integer            :: ErrorSpaceDim


 if (iargc().ne.2) then
   write(stderr,*) 'Usage: assim <init file> <time index> '
   ERROR_STOP
 end if

#ifdef ASSIM_PARALLEL
   call parallInit()
#endif

 call getarg(1,str); call init(str)
 call getarg(2,str); read(str,*) ntime  
 call getInitValue(initfname,'ErrorSpace.dimension',ErrorSpaceDim,default=0) 


 allocate(xf(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel), &
      xa(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel), &
      Sa(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel,ErrorSpaceDim), &
      Sf(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel,ErrorSpaceDim))

 !write(ntimeindex,'(I3.3,A)') ntime,'.'
 call fmtIndex('',ntime,'.',ntimeindex)

 call getInitValue(initfname,'ErrorSpace.type',enstype,default=1)
 if (enstype == 1) then
   ! Sf are error modes and xf is the forecast
   call loadVectorSpace('ErrorSpace.init',ModMLParallel,Sf)
   call loadVector('Forecast'//trim(ntimeindex)//'value',ModMLParallel,xf)
!$omp parallel
   call assim(ntime,Sf,Sa,xf,xa)
!$omp end parallel
 elseif (enstype == 2) then
   ! Sf are ensemble members
   call loadEnsemble('ErrorSpace.init',ModMLParallel,Sf)

!$omp parallel
   call assim(ntime,Sf,Sa)   
!$omp end parallel

 elseif (enstype == 3) then
   ! Sf are ensemble members to be recentred
   call loadEnsemble('ErrorSpace.init',ModMLParallel,Sf)
   allocate(mE(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel))
   ! mE: ensemble mean
   mE = sum(Sf,2) / size(Sf,2)

   call loadVector('Forecast'//trim(ntimeindex)//'value',ModMLParallel,xf)

   ! difference between ensemble mean and the forecast
   mE = mE - xf

   ! recenter ensemble
   do k=1,size(Sf,2)      
     Sf(:,k) = Sf(:,k) - mE
   end do

!$omp parallel
   call assim(ntime,Sf,Sa)   
!$omp end parallel

 else
   write(stderr,*) 'Unknown enstype ',enstype,'. It should be 1, 2 or 3'
 end if


 call done()
#ifdef ASSIM_PARALLEL
 call parallDone()
#endif


 deallocate(xf,xa,Sa,Sf)

 contains


end program assimtest

