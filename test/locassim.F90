!
!  OAK, Ocean Assimilation Kit
!  Copyright(c) 2002-2015 Alexander Barth and Luc Vandenblucke
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

program locassim
 use matoper
 use covariance
 use rrsqrt
 use ufileformat
 use initfile
 use assimilation

 implicit none

 real, allocatable  :: xf(:), xa(:), &
      Sf(:,:), Sa(:,:), &
      Ef(:,:), Ea(:,:), &
      Hc(:,:), yo(:)
 
 class(DiagCovar), allocatable :: R
 type(SparseMatrix) :: H

 integer            :: iargc, ntime, enstype
 character(len=124) :: str
 character(len=124) :: ntimeindex

 if (iargc().ne.2) then
   write(stderr,*) 'Usage: assim <init file> <time index> '
   ERROR_STOP
 end if

#ifdef ASSIM_PARALLEL
   call parallInit()
#endif

 call getarg(1,str); call init(str)
 call getarg(2,str); read(str,*) ntime  

 allocate(xf(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel), &
      xa(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel), &
      Ea(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel,ErrorSpaceDim), &
      Ef(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel,ErrorSpaceDim), &
      Sa(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel,ErrorSpaceDim-1), &
      Sf(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel,ErrorSpaceDim-1))

 !write(ntimeindex,'(I3.3,A)') ntime,'.'
 call fmtIndex('',ntime,'.',ntimeindex)

 call getInitValue(initfname,'ErrorSpace.type',enstype,default=1)
 if (enstype == 1) then
   ! Sf are error modes and xf is the forecast
   call loadVectorSpace('ErrorSpace.init',ModMLParallel,Sf)
   call loadVector('Forecast'//trim(ntimeindex)//'value',ModMLParallel,xf)
!$omp parallel
!!!   call assim(ntime,Sf,Sa,xf,xa)
!$omp end parallel
 else
   ! Sf are ensemble members
   call loadEnsemble('ErrorSpace.init',ModMLParallel,Ef)
   call ens2sqrt(Ef,xf,Sf)

   call fmtIndex('',ntime,'.',infix)
   call MemoryLayout('Obs'//trim(infix),ObsML,rmLPObs)
   m = ObsML%effsize
   k = size(Sf,2)

  allocate(yo(m),invsqrtR(m),Hxf(m),Hxa(m),HSf(m,k),HSa(m,k), &
       yo_Hxf(m), yo_Hxa(m), innov_projection(m), Hshift(m))

  call loadObsTime(ntime,mjd,error)
  call loadObs(ntime,ObsML,yo,invsqrtR)    
   

   

   call locensanalysis(xf,Sf,H,yo,R,lpoints,hc,xa,Sa)


   

   

!$omp parallel
   call assim(ntime,Sf,Sa)   
!$omp end parallel
 end if

 call done()
#ifdef ASSIM_PARALLEL
 call parallDone()
#endif


 deallocate(xf,xa,Sa,Sf)


end program locassim


