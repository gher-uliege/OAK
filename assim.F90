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
      Sf(:,:), Sa(:,:), &
      SfT(:,:), SaT(:,:), W(:)       
 real, pointer      :: yo(:), invsqrtR(:)
 integer            :: mjd, iargc, ntime, startntime, endntime,i
 real               :: seconds
 type(SparseMatrix) :: H
 character(len=124) :: str
 character(len=4) :: ntimeindex
  real,allocatable :: temp(:)


 if (iargc().ne.2.and.iargc().ne.3) then
   write(stderr,*) 'Usage: assim <init file> <time index> '
   write(stderr,*) 'Usage: assim <init file> <start time index> <end time index> '
   stop
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

   call saveStateVector('Analysis'//ntimeindex//'value',xa)

 end do

 call done()

#ifdef ASSIM_PARALLEL
 call parallDone()
#endif

 !  deallocate(xf,xa,W,Sa,Sf)
end program assimtest

!!$
!!$ do ntime=startntime,endntime
!!$  write(ntimeindex,'(I3.3,A)') ntime,'.'
!!$
!!$!  call loadStateVector('Forecast'//ntimeindex//'value',xf)
!!$!  call loadStateVector('Model.norm',W)
!!$
!!$!$omp parallel
!!$  do i=1,10
!!$
!!$!$omp critical (writeStdout)
!!$   write(6,*) 'ntime ', ntime
!!$!$omp end critical (writeStdout)
!!$  
!!$  call loadStateVector('Forecast'//ntimeindex//'value',xf)
!!$
!!$!$omp critical (writeStdout)
!!$   write(6,*) 'sum(xf) ',i,sum(xf),omp_get_thread_num()
!!$!$omp end critical (writeStdout)
!!$
!!$  end do
!!$!  call assim(ntime,xf,Sf,xa,Sa)
!!$!$omp end parallel
!!$
!!$  end do
