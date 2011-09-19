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

  implicit none

  real, allocatable  :: variance(:), covariance(:), &
                        Sf(:,:)

  integer            :: iargc,i,j,k,v,index
  character(len=124) :: str
  type(MemLayout) :: ML
  logical :: valid

 if (iargc().ne.1) then
   write(stderr,*) 'Usage: covar <init file>  '
   write(stderr,*) 'Usage: assim <init file> <start time index> <end time index> '
   stop
 end if

 call getarg(1,str); initfname = str

  call MemoryLayout('Model.',ML)
  call getInitValue(initfname,'ErrorSpace.dimension',ErrorSpaceDim,default=0)

 allocate(variance(ML%effsize), &
    covariance(ML%effsize), &
    Sf(ML%effsize,ErrorSpaceDim))


 call getInitValue(initfname,'Point.i',i)
 call getInitValue(initfname,'Point.j',j)
 call getInitValue(initfname,'Point.k',k)
 call getInitValue(initfname,'Point.v',v)

! from the three spatial indexes (i,j,k) and the variable index v 
! get the index within the state vector

 index = sub2ind(ML,v,i,j,k,valid)

 if (.not.valid) then
   write(stderr,*) 'It seems to me that the point (i,j,k) = ',i,j,k,' of the ',v,'th variable is a land point.'
   write(stderr,*) 'check the land mask (Model.mask)'
   stop
 end if

! load the error modes of the covariance matrix

 call loadVectorSpace('ErrorSpace.init',ML,Sf)

 variance = sum(Sf**2,2)

  if (presentInitValue(initfname,'Variance')) then 
    call saveVector('Variance',ML,variance)
  end if

 covariance = sum(spread(Sf(index,:),1,ML%effsize) * Sf,2)

  if (presentInitValue(initfname,'Covariance')) then 
    call saveVector('Covariance',ML,covariance)
  end if

  if (presentInitValue(initfname,'Correlation')) then 
    where (variance.lt.0) variance = 0.
    covariance = covariance/sqrt(variance*variance(index))

    call saveVector('Correlation',ML,covariance)
  end if

  deallocate(variance,covariance,Sf)
end program 

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
