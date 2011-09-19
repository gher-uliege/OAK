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


program genobsoper
 use matoper
 use rrsqrt
 use ufileformat
 use date
 use initfile
 use assimilation

 implicit none

 integer :: iargc,ntime,startntime,endntime
 integer, pointer   :: Hindex(:,:)
 real, pointer      :: Hop(:,:), Hcoeff(:)
 character(len=124) :: str,path,prefix
 type(MemLayout)    :: ObsML
 type(sparseMatrix) :: H

 if (iargc().ne.2.and.iargc().ne.3) then
   write(stderr,*) 'Usage: genobsoper <init file> <start time index>  [<end time index>] '
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
   write(prefix,'(A,I3.3,A)') 'Obs',ntime,'.'
   call MemoryLayout(prefix,ObsML)

   call genObservationOper(ntime,ObsML,Hindex,Hcoeff)

   allocate(Hop(9,size(Hcoeff)))
   Hop(1:8,:) = Hindex
   Hop(9,:) = Hcoeff

   write(prefix,'(A,I3.3,A)') 'Diag',ntime,'.'
   call getInitValue(initfname,trim(prefix)//'path',path,default='')
   call getInitValue(initfname,trim(prefix)//'H',str)

   call usave(trim(path)//str,Hop,9999.)

!#ifdef .true.
   call packSparseMatrix(Hindex,Hcoeff,ObsML,ModML,H)
   call usave(trim(path)//trim(str)//'.i',real(H%i(1:H%nz)),9999.)
   call usave(trim(path)//trim(str)//'.j',real(H%j(1:H%nz)),9999.)
   call usave(trim(path)//trim(str)//'.s',H%s(1:H%nz),9999.)
!#endif

 end do
end program genobsoper
