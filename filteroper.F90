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


#define DEBUG

program filteroper
 use ufileformat
 use assimilation
 use initfile

 implicit none


 real, pointer, dimension(:) :: x,y,z,lenx,leny,lenz
 real, pointer, dimension(:,:) :: Cop
 integer :: iargc

 integer, pointer     :: tmpRindex(:,:)
 real, pointer        :: tmpRcoeff(:)
 type(MemLayout) :: ML

 integer :: m,i,j,i1,j1,k1,i2,j2,k2, linindex1, linindex2,nz,nz1,nzmax, &
      status,istat
 character(len=128)   :: prefix,str

 real                 :: valex, sumcoeff

 real :: corrlen, logcorr, minlogcorr, alpha,alphax,alphay,alphaz
 real :: mincorr = 1e-3
 real :: EarthRadius = 6400000
 real :: pi = 3.141

 integer :: disti,distj,distk
 integer :: maxdisti = 5
 integer :: maxdistj = 5
 integer :: maxdistk = 0

 if (iargc().ne.1) then
   write(stderr,*) 'Usage: filteroper <initfile>   '
   write(stderr,*)     
   stop
 end if

 call getarg(1,str); initfname=str

 call MemoryLayout('Model.',ML)

 allocate(x(ML%effsize),y(ML%effsize),z(ML%effsize),&
      lenx(ML%effsize),leny(ML%effsize),lenz(ML%effsize))

 call loadVector('Correlation.lenx',ML,lenx)
 call loadVector('Correlation.leny',ML,leny)
 call loadVector('Correlation.lenz',ML,lenz)

 prefix = 'Model.'

 call loadVector(trim(prefix)//'gridX',ML,x)
 call loadVector(trim(prefix)//'gridY',ML,y)
 call loadVector(trim(prefix)//'gridZ',ML,z)


 minlogcorr = log(mincorr)

 status = 0
 nz = 0

 nzmax = 0
 do m=1,ML%nvar
   disti = maxdisti 
   distj = maxdistj 
   distk = maxdistk 

   if (all(lenx(ML%startindexsea(m):ML%endindexsea(m)).eq.0)) disti = 0
   if (all(leny(ML%startindexsea(m):ML%endindexsea(m)).eq.0)) distj = 0
   if (all(lenz(ML%startindexsea(m):ML%endindexsea(m)).eq.0)) distk = 0

   nzmax = nzmax + ML%varsizesea(m) * &
        min((2*disti+1),ML%varshape(1,m)) * &
        min((2*distj+1),ML%varshape(2,m)) * &
        min((2*distk+1),ML%varshape(3,m))
 end do

! nzmax = 3*ML%effsize
 write(stdout,*) 'nz_max ',nzmax

#   ifdef PROFILE
 call cpu_time(cputime(2))
#   endif

 write(stdout,*) 'Allocating ',(4*9*nzmax)/1000000,' MB. festhalten !!!'
 allocate(tmpRindex(8,nzmax),tmpRcoeff(nzmax))

#   ifdef PROFILE
 call cpu_time(cputime(3))
#   endif

 do m=1,ML%nvar
   write(stdout,*) 'var shape ',ML%varshape(1:3,m)

   do k2=1,ML%varshape(3,m)
!!$     do j2=50,50
!!$       do i2=50,50
       do j2=1,ML%varshape(2,m)
         do i2=1,ML%varshape(1,m)
         linindex2 = ML%StartIndex(m) + i2-1 + ML%varshape(1,m) * (j2-1 + ML%varshape(2,m) * (k2-1))         
         j = ML%SeaIndex(linindex2)

         if (j.ne.-1) then

           nz1 = nz
           sumCoeff = 0

           if (lenx(j).eq.0) then
             disti = 0
             alphax = 0
           else
             alphax = lenx(j)**(-2)
             disti = maxdisti
           end if

           if (leny(j).eq.0) then
             distj = 0
             alphay = 0
           else
             alphay = leny(j)**(-2)
             distj = maxdistj
           end if

           if (lenz(j).eq.0) then
             distk = 0
             alphaz = 0
           else
             alphaz = lenz(j)**(-2)
             distk = maxdistk
           end if


           do k1=max(1,k2-distk),min(ML%varshape(3,m),k2+distk)
             do j1=max(1,j2-distj),min(ML%varshape(2,m),j2+distj)
               do i1=max(1,i2-disti),min(ML%varshape(1,m),i2+disti)

                 linindex1 = ML%StartIndex(m) + i1-1 + ML%varshape(1,m) * (j1-1 + ML%varshape(2,m) * (k1-1))         
                 i = ML%SeaIndex(linindex1)

                 if (i.ne.-1) then

                   logcorr =  -(alphax * (x(i)-x(j))**2 + alphay * (y(i)-y(j))**2 + alphaz * (z(i)-z(j))**2)

                   if (logcorr.ge.minlogcorr) then
#                    ifdef DEBUG
                     if (nz.ge.size(tmpRcoeff)) then
                       write(stdout,*) 'genObservationOper: ERROR: ', &
                            'buffer variable too small!!! '
                       call flush(stdout,istat)
                       stop
                     end if
#                    endif

                     nz = nz+1
                     tmpRindex(1,nz) = m
                     tmpRindex(2,nz) = i2
                     tmpRindex(3,nz) = j2
                     tmpRindex(4,nz) = k2

                     tmpRindex(5,nz) = m
                     tmpRindex(6,nz) = i1
                     tmpRindex(7,nz) = j1
                     tmpRindex(8,nz) = k1

                     tmpRcoeff(nz) = exp(logcorr)
                     sumCoeff = sumCoeff+tmpRcoeff(nz)               
                   end if
                 end if
               end do
             end do
           end do

           tmpRcoeff(nz1+1:nz) = tmpRcoeff(nz1+1:nz)/sumCoeff
         end if
       end do
     end do
   end do
 end do

#   ifdef PROFILE
 call cpu_time(cputime(4))
#   endif

 write(stdout,*) 'nz ',nz
 write(stdout,*) 'dens ',sqrt(1.*nz/ML%effsize),ML%effsize
 write(stdout,*) 'nz_max ',(2*maxdisti+1)*(2*maxdistj+1)*ML%effsize


#   ifdef PROFILE
 call cpu_time(cputime(5))

 write(stdlog,*) 'Profiling: genObservationCorr'
 write(stdlog,*) 'load data  ',cputime(2)-cputime(1)
 write(stdlog,*) 'allocation ',cputime(3)-cputime(2)
 write(stdlog,*) 'main loop  ',cputime(4)-cputime(3)
 write(stdlog,*) 'copy data  ',cputime(5)-cputime(4)
 call flush(stdlog,istat)
#   endif




 allocate(Cop(9,nz))
 Cop(1:8,:) = tmpRindex(:,1:nz)
 Cop(9,:) = tmpRcoeff(1:nz)
 deallocate(tmpRindex,tmpRcoeff)

 call getInitValue(initfname,'Filter',str,default='filter.u')
 call usave(str,Cop,9999.)      
 deallocate(Cop)





end program filteroper
