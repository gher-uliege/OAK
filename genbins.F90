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


program generatebins
 use matoper
 use rrsqrt
 use ufileformat
 use date
 use initfile
 use assimilation

 implicit none

 integer :: iargc,ntime,startntime,endntime
 integer, pointer   :: Hindex(:,:)
 real, pointer      :: Hop(:,:), Hcoeff(:),yo(:),invsqrtR(:),Hshift(:)
 integer, pointer :: binlenx(:), binleny(:), binlenz(:)
 character(len=256) :: prefix, Oprefix, Dprefix

 character(len=124) :: str,path
 type(MemLayout) :: ObsML,BinsML
 type(SparseMatrix) :: Hbins,Hobs

 if (iargc().ne.2.and.iargc().ne.3) then
   write(stderr,*) 'Usage: genbins <init file> <start time index>  [<end time index>] '
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
   call fmtIndex('Diag',ntime,'.',Dprefix)
   call fmtIndex('Obs',ntime,'.',Oprefix)

   call MemoryLayout(Oprefix,ObsML)
   allocate(yo(ObsML%effsize),Hshift(ObsML%effsize),invsqrtR(ObsML%effsize))
   invsqrtR = 1.

   call loadVector(trim(Oprefix)//'value',ObsML,yo)
   call loadObservationOper(ntime,ObsML,Hobs,Hshift,invsqrtR)

   ! points out of model grid
   where (invsqrtR.eq.0.) ObsML%mask(ObsML%invindex(:)) = 0. 

   call getInitValue(initfname,trim(Oprefix)//'binlenx',binlenx)
   call getInitValue(initfname,trim(Oprefix)//'binleny',binleny)
   call getInitValue(initfname,trim(Oprefix)//'binlenz',binlenz)

   call genbins(ObsML,binlenx,binleny,binlenz,BinsML,Hindex,Hcoeff)
   call packSparseMatrix(Hindex,Hcoeff,BinsML,ObsML,Hbins)
   call saveVector(trim(Dprefix)//'yo',BinsML,Hbins.x.yo)
   call saveSparseMatrix(trim(Dprefix)//'H',BinsML,ModML,Hbins.x.Hobs)
 end do

contains

 subroutine genbins(ObsML,binlenx,binleny,binlenz,BinsML,Hindex,Hcoeff)
  use assimilation
  implicit none

  type(MemLayout), intent(in) :: ObsML
  integer, intent(in) :: binlenx(:),binleny(:),binlenz(:)
  type(MemLayout), intent(out) :: BinsML
  integer, pointer :: Hindex(:,:)
  real, pointer :: Hcoeff(:)

  integer :: m,mmax,i,j,k,ip,jp,kp,linindex,linindexpacked,linindexp

  mmax =  ObsML%nvar
  BinsML%nvar = mmax

  allocate(BinsML%ndim(mmax),BinsML%varshape(3,mmax),BinsML%varsize(mmax),BinsML%varsizesea(mmax), &
       BinsML%startIndex(mmax),BinsML%endIndex(mmax),BinsML%startIndexSea(mmax),BinsML%endIndexSea(mmax), &
       Hindex(8,ObsML%totsizesea),Hcoeff(ObsML%totsizesea))

  BinsML%ndim=ObsML%ndim
  BinsML%varshape = 1

  BinsML%startIndex(1) = 1

  do m=1,mmax
    BinsML%varshape(1,m) = ceiling(ObsML%varshape(1,m)/float(binlenx(m)))
    BinsML%varshape(2,m) = ceiling(ObsML%varshape(2,m)/float(binleny(m)))
    BinsML%varshape(3,m) = ceiling(ObsML%varshape(3,m)/float(binlenz(m)))
    BinsML%varsize(m) = product(BinsML%varshape(:,m))
    BinsML%endIndex(m) = BinsML%startIndex(m) + BinsML%varsize(m)-1
    if (m.ne.mmax) BinsML%StartIndex(m+1) = BinsML%EndIndex(m)+1
  end do

  BinsML%totsize = sum(BinsML%varsize)
  allocate(BinsML%mask(BinsML%totsize),BinsML%seaindex(BinsML%totsize))

  BinsML%mask = 0
  linindex = 1
  linindexpacked = 1
  BinsML%startIndexsea(1) = 1

  ! make the mask of the subsampled grid and count elements in each bin

  do m=1,mmax
    do k=1,ObsML%varshape(3,m)
      do j=1,ObsML%varshape(2,m)
        do i=1,ObsML%varshape(1,m)
          if (ObsML%mask(linindex).eq.1) then
            ! transform (i,j,k) -> to the subsampled (binned) grid indices (ip,jp,kp)
            ip = (i-1)/binlenx(m)+1
            jp = (j-1)/binleny(m)+1
            kp = (k-1)/binlenz(m)+1

            linindexp = BinsML%StartIndex(m) + ip-1 + BinsML%varshape(1,m) * (jp-1 + BinsML%varshape(2,m) * (kp-1))         
            ! how many elements are yet in the bin (ip,jp,kp) ?
            BinsML%mask(linindexp) = BinsML%mask(linindexp)+1

            linindexpacked = linindexpacked+1
          end if
          linindex = linindex+1
        end do
      end do
    end do

    BinsML%varsizesea(m) = count(BinsML%mask(BinsML%StartIndex(m):BinsML%endIndex(m)).ne.0)
    BinsML%endIndexsea(m) = BinsML%startIndexsea(m) + BinsML%varsizesea(m)-1
    if (m.ne.mmax) BinsML%StartIndexsea(m+1) = BinsML%EndIndexsea(m)+1
  end do

  BinsML%totsizeSea = BinsML%EndIndexsea(mmax)

  ! Hindex and Hcoeff

  linindex = 1
  linindexpacked = 1

  do m=1,mmax
    do k=1,ObsML%varshape(3,m)
      do j=1,ObsML%varshape(2,m)
        do i=1,ObsML%varshape(1,m)
          if (ObsML%mask(linindex).eq.1) then
            ! transform (i,j,k) -> to the subsampled (binned) grid indices (ip,jp,kp)
            ip = (i-1)/binlenx(m)+1
            jp = (j-1)/binleny(m)+1
            kp = (k-1)/binlenz(m)+1

            Hindex(:,linindexpacked) = (/ m,ip,jp,kp,m,i,j,k /)

            linindexp = BinsML%StartIndex(m) + ip-1 + BinsML%varshape(1,m) * (jp-1 + BinsML%varshape(2,m) * (kp-1))         
            Hcoeff(linindexpacked) = 1./BinsML%mask(linindexp)

            linindexpacked = linindexpacked+1
          end if
          linindex = linindex+1
        end do
      end do
    end do
  end do

  ! mask with only 1 and 0

  where (BinsML%mask.gt.1) BinsML%mask=1

  ! seaindex and invindex

  allocate(BinsML%invindex(BinsML%totsizesea))

  BinsML%SeaIndex = -1
  j=1
  do i=1,BinsML%totsize
    if (.not.BinsML%removeLandPoints.or.BinsML%Mask(i).eq.1) then
      BinsML%SeaIndex(i) = j
      BinsML%invindex(j) = i
      j=j+1
    end if
  end do

  BinsML%removeLandPoints = .true.
  BinsML%effsize = BinsML%totsizeSea


 end subroutine genbins

end program generatebins
