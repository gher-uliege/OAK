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

!_______________________________________________________
!


program oper
use matoper
use rrsqrt
implicit none

integer, parameter :: jmax = 5,	mmax = 3, kmax = 2
integer :: info
real :: xf(jmax), xa(jmax), yo(mmax), Hxf(mmax)
real :: Sf(jmax,kmax), Sa(jmax,kmax), invsqrtR(mmax), HSf(mmax,kmax)
real :: Ef(jmax,kmax+1), Ea(jmax,kmax+1), HEf(mmax,kmax+1), ampl(kmax+1)
real :: H(mmax,jmax), tmp(2,2) 
real :: SfT(kmax,jmax), SaT(kmax,jmax), &
  U(jmax,jmax), lam(kmax), VT(kmax,kmax), valex
real,pointer :: S(:,:), Sr(:,:), SrT(:,:)
  real :: M3x3(3,3),tmp1000(1000)

  integer :: i
  real :: seconds
  integer, pointer        :: variable(:)
  real,    pointer        :: x(:),y(:),z(:),observation(:), &
                             invsqrtRo(:)
  logical,    pointer     :: out(:)
  type(SparseMatrix) :: Hp, H2
  integer, pointer :: Hindex(:,:)
  real, pointer :: Hcoeff(:), Hr(:,:)
  character(len=124), pointer :: filenames(:)

Hp%m = mmax; Hp%n = jmax; Hp%nz = 4;
allocate(Hp%i(4),Hp%j(4),Hp%s(4))
Hp%i(1)=1; Hp%j(1)=1; Hp%s(1)=1.1;
Hp%i(2)=2; Hp%j(2)=2; Hp%s(2)=1.2;
Hp%i(3)=3; Hp%j(3)=3; Hp%s(3)=1.3;
Hp%i(4)=2; Hp%j(4)=5; Hp%s(4)=.3;


H = reshape((/1.1,0. ,0. ,  &
              0. ,1.2,0. ,  &
              0. ,0. ,1.3,  &
              0. ,0. ,0. ,  &
              0. ,0.3,0.  /),(/mmax,jmax/))

Ef = reshape((/ (i,i=1,jmax*(kmax+1)) /),(/jmax,kmax+1/))
Ef = sin(1.*Ef)


HEf = H.x.Ef



xf = (/5,7,3,8,1/)
yo = (/2,3,5/)
Hxf = H.x.xf
Hxf = matmul(H,xf)
Sf = reshape((/1,2,3,4,5,  6,7,8,9,10/),(/jmax,kmax/))
SfT = transpose(Sf)
HSf = H.x.Sf
invsqrtR = (/1.,2.,3./)**(-1)
M3x3 = HSf.xt.H

!call ensanalysis(Ef,HEf,yo,invsqrtR,Ea, ampl)

!write(6,*) 'Ea ',Ea
!write(6,*) 'mean(Ea) ',sum(Ea,2)/(kmax+1)

write(6,*) 'enscov(Ef) ',enscov(Ef)
!write(6,*) 'enscov(Ea) ',enscov(Ea)


call ens2sqrt(Ef,xf,Sf) 

!write(6,*) '(Sf.xt.Sf) ',(Sf.xt.Sf)
write(6,*) '(Sf.xt.Sf) ',matmul(Sf,transpose(Sf))
write(6,*) 'xf ',xf

call sqrt2ens(xf,Sf,Ea) 

write(6,*) 'enscov(Ef) ',enscov(Ea)
write(6,*) 'xf ',sum(Ea,2)/(kmax+1)

stop



end program 

!  0.8090215      0.8994658      0.4245623       1.161173
