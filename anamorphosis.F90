!
!  OAK, Ocean Assimilation Kit
!  Copyright(c) 2002-2012 Alexander Barth and Luc Vandenblucke
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
#define stddebug 6

!#define LINEAR_ANAMORPHOSIS
!#define LN_ANAMORPHOSIS


module anamorphosis

 integer, parameter :: maxLen = 256;
 real, save, pointer :: anamorphosisx(:), anamorphosisy(:)

! Thermal expansion coefficient and
! Saline contraction coefficient 
! for T=15 and S=35

  real, parameter :: alpha = 2.1474e-04
  real, parameter :: beta = 7.5137e-04
!  real, parameter :: minimumx = 1e-15
!  real, parameter :: minimumx = 1e-10
  real, parameter :: minimumx = 1e-8


 type AnamTransVar
   ! type: 
   ! 1 = idenity
   ! 2 = log
   ! 3 = mapping x -> y

   integer :: type
   real, pointer :: x(:)
   real, pointer :: y(:)
 end type AnamTransVar

 type AnamorphosisTrans
   ! one anam structure for each variable
   type(AnamTransVar), pointer :: anam(:)   
 end type AnamorphosisTrans
  
 type(AnamorphosisTrans), save :: AnamTrans

contains
 subroutine initAnamorphosis(initfname)
  use ufileformat
  use initfile
  implicit none
  character(len=*) :: initfname
  character(len=maxLen), pointer :: xname(:), yname(:)
  integer :: i
  real :: valex
  character(len=30) :: infoformat = '(A50,2E14.5)'
  character(len=maxLen) :: path

  if (.not.presentInitValue(initfname,'Anamorphosis.x')) then
    write(stddebug,'("== no anamorphosis transform loaded ==")')    
    return
  end if

# ifdef DEBUG
  write(stddebug,'("== load anamorphosis transform ==")')
# endif

  call getInitValue(initfname,'Anamorphosis.path',path,default='')
  call getInitValue(initfname,'Anamorphosis.x',xname)
  call getInitValue(initfname,'Anamorphosis.y',yname)

  allocate(AnamTrans%anam(size(xname)))

  do i=1,size(xname)
    if (xname(i) == '_id') then
      AnamTrans%anam(i)%type = 1
      write(stddebug,*) 'variable ',i,' id '
    elseif (xname(i) == '_log') then
      write(stddebug,*) 'variable ',i,' log '
      AnamTrans%anam(i)%type = 2
    else
      write(stddebug,*) 'variable ',i,' custom '
      AnamTrans%anam(i)%type = 3
      call uload(trim(path)//xname(i),AnamTrans%anam(i)%x,valex)
      call uload(trim(path)//yname(i),AnamTrans%anam(i)%y,valex)

#ifdef DEBUG
      write(stddebug,'(A50,A14,A14)') 'loaded variable','min','max'
      write(stddebug,infoformat) trim(path)//trim(xname(i)),       &
           minval(AnamTrans%anam(i)%x),maxval(AnamTrans%anam(i)%x)
      write(stddebug,infoformat) trim(path)//trim(yname(i)),       &
           minval(AnamTrans%anam(i)%y),maxval(AnamTrans%anam(i)%y)
#endif


    end if      
  end do

  deallocate(xname,yname)
 end subroutine initAnamorphosis


 subroutine doneAnamorphosis
  implicit none
  deallocate(AnamTrans%anam)
 end subroutine doneAnamorphosis

 subroutine initAnamorphosis_old(initfname)
  use ufileformat
  use initfile
  implicit none
  character(len=*) :: initfname
  character(len=256) :: str
  real :: valex

  call getInitValue(initfname,'Anamorphosis.X',str)
#ifdef DEBUG
  write(stddebug,*) 'load Anamorphosis.X',str
#endif
  call uload(str,anamorphosisx,valex)

  call getInitValue(initfname,'Anamorphosis.Y',str)
#ifdef DEBUG
  write(stddebug,*) 'load Anamorphosis.Y',str
#endif
  call uload(str,anamorphosisy,valex)

 end subroutine 


 subroutine TStransform(mask,z,T,S,X,Y)
  implicit none
  integer, dimension(:,:,:), intent(in) :: mask
  real, dimension(:,:,:), intent(in) :: z,T,S
  real, dimension(:,:,:), intent(out) :: X,Y

  integer :: i,j,k,imax,jmax,kmax,nbout
  real :: dT,dS,dZ
  logical :: out

  imax = size(mask,1)
  jmax = size(mask,2)
  kmax = size(mask,3)
  nbout = 0

  !  write(6,*) 'min dz ',minval(Z(:,:,
  do k=1,kmax
    do j=1,jmax
      do i=1,imax
        if (mask(i,j,k).eq.1) then
          !      if (mask(i,j,k-1).eq.0.or.k.eq.1) then
          if (k.eq.1.or.mask(i,j,k-1).eq.0) then
            X(i,j,k) = T(i,j,k)
            Y(i,j,k) = S(i,j,k)
          else 
            dZ = Z(i,j,k)-Z(i,j,k-1)
            dT = alpha * (T(i,j,k)-T(i,j,k-1))/dZ
            dS = beta * (S(i,j,k)-S(i,j,k-1))/dZ

            X(i,j,k) = anamorphfun(dT - dS,out)
            Y(i,j,k) = dT + dS

            if (out) nbout = nbout+1
          end if
        end if
      end do
    end do
  end do

  if (nbout.ne.0) write(stdout,*) 'nb of out of domain ',nbout
 end subroutine TStransform

 !___________________________________________________________________
 !

 subroutine invTStransform(mask,z,X,Y,T,S)
  implicit none
  integer, dimension(:,:,:), intent(in) :: mask
  real, dimension(:,:,:), intent(in) :: z,X,Y
  real, dimension(:,:,:), intent(out) :: T,S

  integer :: i,j,k,imax,jmax,kmax,nbout
  real :: dT,dS,dZ,sumXY,diffXY,X2
  logical :: out


  imax = size(mask,1)
  jmax = size(mask,2)
  kmax = size(mask,3)
  nbout = 0

  do k=1,kmax
    do j=1,jmax
      do i=1,imax
        if (mask(i,j,k).eq.1) then
          !      if (mask(i,j,k-1).eq.0.or.k.eq.1) then
          if (k.eq.1.or.mask(i,j,k-1).eq.0) then
            T(i,j,k) = X(i,j,k)
            S(i,j,k) = Y(i,j,k)
          else 
            dZ = Z(i,j,k)-Z(i,j,k-1)
            X2 = invanamorphfun(X(i,j,k),out)
            dT = (X2 + Y(i,j,k))/2.
            dS = (-X2 + Y(i,j,k))/2.

            T(i,j,k) = T(i,j,k-1) + dT*dZ/alpha
            S(i,j,k) = S(i,j,k-1) + dS*dZ/beta

            if (out) nbout = nbout+1
          end if
        end if
      end do
    end do
  end do

  if (nbout.ne.0) write(stdout,*) 'nb of out of domain ',nbout

 end subroutine invTStransform



 function anamorphfun(x,out) result(y)
  implicit none
  real, intent(in) :: x
  logical, intent(out) :: out
  real :: y

#if defined(LINEAR_ANAMORPHOSIS)
  y = x
  out = .false.
#elif defined(LN_ANAMORPHOSIS)
  y = log(max(x,minimumx));
  out = .false.
#else
  y = interp1(anamorphosisy,anamorphosisx,x,out)
  !  write(6,*) 'x,y,',x,y,out
  !  write(6,*) 'x,y,',1.,interp1(anamorphosisy,anamorphosisx,1.,out),out
  !  stop

#endif
 end function anamorphfun

 function invanamorphfun(x,out) result(y)
  implicit none
  real, intent(in) :: x
  logical, intent(out) :: out
  real :: y

#if defined(LINEAR_ANAMORPHOSIS)
  y = x
  out = .false.
#elif defined(LN_ANAMORPHOSIS)
  y = exp(x);
  out = .false.
#else
  y = interp1(anamorphosisx,anamorphosisy,x,out)
#endif

 end function invanamorphfun


!!$  function INTERP1(X,Y,XI,out) result(yi)
!!$   implicit none
!!$   real, intent(in) :: x(:),y(:),xi(:)
!!$   logical, optional, intent(out) :: out(size(xi))
!!$   real :: yi(size(xi))
!!$
!!$   real :: alpha
!!$   integer :: j,k,kp
!!$
!!$   
!!$
!!$   do j=1,size(xi)
!!$
!!$     k = -1
!!$
!!$     do kp = 1,size(x)-1
!!$       if (x(kp).le.xi(j).and.xi(j).lt.x(kp+1)) then
!!$         k = kp
!!$         exit
!!$       end if
!!$     end do
!!$
!!$     if (k.ne.-1) then
!!$       alpha = (xi(j)-x(k))/(x(k+1)-x(k))
!!$       yi(j) = alpha * y(k) + (1-alpha) * y(k+1)
!!$     else
!!$       if (xi(j).lt.x(1)) then
!!$         yi = y(1)
!!$       else
!!$         yi = y(size(x))
!!$       end if
!!$     end if
!!$
!!$     if (present(out)) then
!!$       out(j)=k.eq.-1
!!$     end if
!!$   end do
!!$
!!$  end function INTERP1



 function interp1(X,Y,XI,out) result(yi)
  implicit none
  real, intent(in) :: x(:),y(:),xi
  logical, optional, intent(out) :: out
  real :: yi

  real :: alpha
  integer :: j,k,kp

  yi = xi
  k = -1

  do kp = 1,size(x)-1
    if (x(kp).le.xi.and.xi.lt.x(kp+1)) then
      k = kp
      exit
    end if
  end do

  if (k.ne.-1) then
    alpha = (xi-x(k))/(x(k+1)-x(k))
    yi = (1-alpha) * y(k) + alpha * y(k+1)
  else
    if (xi.lt.x(1)) then
      yi = y(1)
    else
      yi = y(size(x))
    end if
  end if

  if (present(out)) then
    out=k.eq.-1
  end if


 end function interp1
end module

