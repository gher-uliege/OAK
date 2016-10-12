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
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, 
!  Boston, MA  02110-1301, USA.
!


module covariance

use matoper

integer, parameter :: locensanalysis_SSt = 1
integer, parameter :: locensanalysis_Pc = 2


! Abstract covariance matrix

type, abstract :: Covar
  ! size of the covariance matrix
  integer         :: n
 contains
  procedure(Covar_mtimes_vec), deferred :: mtimes_vec

  procedure :: print => Covar_print
  procedure :: full => Covar_full
  procedure :: mtimes_mat => Covar_mtimes_mat
  procedure :: mldivide_vec => Covar_mldivide_vec
  procedure :: mldivide_mat => Covar_mldivide_mat

  generic, public :: mtimes => mtimes_vec, mtimes_mat
  generic, public :: mldivide => mldivide_vec, mldivide_mat
end type Covar

abstract interface
  function Covar_mtimes_vec(this,x) result(M)
   import Covar
   class(Covar), intent(in) :: this
   real, intent(in) :: x(:)
   real :: M(size(x,1))
  end function Covar_mtimes_vec
end interface


type, extends(Covar) :: DiagCovar
  real, pointer :: diag(:)
   contains
        procedure :: init => DiagCovar_init
        procedure :: full => DiagCovar_full
        procedure :: mtimes_vec => DiagCovar_mtimes_vec
        procedure :: mldivide_vec => DiagCovar_mldivide_vec
end type DiagCovar

#ifdef SUPPORT_CON
interface DiagCovar
  module procedure newDiagCovar
end interface DiagCovar
#endif


type, extends(Covar) :: SMWCovar
  real, allocatable :: C(:), B(:,:), D(:,:)
  contains
        procedure :: initialize => SMWCovar_init
        procedure :: full => SMWCovar_full
        procedure :: mtimes_vec => SMWCovar_mtimes_vec
        procedure :: mldivide_vec => SMWCovar_mldivide_vec
end type SMWCovar


interface
  subroutine locpoints_(i,nnz,j,w,onlyj)
   integer, intent(in) :: i
   integer, intent(out) :: nnz,j(:)
   real, intent(out) :: w(:)
   integer, optional, intent(in) :: onlyj(:)  
  end subroutine locpoints_
end interface

type, extends(Covar) :: LocCovar
  real          :: len
  real, pointer :: x(:,:)
  real, pointer :: S(:,:)
  procedure(locpoints_), pointer, nopass:: lpoints
   contains
        procedure :: init => LocCovar_init
        procedure :: mtimes_vec => LocCovar_mtimes_vec
        procedure :: mtimes_max => LocCovar_mtimes_mat
 end type LocCovar

#ifdef SUPPORT_CON
interface LocCovar
  module procedure newLocCovar
end interface LocCovar
#endif

type, extends(Covar) :: ConsCovar
  real, pointer :: h(:,:)
!  class(LocCovar), pointer :: C
  class(Covar), pointer :: C
   contains
        procedure :: initialize => ConsCovar_init
        procedure :: mtimes_vec => ConsCovar_mtimes_vec
        procedure :: mtimes_mat => ConsCovar_mtimes_mat
end type ConsCovar

interface matmul
  module procedure              &
        Covar_matmul_vec,       &
        Covar_matmul_mat
end interface matmul

interface operator(.x.)
  module procedure              &
        Covar_matmul_vec,       &
        Covar_matmul_mat
end interface 

!interface operator(*)
!  module procedure              &
!        Covar_matmul_vec,       &
!        Covar_matmul_mat
!end interface 


 contains

!---------------------------------------------------------------------
! Covar: abstract covariance matrix
!---------------------------------------------------------------------

  subroutine Covar_print(this)
   implicit none
   class(Covar), intent(in) :: this
   print *, 'Covar: ', this%n,this%n
  end subroutine Covar_print

!---------------------------------------------------------------------

  function Covar_full(this) result(M)
   implicit none
   class(Covar), intent(in) :: this
   real :: M(this%n,this%n)
   integer :: i
   real :: e(this%n)

   ! apply Covar to identity matrix
   e = 0
   do i = 1,this%n
     e(i) = 1.
     M(:,i) = this%mtimes_vec(e)
     e(i) = 0.
   end do
  end function Covar_full

!---------------------------------------------------------------------
 
  function Covar_mtimes_mat(this,A) result(Q)
   implicit none
   class(Covar), intent(in) :: this
   real, intent(in)         :: A(:,:)
   real                     :: Q(size(A,1),size(A,2))
   integer :: i

   do i = 1,size(A,2)
     Q(:,i) = this%mtimes_vec(A(:,i))
   end do
  end function Covar_mtimes_mat

!---------------------------------------------------------------------

  function Covar_mldivide_mat(this,A) result(Q)
   implicit none
   class(Covar) :: this
   real, intent(in) :: A(:,:)
   real :: Q(size(A,1),size(A,2))
   integer :: i

   do i = 1,size(A,2)
     !Q(:,i) = Covar_mldivide_vec(this,A(:,i))
     Q(:,i) = this%mldivide_vec(A(:,i))
   end do
  end function Covar_mldivide_mat

!---------------------------------------------------------------------

  function Covar_mldivide_vec(this,v) result(Q)
   use matoper
   implicit none
   class(Covar) :: this
   real, intent(in) :: v(:)
   real :: Q(size(v))

   Q = pcg(fun,v)

   contains
    function fun(x) result(y)
     implicit none
     real, intent(in) :: x(:)
     real :: y(size(x))
     
     y = matmul(this,x)
    end function fun

  end function Covar_mldivide_vec

  !---------------------------------------------------------------------
  ! for generic function matmul and .x.

  function Covar_matmul_mat(this,x) result (Px)
   implicit none
   class(Covar), intent(in) :: this
   real, intent(in) :: x(:,:)
   real ::  Px(size(x,1),size(x,2))

   Px = this%mtimes(x)
  end function 

  function Covar_matmul_vec(this,x) result (Px)
   implicit none
   class(Covar), intent(in) :: this
   real, intent(in) :: x(:)
   real ::  Px(size(x,1))

   Px = this%mtimes(x)
  end function 

!---------------------------------------------------------------------
! DiagCovar: diagonal covariance matrix
!---------------------------------------------------------------------

  function newDiagCovar(diag) result(this)
    type(DiagCovar) :: this
    real, target :: diag(:)

    this%n = size(diag)
    this%diag => diag 
  end function newDiagCovar

!---------------------------------------------------------------------

  subroutine DiagCovar_init(this,C)
    class(DiagCovar) :: this
    real, target :: C(:)

    this%n = size(C)
    this%diag => C 
   end subroutine DiagCovar_init

!---------------------------------------------------------------------
 
  function DiagCovar_full(this) result(M)
   implicit none
   class(DiagCovar), intent(in) :: this
   real :: M(this%n,this%n)

   M = diag(this%diag)
  end function

!---------------------------------------------------------------------

  function DiagCovar_mtimes_vec(this,x) result(C)
   implicit none
   class(DiagCovar), intent(in) :: this
   real, intent(in) :: x(:)
   real :: C(size(x,1))

   C = this%diag*x    
  end function DiagCovar_mtimes_vec

!---------------------------------------------------------------------
 
  function DiagCovar_mldivide_vec(this,v) result(Q)
   implicit none
   class(DiagCovar) :: this
   real, intent(in) :: v(:)
   real :: Q(size(v))

   Q = v / this%diag
  end function DiagCovar_mldivide_vec


!---------------------------------------------------------------------
! SMWCovar: covariance matrix C +  B*B'
! Create the covariance matrix CSMW = C +  B*B' than can be inverted efficiently
! using the Sherman-Morrison-Woodbury formula, if size(B,2) is much smaller than
! size(C,1):
!
! inv(C +  B*B') = inv(C) -  inv(C)*B*inv(B'*inv(C)*B +  I)*B'*inv(C) 
!---------------------------------------------------------------------

  subroutine SMWCovar_init(this,C,B)
   use matoper
   implicit none   
   class(SMWCovar) :: this
   real, intent(in) :: B(:,:), C(:)

   allocate(this%C(size(C)), &
        this%B(size(B,1),size(B,2)), &
        this%D(size(B,2),size(B,2)))

   this%n = size(C)
   this%C = C
   this%B = B
   this%D = inv(matmul(transpose(B), (1/C) .dx. B) + eye(size(B,2)))
  end subroutine SMWCovar_init


!---------------------------------------------------------------------
 
  function SMWCovar_full(this) result(M)
   implicit none
   class(SMWCovar), intent(in) :: this
   real :: M(this%n,this%n)

   M = diag(this%C) + matmul(this%B,transpose(this%B))
  end function
  
!---------------------------------------------------------------------

  function SMWCovar_mtimes_vec(this,x) result(C)
    class(SMWCovar), intent(in) :: this
    real, intent(in) :: x(:)
    real :: C(size(x,1))

    C = this%C*x + matmul(this%B,matmul(transpose(this%B),x))
   end function SMWCovar_mtimes_vec

!---------------------------------------------------------------------
 
  function SMWCovar_mldivide_vec(this,v) result(Q)
   implicit none
   class(SMWCovar) :: this
   real, intent(in) :: v(:)
   real :: Q(size(v))

   Q = v / this%C

   Q = Q - matmul(this%B,matmul(this%D,(matmul(transpose(this%B),Q))))/this%C
  end function SMWCovar_mldivide_vec


!---------------------------------------------------------------------
! LocCovar: localized ensemble covariance matrix
!---------------------------------------------------------------------

  function locfun(r) result(fun)
   real :: r,fun

   ! optim
   if (r <= 1.) then
     !fun = -r**5/4 + r**4/2 + 5*r**3/8 - 5./3*r**2 + 1
     fun = (((-r/4. + 1./2.) * r + 5./8.) * r - 5./3.) * r**2 + 1.

   elseif (r <= 2.) then
     !fun = r**5/12 - r**4/2 + 5*r**3/8 + 5./3*r**2 - 5*r + 4 - &
     ! 2    ./(3*r)
     fun = ((((r/12. - 1./2.) * r + 5./8.) * r + 5./3.) * r - 5.) * r + 4 - &
          2./(3*r)
   else
     fun = 0
   end if

!   write(6,*) 'locfun ',fun,r

  end function locfun



  ! cartesian distance
  function cdist(x0,x1) result(dist)
   real, intent(in) :: x0(:), x1(:)
   real :: dist

   dist = sqrt(sum((x0-x1)**2))
   
  end function cdist

  ! return all grid points with non-zero weight
  ! relative to grid point i

  ! j,w should be large enought

  subroutine locpoints(i,x,len,nnz,j,w,onlyj)
   integer, intent(in) :: i
   real, intent(in) :: len, x(:,:)
   integer, intent(out) :: nnz, j(:)
   real, intent(out) :: w(:)
   integer, optional, intent(in) :: onlyj(:)  

   integer :: k
   real :: dist, tmp

   nnz = 0

   if (present(onlyj)) then
     do k = 1,size(onlyj)
       dist = cdist(x(i,:),x(onlyj(k),:))

       tmp = locfun(dist/len)

       if (tmp /= 0.) then
         nnz = nnz + 1
         j(nnz) = onlyj(k)
         w(nnz) = tmp
       end if
     end do

   else
     do k = 1,size(x,1)
       dist = cdist(x(i,:),x(k,:))

       tmp = locfun(dist/len)

       if (tmp /= 0.) then
         nnz = nnz + 1
         j(nnz) = k
         w(nnz) = tmp
       end if
     end do
   end if
  end subroutine locpoints

  subroutine LocCovar_init(this, S, lp)
   class(LocCovar) :: this
   real, pointer :: S(:,:)
   procedure(locpoints_) :: lp

   this%n = size(S,1)
   this%lpoints => lp
   this%S => S
  end subroutine LocCovar_init

  function newLocCovar(S, lp) result(LC)
   type(LocCovar) :: LC
   real, pointer :: S(:,:)
   procedure(locpoints_) :: lp

   LC%n = size(S,1)
   LC%lpoints => lp
   LC%S => S
  end function newLocCovar


  function LocCovar_mtimes_vec(this,x) result (Px)
   class(LocCovar), intent(in) :: this
   real, intent(in) :: x(:)
   real :: Px(size(x,1))

   integer :: j(this%n), i, nnz
   real :: rho(this%n), p(this%n), tmp(this%n)

   Px = 0

   do i = 1,this%n
    call this%lpoints(i,nnz,j,rho)
    p(1:nnz) = matmul(this%S(i,:),transpose(this%S(j(1:nnz),:)))    
    tmp(1:nnz) = rho(1:nnz) * p(1:nnz)
    Px(i) = Px(i) + sum(tmp(1:nnz) * x(j(1:nnz)))
!      if (mod(i,1) == 0) write(6,*) 'i,n,',i,this%n,nnz, Px(i)
!      if (mod(i,5000) == 0) stop
  end do
 end function LocCovar_mtimes_vec

 function LocCovar_mtimes_mat(this,x) result (Px)
   class(LocCovar), intent(in) :: this
   real, intent(in) :: x(:,:)
   real ::  Px(size(x,1),size(x,2))

   integer :: j(this%n), i, k, nnz
   real :: rho(this%n), p(this%n), tmp(this%n)
   Px = 0

   do i = 1,this%n
     call this%lpoints(i,nnz,j,rho)
     p(1:nnz) = matmul(this%S(i,:),transpose(this%S(j(1:nnz),:)))    
     tmp(1:nnz) = rho(1:nnz) * p(1:nnz)
     
     do k = 1,size(x,2)
       Px(i,k) = Px(i,k) + sum(tmp(1:nnz) * x(j(1:nnz),k))
       if (mod(i,1000) == 0) write(6,*) 'i,n,',i,this%n,nnz, Px(i,k)
     end do
!      if (mod(i,1000) == 0) stop
   end do


  end function LocCovar_mtimes_mat


  ! compute H * P * H' y : matmul(matul(this, transpose(H),H) y
  ! where H is sparse

  function loccovar_project(this,H,y) result (HPHy)
   class(LocCovar), intent(in) :: this
   type(SparseMatrix), intent(in) :: H
   real, intent(in) :: y(:)
   
   real ::  PHy(this%n)
   integer :: j(this%n), i, nnz, Hjind
   real :: rho(this%n), p(this%n), tmp(this%n)
   real :: Hty(this%n), HPHy(size(y))
   integer :: unique_Hj(H%nz), unique_Hjnz
   
   Hty = y.x.H
!   Cy = H.x.(LC.x.(Hty)) 

   unique_Hj = H%j
!   write(6,*) ' unique_Hj ',H%nz,unique_Hj
   call unique(unique_Hj,unique_Hjnz)
   PHy = 0

!   do i = 1,this%n
   do Hjind = 1,unique_Hjnz
     i = unique_Hj(Hjind)
     
     call this%lpoints(i,nnz,j,rho,onlyj = unique_Hj)
     p(1:nnz) = matmul(this%S(i,:),transpose(this%S(j(1:nnz),:)))    
     tmp(1:nnz) = rho(1:nnz) * p(1:nnz)
     
     PHy(i) = PHy(i) + sum(tmp(1:nnz) * Hty(j(1:nnz)))
   end do

   HPHy = H.x.PHy
  end function loccovar_project


!---------------------------------------------------------------------
! CovarCovar: localized ensemble covariance matrix with global constrains
!---------------------------------------------------------------------

  subroutine ConsCovar_init(this,C,h)
   class(ConsCovar) :: this
!   class(LocCovar), pointer :: C
   class(Covar), target :: C
   real, pointer :: h(:,:)

   this%n = C%n
   this%C => C
   this%h => h
  end subroutine ConsCovar_init

!---------------------------------------------------------------------

  function newConsCovar(C, H) result(this)
   type(ConsCovar) :: this
   class(Covar), target :: C
   real, pointer :: H(:,:)

   this%n = C%n
   this%C => C
   this%h => h
  end function newConsCovar

!---------------------------------------------------------------------

  function ConsCovar_mtimes_vec(this,x) result (Px)
   class(ConsCovar), intent(in) :: this
   real, intent(in) :: x(:)
   real ::  Px(size(x,1)), x2(this%n) 

   x2 = x - matmul(this%h,matmul(transpose(this%h),x))
   Px = this%C .x. x2
   Px = Px - matmul(this%h,matmul(transpose(this%h),Px))
  end function ConsCovar_mtimes_vec

!---------------------------------------------------------------------

  function ConsCovar_mtimes_mat(this,A) result (Px)
   implicit none
   class(ConsCovar), intent(in) :: this
   real, intent(in) :: A(:,:)
   real ::  Px(size(A,1),size(A,2)), x2(this%n,size(A,2)) 

   x2 = A - matmul(this%h,matmul(transpose(this%h),A))
   Px = this%C .x. x2
   Px = Px - matmul(this%h,matmul(transpose(this%h),Px))
  end function ConsCovar_mtimes_mat

!---------------------------------------------------------------------

! method: 1 -> SS^T
! method: 2 -> Pc

 subroutine locensanalysis(xf,S,Hs,yo,R,lpoints,Hc,xa,Sa,method)
  use matoper
  implicit none
  real, intent(in) :: xf(:), yo(:)
  class(Covar), intent(in) :: R
  real, pointer, intent(in) :: S(:,:), Hc(:,:)
  type(SparseMatrix), intent(in) :: Hs
  real, intent(inout) :: xa(:)
  real, intent(inout), optional :: Sa(:,:)
  integer, intent(in), optional :: method
  procedure(locpoints_) :: lpoints

  type(LocCovar) :: LC
  type(ConsCovar) :: Pc
  integer :: i,m,n,Nens
  integer :: method_ = locensanalysis_SSt

  real, allocatable :: d(:)
  real, allocatable :: Sp(:,:), Sigma(:), KSp(:,:), PaS(:,:)
  real, allocatable :: sqrtPaS(:,:),S2(:,:), A(:,:), LCHc(:,:)

  if (present(method)) method_ = method

  n = size(xf)
  m = size(yo)
  Nens = size(S,2)

  if (size(Hc,1) /= n) then
    write(0,*) 'locensanalysis: parameter Hc has wrong size'
    stop
  end if

  if (size(xa) /= n) then
    write(0,*) 'locensanalysis: parameter xa has wrong size'
    stop
  end if

  if (present(Sa)) then
    if (size(Sa,1) /= n .or. size(Sa,2) /= Nens) then
      write(0,*) 'locensanalysis: parameter Sa has wrong size'
      stop
    end if
  end if

  ! local covariance
  LC = newLocCovar(S,lpoints)

  ! with conservation 
  Pc = newConsCovar(LC,Hc)

  !write(6,*)'locensanalysis:',__LINE__
  allocate(d(m),LCHc(Pc%n,size(Hc,2)))

  d = yo - (Hs.x.xf)

  ! covariance between conserved quantify and every other grid point
  ! using localized covariance LC
  LCHc = LC.x.Hc

  ! the function K is the Kalman gain
  ! it uses LCHc and Pc
  xa = xf + K(d,n)

  if (present(Sa)) then

    allocate(Sigma(Nens)) 
    allocate(Sp(n,Nens))
    allocate(PaS(Nens,Nens))
    allocate(S2(Nens,Nens))
    allocate(KSp(m,Nens))
    allocate(sqrtPaS(Nens,Nens))

    do i = 1,Nens
      Sp(:,i) = S(:,i) - K(Hs.x.S(:,i),n)
    end do
    
    do i = 1,Nens
      KSp(:,i) = iC(Hs.x.(Pc.x.Sp(:,i)))
    end do

    S2 = matmul(transpose(Sp),Sp)

    if (method_ == locensanalysis_SSt) then
      PaS = matmul(S2,S2) + (KSp.tx.(R.x.KSp))
    elseif (method_ == locensanalysis_Pc) then
      allocate(A(n,Nens))
      A = Sp - (Hs.tx.KSp)
      PaS = (A.tx.(Pc.x.A)) + (KSp.tx.(R.x.KSp))
      deallocate(A)
    else
      write(0,*) 'unknown method ',method_
      stop
    end if
    PaS = (PaS + transpose(PaS)) / 2
    
    sqrtPaS = sqrtm(PaS)
    
    Sa = Sp.x.(inv(S2).x.sqrtPaS)


    deallocate(Sigma) 
    deallocate(Sp)
    deallocate(PaS)
    deallocate(S2)
    deallocate(KSp)
    deallocate(sqrtPaS)
  
  end if

  deallocate(d)
contains


 ! computes (H P H' + R) * y
 ! where P is a local covariance matrix with constrain

 function fun_Cx(y) result(Cy)
  real, intent(in) :: y(:)
  real :: Hty(Pc%n)
  real :: Cy(size(y))

!  if (.true.) then  
  if (.false.) then  
    Cy = locenscovx(Pc,Hs,R,y)
  else
    ! optimized version
   Hty = y.x.Hs
!   Cy = (H .x. (Pc.x.Hty))

   Cy = loccovar_project(LC,Hs,y)
   Cy = Cy - (Hs.x.(Hc.x.(LCHc.tx.Hty)))
   Cy = Cy - (Hs.x.(LCHc.x.(Hc.tx.Hty)))
   Cy = Cy + (Hs.x.(Hc.x.(Hc.tx.(LCHc.x.(Hc.tx.Hty)))))
   
   Cy = Cy + (R.x.y)
  end if
 end function fun_Cx
 
 function iC(x) result(y)
  real, intent(in) :: x(:)
  real :: y(size(x))

  y = pcg(fun_Cx,x)
 end function iC

 ! Kalman gain
 function K(d,n) result(x)
  real, intent(in) :: d(:)
  integer :: n
  real :: x(n)

  x = Pc.x.(iC(d).x.Hs)
 end function K
 

 function locenscovx(Pc,H,R,y) result(Cy)
  use matoper
  class(ConsCovar), intent(in) :: Pc
  type(SparseMatrix) :: H
  class(Covar) :: R
  real :: y(:)
  real :: Cy(size(y,1))
  real :: Hty(Pc%n)

  Hty = y.x.H
  Cy = (H .x. (Pc.x.Hty))  + (R.x.y)
 end function locenscovx

 end subroutine locensanalysis

 end module covariance
