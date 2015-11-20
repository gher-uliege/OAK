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


type Covar
  integer         :: n
  real, pointer   :: CM(:,:)
   contains
!        procedure :: init => Covar_init
        procedure :: print => Covar_print        
        procedure :: mulmat => Covar_mulmat
        procedure :: mulvec => Covar_mulvec
end type Covar



type, extends(Covar) :: DiagCovar
  real, pointer :: diag(:)
   contains
        procedure :: init => DiagCovar_init
        procedure :: print => DiagCovar_print
        procedure :: mulmat => DiagCovar_mulmat
        procedure :: mulvec => DiagCovar_mulvec
end type DiagCovar

#ifdef SUPPORT_CON
interface DiagCovar
  module procedure newDiagCovar
end interface DiagCovar
#endif

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
        procedure :: mulvec => loccovar_smult_vec
        procedure :: mulmat => loccovar_smult_mat
        procedure :: print => lc_print        
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
        procedure :: initialize => consInitialize
        procedure :: mulvec => conscovar_smult_vec
        procedure :: mulmat => conscovar_smult_mat
end type ConsCovar




interface operator(.x.)
  module procedure              &
   covar_mult_mat, &
   covar_mult_vec
end interface 

interface operator(*)
  module procedure              &
   covar_mult_mat, &
   covar_mult_vec
end interface 


 contains


  subroutine Covar_init(this,C)
    class(Covar) :: this
    real, target :: C(:,:)

    this%n = size(C,1)
    this%CM => C
   end subroutine Covar_init

  subroutine Covar_print(this)
    class(Covar), intent(in) :: this
    print *, 'Covar: ', this%n,this%n
   end subroutine Covar_print

  function Covar_mulmat(this,x) result(C)
    class(Covar), intent(in) :: this
    real, intent(in) :: x(:,:)
    real :: C(size(x,1),size(x,2))
    C = matmul(this%CM,x)
   end function Covar_mulmat


  function Covar_mulvec(this,x) result(C)
    class(Covar), intent(in) :: this
    real, intent(in) :: x(:)
    real :: C(size(x,1))

    C = matmul(this%CM,x)
   end function Covar_mulvec


  function Covar_mult_mat(this,x) result (Px)
   class(Covar), intent(in) :: this
   real, intent(in) :: x(:,:)
   real ::  Px(size(x,1),size(x,2))

   Px = this%mulmat(x)
  end function 

  function Covar_mult_vec(this,x) result (Px)
   class(Covar), intent(in) :: this
   real, intent(in) :: x(:)
   real ::  Px(size(x,1))

   Px = this%mulvec(x)
  end function 


  function newDiagCovar(diag) result(this)
    type(DiagCovar) :: this
    real, target :: diag(:)

    this%n = size(diag)
    this%diag => diag 
  end function newDiagCovar


  subroutine DiagCovar_init(this,C)
    class(DiagCovar) :: this
    real, target :: C(:)

    this%n = size(C)
    this%diag => C 
   end subroutine DiagCovar_init

  subroutine DiagCovar_print(this)
    class(DiagCovar), intent(in) :: this
    print *, 'Diag Covar: ', this%n,this%n
   end subroutine DiagCovar_print

  function DiagCovar_mulvec(this,x) result(C)
    class(DiagCovar), intent(in) :: this
    real, intent(in) :: x(:)
    real :: C(size(x,1))
    integer :: i,j

    C = this%diag*x    
   end function DiagCovar_mulvec

  function DiagCovar_mulmat(this,x) result(C)
    class(DiagCovar), intent(in) :: this
    real, intent(in) :: x(:,:)
    real :: C(size(x,1),size(x,2))
    integer :: i,j

    do j=1,size(x,2)
      do i=1,size(x,1)
        C(i,j) = this%diag(i)*x(i,j)
      end do
    end do
   end function DiagCovar_mulmat

!---------------------------------
! LocCovar
!  



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

  subroutine lc_print(this)
    class(LocCovar), intent(in) :: this
    print *, 'lc : r = ', this%n
   end subroutine lc_print



  function loccovar_smult_vec(this,x) result (Px)
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
 end function loccovar_smult_vec

 function loccovar_smult_mat(this,x) result (Px)
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


  end function loccovar_smult_mat


  ! compute H * P * H' y : matmul(matul(this, transpose(H),H) y
  ! where H is sparse

  function loccovar_project(this,H,y) result (HPHy)
   class(LocCovar), intent(in) :: this
   type(SparseMatrix), intent(in) :: H
   real, intent(in) :: y(:)
   
   real ::  PHy(this%n)
   integer :: j(this%n), i, k, nnz, Hjind
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


  subroutine consInitialize(self,C,h)
   class(ConsCovar) :: self
!   class(LocCovar), pointer :: C
   class(Covar), target :: C
   real, pointer :: h(:,:)

   self%n = C%n
   self%C => C
   self%h => h
  end subroutine consInitialize

  function newConsCovar(C, H) result(self)
   type(ConsCovar) :: self
   class(Covar), target :: C
   real, pointer :: H(:,:)

   self%n = C%n
   self%C => C
   self%h => h
  end function newConsCovar


  function conscovar_smult_vec(this,x) result (Px)
   class(ConsCovar), intent(in) :: this
   real, intent(in) :: x(:)
   real ::  Px(size(x,1)), x2(this%n) 

   x2 = x - matmul(this%h,matmul(transpose(this%h),x))
   Px = this%C .x. x2
   Px = Px - matmul(this%h,matmul(transpose(this%h),Px))
  end function conscovar_smult_vec

  function conscovar_smult_mat(this,x) result (Px)
   class(ConsCovar), intent(in) :: this
   real, intent(in) :: x(:,:)
   real ::  Px(size(x,1),size(x,2)), x2(this%n,size(x,2)) 

   x2 = x - matmul(this%h,matmul(transpose(this%h),x))
   Px = this%C .x. x2
   Px = Px - matmul(this%h,matmul(transpose(this%h),Px))
  end function conscovar_smult_mat
 

! method: 1 -> SS^T
! method: 2 -> Pc

 subroutine locensanalysis(xf,S,Hs,yo,R,lpoints,Hc,xa,Sa,method)
  use matoper
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
