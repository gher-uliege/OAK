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


module matoper2

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


interface assert
  module procedure assert_bool
  module procedure assert_scal
  module procedure assert_vec
  module procedure assert_mat
end interface 



 contains
  subroutine assert_bool(cond,msg)
   logical, intent(in) :: cond
   character(len=*) :: msg
   
   if (cond) then
     write(6,*) msg, ': OK '
   else
     write(6,*) msg, ': FAIL '
     stop
   end if
  end subroutine assert_bool


  subroutine assert_scal(found,expected,tol,msg)
   real, intent(in) :: found, expected, tol
   character(len=*) :: msg
   
   real :: maxdiff
   maxdiff = abs(found - expected)
   
   if (maxdiff < tol) then
     write(6,*) msg, ': OK '
   else
     write(6,*) msg, ': FAIL ', maxdiff
     write(6,*) 'found ',found
     write(6,*) 'expected ',expected

     stop
   end if


  end subroutine assert_scal

  subroutine assert_vec(found,expected,tol,msg)
   real, intent(in) :: found(:), expected(:), tol
   character(len=*) :: msg
   
   real :: maxdiff
   maxdiff = maxval(abs(found - expected))
   
   if (maxdiff < tol) then
     write(6,*) msg, ': OK '
   else
     write(6,*) msg, ': FAIL ', maxdiff
     write(6,*) 'found ',found
     write(6,*) 'expected ',expected

     stop
   end if


  end subroutine assert_vec

  subroutine assert_mat(found,expected,tol,msg)
   real, intent(in) :: found(:,:), expected(:,:), tol
   character(len=*) :: msg
   
   real :: maxdiff
   maxdiff = maxval(abs(found - expected))
   
   if (maxdiff < tol) then
     write(6,*) msg, ': OK '
   else
     write(6,*) msg, ': FAIL ', maxdiff
!     write(6,*) 'found ',found
!     write(6,*) 'expected ',expected

     stop
   end if


  end subroutine assert_mat


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


  subroutine test_locfun()   
   call assert(locfun(0.), 1., 1e-7, 'locfun (1)')
   call assert(locfun(0.4), 0.783573333333333, 1e-6, 'locfun (2)')
   call assert(locfun(1.5), 0.0164930555555556, 1e-6, 'locfun (3)')
   call assert(locfun(2.5), 0., 1e-7, 'locfun (4)')
  end subroutine test_locfun

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
      end do    
    end do
  end function 


  ! compute H * P * H' : matmul(matul(this, transpose(H),H)
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


  function pcg(fun,b,x0,tol,maxit,pc,nit) result(x)

   interface 
     function fun(x) result(y)
     real, intent(in) :: x(:)
     real :: y(size(x))
    end function fun
  end interface

   interface 
     function pc(x) result(y)
     real, intent(in) :: x(:)
     real :: y(size(x))
    end function pc
  end interface

  optional pc
!   class(Covar), intent(in) :: A
   real, intent(in) :: b(:)
   real :: x(size(b))
   real, intent(in), optional :: x0(:)
   real, intent(in), optional :: tol
   integer, intent(in), optional :: maxit
   integer, intent(out), optional :: nit

!   class(Covar) :: pc
   real :: tol_, zr_old, zr_new
   real, pointer :: alpha(:), beta(:)
   integer :: maxit_
   integer :: n, k
   real :: tol2
   real, dimension(size(x)) :: Ap, p, r, r_old, z

   n = size(b)
   ! default parameters
   maxit_ = min(n,100)
   tol_ = 1e-6
   
   if (present(tol)) tol_ = tol
   if (present(maxit)) maxit_ = maxit


   allocate(alpha(maxit_+1),beta(maxit_+1))
   
   ! initial guess
   if (present(x0)) then
     x = x0
   else
     ! random initial vector
     x = reshape(randn(n,1),(/ n /))
   end if   

   tol2 = tol_**2

   ! gradient at initial guess
   r = b - fun(x)
   
   ! quick exit
   if (sum(r**2) < tol2) then
     if (present(nit)) nit = k
     return
   endif

   
   ! apply preconditioner
   
   if (present(pc)) then
     z = pc(r)
   else
     z = r
   endif

   
   ! first search direction == gradient
   p = z

   ! compute: r' * inv(M) * z (we will need this product at several
   ! occasions)
   
   zr_old = sum(r*z)

   ! r_old: residual at previous iteration
   r_old = r
   
   do k=1,maxit_    
!     write(6,*) ' k',k,sum(r*r),maxit_,tol2
     ! compute A*p
     Ap = fun(p)
     !maxdiff(A*p,Ap)
     
     ! how far do we need to go in direction p?
     ! alpha is determined by linesearch
     
     ! alpha z'*r / (p' * A * p)
     alpha(k) = zr_old / ( sum(p * Ap))
     
     ! get new estimate of x
     x = x + alpha(k)*p
     
     ! recompute gradient at new x. Could be done by
     ! r = b-fun(x)
     ! but this does require an new call to fun
     r = r - alpha(k)*Ap
     
     ! apply pre-conditionner
     if (present(pc)) then
       z = pc(r)
     else
       z = r
     endif
    
     
     zr_new = sum(r*z)
     
     if (sum(r*r) < tol2) then
       if (present(nit)) nit = k
       exit
     endif
     
     !Fletcher-Reeves
     beta(k+1) = zr_new / zr_old
     !Polak-Ribiere
     !beta(k+1) = r'*(r-r_old) / zr_old
     !Hestenes-Stiefel
     !beta(k+1) = r'*(r-r_old) / (p'*(r-r_old))
     !beta(k+1) = r'*(r-r_old) / (r_old'*r_old)
     
     
     ! norm(p)
     p = z + beta(k+1)*p
     zr_old = zr_new
     r_old = r
   enddo


  end function pcg


  subroutine test_pcg()
   real :: x(3)
   integer :: nit
   
   x = pcg(test_fun,(/ 1.,2.,3. /),nit=nit)
   call assert(x,(/ 2.,4.,6. /),1e-6,'conjugate gradient')

   contains

    function test_fun(x) result(y)
     real, intent(in) :: x(:)
     real :: y(size(x))
     y = x/2.
    end function test_fun

    
  end subroutine test_pcg



 subroutine test_chol()
  real :: A(3,3), U(3,3)

  A = reshape([2.,1.,1., 1.,2.,1., 1.,1.,2.],[3,3])
  U = chol(A)
  call assert(matmul(transpose(U),U), &
      A, &
      1e-6,'Cholesky factorization')
 end subroutine test_chol


 subroutine test_sqrtm()
  real :: A(3,3)
  real :: S(size(A,1),size(A,1))

  A = reshape([2.,1.,1., 1.,2.,1., 1.,1.,2.],[3,3])
  S = sqrtm(A)

  call assert(matmul(S,S), &
      A, &
      1e-6,'sqrtm factorization')
 end subroutine test_sqrtm


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

  real, allocatable :: tmp(:), d(:)
  real, allocatable :: Sp(:,:), Sigma(:), KSp(:,:), PaS(:,:)
  real, allocatable :: sqrtPaS(:,:),S2(:,:), A(:,:)

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

  allocate(tmp(m),d(m))

  d = yo - (Hs.x.xf)

  tmp = pcg(fun_Cx,yo - (Hs.x.xf))

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
      write(0,*) 'method'
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

  deallocate(tmp,d)
contains



 function fun_Cx(x) result(y)
  real, intent(in) :: x(:)
  real :: y(size(x))

!  write(6,*) 'call fun_Cx'
  y = locenscovx(Pc,Hs,R,x)
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


 subroutine test_covar
  use matoper
  type(DiagCovar) :: C
  !class(DiagCovar), allocatable :: C
  real :: mdiag(2) = (/ 2.,3. /)

  !allocate(C)
  C = newDiagCovar(mdiag)
  !C = DiagCovar(mdiag)
  
  call innersub(C)

  contains 
   subroutine innersub(B)
    class(Covar) :: B
    real :: C(2,2)

    C = 2. *  eye(2)

    call B%print()

    call assert(B%mulmat(C), &
       reshape([4.,0.,0.,6.],[2,2]), &
       1e-7, 'DiagCovar*Matrix (method)')

    call assert(B.x.(C), &
       reshape([4.,0.,0.,6.],[2,2]), &
       1e-7, 'DiagCovar*Matrix (.x. operator)')

    call assert(B * (C), &
       reshape([4.,0.,0.,6.],[2,2]), &
       1e-7, 'DiagCovar*Matrix (* operator)')


   end subroutine innersub
  end subroutine test_covar


  ! sort vector A in-place
  ! optional argument ind is the sort index
  ! such that
  ! sortedA = A
  ! call sort(sortedA,ind)
  ! all(sortedA,A(ind)) is true

#define DATA_TYPE integer 

  subroutine sort(A,ind)
   DATA_TYPE, intent(inout), dimension(:) :: A
   integer, intent(out), dimension(size(A)), optional :: ind
   
   integer :: sort_index(size(A)), i
   
   sort_index = [(i, i=1,size(A))]
   
   call sort_(A,sort_index)
   if (present(ind)) ind = sort_index

  contains
   recursive subroutine sort_(A,ind)
    DATA_TYPE, intent(inout), dimension(:) :: A
    integer, intent(inout), dimension(:) :: ind
    integer :: iq

    if (size(A) > 1) then
      call sort_partition(A, ind, iq)
      call sort_(A(:iq-1),ind(:iq-1))
      call sort_(A(iq:),ind(iq:))
    end if
   end subroutine sort_

   subroutine sort_partition(A, ind, marker)
    DATA_TYPE, intent(inout), dimension(:) :: A
    integer, intent(inout), dimension(:) :: ind
    integer, intent(out) :: marker
    integer :: i, j, tempind
    DATA_TYPE :: temp
    DATA_TYPE :: x      ! pivot point

    x = A(1)
    i = 0
    j = size(A) + 1
    
    do
      j = j-1
      do
        if (A(j) <= x) exit
        j = j-1
      end do
      
      i = i+1
      do
        if (A(i) >= x) exit
        i = i+1
      end do
      

      if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp

        tempind = ind(i)
        ind(i) = ind(j)
        ind(j) = tempind        
      else if (i == j) then
        marker = i+1
        return
      else
        marker = i
        return
      end if
    end do
    
   end subroutine sort_partition

  end subroutine sort


  subroutine test_sort
   implicit none
   integer, parameter :: n = 10
   integer :: ind(n)
   integer, dimension(1:n) :: A = &  
        (/0, 50, 20, 25, 90, 10, 5, 20, 99, 75/), sortedA
   
   sortedA = A
   call sort(sortedA,ind)
   call assert(all(sortedA(1:n-1) <= sortedA(2:n)),'quick sort (1)')
   call assert(all(sortedA == A(ind)),'quick sort (2)')
  end subroutine test_sort


  subroutine unique(A,n,ind)
   DATA_TYPE, intent(inout) :: A(:)
   integer, intent(out) :: n
   integer, intent(out), dimension(size(A)), optional :: ind
   
   integer :: unique_index(size(A)), i
   
   unique_index = [(i, i=1,size(A))]
   
   call unique_(A,n,unique_index)
   if (present(ind)) ind = unique_index

  contains
   subroutine unique_(A,n,ind)
   DATA_TYPE, intent(inout) :: A(:)
   integer, intent(out) :: n
   integer, intent(out), dimension(size(A)), optional :: ind

   integer :: i
   
   call sort(A,ind)
   n = 1
   do i = 1,size(A)-1
     A(n) = A(i)
     ind(n) = ind(i)

     if (A(i) /= A(i+1)) then
       n = n+1
     end if
   end do

   A(n) = A(size(A))
   ind(n) = ind(size(A))

  end subroutine unique_
 end subroutine unique
  subroutine test_unique
   implicit none
   integer, parameter :: n = 10
   integer :: ind(n), i, nu
   integer, dimension(1:n) :: A = &  
        (/0, 20, 20, 25, 90, 10, 5, 20, 90, 75/), sortedA, uA, c
   
   uA = A
   call unique(uA,nu,ind)

   ! all elements in A must be one time in uA
   do i = 1,n
     if (count(A(i) == uA(1:nu)) /= 1) then
       write(6,*) 'unique (1): FAILED'
       stop
     end if
   end do

   write(6,*) 'unique (1): OK'

   call assert(all(uA(1:nu) == A(ind(1:nu))),'unique (2)')

  end subroutine test_unique

  


 end module matoper2
