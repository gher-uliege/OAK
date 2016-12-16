! include the fortran preprocessor definitions
#include "../ppdef.h"

! cd ~/Assim/OAK-nonDiagR &&  make test/test_covariance && test/test_covariance

#define PROFILE
module test_suite

contains

 function testens(n,Nens) result(E)
  integer :: n,Nens
  integer :: i,j 
  real :: E(n,Nens)
  real, parameter :: pi = 3.141592654

  do j = 1,Nens
    do i = 1,n
      E(i,j) = sin( pi * j * real(i-1) / real(n-1)) / j
    end do
  end do
 end function testens

 subroutine run_test(sz)
  use matoper
  use covariance

  integer, intent(in) :: sz(:)

  integer :: n
  real :: len = 3
  integer :: Nens

  type(LocCovar) :: LC
  type(ConsCovar) :: Pc
  type(DiagCovar) :: Rc
  real, pointer :: x(:,:)
  integer :: i,j,k,l,m,nnz
  real, pointer :: v(:), v2(:), S(:,:),w(:),Hc(:,:),P(:,:),Pr(:,:)
  real, allocatable :: H(:,:), R(:,:), xf(:), xa(:), yo(:), xa2(:), tmp(:), d(:)
  real, allocatable :: Sp(:,:), Sigma(:), KSp(:,:), PaS(:,:), A(:,:), Sa2(:,:)
  real, allocatable :: Sa(:,:), diagR(:), Pa(:,:), S2(:,:)
  real, allocatable :: sqrtPaS(:,:)
  integer, pointer :: jj(:)
  type(SparseMatrix) :: Hs
  write(6,*) 'Running test with a domain ',sz


  n = product(sz)
  ! substract 1 because of the constrain and one more 
  ! because the last column of testens would be all zero otherwise
  ! sin(pi * (i-1)) = 0

  Nens = min(n-2,20)

  allocate(x(n,size(sz)))
  allocate(v(n))
  allocate(v2(n))
  allocate(S(n,Nens))
  allocate(w(n),jj(n))


  v = 0
  v(1) = 1
  S = 1

  S = randn(n,Nens)
  S = testens(n,Nens)

  allocate(Hc(n,1))
  Hc = 1
  Hc = Hc/sqrt(sum(Hc**2))

  ! apply constaint on S
  do i = 1,Nens
    do j = 1,size(Hc,2)
      S(:,i) = S(:,i) - Hc(:,j) * sum(Hc(:,j) * S(:,i))
    end do
  end do

  !write(6,*) 'S',S(5,:)
  do i = 1,size(Hc,2)
    call assert(maxval(abs(matmul(transpose(Hc),S))), &
         0.,6e-5,'verifying constraint (Sf)')
  end do

  l = 0

  if (size(sz) == 2) then
    do j = 1,sz(2)
      do i = 1,sz(1)
        l = l+1
        x(l,1) = i
        x(l,2) = j
      end do
    end do
  else
    do k = 1,sz(3)
      do j = 1,sz(2)
        do i = 1,sz(1)
          l = l+1
          x(l,1) = i
          x(l,2) = j
          x(l,3) = k
        end do
      end do
    end do
  end if

  call locpoints(1,x,len,nnz,jj,w)
!  write(6,*) 'j, w',jj
!  write(6,*) 'j, w',w

  call locpoints(1,x,len,nnz,jj,w,[2])
!  write(6,*) 'w',jj(1:nnz)
!  write(6,*) 'j',w(1:nnz)

  allocate(P(n,n))
  allocate(Pr(n,n))

  ! full P
  P = matmul(S,transpose(S))

  ! local P
  do j=1,n
    do i=1,n
      P(i,j) = P(i,j) * locfun(cdist(x(i,:),x(j,:))/len)
    end do
  end do


  ! local covariance
  LC = newLocCovar(S,lpoints)
  call assert(LC.x.v,matmul(P,v),5e-5,'local covariance')

  ! with conservation 
  Pc = newConsCovar(LC,Hc)

  Pr = eye(n) - matmul(Hc,transpose(Hc))
  P = matmul(Pr,matmul(P,Pr))

  call assert(Pc.x.v,matmul(P,v),5e-5,'constrained local covariance')


  m = 3
  allocate(H(m,n),R(m,m),xf(n),xa(n),xa2(n),yo(m))
  allocate(Pa(n,n))

  R = eye(m)
  H = 0
  H(1:m,1:m) = eye(m)
  yo = 1
  xf = 0

  xa = xf + matmul(matmul(P,transpose(H)), &
       matmul(inv(matmul(matmul(H,P),transpose(H)) + R), &
       (yo - matmul(H,xf))))

  ! write(6,*) 'xa',xa

  Pa = P - ((P.xt.H).x.(inv(((H.x.P).xt.H) + R).x.(H.x.P)))

  xa2 = xf + ((P.xt.H).x.(inv(((H.x.P).xt.H) + R).x.(yo - (H.x.xf))))
  call assert(xa,xa2,1e-6,'analysis')

  Hs = sparse(H)

  call assert(((H.x.P).xt.H), &
       ((Hs.x.P).xt.Hs), &
       1e-6,'covariance at observations')

  xa2 = xf + ((P.xt.Hs).x.(inv(((Hs.x.P).xt.Hs) + R).x.(yo - (Hs.x.xf))))
  call assert(xa,xa2,1e-6,'analysis using sparse matrix')

  allocate(diagR(m))
  diagR = 1
  Rc = newDiagCovar(diagR)

  allocate(tmp(m),d(m))

  d = yo - (Hs.x.xf)

  ! call assert( &
  !      (((Hs.x.P).xt.Hs) + R).x.(yo - (Hs.x.xf)), &
  !      (((Hs.x.P).xt.Hs).x.d) + (R.x.d), &
  !      1e-6,'Cx')

  call assert( &
       matmul(transpose(H),d), &
       H.tx.d, &
       1e-6,'sparse mat*vector')

  call assert( &
       matmul(transpose(H),d), &
       d.x.Hs, &
       1e-6,'vector*sparse mat')

  tmp = fun_Cx(d)
  call assert(tmp, &
       (((Hs.x.P).xt.Hs) + R).x.(yo - (Hs.x.xf)), &
       5e-5,'Cx')

  tmp = fun_Cx2(d)
  call assert(tmp, &
       (((Hs.x.P).xt.Hs) + R).x.(yo - (Hs.x.xf)), &
       5e-5,'Cx (optimized)')

  tmp = pcg(fun_Cx,yo - (Hs.x.xf))
  ! tmp = fun_Cx(yo - (Hs.x.xf))
  call assert(tmp, &
       inv(((Hs.x.P).xt.Hs) + R).x.d, &
       2e-6,'solving system using pcg')

  call assert(iC(d), &
       inv(((Hs.x.P).xt.Hs) + R).x.d, &
       2e-6,'solving system using pcg (2)')


  xa2 = xf + ((P.xt.Hs).x.tmp) 
  call assert(xa,xa2,2e-5,'analysis using sparse matrix, pcg')

  xa2 = xf + (Pc.x.(tmp.x.Hs)) 
  call assert(xa,xa2,2e-5,'analysis using sparse matrix, pcg, ConsCovar')

  xa2 = xf + KG(d,n)
  call assert(xa,xa2,2e-5,'analysis using KG')

  allocate(Sp(n,Nens),S2(Nens,Nens),Sigma(Nens)) 
  allocate(PaS(Nens,Nens))
  allocate(A(n,Nens))
  allocate(KSp(m,Nens))
  allocate(sqrtPaS(Nens,Nens))

  do i = 1,Nens
    Sp(:,i) = S(:,i) - KG(Hs.x.S(:,i),n);
  end do

  do i = 1,Nens
    KSp(:,i) = iC(Hs.x.(Pc.x.Sp(:,i)));
  end do

  S2 = matmul(transpose(Sp),Sp)
  PaS = matmul(S2,S2) + (KSp.tx.(R.x.KSp))
    
  PaS = (PaS + transpose(PaS)) / 2
    
  sqrtPaS = sqrtm(PaS)
 
  allocate(Sa(n,Nens)) 
  Sa = Sp.x.(inv(S2).x.sqrtPaS)

!  write(6,*) 'here ',maxval(matmul(sqrtPaS,sqrtPaS) - PaS),maxval(PaS)
  call assert(matmul(sqrtPaS,sqrtPaS), &
       PaS, &
       1e-5 * maxval(abs(PaS)),'sqrtm factorization')

  allocate(Sa2(n,Nens)) 

  write(6,*) 'method ',locensanalysis_SSt
  call locensanalysis(xf,S,Hs,yo,Rc,lpoints,Hc,xa2,Sa2,method=locensanalysis_SSt)
    
  call assert(xa,xa2,2e-5,'analysis using locensanalysis (xa)')
  call assert(Sa.xt.Sa,Sa2.xt.Sa2,1e-4,'analysis using locensanalysis (Sa)')
    
  do i = 1,size(Hc,2)      
    call assert(sum(Hc(:,i) * (xf-xa)),0.,6e-5,'verifying constraint (xa)')
    call assert(maxval(abs(matmul(transpose(Hc),Sa))), &
         0.,6e-5,'verifying constraint (Sa)')
  end do

  write(6,*) 'method ',locensanalysis_Pc
  call locensanalysis(xf,S,Hs,yo,Rc,lpoints,Hc,xa2,Sa2,method=locensanalysis_Pc)
    
  call assert(xa,xa2,2e-5,'analysis using locensanalysis (xa)')

  do i = 1,size(Hc,2)      
    call assert(sum(Hc(:,i) * (xf-xa)),0.,6e-5,'verifying constraint (xa)')
    call assert(maxval(abs(matmul(transpose(Hc),Sa))), &
         0.,6e-5,'verifying constraint (Sa)')
  end do


  deallocate(x)
  deallocate(v)
  deallocate(v2)
  deallocate(S)
  deallocate(w,jj)

  deallocate(P)
  deallocate(Pr)
  deallocate(Hc)
  deallocate(H,R,xf,xa,xa2,yo)
  deallocate(tmp,d)


  deallocate(Sp)
  deallocate(S2,Sigma) 
  deallocate(PaS)
  deallocate(A)
  deallocate(KSp)
  deallocate(sqrtPaS)

 contains



  function fun_Cx(x) result(y)
   real, intent(in) :: x(:)
   real :: y(size(x))
   y = locenscovx(Pc,Hs,Rc,x)
  end function fun_Cx

  ! optimized version
  function fun_Cx2(x) result(y)
   real, intent(in) :: x(:)
   real :: y(size(x))

   y = locenscovx2(Pc,Hs,Rc,x)
  end function fun_Cx2


  function iC(x) result(y)
   real, intent(in) :: x(:)
   real :: y(size(x))

   y = pcg(fun_Cx,x)
  end function iC

  ! Kalman gain
  function KG(d,n) result(x)
   real, intent(in) :: d(:)
   integer :: n
   real :: x(n)

   !  x = ((P.xt.Hs).x.iC(d))  
   x = Pc.x.(iC(d).x.Hs)
  end function KG


  function locenscovx(Pc,H,Rc,y) result(Cy)
   use matoper
   class(ConsCovar), intent(in) :: Pc
   type(SparseMatrix) :: H
   class(Covar) :: Rc
   real :: y(:)
   real :: Cy(size(y,1))
   real :: Hty(Pc%n)

   Hty = y.x.H
   Cy = (H .x. (Pc.x.Hty))  + (Rc.x.y)
  end function locenscovx

  function locenscovx2(Pc,Hs,Rc,y) result(Cy)
   use matoper
   class(ConsCovar), intent(in) :: Pc
   type(SparseMatrix) :: Hs
   class(Covar) :: Rc
   real :: y(:)
   real :: Cy(size(y,1))
   real :: Hty(Pc%n), A(Pc%n,size(Hc,2))

   Hty = y.x.Hs
!   Cy = (H .x. (Pc.x.Hty))

   A = LC.x.Hc
   Cy = loccovar_project(LC,Hs,y)

   Cy = Cy - (Hs.x.(Hc.x.(A.tx.Hty)))
   Cy = Cy - (Hs.x.(A.x.(Hc.tx.Hty)))
   Cy = Cy + (Hs.x.(Hc.x.(Hc.tx.(A.x.(Hc.tx.Hty)))))

   
   Cy = Cy + (Rc.x.y)

  end function locenscovx2


  subroutine lpoints(i,nnz,j,w,onlyj)
   integer, intent(in) :: i
   integer, intent(out) :: nnz,j(:)
   real, intent(out) :: w(:)
   integer, optional, intent(in) :: onlyj(:)  

   call locpoints(i,x,len,nnz,j,w,onlyj)  
  end subroutine lpoints

 end subroutine run_test



 subroutine run_test_large(sz,computeSa)
  use covariance
  use ndgrid, only: regulargrid, init_regulargrid, near_regulargrid
  integer, intent(in) :: sz(:)
  logical :: computeSa

  integer :: n
  real :: len = 3
  integer :: Nens = 20

  class(ConsCovar), allocatable :: Pc
  class(DiagCovar), allocatable :: Rc
  real, pointer :: x(:,:)
  integer :: i,j,k,l,m
  real, pointer :: S(:,:), Hc(:,:)
  real, allocatable :: xf(:), xa(:), yo(:), xa2(:)
  real, allocatable :: Sa(:,:), diagR(:), Sa2(:,:)
  type(SparseMatrix) :: Hs
  type(regulargrid) :: g

# ifdef PROFILE
  real(8) :: cputime(2)
# endif

  write(6,*) 'Running test with a domain ',sz
  n = product(sz)

  allocate(x(n,size(sz)))
  allocate(S(n,Nens))

  S = 1
!  S = randn(n,Nens)
  S = testens(n,Nens)

  allocate(Pc)
  allocate(Hc(n,1))
  Hc = 1
  Hc = Hc/sqrt(sum(Hc**2))

  ! apply constaint on S
  do i = 1,Nens
    do j = 1,size(Hc,2)
      S(:,i) = S(:,i) - Hc(:,j) * sum(Hc(:,j) * S(:,i))
    end do
  end do

  do i = 1,size(Hc,2)
    call assert(maxval(abs(matmul(transpose(Hc),S))), &
         0.,6e-5,'verifying constraint (Sf)')
  end do

  l = 0

  if (size(sz) == 2) then
    do j = 1,sz(2)
      do i = 1,sz(1)
        l = l+1
        x(l,1) = i
        x(l,2) = j
      end do
    end do
  else
    do k = 1,sz(3)
      do j = 1,sz(2)
        do i = 1,sz(1)
          l = l+1
          x(l,1) = i
          x(l,2) = j
          x(l,3) = k
        end do
      end do
    end do
  end if

  call init_regulargrid(g,sz,0.*sz,0.*sz+1)

  m = n/2

  allocate(Rc)
  allocate(diagR(m))
  diagR = 1
  call Rc%init(diagR)

  allocate(Hs%i(m),Hs%j(m),Hs%s(m))
  Hs%m = m
  Hs%n = n
  Hs%i = (/ (i,i=1,m) /)
  Hs%j = (/ (i,i=1,m) /)
  Hs%s = 1
  Hs%nz = m

  allocate(xf(n),xa(n), xa2(n), yo(m))
  xf = 0
  yo = 1

  xf = (/ (2*i,i=1,n) /)

  call assert(Hs .x. xf, (/ (2.*i,i=1,m) /) ,6e-5,'verifying obsoperator')

  allocate(Sa(n,Nens),Sa2(n,Nens)) 

  if (computeSa) then
# ifdef PROFILE
  call cpu_time(cputime(1))
# endif    
    call locensanalysis(xf,S,Hs,yo,Rc,lpoints,Hc,xa,Sa)
# ifdef PROFILE
  call cpu_time(cputime(2))
  write(stdout,*) 'locensanalysis CPU time  ',cputime(2)-cputime(1)
# endif    



! # ifdef PROFILE
!   call cpu_time(cputime(1))
! # endif    
!     call locensanalysis(xf,S,Hs,yo,Rc,lpoints_regulargrid,Hc,xa2,Sa2)

! # ifdef PROFILE
!   call cpu_time(cputime(2))
!   write(stdout,*) 'locensanalysis rg CPU time  ',cputime(2)-cputime(1)
! # endif    

!       call assert(xa,xa2,6e-5,'checking xa ')
!       call assert(Sa,Sa2,6e-5,'checking Sa ')

  else
    call locensanalysis(xf,S,Hs,yo,Rc,lpoints,Hc,xa)
  end if


  do i = 1,size(Hc,2)      
    call assert(sum(Hc(:,i) * (xf-xa)),0.,6e-5,'verifying constraint (xa)')

    if (computeSa) then
      call assert(maxval(abs(matmul(transpose(Hc),Sa))), &
           0.,6e-5,'verifying constraint (Sa)')
    end if
  end do



  deallocate(x)
  deallocate(S)
  deallocate(Hc)
  deallocate(xf,xa,yo,Sa)


 contains


  subroutine lpoints(i,nnz,j,w,onlyj)
   integer, intent(in) :: i
   integer, intent(out) :: nnz,j(:)
   real, intent(out) :: w(:)
   integer, optional, intent(in) :: onlyj(:)  

   call locpoints(i,x,len,nnz,j,w,onlyj)  
  end subroutine lpoints


  subroutine lpoints_regulargrid(i,nnz,j,w,onlyj)
   integer, intent(in) :: i
   integer, intent(out) :: nnz,j(:)
   real, intent(out) :: w(:)
   integer, optional, intent(in) :: onlyj(:)  
   integer :: l

   ! write(6,*) 'i',i
   ! call locpoints(i,x,len,nnz,j,w,onlyj)  
   ! write(6,*) 'locpoints',sum(j(1:nnz)),sum(w(1:nnz))
   ! write(6,*) 'locpoints',nnz
   ! write(6,*) 'len',len
   call near_regulargrid(g,x(i,:),cdist,2*len,j,w,nnz)

   do l = 1,nnz
     w(l) = locfun(w(l)/len)
   end do

   !write(6,*) 'locpoints',j(1:nnz),nnz
   !write(6,*) 'locpoints',sum(j(1:nnz)),sum(w(1:nnz))
   !stop
  end subroutine lpoints_regulargrid

 end subroutine run_test_large

end module test_suite



program test
 use test_suite
 use covariance

 if (kind(real(1.)) <= 4) then
   write(0,*) 'Error: tests require at least double precision'
   stop
 end if

 call test_DiagCovar
 call test_SMWCovar
 call test_DCDCovar
 call test_locfun 

 ! same results as matlab code test_covariance_fortran
 call run_test([3,2]) ! ok

 call run_test([5,5])
 call run_test([5,5,10]) ! ok
! call run_test_large([5,5,10],.false.) ! ok in double precision

 !call run_test_large([5,5,10],.true.) ! ok in double precision
 !call run_test_large([50,5,10],.true.) ok in double precision

! call run_test_large([30,30,20],.false.) ! ok
! call run_test_large([80,80,30],.false.)

contains

!----------------------------------------------------------------------

  subroutine test_locfun()   
   call assert(locfun(0.), 1., 1e-7, 'locfun (1)')
   call assert(locfun(0.4), 0.783573333333333, 1e-6, 'locfun (2)')
   call assert(locfun(1.5), 0.0164930555555556, 1e-6, 'locfun (3)')
   call assert(locfun(2.5), 0., 1e-7, 'locfun (4)')
  end subroutine test_locfun

!----------------------------------------------------------------------


  subroutine test_DiagCovar
   use matoper
   implicit none
   integer, parameter :: m = 10
   real :: C(m), F(m,m)
   type(DiagCovar) :: Cov
   integer :: i

   write(6,*) '# Testing diagonal covariance matrix'
   C = [(mod(i,5)+1., i=1,m)]
   F = diag(C)

   call DiagCovar_init(Cov,C)

   call test_covariance_matrix(m,F,Cov)


  end subroutine test_DiagCovar


  subroutine test_SMWCovar
   use matoper
   implicit none
   integer, parameter :: m = 10, N = 2
   real :: C(m), B(m,N), F(m,m)
   type(SMWCovar) :: Cov
   integer :: i

   write(6,*) '# Testing SMWCovar matrix'

   C = [(mod(i,5)+1., i=1,m)]
   B = reshape([(mod(i,5)+1., i=1,m*N)],[m,N])
   
   call SMWCovar_init(Cov,C,B)

   F = matmul(B,transpose(B)) + diag(C)

   call test_covariance_matrix(m,F,Cov)
  end subroutine test_SMWCovar

!---------------------------------------------------------------------

  subroutine test_DCDCovar
   use matoper
   implicit none
   integer, parameter :: m = 10, N = 2
   real :: C(m), B(m,N), F(m,m), D(m)
   type(SMWCovar) :: innerCov
   type(DCDCovar) :: Cov
   integer :: i

   write(6,*) '# Testing DCDCovar matrix'

   C = [(mod(i,5)+1., i=1,m)]
   B = reshape([(mod(i,5)+1., i=1,m*N)],[m,N])
   D = [(mod(10*i,5)+1., i=1,m)]
   
   call SMWCovar_init(innerCov,C,B)
   call DCDCovar_init(Cov,D,innerCov)

   F = matmul(B,transpose(B)) + diag(C)
   F = matmul(diag(1./D),matmul(F,diag(1./D)))

   call test_covariance_matrix(m,F,Cov)
  end subroutine 

!---------------------------------------------------------------------

  subroutine test_covariance_matrix(m,F,Cov)
   use matoper
   implicit none
   integer, intent(in)        :: m
   class(Covar), intent(in) :: Cov
   real, intent(in)           :: F(m,m)

   real :: y(m), z(m), z_ref(m)
   real :: A(m,m), D(m,m), D_ref(m,m)
   logical :: mask(m)
   integer :: i,i1,i2

   ! compare full matrix

   call assert(F,Cov%full(), 1e-7, 'compare to full matrix')

   ! multiply by a vector

   y = [(mod(i,2)+1., i=1,m)]
   z = matmul(Cov,y)
   z_ref = matmul(F,y)
   call assert(z,z_ref, 1e-7, 'multiply vector (matmul)')

   ! multiply by a vector

   y = [(mod(i,2)+1., i=1,m)]
   z = Cov.x.y
   z_ref = matmul(F,y)
   call assert(z,z_ref, 1e-7, 'multiply vector (.x.)')

   ! multiply by a matrix
   
   A = reshape([(mod(i,10), i=1,m*m)],[m,m])
   D = matmul(Cov,A)
   D_ref = matmul(F,A)
   call assert(D,D_ref, 1e-7, 'multiply matrix (matmul)')

   ! multiply by a matrix
   
   A = reshape([(mod(i,10), i=1,m*m)],[m,m])
   !D = mat_mul_covar(A,Cov)
   !D = matmul(A,Cov)
   D = A.x.Cov
   D_ref = matmul(A,F)
   call assert(D,D_ref, 1e-7, 'multiply matrix (matmul, 2)')

   ! multiply by a matrix
   
   A = reshape([(mod(i,10), i=1,m*m)],[m,m])
   D = Cov.x.A
   D_ref = matmul(F,A)
   call assert(D,D_ref, 1e-7, 'multiply matrix (.x.)')

   ! solve for a vector


   y = [(mod(i,2)+1., i=1,m)]
   z = Cov%mldivide(y)
   z_ref = matmul(inv(F),y)
   call assert(z,z_ref, 1e-7, 'solve for a vector')

   ! solve for a matrix

   A = reshape([(mod(i,10), i=1,m*m)],[m,m])
   D = Cov%mldivide(A)
   D_ref = matmul(inv(F),A)
   call assert(D,D_ref, 1e-7, 'solve for a matrix')

   ! subset using a mask
   mask = .false.
   mask(2:m) = .true.
   call assert(F(2:m,2:m),Covar_full(Cov%pack(mask)), 1e-7, 'pack')

   ! subset using indices
   i1 = 2
   i2 = m
   call assert(F(i1:i2,i1:i2),Covar_full(Cov%sub(i1,i2)), 1e-7, 'subscripts')

   

  end subroutine test_covariance_matrix


end program test


