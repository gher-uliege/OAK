
module test_suite

contains

 subroutine run_test
  use testmod
  integer :: sz(2) = (/5,5/), n
  real :: len = 3
  integer :: Nens = 20

  class(LocCovar), pointer :: LC
  class(ConsCovar), allocatable :: Pc
  real, pointer :: x(:,:)
  integer :: i,j,l,m,nnz
  real, pointer :: v(:), v2(:), S(:,:),w(:),hc(:,:),P(:,:),Pr(:,:)
  real, allocatable :: H(:,:), R(:,:), xf(:), xa(:), yo(:), xa2(:), tmp(:), d(:)
  real, allocatable :: Sp(:,:), U(:,:),Sigma(:), KU(:,:), PaU(:,:), A(:,:), Sa2(:,:)
  real, allocatable :: Sa(:,:)
  real, allocatable :: sqrtPaU(:,:)
  integer, pointer :: jj(:)
  type(SparseMatrix) :: Hs

  n = product(sz)
  allocate(x(n,2))
  allocate(v(n))
  allocate(v2(n))
  allocate(S(n,Nens))
  allocate(w(n),jj(n))


  v = 0
  v(1) = 1
  S = 1

  S = randn(n,Nens)

  l = 0
  do j = 1,sz(2)
    do i = 1,sz(1)
      l = l+1
      x(l,1) = i
      x(l,2) = j
    end do
  end do

  call locpoints(1,x,len,nnz,jj,w)
  ! write(6,*) 'j, w',jj
  ! write(6,*) 'j, w',w

  call test_locfun 
  call test_pcg

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


  allocate(LC)
  call LC%initialize(S,lpoints)
  call assert(LC.x.v,matmul(P,v),1e-6,'local covariance (2)')

  ! with conservation 

  allocate(Pc,hc(n,1))
  hc = 1
  hc = hc/sqrt(sum(hc**2))
  call Pc%initialize(LC,hc)

  Pr = eye(n) - matmul(hc,transpose(hc))
  P = matmul(Pr,matmul(P,Pr))

  call assert(Pc.x.v,matmul(P,v),1e-5,'constrained local covariance')


  m = 3
  allocate(H(m,n),R(m,m),xf(n),xa(n),xa2(n),yo(m))

  R = eye(m)
  H = 0
  H(1:m,1:m) = eye(m)
  yo = 1
  xf = 0

  xa = xf + matmul(matmul(P,transpose(H)), &
       matmul(inv(matmul(matmul(H,P),transpose(H)) + R), &
       (yo - matmul(H,xf))))

  ! write(6,*) 'xa',xa


  xa2 = xf + ((P.xt.H).x.(inv(((H.x.P).xt.H) + R).x.(yo - (H.x.xf))))
  call assert(xa,xa2,1e-6,'analysis')

  Hs = sparse(H)

  call assert(((H.x.P).xt.H), &
       ((Hs.x.P).xt.Hs), &
       1e-6,'covariance at observations')

  xa2 = xf + ((P.xt.Hs).x.(inv(((Hs.x.P).xt.Hs) + R).x.(yo - (Hs.x.xf))))

  call assert(xa,xa2,1e-6,'analysis using sparse matrix')

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
       1e-5,'Cx')

  !stop
  tmp = pcg(fun_Cx,yo - (Hs.x.xf))
  ! tmp = fun_Cx(yo - (Hs.x.xf))
  call assert(tmp, &
       inv(((Hs.x.P).xt.Hs) + R).x.d, &
       1e-6,'solving system using pcg')

  call assert(iC(d), &
       inv(((Hs.x.P).xt.Hs) + R).x.d, &
       1e-6,'solving system using pcg (2)')


  xa2 = xf + ((P.xt.Hs).x.tmp) 
  call assert(xa,xa2,2e-5,'analysis using sparse matrix, pcg')

  xa2 = xf + (Pc.x.(tmp.x.Hs)) 
  call assert(xa,xa2,2e-5,'analysis using sparse matrix, pcg, ConsCovar')

  xa2 = xf + K(d,n)
  call assert(xa,xa2,2e-5,'analysis using K')

  allocate(Sp(n,Nens),U(n,Nens),Sigma(Nens)) 
  allocate(PaU(Nens,Nens))
  allocate(A(n,Nens))
  allocate(KU(m,Nens))
  allocate(sqrtPaU(Nens,Nens))

  do i = 1,Nens
    Sp(:,i) = S(:,i) - K(Hs.x.S(:,i),n);
  end do

  Sigma = svd(Sp,U=U,V=sqrtPaU)

  !  call gesvd('s','n',Sp,Sigma,U,A,info)

  !  write(6,*) 'svd(Sp,U)', svd(Sp,U)


  do i = 1,Nens
    KU(:,i) = iC(Hs.x.(Pc.x.U(:,i)));
  end do

  ! A = U - H' K' U
  A = U - (Hs.tx.KU)

  ! PaU = A' * (Pc * A) + KU' * (R * KU)
  PaU = (A.tx.(Pc.x.A)) + (KU.tx.(R.x.KU))



  call test_chol()


  !write(6,*) 'mv ',maxval(abs(PaU - transpose(PaU)))

  PaU = (PaU + transpose(PaU)) / 2

  !write(6,*) 'mv ',svd(PaU)

  sqrtPaU = chol(PaU);

  allocate(Sa(n,Nens)) 
  Sa = U .x. sqrtPaU;

  call assert(matmul(transpose(sqrtPaU),sqrtPaU), &
       PaU, &
       1e-5,'Cholesky factorization (2)')

  allocate(Sa2(n,Nens)) 
  call locensanalysis(xf,S,Hs,yo,R,lpoints,hc,xa2,Sa2)

  call assert(xa,xa2,2e-5,'analysis using locensanalysis (xa)')
  call assert(Sa,Sa2,4e-5,'analysis using locensanalysis (Sa)')


  deallocate(x)
  deallocate(v)
  deallocate(v2)
  deallocate(S)
  deallocate(w,jj)

  deallocate(P)
  deallocate(Pr)
  deallocate(Pc,hc)
  deallocate(H,R,xf,xa,xa2,yo)
  deallocate(tmp,d)


  deallocate(Sp)
  deallocate(U,Sigma) 
  deallocate(PaU)
  deallocate(A)
  deallocate(KU)
  deallocate(sqrtPaU)

 contains



  function fun_Cx(x) result(y)
   real, intent(in) :: x(:)
   real :: y(size(x))
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

   !  x = ((P.xt.Hs).x.iC(d))  
   x = Pc.x.(iC(d).x.Hs)
  end function K


  function locenscovx(Pc,H,R,y) result(Cy)
   use matoper
   class(ConsCovar), intent(in) :: Pc
   type(SparseMatrix) :: H
   real :: R(:,:)
   real :: y(:)
   real :: Cy(size(y,1))
   real :: Hty(Pc%n)


   Hty = y.x.H
   Cy = (H .x. (Pc.x.Hty))  + (R.x.y)
  end function locenscovx


  subroutine lpoints(i,nnz,j,w)
   integer, intent(in) :: i
   integer, intent(out) :: nnz,j(:)
   real, intent(out) :: w(:)

   call locpoints(i,x,len,nnz,j,w)  
  end subroutine lpoints

 end subroutine run_test

end module test_suite



program test
 use test_suite

 call run_test

end program test


