! cd /home/abarth/Assim/OAK-nonDiagR &&  make test/test_matoper && test/test_matoper


program test_matoper
 use matoper
 implicit none


 !call benchmark_matoper
 call test_factorial
 call test_pcg
 call test_chol
 call test_sqrtm
 call test_quicksort
 call test_mergesort
 call test_unique
 call test_sparse
 call test_randperm


 contains

 subroutine test_factorial
  implicit none

  call assert(factorial(0),1,'test factorial (1)')
  call assert(factorial(6),720,'test factorial (2)')
  call assert(nchoosek(5,3),10,'test nchoosek')
  
  
 end subroutine test_factorial

 subroutine test_chol()
  implicit none
  real :: A(3,3), U(3,3)

  A = reshape([2.,1.,1., 1.,2.,1., 1.,1.,2.],[3,3])
  U = chol(A)
  call assert(matmul(transpose(U),U), &
      A, &
      1e-6,'Cholesky factorization')
 end subroutine test_chol

  !_______________________________________________________
  !

 subroutine test_sqrtm()
  implicit none
  real :: A(3,3)
  real :: S(size(A,1),size(A,1))

  A = reshape([2.,1.,1., 1.,2.,1., 1.,1.,2.],[3,3])
  S = sqrtm(A)

  call assert(matmul(S,S), &
      A, &
      1e-6,'sqrtm factorization')
 end subroutine test_sqrtm


  !_______________________________________________________
  !

  subroutine test_quicksort
   implicit none
   integer, parameter :: n = 10
   integer :: ind(n)
   integer, dimension(1:n) :: A = &  
        (/0, 50, 20, 25, 90, 10, 5, 20, 99, 75/), sortedA
   
   sortedA = A
   call quicksort(sortedA,ind)
   call assert(all(sortedA(1:n-1) <= sortedA(2:n)),'quick sort (1)')
   call assert(all(sortedA == A(ind)),'quick sort (2)')
  end subroutine test_quicksort

  !_______________________________________________________
  !

  subroutine test_mergesort
   implicit none
   integer, parameter :: n = 11
   integer :: ind(n)
   integer, dimension(1:n) :: A = &  
        (/0, 50, 20, 25, 90, 10, 5, 20, 99, 75, 2/), sortedA
   
   sortedA = A
   call sort(sortedA,ind)
   call assert(all(sortedA(1:n-1) <= sortedA(2:n)),'merge sort (1)')
   call assert(all(sortedA == A(ind)),'merge sort (2)')
  end subroutine test_mergesort

  !_______________________________________________________
  !

  subroutine test_unique
   implicit none
   integer, parameter :: n = 10
   integer :: ind(n), ind2(n), i, nu
   integer, dimension(1:n) :: A = &  
        (/0, 20, 20, 25, 90, 10, 5, 20, 90, 75/), uA
   logical :: success
   
   uA = A
   call unique(uA,nu,ind,ind2)

   success = .true.
   ! all elements in A must be one time in uA
   do i = 1,n
     if (count(A(i) == uA(1:nu)) /= 1) then
       success = .false.
       exit
     end if
   end do

   call assert(success,'unique (1)')
   call assert(all(uA(1:nu) == A(ind(1:nu))),'unique (2)')
   call assert(all(A == uA(ind2)),'unique (3)')

  end subroutine test_unique

  !_______________________________________________________
  !
  
  function test_fun(x) result(y)
   implicit none
   real, intent(in) :: x(:)
   real :: y(size(x))
   y = x/2.
  end function test_fun
  
  subroutine test_pcg()
   implicit none
   real :: x(3)
   integer :: nit
   
   x = pcg(test_fun,(/ 1.,2.,3. /),nit=nit)
   call assert(x,(/ 2.,4.,6. /),1e-6,'conjugate gradient')
  end subroutine test_pcg

  !_______________________________________________________
  !

  subroutine benchmark_matoper
   implicit none

 integer, parameter :: r = 100, m=100000
 real(8) :: lambda(r), UT(r,r), yo(m), Hxf(m), invsqrtR(m), HSf(m,r)

 real(8) :: dummy(1,1)
 integer :: info,i,j

 real(8) :: temp(m,r), ampl(r)
 real(8) :: start, finish

#ifndef WITH_OPERATORS
 real(8) :: tmp(m), tmp2(r), tmp3(r)
#endif 

 do i=1,r
   yo(i) = mod(real(i,8),10._8)
   Hxf(i) = mod(real(i,8),5._8)
   invsqrtR(i) = 6._8
 end do

 do j=1,r
   do i=1,m
     HSf(i,j) =  mod(real(i+j,8),10._8)/10._8 - 0.5_8
     temp(i,j) = invsqrtR(i)*HSf(i,j)     
   end do
 end do


  call gesvd('n','a',temp,lambda,dummy,UT,info)
  lambda = (1+lambda**2)**(-1)

  call cpu_time(start)

  do i=1,1000

#ifdef WITH_OPERATORS
    ampl = UT.tx.(lambda.dx.(UT.x.(HSf.tx.(invsqrtR**2*(yo-Hxf)))))
#else
    tmp = invsqrtR**2*(yo-Hxf)
    call DGEMV('T',m,r,1._8,HSf,m,tmp,1,0._8,tmp2,1)
    call DGEMV('N',r,r,1._8,UT,r,tmp2,1,0._8,tmp3,1)
    tmp3 = lambda * tmp3
    call DGEMV('T',r,r,1._8,UT,r,tmp3,1,0._8,ampl,1)
#endif

    ! make loop interation depedent on each other
    yo(1) = yo(1)+ampl(1)/1000._8

  end do

  call cpu_time(finish)

!  write(6,*) 'check sum ',sum(ampl)
  write(6,*) finish-start

 end subroutine benchmark_matoper

 !_______________________________________________________
 !

 subroutine test_sparse
  implicit none
  type(SparseMatrix) :: A,B,C
  integer, parameter :: n = 3
  real, parameter :: tol = 1e-6

  A = speye(3)

  A = sparse([1,2,3],[1,2,3],[1.,1.,1.],n,n)
  B = A + A
  call assert(full(B),full(A) + full(A),tol,'add sparse matrices (1)')

  B = sparse_compress(A + A)
  call assert(full(B),full(A) + full(A),tol,'add sparse matrices (2)')


  A = sparse([1,2,3,2],[1,2,3,1],[1.,2.,1.,3.],n,n)
  B = sparse([1,3,2,2],[1,1,2,1],[1.,2.,1.,3.],n,n)
  call assert(full(sparse_compress(A+B)),full(A) + full(B),tol,'add sparse matrices (3)')

  call assert(full(sparse_compress(A.x.B)),full(A).x.full(B),tol,'mult sparse matrices')

  A = sparse([1,2,3,2],[1,2,3,2],[1.,2.,1.,3.],n,n)
  A = sparse_compress(A)
  call assert(A%nz,3,'sparse compress')
   
!  call spprint(A)


 end subroutine test_sparse

 !_______________________________________________________
 !

 subroutine test_randperm
  implicit none
  
  integer, parameter :: n = 10
  integer :: ind(n), nu
  ind = randperm(n)
  call unique(ind,nu)

  ! check if we have no repetion
  call assert(nu,n,'randperm')
 end subroutine test_randperm
 

end program test_matoper
