

program test_matoper
 use matoper
 implicit none


 !call benchmark_matoper
 call test_pcg
 call test_chol
 call test_sqrtm
 call test_sort
 call test_mergesort
 call test_unique

 contains

 subroutine test_chol()
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

  !_______________________________________________________
  !

  subroutine mergesort (a)
   implicit none
   integer :: a(:)
   integer :: b(size(a))

   integer :: rght, rend, num
   integer :: i,j,m, k, left

   num = size(a)
   k = 1
   do while (k < num) 
     do left=0,num-k-1,2*k
       rght = left + k        
       rend = rght + k
       if (rend > num) rend = num

       m = left+1 
       i = left+1
       j = rght+1

       ! merge
       do while (i-1 < rght .and. j-1 < rend)
         if (a(i) <= a(j)) then         
           b(m) = a(i) 
           i = i+1
         else
           b(m) = a(j) 
           j = j+1
         end if

         m = m+1
       end do

       do while (i-1 < rght) 
         b(m)=a(i) 
         i = i+1
         m = m+1
       end do

       do while (j-1 < rend)
         b(m)=a(j) 
         j = j+1
         m = m+1
       end do

       ! copy over
       do m=left+1,rend
         a(m) = b(m)
       end do
     end do

     k = k*2
   end do
  end subroutine mergesort

  !_______________________________________________________
  !

  subroutine test_mergesort
   implicit none
   integer, parameter :: n = 11
   integer :: ind(n)
   integer, dimension(1:n) :: A = &  
        (/0, 50, 20, 25, 90, 10, 5, 20, 99, 75, 2/), sortedA
   
   sortedA = A
   call mergesort(sortedA)
   call assert(all(sortedA(1:n-1) <= sortedA(2:n)),'merge sort (1)')
!   call assert(all(sortedA == A(ind)),'quick sort (2)')
  end subroutine test_mergesort

  !_______________________________________________________
  !

  subroutine test_unique
   implicit none
   integer, parameter :: n = 10
   integer :: ind(n), ind2(n), i, nu
   integer, dimension(1:n) :: A = &  
        (/0, 20, 20, 25, 90, 10, 5, 20, 90, 75/), uA
   
   uA = A
   call unique(uA,nu,ind,ind2)

   ! all elements in A must be one time in uA
   do i = 1,n
     if (count(A(i) == uA(1:nu)) /= 1) then
       write(6,*) 'unique (1): FAILED'
       stop
     end if
   end do

   write(6,*) 'unique (1): OK'

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


end program test_matoper
