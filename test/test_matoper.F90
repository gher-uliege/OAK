
program test_matoper
 use matoper
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

  do i=1,100

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

  write(6,*) 'check sum ',sum(ampl)
  write(6,*) finish-start




end program test_matoper
