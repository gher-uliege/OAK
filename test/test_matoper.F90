! Copyright(c) 2002-2009 Alexander Barth and Luc Vandenblucke

program test_rrsqrt
 use matoper
 implicit none
 integer, parameter :: r = 4, jmax = 20, m=10
 real :: lambda(r), UT(r,r), yo(m), Hxf(m), invsqrtR(m), HSf(m,r)

 real :: dummy(1,1)
 integer :: info,i,j

 real :: temp(r,r), ampl(r), ampl2(r)
 real :: sd(m), tmp2(r), tmp3(r), tmp4(r)

#ifdef WITH_OPERATORS

#else
 
#endif 

 do j=1,jmax
   do i=1,r
     temp(i,j) = i+j
   end do
 end do

  call gesvd('n','a',temp,lambda,dummy,UT,info)

!#ifdef WITH_OPERATORS
  ampl = (UT.tx.(lambda.dx.(UT.x.(HSf.tx.(invsqrtR**2*(yo-Hxf))))))  
  write(6,*) 'ampl ',ampl
!#else

  sd = invsqrtR**2*(yo-Hxf)
  call DGEMV('T',m,r,1.,HSf,m,sd,1,0.,tmp2,1)
  call DGEMV('N',r,r,1.,UT,r,tmp2,1,0.,tmp3,1)
  tmp3 = lambda * tmp3
  call DGEMV('T',r,r,1.,UT,r,tmp3,1,0.,ampl2,1)

  write(6,*) 'ampl2 ',ampl2
!#endif

end program test_rrsqrt
