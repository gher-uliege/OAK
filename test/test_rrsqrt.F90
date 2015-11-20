! Copyright(c) 2002-2009 Alexander Barth and Luc Vandenblucke

program test_rrsqrt
 use matoper
 implicit none
 integer, parameter :: imax = 4, jmax = 20
 real :: lambda(imax), UT(imax,imax)
 real :: dummy(1,1)
 integer :: info,r,i,j

 real :: temp(imax,jmax)


 do j=1,jmax
   do i=1,imax
     temp(i,j) = i+j
   end do
 end do

 call gesvd('n','a',temp,lambda,dummy,UT,info)

 write(*,*) temp.tx.temp
! write(*,*) UT.x.(diag(lambda**2).xt.UT)
 write(*,*) UT.tx.(diag(lambda**2).x.UT)

end program test_rrsqrt
