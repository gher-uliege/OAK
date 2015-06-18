program toymodel

implicit none
integer, parameter :: n = 1000
integer, parameter :: Ntime = 100

real :: x(n)
integer :: i

! initialize
x = 0
x(1) = 1

! time loop
do i = 1,Ntime
  x(2:n) = x(1:n-1)
  x(1) = x(n)


end do

end program toymodel
