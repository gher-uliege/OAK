! example program how to use ndgrid
! for a two dimensional grid

program test_ndgrid
 use ndgrid
 implicit none

 ! size of the domain
 integer, parameter :: m=10,n=20

 ! x and y coordinates 
 ! f field to interpolate
 real, dimension(m,n) :: x,y,f

 ! land-sea mask (land or invalid points equal to .true.)
 logical, dimension(m,n) :: mask

 ! for bi-linear interpolation we have maximal 2**n interpolation coefficients
 ! for n=2, 2**n is 4
 integer :: i,j, indexes(2,4), nbp
 type(grid) :: g

 ! xi,yi location to interpolate
 real :: xi = 3.1, yi = 6.3, fi, coeff(4)
 logical :: out
 
! initialize coordinate, field and mask
 do j=1,n
   do i=1,m
     x(i,j) = i
     y(i,j) = j
     f(i,j) = fun(x(i,j),y(i,j))
     mask(i,j) = .false.
   end do
 end do

 ! initialize grid
 call initgrid(g,x,y,mask)

! get interpolated value
 call interp(g,f,(/xi,yi/),fi,out)

 write(6,*) 'fi ',fi
 write(6,*) 'ref ',fun(xi,yi) ! should be the same as fi
 write(6,*) 'out ',out        ! should be .false.

! get interpolation coefficients  
 call cinterp(g,(/xi,yi/),indexes,coeff,nbp)
 fi = 0
 do i=1,nbp
   fi = fi + coeff(i) * f(indexes(1,i),indexes(2,i))
 end do

 write(6,*) 'fi ',fi ! the same as before
contains 

function fun(x,y)
  implicit none
  real :: x,y,fun
  fun = 2*x+4*y
end function fun

end program test_ndgrid

! LocalWords:  ndgrid yi fi
