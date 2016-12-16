! example program how to use ndgrid
! for a two dimensional grid

! cd ~/projects/Fortran/OAK &&  make test/test_ndgrid && test/test_ndgrid

program test_ndgrid
 use matoper
 use ndgrid
 implicit none

 call test_regular_ndgrid_2d(10,20)
 ! 2D
 call test_ndgrid_2d(10,20)
 call test_ndgrid_2d(1,3)
 call test_ndgrid_2d(3,1)

 ! nD
 call test_ndgrid_nd((/10/))
 call test_ndgrid_nd((/10,20/))
 call test_ndgrid_nd((/2,2,2/))
 call test_ndgrid_nd((/2,2,2,2/))

 ! degenerated cases

 call test_ndgrid_nd((/1,2/))
 call test_ndgrid_nd((/2,1,2/))
 call test_ndgrid_nd((/10,1,20/))
 call test_ndgrid_nd((/10,3,1,20/))
 call test_ndgrid_nd((/10,3,1/))
 call test_ndgrid_nd((/1,10,3,1/))
 call test_ndgrid_nd((/1,1,1,10/))
 call test_ndgrid_nd((/1,1,1,10,1/))


contains 

 function fun(x,y)
  implicit none
  real :: x,y,fun
  fun = 2*x+4*y
 end function fun

 subroutine test_ndgrid_2d(m,n)

  ! size of the domain
  integer, intent(in) :: m,n

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
  real :: xi, yi, fi, coeff(4), ref
  logical :: out

  ! initialize coordinate, field and mask
  do j=1,n
    do i=1,m
      x(i,j) = i+1
      y(i,j) = j+2
      f(i,j) = fun(x(i,j),y(i,j))
      mask(i,j) = .false.
    end do
  end do

  ! (xi,yi) is the middle of the domain

  xi = sum(x)/(n*m)
  yi = sum(y)/(n*m)

  ! initialize grid
  call initgrid(g,x,y,mask)

  ! get interpolated value
  call interp(g,f,(/xi,yi/),fi,out)

  ref = fun(xi,yi)

  ! write(6,*) 'fi ',fi
  ! write(6,*) 'ref ',fun(xi,yi) ! should be the same as fi
  ! write(6,*) 'out ',out        ! should be .false.

  if (abs(ref-fi) > 1e-6 .or. out) then
    write(6,*) 'domain ',m,n,'FAILED'
    stop
  end if


  ! get interpolation coefficients  
  call cinterp(g,(/xi,yi/),indexes,coeff,nbp)
  fi = 0
  do i=1,nbp
    fi = fi + coeff(i) * f(indexes(1,i),indexes(2,i))
  end do

  call assert(.not.out,'interp domain '//trim(fmtsize([m,n])) // ' out')
  call assert(ref,fi,1e-6,'interp domain '//trim(fmtsize([m,n])) // ' interpolation')

 end subroutine test_ndgrid_2d


!-----------------------
! regular grid

 subroutine test_regular_ndgrid_2d(m,n)

  ! size of the domain
  integer, intent(in) :: m,n

  ! x and y coordinates 
  ! f field to interpolate
  real, dimension(m,n) :: x,y,f

  ! land-sea mask (land or invalid points equal to .true.)
  logical, dimension(m,n) :: mask

  ! for bi-linear interpolation we have maximal 2**n interpolation coefficients
  ! for n=2, 2**n is 4
  integer :: i,j, indexes(2,4), nbp
  type(regulargrid) :: g

  ! xi,yi location to interpolate
  real :: xi, yi, fi, coeff(4), ref
  logical :: out

  ! initialize coordinate, field and mask
  do j=1,n
    do i=1,m
      x(i,j) = i+1
      y(i,j) = j+2
      f(i,j) = fun(x(i,j),y(i,j))
      mask(i,j) = .false.
    end do
  end do

  ! (xi,yi) is the middle of the domain

  xi = sum(x)/(n*m)
  yi = sum(y)/(n*m)

  ! initialize grid
  call init_regulargrid(g,[m,n],[1.,2.],[1.,1.],reshape(mask,[m*n]))

  ! get interpolated value
  call interp(g,f,(/xi,yi/),fi,out)

  ref = fun(xi,yi)

  ! write(6,*) 'fi ',fi
  ! write(6,*) 'ref ',fun(xi,yi) ! should be the same as fi
  ! write(6,*) 'out ',out        ! should be .false.

  if (abs(ref-fi) > 1e-6 .or. out) then
    write(6,*) 'domain ',m,n,'FAILED'
    stop
  end if


  ! get interpolation coefficients  
  call cinterp(g,(/xi,yi/),indexes,coeff,nbp)
  fi = 0
  do i=1,nbp
    fi = fi + coeff(i) * f(indexes(1,i),indexes(2,i))
  end do

  call assert(.not.out,'interp domain '//trim(fmtsize([m,n])) // ' out')
  call assert(ref,fi,1e-6,'interp domain '//trim(fmtsize([m,n])) // ' interpolation')

 end subroutine 




 function fun_nd(x)
  implicit none
  real :: x(:),fun_nd
  integer :: i

  fun_nd = sum((/(2 * i * x(i),i=1,size(x))/))
 end function fun_nd


 ! output a string from vector of integer size
 ! e.g. is size = [12,13,10] the output will be 
 ! '12x13x10'

 function fmtsize(sz) result(s)
  implicit none
  integer, intent(in) :: sz(:)
  character(len=80) :: s, tmp
  integer :: i

  do i = 1,size(sz)
    write(tmp,*) sz(i)
    if (i == 1) then
      s = adjustl(tmp)
    else
      s = trim(s)//'x'//adjustl(tmp)
    end if
  end do

 end function fmtsize

 subroutine test_ndgrid_nd(sz)
  implicit none

  ! size of the domain
  integer, intent(in) :: sz(:)

  ! x(:,1),x(:,2),.. coordinates 
  ! f field to interpolate
  real, dimension(product(sz),size(sz)) :: x
  real, dimension(product(sz)) :: f

  ! number of dimensions
  integer :: n
  ! land-sea mask (land or invalid points equal to .true.)
  logical, dimension(product(sz)) :: mask

  ! for bi-linear interpolation we have maximal 2**n interpolation coefficients
  ! for n=2, 2**n is 4
  integer :: i,j, indexes(size(sz),2**size(sz)), nbp, ioffset(size(sz)), linindex
  type(grid) :: g

  ! xi,yi location to interpolate
  real :: xi(size(sz)), fi, coeff(2**size(sz)), ref, l
  logical :: out

  ! number of dimensions
  n = size(sz)

  write(6,*) 'Testing domain ',trim(fmtsize(sz))

  ! initialize coordinate, field and mask

  ioffset(1) = 1
  do i=2,n
    ioffset(i) = ioffset(i-1)*sz(i-1)
  end do

  mask = .false.
  do j=1,product(sz)
    linindex = j-1

    do i=n,1,-1
      x(j,i) = linindex/ioffset(i)
      linindex = linindex - x(j,i)*ioffset(i)
    end do

    if (sum(x(j,:) * ioffset) .ne. j-1) then
      write(0,*) 'fail', j
      write(0,*) 'ioffset', ioffset
      write(0,*) 'x(j,:)', x(j,:)
      write(0,*) 'sz', sz

      stop
    end if

    f(j) = fun_nd(x(j,:))
  end do

  ! (xi,yi) is the middle of the domain

  xi = sum(x,1)/product(sz)

  ! initialize grid
  call initgrid(g,n,sz,mask)
  do i=1,n
    call setCoord(g,i,x(:,i))
  end do

  ! get interpolated value
  call interp(g,f,xi,fi,out)

  ref = fun_nd(xi)

  call assert(.not.out,'interp domain '//trim(fmtsize(sz)) // ' out')
  call assert(ref,fi,1e-6,'interp domain '//trim(fmtsize(sz)) // ' interpolation')

  ! get interpolation coefficients  
  call cinterp(g,xi,indexes,coeff,nbp)
  fi = 0
  do i=1,nbp
    fi = fi + coeff(i) * f(sum((indexes(:,i)-1) * ioffset)+1)
  end do

  call assert(ref,fi,1e-6,'cinterp domain '//trim(fmtsize(sz)) // ' interpolation')


 end subroutine test_ndgrid_nd

end program test_ndgrid

! LocalWords:  ndgrid yi fi
