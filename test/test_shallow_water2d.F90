! gfortran -Wall  -g -O2 -I/usr/include  -fdefault-real-8 -o test_shallow_water2d  shallow_water2d.F90 test_shallow_water2d.F90 -L/usr/lib -lnetcdff -lnetcdf  && time ./test_shallow_water2d

program test_shallow_water2d
 use shallow_water2d
 implicit none
 integer, parameter :: m = 100, n = 100, Ntime = 2
 real :: x, y, h(m,n), zeta(m,n,Ntime), dx, dy
 real :: U(m-1,n,Ntime), V(m,n-1,Ntime)
 logical :: mask(m,n) 
 type(domain) :: dom
 real :: dt,g,f, ew
 integer :: timecounter, timeindex, i,j
 character(len=*), parameter :: fname = 'example.nc'
 
 ! setup mask
 mask = .false.
 mask(2:m-1,2:m-1) = .true.
 mask(10:40,50:60) = .false.
 mask(70,20) = .false.

 ! grid spacing
 dx = 1000
 dy = 1000

 do j = 1,n
   do i = 1,m
      x = dx*(i-1)
      y = dy*(j-1)

      ew = (erf((x - 70e3)/10e3)+1)/2
      ! depth of domain
      !h(i,j) = 100
      h(i,j) = 5000* (exp(-(y - 49500)**2 / 20000**2)*ew  - 2*ew+1)
      
      zeta(i,j,1) = exp(-(x/(20*dx))**2 - (y/(20*dy))**2)
      !zeta(i,j,1) = exp(-x/(20*dx))
      
    end do
  end do

  
  mask = h > 0
  ! close east and west
  mask([1,m],:) = .false. 
  ! close north and south
  mask(:,[1,n]) = .false.
  where (h < 2) h = 2

  where (.not.mask) zeta(:,:,1) = 0

  call init_domain(dom,dx,dy,mask,h)

! zeta = 0
! zeta(4,4,1) = 1
 U = 0
 V = 0
 dt = 2
 g = 9.81
 f = 1e-4
 timeindex = 1

 ! save initial conditions
 call shallow_water2d_save(dom,timeindex,zeta(:,:,1),U(:,:,1),V(:,:,1),fname)
 timeindex = timeindex + 1

 ! time loop
 do timecounter = 1,3000
   call shallow_water2d_step(dom,timecounter,zeta(:,:,1),U(:,:,1),V(:,:,1), &
        dt,g,f,zeta(:,:,2),U(:,:,2),V(:,:,2))

   zeta(:,:,1) = zeta(:,:,2)
   U(:,:,1) = U(:,:,2)
   V(:,:,1) = V(:,:,2)

   if (mod(timecounter,10) == 0) then
     call diag(dom,timecounter,zeta(:,:,1),U(:,:,1),V(:,:,1),g)
   end if

   if (mod(timecounter,100) == 0) then
     call shallow_water2d_save(dom,timeindex,zeta(:,:,1),U(:,:,1),V(:,:,1),fname)
     timeindex = timeindex + 1
   end if
 end do
end program test_shallow_water2d
