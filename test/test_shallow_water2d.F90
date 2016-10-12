! gfortran -Wall  -g -O2 -I/usr/include  -fdefault-real-8 -o test_shallow_water2d  shallow_water2d.F90 test_shallow_water2d.F90 -L/usr/lib -lnetcdff -lnetcdf  && time ./test_shallow_water2d

program test_shallow_water2d
 use shallow_water2d
 implicit none
 ! size of time and number of time 
 integer, parameter :: m = 100, n = 100
 ! number time steps in memory
 integer, parameter :: Ntime = 2
 ! grid spacing in meters
 real :: dx, dy
 ! position of initial perturbation
 real :: xc, yc
 ! size of initial perturbation
 real :: lenx, leny, x, y, h(m,n)
 ! surface elevation
 real :: zeta(m,n,Ntime)
 ! depth-integrated currents
 real :: U(m-1,n,Ntime), V(m,n-1,Ntime)
 ! land-sea mask (see is .true., land is .false.)
 logical :: mask(m,n) 
 ! all information of the model domain
 type(domain) :: dom
 ! time step
 real :: dt
 ! acceleration due to gravity
 real :: g
 ! Coriolis parameter
 real :: f
 real :: ew, tmp
 ! number of time steps to simulate
 integer :: timecounter
 ! time index in NetCDF file
 integer :: timeindex
 ! spatial indices
 integer :: i,j
 ! output file name
 character(len=*), parameter :: fname = 'example.nc'

 ! grid spacing (in meters)
 dx = 1000
 dy = 1000

 do j = 1,n
   do i = 1,m
     x = dx*(i-1)
     y = dy*(j-1)
     ew = (erf((x - 50e3)/10e3)+1)/2; 
     tmp = (erf((x - 80e3)/10e3)+1)/2;
     ! depth of domain
     !h(i,j) = 100
     !h(i,j) = 5000* (exp(-(y - 49500)**2 / 20000**2)*ew  - 2*ew+1)

     h(i,j) =  5000*(1-ew) + 100*(1-tmp) - 20 + 20*exp(-(y - 49500)**2 / 20000**2)*tmp; 
   end do
 end do

 ! setup mask
 mask = h > 0
 ! close eastern and western boundary
 mask([1,m],:) = .false. 
 ! close northern and southern boundary
 mask(:,[1,n]) = .false.
 ! set minimum depth
 where (h < 2) h = 2


 call init_domain(dom,dx,dy,mask,h)
 
 ! time step and other model parameters
 dt = 2
 dt = 1
 g = 9.81
 f = 1e-4

 ! initial condition
 ! position of initial perturbation and extension (in meters)

 xc = 0
 yc = 0
 lenx = 20e3
 leny = 20e3

 do j = 1,n
   do i = 1,m
     x = dx*(i-1)
     y = dy*(j-1)
     zeta(i,j,1) = exp(-((x-xc)/lenx)**2 - ((y-yc)/leny)**2)
   end do
 end do
 where (.not.mask) zeta(:,:,1) = 0

 ! initially the velocity is zero
 U = 0
 V = 0

 ! time index in NetCDF file
 timeindex = 1

 ! save initial conditions
 call shallow_water2d_save(dom,timeindex,zeta(:,:,1),U(:,:,1),V(:,:,1),fname)
 timeindex = timeindex + 1

 ! time loop
 do timecounter = 1,9000
   call shallow_water2d_step(dom,timecounter,zeta(:,:,1),U(:,:,1),V(:,:,1), &
        dt,g,f,zeta(:,:,2),U(:,:,2),V(:,:,2))

   ! swap time instances
   zeta(:,:,1) = zeta(:,:,2)
   U(:,:,1) = U(:,:,2)
   V(:,:,1) = V(:,:,2)

   if (mod(timecounter,10) == 0) then
     call diag(dom,timecounter,zeta(:,:,1),U(:,:,1),V(:,:,1),g)
   end if

   if (mod(timecounter,30) == 0) then
     call shallow_water2d_save(dom,timeindex,zeta(:,:,1),U(:,:,1),V(:,:,1),fname)
     timeindex = timeindex + 1
   end if
 end do
end program test_shallow_water2d
