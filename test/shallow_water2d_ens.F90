! gfortran -Wall  -g -O2 -I/usr/include  -fdefault-real-8 -o shallow_water2d_ens  shallow_water2d.F90 shallow_water2d_ens.F90 -L/usr/lib -lnetcdff -lnetcdf  && time ./shallow_water2d_ens

program test_shallow_water2d
 use shallow_water2d
 use matoper, only: randn
 implicit none
 ! size of time and number of time 
 integer, parameter :: m = 100, n = 100
 ! size of ensemble
 integer, parameter :: Nens = 1
 ! number time steps in memory
 integer, parameter :: Ntime = 2
 ! grid spacing in meters
 real :: dx, dy
 ! position of initial perturbation
 real :: xc, yc
 real :: zetac
 ! size of initial perturbation
 real :: lenx, leny, x, y, h(m,n)
 ! surface elevation
 real :: zeta(m,n,Ntime,Nens)
 ! depth-integrated currents
 real :: U(m-1,n,Ntime,Nens), V(m,n-1,Ntime,Nens)
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

 integer :: memberindex
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

 ! time index in NetCDF file
 timeindex = 1


 ! initial condition

 ! initially the velocity is zero
 U = 0
 V = 0

 ! ensemble loop
 do  memberindex = 1,Nens

   ! position of initial perturbation and extension (in meters)

   xc = 10e3 + 10e3 * randn()
   yc = 50e3 + 10e3 * randn()
   zetac = randn()
   lenx = 20e3
   leny = 20e3

   write(6,*) 'xc,yc,zetac',xc,yc,zetac

   do j = 1,n
     do i = 1,m
       x = dx*(i-1)
       y = dy*(j-1)
       zeta(i,j,1,memberindex) = zetac * exp(-((x-xc)/lenx)**2 - ((y-yc)/leny)**2)
     end do
   end do
   where (.not.mask) zeta(:,:,1,memberindex) = 0

   ! save initial conditions
   call shallow_water2d_save(dom,timeindex,zeta(:,:,1,memberindex), &
        U(:,:,1,memberindex),V(:,:,1,memberindex),fname,memberindex = memberindex)

 end do



 timeindex = timeindex + 1

 ! time loop
 time: do timecounter = 1,90
   ! ensemble loop
   ensemble: do memberindex = 1,Nens

     call shallow_water2d_step(dom,timecounter, & 
          zeta(:,:,1,memberindex),U(:,:,1,memberindex),V(:,:,1,memberindex), &
          dt,g,f,zeta(:,:,2,memberindex),U(:,:,2,memberindex),V(:,:,2,memberindex))

     ! swap time instances
     zeta(:,:,1,memberindex) = zeta(:,:,2,memberindex)
     U(:,:,1,memberindex) = U(:,:,2,memberindex)
     V(:,:,1,memberindex) = V(:,:,2,memberindex)

     !     if (mod(timecounter,10) == 0) then
     !       call diag(dom,timecounter,zeta(:,:,1,memberindex), &
     !            U(:,:,1,memberindex),V(:,:,1,memberindex),g)
     !     end if

     if (mod(timecounter,30) == 0) then
       call shallow_water2d_save(dom,timeindex,zeta(:,:,1,memberindex), &
            U(:,:,1,memberindex),V(:,:,1,memberindex),fname,memberindex = memberindex)

       if (memberindex == Nens) timeindex = timeindex + 1
     end if
   end do ensemble
 end do time
end program test_shallow_water2d
