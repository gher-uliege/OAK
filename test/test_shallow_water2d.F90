! gfortran -Wall  -g -O2 -I/usr/include  -fdefault-real-8 -o test_shallow_water2d  shallow_water2d.F90 test_shallow_water2d.F90 -L/usr/lib -lnetcdff -lnetcdf  && time ./test_shallow_water2d

program test_shallow_water2d
 use shallow_water2d
 implicit none
 integer, parameter :: m = 100, n = 100, Ntime = 2
 real :: h(m,n), zeta(m,n,Ntime), dx, dy
 real :: U(m-1,n,Ntime), V(m,n-1,Ntime)
 logical :: mask(m,n) 
 type(domain) :: dom
 real :: dt,g,f
 integer :: timecounter, timeindex 
 character(len=*), parameter :: fname = 'example.nc'
 
 mask = .false.
 mask(2:m-1,2:m-1) = .true.
 mask(10:40,50:60) = .false.
 mask(70,20) = .false.

 dx = 1000
 dy = 1000

 h = 100
 call init_domain(dom,dx,dy,mask,h)

 zeta = 0
 zeta(4,4,1) = 1
 U = 0
 V = 0
 dt = 2
 g = 9.81
 f = 1e-4
 timeindex = 1

 do timecounter = 1,1000
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
