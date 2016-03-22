! cd /home/abarth/Assim/OAK_classic/test; gfortran -Wall  -g -O2 -I/usr/include -I$HOME/Assim/OAK_classic  -fdefault-real-8 -o shallow_water2d_ens  shallow_water2d.F90 shallow_water2d_ens.F90 -L/usr/lib -lnetcdff -lnetcdf -L.. -loak -llapack -lblas  && time ./shallow_water2d_ens

! cd /home/abarth/Assim/OAK_classic/; make; cd /home/abarth/Assim/OAK_classic/test; gfortran -Wall  -g -O2 -I/usr/include -I$HOME/Assim/OAK_classic  -fdefault-real-8 -o shallow_water2d_ens  shallow_water2d.F90 shallow_water2d_ens.F90 -L/usr/lib -lnetcdff -lnetcdf -L.. -loak -llapack -lblas  && time ./shallow_water2d_ens


!#define EXTRACT_OBS
#define OAK
!#define ENS

#ifdef OAK
#define ENS
#endif

program test_shallow_water2d
 use shallow_water2d
 use matoper, only: randn
#ifdef EXTRACT_OBS
 use ndgrid, only: initgrid, grid, interp
 use netcdf
 use initfile
#endif
#ifdef OAK
 use oak
 use assimilation
#endif

 implicit none
 ! size of time
 integer, parameter :: m = 100, n = 100
 ! size of time steps
 integer, parameter :: Nsteps = 3000
! integer, parameter :: Nsteps = 90
#ifdef ENS
 ! size of ensemble
 integer, parameter :: Nens = 100
#else
 integer, parameter :: Nens = 1
#endif
 ! grid spacing in meters
 real :: dx, dy
 ! position of initial perturbation
 real :: xc, yc
 real :: zetac
 ! size of initial perturbation
 real :: lenx, leny, x, y, h(m,n)
 ! surface elevation
 real :: zeta(m,n,Nens)
 ! depth-integrated currents
 real :: U(m-1,n,Nens), V(m,n-1,Nens)
 ! surface elevation (next time step)
 real :: zeta_next(m,n)
 ! depth-integrated currents (next time step)
 real :: U_next(m-1,n), V_next(m,n-1)
 ! land-sea mask (see is .true., land is .false.)
 logical :: mask(m,n) 
 ! all information of the model domain
 type(domain) :: dom
 ! time step
 real :: dt
 ! time (seconds since start)
 real(8) :: time2
 ! acceleration due to gravity
 real :: g
 ! Coriolis parameter
 real :: f
 real :: ew, tmp, tmp2
 ! number of time steps to simulate
 integer :: timecounter
 ! time index in NetCDF file
 integer :: timeindex

 integer :: memberindex
 ! spatial indices
 integer :: i,j
 ! output file name
#ifdef EXTRACT_OBS
 character(len=*), parameter :: fname = 'truth.nc'
#elif defined(OAK)
 character(len=*), parameter :: fname = 'assim_ensemble.nc'
#elif defined(ENS)
 character(len=*), parameter :: fname = 'free_ensemble.nc'
#else
 character(len=*), parameter :: fname = 'free.nc'
#endif


#ifdef EXTRACT_OBS
 character(len=*), parameter :: obsloc = 'obsloc.init'
 integer :: nobs
 integer :: l, obsindex
 type(grid) :: modgrid 
 real :: timeobs(Nsteps)
 real, pointer :: xobs(:),yobs(:),zetaobs(:,:)
 logical :: out
 real(8) :: fillvalue = NF90_FILL_DOUBLE
#endif

#ifdef OAK
  type(oakconfig) :: config
!  integer, allocatable :: subdomain(:)
#endif

  ! set seed for random numbers generator
  call init_random_seed()


#ifdef EXTRACT_OBS
 obsindex = 0
 call getInitValue(obsloc,'Observations.gridX',xobs)
 call getInitValue(obsloc,'Observations.gridY',yobs)

 nobs = size(xobs)
 if (nobs > 5) then
   write(6,*) 'You can only choose 5 observations.'
   stop
 end if
 allocate(zetaobs(nobs,Nsteps))
#endif

 ! grid spacing (in meters)
 dx = 1000
 dy = 1000

 do j = 1,n
   do i = 1,m
     x = dx*(i-1)
     y = dy*(j-1)
     ew = (erf((x - 30e3)/10e3)+1)/2; 
     tmp = (erf((x - 78e3)/10e3)+1)/2;
     tmp = (erf((x - 85e3)/10e3)+1)/2;
     tmp2 = (erf((x - 99e3)/2e3)+1)/2;
     ! depth of domain
     !h(i,j) = 100
     !h(i,j) = 5000* (exp(-(y - 49500)**2 / 20000**2)*ew  - 2*ew+1)

     h(i,j) =  6000*(1-ew) - 100*tmp - 50*tmp2 + 80  & 
          + 40*exp(-(y - 70e3)**2 / 5e3**2)*tmp &
          + 40*exp(-(y - 30e3)**2 / 10e3**2)*tmp; 
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
#ifdef EXTRACT_OBS
 call initgrid(modgrid,dom%x,dom%y,.not.mask)
#endif

#ifdef OAK
 call oak_init(config,'shallow_water2d.init',Nens=Nens)
#endif

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

#ifdef FREE
   xc = 22630.742250000112        
   yc = 50925.772325663514       
   zetac = 1.8713348915660335
#ifif defined(EXTRACT_OBS)
   xc =  5282.5894368342
   yc =  34360.9243460298
   zetac = 1.10340012266332  
#else
   xc = 20e3 + 10e3 * randn()
   yc = 50e3 + 30e3 * randn()
   zetac = randn()
#endif

   lenx = 30e3
   leny = 60e3

   write(6,*) 'xc,yc,zetac',xc,yc,zetac


   do j = 1,n
     do i = 1,m
       x = dx*(i-1)
       y = dy*(j-1)
       zeta(i,j,memberindex) = zetac * exp(-((x-xc)/lenx)**2 - ((y-yc)/leny)**2)
     end do
   end do
   where (.not.mask) zeta(:,:,memberindex) = 0

   ! save initial conditions
   call shallow_water2d_save(dom,timeindex,zeta(:,:,memberindex), &
        U(:,:,memberindex),V(:,:,memberindex),fname,memberindex = memberindex)

 end do



 timeindex = timeindex + 1

 ! time loop
 time: do timecounter = 1,Nsteps
   time2 = dt * timecounter

   ! ensemble loop
   ensemble: do memberindex = 1,Nens

     call shallow_water2d_step(dom,timecounter, & 
          zeta(:,:,memberindex),U(:,:,memberindex),V(:,:,memberindex), &
          dt,g,f,zeta_next,U_next,V_next)


     ! swap time instances
     zeta(:,:,memberindex) = zeta_next
     U(:,:,memberindex) = U_next
     V(:,:,memberindex) = V_next

#    ifdef DIAG     
     if (mod(timecounter,10) == 0) then
       call diag(dom,timecounter,zeta(:,:,memberindex), &
            U(:,:,memberindex),V(:,:,memberindex),g)
     end if
#    endif

     if (mod(timecounter,30) == 0) then
#ifdef EXTRACT_OBS
       obsindex = obsindex+1
       timeobs(obsindex) = time2

       do l = 1,nobs
         call interp(modgrid,zeta(:,:,memberindex),[xobs(l),yobs(l)],zetaobs(l,obsindex),out)
         if (out) zetaobs(l,obsindex) = fillvalue
         write(6,*) 'zetaobs',zetaobs(l,obsindex)
       end do
#endif

       call shallow_water2d_save(dom,timeindex,zeta(:,:,memberindex), &
            U(:,:,memberindex),V(:,:,memberindex),fname,memberindex = memberindex)

       if (memberindex == Nens) timeindex = timeindex + 1
     end if

#    ifdef OAK
     call oak_assim_ens(config,time2,zeta,U,V)
#    endif

   end do ensemble
 end do time

 write(6,*) '2D-shallow water model finished.'

#ifdef EXTRACT_OBS
! write(6,*) 'zetaobs',zetaobs(:,1:obsindex)

 call writeobs('observations.nc',xobs,yobs,timeobs(1:obsindex), &
      zetaobs(:,1:obsindex))

 deallocate(xobs,yobs,zetaobs)
#endif

contains 

#ifdef EXTRACT_OBS
 subroutine writeobs(fname,xobs,yobs,timeobs,zetaobs)
  implicit none
  real, intent(in) :: xobs(:), yobs(:), timeobs(:), zetaobs(:,:)
  character(len=*), intent(in) :: fname

  integer :: ncid, dimid_nobs, dimid_time, varid_zetaobs, varid_xobs, & 
       varid_yobs, varid_timeobs

  call check(nf90_create(fname,ior(nf90_clobber,nf90_netcdf4),ncid))
  call check(nf90_def_dim(ncid, 'nobs',   size(zetaobs,1), dimid_nobs))
  call check(nf90_def_dim(ncid, 'time',   size(zetaobs,2), dimid_time))

  varid_xobs = def_var(ncid,'xobs',[dimid_nobs], &
       standard_name='projection_x_coordinate', &
       units='m')

  varid_yobs = def_var(ncid,'yobs',[dimid_nobs], &
       standard_name='projection_y_coordinate', &
       units='m')

  varid_timeobs = def_var(ncid,'timeobs',[dimid_time], &
       standard_name='time', &
       units='s')

  varid_zetaobs = def_var(ncid,'zetaobs',[dimid_nobs, dimid_time], &
       standard_name='sea_surface_elevation', &
       units='m', &
       fillvalue=fillvalue)
  
  call check(nf90_enddef(ncid))

  call check(nf90_put_var(ncid,varid_xobs,xobs))
  call check(nf90_put_var(ncid,varid_yobs,yobs))
  call check(nf90_put_var(ncid,varid_timeobs,timeobs))
  call check(nf90_put_var(ncid,varid_zetaobs,zetaobs))

  call check(nf90_close(ncid))

 end subroutine writeobs

 function def_var(ncid,varname,dimids,standard_name,units,fillvalue) result(varid)
   implicit none
   integer, intent(in)          :: ncid
   character(len=*), intent(in) :: varname
   integer, intent(in)          :: dimids(:)
   character(len=*), intent(in), optional :: standard_name, units
   real, optional :: fillvalue
   integer :: varid

   call check(nf90_def_var(ncid, varname, nf90_double, &
         dimids, varid))

   if (present(standard_name)) then
     call check(nf90_put_att(ncid, varid, &
          'standard_name', standard_name))
   end if
   if (present(units)) then
     call check(nf90_put_att(ncid, varid, 'units', units))
   end if
   if (present(fillvalue)) then
     call check(nf90_put_att(ncid, varid, '_Fillvalueue', fillvalue))   
   end if
   end function def_var

 subroutine check(status)
  implicit none
  integer, intent ( in) :: status

  if(status /= nf90_noerr) then
    write(6,*) 'NetCDF error: ',trim(nf90_strerror(status))
    stop "Stopped"
  end if
 end subroutine check
#endif


 subroutine init_random_seed()
  implicit none
  integer :: i, n
  integer, dimension(:), allocatable :: seed

  call random_seed(size = n)
  allocate(seed(n))  
  ! fixed arbitrary seed
  seed = 45 * [(i - 1, i = 1, n)]
  call random_seed(put = seed)

  deallocate(seed)
 end subroutine init_random_seed


end program test_shallow_water2d
