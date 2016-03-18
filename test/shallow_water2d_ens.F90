! cd /home/abarth/Assim/OAK_classic/test; gfortran -Wall  -g -O2 -I/usr/include -I$HOME/Assim/OAK_classic  -fdefault-real-8 -o shallow_water2d_ens  shallow_water2d.F90 shallow_water2d_ens.F90 -L/usr/lib -lnetcdff -lnetcdf -L.. -loak -llapack -lblas  && time ./shallow_water2d_ens

! cd /home/abarth/Assim/OAK_classic/; make; cd /home/abarth/Assim/OAK_classic/test; gfortran -Wall  -g -O2 -I/usr/include -I$HOME/Assim/OAK_classic  -fdefault-real-8 -o shallow_water2d_ens  shallow_water2d.F90 shallow_water2d_ens.F90 -L/usr/lib -lnetcdff -lnetcdf -L.. -loak -llapack -lblas  && time ./shallow_water2d_ens


!#define EXTRACT_OBS
#define OAK

program test_shallow_water2d
 use shallow_water2d
 use matoper, only: randn
#ifdef EXTRACT_OBS
 use ndgrid, only: initgrid, grid, interp
 use netcdf
#endif
#ifdef OAK
 use oak
 use assimilation
#endif

 implicit none
 ! size of time
 integer, parameter :: m = 100, n = 100
 ! size of time steps
 integer, parameter :: Nsteps = 18
 ! size of ensemble
 integer, parameter :: Nens = 10
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
 ! time (seconds since start)
 real(8) :: time2
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

#ifdef EXTRACT_OBS
 integer, parameter :: nobs = 2
 integer :: l, obsindex
 type(grid) :: modgrid 
 real :: timeobs(Nsteps), xobs(nobs),yobs(nobs),zetaobs(nobs,Nsteps)
 logical :: out
 real(8) :: fillvalue = NF90_FILL_DOUBLE
#endif

#ifdef OAK
  type(oakconfig) :: config
!  integer, allocatable :: subdomain(:)
#endif

#ifdef EXTRACT_OBS
 obsindex = 0
 xobs = [50e3, 50e3]
 yobs = [40e3, 60e3]
#endif

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
#ifdef EXTRACT_OBS
 call initgrid(modgrid,dom%x,dom%y,.not.mask)
#endif

#ifdef OAK
 call oak_init(config,'shallow_water2d.init')
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

   xc = 10e3 + 10e3 * randn()
   yc = 50e3 + 10e3 * randn()
   zetac = randn()
   !xc = 2458.0188865562832        
   !yc = 54634.573638063470       
   !zetac = 0.55745057948576626     

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
 time: do timecounter = 1,Nsteps
   time2 = dt * timecounter

   ! ensemble loop
   ensemble: do memberindex = 1,Nens

     call shallow_water2d_step(dom,timecounter, & 
          zeta(:,:,1,memberindex),U(:,:,1,memberindex),V(:,:,1,memberindex), &
          dt,g,f,zeta(:,:,2,memberindex),U(:,:,2,memberindex),V(:,:,2,memberindex))


     ! swap time instances
     zeta(:,:,1,memberindex) = zeta(:,:,2,memberindex)
     U(:,:,1,memberindex) = U(:,:,2,memberindex)
     V(:,:,1,memberindex) = V(:,:,2,memberindex)

#    ifdef DIAG     
     if (mod(timecounter,10) == 0) then
       call diag(dom,timecounter,zeta(:,:,1,memberindex), &
            U(:,:,1,memberindex),V(:,:,1,memberindex),g)
     end if
#    endif

     if (mod(timecounter,30) == 0) then
#ifdef EXTRACT_OBS
       obsindex = obsindex+1
       timeobs(obsindex) = time2

       do l = 1,nobs
         call interp(modgrid,zeta(:,:,1,memberindex),[xobs(l),yobs(l)],zetaobs(l,obsindex),out)
         if (out) zetaobs(l,obsindex) = fillvalue
         write(6,*) 'zetaobs',zetaobs(l,obsindex)
       end do
#endif

!       call shallow_water2d_save(dom,timeindex,zeta(:,:,1,memberindex), &
!            U(:,:,1,memberindex),V(:,:,1,memberindex),fname,memberindex = memberindex)

       if (memberindex == Nens) timeindex = timeindex + 1
     end if

#    ifdef OAK
     call oak_assim_ens(config,time2,zeta(:,:,1,:),U(:,:,1,:),V(:,:,1,:))
#    endif

   end do ensemble
 end do time

 write(6,*) 'zeta(50,50,1,1)',zeta(50,50,1,1)
#ifdef EXTRACT_OBS
 write(6,*) 'zetaobs',zetaobs(:,1:obsindex)

 call writeobs('observations.nc',xobs,yobs,timeobs(1:obsindex), &
      zetaobs(:,1:obsindex))
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


end program test_shallow_water2d
