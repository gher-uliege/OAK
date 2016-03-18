
! module for linear shallow water model in two dimensions

module shallow_water2d

 ! stucture describing the model domain
 ! we use a C-grid




 !  +----------------+----------------+---........---+----------------+
 !  |                |                |              |                |
 !  |                |                |              |                |
 !  |      rho       u      rho       u              u      rho       |
 !  |     (1,m)    (1,m)   (2,m)    (2,m)         (n-1,m)  (n,m)      |
 !  |                |                |              |                |
 !  |                |                |              |                |
 !  +-------v--------+-------v--------+---........---+-------v--------+
 !  |    (1,m-1)     |    (2,m-1)     |              |    (n,m-1)     |
 !  .                .                .              .                .
 !  .                .                .              .                .
 !  .                .                .              .                .
 !  .                .                .              .                .
 !  |                |                |              |                |
 !  +-------v--------+-------v--------+---........---+-------v--------+
 !  |     (1,2)      |     (2,2)      |              |     (n,2)      |
 !  |                |                |              |                |
 !  |                |                |              |                |
 !  |      rho       u      rho       u              u      rho       |
 !  |     (1,2)    (1,2)   (2,2)    (2,2)         (n-1,2)  (n,2)      |
 !  |                |                |              |                |
 !  |                |                |              |                |
 !  +-------v--------+-------v--------+---........---+-------v--------+
 !  |     (1,1)      |     (2,1)      |              |     (n,1)      |
 !  |                |                |              |                |
 !  |                |                |              |                |
 !  |      rho       u      rho       u              u      rho       |
 !  |     (1,1)    (1,1)   (2,1)    (2,1)         (n-1,1)  (n,1)      |
 !  |                |                |              |                |
 !  |                |                |              |                |
 !  +----------------+----------------+---........---+----------------+

 type domain
   ! x,y are the horizontal coordinates
   real, allocatable    :: x(:,:), y(:,:)
   ! depth of the ocean at rho, u and v points
   real, allocatable    :: h(:,:), h_u(:,:), h_v(:,:)
   ! mask at rho, u and v points
   ! .true. is sea and .false. is land
   logical, allocatable :: mask(:,:), mask_u(:,:), mask_v(:,:)
   real, allocatable    :: pm(:,:), pm_u(:,:), pm_v(:,:)
   real, allocatable    :: pn(:,:), pn_u(:,:), pn_v(:,:)  
   real, allocatable    :: dS(:,:), dS_u(:,:), dS_v(:,:)  
 end type domain

 type parameters

 end type parameters

 type state
   real, allocatable    :: zeta(:,:), U(:,:), V(:,:)
 end type state

contains

 subroutine init_domain(dom,dx,dy,mask,h)
  implicit none
  type(domain), intent(out) :: dom
  real, intent(in)          :: dx, dy, h(:,:)
  logical, intent(in)       :: mask(:,:)

  integer :: m,n,i,j
  m = size(h,1)
  n = size(h,2)

  allocate( &
       dom%x(m,n),dom%y(m,n), &
       dom%h(m,n),dom%mask(m,n), &
       dom%pm(m,n),dom%pn(m,n), &
       dom%h_u(m-1,n),dom%mask_u(m-1,n), &
       dom%pm_u(m-1,n),dom%pn_u(m-1,n), &
       dom%h_v(m,n-1),dom%mask_v(m,n-1), &
       dom%pm_v(m,n-1),dom%pn_v(m,n-1), &
       dom%dS(m,n) &
       )

  dom%h = h
  dom%mask = mask

  ! inverse of the resolution
  dom%pm = 1./dx
  dom%pn = 1./dy

  dom%dS = 1./(dom%pm*dom%pn)

  ! stagger u
  dom%h_u    = (dom%h(1:m-1,:)      +   dom%h(2:m,:))/2
  dom%pm_u   = (dom%pm(1:m-1,:)     +   dom%pm(2:m,:))/2
  dom%pn_u   = (dom%pn(1:m-1,:)     +   dom%pn(2:m,:))/2
  dom%mask_u =  dom%mask(1:m-1,:) .and. dom%mask(2:m,:)

  ! stagger v
  dom%h_v    = (dom%h(:,1:n-1)      +   dom%h(:,2:n))/2
  dom%pm_v   = (dom%pm(:,1:n-1)     +   dom%pm(:,2:n))/2
  dom%pn_v   = (dom%pn(:,1:n-1)     +   dom%pn(:,2:n))/2
  dom%mask_v =  dom%mask(:,1:n-1) .and. dom%mask(:,2:n)

  ! grid
  do j = 1,n
    do i = 1,m
      dom%x(i,j) = dx*i
      dom%y(i,j) = dy*j
    end do
  end do

 end subroutine init_domain

 subroutine shallow_water2d_step(dom,timecounter,zetai,Ui,Vi,dt,g,f,zeta,U,V)
  implicit none
  type(domain), intent(in) :: dom
  integer, intent(in) :: timecounter
  real, intent(in) :: zetai(:,:), Ui(:,:), Vi(:,:)
  real, intent(out) :: zeta(:,:), U(:,:), V(:,:)
  real, intent(in)   :: dt,g,f

  integer :: m,n
  real :: U2(size(zeta,1)-1,size(zeta,2))
  real :: V2(size(zeta,1),size(zeta,2)-1)
  real :: c(size(zeta,1),size(zeta,2)),dt_max
  real :: t, Up, Vp
  integer :: s, i, j

  ! linear shallow water model
  ! Input:
  !   domain: structure define the model domain (mask, h, pm and pn at density points)
  !   zetai: initial elevation
  !   Ui,Vi: initial transport
  !   N: number of time steps
  !   dt: time steps
  !   g: acceleration due to gravity
  !   f: Coriolis frequency
  !   plot_every: number of time between plots (if negative plots are disabled)
  !   save_every: number of time between plots (if negative plots are disabled)
  !   zetabc: callback function for boundary conditions
  ! 
  ! Output:
  !   zeta: elevation after N time steps
  !   U,V: transport after N time steps

  !  subroutine shallow_water2d(dom,zetai,Ui,Vi,N,dt,g,f,plot_every,save_every,zetabc,output)
  !  implicit none

  c = sqrt(g*dom%h)

  m = size(zetai,1)
  n = size(zetai,2)

  dt_max = 0.2 * minval([minval(1./(sqrt(g*dom%h) * dom%pm)), &
       minval(1./(sqrt(g*dom%h) * dom%pn))])

  ! fprintf('Number of time steps:     %10.2f\n',N)
  ! fprintf('Integration length:       %10.2f days\n',N*dt/24/60/60)
  ! fprintf('Time step:                %10.2f s\n',dt)
  ! fprintf('Time step (max):          %10.2f s\n',dt_max)
  ! fprintf('Dimension in x-direction: %7.0f \n',m)
  ! fprintf('Dimension in y-direction: %7.0f \n',n)



  ! zetam = zeros(m,n)
  ! Um = zeros(m-1,n)
  ! Vm = zeros(m,n-1)

  ! zetav = zeros(m,n)
  ! Uv = zeros(m-1,n)
  ! Vv = zeros(m,n-1)


  zeta = zetai
  U = Ui
  V = Vi

  ! if strcmp(output(end-2:end),'.nc')
  !   nc = netcdf(output,'c')
  !   nc('xi_rho') = m
  !   nc('eta_rho') = n
  !   nc('time') = 0
  !   nc{'ocean_time'} = ncfloat('time')
  !   nc{'zeta'} = ncfloat('time','eta_rho','xi_rho')
  !   ihis = 0
  ! end if

  ! fprintf('   Step    Pot.Energy   Kin.Energy   Tot.Energy      Volumn\n')

  ! do timecounter=1:N
  t = dt*timecounter

  U2 = U / dom%pn_u
  V2 = V / dom%pm_v

  do I = 2,m-1
    do J = 2,n-1

      zeta(I,J) = zeta(I,J) &
           - dt * (U2(I,J) - U2(I-1,J) + V2(I,J) - V2(I,J-1)) &
           *dom%pm(I,J)*dom%pn(I,J)
    end do
  end do

  ! !  vol1 = sum(sum(zeta(2:m-1,2:n-1)))
  ! ! fprintf('%7.0f %12.2f %12.2f %12.2f \n',timecounter,vol0,vol1,vol0-vol1)

  !   zeta = zetabc(zeta,t)

  do s=0,1

    if (mod(timecounter+s,2) == 0) then
      do I = 1,m-1
        do J = 2,n-1

          Vp = (V(I,J-1)+V(I,J)+ V(I+1,J-1)+V(I+1,J))/4

          U(I,J) = U(I,J) + f*dt * Vp & 
               - dt * g*dom%h_u(I,J) * (zeta(I+1,J)-zeta(I,J)) * dom%pm_u(I,J)
        end do
      end do

      where (.not.dom%mask_u) U = 0
    else
      do I = 2,m-1
        do J = 1,n-1      
          Up = (U(I-1,J)+U(I,J)+ U(I-1,J+1)  +U(I,J+1))/4

          V(I,J) = V(I,J) - f*dt * Up &
               - dt * g*dom%h_v(I,J) * (zeta(I,J+1)-zeta(I,J)) * dom%pn_v(I,J)
        end do
      end do

      where (.not.dom%mask_v) V = 0
    end if
  end do

 end subroutine shallow_water2d_step

 subroutine diag(dom,timecounter,zeta,U,V,g)
  implicit none
  type(domain), intent(in) :: dom
  integer, intent(in) :: timecounter
  real, intent(in) :: zeta(:,:), U(:,:), V(:,:)
  real, intent(in)   :: g

  real :: Epot, Ekin, Etot, Vol
  integer :: m,n

  m = size(dom%h,1)
  n = size(dom%h,2) 
  Epot = sum(g * zeta**2)
  Ekin = sum(U**2 / dom%h_u) + sum(V**2 / dom%h_v)
  Etot = Epot + Ekin
  Vol = sum(zeta(2:m-1,2:n-1) * dom%dS(2:m-1,2:n-1))

  write(6,'(I10,4F12.2)') timecounter,Epot,Ekin,Etot,Vol
 end subroutine diag

 subroutine shallow_water2d_save(dom,timeindex,zeta,U,V,fname,memberindex)
  use netcdf
  implicit none
  type(domain), intent(in) :: dom
  real, intent(in) :: zeta(:,:), U(:,:), V(:,:)
  integer, intent(in) :: timeindex
  character(len=*), intent(in) :: fname
  integer, intent(in),  optional  :: memberindex

  real(8) :: fillval = NF90_FILL_DOUBLE

  integer :: ncid, status, &
       dimid_xi, dimid_xi_u, dimid_xi_v,  &
       dimid_eta, dimid_eta_u, dimid_eta_v,  &
       dimid_time, dimid_member, & 
       varid_zeta, varid_ubar, varid_vbar, &
       varid_h, varid_mask, varid_mask_u, varid_mask_v

  logical :: createfile
  integer, allocatable :: dimids(:), startindex(:)

  write(6,*)' save'

  if (present(memberindex)) then
    allocate(dimids(4),startindex(4))
    startindex = [1,1,timeindex,memberindex]
  else
    allocate(dimids(3),startindex(3))
    startindex = [1,1,timeindex]
  end if
  
  if (timeindex == 1) then
    if (present(memberindex)) then
      createfile = memberindex == 1
    else
      createfile = .true.
    end if
  else
    createfile = .false.
  end if

  if (createfile) then
    ! create new file

    status = nf90_create(fname,ior(nf90_clobber,nf90_netcdf4),ncid)
    call check(status)

    ! define the dimensions

    call check(nf90_def_dim(ncid, 'xi_rho', size(zeta,1),   dimid_xi))
    call check(nf90_def_dim(ncid, 'xi_u',   size(zeta,1)-1, dimid_xi_u))
    call check(nf90_def_dim(ncid, 'xi_v',   size(zeta,1),   dimid_xi_v))

    call check(nf90_def_dim(ncid, 'eta_rho', size(zeta,1),   dimid_eta))
    call check(nf90_def_dim(ncid, 'eta_u',   size(zeta,1),   dimid_eta_u))
    call check(nf90_def_dim(ncid, 'eta_v',   size(zeta,1)-1, dimid_eta_v))

    call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid_time))
    if (present(memberindex)) then
      call check(nf90_def_dim(ncid, 'member', NF90_UNLIMITED, dimid_member))
    end if

    ! define variable
    call check(nf90_def_var(ncid, 'h', nf90_double, &
         [dimid_xi, dimid_eta], varid_h))

    call check(nf90_def_var(ncid, 'mask', nf90_double, &
         [dimid_xi, dimid_eta], varid_mask))

    call check(nf90_def_var(ncid, 'mask_u', nf90_double, &
         [dimid_xi_u, dimid_eta_u], varid_mask_u))

    call check(nf90_def_var(ncid, 'mask_v', nf90_double, &
         [dimid_xi_v, dimid_eta_v], varid_mask_v))

    varid_zeta = def_var('zeta',[dimid_xi, dimid_eta], &
         'sea_surface_elevation','m')
    varid_ubar = def_var('ubar',[dimid_xi_u, dimid_eta_u], &
         'barotropic_eastward_sea_water_velocity','m s-1')
    varid_vbar = def_var('vbar',[dimid_xi_v, dimid_eta_v], &
         'barotropic_northward_sea_water_velocity','m s-1')

    call check(nf90_enddef(ncid))

    ! put grid
    call check(nf90_put_var(ncid,varid_h,dom%h))
    call check(nf90_put_var(ncid,varid_mask,merge(1,0,dom%mask)))
    call check(nf90_put_var(ncid,varid_mask_u,merge(1,0,dom%mask_u)))
    call check(nf90_put_var(ncid,varid_mask_v,merge(1,0,dom%mask_v)))


  else
    call check(nf90_open(fname,NF90_WRITE,ncid))
    call check(nf90_inq_varid(ncid, 'zeta', varid_zeta))
    call check(nf90_inq_varid(ncid, 'mask', varid_mask))
    call check(nf90_inq_varid(ncid, 'mask_u', varid_mask_u))
    call check(nf90_inq_varid(ncid, 'mask_v', varid_mask_v))
    call check(nf90_inq_varid(ncid, 'ubar', varid_ubar))
    call check(nf90_inq_varid(ncid, 'vbar', varid_vbar))
  end if

  call check(nf90_put_var(ncid,varid_zeta, &
       merge(real(zeta,8),fillval,dom%mask),start=startindex))
  call check(nf90_put_var(ncid,varid_ubar, &
       merge(real(U/dom%h_u,8),fillval,dom%mask_u),start=startindex))
  call check(nf90_put_var(ncid,varid_vbar, &
       merge(real(V/dom%h_v,8),fillval,dom%mask_v),start=startindex))

  call check(nf90_close(ncid))

  deallocate(dimids,startindex)

 contains
  function def_var(varname,hdims,stdname,units) result(varid)
   implicit none
   character(len=*), intent(in) :: varname
   integer, intent(in)          :: hdims(2)
   character(len=*), intent(in) :: stdname, units
   integer :: varid

   dimids(1:2) = hdims
   dimids(3) = dimid_time
   
   if (present(memberindex)) dimids(4) = dimid_member
   
   call check(nf90_def_var(ncid, varname, nf90_double, &
         dimids, varid))
    call check(nf90_put_att(ncid, varid, &
         'standard_name', stdname))
    call check(nf90_put_att(ncid, varid, 'units', units))
    call check(nf90_put_att(ncid, varid, '_FillValue', fillval))   
   end function def_var

  subroutine check(status)
   implicit none
   integer, intent ( in) :: status

   if(status /= nf90_noerr) then
     write(6,*) 'NetCDF error: ',trim(nf90_strerror(status))
     stop "Stopped"
   end if
  end subroutine check

 end subroutine shallow_water2d_save


end module shallow_water2d


