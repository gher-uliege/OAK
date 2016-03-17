! gfortran -o shallow_water2d shallow_water2d.F90

module shallow_water2d

type domain
  real, allocatable    :: h(:,:), h_u(:,:), h_v(:,:)
  logical, allocatable :: mask(:,:), mask_u(:,:), mask_v(:,:)
  real, allocatable    :: pm(:,:), pm_u(:,:), pm_v(:,:)
  real, allocatable    :: pn(:,:), pn_u(:,:), pn_v(:,:)  
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

  integer :: m,n
  m = size(h,1)
  n = size(h,2)
  
  allocate( &
       dom%h(m,n),dom%mask(m,n), &
       dom%pm(m,n),dom%pn(m,n), &
       dom%h_u(m-1,n),dom%mask_u(m-1,n), &
       dom%pm_u(m-1,n),dom%pn_u(m-1,n), &
       dom%h_v(m,n-1),dom%mask_v(m,n-1), &
       dom%pm_v(m,n-1),dom%pn_v(m,n-1))

  dom%h = h
  dom%mask = mask

  ! inverse of the resolution
  dom%pm = 1./dx
  dom%pn = 1./dy

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

       


 end subroutine init_domain

subroutine shallow_water2d_step(dom,zetai,Ui,Vi,timecounter,dt,g,f,zeta,U,V)
 implicit none
 type(domain), intent(in) :: dom
 real, intent(in) :: zetai(:,:), Ui(:,:), Vi(:,:)
 real, intent(out) :: zeta(:,:), U(:,:), V(:,:)
 integer, intent(in) :: timecounter
 real, intent(in)   :: dt,g,f
 
 integer :: m,n
 real :: dS(size(zeta,1),size(zeta,2))
 real :: U2(size(zeta,1)-1,size(zeta,2))
 real :: V2(size(zeta,1),size(zeta,2)-1)
 real :: c(size(zeta,1),size(zeta,2)),dt_max
 real :: t, Up, Vp, Epot, Ekin, Etot, Vol
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

 dS = 1./(dom%pm*dom%pn)

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

          Vp = (V(I,J-1)+V(I,J)+ V(I+1,J-1)+V(I+1,J))/4;
      
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

  

!   zetam = zetam + zeta
!   Um = Um + U
!   Vm = Vm + V

!   zetav = zetav + zeta.^2
!   Uv = Uv + U.^2
!   Vv = Vv + V.^2
  
!   if mod(timecounter,plot_every) == 0 && plot_every > 0
!     ! close is necessary if matlab is run with -nodisplay (otherwise I get a Segmentation fault)
!     !close
!     !clf
!     hold on
!     sc = 1e3
!     !pcolor(x/sc,y/sc,zeta) shading('flat'),
!     tmp = zeta
!     tmp(mask == 0) = NaN
!     !pcolor(domain.lon_rho,domain.lat_rho,tmp) shading('flat'),
!     pcolor(tmp) shading('flat'),
!     ca = [-.3 1]
!     caxis(ca)
!     colorbar
!     title(['time step ' num2str(timecounter)]) 
    
!     if 0
!       Uc = U(2:end,:)+U(1:end-1,:) 
!       Vc = V(:,2:end)+V(:,1:end-1)
!       h = 4
!       s = 1
!       quiver(x(1:h:end,1:h:end)/sc,y(1:h:end,1:h:end)/sc,s*Uc(1:h:end,1:h:end),s*Vc(1:h:end,1:h:end))
!       axis([-L L -L L]/sc)
!     end
!     drawnow
!     hold off

!     !print('-dpng',sprintf('shallow_water_%06d.png',timecounter/plot_every-1))
!   end if

  if (mod(timecounter,1) == 0) then
    Epot = sum(g * zeta**2)
    Ekin = sum(U**2 / dom%h_u) + sum(V**2 / dom%h_v)
    Etot = Epot + Ekin
    Vol = sum(zeta(2:m-1,2:n-1) * dS(2:m-1,2:n-1))

    write(6,'(I10,4F12.2)') timecounter,Epot,Ekin,Etot,Vol
   end if
! end do

! Uc = U(2:m-1,:)+U(1:m-2,:) 
! Vc = V(:,2:n-1)+V(:,1:n-2)

! [zetam,zetav] = stats(zetam,zetav,N)
! [Um,Uv] = stats(Um,Uv,N)
! [Vm,Vv] = stats(Vm,Vv,N)

! if strcmp(output(end-2:end),'.nc')
!   close(nc)
! end if

end subroutine shallow_water2d_step
end module


program test_shallow_water2d
 use shallow_water2d
 implicit none
 integer, parameter :: m = 100, n = 100, Ntime = 2
 real :: h(m,n), zeta(m,n,Ntime), dx, dy
 real :: U(m-1,n,Ntime), V(m,n-1,Ntime)
 logical :: mask(m,n) 
 type(domain) :: dom
 real :: dt,g,f
 integer :: timecounter, i
 
 mask = .false.
 mask(2:m-1,2:m-1) = .true.;
 mask(10:40,50:60) = .false.;
 mask(70,20) = .false.;

 dx = 1000;
 dy = 1000;

 h = 100
 call init_domain(dom,dx,dy,mask,h)

 zeta = 0
 zeta(4,4,1) = 1
 U = 0
 V = 0
 dt = 2;
 g = 9.81;
 f = 1e-4;

 do timecounter = 1,10
   call shallow_water2d_step(dom,zeta(:,:,1),U(:,:,1),V(:,:,1), &
        timecounter,dt,g,f,zeta(:,:,2),U(:,:,2),V(:,:,2))

   zeta(:,:,1) = zeta(:,:,2)
   U(:,:,1) = U(:,:,2)
   V(:,:,1) = V(:,:,2)
 end do
end program test_shallow_water2d
