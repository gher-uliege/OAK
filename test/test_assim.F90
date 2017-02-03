! Copyright(c) 2002-2016 Alexander Barth and Luc Vandenblucke

program test_assim
 use assimilation
 use matoper
 use ufileformat
 use initfile
 implicit none

 integer, parameter :: imax = 3, jmax = 3

 ! state vector size
 integer, parameter :: n = imax*jmax*2

 ! Ensemble size
 integer, parameter :: Nens = 10

 ! number of observations
 integer, parameter :: m = 1

 real :: x(imax,jmax), y(imax,jmax), z(imax,jmax)
 integer :: mask(imax,jmax)

 real :: xobs(m), yobs(m), zobs(m), yo(m)
 integer :: maskobs(m)

 real :: R(m,m)
 integer :: i,j

 real :: Ef(n,Nens)  ! Forecast ensemble
 real :: xf(n)          ! Initial state (ensemble mean)
 real :: H(m,n)         ! Observation operator
 real :: HEf(m,Nens) ! Observed forecast ensemble
 real :: Ea(n,Nens)  ! Analysis ensemble 
 real :: xa(n)          ! Analysis state (ensemble mean)
 real :: Eap(n,Nens) ! Analysis ensemble perturbations
 real :: Efp(n,Nens) ! Forecast ensemble perturbations

 ! variables for verification
 real :: Pf(n,n)        ! Forecast Ensemble covariance
 real :: K(n,m)         ! Kalman gain
 real :: Pa_check(n,n)  ! Analysis ensemble covariance
 real :: xa_check(n)    ! Analysis state (ensemble mean)

 type(MemLayout) :: ObsML

 ! tolerance for checking
 real, parameter :: tol = 1e-5

 character(len=MaxFNameLength) :: path
 character(len=MaxFNameLength), pointer :: masknames(:), xnames(:), ynames(:), &
      znames(:), names(:)

 integer :: Rtype ! 0=diagonal R, 1=SMWCovar
 integer, parameter :: RCdim = 10
 real :: RD(m), RC(m,RCdim)


 mask = 1
 z = 0
 do j = 1,jmax
   do i = 1,imax
     x(i,j) = i
     y(i,j) = j
   end do
 end do

 do Rtype = 1,2

   if (Rtype == 1) then
     initfname = 'test/test_assim.init'
     write (6,*) '= Test-case with diagonal R ='
   else
     initfname = 'test/test_assim_SMWCovar.init'
     write (6,*) '= Test-case with non-diagonal R (SMWCovar) ='
   end if

   call getInitValue(initfname,'Model.mask',masknames)
   call getInitValue(initfname,'Model.gridX',xnames)
   call getInitValue(initfname,'Model.gridY',ynames)
   call getInitValue(initfname,'Model.gridZ',znames)
   call getInitValue(initfname,'Model.path',path)

   do i = 1,size(masknames)
     call usave(trim(path)//masknames(i),real(mask),DefaultValex)
     call usave(trim(path)//xnames(i),x,DefaultValex)
     call usave(trim(path)//ynames(i),y,DefaultValex)
     call usave(trim(path)//znames(i),z,DefaultValex)
   end do

   deallocate(masknames,xnames,ynames,znames)

   call MemoryLayout('Model.',ModML)

   ! Ef = randn(n,Nens)
   Ef = reshape([(sin(3.*i),i=1,n*Nens)],[n,Nens])

   xobs = x(2,2)
   yobs = y(2,2)
   zobs = z(2,2)
   maskobs = mask(2,2)
   yo = 1


   call getInitValue(initfname,'Obs001.mask',masknames)
   call getInitValue(initfname,'Obs001.gridX',xnames)
   call getInitValue(initfname,'Obs001.gridY',ynames)
   call getInitValue(initfname,'Obs001.gridZ',znames)
   call getInitValue(initfname,'Obs001.path',path)

   i = 1
   call usave(trim(path)//masknames(i),real(maskobs),DefaultValex)
   call usave(trim(path)//xnames(i),xobs,DefaultValex)
   call usave(trim(path)//ynames(i),yobs,DefaultValex)
   call usave(trim(path)//znames(i),zobs,DefaultValex) 

   call MemoryLayout('Obs001.',ObsML)

   ! set observed value
   call saveVector('Obs001.value',ObsML,yo)

   ! obs. error covariance
   if (Rtype == 1) then
     R = 2*eye(m)
     call saveVector('Obs001.rmse',ObsML,sqrt(diag(R)))
   else
     RD = 2
     RC = reshape([(sin(3.*i),i=1,m*RCdim)],[m,RCdim])

     R = diag(RD) + (RC.xt.RC)

     ! set RMSE
     call saveVector('Obs001.SMWCovar.SqrtDiag',ObsML,sqrt(RD))
     call saveEnsemble('Obs001.SMWCovar.ErrorSpace',ObsML,RC)
   end if

   H = 0
   H(1,sub2ind([imax,jmax],[2,2])) = 1


   ! Compute initial ensemble mean state
   xf = sum(Ef,2)/Nens
   Efp = Ef - spread(xf,2,Nens)

   ! Observed ensemble
   HEf = matmul(H,Ef)

   ! For verification, compute ensemble covariance and Kalman gain
   Pf = matmul(Efp,transpose(Efp)) / (Nens-1)

   K = matmul( &
        matmul(Pf, transpose(H)), &
        inv(matmul(H,matmul(Pf,transpose(H))) + R))


   ! Analyzed covariance and mean state
   Pa_check = Pf - matmul(K,matmul(H,Pf))
   xa_check = xf + matmul(K,(yo - matmul(H,xf)))

   ! Perform analysis step
   call init(initfname)
   call assim(1,Ef,Ea)

   xa = sum(Ea,2)/Nens
   Eap = Ea - spread(xa,2,Nens)

   ! check results
   call assert(xa,xa_check,tol, &
        'analysis ensemble mean')

   call assert(matmul(Eap,transpose(Eap)) / (Nens-1.),Pa_check, tol, &
        'analysis ensemble covariance')

   call MemoryLayoutDone(ObsML)
   call done()
 end do
end program test_assim

! Local Variables:
! compile-command: "cd ~/Assim/OAK-nonDiagR && make test/test_assim && test/test_assim"
! End:
