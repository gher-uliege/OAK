! Copyright(c) 2002-2016 Alexander Barth and Luc Vandenblucke

! cd ~/Assim/OAK-nonDiagR &&  make test/test_rrsqrt && test/test_rrsqrt

module testing_rrsqrt

end module testing_rrsqrt


program test_rrsqrt
 use matoper
 use testing_rrsqrt
 implicit none
 ! tolerance for checking
 real, parameter :: tol = 1e-8

 call test_analysis
 call test_analysis_covar

contains

 subroutine test_analysis
  use matoper
  use rrsqrt
  implicit none

  ! number of observations
  integer, parameter :: m = 5
  
  ! state vector size
  integer, parameter :: n = 10
  
  ! ensemble size
  integer, parameter :: dim_ens = 3

  real :: y(m)           ! Vector of observations
  real :: Sf(n,dim_ens)  ! 
  real :: Ef(n,dim_ens)  ! Forecast ensemble
  real :: xf(n)          ! Initial state (ensemble mean)
  real :: H(m,n)         ! Observation operator
  real :: HEf(m,dim_ens) ! Observed forecast ensemble
  real :: HSf(m,dim_ens) ! Observed S
  real :: Hxf(m)         ! Observed forecast
  real :: Sa(n,dim_ens)  ! 
  real :: xa(n)          ! Analysis state (ensemble mean)
  real :: Efp(n,dim_ens) ! Forecast ensemble perturbations

  real :: invsqrtR(m)
  ! Observation error covariance matrix
  real :: R(m,m)

  ! variables for verification
  real :: Pf(n,n)        ! Forecast Ensemble covariance
  real :: K(n,m)         ! Kalman gain
  real :: Pa_check(n,n)  ! Analysis ensemble covariance
  real :: xa_check(n)    ! Analysis state (ensemble mean)
  
  ! tolerance for checking
  real, parameter :: tol = 1e-8
  

  ! index 
  integer :: i

  ! Initialize arrays
  y = [(i,i=1,m)]
  Ef = reshape([(sin(3.*i),i=1,n*dim_ens)],[n,dim_ens])
  H = reshape([(i,i=1,m*n)],[m,n])
  
  ! R: diagonal matrix with diagonal elements equal to 2
  R = 2*eye(m)
  invsqrtR = 1./sqrt(2.)

  ! Compute initial ensemble mean state
  xf = sum(Ef,2)/dim_ens
  Efp = Ef - spread(xf,2,dim_ens)
  Sf = Efp / sqrt(dim_ens-1.)

  ! Observed part
  HEf = matmul(H,Ef)
  HSf = matmul(H,Sf)
  Hxf = matmul(H,xf)
  
  ! For verification, compute ensemble covariance and Kalman gain
  Pf = matmul(Efp,transpose(Efp)) / (dim_ens-1)

  K = matmul( &
       matmul(Pf, transpose(H)), &
       inv(matmul(H,matmul(Pf,transpose(H))) + R))
  
  ! Analyzed covariance and mean state
  Pa_check = Pf - matmul(K,matmul(H,Pf))
  xa_check = xf + matmul(K,(y - matmul(H,xf)))

  ! Perform analysis step
  call analysis(xf,Hxf,y,Sf,HSf,invsqrtR, xa,Sa)
  
  ! check results
  call assert(xa,xa_check,tol,'analysis ensemble mean')

!  write(6,*) 'xa',xa
!  write(6,*) 'xa',xa_check
  
  call assert(matmul(Sa,transpose(Sa)),Pa_check, tol, &
       'analysis ensemble variance')

  
 end subroutine test_analysis


 !_______________________________________________________
 !

 subroutine testing_analysis_covar(xf,Sf,H,yo,CovarR)
  use covariance
  use matoper
  use rrsqrt
  implicit none

  real, intent(in) :: xf(:),  H(:,:), yo(:), &
       Sf(:,:)
  class(Covar), intent(in) :: CovarR

  ! Observation error covariance matrix
  real :: R(size(yo),size(yo))
  real :: Pf(size(xf),size(xf))
  real :: K(size(xf),size(yo))
  real :: xa(size(xf))
  real :: Sa(size(xf),size(Sf,2))
  real :: Hxf(size(yo))
  real :: HSf(size(yo),size(Sf,2))

  real :: Pa_check(size(xf),size(xf))
  real :: xa_check(size(xf))

  ! Observed part
  HSf = matmul(H,Sf)
  Hxf = matmul(H,xf)
  R = CovarR%full()

  ! For verification, compute ensemble covariance and Kalman gain
  Pf = matmul(Sf,transpose(Sf))

  K = matmul( &
       matmul(Pf, transpose(H)), &
       inv(matmul(H,matmul(Pf,transpose(H))) + R))
  
  ! Analyzed covariance and mean state
  Pa_check = Pf - matmul(K,matmul(H,Pf))
  xa_check = xf + matmul(K,(yo - matmul(H,xf)))

  call analysis_covar(xf,Hxf,yo,Sf,HSf, CovarR, xa,Sa)

  ! check results
  call assert(xa,xa_check,tol,'analysis ensemble mean')

!  write(6,*) 'xa',xa
!  write(6,*) 'xa',xa_check
  
  call assert(matmul(Sa,transpose(Sa)),Pa_check, tol, &
       'analysis ensemble variance')


 end subroutine testing_analysis_covar

 !_______________________________________________________
 !

 subroutine testing_local_analysis_covar(xf,Sf,H,yo,CovarR)
  use covariance
  use matoper
  use rrsqrt
  implicit none

  real, intent(in) :: xf(:),  H(:,:), yo(:), &
       Sf(:,:)
  class(Covar), intent(in) :: CovarR

  ! Observation error covariance matrix
  real :: R(size(yo),size(yo))
  real :: Pf(size(xf),size(xf))
  real :: K(size(xf),size(yo))
  real :: xa(size(xf))
  real :: Sa(size(xf),size(Sf,2))
  real :: Hxf(size(yo))
  real :: HSf(size(yo),size(Sf,2))

  real :: Pa_check(size(xf),size(xf))
  real :: xa_check(size(xf))

  ! local assimilation
  ! allocate maximum size
  integer :: zoneSize(size(xf)) 
  real :: weight(size(yo))
  logical :: relevantObs(size(yo))
  real, allocatable :: Rloc(:,:)
  
  integer :: i, j, l
  ! number of local observations
  integer :: mloc
  ! index of local observations
  integer, allocatable :: iloc(:)
  ! obs. oper for local observations
  real, allocatable :: Hloc(:,:)

  ! Observed part
  HSf = matmul(H,Sf)
  Hxf = matmul(H,xf)
  R = CovarR%full()

  ! For verification, compute ensemble covariance and Kalman gain
  Pf = matmul(Sf,transpose(Sf))

  K = matmul( &
       matmul(Pf, transpose(H)), &
       inv(matmul(H,matmul(Pf,transpose(H))) + R))
  
  ! Analyzed covariance and mean state
  Pa_check = Pf - matmul(K,matmul(H,Pf))
  xa_check = xf + matmul(K,(yo - matmul(H,xf)))

  call analysis(xf,Hxf,yo,Sf,HSf, 1/sqrt(CovarR%diag()), xa,Sa)

  ! single zone
  zoneSize(1) = size(xf)
  call locAnalysis(zoneSize(1:1),selectAllObservations,xf,Hxf,yo,Sf,HSf, 1/sqrt(CovarR%diag()), xa,Sa)

  ! check results
  call assert(xa,xa_check,tol,'analysis ensemble mean')

!  write(6,*) 'xa',xa
!  write(6,*) 'xa',xa_check
  
  call assert(matmul(Sa,transpose(Sa)),Pa_check, tol, &
       'analysis ensemble variance')

  ! every grid point a separate zone
  zoneSize = [(1,i=1,size(xf))]
  call locAnalysis(zoneSize,selectAllObservations,xf,Hxf,yo,Sf,HSf, 1/sqrt(CovarR%diag()), xa,Sa)
  
  ! check results
  call assert(xa,xa_check,tol,'analysis ensemble mean')
  call assert(matmul(Sa,transpose(Sa)),Pa_check, tol, &
       'analysis ensemble variance')

  ! every grid point a separate zone
  zoneSize = [(1,i=1,size(xf))]
!  call locAnalysis(zoneSize,selectObservations,xf,Hxf,yo,Sf,HSf, 1/sqrt(CovarR%diag()), xa,Sa)
  call locAnalysis_covar(zoneSize,selectObservations,xf,Hxf,yo,Sf,HSf, &
       CovarR,xa,Sa)
  
  ! loop over zones
  j = 1
  do i = 1,size(zoneSize)
    call selectObservations(i,weight,relevantObs)
    l = j+zoneSize(i)-1

    mloc = count(relevantObs)
    allocate(iloc(mloc),Hloc(mloc,size(xf)),Rloc(mloc,mloc))
    ! local indices    
    iloc = pack([(i,i=1,size(yo))],relevantObs)

    Rloc = R(iloc,iloc) * diag(1./(weight(iloc)**2))
    Hloc = H(iloc,:)

    ! ifort 14.0.1 produces wrong results with -O3
    ! xa_check(j:l) = &
    !      xf(j:l) & 
    !      + matmul(                           &
    !          ! Kalmain gain                    &
    !          matmul( &
    !            matmul(Pf(j:l,:), transpose(Hloc)), &
    !            inv(matmul(matmul(Hloc,Pf),transpose(Hloc)) + Rloc)), &
    !          ! innovation vector                               &
    !          yo(iloc) - matmul(Hloc,xf)) 

    xa_check(j:l) = xf(j:l) +  &
         ((Pf(j:l,:).xt.Hloc).x.(inv(((Hloc.x.Pf).xt.Hloc) + Rloc).x. &
         (yo(iloc) - (Hloc.x.xf))))
  
    
    deallocate(iloc,Hloc,Rloc)
    j = j+zoneSize(i)
  end do

  ! check results
  call assert(xa,xa_check,tol,'analysis ensemble mean')
 end subroutine testing_local_analysis_covar

 !_______________________________________________________  
 !

 subroutine selectAllObservations(i,c,relevantObs)
  integer, intent(in) :: i
  real, intent(out) :: c(:)
  logical, intent(out) :: relevantObs(:)
  
  c = 1.
  relevantObs = .true.
 end subroutine selectAllObservations

 !_______________________________________________________  
 !

 subroutine selectObservations(i,c,relevantObs)
  use covariance
  implicit none
  integer, intent(in) :: i
  real, intent(out) :: c(:)
  logical, intent(out) :: relevantObs(:)
 
  ! need to be consitent with test_analysis_covar
  real :: xmod(10), xobs(size(c)), dist
  integer :: j
  xmod = [((j-1.)/(size(xmod)-1.),j=1,size(xmod))]
  xobs = [((j-1.)/(size(xobs)-1.),j=1,size(xobs))]

  ! every grid point is it own zone
  c = locfun(abs(xmod(i) - xobs) / 0.21)
  relevantObs = c /= 0
!  write(6,*) 'i ',i,c
 end subroutine selectObservations

 !_______________________________________________________  
 !


 subroutine test_analysis_covar
  use matoper
  use rrsqrt
  implicit none

  ! number of observations
  integer, parameter :: m = 5
  
  ! state vector size
  integer, parameter :: n = 10
  
  ! ensemble size
  integer, parameter :: dim_ens = 3

  real :: y(m)           ! Vector of observations
  real :: Sf(n,dim_ens)  ! 
  real :: Ef(n,dim_ens)  ! Forecast ensemble
  real :: xf(n)          ! Initial state (ensemble mean)
  real :: H(m,n)         ! Observation operator
  real :: Efp(n,dim_ens) ! Forecast ensemble perturbations

  ! Observation error covariance matrix
  type(DiagCovar) :: DiagCovarR
  type(SMWCovar) :: SMWCovarR

    

  ! index 
  integer :: i

  ! Initialize arrays
  y = [(i,i=1,m)]
  Ef = reshape([(sin(3.*i),i=1,n*dim_ens)],[n,dim_ens])
  H = reshape([(i,i=1,m*n)],[m,n])
  
  ! Compute initial ensemble mean state
  xf = sum(Ef,2)/dim_ens
  Efp = Ef - spread(xf,2,dim_ens)
  Sf = Efp / sqrt(dim_ens-1.)


  ! Testing analysis step

  write(6,*) '= Diagonal error observation covariance = '
  ! R: diagonal matrix with diagonal elements equal to 2
  call DiagCovarR%init([(2.,i=1,m)])

  call testing_analysis_covar(xf,Sf,H,y,DiagCovarR)

  write(6,*) '= SMWCovar error observation covariance = '
  ! R: diagonal matrix with diagonal elements equal to 2
  
  call SMWCovarR%init( &
       [(2.,i=1,m)], &
       reshape([(sin(3.*i),i=1,m*dim_ens)],[m,dim_ens]) &
       )

  call testing_analysis_covar(xf,Sf,H,y,SMWCovarR)
  
  write(6,*) '= Diagonal error observation covariance (local) = '  
  call testing_local_analysis_covar(xf,Sf,H,y,DiagCovarR)

 end subroutine test_analysis_covar

end program test_rrsqrt
