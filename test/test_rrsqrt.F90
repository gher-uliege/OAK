! Copyright(c) 2002-2009 Alexander Barth and Luc Vandenblucke

program test_rrsqrt
 use matoper
 implicit none

 call test_analysis

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
  real :: Ea(n,dim_ens)  ! Analysis ensemble 
  real :: xa(n)          ! Analysis state (ensemble mean)
  real :: Eap(n,dim_ens) ! Analysis ensemble perturbations
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
  
  ! if debug is true, then internal checks are activated
  logical :: debug = .true. 


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

end program test_rrsqrt
