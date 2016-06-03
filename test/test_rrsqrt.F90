! Copyright(c) 2002-2016 Alexander Barth and Luc Vandenblucke

! cd ~/Assim/OAK-nonDiagR &&  make test/test_rrsqrt && test/test_rrsqrt

module testing_rrsqrt

end module testing_rrsqrt


program test_rrsqrt
 use matoper
 use testing_rrsqrt
 implicit none
 ! tolerance for checking
 real :: tol

 if (kind(tol) == 4) then
   tol = 5e-4
 else
   tol = 1e-8
 end if
 
 call test_analysis

contains

 subroutine testing_analysis(xf,Sf,H,yo,CovarR)
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

  call analysis(xf,Hxf,yo,Sf,HSf, CovarR, xa,Sa)

  ! check results
  call assert(xa,xa_check,tol,'analysis ensemble mean')

!  write(6,*) 'xa',xa
!  write(6,*) 'xa',xa_check
  
  call assert(matmul(Sa,transpose(Sa)),Pa_check, tol, &
       'analysis ensemble variance')


 end subroutine testing_analysis

 !_______________________________________________________
 !

 subroutine testing_local_analysis(xf,Sf,H,yo,CovarR)
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
  real :: Pa(size(xf),size(xf))
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
  real, allocatable :: Rloc(:,:), invRloc(:,:)
  
  integer :: i, j, l, ilocobs
  ! number of local observations
  integer :: mloc
  ! index of local observations
  integer, allocatable :: iloc(:)
  ! obs. oper for local observations
  real, allocatable :: Hloc(:,:)

  logical :: localise_obs

  ! Observed part
  HSf = matmul(H,Sf)
  Hxf = matmul(H,xf)
  R = CovarR%full()

  ! For verification, compute ensemble covariance and Kalman gain
  Pf = matmul(Sf,transpose(Sf))

!  call assert(abs(det(Pf)) > tol,'det(Pf) is different from zero')

  K = matmul( &
       matmul(Pf, transpose(H)), &
       inv(matmul(H,matmul(Pf,transpose(H))) + R))
  
  ! Analyzed covariance and mean state
  Pa_check = Pf - matmul(K,matmul(H,Pf))
  xa_check = xf + matmul(K,(yo - matmul(H,xf)))

  call analysis(xf,Hxf,yo,Sf,HSf, CovarR, xa,Sa)

  ! single zone
  zoneSize(1) = size(xf)
  call locAnalysis(zoneSize(1:1),selectAllObservations,xf,Hxf,yo,Sf,HSf, CovarR, xa,Sa)

  ! check results
  call assert(xa,xa_check,tol,'analysis ensemble mean')
  
  call assert(matmul(Sa,transpose(Sa)),Pa_check, tol, &
       'analysis ensemble variance')

  ! every grid point a separate zone
  zoneSize = [(1,i=1,size(xf))]
  call locAnalysis(zoneSize,selectAllObservations,xf,Hxf,yo,Sf,HSf, CovarR, xa,Sa)
  
  ! check results
  call assert(xa,xa_check,tol,'analysis ensemble mean')
  call assert(matmul(Sa,transpose(Sa)),Pa_check, tol, &
       'analysis ensemble variance')


  
  do ilocobs = 1,2
    if (ilocobs == 1) then
      write(6,*) '== without observation localisation == '
      localise_obs = .false.
    else
      write(6,*) '== with observation localisation == '
      localise_obs = .true.
    end if

    ! every grid point a separate zone
    zoneSize = [(1,i=1,size(xf))]
    call locAnalysis(zoneSize,selectObservations,xf,Hxf,yo,Sf,HSf, &
         CovarR,xa,Sa,localise_obs = localise_obs)
  
    ! loop over zones
    j = 1
    do i = 1,size(zoneSize)
      call selectObservations(i,weight,relevantObs)
      l = j+zoneSize(i)-1
      
      if (.not.localise_obs) then
        relevantObs = .true.
      end if
      
      mloc = count(relevantObs)
      allocate(iloc(mloc),Hloc(mloc,size(xf)),Rloc(mloc,mloc),invRloc(mloc,mloc))
      ! local indices    
      iloc = pack([(i,i=1,size(yo))],relevantObs)
      !    Rloc = R(iloc,iloc) * diag(1./(weight(iloc)**2))
      
      Rloc = inv(inv(R(iloc,iloc)) * spread(weight(iloc),dim=2,ncopies=mloc) * spread(weight(iloc),dim=1,ncopies=mloc))
      
      invRloc = inv(R(iloc,iloc)) * spread(weight(iloc),dim=2,ncopies=mloc) * spread(weight(iloc),dim=1,ncopies=mloc)
      
      !Rloc = (diag(1./weight(iloc)).x.R(iloc,iloc)).x.diag(1./weight(iloc))
      !Rloc = R(iloc,iloc)
      
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


      Pa = inv(inv(Pf) + (Hloc.tx.(invRloc.x.Hloc)))

      !    Pa_check = Pf - ((Pf.xt.Hloc).x.(inv(((Hloc.x.Pf).xt.Hloc) + Rloc).x.(Hloc.x.Pf)))
      !    call assert(Pa,Pa_check,tol,'analysis Pa (non-diag)')

      xa_check(j:l) = xf(j:l) +  &
           (Pa(j:l,:).x.(Hloc.tx.(invRloc .x.(yo(iloc) - (Hloc.x.xf)))))

      deallocate(iloc,Hloc,Rloc,invRloc)
      j = j+zoneSize(i)
    end do

    ! check results
    call assert(xa,xa_check,tol,'analysis ensemble mean (non-diag)')
    !write(6,*) 'xa ',xa
    !write(6,*) 'xa_check ',xa_check
  end do
 end subroutine testing_local_analysis

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
 
  ! need to be consitent with test_analysis
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


 subroutine test_analysis
  use matoper
  use rrsqrt
  implicit none

  ! number of observations
  integer, parameter :: m = 5
  
  ! state vector size
  integer, parameter :: n = 10
  
  ! ensemble size (make a full-rank Pf matrix)
  integer, parameter :: dim_ens = 12

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
  Ef = reshape([(sin(3.*i*i),i=1,n*dim_ens)],[n,dim_ens])
  H = reshape([(i,i=1,m*n)],[m,n])
  
  ! Compute initial ensemble mean state
  xf = sum(Ef,2)/dim_ens
  Efp = Ef - spread(xf,2,dim_ens)
  Sf = Efp / sqrt(dim_ens-1.)


  ! Testing analysis step

  write(6,*) '= Diagonal error observation covariance = '
  ! R: diagonal matrix with diagonal elements equal to 2
  call DiagCovarR%init([(2.,i=1,m)])

  call testing_analysis(xf,Sf,H,y,DiagCovarR)

  write(6,*) '= SMWCovar error observation covariance = '
  ! R: diagonal matrix with diagonal elements equal to 2
  
  call SMWCovarR%init( &
       [(2.,i=1,m)], &
       reshape([(sin(3.*i),i=1,m*dim_ens)],[m,dim_ens]) &
       )

  call testing_analysis(xf,Sf,H,y,SMWCovarR)
  
  write(6,*) '= Diagonal observation error covariance (local) = '  
  call testing_local_analysis(xf,Sf,H,y,DiagCovarR)

  write(6,*) '= SMW observation error covariance (local) = '  
  call testing_local_analysis(xf,Sf,H,y,SMWCovarR)

 end subroutine test_analysis

end program test_rrsqrt
