! Copyright(c) 2002-2009 Alexander Barth and Luc Vandenblucke

! include the fortran preprocessor definitions
#include "ppdef.h"

!#define OPTIM_ZONE_OBS
#define ROTATE_ENSEMBLE
#define COPY_ERRORSPACE

module rrsqrt

contains
 subroutine reduceErrorSpace(S,W,Sr)
  use matoper

  implicit none
  real, intent(in) :: S(:,:), W(size(S,1))
  real, intent(out) :: Sr(:,:)

  real :: lam(size(Sr,2)), &
       V(size(S,2),size(Sr,2))
  integer info
  integer m, n, r

  ! 1 <= i <= m
  ! 1 <= j,jp <= n
  ! 1 <= k <= m
  integer i,j,jp,k
 
  ! upper part of the time covariance 
  real :: P((size(S,2)*(size(S,2)+1))/2)

  integer idummy

  ! working space for sspevx
  real :: work(max(8*size(S,2),size(S,1)))
  integer :: iwork(5*size(S,2)), ifail(size(S,2))

  ! BLAS scalar dot product
  real sdot

  ! LAPACK Machine precision routine
  real slamch

  m = size(S,1)
  n = size(S,2)
  r = size(Sr,2)

  lam = 0

  do jp=1,n
    do j=1,jp
      P(j + (jp-1)*jp/2) = 0

      do i=1,m
        P(j + (jp-1)*jp/2) = P(j + (jp-1)*jp/2) + S(i,j)*W(i)*S(i,jp)
      end do
    end do
  end do

  call sspevx('V','I','U',n,P,-1.,-1.,n-r+1,n,2*SLAMCH('S'),idummy,lam,V,n,work,iwork,ifail,info)

  Sr = (1/sqrt(W)).dx.(S.x.V)

 end subroutine reduceErrorSpace

 !_______________________________________________________
 !
 ! compute the analysis correction (xa-xf) associated to the observation,
 ! their error and the model error space:
 !
 ! K = Pf H^T [H Pf H^T + R]^{-1}
 ! xa_xf = K (yo-Hxf)
 ! Pa = Pf - K H Pf
 !
 ! Pf = Sf Sf^T 
 ! Pa = Sa Sa^T 


 !_______________________________________________________  
 !


 subroutine analysisIncrement(Hxf,yo,Sf,HSf,invsqrtR, xa_xf, Sa, amplitudes)
  use ufileformat
  use matoper

  real, intent(in) :: Hxf(:), yo(:), &
       Sf(:,:), HSf(:,:), invsqrtR(:)
  real, intent(out) :: xa_xf(:)
  real, intent(out), optional :: Sa(:,:)
  real, intent(out), optional :: amplitudes(size(Sf,2))

  real :: lambda(size(Sf,2)), UT(size(Sf,2),size(Sf,2)), &
       ampl(size(Sf,2)), isR(size(invsqrtR))
  real :: dummy(1,1)
  integer :: info,r
#ifdef ROTATE_ENSEMBLE
  real, dimension(size(Sf,2)) :: v,w
  real :: rotationMatrix(size(Sf,2),size(Sf,2))
#endif

#ifdef ALLOCATE_LOCAL_VARS
  real, allocatable :: temp(:,:) 
  integer :: i,j
#endif

  r = size(Sf,2)
  lambda = 0
#ifdef ALLOCATE_LOCAL_VARS
  allocate(temp(size(HSf,1),size(HSf,2)))
  !    temp = invsqrtR.dx.HSf

  do j=1,size(HSf,2)
    do i=1,size(HSf,1)
      temp(i,j) = invsqrtR(i)*HSf(i,j)
    end do
  end do

!  write (6,*) 'sum svd ',sum(temp),sum(HSf),sum(invsqrtR)

  call gesvd('n','a',temp,lambda,dummy,UT,info)

! decomposition satisfies:
! temp.tx.temp = UT.tx.(diag(lambda**2).x.UT)

  deallocate(temp)
#else
  call gesvd('n','a',invsqrtR.dx.HSf,lambda,dummy,UT,info)
#endif
  lambda = (1+lambda**2)**(-1)

  ampl = (UT.tx.(lambda.dx.(UT.x.(HSf.tx.(invsqrtR**2*(yo-Hxf))))))

! check if any element of ampl is nan
  if (any(ampl /= ampl)) then
    write (0,*) 'error in ',__FILE__,__LINE__
    write (0,*) 'ampl ',ampl
    stop
  end if
  xa_xf = Sf.x.ampl
!    write (0,*) 'ampl ',sum(ampl)


#ifndef ROTATE_ENSEMBLE
  if (present(Sa)) Sa = Sf.x.(UT.tx.(sqrt(lambda).dx.UT))
#else
  if (present(Sa)) then
!
! rotate Sa such that the sum of all columns is equal to the sum of all colums of Sf
! If Sf is constructed directly from an emsemble by substracting the ensemble mean,
! then the sum of the colums of Sf is zero. In the following algorithm this property is
! maintained. An ensemble can thus be easely constructed from Sa by adding the new ensemble 
! mean xa.
!
! If no ensembles are used, the this rotation is not necessary.
!

    w = 1./sqrt(1.*r)
    v = UT.tx.(sum(UT,2)/sqrt(lambda))
    v = normate(v)

    Sa = Sf.x.(UT.tx.(sqrt(lambda).dx.(UT.x.RotateVector(w,v))))


  end if
#endif

  if (present(amplitudes)) amplitudes = ampl
 end subroutine analysisIncrement


 !_______________________________________________________
 !

 subroutine analysis(xf,Hxf,yo,Sf,HSf,invsqrtR, xa,Sa, amplitudes)
  use matoper

  real, intent(in) :: xf(:),  Hxf(:), yo(:), &
       Sf(:,:), HSf(:,:), invsqrtR(:)
  real, intent(out) :: xa(:)
  real, intent(out), optional :: Sa(:,:), amplitudes(size(Sf,2))

  call analysisincrement(Hxf,yo,Sf,HSf,invsqrtR, xa, Sa, amplitudes)
  xa = xf+xa
 end subroutine analysis

 !_______________________________________________________
 !

 subroutine biasedAnalysis(gamma,xf,bf,Hxf,Hbf,yo,Sf,HSf,invsqrtR, xa,ba,Sa,amplitudes)
  use matoper
  implicit none
  real, intent(in) :: gamma
  real, intent(in) :: xf(:), bf(:), Hxf(:), yo(:), &
       Sf(:,:), HSf(:,:), invsqrtR(:), Hbf(:)
  real, intent(out) :: xa(:), ba(:), Sa(:,:)
  real, intent(out), optional :: amplitudes(size(Sf,2))

  real :: ampl(size(Sf,2))

#   ifndef ALLOCATE_LOCAL_VARS
  real :: unbiasedxf(size(xf)), unbiasedHxf(size(Hxf)) 
#   else
  real, pointer :: unbiasedxf(:), unbiasedHxf(:) 
  allocate(unbiasedxf(size(xf)),unbiasedHxf(size(Hxf)))
#   endif

  call analysisIncrement(Hxf,yo,Sf,HSf,invsqrtR,ba,amplitudes=ampl)
  ba = bf - gamma * ba

  unbiasedxf = xf-ba
  unbiasedHxf = Hxf-(Hbf - gamma * (HSf.x.ampl))

  call analysisincrement(unbiasedHxf,yo,Sf,HSf,sqrt(1-gamma)*invsqrtR,xa)
  xa = unbiasedxf + xa
  Sa = Sf

#   ifdef ALLOCATE_LOCAL_VARS
  deallocate(unbiasedxf,unbiasedHxf)
#   endif

  if (present(amplitudes)) amplitudes = ampl
 end subroutine biasedAnalysis

 !_______________________________________________________
 !

 subroutine locAnalysisIncrement(zoneSize,selectObservations,Hxf,yo,Sf,HSf,invsqrtR,xa_xf, Sa, amplitudes) 
  use matoper
# ifdef ASSIM_PARALLEL
  use parall
# endif 
  implicit none

  ! interface of the function computing the horizontal correlation between 
  ! the ith component of xf and the jth component of yo

  interface 
    subroutine selectObservations(i,c,relevantObs)
     integer, intent(in) :: i
     real, intent(out) :: c(:)
     logical, intent(out) :: relevantObs(:)
    end subroutine selectObservations
  end interface

! size of each zone
! if a zone is a vertical water column, then all zoneSize 
! are equal to the number of vertical layer

  integer, intent(in) :: zoneSize(:)


  real, intent(in) :: Hxf(:), yo(:), &
       Sf(:,:), HSf(:,:), invsqrtR(:)
  real, intent(out) :: xa_xf(:)
  real, intent(out), optional :: Sa(:,:), amplitudes(size(Sf,2),size(zoneSize))

  real, dimension(size(yo)) :: invsqrtRzone,weight

  logical :: noRelevantObs,relevantObs(size(yo))
  integer :: zi, zi1,zi2,i,j, NZones, nbselectedZones,nbObservations

  integer, dimension(size(zoneSize)) :: startIndex,endIndex

# ifdef OPTIM_ZONE_OBS
  real, dimension(size(yo)) :: yozone
  real, dimension(size(yo),size(Sf,2)) :: HSfzone
  integer :: nObs
# endif

  integer :: baseIndex = 1, i1,i2
  real,allocatable :: temp(:)

!$omp single
  xa_xf = 0
  if (present(amplitudes)) amplitudes = 0
  if (present(Sa)) Sa = Sf
!$omp end single

  NZones = size(zoneSize)
  startIndex(1) = 1
  endIndex(1) = zoneSize(1)

  do zi=2,NZones
    startIndex(zi) = endIndex(zi-1)+1
    endIndex(zi)   = endIndex(zi-1)+zoneSize(zi)
  end do


  nbselectedZones = 0

# ifdef ASSIM_PARALLEL
  zi1 = startZIndex(procnum)
  zi2 =   endZIndex(procnum)

  ! if every vector contains only the local subset
  baseIndex = -startIndex(zi1)+1
# else
  zi1 = 1
  zi2 = Nzones

  ! if every vector is a global vector
  baseIndex = 0
# endif

  ! loop over each zone    
!$omp do schedule(dynamic)
  zonesLoop: do zi=zi1,zi2
    i1 = startIndex(zi) + baseIndex
    i2 =   endIndex(zi) + baseIndex

!    write(stdout,*) 'zi ',zi

    call selectObservations(startIndex(zi),weight,relevantObs)

    nbObservations = count(relevantObs)
    if (nbObservations.eq.0) cycle zonesLoop



#   ifndef OPTIM_ZONE_OBS
    invsqrtRzone = invsqrtR * weight
    if (present(amplitudes)) then
      call analysisIncrement(Hxf,yo,Sf(i1:i2,:),HSf,invsqrtRzone,  &
           xa_xf(i1:i2),Sa(i1:i2,:),amplitudes(:,zi))
    else
      call analysisIncrement(Hxf,yo,Sf(i1:i2,:),HSf,invsqrtRzone,xa_xf(i1:i2),Sa(i1:i2,:))
    end if

#   else

    if (nbObservations.eq.size(yo)) then
      invsqrtRzone = invsqrtR * weight
      call analysisIncrement(Hxf,yo,Sfzone,HSf,invsqrtRzone,xa_xf(i1:i2),Sa(i1:i2,:))
    else
      ! selection only the relevant observations
      nObs = 0
      do j=1,size(yo)
        if (relevantObs(j)) then
          nObs = nObs+1
          yozone(nObs) = yo(j)
          HSfzone(nObs,:) = HSf(j,:)
          invsqrtRzone(nObs) = invsqrtR(j) * weight(j)
        end if
      end do

      call analysisIncrement(Hxf,yozone(1:nObs),Sfzone,HSfzone(1:nObs,:),invsqrtRzone(1:nObs),xa_xf(i1:i2),Sa(i1:i2,:))
      !        write(stdout,*) 'nObs ',nObs,sum(yozone),sum(yo),sum(invsqrtRzone),sum(invsqrtR * weight),sum(HSfzone),sum(HSf)
    end if
#   endif


!$omp critical
    nbselectedZones = nbselectedZones+1
!$omp end critical
  end do zonesLoop

  write(stdout,*) 'nbselectedZones ',nbselectedZones,NZones,sum(xa_xf)

!$omp end master
!$omp barrier

 end subroutine locAnalysisIncrement



 !_______________________________________________________
 !
 subroutine locAnalysis(zoneSize,selectObservations,xf,Hxf,yo,Sf,HSf,invsqrtR, xa,Sa, amplitudes)
  use matoper
  implicit none

  ! interface of the function computing the horizontal correlation between 
  ! the ith component of xf and the jth component of yo

  interface 
    subroutine selectObservations(i,c,relevantObs)
     integer, intent(in) :: i
     real, intent(out) :: c(:)
     logical, intent(out) :: relevantObs(:)
    end subroutine selectObservations
  end interface

  integer :: zoneSize(:)

  real, intent(in) :: xf(:),  Hxf(:), yo(:), &
       Sf(:,:), HSf(:,:), invsqrtR(:)
  real, intent(out) :: xa(:)
  real, intent(out), optional :: Sa(:,:), amplitudes(size(Sf,2),size(zoneSize))

  call locAnalysisIncrement(zoneSize,selectObservations,Hxf,yo,Sf,HSf,invsqrtR,xa,Sa,amplitudes)


!$omp master
  xa = xf + xa
!$omp end master

 end subroutine locAnalysis
 !_______________________________________________________
 !
 ! biasedLocAnalysis:
 !
 ! 
 !
 !_______________________________________________________
 !
 subroutine biasedLocAnalysis(zoneSize,selectObservations, &
      gamma,H,Hshift,xf,bf,Hxf,Hbf,yo,Sf,HSf,invsqrtR, &
      xa,ba,Sa, amplitudes)
  use matoper
  implicit none

  ! interface of the function computing the horizontal correlation between 
  ! the ith component of xf and the jth component of yo

  interface 
    subroutine selectObservations(i,c,relevantObs)
     integer, intent(in) :: i
     real, intent(out) :: c(:)
     logical, intent(out) :: relevantObs(:)
    end subroutine selectObservations
  end interface

  integer :: zoneSize(:)

  real, intent(in) :: gamma

  ! opervation operator

  type(SparseMatrix), intent(in) :: H
  real, intent(in) :: Hshift(:)

  real, intent(in) :: xf(:), bf(:), Hxf(:), yo(:), &
       Sf(:,:), HSf(:,:), invsqrtR(:), Hbf(:)
  real, intent(inout) :: xa(:), ba(:)
  real, intent(out), optional :: Sa(:,:), amplitudes(size(Sf,2))

  ! every thread as a local "unbiasedxf" and "unbiasedHxf"
#   ifndef ALLOCATE_LOCAL_VARS
  real :: unbiasedxf(size(xf)), unbiasedHxf(size(Hxf)) 
#   else
  real, pointer :: unbiasedxf(:), unbiasedHxf(:) 
  allocate(unbiasedxf(size(xf)),unbiasedHxf(size(Hxf)))
#   endif

  call locAnalysisIncrement(zoneSize,selectObservations,Hxf,yo,Sf,HSf,invsqrtR,ba)
!$omp master
  ba = bf - gamma * ba
!$omp end master
!$omp barrier

  ! every thread compute its local copy

  unbiasedxf = xf-ba
  unbiasedHxf = (H.x.unbiasedxf) + Hshift

  call locAnalysisIncrement(zoneSize,selectObservations,unbiasedHxf,yo,Sf,HSf,sqrt(1-gamma)*invsqrtR,xa)

!$omp master
  xa = unbiasedxf + xa
  if (present(Sa)) Sa = Sf
!$omp end master

#   ifdef ALLOCATE_LOCAL_VARS
  deallocate(unbiasedxf,unbiasedHxf)
#   endif

 end subroutine biasedLocAnalysis

 !_______________________________________________________
 !


 !_______________________________________________________
 !
 ! The Mahalanobis length gives a measure of covariance and error 
 ! comptability. 
 !
 ! MahalanobisLength = [(yo - Hx)^T (HS (HS)^T + R)^-1 (yo - Hx)]^(1/2)
 !
 ! formula used:  
 !
 ! (yo - Hx)^T (HS (HS)^T + R)^-1 (yo - Hx)
 !   = a^T a - b^T b
 !
 ! where:
 !   a = R^(-1/2) (yo - Hx)
 !   b = (Lambda^2 + I)^(-1/2) U^T HS^T R^(-1) (yo - Hx)
 !
 ! given the svd transformation:
 !   R^{-1/2} HS = V Lambda U^T
 !
 !
 ! See also:
 ! Mardia,K.V.,Kent,J.T.,Bibby,J.M.,1979.
 ! Multivariate Analysis.Academic Press,Reading MA.  
 ! 
 ! Thacker, W.C., Data-model-error compatibility, 
 ! Ocean Modelling 5 (2003), 233-247


 real function MahalanobisLength(yo_Hx,HS,invsqrtR) 
  use matoper
  implicit none
  real, intent(in) :: yo_Hx(:), HS(:,:), invsqrtR(:)
  real :: UT(size(HS,2),size(HS,2))
  real, dimension(size(yo_Hx)) :: a
  real, dimension(size(HS,2)) :: lambda,b

  real :: dummy(1,1)
  integer :: info,istat

  a = invsqrtR * yo_Hx

  call gesvd('n','a',invsqrtR.dx.HS,lambda,dummy,UT,info)
  !    write(stdout,*) 'lambda ',lambda
  !    write(stdout,*) '1+lambda**2 ',1+lambda**2

  lambda = sqrt(1+lambda**2)**(-1)

  b = lambda.dx.(UT.x.(HS.tx.(invsqrtR**2*(yo_Hx))))

  !    write(stdout,*) 'MahalanobisLength 1 ',sum(a**2)
  !    write(stdout,*) 'MahalanobisLength 2 ',sum(b**2),sum(UT),lambda
  !    write(stdout,*) 'MahalanobisLength 2 ',sum(invsqrtR**2*(yo_Hx)),&
  !         sum((HS.tx.(invsqrtR**2*(yo_Hx)))), &
  !         sum((UT.x.(HS.tx.(invsqrtR**2*(yo_Hx)))))
  !    call flush(99,istat)
  MahalanobisLength = sqrt(sum(a**2) - sum(b**2))
 end function MahalanobisLength

 !_______________________________________________________
 !

 subroutine ensAnalysis(Ef,HEf,yo,invsqrtR,Ea, amplitudes)
  implicit none

  real, intent(in) :: yo(:), invsqrtR(:)

! PERFORMANCE BUG
! avoid copies of Ef !!!!
!

#ifdef COPY_ERRORSPACE
  ! for elegance
      real, intent(in) :: Ef(:,:),  HEf(:,:)
#else
  ! for efficiency
   real, intent(inout) :: Ef(:,:),  HEf(:,:)
#endif

  real, intent(out) :: Ea(:,:)
  real, intent(out), optional :: amplitudes(size(Ef,2))


  ! N: ensemble size
  integer :: N,i
  real, dimension(size(Ef,1)) :: xf,xa
  real, dimension(size(HEf,1)) :: Hxf

#ifdef COPY_ERRORSPACE
  real :: Sf(size(Ef,1),size(Ef,2)), HSf(size(HEf,1),size(HEf,2))
#endif

  N = size(Ef,2)
  xf = sum(Ef,2)/N
  Hxf = sum(HEf,2)/N

#ifdef COPY_ERRORSPACE
  do i=1,N
    Sf(:,i) = (Ef(:,i) - xf)/sqrt(N-1.)
    HSf(:,i) = (HEf(:,i) - Hxf)/sqrt(N-1.)
  end do

  call  analysis(xf,Hxf,yo,Sf,HSf,invsqrtR, xa,Ea, amplitudes)  
#else
  do i=1,N
    Ef(:,i) = (Ef(:,i) - xf)/sqrt(N-1.)
    HEf(:,i) = (HEf(:,i) - Hxf)/sqrt(N-1.)
  end do

  call  analysis(xf,Hxf,yo,Ef,HEf,invsqrtR, xa,Ea, amplitudes)
#endif


  do i=1,N
    Ea(:,i) = xa + sqrt(N-1.) * Ea(:,i)
  end do

 end subroutine ensanalysis
 !_______________________________________________________
 !

 subroutine analysisAnamorph(xf,Hxf,yo,Sf,HSf,invsqrtR,  &
  anamorph,invanamorph, &
       xa,Sa, amplitudes)
  use matoper
  use ufileformat
  implicit none

  real, intent(in) :: xf(:),  Hxf(:), yo(:), &
       Sf(:,:), HSf(:,:), invsqrtR(:)
  real, intent(out) :: xa(:)

  interface 
    subroutine anamorph(x)
     real, intent(inout) :: x(:)
    end subroutine anamorph

    subroutine invanamorph(x)
     real, intent(inout) :: x(:)
    end subroutine invanamorph
  end interface

  integer :: i,N
  real, intent(out), optional :: Sa(:,:), amplitudes(size(Sf,2)+1)    

! transformation matix from RRSQRT to ensemble

  real :: Omega(size(Sf,2)+1,size(Sf,2))

! ensemble
  real, allocatable :: E(:,:),HEf(:,:)

  allocate(E(size(xf),size(Sf,2)+1),HEf(size(yo),size(Sf,2)+1))

  N = size(Sf,2)+1
  call sqrt2ens(xf,Sf,E,Omega)

! ensemble at observation locations
  HEf = spread(Hxf,2,N) + (HSf.xt.Omega)

!  call usave('/u/abarth/Assim/Data2/Ef.u',E,0.)

  do i=1,N
    call anamorph(E(:,i))
  end do

!  call usave('/u/abarth/Assim/Data2/Ef2.u',E,0.)

  !call ensanalysis(Ef,HEf,yo,invsqrtR,Ea, amplitudes)
  call ensanalysis(E,HEf,yo,invsqrtR,E,amplitudes)

!  call usave('/u/abarth/Assim/Data2/Ea2.u',E,0.)

  do i=1,N
    call invanamorph(E(:,i))
  end do

!  call usave('/u/abarth/Assim/Data2/E.u',E,0.)

  call  ens2sqrt(E,xa,Sa) 


 end subroutine analysisAnamorph

 !_______________________________________________________
 !

 ! create orthogonal matrix that rotates w onto v
 !  RotateVector(w,v) * w = v

 function RotateVector(w,v) result(Omega)
  use matoper
  implicit none
  real, intent(in) :: w(:),v(:)
  real :: Omega(size(v),size(v))

  Omega = (v.xt.w) + (perpSpace(v).xt.perpSpace(w))
 end function RotateVector

 !_______________________________________________________
 !
 ! normate vectors

 function normate(v) result(w)
  implicit none
  real, intent(in) :: v(:)
  real :: w(size(v))

  w = v / sqrt(sum(v**2))
 end function normate

 !_______________________________________________________
 !
 ! ensemble covariance
 !
 !

 function enscov(E) result(P)
  use matoper
  implicit none
  real, intent(in) :: E(:,:)
  real :: P(size(E,1),size(E,1))

  integer :: N
  real :: mean(size(E,1))

  N = size(E,2);
  mean = sum(E,2)/N;

  P = (E - spread(mean,2,N)).xt.(E - spread(mean,2,N))

  P = P/(N-1.)
 end function enscov

 !
 ! convert an ensemble to RRSQRT representation
 !
 ! Pham 2001

 subroutine ens2sqrt(E,mean,S) 
  use matoper
  implicit none
  real, intent(in) :: E(:,:)
  real, intent(out) :: mean(size(E,1)), S(size(E,1),size(E,2)-1)

  integer :: N,i
  real :: alpha
  real :: shift(size(E,1))

  N = size(E,2);
  mean = sum(E,2)/N;

  alpha = (-1 + sqrt(1.*N))/(N-1);
  shift = (1-alpha) * mean + alpha * E(:,N);

  do i=1,N-1
    S(:,i) = (E(:,i)-shift)/sqrt(N-1.);
  end do

 end subroutine ens2sqrt


 !
 ! convert RRSQRT representation to an ensemble 
 !
 ! Pham 2001


 subroutine sqrt2ens(mean,S,E,Omega)
  use matoper
  implicit none
  real, intent(in) :: mean(:), S(:,:)
  real, intent(out) ::  E(size(S,1),size(S,2)+1)
  real, intent(out), optional :: Omega(size(S,2)+1,size(S,2))
  integer :: N,i
  real :: alpha
  real :: Om(size(S,2)+1,size(S,2)),w(size(S,2)+1)


  N = size(S,2)+1;

  w = 1./sqrt(1.*N);

  Om =  sqrt(N-1.)* (perpSpace(w).x.randOrthMatrix(N-1)) 

  E = S.xt.Om

  do i=1,N
    E(:,i) = mean + E(:,i)
  end do

  if (present(Omega)) Omega = Om
 end subroutine sqrt2ens

end module rrsqrt

