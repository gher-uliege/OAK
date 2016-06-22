! Copyright(c) 2002-2016 Alexander Barth and Luc Vandenblucke

! cd /home/abarth/Assim/OAK-nonDiagR &&  make test/test_bootstrap && test/test_bootstrap

program test_boostrap
use rrsqrt
implicit none

call test_analysis

contains

! [varK] = kalman_sens_randperm(X,HX,R,Nsam)
! Kalman gain sensitivity by splitting the ensemble
!
! Input:
!   X: ensemble (n x Nens)
!   HX: observed part of ensemble (m x Nens)
!   R: observation error covariance
!   Nsam: number of samples (100 default)
!
! Output:
!   varK: variance of Kalman gain

function kalman_inc_randperm(X,HX,R,Nsam,yo) result(varinc)
 use covariance
 implicit none
 real, intent(in) :: X(:,:), HX(:,:), yo(:)
 real :: varinc(size(X,1))
 integer :: Nsam
 class(covar) :: R

 real, dimension(size(X,1)) :: K1d,K2d, deltainc
 integer :: n,Nens, ind(size(X,2)), i, N2

 n = size(X,1)
 Nens = size(X,2)

 varinc = 0

 N2 = Nens/2

 do i=1,Nsam
  ind = randperm(Nens)

  call ensAnalysis(X(:,ind(1:N2)),HX(:,ind(1:N2)),yo, R,xa=K1d)
  call ensAnalysis(X(:,ind(N2+1:)),HX(:,ind(N2+1:)),yo, R,xa=K2d)
  
  deltainc = K1d - K2d
  varinc = varinc + deltainc**2
end do

  
varinc = varinc/Nsam
end function kalman_inc_randperm

subroutine ensstat(E,mean,var)
 implicit none
 real, intent(in) :: E(:,:)
 real, optional, intent(out) :: mean(size(E,1))
 real, optional, intent(out) :: var(size(E,1))

 real :: mean_(size(E,1))
 integer :: i

 ! two-pass (for mean and variance) as one-pass can loose precision
 mean_ = sum(E,2) / size(E,2)

 if (present(var)) then
   var = 0

   do i=1,size(E,2)
     var = var + (E(:,i) - mean_)**2
   end do
 end if
 if (present(mean)) mean = mean_
end subroutine ensstat


subroutine bootstrap_analysis(Ef,HEf,yo,R,Nsam,filter,varincf,errredf,varxf,envel)
 use covariance
 implicit none
 real, intent(in) :: Ef(:,:), yo(:), HEf(:,:)
 class(covar), intent(in) :: R
 integer, intent(in) :: Nsam

 real, intent(out), dimension(size(Ef,1)) :: varincf, errredf, varxf, envel

 real, dimension(size(Ef,1)) :: varinc, errred, varxa, incf, xa
 real, dimension(size(HEf,1)) :: d
 real :: Ea(size(Ef,1),size(Ef,2))
 real :: alpha = 1
  interface 
    function filter(x) result(xfilter)
     real, intent(in) :: x(:)
     real :: xfilter(size(x))
    end function filter
  end interface

  d = yo - sum(HEf,2)/size(HEf,2)

  call ensanalysis(Ef,HEf,yo,R,Ea,xa)
  call ensstat(Ea,varxa)
  call ensstat(Ef,varxf)

  errred = varxf - varxa
  
!  [inc,S,Sa] = ensKinc(Ef,HEf,R,d)

  varinc = kalman_inc_randperm(Ef,HEf,R,Nsam,yo)

  varincf = filter(varinc)
  errredf = filter(errred)

  alpha = 1
  envel = exp( - alpha * (varincf / errredf)**2 * (varxf / errredf))

  envel = filter(envel)
  envel = envel/maxval(envel)
end subroutine bootstrap_analysis



 subroutine test_analysis
  use matoper
  use rrsqrt
  implicit none

  ! number of observations
  integer, parameter :: m = 5
  
  ! state vector size
  integer, parameter :: n = 10
  
  ! ensemble size
  integer, parameter :: dim_ens = 4

  real :: yo(m)           ! Vector of observations
  real :: Sf(n,dim_ens)  ! 
  real :: Ef(n,dim_ens)  ! Forecast ensemble
  real :: HEf(m,dim_ens) ! Obs. forecast ensemble
  real :: Ea(n,dim_ens)  ! Analysis ensemble
  real :: xf(n)          ! Initial state (ensemble mean)
  real :: H(m,n)         ! Observation operator
  real :: Efp(n,dim_ens) ! Forecast ensemble perturbations
  real :: xa(n)
  real :: Sa(n,dim_ens)

  real :: varinc(n)
  ! permutation vector
  integer :: ind(dim_ens),ind1(dim_ens/2),ind2(dim_ens/2)

  ! Observation error covariance matrix
  type(DiagCovar) :: R
  type(SMWCovar) :: SMWCovarR
   
  ! index 
  integer :: i,Nsam = 100


  ! Initialize arrays
  yo = [(i,i=1,m)]
  Ef = reshape([(sin(3.*i*i),i=1,n*dim_ens)],[n,dim_ens])
  H = reshape([(i,i=1,m*n)],[m,n])


  
  ! Compute initial ensemble mean state
  xf = sum(Ef,2)/dim_ens
  Efp = Ef - spread(xf,2,dim_ens)
  Sf = Efp / sqrt(dim_ens-1.)

  HEf = H.x.Ef
  ! Testing analysis step

  write(6,*) '# bootstrap analysis'
  ! R: diagonal matrix with diagonal elements equal to 2
  call R%init([(2.,i=1,m)])

  varinc = kalman_inc_randperm(Ef,HEf,R,Nsam,yo)


 end subroutine test_analysis

end program test_boostrap

