! cd /home/abarth/Assim/OAK-nonDiagR &&  make test/test_nondiag && test/test_nondiag

#include "../ppdef.h"

#ifdef HAS_CHOLMOD

module mod_spline
 use covariance
 use ndgrid, only: basegrid, grid


 type, extends(grid) :: spline
   ! nelem number of grid points, e.g. 110 for a 10 by 11 grid
   ! nelem should be equal to product(gshape)
   integer :: nelem

   ! inverse of local resolution
   real, allocatable :: pm(:,:)

   ! TODO: use data of parent class
   real, allocatable :: x(:,:)

   ! staggered inverse of local resolution
   ! spm(:,i,j) = inverse of local resolution of dimension i 
   ! staggered in direction j
   real, allocatable :: spm(:,:,:)

   ! snu(:,i) "diffusion coefficient"  along dimension i including metric, staggered
   real, allocatable :: snu(:,:)
   logical, allocatable :: smasked(:,:)

 end type spline



contains

 !_______________________________________________________
 !
 ! shift the field by 1 grid point in the direction dim

 function shift(gshape,masked,f,dim) result(sf)
  use matoper
  implicit none
  integer, intent(in) :: gshape(:)
  logical, intent(in) :: masked(:)
  real, intent(in) :: f(:)
  integer, intent(in) :: dim
  real :: sf(size(masked))

  integer :: l1,l2,subs1(size(gshape)),subs2(size(gshape))

  sf = 0

  do l1 = 1,product(gshape)
    ! get subscripts
    ! sub1 and l1 corresponds to (i,j) if dim = 1
    subs1 = ind2sub(gshape,l1)

    if (subs1(dim) /= 1) then
      ! sub2 and l2 corresponds to (i-1,j) if dim = 1
      subs2 = subs1
      subs2(dim) = subs2(dim)-1
      l2 = sub2ind(gshape,subs2)

      sf(l1) = f(l2)
    end if
  end do

 end function shift

 !_______________________________________________________
 !

 function shift_op(gshape,masked,dim) result(S)
  use matoper
  implicit none
  integer, intent(in) :: gshape(:)
  logical, intent(in) :: masked(:)
  integer, intent(in) :: dim
  type(SparseMatrix) :: S

  integer :: l1,l2,subs1(size(gshape)),subs2(size(gshape)),maxnz,nelem

  nelem = product(gshape)
  S%nz = 0
  S%m = nelem
  S%n = nelem
  maxnz = nelem
  allocate(S%i(maxnz),S%j(maxnz),S%s(maxnz))
  S%s = 0

  do l1 = 1,nelem    
    ! get subscripts
    ! sub1 and l1 corresponds to (i,j) if dim = 1
    subs1 = ind2sub(gshape,l1)

    if (subs1(dim) /= 1) then
      ! sub2 and l2 corresponds to (i-1,j) if dim = 1
      subs2 = subs1
      subs2(dim) = subs2(dim)-1
      l2 = sub2ind(gshape,subs2)

      S%nz = S%nz+1
      S%i(S%nz) = l1
      S%j(S%nz) = l2
      S%s(S%nz) = 1.
    end if
  end do

  !  write(6,*) 'S%s',S%s(1:S%nz)
 end function shift_op

 !_______________________________________________________
 !


 function diff(gshape,masked,f,dim) result(df)
  use matoper
  implicit none
  integer, intent(in) :: gshape(:)
  logical, intent(in) :: masked(:)
  real, intent(in) :: f(:)
  integer, intent(in) :: dim
  real :: df(size(masked))

  integer :: l1,l2,subs1(size(gshape)),subs2(size(gshape))

  df = 0

  do l1 = 1,product(gshape)
    ! get subscripts
    ! sub1 and l1 corresponds to (i,j) if dim = 1
    subs1 = ind2sub(gshape,l1)

    if (subs1(dim) /= gshape(dim)) then
      ! sub2 and l2 corresponds to (i+1,j) if dim = 1
      subs2 = subs1
      subs2(dim) = subs2(dim)+1
      l2 = sub2ind(gshape,subs2)

      if (.not.masked(l1).and..not.masked(l2)) then
        df(l1) = f(l2) - f(l1)
      end if
    end if
  end do

 end function diff

 !_______________________________________________________
 !

 function diff_op(gshape,masked,dim) result(S)
  use matoper
  implicit none
  integer, intent(in) :: gshape(:)
  logical, intent(in) :: masked(:)
  integer, intent(in) :: dim
  type(SparseMatrix) :: S

  integer :: l1,l2,subs1(size(gshape)),subs2(size(gshape)),maxnz,nelem

  nelem = product(gshape)
  S%nz = 0
  S%m = nelem
  S%n = nelem
  maxnz = 2*nelem
  allocate(S%i(maxnz),S%j(maxnz),S%s(maxnz))
  S%s = 0

  do l1 = 1,nelem    
    ! get subscripts
    ! sub1 and l1 corresponds to (i,j) if dim = 1
    subs1 = ind2sub(gshape,l1)

    if (subs1(dim) /= gshape(dim)) then
      ! sub2 and l2 corresponds to (i+1,j) if dim = 1
      subs2 = subs1
      subs2(dim) = subs2(dim)+1
      l2 = sub2ind(gshape,subs2)

      if (.not.masked(l1).and..not.masked(l2)) then
        S%nz = S%nz+1
        S%i(S%nz) = l1
        S%j(S%nz) = l2
        S%s(S%nz) = 1.

        S%nz = S%nz+1
        S%i(S%nz) = l1
        S%j(S%nz) = l1
        S%s(S%nz) = -1.
      end if
    end if

  end do

  !  write(6,*) 'S%s',S%s(1:S%nz)
 end function diff_op


 !_______________________________________________________
 !

 function grad(gshape,masked,pm,f,dim) result(df)
  use matoper
  implicit none
  integer, intent(in) :: gshape(:)
  logical, intent(in) :: masked(:)
  real, intent(in) :: pm(:,:)
  real, intent(in) :: f(:)
  integer, intent(in) :: dim
  real :: df(size(masked))

  df = stagger(gshape,masked,pm(:,dim),dim) * diff(gshape,masked,f,dim)
 end function grad

 !_______________________________________________________
 !


 function grad_op(gshape,masked,pm,dim) result(S)
  use matoper
  implicit none
  integer, intent(in) :: gshape(:)
  logical, intent(in) :: masked(:)
  real, intent(in) :: pm(:,:)
  integer, intent(in) :: dim
  type(SparseMatrix) :: S

  S = spdiag(stagger(gshape,masked,pm(:,dim),dim)) .x. diff_op(gshape,masked,dim)
 end function grad_op

 !_______________________________________________________
 !

 function laplacian(conf,f) result(Lf)
  use matoper
  implicit none
  type(spline), intent(in) :: conf
  real, intent(in) :: f(:)

  real :: Lf(conf%nelem)


  integer :: i
  real, dimension(conf%nelem) :: df, ddf

  ! compte laplacian

  Lf = 0
  do i = 1,conf%n
    df = diff(conf%gshape,conf%masked,f,i)

!    write(6,*) 'df',df,conf%masked

    !where (smasked(:,i)) df = snu(:,i) * df
    df = conf%snu(:,i) * df

    ddf = diff(conf%gshape,conf%smasked(:,i),df,i)
    Lf = Lf + shift(conf%gshape,conf%masked,ddf,i)
  end do

  ! product(pm,2) is the inverse of the volume of a grid cell
  Lf = Lf * product(conf%pm,2)

 ! write(6,*) 'Lf',Lf
 end function laplacian

 !_______________________________________________________
 !

 function laplacian_op(conf) result(Lap)
  use matoper
  implicit none
  type(spline), intent(in) :: conf
  type(SparseMatrix) :: Lap,S,D2
  real :: Lf(conf%nelem)
  integer :: i

  ! compute laplacian

  Lap = spzeros(conf%nelem,conf%nelem)

  Lf = 0
  do i = 1,conf%n
    S = diff_op(conf%gshape,conf%masked,i)
    S = spdiag(conf%snu(:,i)) .x. S
    D2 = diff_op(conf%gshape,conf%smasked(:,i),i) .x. S

    ! add sparse matrices
    Lap = Lap + (shift_op(conf%gshape,conf%masked,i) .x. D2)
  end do

  ! product(pm,2) is the inverse of the volume of a grid cell

  Lap = spdiag(product(conf%pm,2)) .x. Lap
 end function laplacian_op

 !_______________________________________________________
 !

 function divand_kernel_coeff(n,alpha) result(mu)
  use matoper
  implicit none
  integer, intent(in) :: n
  real, intent(in) :: alpha(:)
  real :: mu
  real, parameter :: pi = 4. * atan(1.)

  integer :: m,k
  real :: nu

  m = size(alpha)-1
  
  ! check if alpha is a sequence of binomial coefficients
  do k=0,m
    if (nchoosek(m,k) /= alpha(k+1)) then
      write(6,*) 'unsupported sequence of alpha'
      stop
    end if
  end do

  nu = m-n/2.
  mu = (4*pi)**(n/2.) * gamma(real(m)) / gamma(nu)
 end function divand_kernel_coeff

 !_______________________________________________________
 !
 ! compute 
 ! inv(B)
 ! where B is the background covariance matrix in DIVA

 function invB_op(conf,alpha,len) result(iB)
  use matoper
  implicit none
  type(spline), intent(in) :: conf
  real, intent(in) :: alpha(:)
  real, intent(in) :: len
  ! WE: diagonal matrix with elements equal to the square root of the "volume" 
  ! of each grid cell
  type(SparseMatrix) :: iB,Lap,gradient, WE, WEs, Lapi, Lap_scaled

  real :: coeff, d(conf%nelem), leni
  integer :: i, j

  if (size(alpha) < 1) then
    write(6,*) 'error alpha to small', alpha
  end if

  ! TODO effective dimension
  ! d: inverse of the volume of each grid cell
  d = product(conf%pm,2)
  WE = spdiag(1./sqrt(d))

  ! laplacian
  Lap = laplacian_op(conf)
  ! laplacian to a power
  Lapi = speye(conf%nelem)
  ! len to a power
  leni = 1.

  if (.true.) then
  iB = spzeros(conf%nelem,conf%nelem)

  do j = 1,size(alpha)
    if (mod(j,2) == 1) then
      ! i is odd
 
     ! constrain on laplacian
      ! scale by grid cell
      Lap_scaled = WE .x. Lapi
      iB = iB + ((leni * alpha(j)) .x. (transpose(Lap_scaled).x.Lap_scaled))
      !write(6,*) 'j',j,'lap','iB%nz',iB%nz

    else
      ! i is even

      ! constrain on gradient
      do i = 1,conf%n
        WEs = spdiag(1./sqrt(product(conf%spm(:,:,i),2)))
        gradient = grad_op(conf%gshape,conf%masked,conf%pm,i)
!        write(6,*) 'j',j,'grad',__LINE__,'gradient%nz',gradient%nz
        gradient = gradient.x.Lapi
!        gradient = gradient.x.(speye(conf%nelem))
!        write(6,*) 'j',j,'grad',__LINE__,'gradient%nz',gradient%nz
        ! scale by grid cell
        gradient = WEs .x. gradient
!        write(6,*) 'j',j,'grad',__LINE__,'gradient%nz',gradient%nz
        iB = iB + ((leni * alpha(j)) .x. (transpose(gradient).x.gradient))
      end do

      Lapi = Lapi.x.Lap
    end if

    leni = leni * len
  end do
  else
  ! contrain on value
  iB = spdiag(alpha(1)/d)


  ! contrain on gradient
  do i = 1,conf%n
    WEs = spdiag(1./sqrt(product(conf%spm(:,:,i),2)))
    gradient = grad_op(conf%gshape,conf%masked,conf%pm,i)
    ! scale by grid cell
    gradient = WEs .x. gradient
    iB = iB + (alpha(2) .x. (transpose(gradient).x.gradient))
  end do

  ! contrain on laplacian
  ! scale by grid cell
  Lap = WE .x. Lap
  iB = iB + (alpha(3) .x. (transpose(Lap).x.Lap))
  end if

!  write(6,*) 'iB%nz',iB%nz

  ! normalize iB
  coeff = divand_kernel_coeff(conf%n,alpha)
  iB = (1./coeff) .x. iB



!  write(6,*) 'coeff',coeff,conf%n,alpha

 end function invB_op

 !_______________________________________________________
 !
 ! compute 
 ! x' inv(B) x
 ! where B is the background covariance matrix in DIVA

 function invB(conf,alpha,len,x) result(xiBx)
  use matoper
  implicit none
  type(spline), intent(in) :: conf
  real, intent(in) :: alpha(:), x(:)
  real, intent(in) :: len
  real :: xiBx, tmp(conf%nelem), Lapx(conf%nelem)

  real :: coeff, d(conf%nelem), leni, ds(conf%nelem)
  integer :: i, j

  if (size(alpha) < 1) then
    write(6,*) 'error alpha to small', alpha
  end if

  ! d: inverse of the volume of each grid cell
  d = product(conf%pm,2)
  Lapx = x
  ! len to a power
  leni = 1.

  if (.true.) then
    xiBx = 0.

    do j = 1,size(alpha)
      if (mod(j,2) == 1) then
      ! i is odd
 
        ! constrain on laplacian
        ! scale by grid cell
        xiBx = xiBx + leni * alpha(j) * sum(Lapx**2/d)
      else
        ! i is even
        
        ! constrain on gradient
        do i = 1,conf%n
          tmp = grad(conf%gshape,conf%masked,conf%pm,Lapx,i)
          ds = product(conf%spm(:,:,i),2)
          
          where (.not.conf%smasked(:,i))
            tmp = tmp / sqrt(ds)
          end where
          
          xiBx = xiBx + leni * alpha(j) * sum(tmp**2)
        end do

        Lapx = laplacian(conf,Lapx)
      end if

      leni = leni * len
    end do


  else

  ! contrain on value
  xiBx = alpha(1) * sum(x**2/d)

  ! contrain on gradient
  do i = 1,conf%n
    tmp = grad(conf%gshape,conf%masked,conf%pm,x,i)
    ds = product(conf%spm(:,:,i),2)

    where (.not.conf%smasked(:,i))
      tmp = tmp / sqrt(ds)
    end where
    
    xiBx = xiBx + alpha(2) * sum(tmp**2)
  end do

  ! contrain on laplacian
  tmp = laplacian(conf,x)
  xiBx = xiBx + alpha(3) * sum(tmp**2/d)
  end if

  ! normalize iB
  coeff = divand_kernel_coeff(conf%n,alpha)
  xiBx = xiBx/coeff

 end function invB

 !_______________________________________________________
 !

 function stagger_masked(gshape,masked,dim) result(smasked)
  use matoper
  implicit none

  integer, intent(in) :: dim,gshape(:)
  logical, intent(in) :: masked(:)
  logical :: smasked(size(masked))

  integer :: l1,l2,subs1(size(gshape)),subs2(size(gshape))

  smasked = .true.

  do l1 = 1,product(gshape)
    ! get subscripts
    ! sub1 and l1 corresponds to (i,j) if dim = 1
    subs1 = ind2sub(gshape,l1)

    if (subs1(dim) /= gshape(dim)) then
      ! sub2 and l2 corresponds to (i+1,j) if dim = 1
      subs2 = subs1
      subs2(dim) = subs2(dim)+1
      l2 = sub2ind(gshape,subs2)

      smasked(l1) = masked(l1).or.masked(l2)
    end if
  end do

 end function stagger_masked


 !_______________________________________________________
 !

 function stagger(gshape,masked,f,dim) result(sf)
  use matoper
  implicit none

  integer, intent(in) :: dim,gshape(:)
  logical, intent(in) :: masked(:)
  real, intent(in) :: f(:)
  real :: sf(size(f))

  integer :: l1,l2,subs1(size(gshape)),subs2(size(gshape))

  sf = 0

  do l1 = 1,product(gshape)
    ! get subscripts
    ! sub1 and l1 corresponds to (i,j) if dim = 1
    subs1 = ind2sub(gshape,l1)

    if (subs1(dim) /= gshape(dim)) then
      ! sub2 and l2 corresponds to (i+1,j) if dim = 1
      subs2 = subs1
      subs2(dim) = subs2(dim)+1
      l2 = sub2ind(gshape,subs2)

      if (.not.masked(l1).and..not.masked(l2)) then
        sf(l1) = (f(l2) + f(l1))/2.
      end if
    end if
  end do

 end function stagger

 !_______________________________________________________
 !

 subroutine initspline(conf,gshape,masked,pm,x)
  implicit none
  type(spline) :: conf
  integer, intent(in) :: gshape(:)
  logical, intent(in) :: masked(:)
  real, intent(in) :: pm(:,:), x(:,:)

  integer :: nelem,n,i,j

  n = size(gshape)
  nelem = product(gshape)
  allocate(conf%gshape(n),conf%masked(nelem),conf%pm(nelem,n), &
       conf%x(nelem,n), conf%spm(nelem,n,n), conf%snu(nelem,n), &
       conf%smasked(nelem,n))

  conf%n = n
  conf%nelem = nelem
  conf%gshape = gshape
  conf%masked = masked
  conf%pm = pm
  conf%x = x

  ! precompute
  conf%snu = 1

  do i = 1,n
    ! staggered mask
    conf%smasked(:,i) = stagger_masked(gshape,masked,i)

    ! staggered grid metric
    do j = 1,n
      conf%spm(:,j,i) = stagger(gshape,masked,pm(:,j),i)
    end do

    ! coefficient for laplacian
    do j = 1,n
      if (j == i) then
        conf%snu(:,i) = conf%snu(:,i) * conf%spm(:,j,i)
      else
        conf%snu(:,i) = conf%snu(:,i) / conf%spm(:,j,i)
      end if
    end do

    !    write(6,*) 'conf%snu(:,i)',pack(conf%snu(:,i),.not.conf%smasked(:,i))
  end do

  where (conf%smasked)
    conf%snu = 0
  end where
  !  write(6,*) 'conf%snu(:,i)',conf%snu

 end subroutine initspline


 !_______________________________________________________
 !
 ! should use init_regulargrid
 subroutine initspline_rectdom(conf,gshape,x0,x1,masked)
  use matoper
  implicit none
  type(spline) :: conf
  real, intent(in) :: x0(:), x1(:) ! start and end points
  integer, intent(in) :: gshape(:)
  logical, optional :: masked(:)

  real, allocatable :: pm(:,:), x(:,:)
  logical, allocatable :: masked_(:)

  integer :: i,j,nelem, n, subs(size(gshape))
  n = size(gshape)
  nelem = product(gshape)

  allocate(masked_(nelem),x(nelem,n),pm(nelem,n))

  if (present(masked)) then
    masked_ = masked
  else
    masked_ = .false.
  end if

  do j = 1,nelem
    subs = ind2sub(gshape,j)
    x(j,:) = x0 + (x1-x0) * (subs-1.)/ (gshape-1.) 
  end do

  do i = 1,n
    pm(:,i) = (gshape(i)-1.) / (x1(i)-x0(i)) 
  end do

  call initspline(conf,gshape,masked_,pm,x)

  deallocate(masked_,x,pm)

 end subroutine initspline_rectdom
 !_______________________________________________________
 !

 subroutine donespline(conf)
  implicit none
  type(spline) :: conf

  deallocate(conf%gshape,conf%masked,conf%pm, &
       conf%x, conf%spm, conf%snu, &
       conf%smasked)
 end subroutine donespline

 !_______________________________________________________
 !


end module mod_spline
#endif




!_______________________________________________________
!


program test_nondiag
#ifdef HAS_CHOLMOD

 use mod_spline
 use matoper
  use ndgrid, only: cellgrid, setupgrid, near
 implicit none
 type(SparseMatrix) :: iB

 ! for Rfun and PreCondfun
 type(spline) :: conf
 real :: len, lenCov
 type(cellgrid) :: cg
 type(SparseSolver) :: solver

 ! tolerance for checking
 real :: tol

 if (kind(tol) == 4) then
   tol = 0.1
 else
   tol = 1e-7
 end if

 call test_kernel
 call test_diff
 call test_2d
 call test_2d_mask
 call test_2d_small
 call test_3d
 call test_grad2d

 call test_loc_cov
 ! call test_loc_cov_large
contains

 subroutine test_kernel
  implicit none
  real, parameter :: pi = 4. * atan(1.)
  real :: coeff

  coeff = divand_kernel_coeff(2,[1.,2.,1.])
  call assert(coeff,4*pi,tol,'kernel coeff (2d) 121')

  coeff = divand_kernel_coeff(2,[1.,3.,3.,1.])
  call assert(coeff,8*pi,tol,'kernel coeff (2d) 1331')

  coeff = divand_kernel_coeff(3,[1.,2.,1.])
  call assert(coeff,8*pi,tol,'kernel coeff (3d) 121')

  coeff = divand_kernel_coeff(3,[1.,3.,3.,1.])
  call assert(coeff,32*pi,tol,'kernel coeff (3d) 1331')


 end subroutine test_kernel


 function grad2d(masked,f,pm,dim) result(df)
  implicit none

  logical, intent(in) :: masked(:,:)
  real, intent(in) :: f(:,:),pm(:,:,:)
  integer, intent(in) :: dim
  real :: df(size(f,1),size(f,2))

  integer :: i,j

  df = 0
  if (dim == 1) then
    do j = 1,size(df,2)
      do i = 1,size(df,1)-1
        if (.not.masked(i+1,j).and..not.masked(i,j)) then
          df(i,j) = (f(i+1,j) - f(i,j)) * (pm(i+1,j,1) + pm(i,j,1))/2.
        end if
      end do
    end do
  end if

 end function grad2d

 subroutine test_grad2d
  use matoper
  implicit none
  integer, parameter :: gshape(2) = [3,4]
  real :: x(gshape(1),gshape(2),2), pm(gshape(1),gshape(2),2)
  real :: df(gshape(1),gshape(2))
  logical :: masked(gshape(1),gshape(2))

  integer :: i,j

  do j = 1,gshape(2)
    do i = 1,gshape(1)
      x(i,j,1) = i
      x(i,j,2) = j
    end do
  end do

  pm = 1
  masked = .false.
  df = grad2d(masked,3*x(:,:,1) + 2*x(:,:,2),pm,1)

  call assert(maxval(abs(df(1:gshape(1)-1,:) - 3)),0.,1e-10,'grad x')


 end subroutine test_grad2d


 !_______________________________________________________
 !

 subroutine test_diff
  use matoper
  implicit none
  integer, parameter :: gshape(2) = [3,4]
  real :: x(gshape(1)*gshape(2),2), pm(gshape(1)*gshape(2),2)
  real :: df(gshape(1)*gshape(2)), df2(gshape(1),gshape(2))
  logical :: masked(gshape(1)*gshape(2))
  type(spline) :: conf
  type(SparseMatrix) :: Sx,Sy,Lap,Lap2

  integer :: i,j,l

  l = 1
  do j = 1,gshape(2)
    do i = 1,gshape(1)
      x(l,1) = 2.*i
      x(l,2) = 3.*j
      l = l+1
    end do
  end do

  pm(:,1) = 1./2.
  pm(:,2) = 1./3.
  masked = .false.

  !call initspline(conf,gshape,masked,pm,x)
  call initspline_rectdom(conf,gshape,[2.,3.],[2.*gshape(1), 3.*gshape(2)])
  call assert(conf%x,x,tol,'rectdom x')
  call assert(conf%pm,pm,tol,'rectdom pm')

  df = grad(gshape,masked,pm,3*x(:,1) + 2*x(:,2),1)
  df2 = reshape(df,gshape)    
  call assert(maxval(abs(df2(1:gshape(1)-1,:) - 3)),0.,tol,'gradv x')

  df = grad(gshape,masked,pm,3*x(:,1) + 2*x(:,2),2)
  df2 = reshape(df,gshape)    
  call assert(maxval(abs(df2(:,1:gshape(2)-1) - 2)),0.,tol,'gradv y')



  df = laplacian(conf,3*x(:,1) + 2*x(:,2))
  df2 = reshape(df,gshape)    
  call assert(maxval(abs(df2(2:gshape(1)-1,2:gshape(2)-1))),0.,tol,'laplacian (1)')

  df = laplacian(conf,3*x(:,1)**2 + 2*x(:,2))
  df2 = reshape(df,gshape)    
  call assert(maxval(abs(df2(2:gshape(1)-1,2:gshape(2)-1) - 6)),0.,tol,'laplacian (2)')

  Sx = grad_op(gshape,masked,pm,1)
  df2 = reshape(Sx .x. (3*x(:,1) + 2*x(:,2)) ,gshape)    
  call assert(maxval(abs(df2(1:gshape(1)-1,:) - 3)),0.,tol,'gradv op x')

  Sy = grad_op(gshape,masked,pm,2)
  df2 = reshape(Sy .x. (3*x(:,1) + 2*x(:,2)) ,gshape)    
  call assert(maxval(abs(df2(:,1:gshape(2)-1) - 2)),0.,tol,'gradv op x')

  Lap = laplacian_op(conf)
  df2 = reshape(Lap .x. (3*x(:,1)**2 + 2*x(:,2)) ,gshape)    
  call assert(maxval(abs(df2(2:gshape(1)-1,2:gshape(2)-1) - 6)),0.,tol,'laplacian op (1)')

  Lap2 = Lap.x.Lap
  call assert(full(Lap2),matmul(full(Lap),full(Lap)),tol,'laplacian op (2)')

  Lap2 = transpose(Lap).x.Lap
  call assert(full(Lap2),matmul(transpose(full(Lap)),full(Lap)),tol,'laplacian op (3)')

  call donespline(conf)
 end subroutine test_diff



 !_______________________________________________________
 !

 function fun(x) result(iBx)
  implicit none
  real, intent(in) :: x(:)
  real :: iBx(size(x))

  iBx = iB.x.x
 end function fun

 !_______________________________________________________
 ! 

 subroutine test_symmetry(fi)
  use matoper
  implicit none
  real, intent(in) :: fi(:,:)
  integer :: gshape(2) 

  gshape = shape(fi)

  call assert(fi,fi(gshape(1):1:-1,:),tol,'check symmetry - invert x')
  call assert(fi,fi(:,gshape(2):1:-1),tol,'check symmetry - invert y')
  call assert(fi,transpose(fi),tol,'check symmetry - x <-> y')

 end subroutine test_symmetry

 !_______________________________________________________
 ! 

 subroutine test_symmetry3(fi)
  use matoper
  implicit none
  real, intent(in) :: fi(:,:,:)
  integer :: gshape(3) 

  gshape = shape(fi)

  call assert(fi,fi(gshape(1):1:-1,:,:),tol,'check symmetry - invert x')
  call assert(fi,fi(:,gshape(2):1:-1,:),tol,'check symmetry - invert y')
  call assert(fi,fi(:,:,gshape(3):1:-1),tol,'check symmetry - invert z')


 end subroutine test_symmetry3

 !_______________________________________________________
 ! 

 subroutine test_2d_small
  use matoper
  use ufileformat
  implicit none
  integer, parameter :: gshape(2) = [3,3]
  real, parameter :: x0(2) = [-5.,-5.], x1(2) = [5.,5.]
  real, parameter :: alpha(3) = [1,2,1]
  type(spline) :: conf
  real :: fi(gshape(1),gshape(2))
  real, dimension(gshape(1)*gshape(2)) :: x, Bx, r, kernel
  integer :: i
  ! for pcg
  !  integer :: maxit = 10000, nit
  !  real :: relres
  real :: xiBx, len

  write(6,*) '= test_2d_small ='

  len = 1
  call initspline_rectdom(conf,gshape,x0,x1)
  x = 0
  !x((size(x)+1)/2) = 1
  ! middle of domain
  x(sub2ind(gshape,(gshape+1)/2)) = 1

  iB = invB_op(conf,alpha,len)

  xiBx = invB(conf,alpha,len,x)
  call assert(xiBx,x.x.(iB.x.x),tol,'check invB')

  !  Bx = pcg(fun,x,x,tol=tol,maxit=maxit,nit=nit,relres=relres)
  Bx = symsolve(iB,x)

  call assert(iB.x.Bx,x,tol,'check inv(B)*x')
  !  write(6,*) 'x',x
  !  write(6,*) 'nit,relres',nit,relres
  !  write(6,*) 'Bx',Bx

  

  fi = reshape(Bx,gshape)
  call test_symmetry(fi)

 end subroutine test_2d_small


 !_______________________________________________________
 ! 

 subroutine test_2d
  use matoper
  use ufileformat
  implicit none
  integer, parameter :: gshape(2) = [191,191]
  real, parameter :: x0(2) = [-5.,-5.], x1(2) = [5.,5.]
  real, parameter :: alpha(3) = [1,2,1]
  type(spline) :: conf
  real :: fi(gshape(1),gshape(2))
  real, dimension(gshape(1)*gshape(2)) :: x, Bx, r, kernel
  integer :: i
  ! for pcg
  !  integer :: maxit = 10000, nit
  !  real :: relres
  real :: xiBx
  real :: len = 1
  real :: kernel_tol

  write(6,*) '= test_2d ='

  call initspline_rectdom(conf,gshape,x0,x1)
  x = 0
  !x((size(x)+1)/2) = 1
  ! middle of domain
  x(sub2ind(gshape,(gshape+1)/2)) = 1

  iB = invB_op(conf,alpha,len)

  xiBx = invB(conf,alpha,len,x)
  call assert(xiBx,x.x.(iB.x.x),tol,'check invB')

  !  Bx = pcg(fun,x,x,tol=tol,maxit=maxit,nit=nit,relres=relres)
  Bx = symsolve(iB,x)

  call assert(iB.x.Bx,x,tol,'check inv(B)*x')
  !  write(6,*) 'x',x
  !  write(6,*) 'nit,relres',nit,relres
  !  write(6,*) 'Bx',Bx

  

  fi = reshape(Bx,gshape)
  call test_symmetry(fi)


  r = sqrt(sum(conf%x**2,2))

  ! where (r /= 0)
  !   kernel = r * bessel_k1(r)  ! not yet implemented
  ! elsewhere
  !   kernel = 1
  ! end where

  ! octave code to generate the sequence of numbers
  ! [x,y] = ndgrid(linspace(-5,5,191)); 
  ! r = sqrt(x.^2 + y.^2); 
  ! ff = r.* besselk(1,r) ; 
  ! ff(isnan(ff)) = 1; 
  ! ff(96,96 + [1:50])

  !call usave('fi.dat',fi,-9999.)

  if (kind(tol) == 4) then
    kernel_tol = 0.1
  else
    kernel_tol = 0.004
  end if

  call assert(fi(96,[(96+i,i=1,50)]), &
       [0.99507,0.98409,0.96919,0.95146,0.93163,0.91024,0.88769,0.86431, &
        0.84035,0.81603,0.79152,0.76698,0.74251,0.71822,0.69419,0.67050, &
        0.64719,0.62431,0.60191, 0.58000,0.55862,0.53778,0.51750,0.49778, &
        0.47863,0.46004,0.44203, 0.42459,0.40771,0.39139,0.37562,0.36038, &
        0.34568,0.33150,0.31782,0.30465,0.29195,0.27973, 0.26797,0.25666, &
        0.24577,0.23531,0.22526,0.21560,0.20633,0.19742,0.18887,0.18067, &
        0.17280,0.16525],kernel_tol,'check theoretical 2d-121-kernel')


!  call assert(Bx,kernel,0.006,'check theoretical 2d-121-kernel')



 end subroutine test_2d


 !_______________________________________________________
 ! 

 subroutine test_2d_mask
  use matoper
  use ufileformat
  implicit none
  integer, parameter :: gshape(2) = [191,191]
  real, parameter :: x0(2) = [-3.,-3.], x1(2) = [3.,3.]
  real, parameter :: alpha(3) = [1,2,1]
  type(spline) :: conf
  real :: x(gshape(1)*gshape(2)), Bx(gshape(1)*gshape(2)), fi(gshape(1),gshape(2))
  logical :: masked(gshape(1),gshape(2))
  ! for pcg
  !  integer :: maxit = 10000, nit
  !  real :: relres
  real :: xiBx
  real :: len = 1

  write(6,*) '= test_2d_mask ='

  masked = .false.
  masked(80,:) = .true.
  
  call initspline_rectdom(conf,gshape,x0,x1,masked=reshape(masked,[product(gshape)]))
  x = 0
  !x((size(x)+1)/2) = 1
  ! middle of domain
  x(sub2ind(gshape,(gshape+1)/2)) = 1

  !  call test_symmetry(reshape(x,gshape))

  fi = reshape(laplacian(conf,x),gshape)
  !  write(6,*) 'fi ',fi
  call test_symmetry(fi(2:gshape(1)-1,2:gshape(1)-1))

  iB = invB_op(conf,alpha,len)

  call test_symmetry(reshape(iB.x.x,gshape))

  xiBx = invB(conf,alpha,len,x)
  call assert(xiBx,x.x.(iB.x.x),tol,'check invB with masked')

  !  Bx = pcg(fun,x,x,tol=tol,maxit=maxit,nit=nit,relres=relres)
  Bx = symsolve(iB,x)

  !write(6,*) 'iB.x.Bx- x',maxval(abs((iB.x.Bx) - x)),maxval(x)
  call assert(iB.x.Bx,x,tol,'check inv(B)*x with masked')
  !  write(6,*) 'x',x
  !  write(6,*) 'nit,relres',nit,relres
  !  write(6,*) 'Bx',Bx


  fi = reshape(Bx,gshape)

  where (masked) fi = -9999.
  !call usave('fi_mask.dat',fi,-9999.)


 end subroutine 

 !_______________________________________________________
 ! 

 subroutine test_3d
  use matoper
  use ufileformat
  implicit none
  !integer, parameter :: gshape(3) = [4,4,4]
  !integer, parameter :: gshape(3) = [101,101,48]
  integer, parameter :: gshape(3) = [21,21,9]
  real, parameter :: x0(3) = [-5.,-5.,-5.], x1(3) = [5.,5.,5.]
  real, parameter :: alpha(4) = [1,3,3,1]
  !real, parameter :: alpha(3) = [1,2,1]
  type(spline) :: conf
  real :: fi(gshape(1),gshape(2),gshape(3))
  real, dimension(gshape(1)*gshape(2)*gshape(3)) :: x, Bx, r, kernel
  integer :: i
  ! for pcg
  !  integer :: maxit = 10000, nit
  !  real :: relres
  real :: xiBx
  real :: len = 1

  write(6,*) '= test_3d ='

  call initspline_rectdom(conf,gshape,x0,x1)
  x = 0
  !x((size(x)+1)/2) = 1
  ! middle of domain
  x(sub2ind(gshape,(gshape+1)/2)) = 1

  iB = invB_op(conf,alpha,len)

  xiBx = invB(conf,alpha,len,x)

  call assert(xiBx,x.x.(iB.x.x),tol,'check invB')

  !  Bx = pcg(fun,x,x,tol=tol,maxit=maxit,nit=nit,relres=relres)
  Bx = symsolve(iB,x)

  call assert(iB.x.Bx,x,tol,'check inv(B)*x')
  !  write(6,*) 'x',x
  !  write(6,*) 'nit,relres',nit,relres
  !  write(6,*) 'Bx',Bx

  

  fi = reshape(Bx,gshape)
  call test_symmetry3(fi)

  r = sqrt(sum(conf%x**2,2))

  ! where (r /= 0)
  !   kernel = r * bessel_k1(r)  ! not yet implemented
  ! elsewhere
  !   kernel = 1
  ! end where

  ! octave code to generate the sequence of numbers
  ! [x,y] = ndgrid(linspace(-5,5,191)); 
  ! r = sqrt(x.^2 + y.^2); 
  ! ff = r.* besselk(1,r) ; 
  ! ff(isnan(ff)) = 1; 
  ! ff(96,96 + [1:50])

  ! call assert(fi(96,[(96+i,i=1,50)]), &
  !      [0.99507,0.98409,0.96919,0.95146,0.93163,0.91024,0.88769,0.86431, &
  !       0.84035,0.81603,0.79152,0.76698,0.74251,0.71822,0.69419,0.67050, &
  !       0.64719,0.62431,0.60191, 0.58000,0.55862,0.53778,0.51750,0.49778, &
  !       0.47863,0.46004,0.44203, 0.42459,0.40771,0.39139,0.37562,0.36038, &
  !       0.34568,0.33150,0.31782,0.30465,0.29195,0.27973, 0.26797,0.25666, &
  !       0.24577,0.23531,0.22526,0.21560,0.20633,0.19742,0.18887,0.18067, &
  !       0.17280,0.16525],0.004,'check theoretical 2d-121-kernel')


!  call assert(Bx,kernel,0.006,'check theoretical 2d-121-kernel')

!  call usave('fi.dat',fi,-9999.)


 end subroutine test_3d


 !_______________________________________________________
 ! 

 function Rfun(y) result(Ry)
  implicit none
  real, intent(in) :: y(:)
  real             :: Ry(size(y))
  
  real, dimension(size(y)) :: distnear
  integer, dimension(size(y)) :: ind
  integer :: i, j, k, l, nnz, status
  
  Ry = 0
  do j = 1,conf%nelem
    !call near(cg,conf%x(j,:),conf%x,cdist,2*lenCov,ind,distnear,nnz)
    call near_rectdom(cg,conf%x(j,:),conf%x,cdist,2*lenCov,ind,distnear,nnz)
    do k = 1,nnz
      Ry(j) = Ry(j) + locfun(distnear(k)/lenCov) * y(ind(k))
    end do
  end do


!  Ry = iB .x. y

 end function Rfun


 function PreCondfun(y) result(Ry)
  implicit none
  real, intent(in) :: y(:)
  real             :: Ry(size(y))

  integer          :: status
!  Ry = solver%solve(y,status)
!  Ry = 10*y
  Ry = iB.x.y
 end function PreCondfun

 !_______________________________________________________
 ! 
 ! return list of points with a maximum distance maxdist
 ! for point on a regular grid

 subroutine near_rectdom(cg,x,xpos,distfun,maxdist,ind,dist,n)
  use ndgrid, only: distind
  implicit none
  type(cellgrid), intent(in) :: cg
  real, intent(in) :: x(:),xpos(:,:)
  procedure(distind) :: distfun
  real, intent(in) :: maxdist
  real, intent(out) :: dist(:)
  integer, intent(out) :: ind(:), n

  real :: x0(conf%n), dx(conf%n)
  integer :: i,j,k
  integer :: index0(conf%n), index1(conf%n), subshape(conf%n)

  ! should use near_regulargrid
  x0 = conf%x(1,:)
  dx = (conf%x(conf%nelem,:) - x0) / (conf%gshape-1)


  index0 = max(floor((x-x0 - maxdist)/dx) +1,1)
  index1 = min(ceiling((x-x0 + maxdist)/dx) +1,conf%gshape)
  subshape = index1-index0+1
  n = 0
  
  do k = 1,product(subshape)
    i = sub2ind(conf%gshape,ind2sub(subshape,k) + index0-1)
  
    n = n+1
    dist(n) = cdist(conf%x(i,:),x)
    ind(n) = i
  end do
  
 end subroutine near_rectdom


 subroutine test_loc_cov
  use covariance
  use matoper
  use ufileformat
  implicit none
!  integer, parameter :: gshape(2) = [5,5]
  integer, parameter :: gshape(2) = [11,11]
!  integer, parameter :: gshape(2) = [21,21]
  real, parameter :: x0(2) = [-5.,-5.], x1(2) = [5.,5.]
  real, parameter :: alpha(3) = [1,2,1]
  ! start and end index of sub-region
  integer :: index0(2), index1(2), subshape(2)
  real :: dx(2)
!  type(spline) :: conf

  real :: fi(gshape(1),gshape(2))
  real, dimension(gshape(1)*gshape(2)) :: x, Bx, r, kernel, Cx, Cx2, distnear, x2
  integer, dimension(gshape(1)*gshape(2)) :: ind
  integer :: i, j, k, l, nnz, status
  real :: xiBx, dist
  ! for pcg
  integer :: maxit = 10000, nit
  real :: relres

  write(6,*) '= test_loc_cov ='
  
  call initspline_rectdom(conf,gshape,x0,x1)
  lenCov = 2.
  len = lenCov/2.3

  dx = (x1-x0)/(gshape-1)
  x = 0
  !x((size(x)+1)/2) = 1
  ! middle of domain
  x(sub2ind(gshape,(gshape+1)/2)) = 1

  Cx = 0
  do j = 1,conf%nelem
    do i = 1,conf%nelem
      dist = cdist(conf%x(i,:),conf%x(j,:))
      Cx(j) = Cx(j) + locfun(dist/lenCov) * x(i)
    end do
  end do

  cg = setupgrid(conf%x,[lenCov/10.,lenCov/10.])
  Cx2 = 0
  do j = 1,conf%nelem
!    call near(cg,conf%x(j,:),conf%x,cdist,2*lenCov,ind,distnear,nnz)
    call near_rectdom(cg,conf%x(j,:),conf%x,cdist,2*lenCov,ind,distnear,nnz)
    do k = 1,nnz
      Cx2(j) = Cx2(j) + locfun(distnear(k)/lenCov) * x(ind(k))
    end do


! #ifdef OLD
!     call near(cg,conf%x(j,:),conf%x,cdist,2*lenCov,ind,distnear,nnz)
!     do k = 1,nnz
!       Cx2(j) = Cx2(j) + locfun(distnear(k)/lenCov) * x(ind(k))
!     end do
! #else
!     ! i range = (conf%x(j,:)-x0 - 2*lenCov)/dx +1 ...  (conf%x(j,:)-x0 + 2*lenCov)/dx +1

!     index0 = max(floor((conf%x(j,:)-x0 - 2*lenCov)/dx) +1,1)
!     index1 = min(ceiling((conf%x(j,:)-x0 + 2*lenCov)/dx) +1,conf%gshape)
!     subshape = index1-index0+1

!     do k = 1,product(index1-index0+1)
!       i = sub2ind(conf%gshape,ind2sub(subshape,k) + index0-1)

!       Cx2(j) = Cx2(j) + locfun(cdist(conf%x(i,:),conf%x(j,:))/lenCov) * x(i)
!     end do
! #endif
  end do

  call assert(Cx,Cx2,tol,'Cx')
  call assert(Cx,Rfun(x),tol,'Cx (2)')


  iB = invB_op(conf,alpha,len)
  call solver%init(iB)
  x2 = solver%solve(x,status)

  Bx = pcg(Rfun,x,x,tol=1e-4,maxit=maxit,nit=nit,relres=relres)
  write(6,*) 'nit',nit
  call assert(Rfun(Bx),x,1e-4,'check inv(B)*x without preconditioner')

  Bx = pcg(Rfun,x,x,tol=1e-4,maxit=maxit,nit=nit,relres=relres,pc=PreCondfun)
  write(6,*) 'nit',nit

  call assert(Rfun(Bx),x,1e-4,'check inv(B)*x with preconditioner')


  call solver%done()
  call donespline(conf)

 end subroutine test_loc_cov

 
 !_______________________________________________________
 ! 

 subroutine test_loc_cov_large
  use covariance
  use matoper
  use ufileformat
  implicit none
!  integer, parameter :: gshape(2) = [5,5]
!  integer, parameter :: gshape(2) = [11,11]
  integer, parameter :: gshape(2) = [31,31]
  real, parameter :: x0(2) = [-5.,-5.], x1(2) = [5.,5.]
  real, parameter :: alpha(3) = [1,2,1]
!  type(spline) :: conf

  real :: fi(gshape(1),gshape(2))
  real, dimension(gshape(1)*gshape(2)) :: x, Bx, r, kernel, Cx, Cx2, distnear, x2
  integer, dimension(gshape(1)*gshape(2)) :: ind
  integer :: i, j, k, l, nnz, status
  real :: xiBx, dist
  ! for pcg
  integer :: maxit = 10000, nit
  real :: relres

  write(6,*) '= test_loc_cov_large ='
  
  call initspline_rectdom(conf,gshape,x0,x1)
  lenCov = 2.
  len = lenCov/2.3

  x = 0
  !x((size(x)+1)/2) = 1
  ! middle of domain
  x(sub2ind(gshape,(gshape+1)/2)) = 1


  ! Cx = 0
  ! do j = 1,conf%nelem
  !   do i = 1,conf%nelem
  !     dist = cdist(conf%x(i,:),conf%x(j,:))
  !     Cx(j) = Cx(j) + locfun(dist/lenCov) * x(i)
  !   end do
  ! end do

  cg = setupgrid(conf%x,[lenCov/10.,lenCov/10.])
  ! Cx2 = 0
  ! do j = 1,conf%nelem
  !   call near(cg,conf%x(j,:),conf%x,cdist,2*lenCov,ind,distnear,nnz)
  !   do k = 1,nnz
  !     Cx2(j) = Cx2(j) + locfun(distnear(k)/lenCov) * x(ind(k))
  !   end do
  ! end do


  ! call assert(Cx,Cx2,tol,'Cx')
  ! call assert(Cx,Rfun(x),tol,'Cx (2)')


  iB = invB_op(conf,alpha,len)

  Bx = pcg(Rfun,x,x,tol=1e-4,maxit=maxit,nit=nit,relres=relres,pc=PreCondfun)
  write(6,*) 'nit',nit

  call assert(Rfun(Bx),x,1e-4,'check inv(B)*x with preconditioner')
  call donespline(conf)

 end subroutine test_loc_cov_large
#endif

end program test_nondiag
