! cd /home/abarth/Assim/OAK-nonDiagR &&  make test/test_nondiag && test/test_nondiag

module mod_spline
 use covariance

 type spline
   ! ndim number of dimensions (e.g. 2 for lon/lat)
   ! n number of grid points, e.g. 110 for a 10 by 11 grid
   integer :: ndim, n

   ! size of the domain, n should be equal to product(sz)
   integer, allocatable :: sz(:)

   ! mask is true for "within domain"/"sea" grid points and false 
   ! for "outside the domain"/"land"
   logical, allocatable :: mask(:)

   ! inverse of local resolution
   real, allocatable :: pm(:,:)
   real, allocatable :: x(:,:)

   ! staggered inverse of local resolution
   ! spm(:,i,j) = inverse of local resolution of dimension i 
   ! staggered in direction j
   real, allocatable :: spm(:,:,:)

   ! snu(:,i) "diffusion coefficient"  along dimension i including metric, staggered
   real, allocatable :: snu(:,:)
   logical, allocatable :: smask(:,:)

 end type spline



contains

 !_______________________________________________________
 !
 ! shift the field by 1 grid point in the direction dim

 function shift(sz,mask,f,dim) result(sf)
  use matoper
  implicit none
  integer, intent(in) :: sz(:)
  logical, intent(in) :: mask(:)
  real, intent(in) :: f(:)
  integer, intent(in) :: dim
  real :: sf(size(mask))

  integer :: l1,l2,subs1(size(sz)),subs2(size(sz))

  sf = 0

  do l1 = 1,product(sz)
    ! get subscripts
    ! sub1 and l1 corresponds to (i,j) if dim = 1
    subs1 = ind2sub(sz,l1)

    if (subs1(dim) /= 1) then
      ! sub2 and l2 corresponds to (i-1,j) if dim = 1
      subs2 = subs1
      subs2(dim) = subs2(dim)-1
      l2 = sub2ind(sz,subs2)

      sf(l1) = f(l2)
    end if
  end do

 end function shift

 !_______________________________________________________
 !

 function shift_op(sz,mask,dim) result(S)
  use matoper
  implicit none
  integer, intent(in) :: sz(:)
  logical, intent(in) :: mask(:)
  integer, intent(in) :: dim
  type(SparseMatrix) :: S

  integer :: l1,l2,subs1(size(sz)),subs2(size(sz)),maxnz,n

  n = product(sz)
  S%nz = 0
  S%m = n
  S%n = n
  maxnz = n
  allocate(S%i(maxnz),S%j(maxnz),S%s(maxnz))
  S%s = 0

  do l1 = 1,n    
    ! get subscripts
    ! sub1 and l1 corresponds to (i,j) if dim = 1
    subs1 = ind2sub(sz,l1)

    if (subs1(dim) /= 1) then
      ! sub2 and l2 corresponds to (i-1,j) if dim = 1
      subs2 = subs1
      subs2(dim) = subs2(dim)-1
      l2 = sub2ind(sz,subs2)

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


 function diff(sz,mask,f,dim) result(df)
  use matoper
  implicit none
  integer, intent(in) :: sz(:)
  logical, intent(in) :: mask(:)
  real, intent(in) :: f(:)
  integer, intent(in) :: dim
  real :: df(size(mask))

  integer :: l1,l2,subs1(size(sz)),subs2(size(sz))

  df = 0

  do l1 = 1,product(sz)
    ! get subscripts
    ! sub1 and l1 corresponds to (i,j) if dim = 1
    subs1 = ind2sub(sz,l1)

    if (subs1(dim) /= sz(dim)) then
      ! sub2 and l2 corresponds to (i+1,j) if dim = 1
      subs2 = subs1
      subs2(dim) = subs2(dim)+1
      l2 = sub2ind(sz,subs2)

      if (mask(l1).and.mask(l2)) then
        df(l1) = f(l2) - f(l1)
      end if
    end if
  end do

 end function diff

 !_______________________________________________________
 !

 function diff_op(sz,mask,dim) result(S)
  use matoper
  implicit none
  integer, intent(in) :: sz(:)
  logical, intent(in) :: mask(:)
  integer, intent(in) :: dim
  type(SparseMatrix) :: S

  integer :: l1,l2,subs1(size(sz)),subs2(size(sz)),maxnz,n

  n = product(sz)
  S%nz = 0
  S%m = n
  S%n = n
  maxnz = 2*n
  allocate(S%i(maxnz),S%j(maxnz),S%s(maxnz))
  S%s = 0

  do l1 = 1,n    
    ! get subscripts
    ! sub1 and l1 corresponds to (i,j) if dim = 1
    subs1 = ind2sub(sz,l1)

    if (subs1(dim) /= sz(dim)) then
      ! sub2 and l2 corresponds to (i+1,j) if dim = 1
      subs2 = subs1
      subs2(dim) = subs2(dim)+1
      l2 = sub2ind(sz,subs2)

      if (mask(l1).and.mask(l2)) then
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

 function grad(sz,mask,pm,f,dim) result(df)
  use matoper
  implicit none
  integer, intent(in) :: sz(:)
  logical, intent(in) :: mask(:)
  real, intent(in) :: pm(:,:)
  real, intent(in) :: f(:)
  integer, intent(in) :: dim
  real :: df(size(mask))

  df = stagger(sz,mask,pm(:,dim),dim) * diff(sz,mask,f,dim)
 end function grad

 !_______________________________________________________
 !


 function grad_op(sz,mask,pm,dim) result(S)
  use matoper
  implicit none
  integer, intent(in) :: sz(:)
  logical, intent(in) :: mask(:)
  real, intent(in) :: pm(:,:)
  integer, intent(in) :: dim
  type(SparseMatrix) :: S

  S = spdiag(stagger(sz,mask,pm(:,dim),dim)) .x. diff_op(sz,mask,dim)
 end function grad_op

 !_______________________________________________________
 !

 function laplacian(conf,f) result(Lf)
  use matoper
  implicit none
  type(spline), intent(in) :: conf
  real, intent(in) :: f(:)

  real :: Lf(conf%n)


  integer :: i
  real, dimension(conf%n) :: df, ddf

  ! compte laplacian

  Lf = 0
  do i = 1,conf%ndim
    df = diff(conf%sz,conf%mask,f,i)

    !where (smask(:,i)) df = snu(:,i) * df
    df = conf%snu(:,i) * df

    ddf = diff(conf%sz,conf%smask(:,i),df,i)
    Lf = Lf + shift(conf%sz,conf%mask,ddf,i)
  end do

  ! product(pm,2) is the inverse of the volume of a grid cell
  Lf = Lf * product(conf%pm,2)

 end function laplacian

 !_______________________________________________________
 !

 function laplacian_op(conf) result(Lap)
  use matoper
  implicit none
  type(spline), intent(in) :: conf
  type(SparseMatrix) :: Lap,S,D2
  real :: Lf(conf%n)
  integer :: i

  ! compute laplacian

  Lap = spzeros(conf%n,conf%n)

  Lf = 0
  do i = 1,conf%ndim
    S = diff_op(conf%sz,conf%mask,i)
    S = spdiag(conf%snu(:,i)) .x. S
    D2 = diff_op(conf%sz,conf%smask(:,i),i) .x. S

    ! add sparse matrices
    Lap = Lap + (shift_op(conf%sz,conf%mask,i) .x. D2)
  end do

  ! product(pm,2) is the inverse of the volume of a grid cell

  Lap = spdiag(product(conf%pm,2)) .x. Lap
 end function laplacian_op

 !_______________________________________________________
 !

 function divand_kernel_coeff(ndim,alpha) result(mu)
  use matoper
  implicit none
  integer, intent(in) :: ndim
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

  nu = m-ndim/2.
  mu = (4*pi)**(ndim/2.) * gamma(real(m)) / gamma(nu)
 end function divand_kernel_coeff

 !_______________________________________________________
 !
 ! compute 
 ! inv(B)
 ! where B is the background covariance matrix in DIVA

 function invB_op(conf,alpha) result(iB)
  use matoper
  implicit none
  type(spline), intent(in) :: conf
  real, intent(in) :: alpha(:)
  ! WE: diagonal matrix with elements equal to the square root of the "volume" 
  ! of each grid cell
  type(SparseMatrix) :: iB,Lap,gradient, WE, WEs, Lapi, Lap_scaled

  real :: coeff, d(conf%n)
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
  Lapi = speye(conf%n)

  if (.true.) then
  iB = spzeros(conf%n,conf%n)

  do j = 1,size(alpha)
    if (mod(j,2) == 1) then
      ! i is odd
 
     ! constrain on laplacian
      ! scale by grid cell
      Lap_scaled = WE .x. Lapi
      iB = iB + (alpha(j) .x. (transpose(Lap_scaled).x.Lap_scaled))
      !write(6,*) 'j',j,'lap','iB%nz',iB%nz

    else
      ! i is even

      ! constrain on gradient
      do i = 1,conf%ndim
        WEs = spdiag(1./sqrt(product(conf%spm(:,:,i),2)))
        gradient = grad_op(conf%sz,conf%mask,conf%pm,i)
!        write(6,*) 'j',j,'grad',__LINE__,'gradient%nz',gradient%nz
        gradient = gradient.x.Lapi
!        gradient = gradient.x.(speye(conf%n))
!        write(6,*) 'j',j,'grad',__LINE__,'gradient%nz',gradient%nz
        ! scale by grid cell
        gradient = WEs .x. gradient
!        write(6,*) 'j',j,'grad',__LINE__,'gradient%nz',gradient%nz
        iB = iB + (alpha(j) .x. (transpose(gradient).x.gradient))
      end do

      Lapi = Lapi.x.Lap
    end if
  end do
  else
  ! contrain on value
  iB = spdiag(alpha(1)/d)


  ! contrain on gradient
  do i = 1,conf%ndim
    WEs = spdiag(1./sqrt(product(conf%spm(:,:,i),2)))
    gradient = grad_op(conf%sz,conf%mask,conf%pm,i)
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
  coeff = divand_kernel_coeff(conf%ndim,alpha)
  iB = (1./coeff) .x. iB



!  write(6,*) 'coeff',coeff,conf%ndim,alpha

 end function invB_op

 !_______________________________________________________
 !
 ! compute 
 ! x' inv(B) x
 ! where B is the background covariance matrix in DIVA

 function invB(conf,alpha,x) result(xiBx)
  use matoper
  implicit none
  type(spline), intent(in) :: conf
  real, intent(in) :: alpha(:), x(:)
  real :: xiBx, tmp(conf%n), Lapx(conf%n)

  real :: coeff, d(conf%n), ds(conf%n)
  integer :: i, j

  if (size(alpha) < 1) then
    write(6,*) 'error alpha to small', alpha
  end if

  ! d: inverse of the volume of each grid cell
  d = product(conf%pm,2)
  Lapx = x

  if (.true.) then
    xiBx = 0.

    do j = 1,size(alpha)
      if (mod(j,2) == 1) then
      ! i is odd
 
        ! constrain on laplacian
        ! scale by grid cell
        xiBx = xiBx + alpha(j) * sum(Lapx**2/d)
      else
        ! i is even
        
        ! constrain on gradient
        do i = 1,conf%ndim
          tmp = grad(conf%sz,conf%mask,conf%pm,Lapx,i)
          ds = product(conf%spm(:,:,i),2)
          
          where (conf%smask(:,i))
            tmp = tmp / sqrt(ds)
          end where
          
          xiBx = xiBx + alpha(j) * sum(tmp**2)
        end do

        Lapx = laplacian(conf,Lapx)
      end if
    end do


  else

  ! contrain on value
  xiBx = alpha(1) * sum(x**2/d)

  ! contrain on gradient
  do i = 1,conf%ndim
    tmp = grad(conf%sz,conf%mask,conf%pm,x,i)
    ds = product(conf%spm(:,:,i),2)

    where (conf%smask(:,i))
      tmp = tmp / sqrt(ds)
    end where
    
    xiBx = xiBx + alpha(2) * sum(tmp**2)
  end do

  ! contrain on laplacian
  tmp = laplacian(conf,x)
  xiBx = xiBx + alpha(3) * sum(tmp**2/d)
  end if

  ! normalize iB
  coeff = divand_kernel_coeff(conf%ndim,alpha)
  xiBx = xiBx/coeff

 end function invB

 !_______________________________________________________
 !

 function stagger_mask(sz,mask,dim) result(smask)
  use matoper
  implicit none

  integer, intent(in) :: dim,sz(:)
  logical, intent(in) :: mask(:)
  logical :: smask(size(mask))

  integer :: l1,l2,subs1(size(sz)),subs2(size(sz))

  smask = .false.

  do l1 = 1,product(sz)
    ! get subscripts
    ! sub1 and l1 corresponds to (i,j) if dim = 1
    subs1 = ind2sub(sz,l1)

    if (subs1(dim) /= sz(dim)) then
      ! sub2 and l2 corresponds to (i+1,j) if dim = 1
      subs2 = subs1
      subs2(dim) = subs2(dim)+1
      l2 = sub2ind(sz,subs2)

      smask(l1) = mask(l1).and.mask(l2)
    end if
  end do

 end function stagger_mask


 !_______________________________________________________
 !

 function stagger(sz,mask,f,dim) result(sf)
  use matoper
  implicit none

  integer, intent(in) :: dim,sz(:)
  logical, intent(in) :: mask(:)
  real, intent(in) :: f(:)
  real :: sf(size(f))

  integer :: l1,l2,subs1(size(sz)),subs2(size(sz))

  sf = 0

  do l1 = 1,product(sz)
    ! get subscripts
    ! sub1 and l1 corresponds to (i,j) if dim = 1
    subs1 = ind2sub(sz,l1)

    if (subs1(dim) /= sz(dim)) then
      ! sub2 and l2 corresponds to (i+1,j) if dim = 1
      subs2 = subs1
      subs2(dim) = subs2(dim)+1
      l2 = sub2ind(sz,subs2)

      if (mask(l1).and.mask(l2)) then
        sf(l1) = (f(l2) + f(l1))/2.
      end if
    end if
  end do

 end function stagger

 !_______________________________________________________
 !

 subroutine initspline(conf,sz,mask,pm,x)
  implicit none
  type(spline) :: conf
  integer, intent(in) :: sz(:)
  logical, intent(in) :: mask(:)
  real, intent(in) :: pm(:,:), x(:,:)

  integer :: n,ndim,i,j

  ndim = size(sz)
  n = product(sz)
  allocate(conf%sz(ndim),conf%mask(n),conf%pm(n,ndim), &
       conf%x(n,ndim), conf%spm(n,ndim,ndim), conf%snu(n,ndim), &
       conf%smask(n,ndim))

  conf%ndim = ndim
  conf%n = n
  conf%sz = sz
  conf%mask = mask
  conf%pm = pm
  conf%x = x

  ! precompute
  conf%snu = 1

  do i = 1,ndim
    ! staggered mask
    conf%smask(:,i) = stagger_mask(sz,mask,i)

    ! staggered grid metric
    do j = 1,ndim
      conf%spm(:,j,i) = stagger(sz,mask,pm(:,j),i)
    end do

    ! coefficient for laplacian
    do j = 1,ndim
      if (j == i) then
        conf%snu(:,i) = conf%snu(:,i) * conf%spm(:,j,i)
      else
        conf%snu(:,i) = conf%snu(:,i) / conf%spm(:,j,i)
      end if
    end do

    !    write(6,*) 'conf%snu(:,i)',pack(conf%snu(:,i),.not.conf%smask(:,i))
  end do

  where (.not.conf%smask)
    conf%snu = 0
  end where
  !  write(6,*) 'conf%snu(:,i)',conf%snu

 end subroutine initspline


 !_______________________________________________________
 !

 subroutine initspline_rectdom(conf,sz,x0,x1,mask)
  use matoper
  implicit none
  type(spline) :: conf
  real, intent(in) :: x0(:), x1(:) ! start and end points
  integer, intent(in) :: sz(:)
  logical, optional :: mask(:)

  real, allocatable :: pm(:,:), x(:,:)
  logical, allocatable :: mask_(:)

  integer :: i,j,n, ndim, subs(size(sz))
  ndim = size(sz)
  n = product(sz)

  allocate(mask_(n),x(n,ndim),pm(n,ndim))

  if (present(mask)) then
    mask_ = mask
  else
    mask_ = .true.
  end if

  do j = 1,n
    subs = ind2sub(sz,j)
    x(j,:) = x0 + (x1-x0) * (subs-1.)/ (sz-1.) 
  end do

  do i = 1,ndim
    pm(:,i) = (sz(i)-1.) / (x1(i)-x0(i)) 
  end do

  call initspline(conf,sz,mask_,pm,x)

  deallocate(mask_,x,pm)

 end subroutine initspline_rectdom
 !_______________________________________________________
 !

 subroutine donespline(conf)
  implicit none
  type(spline) :: conf

  deallocate(conf%sz,conf%mask,conf%pm, &
       conf%x, conf%spm, conf%snu, &
       conf%smask)
 end subroutine donespline

 !_______________________________________________________
 !


end module mod_spline




!_______________________________________________________
!

program test_nondiag
 use mod_spline
 use matoper
  use ndgrid, only: cellgrid, setupgrid, near
 implicit none
 type(SparseMatrix) :: iB

 ! for Rfun and PreCondfun
 type(spline) :: conf
 real :: len
 type(cellgrid) :: cg
 type(SparseSolver) :: solver

 real :: tol = 1e-7

 call test_kernel
 call test_diff()
 call test_2d()
 call test_2d_mask()
 call test_2d_small()
 call test_3d()
 call test_grad2d()

! call test_loc_cov
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


 function grad2d(mask,f,pm,dim) result(df)
  implicit none

  logical, intent(in) :: mask(:,:)
  real, intent(in) :: f(:,:),pm(:,:,:)
  integer, intent(in) :: dim
  real :: df(size(f,1),size(f,2))

  integer :: i,j

  df = 0
  if (dim == 1) then
    do j = 1,size(df,2)
      do i = 1,size(df,1)-1
        if (mask(i+1,j).and.mask(i,j)) then
          df(i,j) = (f(i+1,j) - f(i,j)) * (pm(i+1,j,1) + pm(i,j,1))/2.
        end if
      end do
    end do
  end if

 end function grad2d

 subroutine test_grad2d
  use matoper
  implicit none
  integer, parameter :: sz(2) = [3,4]
  real :: x(sz(1),sz(2),2), pm(sz(1),sz(2),2)
  real :: df(sz(1),sz(2))
  logical :: mask(sz(1),sz(2))

  integer :: i,j

  do j = 1,sz(2)
    do i = 1,sz(1)
      x(i,j,1) = i
      x(i,j,2) = j
    end do
  end do

  pm = 1
  mask = .true.
  df = grad2d(mask,3*x(:,:,1) + 2*x(:,:,2),pm,1)

  call assert(maxval(abs(df(1:sz(1)-1,:) - 3)),0.,1e-10,'grad x')


 end subroutine test_grad2d


 !_______________________________________________________
 !

 subroutine test_diff
  use matoper
  implicit none
  integer, parameter :: sz(2) = [3,4]
  real :: x(sz(1)*sz(2),2), pm(sz(1)*sz(2),2)
  real :: df(sz(1)*sz(2)), df2(sz(1),sz(2))
  logical :: mask(sz(1)*sz(2))
  type(spline) :: conf
  type(SparseMatrix) :: Sx,Sy,Lap,Lap2

  integer :: i,j,l

  l = 1
  do j = 1,sz(2)
    do i = 1,sz(1)
      x(l,1) = 2.*i
      x(l,2) = 3.*j
      l = l+1
    end do
  end do

  pm(:,1) = 1./2.
  pm(:,2) = 1./3.
  mask = .true.

  !call initspline(conf,sz,mask,pm,x)
  call initspline_rectdom(conf,sz,[2.,3.],[2.*sz(1), 3.*sz(2)])
  call assert(conf%x,x,1e-10,'rectdom x')
  call assert(conf%pm,pm,1e-10,'rectdom pm')

  df = grad(sz,mask,pm,3*x(:,1) + 2*x(:,2),1)
  df2 = reshape(df,sz)    
  call assert(maxval(abs(df2(1:sz(1)-1,:) - 3)),0.,1e-10,'gradv x')

  df = grad(sz,mask,pm,3*x(:,1) + 2*x(:,2),2)
  df2 = reshape(df,sz)    
  call assert(maxval(abs(df2(:,1:sz(2)-1) - 2)),0.,1e-10,'gradv y')



  df = laplacian(conf,3*x(:,1) + 2*x(:,2))
  df2 = reshape(df,sz)    
  call assert(maxval(abs(df2(2:sz(1)-1,2:sz(2)-1))),0.,1e-10,'laplacian (1)')

  df = laplacian(conf,3*x(:,1)**2 + 2*x(:,2))
  df2 = reshape(df,sz)    
  call assert(maxval(abs(df2(2:sz(1)-1,2:sz(2)-1) - 6)),0.,1e-10,'laplacian (2)')

  Sx = grad_op(sz,mask,pm,1)
  df2 = reshape(Sx .x. (3*x(:,1) + 2*x(:,2)) ,sz)    
  call assert(maxval(abs(df2(1:sz(1)-1,:) - 3)),0.,1e-10,'gradv op x')

  Sy = grad_op(sz,mask,pm,2)
  df2 = reshape(Sy .x. (3*x(:,1) + 2*x(:,2)) ,sz)    
  call assert(maxval(abs(df2(:,1:sz(2)-1) - 2)),0.,1e-10,'gradv op x')

  Lap = laplacian_op(conf)
  df2 = reshape(Lap .x. (3*x(:,1)**2 + 2*x(:,2)) ,sz)    
  call assert(maxval(abs(df2(2:sz(1)-1,2:sz(2)-1) - 6)),0.,1e-10,'laplacian op (1)')

  Lap2 = Lap.x.Lap
  call assert(full(Lap2),matmul(full(Lap),full(Lap)),1e-10,'laplacian op (2)')

  Lap2 = transpose(Lap).x.Lap
  call assert(full(Lap2),matmul(transpose(full(Lap)),full(Lap)),1e-10,'laplacian op (3)')

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
  integer :: sz(2) 

  sz = shape(fi)

  call assert(fi,fi(sz(1):1:-1,:),tol,'check symmetry - invert x')
  call assert(fi,fi(:,sz(2):1:-1),tol,'check symmetry - invert y')
  call assert(fi,transpose(fi),tol,'check symmetry - x <-> y')

 end subroutine test_symmetry

 !_______________________________________________________
 ! 

 subroutine test_symmetry3(fi)
  use matoper
  implicit none
  real, intent(in) :: fi(:,:,:)
  integer :: sz(3) 

  sz = shape(fi)

  call assert(fi,fi(sz(1):1:-1,:,:),tol,'check symmetry - invert x')
  call assert(fi,fi(:,sz(2):1:-1,:),tol,'check symmetry - invert y')
  call assert(fi,fi(:,:,sz(3):1:-1),tol,'check symmetry - invert z')


 end subroutine test_symmetry3

 !_______________________________________________________
 ! 

 subroutine test_2d_small
  use matoper
  use ufileformat
  implicit none
  integer, parameter :: sz(2) = [3,3]
  real, parameter :: x0(2) = [-5.,-5.], x1(2) = [5.,5.]
  real, parameter :: alpha(3) = [1,2,1]
  type(spline) :: conf
  real :: fi(sz(1),sz(2))
  real, dimension(sz(1)*sz(2)) :: x, Bx, r, kernel
  integer :: i
  ! for pcg
  !  integer :: maxit = 10000, nit
  !  real :: relres
  real :: xiBx

  write(6,*) '= test_2d_small ='

  call initspline_rectdom(conf,sz,x0,x1)
  x = 0
  !x((size(x)+1)/2) = 1
  ! middle of domain
  x(sub2ind(sz,(sz+1)/2)) = 1

  iB = invB_op(conf,alpha)

  xiBx = invB(conf,alpha,x)
  call assert(xiBx,x.x.(iB.x.x),tol,'check invB')

  !  Bx = pcg(fun,x,x,tol=tol,maxit=maxit,nit=nit,relres=relres)
  Bx = symsolve(iB,x)

  call assert(iB.x.Bx,x,tol,'check inv(B)*x')
  !  write(6,*) 'x',x
  !  write(6,*) 'nit,relres',nit,relres
  !  write(6,*) 'Bx',Bx

  

  fi = reshape(Bx,sz)
  call test_symmetry(fi)

 end subroutine test_2d_small


 !_______________________________________________________
 ! 

 subroutine test_2d
  use matoper
  use ufileformat
  implicit none
  integer, parameter :: sz(2) = [191,191]
  real, parameter :: x0(2) = [-5.,-5.], x1(2) = [5.,5.]
  real, parameter :: alpha(3) = [1,2,1]
  type(spline) :: conf
  real :: fi(sz(1),sz(2))
  real, dimension(sz(1)*sz(2)) :: x, Bx, r, kernel
  integer :: i
  ! for pcg
  !  integer :: maxit = 10000, nit
  !  real :: relres
  real :: xiBx
  real :: len = 1

  write(6,*) '= test_2d ='

  call initspline_rectdom(conf,sz,x0,x1)
  x = 0
  !x((size(x)+1)/2) = 1
  ! middle of domain
  x(sub2ind(sz,(sz+1)/2)) = 1

  iB = invB_op(conf,alpha)

  xiBx = invB(conf,alpha,x)
  call assert(xiBx,x.x.(iB.x.x),tol,'check invB')

  !  Bx = pcg(fun,x,x,tol=tol,maxit=maxit,nit=nit,relres=relres)
  Bx = symsolve(iB,x)

  call assert(iB.x.Bx,x,tol,'check inv(B)*x')
  !  write(6,*) 'x',x
  !  write(6,*) 'nit,relres',nit,relres
  !  write(6,*) 'Bx',Bx

  

  fi = reshape(Bx,sz)
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

  call assert(fi(96,[(96+i,i=1,50)]), &
       [0.99507,0.98409,0.96919,0.95146,0.93163,0.91024,0.88769,0.86431, &
        0.84035,0.81603,0.79152,0.76698,0.74251,0.71822,0.69419,0.67050, &
        0.64719,0.62431,0.60191, 0.58000,0.55862,0.53778,0.51750,0.49778, &
        0.47863,0.46004,0.44203, 0.42459,0.40771,0.39139,0.37562,0.36038, &
        0.34568,0.33150,0.31782,0.30465,0.29195,0.27973, 0.26797,0.25666, &
        0.24577,0.23531,0.22526,0.21560,0.20633,0.19742,0.18887,0.18067, &
        0.17280,0.16525],0.004,'check theoretical 2d-121-kernel')


!  call assert(Bx,kernel,0.006,'check theoretical 2d-121-kernel')

  call usave('fi.dat',fi,-9999.)


 end subroutine test_2d


 !_______________________________________________________
 ! 

 subroutine test_2d_mask
  use matoper
  use ufileformat
  implicit none
  integer, parameter :: sz(2) = [191,191]
  real, parameter :: x0(2) = [-3.,-3.], x1(2) = [3.,3.]
  real, parameter :: alpha(3) = [1,2,1]
  type(spline) :: conf
  real :: x(sz(1)*sz(2)), Bx(sz(1)*sz(2)), fi(sz(1),sz(2))
  logical :: mask(sz(1),sz(2))
  ! for pcg
  !  integer :: maxit = 10000, nit
  !  real :: relres
  real :: xiBx

  write(6,*) '= test_2d_mask ='

  mask = .true.
  mask(80,:) = .false.
  
  call initspline_rectdom(conf,sz,x0,x1,mask=reshape(mask,[product(sz)]))
  x = 0
  !x((size(x)+1)/2) = 1
  ! middle of domain
  x(sub2ind(sz,(sz+1)/2)) = 1

  !  call test_symmetry(reshape(x,sz))

  fi = reshape(laplacian(conf,x),sz)
  !  write(6,*) 'fi ',fi
  call test_symmetry(fi(2:sz(1)-1,2:sz(1)-1))

  iB = invB_op(conf,alpha)

  call test_symmetry(reshape(iB.x.x,sz))

  xiBx = invB(conf,alpha,x)
  call assert(xiBx,x.x.(iB.x.x),tol,'check invB with mask')

  !  Bx = pcg(fun,x,x,tol=tol,maxit=maxit,nit=nit,relres=relres)
  Bx = symsolve(iB,x)

  call assert(iB.x.Bx,x,tol,'check inv(B)*x with mask')
  !  write(6,*) 'x',x
  !  write(6,*) 'nit,relres',nit,relres
  !  write(6,*) 'Bx',Bx


  fi = reshape(Bx,sz)

  where (.not.mask) fi = -9999.
  call usave('fi_mask.dat',fi,-9999.)


 end subroutine 

 !_______________________________________________________
 ! 

 subroutine test_3d
  use matoper
  use ufileformat
  implicit none
  !integer, parameter :: sz(3) = [4,4,4]
  !integer, parameter :: sz(3) = [101,101,48]
  integer, parameter :: sz(3) = [21,21,9]
  real, parameter :: x0(3) = [-5.,-5.,-5.], x1(3) = [5.,5.,5.]
  real, parameter :: alpha(4) = [1,3,3,1]
  !real, parameter :: alpha(3) = [1,2,1]
  type(spline) :: conf
  real :: fi(sz(1),sz(2),sz(3))
  real, dimension(sz(1)*sz(2)*sz(3)) :: x, Bx, r, kernel
  integer :: i
  ! for pcg
  !  integer :: maxit = 10000, nit
  !  real :: relres
  real :: xiBx

  write(6,*) '= test_3d ='

  call initspline_rectdom(conf,sz,x0,x1)
  x = 0
  !x((size(x)+1)/2) = 1
  ! middle of domain
  x(sub2ind(sz,(sz+1)/2)) = 1

  iB = invB_op(conf,alpha)

  xiBx = invB(conf,alpha,x)

  call assert(xiBx,x.x.(iB.x.x),tol,'check invB')

  !  Bx = pcg(fun,x,x,tol=tol,maxit=maxit,nit=nit,relres=relres)
  Bx = symsolve(iB,x)

  call assert(iB.x.Bx,x,tol,'check inv(B)*x')
  !  write(6,*) 'x',x
  !  write(6,*) 'nit,relres',nit,relres
  !  write(6,*) 'Bx',Bx

  

  fi = reshape(Bx,sz)
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
  do j = 1,conf%n
    call near(cg,conf%x(j,:),conf%x,cdist,2*len,ind,distnear,nnz)
    do k = 1,nnz
      Ry(j) = Ry(j) + locfun(distnear(k)/len) * y(ind(k))
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

 subroutine test_loc_cov
  use covariance
  use matoper
  use ufileformat
  implicit none
  integer, parameter :: sz(2) = [5,5]
!  integer, parameter :: sz(2) = [11,11]
  real, parameter :: x0(2) = [-5.,-5.], x1(2) = [5.,5.]
  real, parameter :: alpha(3) = [1,2,1]
!  type(spline) :: conf

  real :: fi(sz(1),sz(2))
  real, dimension(sz(1)*sz(2)) :: x, Bx, r, kernel, Cx, Cx2, distnear, x2
  integer, dimension(sz(1)*sz(2)) :: ind
  integer :: i, j, k, l, nnz, status
  real :: xiBx, dist
  ! for pcg
  integer :: maxit = 10000, nit
  real :: relres
  real :: lenSpline

  write(6,*) '= test_loc_cov ='
  
  call initspline_rectdom(conf,sz,x0,x1)
  len = 1.
  lenSpline = len/2.3

  x = 0
  !x((size(x)+1)/2) = 1
  ! middle of domain
  x(sub2ind(sz,(sz+1)/2)) = 1


  Cx = 0
  do j = 1,conf%n
    do i = 1,conf%n
      dist = cdist(conf%x(i,:),conf%x(j,:))
      Cx(j) = Cx(j) + locfun(dist/len) * x(i)
    end do
  end do


  cg = setupgrid(conf%x,[len/10.,len/10.])
  Cx2 = 0
  do j = 1,conf%n
    call near(cg,conf%x(j,:),conf%x,cdist,2*len,ind,distnear,nnz)
    do k = 1,nnz
      Cx2(j) = Cx2(j) + locfun(distnear(k)/len) * x(ind(k))
    end do
  end do

  call assert(Cx,Cx2,tol,'Cx')

!  call assert(Cx,Rfun(x),tol,'Cx (2)')


  iB = invB_op(conf,alpha)
  call solver%init(iB)
  x2 = solver%solve(x,status)

  Bx = pcg(Rfun,x,x,tol=1e-4,maxit=maxit,nit=nit,relres=relres)
  call assert(Rfun(Bx),x,1e-5,'check inv(B)*x without preconditioner')
  write(6,*) 'nit',nit

  Bx = pcg(Rfun,x,x,tol=1e-4,maxit=maxit,nit=nit,relres=relres,pc=PreCondfun)
  call assert(Rfun(Bx),x,1e-5,'check inv(B)*x with preconditioner')
  write(6,*) 'nit',nit


  call solver%done()

 end subroutine test_loc_cov
end program test_nondiag
