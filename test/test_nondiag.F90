! cd /home/abarth/Assim/OAK-nonDiagR &&  make test/test_nondiag && test/test_nondiag

module mod_spline

 type spline
   ! ndim number of dimensions (e.g. 2 for lon/lat)
   ! n number of grid points
   integer :: ndim, n

   ! size of the domain
   integer, allocatable :: sz(:)

   ! mask is true for "within domain"/"sea" grid points and false 
   ! for "outside the domain"/"land"
   logical, allocatable :: mask(:)

   ! inverse of local resolution
   real, allocatable :: pm(:,:)
   real, allocatable :: x(:,:)

   ! computed
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

  Lap = spzero(conf%n,conf%n)

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
 ! compute 
 ! inv(B)
 ! where B is the background covariance matrix in DIVA

 function invB_op(conf,alpha) result(iB)
  use matoper
  implicit none
  type(spline), intent(in) :: conf
  real, intent(in) :: alpha(:)
  type(SparseMatrix) :: iB,Lap,gradient

  integer :: i

  if (size(alpha) /= 3) then
    write(6,*) 'error alpha to small', alpha
  end if

  ! contrain on value
  iB = spdiag([(alpha(1),i=1,conf%n)])

  ! contrain on gradient
  do i = 1,conf%ndim
    gradient = grad_op(conf%sz,conf%mask,conf%pm,i)
    iB = iB + (alpha(2) .x. (transpose(gradient).x.gradient))
  end do

  ! contrain on laplacian
  Lap = laplacian_op(conf)
  iB = iB + (alpha(3) .x. (transpose(Lap).x.Lap))


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
  real :: xiBx, tmp(conf%n)

  integer :: i

  if (size(alpha) /= 3) then
    write(6,*) 'error alpha to small', alpha
  end if

  ! contrain on value
  xiBx = alpha(1) * sum(x**2)

  ! contrain on gradient
  do i = 1,conf%ndim
    tmp = grad(conf%sz,conf%mask,conf%pm,x,i)
    xiBx = xiBx + alpha(2) * sum(tmp**2)
  end do

  ! contrain on laplacian
  tmp = laplacian(conf,x)
  xiBx = xiBx + alpha(3) * sum(tmp**2)
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

 subroutine initspline_rectdom(conf,sz,x0,x1)
  use matoper
  implicit none
  type(spline) :: conf
  real, intent(in) :: x0(:), x1(:) ! start and end points
  integer, intent(in) :: sz(:)

  real, allocatable :: pm(:,:), x(:,:)
  logical, allocatable :: mask(:)

  integer :: i,j,n, ndim, subs(size(sz))
  ndim = size(sz)
  n = product(sz)

  allocate(mask(n),x(n,ndim),pm(n,ndim))

  mask = .true.
  do j = 1,n
    subs = ind2sub(sz,j)
    x(j,:) = x0 + (x1-x0) * (subs-1.)/ (sz-1.) 
  end do

  do i = 1,ndim
    pm(:,i) = (sz(i)-1.) / (x1(i)-x0(i)) 
  end do

  call initspline(conf,sz,mask,pm,x)

  deallocate(mask,x,pm)

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
 implicit none
 type(SparseMatrix) :: iB
 real :: tol = 1e-7

 call test_diff()
 call test_2d()
 call test_grad2d()
contains

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

 subroutine test_2d
  use matoper
  use ufileformat
  implicit none
  integer, parameter :: sz(2) = [191,191]
  real, parameter :: x0(2) = [-3.,-3.], x1(2) = [3.,3.]
  real, parameter :: alpha(3) = [1,2,1]
  type(spline) :: conf
  real :: x(sz(1)*sz(2)), Bx(sz(1)*sz(2)), fi(sz(1),sz(2))
  ! for pcg
  !  integer :: maxit = 10000, nit
  !  real :: relres
  real :: xiBx

  call initspline_rectdom(conf,sz,x0,x1)
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
  call assert(xiBx,x.x.(iB.x.x),tol,'check invB')

  !  Bx = pcg(fun,x,x,tol=tol,maxit=maxit,nit=nit,relres=relres)
  Bx = symsolve(iB,x)

  call assert(iB.x.Bx,x,tol,'check inv(B)*x')
  !  write(6,*) 'x',x
  !  write(6,*) 'nit,relres',nit,relres
  !  write(6,*) 'Bx',Bx


  fi = reshape(Bx,sz)
  call test_symmetry(fi)

  call usave('fi.dat',fi,-9999.)


 end subroutine test_2d

end program test_nondiag
