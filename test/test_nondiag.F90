module spline

 type config
   integer :: ndim, n
   integer, allocatable :: sz(:)
   logical, allocatable :: mask(:)
   real, allocatable :: pm(:,:)
   real, allocatable :: x(:,:)

   ! computed
   real, allocatable :: spm(:,:,:)
   ! snu(:,i) "diffusion coefficient"  along dimension i including metric, staggered
   real, allocatable :: snu(:,:)
   logical, allocatable :: smask(:,:)

 end type config

contains

 !_______________________________________________________
 !


 function diff(conf,f,dim) result(df)
  use matoper
  implicit none
  type(config), intent(in) :: conf
  real, intent(in) :: f(:)
  integer, intent(in) :: dim
  real :: df(conf%n)

  integer :: l1,l2,subs1(conf%ndim),subs2(conf%ndim)

  df = 0

  do l1 = 1,conf%n
    ! get subscripts
    ! sub1 and l1 corresponds to (i,j) if dim = 1
    subs1 = ind2sub(conf%sz,l1)

    if (subs1(dim) /= conf%sz(dim)) then
      ! sub2 and l2 corresponds to (i+1,j) if dim = 1
      subs2 = subs1
      subs2(dim) = subs2(dim)+1
      l2 = sub2ind(conf%sz,subs2)

      if (conf%mask(l1).and.conf%mask(l2)) then
        df(l1) = (f(l2) - f(l1)) * (conf%pm(l2,1) + conf%pm(l1,1))/2.
      end if
    end if
  end do

 end function diff

 !_______________________________________________________
 !


 function diff_op(conf,dim) result(S)
  use matoper
  implicit none
  type(config), intent(in) :: conf
  integer, intent(in) :: dim
  type(SparseMatrix) :: S

  integer :: l1,l2,subs1(conf%ndim),subs2(conf%ndim),maxnz

  S%nz = 0
  S%m = conf%n
  S%n = conf%n
  maxnz = 2*conf%n
  allocate(S%i(maxnz),S%j(maxnz),S%s(maxnz))

  do l1 = 1,conf%n
    ! get subscripts
    ! sub1 and l1 corresponds to (i,j) if dim = 1
    subs1 = ind2sub(conf%sz,l1)

    if (subs1(dim) /= conf%sz(dim)) then
      ! sub2 and l2 corresponds to (i+1,j) if dim = 1
      subs2 = subs1
      subs2(dim) = subs2(dim)+1
      l2 = sub2ind(conf%sz,subs2)

      if (conf%mask(l1).and.conf%mask(l2)) then
        S%nz = S%nz+1
        S%i(S%nz) = l1
        S%j(S%nz) = l2
        S%s(S%nz) = (conf%pm(l2,1) + conf%pm(l1,1))/2

        S%nz = S%nz+1
        S%i(S%nz) = l1
        S%j(S%nz) = l1
        S%s(S%nz) = -(conf%pm(l2,1) + conf%pm(l1,1))/2
      end if
    end if

  end do

 end function 

 !_______________________________________________________
 !

 function laplacian(conf,f) result(Lf)
  use matoper
  implicit none
  type(config), intent(in) :: conf
  real, intent(in) :: f(:)

  real :: Lf(conf%n)


  integer :: i
  real, dimension(conf%n) :: df, ddf

  ! compte laplacian

  Lf = 0
  do i = 1,conf%ndim
    df = diff(conf,f,i)

    !where (smask(:,i)) df = snu(:,i) * df
    df = conf%snu(:,i) * df

    ddf = diff(conf,df,i)
    Lf = Lf + ddf
  end do

  ! product(pm,2) is the inverse of the volume of a grid cell
  Lf = Lf * product(conf%pm,2)

 end function laplacian

 !_______________________________________________________
 !
 ! compute 
 ! inv(B) * x
 ! where B is the background covariance matrix in DIVA

 ! function invBx(conf,x,alpha) result(iBx)
 !  use matoper
 !  implicit none
 !  type(config), intent(in) :: conf
 !  real, intent(in) :: x(:), alpha(:)
 !  real :: iBx(conf%n)

 !  integer :: i
 !  iBx = x

 !  do i = 1,conf%ndim
 !  end do


 ! end function invBx

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

 subroutine initconfig(conf,sz,mask,pm,x)
  implicit none
  type(config) :: conf
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
      if (j==i) then
        conf%snu(:,i) = conf%snu(:,i) * conf%spm(:,j,i)
      else
        conf%snu(:,i) = conf%snu(:,i) / conf%spm(:,j,i)
      end if
    end do
  end do

 end subroutine initconfig

 !_______________________________________________________
 !

 subroutine doneconfig(conf)
  implicit none
  type(config) :: conf

  deallocate(conf%sz,conf%mask,conf%pm, &
       conf%x, conf%spm, conf%snu, &
       conf%smask)
 end subroutine doneconfig

 !_______________________________________________________
 !


end module spline




!_______________________________________________________
!

program test_nondiag
 use spline
 implicit none

 call test_2d()
 call test_grad2d()
 call test_diff()
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
  write(6,*) 'mv',maxval(abs(df(1:sz(1)-1,:) - 3))


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
  type(config) :: conf
  type(SparseMatrix) :: Sx,Sy

  integer :: i,j,l

  l = 1
  do j = 1,sz(2)
    do i = 1,sz(1)
      x(l,1) = i
      x(l,2) = j
      l = l+1
    end do
  end do

  pm = 1
  mask = .true.

  call initconfig(conf,sz,mask,pm,x)

  df = diff(conf,3*x(:,1) + 2*x(:,2),1)
  df2 = reshape(df,sz)    
  call assert(maxval(abs(df2(1:sz(1)-1,:) - 3)),0.,1e-10,'gradv x')

  df = diff(conf,3*x(:,1) + 2*x(:,2),2)
  df2 = reshape(df,sz)    
  call assert(maxval(abs(df2(:,1:sz(2)-1) - 2)),0.,1e-10,'gradv y')



  df = laplacian(conf,3*x(:,1) + 2*x(:,2))
  df2 = reshape(df,sz)    
  call assert(maxval(abs(df2(1:sz(2)-2,1:sz(2)-2))),0.,1e-10,'laplacian (1)')


  df = laplacian(conf,3*x(:,1)**2 + 2*x(:,2))
  df2 = reshape(df,sz)    
  call assert(maxval(abs(df2(1:sz(1)-2,1:sz(2)-2) - 6)),0.,1e-10,'laplacian (2)')

  Sx = diff_op(conf,1)
  df2 = reshape(Sx .x. (3*x(:,1) + 2*x(:,2)) ,sz)    
  call assert(maxval(abs(df2(1:sz(1)-1,:) - 3)),0.,1e-10,'gradv op x')

  Sy = diff_op(conf,2)
  df2 = reshape(Sy .x. (3*x(:,1) + 2*x(:,2)) ,sz)    
  call assert(maxval(abs(df2(:,1:sz(2)-1) - 2)),0.,1e-10,'gradv op x')

  call doneconfig(conf)
 end subroutine test_diff



 subroutine test_2d
  implicit none

!  integer, parameter :: size(2) = [10,11]
!  real :: x(size(1),size(2),2), pm(size(1),size(2),2), pn(size(1),size(2),2)

 end subroutine test_2d

end program test_nondiag
