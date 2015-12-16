program test_nondiag
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

    call assert_scal(maxval(abs(df(1:sz(1)-1,:) - 3)),0.,1e-10,'grad x')
    write(6,*) 'mv',maxval(abs(df(1:sz(1)-1,:) - 3))
    

   end subroutine test_grad2d

  !_______________________________________________________
  !


 function diff(sz,mask,f,pm,dim) result(df)
  use matoper
  implicit none

  logical, intent(in) :: mask(:)
  real, intent(in) :: f(:),pm(:,:)
  integer, intent(in) :: dim,sz(:)
  real :: df(size(f,1))

  integer :: i,j,l1,l2,subs1(size(sz)),subs2(size(sz))

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
        df(l1) = (f(l2) - f(l1)) * (pm(l2,1) + pm(l1,1))/2.
      end if
    end if
  end do

 end function diff

  !_______________________________________________________
  !

 function laplacian(sz,mask,f,pm) result(Lf)
  use matoper
  implicit none

  logical, intent(in) :: mask(:)
  real, intent(in) :: f(:),pm(:,:)
  integer, intent(in) :: sz(:)
  real :: Lf(size(f,1))
  logical :: smask(size(mask),size(sz))


  integer :: i,j,l1,l2,subs1(size(sz)),subs2(size(sz)), n
  real, dimension(size(mask)) :: df, ddf
  ! staggered version of pm (# grid point,dimension,staggering dimension)
  real, dimension(size(mask),size(sz),size(sz)) :: spm
  ! snu(:,i) "diffusion coefficient"  along dimension i including metric, staggered
  real, dimension(size(mask),size(sz)) :: snu

  ! number of dimensions
  n = size(sz)
  ! precompute
  snu = 1

  do i = 1,n
    smask(:,i) = stagger_mask(sz,mask,i)
  
    do j = 1,n
      spm(:,j,i) = stagger(sz,mask,pm(:,j),i)
    end do

    do j = 1,n
      if (j==i) then
        snu(:,i) = snu(:,i) * spm(:,j,i)
      else
        snu(:,i) = snu(:,i) / spm(:,j,i)
      end if
    end do
  end do

  Lf = 0
  do i = 1,n
    df = diff(sz,mask,f,pm,i)

    !where (smask(:,i)) df = snu(:,i) * df
    df = snu(:,i) * df

    ddf = diff(sz,smask(:,i),df,pm,i)
    Lf = Lf + ddf
  end do

  ! product(pm,2) is the inverse of the volume of a grid cell
  Lf = Lf * product(pm,2)

 end function laplacian

  !_______________________________________________________
  !



 function stagger_mask(sz,mask,dim) result(smask)
  use matoper
  implicit none

  integer, intent(in) :: dim,sz(:)
  logical, intent(in) :: mask(:)
  logical :: smask(size(mask))

  integer :: i,j,l1,l2,subs1(size(sz)),subs2(size(sz))

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

  integer :: i,j,l1,l2,subs1(size(sz)),subs2(size(sz))

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

 end function 

  !_______________________________________________________
  !
   
   subroutine test_diff
    use matoper
    integer, parameter :: sz(2) = [3,4]
    real :: x(sz(1)*sz(2),2), pm(sz(1)*sz(2),2)
    real :: df(sz(1)*sz(2)), df2(sz(1),sz(2))
    logical :: mask(sz(1)*sz(2))
    
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
    df = diff(sz,mask,3*x(:,1) + 2*x(:,2),pm,1)
    df2 = reshape(df,sz)    
    call assert(maxval(abs(df2(1:sz(1)-1,:) - 3)),0.,1e-10,'gradv x')

    df = diff(sz,mask,3*x(:,1) + 2*x(:,2),pm,2)
    df2 = reshape(df,sz)    
    call assert(maxval(abs(df2(:,1:sz(2)-1) - 2)),0.,1e-10,'gradv y')



    df = laplacian(sz,mask,3*x(:,1) + 2*x(:,2),pm)
    df2 = reshape(df,sz)    
    call assert(maxval(abs(df2(1:sz(2)-2,1:sz(2)-2))),0.,1e-10,'laplacian (1)')


    df = laplacian(sz,mask,3*x(:,1)**2 + 2*x(:,2),pm)
    df2 = reshape(df,sz)    
    call assert(maxval(abs(df2(1:sz(1)-2,1:sz(2)-2) - 6)),0.,1e-10,'laplacian (2)')
    

   end subroutine test_diff
    


   subroutine test_2d
    implicit none

    integer, parameter :: size(2) = [10,11]
    real :: x(size(1),size(2),2), pm(size(1),size(2),2), pn(size(1),size(2),2)

   end subroutine test_2d

  end program test_nondiag
