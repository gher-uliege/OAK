program test_nondiag
 implicit none

 call test_2d()
 call test_grad2d()
 call test_grad2dv()
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


 function grad2dv(sz,mask,f,pm,dim) result(df)
  use matoper
  implicit none

  logical, intent(in) :: mask(:)
  real, intent(in) :: f(:),pm(:,:)
  integer, intent(in) :: dim,sz(:)
  real :: df(size(f,1))

  integer :: i,j,l1,l2,subs(size(sz))

  df = 0
  if (dim == 1) then
    do j = 1,sz(2)
      do i = 1,sz(1)-1
        l1 = sub2ind(sz,[i,j])
        l2 = sub2ind(sz,[i+1,j])
        
        if (mask(l1).and.mask(l2)) then
          df(l1) = (f(l2) - f(l1)) * (pm(l2,1) + pm(l1,1))/2.
        end if
      end do
    end do
  end if

 end function grad2dv
   
   subroutine test_grad2dv
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
    df = grad2dv(sz,mask,3*x(:,1) + 2*x(:,2),pm,1)
    df2 = reshape(df,sz)
    
    call assert(maxval(abs(df2(1:sz(1)-1,:) - 3)),0.,1e-10,'gradv x')
    

   end subroutine test_grad2dv
    


   subroutine test_2d
    implicit none

    integer, parameter :: size(2) = [10,11]
    real :: x(size(1),size(2),2), pm(size(1),size(2),2), pn(size(1),size(2),2)

   end subroutine test_2d

  end program test_nondiag
