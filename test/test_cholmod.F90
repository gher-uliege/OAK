! cd /home/abarth/Assim/OAK-nonDiagR &&  make test/test_cholmod && test/test_cholmod


program test_cholmod 
 use matoper
 use iso_c_binding, only: c_int, c_float, c_loc, c_double, c_size_t

 ! int solve(INDEX n,INDEX nz, int* Si, int* Sj,double* Ss, double* bb, double* xx) {

 interface
   integer(c_int) function  solve_cholmod(n,nz,Si,Sj,Ss,bb,xx) bind(c,name="solve")
    use, intrinsic :: iso_c_binding, only: c_int, c_size_t, c_double, c_ptr
    implicit none
    integer(c_size_t), value :: n,nz
    type(c_ptr), value :: Si, Sj, Ss, bb, xx
   end function solve_cholmod
 end interface

 real, parameter :: tol = 1e-7


 call simple_test()
 call sparse_diag_test()
 call sparse_test()

contains

 subroutine simple_test()
  implicit none
  integer(c_size_t), parameter :: n = 4, nz = 4

  type(SparseMatrix) :: S
  real :: b(n)

  integer(c_int), target :: Si(nz) = [1,2,3,4]-1
  integer(c_int), target :: Sj(nz) = [1,2,3,4]-1
  real(c_double), target :: Ss(nz) = [2,2,2,2]
  integer(c_int) :: status

  real(c_double), target :: bb(n) = [1,2,3,4]
  real(c_double), target :: xx(n) = [1,2,3,4]

  status = solve_cholmod(n,nz,c_loc(Si),c_loc(Sj),c_loc(Ss),c_loc(bb),c_loc(xx))

  call assert(xx,bb/2.,tol,'simple cholmod test')

  !S  = spdiag([1.,2.,3.,4.])
  !b = 2

 end subroutine simple_test



 subroutine sparse_diag_test()
  implicit none
  integer, parameter :: m = 4
  type(SparseMatrix) :: S
  real :: b(m)

  integer(c_size_t) :: n, nz
  integer(c_int), target, allocatable :: Si(:), Sj(:)
  real(c_double), target, allocatable :: Ss(:), bb(:), xx(:)
  integer(c_int) :: status
  

  S  = spdiag([2.,2.,2.,2.])
  b = [1,2,3,4]

  allocate(Si(S%nz),Sj(S%nz),Ss(S%nz),bb(m),xx(m))
  Si = S%i-1
  Sj = S%j-1
  Ss = S%s
  bb = b

  n = S%m
  nz = S%nz

  status = solve_cholmod(n,nz,c_loc(Si),c_loc(Sj),c_loc(Ss),c_loc(bb),c_loc(xx))

  call assert(xx,b/2.,tol,'sparse type diag cholmod test')

  deallocate(Si,Sj,Ss)

 end subroutine sparse_diag_test


 subroutine sparse_test()
  implicit none
  integer, parameter :: m = 4
  type(SparseMatrix) :: S
  real :: b(m)

  integer(c_size_t) :: n, nz
  integer(c_int), target, allocatable :: Si(:), Sj(:)
  real(c_double), target, allocatable :: Ss(:), bb(:), xx(:)
  integer(c_int) :: status

  integer :: l

  !S  = spdiag([1.,1.,1.,1.])
  S%nz = 6
  allocate(S%i(S%nz),S%j(S%nz),S%s(S%nz))
  S%m = 4
  S%n = 4
  S%i = [1.,2.,3.,4.,  1.,  2.]
  S%j = [1.,2.,3.,4.,  2.,  1.]
  S%s = [1.,1.,1.,1., 0.1, 0.1]

  b = 1

  nz = 0
  allocate(Si(S%nz),Sj(S%nz),Ss(S%nz))
  allocate(bb(m),xx(m))
!  allocate(Si(5),Sj(5),Ss(5))
!  allocate(Si(4),Sj(4),Ss(4))

  do l = 1,size(Ss)   
!  do l = 1,4
    ! only upper part
    if (S%j(l) >= S%i(l)) then
      nz = nz+1
      Si(nz) = S%i(l)-1
      Sj(nz) = S%j(l)-1
      Ss(nz) = S%s(l)
    end if
  end do

  write(6,*) 'nz',nz
  bb = b

  n = S%m

  status = solve_cholmod(n,nz,c_loc(Si),c_loc(Sj),c_loc(Ss),c_loc(bb),c_loc(xx))

  write(6,*) 'xx',xx

  call assert(xx,[0.909090909090909,0.909090909090909,1.,1.],tol,'sparse type cholmod test')

  deallocate(Si,Sj,Ss)

 end subroutine sparse_test

end program test_cholmod
