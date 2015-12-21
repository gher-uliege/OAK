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



 subroutine sparse_test()
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

  call assert(xx,b/2.,tol,'sparse type cholmod test')

  deallocate(Si,Sj,Ss)

 end subroutine sparse_test

end program test_cholmod
