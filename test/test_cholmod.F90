! cd /home/abarth/Assim/OAK-nonDiagR &&  make test/test_cholmod && test/test_cholmod


program test_cholmod 
 use matoper
 use iso_c_binding, only: c_int, c_double, c_size_t, c_loc

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


  integer(c_int), target :: Si(nz) = [1,2,3,4]-1
  integer(c_int), target :: Sj(nz) = [1,2,3,4]-1
  real(c_double), target :: Ss(nz) = [2,2,2,2]
  integer(c_int) :: status

  real(c_double), target :: bb(n) = [1,2,3,4]
  real(c_double), target :: xx(n) = [1,2,3,4]

  status = solve_cholmod(n,nz,c_loc(Si),c_loc(Sj),c_loc(Ss),c_loc(bb),c_loc(xx))

  call assert(xx,bb/2.,tol,'simple cholmod test')

 end subroutine simple_test

 subroutine sparse_test()
  implicit none
  integer, parameter :: m = 4
  type(SparseMatrix) :: S
  real :: b(m), x(m)

  integer :: status

  S%nz = 6
  allocate(S%i(S%nz),S%j(S%nz),S%s(S%nz))
  S%m = 4
  S%n = 4
  S%i = [1 ,2 ,3 ,4 ,  1 ,  2 ]
  S%j = [1 ,2 ,3 ,4 ,  2 ,  1 ]
  S%s = [1.,1.,1.,1., 0.1, 0.1]

  b = 1
   
  x = symsolve(S,b,status)

  call assert(real(status),0.,tol,'solver status')
  call assert(x,[0.909090909090909,0.909090909090909,1.,1.],tol,'sparse type cholmod test')

  call assert(b,S.x.x,tol,'sparse type cholmod residual test')

 end subroutine sparse_test

end program test_cholmod
