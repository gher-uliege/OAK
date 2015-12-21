! cd /home/abarth/Assim/OAK-nonDiagR &&  make test/test_cholmod && test/test_cholmod


program test_cholmod 
 use matoper


 real, parameter :: tol = 1e-7

 call sparse_test()

contains


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
