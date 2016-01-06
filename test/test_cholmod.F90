! cd /home/abarth/Assim/OAK-nonDiagR &&  make test/test_cholmod && test/test_cholmod


module cholmod
 use, intrinsic :: iso_c_binding, only: c_ptr
 use matoper

  type SparseSolver
    type(c_ptr), pointer :: cholmod_common, A, L
    contains
     procedure :: init => SparseSolver_init
     procedure :: solve => SparseSolver_solve
     procedure :: done => SparseSolver_done

  end type SparseSolver

  interface
    integer(c_int) function cholmod_start(cholmod_common) bind(c, name='cholmod_wrapper_start')
    use, intrinsic :: iso_c_binding, only: c_int, c_ptr
    implicit none
    type (c_ptr), value :: cholmod_common
   end function 

    integer(c_int) function cholmod_finish(cholmod_common) bind(c, name='cholmod_wrapper_finish')
    use, intrinsic :: iso_c_binding, only: c_int, c_ptr
    implicit none
    type (c_ptr), value :: cholmod_common
   end function 


   integer(c_int) function cholmod_matrix(cholmod_common,n,nz,Si,Sj,Ss,A) bind(c,name="cholmod_wrapper_matrix")
    use, intrinsic :: iso_c_binding, only: c_int, c_size_t, c_ptr
    implicit none
    type (c_ptr), value :: cholmod_common
    integer(c_size_t), value :: n,nz
    type(c_ptr), value :: Si, Sj, Ss
    type (c_ptr), value :: A ! cholmod_sparse
   end function cholmod_matrix

   integer(c_int) function cholmod_factorize(cholmod_common,A,L) bind(c,name="cholmod_wrapper_factorize")
    use, intrinsic :: iso_c_binding, only: c_int, c_ptr
    implicit none
    type (c_ptr), value :: cholmod_common
    type (c_ptr), value :: A ! cholmod_sparse
    type (c_ptr), value :: L ! cholmod_factor
   end function cholmod_factorize


   integer(c_int) function cholmod_solve(cholmod_common,A,L,b,x) bind(c,name="cholmod_wrapper_solve")
    use, intrinsic :: iso_c_binding, only: c_int, c_size_t, c_double, c_ptr
    implicit none
    type(c_ptr), value :: cholmod_common
    type(c_ptr), value :: A ! cholmod_sparse
    type(c_ptr), value :: L ! cholmod_factor
    type(c_ptr), value :: b, x
   end function cholmod_solve

   integer(c_int) function cholmod_free(cholmod_common,A,L) bind(c,name="cholmod_wrapper_free")
    use, intrinsic :: iso_c_binding, only: c_int, c_ptr
    implicit none
    type (c_ptr), value :: cholmod_common
    type (c_ptr), value :: A ! cholmod_sparse
    type (c_ptr), value :: L ! cholmod_factor
   end function cholmod_free

  end interface


  contains

 subroutine SparseSolver_init(this,S) 
  use, intrinsic :: iso_c_binding, only: c_int, c_size_t, c_double, c_loc
  implicit none
  class(SparseSolver) :: this
  type(SparseMatrix) :: S

  ! variables for c wrapper
  integer(c_size_t) :: n, nz
  integer(c_int), target, allocatable :: Si(:), Sj(:)
  real(c_double), target, allocatable :: Ss(:)
  integer(c_int) :: status

  integer :: l

  nz = 0
  allocate(this%cholmod_common,this%A,this%L)
  allocate(Si(S%nz),Sj(S%nz),Ss(S%nz))

  do l = 1,size(Ss)   
    ! only upper part
    if (S%j(l) >= S%i(l)) then
      nz = nz+1
      Si(nz) = S%i(l)-1
      Sj(nz) = S%j(l)-1
      Ss(nz) = S%s(l)
    end if
  end do


  n = S%m

  status = cholmod_start(c_loc(this%cholmod_common))

  status = cholmod_matrix(c_loc(this%cholmod_common),n,nz,c_loc(Si),c_loc(Sj), &
       c_loc(Ss), c_loc(this%A))

  status = cholmod_factorize(c_loc(this%cholmod_common),c_loc(this%A),c_loc(this%L))

  
 end subroutine SparseSolver_init

 function SparseSolver_solve(this,b,stat) result(x)
  use, intrinsic :: iso_c_binding, only: c_int, c_double, c_loc
  implicit none
  class(SparseSolver) :: this
  real, intent(in) :: b(:)
  integer, intent(out), optional :: stat
  real :: x(size(b))

  real(c_double), target, allocatable :: bb(:), xx(:)
  integer(c_int) :: status

  allocate(bb(size(b)),xx(size(b)))
  bb = b

  status = cholmod_solve(c_loc(this%cholmod_common),c_loc(this%A), &
       c_loc(this%L),c_loc(bb), & 
       c_loc(xx))

  x = xx
  if (present(stat)) stat = status
 end function SparseSolver_solve


 subroutine SparseSolver_done(this) 
  use, intrinsic :: iso_c_binding, only: c_int, c_loc
  implicit none
  class(SparseSolver) :: this
  integer(c_int) :: status

  status = cholmod_free(c_loc(this%cholmod_common),c_loc(this%A),c_loc(this%L))
  status = cholmod_finish(c_loc(this%cholmod_common))

  deallocate(this%cholmod_common,this%A,this%L)
 end subroutine SparseSolver_done
   
end module cholmod

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
 call factorize_test()
 call sparse_solver_test()
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


 subroutine factorize_test()
  use cholmod
  use iso_c_binding, only: c_int, c_double, c_size_t, c_loc
  implicit none
  

  type(c_ptr), target :: cholmod_common, A, L
  integer(c_int) :: status
  integer(c_size_t), parameter :: n = 4, nz = 4


  integer(c_int), target :: Si(nz) = [1,2,3,4]-1
  integer(c_int), target :: Sj(nz) = [1,2,3,4]-1
  real(c_double), target :: Ss(nz) = [2,2,2,2]

  real(c_double), target :: bb(n) = [1,2,3,4]
  real(c_double), target :: b2(n) = [4,5,2,-4]
  real(c_double), target :: xx(n) = [1,2,3,4]

  
  status = cholmod_start(c_loc(cholmod_common))
!  write(6,*) 'cholmod_common ',cholmod_common

  status = cholmod_matrix(c_loc(cholmod_common),n,nz,c_loc(Si),c_loc(Sj), &
       c_loc(Ss), c_loc(A))

  status = cholmod_factorize(c_loc(cholmod_common),c_loc(A),c_loc(L))


  status = cholmod_solve(c_loc(cholmod_common),c_loc(A),c_loc(L),c_loc(bb), & 
       c_loc(xx))

  call assert(xx,bb/2.,tol,'factorize cholmod test (1)')

  status = cholmod_solve(c_loc(cholmod_common),c_loc(A),c_loc(L),c_loc(b2), & 
       c_loc(xx))

  call assert(xx,b2/2.,tol,'factorize cholmod test (2)')


  status = cholmod_free(c_loc(cholmod_common),c_loc(A),c_loc(L))

  status = cholmod_finish(c_loc(cholmod_common))

 end subroutine factorize_test



 subroutine sparse_solver_test()
  use cholmod
  implicit none
  integer, parameter :: m = 4
  type(SparseMatrix) :: S
  type(SparseSolver) :: solver
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

  call solver%init(S)
  x = solver%solve(b,status)
  call solver%done()

!  call assert(real(status),0.,tol,'solver status')
  call assert(x,[0.909090909090909,0.909090909090909,1.,1.],tol,'sparse type cholmod test')

  call assert(b,S.x.x,tol,'sparse type cholmod residual test')

 end subroutine sparse_solver_test


end program test_cholmod
