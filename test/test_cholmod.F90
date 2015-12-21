! cd /home/abarth/Assim/OAK-nonDiagR &&  make test/test_cholmod && test/test_cholmod


program test_cholmod 
 use matoper
 USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT, C_LOC, c_double, c_size_t

! int solve(INDEX n,INDEX nz, int* Si, int* Sj,double* Ss, double* bb, double* xx) {

 interface
  integer(c_int) function  solve_choldmod(n,nz,Si,Sj,Ss,bb,xx) bind(c,name="solve")
    use, intrinsic :: iso_c_binding, only: c_int, c_size_t, c_double, c_ptr
    implicit none
    integer(c_size_t), value :: n,nz
    type(c_ptr), value :: Si, Sj, Ss, bb, xx
   end function solve_choldmod
end interface

integer(c_size_t), parameter :: n = 4, nz = 4

type(SparseMatrix) :: S
real :: b(n)

integer(c_int), target :: Si(nz) = [1,2,3,4]-1
integer(c_int), target :: Sj(nz) = [1,2,3,4]-1
real(c_double), target :: Ss(nz) = [2,2,2,2]
integer(c_int) :: status

real(c_double), target :: bb(n) = [1,2,3,4]
real(c_double), target :: xx(n) = [1,2,3,4]

status = solve_choldmod(n,nz,c_loc(Si),c_loc(Sj),c_loc(Ss),c_loc(bb),c_loc(xx))

write(6,*) 'xx ',xx
!S  = spdiag([1.,2.,3.,4.])
!b = 2


 
end program test_cholmod
