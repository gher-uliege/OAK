!
!  OAK, Ocean Assimilation Kit
!  Copyright(c) 2002-2016 Alexander Barth and Luc Vandenblucke
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!

! include the fortran preprocessor definitions
#include "ppdef.h"

!#define PROFILE

module matoper
use, intrinsic :: iso_c_binding, only: c_ptr

! sparse matrix in trippled format

type SparseMatrix
  ! nz: non-zero elements
  ! m: number of rows
  ! n: number of columns
  ! i,j: indices of non-zero elements
  ! s: value of non-zero elements
  integer          :: nz,m,n
  integer, allocatable :: i(:),j(:)
  real, allocatable    :: s(:)
end type SparseMatrix


! matrix multiplication
! A B

interface operator(.x.)
  module procedure              &
! single precision
          smat_mult_smat,       &
    ssparsemat_mult_smat,       &
          smat_mult_ssparsemat, &
          smat_mult_svec,       &
    ssparsemat_mult_svec,       &
    ssparsemat_mult_ssparsemat, &
          svec_mult_smat,       &
          svec_mult_ssparsemat, &
          svec_mult_svec,       &
         sscal_mult_ssparsemat, &
! double precision
          dmat_mult_dmat,       &
    ssparsemat_mult_dmat,       &
          dmat_mult_ssparsemat, &
          dmat_mult_dvec,       &
    ssparsemat_mult_dvec,       &
          dvec_mult_dmat,       &
          dvec_mult_ssparsemat, &
          dvec_mult_dvec,       &
          dscal_mult_ssparsemat
end interface




! A B'

interface operator(.xt.)
  module procedure &
! single precision
    smat_mult_smatT, &
    smat_mult_ssparsematT, &
    svec_mult_ssparsematT, &
    svec_mult_svecT, &
! double precision
    dmat_mult_dmatT, &
    dmat_mult_ssparsematT, &
    dvec_mult_ssparsematT, &
    dvec_mult_dvecT
end interface

! A' B

interface operator(.tx.)
  module procedure &
! single precision
    smatT_mult_smat, &
    smatT_mult_svec, &
    ssparsematT_mult_smat, &
    ssparsematT_mult_svec, &
! double precision
    dmatT_mult_dmat, &
    dmatT_mult_dvec, &
    dsparsematT_mult_dmat, &
    dsparsematT_mult_dvec
end interface

! diag(A) B

interface operator(.dx.)
  module procedure & 
! single precision
    sdiag_mult_smat, &
    sdiag_mult_svec, &
! double precision
    ddiag_mult_dmat, &
    ddiag_mult_dvec
end interface

! A diag(B)

interface operator(.xd.)
  module procedure  &
! single precision
     smat_mult_sdiag, &
! double precision
     dmat_mult_ddiag
end interface

interface operator(+)
  module procedure           &
       ssparsemat_add_ssparsemat
end interface 


interface randn
  module procedure &
       randn_scal, &
       randn_vec, &
       randn_mat
end interface randn

interface inv
  module procedure &
    sinv, &
    dinv
end interface

interface det
  module procedure &
    sdet, &
    ddet
end interface

interface sparse
  module procedure &
       sparse_mat, &
       sparse_vec
end interface sparse

interface transpose
  module procedure &
       sptranspose
end interface transpose

! gesvd obsolet singular value decomposition interface, use svd
interface gesvd
  module procedure &
    gesvd_sf90, &
    gesvd_df90
end interface
 
interface svd
  module procedure &
    svd_sf90, &
    svd_df90
end interface
  
interface ssvd
  module procedure &
    ssvd_singleprec, &
    ssvd_doubleprec
end interface
  

interface symeig
  module procedure &
    symeig_sf90, &
    symeig_df90
end interface


interface chol
  module procedure &
    schol, &
    dchol
end interface

interface sqrtm
  module procedure &
    ssqrtm, &
    dsqrtm
end interface


interface eye
  module procedure &
!    seye, &
    deye
end interface

interface trace
  module procedure &
    strace, &
    dtrace
end interface

interface diag
  module procedure &
    sdiag,     &
    ddiag,     &
    cdiag,     &
    sdiag_mat, &
    ddiag_mat
end interface

!----------------------------------------------------------
!
! Sparse matrix solver based on CHOLMOD to solve systems
! A x = b
! where A is sparse, symmetric and positive defined
!----------------------------------------------------------

#ifdef HAS_CHOLMOD
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

#endif

!----------------------------------------------------------
!


interface assert
  module procedure assert_bool
  module procedure assert_scal
  module procedure assert_int
  module procedure assert_vec
  module procedure assert_mat
  module procedure assert_array3
end interface assert


contains

!----------------------------------------------------------
!
! default precision
!
!----------------------------------------------------------


!_______________________________________________________
!
! create the identity matrix
!

 function seye(n) result(E)
  implicit none
  integer, intent(in) :: n
  real(4) :: E(n,n)
  integer :: i

  E = 0.

  do i=1,n
    E(i,i) = 1.
  end do

 end function 

!_______________________________________________________
!
! create a matrix filled with 1
!

 function ones(n,m) result(E)
  implicit none
  integer, intent(in) :: n,m
  real :: E(n,m)

  E = 1.
 end function 

!_______________________________________________________
!
! create a matrix filled with 1
!

 function zeros(n,m) result(E)
  implicit none
  integer, intent(in) :: n,m
  real :: E(n,m)

  E = 0.
 end function 


!_______________________________________________________
!
! create a matrix filled with uniformly distributed 
! random number between 0 and 1
!

 function rand(n,m) result(E)
  implicit none
  integer, intent(in) :: n,m
  real :: E(n,m)

  call random_number(E)
 end function 

  !_______________________________________________________
  !

  function randn_scal() result(E)

    ! Return a normally distributed random scalar having
    ! zero mean and variance one.
    
    implicit none

    ! Output
    real :: E               ! random vector

    ! Local variable
    real :: tmp

    ! Local variable
    integer :: i

    call random_number(E)
    do i=1,11
       call random_number(tmp)
       E = E + tmp
    end do

    E = E-6
  end function randn_scal

  !_______________________________________________________
  !

  function randn_mat(n,m) result(E)

    ! Return a matrix with normally distributed random elements having
    ! zero mean and variance one.
    ! Better would be to implement the "Box-Muller transformation".

    implicit none

    ! Inputs 
    integer, intent(in) :: n,m ! the size of the random matrix

    ! Output
    real :: E(n,m)             ! random matrix

    ! Local variable
    integer :: i
    real :: tmp(n,m)   

    ! Sum 12 normally distributed random variables and substract 6
    call random_number(E)

    do i=1,11
      call random_number(tmp)
       E = E + tmp
    end do

    E = E-6
  end function randn_mat

  !_______________________________________________________
  !

  function randn_vec(n) result(E)

    ! Return a vector with normally distributed random elements having
    ! zero mean and variance one.
    
    implicit none

    ! Input
    integer, intent(in) :: n   ! length of the random vector

    ! Output
    real :: E(n)               ! random vector

    ! Local variable
    real :: tmp(n)

    ! Local variable
    integer :: i

    call random_number(E)
    do i=1,11
       call random_number(tmp)
       E = E + tmp;
    end do

    E = E-6
  end function randn_vec

  !_______________________________________________________
  !
  ! random permutation
  ! Knuth shuffles
  ! https://en.wikipedia.org/w/index.php?title=Random_permutation&oldid=716663005

  function randperm(n) result(ind)
   implicit none
   integer, intent(in) :: n ! length of the permutation vector

   integer :: ind(n), i, j, tmp
   real :: r
   ind = [(i,i=1,n)]
   do i = 1,n-1
     call random_number(r)
     j = floor(r*(n-i))
     ! swap ind(i) and ind(i+j)
     tmp = ind(i)
     ind(i) = ind(i+j)
     ind(i+j) = tmp
   end do
  end function 
!_______________________________________________________
!

function randUnitVector(n) result(w)
implicit none
integer, intent(in) :: n
real :: w(n)

w = reshape(randn(n,1),(/ n /));
w = w/sqrt(sum(w*w));

end function

!_______________________________________________________
!
! generates n-1 orthonormal vectors perpendicular to w
!
! Hoteit et al, 2002
! https://doi.org/10.1016/S0924-7963(02)00129-X

function perpSpace(w) result(H)
implicit none
real, intent(in) :: w(:)
real :: H(size(w),size(w)-1)

real :: alpha
integer :: n,i,j

n = size(w);

H = 0
alpha = - 1/(abs(w(n))+1);

do j=1,n-1
  do i=1,n-1
    H(i,j) = alpha * w(i)*w(j);
    if (i.eq.j) then
      H(i,j) = H(i,j)+1;
    end if
  end do
end do

do j=1,n-1
    H(n,j) = alpha * (w(n)+sign(1.,w(n)))*w(j);
end do

end function

!_______________________________________________________
!


function randOrthMatrix(n) result(Omega)
implicit none
integer, intent(in) :: n
real :: Omega(n,n)

integer :: i
real :: w(n)

Omega = 0

Omega(1:1,1) = randUnitVector(1);

do i=2,n
  w(1:i) = randUnitVector(i);

  Omega(1:i,1:i-1) =  (perpSpace(w(1:i))).x.(Omega(1:i-1,1:i-1));
  Omega(1:i,i) =  w(1:i);
end do

end function

!_______________________________________________________
!
! create a diagonal complex matrix
!

 function cdiag(d) result(A)
  implicit none
  complex, intent(in) :: d(:)
  complex :: A(size(d),size(d))
  integer :: i

  A = 0.

  do i=1,size(d)
    A(i,i) = d(i)
  end do

 end function cdiag

!_______________________________________________________
!
! create a real diagonal sparse
!

 function spdiag(d) result(S)
  implicit none
  real, intent(in) :: d(:)
  integer :: i
  type(SparseMatrix) :: S

  S%m = size(d,1)
  S%n = size(d,1)
  S%nz = size(d,1)
  allocate(S%i(S%nz),S%j(S%nz),S%s(S%nz))

  do i=1,S%m
    S%i(i) = i
    S%j(i) = i
    S%s(i) = d(i)
  end do
  
 end function spdiag

!_______________________________________________________
!
! create an identidy sparse matrix
!

 function speye(n) result(S)
  implicit none
  integer, intent(in) :: n
  type(SparseMatrix) :: S
  integer :: i
  
  S = spdiag([(1.,i=1,n)])
 end function speye

!_______________________________________________________
!
! zero sparse matrix
!

 function spzeros(m,n) result(S)
  implicit none
  integer, intent(in) :: m,n
  type(SparseMatrix) :: S

  S%m = m
  S%n = n
  S%nz = 0
  allocate(S%i(S%nz),S%j(S%nz),S%s(S%nz))
 end function spzeros

!_______________________________________________________
!
! print sparse matrix
!

 subroutine spprint(S)
  implicit none
  type(SparseMatrix), intent(in) :: S
  type(SparseMatrix) :: T

  integer :: i
  T = sparse_compress(S)

  do i = 1,T%nz
    write(6,*) 'i,j,s',T%i(i),T%j(i),T%s(i)
  end do
 end subroutine spprint

!_______________________________________________________
!
! transpose a sparse matrix
!

 function sptranspose(S) result(ST)
  implicit none
  type(SparseMatrix), intent(in) :: S
  type(SparseMatrix) :: ST

  ST%m = S%n
  ST%n = S%m
  ST%nz = S%nz
  allocate(ST%i(S%nz),ST%j(S%nz),ST%s(S%nz))

  ST%i = S%j(1:S%nz)
  ST%j = S%i(1:S%nz)
  ST%s = S%s(1:S%nz)
 end function sptranspose

!_______________________________________________________
!
! S = SPARSE(X) converts a full matrix to sparse form by
!    squeezing out any zero elements.

function sparse_mat(X) result(S)
implicit none
real, intent(in)   :: X(:,:)
type(SparseMatrix) :: S

integer :: i,j,k

S%m = size(X,1)
S%n = size(X,2)
S%nz = count(X.ne.0)

allocate(S%i(S%nz),S%j(S%nz),S%s(S%nz))
k = 1

 do j=1,S%n
   do i=1,S%m
     if (X(i,j).ne.0) then
       S%i(k) = i
       S%j(k) = j
       S%s(k) = X(i,j)
       k=k+1
     end if
   end do
 end do

end function sparse_mat

!_______________________________________________________
!

function sparse_vec(i,j,s,m,n) result(B)
 implicit none
 integer, intent(in) :: i(:),j(:) ! indices of non-zero elements
 real, intent(in) :: s(:)         ! value of non-zero elements
 integer, intent(in), optional :: m,n ! size of matrix

 type(SparseMatrix) :: B

 if (present(m)) then
   B%m = m
 else
   B%m = maxval(i)
 end if

 if (present(n)) then
   B%n = n
 else
   B%n = maxval(j)
 end if

 B%nz = size(i)

 allocate(B%i(B%nz),B%j(B%nz),B%s(B%nz))
 B%i = i
 B%j = j
 B%s = s
end function sparse_vec

!_______________________________________________________
!
! S = full(X) converts a sparse matrix to a dens matrix

function full(S) result(X)
 implicit none
 type(SparseMatrix),intent(in) :: S
 real   :: X(S%m,S%n)
 integer :: k
 X = 0

 do k=1,S%nz
   X(S%i(k),S%j(k)) = X(S%i(k),S%j(k)) + S%s(k)
 end do
end function full





!_______________________________________________________
!
! operator:  .x.
!_______________________________________________________
!

function ssparsemat_mult_ssparsemat_old(A,B) result(C)
 implicit none
 type(SparseMatrix), intent(in) :: A
 type(SparseMatrix), intent(in) :: B
 type(SparseMatrix) :: C
 integer :: k,l,nz,l1,l2,llower,lupper

 if (A%n.ne.B%m) then
   write(stderr,*) 'ssparsemat_mult_ssparsemat: size not conform: A.x.B '
   write(stderr,*) 'shape(A) ',A%m,A%n
   write(stderr,*) 'shape(B) ',B%m,B%n
   ERROR_STOP
 end if

 C%m = A%m
 C%n = B%n

 ! count
!!$  nz = 0
!!$  do k=1,A%nz
!!$    do l=1,B%nz
!!$      if (A%j(k).eq.B%i(l)) nz = nz+1
!!$    end do
!!$  end do

 ! borne sup.
 nz = 10*(A%nz+B%nz)
! nz = 2714888
!  write(6,*) 'nz ',nz, A%nz,B%nz
 allocate(C%i(nz),C%j(nz),C%s(nz))
 nz=0

! if (all(B%i(2:B%nz).ge.B%i(:B%nz-1))) then
 if (.false.) then
   ! optimized version
   !write(stdout,*) 'optimised version.'

   search: do k=1,A%nz
     l1 = 1
     l2 = B%nz
     ! quick cycle if possible
     if (A%j(k).lt.B%i(l1).or.A%j(k).gt.B%i(l2)) cycle search
!     if (all(A%j(k).ne.B%i)) cycle

     ! dichotomic search for llower
     l1 = 1
     l2 = B%nz
     llower = (l1+l2)/2

     do while ( &
     ! stop criteria
       .not.(B%i(llower).eq.A%j(k).and.(llower.eq.1.or.B%i(llower-1).ne.A%j(k))))

       if (B%i(llower).ge.A%j(k)) then
         l2 = llower-1
       else
         l1 = llower+1
       end if
       llower = (l1+l2)/2

       if (l2.lt.l1) then
          cycle search
       end if
     end do

     ! dichotomic search for lupper
     l1 = 1
     l2 = B%nz
     lupper = (l1+l2)/2

     do while ( &
     ! stop criteria
       .not.(B%i(lupper).eq.A%j(k).and.(lupper.eq.B%nz.or.B%i(lupper+1).ne.A%j(k))))

       if (B%i(lupper).le.A%j(k)) then
         l1 = lupper+1
       else
         l2 = lupper-1
       end if
       lupper = (l1+l2)/2
     end do

     do l=llower,lupper
       nz = nz+1

       if (nz.gt.size(C%i)) then
         write(stderr,*) 'Error: sorry buffer to small'
         ERROR_STOP
       end if
       C%i(nz) = A%i(k)
       C%j(nz) = B%j(l)
       C%s(nz) = A%s(k)*B%s(l)
     end do

!         if (mod(k,100).eq.0) then
!           write(stderr,*) 'nz ',k,A%nz
!         end if

   end do search

 else
   ! general version: take a coffee
  ! write(stdout,*) 'Warning: unoptimised version.'
   do k=1,A%nz
     do l=1,B%nz
       if (A%j(k).eq.B%i(l)) then

         if (nz == size(C%s)) then
           write(stderr,*) 'Error: sorry buffer to small'
           ERROR_STOP
         end if

         nz = nz+1
         C%i(nz) = A%i(k)
         C%j(nz) = B%j(l)
         C%s(nz) = A%s(k)*B%s(l)


         !if (mod(nz,100).eq.0) then
         !  write(stderr,*) 'nz ',nz
         !end if
       end if
     end do
   end do
 end if
 C%nz = nz

end function 


!_______________________________________________________
!

function ssparsemat_mult_ssparsemat(A,B) result(C)
 implicit none
 type(SparseMatrix), intent(in) :: A
 type(SparseMatrix), intent(in) :: B
 type(SparseMatrix) :: C
 integer :: k,l,nz,l1,l2,llower,lupper

 integer :: Bi(B%nz),Bj(B%nz), Bindex(B%nz)
 real :: Bs(B%nz)

 if (A%n /= B%m) then
   write(stderr,*) 'ssparsemat_mult_ssparsemat: size not conform: A.x.B '
   write(stderr,*) 'shape(A) ',A%m,A%n
   write(stderr,*) 'shape(B) ',B%m,B%n
   ERROR_STOP
 end if

 C%m = A%m
 C%n = B%n

 ! count
!!$  nz = 0
!!$  do k=1,A%nz
!!$    do l=1,B%nz
!!$      if (A%j(k).eq.B%i(l)) nz = nz+1
!!$    end do
!!$  end do

 ! sort along index i of B

 Bi = B%i(1:B%nz)
 call sort(Bi,Bindex)
 Bj = B%j(Bindex(1:B%nz))
 Bs = B%s(Bindex(1:B%nz))

 ! count non-zero elements
 nz = 0
 do k=1,A%nz
   call find(Bi,A%j(k),llower,lupper)

   if (llower /= -1) then
     nz = nz + lupper-llower+1
   end if
 end do

 !write(6,*) 'nz ',nz,A%nz,B%nz,A%n,A%m
 allocate(C%i(nz),C%j(nz),C%s(nz))
 nz=0

 search: do k=1,A%nz

   call find(Bi,A%j(k),llower,lupper)

   if (llower == -1) cycle search

   do l=llower,lupper
     nz = nz+1

     if (nz > size(C%i)) then
       write(stderr,*) 'Error: sorry buffer to small'
       ERROR_STOP
     end if
     C%i(nz) = A%i(k)
     C%j(nz) = Bj(l)
     C%s(nz) = A%s(k)*Bs(l)
   end do

 end do search

 !write(6,*) 'nz2 ',nz,A%nz,B%nz,A%n,A%m

 C%nz = nz

 C = sparse_compress(C)

 contains 

  ! find all elementes equal to elem in sorted vector vec such that
  ! vec(l) == elem for all l in [llower,lupper]
  ! vec(l) /= elem for all l not in [llower,lupper]

  subroutine find(vec,elem,llower,lupper)
   implicit none
   integer, intent(in) :: vec(:)
   integer, intent(in) :: elem
   integer, intent(out) :: llower,lupper
   
   integer :: l1,l2

   l1 = 1
   l2 = size(vec)

   ! quick cycle if possible
   if (elem < vec(l1).or.elem > vec(l2)) then
     llower = -1
     lupper = -1
     return
   end if

   !     if (all(elem /= vec)) cycle

   ! dichotomic search for llower
   ! vec(l)  == elem for all l, llower <= l <= lupper

   l1 = 1
   l2 = size(vec)

   llower = (l1+l2)/2

   do while (.true.)
     ! stop criteria
     ! beware to avoid the evaluation of vec(llower-1) unless we are sure is 
     ! not 1 (Fortran operators are neither short-circuit nor eager: the 
     ! language specification allows the compiler to select the method for 
     ! optimization).

     if (vec(llower) == elem) then
       if (llower == 1) exit
       if (vec(llower-1) /= elem) exit
     end if

     if (vec(llower) >= elem) then
       l2 = llower-1
     else
       l1 = llower+1
     end if
     llower = (l1+l2)/2

     if (l2 < l1) then
       ! not found
       llower = -1
       lupper = -1
       return
     end if
   end do

   ! dichotomic search for lupper
   l1 = 1
   l2 = size(vec)
   lupper = (l1+l2)/2

   do while (.true.)
     ! stop criteria
     if (vec(lupper) == elem) then
       if (lupper == B%nz) exit
       if (vec(lupper+1) /= elem) exit
     end if

     ! change lupper
     if (vec(lupper) <= elem) then
       l1 = lupper+1
     else
       l2 = lupper-1
     end if
     lupper = (l1+l2)/2
   end do


  end subroutine find

end function ssparsemat_mult_ssparsemat


!_______________________________________________________
!

function ssparsemat_add_ssparsemat(A,B) result(C)
 implicit none
 type(SparseMatrix), intent(in) :: A
 type(SparseMatrix), intent(in) :: B
 type(SparseMatrix) :: C

 if (A%m /= B%m .or. A%n /= B%n) then
   write(stderr,*) 'ssparsemat_add_ssparsemat: size not conform: A + B '
   write(stderr,*) 'shape(A) ',A%m,A%n
   write(stderr,*) 'shape(B) ',B%m,B%n
   ERROR_STOP
 end if

 C%m = A%m
 C%n = A%n
 C%nz = A%nz + B%nz
 allocate(C%i(C%nz),C%j(C%nz),C%s(C%nz))
 ! could be optimized
 C%i = [A%i(1:A%nz), B%i(1:B%nz)]
 C%j = [A%j(1:A%nz), B%j(1:B%nz)]
 C%s = [A%s(1:A%nz), B%s(1:B%nz)]

 C = sparse_compress(C)

end function ssparsemat_add_ssparsemat

!_______________________________________________________
!

function sparse_compress(A) result(B)
 implicit none
 type(SparseMatrix), intent(in) :: A
 type(SparseMatrix) :: B

 integer :: l(A%nz), sort_index(A%nz), i, j, nz, lprevious

 ! sort index
 l = A%n * (A%i(1:A%nz)-1) + A%j(1:A%nz)-1
 
 ! sort first by i then by j
 call sort(l,sort_index)

 nz = A%nz
 B%m = A%m
 B%n = A%n
 allocate(B%i(nz),B%j(nz),B%s(nz))

 ! counter of non-zeros elements in B
 j = 0
 ! linear sort index of current element
 lprevious = -1

 do i=1,A%nz

   if (l(i) /= lprevious) then
     ! new non-zero element
     j = j + 1
     B%i(j) = A%i(sort_index(i))
     B%j(j) = A%j(sort_index(i))
     B%s(j) = A%s(sort_index(i))
     lprevious = l(i)
   else
     ! add to previous element
     B%s(j) = B%s(j) + A%s(sort_index(i))
   end if
 end do

 ! set non-zero elements
 B%nz = j

! write(6,*) 'B',B%nz,count(B%s(1:B%nz) /= 0)

end function sparse_compress
!_______________________________________________________
!
! Return the factorial of N where N is a positive integer.
 
function factorial(n) result (res) 
 implicit none
 integer, intent (in) :: n
 integer :: res
 integer :: i
 
 res = product ((/(i, i = 1, n)/))
 
end function factorial
 
!_______________________________________________________
!
! Compute the binomial coefficient

function nchoosek(n, k) result (res) 
 implicit none
 integer, intent (in) :: n
 integer, intent (in) :: k
 integer :: res
 
 res = factorial (n) / (factorial (k) * factorial (n - k))
 
end function nchoosek


!_______________________________________________________
!

  subroutine permute(indeces,x,y)
   implicit none
   integer, intent(in) :: indeces(:)
   real, intent(in) :: x(size(indeces))
   real, intent(out) :: y(size(indeces))

   integer :: i
   real :: tmp(size(indeces))

   do i=1,size(indeces)
     tmp(i) = x(indeces(i))    
   end do

   do i=1,size(indeces)
     y(i) = tmp(i)
   end do

  end subroutine permute

!--------------------------------------------

  subroutine ipermute(indeces,x,y)
   implicit none
   integer, intent(in) :: indeces(:)
   real, intent(in) :: x(size(indeces))
   real, intent(out) :: y(size(indeces))

   integer :: i
   real :: tmp(size(indeces))

   do i=1,size(indeces)
     tmp(indeces(i)) = x(i)    
   end do

   do i=1,size(indeces)
     y(i) = tmp(i)
   end do


  end subroutine ipermute



!--------------------------------------------
!
! single precision
!
!--------------------------------------------



#define REAL_TYPE real(kind=4)

! blas and lapack subroutines

#define gemv_TYPE sgemv
#define getrf_TYPE sgetrf
#define getri_TYPE sgetri
#define gesvd_TYPE sgesvd
#define spevx_TYPE sspevx
#define copy_TYPE scopy
#define lamch_TYPE slamch
#define dot_TYPE sdot
#define gemm_TYPE sgemm
#define syevx_TYPE ssyevx
#define syev_TYPE ssyev
#define spotrf_TYPE spotrf

#define diag_TYPE sdiag
#define diag_TYPE_mat sdiag_mat
#define trace_TYPE strace

! operators

#define mat_mult_mat_TYPE smat_mult_smat
#define mat_mult_vec_TYPE smat_mult_svec
#define vec_mult_mat_TYPE svec_mult_smat
#define vec_mult_vec_TYPE svec_mult_svec

#define ssparsemat_mult_vec_TYPE ssparsemat_mult_svec
#define ssparsemat_mult_mat_TYPE ssparsemat_mult_smat
#define mat_mult_ssparsemat_TYPE smat_mult_ssparsemat
#define vec_mult_ssparsemat_TYPE svec_mult_ssparsemat
#define scal_mult_ssparsemat_TYPE sscal_mult_ssparsemat

#define mat_mult_matT_TYPE smat_mult_smatT
#define mat_mult_ssparsematT_TYPE smat_mult_ssparsematT
#define vec_mult_ssparsematT_TYPE svec_mult_ssparsematT
#define vec_mult_vecT_TYPE svec_mult_svecT
#define matT_mult_mat_TYPE smatT_mult_smat
#define matT_mult_vec_TYPE smatT_mult_svec
#define ssparsematT_mult_mat_TYPE ssparsematT_mult_smat
#define ssparsematT_mult_vec_TYPE ssparsematT_mult_svec

#define diag_mult_mat_TYPE sdiag_mult_smat
#define diag_mult_vec_TYPE sdiag_mult_svec
#define mat_mult_diag_TYPE smat_mult_sdiag

! subroutines

#define inv_TYPE sinv
#define det_TYPE sdet
#define svd_TYPE svd_sf90
#define ssvd_TYPE ssvd_singleprec
#define gesvd_f90_TYPE gesvd_sf90 
#define symeig_TYPE symeig_sf90
#define chol_TYPE schol
#define sqrtm_TYPE ssqrtm

#include "matoper_inc.F90"




! undefine all macros

#undef REAL_TYPE

! blas and lapack subroutines

#undef gemv_TYPE
#undef getrf_TYPE
#undef getri_TYPE
#undef gesvd_TYPE
#undef spevx_TYPE
#undef copy_TYPE
#undef ssvd_TYPE
#undef lamch_TYPE
#undef dot_TYPE
#undef gemm_TYPE
#undef syevx_TYPE
#undef syev_TYPE
#undef spotrf_TYPE

#undef diag_TYPE
#undef diag_TYPE_mat
#undef trace_TYPE

! operators


#undef mat_mult_mat_TYPE
#undef mat_mult_vec_TYPE
#undef vec_mult_mat_TYPE
#undef vec_mult_vec_TYPE

#undef ssparsemat_mult_vec_TYPE
#undef ssparsemat_mult_mat_TYPE
#undef mat_mult_ssparsemat_TYPE
#undef vec_mult_ssparsemat_TYPE
#undef scal_mult_ssparsemat_TYPE

#undef mat_mult_matT_TYPE
#undef mat_mult_ssparsematT_TYPE
#undef vec_mult_ssparsematT_TYPE 
#undef vec_mult_vecT_TYPE
#undef matT_mult_mat_TYPE
#undef matT_mult_vec_TYPE
#undef ssparsematT_mult_mat_TYPE
#undef ssparsematT_mult_vec_TYPE

#undef diag_mult_mat_TYPE
#undef diag_mult_vec_TYPE
#undef mat_mult_diag_TYPE

! subroutines

#undef inv_TYPE
#undef det_TYPE
#undef svd_TYPE
#undef ssvd_TYPE
#undef gesvd_f90_TYPE
#undef symeig_TYPE
#undef chol_TYPE
#undef sqrtm_TYPE





!--------------------------------------------
!
! double precision
!
!--------------------------------------------


!_______________________________________________________
!
! create the identity matrix
!

 function deye(n) result(E)
  implicit none
  integer, intent(in) :: n
  real(8) :: E(n,n)
  integer :: i

  E = 0.

  do i=1,n
    E(i,i) = 1.
  end do

 end function 

! define macros for inclusion of matoper_inc.F90

#define REAL_TYPE real(kind=8)

! blas and lapack subroutines

#define gemv_TYPE DGEMV
#define getrf_TYPE dgetrf
#define getri_TYPE dgetri
#define gesvd_TYPE dgesvd
#define spevx_TYPE dspevx
#define copy_TYPE dcopy
#define lamch_TYPE dlamch
#define dot_TYPE ddot
#define gemm_TYPE dGEMM
#define syevx_TYPE dsyevx
#define syev_TYPE dsyev
#define spotrf_TYPE dpotrf

#define diag_TYPE ddiag
#define diag_TYPE_mat ddiag_mat
#define trace_TYPE dtrace

! operators

#define mat_mult_mat_TYPE dmat_mult_dmat
#define mat_mult_vec_TYPE dmat_mult_dvec
#define vec_mult_mat_TYPE dvec_mult_dmat
#define vec_mult_vec_TYPE dvec_mult_dvec

#define ssparsemat_mult_vec_TYPE ssparsemat_mult_dvec
#define ssparsemat_mult_mat_TYPE ssparsemat_mult_dmat
#define mat_mult_ssparsemat_TYPE dmat_mult_ssparsemat
#define vec_mult_ssparsemat_TYPE dvec_mult_ssparsemat
#define scal_mult_ssparsemat_TYPE dscal_mult_ssparsemat

#define mat_mult_matT_TYPE dmat_mult_dmatT
#define mat_mult_ssparsematT_TYPE dmat_mult_ssparsematT
#define vec_mult_ssparsematT_TYPE dvec_mult_ssparsematT
#define vec_mult_vecT_TYPE dvec_mult_dvecT
#define matT_mult_mat_TYPE dmatT_mult_dmat
#define matT_mult_vec_TYPE dmatT_mult_dvec
#define ssparsematT_mult_mat_TYPE dsparsematT_mult_dmat
#define ssparsematT_mult_vec_TYPE dsparsematT_mult_dvec

#define diag_mult_mat_TYPE ddiag_mult_dmat
#define diag_mult_vec_TYPE ddiag_mult_dvec
#define mat_mult_diag_TYPE dmat_mult_ddiag

! subroutines

#define inv_TYPE dinv
#define det_TYPE ddet
#define svd_TYPE svd_df90
#define ssvd_TYPE ssvd_doubleprec
#define gesvd_f90_TYPE gesvd_df90 
#define symeig_TYPE symeig_df90
#define chol_TYPE dchol
#define sqrtm_TYPE dsqrtm


#include "matoper_inc.F90"



  !_______________________________________________________
  !

  subroutine assert_message(cond,msg)
    
    ! Writes an error message if the specified condition is no true

    implicit none
    
    ! Inputs
    logical, intent(in) :: cond  ! condition to check   
    character(len=*) :: msg      ! message to print while checking
    character(len=7), parameter :: fmt = '(A40,A)'

    character(len=5) :: color_reset = '', color_ok = '', color_fail = ''

#ifdef COLOR
    color_ok    = achar(27)//'[32m'
    color_fail  = achar(27)//'[31m'
    color_reset = achar(27)//'[0m'
#endif

    if (cond) then
      write(6,fmt) msg, ': ' // color_ok   // '  OK  ' // color_reset
    else
      write(6,fmt) msg, ': ' // color_fail // ' FAIL ' // color_reset
    end if    
   end subroutine assert_message

  !_______________________________________________________
  !

  subroutine assert_bool(success,msg)
    
    ! Produce an error if the specified condition is not true

    implicit none
    
    ! Inputs
    logical, intent(in) :: success  ! condition to check   
    character(len=*)    :: msg      ! message to print while checking

    call assert_message(success,msg)
    if (.not.success) then
       ERROR_STOP
    end if
  end subroutine assert_bool

  !_______________________________________________________
  !

  subroutine assert_scal(found,expected,tol,msg)
    
    ! Produce an error if the scalar found is not the same as expected but
    ! equality comparison for numeric data uses a tolerance tol.

    implicit none
    
    ! Inputs
    real, intent(in) :: found    ! found value
    real, intent(in) :: expected ! expected value
    real, intent(in) :: tol      ! tolerance for comparision
    character(len=*) :: msg      ! message to print while checking

    ! Local variable
    real :: maxdiff
    maxdiff = abs(found - expected)

    call assert_message(maxdiff < tol,msg)

    ! maxdiff can be NaN
    if (.not.(maxdiff < tol)) then
       write(6,*) 'found ',found
       write(6,*) 'expected ',expected
       ERROR_STOP
    end if


  end subroutine assert_scal

  !_______________________________________________________
  !

  subroutine assert_int(found,expected,msg)    

    ! Produce an error if the integer found is not the same as expected

    implicit none
    
    ! Inputs
    integer, intent(in) :: found    ! found value
    integer, intent(in) :: expected ! expected value
    character(len=*)    :: msg      ! message to print while checking

    ! Local variable
    real :: maxdiff

    call assert_message(found == expected,msg)

    if (found /= expected) then
       write(6,*) 'found ',found
       write(6,*) 'expected ',expected
       ERROR_STOP
    end if
  end subroutine 

  !_______________________________________________________
  !

  subroutine assert_vec(found,expected,tol,msg)
    
    ! Produce an error if the vector found is not the same as expected but
    ! equality comparison for numeric data uses a tolerance tol.

    implicit none
    
    ! Inputs
    real, intent(in) :: found(:)    ! vector found 
    real, intent(in) :: expected(:) ! vector expected
    real, intent(in) :: tol         ! tolerance for comparision
    character(len=*) :: msg         ! message to print while checking

    ! Local variable
    real :: maxdiff
    maxdiff = maxval(abs(found - expected))

    call assert_message(maxdiff < tol,msg)

    if (.not.(maxdiff < tol)) then
       write(6,*) 'found ',found
       write(6,*) 'expected ',expected
       ERROR_STOP
    end if


  end subroutine assert_vec

  !_______________________________________________________
  !

  subroutine assert_mat(found,expected,tol,msg)
    
    ! Produce an error if the matrix found is not the same as expected but
    ! equality comparison for numeric data uses a tolerance tol.

    implicit none
    
    ! Inputs
    real, intent(in) :: found(:,:)     ! matrix found
    real, intent(in) :: expected(:,:)  ! matrix expected
    real, intent(in) :: tol            ! tolerance for comparision
    character(len=*) :: msg            ! message to print while checking

    ! Local variable
    real :: maxdiff
    maxdiff = maxval(abs(found - expected))

    call assert_message(maxdiff < tol,msg)

    if (.not.(maxdiff < tol)) then
       ! often too large to print
       !     write(6,*) 'found ',found
       !     write(6,*) 'expected ',expected
       ERROR_STOP
    end if


  end subroutine assert_mat

  !_______________________________________________________
  !

  subroutine assert_array3(found,expected,tol,msg)
    
    ! Produce an error if the matrix found is not the same as expected but
    ! equality comparison for numeric data uses a tolerance tol.

    implicit none
    
    ! Inputs
    real, intent(in) :: found(:,:,:)     ! matrix found
    real, intent(in) :: expected(:,:,:)  ! matrix expected
    real, intent(in) :: tol              ! tolerance for comparision
    character(len=*) :: msg              ! message to print while checking

    ! Local variable
    real :: maxdiff
    maxdiff = maxval(abs(found - expected))

    call assert_message(maxdiff < tol,msg)

    if (.not.(maxdiff < tol)) then
       ! often too large to print
       !     write(6,*) 'found ',found
       !     write(6,*) 'expected ',expected
       ERROR_STOP
    end if
  end subroutine assert_array3


  ! Redefine DATA_TYPE to the type needed for the subroutines sort and unique
#define DATA_TYPE integer 

  !_______________________________________________________
  !
  ! merge sort (non-recursive)
  !
  ! Sort all elements of vector A inplace.
  ! The optional argument ind is the sort index such that
  ! sortedA = A
  ! call sort(sortedA,ind)
  ! A(ind) is the same as sortedA

  subroutine sort (a,ind)
   implicit none
   DATA_TYPE, intent(inout) :: a(:)
   integer, intent(out), optional :: ind(:)

   DATA_TYPE :: b(size(a))
   integer :: inda(size(a)), indb(size(a))
   integer :: right, rend, num
   integer :: i,j,m, k, left

   num = size(a)
   inda = [(i, i=1,num)]

   ! k is the window size starting with 1, 2, 4, 8, ...
   k = 1
   do while (k < num) 
     do left=1,num-k,2*k

       ! divide array a into two subarrays a(left:right-1), a(right:rend)
       ! every subarray can be assumed already sorted (since we start
       ! bottom-up with k=1)

       right = left + k
       rend = min(right + k-1,num)

       ! merge two subarrays
       ! m: index in merged array
       m = left
       ! i,j: indices of subarrays
       i = left
       j = right

       ! put the smallest value of seach subarray in b until one subarray is 
       ! finished
       do while (i <= right-1 .and. j <= rend)
         if (a(i) <= a(j)) then         
           b(m) = a(i)
           indb(m) = inda(i)
           i = i+1
         else
           b(m) = a(j) 
           indb(m) = inda(j)
           j = j+1
         end if

         m = m+1
       end do

       ! put the remainer of first subarray in b (if any)
       do while (i <= right-1) 
         b(m)=a(i) 
         indb(m) = inda(i)
         i = i+1
         m = m+1
       end do

       ! put the remainer of second subarray in b (if any)
       do while (j <= rend)
         b(m)=a(j) 
         indb(m) = inda(j)
         j = j+1
         m = m+1
       end do

       ! copy over
       a(left:rend) = b(left:rend)
       inda(left:rend) = indb(left:rend)
     end do

     k = k*2
   end do

   if (present(ind)) ind = inda
  end subroutine sort

  !_______________________________________________________
  !
  ! quick sort (recursive)
  !
  ! Sort all elements of vector A inplace.
  ! The optional argument ind is the sort index such that
  ! sortedA = A
  ! call sort(sortedA,ind)
  ! A(ind) is the same as sortedA

  subroutine quicksort(A,ind)   
    implicit none
    
    ! Input
    DATA_TYPE, intent(inout), dimension(:) :: A               ! vector to sort
    
    ! Optional output
    integer, intent(out), dimension(size(A)), optional :: ind ! sort index

    ! Local variable
    integer :: sort_index(size(A))
    integer :: i

    sort_index = [(i, i=1,size(A))]

    call sort_(A,sort_index)
    if (present(ind)) ind = sort_index

  contains
    recursive subroutine sort_(A,ind)
      implicit none

      DATA_TYPE, intent(inout), dimension(:) :: A
      integer, intent(inout), dimension(:) :: ind
      integer :: iq

      if (size(A) > 1) then
         call sort_partition(A, ind, iq)
         call sort_(A(:iq-1),ind(:iq-1))
         call sort_(A(iq:),ind(iq:))
      end if
    end subroutine sort_

    subroutine sort_partition(A, ind, marker)
      implicit none

      DATA_TYPE, intent(inout), dimension(:) :: A
      integer, intent(inout), dimension(:) :: ind
      integer, intent(out) :: marker
      integer :: i, j, tempind
      DATA_TYPE :: temp
      DATA_TYPE :: x      ! pivot point

      x = A(1)
      i = 0
      j = size(A) + 1

      do
         j = j-1
         do
            if (A(j) <= x) exit
            j = j-1
         end do

         i = i+1
         do
            if (A(i) >= x) exit
            i = i+1
         end do


         if (i < j) then
            ! exchange A(i) and A(j)
            temp = A(i)
            A(i) = A(j)
            A(j) = temp

            tempind = ind(i)
            ind(i) = ind(j)
            ind(j) = tempind        
         else if (i == j) then
            marker = i+1
            return
         else
            marker = i
            return
         end if
      end do

    end subroutine sort_partition

   end subroutine quicksort

  !_______________________________________________________
  !

  subroutine unique(A,n,ind,ind2)
    
    ! Return all unique elements of A (inplace)
    ! n is the number of unique elements

    implicit none

    ! Input and output
    DATA_TYPE, intent(inout) :: A(:)   ! vector, output unique elements for A

    ! Output
    integer, intent(out) :: n           ! number of unique elements

    ! Optional outputs:
    ! indices such that
    ! ind is a vector of indices for all unique elements in A
    ! Anew(1:n) = Aold(ind)
    ! Aold = Anew(ind2)   
    integer, intent(out), dimension(size(A)), optional :: ind, ind2

    ! Local variables:
    integer :: unique_index(size(A)), unique_index2(size(A)), i

    unique_index = [(i, i=1,size(A))]

    call unique_(A,n,unique_index,unique_index2)
    if (present(ind)) ind = unique_index
    if (present(ind2)) ind2 = unique_index2

  contains
    subroutine unique_(A,n,ind,ind2)
      implicit none

      DATA_TYPE, intent(inout) :: A(:)
      integer, intent(out) :: n
      integer, intent(out), dimension(size(A)) :: ind
      integer, intent(out), dimension(size(A)) :: ind2

      integer :: i
      integer, dimension(size(A)) :: sort_ind

      ind2 = 1
      call sort(A,sort_ind)

      n = 1
      do i = 1,size(A)-1
         ind2(sort_ind(i)) = n

         A(n) = A(i)
         ind(n) = sort_ind(i)

         if (A(i) /= A(i+1)) then
            n = n+1
         end if
      end do

      A(n) = A(size(A))
      ind(n) = sort_ind(size(A))
      ind2(sort_ind(size(A))) = n

    end subroutine unique_
  end subroutine unique

  !_______________________________________________________
  !


 function pcg(fun,b,x0,tol,maxit,pc,nit,relres) result(x)

  interface 
    function fun(x) result(y)
     real, intent(in) :: x(:)
     real :: y(size(x))
    end function fun
  end interface

  interface 
    function pc(x) result(y)
     real, intent(in) :: x(:)
     real :: y(size(x))
    end function pc
  end interface

  optional pc
  !   class(Covar), intent(in) :: A
  real, intent(in) :: b(:)
  real :: x(size(b))
  real, intent(in), optional :: x0(:)
  real, intent(in), optional :: tol
  integer, intent(in), optional :: maxit
  integer, intent(out), optional :: nit
  ! relative residual 
  ! |fun(x) - b| / |b|
  real, intent(out), optional :: relres

  !   class(Covar) :: pc
  real :: tol_, zr_old, zr_new
  real, pointer :: alpha(:), beta(:)
  integer :: maxit_
  integer :: n, k
  real :: tol2
  real, dimension(size(x)) :: Ap, p, r, r_old, z
# ifdef PROFILE
  real(8) :: cputime(2)
# endif

  n = size(b)
  ! default parameters
  maxit_ = min(n,100)
  tol_ = 1e-6

  if (present(tol)) tol_ = tol
  if (present(maxit)) maxit_ = maxit


  allocate(alpha(maxit_+1),beta(maxit_+1))

  ! initial guess
  if (present(x0)) then
    x = x0
  else
    ! random initial vector
    x = reshape(randn(n,1),(/ n /))
  end if

  tol2 = tol_**2

  ! gradient at initial guess
  r = b - fun(x)

  ! quick exit
  if (sum(r**2) < tol2) then
    if (present(nit)) nit = 0
    if (present(relres)) relres = sqrt(sum(r**2)/sum(b**2))
    return
  endif


  ! apply preconditioner

  if (present(pc)) then
    z = pc(r)
  else
    z = r
  endif


  ! first search direction == gradient
  p = z

  ! compute: r' * inv(M) * z (we will need this product at several
  ! occasions)

  zr_old = sum(r*z)

  ! r_old: residual at previous iteration
  r_old = r

  do k=1,maxit_    
    !     write(6,*) ' k',k,sum(r*r),maxit_,tol2

#   ifdef PROFILE
    call cpu_time(cputime(1))
#   endif    

    ! compute A*p
    Ap = fun(p)

#   ifdef PROFILE
    call cpu_time(cputime(2))
    write(stdout,*) 'pcg fun',cputime(2)-cputime(1)
#   endif

    !maxdiff(A*p,Ap)

    ! how far do we need to go in direction p?
    ! alpha is determined by linesearch

    ! alpha z'*r / (p' * A * p)
    alpha(k) = zr_old / ( sum(p * Ap))

    ! get new estimate of x
    x = x + alpha(k)*p

    ! recompute gradient at new x. Could be done by
    ! r = b-fun(x)
    ! but this does require an new call to fun
    r = r - alpha(k)*Ap

    ! apply pre-conditionner
    if (present(pc)) then
      z = pc(r)
    else
      z = r
    endif


    zr_new = sum(r*z)

    if (sum(r*r) < tol2) then
      if (present(nit)) nit = k
      exit
    endif

    !Fletcher-Reeves
    beta(k+1) = zr_new / zr_old
    !Polak-Ribiere
    !beta(k+1) = r'*(r-r_old) / zr_old
    !Hestenes-Stiefel
    !beta(k+1) = r'*(r-r_old) / (p'*(r-r_old))
    !beta(k+1) = r'*(r-r_old) / (r_old'*r_old)


    ! norm(p)
    p = z + beta(k+1)*p
    zr_old = zr_new
    r_old = r
  enddo

  if (present(relres)) relres = sqrt(sum((fun(x) - b)**2)/sum(b**2))


 end function pcg

  !_______________________________________________________
  !
  ! solves the system S x = b
  ! A is a sparse symetric positive defined matrix


#ifdef HAS_CHOLMOD

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


  !_______________________________________________________
  !

 function symsolve(S,b,status) result(x)
  implicit none
  type(SparseMatrix), intent(in) :: S
  integer, intent(out), optional :: status ! negative if error, otherwise 0
  real, intent(in) :: b(:)
  real :: x(size(b))

  integer :: stat
  type(SparseSolver) :: solver

  call solver%init(S)
  x = solver%solve(b,stat)
  call solver%done()


  if (present(status)) status = stat
 end function symsolve

#endif

  !_______________________________________________________
  !

  function ind2sub(dims,ind) result(sub)
    
    ! Create subscripts sub(1), sub(2), ... from a linear index ind.
    ! The linear index is interpreted in column-major order. All indices are 
    ! 1-based.
    ! https://en.wikipedia.org/wiki/Row-major_order
    ! This function is the same as ind2sub in matlab/octave.

    implicit none

    ! Input
    integer, intent(in) :: dims(:)  ! size of an array
    integer, intent(in) :: ind      ! global linear index

    ! Output
    integer :: sub(size(dims))      ! subscript

    ! Local variables
    integer :: j,tmp,offset(size(dims))

    offset(1) = 1
    do j = 2,size(dims)
       offset(j) = offset(j-1) * dims(j-1)
    end do

    ! tmp is here 0-based
    tmp = ind - 1

    do j = size(dims), 1, -1
       ! integer division
       sub(j) = tmp / offset(j)
       tmp = tmp - sub(j) * offset(j)
    end do

    ! make sub 1-based
    sub = sub+1
  end function ind2sub

  !_______________________________________________________
  !

  function sub2ind(dims,sub) result(ind)
    
    ! Create a linear index ind from subscripts sub(1), sub(2), ....
    ! The linear index is interpreted in column-major order. All indices are 1-based.
    ! https://en.wikipedia.org/wiki/Row-major_order
    ! This function is the same as ind2sub in matlab/octave.

    implicit none

    ! Input
    integer, intent(in) :: dims(:)  ! size of an array
    integer, intent(in) :: sub(size(dims))      ! subscript

    ! Output
    integer             :: ind      ! global linear index

    ! Local variables
    integer :: j,offset(size(dims))

    ind = 0
    offset(1) = 1
    do j = 2,size(dims)
       offset(j) = offset(j-1) * dims(j-1)
    end do

    ind = sum((sub-1) * offset)+1
  end function 



!!$
!!$!_______________________________________________________
!!$!
!!$
!!$subroutine gesvd_sf90(JOBU,JOBVT,A,S,U,VT,INFO)
!!$implicit none
!!$character, intent(in) :: jobu,jobvt
!!$real(kind=4), intent(in) :: A(:,:)
!!$real(kind=4), intent(out) :: S(:), U(:,:), VT(:,:)
!!$integer, intent(out) :: info
!!$
!!$real(kind=4) :: dummy,rlwork
!!$integer :: lwork
!!$real(kind=4), allocatable :: work(:)
!!$
!!$#ifndef ALLOCATE_LOCAL_VARS
!!$  real(kind=4) :: copyA(size(A,1),size(A,2))
!!$#else
!!$  real(kind=4), pointer :: copyA(:,:)
!!$  allocate(copyA(size(A,1),size(A,2)))
!!$#endif
!!$
!!$
!!$copyA = A
!!$
!!$call SGESVD( JOBU, JOBVT,size(A,1),size(A,2), copyA, size(A,1), &
!!$   S, U,size(U,1), VT,size(VT,1), &
!!$     rlWORK, -1, INFO)
!!$lwork = rlwork+0.5
!!$allocate(work(lwork))
!!$
!!$call SGESVD( JOBU, JOBVT,size(A,1),size(A,2), copyA, size(A,1), &
!!$   S, U,size(U,1), VT,size(VT,1), &
!!$     WORK, lwork, INFO)
!!$deallocate(work)
!!$#ifdef ALLOCATE_LOCAL_VARS
!!$  deallocate(copyA)
!!$#endif
!!$
!!$end subroutine 
!!$
!!$!_______________________________________________________
!!$!
!!$
!!$subroutine gesvd_df90(JOBU,JOBVT,A,S,U,VT,INFO)
!!$implicit none
!!$character, intent(in) :: jobu,jobvt
!!$real(kind=8), intent(in) :: A(:,:)
!!$real(kind=8), intent(out) :: S(:), U(:,:), VT(:,:)
!!$integer, intent(out) :: info
!!$
!!$real(kind=8) :: dummy,rlwork
!!$integer :: lwork
!!$real(kind=8), allocatable :: work(:)
!!$
!!$#ifndef ALLOCATE_LOCAL_VARS
!!$  real(kind=4) :: copyA(size(A,1),size(A,2))
!!$#else
!!$  real(kind=4), pointer :: copyA(:,:)
!!$  allocate(copyA(size(A,1),size(A,2)))
!!$#endif
!!$
!!$copyA = A
!!$
!!$call DGESVD( JOBU, JOBVT,size(A,1),size(A,2), copyA, size(A,1), &
!!$   S, U,size(U,1), VT,size(VT,1), &
!!$     rlWORK, -1, INFO)
!!$
!!$lwork = rlwork+0.5
!!$allocate(work(lwork))
!!$
!!$call DGESVD( JOBU, JOBVT,size(A,1),size(A,2), copyA, size(A,1), &
!!$   S, U,size(U,1), VT,size(VT,1), &
!!$     WORK, lwork, INFO)
!!$deallocate(work)
!!$
!!$#ifdef ALLOCATE_LOCAL_VARS
!!$  deallocate(copyA)
!!$#endif
!!$
!!$
!!$end subroutine 
!!$
!!$
!!$
!!$!_______________________________________________________
!!$!
!!$! computes eigenvalue/-vector of a symetric matrix
!!$!
!!$! A = V' diag(E) V 
!!$
!!$subroutine symeig_sf90(A,E,V,nbiggest,nsmallest,indices,INFO)
!!$implicit none
!!$real(kind=4), intent(in)  :: A(:,:)
!!$real(4),         intent(out) :: E(size(A,1))
!!$integer, optional, intent(in) :: nbiggest, nsmallest,indices(2)
!!$real(4), optional, target, intent(out) :: V(:,:)
!!$integer, optional, intent(out) :: info
!!$
!!$character :: jobz
!!$real(kind=4), pointer :: pV(:,:)
!!$real(kind=4) :: rlwork, tmp
!!$integer :: lwork, myinfo, N, iwork(5*size(A,1)), ifail(size(A,1)), i,j
!!$real(kind=4), allocatable :: work(:)
!!$
!!$! LAPACK Machine precision routine
!!$
!!$real :: slamch
!!$integer :: r,idummy,ind(2)
!!$
!!$
!!$#ifndef ALLOCATE_LOCAL_VARS
!!$  real(kind=4) :: copyA(size(A,1),size(A,2))
!!$#else
!!$  real(kind=4), pointer :: copyA(:,:)
!!$  allocate(copyA(size(A,1),size(A,2)))
!!$#endif
!!$
!!$
!!$
!!$jobz='n'
!!$N = size(A,1)
!!$r = n
!!$
!!$if (present(V)) then
!!$  jobz='v'
!!$  pV => V
!!$else
!!$  allocate(pV(1,1))
!!$end if
!!$
!!$ind = (/ 1,n /)
!!$
!!$if (present(nbiggest))  ind = (/ n-nbiggest+1,n /)
!!$if (present(nsmallest)) ind = (/ 1,nsmallest /)
!!$if (present(indices))   ind = indices
!!$
!!$! protect content of A
!!$
!!$copyA = A
!!$
!!$! determine the optimal size of work
!!$
!!$call SSYEVX(JOBZ,'I','U',n,copyA,n,-1.,-1.,ind(1),ind(2),     &
!!$     2*SLAMCH('S'),idummy,E,pV,n,rlWORK,-1, IWORK,   &
!!$     IFAIL, myINFO )
!!$
!!$lwork = rlwork+0.5
!!$allocate(work(lwork))
!!$
!!$call SSYEVX(JOBZ,'I','U',n,copyA,n,-1.,-1.,ind(1),ind(2),     &
!!$     2*SLAMCH('S'),idummy,E,pV,n, WORK, LWORK, IWORK,   &
!!$     IFAIL, myINFO )
!!$
!!$  if (present(nbiggest)) then
!!$! sort in descending order
!!$    do i=1,nbiggest/2
!!$      tmp = E(i)
!!$      E(i) = E(nbiggest-i+1)
!!$      E(nbiggest-i+1) = tmp
!!$
!!$      if (present(V)) then
!!$       do j=1,n
!!$        tmp = V(j,i) 
!!$        V(j,i) = V(j,nbiggest-i+1)
!!$        V(j,nbiggest-i+1) = tmp
!!$       end do
!!$      end if
!!$    end do
!!$  end if
!!$
!!$deallocate(work)
!!$#ifdef ALLOCATE_LOCAL_VARS
!!$  deallocate(copyA)
!!$#endif
!!$
!!$
!!$if (.not.present(V)) deallocate(pV)
!!$if (present(info)) info = myinfo
!!$end subroutine
!!$
!!$!_______________________________________________________
!!$!


end module
