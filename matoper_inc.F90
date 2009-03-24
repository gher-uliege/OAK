! Copyright(c) 2002-2009 Alexander Barth and Luc Vandenblucke


!_______________________________________________________
!
! create a diagonal matrix
!

 function diag_TYPE(d) result(A)
  implicit none
  REAL_TYPE, intent(in) :: d(:)
  REAL_TYPE :: A(size(d),size(d))
  integer :: i

  A = 0.

  do i=1,size(d)
    A(i,i) = d(i)
  end do

 end function 

!_______________________________________________________
!
! computes the trace of a matrix
!

 function trace_TYPE(A) result(t)
  implicit none
  REAL_TYPE, intent(in) :: A(:,:)
  REAL_TYPE :: t
  integer :: i

  t = 0.

  do i=1,min(size(A,1),size(A,2))
    t = t + A(i,i)
  end do

 end function 



!_______________________________________________________
!
! operator:  .x.
!_______________________________________________________
!


function mat_mult_mat_TYPE(A,B) result(C)
implicit none
REAL_TYPE, intent(in) :: A(:,:), B(:,:)
REAL_TYPE :: C(size(A,1),size(B,2))

if (size(A,2).ne.size(B,1)) then
  write(stderr,*) 'smat_mult_smat: size not conform: A.x.B '
  write(stderr,*) 'shape(A) ',shape(A)
  write(stderr,*) 'shape(B) ',shape(B)
  call abort()
  stop
end if


call gemm_TYPE ('n','n',size(A,1),size(B,2),size(A,2),1._4, &
  A, size(A,1), B, size(B,1), &
  0._4, C, size(A,1))

end function

!_______________________________________________________
!

function mat_mult_vec_TYPE(A,B) result(C)
implicit none
REAL_TYPE, intent(in) :: A(:,:), B(:)
REAL_TYPE :: C(size(A,1))

if (size(A,2).ne.size(B)) then
  write(stderr,*) 'smat_mult_svec: size not conform: A.x.B '
  write(stderr,*) 'shape(A) ',shape(A)
  write(stderr,*) 'shape(B) ',shape(B)
  stop
end if

call gemv_TYPE('n',size(A,1),size(A,2),1._4,A,size(A,1), &
     B,1, 0._4,C,1)

end function

!_______________________________________________________
!

!_______________________________________________________
!

function vec_mult_ssparsemat_TYPE(A,B) result(C)
implicit none
REAL_TYPE, intent(in)       :: A(:)
type(SparseMatrix), intent(in) :: B
REAL_TYPE :: C(B%n)
integer :: k,l

if (B%m.ne.size(A)) then
  write(stderr,*) 'svec_mult_ssparsemat: size not conform: A.x.B '
  write(stderr,*) 'shape(A) ',shape(A)
  write(stderr,*) 'shape(B) ',B%m,B%n
  stop
end if

C = 0

  do k=1,B%nz
    C(B%j(k)) = C(B%j(k)) + A(B%i(k)) * B%s(k)
  end do

end function

!_______________________________________________________
!


!_______________________________________________________
!

function vec_mult_mat_TYPE(A,B) result(C)
implicit none
REAL_TYPE, intent(in) :: A(:), B(:,:)
REAL_TYPE :: C(size(B,2))

! B' A

call gemv_TYPE('t',size(B,1),size(B,2),1._4,B,size(B,1), &
     A,1, 0._4,C,1)

end function

!_______________________________________________________
!

function vec_mult_vec_TYPE(A,B) result(C)
implicit none
REAL_TYPE, intent(in) :: A(:), B(size(A))
REAL_TYPE :: C, dot_TYPE

if (size(A).ne.size(B)) then
  write(stderr,*) 'svec_mult_svec: size not conform: A.x.B '
  write(stderr,*) 'shape(A) ',shape(A)
  write(stderr,*) 'shape(B) ',shape(B)
  stop
end if

C = dot_TYPE(size(A), A, 1, B, 1)
end function




!_______________________________________________________
!

function ssparsemat_mult_vec_TYPE(A,B) result(C)
implicit none
type(SparseMatrix), intent(in) :: A
REAL_TYPE, intent(in) :: B(:)
REAL_TYPE :: C(A%m)
integer :: k

if (A%n.ne.size(B)) then
  write(stderr,*) 'ssparcemat_mult_svec: size not conform: A.x.B '
  write(stderr,*) 'shape(A) ',A%m,A%n
  write(stderr,*) 'shape(B) ',shape(B)
  stop
end if

C = 0

do k=1,A%nz
  C(A%i(k)) = C(A%i(k)) + A%s(k) * B(A%j(k))
end do

end function

!_______________________________________________________
!

function ssparsemat_mult_mat_TYPE(A,B) result(C)
implicit none
type(SparseMatrix), intent(in) :: A
REAL_TYPE, intent(in) :: B(:,:)
REAL_TYPE :: C(A%m,size(B,2))
integer :: k,l

if (A%n.ne.size(B,1)) then
  write(stderr,*) 'ssparcemat_mult_smat: size not conform: A.x.B '
  write(stderr,*) 'shape(A) ',A%m,A%n
  write(stderr,*) 'shape(B) ',shape(B)
  stop
end if

C = 0

do l=1,size(B,2)
  do k=1,A%nz
    C(A%i(k),l) = C(A%i(k),l) + A%s(k) * B(A%j(k),l)
  end do
end do

end function

!_______________________________________________________
!

function mat_mult_ssparsemat_TYPE(A,B) result(C)
implicit none
REAL_TYPE, intent(in)       :: A(:,:)
type(SparseMatrix), intent(in) :: B
REAL_TYPE :: C(size(A,1),B%n)
integer :: k,l

if (B%m.ne.size(A,2)) then
  write(stderr,*) 'smat_mult_ssparsemat: size not conform: A.x.B '
  write(stderr,*) 'shape(A) ',shape(A)
  write(stderr,*) 'shape(B) ',B%m,B%n
  stop
end if

C = 0

do k=1,B%nz
  do l=1,size(A,1)
    C(l,B%j(k)) = C(l,B%j(k)) + A(l,B%i(k)) * B%s(k)
  end do
end do

end function



!_______________________________________________________
!



!_______________________________________________________
!
! operator: .xt.
!_______________________________________________________
!

function mat_mult_matT_TYPE(A,B) result(C)
implicit none
REAL_TYPE, intent(in) :: A(:,:), B(:,:)
REAL_TYPE :: C(size(A,1),size(B,1))

call gemm_TYPE('n','t',size(A,1),size(B,1),size(A,2),1._4, &
  A, size(A,1), B, size(B,1), &
  0._4, C, size(A,1))

end function

!_______________________________________________________
!

function mat_mult_ssparsematT_TYPE(A,B) result(C)
implicit none
REAL_TYPE, intent(in)       :: A(:,:)
type(SparseMatrix), intent(in) :: B
REAL_TYPE :: C(size(A,1),B%m)
integer :: k,l

if (B%n.ne.size(A,2)) then
  write(stderr,*) 'smat_mult_ssparsematT: size not conform: A.xt.B '
  write(stderr,*) 'shape(A) ',shape(A)
  write(stderr,*) 'shape(B) ',B%m,B%n
  stop
end if

C = 0

do k=1,B%nz
  do l=1,size(A,1)
    C(l,B%i(k)) = C(l,B%i(k)) + A(l,B%j(k)) * B%s(k)
  end do
end do

end function

!_______________________________________________________
!

!_______________________________________________________
!

function vec_mult_vecT_TYPE(A,B) result(C)
implicit none
REAL_TYPE, intent(in) :: A(:), B(:)
REAL_TYPE :: C(size(A),size(B))

integer :: i,j

do j=1,size(B)
do i=1,size(A)
  C(i,j) = A(i)*B(j)
end do
end do
end function


!_______________________________________________________
!
! operator: .tx.
!_______________________________________________________
!

function matT_mult_mat_TYPE(A,B) result(C)
implicit none
REAL_TYPE, intent(in) :: A(:,:), B(:,:)
REAL_TYPE :: C(size(A,2),size(B,2))

call gemm_TYPE ('t','n',size(A,2),size(B,2),size(A,1),1._4, &
  A, size(A,1), B, size(B,1), &
  0._4, C, size(A,2))

end function

!_______________________________________________________
!

function matT_mult_vec_TYPE(A,B) result(C)
implicit none
REAL_TYPE, intent(in) :: A(:,:), B(:)
REAL_TYPE :: C(size(A,2))

C = vec_mult_mat_TYPE(B,A)

end function

!_______________________________________________________
!
! operator: .dx.
!_______________________________________________________
!

function diag_mult_mat_TYPE(A,B) result(C)
implicit none
REAL_TYPE, intent(in) :: A(:), B(:,:)
REAL_TYPE :: C(size(B,1),size(B,2))

integer :: i,j

do j=1,size(B,2)
  do i=1,size(B,1)
    C(i,j) = A(i)*B(i,j)
  end do
end do
end function

!_______________________________________________________
!

function diag_mult_vec_TYPE(A,B) result(C)
implicit none
REAL_TYPE, intent(in) :: A(:), B(:)
REAL_TYPE :: C(size(B,1))

C = A*B
end function

!_______________________________________________________
!
! operator: .xd.
!_______________________________________________________
!

function mat_mult_diag_TYPE(A,B) result(C)
implicit none
REAL_TYPE, intent(in) :: A(:,:), B(:)
REAL_TYPE :: C(size(A,1),size(A,2))

integer :: i,j

do j=1,size(A,2)
  do i=1,size(A,1)
    C(i,j) = A(i,j)*B(j)
  end do
end do
end function

!_______________________________________________________
!
! inv
!_______________________________________________________
!

!
! PGI compiler bug: size(B,1) == 0
!

function inv_TYPE(A,det) result(B)
 implicit none
 REAL_TYPE, intent(in) :: A(:,:)
 REAL_TYPE, optional, intent(out) :: det

 REAL_TYPE :: B(size(A,1),size(A,2))

 integer :: IPIV(min(size(A,1),size(A,2))), info, lwork,i
 REAL_TYPE :: worksize
 REAL_TYPE, allocatable :: work(:)


! write(6,*) 'shape ',shape(B),size(B,1),size(A,1)
 B = A
 call getrf_TYPE(size(A,1),size(A,2), B, size(A,1), IPIV, INFO )

 if (present(det)) then
   det = 1.
   do i=1,size(B,1)
     if (ipiv(i).ne.i) then
       det = -det * b(i,i)
     else
       det = det * b(i,i)
     end if
   end do
 end if

 lwork = -1
 call getri_TYPE(size(A,1), B, size(A,1), IPIV, worksize, -1, INFO )
 lwork = worksize+.5
 allocate(work(lwork))

 call getri_TYPE(size(A,1), B, size(A,1), IPIV, work, lwork, INFO )
 deallocate(work)


end function 


!_______________________________________________________
!
! determinant
!_______________________________________________________
!

function det_TYPE(A) result(det)
 implicit none
 REAL_TYPE, intent(in) :: A(:,:)

 REAL_TYPE :: B(size(A,1),size(A,2))
 integer :: IPIV(min(size(A,1),size(A,2))), info, lwork,i
 REAL_TYPE :: worksize
 REAL_TYPE :: det
 REAL_TYPE, allocatable :: work(:)

 B = A
 call getrf_TYPE(size(B,1),size(B,2), B, size(B,1), IPIV, INFO )
 det = 0.
 if (info.ne.0) then
   return
 end if
 det = 1.
 do i=1,size(B,1)
   if (ipiv(i).ne.i) then
     det = -det * b(i,i)
   else
     det = det * b(i,i)
   end if

 end do

end function 




! computes all singular vectors and aingular values

function svd_TYPE(A,U,V,VT,work,INFO) result(S)
 implicit none
 REAL_TYPE, intent(in) :: A(:,:)
 REAL_TYPE  :: S(min(size(A,1),size(A,2)))
 REAL_TYPE, optional, target, intent(out) :: U(:,:), VT(:,:),  V(:,:), work(:)
 integer, optional, intent(out) :: info

 character :: jobu,jobvt
 REAL_TYPE, pointer :: pU(:,:), pVT(:,:), pwork(:)
 REAL_TYPE, target :: dummy(1,1)
 REAL_TYPE :: rlwork
 integer :: lwork, myinfo,ldU,ldVT

#ifndef ALLOCATE_LOCAL_VARS
 REAL_TYPE :: copyA(size(A,1),size(A,2))
#else
 REAL_TYPE, pointer :: copyA(:,:)
 allocate(copyA(size(A,1),size(A,2)))
#endif

 jobu='n'
 jobvt='n'
 pU => dummy
 pVT => dummy

 if (present(U)) then
   jobu='a'
   pU => U
 end if

 if (present(VT)) then
   jobvt='a'
   pVT => VT 
 else
   if (present(V)) then
     jobvt='a'
     pVT => V   
   end if
 end if

 copyA = A

 if (present(work)) then
   pwork => work
 else
   !  write(6,*) 'Normalise X1',size(pU,1),size(A,1)
   !  write(6,*) 'Normalise A', shape(    A)


   call gesvd_TYPE( JOBU, JOBVT,size(A,1),size(A,2), copyA, size(A,1), &
        S, pU,size(pU,1), pVT,size(pVT,1), &
        rlWORK, -1, myINFO)

!!$  write(6,*) 'Normalise X2', rlwork,shape(copyA),jobu,jobvt
!!$  write(6,*) 'Normalise pVT', shape( pVT   )
!!$  write(6,*) 'Normalise Pu', shape( pU)
!!$  write(6,*) 'Normalise S ', shape( S)
!!$  write(6,*) 'Normalise S ', "GESVD_TYPE"
!!$  write(6,*) 'Normalise S ', kind(S),kind(copyA),kind(work)

   lwork = rlwork+0.5
   allocate(pwork(lwork))
 end if



 call gesvd_TYPE( JOBU, JOBVT,size(A,1),size(A,2), copyA, size(A,1), &
      S, pU,size(pU,1), pVT,size(pVT,1), &
      pWORK, size(pWork), myINFO)

 if (.not.present(work)) deallocate(pwork)

!  write(6,*) 'Normalise X3'

#ifdef ALLOCATE_LOCAL_VARS
deallocate(copyA)
#endif

if (present(V))    V = transpose(pVT);
if (present(info)) info = myinfo;

end function



!_______________________________________________________
!


! computes the first singular vectors and singular values

!_______________________________________________________
!


subroutine ssvd_TYPE(job,X,U,S,V,info)
implicit none

! if job = 'N' or 'n' then X = U S V'
! if job = 'S' or 's' then X = U S^{-1/2} V'

character, intent(in) :: job

! m max space index (m => n)
! n max time index
! r number of eofs wanted (r <= n)

integer :: m, n, r
REAL_TYPE, intent(in) :: X(:,:)
REAL_TYPE, intent(out) :: U(:,:), S(:), V(:,:)

! if info = 0 every thing is fine
! if info <> 0 convergence failed

integer, intent(out) :: info

! 1 <= i <= m
! 1 <= j,jp <= n
! 1 <= k <= m

integer i,j,jp,k

! upper part of the time covariance 

REAL_TYPE, allocatable :: P(:)
integer idummy

! working space for sspevx

REAL_TYPE, allocatable :: work(:)
integer, allocatable :: iwork(:), ifail(:)

! BLAS scalar dot product

REAL_TYPE dot_TYPE

! LAPACK Machine precision routine

REAL_TYPE lamch_TYPE


m = size(X,1)
n = size(X,2)
r = size(U,2)

allocate(P((n*(n+1))/2),work(max(8*n,m)),iwork(5*n), ifail(n))

P(1) = lamch_TYPE('S')

do jp=1,n
  do j=1,jp
    P(j + (jp-1)*jp/2) = dot_TYPE(m,X(1,j),1,X(1,jp),1)
  enddo
enddo

call spevx_TYPE('V','I','U',n,P,-1.,-1.,n-r+1,n,2*lamch_TYPE('S'),idummy,S,V,n,work,iwork,ifail,info)

do k=1,r
  do i=1,m
    U(i,k) = dot_TYPE(n,X(i,1),m,V(1,k),1) / sqrt(S(k))
  enddo
enddo

if ((job.eq.'N').or.(job.eq.'n')) then
  do k=1,r
    S(k) = sqrt(S(k))
  enddo
endif

! sort in descending order

do k=1,r/2
  work(1) = S(k)
  S(k) = S(r-k+1)
  S(r-k+1) = work(1)

  call copy_TYPE(m,U(1,k),1,work,1)
  call copy_TYPE(m,U(1,r-k+1),1,U(1,k),1)
  call copy_TYPE(m,work,1,U(1,r-k+1),1)

  call copy_TYPE(n,V(1,k),1,work,1)
  call copy_TYPE(n,V(1,r-k+1),1,V(1,k),1)
  call copy_TYPE(n,work,1,V(1,r-k+1),1)
enddo

deallocate(P,work,iwork, ifail)

end subroutine





!_______________________________________________________
!

subroutine gesvd_f90_TYPE(JOBU,JOBVT,A,S,U,VT,INFO)
implicit none
character, intent(in) :: jobu,jobvt
REAL_TYPE, intent(in) :: A(:,:)
REAL_TYPE, intent(out) :: S(:), U(:,:), VT(:,:)
integer, intent(out) :: info

REAL_TYPE :: dummy,rlwork
integer :: lwork
REAL_TYPE, allocatable :: work(:)

#ifndef ALLOCATE_LOCAL_VARS
  REAL_TYPE :: copyA(size(A,1),size(A,2))
#else
  REAL_TYPE, pointer :: copyA(:,:)
  allocate(copyA(size(A,1),size(A,2)))
#endif


copyA = A

call gesvd_TYPE( JOBU, JOBVT,size(A,1),size(A,2), copyA, size(A,1), &
   S, U,size(U,1), VT,size(VT,1), &
     rlWORK, -1, INFO)
lwork = rlwork+0.5
allocate(work(lwork))

call gesvd_TYPE( JOBU, JOBVT,size(A,1),size(A,2), copyA, size(A,1), &
   S, U,size(U,1), VT,size(VT,1), &
     WORK, lwork, INFO)
deallocate(work)
#ifdef ALLOCATE_LOCAL_VARS
  deallocate(copyA)
#endif

end subroutine 


!_______________________________________________________
!
! computes eigenvalue/-vector of a symetric matrix
!
! A = V' diag(E) V 

subroutine symeig_TYPE(A,E,V,nbiggest,nsmallest,indices,INFO)
implicit none
REAL_TYPE, intent(in)  :: A(:,:)
REAL_TYPE,         intent(out) :: E(size(A,1))
integer, optional, intent(in) :: nbiggest, nsmallest,indices(2)
REAL_TYPE, optional, target, intent(out) :: V(:,:)
integer, optional, intent(out) :: info

character :: jobz
REAL_TYPE, pointer :: pV(:,:)
REAL_TYPE :: rlwork, tmp
integer :: lwork, myinfo, N, iwork(5*size(A,1)), ifail(size(A,1)), i,j
REAL_TYPE, allocatable :: work(:)

! LAPACK Machine precision routine

real :: slamch
integer :: r,idummy,ind(2)


#ifndef ALLOCATE_LOCAL_VARS
  REAL_TYPE :: copyA(size(A,1),size(A,2))
#else
  REAL_TYPE, pointer :: copyA(:,:)
  allocate(copyA(size(A,1),size(A,2)))
#endif



jobz='n'
N = size(A,1)
r = n

if (present(V)) then
  jobz='v'
  pV => V
else
  allocate(pV(1,1))
end if

ind = (/ 1,n /)

if (present(nbiggest))  ind = (/ n-nbiggest+1,n /)
if (present(nsmallest)) ind = (/ 1,nsmallest /)
if (present(indices))   ind = indices

! protect content of A

copyA = A

! determine the optimal size of work

call syevx_TYPE(JOBZ,'I','U',n,copyA,n,-1.,-1.,ind(1),ind(2),     &
     2*SLAMCH('S'),idummy,E,pV,n,rlWORK,-1, IWORK,   &
     IFAIL, myINFO )

lwork = rlwork+0.5
allocate(work(lwork))

call syevx_TYPE(JOBZ,'I','U',n,copyA,n,-1.,-1.,ind(1),ind(2),     &
     2*SLAMCH('S'),idummy,E,pV,n, WORK, LWORK, IWORK,   &
     IFAIL, myINFO )

  if (present(nbiggest)) then
! sort in descending order
    do i=1,nbiggest/2
      tmp = E(i)
      E(i) = E(nbiggest-i+1)
      E(nbiggest-i+1) = tmp

      if (present(V)) then
       do j=1,n
        tmp = V(j,i) 
        V(j,i) = V(j,nbiggest-i+1)
        V(j,nbiggest-i+1) = tmp
       end do
      end if
    end do
  end if

deallocate(work)
#ifdef ALLOCATE_LOCAL_VARS
  deallocate(copyA)
#endif


if (.not.present(V)) deallocate(pV)
if (present(info)) info = myinfo
end subroutine

!_______________________________________________________
!


