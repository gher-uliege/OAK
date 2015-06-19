program toymodel

#ifdef MPI
 use mpi
#endif

implicit none
integer :: n
integer, parameter :: Ntime = 1000, halo = 1

real, allocatable :: x(:),xinit(:)
integer :: i
integer :: j0, j1, k0, k1, nl

#ifdef MPI
integer :: comm, rank, nprocs, ierr, tag
integer :: precision = MPI_DOUBLE_PRECISION
#endif

n = 1000

#ifdef MPI
 comm = mpi_comm_world
 call mpi_init(ierr)
 call mpi_comm_size(comm, nprocs, ierr)
 call mpi_comm_rank(comm, rank, ierr)

 ! indexes in global vector
 j0 = (rank * n)/ncpus + 1
 j1 = ((rank+1) * n)/ncpus

 ! indexes in local vector
 nl = j1 - j0 + 1
 k0 = halo+1
 k1 = halo+nl

 tag = 1
 write(6,*) 'size rank master',nprocs,rank,comm,j0,j1,k0,k1
#else
 j0 = 1
 j1 = n

 k0 = halo+1
 k1 = halo+n
 nl = n
#endif

! initialize j1-j0+1 grid points and two halo points
!allocate(x(j0-halo:j1+halo),xinit(j0-halo:j1+halo))

allocate(x(1:nl+2*halo),xinit(1:nl+2*halo))
xinit = 0
if (n/2 <= j1+halo) then
  xinit(n/2 - (j0-halo)) = 1
end if

x = xinit

! time loop
do i = 1,Ntime  
  ! single time step
  call step(i,x)

  ! boundary conditions
  call bc(i,x)
end do

write(6,*) 'sum(x) ',sum(x), maxval(abs(x - xinit))
if (maxval(abs(x - xinit)) > 1e-6) then
  write(6,*) 'fail'
else
  write(6,*) 'OK'
end if
deallocate(x,xinit)

#ifdef MPI
call mpi_finalize(ierr)
#endif

contains

 subroutine step(i,x)
  implicit none
  integer, intent(in) :: i
  real, intent(inout) :: x(:)

  x(k0:k1) = x(k0+1:k1+1)

 end subroutine step


 subroutine bc(i,x)
  integer, intent(in) :: i
  real, intent(inout) :: x(:)

#ifndef MPI

!  write(6,*) 'x',shape(x),j0,j0+halo-1
!  write(6,*) 'x',shape(x),j1+1,j1+halo
  x(k1+1:k1+halo) = x(k0:k0+halo-1)
  x(k0-halo:k0-1) = x(k1-halo+1:k1)

#else

!  Send to the "left"
!
 write(6,*) 'send to',modulo(rank-1,ncpus)
 call mpi_send(x(k0:k0+halo-1), halo, precision, modulo(rank-1,ncpus), tag, &
        comm, ierr )

 call mpi_recv(x(k1+1:k1+halo), halo, precision, modulo(rank+1,ncpus), tag, &
        comm, status, ierr )

 tag = tag+1

!  Send to the "right"
!
 write(6,*) 'send to',modulo(rank+1,ncpus)
 call mpi_send(x(k1-halo+1:k1), halo, precision, modulo(rank+1,ncpus), tag, &
        comm, ierr )

 call mpi_recv(x(k0-halo:k0-1), halo, precision, modulo(rank-1,ncpus), tag, &
        comm, status, ierr )

 tag = tag+1

#endif
  
 end subroutine bc
end program toymodel
