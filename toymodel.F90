program toymodel

#ifdef MPI
 use mpi
#endif

implicit none
integer :: n, Ntime
integer, parameter :: halo = 1

real, allocatable :: x(:),xinit(:)
integer :: i
integer :: j0, j1, k0, k1, nl

#ifdef MPI
integer :: comm, rank, nprocs, ierr, tag
!integer :: precision = MPI_DOUBLE_PRECISION
integer :: precision = MPI_REAL
integer status(MPI_STATUS_SIZE)
#endif

n = 5
Ntime = n

#ifdef MPI
 comm = mpi_comm_world
 call mpi_init(ierr)
 call mpi_comm_size(comm, nprocs, ierr)
 call mpi_comm_rank(comm, rank, ierr)

 ! indexes in global vector
 j0 = (rank * n)/nprocs + 1
 j1 = ((rank+1) * n)/nprocs

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

do i=k0,k1
  xinit(i) = i-k0+j0
end do

 write(6,*) 'k0,k1,halo',k0,k1,halo

call bc(i,xinit)
write(6,*) 'xinit ',xinit, 'rank',rank

x = xinit

! time loop
do i = 1,Ntime  
  ! single time step
  call step(i,x)

  ! boundary conditions
  call bc(i,x)
end do

write(6,*) 'x ',x, 'rank',rank
write(6,*) 'xinit ',xinit, 'rank',rank

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
  implicit none
  integer, intent(in) :: i
  real, intent(inout) :: x(:)

  write(6,*) 'k0,k1,halo',k0,k1,halo
  write(6,*) 'nprocs',nprocs, rank

#ifndef MPI

!  write(6,*) 'x',shape(x),j0,j0+halo-1
!  write(6,*) 'x',shape(x),j1+1,j1+halo
  x(k1+1:k1+halo) = x(k0:k0+halo-1)
  x(k0-halo:k0-1) = x(k1-halo+1:k1)

#else

!  Send to the "left"
!
 write(6,*) 'send to',modulo(rank-1,nprocs)
 write(6,*) 'send to',k0,k0+halo-1


 call mpi_send(x(k0:k0+halo-1), halo, precision, modulo(rank-1,nprocs), tag, &
        comm, ierr )

 write(6,*) 'get from',k1+1,k1+halo, shape(x)
  write(6,*) 'nprocs',nprocs
 
 call mpi_recv(x(k1+1:k1+halo), halo, precision, modulo(rank+1,nprocs), tag, &
        comm, status, ierr )

  write(6,*) 'status, ierr',status, ierr, precision
  write(6,*) 'nprocs',nprocs
 tag = tag+1

!  Send to the "right"
!
 write(6,*) 'k0,k1,halo',k0,k1,halo
  write(6,*) 'nprocs',nprocs, rank
 write(6,*) 'send to',modulo(rank+1,nprocs),x(k1-halo+1:k1),k0-halo,k0-1
 call mpi_send(x(k1-halo+1:k1), halo, precision, modulo(rank+1,nprocs), tag, &
        comm, ierr )

 write(6,*) 'get from',modulo(rank-1,nprocs),k0-halo,k0-1
 call mpi_recv(x(k0-halo:k0-1), halo, precision, modulo(rank-1,nprocs), tag, &
        comm, status, ierr )

 tag = tag+1

#endif
  
 end subroutine bc
end program toymodel
