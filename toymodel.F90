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
 real :: sumx

 integer :: rank
#ifdef MPI
 integer :: comm, nprocs, ierr, tag
 !integer :: precision = MPI_DOUBLE_PRECISION
 integer :: precision = MPI_REAL
 integer :: status(MPI_STATUS_SIZE)
#endif

 n = 100
 Ntime = n

#ifdef MPI

 call mpi_init(ierr)

#ifdef OAK
 call init(mpi_comm_world,comm)
#else
 comm = mpi_comm_world
#endif

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
 call MPI_Barrier (comm, ierr )
#else
 rank = 0
 j0 = 1
 j1 = n

 k0 = halo+1
 k1 = halo+n
 nl = n
#endif

 ! initialize j1-j0+1 grid points and two halo points

 allocate(x(1:nl+2*halo),xinit(1:nl+2*halo))
 xinit = 0

 do i=k0,k1
   xinit(i) = i-k0+j0
 end do

 call bc(i,xinit)

 x = xinit

 ! time loop
 do i = 1,Ntime  
   ! single time step
   call step(i,x)

   ! boundary conditions
   call bc(i,x)
 end do

! write(6,*) 'x ',x, 'rank',rank
! write(6,*) 'xinit ',xinit, 'rank',rank

# ifdef MPI
 call MPI_Reduce(sum(x(k0:k1)), sumx, 1, precision, MPI_SUM, 0, comm, ierr)
#else
 sumx = sum(x(k0:k1))
#endif 


 if (rank == 0) then
   write(6,*) 'sum(x) ',sumx
 end if

! write(6,*) 'sum(x) ',sum(x), maxval(abs(x - xinit))
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

#ifndef MPI
  x(k1+1:k1+halo) = x(k0:k0+halo-1)
  x(k0-halo:k0-1) = x(k1-halo+1:k1)
#else
  !  Send to the "left"

  call mpi_send(x(k0:k0+halo-1), halo, precision, modulo(rank-1,nprocs), tag, &
       comm, ierr )

  call mpi_recv(x(k1+1:k1+halo), halo, precision, modulo(rank+1,nprocs), tag, &
       comm, status, ierr )

  tag = tag+1

  !  Send to the "right"

  call mpi_send(x(k1-halo+1:k1), halo, precision, modulo(rank+1,nprocs), tag, &
       comm, ierr )

  call mpi_recv(x(k0-halo:k0-1), halo, precision, modulo(rank-1,nprocs), tag, &
       comm, status, ierr )

  tag = tag+1
#endif

 end subroutine bc

 
 subroutine init(comm,comm_ensmember)
  implicit none

  integer, intent(in) :: comm
  integer, intent(out) :: comm_ensmember

     integer NPROCS
   parameter(NPROCS=8)
   integer rank, new_rank, sendbuf, recvbuf, numtasks, nprocs_check
   integer ranks1(4), ranks2(4), ierr
   integer orig_group, new_group, new_comm
   data ranks1 /0, 1, 2, 3/, ranks2 /4, 5, 6, 7/

   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs_check, ierr)

   if (nprocs_check .ne. NPROCS) then
     print *, 'Must specify NPROCS= ',NPROCS,' Terminating.'
     call MPI_FINALIZE(ierr)
     stop
   endif

   sendbuf = rank

!  Extract the original group handle
   call MPI_COMM_GROUP(MPI_COMM_WORLD, orig_group, ierr)

!  Divide tasks into two distinct groups based upon rank
   if (rank .lt. NPROCS/2) then
      call MPI_GROUP_INCL(orig_group, NPROCS/2, ranks1,  &
                   new_group, ierr)
   else 
      call MPI_GROUP_INCL(orig_group, NPROCS/2, ranks2,  &
                   new_group, ierr)
   endif

   call MPI_COMM_CREATE(MPI_COMM_WORLD, new_group,  &
                    comm_ensmember, ierr)
   call MPI_ALLREDUCE(sendbuf, recvbuf, 1, MPI_INTEGER, &
                    MPI_SUM, comm_ensmember, ierr)

   call MPI_GROUP_RANK(new_group, new_rank, ierr)
   print *, 'rank= ',rank,' newrank= ',new_rank,' recvbuf= ', &
       recvbuf


   call MPI_Barrier (MPI_COMM_WORLD, ierr )


!   stop
 end subroutine init

end program toymodel
