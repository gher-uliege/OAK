
program toymodel
#ifdef MODEL_PARALLEL
 use mpi
#endif
#ifdef OAK
 use oak
#endif

 implicit none
 integer :: n, Ntime
 integer, parameter :: halo = 1

 real, allocatable :: x(:),xinit(:)
 integer :: i
 integer :: j0, j1, k0, k1, nl
 real :: sumx

 integer :: rank
#ifdef MODEL_PARALLEL
 integer :: comm, nprocs, ierr, tag
 integer :: precision = mpi_double_precision
 !integer :: precision = mpi_real
 integer :: status(mpi_status_size)
#endif

#ifdef OAK
  type(oakconfig) :: config
#endif
 n = 8
 Ntime = n

#ifdef MODEL_PARALLEL

 call mpi_init(ierr)

#ifdef OAK
 call oak_init(config,mpi_comm_world,comm,fname='test_assim.init')
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
! write(6,*) 'size rank master',nprocs,rank,comm,j0,j1,k0,k1
 call mpi_barrier(comm, ierr)
#else
 ! model does not run in parallel
 nprocs = 1
 rank = 0
 
 j0 = 1
 j1 = n

 k0 = halo+1
 k1 = halo+n
 nl = n

#ifdef OAK
 call oak_init(config)
#endif

#endif


#ifdef OAK
 call oak_domain(config,nl,partition=[(i,i=j0,j1)])

 do i = 1,nprocs
   config%dom(((i-1) * n)/nprocs + 1:(i * n)/nprocs) = i   
 end do
 write(6,*) 'config%dom ',config%dom
#endif

 ! initialize j1-j0+1 grid points and two halo points

 allocate(x(1:nl+2*halo),xinit(1:nl+2*halo))
 xinit = 0

 do i=k0,k1
   xinit(i) = i-k0+j0
 end do

 call bc(i,xinit)

 x = xinit
! write(6,*) 'xinit ',xinit, 'rank',rank


 ! time loop
 do i = 1,Ntime  
   ! single time step
   call step(i,x)

   ! boundary conditions
   call bc(i,x)

#ifdef OAK
   call oak_assim(config,i,x(k0:k1))
#endif
 end do

! write(6,*) 'x ',x, 'rank',rank
! write(6,*) 'xinit ',xinit, 'rank',rank

 call check_results(x)

 deallocate(x,xinit)

#ifdef OAK
 call oak_done(config)
#endif
#ifdef MODEL_PARALLEL
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

#ifndef MODEL_PARALLEL
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


 subroutine check_results(x)
  implicit none
  real, intent(in) :: x(:)
  real :: sumx_ref

# ifdef MODEL_PARALLEL
 call mpi_reduce(sum(x(k0:k1)), sumx, 1, precision, mpi_sum, 0, comm, ierr)
#else
 sumx = sum(x(k0:k1))
#endif 

 if (rank == 0) then
   sumx_ref = n*(n+1)/2
   if (abs(sumx - sumx_ref) > 1e-6) then
     write(6,*) 'sum: FAIL'
     write(6,*) 'got: ',sumx
     write(6,*) 'explected: ',sumx_ref
   else
     write(6,*) 'sum:  OK'
   end if

 end if

! write(6,*) 'sum(x) ',sum(x), maxval(abs(x - xinit))
 if (maxval(abs(x - xinit)) > 1e-6) then
   write(6,*) 'sub-domain check: FAIL'
 else
   write(6,*) 'sub-domain check:  OK'
 end if

 end subroutine 


end program toymodel
