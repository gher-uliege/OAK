#ifdef OAK

module oak

 integer :: comm_all
 ! communicator of data assimilation
 integer :: comm_da
 logical :: model_uses_mpi
 contains

 subroutine oak_init(comm,comm_ensmember)
  use mpi
  implicit none
  
  integer, intent(in), optional :: comm
  integer, intent(out), optional :: comm_ensmember
  integer :: Nens = 2, nprocs
  integer :: ierr

  if (present(comm)) then
    ! model uses also MPI, split communicators
    call oak_split_comm(comm,Nens,comm_ensmember,.true.)
    comm_all = comm
    model_uses_mpi = .true.
  else
    ! modes does not use MPI, we can use MPI alone
    comm_all = mpi_comm_world
    call mpi_init(ierr)
    model_uses_mpi = .false.
  end if

  call mpi_comm_size(comm_all, nprocs, ierr)
  call oak_split_comm(comm_all,nprocs/Nens,comm_da,.false.)
  
 end subroutine oak_init


 subroutine oak_done()
  use mpi
  implicit none
  integer :: ierr

  if (.not.model_uses_mpi) then
    call mpi_finalize(ierr)
  end if
 end subroutine oak_done

!------------------------------------------------------------------------------
!
! split communicator in n sub-groups, either by having continous ranks (mode = 
! true) or by having ranks separated by size/n (where size is the group of the 
! communicator comm)
!

 subroutine oak_split_comm(comm,N,comm_ensmember,mode)
  implicit none

  integer, intent(in) :: comm
  ! number of ensemble members
  integer, intent(in) :: N
  logical, intent(in) :: mode
  integer, intent(out) :: comm_ensmember
  
  integer :: rank, nprocs
  integer :: ierr
  integer :: group, group_ensmember
  integer :: nprocs_ensmember, i
  integer, allocatable :: ranks(:)
  
  call mpi_comm_rank(comm, rank, ierr)
  call mpi_comm_size(comm, nprocs, ierr)
  
  if (modulo(nprocs,N) /= 0) then
    print *, 'Error: number of processes (',nprocs,') must be divisible by number of ensemble members (',N,')'
    
    call mpi_finalize(ierr)
    stop
  endif
  
  nprocs_ensmember = nprocs/N
  !   write(6,*) 'nprocs_ensmember',nprocs_ensmember

  ! ranks belonging to the same ensemble member group   

  allocate(ranks(nprocs_ensmember))

  if (mode) then
    ranks = [(i,i=0,nprocs_ensmember-1)] + &
         nprocs_ensmember*(rank/nprocs_ensmember)
  else
! CONTINUE HERE
    ranks = [((nprocs/nprocs_ensmember) * i,i=0,nprocs_ensmember-1)] + &
         modulo(rank,N)
  end if

  write(6,*) 'ranks ',rank,':',ranks

!  Extract the original group handle
  call mpi_comm_group(comm, group, ierr)
  
!  Divide group into N distinct groups based upon rank
  call mpi_group_incl(group, nprocs_ensmember, ranks,  &
       group_ensmember, ierr)
  
  call mpi_comm_create(comm, group_ensmember,  &
       comm_ensmember, ierr)
  
  deallocate(ranks)
  
  call mpi_barrier(comm, ierr )
 end subroutine oak_split_comm

 subroutine oak_assim(i,x)
  use mpi
  implicit none
  integer, intent(in) :: i
  real, intent(inout) :: x(:)

  integer :: n = 100
  

  ! x(space) -> x(space/Nens,Nens)
!  x(2:n/Nens)
  
  
  
 end subroutine oak_assim



end module oak
#endif


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
 !integer :: precision = mpi_double_precision
 integer :: precision = mpi_real
 integer :: status(mpi_status_size)
#endif

 n = 100
 Ntime = n

#ifdef MODEL_PARALLEL

 call mpi_init(ierr)

#ifdef OAK
 call oak_init(mpi_comm_world,comm)
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
 rank = 0
 j0 = 1
 j1 = n

 k0 = halo+1
 k1 = halo+n
 nl = n

#ifdef OAK
 call oak_init()
#endif

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

#ifdef OAK
   call oak_assim(i,x)
#endif
 end do

! write(6,*) 'x ',x, 'rank',rank
! write(6,*) 'xinit ',xinit, 'rank',rank

 call check_results(x)

 deallocate(x,xinit)

#ifdef OAK
 call oak_done()
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


# ifdef MODEL_PARALLEL
 call mpi_reduce(sum(x(k0:k1)), sumx, 1, precision, mpi_sum, 0, comm, ierr)
#else
 sumx = sum(x(k0:k1))
#endif 

 if (rank == 0) then
   if (abs(sumx - n*(n+1)/2) > 1e-6) then
     write(6,*) 'sum: FAIL'
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
