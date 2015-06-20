#ifdef OAK

module oak

 integer :: comm_all
 ! communicator between ensemble members of the same sub-domain
 integer :: comm_da
 logical :: model_uses_mpi
 ! ensemble size
 integer :: Nens
! integer :: n
 integer, allocatable :: startIndex(:),endIndex(:)

 ! size of local distributed state vector
 integer, allocatable :: locsize(:)

 integer :: schemetype = 1

 contains

 subroutine oak_init(comm,comm_ensmember)
  use mpi
  implicit none
  
  ! size of model 
  integer, intent(in), optional :: comm
  integer, intent(out), optional :: comm_ensmember
  integer :: nprocs
  integer :: ierr

  Nens = 2
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



 subroutine oak_domain(nl)
  use mpi
  implicit none
  
  ! size of model 
  integer, intent(in) :: nl
  integer, allocatable :: itemp(:)
  integer :: i

  allocate(startIndex(Nens),endIndex(Nens),locsize(Nens),itemp(Nens+1))
  

  ! indices for distributed state vector
  itemp = [ (1 + ( (i-1)*nl)/Nens, i=1,Nens+1)]
  !write(6,*) 'itemp ',itemp

  startIndex = itemp(1:Nens)
  endIndex = itemp(2:Nens+1)-1

  ! local size for assimilation
  locsize = endIndex - startIndex+1

  deallocate(itemp)

 end subroutine 


 subroutine oak_done()
  use mpi
  implicit none
  integer :: ierr

  if (.not.model_uses_mpi) then
    call mpi_finalize(ierr)
  end if

  deallocate(startIndex,endIndex,locsize)
 end subroutine oak_done

!------------------------------------------------------------------------------
!
! split communicator in N sub-groups, either by having continous ranks (mode = 
! true) or by having ranks separated by size/N (where size is the group of the 
! communicator comm)
!

 subroutine oak_split_comm(comm,N,new_comm,mode)
  implicit none

  integer, intent(in) :: comm
  ! number of sub-groups
  integer, intent(in) :: N
  logical, intent(in) :: mode
  integer, intent(out) :: new_comm
  
  integer :: rank, nprocs
  integer :: ierr
  integer :: group, new_group
  integer :: new_group_nprocs, i
  integer, allocatable :: ranks(:)
  
  call mpi_comm_rank(comm, rank, ierr)
  call mpi_comm_size(comm, nprocs, ierr)
  
  if (modulo(nprocs,N) /= 0) then
    print *, 'Error: number of processes (',nprocs, &
         ') must be divisible by the number ',N    
    call mpi_finalize(ierr)
    stop
  endif
  
  new_group_nprocs = nprocs/N
  !   write(6,*) 'new_group_nprocs',new_group_nprocs

  ! ranks belonging to the same ensemble member group   

  allocate(ranks(new_group_nprocs))

  if (mode) then
    ranks = [(i,i=0,new_group_nprocs-1)] + &
         new_group_nprocs*(rank/new_group_nprocs)
  else
    ranks = [((nprocs/new_group_nprocs) * i,i=0,new_group_nprocs-1)] + &
         modulo(rank,N)
  end if

  !write(6,*) 'ranks ',rank,':',ranks

!  Extract the original group handle
  call mpi_comm_group(comm, group, ierr)
  
!  Divide group into N distinct groups based upon rank
  call mpi_group_incl(group, new_group_nprocs, ranks,  &
       new_group, ierr)
  
  call mpi_comm_create(comm, new_group,  &
       new_comm, ierr)
  
  deallocate(ranks)
  
  call mpi_barrier(comm, ierr )
 end subroutine oak_split_comm

!-------------------------------------------------------------

 subroutine oak_assim(ntime,x)
  use mpi
  implicit none
  integer, intent(in) :: ntime
  real, intent(inout) :: x(:)
  integer :: ierr, i, obs_dim, r, myrank
  real, allocatable :: xloc(:,:) 
  integer, dimension(Nens) :: recvcounts, sdispls, rdispls

  real, allocatable :: Hx(:),HE(:,:), d(:), meanHx(:), HSf(:)
  real, pointer :: yo(:), invsqrtR(:)
  real :: scale

 write(6,*) 'x',x
  
  call mpi_comm_rank(comm_da, myrank, ierr)

  ! x(space) -> x(space/Nens,Nens)
!  x(2:n/Nens)
 
 ! collect the local part of the state vector

 recvcounts = locsize(myrank+1)

! allocate(xloc(sum(recvcounts)))
 allocate(xloc(locsize(myrank+1),Nens))

 write(6,*) 'model ',myrank,'counts ',locsize,recvcounts

! The type signature associated with sendcount[j], sendtype at process i must be equal to the type signature associated with recvcount[i], recvtype at process j. This implies that the amount of data sent must be equal to the amount of data received, pairwise between every pair of processes.

 sdispls = startIndex-1
 rdispls(1) = 0
 do i=2,Nens
   rdispls(i) = rdispls(i-1) + recvcounts(i-1)
 end do

 !write(6,*) 'model ',myrank,'disp ',sdispls,rdispls

 call mpi_alltoallv(x, locsize, sdispls, mpi_real, &
      xloc, recvcounts, rdispls, mpi_real, comm_da, ierr)


 write(6,*) 'xloc1',xloc(:,1)
 write(6,*) 'xloc2',xloc(:,2)


 !write(6,*) 'model ',myrank,'has loc ',xloc

 ! analysis

 if (schemetype == 1) then
   ! global
   
 else

 end if
 

 call mpi_alltoallv(xloc, recvcounts, rdispls, mpi_real, &
      x, locsize, sdispls, mpi_real, comm_da, ierr)

 !write(6,*) 'model ',myrank,'has global',x
  
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

 n = 8
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


#ifdef OAK
 call oak_domain(nl)
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
   call oak_assim(i,x(k0:k1))
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
