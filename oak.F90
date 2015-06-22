module oak

 type oakconfig
   integer :: comm_all
   ! communicator between ensemble members of the same sub-domain
   integer :: comm_da
   integer :: comm_ensmember
   logical :: model_uses_mpi
   ! ensemble size
   integer :: Nens
   ! integer :: n
   
   integer :: mpi_precision

   integer, allocatable :: startIndex(:),endIndex(:)
   
   ! size of local distributed state vector
   integer, allocatable :: locsize(:), partition(:)

   real, allocatable :: gridx(:), gridy(:), gridz(:), gridt(:)
   
   integer :: schemetype = 1   
 end type oakconfig

 integer :: comm_all
 ! communicator between ensemble members of the same sub-domain
 integer :: comm_da
 integer :: comm_ensmember
 logical :: model_uses_mpi
 ! ensemble size
 integer :: Nens
! integer :: n

 integer :: mpi_precision

 integer, allocatable :: startIndex(:),endIndex(:)

 ! size of local distributed state vector
 integer, allocatable :: locsize(:), partition(:)

 real, allocatable :: gridx(:), gridy(:), gridz(:), gridt(:)

 integer :: schemetype = 1

 contains

 subroutine oak_init(config,comm,comm_ensmember_)
  use mpi
  implicit none
  
  type(oakconfig), intent(out) :: config
  integer, intent(in), optional :: comm
  integer, intent(out), optional :: comm_ensmember_
  integer :: nprocs
  integer :: ierr

  Nens = 2
  if (present(comm)) then
    ! model uses also MPI, split communicators
    call oak_split_comm(comm,Nens,comm_ensmember_,.true.)
    comm_all = comm
    model_uses_mpi = .true.

    !write(6,*) 'comm_ensmember_',comm_ensmember_
    comm_ensmember = comm_ensmember_
  else
    ! modes does not use MPI, we can use MPI alone
    comm_all = mpi_comm_world
    ! comm_ensmember should never be used
    comm_ensmember = -1
    call mpi_init(ierr)
    model_uses_mpi = .false.
  end if

  call mpi_comm_size(comm_all, nprocs, ierr)
  call oak_split_comm(comm_all,nprocs/Nens,comm_da,.false.)

  if (kind(gridx) == 4) then
    mpi_precision = mpi_real
  elseif (kind(gridx) == 8) then
    mpi_precision = mpi_double_precision
  else
    write(6,*) 'unknown precision'
    stop
  end if
 end subroutine oak_init



 subroutine oak_domain(config,nl,part,gridx_,gridy_,gridz_,gridt_)
  use mpi
  implicit none

  type(oakconfig), intent(inout) :: config  
  ! size of model 
  integer, intent(in) :: nl
  integer, intent(in), optional :: part(:)
  real, intent(in), optional :: gridx_(:)
  real, intent(in), optional :: gridy_(:)
  real, intent(in), optional :: gridz_(:)
  real, intent(in), optional :: gridt_(:)
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

  if (present(part)) then
    allocate(partition(size(part)))
    partition = part
  end if

  if (present(gridx_)) then
    allocate(gridx(size(gridx_)))
    gridx = gridx_
  end if

  if (present(gridy_)) then
    allocate(gridy(size(gridy_)))
    gridy = gridy_
  end if

  if (present(gridz_)) then
    allocate(gridz(size(gridz_)))
    gridz = gridz_
  end if

  if (present(gridt_)) then
    allocate(gridt(size(gridt_)))
    gridt = gridt_
  end if

 end subroutine 


 subroutine oak_done(config)
  use mpi
  implicit none
  type(oakconfig), intent(inout) :: config
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


!  Extract the original group handle
  call mpi_comm_group(comm, group, ierr)
  
!  Divide group into N distinct groups based upon rank
  call mpi_group_incl(group, new_group_nprocs, ranks,  &
       new_group, ierr)
  
  call mpi_comm_create(comm, new_group,  &
       new_comm, ierr)
  
  deallocate(ranks)
  !write(6,*) 'ranks ',rank,':',ranks, new_comm
  
  call mpi_barrier(comm, ierr )
 end subroutine oak_split_comm


!-------------------------------------------------------------

subroutine oak_obsoper(config,x,Hx,comm_ensmember)
 use mpi
 implicit none
 type(oakconfig), intent(inout) :: config
 real, intent(in) :: x(:)
 real, intent(out) :: Hx(:)
 integer, intent(in) :: comm_ensmember
 integer :: rank, m = 1, ierr
 ! use communicator for each ensemble member


 if (model_uses_mpi) then
   call mpi_comm_rank(comm_ensmember, rank, ierr)
   if (rank == 0) then
     Hx = x(1) 
   end if
   
   call mpi_bcast(Hx, m, mpi_precision, 0, comm_ensmember, ierr)
 else
   Hx = x(1)
 end if
end subroutine oak_obsoper


!-------------------------------------------------------------

subroutine oak_load_obs(ntime,yo,invsqrtR)
 implicit none
 integer, intent(in) :: ntime
 real, intent(out), pointer :: yo(:), invsqrtR(:)

 allocate(yo(1),invsqrtR(1))
 yo = 1
 invsqrtR = 1  
end subroutine oak_load_obs


!-------------------------------------------------------------

subroutine selectObs(i,w) 
 implicit none
 integer, intent(in) :: i
 real, intent(out) :: w(:)
 
 
end subroutine selectObs

!-------------------------------------------------------------

 subroutine oak_assim(config,ntime,x)
  use mpi
  use sangoma_ensemble_analysis
  implicit none
  type(oakconfig), intent(inout) :: config
  integer, intent(in) :: ntime
  real, intent(inout) :: x(:)
  integer :: ierr, i, obs_dim, r, myrank
  real, allocatable :: xloc(:,:) 
  integer, dimension(Nens) :: recvcounts, sdispls, rdispls

  real, allocatable :: Hx(:),HEf(:,:), d(:), meanHx(:), HSf(:), diagR(:), Ef(:,:), Ea(:,:)
  real, pointer :: yo(:), invsqrtR(:)
  real :: scale
  character(len=64) :: method = 'ETKF'

  ! load observations
  call oak_load_obs(ntime,yo,invsqrtR)

  obs_dim = size(yo)
  allocate(Hx(obs_dim),HEf(obs_dim,Nens),d(obs_dim),HSf(obs_dim),meanHx(obs_dim),diagR(obs_dim))

  ! extract observations
  call oak_obsoper(config,x,Hx,comm_ensmember)

  call mpi_allgather(Hx, obs_dim, mpi_precision, HEf, obs_dim, mpi_precision, comm_da, ierr)

!  write(6,*) 'Hx',Hx

!  write(6,*) 'x',x
  
  call mpi_comm_rank(comm_da, myrank, ierr)
   
  ! collect the local part of the state vector
  
  recvcounts = locsize(myrank+1)
  
  ! allocate(xloc(sum(recvcounts)))
  allocate(xloc(locsize(myrank+1),Nens))

!  write(6,*) 'model ',myrank,'counts ',locsize,recvcounts

  ! The type signature associated with sendcount[j], sendtype at process i must be equal to the type signature associated with recvcount[i], recvtype at process j. This implies that the amount of data sent must be equal to the amount of data received, pairwise between every pair of processes.

  sdispls = startIndex-1
  rdispls(1) = 0
  do i=2,Nens
    rdispls(i) = rdispls(i-1) + recvcounts(i-1)
  end do
  
  !write(6,*) 'model ',myrank,'disp ',sdispls,rdispls
  
  call mpi_alltoallv(x, locsize, sdispls, mpi_precision, &
       xloc, recvcounts, rdispls, mpi_precision, comm_da, ierr)


  write(6,*) 'xloc1',xloc(:,1)
  write(6,*) 'xloc2',xloc(:,2)


  !write(6,*) 'model ',myrank,'has loc ',xloc
  
  ! analysis
  
  if (schemetype == 1) then
    ! global
    
  else
!    call local_ensemble_analysis(Ef,HEf,yo,diagR,partition,selectObs,method,Ea)
  end if
  

  call mpi_alltoallv(xloc, recvcounts, rdispls, mpi_precision, &
       x, locsize, sdispls, mpi_precision, comm_da, ierr)
  
 !write(6,*) 'model ',myrank,'has global',x
  
 end subroutine oak_assim



end module oak
