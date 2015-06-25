module oak
 use ndgrid
 use assimilation

 ! integer, parameter :: maxLen = 256

 type oakconfig
   character(len=maxLen) :: initfname

   ! communicator model + assimilation
   integer :: comm_all
   ! communicator between ensemble members (of the same sub-domain)
   integer :: comm_da
   ! communicator between sub-domains (of the same ensemble member)
   integer :: comm_ensmember
   ! flag wheter the model uses MPI
   logical :: model_uses_mpi
   ! ensemble size
   integer :: Nens

   ! integer :: n

   integer, allocatable :: startIndex(:),endIndex(:)

   ! size of local distributed state vector
   integer, allocatable :: locsize(:), partition(:)

!   logical, allocatable :: mask(:)
   real, allocatable :: gridx(:), gridy(:), gridz(:), gridt(:)

!   type(MemLayout) :: ModML
!   type(grid), allocatable :: ModelGrid(:)
   integer :: schemetype = 1
   integer :: obsntime = 1   
 end type oakconfig

contains

 subroutine oak_init(config,comm,comm_ensmember_,fname)
  use mpi
  use initfile
  implicit none

  type(oakconfig), intent(out) :: config
  integer, intent(in), optional :: comm
  integer, intent(out), optional :: comm_ensmember_
  character(len=maxLen), optional :: fname

  integer :: nprocs, nprocs_ensmember, vmax
  integer :: ierr

  config%initfname = ''

  if (present(fname)) then
    ! remove this one
    initfname = fname
    config%initfname = fname
    call oak_load_model_grid(config)

    call getInitValue(initfname,'ErrorSpace.dimension',config%Nens,default=0)

  end if

  if (present(comm)) then
    ! model uses also MPI, split communicators
    call oak_split_comm(comm,config%Nens,comm_ensmember_,.true.)
    config%comm_all = comm
    config%model_uses_mpi = .true.

    !write(6,*) 'comm_ensmember_',comm_ensmember_
    config%comm_ensmember = comm_ensmember_
    call mpi_comm_size(config%comm_ensmember, nprocs_ensmember, ierr)
  else
    ! modes does not use MPI, we can use MPI alone
    config%comm_all = mpi_comm_world
    ! comm_ensmember should never be used
    config%comm_ensmember = -1

    call mpi_init(ierr)
    config%model_uses_mpi = .false.
    nprocs_ensmember = 1
  end if

  call mpi_comm_size(config%comm_all, nprocs, ierr)
  call oak_split_comm(config%comm_all,nprocs/config%Nens,config%comm_da,.false.)



  ! print diagnostic information
  write(6,*) 'OAK is initialized'
  write(6,*) 'Total number of processes',nprocs
  write(6,*) 'Processes per ensemble member',nprocs_ensmember

 end subroutine oak_init

 !--------------------------------------------------------------------

 subroutine oak_load_model_grid(config)
  use mpi
  use initfile
  use ufileformat
  use parall
  implicit none
  type(oakconfig), intent(inout) :: config


  call parallInit(communicator=config%comm_da)
  call init(config%initfname)

  ! if (config%initfname /= '') then
  !   call getInitValue(initfname,'schemetype',schemetype,default=GlobalScheme)
  !   ! variables for local assimilation

  !   !call load_grid()

  !   if (schemetype.eq.LocalScheme) then
  !   !  call load_partition()
  !   end if


  ! end if

 contains

  subroutine load_grid()
   implicit none
   character(len=MaxFNameLength), pointer   ::  &
        filenamesX(:),filenamesY(:),filenamesZ(:),    &
        filenamesT(:)
   character(len=MaxFNameLength)            :: path
   integer :: vmax
   integer :: v,n


   ! Models Memory Layout 
   call MemoryLayout('Model.',ModML)

   vmax = ModML%nvar

   allocate(ModelGrid(vmax),hres(vmax))
   hres = 0

   call getInitValue(initfname,'Model.path',path,default='')

   !
   ! define model grid
   !

   call getInitValue(config%initfname,'Model.gridX',filenamesX)
   call getInitValue(config%initfname,'Model.gridY',filenamesY)
   call getInitValue(config%initfname,'Model.gridZ',filenamesZ)

   write(6,*) 'filenamesX',filenamesX(1)
   write(6,*) 'filenamesY',filenamesY(1)
   write(6,*) 'filenamesZ',filenamesZ(1),vmax,trim(path)
   do v=1,vmax
     n = ModML%ndim(v)

     ! initialisze the model grid structure ModelGrid(v)
     call initgrid(ModelGrid(v),n,ModML%varshape(1:n,v), &
          ModML%Mask(ModML%StartIndex(v):ModML%EndIndex(v)).eq.0)

     ! set the coordinates of the model grid
     call setCoord(ModelGrid(v),1,trim(path)//filenamesX(v))
     call setCoord(ModelGrid(v),2,trim(path)//filenamesY(v))

     if (n > 2) then
       call setCoord(ModelGrid(v),3,trim(path)//filenamesZ(v))


       if (n > 3) then
         call getInitValue(initfname,'Model.gridT',filenamesT)
         call setCoord(ModelGrid(v),4,trim(path)//filenamesT(v))
         deallocate(filenamesT)
       end if
     end if


   end do

   ! legacy for obs* routines
   !ModML = ModML
   !ModelGrid = ModelGrid

   deallocate(filenamesX,filenamesY,filenamesZ)
  end subroutine load_grid


  subroutine load_partition()
   implicit none
   real, pointer :: tmp(:)

   allocate(tmp(ModML%effsize),partition(ModML%effsize), &
        hCorrLengthToObs(ModML%effsize),hMaxCorrLengthToObs(ModML%effsize))

   call loadVector('Zones.partition',ModML,tmp)
   ! convertion real -> integer
   partition = nint(tmp)
   deallocate(tmp)

   call loadVector('Zones.corrLength',ModML,hCorrLengthToObs)
   call loadVector('Zones.maxLength',ModML,hMaxCorrLengthToObs)

   ! optimization: permute element in the vector so that element in the same
   ! sub-zone are close in memory

  end subroutine load_partition
 end subroutine oak_load_model_grid

 subroutine oak_domain(config,nl,partition,gridx,gridy,gridz,gridt)
  use mpi
  implicit none

  type(oakconfig), intent(inout) :: config  
  ! size of model 
  integer, intent(in) :: nl
  integer, intent(in), optional :: partition(:)
  real, intent(in), optional :: gridx(:)
  real, intent(in), optional :: gridy(:)
  real, intent(in), optional :: gridz(:)
  real, intent(in), optional :: gridt(:)
  integer, allocatable :: itemp(:)

  integer :: i

  allocate(config%startIndex(config%Nens),config%endIndex(config%Nens),config%locsize(config%Nens),itemp(config%Nens+1))


  ! indices for distributed state vector
  itemp = [ (1 + ( (i-1)*nl)/config%Nens, i=1,config%Nens+1)]
  !write(6,*) 'itemp ',itemp

  config%startIndex = itemp(1:config%Nens)
  config%endIndex = itemp(2:config%Nens+1)-1

  ! local size for assimilation
  config%locsize = config%endIndex - config%startIndex+1

  deallocate(itemp)

  if (present(partition)) then
    allocate(config%partition(size(partition)))
    config%partition = partition
  end if

  if (present(gridx)) then
    allocate(config%gridx(size(gridx)))
    config%gridx = gridx
  end if

  if (present(gridy)) then
    allocate(config%gridy(size(gridy)))
    config%gridy = gridy
  end if

  if (present(gridz)) then
    allocate(config%gridz(size(gridz)))
    config%gridz = gridz
  end if

  if (present(gridt)) then
    allocate(config%gridt(size(gridt)))
    config%gridt = gridt
  end if

 end subroutine oak_domain


 subroutine oak_done(config)
  use mpi
  use ndgrid
  implicit none
  type(oakconfig), intent(inout) :: config
  integer :: ierr, v

  if (.not.config%model_uses_mpi) then
    call mpi_finalize(ierr)
  end if

  deallocate(config%startIndex,config%endIndex,config%locsize)

  if (allocated(config%gridx)) deallocate(config%gridx)
  if (allocated(config%gridy)) deallocate(config%gridy)
  if (allocated(config%gridz)) deallocate(config%gridz)
  if (allocated(config%gridt)) deallocate(config%gridt)

  do v=1,ModML%nvar
    ! deallocate the model grid structure ModelGrid(v)
    call donegrid(ModelGrid(v))
  end do

  deallocate(ModelGrid)

  call MemoryLayoutDone(ModML)

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

  call mpi_barrier(comm, ierr)
 end subroutine oak_split_comm


 !-------------------------------------------------------------

 subroutine oak_obsoper(config,x,Hx)
  use mpi
  implicit none
  type(oakconfig), intent(inout) :: config
  real, intent(in) :: x(:)
  real, intent(out) :: Hx(:)
  integer :: rank, m = 1, ierr
  ! use communicator for each ensemble member


  if (config%model_uses_mpi) then
    call mpi_comm_rank(config%comm_ensmember, rank, ierr)
    if (rank == 0) then
      Hx = x(1) 
    end if

    call mpi_bcast(Hx, m, DEFAULT_REAL, 0, config%comm_ensmember, ierr)
  else
    Hx = x(1)
  end if
 end subroutine oak_obsoper


 !-------------------------------------------------------------

 subroutine oak_load_obs(config,ntime,ObsML,yo,invsqrtR,H,Hshift, &
      obsX, obsY, obsZ, obsT)
  use matoper
  use initfile
  implicit none
  type(oakconfig), intent(inout) :: config
  integer, intent(in) :: ntime
  type(MemLayout), intent(out) :: ObsML
  real, intent(out), pointer :: yo(:), invsqrtR(:), Hshift(:)
  real, intent(out), pointer, dimension(:) :: obsX, obsY, obsZ, obsT

  type(SparseMatrix), intent(out) :: H

  character(len=256)             :: prefix
  real(8) :: mjd0
  integer :: m

  call fmtIndex('Obs',ntime,'.',prefix)
  call MemoryLayout(prefix,ObsML)

  m = ObsML%effsize

  allocate(yo(m),invsqrtR(m),Hshift(m))
  call loadObs(ntime,ObsML,yo,invsqrtR)    

  call loadObservationOper(ntime,ObsML,H,Hshift,invsqrtR)

!  call loadVector(trim(prefix)//'gridX',ObsML,obsX)
!  call loadVector(trim(prefix)//'gridY',ObsML,obsY)
!  call loadVector(trim(prefix)//'gridZ',ObsML,obsZ)

  if (presentInitValue(initfname,trim(prefix)//'gridT')) then
!    call loadVector(trim(prefix)//'gridT',ObsML,obsT)
  end if

 end subroutine oak_load_obs



 !-------------------------------------------------------------
 !
 ! x is the state vector with all variables concatenated and masked points 
 ! removed

 subroutine oak_assim(config,ntime,x)
  use mpi
  use matoper
  use sangoma_ensemble_analysis
  implicit none
  type(oakconfig), intent(inout) :: config
  integer, intent(in) :: ntime
  real, intent(inout) :: x(:)
  integer :: ierr, i, obs_dim, r, myrank
  integer, dimension(config%Nens) :: recvcounts, sdispls, rdispls

  real, allocatable :: Hx(:),HEf(:,:), d(:), meanHx(:), HSf(:), diagR(:), Ef(:,:), Ea(:,:)
  real, pointer :: yo(:), invsqrtR(:), Hshift(:)
  real :: scale
  character(len=64) :: method = 'ETKF'
  type(MemLayout) :: ObsML
  type(SparseMatrix) :: H
  real(8) :: mjd0
  real, pointer, dimension(:) :: obsX, obsY, obsZ, obsT


  ! fixme ntime
  ! check if observations are available
  call loadObsTime(ntime,mjd0,ierr)
  if (ierr /= 0) then
    return
  end if

  ! load observations
  call oak_load_obs(config,ntime,ObsML,yo,invsqrtR,H,Hshift, &
       obsX, obsY, obsZ, obsT)

  obs_dim = size(yo)
  allocate(Hx(obs_dim),HEf(obs_dim,config%Nens),d(obs_dim),HSf(obs_dim), &
       meanHx(obs_dim),diagR(obs_dim))

  ! extract observations
  call oak_obsoper(config,x,Hx)

  call mpi_allgather(Hx, obs_dim, DEFAULT_REAL, HEf,  &
       obs_dim, DEFAULT_REAL, config%comm_da, ierr)

  !  write(6,*) 'Hx',Hx

  !  write(6,*) 'x',x

  call mpi_comm_rank(config%comm_da, myrank, ierr)

  ! collect the local part of the state vector

  recvcounts = config%locsize(myrank+1)

  ! allocate(xloc(sum(recvcounts)))
  !allocate(Ef(config%locsize(myrank+1),config%Nens))

  allocate( &
       Ea(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel,ErrorSpaceDim), &
       Ef(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel,ErrorSpaceDim))



  !  write(6,*) 'model ',myrank,'counts ',config%locsize,recvcounts

  ! The type signature associated with sendcount[j], sendtype at process i must be equal to the type signature associated with recvcount[i], recvtype at process j. This implies that the amount of data sent must be equal to the amount of data received, pairwise between every pair of processes.

  sdispls = config%startIndex-1
  rdispls(1) = 0
  do i=2,config%Nens
    rdispls(i) = rdispls(i-1) + recvcounts(i-1)
  end do

  !write(6,*) 'model ',myrank,'disp ',sdispls,rdispls

  call mpi_alltoallv(x, config%locsize, sdispls, DEFAULT_REAL, &
       Ef, recvcounts, rdispls, DEFAULT_REAL, config%comm_da, ierr)




  write(6,*) 'Ef1',Ef(:,1)
  write(6,*) 'Ef2',Ef(:,2)

!   call assim(ntime,Ef,Ea)


  !write(6,*) 'model ',myrank,'has loc ',Ef

  ! analysis

  if (config%schemetype == 1) then
    ! global

  else
!    call local_ensemble_analysis(Ef,HEf,yo,diagR,partition,selectObs,method,Ea)
  end if


  call mpi_alltoallv(Ef, recvcounts, rdispls, DEFAULT_REAL, &
       x, config%locsize, sdispls, DEFAULT_REAL, config%comm_da, ierr)

  !write(6,*) 'model ',myrank,'has global',x

  !-------------------------------------------------------------
 contains
  subroutine selectObs(i,w) 
   implicit none
   integer, intent(in) :: i
   real, intent(out) :: w(:)


  end subroutine selectObs
 end subroutine oak_assim



end module oak
