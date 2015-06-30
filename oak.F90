! module for on-line coupling between model and assimilation routines

module oak
 use ndgrid
 use assimilation

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

   ! weights for particle filter
   real, allocatable :: weightf(:), weighta(:)

   integer :: obsntime = 1   

   ! domain index
   ! size state vector
   integer, allocatable :: dom(:)

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


  allocate(config%weightf(config%Nens),config%weighta(config%Nens))

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


  call parallInit(communicator=config%comm_all)
  call init(config%initfname)

  allocate(config%dom(ModML%effsize))

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

  deallocate(config%dom)

  do v=1,ModML%nvar
    ! deallocate the model grid structure ModelGrid(v)
    call donegrid(ModelGrid(v))
  end do

  deallocate(ModelGrid)
  deallocate(config%weightf,config%weighta)

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

 subroutine oak_gather_members(config,xdomain,E)
  use parall
  implicit none  
  type(oakconfig), intent(inout) :: config
  real, intent(in) :: xdomain(:)
  real, intent(out) :: E(:,:)

  integer :: i1,i2,i,ndom, subdomsize, needsize,k, p, rank, & 
       ensmember, tag, ierr, source, dest, source_idom, source_ensmember

  integer, allocatable :: ind(:), ind2(:)
  logical, allocatable :: need(:)
  real,    allocatable :: data(:)
  integer :: status(mpi_status_size)

  tag = 1

  ! number of processes
  call mpi_comm_size(config%comm_all, p, ierr)

  ! number of domains
  ndom = maxval(config%dom)

  rank = procnum - 1

  do source = 0,p-1
    do dest = 0,p-1
      source_idom = mod(source,ndom) + 1
      ! ensemble member index (integer division)
      source_ensmember = source/ndom + 1


      i1 = startIndexZones(startZIndex(dest+1))
      i2 =   endIndexZones(endZIndex(dest+1))

      subdomsize = count(config%dom == source_idom)
      allocate(ind(subdomsize),need(subdomsize))

      ind = pack(invZoneIndex,config%dom == source_idom)
      need = i1 <= ind .and. ind <= i2
      needsize = count(need)
      allocate(ind2(needsize),data(needsize))

      ind2 = pack(ind,need)

      if (count(need) > 0) then
        if (rank == source) then
          call mpi_send(pack(xdomain,need), needsize, DEFAULT_REAL, dest, & 
               tag, config%comm_all, ierr)
        end if

        if (rank == dest) then
          call mpi_recv(data, needsize, & 
               DEFAULT_REAL, source, tag, config%comm_all, status, ierr)

          E(ind2 - i1 + 1,source_ensmember) = data
        end if
      end if

      deallocate(ind,ind2,need,data)

    end do
  end do

 end subroutine oak_gather_members

 !-------------------------------------------------------------

 subroutine oak_spread_members(config,E,xdomain)
  use parall
  implicit none  
  type(oakconfig), intent(inout) :: config
  real, intent(out) :: xdomain(:)
  real, intent(in) :: E(:,:)

  integer :: i1,i2,i,ndom, subdomsize, needsize,k, p, rank, & 
       ensmember, tag, ierr, dest, source, dest_idom, dest_ensmember

  integer, allocatable :: ind(:), ind2(:)
  logical, allocatable :: need(:)
  real,    allocatable :: data(:)
  integer :: status(mpi_status_size)

  tag = 1

  ! number of processes
  call mpi_comm_size(config%comm_all, p, ierr)

  ! number of domains
  ndom = maxval(config%dom)

  rank = procnum - 1

  do dest = 0,p-1
    do source = 0,p-1
      dest_idom = mod(dest,ndom) + 1
      ! ensemble member index (integer division)
      dest_ensmember = dest/ndom + 1


      i1 = startIndexZones(startZIndex(source+1))
      i2 =   endIndexZones(endZIndex(source+1))

      subdomsize = count(config%dom == dest_idom)
      allocate(ind(subdomsize),need(subdomsize))

      ind = pack(invZoneIndex,config%dom == dest_idom)
      need = i1 <= ind .and. ind <= i2
      needsize = count(need)
      allocate(ind2(needsize),data(needsize))

      ind2 = pack(ind,need)

      if (count(need) > 0) then
        if (rank == source) then
          call mpi_send(E(ind2 - i1 + 1,dest_ensmember), needsize, & 
               DEFAULT_REAL, dest, tag, config%comm_all, status, ierr)
        end if

        if (rank == dest) then
          call mpi_recv(data, needsize, DEFAULT_REAL, source, & 
               tag, config%comm_all, ierr)
          xdomain(pack([(i,i=1,size(xdomain))],need)) = data
        end if

      end if

      deallocate(ind,ind2,need,data)

    end do
  end do

 end subroutine oak_spread_members

 !-------------------------------------------------------------

 subroutine oak_gather_master(config,xdomain,E)
  use parall
  implicit none  
  type(oakconfig), intent(inout) :: config
  real, intent(in) :: xdomain(:)
  real, intent(out) :: E(:,:)

  integer :: i1,i2,i,ndom, subdomsize, needsize,k, p, rank, & 
       ensmember, tag, ierr, source, dest, source_idom, source_ensmember

  integer, allocatable :: ind(:), ind2(:)
  logical, allocatable :: need(:)
  real,    allocatable :: data(:)
  integer :: status(mpi_status_size), datasize

  tag = 1

  ! number of processes
  call mpi_comm_size(config%comm_all, p, ierr)

  ! number of domains
  ndom = maxval(config%dom)

  rank = procnum - 1


  call mpi_send(xdomain, size(xdomain), DEFAULT_REAL, 0, & 
       tag, config%comm_all, ierr)

  ! gather all data
  if (rank == 0) then
    ! collect from subdomains to form perturmed state vector
  
    do source = 0,p-1
      source_idom = mod(source,ndom) + 1;
      ! ensemble member index (integer division)
      source_ensmember = source/ndom + 1

      datasize = count(config%dom == source_idom)
      allocate(data(datasize))
      call mpi_recv(data, size(data), & 
           DEFAULT_REAL, source, tag, config%comm_all, status, ierr)

      E(pack(invZoneIndex,config%dom == source_idom),source_ensmember) = data
      deallocate(data)
    end do
  end if

 end subroutine

 !-------------------------------------------------------------

 !-------------------------------------------------------------

 subroutine oak_spread_master(config,E,xdomain)
  use parall
  implicit none  
  type(oakconfig), intent(inout) :: config
  real, intent(in) :: E(:,:)
  real, intent(out) :: xdomain(:)

  integer :: i1,i2,i,ndom, subdomsize, needsize,k, p, rank, & 
       ensmember, tag, ierr, source, dest, source_idom, source_ensmember

  integer, allocatable :: ind(:), ind2(:)
  logical, allocatable :: need(:)
  real,    allocatable :: data(:)
  integer :: status(mpi_status_size), datasize

  tag = 1

  ! number of processes
  call mpi_comm_size(config%comm_all, p, ierr)

  ! number of domains
  ndom = maxval(config%dom)

  rank = procnum - 1


  ! gather all data
  if (rank == 0) then
    ! collect from subdomains to form perturmed state vector
  
    do source = 0,p-1
      source_idom = mod(source,ndom) + 1;
      ! ensemble member index (integer division)
      source_ensmember = source/ndom + 1

      datasize = count(config%dom == source_idom)
      allocate(data(datasize))

      data = E(pack(invZoneIndex,config%dom == source_idom),source_ensmember)

      call mpi_send(data, size(data), & 
           DEFAULT_REAL, source, tag, config%comm_all, ierr)

      deallocate(data)
    end do
  end if

  call mpi_recv(xdomain, size(xdomain), DEFAULT_REAL, 0, & 
       tag, config%comm_all, status, ierr)


 end subroutine

 !-------------------------------------------------------------


 !-------------------------------------------------------------
 !
 ! x is the state vector with all variables concatenated and masked points 
 ! removed

 subroutine oak_assim(config,ntime,x)
  use mpi
  use matoper
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

      write(6,*) 'invZoneIndex ',allocated(invZoneIndex)

  if (schemetype == LocalScheme) then
    allocate( &
         Ea(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel,ErrorSpaceDim), &
         Ef(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel,ErrorSpaceDim))

    call oak_gather_members(config,x,Ef)

    !  Write(6,*) 'model ',myrank,'counts ',config%locsize,recvcounts

    write(6,*) 'x',x
    write(6,*) 'Ef1',Ef(:,1)
    write(6,*) 'Ef2',Ef(:,2)

    Ea = Ef


    !write(6,*) 'model ',myrank,'has loc ',Ef
    
    ! analysis

    ! weights are not used unless  schemetype == EWPFScheme
    call assim(ntime,Ef,Ea,config%weightf,config%weighta)

    call oak_spread_members(config,Ea,x)
  else
    ! collect at master
    allocate( &
         Ea(ModML%effsize,ErrorSpaceDim), &
         Ef(ModML%effsize,ErrorSpaceDim))

    call oak_gather_master(config,x,Ef)
    Ea = Ef
    call oak_spread_master(config,Ea,x)   
  end if

  deallocate(Ea,Ef)
  !-------------------------------------------------------------
 contains
  subroutine selectObs(i,w) 
   implicit none
   integer, intent(in) :: i
   real, intent(out) :: w(:)


  end subroutine selectObs
 end subroutine oak_assim



end module oak
