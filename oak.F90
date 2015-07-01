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

   ! weights for particle filter
   real, allocatable :: weightf(:), weighta(:)

   ! time index of the next observation
   integer :: obsntime = 1   

   ! domain index 
   integer, allocatable :: dom(:)

   ! number of domains
   integer :: ndom

 end type oakconfig

contains

 subroutine oak_init(config,fname,comm,comm_ensmember_)
  use mpi
  use initfile
  use parall, only: parallInit
  use initfile
  implicit none

  type(oakconfig), intent(out) :: config
  character(len=maxLen) :: fname
  integer, intent(in), optional :: comm
  integer, intent(out), optional :: comm_ensmember_

  integer :: nprocs, nprocs_ensmember, vmax
  integer :: ierr

  config%initfname = ''

  initfname = fname
  config%initfname = fname
  call getInitValue(initfname,'ErrorSpace.dimension',config%Nens,default=0)

  call parallInit(communicator=config%comm_all)
  call init(config%initfname)
  allocate(config%dom(ModML%effsize))

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
!  write(6,*) 'OAK is initialized'
!  write(6,*) 'Total number of processes',nprocs
!  write(6,*) 'Processes per ensemble member',nprocs_ensmember

 end subroutine oak_init

 !--------------------------------------------------------------------


 subroutine oak_domain_decomposition(config,dom)
  implicit none
  type(oakconfig), intent(inout) :: config
  integer, intent(in) :: dom(:)
  
  config%dom = dom
  ! number of domains
  config%ndom = maxval(config%dom)
  
 end subroutine oak_domain_decomposition

!-------------------------------------------------------------------------------

 subroutine oak_done(config)
  use mpi
  use ndgrid
  implicit none
  type(oakconfig), intent(inout) :: config
  integer :: ierr, v

  if (.not.config%model_uses_mpi) then
    call mpi_finalize(ierr)
  end if


  if (allocated(config%dom)) deallocate(config%dom)

  do v=1,ModML%nvar
    ! deallocate the model grid structure ModelGrid(v)
    call donegrid(ModelGrid(v))
  end do

  deallocate(ModelGrid)

  if (allocated(config%weightf)) deallocate(config%weightf)
  if (allocated(config%weighta)) deallocate(config%weighta)
!  deallocate(config%weightf,config%weighta)

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


  rank = procnum - 1

  do source = 0,p-1
    do dest = 0,p-1
      source_idom = mod(source,config%ndom) + 1
      ! ensemble member index (integer division)
      source_ensmember = source/config%ndom + 1


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
      dest_idom = mod(dest,config%ndom) + 1
      ! ensemble member index (integer division)
      dest_ensmember = dest/config%ndom + 1


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
      source_idom = mod(source,config%ndom) + 1;
      ! ensemble member index (integer division)
      source_ensmember = source/config%ndom + 1

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
      source_idom = mod(source,config%ndom) + 1;
      ! ensemble member index (integer division)
      source_ensmember = source/config%ndom + 1

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

 subroutine oak_assim(config,time,x)
  use mpi
  use parall
  use matoper
  use initfile
  implicit none
  type(oakconfig), intent(inout) :: config
  real(8), intent(in) :: time
  real, intent(inout) :: x(:)
  integer :: ierr, i, obs_dim, r, myrank
  integer, dimension(config%Nens) :: recvcounts, sdispls, rdispls

  real, allocatable :: Hx(:),HEf(:,:), d(:), meanHx(:), HSf(:), diagR(:), Ef(:,:), Ea(:,:)
  real, pointer :: yo(:), invsqrtR(:), Hshift(:)
  real :: scale
  character(len=64) :: method = 'ETKF'
  character(len=256)             :: infix
  type(MemLayout) :: ObsML
  type(SparseMatrix) :: H
  ! force double precision of time
  real(8) :: obstime, obstime_next, model_dt
  real, pointer, dimension(:) :: obsX, obsY, obsZ, obsT
  integer :: ntime, obsVec, dt_obs


  if (schemetype == EWPFScheme) then
    ! time of the next observation
    call loadObsTime(config%obsntime,obstime,ierr)    
    call getInitValue(initfname,'Config.dt',model_dt)

    allocate( &
         Ea(ModML%effsize,ErrorSpaceDim), &
         Ef(ModML%effsize,ErrorSpaceDim))

    call oak_gather_master(config,x,Ef)

    if (procnum == 1) then
      if (time < obstime) then
        ! proposal step

        ntime = time / model_dt
        call loadObsTime(config%obsntime+1,obstime_next,ierr)    

        if (ierr == 0) then
          obsVec = obstime_next / model_dt
          dt_obs = (obstime_next - obstime) / model_dt

          call fmtIndex('',config%obsntime+1,'.',infix)
          call MemoryLayout('Obs'//trim(infix),ObsML,.true.)
          allocate(yo(ObsML%effsize),invsqrtR(ObsML%effsize),Hshift(ObsML%effsize))

          call loadObs(config%obsntime+1,ObsML,yo,invsqrtR)    
          call loadObservationOper(config%obsntime+1,ObsML,H,Hshift,invsqrtR)

          call ewpf_proposal_step(ntime,obsVec,dt_obs,Ef,Ea,config%weightf,yo,invsqrtR,H)
          deallocate(yo,invsqrtR,Hshift)
        end if
      else
        ! analysis step
        call assim(config%obsntime,Ef,Ea, & 
             weightf=config%weightf,weighta=config%weighta)
        config%weightf = config%weighta
      end if
    end if

    call oak_spread_master(config,Ea,x)   

    deallocate(Ea,Ef)
  else

    ! fixme ntime
    ! check if observations are available
    call loadObsTime(config%obsntime,obstime,ierr)
    if (ierr /= 0) then
      return
    end if

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
      
      call assim(config%obsntime,Ef,Ea)
      
      call oak_spread_members(config,Ea,x)
    else
    ! collect at master
      allocate( &
           Ea(ModML%effsize,ErrorSpaceDim), &
           Ef(ModML%effsize,ErrorSpaceDim))
      
      call oak_gather_master(config,x,Ef)
      Ea = Ef
      if (procnum == 1) then
        write(6,*) 'model ',procnum,'has loc ',shape(Ef)
        ! weights are not used unless  schemetype == EWPFScheme
        call assim(config%obsntime,Ef,Ea,weightf=config%weightf,weighta=config%weighta)
        ! for testing
        Ea = Ef
      end if
      call oak_spread_master(config,Ea,x)   
    end if


    config%obsntime = config%obsntime +1
    deallocate(Ea,Ef)
  end if
  !-------------------------------------------------------------
 contains
  subroutine selectObs(i,w) 
   implicit none
   integer, intent(in) :: i
   real, intent(out) :: w(:)


  end subroutine selectObs
 end subroutine oak_assim



end module oak
