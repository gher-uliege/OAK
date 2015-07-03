!
!  OAK, Ocean Assimilation Kit
!  Copyright(c) 2015 Alexander Barth and Luc Vandenblucke
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!


! include the fortran preprocessor definitions
#include "ppdef.h"

! module for on-line coupling between model and assimilation routines

module oak
 use assimilation

 type oakconfig
   ! communicator model + assimilation
   integer :: comm_all

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

   real(8) :: starttime, obstime_previous
   
 end type oakconfig

contains

 subroutine oak_init(config,fname,comm,comm_ensmember_,starttime)
  use mpi
  use initfile
  use parall, only: parallInit
  use initfile
  implicit none

  type(oakconfig), intent(out)   :: config
  character(len=*)               :: fname
  integer, intent(in), optional  :: comm
  integer, intent(out), optional :: comm_ensmember_
  real(8), intent(in), optional :: starttime

  integer :: nprocs, nprocs_ensmember
  integer :: ierr

  initfname = fname
  config%starttime = 0
  if (present(starttime)) config%starttime = starttime
  config%obstime_previous = config%starttime

  call getInitValue(initfname,'ErrorSpace.dimension',config%Nens,default=0)

  call parallInit(communicator=config%comm_all)
  call init(initfname)
  allocate(config%dom(ModML%effsize))

  if (present(comm)) then
    ! model uses also MPI, split communicators
    call oak_split_comm(comm,config%Nens,comm_ensmember_,.true.)
    config%comm_all = comm
    config%model_uses_mpi = .true.

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


  allocate(config%weightf(config%Nens),config%weighta(config%Nens))
  config%weightf = 1./config%Nens
  
  ! print diagnostic information
  !  write(6,*) 'OAK is initialized'
  !  write(6,*) 'Total number of processes',nprocs
  !  write(6,*) 'Processes per ensemble member',nprocs_ensmember

 end subroutine oak_init

 !--------------------------------------------------------------------


 subroutine oak_domain_decomposition(config,dom)
  implicit none
  type(oakconfig), intent(inout) :: config
  integer, intent(in)            :: dom(:)

  config%dom = dom
  ! number of domains
  config%ndom = maxval(config%dom)

 end subroutine oak_domain_decomposition

 !-------------------------------------------------------------------------------

 subroutine oak_done(config)
  use mpi
  use ndgrid
  use assimilation
  implicit none
  type(oakconfig), intent(inout) :: config
  integer :: ierr, v

  if (.not.config%model_uses_mpi) then
    call mpi_finalize(ierr)
  end if


  if (allocated(config%dom)) deallocate(config%dom)
  if (allocated(config%weightf)) deallocate(config%weightf)
  if (allocated(config%weighta)) deallocate(config%weighta)

  call done()
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

 subroutine oak_gather_members(config,xdomain,E)
  use parall
  implicit none  
  type(oakconfig), intent(inout) :: config
  real, intent(in) :: xdomain(:)
  real, intent(out) :: E(:,:)

  integer :: i1,i2, subdomsize, needsize, p, rank, & 
       tag, ierr, source, dest, source_idom, source_ensmember

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

  integer :: i1,i2,i,ndom, subdomsize, needsize, p, rank, & 
       tag, ierr, dest, source, dest_idom, dest_ensmember

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
               DEFAULT_REAL, dest, tag, config%comm_all, ierr)
        end if

        if (rank == dest) then
          call mpi_recv(data, needsize, DEFAULT_REAL, source, & 
               tag, config%comm_all, status, ierr)
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

  integer :: ndom, p, rank, & 
       tag, ierr, source, source_idom, source_ensmember

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

 end subroutine oak_gather_master

 !-------------------------------------------------------------

 !-------------------------------------------------------------

 subroutine oak_spread_master(config,E,xdomain)
  use parall
  implicit none  
  type(oakconfig), intent(inout) :: config
  real, intent(in) :: E(:,:)
  real, intent(out) :: xdomain(:)

  integer :: p, rank, & 
       tag, ierr, source, source_idom, source_ensmember

  real,    allocatable :: data(:)
  integer :: status(mpi_status_size), datasize

  tag = 1

  ! number of processes
  call mpi_comm_size(config%comm_all, p, ierr)

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

 end subroutine oak_spread_master

 !-------------------------------------------------------------

 subroutine oak_perturb(config,x)
  use mpi
  use matoper
  use assimilation
  use parall
  implicit none
  type(oakconfig), intent(inout) :: config
  real, intent(inout) :: x(:)

  real, allocatable :: E(:,:), meanx(:)
!  write(6,*) 'x',x
!  x = x + reshape(randn(size(x),1),[size(x)])
!  write(6,*) 'x',x

    allocate( &
         E(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel, &
         ErrorSpaceDim), &
         meanx(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel))
    
    !write(6,*) 'orig',x

  if (schemetype == LocalScheme) then    
    call loadVectorSpace('ErrorSpace.init',ModMLParallel,E)
    call oak_spread_members(config,E,x)
  else
    if (procnum == 1) then
      !write(6,*) 'load vector space'
      call loadVectorSpace('ErrorSpace.init',ModML,E,meanx)
      E = sqrt(1.*ErrorSpaceDim - ASSIM_SCALING) * E + spread(meanx,2,ErrorSpaceDim)
      !write(6,*) 'ensemble',E
    end if
    
    call oak_spread_master(config,E,x)
  end if

!  write(6,*) 'load vector space',x
  
  
  deallocate(E,meanx)
 end subroutine oak_perturb
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
  integer :: ierr
  real, allocatable :: Ef(:,:), Ea(:,:)
  real, pointer :: yo(:), invsqrtR(:), Hshift(:)
  character(len=256)             :: infix
  type(MemLayout) :: ObsML
  type(SparseMatrix) :: H
  ! force double precision of time
  real(8) :: obstime, obstime_next, model_dt
  integer :: ntime, obsVec, dt_obs

  ! time of the next observation
  call loadObsTime(config%obsntime,obstime,ierr)    
  if (ierr /= 0) then
    ! no more observations are available
    write(6,*) 'no more obs',procnum
    return
  end if
  

  if (schemetype == EWPFScheme) then

    call getInitValue(initfname,'Config.dt',model_dt)

    allocate( &
         Ea(ModML%effsize,ErrorSpaceDim), &
         Ef(ModML%effsize,ErrorSpaceDim))

    !write(6,*) 'diff f',x
    call oak_gather_master(config,x,Ef)

    if (time < obstime) then
      ! proposal step
      if (procnum == 1) then
        write(6,*) 'proposal step',time
        dbg(Ef)

        ntime = nint(time / model_dt)
        call loadObsTime(config%obsntime+1,obstime_next,ierr)    

        if (ierr == 0) then
          obsVec = nint(obstime / model_dt)
          dt_obs = nint((obstime - config%obstime_previous) / model_dt)

          call fmtIndex('',config%obsntime+1,'.',infix)
          call MemoryLayout('Obs'//trim(infix),ObsML,.true.)
          allocate(yo(ObsML%effsize),invsqrtR(ObsML%effsize),Hshift(ObsML%effsize))

          call loadObs(config%obsntime+1,ObsML,yo,invsqrtR)    
          call loadObservationOper(config%obsntime+1,ObsML,H,Hshift,invsqrtR)

          !write(6,*) 'weightf ',procnum,config%weightf
          call ewpf_proposal_step(ntime,obsVec,dt_obs,Ef,Ea,config%weightf,yo,invsqrtR,H)
          !write(6,*) 'weightf ',procnum,config%weightf


          deallocate(yo,invsqrtR,Hshift)
        end if
      end if
    else
      ! analysis step
      if (procnum == 1) then
        write(6,*) 'analysis step',time

        call assim(config%obsntime,Ef,Ea, & 
             weightf=config%weightf,weighta=config%weighta)
        config%weightf = config%weighta

        dbg(Ea)
        dbg(config%weightf)

      end if

      config%obsntime = config%obsntime +1    
      config%obstime_previous = obstime
    end if

    call oak_spread_master(config,Ea,x)   

    !write(6,*) 'diff a',x
    dbg(x)
    deallocate(Ea,Ef)
  else

    if (schemetype == LocalScheme) then
      allocate( &
           Ea(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel,ErrorSpaceDim), &
           Ef(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel,ErrorSpaceDim))

      call oak_gather_members(config,x,Ef)

      !Ea = Ef
      ! analysis
      write(6,*) 'Ef ',Ef
      call assim(config%obsntime,Ef,Ea)
      write(6,*) 'Ea ',Ea
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
    config%obstime_previous = obstime
    deallocate(Ea,Ef)
  end if
 end subroutine oak_assim

end module oak
