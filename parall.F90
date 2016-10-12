!
!  OAK, Ocean Assimilation Kit
!  Copyright(c) 2002-2012 Alexander Barth and Luc Vandenblucke
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

module parall
# ifdef MPI  
  use mpi
# endif
# ifdef PVM
  include 'fpvm3.h'
# endif

#ifdef ASSIM_PARALLEL
 ! procnum: process number (one-based) (1 <= procnum <= nbprocs)
 ! nbprocs: number of processes

 integer, save :: procnum=1, nbprocs=1
 integer, save, allocatable :: procid(:), procSpeed(:), cumulProcSpeed(:)

 ! communicator
 integer :: comm

 ! indices for parallelisation (zones)
 ! vector of nbprocs integer
 integer, allocatable :: startZIndex(:), endZIndex(:)
#endif



contains

!_______________________________________________________


 subroutine parallInit(num,nb,speed,communicator)
  implicit none
  integer, intent(in), optional :: num,nb
  integer, intent(in), optional :: speed
  integer, intent(in), optional :: communicator

  integer :: i,inum,istat,info,ierr
  logical :: flag

#ifdef MPI  
  comm = mpi_comm_world
  if (present(communicator)) comm = communicator

  call mpi_initialized(flag, ierr)
  if (.not.flag) call mpi_init(ierr)

  call mpi_comm_rank(comm, procnum, ierr)
  call mpi_comm_size(comm, nbprocs, ierr)

  procnum = procnum+1

  allocate(procid(nbprocs),procSpeed(nbprocs),cumulProcSpeed(nbprocs+1))

  procSpeed = 1

  do i=1,nbprocs
    procid(i) = i-1
  end do
#endif


#ifdef PVM


  call pvmfjoingroup('worker', inum )
  call pvmfbarrier('worker',nbprocs,info)

  allocate(procid(nbprocs),procSpeed(nbprocs),cumulProcSpeed(nbprocs+1))

  do i=1,nbprocs
    if (i.eq.procnum) then
      ! send to all workers my models tid
      call pvmfmytid(procid(i))

      if (present(speed)) then
        procSpeed(i) = speed
      else
        procSpeed(i) = 1
      end if

#     ifdef DEBUG
      write(stdout,*) 'start sending broadcast procid ',i,procid(i),procSpeed(i)
      call flush(stdout,istat)
#     endif

      call pvmfinitsend(pvmdefault,info)
      call pvmfpack(integer4,procid(i), 1, 1, info )
      call pvmfpack(integer4,procSpeed(i), 1, 1, info )
      call pvmfbcast('worker',5555 + i, info )

#     ifdef DEBUG
      write(stdout,*) 'end sending broadcast procid '
      call flush(stdout,istat)
#     endif
    else
      ! recev the the models tid
#     ifdef DEBUG
      write(stdout,*) 'start reciving broadcast procid '
      call flush(stdout,istat)
#     endif

      call pvmfrecv(-1, 5555+i, info)
      call pvmfunpack(integer4,procid(i), 1, 1, info)
      call pvmfunpack(integer4,procSpeed(i), 1, 1, info )

#     ifdef DEBUG
      write(stdout,*) 'end recieving broadcast procid ',i,procid(i),procSpeed(i)
      call flush(stdout,istat)
#     endif 
    end if
  end do

#endif

#ifdef ASSIM_PARALLEL

! cumulative sum of the speed of each processor
! used for splitting a do-loop acording to the speed of each processor

  cumulProcSpeed(1) = 0  
  do i=1,nbprocs
    cumulProcSpeed(i+1) = cumulProcSpeed(i) + procSpeed(i)
  end do  

# ifdef DEBUG
  write(stdout,*) 'all procid ',procid(:),procSpeed(:),cumulProcSpeed(:)
  call flush(stdout,istat)
# endif

#endif

 end subroutine parallInit
            
!_______________________________________________________
!
! if nzones iterations must be carried out, then split the
! number of iterations among nbprocs processes taking the speed
! into accound
! A process i, must do the iterations 
! from startZIndex(i) to endZIndex(i)

! startZIndex and endZIndex are global variables

subroutine parallPartion(nzones)
  implicit none

!  integer, intent(in) :: zoneSize(:)
!  integer, intent(in) :: zoneIndex(:)
  integer, intent(in) :: nzones

#ifdef ASSIM_PARALLEL

  allocate(startZIndex(nbprocs),endZIndex(nbprocs))
  startZIndex =(nzones*cumulProcSpeed(1:nbprocs))/cumulProcSpeed(nbprocs+1) + 1
  endZIndex =  (nzones*cumulProcSpeed(2:nbprocs+1))/cumulProcSpeed(nbprocs+1)

# ifdef DEBUG
  write(stdout,*) 'partitioning start',startZIndex
  write(stdout,*) 'partitioning end  ',endZIndex
# endif

#endif

end subroutine


 !_______________________________________________________

!
! change of interface of parallSyncronise
! previously the first 2 arguments where scalars of the start and end indices of x
! now they are vectors

 subroutine parallSyncronise(i1,i2,x,tag)
  implicit none

  integer, intent(in) :: i1(:),i2(:)
  real, intent(inout) :: x(:)
  integer, intent(in) :: tag


  integer :: otheri1,otheri2,j,info,ierr,istat
  integer, allocatable :: rcount(:),rdispls(:)
  integer :: myi1,myi2


#ifdef MPI  

  allocate(rcount(nbprocs),rdispls(nbprocs))

  ! revieve count and displacements
  rcount = i2 - i1 + 1
  rdispls = i1-1

# ifdef DEBUG
  write(stdout,*) 'rdispls ',rdispls
  write(stdout,*) 'rcount ',rcount
# endif

! only master get the complete x
  call mpi_gatherv(x(i1(procnum):i2(procnum)),rcount(procnum),DEFAULT_REAL,x,rcount, &
     rdispls,DEFAULT_REAL,0, mpi_comm_world, ierr)


  deallocate(rcount,rdispls)
#endif


#ifdef PVM

  myi1 = i1(procnum)
  myi2 = i2(procnum)

  do j=1,nbprocs
    if (j.eq.procnum) then
    ! my turn to send
#     ifdef DEBUG
      write(stdout,*) 'proc ',j,'start broadcast. checksum ',sum(x(myi1:myi2)),myi1,myi2
      call flush(stdout,istat)     
#     endif

      call pvmfinitsend (pvminplace,info)
      call pvmfpack(integer4,myi1,1,1,info)
      call pvmfpack(integer4,myi2,1,1,info)
      call pvmfpack(real4,x(myi1:myi2),myi2-myi1+1,1,info)
      call pvmfbcast('worker', tag + j, info )

#     ifdef DEBUG
      write(stdout,*) 'proc ',j,'finish broadcast'
      call flush(stdout,istat)     
#     endif
    else
      ! reciev from j

#     ifdef DEBUG
      write(stdout,*) 'proc ',procnum,'start recieving from ',j
      call flush(stdout,istat)     
#     endif

      call pvmfrecv(-1,tag+j, info)
      call pvmfunpack(integer4,otheri1,1,1,info)
      call pvmfunpack(integer4,otheri2,1,1,info)
      call pvmfunpack(real4,x(otheri1:otheri2),otheri2-otheri1+1,1,info)

#     ifdef DEBUG
      write(stdout,*) 'proc ',procnum,'finish recieving from ',j,' checksum ',sum(x(otheri1:otheri2)),otheri1,otheri2
      call flush(stdout,istat)     
#     endif
    end if
  end do

#     ifdef DEBUG
      write(stdout,*) ' checksum ',sum(x)
      call flush(stdout,istat)     
#     endif

#endif

 end subroutine

 !_______________________________________________________


 subroutine parallSyncronise2(myi1,myi2,myx,allx,tag,whoRecieves)
  implicit none

  integer, intent(in) :: myi1,myi2
  real, intent(in) :: myx(:)
  real, intent(out) :: allx(:)
  integer, intent(in) :: tag
  integer, intent(in), optional :: whoRecieves(:)


  integer :: otheri1,otheri2,j,info,k,istat

#ifdef PVM
  write(stdout,*) 'parallSyncronise2'

  if (present(whoRecieves)) then
    do j=1,nbprocs
      ! my turn to send
      if (j.eq.procnum) then
        ! send to all expecpt to my self
        do k=1,size(whoRecieves)
          if (j.ne.procnum) then

            call pvmfinitsend (pvminplace,info)
            call pvmfpack(integer4,myi1,1,1,info)
            call pvmfpack(integer4,myi2,1,1,info)
            call pvmfpack(real4,myx,myi2-myi1+1,1,info)
            call pvmfsend(procid(whoRecieves(k)),tag+k,info)  
          end if
        end do
      else
        ! reciev from j
        if (any(whoRecieves.eq.procnum)) then
          CALL PVMFRECV(procid(j),tag+procnum, info)
          call pvmfunpack(integer4,otheri1,1,1,info)
          call pvmfunpack(integer4,otheri2,1,1,info)
          call pvmfunpack(real4,allx(otheri1:otheri2),otheri2-otheri1+1,1,info)
        end if
      end if
    end do

  else

    ! syncronise all

    do j=1,nbprocs
      if (j.eq.procnum) then
        ! my turn to send
#     ifdef DEBUG
        write(stdout,*) 'proc ',j,'start broadcast'
        call flush(stdout,istat)     
#     endif

        call pvmfinitsend (pvminplace,info)
        call pvmfpack(integer4,myi1,1,1,info)
        call pvmfpack(integer4,myi2,1,1,info)
        call pvmfpack(real4,myx,myi2-myi1+1,1,info)
        call pvmfbcast('worker', tag + j, info )

#     ifdef DEBUG
        write(stdout,*) 'proc ',j,'finish broadcast'
        call flush(stdout,istat)     
#     endif
      else
        ! reciev from j

#     ifdef DEBUG
        write(stdout,*) 'proc ',procnum,'start recieving from ',j
        call flush(stdout,istat)     
#     endif

        call pvmfrecv(-1,tag+j, info)
        call pvmfunpack(integer4,otheri1,1,1,info)
        call pvmfunpack(integer4,otheri2,1,1,info)
        call pvmfunpack(real4,allx(otheri1:otheri2),otheri2-otheri1+1,1,info)

#     ifdef DEBUG
        write(stdout,*) 'proc ',procnum,'finish recieving from ',j
        call flush(stdout,istat)     
#     endif
      end if
    end do
  end if
#else
  write(stderr,*) 'Error: not implemented',__FILE__,__LINE__
  allx = 0
#endif

 end subroutine parallSyncronise2

 !_______________________________________________________



 subroutine parallJoinParts(localSVsize,localxf,xf)
  implicit none
  integer, intent(in) :: localSVsize(:)
  real, intent(in) :: localxf(:)
  real, intent(out) :: xf(:)

#ifdef PVM

  integer :: i,info,ptr,istat

  if (procnum.eq.1) then
    xf(1:localSVsize(1)) = localxf
    ptr = localSVsize(1)+1

    ! collect local state vector form other grid
    do i=2,nbprocs
#     ifdef DEBUG
      write(stdout,*) 'parallJoinParts 251: waiting for localxf from',i
      call flush(stdout,istat)
#     endif

      call pvmfrecv(procid(i), 6555+i, info)
      call pvmfunpack(real4,xf(ptr),localSVsize(i), 1, info)

#     ifdef DEBUG
      write(stdout,*) 'parallJoinParts 259: get localxf from',i, &
                 'checksum ',sum(xf(ptr:ptr+localSVsize(i)))
      call flush(stdout,istat)
#     endif

      ptr = ptr+localSVsize(i)
    end do
  else
    ! send local state vector 
#   ifdef DEBUG
    write(stdout,*) 'parallJoinParts 269: send localxf of',i, &
             'checksum ',sum(localxf)
    call flush(stdout,istat)
#   endif

    call pvmfinitsend(pvminplace,info)
    call pvmfpack(real4,localxf,localSVsize(procnum),1,info)
    call pvmfsend(procid(1),6555+procnum,info)        

#   ifdef DEBUG
    write(stdout,*) 'parallJoinParts 279: end send localxf of',i, &
             'checksum ',sum(localxf)
    call flush(stdout,istat)
#   endif
  end if
#else
  write(stderr,*) 'Error: not implemented',__FILE__,__LINE__
  xf = 0
#endif

 end subroutine parallJoinParts

 !_______________________________________________________

 subroutine parallSplitParts(localSVsize,localxa,xa)
  implicit none
  integer, intent(in) :: localSVsize(:)
  real, intent(out) :: localxa(:)
  real, intent(in) :: xa(:)

  integer :: i,info,ptr,istat

#ifdef PVM

  if (procnum.eq.1) then
    localxa = xa(1:localSVsize(1))
    ptr = localSVsize(1)+1

    do i=2,nbprocs
#     ifdef DEBUG
      write(stdout,*) 'parallSplitParts 302: send localxa to',i, &
                 'checksum ',sum(xa(ptr:ptr+localSVsize(i)))
      call flush(stdout,istat)
#     endif
      call pvmfinitsend (pvminplace,info)
      call pvmfpack(real4,xa(ptr),localSVsize(i),1,info)
      call pvmfsend(procid(i),7555+i,info)        

#     ifdef DEBUG
      write(stdout,*) 'parallSplitParts 311: end send localxa to',i, &
                 'checksum ',sum(xa(ptr:ptr+localSVsize(i)))
      call flush(stdout,istat)
#     endif

      ptr = ptr+localSVsize(i)
    end do
  else
#   ifdef DEBUG
    write(stdout,*) 'parallSplitParts: wait localxa from',procid(1)
    call flush(stdout,istat)
#   endif

    CALL PVMFRECV(procid(1), 7555+procnum, info)
    CALL PVMFUNPACK(real4,localxa,localSVsize(procnum), 1, info)

#   ifdef DEBUG
    write(stdout,*) 'parallSplitParts: got localxa from',procid(1), &
             'checksum ',sum(localxa)
    call flush(stdout,istat) 
#   endif
  end if

#else
  write(stderr,*) 'Error: not implemented',__FILE__,__LINE__
  localxa = 0
#endif

 end subroutine parallSplitParts

 !_______________________________________________________
 !
 ! gather pieces of a distributed state vector at process 
 ! with procnum 1
 ! Input:
 !  xf: subset of state vector
 !  startIndexZones,endIndexZones: start and end indeces of each zone
 !    vector of Nzones integers where Nzones if the number of zones
 !  vectype: 0 for state vector (default) and 1 for vector of a size equal to the 
 !     number of zones
 ! Output:
 !  xt: gathered total vector

 subroutine parallGather(xf,xt,startIndexZones,endIndexZones, vectype)
  implicit none
  real, intent(in) :: xf(:)
  real, intent(out):: xt(:)
  integer, intent(in) :: startIndexZones(:),endIndexZones(:)
  integer, intent(in), optional :: vectype
  integer, allocatable :: rcount(:),rdispls(:)
  integer          :: ierr,k,j1,j2,baseIndex
  
  integer :: type
  real :: dummy(1)

  type = 0
  if (present(vectype)) type = vectype

#ifdef ASSIM_PARALLEL

  if (type == 0) then
    j1 = startIndexZones(startZIndex(procnum))
    j2 =   endIndexZones(  endZIndex(procnum))
  else
    j1 = startZIndex(procnum)
    j2 =   endZIndex(procnum)
  end if


  baseIndex = -j1+1

  allocate(rcount(nbprocs),rdispls(nbprocs))

  if (type == 0) then
    ! revieve count and displacements
    rcount = endIndexZones(endZIndex) - startIndexZones(startZIndex) + 1
    rdispls = startIndexZones(startZIndex) -1
  else
    ! revieve count and displacements
    rcount = endZIndex - startZIndex + 1
    rdispls = startZIndex -1
  end if

# ifdef DEBUG
!  write(stdout,*) 'rdispls ',rdispls,type
!  write(stdout,*) 'rcount ',rcount
# endif


  ! only master get the complete x

  if (procnum == 1) then
    call mpi_gatherv(xf,rcount(procnum),DEFAULT_REAL,xt,rcount,rdispls,DEFAULT_REAL,0, mpi_comm_world, ierr)
  else
    call mpi_gatherv(xf,rcount(procnum),DEFAULT_REAL,dummy,rcount,rdispls,DEFAULT_REAL,0, mpi_comm_world, ierr)
  end if
    
  deallocate(rcount,rdispls)

#else
  xt = xf
#endif
 end subroutine



 !_______________________________________________________
 !
 ! scatter pieces of a distributed state vector at process 
 ! with procnum 1
 ! Input:
 !  xt: scattered total vector
 !  startIndexZones,endIndexZones: start and end indeces of each zone
 !    vector of Nzones integers where Nzones if the number of zones
 ! Output:
 !  xf: subset of state vector

 subroutine parallScatter(xt,xp,startIndexZones,endIndexZones)
  implicit none
  real, intent(in):: xt(:)
  real, intent(out) :: xp(:)
  integer, intent(in) :: startIndexZones(:),endIndexZones(:)
  integer, allocatable :: rcount(:),rdispls(:)
  integer          :: ierr,k,j1,j2,baseIndex

  real :: dummy(1)

#ifdef ASSIM_PARALLEL
  j1 = startIndexZones(startZIndex(procnum))
  j2 =   endIndexZones(  endZIndex(procnum))

  baseIndex = -j1+1

  allocate(rcount(nbprocs),rdispls(nbprocs))

  ! revieve count and displacements
  rcount = endIndexZones(endZIndex) - startIndexZones(startZIndex) + 1
  rdispls = startIndexZones(startZIndex) -1

# ifdef DEBUG
!  write(stdout,*) 'rdispls ',rdispls
!  write(stdout,*) 'rcount ',rcount
# endif

  ! only master get the complete x

  if (procnum == 1) then
    call mpi_scatterv(xt,rcount,rdispls,DEFAULT_REAL,xp,rcount(procnum),DEFAULT_REAL,0, mpi_comm_world, ierr)
  else
    call mpi_scatterv(dummy,rcount,rdispls,DEFAULT_REAL,xp,rcount(procnum),DEFAULT_REAL,0, mpi_comm_world, ierr)
  end if
    
  deallocate(rcount,rdispls)

#else
  xp = xt
#endif
 end subroutine


!_______________________________________________________

! Finalize parallel job

 subroutine parallDone
#ifdef MPI  
  integer :: ierr

  call mpi_finalize(ierr)
#endif

#ifdef ASSIM_PARALLEL
  deallocate(procid,procSpeed,cumulProcSpeed)
  deallocate(startZIndex, endZIndex)
#endif

 end subroutine


end module parall
