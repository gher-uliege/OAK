! include the fortran preprocessor definitions
#include "ppdef.h"

!#define DEBUG
!#define MPI

module parall
# ifdef MPI  
  include 'mpif.h'
# endif
# ifdef PVM
  include 'fpvm3.h'
# endif

#ifdef ASSIM_PARALLEL
 integer, save :: procnum=1, nbprocs=1
 integer, save, allocatable :: procid(:), procSpeed(:), cumulProcSpeed(:)


! indices for parallelisation

! indices for zones
 integer, allocatable :: startZIndex(:), endZIndex(:)

! indices for state vector
!  integer, allocatable ::  startIndex(:),endIndex(:)

!  integer :: mystartIndex, myendIndex
#endif



contains


 subroutine parallInit(num,nb,speed)
  implicit none
  integer, intent(in), optional :: num,nb
  integer, intent(in), optional :: speed

  integer :: i,inum,istat,info,ierr


#ifdef MPI  
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, procnum, ierr)
  call mpi_comm_size(mpi_comm_world, nbprocs, ierr)

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

  allocate(startZIndex(nbprocs),endZIndex(nbprocs))
  startZIndex =(nzones*cumulprocspeed(1:nbprocs))/cumulprocspeed(nbprocs+1) + 1
  endZIndex =  (nzones*cumulprocspeed(2:nbprocs+1))/cumulprocspeed(nbprocs+1)

# ifdef DEBUG
  write(stdout,*) 'partitioning start',startZIndex
  write(stdout,*) 'partitioning end  ',endZIndex
# endif

end subroutine



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
  call mpi_gatherv(x(i1(procnum):i2(procnum)),rcount(procnum),mpi_real,x,rcount,rdispls,mpi_real,0, mpi_comm_world, ierr)


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


 subroutine parallSyncronise2(myi1,myi2,myx,allx,tag,whoRecieves)
  implicit none

  integer, intent(in) :: myi1,myi2
  real, intent(in) :: myx(:)
  real, intent(out) :: allx(:)
  integer, intent(in) :: tag
  integer, intent(in), optional :: whoRecieves(:)


  integer :: otheri1,otheri2,j,info,k,istat

  write(stdout,*) 'parallSyncronise2'

# ifdef PVM


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

#endif

 end subroutine parallSyncronise2



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

#endif

 end subroutine parallJoinParts


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

#endif


 end subroutine parallSplitParts

! Finalize parallel job

 subroutine parallDone
#ifdef MPI  
  integer :: ierr

  call mpi_finalize(ierr)
#endif

  deallocate(procid,procSpeed,cumulProcSpeed)
  deallocate(startZIndex, endZIndex)
 end subroutine


end module parall
