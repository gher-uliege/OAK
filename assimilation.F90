!
!  OAK, Ocean Assimilation Kit
!  Copyright(c) 2002-2015 Alexander Barth and Luc Vandenblucke
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

#define DEBUG
#define PROFILE

module assimilation
 use ndgrid

 interface loadStateVector
   module procedure loadStateVector_byinitval, loadStateVector_byfilenames
 end interface

 interface saveStateVector 
   module procedure saveStateVector_byinitval, saveStateVector_byfilenames
 end interface

 interface loadVector
   module procedure loadVector_byinitval, loadVector_byfilenames
 end interface

 interface saveVector
   module procedure saveVector_byinitval, saveVector_byfilenames
 end interface

 interface MemoryLayout
   module procedure MemoryLayout_byinitval, MemoryLayout_byfilenames 
 end interface

 interface loadSparseMatrix
   module procedure loadSparseMatrix_byinitval, loadSparseMatrix_byfilenames
 end interface

 interface saveSparseMatrix
   module procedure saveSparseMatrix_byinitval, saveSparseMatrix_byfilenames
 end interface

 logical, parameter :: removeLandPoints = .true.
! logical, parameter :: removeLandPoints = .false.
 integer, parameter :: maxLen = 256
 
 character(len=maxLen) :: initfname, localInitfname, globalInitfname

 type(grid), allocatable :: ModelGrid(:)

 real, allocatable          :: maxCorrection(:)

 ! horizontal resolution of each variable (approximation)

 real, allocatable         :: hres(:)
 integer                   :: StateVectorSize, StateVectorSizeSea
! integer                   :: ErrorSpaceDim

 ! possible values of "runtype" are:
 ! FreeRun     = 0 = do nothing, i.e. a pure run of the model 
 ! FreeCompRun = 1 = still do not assimilate, but compare model to observations
 ! AssimRun    = 2 = assimilate observations

 integer, target :: runtype
 integer, parameter :: &
      FreeRun     = 0,    &
      FreeCompRun = 1,    &
      AssimRun    = 2

 ! How to propagate the model error in time ?

 integer :: moderrtype
 integer, parameter :: &
      ConstModErr    = 0, &
      FFModErr       = 1, &
      SpaceFFModErr  = 2, &    
      LagrModErr     = 3, &
      EnsembleModErr = 4

 ! Is a bias expected ?
 ! 0 = standard bias-blind analysis
 ! 1 = A fraction of the error (gamma) is a systematic error 
 !     and the rest (1-gamma) is random

 integer :: biastype
 integer, parameter :: &
      NoBias            = 0, &           ! (default)
      ErrorFractionBias = 1    

 integer :: Anamorphosistype
 integer, parameter :: &
      NoAnamorphosis = 0, &              ! (default)
      TSAnamorphosis = 1    


 integer :: metrictype
 integer, parameter :: &
      CartesianMetric  = 0, &           
      SphericalMetric  = 1, &                  ! (default)
      SphericalMetricApprox  = 2 


 ! value for the _FillValue attribute
 real :: FillValue = 9999.

 ! fortran unit of logfile
 ! the logfile contains simple diagnostics such as rmse with observations

 integer :: stdlog

 ! fortran unit of debugfile
 ! the debugfile contains debugging information such as 
 ! input/output, memory layout structure,... 

 integer :: stddebug

 logical :: rmLPObs = .true.
 !  logical :: rmLPObs = .false.

 ! bias of the model's forecast
 ! a priorie and a posteriorie to the observations

 real, pointer :: biasf(:), biasa(:)
 real :: biasgamma  

 !
 ! the structure "MemLayout" described the arrangement in memory of a vector
 ! composed by a certain number of variables. For example the state vector composed 
 ! by elevation, temperature, salinity,... or the observation vector composed by all 
 ! observation of a given time: e.g. temperature and salinity profile
 ! 
 ! nvar = number of variables
 ! totsize: size of the vector including land points
 ! totsizesea: size of the vector excluding land points
 ! removeLandPoints: true if the land points have to be removed
 !
 ! effsize = totsize if removeLandPoints = false.
 !         = totsizesea if removeLandPoints = true
 ! ndim: dimension of each variable
 ! varshape: the shape of each variable
 ! varsize: the size of each variable including land points
 ! varsizesea: the size of each variable excluding land points
 ! mask: 1 if sea point and 0 if land point
 ! startIndex,endIndex: start- and end-index of a variable in vector containing 
 !   also land points
 ! startIndexSea,endIndexSea: start- and end-index of a variable in vector containing 
 !   only sea points
 ! seaindex: vector of indexes transforming the space with masked points into
 !   space only with unmasked points
 ! invindex: vector of indexes transforming space only with unmasked points into 
 !   the space with masked points

 type MemLayout
   integer :: nvar, totsize, totsizesea, effsize
   integer, pointer :: ndim(:), varshape(:,:), varsize(:), varsizesea(:), &
        mask(:),startIndex(:),endIndex(:),startIndexSea(:),endIndexSea(:), &
        seaindex(:),invindex(:)

   ! names of individual variables
   character(len=maxLen), pointer :: varnames(:)

   logical :: removeLandPoints 
   logical :: permute

   ! in single CPU version those have sensible values
   logical :: distributed
   integer :: startIndexParallel, endIndexParallel
 end type MemLayout

! entire model state vector

 type(MemLayout), target ::  ModML

! model state vector for specific node

 type(MemLayout), target ::  ModMLParallel

#if defined(NEST) & defined(ASSIM)
 type(MemLayout), pointer :: LocalModML
#endif


! local or global assimilation ?

  integer :: schemetype
  integer, parameter :: &
      GlobalScheme  = 0, &           ! (default)
      LocalScheme   = 1, &    
      CLocalScheme   = 2,  &
      EWPFScheme   = 3

 ! type of localization
 ! 1: horizontal distance (default)
 ! 2: only gridZ distance
 ! 3: only gridT distance
  integer :: loctype

! partition for local assimilation

! zoneIndex: permutation vector for state vector such that all 
! zones are continuous elements


  integer, allocatable :: partition(:),zoneSize(:),zoneIndex(:),invZoneIndex(:)


! startIndex and endIndex for local zones

  integer, allocatable :: startIndexZones(:),endIndexZones(:)


! 
! variables read by callback function selectObservations
!

  real, allocatable :: obsGridX(:),obsGridY(:),obsGridZ(:),obsGridT(:)
  real, allocatable :: hCorrLengthToObs(:), hMaxCorrLengthToObs(:)

 !_______________________________________________________
 !

contains

 !_______________________________________________________
 !

 subroutine init(fname)
  implicit none
  character(len=*), intent(in) :: fname

  call localInit(fname)
  call globalInit(fname)

 end subroutine init

 !_______________________________________________________
 !
 !_______________________________________________________
 !

 subroutine globalInit(fname)
  use initfile
  use ufileformat
  use anamorphosis, only: initAnamorphosis
# ifdef ASSIM_PARALLEL
  use parall
# endif
  use matoper
  use covariance
  implicit none
  character(len=*), intent(in) :: fname
  integer                      :: v,vmax,n
  character(len=MaxFNameLength), pointer   :: filenamesX(:),filenamesY(:),filenamesZ(:),    &
                                              filenamesT(:)
  character(len=MaxFNameLength)            :: path
  real, pointer                :: maxCorr(:), tmp(:), tmp_hres(:)
  integer                      :: NZones, zi, istat
  
! for paritioning
  integer, allocatable :: tmpi(:)

  initfname = fname

  call getInitValue(initfname,'runtype',runtype,default=AssimRun)
  call getInitValue(initfname,'metrictype',metrictype,default=SphericalMetric)

  call getInitValue(initfname,'Config.FillValue',FillValue,default=9999.)

  if (runtype.eq.FreeRun) return

  call getInitValue(initfname,'moderrtype',moderrtype,default=ConstModErr)
  call getInitValue(initfname,'biastype',biastype,default=NoBias)
  call getInitValue(initfname,'schemetype',schemetype,default=GlobalScheme)
  call getInitValue(initfname,'loctype',loctype,default=1)
  call getInitValue(initfname,'anamorphosistype',anamorphosistype,default=NoAnamorphosis)

! # ifdef ASSIM_PARALLEL
!   if (schemetype /= LocalScheme) then
!     write(stderr,*) 'Error: for parallel version schemetype should be 1 (local assimilation)'
!     ERROR_STOP
!   end if
! # endif
  
  call initAnamorphosis(fname,stddebug)


  if (biastype.eq.ErrorFractionBias) then
    call getInitValue(initfname,'Bias.gamma',biasgamma)       
  end if

  call getInitValue(initfname,'Model.path',path,default='')

  ! Models Memory Layout 

  call MemoryLayout('Model.',ModML)
  vmax = ModML%nvar
  StateVectorSize = ModML%EndIndex(vmax)
  StateVectorSizeSea = ModML%totsizesea

  allocate( &
       ModelGrid(vmax), &
       hres(vmax), &
       biasf(StateVectorSizeSea), &
       biasa(StateVectorSizeSea))

!
! define model grid
!

  call getInitValue(initfname,'Model.gridX',filenamesX)  

  do v=1,vmax

    n = ModML%ndim(v)

    ! initialisze the model grid structure ModelGrid(v)
    call initgrid(ModelGrid(v),n,ModML%varshape(1:n,v), &
         ModML%Mask(ModML%StartIndex(v):ModML%EndIndex(v)).eq.0)

    ! set the coordinates of the model grid
    call setCoord(ModelGrid(v),1,trim(path)//filenamesX(v))

    if (n > 1) then
      if (.not.associated(filenamesY)) then
        call getInitValue(initfname,'Model.gridY',filenamesY)
      end if

      call setCoord(ModelGrid(v),2,trim(path)//filenamesY(v))

      if (n > 2) then
        if (.not.associated(filenamesZ)) then
          call getInitValue(initfname,'Model.gridZ',filenamesZ)
        end if

        call setCoord(ModelGrid(v),3,trim(path)//filenamesZ(v))

        if (n > 3) then
          if (.not.associated(filenamesT)) then
            call getInitValue(initfname,'Model.gridT',filenamesT)
          end if

          call setCoord(ModelGrid(v),4,trim(path)//filenamesT(v))
          deallocate(filenamesT)
          
          if (n > 4) then
            write(stderr,*) 'The dimension of variable ',trim(ModML%varnames(v)),' is ',n
            write(stderr,*) 'Error: Only 4-d grids are supported for now. '
            ERROR_STOP 
          end if
        end if
      end if
    end if

  end do

  ! typical horizontal resolution to prioritize model grids in genObsOperator
  ! in case of overlapping grids
  if (presentInitValue(initfname,'Model.hres')) then
    call getInitValue(initfname,'Model.hres',tmp_hres)
    hres = tmp_hres
    deallocate(tmp_hres)
  else
    ! grids should not overlap in this case
    hres = 0
  end if

  if (associated(filenamesX)) deallocate(filenamesX)
  if (associated(filenamesY)) deallocate(filenamesY)
  if (associated(filenamesZ)) deallocate(filenamesZ)
  if (associated(filenamesT)) deallocate(filenamesT)

!  call getInitValue(initfname,'ErrorSpace.dimension',ErrorSpaceDim,default=0)

  biasf = 0.


  ModMLParallel = ModML

! variables for local assimilation

  if (schemetype == LocalScheme .or. schemetype == CLocalScheme) then
    allocate(tmp(ModML%effsize),partition(ModML%effsize), &
         hCorrLengthToObs(ModML%effsize),hMaxCorrLengthToObs(ModML%effsize))
    call loadVector('Zones.partition',ModML,tmp)
    ! convertion real -> integer
    partition = tmp+.5
    deallocate(tmp)

    ! deal with gaps in partition
    allocate(tmpi(ModML%effsize))
    tmpi = partition
    call unique(tmpi,n,ind2 = partition)
    deallocate(tmpi)

    call loadVector('Zones.corrLength',ModML,hCorrLengthToObs)
    call loadVector('Zones.maxLength',ModML,hMaxCorrLengthToObs)

    allocate(zoneSize(maxval(partition)),zoneIndex(StateVectorSizeSea),invZoneIndex(StateVectorSizeSea))

    call initPartition(partition,zoneSize,zoneIndex)

    NZones = size(zoneSize)
    allocate(startIndexZones(NZones),endIndexZones(NZones))

    startIndexZones(1) = 1
    endIndexZones(1) = zoneSize(1)

    do zi=2,size(zoneSize)
      startIndexZones(zi) = endIndexZones(zi-1)+1
      endIndexZones(zi)   = endIndexZones(zi-1)+zoneSize(zi)
    end do

#   ifdef ASSIM_PARALLEL
    call parallPartion(NZones)
    ModMLParallel%distributed = .true.
    ModMLParallel%startIndexParallel = startIndexZones(startZIndex(procnum))
    ModMLParallel%endIndexParallel   =   endIndexZones(endZIndex(procnum))

#   ifdef DEBUG
    write(stddebug,*) 'indexes for me',procnum,ModMLParallel%startIndexParallel,ModMLParallel%endIndexParallel
    call flush(stddebug,istat)    
#   endif
#else
    ModMLParallel%distributed = .false.
    ModMLParallel%startIndexParallel = 1
    ModMLParallel%endIndexParallel   =  ModMLParallel%effsize
#   endif     
  else
    ! no permutation
    allocate(invZoneIndex(StateVectorSizeSea))
    invZoneIndex = [(zi,zi=1,StateVectorSizeSea)]
  end if

  ModMLParallel%permute = schemetype == LocalScheme .or. schemetype == CLocalScheme

  if (schemetype == CLocalScheme) then
     ModML%permute = .true.
  end if

! Maximum correction

  allocate(tmp(ModML%effsize))
  allocate(maxCorrection(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel));

  if (presentInitValue(initfname,'Model.maxCorrection')) then
    call getInitValue(initfname,'Model.maxCorrection',maxCorr)   
    do v=1,ModML%nvar
      tmp(ModML%startIndexSea(v):ModML%endIndexSea(v)) = maxCorr(v);
    end do
    deallocate(maxCorr)
  else
    tmp = huge(tmp(1))
  end if

  if (ModMLParallel%permute) call permute(zoneIndex,tmp,tmp)

  maxCorrection = tmp(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel);


  

  deallocate(tmp)

# ifdef DEBUG
  write(stddebug,*) '= Assimilation parameters ='
  write(stddebug,'(20A,I3)') 'run type',runtype
  write(stddebug,'(20A,I3)') 'metric',metrictype
  write(stddebug,'(20A,I3)') 'scheme',schemetype
  write(stddebug,'(20A,I3)') 'localization',loctype
  write(stddebug,'(20A,I3)') 'anamorphosis',anamorphosistype
# endif  


 end subroutine globalInit

 !_______________________________________________________
 !
 !_______________________________________________________
 !

 subroutine localInit(fname)
  use initfile
  use ufileformat
#ifdef ASSIM_PARALLEL
  use parall, only: procnum
#endif
  implicit none
  character(len=*), intent(in) :: fname
  character(len=maxLen)            :: str,postfix=''

  initfname = fname

#ifdef ASSIM_PARALLEL
  write(postfix,'(A,I5.5)') '-',procnum
#endif

! open log file with unit stdlog

  if (presentInitValue(initfname,'logfile')) then
    call getInitValue(initfname,'logfile',str)
    stdlog = 912391
    open(stdlog,file=trim(str)//postfix,status='unknown',position='append')
  else
    stdlog = stdout
  end if

  stddebug = 6

#ifdef DEBUG
! open debug file with unit stddebug

  if (presentInitValue(initfname,'debugfile')) then
    call getInitValue(initfname,'debugfile',str)
    stddebug = 92392
    open(stddebug,file=trim(str)//postfix,status='unknown',position='append')
  else
    stddebug = stdout
  end if
#endif
 
 end subroutine localInit

 !_______________________________________________________
 !

 subroutine done
  implicit none
  integer :: v

  do v=1,ModML%nvar
    ! initialisze the model grid structure ModelGrid(v)
    call donegrid(ModelGrid(v))
  end do

  call MemoryLayoutDone(ModML)

  deallocate(biasf,biasa,maxCorrection)
  deallocate(ModelGrid,hres)

   if (schemetype.eq.LocalScheme) then
     deallocate(partition)
     deallocate(hCorrLengthToObs)
     deallocate(hMaxCorrLengthToObs)
     deallocate(zoneSize)
     deallocate(zoneIndex)
     deallocate(invZoneIndex)
     deallocate(startIndexZones,endIndexZones)
   end if
  close(stdlog)
#ifdef DEBUG
  close(stddebug)
#endif
  end subroutine done

!
! in partition 
!

 subroutine initPartition(partition,zoneSize,zoneIndex)
 use ufileformat
 implicit none

 integer :: partition(StateVectorSizeSea)
 integer :: zoneSize(:)

 integer, dimension(size(zoneSize)) :: startIndex,cIndex

 integer :: zoneIndex(StateVectorSizeSea)

 integer :: j,NZones,zi
 integer :: istat


#ifdef DEBUG
 write(stddebug,*) 'Start the partitioning of the state vector'
 call flush(stddebug,istat)
#endif

 NZones = maxval(partition)

! size of each zone

  zoneSize = 0
  do j=1,StateVectorSizeSea
    zi = partition(j)

    if (zi < 1) then
      write(stderr,*) 'Error: all values in partition vector must be greater than one ',j,zi
      ERROR_STOP
    end if

    zoneSize(zi) = zoneSize(zi)+1
  end do

  if (any(zoneSize.eq.0)) then
    write(stderr,*) 'Error: all number in partition file have to be integer, starting with 1 and must follow without gap'
    ERROR_STOP
  end if

! start index of each zone in vector zoneIndex

  startIndex(1) = 1
  do zi=1,NZones-1
    startIndex(zi+1) = startIndex(zi)+zoneSize(zi)
  end do

  cIndex = startIndex

  do j=1,StateVectorSizeSea
    zi = partition(j)
    zoneIndex(cIndex(zi)) = j
    invZoneIndex(j) = cIndex(zi)

    cIndex(zi) = cIndex(zi)+1
  end do

#ifdef DEBUG   
   write(stddebug,*) 'End the partitioning of the state vector'
   call flush(stddebug,istat)
#endif
 
 end subroutine

 !_______________________________________________________
 !
subroutine fmtIndex(str1,index,str2,str)
 implicit none
 integer, intent(in) :: index
 character(len=*), intent(in) :: str1, str2
 character(len=*), intent(out) :: str
 character(len=128) :: tmp

 if (index < 999) then
   write(str,'(A,I3.3,A)') str1,index,str2
 else
   write(tmp,*) index

   str = trim(str1) // trim(adjustl(tmp)) // trim(str2)
 end if

end subroutine fmtIndex

 !_______________________________________________________
 !

 subroutine checkVar(var,name)
  use ufileformat
  implicit none
  
  real, intent(in) :: var(:)
  integer, save :: i = 0
  character(len=*) :: name
  character(len=256) :: filename

#ifdef ASSIM_PARALLEL
  write(filename,'("debugM",I3.3)') i
#else
  write(filename,'("debugS",I3.3)') i
#endif

  call usave(trim(filename)//name,var,9999.)
  write(6,*) 'write ',trim(filename)//name
  i = i+1

 end subroutine checkVar

 !_______________________________________________________
 !
 ! high level generic input and output routines
 !
 !
 !_______________________________________________________
 !
!#define ASSIM_LOAD_OPTIM
#ifdef ASSIM_LOAD_OPTIM

 subroutine loadVector_byfilenames(path,filenames,ML,vector)
  use initfile
  use ufileformat
  use matoper
  use parall
  implicit none

  character(len=*), intent(in) :: path
  character(len=*), intent(in) :: filenames(:)
  type(MemLayout),  intent(in) :: ML
  real, intent(out) :: vector(:)
  character(len=30) :: infoformat = '(A50,2E14.5)'

  real, allocatable :: x(:),xt(:)

  integer :: v,i,j0,j1
  integer :: istat
  real :: valex

  allocate(x(ML%totsize))

  if (size(filenames).ne.ML%nvar) then
    write(stderr,*) 'ERROR: Vector to load should be composed by ',ML%nvar, &
         ' variables, but ',size(filenames), ' filenames are found: '
    write(stderr,*) (trim(filenames(v))//' ',v=1,size(filenames))
    call flush(stderr,istat)
    ERROR_STOP
  end if


  if (procnum == 1 .or..not.ML%distributed) then 
    ! load variables

#ifdef DEBUG
    write(stddebug,'(A50,A14,A14)') 'loaded variable','min','max'
#endif

    do v=1,size(filenames)
      i = ML%ndim(v)
      j0 = ML%StartIndex(v)
      j1 = ML%EndIndex(v)

      !    call ureadfull(trim(path)//filenames(v),x(ML%StartIndex(v)),valex,check_numel=ML%varsize(v))

      call ureadfull_srepmat(trim(path)//filenames(v),x(j0),valex,      &
           check_numel = ML%varsize(v),                                                 &
           force_shape = ML%varshape(1:i,v))

#ifdef DEBUG
      write(stddebug,infoformat) trim(path)//trim(filenames(v)),       &
           minval(x(j0:j1),x(j0:j1) /= valex),  &
           maxval(x(j0:j1),x(j0:j1) /= valex)
      call flush(stddebug,istat)
#endif

    end do

!    write(stdout ,*) 'load vector shape(vector) ',procnum, shape(vector), lbound(vector), ubound(vector), 'new'

    ! remove masked points if requested
    allocate(xt(ML%effsize))

    if (ML%removeLandPoints) then
      xt = pack(x,ML%Mask.eq.1)
    else
      where (ML%Mask.eq.0) x=0
      xt = x
    end if

    ! permute state vector
    if (ML%permute) call permute(zoneIndex,xt,xt)
  end if

  ! distribute 
  if (ML%distributed) then
    call parallScatter(xt,vector,startIndexZones,endIndexZones)
  else
    vector = xt
  end if

# ifdef DEBUG
  write(stddebug, *) 
# endif

  deallocate(x)
 end subroutine loadVector_byfilenames

#else

 subroutine loadVector_byfilenames(path,filenames,ML,vector)
  use initfile
  use ufileformat
  use matoper
  use parall
  implicit none

  character(len=*), intent(in) :: path
  character(len=*), intent(in) :: filenames(:)
  type(MemLayout),  intent(in) :: ML
  real, intent(out) :: vector(:)
  character(len=30) :: infoformat = '(A50,2E14.5)'

  real, allocatable :: x(:)

  integer :: v,k,i,l,istat,j0,j1
  real :: valex

  allocate(x(ML%totsize))

  if (size(filenames).ne.ML%nvar) then
    write(stderr,*) 'ERROR: Vector to load should be composed by ',ML%nvar, &
         ' variables, but ',size(filenames), ' filenames are found: '
    write(stderr,*) (trim(filenames(v))//' ',v=1,size(filenames))
    call flush(stderr,istat)
    ERROR_STOP
  end if


#ifdef DEBUG
  write(stddebug,'(A50,A14,A14)') 'loaded variable','min','max'
#endif

  do v=1,size(filenames)
    i = ML%ndim(v)
    j0 = ML%StartIndex(v)
    j1 = ML%EndIndex(v)

!    call ureadfull(trim(path)//filenames(v),x(ML%StartIndex(v)),valex,check_numel=ML%varsize(v))

    call ureadfull_srepmat(trim(path)//filenames(v),x(j0),valex,      &
       check_numel = ML%varsize(v),                                                 &
       force_shape = ML%varshape(1:i,v))

#ifdef DEBUG
    write(stddebug,infoformat) trim(path)//trim(filenames(v)),       &
         minval(x(j0:j1),x(j0:j1) /= valex),  &
         maxval(x(j0:j1),x(j0:j1) /= valex)
    call flush(stddebug,istat)
#endif

  end do

# ifdef ASSIM_PARALLEL
  !  write(stdout,*) 'distributed', ML%distributed
  !  write(stdout,*) 'rmland', ML%removeLandPoints
  
!!$  if (ML%distributed) then
!!$    vector = 99999.
!!$
!!$    do l=ML%startIndexParallel,ML%endIndexParallel
!!$
!!$     ! select only points corresponding to sea values and
!!$     ! map index l to look only for data belonging to the zone
!!$
!!$      vector(l) = x(ML%invindex(zoneIndex(l)))
!!$    end do
!!$
!!$  else
!!$
!!$    do l=1,size(vector)
!!$      vector(l) = x(ML%invindex(l))
!!$    end do
!!$
!!$  end if
!!$
!!$  !write(stdout,*)  'i ',ML%startIndexParallel,ML%endIndexParallel,size(vector)
!!$
!!$  if (ML%permute .and. .not. ML%distributed) call permute(zoneIndex,vector,vector)    
!!$


!  write(stdout ,*) 'load vector shape(vector) ',procnum, shape(vector), lbound(vector), ubound(vector) 

  if (ML%permute) then
    do l=ML%startIndexParallel,ML%endIndexParallel
     ! select only points corresponding to sea values and
     ! map index l to look only for data belonging to the zone      

     ! FIX ME
      if (size(vector) == ML%effsize) then
        vector(l) = x(ML%invindex(zoneIndex(l)))
      else
        vector(l-ML%startIndexParallel+1) = x(ML%invindex(zoneIndex(l)))
      end if
    end do

  else

    write(stddebug,*) 'files ',filenames
    write(stddebug,*) 'shape ',procnum,shape(vector),lbound(vector),ubound(vector),shape(x), &
         ML%startIndexParallel,ML%endIndexParallel
    do l=ML%startIndexParallel,ML%endIndexParallel
     ! select only points corresponding to sea values
      vector(l) = x(ML%invindex(l))
    end do
  end if

# else  
  if (ML%removeLandPoints) then
    vector = pack(x,ML%Mask.eq.1)
  else
    where (ML%Mask.eq.0) x=0
    vector = x
  end if

  if (ML%permute) call permute(zoneIndex,vector,vector)

# endif

# ifdef DEBUG
  write(stddebug, *) 
# endif

  deallocate(x)
 end subroutine loadVector_byfilenames

#endif
 !_______________________________________________________
 !

 subroutine loadVector_byinitval(str,ML,vector)
  use initfile
  use ufileformat
  character(len=*), intent(in) :: str
  type(MemLayout),  intent(in) :: ML
  real, intent(out) :: vector(:)

  character(len=maxLen), pointer :: filenames(:)
  character(len=maxLen) :: path
  real, pointer :: constVal(:)
  integer :: i

# ifdef DEBUG
  write(stddebug,'("== load Vector (",A,") ==")') str
# endif

  call getInitValue(initfname,str(1:index(str,'.',.true.))//'path',path,default='')

  if (presentInitValue(initfname,str)) then
    call getInitValue(initfname,str,filenames)

    call loadVector_byfilenames(path,filenames,ML,vector)

    deallocate(filenames)
  else
    call getInitValue(initfname,trim(str)//'.const',constVal)
    if (size(constVal).ne.ML%nvar) then
      write(stderr,*) 'ERROR: Vector "',trim(str)//'.const','" should be composed by ',ML%nvar, &
         ' variables, but ',size(constVal), ' constans are found: ',constVal(:)
      call flush(stderr,istat)
      stop
    end if

    do i=1,ML%nvar
       vector(ML%startIndexSea(i):ML%endIndexSea(i)) = constVal(i)
    end do    

    deallocate(constVal)
  end if



 end subroutine loadVector_byinitval

 !_______________________________________________________
 !


 !_______________________________________________________
 !
 ! where mask is false, the corresponding value will be 
 ! a exclusion point
 !_______________________________________________________
 !

 subroutine saveVector_byfilenames(path,filenames,ML,vector,mask)
  use initfile
  use ufileformat
  use parall
  use matoper
  
  implicit none
  character(len=*), intent(in) :: path
  character(len=*), intent(in) :: filenames(:)
  type(MemLayout),  intent(in) :: ML
  real, intent(in), target :: vector(:)
  logical, intent(in), optional :: mask(:)

  real, pointer :: x(:),vec(:),tmp1(:),xt(:)
  integer :: v,i,j,j0,j1
  integer :: istat
  character(len=30) :: infoformat = '(A50,2E14.5)'

  if (size(filenames).ne.ML%nvar) then
    write(stderr,*) 'ERROR: Vector to save should be composed by ', &
         ML%nvar, ' variables, but ',size(filenames), &
         ' filenames are found: '
    write(stderr,*) (trim(filenames(v))//' ',v=1,size(filenames))
    call flush(stderr,istat)
    stop
  end if

  ! add additional exclusion points if mask is present

  if (present(mask)) then
    allocate(vec(size(vector)))
    where (mask)
      vec = vector
    elsewhere
      vec = FillValue
    end where
  else
! This does not work in PGI, but in g95
! bug in PGI compiler?

!    vec => vector

    allocate(vec(size(vector)))
    vec = vector
  end if

#ifdef DEBUG
  write(stddebug,'(A50,A14,A14)') 'saved variable','min','max'
#endif

  if (ML%distributed) then
    allocate(xt(ML%effsize))
    call parallGather(vec,xt,startIndexZones,endIndexZones)
 
    if (procnum.eq.1) then

      if (ML%removeLandPoints) then
        do v=1,size(filenames)

          allocate(tmp1(ML%varsize(v)) )

          j = 0

          do i=1,ML%varsize(v)
            if (ML%Mask(ML%startIndex(v) + i-1).eq.0) then
              tmp1(i) = FillValue
            else
              tmp1(i) = xt(invZoneIndex(ML%startIndexSea(v) + j))
              j = j+1
            end if
          end do


#ifdef DEBUG
          write(stddebug,infoformat) trim(path)//trim(filenames(v)),  &
               minval(tmp1,tmp1 /= FillValue),                     &
               maxval(tmp1,tmp1 /= FillValue)
#endif

          call usave(trim(path)//filenames(v), &
               tmp1,FillValue,4,ML%varshape(1:4,v),.false.)

          deallocate(tmp1)
        end do
      end if
    end if

    deallocate(xt)

  else

    if (procnum.eq.1) then
    if (ML%permute) call ipermute(zoneIndex,vec,vec)

    if (ML%removeLandPoints) then
      do v=1,size(filenames)
        j0 = ML%StartIndexSea(v)
        j1 = ML%EndIndexSea(v)


#ifdef DEBUG
        write(stddebug,infoformat) trim(path)//trim(filenames(v)),        &
             minval(vec(j0:j1)),  &
             maxval(vec(j0:j1))
#endif

        call usave(trim(path)//filenames(v), &
             reshape( &
             unpack(vec(j0:j1), &
             ML%Mask(ML%StartIndex(v):ML%EndIndex(v)).eq.1,FillValue),  &
                                !                     ML%varshape(1:ML%ndim(v),v)), &
             !(/ ML%varshape(1,v),ML%varshape(2,v),ML%varshape(3,v) /)       &
             ! ML%varshape(1:ML%ndim(v),v) &
             ML%varshape(1:4,v) &
             ),  &
             FillValue);


      end do
    else
      allocate(x(ML%totsize))
      where (ML%mask.eq.1)
        x = vec
      elsewhere
        x = FillValue
      end where

      do v=1,size(filenames)
        j0 = ML%StartIndex(v)
        j1 = ML%EndIndex(v)

#ifdef DEBUG
        write(stddebug,infoformat) trim(path)//trim(filenames(v)),  &
               minval(x(j0:j1),x(j0:j1) /= FillValue),              &
               maxval(x(j0:j1),x(j0:j1) /= FillValue)
#endif

        call usave(trim(path)//filenames(v), &
             reshape(x(j0:j1), &
                                !                     ML%varshape(1:ML%ndim(v),v)), &
             !(/ ML%varshape(1,v),ML%varshape(2,v),ML%varshape(3,v) /)), &
             ML%varshape(1:4,v)), &
             FillValue);



      end do

      deallocate(x)
      end if
    end if
  end if

#ifdef DEBUG
  write(stddebug,*) 
#endif

! Due to bug in PGI
!  if (present(mask)) deallocate(vec)
  deallocate(vec)
 end subroutine saveVector_byfilenames

 !_______________________________________________________
 !

 subroutine saveVector_byinitval(str,ML,vector,mask)
  use initfile
  use ufileformat
!  use parall
  character(len=*), intent(in) :: str
  type(MemLayout),  intent(in) :: ML
  real, intent(in) :: vector(:)
  logical, intent(in), optional :: mask(:)

  character(len=maxLen), pointer :: filenames(:)
  character(len=maxLen) :: path

!  write(stdout,*) 'minval(vector) ',__LINE__,minval(vector),size(vector),procnum

# ifdef DEBUG
  write(stddebug,'("== save Vector (",A,") ==")') str
# endif

  call getInitValue(initfname,str(1:index(str,'.',.true.))//'path',path,default='')
  call getInitValue(initfname,str,filenames)
  call saveVector_byfilenames(path,filenames,ML,vector,mask)
  deallocate(filenames)
 end subroutine saveVector_byinitval

 !_______________________________________________________
 !

 !_______________________________________________________
 !

 subroutine loadVectorSpace(str,ML,S,mean)
  use initfile
  use ufileformat
  use matoper
  use parall

  implicit none
  character(len=*), intent(in) :: str
  type(MemLayout),  intent(in) :: ML
  real, intent(out)            :: S(:,:)
  real, intent(out), optional  :: mean(:)

  integer                      :: v,k,dim,enstype
  character(len=maxLen), pointer   :: filenames(:), formats(:)
  character(len=maxLen)            :: prefix,path
  real                         :: scale
  real, pointer                :: spaceScale(:)
  logical                      :: doSpaceScaling = .false.
  real, allocatable            :: ensembleMean(:)
# ifdef PROFILE
  real(8) :: cputime(2)
# endif

# ifdef DEBUG
  write(stddebug,'("== load Vector Space (",A,") ==")') str
# endif

# ifdef PROFILE
  call cpu_time(cputime(1))
# endif    

  prefix = str(1:index(str,'.',.true.))
  call getInitValue(initfname,str,formats)
  call getInitValue(initfname,trim(prefix)//'path',path,default='')
  call getInitValue(initfname,trim(prefix)//'dimension',dim)
  call getInitValue(initfname,trim(prefix)//'scale',scale,default=1.)

! load space dependent scale is present 

  if (presentInitValue(initfname,trim(prefix)//'spaceScale')) then
    doSpaceScaling = .true.
    allocate(spaceScale(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel))
    call loadVector(trim(prefix)//'spaceScale',ML,spaceScale)
  end if

!!$omp parallel private(v,filenames)
  allocate(filenames(size(formats)))
!!$omp do
  do k=1,dim
    do v=1,size(filenames)
      write(filenames(v),formats(v)) k
    end do

    call loadVector(path,filenames,ML,S(:,k))
  end do
!!$omp end do
  deallocate(filenames)
!!$omp end parallel

! type = 1 vectors are anomalies
! type = 2 vectors are ensemble members (with mean)
  
  call getInitValue(initfname,trim(prefix)//'type',enstype,default=1)
  if (enstype.eq.2) then
! compute the ensemble mean
    allocate(ensembleMean(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel))

    ensembleMean = sum(S,2)/dim

    if (present(mean)) mean = ensembleMean;

#   ifdef DEBUG   
#    if ASSIM_SCALING == 0
       write(stddebug,*) 'remove ensemble mean and scale each member by 1/sqrt(dim)'
#    else
       write(stddebug,'("remove ensemble mean and scale each member by 1/sqrt(dim-",I2,")")') ASSIM_SCALING
#    endif
#   endif
  end if

! simple post processing of the ensemble
  do k=1,dim
    if (enstype.eq.2)  S(:,k) = (S(:,k)-ensembleMean)/sqrt(1.*dim - ASSIM_SCALING)
    if (doSpaceScaling) S(:,k) = spaceScale * S(:,k)
    if (scale.ne.1) S(:,k) = scale * S(:,k)
  end do

# ifdef PROFILE
  call cpu_time(cputime(2))

  write(stddebug,*) 'Profiling: loadVectorSpace',procnum
  write(stddebug,*) 'load data  ',cputime(2)-cputime(1)
  call flush(stddebug,istat)
# endif

  deallocate(formats)
  if (doSpaceScaling) deallocate(spaceScale)  
  if (enstype.eq.2)  deallocate(ensembleMean)
 end subroutine 

 !_______________________________________________________
 !

 subroutine loadEnsemble(str,ML,S)
  use initfile
  use ufileformat
  use matoper
  use parall
  implicit none
  character(len=*), intent(in) :: str
  type(MemLayout),  intent(in) :: ML
  real, intent(out)            :: S(:,:)

  integer                      :: v,k,dim
  character(len=maxLen), pointer   :: filenames(:), formats(:)
  character(len=maxLen)            :: prefix,path
  real                             :: scale
  real, allocatable                :: ensembleMean(:),ensembleAnom(:)
  real, pointer                    :: spaceScale(:)
  logical                          :: doSpaceScaling = .false.

# ifdef PROFILE
  real(8) :: cputime(2)
# endif

# ifdef DEBUG
  write(stddebug,'("== load ensemble (",A,") ==")') str
# endif

# ifdef PROFILE
  call cpu_time(cputime(1))
# endif    

  prefix = str(1:index(str,'.',.true.))
  call getInitValue(initfname,str,formats)
  call getInitValue(initfname,trim(prefix)//'path',path,default='')
  call getInitValue(initfname,trim(prefix)//'dimension',dim)


!!$omp parallel private(v,filenames)
  allocate(filenames(size(formats)))
!!$omp do
  do k=1,dim
    do v=1,size(filenames)
      write(filenames(v),formats(v)) k
    end do

    call loadVector(path,filenames,ML,S(:,k))
  end do
!!$omp end do
  deallocate(filenames)
!!$omp end parallel

# ifdef PROFILE
  call cpu_time(cputime(2))

  write(stddebug,*) 'Profiling: loadVectorSpace',procnum
  write(stddebug,*) 'load data  ',cputime(2)-cputime(1)
  call flush(stddebug,istat)
# endif

  ! scaling of ensemble
  call getInitValue(initfname,trim(prefix)//'scale',scale,default=1.)

! load space dependent scale is present 
  if (presentInitValue(initfname,trim(prefix)//'spaceScale')) then
    doSpaceScaling = .true.
    allocate(spaceScale(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel))
    call loadVector(trim(prefix)//'spaceScale',ML,spaceScale)
  end if

  if (scale /= 1 .or. doSpaceScaling) then
    write(stddebug,*) 'Apply scaling to ensemble ',scale

    allocate( &
      ensembleMean(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel), &
      ensembleAnom(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel) &     
      )

    ensembleMean = sum(S,2)/dim

    ! simple post processing of the ensemble
    do k=1,dim
      ensembleAnom = S(:,k)-ensembleMean

      if (scale /= 1) ensembleAnom = scale * ensembleAnom
      if (doSpaceScaling) ensembleAnom = spaceScale * ensembleAnom

      S(:,k) = ensembleMean + ensembleAnom
    end do

    deallocate(ensembleMean,ensembleAnom)
  end if

 end subroutine loadEnsemble

 !_______________________________________________________
 !


 !_______________________________________________________
 !

 subroutine saveEnsemble(str,ML,S)
  use initfile
  use ufileformat
  character(len=*), intent(in) :: str
  type(MemLayout),  intent(in) :: ML
  real, intent(in)             :: S(:,:)

  integer                      :: v,k
  character(len=maxLen), pointer   :: filenames(:), formats(:)
  character(len=maxLen)            :: prefix,path

  prefix = str(1:index(str,'.',.true.))
  call getInitValue(initfname,str,formats)
  call getInitValue(initfname,trim(prefix)//'path',path,default='')

  allocate(filenames(size(formats)))

  do k=1,size(S,2)
    do v=1,size(filenames)
      write(filenames(v),formats(v)) k
    end do

    call saveVector(path,filenames,ML,S(:,k))
  end do

  deallocate(formats,filenames)
 end subroutine saveEnsemble


 !_______________________________________________________
 !
 !_______________________________________________________
 !

 subroutine loadErrorSpace(str,S,mean)
  use initfile
  use ufileformat
  implicit none
  character(len=*), intent(in) :: str
  real, intent(out)            :: S(:,:)
  real, intent(out), optional  :: mean(:)

  call loadVectorSpace(str,ModMLParallel,S,mean)
 end subroutine loadErrorSpace

 !_______________________________________________________
 !

 subroutine saveErrorSpace(str,S)
  use initfile
  use ufileformat
  character(len=*), intent(in) :: str
  real, intent(in)             :: S(:,:)

  call saveEnsemble(str,ModMLParallel,S)
 end subroutine saveErrorSpace


 !_______________________________________________________
 !
 !_______________________________________________________
 !

 subroutine loadSparseMatrix_byfilenames(path,filename,ML1,ML2,H,valid1,valid2)
  use ufileformat
  use matoper
  implicit none

  character(*), intent(in) :: filename,path
  type(MemLayout), intent(in) :: ML1,ML2
  type(SparseMatrix), intent(out)  :: H
  logical, optional, intent(out), dimension(:) :: valid1, valid2

  integer :: istat
  integer :: i,j
  real :: valex
  real, pointer :: Hop(:,:)
  integer, allocatable :: Hindex(:,:)
  real, allocatable :: Hcoeff(:)

  call uload(trim(path)//filename,Hop,valex,check_vshape=[9,-1])

  allocate(Hcoeff(size(Hop,2)),Hindex(10,size(Hop,2)))

! bug in compiler ?
  do j=1,size(Hop,2)
    do i=1,10
      Hindex(i,j) = Hop(i,j)
    end do
    Hcoeff(j) = Hop(11,j)
  end do
  !Hindex = Hop(:10,:)
  !Hcoeff = Hop(11,:)

  deallocate(Hop)

  call packSparseMatrix(Hindex,Hcoeff,ML1,ML2,H,valid1,valid2)
  deallocate(Hindex,Hcoeff)

# ifdef DEBUG
  write(stddebug,*) 
  write(stddebug,'(A)') '== Sparse operator =='
  call writeInfoSparseMatrix(stddebug,H)
  call flush(stddebug,istat)
# endif

 end subroutine loadSparseMatrix_byfilenames

 !_______________________________________________________
 !

 subroutine loadSparseMatrix_byinitval(str,ML1,ML2,H,valid1,valid2)
  use initfile
  use matoper
  implicit none
  character(len=*), intent(in) :: str
  type(MemLayout), intent(in) :: ML1,ML2
  type(SparseMatrix), intent(out)  :: H
  logical, optional, intent(out), dimension(:) :: valid1, valid2

  character(256) :: filename,path
  
  call getInitValue(initfname,str(1:index(str,'.',.true.))//'path',path,default='')
  call getInitValue(initfname,str,filename)

  call loadSparseMatrix_byfilenames(path,filename,ML1,ML2,H,valid1,valid2)
 end subroutine loadSparseMatrix_byinitval

 !_______________________________________________________
 !

 subroutine saveSparseMatrix_byfilenames(path,filename,ML1,ML2,H,valid1,valid2)
  use ufileformat
  use matoper
  implicit none

  character(*), intent(in) :: filename,path
  type(MemLayout), intent(in) :: ML1,ML2
  type(SparseMatrix), intent(in)  :: H
  logical, optional, intent(in) :: valid1(:),valid2(:)

  integer :: i,j,nz,k
  integer :: istat
  real :: valex=9999.
  real, allocatable    :: Hop(:,:),Hcoeff(:)
  integer, allocatable :: Hindex(:,:)


  nz = H%nz
  if (present(valid1)) nz = nz + count(.not.valid1)
  if (present(valid2)) nz = nz + count(.not.valid2)

  allocate(Hindex(10,nz),Hcoeff(nz),Hop(11,nz))


  call unpackSparseMatrix(Hindex(:,:H%nz),Hcoeff(:H%nz),ML1,ML2,H)

  ! out of grid values

  k = H%nz

  if (present(valid1)) then
    do i=1,size(valid1)
      if (.not.valid1(i)) then
        k=k+1
        call ind2sub(ML1,i,Hindex(1,k),Hindex(2,k),Hindex(3,k),Hindex(4,k),Hindex(5,k))
        Hindex(6:10,k) = -1
        Hcoeff(k) = 0.
      end if
    end do
  end if

  if (present(valid2)) then
    do i=1,size(valid2)
      if (.not.valid2(i)) then
        k=k+1
        Hindex(1:5,k) = -1
        call ind2sub(ML2,i,Hindex(6,k),Hindex(7,k),Hindex(8,k),Hindex(9,k),Hindex(10,k))
        Hcoeff(k) = 0.
      end if
    end do
  end if


! bug in compiler ?
  do j=1,size(Hcoeff)
    do i=1,10
      Hop(i,j) = Hindex(i,j)
    end do
    Hop(11,j) = Hcoeff(j) 
  end do
  !Hop(1:10,:) = Hindex
  !Hop(11,:)   = Hcoeff

  call usave(trim(path)//filename,Hop,valex)

  deallocate(Hindex,Hcoeff,Hop)
 end subroutine saveSparseMatrix_byfilenames

 !_______________________________________________________
 !

 subroutine saveSparseMatrix_byinitval(str,ML1,ML2,H,valid1,valid2)
  use initfile
  use matoper
  implicit none
  character(len=*), intent(in) :: str
  type(MemLayout), intent(in) :: ML1,ML2
  type(SparseMatrix), intent(in)  :: H
  logical, optional, intent(in) :: valid1(:),valid2(:)

  character(256) :: filename,path=''

  !    call getInitValue(initfname,str(1:index(str,'.',.true.))//'path',path,default='')
  call getInitValue(initfname,str,filename)
  call saveSparseMatrix_byfilenames(path,filename,ML1,ML2,H,valid1,valid2)

 end subroutine saveSparseMatrix_byinitval
 !_______________________________________________________
 !

 ! joins a list of integer where all elements are separated by sep 
 ! (or 'x' is the default value)

 function sizeformat(sz,sep) result(s)
  implicit none
  integer, intent(in) :: sz(:)
  character, optional :: sep
  character(len=50) :: s,str
  character :: sepa
  integer :: i

  if (present(sep)) then
    sepa = sep
  else
    sepa = 'x'
  end if

  s = ''
  do i=1,size(sz)
    write (str,*) sz(i)

    s = trim(s) // trim(ADJUSTL(str)) // sepa
  end do

  s = s(1:len_trim(s)-1)

 end function sizeformat

 !_______________________________________________________
 !

 subroutine MemoryLayout_byfilenames(path,filenames,la,packed)
  use ufileformat
  implicit none
  character(len=*), intent(in) :: path
  character(len=*), intent(in) :: filenames(:)

  type(MemLayout), intent(out) :: la
  logical, intent(in), optional :: packed  
  real, pointer        :: tmp(:)

  integer :: m,mmax,prec,i,j
  integer :: istat
  real :: valex
  logical :: isdegen
  character(len=MaxLen) :: ssize

  if (present(packed)) then
    la%removeLandPoints = packed
  else
    la%removeLandPoints = removeLandPoints
  end if

  mmax = size(filenames)
  la%nvar = mmax
  
  allocate( &
       la%StartIndex(mmax),la%EndIndex(mmax),la%VarSize(mmax),&
       la%StartIndexSea(mmax),la%EndIndexSea(mmax),la%VarSizeSea(mmax),&
       la%ndim(mmax), la%varshape(5,mmax))

  la%varshape = 1
  la%StartIndex(1) = 1

  do m=1,mmax
    call uinquire(trim(path)//filenames(m),valex,prec,la%ndim(m), &
         la%varshape(:,m), isdegen)


    ! make code more general

    if (la%ndim(m) == 4) then
      if (la%varshape(4,m).eq.1) la%ndim(m) = 3

    elseif (la%ndim(m) == 3) then

      if (la%varshape(3,m).eq.1) la%ndim(m) = 2
    end if
    !write(6,*) 'la%varshape(:,m), la%ndim(m)',la%varshape(:,m), la%ndim(m)

    la%VarSize(m) = product(la%varshape(1:la%ndim(m),m))
    la%EndIndex(m) = la%StartIndex(m) + la%VarSize(m)-1 
    if (m.ne.mmax) la%StartIndex(m+1) = la%EndIndex(m)+1
  end do

  la%totsize = la%EndIndex(mmax)

  allocate(la%Mask(la%totsize),la%SeaIndex(la%totsize))

  la%StartIndexSea(1) = 1

  do m=1, mmax
    call uload(trim(path)//filenames(m),tmp,valex)

    where (tmp.eq.valex)
      tmp = 0
    elsewhere
      tmp = 1
    end where

    la%Mask(la%StartIndex(m):la%EndIndex(m)) = tmp

    la%VarSizeSea(m) = count(tmp.eq.1)
    la%EndIndexSea(m) = la%StartIndexSea(m) + la%VarSizeSea(m)-1
    if (m.ne.mmax) la%StartIndexSea(m+1) = la%EndIndexSea(m)+1

    deallocate(tmp)
  end do

  la%totsizeSea = count(la%Mask.eq.1)
  if (la%removeLandPoints) then
    la%effsize = la%totsizeSea
  else
    la%effsize = la%totsize
  end if

  allocate(la%invindex(la%effsize))

  la%SeaIndex = -1
  j=1
  do i=1,la%totsize
    if (.not.la%removeLandPoints.or.la%Mask(i).eq.1) then
      la%SeaIndex(i) = j
      la%invindex(j) = i
      j=j+1
    end if
  end do

  la%distributed = .false.
  la%permute = .false.
  la%startIndexParallel = 1
  la%endIndexParallel = la%effsize

# ifdef DEBUG
  write(stddebug,'(A20,A15)') 'path:                ',trim(path)
  write(stddebug,'(A20,L15)') 'removeLandPoints:    ',la%removeLandPoints
  write(stddebug,'(A20,I15)') 'totsize:             ',la%totsize
  write(stddebug,'(A20,I15)') 'totsizeSea:          ',la%totsizeSea
  write(stddebug,'(A20,I15)') 'effsize:             ',la%effsize


  write(stddebug,'(A)') 
  write(stddebug,'(A2," ",A20,A20,2A30)') '#','Mask','Shape','....... all elements .......','...... only non-masked ......'
  write(stddebug,'(43X,6A10)') 'Size','Start','End','Size','Start','End'

  do m=1,mmax
    ssize = sizeformat(la%varshape(1:la%ndim(m),m))
    write(stddebug,'(I2," ",A20,A20,6I10)') m,trim(filenames(m)),trim(ssize), &
         la%varsize(m), &
         la%StartIndex(m), &
         la%EndIndex(m), &
         la%varsizesea(m), &
         la%StartIndexSea(m), &
         la%EndIndexSea(m)
  end do
  write(stddebug,*)

  call flush(stddebug,istat)
# endif
 end subroutine MemoryLayout_byfilenames

 !_______________________________________________________
 !

 subroutine MemoryLayout_byinitval(prefix,la,packed)
  use initfile
  use ufileformat
  implicit none
  character(len=*) :: prefix
  type(MemLayout), intent(out) :: la
  logical, intent(in), optional :: packed  

  character(len=maxLen), pointer :: filenames(:)
  character(len=maxLen) :: path

# ifdef DEBUG
  write(stddebug,'("== Memory Layout (",A,") ==")') prefix
# endif

  call getInitValue(initfname,trim(prefix)//'path',path,default='')
  call getInitValue(initfname,trim(prefix)//'mask',filenames)
  call MemoryLayout_byfilenames(path,filenames,la,packed)

  call getInitValue(initfname,trim(prefix)//'variables',la%varnames)

  deallocate(filenames)
 end subroutine MemoryLayout_byinitval


 !_______________________________________________________
 !

 function sub2ind(ML,v,i,j,k,n,valid) result(index)
  implicit none
  type(MemLayout) :: ML
  integer :: index
  integer, intent(in) :: v,i,j,k,n
  logical, optional,intent(out) :: valid

  logical :: val


  index = -1
  val = 1 <= v.and.v <= ML%nvar

  if (val) then
    val = &
         1.le.i.and.i.le.ML%varshape(1,v).and. &
         1.le.j.and.j.le.ML%varshape(2,v).and. &
         1.le.k.and.k.le.ML%varshape(3,v).and. &
         1.le.n.and.n.le.ML%varshape(4,v)

    if (val) then
      index = ML%StartIndex(v) + i-1 + & 
           ML%varshape(1,v) * (j-1 + ML%varshape(2,v) * (k-1 + ML%varshape(3,v) * (n-1)))

      index = ML%SeaIndex(index)
      val =  index.ne.-1
    end if


  end if

  if (present(valid)) valid=val

 end function sub2ind

 !_______________________________________________________
 !
 ! index includes only sea points
 ! Fix me:
 ! make it work for distributed and permuted vectors

 subroutine ind2sub(ML,index,v,i,j,k,n)
  implicit none
  type(MemLayout) :: ML
  integer, intent(in) :: index
  integer, intent(out) :: v,i,j,k,n
  integer :: linindex


!  if (ML%permute) then
!    write(stderr,*) 'ind2sub: not implemented '
!  end if
! currently: should be "unpermuted" at calling level

  linindex = ML%invindex(index)

  do v=1,ML%nvar
    if (ML%startIndex(v).le.linindex.and.linindex.le.ML%endIndex(v)) then
      linindex = linindex - ML%startIndex(v) 
      ! i,j,k,n zero-based
      n = linindex/(ML%varshape(1,v)*ML%varshape(2,v)*ML%varshape(3,v))
      linindex = linindex - n*(ML%varshape(1,v)*ML%varshape(2,v)*ML%varshape(3,v))

      k = linindex/(ML%varshape(1,v)*ML%varshape(2,v))
      linindex = linindex - k*(ML%varshape(1,v)*ML%varshape(2,v))

      j = linindex/ML%varshape(1,v)
      linindex = linindex - j*(ML%varshape(1,v))

      i = linindex

      ! i,j,k one-based
      i=i+1; j=j+1; k=k+1; n=n+1;
      return
    end if
  end do
  v = -1
 end subroutine ind2sub


 !_______________________________________________________
 !
 ! deallocate pointers in type MemLayout
 !

 subroutine MemoryLayoutDone(la)
  use initfile
  use ufileformat
  implicit none
  type(MemLayout), intent(out) :: la

  deallocate( &
       la%StartIndex,la%EndIndex,la%VarSize,&
       la%StartIndexSea,la%EndIndexSea,la%VarSizeSea,&
       la%varnames, &
       la%ndim, la%varshape, &
       la%Mask, la%SeaIndex, la%invindex)

# ifdef DEBUG
  write(stddebug,*) 'Memory Layout: deallocate '
# endif
 end subroutine MemoryLayoutDone

 !_______________________________________________________
 !
 !
 ! specialised input and output routines
 ! depreciated

 subroutine loadStateVector_byfilenames(path,filenames,StateVector)
  use initfile
  use ufileformat
  character(len=*), intent(in) :: path
  character(len=*), intent(in) :: filenames(:)
  real, intent(out) :: StateVector(:)

  call loadVector_byfilenames(path,filenames,ModML,StateVector)
 end subroutine loadStateVector_byfilenames

 !_______________________________________________________
 !
 ! depreciated

 subroutine loadStateVector_byinitval(str,StateVector)
  use initfile
  use ufileformat
  use matoper
  character(len=*), intent(in) :: str
  real, intent(out) :: StateVector(:)

  call loadVector_byinitval(str,ModMLParallel,StateVector)
 end subroutine loadStateVector_byinitval

 !_______________________________________________________
 !
 ! depreciated

 subroutine saveStateVector_byfilenames(path,filenames,StateVector)
  use initfile
  use ufileformat
  character(len=*), intent(in) :: path
  character(len=*), intent(in) :: filenames(:)
  real, intent(in) :: StateVector(:)

  call saveVector_byfilenames(path,filenames,ModMLParallel,StateVector)
 end subroutine saveStateVector_byfilenames

 !_______________________________________________________
 !
 ! depreciated

 subroutine saveStateVector_byinitval(str,StateVector)
  use initfile
  use ufileformat
  character(len=*), intent(in) :: str
  real, intent(in) :: StateVector(:)

  call saveVector_byinitval(str,ModMLParallel,StateVector)
 end subroutine saveStateVector_byinitval

 !_______________________________________________________
 !
 ! load the modified julian day of 'ntime'th observation
 ! error is different from zero if this observation does not exist 
 !_______________________________________________________
 !

 subroutine loadObsTime(ntime,mjd0,error)
  use initfile
  use date
  use ufileformat
  implicit none
  integer, intent(in)  :: ntime
  real(8), intent(out) :: mjd0
  integer, intent(out), optional :: error

  character(len=maxLen) :: prefix,str
  integer :: day,month,year,h,min
  integer :: istat
  real :: s,seconds

  !write(prefix,'(A,I3.3,A)') 'Obs',ntime,'.'
  call fmtIndex('Obs',ntime,'.',prefix)

  if (presentInitValue(initfname,trim(prefix)//'time')) then

    !
    !  ObsXXX.time
    !

    call getInitValue(initfname,trim(prefix)//'time',str)


    
    if (index(str,'T') /= 0) then
      ! ISO 8601 time format: YYYY-MM-DDThh:mm:ss

      read(str(1:4),*) year
      read(str(6:7),*) month
      read(str(9:10),*) day

      read(str(12:13),*) h
      read(str(15:16),*) min
      read(str(18:),*) s

      seconds = 60 * (60*h + min) + s
      mjd0 = mjd(year,month,day,seconds)

    elseif (index(str,':') /= 0) then      
      ! time: hh:mm:ss

      read(str(1:2),*) h
      read(str(4:5),*) min
      read(str(7:),*) s

      !
      !  ObsXXX.date
      !

      call getInitValue(initfname,trim(prefix)//'date',str)

      ! date: dd/mm/yyyy 
      read(str(1:2),*) day
      read(str(4:5),*) month
      read(str(7:),*) year

      seconds = 60 * (60*h + min) + s
      mjd0 = mjd(year,month,day,seconds)

    else
      ! time is directly in mjd
      read(str,*) mjd0

    end if


    !  write(stdout,*) str, seconds

    
    if (present(error)) error = 0

  else
    ! no observation with number ntime
    if (present(error)) then
      error = -1
      return
    else
      write(stderr,*) 'No observation number ',ntime
      call flush(stderr,istat)
      stop 'NO OBS'
    end if
  end if

 end subroutine loadObsTime


 !_______________________________________________________
 !
 ! load the measurments and theirs associated error
 ! 
 !_______________________________________________________
 !


 subroutine loadObs(ntime,ObsML,observation,invsqrtR)
  use initfile
  use date
  use ufileformat
  implicit none
  integer, intent(in)  :: ntime
  type(MemLayout), intent(in) :: ObsML
  real,    intent(out) :: observation(:), invsqrtR(:)

  character(len=maxLen)          :: path        
  character(len=maxLen) :: prefix
  integer :: omax
  integer :: istat
  real, parameter :: min_rmse = 0.01

  !write(prefix,'(A,I3.3,A)') 'Obs',ntime,'.'
  call fmtIndex('Obs',ntime,'.',prefix)

  call getInitValue(initfname,trim(prefix)//'path',path,default='')

  omax = size(ObsML%Mask)

  !    allocate(observation(ObsML%effsize),invsqrtR(ObsML%effsize))

  !
  !  ObsXXX.value
  !

  call loadVector(trim(prefix)//'value',ObsML,observation)

#   ifdef DEBUG
  write(stddebug,*) 'sum(obs) ',sum(observation), shape(observation)
  call flush(stddebug,istat)
#   endif
  !
  !  ObsXXX.rmse
  !

  call loadVector(trim(prefix)//'rmse',ObsML,invsqrtR)

  if (any(invsqrtR <= 0)) then
    write (0,*) 'error in ',__FILE__,__LINE__
    write (0,*) 'zero value found in ',trim(prefix)//'rmse'
    write (0,*) 'are your observations so good?'
    ERROR_STOP
  end if

# ifdef LIMIT_iR  
  if (any(invsqrtR < min_rmse .and. ObsML%mask == 1)) then
    write(stdlog,*) 'loadObs: warning observations with too small errors'
  end if

  where (invsqrtR < min_rmse .and. ObsML%mask == 1)
     invsqrtR = min_rmse
  end where
# endif

# ifdef DEBUG
  write(stddebug,*) 'assim invsqrtR',sum(invsqrtR),count(invsqrtR.eq.0) 
# endif

  if (ObsML%removeLandPoints) then
    invsqrtR = 1./invsqrtR
  else
    where (ObsML%mask.eq.1) 
      invsqrtR = 1./invsqrtR
    elsewhere
      invsqrtR = 0
    end where
  end if

 end subroutine loadObs


 !_______________________________________________________
 !
 ! load the measurments correlation
 ! 
 !_______________________________________________________
 !

 subroutine loadObservationCorr(ntime,ObsML,C)
  use initfile
!  use grids
  use ufileformat
  use matoper
  implicit none
  integer, intent(in)  :: ntime
  type(MemLayout), intent(in) :: ObsML
  type(SparseMatrix), intent(out)  :: C

  character(len=maxLen)          :: path

  character(len=maxLen) :: prefix,str,Dprefix
  integer :: mmax
  integer :: istat

  integer :: error
  real :: valex

  real, pointer :: Cop(:,:)
  integer, pointer :: Cindex(:,:)
  real, pointer :: Ccoeff(:)
  logical, pointer :: validobs(:)


  !write(prefix,'(A,I3.3,A)') 'Obs',ntime,'.'
  call fmtIndex('Obs',ntime,'.',prefix)

  call getInitValue(initfname,trim(prefix)//'path',path,default='')

  mmax = size(ObsML%VarSize)

  call getInitValue(initfname,trim(prefix)//'Correlation',str,error)

  if (error.eq.0) then
#ifdef DEBUG
    write(stddebug,*) 'Load observation correlation: ',trim(str)
#endif
    call uload(trim(path)//str,Cop,valex,check_vshape=[9,-1])
    allocate(Ccoeff(size(Cop,2)),Cindex(10,size(Cop,2)))
    Cindex = Cop(1:10,:)
    Ccoeff = Cop(11,:)
    deallocate(Cop)
  else
#ifdef DEBUG
    write(stddebug,*) 'Generate observation correlation'
#endif
    call genObservationCorr(ntime,ObsML,Cindex,Ccoeff)
  end if

  ! save the obervation operator if desiered

  !write(Dprefix,'(A,I3.3,A)') 'Diag',ntime,'.'
  call fmtIndex('Diag',ntime,'.',Dprefix)


  if (presentInitValue(initfname,trim(Dprefix)//'Correlation')) then 
    call getInitValue(initfname,trim(Dprefix)//'Correlation',str)
    call getInitValue(initfname,trim(Dprefix)//'path',path,default='')
    allocate(Cop(11,size(Ccoeff)))
    Cop(1:10,:) = Cindex
    Cop(11,:) = Ccoeff
    call usave(trim(path)//trim(str),Cop,9999.)      
    deallocate(Cop)
  end if

  allocate(validobs(size(Ccoeff)))
  call packSparseMatrix(Cindex,Ccoeff,ObsML,ObsML,C)

  !    where (.not.validobs) invsqrtR = 0
  deallocate(validobs)


# ifdef DEBUG
  write(stddebug,*) 
  write(stddebug,'(A)') '== Observation Correlation Matrix C =='
  call writeInfoSparseMatrix(stddebug,C)
  call flush(stddebug,istat)
# endif

 end subroutine loadObservationCorr


 !_______________________________________________________
 !
 ! load the observation operator
 ! 
 ! order of preceding of groups of keys in initfile
 !
 ! 1. ObsXXX.operator
 !
 ! 3. ObsXXX.HperObs
 !
 ! 2. ObsXXX.variables
 !    ObsXXX.gridX
 !    ObsXXX.gridY
 !    ObsXXX.gridZ
 !    ObsXXX.gridT (optional)
 !
 !_______________________________________________________
 !

 subroutine loadObservationOper(ntime,ObsML,H,Hshift,invsqrtR)
  use initfile
!  use grids
  use ufileformat
  use matoper
  implicit none
  integer, intent(in)  :: ntime
  type(MemLayout), intent(in) :: ObsML
  type(SparseMatrix), intent(out)  :: H
  real, intent(out) :: Hshift(:)
  real, intent(inout) :: invsqrtR(:)

  type(SparseMatrix)  :: filterMod, tmpH
  character(len=maxLen), pointer :: filenames(:)
  character(len=maxLen)          :: path

  character(len=maxLen) :: prefix,str,Dprefix
  integer :: m,mmax,nz,idummy,nbentries

  integer :: istat
  integer :: i,j
  real :: valex
  logical :: isdegen
  real, pointer :: Hop(:,:)
  integer, pointer :: Hindex(:,:)
  real, pointer :: Hcoeff(:)
  logical, allocatable :: validobs(:)
  real, allocatable :: shiftMod(:)

  !write(prefix,'(A,I3.3,A)') 'Obs',ntime,'.'
  call fmtIndex('Obs',ntime,'.',prefix)

  call getInitValue(initfname,trim(prefix)//'path',path,default='')

  mmax = size(ObsML%VarSize)

  if (presentInitValue(initfname,trim(prefix)//'operator')) then

    call getInitValue(initfname,trim(prefix)//'operator',str)
#ifdef DEBUG
    write(stddebug,*) 'Load observation operator: ',trim(str)
#endif
    call uload(trim(path)//str,Hop,valex,check_vshape=[11,-1])
    
    allocate(Hcoeff(size(Hop,2)),Hindex(10,size(Hop,2)))
    do j=1,size(Hop,2)
      do i=1,10
        Hindex(i,j) = Hop(i,j)
      end do
      Hcoeff(j) = Hop(11,j)
    end do
    !Hindex = Hop(:10,:)
    !Hcoeff = Hop(11,:)
    deallocate(Hop)

  elseif (presentInitValue(initfname,trim(prefix)//'HperObs')) then
#ifdef DEBUG
    write(stddebug,*) 'assemble observation operator'
#endif
    call getInitValue(initfname,trim(prefix)//'HperObs',filenames)

    nz = 0
    do m=1,ObsML%nvar
      call uinquire(trim(path)//filenames(m),valex,idummy,idummy,nbentries,idummy,idummy,isdegen)
      nz = nz + nbentries
    end do

    allocate(Hcoeff(nz),Hindex(10,nz))

    nz=1

    do m=1,ObsML%nvar
      call uload(trim(path)//filenames(m),Hop,valex,check_vshape=[11,-1])

      do i=1,size(Hop,2)
        Hindex(1,nz) = m
        Hindex(2:10,nz) = floor(Hop(2:10,i)+.5)
        Hcoeff(nz) = Hop(11,i)
        nz = nz + 1         
      end do
      deallocate(Hop)
    end do

    deallocate(filenames)
  else
#ifdef DEBUG
    write(stddebug,*) 'Generate observation operator'
#endif
    call genObservationOper(ntime,ObsML,Hindex,Hcoeff)

  end if


  ! save the obervation operator if desiered

  !write(Dprefix,'(A,I3.3,A)') 'Diag',ntime,'.'
  call fmtIndex('Diag',ntime,'.',Dprefix)


  if (presentInitValue(initfname,trim(Dprefix)//'H')) then 
    call getInitValue(initfname,trim(Dprefix)//'H',str)
    call getInitValue(initfname,trim(Dprefix)//'path',path,default='')
    allocate(Hop(11,size(Hcoeff)))
    Hop(1:10,:) = Hindex
    Hop(11,:) = Hcoeff
    call usave(trim(path)//trim(str),Hop,9999.)      
    deallocate(Hop)
  end if

  allocate(validobs(ObsML%effsize))
  call packSparseMatrix(Hindex,Hcoeff,ObsML,ModML,H,validobs)

  where (.not.validobs) invsqrtR = 0
  deallocate(validobs,Hindex,Hcoeff)

  if (presentInitValue(initfname,trim(prefix)//'shift')) then
    call loadVector(trim(prefix)//'shift',ObsML,Hshift)
  else
    Hshift = 0.
  end if

  if (presentInitValue(initfname,trim(prefix)//'shiftMod')) then
    allocate(shiftMod(ModML%effsize))
    call loadVector(trim(prefix)//'shiftMod',ModML,shiftMod)    
    Hshift = (H.x.shiftMod) + Hshift
    deallocate(shiftMod)
  end if

  if (presentInitValue(initfname,trim(prefix)//'filterMod')) then
    call loadSparseMatrix(trim(prefix)//'filterMod',ModML,ModML,filterMod)
    ! bug deallocation of H ?
    tmpH%i => H%i
    tmpH%j => H%j
    tmpH%s => H%s
    H = tmpH.x.filterMod
    deallocate(tmpH%i,tmpH%j,tmpH%s,filterMod%i,filterMod%j,filterMod%s)
  end if


  if (schemetype.eq.LocalScheme .and. allocated(invZoneIndex)) then
    ! permute state vector
    do j=1,H%nz
      H%j(j) = invZoneIndex(H%j(j))
    end do
  endif

#ifdef DEBUG
  write(stddebug,*) 
  write(stddebug,'(A)') '== Observation operator H =='
  call writeInfoSparseMatrix(stddebug,H)
  call flush(stddebug,istat)
#endif
 end subroutine loadObservationOper

 !_______________________________________________________
 !
 ! writes some information of sparse matrix H on unit

 subroutine writeInfoSparseMatrix(unit,H)
  use matoper
  implicit none
  integer, intent(in) :: unit
  type(SparseMatrix), intent(in)  :: H
  character(len=MaxLen) :: ssize

  integer :: maxi, maxj

  maxi = maxval(H%i(1:H%nz))
  maxj = maxval(H%j(1:H%nz))

  ssize = sizeformat((/H%m,H%n/))
  write(unit,'(A,A10)')   'Matrix shape:       ',trim(ssize)
  write(unit,'(A,I10)')   '# non-zero elements:',H%nz
  write(unit,'(A,I10)')   'max(H%i):           ',maxi
  write(unit,'(A,I10)')   'min(H%i):           ',minval(H%i(1:H%nz))
  write(unit,'(A,I10)')   'max(H%j):           ',maxj
  write(unit,'(A,I10)')   'min(H%j):           ',minval(H%j(1:H%nz))
  write(unit,'(A,E10.3)') 'max(H%s):           ',maxval(H%s(1:H%nz))
  write(unit,'(A,E10.3)') 'min(H%s):           ',minval(H%s(1:H%nz)) 
  write(unit,'(A,E10.3)') 'mean(H%s):          ',sum(H%s(1:H%nz))/H%nz

  if (maxi > H%m .or. maxj > H%n) then
    write(stderr,*) 'Index exceeds range '
    ERROR_STOP
  end if  
 end subroutine


 !_______________________________________________________
 !
 ! routines for constructing objects
 !
 !  

 !_______________________________________________________
 !
 ! create the observation operator
 ! 
 !_______________________________________________________
 !

 subroutine genObservationOper(ntime,ObsML,Hindex,Hcoeff)
  use initfile
!  use grids
  use ufileformat
  implicit none
  integer, intent(in)  :: ntime
  integer, pointer     :: Hindex(:,:)
  real, pointer        :: Hcoeff(:)

  integer, allocatable :: tmpHindex(:,:)
  real, allocatable    :: tmpHcoeff(:)
  character(len=maxLen), pointer :: varNames(:)
  character(len=maxLen)    :: path
  real, allocatable, dimension(:) :: obsX, obsY, obsZ, obsT

  character(len=maxLen)   :: prefix
  type(MemLayout), intent(in) :: ObsML

  integer              ::  &
       i,j,k,n, &
       v,tv,m,mmax,omaxSea,tn,nz,linindex, &
       tindexes(4,16), tmpm, tn_test
  integer :: istat
  real                 :: tc(16), minres
  logical              :: ingrid

  !write(prefix,'(A,I3.3,A)') 'Obs',ntime,'.'
  call fmtIndex('Obs',ntime,'.',prefix)

  call getInitValue(initfname,trim(prefix)//'path',path,default='')

  mmax = size(ObsML%VarSize)
  omaxSea = ObsML%effsize


  ! assume that land-point are removed 
  ! ObsML%removeLandPoints == .true.

  ! n highest dimensions (here n=4)
  ! 10 = 2 + 2 n
  ! 16 = 2^n
  allocate(tmpHindex(10,16*omaxSea),tmpHcoeff(16*omaxSea))
  nz = 0

  call getInitValue(initfname,trim(prefix)//'variables',varNames)

  ! load position of observations
  allocate(obsX(ObsML%effsize),obsY(ObsML%effsize),          &
       obsZ(ObsML%effsize),obsT(ObsML%effsize))

  call loadVector(trim(prefix)//'gridX',ObsML,obsX)
  call loadVector(trim(prefix)//'gridY',ObsML,obsY)
  call loadVector(trim(prefix)//'gridZ',ObsML,obsZ)

  if (presentInitValue(initfname,trim(prefix)//'gridT')) then
    call loadVector(trim(prefix)//'gridT',ObsML,obsT)
  end if


  do m=1, mmax
    ! preliminary check if variable is known
    v = -1 

    do tv=1,size(ModML%varnames)
      if (varNames(m).eq.ModML%varnames(tv)) then
        v = tv
        exit
      end if
    end do

    ! unknown variable
    if (v.eq.-1) then
      write(stderr,*) 'genObservationOper: WARNING ',trim(varNames(m)),' is unknown !'
      write(stderr,*) 'known variables are:', &
           (trim(ModML%varnames(tv))//' ',tv=1,size(ModML%varnames))
      write(stderr,*) 'unknown variables are treated as absent.'
      call flush(stderr,istat)
    end if

    ! loop over all observation (only non-masked) in the mth observation
    do linindex = ObsML%StartIndexSea(m),ObsML%EndIndexSea(m)
        call ind2sub(ObsML,linindex,tmpm,i,j,k,n)

        minres = huge(minres)
        v = -1
        tn = 0
        ingrid = .false.

        do tv=1,size(ModML%varnames)
          if (varNames(m).eq.ModML%varnames(tv).and.minres.ge.hres(tv)) then
            ! known variable
            v = tv

            ! compute interpolation coefficients
            
            tindexes = 1

            if (ModML%ndim(v).eq.1) then
              call cinterp(ModelGrid(v), (/ obsX(linindex) /), &
                   tindexes(1:1,1:2),tc(1:2),tn_test)
            elseif (ModML%ndim(v).eq.2) then
              call cinterp(ModelGrid(v), (/ obsX(linindex),obsY(linindex) /), &
                   tindexes(1:2,1:4),tc(1:4),tn_test)
            elseif (ModML%ndim(v).eq.3) then
              call cinterp(ModelGrid(v), (/ obsX(linindex),obsY(linindex),obsZ(linindex) /), &
                   tindexes(1:3,1:8),tc(1:8),tn_test)
            elseif (ModML%ndim(v).eq.4) then
              call cinterp(ModelGrid(v), (/ obsX(linindex),obsY(linindex),obsZ(linindex),obsT(linindex) /), &
                   tindexes,tc,tn_test)
            else
              write(stderr,*) 'more than 4 dimensions are not supported'
              ERROR_STOP
            end if


            if (tn_test.ne.0) then
              ! ok variable v is a candidate
              minres = hres(v)
              tn = tn_test
              tmpHindex(6,nz+1:nz+tn) = v
              tmpHindex(7:10,nz+1:nz+tn) = tindexes(:,1:tn)
              tmpHcoeff(nz+1:nz+tn) = tc(1:tn)
              ingrid = .true.
            end if
          end if
        end do

        if (v.eq.-1) then
          ! unknown variable
          tn=1

          tmpHindex(6,nz+1:nz+tn) = -1
          tmpHindex(7:10,nz+1:nz+tn) = 0
          tmpHcoeff(nz+1:nz+tn) = 0
        elseif (.not.ingrid) then
          ! out of domain
          tn=1
          
          tmpHindex(6,nz+1:nz+tn) = v
          tmpHindex(7:10,nz+1:nz+tn) = -1
          tmpHcoeff(nz+1:nz+tn) = 0
        end if

        ! observation part

        tmpHindex(1,nz+1:nz+tn) = m
        tmpHindex(2,nz+1:nz+tn) = i
        tmpHindex(3,nz+1:nz+tn) = j
        tmpHindex(4,nz+1:nz+tn) = k
        tmpHindex(5,nz+1:nz+tn) = n
        !                write(stdout,*) 'n,tc ',n,(tc(1:n))

        nz = nz+tn

#       ifdef DEBUG
            !             if (any(tmpHindex(5,1:nz).eq.0)) then
            !               write(stdout,*) 'genObservationOper: ERROR: '
            !               write(stdout,*) 'i,j,k,nz ',i,j,k,nz
            !               call flush(stdout,istat)
            !             end if

        if (nz.gt.size(tmpHcoeff)) then
          write(stderr,*) 'genObservationOper: ERROR: ', &
               'buffer variable too small!!! '
          call flush(stderr,istat)
        end if
#       endif

    end do
  end do


  allocate(Hindex(10,nz),Hcoeff(nz))
  do j=1,nz
    do i=1,10
      Hindex(i,j) = tmpHindex(i,j)
    end do
    Hcoeff(j) = tmpHcoeff(j)
  end do
!  Hindex = tmpHindex(:,1:nz)
!  Hcoeff = tmpHcoeff(1:nz)

!  write(stddebug,*) 'count(Hindex(1,:).eq.0) ',count(Hindex(1,:).eq.0)
!  write(stddebug,*) 'count(Hindex(5,:).eq.0) ',count(Hindex(5,:).eq.0)

  deallocate(tmpHindex,tmpHcoeff,varNames,obsX,obsY,obsZ)
 end subroutine genObservationOper



 !_______________________________________________________
 !
 ! load the observation correlation matrix
 ! 
 !_______________________________________________________
 !



 subroutine genObservationCorr(ntime,ObsML,Rindex,Rcoeff)
  use initfile
!  use grids
  use ufileformat
  implicit none
  integer, intent(in)  :: ntime
  type(MemLayout),intent(in) :: ObsML

  integer, pointer     :: Rindex(:,:)
  real, pointer        :: Rcoeff(:)

  integer, pointer     :: tmpRindex(:,:)
  real, pointer        :: tmpRcoeff(:)


  real, pointer        :: x(:),y(:),z(:)
  integer :: m,i,j,i1,j1,k1,i2,j2,k2, linindex1, linindex2,nz,nzmax, &
       status
  integer :: istat
  character(len=maxLen)   :: prefix

  real :: corrlen, logcorr, minlogcorr, alphax,alphay,alphaz
  real :: mincorr = 1e-3

  real, parameter :: pi = 3.141592654
  real, parameter :: EarthRadius = 6378137 ! m

  integer :: maxdisti = 10
  integer :: maxdistj = 10
  integer :: maxdistk = 30

# ifdef PROFILE
  real(8) :: cputime(10)
# endif

# ifdef PROFILE
  call cpu_time(cputime(1))
# endif

  !write(prefix,'(A,I3.3,A)') 'Obs',ntime,'.'
  call fmtIndex('Obs',ntime,'.',prefix)

  allocate(x(ObsML%effsize),y(ObsML%effsize),z(ObsML%effsize))

  call loadVector(trim(prefix)//'gridX',ObsML,x)
  call loadVector(trim(prefix)//'gridY',ObsML,y)
  call loadVector(trim(prefix)//'gridZ',ObsML,z)
  call getInitValue(initfname,trim(prefix)//'CorrLen',corrlen)

  alphax =  -(pi*EarthRadius*cos(44.*pi/180.)/180./corrlen)**2
  alphay =  -(pi*EarthRadius/180/corrlen)**2
  alphaz =  0

  minlogcorr = log(mincorr)
  write(stdout,*) 'alphax,alphay,alphaz ',alphax,alphay,alphaz

  status = 0
  nz = 0

  nzmax = 0
  do m=1,ObsML%nvar
    nzmax = nzmax + ObsML%varsizesea(m) * &
         min((2*maxdisti+1),ObsML%varshape(1,m)) * &
         min((2*maxdistj+1),ObsML%varshape(2,m)) * &
         min((2*maxdistk+1),ObsML%varshape(3,m))
  end do

  write(stdout,*) 'nz_max ',nzmax

#   ifdef PROFILE
  call cpu_time(cputime(2))
#   endif

  write(stdout,*) 'Allocating ',(4*9*nzmax)/1000000,' MB. festhalten !!!'
  allocate(tmpRindex(8,nzmax),tmpRcoeff(nzmax))

#   ifdef PROFILE
  call cpu_time(cputime(3))
#   endif

  do m=1,ObsML%nvar
    write(stdout,*) 'var shape ',ObsML%varshape(1:3,m)

    do k2=1,ObsML%varshape(3,m)
!!$     do j2=50,50
!!$       do i2=50,50
      do j2=1,ObsML%varshape(2,m)
        do i2=1,ObsML%varshape(1,m)
          linindex2 = ObsML%StartIndex(m) + i2-1 + ObsML%varshape(1,m) * (j2-1 + ObsML%varshape(2,m) * (k2-1))         
          j = ObsML%SeaIndex(linindex2)

          if (j.ne.-1) then
            do k1=1,ObsML%varshape(3,m)
              do j1=max(1,j2-maxdistj),min(ObsML%varshape(2,m),j2+maxdistj)
                do i1=max(1,i2-maxdisti),min(ObsML%varshape(1,m),i2+maxdisti)

                  linindex1 = ObsML%StartIndex(m) + i1-1 + ObsML%varshape(1,m) * (j1-1 + ObsML%varshape(2,m) * (k1-1))         
                  i = ObsML%SeaIndex(linindex1)

                  if (i.ne.-1) then

                    logcorr =  alphax * (x(i)-x(j))**2 + alphay * (y(i)-y(j))**2 + alphaz * (z(i)-z(j))**2

                    if (logcorr.ge.minlogcorr) then
#                      ifdef DEBUG
                      if (nz.ge.size(tmpRcoeff)) then
                        write(stderr,*) 'genObservationOper: ERROR: ', &
                             'buffer variable too small!!! '
                        call flush(stderr,istat)
                      end if
#                      endif

                      nz = nz+1
                      tmpRindex(1,nz) = m
                      tmpRindex(2,nz) = i1
                      tmpRindex(3,nz) = j1
                      tmpRindex(4,nz) = k1

                      tmpRindex(5,nz) = m
                      tmpRindex(6,nz) = i2
                      tmpRindex(7,nz) = j2
                      tmpRindex(8,nz) = k2

                      tmpRcoeff(nz) = exp(logcorr)
                    end if
                  end if
                end do
              end do
            end do
          end if
        end do
      end do
    end do
  end do

#   ifdef PROFILE
  call cpu_time(cputime(4))
#   endif

  write(stdout,*) 'nz ',nz
  write(stdout,*) 'dens ',sqrt(1.*nz/ObsML%effsize),ObsML%effsize
  write(stdout,*) 'nz_max ',(2*maxdisti+1)*(2*maxdistj+1)*ObsML%effsize

  allocate(Rindex(8,nz),Rcoeff(nz))
  Rindex = tmpRindex(:,1:nz)
  Rcoeff = tmpRcoeff(1:nz)

  deallocate(tmpRindex,tmpRcoeff,x,y,z)

  write(stdout,*) 'Rcoeff ',maxval(Rcoeff),minval(Rcoeff)


#   ifdef PROFILE
  call cpu_time(cputime(5))

  write(stddebug,*) 'Profiling: genObservationCorr'
  write(stddebug,*) 'load data  ',cputime(2)-cputime(1)
  write(stddebug,*) 'allocation ',cputime(3)-cputime(2)
  write(stddebug,*) 'main loop  ',cputime(4)-cputime(3)
  write(stddebug,*) 'copy data  ',cputime(5)-cputime(4)
  call flush(stddebug,istat)
#   endif

 end subroutine genObservationCorr


 !_______________________________________________________

 function obsoper(H,xf) result(Hx)
  use matoper
# ifdef ASSIM_PARALLEL
  use parall
# endif

  implicit none
  type(SparseMatrix), intent(in) :: H
  real, intent(in) :: xf(:)
  real             :: Hx(H%m)

!#define EXACT_OBS_OPER

#ifdef ASSIM_PARALLEL
  integer          :: ierr
#ifdef EXACT_OBS_OPER
  real, allocatable :: xt(:)
#else
  real             :: tmp(H%m)
  integer          :: baseIndex,j1,j2,k
#endif
#endif


#ifndef ASSIM_PARALLEL
  Hx = H.x.xf
#else

  if (size(xf) == StateVectorSizeSea) then
    ! we have already the complete state vector

    Hx = H.x.xf
    return
  end if

#ifdef EXACT_OBS_OPER
  if (procnum.eq.1) allocate(xt(H%n))
  call parallGather(xf,xt,startIndexZones,endIndexZones)
  
  if (procnum == 1) then
    Hx = H.x.xt
    deallocate(xt)
  end if
  call mpi_bcast(Hx,H%m,DEFAULT_REAL,0,comm,ierr)
#else
  tmp = 0

  j1 = startIndexZones(startZIndex(procnum))
  j2 =   endIndexZones(  endZIndex(procnum))

  baseIndex = -j1+1

  do k=1,H%nz
    if (j1 <= H%j(k) .and. H%j(k) <= j2) then
      tmp(H%i(k)) = tmp(H%i(k)) + H%s(k) * xf(H%j(k) + baseIndex)
    end if
  end do

  call mpi_allreduce(tmp, Hx, H%m, DEFAULT_REAL,mpi_sum, comm, ierr)
#endif
#endif

 end function obsoper

!_______________________________________________________
!
! error standard deviation
!_______________________________________________________
!

function stddev(S) result(R)
 implicit none
 real, intent(in) :: S(:,:)

 real :: R(size(S,1))
 integer :: i

 R = 0

 do i=1,size(S,2)
   R = R + S(:,i)**2
 end do

 R = sqrt(R)
end function 


 !_______________________________________________________
 !

 !_______________________________________________________
 !
 ! load the oservations and make the assimilation
 ! save diagnostics if the corresponding entry in the
 ! initfile exists
 !
 ! variables used:
 !
 ! ntime: time index of the observation to assimilate as defined 
 !   by the initilisation file
 ! xfp: forecasted statevector
 ! Sf: forecasted error space (error modes or ensemble)
 ! xap: analysed statevector
 ! Sa: analysed error space (error modes or ensemble)
 !
 ! if xfp and xfp are not present, then Sf is assumed to be an ensemble
 ! xfp is then computed by calculating the mean of all culumns of Sf.
 ! Sa will also be an ensemble in this case 
 !
 ! infix: "XXX." where XXX three digit time index
 ! path: path for diagnostics
 ! output: xa, Sa, biasa (global variable)
 ! str: filename for some diagnostics
 ! obsnames: descriptive names of the observations. For example: SST, SSH, 
 !   CTD_TEM, ...
 ! m: number of scalar observations
 ! n: state vector size
 ! k: error space size
 ! v: index for loop over all observation variables
 ! i1,i2: start- and end-index for a given observation variables
 ! ingrid: number of valid observation inside the model grid, also german firstame
 ! istat: return code of flush
 ! H: obervation operator (sparse matrix)
 ! C: correlation matrix (sparse matrix)
 ! ObsML: memory layout structure of the observations
 ! yo: observation vector
 ! Hxf: model forecast at observation location
 ! Hxa: model analysis at observation location
 ! invsqrtR: inverse of the square root of the diagonal elements of R (error
 !   covariance of observations)
 ! yo_Hxf: difference between yo and Hxf
 ! yo_Hxa: difference between yo and Hxa
 ! Hshift: constant shift of the observation operator
 ! Hbf: model forecasted bias at observation location
 ! HSf: error space of the model forecast at observation location
 ! HSa: error space of the model analysis at observation location
 ! amplitudes: correction (xa-xf) projected on the error space Sf
 ! mjd: modfied julian day of observation ntime
 ! cputime: cpu time for profiling
 ! bindex: index of each timed block 
 !

 subroutine Assim(ntime,Sf,Sa,xfp,xap,Efanam,Eaanam,weightf,weighta)
  use matoper
  use rrsqrt
  use ufileformat
  use initfile
#ifdef ASSIM_PARALLEL
  use parall
#endif
  implicit none
  integer, intent(in)            :: ntime
  real, intent(inout)            :: Sf(:,:)
  real, intent(out)              :: Sa(:,:)
  real, intent(inout), optional, target  :: xfp(:)
  real, intent(out),   optional, target  :: xap(:)
  real, intent(out),   optional  :: Efanam(:,:),Eaanam(:,:)
  real, intent(in),    optional  :: weightf(:)
  real, intent(out),   optional  :: weighta(:)

  ! local variables
  character(len=256)             :: prefix,path, str, str2
  character(len=256)             :: infix
  character(len=256), pointer    :: obsnames(:)
  real, pointer :: xf(:),xa(:)
  real, pointer :: modGrid(:,:)

  integer :: m,n,k,v,i1,i2,ingrid,error
  integer :: istat
  type(SparseMatrix) :: H

  type(MemLayout) :: ObsML

  real, allocatable, dimension(:) :: yo, Hxf, Hxa, invsqrtR, &
       yo_Hxf, yo_Hxa, innov_projection, Hshift, Hbf
  real, allocatable, dimension(:,:) :: HSf, HSa, locAmplitudes


!!$  real, pointer, dimension(:) :: yo, Hxf, Hxa, invsqrtR, &
!!$       yo_Hxf, yo_Hxa, innov_projection, Hshift, Hbf
!!$  real, pointer, dimension(:,:) :: HSf, HSa

  real :: amplitudes(size(Sf,2)), innov_amplitudes(size(Sf,2)), ensampl(size(Sf,2)+1)
  real :: scaling

  real(8) :: mjd
# ifdef PROFILE
  real(8) :: cputime(10)
  integer :: bindex=1
# endif

  ! multiplicative inflation factor
  real :: inflation
  real :: valex
  

# ifdef _OPENMP
  ! shared local variables among the OpenMP threads
  save :: H,yo,invsqrtR,Hxf,Hxa,HSf,HSa, &
       yo_Hxf, yo_Hxa, innov_projection, Hshift,Hbf, mjd,error,m,n,k, &
       locAmplitudes, xf, xa
# endif

  real, allocatable :: E(:,:)

  !
  ! Initialisation
  !

!$omp master
# ifdef PROFILE
  call cpu_time(cputime(bindex)); bindex=bindex+1
# endif

  n = ModML%effsize
  !write(infix,'(I3.3,A)') ntime,'.'
  call fmtIndex('',ntime,'.',infix)


  !
  ! Observations  
  !


  call MemoryLayout('Obs'//trim(infix),ObsML,rmLPObs)

  !    call loadObservationCorr(ntime,ObsML,C)

  m = ObsML%effsize
  k = size(Sf,2)

  allocate(yo(m),invsqrtR(m),Hxf(m),Hxa(m),HSf(m,k),HSa(m,k), &
       yo_Hxf(m), yo_Hxa(m), innov_projection(m), Hshift(m))

  call loadObsTime(ntime,mjd,error)
  call loadObs(ntime,ObsML,yo,invsqrtR)    

  !
  ! load the obervation matrix. All points out of grid will have a 
  ! weight (invsqrtR) = 0
  !

  call loadObservationOper(ntime,ObsML,H,Hshift,invsqrtR)
  scaling = sqrt(real(size(Sf,2)) - ASSIM_SCALING)


  if (present(xfp)) then
    ! Sf represent error modes
    xf => xfp
    xa => xap

    Hxf = obsoper(H,xf) + Hshift

    ! observed part of error modes
    do k=1,size(Sf,2)
      HSf(:,k) = obsoper(H,Sf(:,k))
    end do
  else
    ! assume Sf is an ensemble, compute mean and ensemble anomalies
    allocate(xf(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel),& 
             xa(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel))

    ! apply observation operator to all ensemble members
    ! to get observed part of ensemble members
    do k=1,size(Sf,2)
      HSf(:,k) = obsoper(H,Sf(:,k)) + Hshift
    end do

    if (presentInitValue(initfname,'Diag'//trim(infix)//'Ef')) & 
      call saveEnsemble('Diag'//trim(infix)//'Ef',ModMLParallel,Sf)

    if (presentInitValue(initfname,'Diag'//trim(infix)//'HEf')) & 
      call saveEnsemble('Diag'//trim(infix)//'HEf',ObsML,HSf)

    ! anamorphosis transform
    do k=1,size(Sf,2)
      call anamtransform(.true.,ModML,Sf(:,k))
    end do

    xf = sum(Sf,2)/size(Sf,2)
    Hxf = sum(HSf,2)/size(Sf,2)

    ! compute errors modes Sf and HSf
    do k=1,size(Sf,2)      
       Sf(:,k) = (Sf(:,k)-xf)/scaling
       HSf(:,k) = (HSf(:,k)-Hxf)/scaling
    end do
  end if

  ! now Sf,HSf are always error modes

  yo_Hxf = yo-Hxf
  innov_amplitudes = inv(HSf.tx.HSf).x.(HSf.tx.yo_Hxf)

  innov_projection = HSf.x.innov_amplitudes


!    write(6,*) 'innov ',yo_Hxf,count(invsqrtR.ne.0)
!    write(stdlog,*) 'innov ',yo_Hxf


  ! initialisation for bias-aware assimilation

  if (biastype.eq.ErrorFractionBias) then
    allocate(Hbf(m))
    Hbf = obsoper(H,biasf) + Hshift
    !Hbf = (H.x.biasf) + Hshift
  end if

  ! initialisation for local assimilation

  if (schemetype.eq.LocalScheme) then
    ! load obsGrid{X,Y,Z,T} used by callback 
    ! subroutine selectObservations

    allocate(obsGridX(ObsML%effsize),obsGridY(ObsML%effsize))
    allocate(obsGridZ(ObsML%effsize),obsGridT(ObsML%effsize))
    allocate(locAmplitudes(size(Sf,2),size(zoneSize)))

    call loadVector('Obs'//trim(infix)//'gridX',ObsML,obsGridX)
    call loadVector('Obs'//trim(infix)//'gridY',ObsML,obsGridY)
    call loadVector('Obs'//trim(infix)//'gridZ',ObsML,obsGridZ)
    if (presentInitValue(initfname,'Obs'//trim(infix)//'gridT')) &
        call loadVector('Obs'//trim(infix)//'gridT',ObsML,obsGridT)
  end if


    if (presentInitValue(initfname,'Diag'//trim(infix)//'xf')) &
         call saveVector('Diag'//trim(infix)//'xf',ModMLParallel,xf)

    if (presentInitValue(initfname,'Diag'//trim(infix)//'test')) &
         call saveVector('Diag'//trim(infix)//'test',ModMLParallel,Sf(:,2))

    if (presentInitValue(initfname,'Diag'//trim(infix)//'Hxf')) &
         call saveVector('Diag'//trim(infix)//'Hxf',ObsML,Hxf,invsqrtR.ne.0.)

    if (presentInitValue(initfname,'Diag'//trim(infix)//'yo')) &
         call saveVector('Diag'//trim(infix)//'yo',ObsML,yo,invsqrtR.ne.0.)

    if (presentInitValue(initfname,'Diag'//trim(infix)//'yo-Hxf')) &
         call saveVector('Diag'//trim(infix)//'yo-Hxf',ObsML,yo_Hxf,invsqrtR.ne.0.)

    if (presentInitValue(initfname,'Diag'//trim(infix)//'Sf')) &
         call saveErrorSpace('Diag'//trim(infix)//'Sf',Sf)

    if (presentInitValue(initfname,'Diag'//trim(infix)//'diagHPfHT')) & 
         call saveVector('Diag'//trim(infix)//'diagHPfHT',ObsML,stddev(HSf),invsqrtR.ne.0.)

    if (presentInitValue(initfname,'Diag'//trim(infix)//'stddevHxf')) &
         call saveVector('Diag'//trim(infix)//'stddevHxf',ObsML,stddev(HSf),invsqrtR.ne.0.)

    if (presentInitValue(initfname,'Diag'//trim(infix)//'diagPf')) &
         call saveVector('Diag'//trim(infix)//'diagPf',ModMLParallel,stddev(Sf))

    if (presentInitValue(initfname,'Diag'//trim(infix)//'stddevxf')) &
         call saveVector('Diag'//trim(infix)//'stddevxf',ModMLParallel,stddev(Sf))

    if (presentInitValue(initfname,'Diag'//trim(infix)//'invsqrtR'))  &
         call saveVector('Diag'//trim(infix)//'invsqrtR',ObsML,invsqrtR)

    if (presentInitValue(initfname,'Diag'//trim(infix)//'innov_amplitudes')) then
      call getInitValue(initfname,'Diag'//trim(infix)//'path',path)
      call getInitValue(initfname,'Diag'//trim(infix)//'innov_amplitudes',str)
      call usave(trim(path)//str,innov_amplitudes,9999.)
    end if

    if (presentInitValue(initfname,'Diag'//trim(infix)//'innov_projection'))  &
         call saveVector('Diag'//trim(infix)//'innov_projection',ObsML,innov_projection,invsqrtR.ne.0)

    if (presentInitValue(initfname,'Diag'//trim(infix)//'meanSf')) then
      call saveVector('Diag'//trim(infix)//'meanSf',ModMLParallel,sum(Sf,2)/size(Sf,2))
    end if

!$omp end master
!$omp barrier


  if (runtype.ne.AssimRun) then
!$omp master
     write(stdlog,*) 'Compare model with observation: ',ntime
!$omp end master
  else
!$omp master
     write(stdlog,*) 'Assimilate observation: ',ntime
!$omp end master

     if (schemetype.eq.LocalScheme) then

        ! local assimilation        
        if (biastype.eq.ErrorFractionBias) then
           call biasedLocAnalysis(zoneSize,selectObservations, &
                biasgamma,H,Hshift,xf,biasf,Hxf,Hbf,yo,Sf,HSf,invsqrtR, &
                xa,biasa,Sa, amplitudes)
        else
           call locanalysis(zoneSize,selectObservations, &
                xf,Hxf,yo,Sf,HSf,invsqrtR, xa,Sa,locAmplitudes)
        end if
     elseif (schemetype.eq.CLocalScheme) then
        allocate(modGrid(ModML%effsize,2))
        call loadVector('Model.gridX',ModML,modGrid(:,1))
        call loadVector('Model.gridY',ModML,modGrid(:,2))

        call locanalysis2(modGrid,xf,Sf,H,yo,invsqrtR, xa,Sa)
     elseif (schemetype.eq.EWPFScheme) then
       if (.not.present(weightf) .or. .not.present(weighta)) then
         write(stdout,*) 'Error: the parameters weigthf and weighta not ', &
              ' specified for EWPFScheme'
       end if
         
       !call getInitValue(initfname,'ErrorSpace.path',path)
       !call getInitValue(initfname,'ErrorSpace.weight',str)
       !call uload(trim(path)//str,weight,valex)
       
       !allocate(weighta(size(weight,1)))
       call ewpf_analysis(xf,Sf,weightf,H,invsqrtR,yo,xa,Sa,weighta)

       if (presentInitValue(initfname,'Diag'//trim(infix)//'weightf')) then
         call getInitValue(initfname,'Diag'//trim(infix)//'path',path)
         call getInitValue(initfname,'Diag'//trim(infix)//'weightf',str)
         call usave(trim(path)//str,weightf,9999.)
       end if

       if (presentInitValue(initfname,'Diag'//trim(infix)//'weighta')) then
         call getInitValue(initfname,'Diag'//trim(infix)//'path',path)
         call getInitValue(initfname,'Diag'//trim(infix)//'weighta',str)
         call usave(trim(path)//str,weighta,9999.)
       end if
     else
!$omp master
        if (biastype.eq.ErrorFractionBias) then
           call biasedanalysis(biasgamma,xf,biasf,Hxf,Hbf,yo,Sf,HSf,invsqrtR, &
                xa,biasa,Sa,amplitudes)
        else
           if (anamorphosistype.eq.TSAnamorphosis) then
              call analysisAnamorph2(xf,Hxf,yo,Sf,HSf,invsqrtR,  &
                   anamorphosisTransform,invanamorphosisTransform, &
                   xa,Sa,ensampl)
           elseif (anamorphosistype.eq.2) then
              call ensAnalysisAnamorph2(yo,Sf,HSf,invsqrtR,  &
                   anamorphosisTransform,invanamorphosisTransform, &
                   Sa,ensampl,Efanam,Eaanam)
           else

!!! FIXME flag of ensemble
              !    write(6,*) 'innov ',yo_Hxf

              !         call ensAnalysis(Sf,H.x.Sf,yo,invsqrtR,Sa,amplitudes)
              call analysis(xf,Hxf,yo,Sf,HSf,invsqrtR, xa,Sa,amplitudes)
           end if
        end if

        !      call analysis_sparseR2(xf,Hxf,yo,Sf,HSf,invsqrtR,C, xa,Sa)

!$omp end master

      end if
!$omp barrier

!$omp master

   ! apply inflation

   call getInitValue(initfname,'inflation.mult',inflation,default = 1.)
   if (inflation /= 1) then
     Sa = inflation * Sa
   end if

   ! saturate correction

    write(stdlog,*) 'max Correction reached ', & 
    count(abs(xa-xf).gt.maxCorrection),'ntimes.'

    where (xa-maxCorrection.gt.xf) xa=xf+maxCorrection
    where (xa.lt.xf-maxCorrection) xa=xf-maxCorrection

    if (.not.present(xfp)) then
      ! error modes to ensemble
      ! Sa, HSa represent an ensemble

      do k=1,size(Sa,2)
        Sa(:,k) = xa + Sa(:,k) * scaling
      end do

      ! inverse anamorphosis transform
      do k=1,size(Sa,2)
        call anamtransform(.false.,ModML,Sa(:,k))
      end do

      ! extract observed part
      do k=1,size(Sf,2)
        HSa(:,k) = obsoper(H,Sa(:,k)) + Hshift
      end do

      if (presentInitValue(initfname,'Diag'//trim(infix)//'Ea')) & 
           call saveEnsemble('Diag'//trim(infix)//'Ea',ModMLParallel,Sa)

      if (presentInitValue(initfname,'Diag'//trim(infix)//'HEa')) & 
           call saveEnsemble('Diag'//trim(infix)//'HEa',ObsML,HSa)


      xa = sum(Sa,2)/size(Sa,2)
      Hxa = sum(HSa,2)/size(Sa,2)

      ! ensemble to error modes

      xa = sum(Sa,2)/size(Sa,2)

      do k=1,size(Sa,2)
        Sa(:,k) = (Sa(:,k)-xa)/scaling
        HSa(:,k) = (HSa(:,k)-Hxa)/scaling
      end do
    else
      Hxa = obsoper(H,xa) + Hshift

      do k=1,size(Sf,2)
        HSa(:,k) = obsoper(H,Sa(:,k))
      end do
    end if

    yo_Hxa = yo-Hxa

  ! all calculations end here
  ! only optional diagonistics follow
  !
  ! begin diagonistics
  !


  !
  ! write some information in the log file
  ! rms take only into accound points inside the grid
  !

    ingrid = count(invsqrtR.ne.0.)

    write(stdlog,*) 'Nb_observations ',size(yo)
    write(stdlog,*) 'Nb_rejected_observations ',count(invsqrtR.eq.0.)
    write(stdlog,*) 'Nb_valid_observations ',ingrid
!    write(stdlog,*) 'amplitudes: ',amplitudes
!    write(stdlog,*) 'ensamplitudes: ',ensampl


    call report(stdlog,trim(infix)//'forecast.',mjd,ingrid,invsqrtR,HSf,yo_Hxf)
    if (runtype.eq.AssimRun) then
      call report(stdlog,trim(infix)//'analysis.',mjd,ingrid,invsqrtR,HSa,yo_Hxa)
    end if

    ! obsnames = name of each variable (only for output)

    if (presentInitValue(initfname,'Obs'//trim(infix)//'names')) then
      call getInitValue(initfname,'Obs'//trim(infix)//'names',obsnames)
    else
      ! default names for log VarXX
      allocate(obsnames(ObsML%nvar))
      do v=1,ObsML%nvar
        write(obsnames(v),'(A,I2.2)') 'Var',v
      end do
    end if

    write(stdlog,*) 'Per variables :'

    do v=1,ObsML%nvar
      !write(prefix,'(I3.3,A,A,A)') ntime,'.',trim(obsnames(v)),'.'
      call fmtIndex('',ntime,'.' // trim(obsnames(v)) // '.' ,prefix)


      i1 = ObsML%startIndexSea(v)
      i2 = ObsML%endIndexSea(v)
      ingrid = count(invsqrtR(i1:i2).ne.0.)

      write(stdlog,*) 'Variable number ',v,trim(obsnames(v))
      write(stdlog,*) '  Shape: ',ObsML%varshape(1:ObsML%ndim(v),v)
      write(stdlog,*) '  Size: ',ObsML%varsize(v)
      write(stdlog,*) '  Sea points: ',ObsML%varsizesea(v)
      write(stdlog,*) '  Sea points out of grid: ',count(invsqrtR(i1:i2).eq.0.)

      if (ObsML%varsizesea(v) > 0) then
        call report(stdlog,trim(prefix)//'forecast.',mjd,ingrid,invsqrtR(i1:i2),HSf(i1:i2,:),yo_Hxf(i1:i2))

        if (runtype.eq.AssimRun) then
          call report(stdlog,trim(prefix)//'analysis.',mjd,ingrid,invsqrtR(i1:i2),HSa(i1:i2,:),yo_Hxa(i1:i2))
        end if
      end if

      write(stdlog,*) 
    end do

    call flush(stdlog,istat)

    ! Save Diagnostics if desired


    ! these diagnostics makes only sens when we assimilate

    if (runtype.eq.AssimRun) then

      if (presentInitValue(initfname,'Diag'//trim(infix)//'amplitudes')) then
        call getInitValue(initfname,'Diag'//trim(infix)//'path',path)
        call getInitValue(initfname,'Diag'//trim(infix)//'amplitudes',str)
        
        if (schemetype.eq.LocalScheme) then

#         ifdef ASSIM_PARALLEL
          do k=1,size(LocAmplitudes,1)
            call parallGather(LocAmplitudes(k,startZIndex(procnum):endZIndex(procnum)),LocAmplitudes(k,:),      & 
                 startIndexZones,endIndexZones,1)
          end do
#         endif

          if (procnum == 1) then
            call usave(trim(path)//str,LocAmplitudes,9999.)
          end if
        else 
          call usave(trim(path)//str,amplitudes,9999.)
        end if
      end if


      if (presentInitValue(initfname,'Diag'//trim(infix)//'meanSa')) then
        call saveVector('Diag'//trim(infix)//'meanSa',ModMLParallel,sum(Sa,2)/size(Sa,2))
      end if

      if (presentInitValue(initfname,'Diag'//trim(infix)//'Hxa'))         &
           call saveVector('Diag'//trim(infix)//'Hxa',ObsML,Hxa,invsqrtR.ne.0.)

      if (presentInitValue(initfname,'Diag'//trim(infix)//'Hxa-Hxf'))     &
           call saveVector('Diag'//trim(infix)//'Hxa-Hxf',ObsML,Hxa-Hxf,invsqrtR.ne.0.)

      if (presentInitValue(initfname,'Diag'//trim(infix)//'yo-Hxa')) &
           call saveVector('Diag'//trim(infix)//'yo-Hxa',ObsML,yo_Hxa,invsqrtR.ne.0.)

      if (presentInitValue(initfname,'Diag'//trim(infix)//'xa'))  &
           call saveVector('Diag'//trim(infix)//'xa',ModMLParallel,xa)

      if (presentInitValue(initfname,'Diag'//trim(infix)//'xa-xf')) &
           call saveVector('Diag'//trim(infix)//'xa-xf',ModMLParallel,xa-xf)

      if (presentInitValue(initfname,'Diag'//trim(infix)//'diagPa')) &
           call saveVector('Diag'//trim(infix)//'diagPa',ModMLParallel,stddev(Sa))

      if (presentInitValue(initfname,'Diag'//trim(infix)//'stddevxa')) &
           call saveVector('Diag'//trim(infix)//'stddevxa',ModMLParallel,stddev(Sa))

      if (presentInitValue(initfname,'Diag'//trim(infix)//'diagHPaHT')) &
           call saveVector('Diag'//trim(infix)//'diagHPaHT',ObsML,stddev(HSa),invsqrtR.ne.0.)

      if (presentInitValue(initfname,'Diag'//trim(infix)//'stddevHxa')) &
           call saveVector('Diag'//trim(infix)//'stddevHxa',ObsML,stddev(HSa),invsqrtR.ne.0.)

      if (presentInitValue(initfname,'Diag'//trim(infix)//'Sa')) &
           call saveErrorSpace('Diag'//trim(infix)//'Sa',Sa)

    end if

    if (biastype.eq.ErrorFractionBias) then
      ! output bias estimation

      if (presentInitValue(initfname,'Diag'//trim(infix)//'biasf'))  &
           call saveVector('Diag'//trim(infix)//'biasf',ModMLParallel,biasf)

      if (presentInitValue(initfname,'Diag'//trim(infix)//'biasa'))  &
           call saveVector('Diag'//trim(infix)//'biasa',ModMLParallel,biasa)
    end if


#   ifdef GZIPDiag
    call getInitValue(initfname,'Diag'//trim(infix)//'path',path)
    write(stddebug,*) 'gzip ',trim(path)
    call system('gzip -f '//trim(path)//'*')
    write(stddebug,*) 'end gzip ',trim(path)
#   endif

    deallocate(obsnames)

  !
  ! end diagonistics
  !

!$omp end master

  end if

!$omp master
 
#  ifdef PROFILE
   call cpu_time(cputime(5))
   write(stddebug,*) 'profiling ',cputime(1:5)
   write(stddebug,*) 'profiling ',cputime(2:5)-cputime(1:4)
#  endif

#  ifdef DEBUG
   write(stddebug,*) 'free memory...'
#  endif

  ! free memory

  deallocate(yo,invsqrtR,Hxf,Hxa,HSf,HSa,yo_Hxf,yo_Hxa,innov_projection,H%i,H%j,H%s,Hshift)
  if (biastype.eq.ErrorFractionBias) deallocate(Hbf)
  if (schemetype.eq.LocalScheme) then
    deallocate(obsGridX,obsGridY,locAmplitudes)
    deallocate(obsGridZ,obsGridT)
  end if

  call MemoryLayoutDone(ObsML)

  if (.not.present(xfp)) then
    ! Sa must be an ensemble for output
    do k=1,size(Sa,2)
      Sa(:,k) = xa + Sa(:,k) * sqrt(real(size(Sa,2)) - ASSIM_SCALING)
    end do

    deallocate(xf,xa)
  end if

# ifdef DEBUG
  write(stddebug,*) 'Exit sub assim'
  call flush(stddebug,istat)
# endif

!$omp end master

 end subroutine Assim



 !_______________________________________________________
 !
 ! write some basic statistics 
 !_______________________________________________________
 !

 subroutine report(unit,prefix,mjd,ingrid,invsqrtR,HS,yo_Hx)
  use rrsqrt
  use matoper
  implicit none

  integer, intent(in) :: unit
  character(len=*), intent(in) :: prefix
  integer, intent(in) :: ingrid
  real(8), intent(in) :: mjd
  real, intent(in) :: invsqrtR(:),HS(:,:),yo_Hx(:)

  character(len=*), parameter :: form = '(A,A,E15.7)'
  real :: innov_projection(size(HS,2))
  real :: MahalanobisLen

!  MahalanobisLen = MahalanobisLength(yo_Hx,HS,invsqrtR)
  MahalanobisLen=0

  write(unit,'(A,A,F15.4)') prefix,'mjd                      ',mjd
  write(unit,form) prefix,'rms_yo-Hx                ',sqrt(sum( (yo_Hx)**2,invsqrtR.ne.0.)/ingrid)
  write(unit,form) prefix,'bias_yo-Hx               ',sum(yo_Hx,invsqrtR.ne.0.)/ingrid

  innov_projection = HS.tx.(invsqrtR**2*(yo_Hx))

  ! [(yo-Hx)^T R^-1 (yo-Hx)]^(1/2)

  write(unit,form) prefix,'rms_invsqrtR_yo-Hx       ',sqrt(sum( (invsqrtR*(yo_Hx))**2)/ingrid)
  write(unit,form) prefix,'bias_invsqrtR_yo-Hx      ',sum(invsqrtR*(yo_Hx))/ingrid
  !write(unit,*) prefix,'projection_coeff         ',innov_projection
  write(unit,form) prefix,'projection_yo-Hx_into_HSf',sqrt(sum(innov_projection**2)/ingrid)


  ! [(yo-Hx)^T (H P H^T + R)^-1 (yo-Hx)]^(1/2) 
  !write(unit,form) prefix,'MahalanobisLength_yo-Hx  ', MahalanobisLen

  !write(unit,form) prefix,'chi2_yo-Hx               ', MahalanobisLen**2

 end subroutine report



 ! distance between point p0 and p1
 ! p0 and p1 are (/ longitude,latitude,... /)

 real function distance(p0,p1)
  implicit none
  real, intent(in) :: p0(:), p1(:)
  real :: a, b, C
  real :: coeff
  real, parameter :: EarthRadius = 6378137 ! m
  real, parameter :: pi = 3.141592653589793238462643383279502884197
  real :: d2r = pi/180.
  
  if (metrictype == CartesianMetric) then
    distance = sqrt(sum((p1 - p0)**2))
    
  elseif (metrictype == SphericalMetricApprox) then

     coeff = pi*EarthRadius/(180.)     
     distance = sqrt((coeff * cos((p0(2)+p1(2))* (pi/360.))*(p1(1)-p0(1)))**2 &
          +(coeff * (p1(2)-p0(2)))**2)

  elseif (metrictype == SphericalMetric) then
     a = p0(2) * d2r
     b = p1(2) * d2r
     C = (p1(1) - p0(1)) * d2r
     
     coeff = sin(b) * sin(a) + cos(b) * cos(a) * cos(C)
     coeff = max(min(coeff,1.),-1.)
     ! distance in radian
     distance = acos(coeff)

     ! distance in km
     distance = EarthRadius * distance

  else
    write(stderr,*) 'Unsupported metric: ',metrictype
    write(stderr,*) 'Supported metrics are CartesianMetric (0) and SphericalMetric (1)'
    ERROR_STOP
  end if

 end function distance



 !_______________________________________________________
 !
 ! horizontal correlation function used by locanalysis
 ! 
 !_______________________________________________________


   subroutine selectObservations(ind,weight,relevantObs)
   ! input:
   !   ind:  index of zone
   ! output:
   !   weigth(:): weight of observation between 0 (no weight) and 1 (full 
   !      weight)
   !   relevantObs(:): true of observation should be used or false otherwise
   implicit none

! index of model state vector

     integer, intent(in) :: ind
     real, intent(out) :: weight(:)
!     logical, intent(out) :: relevantObs(size(weight))
     logical, intent(out) :: relevantObs(:)

     integer v,i,j,k,n,index,l

! x,y=longitude and latitude of the element the "index"th component of 
! model state vector
     real x,y,z,t,x4(4),x3(3),x2(2),x1(1)
     logical out

     real, parameter :: pi = 3.141592653589793238462643383279502884197
     real, parameter :: EarthRadius = 6378137 ! m

     logical :: noRelevantObs


     ! is state vector permuted ? yes -> in local analysis
     index = zoneIndex(ind)

     call ind2sub(ModML,index,v,i,j,k,n)

     if (ModML%ndim(v).eq.1) then
       x1 = getCoord(ModelGrid(v),(/ i /),out)
       x = x1(1)
       y = 0
     elseif (ModML%ndim(v).eq.2) then
       x2 = getCoord(ModelGrid(v),(/ i,j /),out)
       x = x2(1)
       y = x2(2)
     elseif (ModML%ndim(v).eq.3) then
       x3 = getCoord(ModelGrid(v),(/ i,j,k /),out)
       x = x3(1)
       y = x3(2)
       z = x3(3)
     elseif (ModML%ndim(v).eq.4) then
       x4 = getCoord(ModelGrid(v),(/ i,j,k,n /),out)
       x = x4(1)
       y = x4(2)
       z = x4(3)
       t = x4(4)
     else
       write(stderr,*) __FILE__,__LINE__,'the number of dimensions should be ',&
            'between 0 and 4 and got ',ModML%ndim(v)
       ERROR_STOP
     end if


!     write(6,*) 'x, y ',x,y,'-',ModML%ndim(v),index,ind,v,i,j,out

   do l = 1,size(weight)
     ! weight is the distance here

     if (loctype == 1) then
       weight(l) = distance((/ obsGridX(l),obsGridY(l) /),(/ x,y /))
     elseif (loctype == 2) then
       weight(l) = abs(obsGridZ(l) - z)
     else ! loctype == 3
       weight(l) = abs(obsGridT(l) - t)
     end if
  
     relevantObs(l) = weight(l) <= hMaxCorrLengthToObs(index)
   end do

!   write(6,*) 'relevantObs ',count(relevantObs),minval(weight),maxval(weight)

   noRelevantObs = .not.any(relevantObs)
!   if (count(relevantObs)  > 0) then
!     stop
!   end if


   if (.not.noRelevantObs) weight = exp(- (weight/hCorrLengthToObs(index))**2)



   end subroutine 


   subroutine locanalysis2(modGrid,xf,Sf,H,yo,invsqrtR, xa,Sa)
    use matoper
    use covariance
    use initfile
    real, intent(in) :: modGrid(:,:), xf(:), yo(:), invsqrtR(:)
    real, intent(in) :: Sf(:,:)
    type(SparseMatrix), intent(in) :: H
    real, intent(out) :: xa(:)
    real, intent(out), optional :: Sa(:,:)

    class(DiagCovar), allocatable :: Rc
    real, pointer :: Hc(:,:)
    real :: len
    real, pointer :: Sf2(:,:)
    integer :: leni, lenj

    integer :: indexj(ModML%effsize), nnz, i
    real :: w(ModML%effsize)
!#define CELLGRID_SEARCH
#ifdef CELLGRID_SEARCH
    ! distance for internal function lpoints_cellgrid
    real :: dist(ModML%effsize)
    type(cellgrid) :: cg
#endif

    allocate(Rc)
    call Rc%init(1./(invsqrtR**2))

    call getInitValue(initfname,'CLoc.len',len)
    call getInitValue(initfname,'CLoc.leni',leni)
    call getInitValue(initfname,'CLoc.lenj',lenj)

    write(6,*) 'len, leni, lenj ',len, leni, lenj,ModML%permute
    allocate(Hc(ModML%effsize,1))
!!! 'fix me'
    Hc = 1
    ! normalize
    do i = 1,size(Hc,2)      
      Hc(:,i) = Hc(:,i) / sqrt(sum(Hc(:,i)**2))
    end do

    allocate(Sf2(ModML%effsize,size(Sf,2)))
    Sf2 = Sf

    !    write(6,*) 'test ',modgrid(1:50,1)
    !    write(6,*) 'test ',modgrid(1:50,2)
    !    write(6,*) 'test ',modgrid(1:50,3)

    ! should give the same
    !    do i = 1,100
    !    call lpoints(i,nnz,indexj,w)
    !    write(6,*) 'sum(indexj) ',sum(indexj(1:nnz)),sum(w(1:nnz)),nnz
    !    call locpoints2(i,modgrid,len,nnz,indexj,w)
    !    write(6,*) 'sum(indexj) ',sum(indexj(1:nnz)),sum(w(1:nnz)),nnz
    ! end do
    !      stop

    !    call locensanalysis(xf,Sf2,H,yo,Rc,lpoints,Hc,xa,Sa)
    ! test

#ifdef CELLGRID_SEARCH
    cg = setupgrid(modGrid,[leni/10.,lenj/10.])
    call locensanalysis(xf,Sf2,H,yo,Rc,lpoints_cellgrid,Hc,xa)
#else
    call locensanalysis(xf,Sf2,H,yo,Rc,lpoints,Hc,xa)
#endif
    Sa = Sf


    do i = 1,size(Hc,2)      
      write(6,*) 'Hc*(xf-xa) ',i, sum(Hc(:,i) * (xf-xa))
    end do

    deallocate(Rc)
    deallocate(Hc)
    deallocate(Sf2)

   contains
    subroutine lpoints(indexi,nnz,indexj,w,onlyj)
     integer, intent(in) :: indexi
     integer, intent(out) :: nnz,indexj(:)
     real, intent(out) :: w(:)
     integer, optional, intent(in) :: onlyj(:)  

     integer :: pi,pj,pk,pn,pv,pindexi
     integer :: i,j,k,n,v
     integer :: tj
     logical :: valid
     real :: tw

     !      call locpoints(indexi,modGrid,len,nnz,indexj,w,onlyj)  

     !pindexi = indexi
     pindexi = zoneIndex(indexi)
     call ind2sub(ModML,pindexi,pv,pi,pj,pk,pn)
     nnz = 0
     !write(6,*) 'indexi ',indexi,pindexi,pv,pi,pj,pk,pn

     ! state vector is normally permuted for efficiency
     ! in this case the z-dimension vary the fastest
     ! therefore it is the most inner-loop

     do v = 1,ModML%nvar
       do n = 1,1
         do j = max(pj-lenj, 1), min(pj+lenj, ModML%varshape(2,v))
           do i = max(pi-leni, 1), min(pi+leni, ModML%varshape(1,v))
             do k = 1,ModML%varshape(3,v)
               tj = sub2ind(ModML,v,i,j,k,n,valid)
               !write(6,*) 'v,i,j,k,n ',v,i,j,k,n,tj

               if (valid) then
                 tj = invZoneIndex(tj)

                 tw = locfun(distance(modGrid(indexi,:), &
                      modGrid(tj,:))/len)

                 if (tw /= 0.) then                
                   nnz = nnz + 1
                   indexj(nnz) = tj
                   w(nnz) = tw
                 end if
               end if
             end do
           end do
         end do
       end do
     end do

     !      write(6,*) 'indexi ',indexi,nnz, maxval(w),minval(w),count(w > 0.5)
    end subroutine lpoints


    subroutine locpoints2(i,x,len,nnz,j,w)
     integer, intent(in) :: i
     real, intent(in) :: len, x(:,:)
     integer, intent(out) :: nnz, j(:)
     real, intent(out) :: w(:)

     integer :: k
     real :: dist(size(x,1)), tmp

     do k = 1,size(x,1)
       dist(k) = distance(x(i,:),x(k,:))
     end do

     nnz = 0

     do k = 1,size(x,1)
       tmp = locfun(dist(k)/len)

       if (tmp /= 0.) then
         nnz = nnz + 1
         j(nnz) = k
         w(nnz) = tmp
       end if
     end do
    end subroutine locpoints2

#ifdef CELLGRID_SEARCH

    subroutine lpoints_cellgrid(indexi,nnz,indexj,w,onlyj)
     integer, intent(in) :: indexi
     integer, intent(out) :: nnz,indexj(:)
     real, intent(out) :: w(:)
     integer, optional, intent(in) :: onlyj(:)  

     integer :: k

     call near(cg,modGrid(indexi,:),modGrid,distance,2*len,indexj,dist,nnz)
     do k = 1,nnz
       w(k) = locfun(dist(k)/len)
     end do
    end subroutine lpoints_cellgrid
#endif

   end subroutine locanalysis2

 !_______________________________________________________
 !
 ! transform objects to packed storage and back
 ! 
 !_______________________________________________________
 !

 subroutine packVector(ML,statevector,v1,v2,v3,v4,v5,v6,v7,v8)
  use matoper
  implicit none
  type(MemLayout), intent(in) ::  ML
  real, intent(in) :: v1(*)
  real, intent(in), optional :: v2(*),v3(*),v4(*),v5(*),v6(*),v7(*),v8(*)
  real, intent(out) :: StateVector(:)

  call packVariable(ML,StateVector,1,v1)
  if (present(v2)) call packVariable(ML,StateVector,2,v2)
  if (present(v3)) call packVariable(ML,StateVector,3,v3)
  if (present(v4)) call packVariable(ML,StateVector,4,v4)
  if (present(v5)) call packVariable(ML,StateVector,5,v5)
  if (present(v6)) call packVariable(ML,StateVector,6,v6)
  if (present(v7)) call packVariable(ML,StateVector,7,v7)
  if (present(v8)) call packVariable(ML,StateVector,8,v8)

  if (ML%permute) then
    call permute(zoneIndex,StateVector,StateVector)
  end if
 end subroutine packVector

 !_______________________________________________________
 !

 subroutine unpackVector(ML,statevector,v1,v2,v3,v4,v5,v6,v7,v8)
  use matoper
  implicit none
  type(MemLayout), intent(in) ::  ML
  real, intent(in) :: StateVector(:)
  real, intent(out) :: v1(*)
  real, intent(out), optional :: v2(*),v3(*),v4(*),v5(*),v6(*),v7(*),v8(*)

  real :: tmp(size(StateVector))
  integer :: i

  if (ML%permute) then
    call ipermute(zoneIndex,statevector,tmp)
  else
    tmp = statevector
  end if
  
  call unpackVariable(ML,tmp,1,v1)
  if (present(v2)) call unpackVariable(ML,tmp,2,v2)
  if (present(v3)) call unpackVariable(ML,tmp,3,v3)
  if (present(v4)) call unpackVariable(ML,tmp,4,v4)
  if (present(v5)) call unpackVariable(ML,tmp,5,v5)
  if (present(v6)) call unpackVariable(ML,tmp,6,v6)
  if (present(v7)) call unpackVariable(ML,tmp,7,v7)
  if (present(v8)) call unpackVariable(ML,tmp,8,v8)
 end subroutine unpackVector


 !_______________________________________________________
 !
 ! v1, v2, ... are variables with masked elements. Each variable
 ! represent an ensemble with size(E,2) members

 subroutine packEnsemble(ML,E,v1,v2,v3,v4,v5,v6,v7,v8)
  use matoper
  implicit none
  type(MemLayout), intent(in) ::  ML
  real, intent(in) :: v1(*)
  real, intent(in), optional :: v2(*),v3(*),v4(*),v5(*),v6(*),v7(*),v8(*)
  real, intent(out) :: E(:,:)

  ! ensemble member index
  integer :: k

  do k = 1,size(E,2)
    call packVariable(ML,E(:,k),1,v1(1+(k-1)*ML%varsize(1)))
    if (present(v2)) call packVariable(ML,E(:,k),2,v2(1+(k-1)*ML%varsize(2)))
    if (present(v3)) call packVariable(ML,E(:,k),3,v3(1+(k-1)*ML%varsize(3)))
    if (present(v4)) call packVariable(ML,E(:,k),4,v4(1+(k-1)*ML%varsize(4)))
    if (present(v5)) call packVariable(ML,E(:,k),5,v5(1+(k-1)*ML%varsize(5)))
    if (present(v6)) call packVariable(ML,E(:,k),6,v6(1+(k-1)*ML%varsize(6)))
    if (present(v7)) call packVariable(ML,E(:,k),7,v7(1+(k-1)*ML%varsize(7)))
    if (present(v8)) call packVariable(ML,E(:,k),8,v8(1+(k-1)*ML%varsize(8)))
  end do


  if (ML%permute) then
    do k = 1,size(E,2)
      call permute(zoneIndex,E(:,k),E(:,k))
    end do
  end if
  
 end subroutine packEnsemble

 !_______________________________________________________
 !

 subroutine unpackEnsemble(ML,E,v1,v2,v3,v4,v5,v6,v7,v8)
  use matoper
  implicit none
  type(MemLayout), intent(in) ::  ML
  real, intent(in) :: E(:,:)

  real, intent(out) :: v1(*)
  real, intent(out), optional :: v2(*),v3(*),v4(*),v5(*),v6(*),v7(*),v8(*)

  real :: tmp(size(E,1))

  ! ensemble member index
  integer :: k


  do k = 1,size(E,2)
    if (ML%permute) then
      call ipermute(zoneIndex,E(:,k),tmp)
    else
      tmp = E(:,k)
    end if

    call unpackVariable(ML,tmp,1,v1(1+(k-1)*ML%varsize(1)))
    if (present(v2)) call unpackVariable(ML,tmp,2,v2(1+(k-1)*ML%varsize(2)))
    if (present(v3)) call unpackVariable(ML,tmp,3,v3(1+(k-1)*ML%varsize(3)))
    if (present(v4)) call unpackVariable(ML,tmp,4,v4(1+(k-1)*ML%varsize(4)))
    if (present(v5)) call unpackVariable(ML,tmp,5,v5(1+(k-1)*ML%varsize(5)))
    if (present(v6)) call unpackVariable(ML,tmp,6,v6(1+(k-1)*ML%varsize(6)))
    if (present(v7)) call unpackVariable(ML,tmp,7,v7(1+(k-1)*ML%varsize(7)))
    if (present(v8)) call unpackVariable(ML,tmp,8,v8(1+(k-1)*ML%varsize(8)))
  end do
 end subroutine unpackEnsemble
 !_______________________________________________________
 !
 ! extract variable number "v" from state vector
 ! (without permutation)
 !_______________________________________________________
 !

 subroutine unpackVariable(ML,statevector,v,var)
  implicit none
  type(MemLayout), intent(in) ::  ML
  real, intent(in) :: StateVector(:)
  real, intent(out) :: var(*)
  integer :: v

  var(1:ML%varsize(v)) = &
       unpack(StateVector(ML%StartIndexSea(v):ML%EndIndexSea(v)), &
       ML%mask(ML%StartIndex(v):ML%EndIndex(v)) == 1,0.)

 end subroutine unpackVariable

 !_______________________________________________________
 !
 ! put a variable number "v" into the state vector
 ! (without permutation)
 !_______________________________________________________
 !

  subroutine packVariable(ML,statevector,v,var)
  implicit none
  type(MemLayout), intent(in) ::  ML
  real, intent(in) :: var(*)
  integer, intent(in) :: v
  real, intent(inout) :: StateVector(:)

  StateVector(ML%StartIndexSea(v):ML%EndIndexSea(v)) = &
       pack(var(1:ML%varsize(v)),ML%mask(ML%StartIndex(v):ML%EndIndex(v)) == 1)
 end subroutine packVariable




 !_______________________________________________________
 !

 !_______________________________________________________
 !
 ! ML1 -> destination space (e.g. ObsML)
 ! ML2 -> departure space (e.g. ModML)
 ! 
 ! 

 subroutine packSparseMatrix(Hindex,Hcoeff,ML1,ML2,H,valid1,valid2)
  use matoper
  implicit none
  integer, intent(in)  :: Hindex(:,:)
  real, intent(in) :: Hcoeff(:)
  type(MemLayout), intent(in) :: ML1,ML2
  type(SparseMatrix), intent(out)  :: H
  logical, optional, intent(out), dimension(:) :: valid1, valid2
  !    logical, intent(out), dimension(:) :: valid1, valid2

  integer :: i,nz,linindex1,linindex2
  logical :: val1, val2

  if (present(valid1)) valid1 = .true.
  if (present(valid2)) valid2 = .true.

  H%m = ML1%effsize
  H%n = ML2%effsize
  ! allocation if all points where valid
  allocate(H%i(size(Hcoeff)),H%j(size(Hcoeff)),H%s(size(Hcoeff)))

  nz = 0

  do i=1,size(Hcoeff)
    ! space 1: destination
    ! transform [Hindex(1,i) Hindex(2,i) Hindex(3,i) Hindex(4,i)] into the
    ! linear index linindex1 and trapp error in variable val1

    linindex1 = sub2ind(ML1,Hindex(1,i),Hindex(2,i),Hindex(3,i),Hindex(4,i),Hindex(5,i),val1)

    ! space 2: origin
    ! transform [Hindex(6,i) Hindex(7,i) Hindex(8,i) Hindex(9,i)] into the
    ! linear index linindex2 and trapp error in variable val2

    linindex2 = sub2ind(ML2,Hindex(6,i),Hindex(7,i),Hindex(8,i),Hindex(9,i),Hindex(10,i),val2)

    ! return the valid flags is desiered

    if (present(valid1).and.val1) valid1(linindex1) = val2
    if (present(valid2).and.val2) valid2(linindex2) = val1
    !       valid1(i) = val1
    !       valid2(i) = val2


    if (val1.and.val2) then
#ifdef DEBUG
      if (nz.ge.size(H%s)) then
        write(stderr,*) 'packSparseMatrix: ERROR: ', &
             'buffer variable too small!!! '
        call flush(stderr,istat)
      end if
#endif

      ! add an entry in sparse matrix H

      nz = nz+1
      H%i(nz) = linindex1
      H%j(nz) = linindex2
      H%s(nz) = Hcoeff(i)
    end if

  end do

  H%nz = nz

# ifdef DEBUG

  write(stddebug,*) 
  write(stddebug,'(A)') '== packSparseMatrix =='
  call writeInfoSparseMatrix(stddebug,H)
  if (present(valid1))  write(stddebug,'(A,I10)') 'count(valid1):      ',count(valid1)
  if (present(valid2))  write(stddebug,'(A,I10)') 'count(valid2):      ',count(valid2)
  call flush(stddebug,istat)

# endif

 end subroutine packSparseMatrix

 !_______________________________________________________
 !


 subroutine unpackSparseMatrix(Hindex,Hcoeff,ML1,ML2,H)
  use matoper
  implicit none
  type(SparseMatrix), intent(in)  :: H
  integer, intent(out)  :: Hindex(10,H%nz)
  real, intent(out) :: Hcoeff(H%nz)
  type(MemLayout), intent(in) :: ML1,ML2

  integer nz 

  do nz=1,H%nz
    call ind2sub(ML1,H%i(nz),Hindex(1,nz),Hindex(2,nz),Hindex(3,nz),Hindex(4,nz),Hindex(5,nz))
    call ind2sub(ML2,H%j(nz),Hindex(6,nz),Hindex(7,nz),Hindex(8,nz),Hindex(9,nz),Hindex(10,nz))
    Hcoeff(nz) = H%s(nz) 
  end do


 end subroutine unpackSparseMatrix

 !_______________________________________________________
 !
 ! assumes 
 ! - TEM and SAL as variables names
 ! - TEM and SAL of the same grid following each other in the state vector definition
 ! Example
 !                   
 ! Model.variables = ['ETA','TEM','SAL','ETA','TEM','SAL']
 ! Model.gridnum   = [    1,    1,    1,    2,    2,    2]
 ! is OK
 !
 ! but
 ! Model.variables = ['ETA','TEM','SAL','ETA','TEM','SAL']
 ! Model.gridnum   = [    1,    1,    2,    2,    2,    1]
 ! is not OK

 subroutine anamorphosisTransform(statevector)
 use ufileformat
 use anamorphosis
 implicit none
 real, intent(inout) :: statevector(:)

 real, allocatable, dimension(:,:,:) :: X,Y

 integer :: v,vT,vS

 v = 1
! return

  do while (v.lt.ModML%nvar)
   if ((ModML%varnames(v).eq.'TEM'.and.ModML%varnames(v+1).eq.'SAL').or. &
       (ModML%varnames(v).eq.'SAL'.and.ModML%varnames(v+1).eq.'TEM')) then

     if (ModML%varnames(v).eq.'TEM') then
       vT = v
       vS = v+1
     else
       vT = v+1
       vS = v
     end if

     allocate(X(ModML%varshape(1,v),ModML%varshape(2,v),ModML%varshape(3,v)), &
              Y(ModML%varshape(1,v),ModML%varshape(2,v),ModML%varshape(3,v)))


! temporaly disabled

!     call TStransform(ModelGrid3D(v)%mask,ModelGrid3D(v)%z, &
!       unpack3DVariable(ModML,statevector,vT),  &
!       unpack3DVariable(ModML,statevector,vS),  &
!       X,Y)

      call packVariable(ModML,statevector,v,X)
      call packVariable(ModML,statevector,v+1,Y)

#ifdef DEBUG
      write(stddebug,*) 'Temperature variable ',vT,' and salinity variable ',vS
      write(stddebug,*) 'tranformed by anamorphosis'
#endif

     deallocate(X,Y)
     v=v+2
   else
     v=v+1
   end if

   
  end do

 end subroutine


 !_______________________________________________________
 !

 subroutine invanamorphosisTransform(statevector)
 use anamorphosis
 use ufileformat
 implicit none
 real, intent(inout) :: statevector(:)

 real, allocatable, dimension(:,:,:) :: T,S

 integer :: v,vT,vS

 v = 1
! return

  do while (v.lt.ModML%nvar)
   if ((ModML%varnames(v).eq.'TEM'.and.ModML%varnames(v+1).eq.'SAL').or. &
       (ModML%varnames(v).eq.'SAL'.and.ModML%varnames(v+1).eq.'TEM')) then

     if (ModML%varnames(v).eq.'TEM') then
       vT = v
       vS = v+1
     else
       vT = v+1
       vS = v
     end if

     allocate(T(ModML%varshape(1,v),ModML%varshape(2,v),ModML%varshape(3,v)), &
              S(ModML%varshape(1,v),ModML%varshape(2,v),ModML%varshape(3,v)))

! disabled

!     call invTStransform(ModelGrid3D(v)%mask,ModelGrid3D(v)%z, &
!       unpack3DVariable(ModML,statevector,v),  &
!       unpack3DVariable(ModML,statevector,v+1),  &
!       T,S)

      call packVariable(ModML,statevector,vT,T)
      call packVariable(ModML,statevector,vS,S)

#ifdef DEBUG
      write(stddebug,*) 'Temperature variable ',vT,' and salinity variable ',vS
      write(stddebug,*) 'tranformed back by anamorphosis'
#endif

     deallocate(T,S)
     v=v+2
   else
     v=v+1
   end if

   
  end do

 end subroutine

! testing


 subroutine analysisAnamorph2(xf,Hxf,yo,Sf,HSf,invsqrtR,  &
  anamorph,invanamorph, &
       xa,Sa, amplitudes)
  use matoper
  use ufileformat
  use rrsqrt
  implicit none

  real, intent(in) :: xf(:),  Hxf(:), yo(:), &
       Sf(:,:), HSf(:,:), invsqrtR(:)
  real, intent(out) :: xa(:)

  interface 
    subroutine anamorph(x)
     real, intent(inout) :: x(:)
    end subroutine anamorph

    subroutine invanamorph(x)
     real, intent(inout) :: x(:)
    end subroutine invanamorph
  end interface

  integer :: i,N
  real, intent(out), optional :: Sa(:,:), amplitudes(size(Sf,2)+1)    

! transformation matix from RRSQRT to ensemble

  real :: Omega(size(Sf,2)+1,size(Sf,2))

! ensemble
  real, allocatable :: E(:,:),HEf(:,:)

  allocate(E(size(xf),size(Sf,2)+1),HEf(size(yo),size(Sf,2)+1))

  N = size(Sf,2)+1
  call sqrt2ens(xf,Sf,E,Omega)

 
! ensemble at observation locations
  HEf = spread(Hxf,2,N) + (HSf.xt.Omega)

!  call saveEnsemble('Ef1.value',ModML,E)
!  call usave('/u/abarth/Assim/Data2/Ef.u',E,0.)

  do i=1,N
    call anamorph(E(:,i))
  end do

!  call saveEnsemble('Ef2.value',ModML,E)
!  call usave('/u/abarth/Assim/Data2/Ef2.u',E,0.)

  !call ensanalysis(Ef,HEf,yo,invsqrtR,Ea, amplitudes)
  call ensanalysis(E,HEf,yo,invsqrtR,E,amplitudes)

!  call saveEnsemble('Ea2.value',ModML,E)
!  call usave('/u/abarth/Assim/Data2/Ea2.u',E,0.)

  do i=1,N
    call invanamorph(E(:,i))
  end do

!  call saveEnsemble('Ea1.value',ModML,E)
!  call usave('/u/abarth/Assim/Data2/E.u',E,0.)

  call  ens2sqrt(E,xa,Sa) 


 end subroutine

!----------------------------------------------------------------------

subroutine ensAnalysisAnamorph2(yo,Ef,HEf,invsqrtR,  &
     anamorph,invanamorph, &
     Ea, amplitudes,Efanam,Eaanam)
  use matoper
  use ufileformat
  use rrsqrt
  implicit none

  real, intent(in) :: yo(:), &
       Ef(:,:), HEf(:,:), invsqrtR(:)

  interface 
    subroutine anamorph(x)
     real, intent(inout) :: x(:)
    end subroutine anamorph

    subroutine invanamorph(x)
     real, intent(inout) :: x(:)
    end subroutine invanamorph
  end interface

  integer :: i,N
  real, intent(out), optional :: Ea(:,:), amplitudes(size(Ef,2)), Efanam(:,:), Eaanam(:,:)

! ensemble
  real, allocatable :: E(:,:),HE(:,:)

  allocate(E(size(Ef,1),size(Ef,2)),HE(size(yo),size(Ef,2)))

 
  N = size(Ef,2)



  if (present(Efanam).and.present(Eaanam)) then
! this version allows to return the anamorphosed ensemble as diagnostics  
    Efanam = Ef
    HE = HEf

    do i=1,N
      call anamorph(Efanam(:,i))
    end do

    call ensAnalysis(Efanam,HE,yo,invsqrtR,Eaanam,amplitudes)

    Ea = Eaanam 
    do i=1,N
      call invanamorph(Ea(:,i))
    end do
      
  else
! this version reduces memory usage but does not return the anamorphosed ensemble as diagnostics  
 
    Ea = Ef
    HE = HEf

    do i=1,N
      call anamorph(Ea(:,i))
    end do

!  call saveEnsemble('Ef2.value',ModML,E)
!  call usave('/u/abarth/Assim/Data2/Ef2.u',E,0.)

  !call ensanalysis(Ef,HEf,yo,invsqrtR,Ea, amplitudes)
    call ensAnalysis(Ea,HE,yo,invsqrtR,Ea,amplitudes)

!  call saveEnsemble('Ea2.value',ModML,E)
!  call usave('/u/abarth/Assim/Data2/Ea2.u',E,0.)

    do i=1,N
      call invanamorph(Ea(:,i))
    end do

!  call saveEnsemble('Ea1.value',ModML,E)
!  call usave('/u/abarth/Assim/Data2/E.u',E,0.)

   end if

 end subroutine

 subroutine anamtransform(foreward,ML,x)
  use anamorphosis
  implicit none
  logical, intent(in) :: foreward
  real, intent(inout) :: x(:)
  type(MemLayout) :: ML

  integer :: v,i,j,k,n,l,nout=0, ti,tj
  logical :: out

  if (.not.associated(AnamTrans%anam)) then
    ! do nothing
    return
  end if

  !write(0,*) 'init ',x(1),foreward,AnamTrans%anam(1)%type

  ! loop over all elements
  do l=1,size(x)    
    ! get variable index v for element l
    call ind2sub(ML,l,v,i,j,k,n)

    !write(0,*) 'init ',x(l),AnamTrans%anam(v)%type,foreward

    if (AnamTrans%anam(v)%type == 2) then
      if (foreward) then
        x(l) = log(x(l))
      else
        x(l) = exp(x(l))
      end if
    elseif (AnamTrans%anam(v)%type == 3) then
      if (foreward) then
        ti = 1
        tj = 2
      else
        ti = 2
        tj = 1
      end if

      x(l) = interp1(&
           AnamTrans%anam(v)%transform(:,ti), &
           AnamTrans%anam(v)%transform(:,tj), &
           x(l),out)

      if (out) then
        nout = nout + 1
        if (x(l) < AnamTrans%anam(v)%transform(1,ti)) then
          x(l) = AnamTrans%anam(v)%transform(1,ti)
        else
          x(l) = AnamTrans%anam(v)%transform(size(AnamTrans%anam(v)%transform,1),ti)
        end if
      end if
    end if
  end do

  !write(0,*) 'trans ',x(1)

  if (nout > 0) then
    write(stddebug,*) 'Warning: anamorphosis extrapolated ',nout,' times '
  end if
end subroutine anamtransform


subroutine ewpf_proposal_step(ntime,obsVec,dt_obs,X,weight,yo,invsqrtR,H)
 use matoper
 use sangoma_ewpf
 use initfile
 implicit none
 integer, intent(in) :: ntime                   ! current model timestep
 integer,intent(in) :: obsVec                     ! model time step at which we have next 
                                                            ! observations, i.e. next analysis time
 integer,intent(in) :: dt_obs                     ! model timesteps betwe
 real, intent(in) :: yo(:)
 real, intent(inout) :: weight(:),X(:,:)
 real, intent(in) :: invsqrtR(:)
 type(SparseMatrix), intent(in)  :: H

! current model timestep
! integer ::  ntime = 0      

 ! model time step at which we have next 
 ! observations, i.e. next analysis time 
 !integer :: obsVec = 10

 ! model timesteps between last and next
 ! observation setst
 !integer :: dt_obs = 3

 ! double precision for sangoma tools
 real(8) :: weight2(size(weight))
 real(8)  :: X2(size(X,1),size(X,2))

 real :: Qscale
 call getInitValue(initfname,'EWPF.Qscale',Qscale,default=0.001)


 !subroutine proposal_step(Ne,Nx,Ny,weight,x_n,y,ntime,obsVec,dt_obs, &
 !          cb_H, cb_HT, cb_Qhalf, cb_solve_r) bind(C, name="proposal_step_")

! write(6,*) 'weight ',__LINE__,weight
 weight2 = weight
 X2 = X

 call proposal_step(size(X,2),size(X,1),size(yo), &
      weight2,X2, &
      real(yo,8), &
      ntime,obsVec,dt_obs, &
      cb_H, cb_HT, cb_Qhalf, cb_solve_r)

 weight = weight2
 X = X2
! write(6,*) 'weight ',__LINE__,weight
! dbg(X)
contains

 subroutine cb_H(Ne,Nx,Ny,vec_in,vec_out) ! bind(C)
  use, intrinsic :: ISO_C_BINDING
  use sangoma_base, only: REALPREC, INTPREC
  implicit none
  
  integer(INTPREC), intent(in) :: Nx,Ny,Ne             ! state, observation and ensemble dimensions
  real(REALPREC), intent(in), dimension(Nx,Ne) :: vec_in      ! input vector in state space to which
  ! to apply the observation operator h, e.g. h(x)
  real(REALPREC), intent(inout), dimension(Ny,Ne) :: vec_out  ! resulting vector in observation space

  vec_out = H.x.vec_in
 end subroutine cb_H

 subroutine cb_HT(Ne,Nx,Ny,vec_in,vec_out) ! bind(C)
  use, intrinsic :: ISO_C_BINDING
  use sangoma_base, only: REALPREC, INTPREC
  implicit none
  
  integer(INTPREC), intent(in) :: Nx,Ny,Ne            ! state, observation and ensemble dimensions
  real(REALPREC), intent(in), dimension(Ny,Ne) :: vec_in     ! input vector in observation space to which
  ! to apply the observation operator h, e.g. h^T(x)
  real(REALPREC), intent(inout), dimension(Nx,Ne) :: vec_out ! resulting vector in state space

  vec_out = H.tx.vec_in
 end subroutine cb_HT


 subroutine cb_solve_r(Ne,Ny,vec_in,vec_out) ! bind(C)
  use, intrinsic :: ISO_C_BINDING
  use sangoma_base, only: REALPREC, INTPREC
  implicit none
  
  integer(INTPREC), intent(in) :: Ny,Ne               ! observation and ensemble dimensions
  real(REALPREC), intent(in), dimension(Ny,Ne) :: vec_in     ! input vector in observation space 
  ! which to apply the inverse observation error
  ! covariances R, e.g. R^{-1}(d)
  real(REALPREC), intent(inout), dimension(Ny,Ne) :: vec_out ! resulting vector in observation space
  integer :: k

  do k = 1,Ne
    vec_out(:,k) = invsqrtR**2 * vec_in(:,k)
  end do
 end subroutine cb_solve_r

 subroutine cb_Qhalf(Ne,Nx,vec_in,vec_out) ! bind(C)
  use, intrinsic :: ISO_C_BINDING
  use sangoma_base, only: REALPREC, INTPREC
  implicit none
  
  integer(INTPREC), intent(in) :: Nx,Ne               ! state and ensemble dimensions
  real(REALPREC), intent(in), dimension(Nx,Ne) :: vec_in     ! vector in state space to which to apply
                                                                     ! the squarerooted model error covariances 
                                                                     ! Q^{1/2}, e.g. Q^{1/2}(d)
  real(REALPREC), intent(inout), dimension(Nx,Ne) :: vec_out ! resulting vector in state space!!!

!  vec_out = sqrtQ.x.vec_in
  vec_out = sqrt(Qscale) * vec_in
 end subroutine cb_Qhalf

 
end subroutine ewpf_proposal_step

!------------------------------------------------------------------

subroutine ewpf_analysis(xf,Sf,weight,H,invsqrtR, &
     !sqrtQ, &
     yo,xa,Sa,weighta)
 use matoper
 use initfile
 use sangoma_ewpf
 use user_base
 implicit none
 real, intent(in) :: xf(:),Sf(:,:),weight(:),invsqrtR(:),yo(:)
 real, intent(out) :: xa(:),Sa(:,:),weighta(:)
 type(SparseMatrix), intent(in)  :: H
! type(SparseMatrix), intent(in)  :: sqrtQ

 real :: Qscale, tmp

 ! double precision for sangoma tools
 real(8), allocatable :: X(:,:)
 real(8) :: weight_analysis(size(weight))
 integer :: i

 ! parameters keep, ... are always in double precision
 tmp = keep
 call getInitValue(initfname,'EWPF.keep',tmp,default=tmp)
 keep = tmp

 tmp = nstd
 call getInitValue(initfname,'EWPF.nstd',tmp,default=tmp)
 nstd = tmp

 tmp = nmean
 call getInitValue(initfname,'EWPF.nmean',tmp,default=tmp)
 nmean = tmp

 tmp = ufac 
 call getInitValue(initfname,'EWPF.ufac',tmp,default=tmp)
 ufac = tmp

 tmp = efacNum
 call getInitValue(initfname,'EWPF.efacNum',tmp,default=tmp)
 efacNum = tmp

 tmp = freetime
 call getInitValue(initfname,'EWPF.freetime',tmp,default=tmp)
 freetime = tmp 

 tmp = nudgefac
 call getInitValue(initfname,'EWPF.nudgefac',tmp,default=tmp)
 nudgefac = tmp

 call getInitValue(initfname,'EWPF.Qscale',Qscale,default=0.001)

 allocate(X(size(xf,1),size(Sf,2)))

 do i=1,size(Sf,2)
   X(:,i) = xf + Sf(:,i)
 end do

 weight_analysis = weight

 call equal_weight_step(size(Sf,2),size(xf),size(yo), &
      weight_analysis,X,real(yo,8), &
      cb_H, cb_HT, cb_solve_r, cb_solve_hqht_plus_r, cb_Qhalf)

 weighta = weight_analysis 

! write(6,*) 'X', __LINE__,X(4,:)
 xa = sum(X,2) / size(X,2)
 do i=1,size(Sf,2)
   Sa(:,i) = X(:,i) - xa
 end do

 deallocate(X)
contains

 subroutine cb_H(Ne,Nx,Ny,vec_in,vec_out) ! bind(C)
  use, intrinsic :: ISO_C_BINDING
  use sangoma_base, only: REALPREC, INTPREC
  implicit none
  
  integer(INTPREC), intent(in) :: Nx,Ny,Ne             ! state, observation and ensemble dimensions
  real(REALPREC), intent(in), dimension(Nx,Ne) :: vec_in      ! input vector in state space to which
  ! to apply the observation operator h, e.g. h(x)
  real(REALPREC), intent(inout), dimension(Ny,Ne) :: vec_out  ! resulting vector in observation space

  vec_out = H.x.vec_in

 end subroutine cb_H

 subroutine cb_HT(Ne,Nx,Ny,vec_in,vec_out) ! bind(C)
  use, intrinsic :: ISO_C_BINDING
  use sangoma_base, only: REALPREC, INTPREC
  implicit none
  
  integer(INTPREC), intent(in) :: Nx,Ny,Ne            ! state, observation and ensemble dimensions
  real(REALPREC), intent(in), dimension(Ny,Ne) :: vec_in     ! input vector in observation space to which
  ! to apply the observation operator h, e.g. h^T(x)
  real(REALPREC), intent(inout), dimension(Nx,Ne) :: vec_out ! resulting vector in state space

  vec_out = H.tx.vec_in
 end subroutine cb_HT

 function hqht_plus_r(vec_in) result (vec_out)
  implicit none
  real, intent(in) :: vec_in(:)
  ! e.g. (HQH^T+R) vec_in
  real :: vec_out(size(vec_in))
  
  ! (HQH^T+R) vec_in
  ! H (Q (H^T * vec_in)) + R * vec_in

  ! R * vec_in
  vec_out = H.x.(Qscale  * (H.tx.vec_in))

  ! R * vec_in
  vec_out = vec_out + vec_in / (invsqrtR**2)
 end function hqht_plus_r


 subroutine cb_solve_hqht_plus_r(Ne,Ny,vec_in,vec_out) ! bind(C)
  use, intrinsic :: ISO_C_BINDING
  use sangoma_base, only: REALPREC, INTPREC
  implicit none
  
  integer(INTPREC), intent(in) :: Ny,Ne                ! observation and ensemble dimensions
  real(REALPREC), intent(in), dimension(Ny,Ne) :: vec_in      ! vector in observation space to which to
  ! apply the observation error covariances R,
  ! e.g. (HQH^T+R)^{-1}(d)
  real(REALPREC), intent(inout), dimension(Ny,Ne) :: vec_out  ! resulting vector in observation space

  integer :: k
  real :: residual(Ny),relres

  do k = 1,Ne
!    write(6,*) 'start ',k,vec_in(:,k)
!    vec_out(:,k) = hqht_plus_r(vec_in(:,k))
!    write(6,*) 'apply ',k,vec_out(:,k)

    vec_out(:,k) = pcg(hqht_plus_r,real(vec_in(:,k)),relres=relres)

!    residual = hqht_plus_r(vec_out(:,k)) - vec_in(:,k)

!    write(6,*) 'residual ', relres
!    write(6,*) 'residual ', sqrt(sum(residual**2)/sum(vec_in(:,k)**2))

    !write(6,*) 'residual ', (hqht_plus_r(vec_out(:,k)) - vec_in(:,k))
  end do

 end subroutine cb_solve_hqht_plus_R

 subroutine cb_solve_r(Ne,Ny,vec_in,vec_out) ! bind(C)
  use, intrinsic :: ISO_C_BINDING
  use sangoma_base, only: REALPREC, INTPREC
  implicit none
  
  integer(INTPREC), intent(in) :: Ny,Ne               ! observation and ensemble dimensions
  real(REALPREC), intent(in), dimension(Ny,Ne) :: vec_in     ! input vector in observation space 
  ! which to apply the inverse observation error
  ! covariances R, e.g. R^{-1}(d)
  real(REALPREC), intent(inout), dimension(Ny,Ne) :: vec_out ! resulting vector in observation space
  integer :: k

  do k = 1,Ne
    vec_out(:,k) = invsqrtR**2 * vec_in(:,k)
  end do
 end subroutine cb_solve_r

 subroutine cb_Qhalf(Ne,Nx,vec_in,vec_out) ! bind(C)
  use, intrinsic :: ISO_C_BINDING
  use sangoma_base, only: REALPREC, INTPREC
  implicit none
  
  integer(INTPREC), intent(in) :: Nx,Ne               ! state and ensemble dimensions
  real(REALPREC), intent(in), dimension(Nx,Ne) :: vec_in     ! vector in state space to which to apply
                                                                     ! the squarerooted model error covariances 
                                                                     ! Q^{1/2}, e.g. Q^{1/2}(d)
  real(REALPREC), intent(inout), dimension(Nx,Ne) :: vec_out ! resulting vector in state space!!!

  vec_out = sqrt(Qscale) * vec_in

 end subroutine cb_Qhalf


 
end subroutine ewpf_analysis



#ifdef CINTERFACE

!
! C-interface
!

 function c2fstring(cstr) result(fstr)
  use iso_c_binding
  implicit none
  character(kind=c_char) :: cstr(*)
  character(len=maxLen)            :: fstr
  integer :: len=0

  fstr = ""

  do while (cstr(len+1) /= C_NULL_CHAR)
    len = len+1
    if (len > maxLen) then
      stop 'error'
    end if
    fstr(len:len) = cstr(len)
  end do
 end function c2fstring



 subroutine  oak_init ( fname ) bind(C)
  use iso_c_binding
  implicit none
  character(kind=c_char) :: fname(*)
  call init(c2fstring(fname))
 end subroutine oak_init


 subroutine  oak_assim (ntime,n,r,Ef,Ea) bind(C)
  use iso_c_binding
  implicit none
  integer(kind=c_int), value :: ntime
  integer(kind=c_int), value :: n,r
  real(kind=c_double) :: Ef(n,r)
  real(kind=c_double) :: Ea(n,r)
  call assim(int(ntime),Ef,Ea)
 end subroutine oak_assim

 subroutine  oak_done() bind(C)
  use iso_c_binding
  implicit none
  call done()
 end subroutine oak_done

#endif
end module assimilation





