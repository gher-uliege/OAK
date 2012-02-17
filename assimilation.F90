!
!  OAK, Ocean Assimilation Kit
!  Copyright(c) 2002-2011 Alexander Barth and Luc Vandenblucke
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


module assimilation
 use ndgrid, cinterp_nd => cinterp

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
 integer, parameter :: maxLen = 256;
 
 character(len=maxLen) :: initfname, localInitfname, globalInitfname

 type(grid), allocatable :: ModelGrid(:)

 real, allocatable          :: maxCorrection(:)

 ! horizontal resolution of each variable (approximation)

 real, allocatable         :: hres(:)
 integer                   :: StateVectorSize, StateVectorSizeSea, &
      ErrorSpaceDim

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
      SphericalMetric  = 1               ! (default)


 ! value for the _FillValue attribute
 real :: FillValue = 9999.

 ! fortran unit of logfile
 ! the logfile contains simple diagnostics such as rmse with observations

 integer :: stdlog

 ! fortran unit of debugfile
 ! the debugfile contains debugging information such as 
 ! input/output, memory layout structure,... 

#ifdef DEBUG
 integer :: stddebug
#endif

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

  integer schemetype
  integer, parameter :: &
      GlobalScheme  = 0, &           ! (default)
      LocalScheme   = 1    

! partition for local assimilation

! zoneIndex: permutation vector for state vector such that all 
! zones are continuous elements


  integer, allocatable :: partition(:),zoneSize(:),zoneIndex(:),invZoneIndex(:)


! startIndex and endIndex for local zones

  integer, allocatable :: startIndexZones(:),endIndexZones(:)


! 
! variables read by callback function selectObservations
!

  real, allocatable :: obsGridX(:),obsGridY(:)
  real, allocatable :: hCorrLengthToObs(:), hMaxCorrLengthToObs(:)

 !_______________________________________________________
 !

contains

 !_______________________________________________________
 !

 subroutine init(fname)
  use initfile
  use ufileformat
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
  use anamorphosis
# ifdef ASSIM_PARALLEL
  use parall
# endif
  use matoper
  implicit none
  character(len=*), intent(in) :: fname

  integer                      :: v,vmax,n
  real                         :: c0,cx,cy,cz
  integer                      :: imax,jmax,kmax, tmpi, i,j
  character(len=MaxFNameLength), pointer   :: filenames(:),filenamesX(:),filenamesY(:),filenamesZ(:)
  character(len=MaxFNameLength)            :: path, str  
  real, pointer                :: maxCorr(:),tmp(:)
  integer                      :: NZones, zi

  initfname = fname

  call getInitValue(initfname,'runtype',runtype,default=AssimRun)
  call getInitValue(initfname,'metrictype',metrictype,default=SphericalMetric)

  call getInitValue(initfname,'Config.FillValue',FillValue,default=9999.)

  if (runtype.eq.FreeRun) return

  call getInitValue(initfname,'moderrtype',moderrtype,default=ConstModErr)
  call getInitValue(initfname,'biastype',biastype,default=NoBias)
  call getInitValue(initfname,'schemetype',schemetype,default=GlobalScheme)
  call getInitValue(initfname,'anamorphosistype',anamorphosistype,default=NoAnamorphosis)


# ifdef ASSIM_PARALLEL
  if (schemetype /= LocalScheme) then
    write(stderr,*) 'Error: for parallel version schemetype should be 1 (local assimilation)'
    ERROR_STOP
  end if
# endif

  if (anamorphosistype.ne.NoAnamorphosis) then
    call initAnamorphosis(fname)
  end if
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

  hres = 0

  call getInitValue(initfname,'Model.gridX',filenamesX)
  call getInitValue(initfname,'Model.gridY',filenamesY)
  call getInitValue(initfname,'Model.gridZ',filenamesZ)

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
        write(stderr,*) 'The dimension of variable ',trim(ModML%varnames(v)),' is ',n
        write(stderr,*) 'Error: Only 3-d grids are supported for now. '

        ERROR_STOP 


        deallocate(filenamesZ)
        call getInitValue(initfname,'Model.gridT',filenamesZ)
        call setCoord(ModelGrid(v),4,trim(path)//filenamesZ(v))


      end if
    end if
    

!   what to do with hres ?
!    hres(v) = cx**2
  end do

  deallocate(filenamesX,filenamesY,filenamesZ)



  call getInitValue(initfname,'ErrorSpace.dimension',ErrorSpaceDim,default=0)

  biasf = 0.


  ModMLParallel = ModML

! variables for local assimilation

  if (schemetype.eq.LocalScheme) then
    allocate(tmp(ModML%effsize),partition(ModML%effsize), &
         hCorrLengthToObs(ModML%effsize),hMaxCorrLengthToObs(ModML%effsize))
    call loadVector('Zones.partition',ModML,tmp)
    ! convertion real -> integer
    partition = tmp+.5
    deallocate(tmp)

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
#   endif
#else
    ModMLParallel%distributed = .false.
    ModMLParallel%startIndexParallel = 1
    ModMLParallel%endIndexParallel   =  ModMLParallel%effsize
#   endif     
  end if

  ModMLParallel%permute = schemetype.eq.LocalScheme

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

 end subroutine globalInit

 !_______________________________________________________
 !
 !_______________________________________________________
 !

 subroutine localInit(fname)
  use initfile
  use ufileformat
  implicit none
  character(len=*), intent(in) :: fname

  integer                      :: v,vmax
  real                         :: c0,cx,cy,cz
  integer                      :: imax,jmax,kmax, tmpi, i,j
  character(len=maxLen), pointer   :: filenames(:)
  character(len=maxLen)            :: path, str  
  real, pointer                :: maxCorr(:),tmp(:)

  initfname = fname

! open log file with unit stdlog

  if (presentInitValue(initfname,'logfile')) then
    call getInitValue(initfname,'logfile',str)
    stdlog = 912391
    open(stdlog,file=str,status='unknown',position='append')
  else
    stdlog = stdout
  end if

#ifdef DEBUG
! open debug file with unit stddebug

  if (presentInitValue(initfname,'debugfile')) then
    call getInitValue(initfname,'debugfile',str)
    stddebug = 92392
    open(stddebug,file=str,status='unknown',position='append')
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

 integer :: i,j,k,NZones,zi,istat


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


  !write(stdout ,*) 'shape(vector) ',procnum, shape(vector), lbound(vector), ubound(vector) 

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
  integer :: v,istat,i,j,j0,j1
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
               tmp1,FillValue,3,(/ ML%varshape(1,v),ML%varshape(2,v),ML%varshape(3,v) /),.false.)

          deallocate(tmp1)
        end do
      end if
    end if

    deallocate(xt)

  else
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
             (/ ML%varshape(1,v),ML%varshape(2,v),ML%varshape(3,v) /) ),  &
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
             (/ ML%varshape(1,v),ML%varshape(2,v),ML%varshape(3,v) /)), &
             FillValue);



      end do

      deallocate(x)
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

# ifdef DEBUG
  write(stddebug,'("== load Vector Space (",A,") ==")') str
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
    write(stddebug,*) 'remove ensemble mean and scale each member by 1/sqrt(dim)'
#   endif
  end if

! simple post processing of the ensemble
  do k=1,dim
    if (enstype.eq.2)  S(:,k) = (S(:,k)-ensembleMean)/sqrt(1.*dim)
    if (doSpaceScaling) S(:,k) = spaceScale * S(:,k)
    if (scale.ne.1) S(:,k) = scale * S(:,k)
  end do


  deallocate(formats)
  if (doSpaceScaling) deallocate(spaceScale)  
  if (enstype.eq.2)  deallocate(ensembleMean)
 end subroutine loadVectorSpace

 !_______________________________________________________
 !
 !_______________________________________________________
 !

 subroutine saveVectorSpace(str,ML,S)
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
 end subroutine saveVectorSpace


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

  call saveVectorSpace(str,ModMLParallel,S)
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

  integer :: istat,i,j
  real :: valex
  real, pointer :: Hop(:,:)
  integer, allocatable :: Hindex(:,:)
  real, allocatable :: Hcoeff(:)

  call uload(trim(path)//filename,Hop,valex)

  allocate(Hcoeff(size(Hop,2)),Hindex(8,size(Hop,2)))

! bug in compiler ?
  do j=1,size(Hop,2)
    do i=1,8
      Hindex(i,j) = Hop(i,j)
    end do
    Hcoeff(j) = Hop(9,j)
  end do
  !Hindex = Hop(:8,:)
  !Hcoeff = Hop(9,:)

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

  integer :: i,j,istat,nz,n,k
  real :: valex=9999.
  real, allocatable    :: Hop(:,:),Hcoeff(:)
  integer, allocatable :: Hindex(:,:)


  nz = H%nz
  if (present(valid1)) nz = nz + count(.not.valid1)
  if (present(valid2)) nz = nz + count(.not.valid2)

  allocate(Hindex(8,nz),Hcoeff(nz),Hop(9,nz))


  call unpackSparseMatrix(Hindex(:,:H%nz),Hcoeff(:H%nz),ML1,ML2,H)

  ! out of grid values

  k = H%nz

  if (present(valid1)) then
    do i=1,size(valid1)
      if (.not.valid1(i)) then
        k=k+1
        call ind2sub(ML1,i,Hindex(1,k),Hindex(2,k),Hindex(3,k),Hindex(4,k))
        Hindex(5,k)=-1
        Hindex(6,k)=-1
        Hindex(7,k)=-1
        Hindex(8,k)=-1
        Hcoeff(k) = 0.
      end if
    end do
  end if

  if (present(valid2)) then
    do i=1,size(valid2)
      if (.not.valid2(i)) then
        k=k+1
        Hindex(1,k)=-1
        Hindex(2,k)=-1
        Hindex(3,k)=-1
        Hindex(4,k)=-1
        call ind2sub(ML2,i,Hindex(5,k),Hindex(6,k),Hindex(7,k),Hindex(8,k))          
        Hcoeff(k) = 0.
      end if
    end do
  end if


! bug in compiler ?
  do j=1,size(Hcoeff)
    do i=1,8
      Hop(i,j) = Hindex(i,j)
    end do
    Hop(9,j) = Hcoeff(j) 
  end do
  !Hop(1:8,:) = Hindex
  !Hop(9,:)   = Hcoeff

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

  integer :: m,mmax,omax,prec,nbmots,i,j,istat
  real :: valex
  logical :: isdegen

  if (present(packed)) then
    la%removeLandPoints = packed
  else
    la%removeLandPoints = .true.
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
  allocate(la%invindex(la%totsizesea))

  la%SeaIndex = -1
  j=1
  do i=1,la%totsize
    if (.not.la%removeLandPoints.or.la%Mask(i).eq.1) then
      la%SeaIndex(i) = j
      la%invindex(j) = i
      j=j+1
    end if
  end do

  if (la%removeLandPoints) then
    la%effsize = la%totsizeSea
  else
    la%effsize = la%totsize
  end if

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
    write(stddebug,'(I2," ",A20,A20,6I10)') m,trim(filenames(m)),trim(sizeformat(la%varshape(1:la%ndim(m),m))), &
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

 function sub2ind(ML,v,i,j,k,valid) result(index)
  implicit none
  type(MemLayout) :: ML
  integer :: index
  integer, intent(in) :: v,i,j,k
  logical, optional,intent(out) :: valid

  logical :: val
  integer :: linindex

  index = -1
  val = 1 <= v.and.v <= ML%nvar

  if (val) then
    val = &
         1.le.i.and.i.le.ML%varshape(1,v).and. &
         1.le.j.and.j.le.ML%varshape(2,v).and. &
         1.le.k.and.k.le.ML%varshape(3,v)

    if (val) then
      index = ML%StartIndex(v) + i-1 + ML%varshape(1,v) * (j-1 + ML%varshape(2,v) * (k-1))
      index = ML%SeaIndex(index)
      val =  index.ne.-1
    end if


  end if

  if (present(valid)) valid=val

 end function sub2ind

 !_______________________________________________________
 !
 ! index includes only sea points

 subroutine ind2sub(ML,index,v,i,j,k)
  implicit none
  type(MemLayout) :: ML
  integer, intent(in) :: index
  integer, intent(out) :: v,i,j,k
  integer :: linindex
  linindex = ML%invindex(index)

  do v=1,ML%nvar
    if (ML%startIndex(v).le.linindex.and.linindex.le.ML%endIndex(v)) then
      linindex = linindex - ML%startIndex(v) 
      ! i,j,k zero-based
      k = linindex/(ML%varshape(1,v)*ML%varshape(2,v))
      linindex = linindex - k*(ML%varshape(1,v)*ML%varshape(2,v))
      j = linindex/ML%varshape(1,v)
      i = linindex - j*(ML%varshape(1,v))

      ! i,j,k one-based
      i=i+1; j=j+1; k=k+1;
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
 !

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

 subroutine saveStateVector_byfilenames(path,filenames,StateVector)
  use initfile
  use ufileformat
  character(len=*), intent(in) :: path
  character(len=*), intent(in) :: filenames(:)
  real, intent(in) :: StateVector(:)

  write(6,*) 'saveStateVector_byfilenames ',ModMLParallel%distributed
  call saveVector_byfilenames(path,filenames,ModMLParallel,StateVector)
 end subroutine saveStateVector_byfilenames

 !_______________________________________________________
 !

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
  integer :: day,month,year,h,min,istat
  real :: s,seconds

  write(prefix,'(A,I3.3,A)') 'Obs',ntime,'.'


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
  character(len=maxLen) :: prefix,str
  integer :: day,month,year,h,min,omax
  real :: s,seconds

  integer :: prec,imax,jmax,kmax,nbmots,istat
  real :: valex
  logical :: isdegen

  real, parameter :: min_rmse = 0.01

  write(prefix,'(A,I3.3,A)') 'Obs',ntime,'.'
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
  integer :: v,vmax,m,mmax,n,nz,istat

  integer :: error
  real :: valex

  real, pointer :: Cop(:,:)
  integer, pointer :: Cindex(:,:)
  real, pointer :: Ccoeff(:)
  logical, pointer :: validobs(:)


  write(prefix,'(A,I3.3,A)') 'Obs',ntime,'.'
  call getInitValue(initfname,trim(prefix)//'path',path,default='')

  mmax = size(ObsML%VarSize)

  call getInitValue(initfname,trim(prefix)//'Correlation',str,error)

  if (error.eq.0) then
#ifdef DEBUG
    write(stddebug,*) 'Load observation correlation: ',trim(str)
#endif
    call uload(trim(path)//str,Cop,valex)
    allocate(Ccoeff(size(Cop,2)),Cindex(8,size(Cop,2)))
    Cindex = Cop(1:8,:)
    Ccoeff = Cop(9,:)
    deallocate(Cop)
  else
#ifdef DEBUG
    write(stddebug,*) 'Generate observation correlation'
#endif
    call genObservationCorr(ntime,ObsML,Cindex,Ccoeff)
  end if

  ! save the obervation operator if desiered

  write(Dprefix,'(A,I3.3,A)') 'Diag',ntime,'.'
  if (presentInitValue(initfname,trim(Dprefix)//'Correlation')) then 
    call getInitValue(initfname,trim(Dprefix)//'Correlation',str)
    call getInitValue(initfname,trim(Dprefix)//'path',path,default='')
    allocate(Cop(9,size(Ccoeff)))
    Cop(1:8,:) = Cindex
    Cop(9,:) = Ccoeff
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
  integer :: v,vmax,m,mmax,n,nz,idummy,nbentries

  integer :: error,istat,i,j
  real :: valex
  logical :: isdegen
  real, pointer :: Hop(:,:)
  integer, pointer :: Hindex(:,:)
  real, pointer :: Hcoeff(:)
  logical, allocatable :: validobs(:)
  real, allocatable :: shiftMod(:)

  write(prefix,'(A,I3.3,A)') 'Obs',ntime,'.'
  call getInitValue(initfname,trim(prefix)//'path',path,default='')

  mmax = size(ObsML%VarSize)

  if (presentInitValue(initfname,trim(prefix)//'operator')) then

    call getInitValue(initfname,trim(prefix)//'operator',str)
#ifdef DEBUG
    write(stddebug,*) 'Load observation operator: ',trim(str)
#endif
    call uload(trim(path)//str,Hop,valex)
    allocate(Hcoeff(size(Hop,2)),Hindex(8,size(Hop,2)))
    do j=1,size(Hop,2)
      do i=1,8
        Hindex(i,j) = Hop(i,j)
      end do
      Hcoeff(j) = Hop(9,j)
    end do
    !Hindex = Hop(:8,:)
    !Hcoeff = Hop(9,:)
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

    allocate(Hcoeff(nz),Hindex(8,nz))

    nz=1

    do m=1,ObsML%nvar
      call uload(trim(path)//filenames(m),Hop,valex)
      do i=1,size(Hop,2)
        Hindex(1,nz) = m
        Hindex(2:8,nz) = floor(Hop(2:8,i)+.5)
        Hcoeff(nz) = Hop(9,i)
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

  write(Dprefix,'(A,I3.3,A)') 'Diag',ntime,'.'
  if (presentInitValue(initfname,trim(Dprefix)//'H')) then 
    call getInitValue(initfname,trim(Dprefix)//'H',str)
    call getInitValue(initfname,trim(Dprefix)//'path',path,default='')
    allocate(Hop(9,size(Hcoeff)))
    Hop(1:8,:) = Hindex
    Hop(9,:) = Hcoeff
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


  if (schemetype.eq.LocalScheme) then
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

  integer :: maxi, maxj

  maxi = maxval(H%i(1:H%nz))
  maxj = maxval(H%j(1:H%nz))

  write(unit,'(A,A10)')   'Matrix shape:       ',sizeformat((/H%m,H%n/))
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
  real, allocatable, dimension(:) :: obsX, obsY, obsZ

  character(len=maxLen)   :: prefix,str
  type(MemLayout), intent(in) :: ObsML

  integer              :: ti(8),tj(8),tk(8), &
       i,j,k, istat, &
       v,tv,vmax,m,mmax,omaxSea,n,tn,nz,linindex, &
       tindexes(3,8), tmpm
  real                 :: tc(8), valex, minres
  logical              :: isdegen

  write(prefix,'(A,I3.3,A)') 'Obs',ntime,'.'
  call getInitValue(initfname,trim(prefix)//'path',path,default='')

  mmax = size(ObsML%VarSize)
  omaxSea = ObsML%effsize


  ! assume that land-point are removed 
  ! ObsML%removeLandPoints == .true.

  allocate(tmpHindex(8,8*omaxSea),tmpHcoeff(8*omaxSea))
  nz = 0

  call getInitValue(initfname,trim(prefix)//'variables',varNames)

  ! load position of observations
  allocate(obsX(ObsML%effsize),obsY(ObsML%effsize),obsZ(ObsML%effsize))

  call loadVector(trim(prefix)//'gridX',ObsML,obsX)
  call loadVector(trim(prefix)//'gridY',ObsML,obsY)
  call loadVector(trim(prefix)//'gridZ',ObsML,obsZ)

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
        call ind2sub(ObsML,linindex,tmpm,i,j,k)

        minres = huge(minres)
        v = -1
        n = -1

        do tv=1,size(ModML%varnames)
          if (varNames(m).eq.ModML%varnames(tv).and.minres.ge.hres(tv)) then
            ! known variable
            v = tv

            ! compute interpolation coefficients
            
            tindexes=1
            if (ModML%ndim(v).eq.2) then
              call cinterp_nd(ModelGrid(v), (/ obsX(linindex),obsY(linindex) /), &
                       tindexes(1:2,1:4),tc(1:4),tn)
              ti(1:tn) =  tindexes(1,1:tn)
              tj(1:tn) =  tindexes(2,1:tn)
              tk(1:tn) = 1
            else
              call cinterp_nd(ModelGrid(v), (/ obsX(linindex),obsY(linindex),obsZ(linindex) /), &
                   tindexes,tc,tn)
              
              ti(1:tn) =  tindexes(1,1:tn)
              tj(1:tn) =  tindexes(2,1:tn)
              tk(1:tn) =  tindexes(3,1:tn)
            end if


            if (tn.ne.0) then
              ! ok variable v is a candidate
              minres = hres(v)
              n = tn
              tmpHindex(5,nz+1:nz+n) = v
              tmpHindex(6,nz+1:nz+n) = ti(1:n)
              tmpHindex(7,nz+1:nz+n) = tj(1:n)
              tmpHindex(8,nz+1:nz+n) = tk(1:n)

              tmpHcoeff(nz+1:nz+n) = tc(1:n)
            end if
          end if
        end do

        if (v.eq.-1) then
          ! unknown variable
          n=1

          tmpHindex(5,nz+1:nz+n) = -1
          tmpHindex(6,nz+1:nz+n) = 0
          tmpHindex(7,nz+1:nz+n) = 0
          tmpHindex(8,nz+1:nz+n) = 0
          tmpHcoeff(nz+1:nz+n) = 0
        elseif (n.eq.-1) then
          ! out of domain
          n=1
          
          tmpHindex(5,nz+1:nz+n) = v
          tmpHindex(6,nz+1:nz+n) = -1
          tmpHindex(7,nz+1:nz+n) = -1
          tmpHindex(8,nz+1:nz+n) = -1
          tmpHcoeff(nz+1:nz+n) = 0
        end if

        ! observation part

        tmpHindex(1,nz+1:nz+n) = m
        tmpHindex(2,nz+1:nz+n) = i
        tmpHindex(3,nz+1:nz+n) = j
        tmpHindex(4,nz+1:nz+n) = k
        !                write(stdout,*) 'n,tc ',n,(tc(1:n))

        nz = nz+n

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


  allocate(Hindex(8,nz),Hcoeff(nz))
  do j=1,nz
    do i=1,8
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
  character(len=maxLen), pointer :: varNames(:), &
       gridXnames(:),gridYnames(:),gridZnames(:)
  character(len=maxLen)    :: path

  real, pointer        :: x(:),y(:),z(:)
  real, pointer, dimension(:,:,:) :: gridX,gridY,gridZ
  integer :: m,i,j,i1,j1,k1,i2,j2,k2, linindex1, linindex2,nz,nzmax, &
       status,istat
  character(len=maxLen)   :: prefix,str

  real                 :: valex

  real :: corrlen, logcorr, minlogcorr, alpha,alphax,alphay,alphaz
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

  write(prefix,'(A,I3.3,A)') 'Obs',ntime,'.'

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
  real             :: Hx(H%m), tmp(H%m)
  integer          :: ierr,k,j1,j2,baseIndex

!#define EXACT_OBS_OPER


#ifndef ASSIM_PARALLEL

  Hx = H.x.xf

#else

#ifdef EXACT_OBS_OPER
  integer, allocatable :: rcount(:),rdispls(:)
!  real :: xt(H%n)
   real, allocatable :: xt(:)


   if (procnum.eq.1) allocate(xt(H%n))

!  write(stdout,*) 'xt ',__FILE__,__LINE__,size(xt),size(xf),procnum,H%n
!   allocate(xt(H%n))

  call parallGather(xf,xt,startIndexZones,endIndexZones)

  if (procnum == 1) then
    Hx = H.x.xt
    deallocate(xt)
  end if
!    deallocate(xt)

  call mpi_bcast(Hx,H%m,DEFAULT_REAL,0,mpi_comm_world,ierr)

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

!  write(stdout,*) ' allreduce '
  call mpi_allreduce(tmp, Hx, H%m, DEFAULT_REAL,mpi_sum, mpi_comm_world, ierr)

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
 ! xf: forecasted statevector
 ! Sf: forecasted error space
 ! xa: analysed statevector
 ! Sa: analysed error space
 ! ntime: time index of the observation to assimilate as defined 
 !   by the initilisation file
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

 subroutine Assim(ntime,xf,Sf,xa,Sa,Efanam,Eaanam)
  use matoper
  use rrsqrt
  use ufileformat
  use initfile
#ifdef ASSIM_PARALLEL
  use parall
#endif
  implicit none
  integer, intent(in) :: ntime
!  real, intent(in) :: xf(:),Sf(:,:)
! FIX ME
  real, intent(inout) :: xf(:),Sf(:,:)
  real, intent(out) :: xa(:), Sa(:,:)
  real, optional, intent(out) :: Efanam(:,:),Eaanam(:,:)

  character(len=256) :: prefix,path, str
  character(len=4) :: infix
  character(len=256), pointer :: obsnames(:)

  integer :: m,n,k,v,i1,i2,ingrid,error,istat
  type(SparseMatrix) :: H,C
  type(MemLayout) :: ObsML

  real, allocatable, dimension(:) :: yo, Hxf, Hxa, invsqrtR, &
       yo_Hxf, yo_Hxa, innov_projection, Hshift, Hbf
  real, allocatable, dimension(:,:) :: HSf, HSa, locAmplitudes

!!$  real, pointer, dimension(:) :: yo, Hxf, Hxa, invsqrtR, &
!!$       yo_Hxf, yo_Hxa, innov_projection, Hshift, Hbf
!!$  real, pointer, dimension(:,:) :: HSf, HSa

  real :: amplitudes(size(Sf,2)), innov_amplitudes(size(Sf,2)), ensampl(size(Sf,2)+1)

  real(8) :: mjd
# ifdef PROFILE
  real(8) :: cputime(10)
  integer :: bindex=1
# endif




# ifdef _OPENMP
  ! shared local variables among the OpenMP threads
  save :: H,yo,invsqrtR,Hxf,Hxa,HSf,HSa, &
       yo_Hxf, yo_Hxa, innov_projection, Hshift,Hbf, mjd,error,m,n,k, &
       locAmplitudes
# endif

  real, allocatable :: E(:,:)

  !
  ! Initialisation
  !

!$omp master
# ifdef PROFILE
  call cpu_time(cputime(bindex)); bindex=bindex+1
# endif

  write(infix,'(I3.3,A)') ntime,'.'

  call MemoryLayout('Obs'//infix,ObsML,rmLPObs)

  !    call loadObservationCorr(ntime,ObsML,C)

  m = ObsML%effsize
  n = ModML%effsize
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


  Hxf = obsoper(H,xf) + Hshift

  do k=1,size(Sf,2)
    HSf(:,k) = obsoper(H,Sf(:,k))
  end do

!  write (6,*) 'sum HSf ',sum(HSf)

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
    ! load obsGridX and obsGridY used by callback 
    ! subroutine selectObservations

    allocate(obsGridX(ObsML%effsize),obsGridY(ObsML%effsize),locAmplitudes(size(Sf,2),size(zoneSize)))
    call loadVector('Obs'//infix//'gridX',ObsML,obsGridX)
    call loadVector('Obs'//infix//'gridY',ObsML,obsGridY)

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
!$omp end master

    end if

    !      call analysis_sparseR2(xf,Hxf,yo,Sf,HSf,invsqrtR,C, xa,Sa)

!$omp master


   ! saturate correction

    write(stdlog,*) 'max Correction reached ', & 
    count(abs(xa-xf).gt.maxCorrection),'ntimes.'

    where (xa-maxCorrection.gt.xf) xa=xf+maxCorrection
    where (xa.lt.xf-maxCorrection) xa=xf-maxCorrection

    Hxa = obsoper(H,xa) + Hshift

    do k=1,size(Sf,2)
      HSa(:,k) = obsoper(H,Sa(:,k))
    end do

    yo_Hxa = yo-Hxa

!$omp end master

  end if


  ! all calculations end here
  ! only optional diagonistics follow
  !
  ! begin diagonistics
  !

  !
  ! write some information in the log file
  ! rms take only into accound points inside the grid
  !


!$omp master
! for all
! fix me
!  if (procnum.eq.1) then
  if (.true.) then
    ingrid = count(invsqrtR.ne.0.)


    write(stdlog,*) 'Nb_observations ',size(yo)
    write(stdlog,*) 'Nb_rejected_observations ',count(invsqrtR.eq.0.)
    write(stdlog,*) 'Nb_valid_observations ',ingrid
!    write(stdlog,*) 'amplitudes: ',amplitudes
!    write(stdlog,*) 'ensamplitudes: ',ensampl


    call report(stdlog,infix//'forecast.',mjd,ingrid,invsqrtR,HSf,yo_Hxf)
    if (runtype.eq.AssimRun) then
      call report(stdlog,infix//'analysis.',mjd,ingrid,invsqrtR,HSa,yo_Hxa)
    end if

    ! obsnames = name of each variable (only for output)

    if (presentInitValue(initfname,'Obs'//infix//'names')) then
      call getInitValue(initfname,'Obs'//infix//'names',obsnames)
    else
      ! default names for log VarXX
      allocate(obsnames(ObsML%nvar))
      do v=1,ObsML%nvar
        write(obsnames(v),'(A,I2.2)') 'Var',v
      end do
    end if

    write(stdlog,*) 'Per variables :'

    do v=1,ObsML%nvar
      write(prefix,'(I3.3,A,A,A)') ntime,'.',trim(obsnames(v)),'.'

      i1 = ObsML%startIndexSea(v)
      i2 = ObsML%endIndexSea(v)
      ingrid = count(invsqrtR(i1:i2).ne.0.)

      write(stdlog,*) 'Variable number ',v
      write(stdlog,*) '  Shape: ',ObsML%varshape(1:ObsML%ndim(v),v)
      write(stdlog,*) '  Size: ',ObsML%varsize(v)
      write(stdlog,*) '  Sea points: ',ObsML%varsizesea(v)
      write(stdlog,*) '  Sea points out of grid: ',count(invsqrtR(i1:i2).eq.0.)

      call report(stdlog,trim(prefix)//'forecast.',mjd,ingrid,invsqrtR(i1:i2),HSf(i1:i2,:),yo_Hxf(i1:i2))

      if (runtype.eq.AssimRun) then
        call report(stdlog,trim(prefix)//'analysis.',mjd,ingrid,invsqrtR(i1:i2),HSa(i1:i2,:),yo_Hxa(i1:i2))
      end if

      write(stdlog,*) 
    end do

    call flush(stdlog,istat)

    ! Save Diagnostics if desired

    if (presentInitValue(initfname,'Diag'//infix//'xf')) &
         call saveVector('Diag'//infix//'xf',ModMLParallel,xf)

    if (presentInitValue(initfname,'Diag'//infix//'test')) &
         call saveVector('Diag'//infix//'test',ModMLParallel,Sf(:,2))

    if (presentInitValue(initfname,'Diag'//infix//'Hxf')) &
         call saveVector('Diag'//infix//'Hxf',ObsML,Hxf,invsqrtR.ne.0.)

    if (presentInitValue(initfname,'Diag'//infix//'yo')) &
         call saveVector('Diag'//infix//'yo',ObsML,yo,invsqrtR.ne.0.)

    if (presentInitValue(initfname,'Diag'//infix//'yo-Hxf')) &
         call saveVector('Diag'//infix//'yo-Hxf',ObsML,yo_Hxf,invsqrtR.ne.0.)

    if (presentInitValue(initfname,'Diag'//infix//'Sf')) &
         call saveErrorSpace('Diag'//infix//'Sf',Sf)

    if (presentInitValue(initfname,'Diag'//infix//'diagHPfHT')) & 
         call saveVector('Diag'//infix//'diagHPfHT',ObsML,stddev(HSf),invsqrtR.ne.0.)

    if (presentInitValue(initfname,'Diag'//infix//'stddevHxf')) &
         call saveVector('Diag'//infix//'stddevHxf',ObsML,stddev(HSf),invsqrtR.ne.0.)

    if (presentInitValue(initfname,'Diag'//infix//'diagPf')) &
         call saveVector('Diag'//infix//'diagPf',ModMLParallel,stddev(Sf))

    if (presentInitValue(initfname,'Diag'//infix//'stddevxf')) &
         call saveVector('Diag'//infix//'stddevxf',ModMLParallel,stddev(Sf))

    if (presentInitValue(initfname,'Diag'//infix//'invsqrtR'))  &
         call saveVector('Diag'//infix//'invsqrtR',ObsML,invsqrtR)

    if (presentInitValue(initfname,'Diag'//infix//'innov_amplitudes')) then
      call getInitValue(initfname,'Diag'//infix//'path',path)
      call getInitValue(initfname,'Diag'//infix//'innov_amplitudes',str)
      call usave(trim(path)//str,innov_amplitudes,9999.)
    end if

    if (presentInitValue(initfname,'Diag'//infix//'innov_projection'))  &
         call saveVector('Diag'//infix//'innov_projection',ObsML,innov_projection,invsqrtR.ne.0)

    if (presentInitValue(initfname,'Diag'//infix//'meanSf')) then
      call saveVector('Diag'//infix//'meanSf',ModMLParallel,sum(Sf,2)/size(Sf,2))
    end if

    if (presentInitValue(initfname,'Diag'//infix//'Ef')) then
         allocate(E(size(xf),size(Sf,2)))

         do k=1,size(Sf,2)
           E(:,k) = xf + Sf(:,k) * sqrt(real(size(Sf,2)))
         end do

         call saveErrorSpace('Diag'//infix//'Ef',E)
         deallocate(E)
    end if

    ! these diagnostics makes only sens when we assimilate

    if (runtype.eq.AssimRun) then

      if (presentInitValue(initfname,'Diag'//infix//'amplitudes')) then
        call getInitValue(initfname,'Diag'//infix//'path',path)
        call getInitValue(initfname,'Diag'//infix//'amplitudes',str)
        if (schemetype.eq.LocalScheme) then
          call usave(trim(path)//str,LocAmplitudes,9999.)
        else 
          call usave(trim(path)//str,amplitudes,9999.)
        end if
      end if


      if (presentInitValue(initfname,'Diag'//infix//'meanSa')) then
        call saveVector('Diag'//infix//'meanSa',ModMLParallel,sum(Sa,2)/size(Sa,2))
      end if

      if (presentInitValue(initfname,'Diag'//infix//'Hxa'))         &
           call saveVector('Diag'//infix//'Hxa',ObsML,Hxa,invsqrtR.ne.0.)

      if (presentInitValue(initfname,'Diag'//infix//'Hxa-Hxf'))     &
           call saveVector('Diag'//infix//'Hxa-Hxf',ObsML,Hxa-Hxf,invsqrtR.ne.0.)

      if (presentInitValue(initfname,'Diag'//infix//'yo-Hxa')) &
           call saveVector('Diag'//infix//'yo-Hxa',ObsML,yo_Hxa,invsqrtR.ne.0.)

      if (presentInitValue(initfname,'Diag'//infix//'xa'))  &
           call saveVector('Diag'//infix//'xa',ModMLParallel,xa)

      if (presentInitValue(initfname,'Diag'//infix//'xa-xf')) &
           call saveVector('Diag'//infix//'xa-xf',ModMLParallel,xa-xf)

      if (presentInitValue(initfname,'Diag'//infix//'diagPa')) &
           call saveVector('Diag'//infix//'diagPa',ModMLParallel,stddev(Sa))

      if (presentInitValue(initfname,'Diag'//infix//'stddevxa')) &
           call saveVector('Diag'//infix//'stddevxa',ModMLParallel,stddev(Sa))

      if (presentInitValue(initfname,'Diag'//infix//'diagHPaHT')) &
           call saveVector('Diag'//infix//'diagHPaHT',ObsML,stddev(HSa),invsqrtR.ne.0.)

      if (presentInitValue(initfname,'Diag'//infix//'stddevHxa')) &
           call saveVector('Diag'//infix//'stddevHxa',ObsML,stddev(HSa),invsqrtR.ne.0.)

      if (presentInitValue(initfname,'Diag'//infix//'Sa')) &
           call saveErrorSpace('Diag'//infix//'Sa',Sa)

      if (presentInitValue(initfname,'Diag'//infix//'Ea')) then
         allocate(E(size(xf),size(Sa,2)))

         do k=1,size(Sa,2)
            E(:,k) = xa + Sa(:,k) * sqrt(real(size(Sa,2)))
         end do

         call saveErrorSpace('Diag'//infix//'Ea',E)
         deallocate(E)
      end if

    end if

    if (biastype.eq.ErrorFractionBias) then
      ! output bias estimation

      if (presentInitValue(initfname,'Diag'//infix//'biasf'))  &
           call saveVector('Diag'//infix//'biasf',ModMLParallel,biasf)

      if (presentInitValue(initfname,'Diag'//infix//'biasa'))  &
           call saveVector('Diag'//infix//'biasa',ModMLParallel,biasa)
    end if


#   ifdef GZIPDiag
    call getInitValue(initfname,'Diag'//infix//'path',path)
    write(stddebug,*) 'gzip ',trim(path)
    call system('gzip -f '//trim(path)//'*')
    write(stddebug,*) 'end gzip ',trim(path)
#   endif

    deallocate(obsnames)
  end if

  !
  ! end diagonistics
  !

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
  if (schemetype.eq.LocalScheme) deallocate(obsGridX,obsGridY,locAmplitudes)

  call MemoryLayoutDone(ObsML)

#   ifdef DEBUG

  write(stddebug,*) 'Exit sub assim'
  call flush(stddebug,istat)

#   endif

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


 !_______________________________________________________
 !
 ! horizontal correlation function used by locanalysis
 ! 
 !_______________________________________________________


   subroutine selectObservations(ind,weight,relevantObs)
   implicit none

! index of model state vector

     integer, intent(in) :: ind
     real, intent(out) :: weight(:)
!     logical, intent(out) :: relevantObs(size(weight))
     logical, intent(out) :: relevantObs(:)

     integer v,i,j,k,index,l

! x,y=longitude and latitude of the element the "index"th component of 
! model state vector
     real x,y,coeff,x3(3),x2(2)
     logical out;

     real, parameter :: pi = 3.141592653589793238462643383279502884197
     real, parameter :: EarthRadius = 6378137 ! m

     logical :: noRelevantObs


     ! is state vector permuted ? yes -> in local analysis
     index = zoneIndex(ind)

     call ind2sub(ModML,index,v,i,j,k)

     if (ModML%ndim(v).eq.2) then
      x2 = getCoord(ModelGrid(v),(/ i,j /),out);
      x = x2(1);
      y = x2(2);
    else
      x3 = getCoord(ModelGrid(v),(/ i,j,k /),out);
      x = x3(1);
      y = x3(2); 
    end if 

!    write(6,*) 'x, y ',x,y,'-',ModML%ndim(v),index,v,i,j,out

   do l = 1,size(weight)
     ! weight is the distance here
     weight(l) = distance((/ obsGridX(l),obsGridY(l) /),(/ x,y /))
     relevantObs(l) = weight(l) <= hMaxCorrLengthToObs(index)
   end do

!   write(6,*) 'relevantObs ',count(relevantObs),minval(weight),maxval(weight)

   noRelevantObs = .not.any(relevantObs)

   if (.not.noRelevantObs) weight = exp(- (weight/hCorrLengthToObs(index))**2)

   contains 

    ! distance between point p0 and p1
    ! p0 and p1 are (/ longitude,latitude,... /)

    real function distance(p0,p1)
     implicit none
     real, intent(in) :: p0(:), p1(:)
     real :: d2r = pi/180.
     real :: a, b, C

     if (metrictype == CartesianMetric) then
       distance = sqrt(sum((p1 - p0)**2))

     elseif (metrictype == SphericalMetric) then
       ! assume sperical metric

!#define DISTANCE_SIMPLE
#ifdef DISTANCE_SIMPLE
     real :: coeff
     coeff = pi*EarthRadius/(180.)     
     distance = sqrt((coeff * cos((p0(2)+p1(2))* (pi/360.))*(p1(1)-p0(1)))**2 &
          +(coeff * (p1(2)-p0(2)))**2)
     
#else

     a = p0(2) * d2r
     b = p1(2) * d2r
     C = (p1(1) - p0(1)) * d2r
     
     ! distance in radian
     distance = acos(sin(b) * sin(a) + cos(b) * cos(a) * cos(C))

     ! distance in km
     distance = EarthRadius * distance
#endif
     else
       write(stderr,*) 'Unsupported metric: ',metrictype
       write(stderr,*) 'Supported metrics are CartesianMetric (0) and SphericalMetric (1)'
       ERROR_STOP
     end if

    end function distance




   end subroutine 

 !_______________________________________________________
 !
 ! transform objects to packed storage and back
 ! 
 !_______________________________________________________
 !

 subroutine packVector(ML,statevector,v1,v2,v3,v4,v5,v6,v7,v8)
  implicit none
  type(MemLayout), intent(in) ::  ML
  real, intent(in) :: v1(*)
  real, intent(in), optional :: v2(*),v3(*),v4(*),v5(*),v6(*),v7(*),v8(*)
  real, intent(out) :: StateVector(ML%effsize)

  StateVector(ML%StartIndexSea(1):ML%EndIndexSea(1)) = &
       pack(v1(1:ML%varsize(1)),ML%mask(ML%StartIndex(1):ML%EndIndex(1)).eq.1)

  if (.not.present(v2)) return

  StateVector(ML%StartIndexSea(2):ML%EndIndexSea(2)) = &
       pack(v2(1:ML%varsize(2)),ML%mask(ML%StartIndex(2):ML%EndIndex(2)).eq.1)

  if (.not.present(v3)) return

  StateVector(ML%StartIndexSea(3):ML%EndIndexSea(3)) = &
       pack(v3(1:ML%varsize(3)),ML%mask(ML%StartIndex(3):ML%EndIndex(3)).eq.1)

  if (.not.present(v4)) return

  StateVector(ML%StartIndexSea(4):ML%EndIndexSea(4)) = &
       pack(v4(1:ML%varsize(4)),ML%mask(ML%StartIndex(4):ML%EndIndex(4)).eq.1)

  if (.not.present(v5)) return

  StateVector(ML%StartIndexSea(5):ML%EndIndexSea(5)) = &
       pack(v5(1:ML%varsize(5)),ML%mask(ML%StartIndex(5):ML%EndIndex(5)).eq.1)

  if (.not.present(v6)) return

  StateVector(ML%StartIndexSea(6):ML%EndIndexSea(6)) = &
       pack(v6(1:ML%varsize(6)),ML%mask(ML%StartIndex(6):ML%EndIndex(6)).eq.1)

  if (.not.present(v7)) return

  StateVector(ML%StartIndexSea(7):ML%EndIndexSea(7)) = &
       pack(v7(1:ML%varsize(7)),ML%mask(ML%StartIndex(7):ML%EndIndex(7)).eq.1)

  if (.not.present(v8)) return

  StateVector(ML%StartIndexSea(8):ML%EndIndexSea(8)) = &
       pack(v8(1:ML%varsize(8)),ML%mask(ML%StartIndex(8):ML%EndIndex(8)).eq.1)

 end subroutine packVector

 !_______________________________________________________
 !

 subroutine unpackVector(ML,statevector,v1,v2,v3,v4,v5,v6,v7,v8)
  implicit none
  type(MemLayout), intent(in) ::  ML
  real, intent(in) :: StateVector(ML%effsize)

  real, intent(out) :: v1(*)
  real, intent(out), optional :: v2(*),v3(*),v4(*),v5(*),v6(*),v7(*),v8(*)

  integer :: i

  i=1
  v1(1:ML%varsize(i)) = &
       unpack(StateVector(ML%StartIndexSea(i):ML%EndIndexSea(i)), &
       ML%mask(ML%StartIndex(i):ML%EndIndex(i)).eq.1,0.)

  if (.not.present(v2)) return
  i=2
  v2(1:ML%varsize(i)) = &
       unpack(StateVector(ML%StartIndexSea(i):ML%EndIndexSea(i)), &
       ML%mask(ML%StartIndex(i):ML%EndIndex(i)).eq.1,0.)

  if (.not.present(v3)) return
  i=3
  v3(1:ML%varsize(i)) = &
       unpack(StateVector(ML%StartIndexSea(i):ML%EndIndexSea(i)), &
       ML%mask(ML%StartIndex(i):ML%EndIndex(i)).eq.1,0.)

  if (.not.present(v4)) return
  i=4
  v4(1:ML%varsize(i)) = &
       unpack(StateVector(ML%StartIndexSea(i):ML%EndIndexSea(i)), &
       ML%mask(ML%StartIndex(i):ML%EndIndex(i)).eq.1,0.)

  if (.not.present(v5)) return
  i=5
  v5(1:ML%varsize(i)) = &
       unpack(StateVector(ML%StartIndexSea(i):ML%EndIndexSea(i)), &
       ML%mask(ML%StartIndex(i):ML%EndIndex(i)).eq.1,0.)

  if (.not.present(v6)) return
  i=6
  v6(1:ML%varsize(i)) = &
       unpack(StateVector(ML%StartIndexSea(i):ML%EndIndexSea(i)), &
       ML%mask(ML%StartIndex(i):ML%EndIndex(i)).eq.1,0.)

  if (.not.present(v7)) return
  i=7
  v7(1:ML%varsize(i)) = &
       unpack(StateVector(ML%StartIndexSea(i):ML%EndIndexSea(i)), &
       ML%mask(ML%StartIndex(i):ML%EndIndex(i)).eq.1,0.)

  if (.not.present(v8)) return
  i=8
  v8(1:ML%varsize(i)) = &
       unpack(StateVector(ML%StartIndexSea(i):ML%EndIndexSea(i)), &
       ML%mask(ML%StartIndex(i):ML%EndIndex(i)).eq.1,0.)


 end subroutine unpackVector

 !_______________________________________________________
 !
 ! extract variable number "v" from state vector and put 
 ! it as 3D variable
 ! 
 !_______________________________________________________
 !

 function unpack3DVariable(ML,statevector,v) result(var)
  implicit none
  type(MemLayout), intent(in) ::  ML
  real, intent(in) :: StateVector(ML%effsize)
  integer :: v

  real :: var(ML%varshape(1,v),ML%varshape(2,v),ML%varshape(3,v))

  var = reshape( &
         unpack(StateVector(ML%StartIndexSea(v):ML%EndIndexSea(v)), &
           ML%mask(ML%StartIndex(v):ML%EndIndex(v)).eq.1,0.), &
        (/ ML%varshape(1,v),ML%varshape(2,v),ML%varshape(3,v) /) )

 end function

 !_______________________________________________________
 !
 ! put 3D variable number "v" into the state vector
 ! 
 !_______________________________________________________
 !

  subroutine pack3DVariable(ML,var,v,statevector)
  implicit none
  type(MemLayout), intent(in) ::  ML
  real, intent(in) :: var(:,:,:)
  integer, intent(in) :: v
  real, intent(inout) :: StateVector(ML%effsize)

  StateVector(ML%StartIndexSea(v):ML%EndIndexSea(v)) = &
       pack(reshape(var,(/ ML%varsize(v) /)),ML%mask(ML%StartIndex(v):ML%EndIndex(v)).eq.1)

 end subroutine




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

    linindex1 = sub2ind(ML1,Hindex(1,i),Hindex(2,i),Hindex(3,i),Hindex(4,i),val1)

    ! space 2: origin
    ! transform [Hindex(5,i) Hindex(6,i) Hindex(7,i) Hindex(8,i)] into the
    ! linear index linindex2 and trapp error in variable val2

    linindex2 = sub2ind(ML2,Hindex(5,i),Hindex(6,i),Hindex(7,i),Hindex(8,i),val2)

    ! return the valid flags is desiered

    if (present(valid1).and.val1) valid1(linindex1) = val2
    if (present(valid2).and.val2) valid2(linindex2) = val1
    !       valid1(i) = val1
    !       valid2(i) = val2


    if (val1.and.val2) then
#        ifdef DEBUG
      if (nz.ge.size(H%s)) then
        write(stderr,*) 'packSparseMatrix: ERROR: ', &
             'buffer variable too small!!! '
        call flush(stderr,istat)
      end if
#        endif

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
  integer, intent(out)  :: Hindex(8,H%nz)
  real, intent(out) :: Hcoeff(H%nz)
  type(MemLayout), intent(in) :: ML1,ML2

  integer nz 

  do nz=1,H%nz
    call ind2sub(ML1,H%i(nz),Hindex(1,nz),Hindex(2,nz),Hindex(3,nz),Hindex(4,nz))
    call ind2sub(ML2,H%j(nz),Hindex(5,nz),Hindex(6,nz),Hindex(7,nz),Hindex(8,nz))
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

!     call usave('/u/abarth/Assim/Data2/toto.TEM',unpack3DVariable(ModML,statevector,vT),0.)
!     call usave('/u/abarth/Assim/Data2/toto.SAL',unpack3DVariable(ModML,statevector,vS),0.)

! temporaly disabled

!     call TStransform(ModelGrid3D(v)%mask,ModelGrid3D(v)%z, &
!       unpack3DVariable(ModML,statevector,vT),  &
!       unpack3DVariable(ModML,statevector,vS),  &
!       X,Y)

!     call usave('/u/abarth/Assim/Data2/toto.TEM2',X,0.)
!     call usave('/u/abarth/Assim/Data2/toto.SAL2',Y,0.)

      call pack3DVariable(ModML,X,v,statevector)
      call pack3DVariable(ModML,Y,v+1,statevector)

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

!     call usave('/u/abarth/Assim/Data2/lulu.TEM2',unpack3DVariable(ModML,statevector,v),0.)
!     call usave('/u/abarth/Assim/Data2/lulu.SAL2',unpack3DVariable(ModML,statevector,v+1),0.)

! disabled

!     call invTStransform(ModelGrid3D(v)%mask,ModelGrid3D(v)%z, &
!       unpack3DVariable(ModML,statevector,v),  &
!       unpack3DVariable(ModML,statevector,v+1),  &
!       T,S)

      call pack3DVariable(ModML,T,vT,statevector)
      call pack3DVariable(ModML,S,vS,statevector)

!     call usave('/u/abarth/Assim/Data2/lulu.TEM',T,0.)
!     call usave('/u/abarth/Assim/Data2/lulu.SAL',S,0.)

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

!  call saveVectorSpace('Ef1.value',ModML,E)
!  call usave('/u/abarth/Assim/Data2/Ef.u',E,0.)

  do i=1,N
    call anamorph(E(:,i))
  end do

!  call saveVectorSpace('Ef2.value',ModML,E)
!  call usave('/u/abarth/Assim/Data2/Ef2.u',E,0.)

  !call ensanalysis(Ef,HEf,yo,invsqrtR,Ea, amplitudes)
  call ensanalysis(E,HEf,yo,invsqrtR,E,amplitudes)

!  call saveVectorSpace('Ea2.value',ModML,E)
!  call usave('/u/abarth/Assim/Data2/Ea2.u',E,0.)

  do i=1,N
    call invanamorph(E(:,i))
  end do

!  call saveVectorSpace('Ea1.value',ModML,E)
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

!  call saveVectorSpace('Ef2.value',ModML,E)
!  call usave('/u/abarth/Assim/Data2/Ef2.u',E,0.)

  !call ensanalysis(Ef,HEf,yo,invsqrtR,Ea, amplitudes)
    call ensAnalysis(Ea,HE,yo,invsqrtR,Ea,amplitudes)

!  call saveVectorSpace('Ea2.value',ModML,E)
!  call usave('/u/abarth/Assim/Data2/Ea2.u',E,0.)

    do i=1,N
      call invanamorph(Ea(:,i))
    end do

!  call saveVectorSpace('Ea1.value',ModML,E)
!  call usave('/u/abarth/Assim/Data2/E.u',E,0.)

   end if

 end subroutine

end module assimilation





