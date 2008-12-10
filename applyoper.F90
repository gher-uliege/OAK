! include the fortran preprocessor definitions
#include "ppdef.h"

program applyobsoper
 use matoper
 use rrsqrt
 use date
 use ufileformat
 use assimilation
 use initfile

  implicit none

  real, pointer      :: x(:)
  type(SparseMatrix) :: H
  type(MemLayout) :: ML
  character(256) :: str
  integer iargc

  if (iargc().ne.1) then
     write(stderr,*) 'Usage: applyoper <init file>'
     stop
  end if


  call getarg(1,str); initfname = str

  call MemoryLayout('Model.',ML)
  allocate(x(ML%effsize))

!  call loadVector('In.value',ML,x)
  call loadSparseMatrix('Operator',ML,ML,H)
!  call saveSparseMatrix('filter1.u',ML,ML,H)
   call loadVector('In.value',ML,x)

  call saveVector('Out.value',ML,H.x.x)
end program

