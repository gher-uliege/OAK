! include the fortran preprocessor definitions
#include "ppdef.h"


program anamorph
  use matoper
  use rrsqrt
  use ufileformat
  use date
  use initfile
  use assimilation

  implicit none

  real, pointer      :: v1(:),v2(:),w(:)
  character(256) :: str
  integer iargc,TransformDir

  if (iargc().ne.1) then
     write(stderr,*) 'Usage: anamorph <init file>'
     stop
  end if

  call getarg(1,str); initfname = str

  call init(initfname)

!  call MemoryLayout('Model.',ModML)
  call getInitValue(initfname,'TransformDir',TransformDir)

  allocate(v1(ModML%effsize))

  call loadVector('Vec1.value',ModML,v1)

  if (TransformDir.eq.0) then
    call anamorphosisTransform(v1)
  else
    call invanamorphosisTransform(v1)
  end if
  call saveVector('Vec2.value',ModML,v1)

  deallocate(v1)
end program
