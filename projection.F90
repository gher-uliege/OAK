! include the fortran preprocessor definitions
#include "ppdef.h"


program opermul
  use matoper
  use rrsqrt
  use ufileformat
  use date
  use initfile
  use assimilation

  implicit none

  real, pointer      :: v1(:),v2(:),w(:)
  type(MemLayout) :: ML1
  character(256) :: str
  integer iargc

  if (iargc().ne.1) then
     write(stderr,*) 'Usage: projection <init file>'
     stop
  end if

  call getarg(1,str); initfname = str

  call MemoryLayout('',ML1)

  allocate(v1(ML1%effsize),v2(ML1%effsize),w(ML1%effsize))

  call loadVector('Vec1.value',ML1,v1)
  call loadVector('Vec1.value',ML1,v2)
  call loadVector('Norm.value',ML1,w)

  write(6,*) sum(v1*w*v2)

  deallocate(v1,v2,w)
end program
