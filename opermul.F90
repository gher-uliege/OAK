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

  real, pointer      :: x(:)
  type(SparseMatrix) :: H1,H2,prod
  type(MemLayout) :: ML1,ML2,ML3
  character(256) :: str
  integer iargc
  logical, pointer :: valid1(:),valid2(:)

  if (iargc().ne.1) then
     write(stderr,*) 'Usage: opermul <init file>'
     write(stderr,*) 
     write(stderr,*) ' multiplication of to operators'
     write(stderr,*) 
     write(stderr,*) ' Product = Term2 * Term1'
     write(stderr,*) 
     write(stderr,*) ' Term1 : Space1 -> Space2'
     write(stderr,*) ' Term2 : Space2 -> Space3'
     write(stderr,*) 
     write(stderr,*) ' Product : Space1 -> Space3'
     write(stderr,*) 
     write(stderr,*) ' INPUT:'
     write(stderr,*) '   Space1,Space2,Space3,Term1,Term2'
     write(stderr,*) ' OUTPUT:'
     write(stderr,*) '   Product'

     stop
  end if

  call getarg(1,str); initfname = str

  call MemoryLayout('Space1.',ML1)
  call MemoryLayout('Space2.',ML2)
  call MemoryLayout('Space3.',ML3)

  allocate(valid1(ML1%effsize),valid2(ML3%effsize))

  call loadSparseMatrix('Term1',ML2,ML1,H1)
  call loadSparseMatrix('Term2',ML3,ML2,H2)

  prod = H2.x.H1
  deallocate(H1%i,H1%j,H1%s,H2%i,H2%j,H2%s)
  call saveSparseMatrix('Product',ML3,ML1,prod)
!  call saveSparseMatrix('/u/abarth/soft/Ligur2/Obs/2000/','Product',ML3,ML1,prod)

!!$  call loadSparseMatrix('Term1',ML2,ML1,H1,valid1=valid2)
!!$  call loadSparseMatrix('Term2',ML3,ML2,H2,valid2=valid1)
!!$
!!$  call saveSparseMatrix('Product',ML3,ML1,H1.x.H2,valid1,valid2)
!  call saveSparseMatrix('Product',ML1,ML3,H1.x.H2)
end program
