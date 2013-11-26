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

module setgrid

   type cell
     ! arrary of indices pointing to the first dimension of modGrid
     integer, allocatable :: ind(:)
     ! number of elements in ind with a valid value
     integer :: n
   end type cell
   
   type cellgrid
     type(cell), allocatable :: grid(:,:)
     real, allocatable :: xmin(:), dx(:)
     integer, allocatable :: Ni(:)
     integer :: n
   end type cellgrid

contains

  function setupgrid(modGrid) result(cg)
   real, intent(in) :: modGrid(:,:)
   real, allocatable :: xmax(:)
   integer, allocatable :: gridind(:)
   integer :: n=, i, j, l
   type(cellgrid) :: cg
   
   integer :: startsize = 100

   n = size(modGrid,2)
   cg%n = n

   allocate(cg%xmin(n),xmax(n),cg%dx(n),cg%Ni(n),gridind(n))

   cg%dx = 0.2

   do i = 1,n
     cg%xmin(i) = minval(modGrid(:,i))
     xmax(i) = maxval(modGrid(:,i))
   end do

   ! number of intervales in each dimensions
   cg%Ni = ceiling((xmax - cg%xmin) / cg%dx)

   ! the j(1),j(2),...j(n) box is surrounded by the grid points
   ! lower: xmin(i) + dx(i) * (j(i)-1) for i=1,...,n
   ! upper: xmin(i) + dx(i) * j(i) for i=1,...,n
   

   ! allocate grid of cells and say that nothing is in there
   allocate(cg%grid(Ni(1),Ni(2)))
   do j = 1,Ni(2)
     do i = 1,Ni(1)
       allocate(cg%grid(i,j)%ind(startsize))
       cg%grid(i,j)%n = 0
     end do
   end do

   ! loop throught every coordinate and put index into a cell
   do l = 1,size(modGrid,1)
     ! compute grid index
     gridind = floor((modGrid(l,:) - xmin)/ dx) + 1
     call add(cg%grid(gridind(1),gridind(2)),l)
     if (mod(l,100) == 0) write(6,*) 'l ',l,gridind         
     !grid(gridind(1),gridind(2))
   end do

   contains 

   subroutine add(c,l)
    type(cell), intent(inout) :: c
    integer, intent(in) :: l
    integer, allocatable :: tmp(:)
    integer :: sz

    c%n = c%n + 1

    if (c%n > size(c%ind)) then
      ! need to grow c%ind
      sz = size(c%ind)
      allocate(tmp(sz))

      ! make backup
      tmp = c%ind

      deallocate(c%ind)
      allocate(c%ind(2*sz))
      c%ind(1:sz) = tmp

      deallocate(tmp)
    end if

    c%ind(c%n) = l

   end subroutine add


  end function setupgrid

   subroutine near(x,ind,n)
    real, intent(in) :: x(:)
    gridind = floor((x - xmin)/ dx) + 1
    n = grid(gridind(1),gridind(2))%n
    ind(:n) = grid(gridind(1),gridind(2))%ind
   end subroutine near

end module setgrid

program assimtest
 use matoper
 use rrsqrt
 use ufileformat
 use initfile
 use assimilation
 use parall
 use setgrid

 implicit none

 real, allocatable  :: xf(:), xa(:), &
      Sf(:,:), Sa(:,:)

 integer            :: iargc, ntime, enstype
 character(len=124) :: str
 character(len=124) :: ntimeindex


 real, pointer :: modGrid(:,:)

 if (iargc().ne.2) then
   write(stderr,*) 'Usage: assim <init file> <time index> '
   ERROR_STOP
 end if

#ifdef ASSIM_PARALLEL
   call parallInit()
#endif

 call getarg(1,str); call init(str)

   allocate(modGrid(ModML%effsize,2))
   call loadVector('Model.gridX',ModML,modGrid(:,1))
   call loadVector('Model.gridY',ModML,modGrid(:,2))

 call setupgrid(modGrid)

 stop


 call getarg(2,str); read(str,*) ntime  

 allocate(xf(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel), &
      xa(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel), &
      Sa(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel,ErrorSpaceDim), &
      Sf(ModMLParallel%startIndexParallel:ModMLParallel%endIndexParallel,ErrorSpaceDim))

 !write(ntimeindex,'(I3.3,A)') ntime,'.'
 call fmtIndex('',ntime,'.',ntimeindex)

 call getInitValue(initfname,'ErrorSpace.type',enstype,default=1)
 if (enstype == 1) then
   ! Sf are error modes and xf is the forecast
   call loadVectorSpace('ErrorSpace.init',ModMLParallel,Sf)
   call loadVector('Forecast'//trim(ntimeindex)//'value',ModMLParallel,xf)
!$omp parallel
   call assim(ntime,Sf,Sa,xf,xa)
!$omp end parallel
 else
   ! Sf are ensemble members
   call loadEnsemble('ErrorSpace.init',ModMLParallel,Sf)

!$omp parallel
   call assim(ntime,Sf,Sa)   
!$omp end parallel
 end if

 call done()
#ifdef ASSIM_PARALLEL
 call parallDone()
#endif


 deallocate(xf,xa,Sa,Sf)

 contains


end program assimtest

