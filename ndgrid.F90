!
!  OAK, Ocean Assimilation Kit
!  Copyright(c) 2002-2016 Alexander Barth and Luc Vandenblucke
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
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
!

!
! interpolation of arbitrary N-dimentional grids
!
! coordinate can be givem explicetly (as a variable) or implicitly (as a function)
! 
!
!


!
! include the fortran preprocessor definitions
!
#include "ppdef.h"

#define DATABOX_SEARCH
!#define IMPLICIT_COORDINATE

module ndgrid


#ifdef DATABOX_SEARCH
 type databoxND
   ! xmin(i) and xmax(i) are the minimum and maximum of the ith coordinate
   real, pointer :: xmin(:), xmax(:)

   ! type of the databox, possible values are db_cell, db_splitted and db_notsplitted
   integer :: type

   ! zero-based indexes of the databox
   integer, pointer :: imin(:),imax(:)

   ! sub-databoxes
   type(databoxND), pointer :: db(:)   
 end type databoxND
#endif

! definition of a n-dimensional aribtrary grid

 type, abstract :: basegrid
   ! n number of dimensions (e.g. 2 for lon/lat)

   integer :: n

   ! gshape(n): shape of the grid, gshape(i) is the maximum value of the ith one-based index
   integer, pointer :: gshape(:)

   ! .false. if point on grid is valid, .true. if is masked
   logical, pointer :: masked(:)

   ! offset for mask
   integer, pointer :: ioffset_mask(:)

   ! how the cube should be splitted in hyper-tetrahedrons, for 2D just 
   ! triangles:

   !  a-----------b
   !  |\         /|
   !  | \       / |
   !  |  \     /  |
   !  |   \   /   |
   !  |    \ /    |
   !  |     m
   !  |    / \    |
   !  |   /   \   |
   !  |  /     \  |
   !  | /       \ |
   !  |/         \|
   !  d-----------c
   !
   ! all triangles are: 
   ! a b m
   ! b c m
   ! c d m
   ! d a m

   real, pointer :: tetrahedron(:,:,:)

   contains
    procedure(getCoord_interface), deferred, pass(g) :: getCoord    
    procedure(locate_inteface), deferred, pass(g) :: locate
 end type basegrid

 abstract interface
   function getCoord_interface(g,ind,out) result(x)
    import :: basegrid
    implicit none
    class(basegrid), intent(in) :: g
    integer, intent(in) :: ind(:)
    logical, optional, intent(out) :: out
    real :: x(size(ind))
   end function getCoord_interface   

   subroutine locate_inteface(g,x,ind,out) 
    import :: basegrid
    implicit none
    class(basegrid), intent(inout) :: g
    real, intent(in) :: x(:)
    integer, intent(out) :: ind(size(x))
    logical, intent(out) :: out
   end subroutine locate_inteface
 end interface


 ! grid with a constant increment
 type, extends(basegrid) :: regulargrid
   ! coordinate of point (i,j)  is 
   ! x0(1) + i*dx(1), x0(2) + j*dx(2)

   real, allocatable :: x0(:), dx(:) 

  contains
   procedure, pass(g) :: getCoord => getCoord_regulargrid
   procedure, pass(g) :: locate => locate_regulargrid
 end type regulargrid


 type, extends(basegrid) :: grid

   ! ioffset(i,j) offset of the ith coordinate with repect to the jth index
   integer, pointer :: ioffset(:,:)

   ! dependence(n,n): dependence(i,j) is 1 if the jth coordinate dependend on the ith index, ortherwise 0
   integer, pointer :: dependence(:,:)

   ! coordinates
   !real, pointer :: coord(:,:) 
   real, pointer :: data(:)

   ! startindex(n)
   ! the corridnate value relevant for intex i start at data(startindex(i))
   integer, pointer :: startindex(:)
   integer, pointer :: endindex(:)
   
#  ifdef DATABOX_SEARCH
   type(databoxND) :: db
#  endif

  contains
   procedure, pass(g) :: getCoord => getCoord_grid
   procedure, pass(g) :: locate => locate_grid
 end type grid




! possible values of type
  integer, parameter :: db_cell = 0
  integer, parameter :: db_splitted = 1
  integer, parameter :: db_notsplitted = 2

! grid constructor

interface initgrid
 module procedure initgrid_nd, initgrid_1d, initgrid_2d, initgrid_3d, initgrid_4d
end interface

! top-level grid interpolation subroutine
!
! interpgrid_E : with explicit coordinates and databox search
! interpgrid_I : with implicit coordinates and databox search
! interpgrid_EL: with explicit coordinates and user provided search function
! interpgrid_IL: with implicit coordinates and user provided search function

interface interpgrid 
 module procedure interpgrid_E, interpgrid_I,interpgrid_EL,interpgrid_IL
end interface

interface interpgrid_coeff
 module procedure interpgrid_coeff_E, interpgrid_coeff_I,interpgrid_coeff_EL,interpgrid_coeff_IL
end interface

interface locate_databox
 module procedure locate_databox_E, locate_databox_I
end interface

interface init_databox
 module procedure init_databox_E, init_databox_I
end interface

interface InCube
 module procedure InCube_E, InCube_I
end interface


interface setCoord
 module procedure setCoord_fromFile, setCoord_fromVariable
end interface



 ! put discrete points in cell "boxes"

 ! cell containing a series of points maked by their index
 type cell
   ! arrary of indices pointing to the first dimension of xpos
   integer, allocatable :: ind(:)
   ! number of elements in ind with a valid value
   integer :: n
 end type cell

 ! a regular grid of cells
 type cellgrid
   type(cell), allocatable :: gridn(:)
   real, allocatable :: xmin(:), dx(:)
   integer, allocatable :: Ni(:), offset(:)
   integer :: n
 end type cellgrid

 ! interface for the distance between two points of coordinate x and coordinates y
 interface
   real function distind(x,y)
    real, intent(in) :: x(:),y(:)
   end function distind
 end interface



contains

 subroutine init_regulargrid(g,gshape,x0,dx,masked)
  implicit none
  type(regulargrid) :: g
  integer, intent(in) :: gshape(:)
  real, intent(in) :: x0(:), dx(:) ! offset point and increment
  logical, optional :: masked(:)

  integer :: n
  real, allocatable :: pm(:,:), x(:,:)
  logical, allocatable :: masked_(:)

  n = size(gshape)

  call init_basegrid(g,gshape,masked)
  allocate(g%x0(n),g%dx(n))

  g%x0 = x0
  g%dx = dx

 end subroutine init_regulargrid


 ! get coordinates x of point at indices ind
 function getCoord_regulargrid(g,ind,out) result(x)
  implicit none
  class(regulargrid), intent(in) :: g
  integer, intent(in) :: ind(:)
  logical, optional, intent(out) :: out
  real :: x(size(ind))

  if (present(out)) then
    out = any(ind < 1) .or. any(ind > g%gshape)
  end if

  x = g%x0 + ind * g%dx

  ! TODO define out based on mask
 end function getCoord_regulargrid

 ! get all points (or more!) near the location x
 
 subroutine near_regulargrid(g,x,distfun,maxdist,ind,dist,n)
  use matoper
  implicit none
  type(regulargrid), intent(in) :: g
  real, intent(in) :: x(:) !,xpos(:,:)
  procedure(distind) :: distfun
  ! maximum distance to x
  real, intent(in) :: maxdist
  ! actual distance to x
  real, intent(out) :: dist(:)
  ! indices and number of points
  integer, intent(out) :: ind(:), n

  real :: x0(g%n), dx(g%n), xi(g%n)
  ! distance 
  real :: d
  integer :: i,j,k
  integer :: index0(g%n), index1(g%n), index(g%n), subshape(g%n)
  logical :: out

  x0 = g%x0
  dx = g%dx

  index0 = max(floor((x-x0 - maxdist)/dx) +1,1)
  index1 = min(ceiling((x-x0 + maxdist)/dx) +1,g%gshape)
  subshape = index1-index0+1
  n = 0

  !write(6,*) 'index0',index0
  !write(6,*) 'index1',index1

  do k = 1,product(subshape)
    index = ind2sub(subshape,k) + index0-1
    i = sub2ind(g%gshape,index)

    xi = getCoord_regulargrid(g,index,out)
    !write(6,*) 'index ',index
    !write(6,*) 'xi ',xi
    d = cdist(xi,x)
    !write(6,*) 'd ',d

    if (d <= maxdist) then
      n = n+1
      dist(n) = d     
      ind(n) = i
    end if
  end do
  
 end subroutine near_regulargrid



 subroutine locate_regulargrid(g,x,ind,out) 
  implicit none
  class(regulargrid), intent(inout) :: g
  real, intent(in) :: x(:)
  integer, intent(out) :: ind(size(x))
  logical, intent(out) :: out
  
  ind = floor((x - g%x0) / g%dx) 
  ! 0-based
  ind = ind-1
  out = any(ind < 1 .or. ind > g%gshape)
 end subroutine locate_regulargrid

 !
 ! --------------------------------------------------------------------
 !

 ! x(components , corner point indexes)

 ! c(corner points of cube, corner points of cube)

 ! tetrahedron(corner points of cube, corner points of tetrahedron, index of tetrahedron)

 ! tetrahedron(1:2^n, corner point, thetraeders)  = coefficient

 recursive subroutine split(n,subn,fixeddim,selection,tc)
  implicit none
  integer, intent(in) :: n,subn
  logical, intent(in) :: fixeddim(n),selection(2**n)
!  real, intent(in) :: c(2**n,2**n)

  real, intent(out) :: tc(:,:,:)

  real :: cm(2**n)
  integer :: twon,nbth,subnbth,i,j,k,m
  logical :: fd(n),subselection(2**n),l

  twon = 2**n

  if (subn.eq.1) then
    ! we have a straight line
!!$    do i=1,twon
!!$      tc(i,1:2,1) = pack(c(i,:),selection)
!!$    end do

    tc = 0
    j=1
   
    do i=1,twon
      if (selection(i)) then
        tc(i,j,1) = 1 
        j=j+1
      end if
    end do

    return
  end if

  ! middle point

!  cm = sum(c,2,spread(selection,1,twon))/count(spread(selection,1,twon),2)

  cm = 0
  do i=1,twon
    if (selection(i))  cm(i) = 1./count(selection)
  end do

  nbth = factorial(subn)*2**(subn-1)
  subnbth = factorial(subn-1)*2**(subn-2)

  ! add middle point

  tc(:,subn+1,1:nbth) = spread(cm,2,nbth)
  !  do i=1,nbth
  !    tc(:,subn+1,i) = cm(:)
  !  end do


  m=1
  do i=1,n
    if (.not.fixeddim(i)) then
      fd = fixeddim
      fd(i) = .true.

      do j=0,1
!        subselection = selection .and. (corner_indexes(i,:).eq.j)

        do k=1,twon
          if (btest(k-1,i-1)) then
            subselection(k) = selection(k) .and. (1.eq.j)
          else
            subselection(k) = selection(k) .and. (0.eq.j)
          end if
        end do

        ! split the faces of the cube
        call split(n,subn-1,fd,subselection,tc(:,1:subn,m:m+subnbth-1) )

        m=m+subnbth
      end do
    end if
  end do

 end subroutine split



 !
 ! --------------------------------------------------------------------
 !

 function factorial(n) result(f)
  implicit none
  integer, intent(in) :: n
  integer :: f

  integer ::i
  f = 1

  do i=2,n
    f = f*i
  end do

 end function factorial

 !
 ! --------------------------------------------------------------------
 !

 ! Y(1) first data point located at X(:,1)
 ! X( components, corner points )

 subroutine interp_tetrahedron(n,X,xi,out,coeff) 
  use matoper
  implicit none
  integer, intent(in) :: n
  real, intent(in) :: X(n,n+1),xi(n)
  real, intent(out) :: coeff(n+1)
  logical, intent(out) :: out

  ! better force double precision here
  ! precision of the inverse
  real(8), parameter :: tol = 1e-8
  ! local variables

  real(8) :: M(n+1,n+1), c(n+1), M2(n+1,n+1), d(n+1), d2(n+1)
  real(8) :: determ

  integer :: i, j, nnz, nuncon, offset(n+1), ind(n+1), tmp
  integer :: nz(n+1)
  real(8) :: c2(n), X8(n,n+1), U(n+1,n+1), V(n+1,n+1), S(n+1)

  real(8) :: xc(n), err, best, detM2, testc(n+1)
  logical :: cinit

  ! look for c which satisfies:
  !
  ! 1 = c(1) + c(2) + ... + c(n+1)
  ! xi(:) = c(1) * X(:,1) + c(2) * X(:,2) + ... + c(n+1) * X(:,n+1)

  ! In matrix form
  ! 
  ! [ 1  ]  = M * c
  ! [ xi ]
  ! 
  ! where M
  !
  ! M = [  1      1       ...  1      ]
  !     [  X(:,1) X(:,2)  ... X(:,n+1)]
  !

  M = 1
  d = 1

  ! average all X and make vectors relative to this location
  ! for better numerical precision

  xc = sum(X,2)/(n+1)
  M(2:n+1,:) = X - spread(xc,2,n+1)
  d(2:n+1) = xi - xc

  c = matmul(inv(M,determ),d)

  if (abs(determ) > tol) then
    out = .not.all(0d0-tol <= c .and.c <= 1d0+tol)
    coeff = real(c)
    return
  end if

  !write(6,*) 'degenerated case (0)',determ
  !write(6,*) 'degenerated case (0)',X
  !write(6,*) 'degenerated case (0)',n
  !write(6,*) 'degenerated case (0)',xi

  ! degenerated case
  S = svd(M,U,V)

  ! find all indices of elements in S larger than tol
  ! number of non-redundant constraints
  nnz = count(S > tol)

  nz(1:nnz) = pack((/ (i,i=1,n+1) /),S > tol)


  ! number of coefficients uncontrained
  ! nuncon + nnz = n+1
  nuncon = n+1-nnz


  coeff = 0
  ! is coef initialized with a meanful value?
  cinit = .false.;
  best = 0;
  out = .true.;

  ! offset used to convert linear index i to subscribt
  ! similar to ind2sub in matlab/octave

  offset(1) = 1
  do j = 2,nuncon
    offset(j) = offset(j-1) * (n+1)
  end do

  ! force one coeff to zero after the other

  do i = 1, (n+1)**nuncon

    ! indices for vertices to test if the corresponding c can be zero
    ![ind(1),ind(2)] = ind2sub((n+1)*ones(nuncon,1),i);

    ! tmp and ind are here 0-based
    tmp = i - 1;

    do j = nuncon, 1, -1
      ! integer devision
      ind(j) = tmp / offset(j)
      tmp = tmp - ind(j) * offset(j);
    end do
    ! make ind 1-based
    ind = ind+1;

    M2 = 0

    do j = 1, nuncon
      M2(j,ind(j)) = 1;
    end do

    if (any(sum(M2(1:nuncon,:),2) == 0)) then
      ! all elements of ind must be different
      ! ignore the present case
      cycle
    end if

    ! with S,V
    ! ! remove redundant contraints
    !  S * V' * coeff = U' * d
    !  A*c = 0

    M2(nuncon+1:,:) = matmul(diag(S(nz(1:nnz))), transpose(V(:,nz(1:nnz))))

    d2 = 0
    d2(nuncon+1:) = matmul(transpose(U(:,nz(1:nnz))), d)

    testc = matmul(inv(M2,detM2), d2);

    if (abs(detM2) < tol) then
      ! still degenerated, look for other possibilities
      cycle
    end if

    if (.not.all(0-tol <= testc .and. testc <= 1+tol)) then
      ! no, this is a extrapolation!
      ! go to next iteration
      cycle
    end if

    err = maxval(abs(matmul(M,testc) - d));

    ! but, wait check if it is still a solution to our problem
    ! this is necessary if xi is outside a degenerated X, for example
    ! xi = [0 1]';
    ! X = [0    0.5   1; ...
    !      0      0   0];

    if (err < tol) then
      ! abs(det(M2)) is an indication of the surface of the
      ! "triangle", the smaller the better
      if (abs(detM2) < best  .or. .not.cinit) then
        ! we have a even better coeff
        coeff = testc;
        cinit = .true.;
        best = abs(detM2);
        out = .false.;
      end if
    end if
  end  do

 end subroutine interp_tetrahedron
 !
 ! --------------------------------------------------------------------
 !



 subroutine interp_cube(n,tetrahedron,x,xi,out,c)
  implicit none
  integer, intent(in) :: n
  real, intent(in) :: x(n,2**n),xi(n)
  real, intent(in) :: tetrahedron(:,:,:)
  logical, intent(out) :: out
  real, optional, intent(out) :: c(2**n)


  real :: coeff(n+1)
  integer :: l

  !   write(6,*) 'cube ',allocated(tetrahedron)
  !c = 0

  !  write(6,*) 'cube ',allocated(tetrahedron)

  ! search over all tetrahedrons

  do l=1,size(tetrahedron,3)

    call interp_tetrahedron(n,matmul(x,tetrahedron(:,:,l)),xi,out,coeff)

    if (.not.out)  then
      if (present(c)) c = matmul(tetrahedron(:,:,l),coeff)
      return
    end if
  end do

 end subroutine interp_cube


! four different different variants of the interpolation routine are
! possible weather the corrdinate are given explicitly (as a variabble)
! or implicitly (as a function) or if a point is localised with
! the default routine or localisation function provided by the used.
!
! The different variants are obtained by defining or not the preprocessor
! macros IMPLICIT_COORDINATE and USER_LOCATE


! explicit coordinate
! no user locate

#undef IMPLICIT_COORDINATE
#undef USER_LOCATE

#define interpgrid_coeff_VARIANT    interpgrid_coeff_E
#define locate_testall_VARIANT      locate_testall_E
#define InCube_VARIANT              InCube_E
#define interpgrid_VARIANT          interpgrid_E
#define init_databox_VARIANT        init_databox_E
#define init_databox_helper_VARIANT init_databox_helper_E
#define locate_databox_VARIANT      locate_databox_E
#define define_databox_VARIANT      define_databox_E
#define split_databox_VARIANT       split_databox_E
#define search_boundarybox_VARIANT      search_boundarybox_E

#include "ndgrid_inc.F90"

#undef interpgrid_coeff_VARIANT
#undef locate_testall_VARIANT 
#undef InCube_VARIANT 
#undef interpgrid_VARIANT 
#undef init_databox_VARIANT 
#undef init_databox_helper_VARIANT 
#undef locate_databox_VARIANT 
#undef define_databox_VARIANT 
#undef split_databox_VARIANT 
#undef search_boundarybox_VARIANT

! implicit coordinate
! no user locate

#define IMPLICIT_COORDINATE
#undef USER_LOCATE

#define interpgrid_coeff_VARIANT    interpgrid_coeff_I
#define locate_testall_VARIANT      locate_testall_I
#define InCube_VARIANT              InCube_I
#define interpgrid_VARIANT          interpgrid_I
#define init_databox_VARIANT        init_databox_I
#define init_databox_helper_VARIANT init_databox_helper_I
#define locate_databox_VARIANT      locate_databox_I
#define define_databox_VARIANT      define_databox_I
#define split_databox_VARIANT       split_databox_I
#define search_boundarybox_VARIANT      search_boundarybox_I

#include "ndgrid_inc.F90"

#undef interpgrid_coeff_VARIANT
#undef locate_testall_VARIANT 
#undef InCube_VARIANT 
#undef interpgrid_VARIANT 
#undef init_databox_VARIANT 
#undef init_databox_helper_VARIANT 
#undef locate_databox_VARIANT 
#undef define_databox_VARIANT 
#undef split_databox_VARIANT 
#undef search_boundarybox_VARIANT

! explicit coordinate
! user locate

#undef IMPLICIT_COORDINATE
#define USER_LOCATE

#define interpgrid_coeff_VARIANT    interpgrid_coeff_EL
#define locate_testall_VARIANT      locate_testall_EL
#define InCube_VARIANT              InCube_EL
#define interpgrid_VARIANT          interpgrid_EL
#define init_databox_VARIANT        init_databox_EL
#define init_databox_helper_VARIANT init_databox_helper_EL
#define locate_databox_VARIANT      locate_databox_EL
#define define_databox_VARIANT      define_databox_EL
#define split_databox_VARIANT       split_databox_EL
#define search_boundarybox_VARIANT      search_boundarybox_EL

#include "ndgrid_inc.F90"

#undef interpgrid_coeff_VARIANT
#undef locate_testall_VARIANT 
#undef InCube_VARIANT 
#undef interpgrid_VARIANT 
#undef init_databox_VARIANT 
#undef init_databox_helper_VARIANT 
#undef locate_databox_VARIANT 
#undef define_databox_VARIANT 
#undef split_databox_VARIANT 
#undef search_boundarybox_VARIANT

! implicit coordinate
! user locate

#define IMPLICIT_COORDINATE
#define USER_LOCATE

#define interpgrid_coeff_VARIANT    interpgrid_coeff_IL
#define locate_testall_VARIANT      locate_testall_IL
#define InCube_VARIANT              InCube_IL
#define interpgrid_VARIANT          interpgrid_IL
#define init_databox_VARIANT        init_databox_IL
#define init_databox_helper_VARIANT init_databox_helper_IL
#define locate_databox_VARIANT      locate_databox_IL
#define define_databox_VARIANT      define_databox_IL
#define split_databox_VARIANT       split_databox_IL
#define search_boundarybox_VARIANT      search_boundarybox_IL

#include "ndgrid_inc.F90"

#undef interpgrid_coeff_VARIANT
#undef locate_testall_VARIANT 
#undef InCube_VARIANT 
#undef interpgrid_VARIANT 
#undef init_databox_VARIANT 
#undef init_databox_helper_VARIANT 
#undef locate_databox_VARIANT 
#undef define_databox_VARIANT 
#undef split_databox_VARIANT 
#undef search_boundarybox_VARIANT

! coordinate as grid type
! no user locate

#undef IMPLICIT_COORDINATE
#define GRID_COORDINATE
#undef USER_LOCATE

#define interpgrid_coeff_VARIANT    interpgrid_coeff_G
#define locate_testall_VARIANT      locate_testall_G
#define InCube_VARIANT              InCube_G
#define interpgrid_VARIANT          interpgrid_G
#define init_databox_VARIANT        init_databox_G
#define init_databox_helper_VARIANT init_databox_helper_G
#define locate_databox_VARIANT      locate_databox_G
#define define_databox_VARIANT      define_databox_G
#define split_databox_VARIANT       split_databox_G
#define search_boundarybox_VARIANT      search_boundarybox_G

#include "ndgrid_inc.F90"



! routines independent of the variant

 function linindex2index(m,ioffset) result(ind)
  implicit none
  integer, intent(in) :: m,ioffset(:)
  integer :: ind(size(ioffset))

  integer :: d
  ind(1) = m
  do d=size(ioffset) ,2,-1
    ind(d) = ind(1)/ioffset(d)
    ind(1) = ind(1)-ind(d)*ioffset(d)
  end do

 end function linindex2index



 !_______________________________________________________
 !
 ! deallocate data box
 !_______________________________________________________
 !

 recursive subroutine done_databox(db)
  implicit none
  type(databoxND), intent(out) :: db

  integer :: n

  deallocate(db%imin,db%imax,db%xmin,db%xmax)

!  if (associated(db%db)) then
  if (db%type .eq. db_splitted) then
    ! deallocate all sub data boxes
    do n=1,size(db%db)
      call done_databox(db%db(n))
    end do
    deallocate(db%db)
  end if
 end subroutine done_databox


!_______________________________________________________________
!
! methods for type grid
!
!________________________________________________________________
!

 subroutine init_basegrid(g,gshape,masked)
  implicit none
  class(basegrid), intent(out) :: g
  integer, intent(in) :: gshape(:)
  logical, optional :: masked(:)

  integer :: n,i,nbth

  n = size(gshape)
  allocate(g%gshape(n),g%ioffset_mask(n))
  allocate(g%masked(product(gshape)))

  nbth =  factorial(n)*2**(n-1)
  allocate(g%tetrahedron(2**n,n+1,nbth))

  g%n = n
  g%gshape = gshape
  g%ioffset_mask(1) = 1

  do i=2,n
    g%ioffset_mask(i) = g%ioffset_mask(i-1)*(g%gshape(i-1))
  end do
  
  call split(n,n,(/ (.false.,i=1,n) /),(/ (.true.,i=1,2**n) /),g%tetrahedron)
  
  if (present(masked)) then
    g%masked = masked
  else
    g%masked = .false.
  end if

 end subroutine init_basegrid




!
! Initialization
!
! create a grid g of given dimension, shape and mask
! Input:
!   n: dimension of the grid
!   gshape: shape of the grid
!   masked: true if grid point is invalid (collapsed vector)
! Output: 
!   g: grid

 subroutine initgrid_nd(g,n,gshape,masked,dependence)
  implicit none
  type(grid), intent(out) :: g
  integer, intent(in) :: n,gshape(n)
  integer, optional, intent(in) :: dependence(n,n)
  logical, optional :: masked(:)

  integer :: i,j,nbth
  real, dimension(2**n,2**n) :: identity
!  write(6,*) 'n ',n
!  n = size(gshape)

  g%n = n
  allocate(g%gshape(n),g%ioffset(n,n),g%ioffset_mask(n),g%dependence(n,n),g%startindex(n),g%endindex(n))


  g%gshape = gshape
!  write(6,*) 'n ',   g%gshape
  if (present(dependence)) then
    g%dependence = dependence
  else
    g%dependence = 1
  end if

!  write(6,*) '  g%dependence ', g%dependence(:,1),g%gshape

  g%startindex(1) = 1
  g%endindex(1) = product(g%gshape,g%dependence(:,1).eq.1)
  g%ioffset_mask(1) = 1

  do i=2,n
    g%startindex(i) = g%endindex(i-1) +1
    g%endindex(i) = g%endindex(i-1) + product(g%gshape,g%dependence(:,i).eq.1)
    g%ioffset_mask(i) = g%ioffset_mask(i-1)*(g%gshape(i-1))
  end do

!  write(6,*) '  end ',  g%endindex
!  write(6,*) '  start ',  g%startindex

  do j=1,n
    do i=1,n
      if (g%dependence(i,j).eq.1) then
        g%ioffset(i,j) = product(gshape(1:i-1),g%dependence(1:i-1,j).eq.1)
      else
        g%ioffset(i,j) = 0
      end if
    end do
  end do

  allocate(g%data( g%endindex(n) ))

  nbth =  factorial(n)*2**(n-1) 
  allocate(g%tetrahedron(2**n,n+1,nbth))

 ! identity matrix

  identity = 0.
  do i=1,2**n
    identity(i,i) = 1.
  end do

  g%tetrahedron = 0
  call split(n,n,(/ (.false.,i=1,n) /),(/ (.true.,i=1,2**n) /),g%tetrahedron)

  ! initialise mask

  allocate(g%masked(product(gshape)))

  if (present(masked)) then
    g%masked = masked
  else
    g%masked = .false.
  end if

  ! data box

  call init_databox_G(g%db,g%n,g%gshape,g)
 end subroutine initgrid_nd

! macro to create a vector for an array
#define VEC(x) reshape(x, (/ size(x) /))

! create a 1-d grid based on a x and mask

 subroutine initgrid_1d(g,x,masked)
  implicit none
  type(grid), intent(out) :: g
  real, intent(in)        :: x(:)
  logical, intent(in)     :: masked(:)

  call initgrid_nd(g,1,shape(masked),VEC(masked))
  call setCoord(g,1,VEC(x))
 end subroutine initgrid_1d

! create a 2-d grid based on a x, y and mask

 subroutine initgrid_2d(g,x,y,masked)
  implicit none
  type(grid), intent(out) :: g
  real, intent(in)        :: x(:,:), y(:,:)
  logical, intent(in)     :: masked(:,:)

  call initgrid_nd(g,2,shape(masked),VEC(masked))
  call setCoord(g,1,VEC(x))
  call setCoord(g,2,VEC(y))
 end subroutine initgrid_2d

! create a 3-d grid based on a x, y, z and mask

 subroutine initgrid_3d(g,x,y,z,masked)
  implicit none
  type(grid), intent(out) :: g
  real, intent(in)        :: x(:,:,:), y(:,:,:), z(:,:,:)
  logical, intent(in)     :: masked(:,:,:)

  call initgrid_nd(g,3,shape(masked),VEC(masked))
  call setCoord(g,1,VEC(x))
  call setCoord(g,2,VEC(y))
  call setCoord(g,3,VEC(z))
 end subroutine initgrid_3d

! create a 4-d grid based on a x, y, z, t and mask

 subroutine initgrid_4d(g,x,y,z,t,masked)
  implicit none
  type(grid), intent(out) :: g
  real, intent(in)        :: x(:,:,:,:), y(:,:,:,:), z(:,:,:,:), t(:,:,:,:)
  logical, intent(in)     :: masked(:,:,:,:)

  call initgrid_nd(g,3,shape(masked),VEC(masked))
  call setCoord(g,1,VEC(x))
  call setCoord(g,2,VEC(y))
  call setCoord(g,3,VEC(z))
  call setCoord(g,4,VEC(t))
 end subroutine initgrid_4d



!
! Finalization
!

 subroutine donegrid(g)
  implicit none
  type(grid), intent(out) :: g


  call done_databox(g%db)

  deallocate(g%gshape,g%ioffset,g%ioffset_mask,g%dependence,g%startindex,g%endindex)
  deallocate(g%data)
  deallocate(g%tetrahedron)
  deallocate(g%masked)

 end subroutine



 subroutine setCoord_fromVariable(g,nd,x)
  implicit none
  type(grid), intent(inout) :: g
  integer, intent(in) :: nd
  real, intent(in) :: x(*)

  integer :: i,n,totsize

  do i=g%startindex(nd),g%endindex(nd)
    g%data(i) = x(i-g%startindex(nd)+1)
  end do

 end subroutine 


 subroutine setCoord_fromFile(g,nd,fname)
  use ufileformat
  implicit none
  type(grid), intent(inout) :: g
  integer, intent(in) :: nd
  character(*), intent(in) :: fname

  integer :: ndim,gshape(MaxDimensions),prec
  real :: valex
  real, pointer :: x(:,:,:)


  call ureadfull_srepmat(fname,g%data(g%startindex(nd):),valex,force_shape=g%gshape)
   
 end subroutine 



! ind = one-based index

function getCoord_grid(g,ind,out) result(x)
implicit none
  class(grid), intent(in) :: g
  integer, intent(in) :: ind(:)
  logical, optional, intent(out) :: out
  real :: x(size(ind))

  !dbg(ind)
  if (present(out)) then
    x = getCoord0(g,ind-1,out)
  else
    x = getCoord0(g,ind-1)
  end if
end function



! ind = zero-based index

function getCoord0(g,ind,out) result(x)
implicit none
  type(grid), intent(in) :: g
  integer, intent(in) :: ind(:)
  logical, optional, intent(out) :: out
  real :: x(size(ind))

  integer :: i,linindex

!  dbg(ind)
!  dbg(g%gshape)
  if (present(out)) then
    out = any(ind < 0).or.any(ind >= g%gshape)
    if (out) return
  end if

#ifdef DEBUG
  if (any(ind < 0).or.any(ind >= g%gshape)) then
    write(6,*) 'index (0-based) out of bound ',ind,' shape ',g%gshape
    call abort()
  end if
#endif

  do i=1,g%n
    linindex = g%startindex(i) + sum(ind * g%ioffset(:,i))
    x(i) = g%data(linindex)
  end do
  
  
end function

subroutine locate_grid(g,x,ind,out)
implicit none
  class(grid), intent(inout) :: g
  real, intent(in) :: x(:)
  integer, intent(out) :: ind(size(x))
  logical, intent(out) :: out

  call locate_databox_G(g%db,g%n,g%gshape,g%ioffset_mask,g%tetrahedron,g,x,ind,out)


end subroutine

! compute the interpolation coefficient and indices at the location xi 
! of a field defined on the grid g. 
! Input: 
!   g (type(grid)): grid
!   xi: coordinates of the data point
! Output:
!   indeces: (size n by 2**n) indexes of field (1-based)
!   coeff (size 2**n): interpolation coefficient
!   nbp: number of grid points used to compute the interpolated value
!      nbp is zero if the xi is out of grid.


subroutine cinterp(g,xi,indexes,coeff,nbp)
 implicit none
 class(basegrid), intent(inout) :: g
 real, intent(in) :: xi(:)
 integer, intent(out) :: indexes(:,:),nbp
 real, intent(out) :: coeff(:)

 ! local variables

 integer :: lentot,twon,d,linindex,i,l,j,k

 ! ind = zero-based indexes

 integer, dimension(g%n) :: ind
 integer, dimension(g%n,2**g%n) :: pind
 real, dimension(g%n,2**g%n) :: px
 logical, dimension(2**g%n) :: pmasked
 logical :: out

 lentot = product(g%gshape-1)
 twon = 2**g%n
 nbp = 0

 call g%locate(xi,ind,out) 
 !write(6,*) 'locate ',xi,ind,out

 if (.not.out) then

   ! compute all corner points of the hypercube: px

   do j = 1,twon
     do k=1,g%n
       if  (btest(j-1,k-1) .and. g%gshape(k) > 1) then
         indexes(k,j) = ind(k) + 1
       else
         indexes(k,j) = ind(k) 
       end if
     end do

     !write(6,*) 'indexes(i,:,j) ',i,indexes(i,:,j)

     linindex = sum((indexes(:,j)) * g%ioffset_mask) + 1
     pmasked(j) = g%masked(linindex)
     !px(:,j) = getCoord0(g,indexes(:,j))
     px(:,j) = g%getCoord(indexes(:,j)+1)

     !write(6,*) 'point ',px(:,j)

   end do

   out = any(pmasked)
   ! write(6,*) 'pmasked ',pmasked,linindex,g%gshape
   ! write(6,*) 'px ',px
   ! write(6,*) 'out ',out
   
   if (.not.out)  then
     call interp_cube(g%n,g%tetrahedron,px,xi,out,coeff)
     nbp = twon
     ! write(6,*) 'coeff ',coeff
     ! write(6,*) 'out ',out
     ! write(6,*) 'te ',g%tetrahedron


   end if
 end if

! convert zero-based indexes to one-based indexes

  indexes = indexes+1

  !do j = 1,twon
  !  write(6,*) 'indexes [',indexes(:,j),']',coeff(j)
  !end do
  
end subroutine cinterp

! compute the interpolated value fi at the location xi of a field f defined 
! on the grid g. 
! Input: 
!   g (type(grid)): grid of f
!   f: fields to interpolated (collapsed into a vector)
!   xi: location of the interpolated value
! Output:
!   fi: interpolated value
!   out: true where value out of grid

subroutine interp(g,f,xi,fi,out)
 implicit none
 class(basegrid), intent(inout) :: g
 real, intent(in) :: f(*),xi(:)
 real, intent(out) :: fi
 logical, intent(out) :: out

 real :: coeff(2**g%n)
 integer :: i, nbp, indexes(g%n,2**g%n), linindex

 call cinterp(g,xi,indexes,coeff,nbp)
 out = nbp.eq.0

 fi = 0
 do i=1,nbp

   if (any(indexes(:,i).le.0).or.any(indexes(:,i).gt.g%gshape)) then
     write(0,*) 'invalid index ',i,g%n
     write(0,*) 'indexes(:,i) ',indexes(:,i)
     ERROR_STOP
   end if

   linindex = sum((indexes(:,i)-1) * g%ioffset_mask) + 1
   fi=fi + coeff(i)*f(linindex)
 end do
end subroutine interp

 !_______________________________________________________
 !

subroutine interp_c(g,f,xi,fi,out,indexes,coeff,nbp)
 implicit none
 class(grid), intent(inout) :: g
 real, intent(in) :: f(:),xi(:)
 real, intent(out) :: fi
 logical, intent(out) :: out

 real, intent(out) :: coeff(2**g%n)
 integer, intent(out) :: indexes(g%n,2**g%n),  nbp

 integer :: i,linindex

 call cinterp(g,xi,indexes,coeff,nbp)
 out = nbp.eq.0

 fi = 0
 do i=1,nbp

   if (any(indexes(:,i).le.0).or.any(indexes(:,i).gt.g%gshape)) then
     write(0,*) 'invalid index ',i,g%n
     write(0,*) 'indexes(:,i) ',indexes(:,i)
     ERROR_STOP
   end if

   linindex = sum((indexes(:,i)-1) * g%ioffset_mask) + 1
   fi=fi + coeff(i)*f(linindex)
 end do
end subroutine interp_c

 !_______________________________________________________
 !

subroutine interp_field(grid1,field1,grid2,field2,outgrid)
 implicit none
 class(grid), intent(inout) :: grid1, grid2
 real, intent(in) :: field1(:)
 logical, intent(out), optional :: outgrid(:)
 real, intent(out) :: field2(:)

 logical :: out
 integer :: i2,ind(grid2%n),d
 real :: x2(grid2%n)


 do i2=1,product(grid2%gshape)
! do i2=70401,70401
   if (.not.grid2%masked(i2)) then
     !  i2 -> ind

     ind(1) = i2-1
     do d=grid2%n,2,-1
       ind(d) = ind(1)/grid2%ioffset_mask(d)
       ind(1) = ind(1)-ind(d)*grid2%ioffset_mask(d)
     end do
     !ind = ind+1
     x2 = getCoord0(grid2,ind,out)

     if (.not.out) then
       call interp(grid1,field1,x2,field2(i2),out)
!     write(6,*) 'ind ',i2, ind, out, x2,field2(i2)

     end if

     if (present(outgrid)) outgrid(i2) = out
   else
     field2(i2)=0
   end if
 end do

end subroutine interp_field

 !_______________________________________________________
 !

subroutine interp_field_c(grid1,field1,grid2,field2,outgrid,indexes,coefficients,nz)
 implicit none
 class(grid), intent(inout) :: grid1, grid2
 real, intent(in) :: field1(:)
 logical, intent(out), optional :: outgrid(:)
 real, intent(out) :: field2(:)
 integer, intent(out) :: indexes(:,:)
 real, intent(out) :: coefficients(:)

 logical :: out
 integer :: i2,ind(grid2%n),d
 real :: x2(grid2%n)


 real :: tcoeff(2**grid2%n)
 integer :: tindexes(grid2%n,2**grid2%n),  tnbp, n, nz
 integer :: linindex1,linindex2

 nz = 0
 n = grid2%n

 do i2=1,product(grid2%gshape)
! do i2=70401,70401
   if (.not.grid2%masked(i2)) then
     !  i2 -> ind

     ind(1) = i2-1
     do d=n,2,-1
       ind(d) = ind(1)/grid2%ioffset_mask(d)
       ind(1) = ind(1)-ind(d)*grid2%ioffset_mask(d)
     end do

     x2 = getCoord0(grid2,ind,out)

     ! convert zero-based index to one-based index
     ind = ind+1

     if (.not.out) then
!       call interp_c(grid1,field1,x2,field2(i2),out,tindexes,tcoeff,tnbp)
       call cinterp(grid1,x2,tindexes,tcoeff,tnbp)
     end if

       if (out.or.tnbp.eq.0) then
        tnbp=1
        ! indexes of destination grid
        indexes(1:n,nz+1) = ind
        indexes(n+1:2*n,nz+1) = -1
        coefficients(nz+1) = 0.
       else
         !write(6,*) 'shape(indexes) ',shape(indexes) , nz,n,tnbp
         ! indexes of destination grid
         indexes(1:n,nz+1:nz+tnbp) = spread(ind,2,tnbp)
         ! indexes of source grid
         indexes(n+1:2*n,nz+1:nz+tnbp) = tindexes(:,1:tnbp)
         coefficients(nz+1:nz+tnbp) = tcoeff(1:tnbp)  
       end if
!     write(6,*) 'ind ',i2, ind, out, x2,field2(i2)
!     write(6,*) 'nz ',nz
      nz = nz+tnbp
!     end if

     if (present(outgrid)) outgrid(i2) = out
   else
     field2(i2)=0
   end if
 end do

end subroutine interp_field_c

 !_______________________________________________________
 !

subroutine interp_field_coeff(indexes,coefficients,ioffset1,field1,ioffset2,field2,outgrid,nz)
 implicit none
 integer, intent(in) :: ioffset1(:),ioffset2(:)

 real, intent(in) :: field1(:)
 logical, intent(out), optional :: outgrid(:)
 real, intent(out) :: field2(:)
 integer, intent(in) :: indexes(:,:)
 real, intent(in) :: coefficients(:)

 logical :: out
 integer :: i2,ind(size(ioffset2)),d
 real :: x2(size(ioffset2))


 real :: tcoeff(2**size(ioffset2))
 integer :: tindexes(size(ioffset2),2**size(ioffset2)),  tnbp, n, nz,i
 integer :: linindex1,linindex2

 n = size(ioffset2)

if (present(outgrid)) outgrid = .true.
 field2 = 0.
! write(6,*) 'sum c',sum(coefficients),nz
 do i=1,nz
   linindex2 = sum((indexes(1:n,i) -1) *     ioffset2) + 1

!   if (indexes(i,n+1).ne.-1) then
   if (indexes(n+1,i).ne.-1) then
     linindex1 = sum((indexes(n+1:2*n,i) -1) * ioffset1) + 1
!     linindex1 = sum((indexes(i,n+1:2*n) -1) * ioffset1) + 1
     field2(linindex2) = field2(linindex2) + field1(linindex1) * coefficients(i)

     if (present(outgrid)) outgrid(linindex2) = .false.
   else
     if (present(outgrid)) outgrid(linindex2) = .true.
   end if
 end do


end subroutine 


 ! create a grid of cells
 function setupgrid(xpos,dx) result(cg)
  ! xpos(m,n) xpos(i,:) are the n-dimensional coordinates of the i-th point
  real, intent(in) :: xpos(:,:), dx(:)
  real, allocatable :: xmax(:)
  integer, allocatable :: gridind(:)
  integer :: n, i, j, l, Nt
  type(cellgrid) :: cg
  ! fractional index
  real :: fracindex(size(xpos,2))
  integer :: startsize = 100, ci
  real :: eps

  ! small balue to make the upper-bound slightly larger to ensure that 
  ! the largerst value is also inside a cell
  if (kind(dx) == 4) then
    eps = 1e-3 * minval(dx)
  else
    eps = 1e-6 * minval(dx)
  end if

  n = size(xpos,2)
  cg%n = n

  allocate(cg%xmin(n),xmax(n),cg%dx(n),cg%Ni(n),cg%offset(n),gridind(n))

  cg%dx = dx

  do i = 1,n
    cg%xmin(i) = minval(xpos(:,i))
    xmax(i) = maxval(xpos(:,i)) + eps
  end do

  ! number of intervales in each dimensions
  cg%Ni = ceiling((xmax - cg%xmin) / cg%dx)

  cg%offset(1) = 1
  do i = 1,n-1
    cg%offset(i+1) = cg%offset(i) * cg%Ni(i)
  end do

  ! total number of cells
  Nt = product(cg%Ni)

  ! the j(1),j(2),...j(n) box is surrounded by the grid points
  ! lower: xmin(i) + dx(i) * (j(i)-1) for i=1,...,n
  ! upper: xmin(i) + dx(i) * j(i) for i=1,...,n

  allocate(cg%gridn(Nt))
  do i = 1,Nt
    allocate(cg%gridn(i)%ind(startsize))
    cg%gridn(i)%n = 0
  end do

  ! loop throught every coordinate and put index into a cell
  do l = 1,size(xpos,1)
    ! compute grid index
    fracindex = (xpos(l,:) - cg%xmin)/dx
    gridind = floor(fracindex) + 1
    !write(6,*) 'xpos',xpos(l,:),cg%xmin,gridind, (xpos(l,:) - cg%xmin)/ dx

    ci = sum((gridind-1) * cg%offset) + 1    
    call add(cg%gridn(ci),l)

    !if (mod(l,100) == 0) write(6,*) 'l ',l,gridind         
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


 ! get all points (or more!) near the location x
 
 subroutine near(cg,x,xpos,distfun,maxdist,ind,dist,n)
  implicit none
  type(cellgrid), intent(in) :: cg
  real, intent(in) :: x(:),xpos(:,:)
  procedure(distind) :: distfun
  real, intent(in) :: maxdist
  real, intent(out) :: dist(:)
  integer, intent(out) :: ind(:), n

  ! already checked grid indices
  logical :: ischecked(product(cg%Ni))
  integer :: gridind(cg%n)

  n = 0
  ischecked = .false.
  
  gridind = floor((x - cg%xmin)/ cg%dx) + 1
  call append(gridind)

 contains
  recursive subroutine append(gridind)
   integer, intent(in) :: gridind(:)
   integer :: l,nc,ci, j,k,twon
   integer :: gi(cg%n)

   real :: xb(cg%n)
   logical :: close_corners
   !     real, allocatable :: dist(:)

   !   write(6,*) 'gridindex', gridind

   twon = 2**cg%n

   ! check if gridind is valid
   do l = 1,cg%n
     if (gridind(l) < 1 .or. gridind(l) > cg%Ni(l)) return
   end do

   ! cell index
   ci = sum((gridind-1) * cg%offset) + 1

   ! check if already checked
   if (ischecked(ci)) return

   ! it's a fresh one
   ischecked(ci) = .true.

   ! check if any corner points of cell grid are within distance
   close_corners = .false.
   do j = 1,twon
     do k=1,cg%n
       if (btest(j-1,k-1)) then
         xb(k) = cg%xmin(k) + (gridind(k)-1) * cg%dx(k)
       else
         xb(k) = cg%xmin(k) + (gridind(k)) * cg%dx(k)
       end if

       close_corners = close_corners .or. distfun(x,xb) < maxdist
       !write(6,*) 'distfun(x,xb) < maxdist',distfun(x,xb),maxdist,distfun(x,xb) < maxdist
     end do
   end do

   if (.not. close_corners) then
     !write(6,*) 'no close corners'
     return
   end if

   nc = cg%gridn(ci)%n

   if (n+nc > size(dist)) then
     write(0,*) __FILE__,':',__LINE__,'buffer too small'
     stop
   end if

   do l = 1,nc
     dist(n+l) = distfun(x,xpos(cg%gridn(ci)%ind(l),:))
     !     write(6,*) 'gridindex', gridind,cg%Ni,dist(n+l)
   end do

   ! check if any grid point is near x
   if (any(dist(n+1:n+nc) < maxdist)) then
     ! ok add indices to the list

     ind(n+1:n+nc) = cg%gridn(ci)%ind(1:nc)
     n = n+nc
   end if

   ! recursively check neighbors
   
   do l = 1,cg%n
     gi = gridind
     gi(l) = gridind(l) + 1
     call append(gi)
     
     gi = gridind
     gi(l) = gridind(l) - 1
     call append(gi)
   end do
   
  end subroutine append
 end subroutine near


 ! check results with exhaustive search
 subroutine checknear(cg,x,xpos,distfun,maxdist,ind)
  use matoper
  implicit none
  type(cellgrid), intent(in) :: cg
  real, intent(in) :: x(:),xpos(:,:)
  procedure(distind) :: distfun
  real, intent(in) :: maxdist
  integer, intent(in) :: ind(:)

  integer :: found,l
  real :: distl
  logical :: success 

  success = .true.
  found = 0
  do l = 1,size(xpos,1)
    !if (mod(l,1000) == 0) write(6,*) 'l ',l,size(xpos,1)

    distl = distfun(x,xpos(l,:))

    if (distl < maxdist) then
      if (any(ind == l)) then
        found = found+1
      else
        write(6,*) 'not found ',l,distl,maxdist,x,xpos(l,:),'FAIL'
        success = .false.
        exit
        stop
      end if
    end if
  end do
  
  call assert(success,'find all')
  !write(6,*) 'found all ',found,size(ind),'OK'

 end subroutine checknear

 ! return all indices inside the cell to which x also belongs

 subroutine get(cg,x,ind,n)
  type(cellgrid), intent(in) :: cg
  real, intent(in) :: x(:)
  integer, intent(out) :: ind(:), n

  integer :: gridind(cg%n),ci

  gridind = floor((x - cg%xmin)/ cg%dx) + 1

  if (any(gridind < 0) .or. any(gridind > cg%Ni)) then
    ! out of grid
    n = 0
    return
  end if

  ! cell index
  ci = sum((gridind-1) * cg%offset) + 1
  n = cg%gridn(ci)%n

  if (n > size(ind)) then
    ! buffer too small
    n = -n
    return
  end if

  ind(:n) = cg%gridn(ci)%ind(1:n)

 end subroutine get



 ! cartesian distance
 pure real function cdist(x,y)
  real, intent(in) :: x(:),y(:)
  cdist = sqrt(sum((x-y)**2))
 end function cdist

end module ndgrid
