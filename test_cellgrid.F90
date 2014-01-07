module setgrid

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

 ! create a grid of cells
 function setupgrid(xpos,dx) result(cg)
  real, intent(in) :: xpos(:,:), dx(:)
  real, allocatable :: xmax(:)
  integer, allocatable :: gridind(:)
  integer :: n, i, j, l, Nt
  type(cellgrid) :: cg

  integer :: startsize = 100, ci

  n = size(xpos,2)
  cg%n = n

  allocate(cg%xmin(n),xmax(n),cg%dx(n),cg%Ni(n),cg%offset(n),gridind(n))

  cg%dx = dx

  do i = 1,n
    cg%xmin(i) = minval(xpos(:,i))
    xmax(i) = maxval(xpos(:,i))
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
    gridind = floor((xpos(l,:) - cg%xmin)/ dx) + 1

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


 ! get all points (one more!) near the location x
 
 subroutine near(cg,x,xpos,distfun,maxdist,ind,dist,n)
  type(cellgrid), intent(in) :: cg
  real, intent(in) :: x(:),xpos(:,:)
  procedure(distind) :: distfun
  real, intent(in) :: maxdist
  real, intent(out) :: dist(:)
  integer, intent(out) :: ind(:), n

  ! appended grid index
  integer :: appended(2,product(cg%Ni)), nappended
  integer :: gridind(cg%n)

  nappended = 0
  n = 0

  gridind = floor((x - cg%xmin)/ cg%dx) + 1
  call append(gridind)

 contains
  recursive subroutine append(gridind)
   integer, intent(in) :: gridind(:)
   integer :: l,nc,ci
   integer :: gi(cg%n)

   !     real, allocatable :: dist(:)

   !   write(6,*) 'gridindex', gridind
   ! check if gridind is valid
   do l = 1,cg%n
     if (gridind(l) < 1 .or. gridind(l) > cg%Ni(l)) return
   end do

   ! check if already appended
   do l = 1,nappended
     if (all(gridind == appended(:,l))) return
   end do

   ! it's a fresh one

   ! cell index
   ci = sum((gridind-1) * cg%offset) + 1
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

     ! all cell index to the list of already visited cells
     nappended = nappended+1
     appended(:,nappended) = gridind

     ! recursively check neighbors

     do l = 1,cg%n
       gi = gridind
       gi(l) = gridind(l) + 1
       call append(gi)

       gi = gridind
       gi(l) = gridind(l) - 1
       call append(gi)
     end do
   end if
  end subroutine append
 end subroutine near


 ! check results with exhaustive search
 subroutine checknear(cg,x,xpos,distfun,maxdist,ind,dist)
  type(cellgrid), intent(in) :: cg
  real, intent(in) :: x(:),xpos(:,:)
  procedure(distind) :: distfun
  real, intent(in) :: maxdist
  real, intent(in) :: dist(:)
  integer, intent(in) :: ind(:)

  integer :: found,l
  real :: distl

  found = 0
  do l = 1,size(xpos,1)
    !if (mod(l,1000) == 0) write(6,*) 'l ',l,size(xpos,1)

    distl = distfun(x,xpos(l,:))

    if (distl < maxdist) then
      if (any(ind == l)) then
        found = found+1
      else
        write(6,*) 'not found ',l,distl,x
        stop
      end if
    end if
  end do

  write(6,*) 'found all ',found,size(ind)

 end subroutine checknear

 subroutine get(cg,x,ind,n)
  type(cellgrid), intent(in) :: cg
  real, intent(in) :: x(:)
  integer, intent(out) :: ind(:), n

  integer :: gridind(cg%n),ci

  gridind = floor((x - cg%xmin)/ cg%dx) + 1

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


 pure real function cdist(x,y)
  real, intent(in) :: x(:),y(:)
  cdist = sqrt(sum((x-y)**2))
 end function cdist

 subroutine test_near(sz,maxdist)
  integer, intent(in) :: sz(:)
  real, intent(in) :: maxdist

  real, pointer :: xpos(:,:)
  type(cellgrid) :: cg

  integer :: i,j,k,l, nind, Nsz
  real :: x(2)

  integer, allocatable :: ind(:)
  real, allocatable :: dist(:)

  Nsz = product(sz)
  allocate(xpos(Nsz,2),ind(Nsz),dist(Nsz))

  !   call loadVector('Model.gridX',ModML,xpos(:,1))
  !   call loadVector('Model.gridY',ModML,xpos(:,2))
  l = 0

  if (size(sz) == 2) then
    do j = 1,sz(2)
      do i = 1,sz(1)
        l = l+1
        xpos(l,1) = 2*(i-10)
        xpos(l,2) = j
      end do
    end do
  else
    do k = 1,sz(3)
      do j = 1,sz(2)
        do i = 1,sz(1)
          l = l+1
          xpos(l,1) = i
          xpos(l,2) = j
          xpos(l,3) = k
        end do
      end do
    end do
  end if


  cg = setupgrid(xpos,[4.,4.])

  call get(cg,[11.,11.],ind,nind)
  !  write(6,*) 'nind ',nind
  !  write(6,*) 'ind ',ind(:nind)


  x = [2.,2.]

  call near(cg,x,xpos,cdist,maxdist,ind,dist,nind)
  !  write(6,*) 'nind ',nind
  !  write(6,*) 'ind ',ind(:nind)

  call checknear(cg,x,xpos,cdist,maxdist,ind,dist)

 end subroutine test_near

 subroutine test_large()
  use matoper
  use rrsqrt
  use ufileformat
  use initfile
  use assimilation
  use ndgrid

  real :: maxdist = 2e3

  real, pointer :: xpos(:,:)
  type(cellgrid) :: cg

  integer :: i,j,k,l, nind, Nsz
  integer :: found
  real :: x(2)

  integer, allocatable :: ind(:)
  real, allocatable :: dist(:)
  character(len=MaxFNameLength) :: str
  real :: start, finish
  integer :: ntot  = 0
  integer :: nsearch = 10000

  call getarg(1,str); call init(str)
  Nsz = ModML%effsize
  allocate(xpos(Nsz,2),ind(Nsz),dist(Nsz))

  call loadVector('Model.gridX',ModML,xpos(:,1))
  call loadVector('Model.gridY',ModML,xpos(:,2))

  cg = setupgrid(xpos,[0.1,0.1])

  x = [9.,43.]

  write(6,*) 'Optimized search (mean over ',nsearch,' seaches)'

  call cpu_time(start)
  do i = 1,nsearch
    x(1) = 9 + i / real(nsearch)
    call near(cg,x,xpos,distance,maxdist,ind,dist,nind)

    ntot = ntot + nind
  end do
  call cpu_time(finish)
  print '("Time = ",f9.6," seconds.")',(finish-start)/nsearch
!  write(6,*) 'nind ',nind

  write(6,*) 'ntot ',ntot


  write(6,*) 'Non-optimized search'

  call cpu_time(start)
  call checknear(cg,x,xpos,distance,maxdist,ind(1:nind),dist(1:nind))
  call cpu_time(finish)
  print '("Time = ",f9.6," seconds.")',(finish-start)

 end subroutine test_large
end module setgrid


program test_cellgrid
 use setgrid
 integer :: sz(2)

 ! call test_near([20,20],5.)
 ! call test_near([1000,1000],20.)
 call test_large()


end program test_cellgrid
