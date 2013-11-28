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

 interface
   real function distind(x,y)
    real, intent(in) :: x(:),y(:)
   end function distind
 end interface


contains

 function setupgrid(modGrid,dx) result(cg)
  real, intent(in) :: modGrid(:,:), dx(:)
  real, allocatable :: xmax(:)
  integer, allocatable :: gridind(:)
  integer :: n, i, j, l
  type(cellgrid) :: cg

  integer :: startsize = 100

  n = size(modGrid,2)
  cg%n = n

  allocate(cg%xmin(n),xmax(n),cg%dx(n),cg%Ni(n),gridind(n))

  cg%dx = dx

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
  allocate(cg%grid(cg%Ni(1),cg%Ni(2)))
  do j = 1,cg%Ni(2)
    do i = 1,cg%Ni(1)
      allocate(cg%grid(i,j)%ind(startsize))
      cg%grid(i,j)%n = 0
    end do
  end do

  ! loop throught every coordinate and put index into a cell
  do l = 1,size(modGrid,1)
    ! compute grid index
    gridind = floor((modGrid(l,:) - cg%xmin)/ dx) + 1
    call add(cg%grid(gridind(1),gridind(2)),l)
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


 subroutine near(cg,x,modGrid,distfun,maxdist,ind,dist,n)
  type(cellgrid), intent(in) :: cg
  real, intent(in) :: x(:),modGrid(:,:)
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
!  call append(gridind)
  call append([9,1])

 contains
  recursive subroutine append(gridind)
   integer, intent(in) :: gridind(:)
   integer :: l,nc
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
   nc = cg%grid(gridind(1),gridind(2))%n

   if (n+nc > size(dist)) then
     write(0,*) __FILE__,':',__LINE__,'buffer too small'
     stop
   end if

   do l = 1,nc
     dist(n+l) = distfun(x,modGrid(cg%grid(gridind(1),gridind(2))%ind(l),:))
     !     write(6,*) 'gridindex', gridind,cg%Ni,dist(n+l)
   end do

   ! check if any grid point is near x
   if (any(dist(n+1:n+nc) < maxdist)) then
     ! ok add indices to the list
     write(6,*) 'check dist ',any(dist(n+1:n+nc) < maxdist),gridind

     ind(n+1:n+nc) = cg%grid(gridind(1),gridind(2))%ind(1:nc)
     n = n+nc

     ! all cell index to the list of already visited cells
     nappended = nappended+1
     appended(:,nappended) = gridind

     stop
     ! recursively check neighbors
     call append([gridind(1)+1,gridind(2)])
     call append([gridind(1)-1,gridind(2)])
     call append([gridind(1),gridind(2)+1])
     call append([gridind(1),gridind(2)-1])
   end if
  end subroutine append
 end subroutine near


 ! check results with exhaustive search
 subroutine checknear(cg,x,modGrid,distfun,maxdist,ind,dist)
  type(cellgrid), intent(in) :: cg
  real, intent(in) :: x(:),modGrid(:,:)
  procedure(distind) :: distfun
  real, intent(in) :: maxdist
  real, intent(in) :: dist(:)
  integer, intent(in) :: ind(:)

  integer :: found,l
  real :: distl

  found = 0
  do l = 1,size(modGrid,1)
    !if (mod(l,1000) == 0) write(6,*) 'l ',l,size(modGrid,1)

    distl = distfun(x,modGrid(l,:))

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
  integer :: gridind(cg%n)

  gridind = floor((x - cg%xmin)/ cg%dx) + 1
  n = cg%grid(gridind(1),gridind(2))%n

  if (n > size(ind)) then
    ! buffer too small
    n = -n
    return
  end if

  ind(:n) = cg%grid(gridind(1),gridind(2))%ind(1:n)
 end subroutine get


 pure real function distfun(x,y)
  real, intent(in) :: x(:),y(:)
  distfun = sqrt(sum((x-y)**2))
 end function distfun

 subroutine test_near(sz,maxdist)
  integer, intent(in) :: sz(:)
  real, intent(in) :: maxdist

  real, pointer :: modGrid(:,:)
  type(cellgrid) :: cg

  integer :: i,j,k,l, nind, Nsz
  real :: x(2)

  integer, allocatable :: ind(:)
  real, allocatable :: dist(:)

  Nsz = product(sz)
  allocate(modGrid(Nsz,2),ind(Nsz),dist(Nsz))

  !   call loadVector('Model.gridX',ModML,modGrid(:,1))
  !   call loadVector('Model.gridY',ModML,modGrid(:,2))
  l = 0

  if (size(sz) == 2) then
    do j = 1,sz(2)
      do i = 1,sz(1)
        l = l+1
        modGrid(l,1) = 2*(i-10)
        modGrid(l,2) = j
      end do
    end do
  else
    do k = 1,sz(3)
      do j = 1,sz(2)
        do i = 1,sz(1)
          l = l+1
          modGrid(l,1) = i
          modGrid(l,2) = j
          modGrid(l,3) = k
        end do
      end do
    end do
  end if


  cg = setupgrid(modGrid,[4.,4.])

  call get(cg,[11.,11.],ind,nind)
  !  write(6,*) 'nind ',nind
  !  write(6,*) 'ind ',ind(:nind)


  x = [2.,2.]

  call near(cg,x,modGrid,distfun,maxdist,ind,dist,nind)
  !  write(6,*) 'nind ',nind
  !  write(6,*) 'ind ',ind(:nind)

  call checknear(cg,x,modGrid,distfun,maxdist,ind,dist)

 end subroutine test_near

 subroutine test_large()
  use matoper
  use rrsqrt
  use ufileformat
  use initfile
  use assimilation

  real :: maxdist = 2e3

  real, pointer :: modGrid(:,:)
  type(cellgrid) :: cg

  integer :: i,j,k,l, nind, Nsz
  integer :: found
  real :: x(2)

  integer, allocatable :: ind(:)
  real, allocatable :: dist(:)
  character(len=MaxFNameLength) :: str
  real :: start, finish


  call getarg(1,str); call init(str)
  Nsz = ModML%effsize
  allocate(modGrid(Nsz,2),ind(Nsz),dist(Nsz))

  call loadVector('Model.gridX',ModML,modGrid(:,1))
  call loadVector('Model.gridY',ModML,modGrid(:,2))

  cg = setupgrid(modGrid,[0.1,0.1])

  x = [9.,43.]

  call cpu_time(start)
  do i = 1,100
    x(1) = 9 + i / 100.
    call near(cg,x,modGrid,distance,maxdist,ind,dist,nind)
    ! put code to test here
  end do
  call cpu_time(finish)
  print '("Time = ",f9.6," seconds.")',(finish-start)/100
  write(6,*) 'nind ',nind,distance([9.,43.],[10.,44.])

  call cpu_time(start)
  call checknear(cg,x,modGrid,distance,maxdist,ind(1:nind),dist(1:nind))
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
