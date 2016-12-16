! cd ~/Assim/OAK-nonDiagR; make DEBUG= PROFILING= test/test_cellgrid && time test/test_cellgrid 

program test_cellgrid
 integer :: sz(2)

 call test_near([4,4],2.)
 call test_near([20,20],3.)
 call test_near([1000,1000],5.)

#ifdef LARGE
 ! requires ~/matlab/LocEns/ligurian_sea_loc_assim.init
 call test_large()
#endif

contains

 subroutine test_near(sz,maxdist)
  use ndgrid

  integer, intent(in) :: sz(:)
  real, intent(in) :: maxdist

  real, pointer :: xpos(:,:)
  type(cellgrid) :: cg

  integer :: i,j,k,l, nind, Nsz
  real :: x(2)
  real :: start, finish

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
        xpos(l,1) = i
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


!  cg = setupgrid(xpos,[maxdist,maxdist]/2)
  cg = setupgrid(xpos,[maxdist,maxdist])

!  call get(cg,[11.,11.],ind,nind)
  !  write(6,*) 'nind ',nind
  !  write(6,*) 'ind ',ind(:nind)


  x = [2.,2.]

  call cpu_time(start)
  call near(cg,x,xpos,cdist,maxdist,ind,dist,nind)
!  write(6,*) 'nind ',nind
!  write(6,*) 'ind ',ind(:nind)
#ifdef PROFILE
  call cpu_time(finish)
  print '("Time = ",f9.6," seconds.")',(finish-start)
#endif

  call cpu_time(start)
  call checknear(cg,x,xpos,cdist,maxdist,ind(1:nind))
  call cpu_time(finish)
#ifdef PROFILE
  print '("Time = ",f9.6," seconds.")',(finish-start)
#endif
  ! x = sz

  ! call near(cg,x,xpos,cdist,maxdist,ind,dist,nind)
  ! !write(6,*) 'nind ',nind, maxdist
  ! !write(6,*) 'ind ',ind(:nind)

  ! call cpu_time(start)
  ! call checknear(cg,x,xpos,cdist,maxdist,ind(1:nind))
  ! call cpu_time(finish)
  ! print '("Time = ",f9.6," seconds.")',(finish-start)

 end subroutine test_near

#ifdef LARGE
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
  call checknear(cg,x,xpos,distance,maxdist,ind(1:nind))
  call cpu_time(finish)
  print '("Time = ",f9.6," seconds.")',(finish-start)

 end subroutine test_large
#endif

end program test_cellgrid
