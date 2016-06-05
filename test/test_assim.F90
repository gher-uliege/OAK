! Copyright(c) 2002-2016 Alexander Barth and Luc Vandenblucke

! cd ~/Assim/OAK-nonDiagR &&  make test/test_assim && test/test_assim

program test_assim
 use assimilation
 use matoper
 use ufileformat
 use initfile
 implicit none
 
 integer, parameter :: imax = 10, jmax = 15
 
 real :: x(imax,jmax), y(imax,jmax), z(imax,jmax)
 integer :: mask(imax,jmax)

 integer :: i,j,k

 character(len=MaxFNameLength) :: path
 character(len=MaxFNameLength), pointer :: filenames(:)

 mask = 1
 do j = 1,jmax
   do i = 1,imax
     x(i,j) = i
     y(i,j) = j
   end do
 end do

 initfname = 'test/test_assim.init'
 call getInitValue(initfname,'Model.mask',filenames)
 call getInitValue(initfname,'Model.path',path)
 
 do i = 1,size(filenames)
   call usave(trim(path)//filenames(i),real(mask),DefaultValex)
 end do

 call MemoryLayout('Model.',ModML)

 

 
end program test_assim
