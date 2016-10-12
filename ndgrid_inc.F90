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

!
! --------------------------------------------------------------------
!


! x(comp,data points)
! f(data points )

! xi (components)

subroutine interpgrid_coeff_VARIANT(n,gshape,coord,masked,xi,out,indexes,coeff &
#ifdef USER_LOCATE
 ,locate &
#endif
 )

 !  use matoper
 implicit none

 integer, intent(in) :: n
 integer, intent(in) :: gshape(n)
 logical, intent(in) :: masked(:)
 real, intent(in) :: xi(:,:)
 logical, intent(out) :: out(:)
 integer, intent(out) :: indexes(:,:,:)
 real, intent(out) :: coeff(:,:)

#ifdef  IMPLICIT_COORDINATE
 interface 
   function coord(indexes) result(x)
    integer, intent(in) :: indexes(:)
    real :: x(size(indexes))
   end function coord
 end interface
#elif defined(GRID_COORDINATE)
 type(grid), intent(inout) :: coord
# else
 real, intent(in) :: coord(:,:)
#endif

#ifdef USER_LOCATE
 interface 
   subroutine locate(x,indexes,out)
    real, intent(in) :: x(:)
    integer, intent(out) :: indexes(size(x))
    logical, intent(out) :: out
   end subroutine locate
 end interface
#endif

 ! local variables

 integer :: lentot,m,twon,d,linindex,nbth,i,l,j,k

 ! ind = zero-based indexes

 integer, dimension(n) :: ind
 integer, dimension(n) :: ioffset
 integer, dimension(n) :: ioffset2

 integer, dimension(n,2**n) :: pind
 real, dimension(n,2**n) :: px

 real, dimension(2**n,2**n) :: identity
 logical, dimension(2**n) :: pmasked

 real, pointer,save :: tetrahedron(:,:,:)


# ifdef DATABOX_SEARCH
 type(databoxND) :: databox
# endif

 lentot = product(gshape-1)
 m = size(xi,2)
 twon = 2**n
 nbth =  factorial(n)*2**(n-1) 

 allocate(tetrahedron(twon,n+1,nbth))
 !  allocate(ioffset(n),ioffset2(n))
 ! vectors to access the element of f(:)

 ioffset(1) = 1
 ioffset2(1) = 1

 do d=1,n-1
   ioffset(d+1) = ioffset(d)*(gshape(d))
   ioffset2(d+1) = ioffset2(d)*(gshape(d)-1)
 end do

 ! identity matrix

 identity = 0.
 do i=1,twon
   identity(i,i) = 1.
 end do

 tetrahedron = 0
 call split(n,n,(/ (.false.,i=1,n) /),(/ (.true.,i=1,twon) /),tetrahedron)

 !       write(6,*) 'shape(tetrahedron) ',allocated(tetrahedron)

# ifndef USER_LOCATE
# ifdef DATABOX_SEARCH
 write(6,*) 'initialise databox ',gshape
 call init_databox_VARIANT(databox,n,gshape,coord)
 write(6,*) 'databox done'
# endif
# endif

 ! loop over all points where we want to know fi

 do i=1,m
   !      write(6,*) 'al(tetrahedron) ',m,allocated(tetrahedron)
   !   write(6,*) 'xi ',xi(:,i)


#   ifdef USER_LOCATE
   call locate(xi(:,i),ind,out(i))
#   else
#   ifdef DATABOX_SEARCH
   call locate_databox_VARIANT(databox,n,gshape,ioffset,tetrahedron,coord,xi(:,i),ind,out(i))
#   else
   call locate_testall_VARIANT(n,gshape,ioffset,ioffset2,tetrahedron,coord,xi(:,i),ind,out(i))
#   endif
#   endif

   !    write(6,*) 'located in ',ind,out(i)  

   !   write(6,*) 'i,out ',i,out(i)

   if (.not.out(i)) then

     ! compute all corner points of the hypercube: px

     do j = 1,twon
       !        indexes(i,:,j) = corner_indexes(:,j) + ind

       do k=1,n
         if  (btest(j-1,k-1)) then
           indexes(i,k,j) = ind(k) + 1
         else
           indexes(i,k,j) = ind(k) 
         end if
       end do

       !write(6,*) 'indexes(i,:,j) ',i,indexes(i,:,j)

       linindex = sum((indexes(i,:,j)) * ioffset) + 1
       pmasked(j) = masked(linindex)
#      ifdef IMPLICIT_COORDINATE
       px(:,j) = coord(indexes(i,:,j))
#      elif defined(GRID_COORDINATE)
       px(:,j) = getCoord0(coord,indexes(i,:,j))
#      else
       px(:,j) = coord(:,linindex)
#      endif
     end do

     out(i) = any(pmasked)
     if (.not.out(i)) call interp_cube(n,tetrahedron,px,xi(:,i),out(i),coeff(i,:))

     !   write(6,*) 'out ',out
   end if

 end do

 deallocate(tetrahedron)

 write(6,*) 'inside ',count(.not.out)

end subroutine interpgrid_coeff_VARIANT

!
! --------------------------------------------------------------------
!

subroutine locate_testall_VARIANT(n,gshape,ioffset,ioffset2,tetrahedron,coord,xi,ind,out)
 implicit none
 integer, intent(in) :: n
 integer, intent(in) :: gshape(n),ioffset(n),ioffset2(n)
 real, intent(in) :: tetrahedron(:,:,:)
 real, intent(in) :: xi(n)

 integer, dimension(n) :: ind
 logical, intent(out) :: out

#ifdef  IMPLICIT_COORDINATE
 interface 
   function coord(indexes) result(x)
    integer, intent(in) :: indexes(:)
    real :: x(size(indexes))
   end function coord
 end interface
#elif defined(GRID_COORDINATE)
 type(grid), intent(inout) :: coord
# else
 real, intent(in) :: coord(:,:)
#endif

 integer :: l,d,j,linindex

 integer :: pind(n,2**n)
 logical :: pmasked(2**n)
 real, dimension(n,2**n) :: px

 out = .true.

 search_grid: do l=0,product(gshape-1)-1
   ! transform linear index to zero-based vector indexes

   ind(1) = l
   do d=n,2,-1
     ind(d) = ind(1)/ioffset2(d)
     ind(1) = ind(1)-ind(d)*ioffset2(d)
   end do

   out = .not.InCube_VARIANT(n,gshape,ioffset,tetrahedron,coord,xi,ind)

   if (.not.out) then
     exit search_grid
   end if
 end do search_grid

end subroutine locate_testall_VARIANT



!!$
!!$ !
!!$ ! --------------------------------------------------------------------
!!$ ! This routine in even slower

!!$ subroutine locate2(n,gshape,x,xi,ind,out)
!!$  implicit none
!!$  integer, intent(in) :: n
!!$  integer, intent(in) :: gshape(n)
!!$  real, intent(in) :: x(:,:)
!!$  real, intent(in) :: xi(n)
!!$
!!$  integer, dimension(n) :: ind
!!$  logical, intent(out) :: out
!!$
!!$  integer :: l,d,j,linindex
!!$
!!$  integer :: pind(n,2**n)
!!$  logical :: pmasked(2**n)
!!$  real, dimension(n,2**n) :: px
!!$
!!$  out = .false.
!!$
!!$    ! transform linear index to zero-based vector indexes
!!$  ind = 0
!!$
!!$  search_grid: do l=0,product(gshape-1)-1
!!$
!!$
!!$    out = InCube(n,gshape,x,xi,ind)
!!$
!!$    if (out) then
!!$      exit search_grid
!!$    end if
!!$
!!$    ! increment index
!!$    d = 1
!!$    ind(d) = ind(d)+1
!!$    do while (ind(d).eq.gshape(d)) 
!!$      ind(d) = 0
!!$      d =d+1
!!$      ind(d) = ind(d)+1
!!$    end do
!!$   end do search_grid
!!$
!!$  end subroutine locate2
!!$
!!$


!
! --------------------------------------------------------------------
!
! Test if point xi is in the hypercube formed by the 
! corners x(:,1), x(:,2), ... x(:,2**n)

function  InCube_VARIANT(n,gshape,ioffset,tetrahedron,coord,xi,ind) result(inside)
 implicit none
 ! dimension of the grid
 integer, intent(in) :: n
 ! shape and memory offset 
 integer, intent(in) :: gshape(n),ioffset(n)
 real, intent(in) :: tetrahedron(:,:,:)

 real, intent(in) :: xi(n)
 integer, intent(in) :: ind(n)
 logical :: inside

 logical :: out
 integer :: l,linindex,j
 real, dimension(n,2**n) :: px
 integer :: k,pind(n)

#ifdef  IMPLICIT_COORDINATE
 interface 
   function coord(indexes) result(x)
    integer, intent(in) :: indexes(:)
    real :: x(size(indexes))
   end function coord
 end interface
#elif defined(GRID_COORDINATE)
 type(grid), intent(inout) :: coord
# else
 real, intent(in) :: coord(:,:)
#endif

 ! compute all corner points of the hypercube: px

 !  write(6,*) 'corner_indexes ',corner_indexes
 !  write(6,*) 'commul ',ioffset

 do j = 1,2**n
   ! pind = corner_indexes(:,j) + ind

   do k=1,n
     if  (btest(j-1,k-1)) then
       pind(k) = min(ind(k) + 1,gshape(k)-1)
     else
       pind(k) = ind(k) 
     end if
   end do

#  ifdef IMPLICIT_COORDINATE
   px(:,j) = coord( pind  )
#  elif defined(GRID_COORDINATE)
   px(:,j) = getCoord0(coord,pind)
#  else
   linindex = sum( pind * ioffset) + 1
   px(:,j) = coord(:,linindex)
#  endif

 end do

 ! if (n.eq.4) write(6,*) 'px ',minval(px(3,:)),maxval(px(3,:))

 call interp_cube(n,tetrahedron,px,xi,out)
 inside = .not.out

 ! write(6,*) 'inside ',inside
end function InCube_VARIANT



!
! --------------------------------------------------------------------
!
! interpolate the data points f given on the grid x to the points xi
!
!
! Input:
!
! integer :: n          
!   number of dimensions
!
! integer :: gshape(n)  
!   gshape(i) is the maximum value of the ith index
!
! real    :: x(n,p)     
!   coordinate of data points: 
!   where p =  gshape(1)*gshape(2)*...*gshape(n), i.e. the number of 
!   nodes of the grid
!
! real    :: f(p)       
!   f(i) is value of the fields at the point (x(1,i),x(2,i),...,x(n,i))
!
! logical :: masked(p)  
!   Where masked is .false., the corresponding data point f is valid, 
!   where it is .true. the data points is not valid and will not be used
!
! real    :: xi(n,m)
!   coordinates of the points where we want to have a interpolated value
!
! Output
! real :: fi(m)
!   interpolated value
!
! logical :: out(m)
!   inside is .true. where the inteprolation was successful. If it is false, then the point 
!   lies out of the grid
!
! --------------------------------------------------------------------


subroutine interpgrid_VARIANT(n,gshape,coord,f,masked,xi,fi,out &
#ifdef USER_LOCATE
 ,locate &
#endif
 )

 implicit none
 integer, intent(in) :: n 
 integer, intent(in) :: gshape(n)          
 real, intent(in) :: f(:)                  
 logical, intent(in) :: masked(:)              
 real, intent(in) :: xi(:,:)               
 real, intent(out) :: fi(:)                
 logical, intent(out) :: out(:)            

#ifdef  IMPLICIT_COORDINATE
 interface 
   function coord(indexes) result(x)
    integer, intent(in) :: indexes(:)
    real :: x(size(indexes))
   end function coord
 end interface
#elif defined(GRID_COORDINATE)
 type(grid), intent(inout) :: coord
# else
 real, intent(in) :: coord(:,:)
#endif

#ifdef USER_LOCATE
 interface 
   subroutine locate(x,indexes,out)
    real, intent(in) :: x(:)
    integer, intent(out) :: indexes(size(x))
    logical, intent(out) :: out
   end subroutine locate
 end interface
#endif



 ! local variables

 integer, allocatable :: indexes(:,:,:)
 real, allocatable :: coeff(:,:)

 integer :: m,twon,linindex,i,j

 ! ioffset is used to transform a linear index to a 
 ! vector index
 integer, dimension(n) :: ioffset


 ! dimension

 m = size(xi,2)
 twon = 2**n

 ioffset(1) = 1
 do i=1,n-1
   ioffset(i+1) = ioffset(i)*(gshape(i))
 end do

 allocate(indexes(m,n,twon),coeff(m,twon))

 call interpgrid_coeff_VARIANT(n,gshape,coord,masked,xi,out,indexes,coeff &
#ifdef USER_LOCATE
 ,locate &
#endif
 )

 !  write(6,*) 'out ',out

 ! loop over all points where we want to know fi

 do i=1,m
   !write(6,*) 'i,out ',i,out(i)

   if (.not.out(i)) then

     fi(i) = 0.

     ! extract values of the fields at the cube corner points 
     do j = 1,twon
       if (any(indexes(i,:,j).lt.0).or.any(indexes(i,:,j).ge.gshape)) then
         write(0,*) 'invalid index ',i,j,n,shape(out)
         write(0,*) 'indexes(i,:,j) ',indexes(i,:,j)
         ERROR_STOP
       end if
       linindex = sum(indexes(i,:,j) * ioffset) + 1
       fi(i)=fi(i) + coeff(i,j)*f(linindex)
     end do
   end if
 end do

 deallocate(indexes,coeff)

end subroutine interpgrid_VARIANT








#ifdef DATABOX_SEARCH




!_______________________________________________________
!
! subroutine for the databox structure
!
!_______________________________________________________
!

!_______________________________________________________
!
!
! subroutine for initialising the databox structure,


subroutine init_databox_VARIANT(db,n,gshape,coord)
 implicit none
 type(databoxND), intent(out) :: db
 integer, intent(in) :: n
 integer, intent(in) :: gshape(n)

 integer :: i,ioffset(n)

#ifdef  IMPLICIT_COORDINATE
 interface 
   function coord(indexes) result(x)
    integer, intent(in) :: indexes(:)
    real :: x(size(indexes))
   end function coord
 end interface
#elif defined(GRID_COORDINATE)
 type(grid), intent(inout) :: coord
# else
 real, intent(in) :: coord(:,:)
#endif

 ioffset(1) = 1

 do i=1,n-1
   ioffset(i+1) = ioffset(i)*(gshape(i))
 end do


! call init_databox_helper_VARIANT(db,n,ioffset,coord,(/ (0,i=1,n) /),gshape-1)

 call define_databox_VARIANT(db,n,ioffset,coord,(/ (0,i=1,n) /),gshape-1)

end subroutine init_databox_VARIANT

! helper function for recursive calls
!
! devide grid into 2**n subgrids
!
! For n=2
!                                                                       imax
!  +----------------------------------+----------------------------------+
!  |                      subimax(:,2)|                      subimax(:,1)|
!  |                                  |                                  |
!  |                                  |                                  |
!  |                                  |                                  |
!  |                                  |                                  |
!  |                                  |                                  |
!  |             m = 2                |             m = 1                |
!  |                                  |                                  |
!  |                                  |                                  |
!  |                                  |                                  |
!  |                                  |                                  |
!  |                                  |                                  |
!  |subimin(:,2)                      |subimin(:,1)                      |
!  +----------------------------------+----------------------------------+
!  |                      subimax(:,4)|im                    subimax(:,3)|
!  |                                  |                                  |
!  |                                  |                                  |
!  |                                  |                                  |
!  |                                  |                                  |
!  |                                  |                                  |
!  |            m = 4                 |             m = 3                |
!  |                                  |                                  |
!  |                                  |                                  |
!  |                                  |                                  |
!  |                                  |                                  |
!  |                                  |                                  |
!  |subimin(:,4)                      |subimin(:,3)                      |
!  +----------------------------------+----------------------------------+
!imin
!
!
!  Each of these four subgrid further subdivided in four subgrids
!

recursive subroutine init_databox_helper_VARIANT(db,n,ioffset,coord,imin,imax)
 implicit none
 integer, intent(in) :: n
 integer, intent(in) :: ioffset(n)
 type(databoxND), intent(out) :: db
 integer,  intent(in) :: imin(n),imax(n)

#ifdef  IMPLICIT_COORDINATE
 interface 
   function coord(indexes) result(x)
    integer, intent(in) :: indexes(:)
    real :: x(size(indexes))
   end function coord
 end interface
#elif defined(GRID_COORDINATE)
 type(grid), intent(inout) :: coord
# else
 real, intent(in) :: coord(:,:)
#endif

 integer :: im(n),subimin(n,2**n),subimax(n,2**n),m,p,splitdim,tot,ioffset_subgrid(n)
 integer :: ind(n),d,i,l,linindex
 logical :: subs(2**n)
 real :: xc(n)

 allocate(db%imin(n),db%imax(n),db%xmin(n),db%xmax(n))
 ! make sure the associated status is defined
 nullify(db%db)


 db%imin = imin
 db%imax = imax

 ioffset_subgrid(1) = 1

 do i=1,n-1
   ioffset_subgrid(i+1) = ioffset_subgrid(i)*(imax(i)-imin(i)+1)
 end do


 db%xmin = huge(db%xmin)
 db%xmax = -huge(db%xmax)

 ! search maximum and minimum

 tot = product(imax-imin+1)
 do m=0,tot-1
   ind(1) = m
   do d=n,2,-1
     ind(d) = ind(1)/ioffset_subgrid(d)
     ind(1) = ind(1)-ind(d)*ioffset_subgrid(d)
   end do

   !
   !    ind = linindex2index(m,ioffset_subgrid) 

   ind = ind+imin

   !  write(6,*) 'ioffset2 ',ioffset
#  ifdef IMPLICIT_COORDINATE
   xc = coord(ind)
#  elif defined(GRID_COORDINATE)
   xc = getCoord0(coord,ind)
#  else
   linindex = sum(ind*ioffset)+1
   xc = coord(:,linindex)
#  endif



   do d=1,n
     if (xc(d).lt.db%xmin(d)) db%xmin(d) = xc(d)
     if (xc(d).gt.db%xmax(d)) db%xmax(d) = xc(d)
   end do
 end do

 !write(6,*) 'min2 ', db%xmin(1), db%xmin(2),  db%xmax(1), db%xmax(2)
 ! write(6,*) 'min2 ', db%xmin(1)

 !  if (db%imax(1)-db%imin(1).ne.1.or.db%imax(2)-db%imin(2).ne.1)  then
 if (any(db%imax-db%imin.ne.1))  then
   db%type = db_splitted

   ! middle of data box
   im = (db%imin+db%imax)/2

   do m=1,2**n
     do d=1,n
       if (btest(m-1,d-1)) then
         subimin(d,m) = imin(d)
         subimax(d,m) = im(d)
       else
         subimin(d,m) = im(d)
         subimax(d,m) = imax(d)
       end if
     end do
   end do

!!$     write(6,*) 'subimin ', subimin(1,:)
!!$     write(6,*) 'subimin ', subimin(2,:)
!!$     write(6,*) 'subimax ', subimax(1,:)
!!$     write(6,*) 'subimax ', subimax(2,:)
!!$    stop
   subs = all(subimax-subimin.ge.1,1)

   allocate(db%db(count(subs)))
   p=1
   do m=1,size(subimax,2)
     if (subs(m)) then
       call init_databox_helper_VARIANT(db%db(p),n,ioffset,coord,subimin(:,m),subimax(:,m))
       p=p+1
     end if
   end do
 else
   db%type = db_cell
 end if
end subroutine init_databox_helper_VARIANT





!
! define the folowing attributes of databox:
!   imin,imax,xmin,xmax
!
! but the databox is not splitted
!

recursive subroutine define_databox_VARIANT(db,n,ioffset,coord,imin,imax)
 implicit none
 integer, intent(in) :: n
 integer, intent(in) :: ioffset(n)
 type(databoxND), intent(out) :: db
 integer,  intent(in) :: imin(n),imax(n)

#ifdef  IMPLICIT_COORDINATE
 interface 
   function coord(indexes) result(x)
    integer, intent(in) :: indexes(:)
    real :: x(size(indexes))
   end function coord
 end interface
#elif defined(GRID_COORDINATE)
 type(grid), intent(inout) :: coord
# else
 real, intent(in) :: coord(:,:)
#endif

 integer :: m,tot,ioffset_subgrid(n)
 integer :: ind(n),d,i,l,linindex
 real :: xc(n)

 db%type = db_notsplitted

 allocate(db%imin(n),db%imax(n),db%xmin(n),db%xmax(n))
 db%imin = imin
 db%imax = imax

end subroutine


subroutine search_boundarybox_VARIANT(db,n,ioffset,coord)
 implicit none
 integer, intent(in) :: n
 integer, intent(in) :: ioffset(n)
 type(databoxND), intent(inout) :: db

#ifdef  IMPLICIT_COORDINATE
 interface 
   function coord(indexes) result(x)
    integer, intent(in) :: indexes(:)
    real :: x(size(indexes))
   end function coord
 end interface
#elif defined(GRID_COORDINATE)
 type(grid), intent(inout) :: coord
# else
 real, intent(in) :: coord(:,:)
#endif

 integer :: m,tot,ioffset_subgrid(n)
 integer :: ind(n),d,i,l,linindex
 real :: xc(n)

 ioffset_subgrid(1) = 1

 do i=1,n-1
   ioffset_subgrid(i+1) = ioffset_subgrid(i)*(db%imax(i)-db%imin(i)+1)
 end do


 db%xmin = huge(db%xmin)
 db%xmax = -huge(db%xmax)

 ! search maximum and minimum

 tot = product(db%imax-db%imin+1)
 do m=0,tot-1
   ind(1) = m
   do d=n,2,-1
     ind(d) = ind(1)/ioffset_subgrid(d)
     ind(1) = ind(1)-ind(d)*ioffset_subgrid(d)
   end do

   !
   !    ind = linindex2index(m,ioffset_subgrid) 

   ind = ind+db%imin

   !  write(6,*) 'ioffset2 ',ioffset
#  ifdef IMPLICIT_COORDINATE
   xc = coord(ind)
#  elif defined(GRID_COORDINATE)
   xc = getCoord0(coord,ind)
#  else
   linindex = sum(ind*ioffset)+1
   xc = coord(:,linindex)
#  endif


   do d=1,n
     if (xc(d).lt.db%xmin(d)) db%xmin(d) = xc(d)
     if (xc(d).gt.db%xmax(d)) db%xmax(d) = xc(d)
   end do
 end do

end subroutine

!
! split a databox into subboxes assuming imin,imax,xmin,xmax are already defined
! The attribute db(:) is alloctated and initialised by define_databox_VARIANT
! 


recursive subroutine split_databox_VARIANT(db,n,ioffset,coord)
 implicit none
 integer, intent(in) :: n
 integer, intent(in) :: ioffset(n)
 type(databoxND), intent(inout) :: db

#ifdef  IMPLICIT_COORDINATE
 interface 
   function coord(indexes) result(x)
    integer, intent(in) :: indexes(:)
    real :: x(size(indexes))
   end function coord
 end interface
#elif defined(GRID_COORDINATE)
 type(grid), intent(inout) :: coord
# else
 real, intent(in) :: coord(:,:)
#endif

! subimin/subimax are zero-based indices

 integer :: im(n),subimin(n,2**n),subimax(n,2**n),m,p
 integer :: ind(n),d,l
 logical :: subs(2**n)

 call  search_boundarybox_VARIANT(db,n,ioffset,coord)

 ! check if databox has to be splitted

! if (any(db%imax-db%imin.ne.1))  then
 if (any(db%imax-db%imin > 1))  then
   db%type = db_splitted

   ! middle of data box
   im = (db%imin+db%imax)/2

   do m=1,2**n
     do d=1,n
       if (btest(m-1,d-1)) then
         subimin(d,m) = db%imin(d)
         subimax(d,m) = im(d)
       else
         subimin(d,m) = im(d)
         subimax(d,m) = db%imax(d)
       end if

     !write(6,*) 'subbox ',m,d,subimin(d,m),subimax(d,m)

     end do
   end do

!!$     write(6,*) 'subimin ', subimin(1,:)
!!$     write(6,*) 'subimin ', subimin(2,:)
!!$     write(6,*) 'subimax ', subimax(1,:)
!!$     write(6,*) 'subimax ', subimax(2,:)
!!$    stop
!   subs = all(subimax-subimin.ge.1,1)
   subs = any(subimax-subimin >= 1,1)

   !write(6,*) 'subs ', subs
   allocate(db%db(count(subs)))
!   stop

   p=1
   do m=1,size(subimax,2)
     if (subs(m)) then
       call define_databox_VARIANT(db%db(p),n,ioffset,coord,subimin(:,m),subimax(:,m))
       p=p+1
     end if
   end do
 else
   db%type = db_cell
 end if
end subroutine 




!_______________________________________________________
!
! subroutine for localising a point within a databox 
!_______________________________________________________
!

recursive subroutine locate_databox_VARIANT(db,n,gshape,ioffset,tetrahedron,coord,xi,i,out)
 implicit none
 type(databoxND), intent(inout) :: db
 integer, intent(in) :: n
 integer, intent(in) :: gshape(n)
 integer, intent(in) :: ioffset(n)
 real, intent(in) :: tetrahedron(:,:,:)

 real, intent(in) :: xi(n)
 integer, intent(out) :: i(n)
 logical, intent(out) :: out

#ifdef  IMPLICIT_COORDINATE
 interface 
   function coord(indexes) result(x)
    integer, intent(in) :: indexes(:)
    real :: x(size(indexes))
   end function coord
 end interface
#elif defined(GRID_COORDINATE)
 type(grid), intent(inout) :: coord
# else
 real, intent(in) :: coord(:,:)
#endif

 integer :: p

   if (db%type.eq.db_notsplitted) then
      call split_databox_VARIANT(db,n,ioffset,coord)
   end if

   ! write(6,*) 'xi ',xi,out,db%xmin,db%xmax,db%type,size(db%db)

 !  if (xi.lt.db%xmin(1).or.xi.gt.db%xmax(1).or.yi.lt.db%xmin(2).or.yi.gt.db%xmax(2)) then
 if ( any(xi.lt.db%xmin).or.any(xi.gt.db%xmax)) then

   out = .true.
   ! write(6,*) 'xi ',xi
   ! write(6,*) 'not in ',db%xmin,db%xmax
!!$    write(6,*) ' any(xi.lt.db%xmin)', any(xi.lt.db%xmin)
!!$    write(6,*) ' any(xi.gt.db%xmax)', any(xi.gt.db%xmax)

 else
     ! test if databox has already be splitted, if not do it now

!    write(6,*) 'xi1',xi,out
    !write(6,*) 'xi2',xi,out

   if (db%type.ne.db_cell)  then

     ! search in sub boxes

     do p=1,size(db%db)
    !write(6,*) 'xi3',xi,out
       call locate_databox_VARIANT(db%db(p),n,gshape,ioffset,tetrahedron,coord,xi,i,out)       
       if (.not.out) exit
     end do
   else
     ! verify that the data point is realy in the rectangle

     i = db%imin
     !      if (all(i.lt.gshape)) then
     !write(6,*) 'test i ',i,ioffset
     out = .not.InCube_VARIANT(n,gshape,ioffset,tetrahedron,coord,xi,i)
     !write(6,*) 'i ',i,out
     !     end if
   end if
 end if

    !write(6,*) 'xi ',xi
end subroutine locate_databox_VARIANT


#endif
