! ============================================================================
! module numerics: a collection of useful general-purpose routines
! ============================================================================
#define __ERROR__(message) \
call my_error(mod_name,message,FATAL,thisfile,__LINE__)
#define __NOTE__(message) \
call my_error(mod_name,message,NOTE,thisfile,__LINE__)
#define __ASSERT__(x, message) \
if(.NOT.(x))call my_error(mod_name,message,FATAL,thisfile,__LINE__)

module land_numerics_mod

use fms_mod, only: error_mesg, FATAL, NOTE, mpp_pe, stdout
use mpp_mod, only: mpp_npes, mpp_get_current_pelist, mpp_send, mpp_recv, &
     mpp_sync, mpp_sync_self, EVENT_RECV, COMM_TAG_1,  COMM_TAG_2,       &
     COMM_TAG_3,  COMM_TAG_4, COMM_TAG_5,  COMM_TAG_6, COMM_TAG_7,       &
     COMM_TAG_8, COMM_TAG_9,  COMM_TAG_10, COMM_TAG_11,  COMM_TAG_12,    &
     COMM_TAG_13, COMM_TAG_14,  COMM_TAG_15, COMM_TAG_16,  COMM_TAG_17,    &
     COMM_TAG_18, COMM_TAG_19,  COMM_TAG_20


use mpp_domains_mod, only : domain2d, mpp_get_compute_domain, &
     mpp_get_global_domain, domainUG, mpp_get_ug_compute_domain

use land_data_mod, only : log_version

implicit none
private

! ==== public interfaces =====================================================
public :: bisect    ! finds a position of point in array of bounds
public :: lin_int   ! linear interpolation
public :: ludcmp, lubksb ! LU decomposition and back substitution
public :: tridiag   ! tri-diagonal system solver
public :: nearest   ! nearest point search

public :: horiz_remap_type
public :: horiz_remap_new, horiz_remap_del
public :: horiz_remap_print
public :: horiz_remap

public :: rank_descending ! rank the input array in descending order

public :: numerics_init
! ==== end of public interfaces ==============================================


!     Linear interpolation.
interface lin_int
   module procedure lin_int0
   module procedure lin_int1
   module procedure lin_int2
   module procedure lin_int1m
   module procedure lin_int2m
end interface

interface nearest
   module procedure nearest1D, nearest2D, nearestUG
end interface

logical :: module_is_initialized =.FALSE.
! ==== module constants ======================================================
character(len=*), parameter :: mod_name = 'land_numerics'
#include "../shared/version_variable.inc"
character(len=*), parameter :: thisfile = __FILE__

! ==== public type ===========================================================
! this data structure describes the horizontal remapping: that is, the operation
! of copying the data from the source points to the destination points. The source
! points are not necessarily on the same PE as destination points.
type :: horiz_remap_type
   integer :: n = 0 ! number of points that need remapping on this PE
   integer, pointer :: &
       dst_l(:)=>NULL()    ! unstructured grid indices of destination points
   integer, pointer :: &
       src_l(:)=>NULL(), & ! unstructured grid indices of source points
       src_p(:)=>NULL()    ! processor number of source points
   ! data distribution map: for each processor pair that communicate
   ! (unidirectionally), an entry in the srcPE and dstPE arrays holds their
   ! numbers. This map is the same on each of the PEs that participate in
   ! remapping.
   integer :: mapSize = 0
   integer, pointer :: &     ! for each index:
       srcPE(:) => NULL(), & ! PE that provides the data
       dstPE(:) => NULL()    ! PE that requests and then uses the data
end type horiz_remap_type


contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! ============================================================================
! Initializes the numerics module.
subroutine numerics_init()
  module_is_initialized =.TRUE.
  call log_version(version, mod_name, thisfile)
end subroutine numerics_init


! ============================================================================
!  Finds a position of point in array of bounds. Returns i, such that x is
!  between xx(i) and xx(i+1).
!  Usage:
!     value=bisect( xx, x1, periodic )
function bisect(xx, x1, periodic)
  real, intent(in)              :: xx(:)     ! array of boundaries
  real, intent(in)              :: x1        ! point to locate
  logical, intent(in), optional :: periodic  ! if present and true, the data
                                             ! domain is assumed to be periodic
  ! ---- result type
  integer bisect

  character(*),parameter :: mod_name='bisect'

   ! ---- local vars
  real    :: x              ! duplicate of input value
  integer :: low, high, mid
  integer :: n              ! size of the input array
  logical :: ascending      ! if true, the coordinates are in ascending order

  n = size(xx)
  x = x1

  ! bring the point inside bounds of the period
  if (present(periodic)) then
     if (periodic) then
        __ASSERT__(xx(n)-xx(1)/=0,"periodic bisect: period equal to zero")
        x = modulo(x-min(xx(1),xx(n)),abs(xx(n)-xx(1)))+min(xx(1),xx(n))
     endif
  endif

  ! find the coordinates
  if (x >= xx(1).and.x<=xx(n)) then
     low = 1; high = n
     ascending = xx(n) > xx(1)
     do while (high-low > 1)
        mid = (low+high)/2
        if (ascending.eqv.xx(mid) <= x) then
           low = mid
        else
           high = mid
        endif
     enddo
     bisect = low
  else
     bisect = -1
  endif

end function bisect


! ==============================================================================
!     Linearly interpolates 1-D data.
subroutine lin_int0(data, xx, x, res)
  real, intent(in) :: data(:)    ! data to interpolate
  real, intent(in) :: xx(:)      ! coordinates of the data points
  real, intent(in) :: x          ! coordinates to interpolate to
  real, intent(inout) :: res     ! result of interpolation

  character(*),parameter :: mod_name='lin_int0'

  ! ---- local vars
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(xx,x)
  i2 = i1+1
  __ASSERT__(i1>0.AND.i1<size(xx),"Coordinate is out of range")

  f2 = (x-xx(i1))/(xx(i2)-xx(i1))
  f1 = 1.0-f2

  ! update the result
  res = data(i1)*f1+data(i2)*f2

end subroutine lin_int0


! ==============================================================================
!     Linearly interpolates 1-D data.
subroutine lin_int1(data, xx, x, res)

  real, intent(in) :: data(:,:)    ! data to interpolate
  real, intent(in) :: xx(:)        ! coordinates of the data points
  real, intent(in) :: x            ! coordinates to interpolate to
  real, intent(inout) :: res(:)    ! result of interpolation

  character(*),parameter :: mod_name='lin_int0'

  ! ---- local vars
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(xx,x)
  i2 = i1+1

  __ASSERT__(i1>0.AND.i1<size(xx),"Coordinate is out of range")

  f2 = (x-xx(i1))/(xx(i2)-xx(i1))
  f1 = 1.0-f2

  ! update the result
  res = data(:,i1)*f1+data(:,i2)*f2

end subroutine lin_int1


! ==============================================================================
!     Interpolates prescribed over time.
subroutine lin_int2(data, tt, t, res)

  real, intent(in) :: data(:,:,:)  ! data to interpolate
  real, intent(in) :: tt(:)        ! time moments corresponding to data points
  real, intent(in) :: t            ! time to interpolate to
  real, intent(inout) :: res(:,:)  ! result

  character(*),parameter :: mod_name='lin_int2'
  ! ---- local vars
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(tt,t); i2 = i1+1
  __ASSERT__(i1>0.AND.i1<size(tt),"Coordinate is out of range")

  f2 = (t-tt(i1))/(tt(i2)-tt(i1))
  f1 = 1-f2

  ! update the result
  res = data(:,:,i1)*f1+data(:,:,i2)*f2

end subroutine lin_int2


!     Linearly interpolates 1-D data.
subroutine lin_int1m(data, xx, x, res, mask)
  real, intent(in) :: data(:,:)    ! data to interpolate
  real, intent(in) :: xx(:)        ! coordinates of the data points
  real, intent(in) :: x            ! coordinates to interpolate to
  real, intent(inout) :: res(:)    ! result of interpolation
  logical, intent(in) :: mask(:)   ! valid data mask

  character(*),parameter :: mod_name='lin_int1m'
  ! ---- local vars
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(xx,x)
  i2 = i1+1
  __ASSERT__(i1>0.AND.i1<size(xx),"Coordinate is out of range")

  f2 = (x-xx(i1))/(xx(i2)-xx(i1))
  f1 = 1.0-f2

  ! finally, update the result
  where (mask)
     res = data(:,i1)*f1+data(:,i2)*f2
  endwhere

end subroutine lin_int1m


! ==============================================================================
!     Interpolates prescribed over time.
subroutine lin_int2m(data, tt, t, res, mask)
  real, intent(in) :: data(:,:,:)  ! data to interpolate
  real, intent(in) :: tt(:)        ! time moments corresponding to data points
  real, intent(in) :: t            ! time to interpolate to
  real, intent(inout) :: res(:,:)  ! result
  logical, intent(in) :: mask(:,:) ! interpolation mask

  character(*),parameter :: mod_name='lin_int2m'
  ! ---- local vars
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(tt,t); i2 = i1+1
  __ASSERT__(i1>0.AND.i1<size(tt),"Coordinate is out of range")

  f2 = (t-tt(i1))/(tt(i2)-tt(i1))
  f1 = 1-f2

  ! update the result
  where (mask)
     res = data(:,:,i1)*f1+data(:,:,i2)*f2
  endwhere
end subroutine lin_int2m


! ==============================================================================
! given a matrix a(n,n) replaces it by LU decomposition of a row-wise permutation
! of itself indx(n) is an output that records the permutation. This routine is
! used in combination with lubksb to solve linear equations or invert a matrix
! example:
!    call ludcmp(a,indx)
!    call lubksb(a,indx,b1)
!    call lubksb(a,indx,b2)
subroutine ludcmp(a,indx,status)
  real,    intent(inout) :: a(:,:) ! matrix that gets replaced by its LU decomposition
  integer, intent(out)   :: indx(:) ! index of row permutations affected by partial pivoting
  integer, intent(out), optional :: status

  integer, parameter :: TINY = 1.0e-20
  integer :: n ! size of the matrix
  integer :: i,j,k,imax
  real    :: aamax,dum,sum
  real    :: vv(size(a,1)) ! implicit scaling for each row

  n = size(a,1)
  if(present(status))status = 0

  ! find largest element in each row and calculate scaling
  do i = 1,n
     aamax = 0.0
     do j = 1,n
        if(abs(a(i,j)) > aamax)aamax = abs(a(i,j))
     enddo
     if(.not.(aamax /= 0.0)) then
        if(present(status))then
           status = -1; aamax = TINY
        else
           call error_mesg('ludcmp','Matrix is singular', FATAL)
        endif
     endif
     vv(i) = 1.0/aamax
  enddo

  ! loop over the columns of Crout's method
  do j=1,n
     do i = 1,j-1
        sum = a(i,j)
        do k = 1,i-1
           sum = sum-a(i,k)*a(k,j)
        enddo
        a(i,j) = sum
     enddo
     aamax = 0.0 ! initialize the search for the largest pivot element
     do i = j,n
        sum = a(i,j)
        do k = 1,j-1
           sum = sum-a(i,k)*a(k,j)
        enddo
        a(i,j) = sum
        dum = vv(i)*abs(sum) ! figure of merit for the pivot
        if (dum >= aamax) then ! is it better than the best so far?
           imax = i
           aamax = dum
        endif
     enddo
     if (j /= imax) then ! do we need to interchange rows?
        ! Yes, do so
        do k=1,n
           dum=a(imax,k)
           a(imax,k) = a(j,k)
           a(j,k) = dum
        enddo
        vv(imax) = vv(j)
     endif
     indx(j) = imax
     ! if the pivot element is zero, then the matrix is singular (at least to the
     ! precision of the algorithm). For some applications on singular matrices, it
     ! is desirable to substitute TINY for zero
     if(a(j,j)==0.0) a(j,j) = TINY

     if (j/=n)then
        ! Finally, divide by the pivot element
        dum = 1.0/a(j,j)
        do i = j+1,n
           a(i,j) = a(i,j)*dum
        enddo
     endif
  enddo ! loop over the columns
end subroutine ludcmp


! ==============================================================================
! given a LU decomposition of matrix a(n,n), permutation vector indx, and right-
! hand side b, solves the set of linear equations A*X = B
subroutine lubksb(a,indx,b)
  real,    intent(in)    :: a(:,:)  ! LU-decomposed matrix
  integer, intent(in)    :: indx(:) ! permutation vector, as returned by the ludcmp
  real,    intent(inout) :: b(:)    ! right-hand side vector, returns solution

  integer :: i,ii,j,ll,n
  real    :: sum

  n = size(a,1)
  ii = 0
  do i = 1,n
     ll = indx(i)
     sum = b(ll)
     b(ll) = b(i)
     if (ii/=0) then
        do j = ii,i-1
           sum = sum - a(i,j)*b(j)
        enddo
     else if (sum /= 0.0) then
        ii = i
     endif
     b(i) = sum
  enddo
  do i = n, 1, -1
     sum = b(i)
     do j = i+1,n
        sum = sum-a(i,j)*b(j)
     enddo
     b(i) = sum/a(i,i)
  enddo
end subroutine lubksb


! ============================================================================
! given values of the tri-diagonal matrix coefficients, computes a solution
subroutine tridiag(a,b,c,r,u)
  real, intent(in)  :: a(:),b(:),c(:),r(:)
  real, intent(out) :: u(:)

  integer :: j
  real :: bet, gam(size(a))

  ! check that the sizes are the same
  if(size(a)/=size(b).or.size(a)/=size(c).or.size(a)/=size(r)) &
       call error_mesg('tridiag','sizes of input arrays are not equal',FATAL)
  if(size(u)<size(a)) &
       call error_mesg('tridiag','size of the result is insufficient',FATAL)
  ! check that a(1)==0 and c(N)==0
  if(a(1)/=0.or.c(size(a))/=0) &
       call error_mesg('tridiag','a(1) and c(N) must be equal to 0',FATAL)
  ! decomposition and forward substitution
  bet = b(1)
  u(1) = r(1)/bet
  do j = 2,size(a)
     gam(j) = c(j-1)/bet
     bet = b(j)-a(j)*gam(j)
     if(bet==0) &
          call error_mesg('tridiag','system is ill-defined',FATAL)
     u(j) = (r(j)-a(j)*u(j-1))/bet
  enddo
  ! backward substitution
  do j = size(a)-1,1,-1
     u(j) = u(j)-gam(j+1)*u(j+1)
  enddo
end subroutine tridiag

! ============================================================================
! finds nearest point that is not masked out in input data
! NOTE: implemented in very naive and inefficient way
subroutine nearest1D(mask, lon, lat, plon, plat, iout, jout, dist)
  logical, intent(in) :: mask(:,:)  ! mask of valid input points (.true. if valid point)
  real,    intent(in) :: lon(:)     ! longitudes of input grid central points, radian
  real,    intent(in) :: lat(:)     ! latitudes of input grid central points, radian
  real,    intent(in) :: plon, plat ! coordinates of destination point, radian
  integer, intent(out):: iout, jout ! indices of nearest valid (unmasked) point
  real, optional, intent(out):: dist ! distance to the point

  ! ---- local constants
  character(*),parameter :: mod_name='nearest1D'
  ! ---- local vars
  integer :: i,j
  real    :: r,r1

  __ASSERT__(size(mask,1)==size(lon),'sizes of "mask" and "lon" are inconsistent')
  __ASSERT__(size(mask,2)==size(lat),'sizes of "mask" and "lat" are inconsistent')

  r = HUGE(r)  ! some value larger than any possible distance

  do j = 1, size(mask,2)
  do i = 1, size(mask,1)
     if (.not.mask(i,j)) cycle
     r1 = distance(plon,plat,lon(i),lat(j))
     if ( r1 < r ) then
        iout = i
        jout = j
        r = r1
     endif
  enddo
  enddo
  if (present(dist)) dist = r
end subroutine nearest1D

! ============================================================================
! finds nearest point that is not masked out in input data
! this version works with 2D lon and lat fields
subroutine nearest2D(mask, lon, lat, plon, plat, iout, jout, dist)
  logical, intent(in) :: mask(:,:)  ! mask of valid input points (.true. if valid point)
  real,    intent(in) :: lon(:,:)   ! longitudes of input grid central points, radian
  real,    intent(in) :: lat(:,:)   ! latitudes of input grid central points, radian
  real,    intent(in) :: plon, plat ! coordinates of destination point, radian
  integer, intent(out):: iout, jout ! indices of nearest valid (unmasked) point
  real, optional, intent(out):: dist! distance to the point

  ! ---- local constants
  character(*),parameter :: mod_name='nearest2D'
  ! ---- local vars
  integer :: i,j
  real    :: r,r1

  __ASSERT__(ALL(SHAPE(mask)==SHAPE(lon)),'shapes of "mask" and "lon" are different')
  __ASSERT__(ALL(SHAPE(mask)==SHAPE(lat)),'shapes of "mask" and "lat" are different')

  r = HUGE(r)  ! some value larger than any possible distance

  do j = 1, size(mask,2)
  do i = 1, size(mask,1)
     if (.not.mask(i,j)) cycle
     r1 = distance(plon,plat,lon(i,j),lat(i,j))
     if ( r1 < r ) then
        iout = i
        jout = j
        r = r1
     endif
  enddo
  enddo
  if (present(dist)) dist=r
end subroutine nearest2D

subroutine nearestUG(mask, lon, lat, plon, plat, lout, dist)
  logical, intent(in) :: mask(:)  ! mask of valid input points (.true. if valid point)
  real,    intent(in) :: lon(:)   ! longitudes of input grid central points, radian
  real,    intent(in) :: lat(:)   ! latitudes of input grid central points, radian
  real,    intent(in) :: plon, plat ! coordinates of destination point, radian
  integer, intent(out):: lout       ! indices of nearest valid (unmasked) point
  real, optional, intent(out):: dist! distance to the point

  ! ---- local constants
  character(*),parameter :: mod_name='nearestUG'
  ! ---- local vars
  integer :: l
  real    :: r,r1

  __ASSERT__(ALL(SHAPE(mask)==SHAPE(lon)),'shapes of "mask" and "lon" are different')
  __ASSERT__(ALL(SHAPE(mask)==SHAPE(lat)),'shapes of "mask" and "lat" are different')

  r = HUGE(r)  ! some value larger than any possible distance

  do l = 1, size(mask(:))
     if (.not.mask(l)) cycle
     r1 = distance(plon,plat,lon(l),lat(l))
     if ( r1 < r ) then
        lout = l
        r = r1
     endif
  enddo
  if (present(dist)) dist=r
end subroutine nearestUG


! ============================================================================
! private functions that calculates the distance between two points given their
! coordinates
function distance(lon1, lat1, lon2, lat2) ; real distance
  ! calculates distance between points on unit square
  real, intent(in) :: lon1,lat1,lon2,lat2

  real :: x1,y1,z1, x2,y2,z2
  real :: dlon
  dlon = (lon2-lon1)

  z1 = sin(lat1) ;  z2 = sin(lat2)
  y1 = 0.0       ;  y2 = cos(lat2)*sin(dlon)
  x1 = cos(lat1) ;  x2 = cos(lat2)*cos(dlon)

  ! distance = acos(x1*x2 + z1*z2)
  distance = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
end function distance


! ============================================================================
! deallocate memory associated with the horiz_remap_type data structure
subroutine horiz_remap_del(map)
   type(horiz_remap_type), intent(inout) :: map
#define __DEALLOC__(x)\
if (associated(x)) then; deallocate(x); x=>NULL(); endif
   __DEALLOC__(map%dst_l)
   __DEALLOC__(map%src_l)
   __DEALLOC__(map%src_p)
   map%n=0

   __DEALLOC__(map%srcPE)
   __DEALLOC__(map%dstPE)
   map%mapSize=0
#undef __DEALLOC__
end subroutine

! ============================================================================
! prints remapping information
subroutine horiz_remap_print(map, prefix)
   type(horiz_remap_type), intent(in) :: map
   character(*), intent(in) :: prefix

   integer :: k

   do k = 1, map%n
      write(*,100) prefix,&
         map%src_l(k),map%src_p(k),&
         map%dst_l(k),mpp_pe()
   enddo
100 format(a,'(L:',i4.4,' PE:',i4.4,') -> (L:',i4.4,' PE:',i4.4,')')
end subroutine

! ============================================================================
! given the local mask of the points that need filling, the local mask of
! the valid points, local arrays of coordinates, and the domain, returns the
! remapping information that can be used later to fill the data
subroutine horiz_remap_new(invalid, valid, lon, lat, domain, pes, map)
  logical, intent(in) :: invalid(:) ! mask of points to be filled
  logical, intent(in) :: valid  (:) ! mask of valid input points
  real,    intent(in) :: lon(:)   ! longitudes of input grid central points, radian
  real,    intent(in) :: lat(:)   ! latitudes of input grid central points, radian
  type(domainUG), intent(in) :: domain ! our domain
  integer, intent(in) :: pes(:)     ! list of PEs
  type(horiz_remap_type), intent(out) :: map ! remapping information


  ! --- local constants
  character(*), parameter :: mod_name='horiz_remap_new'
  ! --- local vars
  integer :: ntot ! total number of missing points across all PEs
  integer :: ls,le ! boundaries of our compute domain
  integer :: npes ! total number of PEs
  integer :: root_pe ! root PE for this operation
  integer, allocatable :: np(:) ! number of missing points per processor
  integer :: p ! processor iterator
  real   , allocatable :: glon(:), glat(:) ! global arrays of missing point coordinates
  integer, allocatable :: from_pe(:) ! number of PE the missing points belong to
  real   , allocatable :: dist(:) ! distance to the missing points
  integer, allocatable :: ll(:)   ! indices of the nearest points
  integer, allocatable :: lbuf(:)  !send/receive buffers for indices
  real   , allocatable :: dbuf(:) ! send/receive buffer for distances
  integer :: k,m,n1,l
  integer :: k0

  ! get the size of our domain
  call mpp_get_ug_compute_domain(domain, ls, le)
  ! check the input array shapes
  if(size(invalid(:))/=le-ls+1) then
    call my_error(mod_name,'shape of input array "'//'invalid'//'" must be the same as shape of compute domain',FATAL)
  endif
  if(size(valid(:))/=le-ls+1) then
    call my_error(mod_name,'shape of input array "'//'valid'//'" must be the same as shape of compute domain',FATAL)
  endif
  if(size(lon(:))/=le-ls+1) then
    call my_error(mod_name,'shape of input array "'//'lon'//'" must be the same as shape of compute domain',FATAL)
  endif
  if(size(lat(:))/=le-ls+1) then
    call my_error(mod_name,'shape of input array "'//'lat'//'" must be the same as shape of compute domain',FATAL)
  endif

  ! get the number of missing points for this PE
  map%n = count(invalid)

  ! get the number of PEs that communicate
  npes = size(pes)
  ! and the number of the root PE
  root_pe = pes(1)

  ! [x] compute the global number of missing points and assemble the
  ! array of point numbers per PE on root PE
  ! no need to send data to oneself (rab)
  if (mpp_pe()/=root_pe) then
     call mpp_send(map%n,root_pe,tag=COMM_TAG_1)
     call mpp_recv(ntot,root_pe,tag=COMM_TAG_2)
  else
     allocate(np(npes))
     np(1)=map%n
     do p = 2,npes
        call mpp_recv(np(p),pes(p),tag=COMM_TAG_1)
     enddo
     ntot = sum(np)
     do p = 2, npes
        call mpp_send(ntot,pes(p),tag=COMM_TAG_2)
     enddo
  endif

  call mpp_sync_self()
  ! we don't need to do anything if there are no missing points anywhere
  if (ntot==0) return

  ! [x] allocate global buffers
  allocate(glon(ntot),glat(ntot),from_pe(ntot),dist(ntot),ll(ntot))

  ! allocate buffers for missing point indices and processors
  allocate(map%dst_l(map%n))
  allocate(map%src_l(map%n), map%src_p(map%n))
  ! and fill the coordinates of missing points for this PE
  k = 1
  do l=1,size(invalid(:))
     if (invalid(l)) then
        glon(k)      = lon(l); glat(k)      = lat(l)
        map%dst_l(k) = l+ls-1 
        k = k+1
     endif
  enddo

  ! [x] send the array of point coordinates to root PE and get the global
  ! arrays of point coordinates in return
  if (mpp_pe()/=root_pe) then
     ! non-root PE sends its data to the root
     call mpp_send(map%n,root_pe, tag=COMM_TAG_3)
     if (map%n>0) then
        call mpp_send(glon(1), plen=map%n, to_pe=root_pe, tag=COMM_TAG_4)
        call mpp_send(glat(1), plen=map%n, to_pe=root_pe, tag=COMM_TAG_5)
     endif
     ! and then receives the global data in response
     call mpp_recv(glon(1),glen=ntot,from_pe=root_pe, tag=COMM_TAG_6)
     call mpp_recv(glat(1),glen=ntot,from_pe=root_pe, tag=COMM_TAG_7)
  else
     ! root PE receives data from all PEs and assembles global coordinate arrays
     ! in the order of PEs in the list, except that it puts its own data first.
     from_pe(1:map%n) = root_pe
     k=map%n+1
     do p = 1,npes
        if (pes(p)==root_pe) cycle
        call mpp_recv(n1,pes(p), tag=COMM_TAG_3)
        if (n1>0) then
           call mpp_recv(glon(k),glen=n1,from_pe=pes(p), tag=COMM_TAG_4)
           call mpp_recv(glat(k),glen=n1,from_pe=pes(p), tag=COMM_TAG_5)
           from_pe(k:k+n1-1)=pes(p)
        endif
        k = k+n1
     enddo
     ! then it distributes the resulting array among PEs
     do p = 1,npes
        if (pes(p)==root_pe) cycle
        call mpp_send(glon(1),plen=ntot,to_pe=pes(p), tag=COMM_TAG_6)
        call mpp_send(glat(1),plen=ntot,to_pe=pes(p), tag=COMM_TAG_7)
     enddo
  endif
  call mpp_sync_self()

  ! [x] find the nearest points in the domain
  do k = 1, ntot
     call nearest(valid,lon,lat,glon(k),glat(k),ll(k),dist=dist(k))
  enddo
  ! convert local domain indices to global
  ll(:) = ll(:)+ls-1

  ! [5] send the data to root PE and let it calculate the points corresponding to
  ! the global minimum distance
  if (mpp_pe()/=root_pe) then
     ! non-root PE just sends the data
     call mpp_send(ll(1)  ,plen=ntot,to_pe=root_pe, tag=COMM_TAG_8)
     call mpp_send(dist(1),plen=ntot,to_pe=root_pe, tag=COMM_TAG_9)
     ! and receives the updated data in response
     if(map%n>0) then
        ! receive the nearest point locations and PEs
        call mpp_recv(map%src_l(1),glen=map%n,from_pe=root_pe, tag=COMM_TAG_10)
        call mpp_recv(map%src_p(1),glen=map%n,from_pe=root_pe, tag=COMM_TAG_11)
     endif
     ! receive communication map
     call mpp_recv(map%mapSize,glen=1,from_pe=root_pe, tag=COMM_TAG_12)
     if (map%mapSize>0) then
        allocate (map%srcPE(map%mapSize),map%dstPE(map%mapSize))
        call mpp_recv(map%srcPE(1),glen=map%mapSize,from_pe=root_pe, tag=COMM_TAG_13)
        call mpp_recv(map%dstPE(1),glen=map%mapSize,from_pe=root_pe, tag=COMM_TAG_14)
     endif
  else
     ! root PE does the bulk of processing: it assembles all the data
     ! and sends the relevant parts back to the processors that need them

     ! receive data about domain-specific nearest points from PEs and select
     ! the globally nearest point among them
     allocate(lbuf(ntot),dbuf(ntot))
     ! note that arrays ll, and dist are initially filled with the
     ! nearest points information for the root PE own domain
     from_pe(:) = root_pe
     do p = 1,npes
        if (pes(p)==root_pe) cycle
        call mpp_recv(lbuf(1),glen=ntot,from_pe=pes(p), tag=COMM_TAG_8)
        call mpp_recv(dbuf(1),glen=ntot,from_pe=pes(p), tag=COMM_TAG_9)
        do k = 1,ntot
           ! to avoid dependence on the order of operations, give preference
           ! to the lowest leftmost point among the equidistant points
           if (dbuf(k)<dist(k).or.(&
               dbuf(k)==dist(k).and.lbuf(k)<ll(k))) then
              ll(k)=lbuf(k); dist(k)=dbuf(k); from_pe(k)=pes(p)
           endif
        enddo
     enddo

     ! release buffers
     deallocate (lbuf,dbuf)

     ! create a communication map: arrays srcPE and dstPE listing all pairs that
     ! communicate
     allocate(map%srcPE(ntot),map%dstPE(ntot)) ! allocate max possible number
     k0=1; m=1
     do p = 1, npes
        do k = sum(np(1:p-1))+1, sum(np(1:p))
           if (from_pe(k) == pes(p)) cycle ! skip communications to itself
           if (ANY(map%srcPE(k0:m-1)==from_pe(k))) cycle ! skip src->dst pair that already exists
           ! add current pair to the communication map
           map%srcPE(m)=from_pe(k); map%dstPE(m)=pes(p); m=m+1
        enddo
        k0=m
     enddo

     ! actual number of elements in comm. map is m-1
     map%mapSize=m-1

     ! simply assign the results for the root PE
     if (map%n>0) then
        map%src_l(:) = ll(1:map%n)
        map%src_p(:) = from_pe(1:map%n)
     endif
     ! distribute the results among processors
     k = map%n+1
     do p = 1,npes
        if (pes(p)==root_pe) cycle
        if (np(p)>0) then
           ! send nearest point location
           call mpp_send(ll(k),plen=np(p),to_pe=pes(p), tag=COMM_TAG_10)
           call mpp_send(from_pe(k),plen=np(p),to_pe=pes(p), tag=COMM_TAG_11)
        endif
        ! broadcast comm. map
        call mpp_send(map%mapSize,plen=1,to_pe=pes(p), tag=COMM_TAG_12)
        if (map%mapSize>0) then
           call mpp_send(map%srcPE(1),plen=map%mapSize,to_pe=pes(p), tag=COMM_TAG_13)
           call mpp_send(map%dstPE(1),plen=map%mapSize,to_pe=pes(p), tag=COMM_TAG_14)
        endif
        call mpp_sync_self()
        k = k+np(p)
     enddo

  endif

  call mpp_sync_self()

  deallocate(glon,glat,from_pe,dist,ll)

  ! note that many communications in this routine can be sped up if the data
  ! are combined.
  ! For example instead of sending ll and from_pe one can encode them
  ! in a single integer array [ a(i*3-2)=ii(i), a(i*3-1)=jj(i), a(i*3)=from_pe(i) ]
  ! and send that array.

end subroutine

! ============================================================================
subroutine horiz_remap(map,domain,d)
  type(horiz_remap_type), intent(in)    :: map
  type(domainUG)        , intent(in)    :: domain
  real                  , intent(inout) :: d(:,:) ! field to fill

  character(*),parameter :: mod_name='horiz_remap'
  ! ---- local vars
  integer :: i,j,k,n,ls,le
  integer :: is,ie,js,je
  integer, allocatable :: ll(:)
  real   , allocatable :: buf(:,:)
  logical :: ltmp

  call mpp_get_ug_compute_domain(domain, ls, le)
  ltmp = size(d,1)==le-ls+1
  __ASSERT__(ltmp,'shape of data must be the same as shape of compute domain')

  ! handle the local points
  do i = 1, map%n
     if (map%src_p(i)==mpp_pe()) then
       d(map%dst_l(i)-ls+1,:) = d(map%src_l(i)-ls+1,:)
     endif
  enddo

  ! exchange information with other processors
  do k = 1, map%mapSize
     if (map%srcPE(k)==mpp_pe()) then
        ! get the size of the data from the other PE
        call mpp_recv(n,map%dstPE(k), tag=COMM_TAG_15)
        allocate(ll(n),buf(n,size(d,2)))
        ! get the indices
        call mpp_recv(ll(1),glen=n,from_pe=map%dstPE(k), tag=COMM_TAG_16)
        ! fill the buffer
        do i = 1,n
           if(ll(i)<ls.or.ll(i)>le) call error_mesg('distr_fill','requested index l outside of domain', FATAL)
           buf(i,:) = d(ll(i)-ls+1,:)
        enddo
        ! send the buffer
        call mpp_send(buf(1,1),plen=size(buf),to_pe=map%dstPE(k), tag=COMM_TAG_17)
        call mpp_sync_self()
        deallocate (ll,buf)
     else if (map%dstPE(k)==mpp_pe()) then
        ! send data request
        n = count(map%src_p(:)==map%srcPE(k))
        ! alloacate and fill arrays of requested indices ii and jj
        allocate(ll(n),buf(n,size(d,2)))
        j = 1
        do i = 1, map%n
           if (map%src_p(i)==map%srcPE(k)) then
              ll(j) = map%src_l(i); j = j+1
           endif
        enddo
        ! send the data request
        call mpp_send(n,map%srcPE(k), tag=COMM_TAG_15)
        call mpp_send(ll(1),plen=n,to_pe=map%srcPE(k), tag=COMM_TAG_16)

        ! get the response
        call mpp_recv(buf(1,1),glen=size(buf),from_pe=map%srcPE(k), tag=COMM_TAG_17)
        ! fill the data
        j = 1
        do i = 1,map%n
           if (map%src_p(i)==map%srcPE(k)) then
              d(map%dst_l(i)-ls+1,:) = buf(j,:) ; j = j+1
           endif
        enddo
        call mpp_sync_self()
        deallocate (ll,buf)
     endif
  enddo

end subroutine

! ======================================================================
! ranks array x in descending order: on return, idx() contains indices
! of elements of array x in descending order of x values
subroutine rank_descending(x,idx)
   real,    intent(in)  :: x(:)
   integer, intent(out) :: idx(:)

   integer :: i,n
   integer, allocatable :: t(:)

   n = size(x)
   do i = 1,n
      idx(i) = i
   enddo

   allocate(t((n+1)/2))
   call mergerank(x,idx,n,t)
   deallocate(t)
end subroutine


! =====================================================================
! based on:
! http://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
subroutine merge(x,a,na,b,nb,c,nc)
   integer, intent(in) :: na,nb,nc ! Normal usage: NA+NB = NC
   real, intent(in)       :: x(*)
   integer, intent(in)    :: a(na)    ! B overlays C(NA+1:NC)
   integer, intent(in)    :: b(nb)
   integer, intent(inout) :: c(nc)

   integer :: i,j,k

   i = 1; j = 1; k = 1;
   do while(i <= na .and. j <= nb)
      if (x(a(i)) >= x(b(j))) then
         c(k) = a(i) ; i = i+1
      else
         c(k) = b(j) ; j = j+1
      endif
      k = k + 1
   enddo
   do while (i <= na)
      c(k) = a(i) ; i = i + 1 ; k = k + 1
   enddo
end subroutine merge

recursive subroutine mergerank(x,a,n,t)
  integer, intent(in) :: n
  real,    intent(in) :: x(*)
  integer, dimension(n), intent(inout) :: a
  integer, dimension((n+1)/2), intent (out) :: t

  integer :: na,nb
  integer :: v

  if (n < 2) return
  if (n == 2) then
     if ( x(a(1)) < x(a(2)) ) then
        v = a(1) ; a(1) = a(2) ; a(2) = v
     endif
     return
  endif
  na=(n+1)/2
  nb=n-na

  call mergerank(x,a,na,t)
  call mergerank(x,a(na+1),nb,t)

  if (x(a(na)) < x(a(na+1))) then
     t(1:na)=a(1:na)
     call merge(x,t,na,a(na+1),nb,a,n)
  endif
end subroutine mergerank

! ==============================================================================
! Reports error, including file name and line.
subroutine my_error(mod_name, message, mode, file, line)

  character(len=*), intent(in) :: mod_name
  character(len=*), intent(in) :: message
  integer,          intent(in) :: mode
  character(len=*), intent(in), optional :: file
  integer,          intent(in), optional :: line

  ! ---- local vars
  character(len=512) :: mesg
  if(present(file)) then ! assume that file and line are either both present or not
     write(mesg,'("File ",a," Line ",i4.4," :: ",a)')&
       file, line, trim(message)
  else
    mesg = trim(message)
  endif
  call error_mesg(mod_name, mesg, mode)
end subroutine


end module land_numerics_mod
