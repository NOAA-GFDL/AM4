 
!-----------------------------------------------------------------------
! local data

  real*8 :: sph_out(2,nyo)
  real*8 :: lonb_in(nxi+1), sph_in(nyi+1)

  real*8  :: blon, fac, dtr, hpi, tpi, eps, critical
  integer :: i, j, k, m, n, np, npass, m2, n2, iter

  real*8 :: dlon_in  (nxi), dsph_in  (nyi),  &
            dlon_out (nxo), dsph_out (nyo), &
            facj (2,nyo), faci (2,nxo)
  integer :: jlat (2,nyo), ilon (2,nxo)

  real*8, dimension(nxi, nyi) :: area_in
  real*8, allocatable :: wt(:,:)

  integer :: id, jd, is, ie, js, je
  real*8  :: fis, fie, fjs, fje 

  real*8  :: arsum, asum
  real*8, dimension(nz) :: dwtsum, wtsum   
  character(len=64) :: mesg

  real*8 :: sum_in, sum_out
  real*8 :: sum_area_in, sum_area_out

  integer :: num_iters = 4 
  logical :: do_conservation_check = .false.
!-----------------------------------------------------------------------

   eps = epsilon(wtsum)
   critical = max(eps,min(1.0-eps,critical_frac))
   error_status = 0

   hpi = acos(0.0)
   tpi = 4.*hpi
   dtr = hpi/90.

   if (verbose > 0) then
     print *, 'nlon_in, nlat_in = ',nxi, nyi
     print *, 'nlon_out, nlat_out = ',nxo, nyo
     print *, 'data_in(1,1,1) = ', data_in(1,1,1)
     print *, 'missing_value = ', missing_value
   endif
   if (verbose > 1) then
     do_conservation_check = .true.
   endif

!  check that lon and lat are increasing (skip monotonicity check for now)
   if (blon_in(1,1) > blon_in(1,nxi)) then
     write(*,'(a)') "ERROR: Input grid must go east to west"
     error_status = 1
     return
   endif
   if (blat_in(1,1) > blat_in(1,nyi)) then
     write(*,'(a)') "ERROR: Input grid must go south to north"
     error_status = 1
     return
   endif

!-----------------------------------------------------------------------
!  --- set-up for area of input grid boxes ---

   do i = 1,nxi
       lonb_in(i) = blon_in(1,i)
   enddo
       lonb_in(nxi+1) = blon_in(2,nxi)

   do j = 1, nyi
       sph_in(j) = sin(blat_in(1,j)*dtr)
   enddo
       sph_in(nyi+1) = sin(blat_in(2,nyi)*dtr)

   do j = 1, nyi
       dsph_in(j) = abs(sph_in(j+1)-sph_in(j))
   enddo

   do i = 1,nxi
       dlon_in(i) = abs(blon_in(2,i)-blon_in(1,i))*dtr
   enddo

!  --- area of input grid boxes ---
   do j = 1,nyi
   do i = 1,nxi
       area_in(i,j) = dlon_in(i) * dsph_in(j)
   enddo
   enddo
   if (do_conservation_check) then
     sum_area_in = sum(area_in)
     print *, 'sum area_in = ',sum_area_in
   endif

!-----------------------------------------------------------------------
!  --- set-up for area of output grid boxes ---

   do n = 1, nyo
      sph_out(1,n) = sin(blat_out(1,n)*dtr)
      sph_out(2,n) = sin(blat_out(2,n)*dtr)
   enddo

   do n = 1, nyo
       dsph_out(n) = abs(sph_out(2,n)-sph_out(1,n))
   enddo

   do m = 1,nxo
       dlon_out(m) = abs(blon_out(2,m)-blon_out(1,m))*dtr
   enddo

   if (do_conservation_check) then
     sum_area_out = 0.0
     do n = 1,nyo
     do m = 1,nxo
       sum_area_out = sum_area_out + dlon_out(m) * dsph_out(n)
     enddo
     enddo
     print *, 'sum area_out = ',sum_area_out
   endif
!***********************************************************************

!------ set up latitudinal indexing ------
!------ make sure output grid goes south to north ------

 jlat = 0

 do n2 = 1, 2           ! looping on grid box edges
   do n = 1, nyo   ! looping on output latitudes
     eps = 0.0
     do iter=1,num_iters
       ! find indices from input latitudes
       do j = 1, nyi
         if ( sph_in(j)-sph_out(n2,n) <= eps .and. sph_out(n2,n)-sph_in(j+1) <= eps ) then
           jlat(n2,n) = j
           ! weight with sin(lat) to exactly conserve area-integral
           fac = (sph_out(n2,n)-sph_in(j))/(sph_in(j+1)-sph_in(j))
           if (n2 == 1) facj(n2,n) = 1.0 - fac
           if (n2 == 2) facj(n2,n) = fac
           if (verbose > 3) then
             write(*,601) n,n2,j,sph_out(n2,n),sph_in(j),sph_in(j+1),facj(n2,n)
        601  format (2x,'n',4x,'n2',5x,'j',4x,'sph',4x,'sph_in(j)',4x,'sph_in(j+1)',5x,'facj', /,(i4,2i7,4f10.5))
           endif
           exit
         endif
       enddo
       if ( jlat(n2,n) /= 0 ) exit
       ! did not find this output grid edge in the input grid
       ! increase tolerance for multiple passes
       eps  = epsilon(sph_out)*real(10**iter)
     enddo
     ! no match
     if ( jlat(n2,n) == 0 ) then
       write (mesg,710) n,sph_out(n2,n)
  710  format (': n,sph=',i3,f14.7,40x)
       write(*,'(a)') 'No latitude index found'//trim(mesg)
       error_status = 1
       return
     endif
   enddo
 enddo

!------ set up longitudinal indexing ------

 ilon = 0

 do m2 = 1, 2           ! looping on grid box edges
   do m = 1, nxo   ! looping on output longitudes
     blon = blon_out(m2,m)
     if ( blon < lonb_in(1)     ) blon = blon + 360.
     if ( blon > lonb_in(nxi+1) ) blon = blon - 360.
     eps = 0.0
     do iter=1,num_iters
       ! find indices from input longitudes
       do i = 1, nxi
         if ( (lonb_in(i)-blon) <= eps .and.  (blon-lonb_in(i+1)) <= eps ) then
           ilon(m2,m) = i
           fac = (blon-lonb_in(i))/(lonb_in(i+1)-lonb_in(i))
           if (m2 == 1) faci(m2,m) = 1.0 - fac
           if (m2 == 2) faci(m2,m) = fac
           if (verbose > 3) then
            !print *, 'm,m2,i,blon,lonb_in(i),lonb_in(i+1),eps,fac=',m,m2,i,blon,lonb_in(i),lonb_in(i+1),eps,fac
             write(*,701) m,m2,i,blon,lonb_in(i),lonb_in(i+1),faci(m2,m)
        701  format (2x,'m',4x,'m2',5x,'i',4x,'blon',4x,'lonb_in(i)',4x,'lonb_in(i+1)',5x,'faci', /,(i4,2i7,4f10.5))
           endif
           exit
         endif
       enddo
       if ( ilon(m2,m) /= 0 ) exit
       ! did not find this output grid edge in the input grid
       ! increase tolerance for multiple passes
       eps  = epsilon(blon)*real(10**iter)
     enddo
     ! no match
     if ( ilon(m2,m) == 0 ) then
       if (verbose > 2) then
          print *, 'blon_out,blon,lonb_in,eps=',blon_out(m2,m),blon,lonb_in(1),lonb_in(nxi+1),eps
       endif
       write(*,'(a)') 'No longitude index found'
       error_status = 1
       return
     endif
   enddo
 enddo

!-----------------------------------------------------------------------
! this output may be quite lengthy and is not recommended
! when using more than one processor
  if (verbose > 2) then
      write (*,801) (i,ilon(1,i),ilon(2,i),faci(1,i),faci(2,i),i=1,nxo)
      write (*,802) (j,jlat(1,j),jlat(2,j),facj(1,j),facj(2,j),j=1,nyo)
 801  format (/,2x,'i',4x,'is',5x,'ie',4x,'facis',4x,'facie', /,(i4,2i7,2f10.5))
 802  format (/,2x,'j',4x,'js',5x,'je',4x,'facjs',4x,'facje', /,(i4,2i7,2f10.5))
  endif
!-----------------------------------------------------------------------

  if (do_conservation_check) then
     sum_in = 0.0
     sum_out = 0.0
    !if (present(missing_value)) then
        do k = 1, size(data_in,3)
           sum_in  = sum_in + sum( data_in(:,:,k)*area_in, mask = data_in(:,:,k)/=missing_value )
        enddo
    !else
    !   do k = 1, size(data_in,3)
    !      sum_in  = sum_in + sum( data_in(:,:,k)*area_in )
    !   enddo
    !endif
  endif

!---- loop through output grid boxes ----

   do n = 1, nyo

    ! latitude window
    ! setup ascending latitude indices and weights
    js = jlat(1,n)
    je = jlat(2,n)
    fjs = facj(1,n)
    fje = facj(2,n)

   do m = 1, nxo

    ! longitude window
       is = ilon(1,m)
       ie = ilon(2,m)
      fis = faci(1,m)
      fie = faci(2,m)
      npass = 1
      dwtsum = 0.
       wtsum = 0.
       arsum = 0.

    ! wrap-around on input grid
    ! sum using 2 passes (pass 1: end of input grid)
      if ( ie < is ) then
           ie = nxi
          fie = 1.0
          npass = 2
      endif

      do np = 1, npass

       ! pass 2: beginning of input grid
         if ( np == 2 ) then
              is = 1
             fis = 1.0
              ie = ilon(2,m)
             fie = faci(2,m)
         endif

       ! summing data*weight and weight for single grid point
         id = ie-is+1; jd = je-js+1
         allocate(wt(id,jd))
         wt = area_in(is:ie,js:je)
         wt(1,:)=wt(1,:)*fis;  wt(id,:)=wt(id,:)*fie
         wt(:,1)=wt(:,1)*fjs;  wt(:,jd)=wt(:,jd)*fje
         asum = sum(wt)
         arsum = arsum + asum
         do k = 1, nz
            dwtsum(k) = dwtsum(k) + sum(wt*data_in(is:ie,js:je,k), mask = data_in(is:ie,js:je,k)/=missing_value)
             wtsum(k) =  wtsum(k) + sum(wt                       , mask = data_in(is:ie,js:je,k)/=missing_value)
         end do
         deallocate(wt)

      enddo !np

     !write (*,'(2i4,3e16.8)') m,n,dwtsum(1),wtsum(1),arsum
      do k = 1, size(data_in,3)
         if ( wtsum(k)/arsum >= critical ) then
            data_out(m,n,k) = dwtsum(k) / wtsum(k)
         else
           !if (present(missing_value)) then
               data_out(m,n,k) = missing_value
           !else
           !  !call error_handler ('missing value unexpectedly needed')
           !   error_message = 'missing value unexpectedly needed'
           !   return
           !endif
         endif

         ! integral check
         if (do_conservation_check) then
           !if (present(missing_value)) then
               if (data_out(m,n,k) /= missing_value) then
                  sum_out = sum_out + data_out(m,n,k)*dlon_out(m)*dsph_out(n)
               endif
           !else
           !   sum_out = sum_out + data_out(m,n,k)*dlon_out(m)*dsph_out(n)
           !endif
         endif
      enddo !k

   enddo !m
   enddo !n


  if (do_conservation_check) then
     write(*,'(a,g16.8)') 'Sum  input data = ', sum_in
     write(*,'(a,g16.8)') 'Sum output data = ', sum_out
     write(*,'(a,g16.8)') 'Avg  input data = ', sum_in/sum_area_in
     write(*,'(a,g16.8)') 'Avg output data = ', sum_out/sum_area_out
  endif

