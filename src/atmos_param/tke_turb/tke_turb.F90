module tke_turb_mod

!=======================================================================
!  Turbulence Kinetic Energy parameterization.
!
!  Originally based on GFDL Mellor-Yamada level 2.5 turbulence
!  closure scheme.
!
!  Modified by Chris Golaz
!=======================================================================

 use           mpp_mod, only: input_nml_file

 use           fms_mod, only: file_exist, open_namelist_file,       &
                              error_mesg, FATAL, close_file, note,  &
                              check_nml_error, mpp_pe, mpp_root_pe, &
                              write_version_number, stdlog, stdout, &
                              mpp_chksum

 use     constants_mod, only: rdgas, rvgas, kappa, grav, vonkarm,   &
                              cp_air, hlv, hls, tfreeze

 use monin_obukhov_mod, only: mo_diff

 use sat_vapor_pres_mod,only: compute_qs

 use  diag_manager_mod, only: register_diag_field, send_data

 use  time_manager_mod, only: time_type
 
!---------------------------------------------------------------------
  implicit none
  private
  public :: tke_turb, tke_turb_init, tke_turb_end

  character(len=128) :: version = ''
  character(len=128) :: tagname = ''
  logical            :: module_is_initialized = .false.
 
!---------------------------------------------------------------------
! --- Constants
!---------------------------------------------------------------------

 real, parameter :: p00    = 1000.0e2
 real, parameter :: p00inv = 1./p00
 real, parameter :: d622   = rdgas/rvgas
 real, parameter :: d378   = 1.-d622
 real, parameter :: d608   = d378/d622

 real :: ckm1,  ckm2,  ckm3, ckm4, ckm5, ckm6, ckm7, ckm8
 real :: ckh1,  ckh2,  ckh3, ckh4
 real :: cvfq1, cvfq2, bcq

 real, parameter :: aa1     =  0.92
 real, parameter :: aa2     =  0.74
 real, parameter :: bb1     = 16.6  !bqx (orig: 16.0)
 real, parameter :: bb2     = 10.1  !bqx (orig: 10.0)
 real, parameter :: ccc     =  0.08
 real, parameter :: cc1     =  0.27
 real, parameter :: t00     =  2.7248e2 
 real, parameter :: small   =  1.0e-10

 ! Constants for new scheme

 real, parameter :: cp_air_inv = 1.0/cp_air

 real, parameter :: A1 = 0.92
 real, parameter :: A2 = 0.74
 real, parameter :: B1 = 16.6
 real, parameter :: B2 = 10.1
 real, parameter :: C1 = 0.08

 real, parameter :: alpha1 = A1*(1.0-3.0*C1-6.0*A1/B1)
 real, parameter :: alpha2 = -3.0*A1*A2*((B2-3.0*A2)*(1.0-6.0*A1/B1)-3.0*C1*(B2+6.0*A1))
 real, parameter :: alpha3 = -3.0*A2*(6.0*A1+B2)
 real, parameter :: alpha4 = -9.0*A1*A2
 real, parameter :: alpha5 = A2*(1.0-6.0*A1/B1)

 real, parameter :: be = B1 / (2.0*sqrt(2.0))
 real, parameter :: ke = 0.5421

!---------------------------------------------------------------------
! --- Namelist
!---------------------------------------------------------------------

 real    :: tkemax           =  5.0
 real    :: tkemin           =  0.0
 integer :: tke_option       =  0
 integer :: buoync_option    =  2
 integer :: pbl_depth_option =  0
 real    :: tkecrit          =  0.05
 real    :: parcel_buoy      =  1.0
 real    :: akmax            =  1.0e4
 real    :: akmin_land       =  5.0
 real    :: akmin_sea        =  0.0
 integer :: nk_lim           =  2
! real    :: el0max           =  1.0e6
 real    :: el0max           =  1.0e6
 real    :: el0min           =  0.0
 real    :: alpha_land       =  0.10
 real    :: alpha_sea        =  0.10
 integer :: el_sl            =  1        ! bqx start level to calculate el0
 integer :: leng_option      =  1        ! bqx add another option 2
 integer :: el_sl_option     =  1        ! bqx starting level based on level or pressure above surface
 real    :: pcrit_from_sfc   =  2000.e2  ! bqx critical level (pressure) above surface

 namelist / tke_turb_nml /                            &
         tke_option, buoync_option,                   &
         pbl_depth_option, tkecrit, parcel_buoy,      &
         tkemax,   tkemin,                            &
         akmax,    akmin_land, akmin_sea, nk_lim,     &
         el0max,   el0min,  alpha_land,  alpha_sea,   &
         leng_option, el_sl_option, el_sl,            &!bqx+
         pcrit_from_sfc                                !bqx+ 

!---------------------------------------------------------------------
!--- Diagnostic fields       
!---------------------------------------------------------------------

character(len=10) :: mod_name = 'tke_turb'
real              :: missing_value = -999.99
integer           :: id_h, id_pblh_tke, id_pblh_parcel

!--> h1g, 2015-08-11
integer           :: id_tke_avg_pbl, id_tke_turb, id_dthetav, id_dtheta, id_dqt
!<-- h1g, 2015-08-11

!---------------------------------------------------------------------

 contains

!#######################################################################

 subroutine tke_turb( is, ie, js, je, time, delt, fracland,          &
                      phalf, pfull, zhalf, zfull,                    &
                      tt, qv, qa, ql, qi, um, vm, z0, ustar, bstar,  &
                      tr_tke,                                        & 
                      el0, el, akm, akh, h, tke_avg )

  integer,         intent(in)           :: is,ie,js,je
  type(time_type), intent(in)           :: time
  real,    intent(in)                   :: delt 
  real,    intent(in), dimension(:,:)   :: fracland
  real,    intent(in), dimension(:,:,:) :: phalf, pfull, zhalf, zfull
  real,    intent(in), dimension(:,:,:) :: tt, qv, qa, ql, qi, um, vm
  real,    intent(in), dimension(:,:)   :: z0, ustar, bstar
  real,    intent(inout), dimension(:,:,:) :: tr_tke
  real, intent(out), dimension(:,:)   :: el0
  real, intent(out), dimension(:,:,:) :: el, akm, akh
  real, intent(out), dimension(:,:)   :: h

!---> h1g, 2015-08-11
  real, intent(out), optional, dimension(:,:) :: tke_avg  !averaged TKE within PBL
!<--- h1g, 2015-08-11

  if( present(tke_avg) ) tke_avg = 0.0  ! h1g, 2015-08-11

  if (tke_option == 0) then 

    call tke_turb_legacy( is, ie, js, je, time, delt, fracland,        &
                          phalf, pfull, zhalf, zfull,                  &
                          tt, qv, ql, qi, um, vm, z0, ustar, bstar,    &
                          tr_tke,                                      & 
                          el0, el, akm, akh, h )

  elseif (tke_option == 1) then

   if( present(tke_avg) ) then   ! h1g, 2015-08-11
    call tke_turb_dev( is, ie, js, je, time, delt, fracland,        &
                       phalf, pfull, zhalf, zfull,                  &
                       tt, qv, qa, ql, qi, um, vm,                  & 
                       z0, ustar, bstar,                            &
                       tr_tke,                                      & 
                       el0, el, akm, akh, h, tke_avg=tke_avg )
   else
    call tke_turb_dev( is, ie, js, je, time, delt, fracland,        &
                       phalf, pfull, zhalf, zfull,                  &
                       tt, qv, qa, ql, qi, um, vm,                  & 
                       z0, ustar, bstar,                            &
                       tr_tke,                                      & 
                       el0, el, akm, akh, h )
   endif   ! h1g, 2015-08-11

  end if

 end subroutine tke_turb

!#######################################################################

 subroutine tke_turb_legacy( is, ie, js, je, time, delt, fracland,     &
                             phalf, pfull, zhalf, zfull,               &
                             tt, qv, ql, qi, um, vm, z0, ustar, bstar, &
                             tr_tke,                                   & 
                             el0, el, akm, akh, h )

!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!   is,ie,js,je - i,j indices marking the slab of model working on
!   time        - variable needed for netcdf diagnostics
!   delt        - time step in seconds
!   fracland    - fractional amount of land beneath a grid box
!   phalf       - pressure at half levels
!   pfull       - pressure at full levels
!   zhalf       - height at half levels
!   zfull       - height at full levels
!   tt          - temperature
!   qv          - water vapor
!   ql          - liquid water
!   qi          - ice water
!   um, vm      - wind components
!   z0          - roughness length
!   ustar       - friction velocity (m/sec)
!   bstar       - buoyancy scale (m/sec**2)
!---------------------------------------------------------------------

  integer,         intent(in)           :: is,ie,js,je
  type(time_type), intent(in)           :: time
  real,    intent(in)                   :: delt 
  real,    intent(in), dimension(:,:)   :: fracland
  real,    intent(in), dimension(:,:,:) :: phalf, pfull, zhalf, zfull
  real,    intent(in), dimension(:,:,:) :: tt, qv, ql, qi, um, vm
  real,    intent(in), dimension(:,:)   :: z0, ustar, bstar

!---------------------------------------------------------------------
! Arguments (Intent inout)
!   tr_tke   -  Tracer that stores tke
!---------------------------------------------------------------------

  real,    intent(inout), dimension(:,:,:) :: tr_tke

!---------------------------------------------------------------------
! Arguments (Intent out)
!   el0      -  characteristic length scale
!   el       -  master length scale
!   akm      -  mixing coefficient for momentum
!   akh      -  mixing coefficient for heat and moisture
!   h        -  diagnosed depth of planetary boundary layer (m)
!---------------------------------------------------------------------

  real, intent(out), dimension(:,:)   :: el0
  real, intent(out), dimension(:,:,:) :: el, akm, akh
  real, intent(out), dimension(:,:)   :: h

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------

  logical used
  integer outunit
  integer :: ix, jx, kx, i, j, k
  integer :: kxp, kxm, klim
  real    :: cvfqdt, dvfqdt

  real, dimension(size(um,1),size(um,2)) :: zsfc, x1, x2, akmin,  &
        pblh_tke, pblh_parcel

  real, dimension(size(um,1),size(um,2),size(um,3)-1) ::     &
        dsdzh, shear, buoync, qm2,  qm3, qm4, el2,           &
        aaa,   bbb,   ccc,    ddd,                           &
        xxm1,  xxm2,  xxm3,   xxm4, xxm5

  real, dimension(size(um,1),size(um,2),size(um,3)) ::       &
        thetav, tmp, dsdz, qm, xx1, xx2

  real, dimension(size(um,1),size(um,2),size(um,3)+1) :: tke

!====================================================================

! --- Check to see if tke_turb has been initialized
  if( .not. module_is_initialized ) call error_mesg( ' tke_turb',     &
                                 ' tke_turb_init has not been called',&
                                   FATAL )

! --- Set dimensions etc
  ix  = size( um, 1 )
  jx  = size( um, 2 )
  kx  = size( um, 3 )
  kxp = kx + 1
  kxm = kx - 1

!====================================================================
! --- Copy input tke from tracer array
!====================================================================

  tke(:,:,1) = tkemin
  tke(:,:,2:kx) = tr_tke(:,:,1:kxm)
  tke(:,:,kxp) = bcq * ustar * ustar

!====================================================================
! --- Compute virtual potential temperature
!====================================================================

  tmp(:,:,:) = (pfull(:,:,:)*p00inv)**(-kappa)
  thetav(:,:,:) = tt(:,:,:)*(qv(:,:,:)*d608+1.0)*tmp

!====================================================================
! --- Surface height     
!====================================================================

  zsfc(:,:) = zhalf(:,:,kxp)

!====================================================================
! --- d( )/dz operators: at full levels & at half levels          
!====================================================================

   dsdz(:,:,1:kx)  = 1.0 / ( zhalf(:,:,2:kxp) - zhalf(:,:,1:kx) )
  dsdzh(:,:,1:kxm) = 1.0 / ( zfull(:,:,2:kx)  - zfull(:,:,1:kxm) )

!====================================================================
! --- Wind shear                 
!====================================================================

  xxm1(:,:,1:kxm) = dsdzh(:,:,1:kxm)*( um(:,:,2:kx) - um(:,:,1:kxm) )
  xxm2(:,:,1:kxm) = dsdzh(:,:,1:kxm)*( vm(:,:,2:kx) - vm(:,:,1:kxm) )

  shear = xxm1 * xxm1 + xxm2 * xxm2

!====================================================================
! --- Buoyancy                 
!====================================================================

  xxm1(:,:,1:kxm) = thetav(:,:,2:kx) - thetav(:,:,1:kxm) 
  xxm2(:,:,1:kxm) = 0.5*( thetav(:,:,2:kx) + thetav(:,:,1:kxm) )

  buoync = grav * dsdzh * xxm1 / xxm2

!====================================================================
! --- Some tke stuff
!====================================================================

  do k=1,kx
  do j=1,jx
  do i=1,ix
    xx1(i,j,k) = 2*tke(i,j,k+1)
    if(xx1(i,j,k) > 0.0) then
      qm(i,j,k) = sqrt(xx1(i,j,k))
    else
      qm(i,j,k) = 0.0
    endif
  enddo
  enddo
  enddo

  qm2(:,:,1:kxm)  = xx1(:,:,1:kxm) 
  qm3(:,:,1:kxm)  =  qm(:,:,1:kxm) * qm2(:,:,1:kxm) 
  qm4(:,:,1:kxm)  = qm2(:,:,1:kxm) * qm2(:,:,1:kxm) 

!====================================================================
! --- Characteristic length scale                         
!====================================================================

  xx1(:,:,1:kxm) = qm(:,:,1:kxm)*( pfull(:,:,2:kx) - pfull(:,:,1:kxm) )

  do k = 1, kxm
     xx2(:,:,k) = xx1(:,:,k)  * ( zhalf(:,:,k+1) - zsfc(:,:) )
  end do

  xx1(:,:,kx) =  qm(:,:,kx) * ( phalf(:,:,kxp) - pfull(:,:,kx) )
  xx2(:,:,kx) = xx1(:,:,kx) * z0(:,:)

  x1 = sum( xx1, 3 )
  x2 = sum( xx2, 3 )

!---- should never be equal to zero ----
  if (count(x1 <= 0.0) > 0) call error_mesg( ' tke_turb',  &
                             'divid by zero, x1 <= 0.0', FATAL)
  el0 = x2 / x1
  el0 = el0 * (alpha_land*fracland + alpha_sea*(1.-fracland))

  el0 = min( el0, el0max )
  el0 = max( el0, el0min )

!====================================================================
! --- Master length scale 
!====================================================================

  do k = 1, kxm
     xx1(:,:,k)  = vonkarm * ( zhalf(:,:,k+1) - zsfc(:,:) )
  end do

  x1(:,:) = vonkarm * z0(:,:) 
  xx1(:,:,kx) = x1(:,:)

  do k = 1,kx
    el(:,:,k+1) = xx1(:,:,k) / ( 1.0 + xx1(:,:,k) / el0(:,:) )
  end do
    el(:,:,1)   = el0(:,:)

  el2(:,:,1:kxm) = el(:,:,2:kx) * el(:,:,2:kx)

!====================================================================
! --- Mixing coefficients                     
!====================================================================

  xxm3(:,:,1:kxm) = el2(:,:,1:kxm)*buoync(:,:,1:kxm)
  xxm4(:,:,1:kxm) = el2(:,:,1:kxm)* shear(:,:,1:kxm)
  xxm5(:,:,1:kxm) =  el(:,:,2:kx )*   qm3(:,:,1:kxm)

!-------------------------------------------------------------------
! --- Momentum 
!-------------------------------------------------------------------
 
  xxm1 = xxm5*( ckm1*qm2 + ckm2*xxm3 )
  xxm2 = qm4 + ckm5*qm2*xxm4 + xxm3*( ckm6*xxm4 + ckm7*qm2 + ckm8*xxm3 )

  xxm2 = max( xxm2, 0.2*qm4 )
  xxm2 = max( xxm2, small  )

  akm(:,:,1)    = 0.0
  akm(:,:,2:kx) = xxm1(:,:,1:kxm) / xxm2(:,:,1:kxm)

  akm = max( akm, 0.0 )

!-------------------------------------------------------------------
! --- Heat and moisture 
!-------------------------------------------------------------------

  xxm1(:,:,1:kxm) = ckh1*xxm5(:,:,1:kxm) - ckh2*xxm4(:,:,1:kxm)*akm(:,:,2:kx)
  xxm2(:,:,1:kxm) = qm2(:,:,1:kxm) + ckh3*xxm3(:,:,1:kxm)

  xxm1 = max( xxm1, ckh4*xxm5 )
  xxm2 = max( xxm2, 0.4*qm2   )
  xxm2 = max( xxm2, small     )

  akh(:,:,1)    = 0.0
  akh(:,:,2:kx) = xxm1(:,:,1:kxm) / xxm2(:,:,1:kxm)

!-------------------------------------------------------------------
! --- Bounds 
!-------------------------------------------------------------------

! --- Upper bound
  akm = min( akm, akmax )
  akh = min( akh, akmax )

! --- Lower bound near surface

  akmin = akmin_land*fracland + akmin_sea*(1.-fracland)
  klim = kx - nk_lim + 1
  do  k = klim,kx
    akm(:,:,k) = max( akm(:,:,k), akmin(:,:) )
    akh(:,:,k) = max( akh(:,:,k), akmin(:,:) )
  end do

!====================================================================
! --- Prognosticate turbulence kinetic energy
!====================================================================

  cvfqdt = cvfq1 * delt
  dvfqdt = cvfq2 * delt * 2.0

!-------------------------------------------------------------------
! --- Part of linearized energy disiipation term 
!-------------------------------------------------------------------

  xxm1(:,:,1:kxm) = dvfqdt * qm(:,:,1:kxm) / el(:,:,2:kx)

!-------------------------------------------------------------------
! --- Part of linearized vertical diffusion term
!-------------------------------------------------------------------

  xx1(:,:,1:kx) = el(:,:,2:kxp) * qm(:,:,1:kx)

  xx2(:,:,1)    = 0.5*  xx1(:,:,1)
  xx2(:,:,2:kx) = 0.5*( xx1(:,:,2:kx) + xx1(:,:,1:kxm) )

  xx1 = xx2 * dsdz

!-------------------------------------------------------------------
! --- Implicit time differencing for vertical diffusion 
! --- and energy dissipation term 
!-------------------------------------------------------------------
 
  do k=1,kxm
  do j=1,jx
  do i=1,ix
    aaa(i,j,k) = -cvfqdt * xx1(i,j,k+1) * dsdzh(i,j,k)
    ccc(i,j,k) = -cvfqdt * xx1(i,j,k  ) * dsdzh(i,j,k)
    bbb(i,j,k) =     1.0 - aaa(i,j,k  ) -   ccc(i,j,k) 
    bbb(i,j,k) =           bbb(i,j,k  ) +  xxm1(i,j,k)
    ddd(i,j,k) =           tke(i,j,k+1)
  enddo
  enddo
  enddo

! correction for vertical diffusion of tke surface boundary condition

  do j = 1,jx
  do i = 1,ix
    ddd(i,j,kxm) = ddd(i,j,kxm) - aaa(i,j,kxm) * tke(i,j,kxp)
  enddo
  enddo

! solve tridiagonal system

  call tri_invert( xxm1, ddd, aaa, bbb, ccc ) 

!-------------------------------------------------------------------
! --- Shear and buoyancy terms
!-------------------------------------------------------------------

  xxm2(:,:,1:kxm) =  delt*( akm(:,:,2:kx)* shear(:,:,1:kxm)    &
                          - akh(:,:,2:kx)*buoync(:,:,1:kxm) )

!-------------------------------------------------------------------
! --- Update turbulence kinetic energy
!-------------------------------------------------------------------

  do j=1,jx
  do i=1,ix
    tke(i,j,1) = 0.0
    do k=2,kx
      tke(i,j,k) = xxm1(i,j,k-1) + xxm2(i,j,k-1)
    enddo
  enddo
  enddo

!====================================================================
! --- Bound turbulence kinetic energy
!====================================================================

  tke(:,:,:) = min( tke(:,:,:), tkemax )
  tke(:,:,:) = max( tke(:,:,:), tkemin )

!====================================================================
! --- Compute PBL depth
!====================================================================

  if ( pbl_depth_option == 0 .or. id_pblh_tke > 0) then
    call tke_pbl_depth(zsfc,zhalf,tke,bstar,pblh_tke)
  end if
  if ( pbl_depth_option == 1 .or. id_pblh_parcel > 0) then
    call parcel_pbl_depth(zsfc,zfull,tt,qv,ql,qi,ustar,bstar,pblh_parcel)
  end if

  if ( pbl_depth_option == 0 ) then
    h = pblh_tke
  else if ( pbl_depth_option == 1 ) then
    h = pblh_parcel
  end if

!====================================================================
! --- Copy output tke back to tracer array
!====================================================================

  tr_tke(:,:,:) = tke(:,:,2:kxp)

!====================================================================
! --- Diagnostics
!====================================================================

  if ( id_h > 0 ) then
    used = send_data ( id_h, h, time, is, js )
  end if
  if ( id_pblh_tke > 0 ) then
    used = send_data ( id_pblh_tke, pblh_tke, time, is, js )
  end if
  if ( id_pblh_parcel > 0 ) then
    used = send_data ( id_pblh_parcel, pblh_parcel, time, is, js )
  end if

  end subroutine tke_turb_legacy

!#######################################################################

 subroutine tke_turb_dev( is, ie, js, je, time, delt, fracland,     &
                          phalf, pfull, zhalf, zfull,               &
                          T, qv, qa, ql, qi, um, vm,                &
                          z0, ustar, bstar,                         &
                          tr_tke,                                   & 
                          el0, el, akm, akh, h, tke_avg )

!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!   is,ie,js,je - i,j indices marking the slab of model working on
!   time        - variable needed for netcdf diagnostics
!   delt        - time step in seconds
!   fracland    - fractional amount of land beneath a grid box
!   phalf       - pressure at half levels
!   pfull       - pressure at full levels
!   zhalf       - height at half levels
!   zfull       - height at full levels
!   T           - temperature
!   qv          - water vapor
!   qa          - cloud amount
!   ql          - liquid water
!   qi          - ice water
!   um, vm      - wind components
!   z0          - roughness length
!   ustar       - friction velocity (m/sec)
!   bstar       - buoyancy scale (m/sec**2)
!---------------------------------------------------------------------

  integer,         intent(in)           :: is,ie,js,je
  type(time_type), intent(in)           :: time
  real,    intent(in)                   :: delt 
  real,    intent(in), dimension(:,:)   :: fracland
  real,    intent(in), dimension(:,:,:) :: phalf, pfull, zhalf, zfull
  real,    intent(in), dimension(:,:,:) :: T, qv, qa, ql, qi, um, vm
  real,    intent(in), dimension(:,:)   :: z0, ustar, bstar

!---------------------------------------------------------------------
! Arguments (Intent inout)
!   tr_tke   -  Tracer that stores tke
!---------------------------------------------------------------------

  real,    intent(inout), dimension(:,:,:) :: tr_tke

!---------------------------------------------------------------------
! Arguments (Intent out)
!   el0      -  characteristic length scale
!   el       -  master length scale
!   akm      -  mixing coefficient for momentum
!   akh      -  mixing coefficient for heat and moisture
!   h        -  diagnosed depth of planetary boundary layer (m)
!---------------------------------------------------------------------

  real, intent(out), dimension(:,:)   :: el0
  real, intent(out), dimension(:,:,:) :: el, akm, akh
  real, intent(out), dimension(:,:)   :: h

!---> h1g, 2015-08-11
  real, intent(out), optional, dimension(:,:) :: tke_avg
  integer :: k_pbl
!<--- h1g, 2015-08-11

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------

  logical used
  integer outunit
  integer :: ix, jx, kx, i, j, k
  integer :: kxp, kxm, klim
  real    :: kedt, bedt

  real, dimension(size(um,1),size(um,2)) :: zsfc, x1, x2, akmin,  &
        pblh_tke, pblh_parcel, hleff, Tl, qsl, qslT, gamma,       &
        cqt, cthl, Aclr, Bclr, Acld, Bcld

  real, dimension(size(um,1),size(um,2),size(um,3)-1) ::          &
        dsdzh, shear, buoync

  real, dimension(size(um,1),size(um,2),size(um,3)-1), target ::  &
        xxm1, xxm2, xxm3, xxm4, xxm5

  real, dimension(:,:,:), pointer ::                              &
        aaa, bbb, ccc, ddd, Gh, Sm, Sh

  real, dimension(size(um,1),size(um,2),size(um,3)) ::            &
        dsdz, sl, qt, sv, bhl, bqt

  real, dimension(size(um,1),size(um,2),size(um,3)), target ::    &
        xx1, xx2, xx3

  real, dimension(:,:,:), pointer ::                              &
        tmp

  real, dimension(size(um,1),size(um,2),size(um,3)+1) :: tke, sqrte

  real, dimension(size(um,1),size(um,2))              :: tke_avg_tmp, qa_tmp !h1g, 2015-08-11
  real, dimension(size(um,1),size(um,2))              :: dthetav, dtheta, dqt
  real, dimension(size(um,1),size(um,2),size(um,3))   :: Aclr_3D, Bclr_3D, Acld_3D, Bcld_3D
  real, dimension(size(um,1),size(um,2),size(um,3))   :: bfact  ! bqx+
  real :: cleng !bqx

  cleng = 3. !bqx

!====================================================================

! --- Check to see if tke_turb has been initialized
  if( .not. module_is_initialized ) call error_mesg( ' tke_turb',     &
                                 ' tke_turb_init has not been called',&
                                   FATAL )

! --- Set dimensions etc
  ix  = size( um, 1 )
  jx  = size( um, 2 )
  kx  = size( um, 3 )
  kxp = kx + 1
  kxm = kx - 1

!====================================================================
! --- Copy input tke from tracer array
!====================================================================

  tke(:,:,1) = tkemin
  tke(:,:,2:kx) = tr_tke(:,:,1:kxm)
  tke(:,:,kxp) = bcq * ustar * ustar

!====================================================================
! --- Compute thermodynamic variables
!====================================================================

  do k=1,kx

    ! Effective latent heat
    hleff(:,:) = (min(1.,max(0.,0.05*(T(:,:,k)-tfreeze+20.)))*hlv +            &
                  min(1.,max(0.,0.05*(tfreeze -T(:,:,k)   )))*hls)

    ! Liquid water temperature
    Tl(:,:) = T(:,:,k) - cp_air_inv*hleff(:,:)*(ql(:,:,k)+qi(:,:,k))

    ! Liquid water static energy (sl/cp_air)
    sl(:,:,k) = Tl(:,:) + cp_air_inv*grav*zfull(:,:,k)

    ! Total water
    qt(:,:,k) = qv(:,:,k) + ql(:,:,k) + qi(:,:,k)

    ! Virtual static energy (sv/cp_air)
    sv(:,:,k)                                                                  &
    = sl(:,:,k)                                                                &
      + d608*T(:,:,k)*qt(:,:,k)                                                &
      + (cp_air_inv*hleff(:,:) - T(:,:,k)*(1.0+d608))*(ql(:,:,k) + qi(:,:,k))

    if (buoync_option == 1) then

      ! bhl, bqt
      bhl(:,:,k) = 1.0
      bqt(:,:,k) = d608 * T(:,:,k)

    elseif (buoync_option == 2) then

      ! bhl, bqt
      bhl(:,:,k) = 1.0
      bqt(:,:,k) = d608 * T(:,:,k)

      ! Saturation specific humidity at Tl (qsl) and its derivative (qslT)
      call compute_qs(Tl(:,:), pfull(:,:,k), qsl, dqsdT=qslT)

      ! gamma
      gamma = cp_air_inv * hleff * qslT

      ! cqt, cthl
      cqt = 1.0 / (1.0 + gamma)
      cthl = cqt * gamma * cp_air / hleff * (pfull(:,:,k)*p00inv)**(kappa)
    
      ! Add moist contribution to bhl, bqt
      x1 = qa(:,:,k) * (cp_air_inv*hleff(:,:) - T(:,:,k)*(1.0+d608))
      bhl(:,:,k) = bhl(:,:,k) - cthl * x1
      bqt(:,:,k) = bqt(:,:,k) + cqt  * x1

    elseif (buoync_option == 3) then

      ! Moist buoyancy formulation from  Cuijpers and Duynkerke (1993, JAS)

      ! Saturation specific humidity at T (qsl) and its derivative (qslT)
      call compute_qs(T(:,:,k), pfull(:,:,k), qsl, dqsdT=qslT)

      ! gamma
      gamma = cp_air_inv * hleff * qslT

      ! Clear sky coefficients
      Aclr = 1.0 + d608*qt(:,:,k)
      Bclr = d608*T(:,:,k)

      ! Cloud sky coefficients
      Acld = ( 1.0 - qt(:,:,k) + (1.0+d608)*(qsl + d622*T(:,:,k)*qslT) )  &
             / ( 1.0 + gamma )
      Bcld = cp_air_inv*hleff * Acld/T(:,:,k) - 1.0
      
      ! All sky
      bhl(:,:,k) = (1.0-qa(:,:,k))*Aclr + qa(:,:,k)*Acld
      bqt(:,:,k) = (1.0-qa(:,:,k))*Bclr + qa(:,:,k)*Bcld

!---> h1g, 2015-08-27
    elseif (buoync_option == 4) then

      ! Moist buoyancy formulation from  Cuijpers and Duynkerke (1993, JAS)
      ! bugfix for "buoync_option == 3"

      ! Saturation specific humidity at T (qsl) and its derivative (qslT)
      call compute_qs(T(:,:,k), pfull(:,:,k), qsl, dqsdT=qslT)

      ! gamma
      gamma = cp_air_inv * hleff * qslT

      ! Clear sky coefficients
      Aclr = 1.0 + d608*qt(:,:,k)
      Bclr = d608*T(:,:,k)

      ! Cloud sky coefficients
      Acld = ( 1.0 - qt(:,:,k) + (1.0+d608)*(qsl + T(:,:,k)*qslT) )  &
             / ( 1.0 + gamma )
      Bcld = cp_air_inv*hleff * Acld - T(:,:,k)
      
      ! All sky
      bhl(:,:,k) = (1.0-qa(:,:,k))*Aclr + qa(:,:,k)*Acld
      bqt(:,:,k) = (1.0-qa(:,:,k))*Bclr + qa(:,:,k)*Bcld
!<--- h1g, 2015-08-27
    elseif (buoync_option == 5) then
      call compute_qs(T(:,:,k), pfull(:,:,k), qsl, dqsdT=qslT)
      ! gamma
      gamma = cp_air_inv * hleff * qslT

      ! Clear sky coefficients
      Aclr_3D(:,:,k) = 1.0 + d608*qt(:,:,k)
      Bclr_3D(:,:,k) = d608*T(:,:,k)

      ! Cloud sky coefficients
      Acld_3D(:,:,k) = ( 1.0 - qt(:,:,k) + (1.0+d608)*(qsl + T(:,:,k)*qslT) )  &
             / ( 1.0 + gamma )
      Bcld_3D(:,:,k) = cp_air_inv*hleff * Acld_3D(:,:,k) - T(:,:,k)
    elseif (buoync_option == 6) then !bqx+ option 6

      bfact(:,:,k) = grav / (T(:,:,k) * ( 1.0 + d608 * qv(:,:,k) - ql(:,:,k) - qi(:,:,k)) )
      call compute_qs(T(:,:,k), pfull(:,:,k), qsl, dqsdT=qslT)
      ! gamma
      gamma = cp_air_inv * hleff * qslT

      ! Clear sky coefficients
      Aclr_3D(:,:,k) = (1.0 + d608*qt(:,:,k)) * bfact(:,:,k) * cp_air_inv   
      Bclr_3D(:,:,k) = d608*T(:,:,k) * bfact(:,:,k)          

      ! Cloud sky coefficients

      Acld_3D(:,:,k) = ( ( 1.0 + ( 1.0+d608 )*gamma*cp_air*T(:,:,k)/hleff) / (1.0+gamma ))*bfact(:,:,k)*cp_air_inv 
      Bcld_3D(:,:,k) = hleff * Acld_3D(:,:,k) - bfact(:,:,k) * T(:,:,k) 

    end if

  end do

!---> h1g, 2015-08-31
  if (buoync_option == 5 .or. buoync_option==6) then
    do k=1,kxm
      Aclr(:,:) = 0.5 * ( Aclr_3D(:,:,k) + Aclr_3D(:,:,k+1) )
      Bclr(:,:) = 0.5 * ( Bclr_3D(:,:,k) + Bclr_3D(:,:,k+1) )

      Acld(:,:) = 0.5 * ( Acld_3D(:,:,k) + Acld_3D(:,:,k+1) )
      Bcld(:,:) = 0.5 * ( Bcld_3D(:,:,k) + Bcld_3D(:,:,k+1) )

      ! All sky at level k+1/2
      qa_tmp(:,:) = min( qa(:,:,k), 0.5*(qa(:,:,k)+qa(:,:,k+1)) )
      bhl(:,:,k) = (1.0-qa_tmp(:,:))*Aclr + qa_tmp(:,:)*Acld
      bqt(:,:,k) = (1.0-qa_tmp(:,:))*Bclr + qa_tmp(:,:)*Bcld
    end do
  endif
!<--- h1g, 2015-08-31

!====================================================================
! --- Surface height     
!====================================================================

  zsfc(:,:) = zhalf(:,:,kxp)

!====================================================================
! --- d( )/dz operators: at full levels & at half levels          
!====================================================================

   dsdz(:,:,1:kx)  = 1.0 / ( zhalf(:,:,1:kx)  - zhalf(:,:,2:kxp) )
  dsdzh(:,:,1:kxm) = 1.0 / ( zfull(:,:,1:kxm) - zfull(:,:,2:kx) )

!====================================================================
! --- Wind shear                 
!====================================================================

  xxm1(:,:,1:kxm) = dsdzh(:,:,1:kxm)*( um(:,:,1:kxm) - um(:,:,2:kx) )
  xxm2(:,:,1:kxm) = dsdzh(:,:,1:kxm)*( vm(:,:,1:kxm) - vm(:,:,2:kx) )

  shear = xxm1 * xxm1 + xxm2 * xxm2

!====================================================================
! --- Buoyancy                 
!====================================================================

  if (buoync_option == 0) then

    xxm1(:,:,1:kxm) = sv(:,:,1:kxm) - sv(:,:,2:kx) 
    xxm2(:,:,1:kxm) = 0.5*( sv(:,:,1:kxm) + sv(:,:,2:kx) )
   
    buoync = grav * dsdzh * xxm1 / xxm2

  elseif(buoync_option == 5) then
    xxm1(:,:,1:kxm) = sl(:,:,1:kxm) - sl(:,:,2:kx) 
    xxm2(:,:,1:kxm) = bhl(:,:,1:kxm) 
    xxm3(:,:,1:kxm) = qt(:,:,1:kxm) - qt(:,:,2:kx) 
    xxm4(:,:,1:kxm) = bqt(:,:,1:kxm)  
    xxm5(:,:,1:kxm) = 0.5*( sv(:,:,1:kxm) + sv(:,:,2:kx) )

    buoync = grav * dsdzh / xxm5 * ( xxm1*xxm2 + xxm3*xxm4 )
!bqx+ option 6
  elseif(buoync_option == 6) then
    xxm1(:,:,1:kxm) = sl(:,:,1:kxm) - sl(:,:,2:kx)
    xxm2(:,:,1:kxm) = bhl(:,:,1:kxm)
    xxm3(:,:,1:kxm) = qt(:,:,1:kxm) - qt(:,:,2:kx)
    xxm4(:,:,1:kxm) = bqt(:,:,1:kxm)

    buoync =  dsdzh * ( cp_air * xxm1*xxm2 + xxm3*xxm4 )
! bqx + end
  else
    xxm1(:,:,1:kxm) = sl(:,:,1:kxm) - sl(:,:,2:kx) 
    xxm2(:,:,1:kxm) = 0.5*( bhl(:,:,1:kxm) + bhl(:,:,2:kx) )
    xxm3(:,:,1:kxm) = qt(:,:,1:kxm) - qt(:,:,2:kx) 
    xxm4(:,:,1:kxm) = 0.5*( bqt(:,:,1:kxm) + bqt(:,:,2:kx) )
    xxm5(:,:,1:kxm) = 0.5*( sv(:,:,1:kxm) + sv(:,:,2:kx) )

    buoync = grav * dsdzh / xxm5 * ( xxm1*xxm2 + xxm3*xxm4 )

  end if

!====================================================================
! --- Some tke stuff
!====================================================================

  do k=1,kxp
  do j=1,jx
  do i=1,ix
    if(tke(i,j,k) > 0.0) then
      sqrte(i,j,k) = sqrt(tke(i,j,k))
    else
      sqrte(i,j,k) = 0.0
    endif
  enddo
  enddo
  enddo

!====================================================================
! --- Characteristic length scale                         
!====================================================================
  xx1(:,:,1:kxm) = sqrt (2.0) * sqrte(:,:,2:kx)*( pfull(:,:,2:kx) - pfull(:,:,1:kxm) )  ! bqx updated but no impact on el0

  do k = 1, kxm
     xx2(:,:,k) = xx1(:,:,k)  * ( zhalf(:,:,k+1) - zsfc(:,:) )
  end do

  xx1(:,:,kx) =  sqrt (2.0) * sqrte(:,:,kxp) * ( phalf(:,:,kxp) - pfull(:,:,kx) ) ! bqx updated but no impact on el0

  xx2(:,:,kx) = xx1(:,:,kx) * z0(:,:)

! bqx +
 if (el_sl_option == 1) then
  x1 = sum( xx1(:,:,el_sl:kx), 3 )
  x2 = sum( xx2(:,:,el_sl:kx), 3 )
 elseif (el_sl_option == 2) then  ! bqx+ option 2
  do j=1,jx
  do i=1,ix
   do k=1,kx
    do while ( (phalf (i,j,kxp) - phalf(i,j,k)) >= pcrit_from_sfc) 
     el_sl = k
    enddo
   enddo 
    x1(i,j) = sum( xx1(i,j,el_sl:kx))
    x2(i,j) = sum( xx2(i,j,el_sl:kx))
  enddo
  enddo
 else
    call error_mesg( 'tke_turb',  &
                     'el_sl_option only has 2 options now', FATAL)
 endif 
! bqx end

!---- should never be equal to zero ----
  if (count(x1 <= 0.0) > 0) call error_mesg( ' tke_turb',  &
                             'divid by zero, x1 <= 0.0', FATAL)
   el0 = x2 / x1
   el0 = el0 * (alpha_land*fracland + alpha_sea*(1.-fracland))
!  write(0,*) 'el0=',el0
!  write(0,*) 'pblh_tke=',pblh_tke
  el0 = min( el0, el0max )
  el0 = max( el0, el0min )

!====================================================================
! --- Master length scale 
!====================================================================

  do k = 1, kxm
     xx1(:,:,k)  = vonkarm * ( zhalf(:,:,k+1) - zsfc(:,:) )
  end do

  x1(:,:) = vonkarm * z0(:,:) 
  xx1(:,:,kx) = x1(:,:)

 !bqx add option 2 for leng_option 
  do k = 1,kx
   if (leng_option == 1 ) then
    el(:,:,k+1) = xx1(:,:,k) / ( 1.0 + xx1(:,:,k) / el0(:,:) )
   elseif(leng_option == 2) then 
    el(:,:,k+1) = ( xx1(:,:,k)**(-cleng) + el0(:,:)**(-cleng) )**(-1./cleng)
   else
    call error_mesg( 'tke_turb',  &
                     'leng_option only has 2 options now', FATAL)
   endif
  end do
 !bqx+ end

    el(:,:,1)   = el0(:,:)

!====================================================================
! --- Mixing coefficients                     
!====================================================================

  Gh => xxm1; Sm => xxm2; Sh => xxm3

  Gh(:,:,1:kxm) = - buoync(:,:,1:kxm) * (el(:,:,2:kx)*el(:,:,2:kx)) &
                    / max( 2.0 * tke(:,:,2:kx), small)   ! bqx updated
 
  Gh(:,:,1:kxm) = min( Gh(:,:,1:kxm), 0.0233 )

  Sm = (alpha1+alpha2*Gh) / ((1.0+alpha3*Gh)*(1.0+alpha4*Gh))
  Sh = alpha5 / (1.0+alpha3*Gh)

  akm(:,:,1)    = 0.0

  akm(:,:,2:kx) = Sm * el(:,:,2:kx) * sqrte(:,:,2:kx) * sqrt(2.0)  ! bqx updated

  akh(:,:,1)    = 0.0

  akh(:,:,2:kx) = Sh * el(:,:,2:kx) * sqrte(:,:,2:kx) * sqrt (2.0) ! bqx updated 

  nullify(Gh,Sm,Sh) 

!-------------------------------------------------------------------
! --- Bounds 
!-------------------------------------------------------------------

! --- Upper bound
  akm = max( min( akm, akmax ), 0.0)
  akh = max( min( akh, akmax ), 0.0)

! --- Lower bound near surface

  akmin = akmin_land*fracland + akmin_sea*(1.-fracland)
  klim = kx - nk_lim + 1
  do  k = klim,kx
    akm(:,:,k) = max( akm(:,:,k), akmin(:,:) )
    akh(:,:,k) = max( akh(:,:,k), akmin(:,:) )
  end do

!====================================================================
! --- Prognosticate turbulence kinetic energy
!====================================================================

  aaa => xxm2; bbb => xxm3; ccc => xxm4; ddd => xxm5

  kedt = delt * ke
  bedt = delt / be

!-------------------------------------------------------------------
! --- Part of linearized energy dissipation term 
!-------------------------------------------------------------------
 
  xxm1(:,:,1:kxm) = bedt * sqrte(:,:,2:kx) / el(:,:,2:kx)

!-------------------------------------------------------------------
! --- Part of linearized vertical diffusion term
!-------------------------------------------------------------------

  xx1(:,:,1:kx) = el(:,:,2:kxp) * sqrte(:,:,2:kxp)

  xx2(:,:,1)    = 0.5*  xx1(:,:,1)
  xx2(:,:,2:kx) = 0.5*( xx1(:,:,2:kx) + xx1(:,:,1:kxm) )

  xx1 = xx2 * dsdz

!-------------------------------------------------------------------
! --- Implicit time differencing for vertical diffusion 
! --- and energy dissipation term 
!-------------------------------------------------------------------
 
  do k=1,kxm
  do j=1,jx
  do i=1,ix
    aaa(i,j,k) = -kedt * xx1(i,j,k+1) * dsdzh(i,j,k) 
    ccc(i,j,k) = -kedt * xx1(i,j,k  ) * dsdzh(i,j,k)
    bbb(i,j,k) =   1.0 - aaa(i,j,k  ) -   ccc(i,j,k) 
    bbb(i,j,k) =         bbb(i,j,k  ) +  xxm1(i,j,k)
    ddd(i,j,k) =         tke(i,j,k+1) 
  enddo
  enddo
  enddo

! correction for vertical diffusion of tke surface boundary condition

  do j = 1,jx
  do i = 1,ix
    ddd(i,j,kxm) = ddd(i,j,kxm) - aaa(i,j,kxm) * tke(i,j,kxp)
  enddo
  enddo

! solve tridiagonal system

  call tri_invert( xxm1, ddd, aaa, bbb, ccc ) 

  nullify(aaa,bbb,ccc,ddd)

!-------------------------------------------------------------------
! --- Shear and buoyancy terms
!-------------------------------------------------------------------

  xxm2(:,:,1:kxm) =  delt*( akm(:,:,2:kx)* shear(:,:,1:kxm)    &
                          - akh(:,:,2:kx)*buoync(:,:,1:kxm) )

!-------------------------------------------------------------------
! --- Update turbulence kinetic energy
!-------------------------------------------------------------------

  do j=1,jx
  do i=1,ix
    tke(i,j,1) = 0.0
    do k=2,kx
      tke(i,j,k) = xxm1(i,j,k-1) + xxm2(i,j,k-1)
    enddo
  enddo
  enddo

!====================================================================
! --- Bound turbulence kinetic energy
!====================================================================

  tke(:,:,:) = min( tke(:,:,:), tkemax )
  tke(:,:,:) = max( tke(:,:,:), tkemin )

!====================================================================
! --- Compute PBL depth
!====================================================================

  if ( pbl_depth_option == 0 .or. id_pblh_tke > 0) then
    call tke_pbl_depth(zsfc,zhalf,tke,bstar,pblh_tke)
  end if
  if ( pbl_depth_option == 1 .or. id_pblh_parcel > 0) then
    call parcel_pbl_depth(zsfc,zfull,T,qv,ql,qi,ustar,bstar,pblh_parcel)
  end if

  if ( pbl_depth_option == 0 ) then
    h = pblh_tke
  else if ( pbl_depth_option == 1 ) then
    h = pblh_parcel
  end if

!====================================================================
! --- Copy output tke back to tracer array
!====================================================================

  tr_tke(:,:,:) = tke(:,:,2:kxp)

!--> h1g, 2015-08-11
! calculate averaged TKE within PBL, and moisture and temperature jump across inversion
  dqt     = missing_value
  dthetav = missing_value
  dtheta  = missing_value

  do j=1,jx
    do i=1,ix
      tke_avg_tmp(i,j) = tke(i,j,kxp) * ( phalf(i,j,kxp)-pfull(i,j,kx) )
      k_pbl        = kx
      do k= kx, 2, -1
        if( h(i,j) > zfull(i,j,k)-zsfc(i,j) ) then 
          tke_avg_tmp(i,j) = tke_avg_tmp(i,j) + tke(i,j,k)*( pfull(i,j,k)-pfull(i,j,k-1) )      
          k_pbl        = k-1
        endif
      enddo  
      tke_avg_tmp(i,j) = tke_avg_tmp(i,j)/( phalf(i,j,kxp)-pfull(i,j,k_pbl))

      if( k_pbl<kx .and. k_pbl>1 ) then
        dqt(i,j)       = (qt(i,j,k_pbl-1)-qt(i,j,k_pbl+1))/(zfull(i,j,k_pbl-1)-zfull(i,j,k_pbl+1))
        dqt(i,j)       = dqt(i,j) * 1.e3  ! convert from kg/kg to g/kg

        dthetav(i,j)   = T(i,j,k_pbl-1)*(1.0+d608*qv(i,j,k_pbl-1))*(pfull(i,j,k_pbl-1)*p00inv)**(-kappa) &
                        -T(i,j,k_pbl+1)*(1.0+d608*qv(i,j,k_pbl+1))*(pfull(i,j,k_pbl+1)*p00inv)**(-kappa)
        dthetav(i,j)   = dthetav(i,j)/(zfull(i,j,k_pbl-1)-zfull(i,j,k_pbl+1))

        dtheta(i,j)   = T(i,j,k_pbl-1)*(pfull(i,j,k_pbl-1)*p00inv)**(-kappa) &
                       -T(i,j,k_pbl+1)*(pfull(i,j,k_pbl+1)*p00inv)**(-kappa)
        dtheta(i,j)   = dtheta(i,j)/(zfull(i,j,k_pbl-1)-zfull(i,j,k_pbl+1))
      endif 
    enddo
  enddo
  tke_avg_tmp(:,:) = min(tke_avg_tmp(:,:), tkemax)
  tke_avg_tmp(:,:) = max(tke_avg_tmp(:,:), tkemin)
  if( present(tke_avg))  tke_avg = tke_avg_tmp
!<-- h1g, 2015-08-11

!====================================================================
! --- Diagnostics
!====================================================================

  if ( id_h > 0 ) then
    used = send_data ( id_h, h, time, is, js )
  end if
  if ( id_pblh_tke > 0 ) then
    used = send_data ( id_pblh_tke, pblh_tke, time, is, js )
  end if
  if ( id_pblh_parcel > 0 ) then
    used = send_data ( id_pblh_parcel, pblh_parcel, time, is, js )
  end if

!--> h1g, 2015-08-11
  if ( id_tke_avg_pbl > 0 ) then
    used = send_data ( id_tke_avg_pbl, tke_avg_tmp, time, is, js )
  end if

  if ( id_tke_turb > 0 ) then
    used = send_data ( id_tke_turb, tr_tke, time, is, js, 1 )
  end if

  if ( id_dqt > 0 ) then
    used = send_data ( id_dqt, dqt, time, is, js )
  end if

  if ( id_dthetav > 0 ) then
    used = send_data ( id_dthetav, dthetav, time, is, js )
  end if

  if ( id_dtheta > 0 ) then
    used = send_data ( id_dtheta, dtheta, time, is, js )
  end if

!<-- h1g, 2015-08-11

  end subroutine tke_turb_dev


!#######################################################################

  subroutine tke_turb_init(lonb, latb, axes, time, idim, jdim, kdim)

!=======================================================================
! --- Initialize tke module
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!---------------------------------------------------------------------
!   latb, lonb       - latitudes and longitudes at grid box corners
!   axes, time       - variables needed for netcdf diagnostics
!   idim, jdim, kdim - size of the first 3 dimensions 

 integer,              intent(in) :: idim, jdim, kdim, axes(4)
 type(time_type),      intent(in) :: time
 real, dimension(:,:), intent(in) :: lonb, latb

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------
 integer             :: unit, io, ierr, logunit

!=====================================================================

!---------------------------------------------------------------------
! --- Read namelist
!---------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=tke_turb_nml, iostat=io)
   ierr = check_nml_error(io,'tke_turb_nml')
#else   
  if( file_exist( 'input.nml' ) ) then
! -------------------------------------
   unit = open_namelist_file( )
   ierr = 1
   do while( ierr .ne. 0 )
   READ ( unit,  nml = tke_turb_nml, iostat = io, end = 10 ) 
   ierr = check_nml_error (io, 'tke_turb_nml')
   end do
10 continue
   call close_file( unit )
! -------------------------------------
  end if
#endif

!---------------------------------------------------------------------
! --- Output version
!---------------------------------------------------------------------

  if ( mpp_pe() == mpp_root_pe() ) then
       call write_version_number(version, tagname)
       logunit = stdlog()
       WRITE( logunit, nml = tke_turb_nml ) 
  endif

!---------------------------------------------------------------------
! --- Initialize constants
!---------------------------------------------------------------------

     ckm1 = ( 1.0 - 3.0*ccc )*aa1
     ckm3 =  3.0 * aa1*aa2*    ( bb2 - 3.0*aa2 )
     ckm4 =  9.0 * aa1*aa2*ccc*( bb2 + 4.0*aa1 )
     ckm5 =  6.0 * aa1*aa1
     ckm6 = 18.0 * aa1*aa1*aa2*( bb2 - 3.0*aa2 )
     ckm7 =  3.0 * aa2*        ( bb2 + 7.0*aa1 )
     ckm8 = 27.0 * aa1*aa2*aa2*( bb2 + 4.0*aa1 )
     ckm2 =  ckm3 - ckm4
     ckh1 =  aa2
     ckh2 =  6.0 * aa1*aa2
     ckh3 =  3.0 * aa2*( bb2 + 4.0*aa1 )
     ckh4 =  2.0e-6 * aa2
    cvfq1 = 5.0 * cc1 / 3.0
    cvfq2 = 1.0 / bb1
      bcq = 0.5 * ( bb1**(2.0/3.0) )

!---------------------------------------------------------------------
!--- Register diagnostic fields       
!---------------------------------------------------------------------

  id_h = register_diag_field (mod_name, 'zpbl', axes(1:2),     &
       time, 'diagnosed PBL depth', 'meters',                  &
       missing_value=missing_value )

!---> h1g, 2015-08-11
  id_tke_avg_pbl = register_diag_field (mod_name, 'tke_avg_pbl', axes(1:2),     &
       time, 'diagnosed averaged tke within depth', 'm2/s2',                  &
       missing_value=missing_value )

  id_tke_turb = register_diag_field (mod_name, 'tke_turb', axes(1:3),     &
       time, 'tke diagnosed from tr_tke', 'm2/s2',                  &
       missing_value=missing_value )

  id_dthetav = register_diag_field (mod_name, 'dthetav', axes(1:2),     &
       time, 'virtual potential temperature jump across inversion','K/m', &
       missing_value=missing_value )

  id_dtheta = register_diag_field (mod_name, 'dtheta', axes(1:2),     &
       time, 'potential temperature jump across inversion','K/m', &
       missing_value=missing_value )

  id_dqt = register_diag_field (mod_name, 'dqt', axes(1:2),     &
       time, 'total water content jump across inversion','g/kg/m', &
       missing_value=missing_value )
!<--- h1g, 2015-08-11

!---------------------------------------------------------------------
!--- Done with initialization
!---------------------------------------------------------------------

  module_is_initialized = .true.

!=====================================================================
  end subroutine tke_turb_init

!#######################################################################

  subroutine tke_turb_end

!=======================================================================
! --- Close tke module
!=======================================================================

      module_is_initialized = .false.
 
!=====================================================================

  end subroutine tke_turb_end

!#######################################################################

  subroutine tke_pbl_depth(zsfc,zhalf,tke,bstar,h)

!=======================================================================
! --- Estimate PBL depth from tke
!=======================================================================
 
  real,    intent(in), dimension(:,:)   :: zsfc
  real,    intent(in), dimension(:,:,:) :: zhalf, tke
  real,    intent(in), dimension(:,:)   :: bstar
  real,    intent(out),dimension(:,:)   :: h

!=======================================================================
  real, dimension(size(zhalf,1),size(zhalf,2),size(zhalf,3)):: zhalf_ag
  integer:: i,j,k,ix,jx,nlev
  
  ! --- dimensions
  ix = size( zhalf, 1 )
  jx = size( zhalf, 2 )
  nlev = size( zhalf, 3 ) - 1

  ! --- compute heights relative to surface
  do k = 1,nlev+1
    zhalf_ag(:,:,k) = zhalf(:,:,k) - zsfc(:,:)
  end do

  ! --- Determine pbl depth as the height above ground where tke
  !     first falls beneath a critical value tkecrit.
  do j=1,jx
    do i=1,ix
      if (bstar(i,j).gt.0. .and. tke(i,j,nlev).gt.tkecrit) then
        k = nlev
        do while (k.gt.2 .and. tke(i,j,k-1).gt.tkecrit)
          k = k-1
        enddo
        h(i,j) = zhalf_ag(i,j,k) + &
                 (zhalf_ag(i,j,k-1)-zhalf_ag(i,j,k)) &
                 *(tke(i,j,k)-tkecrit)/(tke(i,j,k)-tke(i,j,k-1))
      else
        h(i,j) = 0.                      
      end if
    enddo
  enddo

!=====================================================================

  end subroutine tke_pbl_depth

!#######################################################################

  subroutine parcel_pbl_depth(zsfc,zfull,t,qv,ql,qi,ustar,bstar,h)

!=======================================================================
! --- Estimate PBL depth by parcel lifting
!=======================================================================
 
  real,    intent(in), dimension(:,:)   :: zsfc
  real,    intent(in), dimension(:,:,:) :: zfull, t, qv, ql, qi
  real,    intent(in), dimension(:,:)   :: ustar, bstar
  real,    intent(out),dimension(:,:)   :: h

!=======================================================================

  real, dimension(size(zfull,1),size(zfull,2),size(zfull,3)) :: &
    zfull_ag, slv, hleff
  real, dimension(size(zfull,1),size(zfull,2)) :: &
    svp, h1, ws, k_t_ref, parcelkick
  integer:: i,j,k,ix,jx,nlev
  
  ! --- dimensions
  ix = size( zfull, 1 )
  jx = size( zfull, 2 )
  nlev = size( zfull, 3 )

  ! --- compute heights relative to surface
  do k = 1,nlev
    zfull_ag(:,:,k) = zfull(:,:,k) - zsfc(:,:)
  end do

  ! --- compute slv
  hleff = (min(1.,max(0.,0.05*(t       -tfreeze+20.)))*hlv + &
           min(1.,max(0.,0.05*(tfreeze -t          )))*hls)
     
  slv = cp_air*t + grav*zfull_ag - hleff*(ql + qi)
  slv = slv*(1+d608*(qv+ql+qi))
  slv = slv / cp_air


  ! --- surface parcel properties
  svp  = slv(:,:,nlev)
  h1   = zfull_ag(:,:,nlev)
  call mo_diff(h1, ustar, bstar, ws, k_t_ref)
  ws = max(small,ws/vonkarm/h1)
  svp  = svp*(1.+(parcel_buoy*ustar*bstar/grav/ws) )
  parcelkick = svp*parcel_buoy*ustar*bstar/grav/ws

  ! --- Determine pbl depth
  do j=1,jx
    do i=1,ix

      if (bstar(i,j).gt.0. .and. svp(i,j).gt.slv(i,j,nlev)) then
        k = nlev
        do while (k.gt.2 .and. svp(i,j).gt.slv(i,j,k-1))
          k = k-1
        enddo
        h(i,j) = zfull_ag(i,j,k) + &
                 (zfull_ag(i,j,k-1)-zfull_ag(i,j,k)) &
                 *(slv(i,j,k)-svp(i,j))/(slv(i,j,k)-slv(i,j,k-1))
      else
        h(i,j) = 0.                      
      end if

    enddo
  enddo

!=====================================================================

  end subroutine parcel_pbl_depth

!#######################################################################

 subroutine tri_invert(x,d,a,b,c)

 !  Solves the tridiagonal system of equations.
 !
 !  The following schematic represents the system of equations solved,
 !  where X is the solution.
 !
 !  | B(1)  A(1)   0     0                .......            0    |  |X(1)|   |D(1)|
 !  | C(2)  B(2)  A(2)   0                .......            0    |  |X(2)|   |D(2)|
 !  |  0    C(3)  B(3)  A(3)  0           .......            0    |  | .. |   | .. |
 !  |  ..........................................                 |  | .. | = | .. |
 !  |  ..........................................                 |  | .. |   | .. |
 !  |                                  C(N-2) B(N-2) A(N-2)  0    |  | .. |   | .. |
 !  |                                    0    C(N-1) B(N-1) A(N-1)|  | .. |   | .. |
 !  |                                    0      0    C(N)   B(N)  |  |X(N)|   |D(N)|

  implicit none

  real, intent(out), dimension(:,:,:) :: x
  real, intent(in),  dimension(:,:,:) :: d
  real, intent(inout), dimension(:,:,:) :: a,b,c

  real, dimension(size(x,1),size(x,2),size(x,3)) :: e,f,g,cc
  real, dimension(size(x,1),size(x,2)) :: bb
  integer :: k

  e(:,:,1) = - a(:,:,1)/b(:,:,1)
  a(:,:,size(x,3)) = 0.0

  do  k= 2,size(x,3)
    g(:,:,k) = 1/(b(:,:,k)+c(:,:,k)*e(:,:,k-1))
    e(:,:,k) = - a(:,:,k)*g(:,:,k)
  end do
  cc = c
  bb = 1.0/b(:,:,1)

  f(:,:,1) =  d(:,:,1)*bb
  do k= 2, size(x,3)
    f(:,:,k) = (d(:,:,k) - cc(:,:,k)*f(:,:,k-1))*g(:,:,k)
  end do

  x(:,:,size(x,3)) = f(:,:,size(x,3))
  do k = size(x,3)-1,1,-1
    x(:,:,k) = e(:,:,k)*x(:,:,k+1)+f(:,:,k)
  end do

  return
  end subroutine tri_invert

 subroutine bulk_Ri(t,qv,u,v,t_ca,q_ca,pfull,zfull,phalf,zhalf,Ri)

! This subroutine computes the height where the bulk Richardson number
! reaches 0.25. This height can be regarded an approximation to the
! PBL height. The algorithm was constructed to enable straight-forward 
! comparisons with radiosondes data.
!
! Reference: Seidel et al. 2012 (doi:10.1029/2012JD018143)

!-----------------------------------------------------------------------
!
!      variables
!
!      -----
!      input
!      -----
!
!      t         temperature (K)
!      qv        water vapor specific humidity (kg vapor/kg air)
!      u         zonal wind (m/s)
!      v         meridional wind (m/s)
!      t_ca      canopy air temperature (K)
!      q_ca      canopy air specific humidity (kg vapor/kg air)
!      zfull     height of full levels (m)
!      zhalf     height of half levels (m)
!      pfull     pressure at full levels (Pa)
!      phalf     pressure at half levels (Pa)
!
!      ------------
!      output
!      ------------
!
!      Ri  bluk Richardson 

! Input/output variables

real,            intent(in),    dimension(:,:,:) :: t, qv, u, v
real,            intent(in),    dimension(:,:)   :: t_ca, q_ca
real,            intent(in),    dimension(:,:,:) :: pfull, zfull
real,            intent(in),    dimension(:,:,:) :: phalf, zhalf

! Local variables

integer                                          :: i, j, k
integer                                          :: nlev, nlat, nlon
real, dimension(size(t,1),size(t,2))             :: zsurf
real, dimension(size(t,1),size(t,2),size(t,3))   :: zfull_ag
real, dimension(size(t,1),size(t,2),size(t,3)+1) :: zhalf_ag

real thetavs, thetavh, zs, us, vs, vv
real, dimension(size(t,3)+1) :: Ri, z

! Constants
! Constants
real, parameter :: eps = 1.0e-8
real, parameter :: Ri_crit = 0.25
real, parameter :: p00 = 1000.0e2

!-----------------------------------------------------------------------

! Initialization
  nlev = size(t,3)
  nlat = size(t,2)
  nlon = size(t,1)

! Compute height above surface
  zsurf(:,:) = zhalf(:,:,nlev+1)
  do k = 1, nlev
     zfull_ag(:,:,k) = zfull(:,:,k) - zsurf(:,:)
     zhalf_ag(:,:,k) = zhalf(:,:,k) - zsurf(:,:)
  end do
  zhalf_ag(:,:,nlev+1) = zhalf(:,:,nlev+1) - zsurf(:,:)

! Horizontal loop
  do j=1,nlat
    do i=1,nlon

!      Extended height
       z(1:nlev) = zfull_ag(i,j,:)
       z(nlev+1) = 2.0         ! Assumed reference (canopy) air height

!      Canopy level properties
       thetavs = t_ca(i,j) * (p00/phalf(i,j,nlev+1))**kappa  &
                           * ( 1.0 + 0.61*q_ca(i,j)/(1.0-q_ca(i,j)) )
       zs = z(nlev+1)
       us = 0.0
       vs = 0.0
       Ri(nlev+1) = 0.0

!      Vertical upward loop
       do k = nlev,1,-1

         thetavh = t(i,j,k) * (p00/pfull(i,j,k))**kappa   &
                            * ( 1.0 + 0.61*qv(i,j,k)/(1.0-qv(i,j,k)) )
         vv = max( (u(i,j,k)-us)**2 + (v(i,j,k)-vs)**2, eps )
         Ri(k) = grav * (thetavh-thetavs) * (z(k)-zs) / (thetavs*vv)

       end do

    end do
  end do

end subroutine bulk_Ri


!#######################################################################

end module tke_turb_mod
