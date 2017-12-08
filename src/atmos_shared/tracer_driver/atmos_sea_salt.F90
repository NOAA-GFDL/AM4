module atmos_sea_salt_mod
! <DESCRIPTION>
!   This module evaluates the change of mass mixing ratio for seasalt
!   particles due to their emission from ocean, and the removal
!   by gravitational settling. The seasalt particles are transported as dry
!   particles.
! </DESCRIPTION>
! <CONTACT EMAIL="Paul.Ginoux@noaa.gov">
!   Paul Ginoux
! </CONTACT>
!-----------------------------------------------------------------------

use        constants_mod, only : PI, GRAV, RDGAS, DENS_H2O, PSTD_MKS, WTMAIR
use              mpp_mod, only : input_nml_file 
use              fms_mod, only : write_version_number, mpp_pe,  mpp_root_pe, &
                                 open_namelist_file, close_file, file_exist, &
                                 check_nml_error, error_mesg,  &
                                 stdlog, stdout, string, lowercase, &
                                 NOTE, FATAL
use     time_manager_mod, only : time_type
use     diag_manager_mod, only : send_data, register_diag_field
use  atmos_cmip_diag_mod, only : register_cmip_diag_field_2d
use   tracer_manager_mod, only : get_number_tracers, get_tracer_index, &
                                 get_tracer_names, set_tracer_atts, & 
                                 query_method, NO_TRACER
use    field_manager_mod, only : parse, MODEL_ATMOS
use atmos_tracer_utilities_mod, only : wet_deposition, dry_deposition

implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_sea_salt_sourcesink, atmos_sea_salt_init, atmos_sea_salt_end, &
        is_seasalt_tracer, param_ak, param_bk

public do_seasalt, n_seasalt_tracers, seasalt_tracers

!---- version number -----
character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'
character(len=7), parameter :: module_name = 'tracers'

! data type to hold the individual seasalt tracer parameters
type :: seasalt_data_type
   character(32) :: name = '' ! name of the tracer
   character(32) :: seasaltscheme = ' ' 
   integer       :: tr   = NO_TRACER ! index of this seasalt tracer in the atmos tracer array
   real          :: ra=1.e-7, rb=1.e-5   ! boundaries of the size distribution bin, m
   real          :: seasaltref = 2.e-6  ! effective radius of the dry seasalt particles, m
   real          :: seasaltden = 2200.0 ! density of dry seasalt particles, kg/m3
   ! diagnostic IDs
   integer       :: id_seasalt_emis = -1, id_seasalt_setl = -1
end type seasalt_data_type

logical :: do_seasalt = .FALSE.
! ---- module data ----
logical :: module_is_initialized = .FALSE.
integer :: n_seasalt_tracers = 0 ! number of seasalt tracers
type(seasalt_data_type), allocatable :: seasalt_tracers(:) ! parameters for specific seasalt tracers
! ---- identification numbers for diagnostic fields ----
integer :: id_seasalt_emis, id_seasalt_ddep
integer :: id_emiss, id_dryss ! cmip

!---------------------------------------------------------------------
!-------- namelist  ---------
character(len=80) :: scheme = " "
real, save :: coef1
real, save :: coef2
real, save :: ch_fine, ch_coarse
logical    :: use_sj_sedimentation_solver = .FALSE.
!---------------------------------------------------------------------
real :: coef_emis1=-999.
real :: coef_emis2=-999.
real :: coef_emis_fine=-999.
real :: coef_emis_coarse=-999.
real :: critical_sea_fraction = 0.0   ! sea-salt aerosol production
                                      ! occurs in grid cells with
                                      ! ocn_flx_fraction .gt. this value
logical :: ulm_ssalt_deposition=.false.  ! Ulm backward compatibility flag
                                         ! to be removed at Verona

logical            :: ssalt_debug = .false.
integer            :: logunit
namelist /ssalt_nml/  scheme, coef_emis1, coef_emis2, &
                      coef_emis_fine, coef_emis_coarse, &
                      critical_sea_fraction, ulm_ssalt_deposition, &
                      use_sj_sedimentation_solver, ssalt_debug

!-----------------------------------------------------------------------
integer, parameter :: nrh= 65   ! number of RH in look-up table
integer, parameter :: nr = 10   ! number of integration points 
                                      ! The difference with nr=100 & nr=5 < 1e-3
real, dimension(nrh) :: rho_table, growth_table
!! Sea salt hygroscopic growth factor from 35 to 99% RH
!! We start at the deliquescence point of sea-salt for RH=37%
!! Any lower RH doesn't affect dry properties
!! Reference: Tang et al., JGR, v102(D19), 23,269-23,275, 1997. 
data growth_table/1.000, 1.000, 1.396, &
       1.413, 1.428, 1.441, 1.454, 1.466, 1.478, 1.490, 1.501, 1.512, &
       1.523, 1.534, 1.545, 1.555, 1.566, 1.577, 1.588, 1.599, 1.610, &
       1.621, 1.632, 1.644, 1.655, 1.667, 1.679, 1.692, 1.704, 1.717, &
       1.730, 1.743, 1.757, 1.771, 1.786, 1.801, 1.816, 1.832, 1.849, &
       1.866, 1.884, 1.903, 1.923, 1.944, 1.966, 1.990, 2.014, 2.041, &
       2.069, 2.100, 2.134, 2.170, 2.210, 2.255, 2.306, 2.363, 2.430, &
       2.509, 2.605, 2.723, 2.880, 3.087, 3.402, 3.919, 5.048/
!! Seal salt density for 65 RH values from 35% to 99% [g/cm3]
data rho_table/2.160, 2.160, 1.490, &
       1.475, 1.463, 1.452, 1.441, 1.432, 1.422, 1.414, 1.406, 1.398, &
       1.390, 1.382, 1.375, 1.368, 1.361, 1.354, 1.347, 1.341, 1.334, &
       1.328, 1.322, 1.315, 1.309, 1.303, 1.297, 1.291, 1.285, 1.279, &
       1.273, 1.267, 1.261, 1.255, 1.249, 1.243, 1.237, 1.231, 1.225, &
       1.219, 1.213, 1.207, 1.201, 1.195, 1.189, 1.183, 1.176, 1.170, &
       1.163, 1.156, 1.150, 1.142, 1.135, 1.128, 1.120, 1.112, 1.103, &
       1.094, 1.084, 1.074, 1.063, 1.051, 1.038, 1.025, 1.011/

real :: betha


contains

!#######################################################################
! this subroutine calculates tendencies for all seasalt tracers, and reports
! total fields, like total seasalt emission and settling
subroutine atmos_sea_salt_sourcesink ( lon, lat, ocn_flx_fraction, pwt, &
       zhalf, pfull, w10m, t, rh, tracer, dsinku, rdt, dt, Time, is,ie,js,je, kbot)

  real, intent(in) :: lon(:,:), lat(:,:) ! geographical coordinates, units?
  real, intent(in) :: ocn_flx_fraction(:,:) ! fraction of land in the grid cell
  real, intent(in) :: w10m(:,:) ! 10-m wind speed, m/s
  real, intent(in) :: pwt(:,:,:) ! mass of each atmos layer, kg/m2
  real, intent(in) :: zhalf(:,:,:) ! z of half-layers, m(?)
  real, intent(in) :: pfull(:,:,:) ! pressure on layers, Pa
  real, intent(in) :: t(:,:,:) ! temperature of atmosphere, degK
  real, intent(in) :: rh(:,:,:) ! relative humidity 
  real, intent(in) :: tracer(:,:,:,:) ! tracer concentrations
  real, intent(in) :: dsinku(:,:,:) ! dry deposition flux at the surface, for diag only
  real, intent(inout) :: rdt(:,:,:,:) ! tendency of tracers, to be updated for seasalt tracers
  real, intent(in) :: dt ! time step
  type(time_type), intent(in) :: Time ! current model time
  integer, intent(in) :: is, ie, js, je ! boundaries of physical window
  integer, intent(in), optional :: kbot(:,:) ! index of bottom level
  ! NOTE that operations are done on physics window; in particular the sizes
  ! of the arrays correspond to the sizes of physics window.

! ---- local vars
  real, dimension(size(tracer,1),size(tracer,2)) :: &
     seasalt_setl, &     ! seasalt sedimentation at the bottom of the atmos
     all_seasalt_setl, & ! total seasalt sedimentation flux at the bottom of the atmos
     seasalt_emis, &     ! seasalt emission flux at the bottom of the atmos
     all_seasalt_emis    ! total seasalt emission flux at the bottom of the atmos
  real, dimension(size(tracer,1),size(tracer,2),size(tracer,3)) :: &
     seasalt_dt           ! calculated seasalt tendency

  integer :: i
  integer :: kd    ! vertical size of our arrays
  integer :: nseasalt ! atmos tracer number that corresponds to the current seasalt tracer
  logical :: used
  
  kd = size(tracer,3)
  
  ! initialize accumulated deposition and emission fields
  all_seasalt_emis(:,:) = 0.0
  all_seasalt_setl(:,:) = 0.0


  do i = 1,n_seasalt_tracers
     nseasalt = seasalt_tracers(i)%tr
     ! calculate sources and sinks for each individual seasalt tracer
     call atmos_seasalt_sourcesink1(ocn_flx_fraction, pwt, &
        seasalt_tracers(i)%seasaltden, seasalt_tracers(i)%seasaltref, &
        seasalt_tracers(i)%ra, seasalt_tracers(i)%rb, &
        seasalt_tracers(i)%seasaltscheme, &
        zhalf, pfull, w10m, t, rh, &
        tracer(:,:,:,nseasalt), seasalt_dt, seasalt_emis, seasalt_setl, dt, &
        is,ie,js,je, kbot)
     ! update seasalt tendencies
     rdt(:,:,:,nseasalt)=rdt(:,:,:,nseasalt)+seasalt_dt(:,:,:)
     
     ! Send the emission data to the diag_manager for output.
     if (seasalt_tracers(i)%id_seasalt_emis > 0 ) then
       used = send_data ( seasalt_tracers(i)%id_seasalt_emis, seasalt_emis, Time, is_in=is,js_in=js )
     endif
     ! Send the settling data to the diag_manager for output.
     if (seasalt_tracers(i)%id_seasalt_setl > 0 ) then
       used = send_data ( seasalt_tracers(i)%id_seasalt_setl, seasalt_setl, Time, is_in=is,js_in=js )
     endif

     ! Accumulate total emission and deposition fluxes for output
     if (id_seasalt_ddep > 0 .or. id_dryss > 0) then
        ! accumulate total seasalt deposition flux
        all_seasalt_setl(:,:) = all_seasalt_setl(:,:) &
               + seasalt_setl(:,:) + pwt(:,:,kd)*dsinku(:,:,nseasalt) ! shouldn't kd be kbot?
     endif
     if (id_seasalt_emis > 0 .or. id_emiss > 0) then
        ! accumulate total seasalt emission flux
        all_seasalt_emis(:,:) = all_seasalt_emis(:,:) + seasalt_emis(:,:) 
     endif
  enddo
  
  if (id_seasalt_ddep > 0) then
     used = send_data (id_seasalt_ddep, all_seasalt_setl(:,:), Time, is_in=is, js_in=js)
  endif
  if (id_seasalt_emis > 0) then
     used = send_data (id_seasalt_emis, all_seasalt_emis(:,:), Time, is_in=is, js_in=js)
  endif

  ! cmip variables
  if (id_dryss > 0) then
     used = send_data (id_dryss, all_seasalt_setl(:,:), Time, is_in=is, js_in=js)
  endif
  if (id_emiss > 0) then
     used = send_data (id_emiss, all_seasalt_emis(:,:), Time, is_in=is, js_in=js)
  endif
  

end subroutine atmos_sea_salt_sourcesink


!#######################################################################
! Given seasalt properties and atmospheric variables, calculate seasalt tendencies
subroutine atmos_seasalt_sourcesink1 ( &
       ocn_flx_fraction, pwt, &
       seasaltden, seasaltref, seasaltra, seasaltrb,seasalt_scheme, &
       zhalf, pfull, w10m, t, rh, &
       seasalt, seasalt_dt, seasalt_emis, seasalt_setl, dt, is,ie,js,je,kbot)

  real, intent(in),  dimension(:,:)   :: ocn_flx_fraction
  real, intent(in) :: seasaltref ! effective radius of the dry seasalt particles, m
  real, intent(in) :: seasaltra  ! lowest radius m
  real, intent(in) :: seasaltrb  ! highest radius
  real, intent(in) :: seasaltden ! density of dry seasalt particles, kg/m3
  real, intent(in),  dimension(:,:)   :: w10m
  character(32),intent(in) :: seasalt_scheme 
  real, intent(in),  dimension(:,:,:) :: pwt, seasalt
  real, intent(in) :: dt
  real, intent(in),  dimension(:,:,:) :: zhalf, pfull, t, rh
  integer, intent(in),  dimension(:,:), optional :: kbot
  integer, intent(in)  :: is, ie, js, je
  real, intent(out), dimension(:,:,:) :: seasalt_dt
  real, intent(out) :: seasalt_emis(:,:) ! seasalt emission
  real, intent(out) :: seasalt_setl(:,:) ! grav. sedimentation flux at the atmos bottom 

  ! ---- local vars
  integer  i, j, k, id, jd, kd, kb, ir, irh
  real, dimension(size(seasalt,3)) :: setl

  real, parameter :: mtcm = 100.  ! meter to cm
  real, parameter :: mtv  = 1.    ! factor conversion for mixing ratio of seasalt
  real, parameter :: ptmb = 0.01  ! pascal to mb
  integer, parameter :: nstep_max = 5  !Maximum number of cyles for settling
  real :: rhb, rcm,step,rwet
  real :: viscosity, free_path, C_c
  real :: ratio_r, rho_wet_seasalt, seasalt_flux
  real :: rho_air
  real :: a1, a2, Bcoef, r, dr, rmid
  real, dimension(size(pfull,3))  :: vdep, seasalt_conc0, seasalt_conc1
  real, dimension(size(pfull,3))  :: dz, air_dens, qn, qn1
  integer :: istep, nstep

  id=size(seasalt,1); jd=size(seasalt,2); kd=size(seasalt,3)

  seasalt_emis(:,:) = 0.0
  seasalt_setl(:,:) = 0.0
  seasalt_dt(:,:,:) = 0.0

!----------- compute seasalt emission ------------
  if (seasalt_scheme .eq. "Martensson") then !  ie, Martensson et al., JGR-Atm, 2003
          do j=1,jd
            do i=1,id
              if (present(kbot)) then
                 kb=kbot(i,j)
              else
                 kb=kd
              endif
              if (ocn_flx_fraction (i,j).gt.critical_sea_fraction) then
                if (seasaltra .lt. 0.01e-6 .and. mpp_pe() == mpp_root_pe()) then
!$OMP critical (SEA_SALT_WRITE_LOG)
                  write (logunit,*) "***WARNING (atmos_sea_salt): lowest radius of seasalt aerosol should be greater than 0.01E-6 m"
!$OMP end critical (SEA_SALT_WRITE_LOG)
                endif
                r=seasaltra
                dr= (seasaltrb - seasaltra)/float(nr)
                seasalt_flux=0.
                do ir=1,nr
                  rmid=r+dr*0.5   ! Dry radius
                  r=r+dr
                  if (rmid .le. 1.4e-6) then
! Martensson et al., JGR-Atm, 2003
                    seasalt_flux = seasalt_flux + &
                       ch_fine*3.84e-4* 4./3.*pi*seasaltden*1e-3*rmid**2.* &
                       (param_ak(rmid)*t(i,j,kb)+param_bk(rmid))*dr/0.4343
                  else
! Monahan (1986)
                    Bcoef=(coef1-alog10(betha*rmid*1.e6))/coef2
                    seasalt_flux = seasalt_flux + &
                       ch_coarse*1.373*4./3.*pi*seasaltden/betha**2*1.e-12* &
                       (1.+0.057*(betha*rmid*1e6)**1.05)*dr*      &
                       10**(1.19*exp(-(Bcoef**2)))
                  endif
                enddo
                seasalt_emis(i,j) = seasalt_flux*ocn_flx_fraction (i,j)*w10m(i,j)**3.41
              endif
            enddo
          enddo
  else
    if (seasalt_scheme .eq. "Smith" ) then
          do j=1,jd
            do i=1,id
              if (ocn_flx_fraction (i,j).gt.critical_sea_fraction) then
!----------------------------------------------------------------
!    Surface emission of sea salt  (Smith et al. (1993))
!------------------------------------------------------------------
                seasalt_flux            = 0.0
                a1=exp(0.155*w10m(i,j)+5.595)
                a2=exp(2.2082*sqrt(w10m(i,j))-3.3986)
                r = seasaltra* 1.e6
                dr= (seasaltrb - seasaltra)/float(nr)* 1.e6
                do ir=1,nr
                  rmid=r+dr*0.5   ! Dry radius
                  r=r+dr
                  seasalt_flux            = seasalt_flux            + &
                     ch_coarse*4.188e-18*rmid**3*seasaltden*betha*( &
                    + coef1*a1*exp(-3.1*(alog(betha*rmid/2.1))**2) &
                    + coef2*a2*exp(-3.3*(alog(betha*rmid/9.2))**2) )
                enddo
                  seasalt_emis(i,j) = seasalt_flux*ocn_flx_fraction (i,j)
              endif
            enddo
          enddo
    else
! Default is :  Monahan et al. (1986)
          r = seasaltra* 1.e6
          dr= (seasaltrb - seasaltra)/float(nr)* 1.e6
          seasalt_flux=0.
          do ir=1,nr
            rmid=r+dr*0.5   ! Dry radius
            r=r+dr
            Bcoef=(coef1-alog10(betha*rmid))/coef2
            seasalt_flux = seasalt_flux + &
               ch_coarse*1.373*4./3.*pi*seasaltden/betha**2*1.e-18* &
               (1.+0.057*(betha*rmid)**1.05)*dr*      &
               10**(1.19*exp(-(Bcoef**2)))
          enddo
          do j=1,jd
            do i=1,id
              if (ocn_flx_fraction (i,j).gt.critical_sea_fraction) then
                seasalt_emis(i,j) = seasalt_flux*(ocn_flx_fraction (i,j))*   &
                                                           w10m(i,j)**3.41
              endif
            enddo
          enddo
    endif
  endif

  seasalt_dt(:,:,kd)=seasalt_dt(:,:,kd)+seasalt_emis(:,:)/pwt(:,:,kd)*mtv

  rcm=seasaltref*mtcm            ! Particles radius in centimeters
!------------------------------------------
!       Solve at the model TOP (layer plev-10)
!------------------------------------------
  do j=1,jd
    do i=1,id
      setl(:)=0.
      air_dens(:)=pfull(i,j,:)/t(i,j,:)/RDGAS
      dz(:) = pwt(i,j,:)/air_dens(:)
      if (present(kbot)) then
         kb=kbot(i,j)
      else
         kb=kd
      endif
!
! Determine the maximum timestep to avoid particles settling more than 1 layer
!
      nstep=1
      do k=1,kb
        rhb=amin1(0.99,rh(i,j,k))
        rhb=amax1(0.001,rhb)
        irh=max0(1,int(rhb*100.-34.))
        rho_wet_seasalt=rho_table(irh)*1000. !Density of wet sea-salt [kg/m3]
        rwet=seasaltref*growth_table(irh) ! Radius of wet sea-salt [m]
        if(ulm_ssalt_deposition) then   ! Default is .false.
           ! Old claculation left inlined to make it easier to remove
           viscosity = 1.458E-6 * t(i,j,k)**1.5/(t(i,j,k)+110.4)
           free_path=6.6e-8*t(i,j,k)/293.15*(PSTD_MKS/pfull(i,j,k))
           C_c=1. + free_path/seasaltref* &            ! Slip correction [none]
             (1.257+0.4*exp(-1.1*seasaltref/free_path))
           vdep(k)=2./9.*C_c*GRAV*rho_wet_seasalt*rwet**2./viscosity
        else
          ! New calculation drops effective radius seasaltref in favor of rwet
          vdep(k)= sedimentation_velocity(t(i,j,k),pfull(i,j,k),rwet,rho_wet_seasalt) ! Settling velocity [m/s]
        endif
      enddo
      if (use_sj_sedimentation_solver) then
        qn(:)=seasalt(i,j,:)
        qn1(1)=qn(1)*dz(1)/(dz(1)+dt*vdep(1))
        do k=2,kb
          qn1(k)=(qn(k)*dz(k)+dt*qn1(k-1)*vdep(k-1)*air_dens(k-1)/air_dens(k))/(dz(k)+dt*vdep(k))
        enddo
        seasalt_dt(i,j,:)=seasalt_dt(i,j,:)+(qn1(:)-qn(:))/dt
        seasalt_setl(i,j) = qn1(kb)*air_dens(kb)/mtv*vdep(kb)

!---> h1g, 2016-04-05
       if( ssalt_debug ) then
!$OMP CRITICAL

!  if (mpp_pe()==mpp_root_pe()) then
!      write(logunit,'("SALT ",2i5," qn(kb)=",e12.4," qn1(kb)=",e12.4," vdep(kb)=",e12.4," air_dens(kb)=",e12.4," dz(kb)=",e12.4," dt=",e12.4," dust_dt=",e12.4," setl=",e12.4)')  &
!                i,j,qn(kb),qn1(kb),vdep(kb),air_dens(kb),dz(kb),dt,seasalt_dt(i,j,kb),seasalt_setl(i,j)
!  endif
!$OMP END CRITICAL  
       end if ! ssalt_debug
!<--- h1g, 2016-04-05

      else
        do k=1,kb
          step = (zhalf(i,j,k)-zhalf(i,j,k+1)) / vdep(k) / 2.
          nstep = max(nstep, int( dt/ step) )
!!! To avoid spending too much time on cycling the settling in case
!!! of very large particles falling through a tiny layer, impose
!!! maximum speed for the selected nstep_max. This is not physically
!!! correct, but as these particles are very large there will be removed
!!! fast enough to not change significantly their lifetime. The proper
!!! way would be to implement semi-lagrangian technique.
          if (nstep.gt.nstep_max) then
            nstep = nstep_max
            vdep(k)=(zhalf(i,j,k)-zhalf(i,j,k+1))*nstep / 2. /dt
          endif
        enddo
        step = dt / nstep
        seasalt_conc1(:) = seasalt(i,j,:)
        do istep = 1, nstep
          seasalt_conc0(:) = seasalt_conc1(:)
          do k=1,kb
            rho_air = pfull(i,j,k)/t(i,j,k)/RDGAS ! Air density [kg/m3]
            if (seasalt_conc0(k).gt.0.) then
!!!           settling flux [kg/m2/s]
              setl(k)=seasalt_conc0(k)*rho_air/mtv*vdep(k)
            else
              setl(k)=0.
            endif
          enddo
          seasalt_setl(i,j)=seasalt_setl(i,j)+setl(kb)*step
          seasalt_conc1(1) = seasalt_conc0(1) - setl(1)/pwt(i,j,1)*mtv * step
          seasalt_conc1(2:kb)= seasalt_conc0(2:kb) &
            + ( setl(1:kb-1) - setl(2:kb) )/pwt(i,j,2:kb)*mtv * step
          where (seasalt_conc1 < 0 ) seasalt_conc1=0.0
        enddo
        seasalt_dt(i,j,:)=seasalt_dt(i,j,:)+ (seasalt_conc1(:)-seasalt(i,j,:))/dt
        seasalt_setl(i,j)=seasalt_setl(i,j)/dt
      endif
    enddo
  enddo

end subroutine atmos_seasalt_sourcesink1


!#######################################################################
! calculates the vertical velocity of seasalt settling
elemental real function sedimentation_velocity(T,p,rwet,rho_wet_seasalt) result(vdep)
   real, intent(in) :: T            ! air temperature, deg K
   real, intent(in) :: p            ! pressure, Pa
   real, intent(in) :: rwet         ! radius of seasalt particles, m
   real, intent(in) :: rho_wet_seasalt ! density of seasalt particles, kg/m3
 
   real :: viscosity, free_path, C_c
   viscosity = 1.458E-6 * T**1.5/(T+110.4)     ! Dynamic viscosity
   free_path = 6.6e-8*T/293.15*(PSTD_MKS/p)
   C_c = 1.0 + free_path/rwet * &              ! Slip correction [none]
               (1.257+0.4*exp(-1.1*rwet/free_path))
   vdep = 2./9.*C_c*GRAV*rho_wet_seasalt*rwet**2/viscosity  ! Settling velocity [m/s]
end function sedimentation_velocity


!######################################################################
! given a tracer index, returns TRUE if this is one of seasalt tracers
function is_seasalt_tracer(tr) result(ret)
  logical :: ret
  integer, intent(in) :: tr

  integer :: i
  ret = .FALSE.
  do i = 1,n_seasalt_tracers
     ret = (seasalt_tracers(i)%tr==tr)
     if (ret) return
  enddo
end function 
!######################################################################
! Coefficient Ak for the parametrization of sea-salt emission 
! as in Table 1 of Martensson et al. (JGR, 2003)
function param_ak(r) result(ak)
real :: r, d, ak
ak=0
d=2.*r
   if (d.le.0.145e-6) then
      ak=-2.576e35*d**4.+5.932e28*d**3.-2.867e21*d**2.-3.003e13*d-2.881e6
    else
      if (d.le.0.419d-6) then 
        ak=-2.452e33*d**4.+2.404e27*d**3.-8.148e20*d**2.+1.183e14*d-6.743e6
      else
        if (d.le.2.8d-6) then
          ak=1.085e29*d**4.-9.841e23*d**3.+3.132e18*d**2.-4.165e12*d+2.181e6
        endif
      endif
    endif
end function
!######################################################################
! Coefficient Bk for the parametrization of sea-salt emission 
! as in Table 1 of Martensson et al. (JGR, 2003)
function param_bk(r) result(bk)
real :: r, d, bk
  bk=0
  d=2.*r
    if (d.le.0.145e-6) then
      bk=7.188e37*d**4.-1.616e31*d**3.+6.791e23*d**2.+1.829e16*d+7.609e8
    else
      if (d.le.0.419d-6) then
        bk=7.368e35*d**4.-7.310e29*d**3.+2.528e23*d**2.-3.787e16*d+2.279e9
      else
        if (d.le.2.8d-6) then
          bk=-2.859e31*d**4.+2.601e26*d**3.-8.297e20*d**2.+1.105e15*d-5.800e8
        endif
      endif
    endif
end function

!######################################################################
!<SUBROUTINE NAME="atmos_seasalt_init">
!<OVERVIEW>
! The constructor routine for the seasalt module.
!</OVERVIEW>
subroutine atmos_sea_salt_init (lonb, latb, axes, Time, mask)
  real,             intent(in) :: lonb(:,:), latb(:,:) ! grid cell boundaries
  type(time_type),  intent(in) :: Time                 ! model time
  integer,          intent(in) :: axes(4)              ! diagnostic axes
  real, optional,   intent(in) :: mask(:,:,:)

  ! ---- local vars
  integer :: outunit, unit, ierr, io
  integer :: n_atm_tracers ! number of prognostic atmos tracers
  integer :: tr ! atmos tracer iterator
  integer :: i  ! running index of seasalt tracers
  character(32)  :: tr_name    ! tracer name
  character(256) :: longname   ! long and meaningful tracer name
  character(32)  :: method     ! method string for parameter retrieval (not used)
  character(256) :: parameters ! parameter string for seasalt tracer
  real    :: value ! temporary storage for parsing input
  
  if (module_is_initialized) return

  call write_version_number (version, tagname)
  logunit = stdlog()
  outunit = stdout()

  ! read namelist.
  if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ssalt_nml, iostat=io)
    ierr = check_nml_error(io,'ssalt_nml')
#else
    unit =  open_namelist_file ( )
    ierr=1; do while (ierr /= 0)
       read (unit, nml=ssalt_nml, iostat=io, end=10)
       ierr = check_nml_error(io, 'ssalt_nml')
    end do
10  call close_file (unit)
#endif
  endif
 
  ! write namelist to the log file
  if (mpp_pe() == mpp_root_pe()) then
     write (logunit, nml=ssalt_nml)
  endif
      if (coef_emis_coarse .le. -990) then
        ch_coarse = 1.0
      else
        ch_coarse = coef_emis_coarse
      endif
      if (coef_emis_fine .le. -990) then
        ch_fine = 1.0
      else
        ch_fine = coef_emis_fine
      endif
      if (scheme .eq. "Smith") then
        if (coef_emis1 .le. -990) then
          coef1 = 1.0
        else
          coef1 = coef_emis1
        endif
        if (coef_emis2 .le. -990) then
          coef2 = 1.0
        else
          coef2 = coef_emis2
        endif
      else
        if (coef_emis1 .le. -990) then
          coef1 = 0.38
        else
          coef1 = coef_emis1
        endif
        if (coef_emis2 .le. -990) then
          coef2 = 0.65
        else
          coef2 = coef_emis2
        endif
      endif

      betha = growth_table(46)  ! Growth factor at 80% RH

  ! find out if there are any seasalt tracers
  call get_number_tracers(MODEL_ATMOS, num_prog=n_atm_tracers)
  n_seasalt_tracers = 0
  do tr = 1, n_atm_tracers
     call get_tracer_names(MODEL_ATMOS,tr,tr_name)
     if (lowercase(tr_name(1:5))=='ssalt') then
        n_seasalt_tracers = n_seasalt_tracers + 1
     endif
  enddo
  
  ! log the number of seasalt tracers
  call error_mesg('atmos_sea_salt_init','number of atmos seasalt tracers ='//string(n_seasalt_tracers), &
                  NOTE)

  ! check if any seasalt is present in the model and set up the flag
  if (n_seasalt_tracers == 0) return ! nothing to do

  ! get the size of our compute domain, for settl flux storage
  
  ! fill the tracer parameters
  allocate(seasalt_tracers(n_seasalt_tracers))
  i = 0
  ierr = 0 
  do tr = 1, n_atm_tracers
     call get_tracer_names(MODEL_ATMOS,tr,name=tr_name,longname=longname)
     if (lowercase(tr_name(1:5)).ne.'ssalt') cycle ! this is not seasalt, we are not interested  

     i = i+1
     seasalt_tracers(i)%name = tr_name
     seasalt_tracers(i)%tr   = tr
     call set_tracer_atts(MODEL_ATMOS,tr_name,longname,'mmr')
     method=' '
     if(query_method('scheme', MODEL_ATMOS, tr, method, parameters)) then
        seasalt_tracers(i)%seasaltscheme=trim(method)
     else
        seasalt_tracers(i)%seasaltscheme='Monahan'
     endif
     method = ''; parameters = ''
     if(query_method('parameters', MODEL_ATMOS, tr, method, parameters)) then
        call parse_and_check(parameters, tr_name, 'ra',      seasalt_tracers(i)%ra,      ierr)
        call parse_and_check(parameters, tr_name, 'rb',      seasalt_tracers(i)%rb,      ierr)
        call parse_and_check(parameters, tr_name, 'ssaltref', seasalt_tracers(i)%seasaltref, ierr)
        call parse_and_check(parameters, tr_name, 'ssaltden', seasalt_tracers(i)%seasaltden, ierr)
     else
        call error_mesg('atmos_sea_salt_init',&
          '"parameters" line is missing from the field table for seasalt tracer "'//trim(tr_name)//'"', NOTE)
        ierr = ierr+1
     endif
     
     ! Register a diagnostic field : total emission of seasalt
     seasalt_tracers(i)%id_seasalt_emis = register_diag_field ( module_name,     &
                     trim(seasalt_tracers(i)%name)//'_emis', axes(1:2),Time,  &
                     trim(seasalt_tracers(i)%name)//'_emis', 'kg/m2/s',       &
                     missing_value=-999.  )
     ! Register a diagnostic field : total settling of seasalt
     seasalt_tracers(i)%id_seasalt_setl = register_diag_field ( module_name,     &
                     trim(seasalt_tracers(i)%name)//'_setl', axes(1:2),Time,  &
                     trim(seasalt_tracers(i)%name)//'_setl', 'kg/m2/s',       &
                     missing_value=-999.  )
  enddo  
  ! print out information about seasalt tracers
  if (mpp_pe()==mpp_root_pe()) then
     call print_table(logunit)
     call print_table(outunit)
  endif
  ! stop the models is there were any errors reading seasalt parameters
  if (ierr>0) &
     call error_mesg('atmos_sea_salt_init', &
         'ERRORS detected reading one or more parameters of seasalt tracers; '// &
         'look for NOTES from atmos_sea_salt_init above',&
         FATAL)
  
  id_seasalt_ddep = register_diag_field ( module_name, &
      'seasalt_ddep', axes(1:2), Time, &
      'total dry deposition and settling of seasalt', 'kg/m2/s')

  id_seasalt_emis = register_diag_field ( module_name, &
      'seasalt_emis', axes(1:2), Time, &
      'total emission of seasalt', 'kg/m2/s')

  ! cmip variables
  id_dryss = register_cmip_diag_field_2d ( module_name, 'dryss', Time, &
                                  'Dry Deposition Rate of Seasalt', 'kg m-2 s-1', &
             standard_name='tendency_of_atmosphere_mass_content_of_seasalt_dry_aerosol_particles_due_to_dry_deposition')

  id_emiss = register_cmip_diag_field_2d ( module_name, 'emiss', Time, &
                                  'Total Emission Rate of Seasalt', 'kg m-2 s-1', &
             standard_name='tendency_of_atmosphere_mass_content_of_seasalt_dry_aerosol_particles_due_to_emission')

 
  do_seasalt = .TRUE.
  module_is_initialized = .TRUE.

 end subroutine atmos_sea_salt_init
!</SUBROUTINE>

!#######################################################################
! read the specified parameter from the input string, and report a
! non-fatal error if something goes wrong
subroutine parse_and_check(parameters, tr_name, par_name, value, nerrors)
  character(*), intent(in) :: &
     parameters, &              ! string to parse
     tr_name,    &              ! tracer name, for reporting only
     par_name                   ! name of the parameter
  real, intent(inout)  :: value ! result of the parsing
  integer, intent(inout) :: nerrors ! error flag

  real :: tmp ! buffer for input value
  
  if ( parse(parameters, par_name, tmp) > 0 ) then 
     value = tmp
  else
     nerrors = nerrors + 1
     call error_mesg('atmos_sea_salt_init', &
       'ERROR reading parameter "'//trim(par_name)//'" for tracer "'//trim(tr_name)//'"',&
       NOTE)
  endif
end subroutine


!######################################################################
subroutine print_table(unit)
   integer, intent(in) :: unit
   
   integer :: i

   write(unit,'(x,121("-"))')
   write(unit,'(3x,99(x,a16))')'seasalt tr. name','emission scheme ','atm. tr. number','ra','rb', &
       'seasaltref','seasaltden','seasaltscheme'
   write(unit,'(x,121("-"))')
   do i = 1,n_seasalt_tracers
      write(unit,'(x,i2,x,a16,1x,a16,99(x,g16.6))')&
         i, trim(seasalt_tracers(i)%name), trim(seasalt_tracers(i)%seasaltscheme), &
         seasalt_tracers(i)%tr, &
         seasalt_tracers(i)%ra, seasalt_tracers(i)%rb, &
         seasalt_tracers(i)%seasaltref, seasalt_tracers(i)%seasaltden
   enddo
   write(unit,'(x,121("-"))')   
end subroutine 

!#######################################################################
!<SUBROUTINE NAME="atmos_sea_salt_end">
!<OVERVIEW>
!  The destructor routine for the seasalt module.
!</OVERVIEW>
subroutine atmos_sea_salt_end
 
    integer :: i
    module_is_initialized = .FALSE.

    if(.not.do_seasalt) return
    
    do_seasalt = .FALSE.
    deallocate(seasalt_tracers)
end subroutine atmos_sea_salt_end
!</SUBROUTINE>

end module atmos_sea_salt_mod
