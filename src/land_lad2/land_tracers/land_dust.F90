module land_dust_mod

#include "../shared/debug.inc"

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use constants_mod, only: PI, rdgas, GRAV, PSTD_MKS, DENS_H2O
use land_constants_mod, only : d608, kBoltz

use fms_mod, only : error_mesg, FATAL, NOTE, file_exist, &
     close_file, check_nml_error, mpp_pe, mpp_root_pe, stdlog, stdout, string, lowercase
use time_manager_mod, only: time_type, time_type_to_real
use diag_manager_mod, only : register_static_field, &
     send_data
use field_manager_mod, only : parse, MODEL_ATMOS, MODEL_LAND
use tracer_manager_mod, only : get_tracer_index, get_tracer_names, query_method, NO_TRACER

use cana_tile_mod, only : canopy_air_mass_for_tracers
use soil_tile_mod, only : soil_ave_wetness
use snow_tile_mod, only : snow_tile_stock_pe
use vegn_tile_mod, only : vegn_tile_LAI, vegn_tile_SAI
use vegn_data_mod, only:  LU_PAST, LU_CROP, LU_SCND, LU_NTRL
use land_tile_mod, only : land_tile_type, land_tile_grnd_T
use land_tile_diag_mod, only : register_tiled_diag_field, send_tile_data
use land_data_mod, only : lnd, log_version
use land_io_mod, only : read_field
use land_tracers_mod, only : ntcana, isphum
use land_tile_diag_mod, only: set_default_diag_filter
use land_debug_mod, only : is_watch_point
use table_printer_mod


implicit none ; private

! ==== public interfaces =======================================================
public :: land_dust_init
public :: land_dust_end
public :: update_land_dust
! ==== end of public interfaces ================================================

! ==== module constants ========================================================
character(len=*), parameter :: module_name = 'land_dust_mod'
#include "../shared/version_variable.inc"
character(len=*), parameter :: diag_name   = 'land_dust'

type :: dust_data_type
   character(32) :: name = ''
   integer       :: tr_cana = NO_TRACER ! index of this dust tracer in canopy air tracer array
   integer       :: tr_atm  = NO_TRACER ! index of this dust tracer in atmos tracer array
   real          :: radius  = 1e-6 ! dust reference radius, m?
   real          :: density = 2500.0 ! dust material density, kg/m3
   logical       :: do_deposition = .FALSE.
   logical       :: do_emission   = .FALSE.
   real          :: alpha_r, alpha_s ! wet deposition parametrs for rain and snow
   real          :: source_fraction = 0.0 ! fraction of the source allocated to this dust tracer

   integer :: & ! diag field ids
     id_emis,      id_ddep,      id_wdep, &
     id_flux_atm,  id_dfdtr, &
     id_con_v_lam, id_con_g_lam, &
     id_con_v,     id_con_g, &
     id_vdep
end type dust_data_type


! ==== module variables ========================================================

!---- namelist -----------------------------------------------------------------
real :: soil_depth = 0.15 ! depth scale for soil wetness averaging, m
real :: c1         = 10.0 ! adjustment factor for bareness calculations
real :: sai_thresh = 0.05 ! stem area index threshold for dust emission
real :: lai_thresh = 0.1  ! leaf area index threshold for dust emission
real :: sliq_thresh= 0.5  ! soil liquid threshold for dust emission
real :: sice_thresh= 0.05 ! soil ice threshold for dust emission
real :: snow_thresh= 1.0  ! snow threshold, kg/m2
real :: u_min      = 2.0  ! units m/s
real :: u_min_crop = 2.0  ! units m/s
real :: u_min_past = 4.0  ! units m/s
real :: frac_bare_crop = 0.25 ! fraction of bare surface for cropfields
real :: frac_bare_past = 0.25 ! fraction of bare surface for pasture
real :: ch         = 3.5e-10 ! dimensional factor [kg s2/m5]
logical :: dependency_soil_moisture = .false.
character(len=256) :: input_file_name = 'INPUT/dust_source.nc'
character(len=64)  :: input_field_name = 'source'
namelist /land_dust_nml/ &
   soil_depth, c1, lai_thresh, sai_thresh, &
   sliq_thresh, sice_thresh, snow_thresh, dependency_soil_moisture, &
   u_min, u_min_crop, u_min_past, frac_bare_crop, frac_bare_past, &
   ch, input_file_name, input_field_name
!---- end of namelist ----------------------------------------------------------


logical :: module_is_initialized = .FALSE.
logical :: do_dust = .FALSE.  ! if true, do dust calculations
integer :: n_dust_tracers = 0 ! number of dust tracers
type(dust_data_type), allocatable :: trdata(:) ! data for specific dust tracers

real, allocatable :: dust_source(:) ! geographically-dependent dust source


real, save    :: dt     ! fast time step, s

!---- diag field IDs -----------------------------------------------------------
integer :: id_soil_wetness, id_soil_iceness, id_dust_emis, id_dust_source, &
           id_ddep_tot, id_wdep_tot, id_fatm_tot, id_cana_dens, id_u_ts, &
           id_bareness


contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

subroutine land_dust_init (mask, id_ug)
  logical,intent(inout) :: mask(:)
  integer,intent(in) :: id_ug !<Unstructured axis id.

  ! ---- local vars
  logical :: used ! return value from send_data
  integer :: i, tr
  integer :: logunit, outunit, unit, io, ierr
  character(32)  :: name ! tracer name
  character(32)  :: method
  character(128) :: parameters
  real    :: value ! temporary storage for parsing input
  type(table_printer_type) :: table

  ! log module version
  call log_version(version, module_name, &
  __FILE__)
  logunit = stdlog()
  outunit = stdout()

  ! read namelist
#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=land_dust_nml, iostat=io)
     ierr = check_nml_error(io, 'land_dust_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=land_dust_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'land_dust_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=land_dust_nml)
  endif

  ! calculate time step
  dt  = time_type_to_real(lnd%dt_fast) ! store in a module variable for convenience

  ! find out if there are any dust tracers
  n_dust_tracers = 0
  do tr = 1, ntcana
     call get_tracer_names(MODEL_LAND, tr, name)
     if (lowercase(name(1:4))=='dust') then
        n_dust_tracers = n_dust_tracers + 1
     endif
  enddo

  ! log number of dust tracers
  call error_mesg('land_dust_init','number of land dust tracers ='//string(n_dust_tracers), NOTE)

  ! check if any dust is present in the model
  do_dust = (n_dust_tracers > 0)
  if (.not.do_dust) return ! nothing to do

  allocate(trdata(n_dust_tracers))
  i = 0
  do tr = 1, ntcana
     if (.not.mask(tr)) cycle ! sip all tracers that were already registered

     call get_tracer_names(MODEL_LAND, tr, name)
     if (lowercase(name(1:4)).ne.'dust') cycle ! this is not dust, we are not interested

     mask(tr) = .FALSE. ! this is our tracer, mark it non-generic

     i = i+1
     trdata(i)%name    = name
     trdata(i)%tr_cana = tr

     trdata(i)%tr_atm = get_tracer_index(MODEL_ATMOS,name)
     if (trdata(i)%tr_atm == NO_TRACER) cycle ! no such tracer in the atmosphere, nothing to do here

     method = ''; parameters = ''
     if (query_method('dry_deposition', MODEL_ATMOS, trdata(i)%tr_atm, method)) then
        trdata(i)%do_deposition = (index(lowercase(method),'lm3')>0)
     endif

     method = ''; parameters = ''
     if (query_method('wet_deposition', MODEL_ATMOS, trdata(i)%tr_atm, method, parameters)) then
        ! try to parse alphar and alphas parameters for any wet deposition method;
        ! we are not intersted in anything else
        if (parse(parameters,'alphar',value) > 0) trdata(i)%alpha_r = value
        if (parse(parameters,'alphas',value) > 0) trdata(i)%alpha_s = value
     endif

     method = ''; parameters = ''
     if(query_method('emission', MODEL_ATMOS, trdata(i)%tr_atm, method, parameters)) then
        if (trim(method)=='lm3') then
           trdata(i)%do_emission = .true.
           if ( parse(parameters, 'source_fraction', value) > 0 ) then
              trdata(i)%source_fraction = value
           endif
           ! TODO: register emission-related diagnostic fields for this tracer here
        endif
     endif
     method = ''; parameters = ''
     if(query_method('parameters', MODEL_ATMOS, trdata(i)%tr_atm, method, parameters)) then
        if ( parse(parameters, 'dustden', value) > 0 ) trdata(i)%density = value
        if ( parse(parameters, 'dustref', value) > 0 ) trdata(i)%radius  = value
     endif
  enddo


  ! log dust tracer information
  call init_with_headers(table, trdata(:)%name)
  call add_row(table, 'cana.tr.number',  trdata(:)%tr_cana)
  call add_row(table, 'atm.tr.number',   trdata(:)%tr_atm)
  call add_row(table, 'radius',          trdata(:)%radius)
  call add_row(table, 'density',         trdata(:)%density)
  call add_row(table, 'do_deposition',   trdata(:)%do_deposition)
  call add_row(table, 'do_emission',     trdata(:)%do_emission)
  call add_row(table, 'source_fraction', trdata(:)%source_fraction)
  call add_row(table, 'alpha_r',         trdata(:)%alpha_r)
  call add_row(table, 'alpha_s',         trdata(:)%alpha_s)
  call print(table,stdlog())
  call print(table,stdout())

  ! read dust source field
  allocate(dust_source(lnd%ls:lnd%le))
  call read_field( input_file_name, input_field_name, &
       lnd%lon, lnd%lat, dust_source, interp='bilinear' )

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('land')

  ! initialize diagnostic fields
  id_soil_wetness = register_tiled_diag_field(diag_name, 'soil_wetness', (/id_ug/),  &
       lnd%time, 'soil wetness for dust emission calculations', 'unitless', missing_value=-1.0)
  id_soil_iceness = register_tiled_diag_field(diag_name, 'soil_iceness', (/id_ug/),  &
       lnd%time, 'soil ice content for dust emission calculations', 'unitless', missing_value=-1.0)
  id_dust_emis = register_tiled_diag_field(diag_name, 'dust_emis', (/id_ug/),  &
       lnd%time, 'dust emission', 'kg/(m2 s)', missing_value=-1.0)
  id_u_ts = register_tiled_diag_field(diag_name, 'u_ts', (/id_ug/),  &
       lnd%time, 'threshold of wind erosion', 'm/s', missing_value=-1.0)
  id_bareness = register_tiled_diag_field(diag_name, 'bareness', (/id_ug/),  &
       lnd%time, 'bareness measure', 'unitless', missing_value=-1.0)
  id_ddep_tot = register_tiled_diag_field(diag_name, 'dust_ddep', (/id_ug/),  &
       lnd%time, 'total dust dry deposition', 'kg/(m2 s)', missing_value=-1.0)
  id_wdep_tot = register_tiled_diag_field(diag_name, 'dust_wdep', (/id_ug/),  &
       lnd%time, 'total dust wet deposition in canopy air', 'kg/(m2 s)', missing_value=-1.0)
  id_fatm_tot = register_tiled_diag_field(diag_name, 'dust_fatm', (/id_ug/),  &
       lnd%time, 'total dust flux to the atmosphere', 'kg/(m2 s)', missing_value=-1.0)
  id_cana_dens = register_tiled_diag_field(diag_name, 'cana_dens', (/id_ug/),  &
       lnd%time, 'density of canopy air', 'kg/m3', missing_value=-1.0)

  id_dust_source = register_static_field ( diag_name, 'dust_source', (/id_ug/), &
       'topographical dust source', missing_value = -1.0 )

  if (id_dust_source .gt. 0) then
      used = send_data(id_dust_source, &
                       dust_source, &
                       lnd%time)
  endif

  do i = 1,n_dust_tracers
     name = trdata(i)%name
     trdata(i)%id_emis = register_tiled_diag_field(diag_name, trim(name)//'_emis', &
       (/id_ug/),  lnd%time, trim(name)//' emission', &
       'kg/(m2 s)', missing_value=-1.0)
     trdata(i)%id_ddep = register_tiled_diag_field(diag_name, trim(name)//'_ddep', &
       (/id_ug/),  lnd%time, trim(name)//' dry deposition', 'kg/(m2 s)', &
       missing_value=-1.0)
     trdata(i)%id_wdep = register_tiled_diag_field(diag_name, trim(name)//'_wdep', &
       (/id_ug/),  lnd%time, trim(name)//' wet deposition in canopy air', 'kg/(m2 s)', &
       missing_value=-1.0)
     trdata(i)%id_flux_atm = &
       register_tiled_diag_field(diag_name, trim(name)//'_flux_atm', &
       (/id_ug/),  lnd%time, trim(name)//' flux to the atmosphere', &
       'kg/(m2 s)', missing_value=-1.0)
     trdata(i)%id_dfdtr = &
       register_tiled_diag_field(diag_name, trim(name)//'_dfdtr', &
       (/id_ug/),  lnd%time,'derivative of '//trim(name)//' flux to the atmosphere', &
       'kg/(m2 s)', missing_value=-1.0)
     trdata(i)%id_con_v_lam = &
       register_tiled_diag_field(diag_name, trim(name)//'_con_v_lam', &
       (/id_ug/),  lnd%time, 'quasi-laminar conductance between canopy and canopy air for '//trim(name), &
       'm/s', missing_value=-1.0)
     trdata(i)%id_con_g_lam = &
       register_tiled_diag_field(diag_name, trim(name)//'_con_g_lam', &
       (/id_ug/),  lnd%time, 'quasi-laminar conductance between ground and canopy air for '//trim(name), &
       'm/s', missing_value=-1.0)
     trdata(i)%id_con_v = &
       register_tiled_diag_field(diag_name, trim(name)//'_con_v', &
       (/id_ug/),  lnd%time, 'total conductance between canopy and canopy air for'//trim(name), &
       'm/s', missing_value=-1.0)
     trdata(i)%id_con_g = &
       register_tiled_diag_field(diag_name, trim(name)//'_con_g', &
       (/id_ug/),  lnd%time, 'total conductance between ground and canopy air for '//trim(name), &
       'm/s', missing_value=-1.0)
     trdata(i)%id_vdep = &
       register_tiled_diag_field(diag_name, trim(name)//'_vdep', &
       (/id_ug/),  lnd%time, 'sedimentation velocity in canopy air for '//trim(name), &
       'm/s', missing_value=-1.0)
  enddo

  module_is_initialized = .TRUE.

end subroutine land_dust_init


! ==============================================================================
subroutine land_dust_end()
  if (allocated(dust_source))   deallocate(dust_source)
  if (allocated(trdata))     deallocate(trdata)
  module_is_initialized = .FALSE.
end subroutine land_dust_end

! ==============================================================================
! calculates the vertical velocity of dust settling
elemental real function sedimentation_velocity(T,p,r,rho) result(vdep)
   real, intent(in) :: T   ! air temperature, deg K
   real, intent(in) :: p   ! pressure, Pa
   real, intent(in) :: r   ! radius of dust particles, m
   real, intent(in) :: rho ! density of dust particles, kg/m3

   real :: viscosity, free_path, C_c
   viscosity = 1.458E-6 * T**1.5/(T+110.4)   ! Dynamic viscosity [kg/(m s)]
   free_path = 6.6e-8*T/293.15*(PSTD_MKS/p)
   C_c = 1.0 + free_path/r * (1.257+0.4*exp(-1.1*r/free_path)) ! slip correction
   vdep = 2./9.*C_c*GRAV*rho*r**2/viscosity  ! Settling velocity [m/s]
end function sedimentation_velocity

! ==============================================================================
elemental real function laminar_conductance(T,p,sphum,ustar,r,rho,vdep) result(cond)
   real, intent(in) :: T     ! air temperature, deg K
   real, intent(in) :: p     ! pressure, Pa
   real, intent(in) :: sphum ! specific humidity
   real, intent(in) :: ustar ! friction velocity, m/s
   real, intent(in) :: r     ! radius of dust particles, m
   real, intent(in) :: rho   ! density of dust particles, kg/m3
   real, intent(in) :: vdep  ! sedimentation velocity, m/s

   ! ---- local vars
   real :: diff      ! diffusion coefficient for the particles, m2/s
   real :: viscosity ! dynamic viscosity, kg/(m s)
   real :: free_path ! air molecule free path, m
   real :: C_c       ! slip correction, unitless
   real :: rho_air   ! air density, kg/m3
   real :: Sc        ! Schmidt number, unitless
   real :: St        ! Stokes number, unitless

   viscosity = 1.458E-6 * T**1.5/(T+110.4)  ! Dynamic viscosity [kg/(m s)]
   free_path = 6.6e-8*T/293.15*(PSTD_MKS/p)
   C_c = 1.0 + free_path/r * (1.257+0.4*exp(-1.1*r/free_path))
   diff = kBoltz*T*C_c/(6*PI*viscosity*r)
   rho_air = p/(rdgas*T *(1+d608*sphum))
   Sc = viscosity/(rho_air*diff)
   St = vdep*ustar**2*rho_air/(grav*viscosity)
   cond = ustar*(Sc**(-2.0/3.0)+10.0**(-3.0/St))
end function laminar_conductance

! ==============================================================================
subroutine update_land_dust(tile, l, tr_flux, dfdtr, &
     precip_l, precip_s, p_surf, ustar, con_g, con_v )
  type(land_tile_type), intent(inout) :: tile
  integer :: l ! unstructure grid cell indices (global)
  real, intent(in) :: tr_flux(:) ! fluxes of tracers
  real, intent(in) :: dfdtr  (:) ! derivatives of tracer fluxes w.r.t. concentrations
  real, intent(in) :: precip_l   ! liquid precipitation, kg/(m2 s)
  real, intent(in) :: precip_s   ! solid precipitation, kg/(m2 s)
  real, intent(in) :: p_surf     ! atmospheric pressure, N/m2
  real, intent(in) :: ustar ! friction velocity, m/s
  real, intent(in) :: con_g ! aerodynamic conductance between canopy air and ground for tracers
  real, intent(in) :: con_v ! aerodynamic conductance between canopy air and canopy for tracers
        ! note that con_v is for entire canopy, *not* for unit leaf area

  ! ---- local constants
  real , parameter :: &
      R_r = 0.001,  & ! radius of cloud-droplets for rain
      R_s = 0.001     ! radius of cloud-droplets for snow
  real, parameter ::  DENS_SNOW = 500. ! Snow density [kg/m3]
  ! WARNING : this is the density defined in atmos_tracer_utilities; it is not
  !           consistent with the snow density in the land model (a namelist
  !           parameter)
  ! ---- local vars
  integer :: tr, idx
  real    :: rho     ! air density
  real    :: emis(n_dust_tracers)
  real    :: dq      ! canopy air dust tendency per time step
  real    :: rwet    ! radius of wet dust particles, m
  real    :: ratio_r ! ratio of wet/dry volumes
  real    :: rho_wet ! wet dust particle density, kg/m3
  real    :: vdep    ! deposition velocity
  real    :: con_v_lam ! quasi-laminar layer conductance for vegetation, m/s
  real    :: con_g_lam ! quasi-laminar layer conductance for ground, m/s
  real    :: cv, cg  ! total conductances for vegetation and ground surface, m/s
  real    :: f_atm   ! flux of the dust to the atmosphere, kg/(m2 s)
  real    :: ddep    ! dry deposition of the dust, kg/(m2 s)
  real    :: wdep    ! wet deposition of the dust, kg/(m2 s)
  real    :: scav    ! rain and snow scavenging coefficent, 1/s
  real    :: wdep_tot, ddep_tot, fatm_tot ! total values of wet dep., dry dep., and atm. flux
  real    :: LAI     ! leaf area index, m2/m2
  real    :: wind10  ! wind at 10 m above displacement height, m/s
  real    :: ustar_s ! friction velocity at the top of quasi-laminar layer

  if (.not.do_dust) return ! no dust, do nothing in this case
  wind10 = 10.0*ustar

  ! calculate dust sources
  call update_dust_source(tile,l,ustar,wind10,emis)

  rho = p_surf/(rdgas*tile%cana%T*(1+d608*tile%cana%tr(isphum))) ! air density
  wdep_tot = 0.0; ddep_tot = 0.0 ; fatm_tot = 0.0 ! initialize total fluxes to zero for accumulation
  if (associated(tile%vegn)) LAI = vegn_tile_LAI(tile%vegn)
  do tr = 1, n_dust_tracers
     idx = trdata(tr)%tr_cana
     if (trdata(tr)%do_deposition) then
        rwet    = trdata(tr)%radius ! Add any particle growth here
        ratio_r = (trdata(tr)%radius/rwet)**3  ! Ratio dry over wet radius cubic power
        rho_wet = ratio_r*trdata(tr)%density+(1.-ratio_r)*DENS_H2O ! Density of wet aerosol [kg/m3]
        ! calculate deposition velocity
        vdep = sedimentation_velocity(tile%cana%T, p_surf, rwet, rho_wet)
        ! calculate the conductance of quasi-laminar layers
        ! is using the same ustar correct? ustar is value calculated for the
        ! momentum flux from the atmos, is it applicable to the laminar
        ! conductances?
        if (associated(tile%vegn)) then
           con_v_lam = LAI*laminar_conductance(tile%vegn%cohorts(1)%Tv, p_surf, tile%cana%tr(isphum), ustar, rwet, rho_wet, vdep)
        else
           con_v_lam = 0.0
        endif
        cv = 0.0
        if (con_v+con_v_lam > 0) cv = con_v*con_v_lam/(con_v+con_v_lam)
! Estimate ustar at the top of the laminar sublayer by reducing exponentially
! its value from the displacement height (tile%land_d in units of meters)
        ustar_s=ustar*exp(-10.*tile%land_d)
        con_g_lam = laminar_conductance(land_tile_grnd_T(tile), p_surf, tile%cana%tr(isphum), ustar, rwet, rho_wet, vdep)
        cg = con_g*con_g_lam/(con_g+con_g_lam) + vdep
     else
        cv = 0 ; cg = 0; con_v_lam = 0 ; con_g_lam = 0 ; vdep = 0
     endif
     scav = 3./4.*canopy_air_mass_for_tracers* &
           (precip_l*trdata(tr)%alpha_r/R_r/DENS_H2O +  &
            precip_s*trdata(tr)%alpha_s/R_s/DENS_SNOW)
     if ( tile%cana%tr(idx) <= 0 ) scav = 0.0

     ! update tracer concentration -- this needs to be done even when do_deposition
     ! is FALSE, because source might be non-zero
     dq = (-tr_flux(idx)-rho*(cv+cg)*tile%cana%tr(idx)-scav*tile%cana%tr(idx)+emis(tr)) &
        / (canopy_air_mass_for_tracers/dt + dfdtr(idx) + rho*(cv+cg)+scav)
     if (is_watch_point()) then
        write(*,*) 'update_land_dust : ', trim(trdata(tr)%name)
        __DEBUG4__(tile%cana%tr(idx), tr_flux(idx), dfdtr(idx), emis(tr))
        __DEBUG3__(rho,cv,cg)
        __DEBUG2__(dq,tile%cana%tr(idx) + dq)
     endif
     tile%cana%tr(idx) = tile%cana%tr(idx) + dq
     ! ---- final values of the fluxes, for diagnostics
     ddep  = rho*(cv+cg)*tile%cana%tr(idx)
     f_atm = tr_flux(idx)+dfdtr(idx)*dq
     wdep  = scav*tile%cana%tr(idx)
     ! ---- diagnostic section
     call send_tile_data(trdata(tr)%id_con_v_lam,  con_v_lam,  tile%diag)
     call send_tile_data(trdata(tr)%id_con_g_lam,  con_g_lam,  tile%diag)
     call send_tile_data(trdata(tr)%id_con_v,      cv,         tile%diag)
     call send_tile_data(trdata(tr)%id_con_g,      cg,         tile%diag)
     call send_tile_data(trdata(tr)%id_vdep,       vdep,       tile%diag)
     call send_tile_data(trdata(tr)%id_emis,       emis(tr),   tile%diag)
     call send_tile_data(trdata(tr)%id_ddep,       ddep,       tile%diag)
     call send_tile_data(trdata(tr)%id_wdep,       wdep,       tile%diag)
     call send_tile_data(trdata(tr)%id_flux_atm,   f_atm,      tile%diag)
     call send_tile_data(trdata(tr)%id_dfdtr,      dfdtr(idx), tile%diag)
     wdep_tot = wdep_tot + wdep
     ddep_tot = ddep_tot + ddep
     fatm_tot = fatm_tot + f_atm
  enddo
  call send_tile_data(id_ddep_tot,  ddep_tot,  tile%diag)
  call send_tile_data(id_wdep_tot,  wdep_tot,  tile%diag)
  call send_tile_data(id_fatm_tot,  fatm_tot,  tile%diag)
  call send_tile_data(id_cana_dens, rho,       tile%diag)
end subroutine update_land_dust


! ==============================================================================
subroutine update_dust_source(tile, l, ustar, wind10, emis)
  type(land_tile_type), intent(inout) :: tile ! it is only "inout" because diagnostics is sent to it
  integer :: l ! unstructure grid indices 
  real, intent(in) :: ustar ! friction velocity, m/s
  real, intent(in) :: wind10 ! wind at 10 m above displacement height, m/s
  real, intent(inout) :: emis(:)

  ! ---- local vars
  real, parameter :: m     = 0.5
  real, parameter :: sigma = 1.0
  real, parameter :: beta = 100.0
  real :: soil_wetness, soil_iceness ! soil properties for dust source calculations
  real :: snow_lmass, snow_fmass ! snow liquid and frozen water mass
  real :: bareness ! barenes factor, unitless
  real :: lambda, drag
  real :: u_ts, u_thresh ! wind erosion threshold, m/s
  real :: dust_emis
  integer :: tr ! tracer index

  u_ts     = 100.0 ! unrealistically big value guaranteed to be above ustar
  bareness = 1.0  ! value for bare ground
  soil_wetness = 0.0 ; soil_iceness = 0.0 ! for glaciers and lakes
  dust_emis = 0.0 ! default value
  u_thresh = u_min
  if (associated(tile%soil)) then
    ! calculate soil average wetness and "iceness"
    call soil_ave_wetness(tile%soil, soil_depth, soil_wetness, soil_iceness)
    ! calculate snow properties
    call snow_tile_stock_pe(tile%snow, snow_lmass, snow_fmass)
    if (associated(tile%vegn)) then
       if (tile%vegn%landuse .eq. LU_PAST ) then
          u_thresh = u_min_past
          bareness = frac_bare_past
       else if (tile%vegn%landuse .eq. LU_CROP ) then
          u_thresh = u_min_crop
          bareness = frac_bare_crop
       else ! NTRL or SCND
          if ((vegn_tile_LAI(tile%vegn)<lai_thresh) .and. &
              (vegn_tile_SAI(tile%vegn)<sai_thresh)) then
             u_thresh=u_min
             bareness = exp(-2.0*vegn_tile_LAI(tile%vegn)/2.0 &
                            -10.*vegn_tile_SAI(tile%vegn))
          else
             u_thresh=u_ts
             bareness = 0.0
          endif
       endif
    endif

    if ((soil_wetness < sliq_thresh).and.(soil_iceness<sice_thresh)) then
       u_ts = u_thresh
       if (dependency_soil_moisture) then
         u_ts=u_ts*(1.2+0.2*log10(soil_wetness+1.e-5))**2
       endif
    endif
    if ((snow_fmass < snow_thresh).and.(wind10 > u_ts)) then
      dust_emis = ch*bareness*dust_source(l)*(wind10-u_ts)*wind10**2
    endif
  endif

  ! distribute dust emission among dust tracers
  do tr = 1,n_dust_tracers
     if (trdata(tr)%do_emission) &
        emis(tr) = dust_emis*trdata(tr)%source_fraction
  enddo

  ! send data to diagnostics
  call send_tile_data(id_soil_wetness, soil_wetness, tile%diag)
  call send_tile_data(id_soil_iceness, soil_iceness, tile%diag)
  call send_tile_data(id_dust_emis,    dust_emis,    tile%diag)
  call send_tile_data(id_u_ts,         u_ts,         tile%diag)
  call send_tile_data(id_bareness,     bareness,     tile%diag)

end subroutine update_dust_source

end module land_dust_mod
