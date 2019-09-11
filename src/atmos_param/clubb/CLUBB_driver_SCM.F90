module clubb_driver_mod

!=======================================================================
! GFDL driver for the CLUBB cloud and turbulence parameterization.
!=======================================================================
!
! Authors: Huan Guo, Chris Golaz.
!
! For a scientific description of clubb, refer to the following
! publications:
!
! Larson, V. E., J.-C. Golaz, and W. R. Cotton, 2002: Small-scale and 
! mesoscale variability in cloudy boundary layers: Joint probability 
! density functions. 
! J. Atmos. Sci., 59, 3519--3539. 
! DOI: 10.1175/1520-0469(2002)059<3519:SSAMVI>2.0.CO;2
!
! Golaz, J.-C., V. E. Larson, and W. R. Cotton, 2002a: A PDF-based model 
! for boundary layer clouds. Part I: Method and model description. 
! J. Atmos. Sci., 59, 3540--3551.
! DOI: 10.1175/1520-0469(2002)059<3540:APBMFB>2.0.CO;2
!
! Golaz, J.-C., V. E. Larson, and W. R. Cotton, 2002b: A PDF-based model 
! for boundary layer clouds. Part II: Model results. 
! J. Atmos. Sci., 59, 3552--3571.
! DOI: 10.1175/1520-0469(2002)059<3552:APBMFB>2.0.CO;2
!
! Larson, V. E. and J.-C. Golaz, 2005: Using probability density 
! functions to derive consistent closure relationships among 
! higher-order moments. 
! Mon. Wea. Rev., 133, 1023--1042.
! DOI: 10.1175/MWR2902.1
!
! Golaz, J.-C., V. E. Larson, J. A. Hansen, D. P. Schanen, and
! B. M. Griffin, 2007: Elucidating model inadequacies in a cloud 
! parameterization by use of an ensemble-based calibration framework. 
! Mon. Wea. Rev., 135, 4077--4096.
! DOI: 10.1175/2007MWR2008.1
!
! Guo, H., J.-C. Golaz, L. J. Donner, V. E. Larson, D. P. Schanen, and 
! B. M. Griffin, 2010: Multi-variate probability density functions with 
! dynamics for cloud droplet activation in large-scale models: single 
! column tests. 
! Geosci. Model Dev., 3, 475--486.
! DOI: 10.5194/gmd-3-475-2010
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
! CLUBB and the GCM use different conventions for their vertical grids.
! Here is a comparison of the grid layouts.
!
!         **** CLUBB ****                        **** GCM ****
!
! ========== zm(gr%nz) ====GP=====       ========== phalf(1) ==========
!
! ---------- zt(gr%nz) ----GP-----       ---------- pfull(1) ----------
!
! ========== zm(gr%nz-1) =========       ========== phalf(2) ==========
!
! ---------- zt(gr%nz-1) ---------       ---------- pfull(2) ----------
!
!              .                                         .
!              .                                         .
!              .                                         .
!
! ========== zm(k+1) ===========         ========== phalf(n-1) ========
!
! ---------- zt(k+1) -----------         ---------- pfull(n-1) --------
!
! ========== zm(k) =============         ========== phalf(n) ==========
!
! ---------- zt(k) -------------         ---------- pfull(n) ----------
!
! ========== zm(k-1) ===========         ========== phalf(n+1) ========
!
! ---------- zt(k-1) -----------         ---------- pfull(n+1) --------
!
!              .                                         .
!              .                                         .
!              .                                         .
!
! ========== zm(2) =============         ========== phalf(kdim) =======
!
! ---------- zt(2) -------------         ---------- pfull(kdim) -------
!
! ========== zm(1) ======GP===== surface ========== phalf(kdim+1) =====
! //////////////////////////////         //////////////////////////////
! ---------- zt(1) ------GP-----
!
!
! CLUBB arrays have vertical dimension of gr%nz
!
! gz%nz = kdim + 1
!
! Conversion between CLUBB and host model levels, applicable
! to both full and half levels:
!
! iz_host = kdim + 2 - iz_clubb
! iz_clubb = kdim + 2 - iz_host
!
!-----------------------------------------------------------------------

! --- Imported from FMS modules ---

use       constants_mod, only: RAD_TO_DEG
use             mpp_mod, only: mpp_pe, mpp_root_pe, stdlog, mpp_chksum,        &
                               mpp_clock_id, mpp_clock_begin, mpp_clock_end,   &
                               CLOCK_MODULE_DRIVER
use    diag_manager_mod, only: register_diag_field, send_data
use    time_manager_mod, only: time_type, get_time, set_time, get_date,        &
                               operator(+), operator(-)
use             fms_mod, only: write_version_number, open_file,                &
                               open_namelist_file, check_nml_error,            &
                               file_exist, error_mesg, close_file,             &
                               read_data, write_data,                          &
                               mpp_error, FATAL, NOTE
use   field_manager_mod, only: MODEL_ATMOS
use  tracer_manager_mod, only: get_number_tracers, get_tracer_index,           &
                               get_tracer_names
use   aerosol_types_mod, only: aerosol_type
use     aer_ccn_act_mod, only: aer_ccn_act_init, aer_ccn_act_end
use   aer_ccn_act_k_mod, only: aer_ccn_act_k
use        ice_nucl_mod, only: ice_nucl_wpdf_init, ice_nucl_wpdf_end
#ifdef CLUBB
use       alt_cloud_mod, only: alt_cloud

! --- Imported from clubb modules ---
use                  clubb_core, only: advance_clubb_core, cleanup_clubb_core, &
                                       setup_clubb_core
use             clubb_precision, only: time_precision
use             constants_clubb, only: cp, em_min, ep1, ep2, grav, kappa, Lv,  &
                                       p0, Rd, rt_tol, thl_tol, w_tol_sqd,     &
                                       zero_threshold
use                  error_code, only: fatal_error, set_clubb_debug_level
use                  grid_class, only: gr, setup_grid_heights, zm2zt, zt2zm
use           parameter_indices, only: nparams
use            parameters_model, only: cloud_frac_min, sclr_dim, T0
use          parameters_tunable, only: read_parameters
use                  stats_subs, only: stats_init, stats_begin_timestep,       &
                                       stats_end_timestep, stats_finalize
use             stats_variables, only: l_stats, l_stats_samp
use               T_in_K_module, only: thlm2T_in_K
use variables_prognostic_module, only: cloud_cover, cloud_frac, edsclrm,       &
                                       edsclrm_forcing, exner,                 &
                                       invrs_rho_ds_zm, invrs_rho_ds_zt,       &
                                       pdf_params, p_in_Pa, rcm, rcm_in_layer, &
                                       RH_crit, rho, rho_ds_zm, rho_ds_zt,     &
                                       rho_zm, rtm, rtm_forcing, rtp2,         &
                                       rtpthlp, sclrm, sclrm_forcing, sclrp2,  &
                                       sclrprtp, sclrpthlp, temp_clubb, thlm,  &
                                       thlm_forcing, thlp2, thv_ds_zm,         &
                                       thv_ds_zt, um, um_forcing, up2, up2,    &
                                       upwp_sfc, upwp, upwp_sfc, vm,           &
                                       vm_forcing, vp2, vpwp, vpwp_sfc,        &
                                       wm_zm, wm_zt, wp2, wp3,                 &
                                       wpedsclrp_sfc, wprcp, wprtp,            &
                                       wprtp_sfc, wpsclrp, wpsclrp_sfc,        &
                                       wpthlp, wpthlp_sfc,                     &
                                       wprtp_forcing,    wpthlp_forcing,       &
                                       rtp2_forcing,     thlp2_forcing,        &
                                       rtpthlp_forcing
#endif
implicit none 
public :: clubb_setup, & 
          clubb_init,  & 
          clubb,       & 
          clubb_end

!--------------------- version number ----------------------------------
character(len=128)   :: version = '$Id$'
character(len=128)   :: tagname = '$Name$'

logical              :: module_is_initialized = .false.
character(len=32)    :: tracer_units, tracer_name
character(len=128)   :: diaglname

character(len=9)     :: mod_name = 'clubb'
real                 :: missing_value = -999.
logical              :: used

! ----- Misc parameters -----

! CLUBB related
real            :: T0_in = 300.0
logical         :: do_clubb_conservation_checks = .true.
logical         :: do_clubb_diss_heat           = .true.

! For the calculation of surface momentum fluxes from u_star
real, parameter :: gust_const = 1.0

! Aerorols related
integer, parameter  :: n_totmass = 4
integer, parameter  :: n_imass   = 12

! For the ice nucleation
integer :: rh_act_opt = 2           ! option as to which rh to use for ice
                                    ! nucleation. 1 => use gridbox mean;
                                    ! 2 => use rh in non-cloudy portion of
                                    ! box when the cloudy part of grid box
                                    ! is less than cf_thresh_nucl; use grid
                                    ! box mean otherwise.

real    :: sea_salt_scale_onl = 1.  ! scaling factor to convert seasalt 
                                    ! aerosol tracer to seasalt available
                                    ! for use as condensation nuclei

real    :: cf_thresh_nucl = 0.98    ! threshold cloud fraction below which
                                    ! the rh of the non-cloudy portion of 
                                    ! the gribox is used in determining
                                    ! ice nucleation (when rh_act_opt is >1)

real    :: var_limit_ice = -999.    ! lower limit to the vertical velocity
                                    ! variance of the pdf used for ice
                                    ! nucleation; default value implies 
                                    ! that no limit is imposed.
LOGICAL :: do_var_limit_ice


! ----- Variables for "clubb_setting_nml" namelist -----

logical   ::                       &
   use_sclr_HOC       = .false.,   &
   use_cloud_cover    = .false.,   &
   do_drp_evap_cf     = .true.,    &
   spread_host_tdcy   = .true.

real ::  cloud_frac_min_in = 0.008  ! minimum cloud fraction

character(len=6)::  saturation_formula_in = "GFDL  " ! "bolton" approx., "flatau" approx, or "GFDL" approx.

logical   ::                        &
   do_BL_gauss         =  .false.,  &
   do_diffK_gauss      =  .false.,  &
   do_quadrature_gauss =  .true.

real ::  Var_w = 0.49  ! variance of vertical velocity (w) m2/s2

real ::  clubb_dt  = 60.0       ! time-step used in clubb [s]

integer ::   debug_level = 2     ! Amount of debugging information

logical :: do_aeromass_clubb_const = .false.
real    :: aeromass_clubb_const = 2.25e-12
real    :: rh_factor = 1.0

real    :: host_dx = 4.0e5
real    :: host_dy = 4.0e5

integer :: icheck_temp = 0

integer ::  do_conv_flag_clubb = 1
real    ::  avg_deltaz = 50.0

integer ::  do_alt_cloud = 0

! ---> h1g, 2012-06-14, add option of do_liquid_only_in_clubb, 
!                       by default .false. (that is, do liquid and ice together)
logical                 ::  do_liquid_only_in_clubb = .false.
! <--- h1g, 2012-06-14

real    ::  qcvar_clubb_min = 0.001,               &
            qcvar_clubb_max = 2.0,                 &
            scale_negative_qcvar_factor = 1.0,     &
            qcvar_factor = 1.0,                    &
            qcvar_factor_const = 1.0,              &
            rcm_min      = 0.0

logical                 :: use_qcvar_in_cloud = .false.
logical                 :: do_scale_negative_qcvar = .false.
logical                 :: use_avg4qcvar = .false.
logical                 :: use_avg_qcvar = .false.


namelist /clubb_setting_nml/                                &
         use_sclr_HOC,                                      &
         use_cloud_cover,                                   &
         do_drp_evap_cf,                                    &
         spread_host_tdcy,                                  &
         cloud_frac_min_in,                                 &
         saturation_formula_in,                             &
         do_quadrature_gauss, do_diffK_gauss, do_BL_gauss,  &
         Var_w,                                             &
         clubb_dt,                                          &
         debug_level,                                       &
         do_aeromass_clubb_const, aeromass_clubb_const,     &
         host_dx, host_dy, icheck_temp,                     &
         do_conv_flag_clubb,                                &
         avg_deltaz,                                        &
         do_alt_cloud,                                      &
         do_liquid_only_in_clubb,                           &
         rh_act_opt, sea_salt_scale_onl,                    &
         var_limit_ice,                                     &
         cf_thresh_nucl,                                    &
         rh_factor,                                         &
         qcvar_clubb_min,  qcvar_clubb_max,                 &
         use_qcvar_in_cloud, do_scale_negative_qcvar,       & 
         use_avg4qcvar, use_avg_qcvar,                      &
         qcvar_factor, scale_negative_qcvar_factor,         &
         qcvar_factor_const,  rcm_min


! ----- Variables related "clubb_stats_setting_nml" namelist -----

type(time_type) :: Time_stats_init
integer         :: unit_stats
integer         :: istats = -99999,  &  ! index for stats column
                   jstats = -99999
integer         :: fstats = -99999      ! index for file name

logical         :: do_stats = .false.   ! whether statistics are computed and output
character(len=10) :: & 
      stats_fmt = "grads"  ! File format for stats; typically GrADS.


real, dimension(10) :: lon_stats = -9999.0,  &  ! lon and lat of grid point selected for statistical output
                       lat_stats = -9999.0      ! in degrees.

character(len=100), dimension(10) :: & 
      fname_prefix = "clubb"       ! Prefix of stats filenames, to be followed by, for example "_zt"
#ifdef CLUBB
real(kind=time_precision) :: & 
      stats_tsamp = 0.0,   & ! Stats sampling interval [s]
      stats_tout = 0.0       ! Stats output interval   [s]
namelist /clubb_stats_setting_nml/ & 
      do_stats, lon_stats, lat_stats, fname_prefix, stats_tsamp, stats_tout, stats_fmt
#endif

 
! ----- Local arrays -----

real, dimension(:),     allocatable :: momentum_heights, thermodynamic_heights

integer :: kdim                              ! number of vertical levels
integer :: iz_clubb                          ! vertical index for clubb levels
integer :: iz_host                           ! vertical index for host model levels
integer :: nt                                ! total no. of tracers
integer :: ntp                               ! total no. of prognostic tracers
integer :: it                                ! index for loop over tracers
integer :: clubb_core_clock                  ! clock index
integer :: nsphum, nql, nqi, nqa, nqn, nqni  ! tracer indices for stratiform clouds
integer :: nupwp, nvpwp, nup2, nvp2, nwp2, & ! tracer indices for clubb variables
           nwprtp, nwpthlp, nrtp2, nthlp2, &
           nrtpthlp, nwp3,                 &
           nrhcrit1, nrhcrit2

! ----- Diagnostic variables -----

integer :: id_sulfate,                   &  ! sulfate aerosol mass concentration
           id_seasalt_sub,               &  ! sub-micron sea salt aerosol mass concentration
           id_seasalt_sup,               &  ! sup-micron sea salt aerosol mass concentration
           id_om,                        &  ! organic matter aerosol mass concentration
           id_enth_clubb_col,            &  ! column enthalpy conservation
! ---> h1g, 2012-06-04, add dissipation heat from kenetic energy
           id_diss_heat_clubb,           &
! <--- h1g, 2012-06-04
           id_wat_clubb_col,             &  ! column water conservation
           id_udt_clubb,                 &  ! clubb u tendency
           id_vdt_clubb,                 &  ! clubb u tendency
           id_tdt_clubb ,                &  ! clubb t tendency
           id_drp_evap_clubb,            &  ! macrophysics drop evaporation
           id_aer_ccn_act_clubb,         &  ! CCN acrivation
           id_Ndrop_act_clubb,           &  ! drop activation
           id_drp_flux_clubb,            &  ! drop flux
           id_qndt_clubb_trsport_only,   &
           id_ice_evap_clubb,            &  ! macrophysics ice crystal evaporation
           id_aer_ice_act_clubb,         &  ! IN activation
           id_Icedrop_act_clubb,         &  ! ice crystal activation
           id_icedrop_flux_clubb,        &  ! ice crystal flux
           id_qnidt_clubb_trsport_only,  &
           id_wm_clubb,                  &  ! wm
           id_omega_clubb,               &  ! omega
           id_wp3_clubb,                 &  ! wp3 (full levels)
           id_qcvar_clubb,               &  ! qcvar_clubb (full levels)
           id_qcvar_cf_clubb,            &  ! qcvar_cf_clubb (full levels)
           id_wp2_clubb,                 &  ! wp2 (half levels)
           id_upwp_clubb,                &  ! upwp (half levels)
           id_vpwp_clubb,                &  ! vpwp (half levels)
           id_up2_clubb,                 &  ! up2 (half levels)
           id_vp2_clubb,                 &  ! vp2 (half levels)
           id_wprtp_clubb,               &  ! wprtp (half levels)
           id_wpthlp_clubb,              &  ! wpthlp (half levels)
           id_rtp2_clubb,                &  ! rtp2 (half levels)
           id_thlp2_clubb,               &  ! thlp2 (half levels)
           id_rtpthlp_clubb,             &  ! rptthlp (half levels)
           id_wp2_vert_avg_clubb,        &  ! vertically averaged wp2 = sum(wp2*rho*dz)/sum(rho*dz)
           id_up2_vert_avg_clubb,        &  ! vertically averaged up2 = sum(up2*rho*dz)/sum(rho*dz)
           id_vp2_vert_avg_clubb,        &  ! vertically averaged vp2 = sum(vp2*rho*dz)/sum(rho*dz)
           id_rtp2_vert_avg_clubb,       &  ! vertically averaged rtp2 = sum(rtp2*rho*dz)/sum(rho*dz)
           id_thlp2_vert_avg_clubb,      &  ! vertically averaged thlp2 = sum(thlp2*rho*dz)/sum(rho*dz)
           id_rcm_clubb,                 &  ! cloud water content in CLUBB
           id_rcp2_clubb,                &  ! cloud water variance on grid box in CLUBB
           id_cf_clubb,                  &  ! cloud fraction in CLUBB 
           id_cf_avg_clubb



integer, dimension(:), allocatable :: id_tracer_clubb
integer, dimension(:), allocatable :: id_tracer_clubb_col

!module variables obtained during initialization from other modules

! ---> h1g, 2012-08-23, add option of applying surface fluxes in host-model
!                       by default .true. (that is, applying surface fluxes in host-model)
logical                 ::  l_host_applies_sfc_fluxes 
! <--- h1g, 2012-08-23
!RSH Obtain from lscloud_driver:
logical :: do_ice_nucl_wpdf    
logical :: do_liq_num

!-----------------------------------------------------------------------

contains

#ifdef CLUBB
  !=====================================================================

  subroutine clubb(is, ie, js, je, lon, lat,                  &
                   Time_next,                                 &
                   dtmain,                                    &
                   phalf, pfull, zhalf, zfull, omega_avg,     &
                   t, q, r, u, v,                             &
                   u_star, b_star, q_star,                    &
                   tdt, qdt, rdt, rdiag, udt, vdt,                   &
                   dcond_ls_liquid, dcond_ls_ice,             &
                   Ndrop_act_clubb, Icedrop_act_clubb,        &
                   totalmass1, imass1,  &
                   diff_t_clubb,                              &
                   qcvar_clubb,                               &
                   tdt_shf,  qdt_lhf ,                        &                   
                   Aerosol, mask,                             &
                   mc_full,                                   &
                   conv_frac_clubb,                           &
                   convective_humidity_ratio_clubb)

  !---------------------------------------------------------------------
  !
  !  is,ie,js,je                        input
  !    starting/ending subdomain i,j indices of data in the physics
  !    window being integrated
  !
  !  lon, lat                           input
  !    longitude, latitude in radians
  !
  !  Time_next                          input
  !    next time, used for diagnostics
  !
  !  dtmain                             input
  !    time step in seconds
  !
  !  phalf                              input
  !    pressure at half levels in pascals
  !
  !  pfull                              input
  !    pressure at full levels in pascals
  !
  !  zhalf                              input
  !    height at half levels in meters
  !
  !  zfull                              input
  !    height at full levels in meters
  !
  !  omega_avg                          input
  !    grid average omega at full levels in pascals per second
  !
  !  t, q                               input
  !    temperature (t) [deg k] and specific humidity
  !    of water vapor (q) [kg/kg] at full model levels,
  !
  !  r                                  input/output
  !    tracer fields at full model levels, unit varies,
  !
  !  u, v                               input
  !    zonal and meridional wind [m/s] at full model levels
  !  
  !  u_star                             input
  !    friction velocity [m/s]
  !
  !  b_star                             input
  !    buoyancy scale    [m/s^2/K]
  !
  !  q_star                             input
  !    moisture scale    dimensionless
  !
  !  tdt, qdt                           input/output
  !    temperature (tdt) [deg k/sec] and specific humidity of water 
  !    vapor (qdt) tendency [1/sec]
  !
  !  rdt                                input/output
  !    tracer tendencies , unit varies, 
  !
  !  udt, vdt                           input/output
  !    zonal and meridional wind tendencies [m/s/s]
  !
  !  dcond_ls_liquid                    output
  !    large-scale liquid condensation rate
  !
  !  dcond_ls_ice                       output
  !    large-scale ice condensation rate
  !
  !  Ndrop_act_clubb                    output
  !    in-cloud (NOT DOMAIN) averaged activated droplet number 
  !    concentration (#/kg)
  !
  !  Icedrop_act_clubb                  output
  !    in-cloud (NOT DOMAIN) averaged activated ice-crystal number 
  !    concentration (#/kg)
  !
  !  diff_t_clubb                       output
  !    tracer eddy diffusivity from clubb
  !
  !  Aerosol                            input, optional
  !    different aerosol species mass concentration
  !
  !  mask                               input, optional
  !    real array indicating the point is above the surface if equal 
  !    to 1.0 and indicating the point is below the surface if equal 
  !    to 0.
  !
  !  mc_full                            input, optional
  !    (total) net convective mass flux [kg/m2/s]
  !    used to adjust environmental vertical velocity 
  !    due to convection:
  !      omega = omega_avg + mc_full*grav
  !
  !  conv_frac_clubb                    input, optional
  !    total convective cloud fraction
  !
  !  convective_humidity_ratio_clubb    input, optional
  !    ratio of the grid average specific humidity 
  !    to environmental specific humidity 
  !
  !---------------------------------------------------------------------

  ! ----- Calling arguments -----

  integer, intent(in)                           ::  is, ie, js, je
  real, intent(in), dimension(:,:)              ::  lon, lat
  type(time_type), intent(in)                   ::  Time_next
  real, intent(in)                              ::  dtmain
  real, intent(in), dimension(:,:,:)            ::  phalf, pfull, zhalf, zfull, omega_avg
  real, intent(in), dimension(:,:,:)            ::  t, q, u, v
  real, intent(in), dimension(:,:,:,:)          ::  r
  real, intent(in), dimension(:,:)              ::  u_star, b_star, q_star
  real, intent(inout), dimension(:,:,:)         ::  tdt, qdt, udt, vdt
  real, intent(inout), dimension(:,:,:,:)       ::  rdt
  real, intent(inout), dimension(:,:,:,ntp+1:)  ::  rdiag
  real, intent(out), dimension(:,:,:)           ::  dcond_ls_liquid
  real, intent(out), dimension(:,:,:)           ::  dcond_ls_ice
  real, intent(out), dimension(:,:,:)           ::  Ndrop_act_clubb
  real, intent(out), dimension(:,:,:)           ::  Icedrop_act_clubb
  real, intent( in), dimension(:,:,:,:)         ::  totalmass1, imass1
  real, intent(out), dimension(:,:,:)           ::  diff_t_clubb
  real, intent(out), optional, dimension(:,:,:) ::  qcvar_clubb   
  type(aerosol_type), intent(in), optional      ::  Aerosol
  real, intent(in), optional, dimension(:,:,:)  ::  mask
  real, intent(in), optional, dimension(:,:,:)  ::  mc_full
  real, intent(in), optional, dimension(:,:,:)  ::  conv_frac_clubb
  real, intent(in), optional, dimension(:,:,:)  ::  convective_humidity_ratio_clubb
  real, intent(in), optional, dimension(:,:)    ::  tdt_shf,  qdt_lhf
  ! ----- Local variables -----

  type(time_type)                             Time

  real, dimension( size(t,3)+1 )                :: diff_t_1d
  real, dimension( size(t,3)+1 )                :: qcvar_clubb_1d, rcp2_zt, cf_avg, rcm_avg, rcp2_avg, qcvar_clubb_avg, qcvar_clubb_tmp, &
                                                   qcvar_cf_clubb_1d
 
  real, dimension( size(t,1), size(t,2), size(t,3) ) :: omega
  real, dimension( size(t,1), size(t,2), size(t,3) ) :: env_qv_scale
  real, dimension( size(t,1), size(t,2), size(t,3) ) :: env_condensate_scale
  real, dimension( size(rdt,1), size(rdt,2), size(rdt,3), size(rdt,4) )      :: rdt_orig

  real, dimension( size(t,3) )                   :: tmp_host
  real, dimension( size(t,3)+1 )                 :: tmp_clubb

  ! Passive scalar concentration due to pure transport [{units vary}/s]
  real, dimension( size(t,3)+1, sclr_dim )       :: sclrm_trsport_only

  ! updated tracers concentration after clubb
  real, dimension( size(t,3)+1 )                 :: trs_clubb

  ! Pressure on momentum levels                    [Pa]
  real, dimension(size(t,3)+1) ::  P_in_Pa_zm
  real, dimension(size(t,3)+1) ::  exner_zm

  !aerosol mass concentration
  real,  dimension(size(t,3)+1, n_totmass) :: aeromass_clubb
  real,  dimension(size(t,3)+1, n_imass)   :: aeroimass_clubb

  real,  dimension(size(t,3)+1)                     ::  cloud_frac_before_clubb
  real,  dimension(size(t,3)+1)                     ::  Ndrop_max
  real,  dimension(size(t,3)+1)                     ::  rcm_before_clubb
  real,  dimension(size(t,3)+1)                     ::  drp_before_clubb
  real,  dimension(size(t,1),size(t,2),size(t,3))   ::  aer_ccn_act_clubb, drp_evap_clubb
  real,  dimension(size(t,1),size(t,2),size(t,3)+1) ::  drp_flux_clubb
  real,  dimension(size(t,3)+1)                     ::  Ncrystal_max
  real,  dimension(size(t,3)+1)                     ::  ice_before_clubb
  real,  dimension(size(t,1),size(t,2),size(t,3))   ::  aer_ice_act_clubb, ice_evap_clubb
  real,  dimension(size(t,1),size(t,2),size(t,3)+1) ::  icedrop_flux_clubb
 
  real    ::  sum_sclrm_trsport
  integer ::  isum_sclrm

  integer ::  isub_clubb, nsub_clubb, i_dtmain, i_clubb_dt
  integer :: err_code    ! valid run?
  integer :: ix, iy, k
  integer :: idim, jdim
  integer :: i_totmass, i_imass

  ! tendencies
  real, dimension(size(t,1),size(t,2),size(t,3)) :: &
              uten_clubb,           &     ! u-component tendency
              vten_clubb,           &     ! u-component tendency
              tten_clubb                  ! air temperature tendency

  !tracer tendencies
  real, dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: rdt_clubb    ! tracer transport tendency

  !tracer transport-only tendencies
  real, dimension(size(r,1),size(r,2),size(r,3)) :: qndt_clubb_trsport_only      !drop number
  real, dimension(size(r,1),size(r,2),size(r,3)) :: qnidt_clubb_trsport_only     !ice number


  ! higher order terms and fluxes from clubb
  !at full levels
  real, dimension(size(t,1), size(t,2), size(t,3)) :: wp3_clubb
  real, dimension(size(t,1), size(t,2), size(t,3)) :: wm_clubb
  real, dimension(size(t,1), size(t,2), size(t,3)) :: rcm_clubb, rcp2_clubb, cf_clubb, qcvar_cf_clubb_3d, cf_avg_clubb_3d

  !at half levels
  real, dimension(size(t,1),size(t,2),size(t,3)+1) :: &
                 wp2_clubb,      &
                 upwp_clubb,     &
                 vpwp_clubb,     &
                 up2_clubb,      &
                 vp2_clubb,      &
                 wprtp_clubb,    &
                 wpthlp_clubb,   &
                 rtp2_clubb,     &
                 thlp2_clubb,    &
                 rtpthlp_clubb

  ! at surface
  real, dimension(size(t,1),size(t,2)) :: wp2_vert_avg_clubb, &
                                          up2_vert_avg_clubb, &
                                          vp2_vert_avg_clubb, &
                                          rtp2_vert_avg_clubb, &
                                          thlp2_vert_avg_clubb
  real                                 :: tmp0D

  ! check mass and energy conservation
  real, dimension(size(t,1),size(t,2))            :: tmp2D_check
  real, dimension(size(t,1),size(t,2),size(t,3))  :: pmass_3d
  !  temperature fix to force exact enthalpy/energy conservation
  real :: tten_clubb_fix,  uv_dissipation

! -->h1g, 2012-06-07,  add Total ice-phase water mixing ratio
  real, dimension(size(t,3)+1)  ::  rfrzm, khzm, khzt
! <--h1g, 2012-06-07

  type(time_type)           :: Time_clubb
  integer                   :: itime_elapsed
  real(kind=time_precision) :: time_elapsed

  integer(kind=8) :: ichk   ! for output of mpp_chksum

  !---------------------------------------------------------------------
  ! ----- Begin Code -----

  ! initialization

  idim = size(t,1)
  jdim = size(t,2)
  kdim = size(t,3)

  tmp_host         = 0.0
  trs_clubb        = 0.0

  uten_clubb       = 0.0
  vten_clubb       = 0.0
  tten_clubb       = 0.0
  rdt_clubb        = 0.0

  wp3_clubb        = 0.0
  wp2_clubb        = 0.0
  upwp_clubb       = 0.0
  vpwp_clubb       = 0.0
  up2_clubb        = 0.0
  vp2_clubb        = 0.0
  wprtp_clubb      = 0.0
  wpthlp_clubb     = 0.0
  rtp2_clubb       = 0.0
  thlp2_clubb      = 0.0
  rtpthlp_clubb    = 0.0

  wp2_vert_avg_clubb = 0.0
  up2_vert_avg_clubb = 0.0
  vp2_vert_avg_clubb = 0.0
  rtp2_vert_avg_clubb = 0.0
  thlp2_vert_avg_clubb = 0.0


  i_dtmain    = dtmain
  i_clubb_dt  = clubb_dt
  Time = Time_next - set_time( i_dtmain )

  if ( mod(i_dtmain, i_clubb_dt) /= 0 ) then
    call error_mesg (  &
      'clubb_driver_mod',  &
      'time step in host atmosphere must be an intergal multiple of time step in clubb',  &
      FATAL )
  endif

  nsub_clubb = i_dtmain/i_clubb_dt
  if ( nsub_clubb < 1 ) then
    call error_mesg (  &
      'clubb_driver_mod',  &
      'time step in host atmosphere must be larger than or equal to the time step in clubb',  &
      FATAL )
  endif

  if ( nqn > 0 ) then 
    aeromass_clubb           = 0.0
    Ndrop_max                = 0.0
    rcm_before_clubb         = 0.0
    cloud_frac_before_clubb  = 0.0
    drp_before_clubb         = 0.0
    aer_ccn_act_clubb        = 0.0
    drp_evap_clubb           = 0.0
    Ndrop_act_clubb          = 0.0
    qndt_clubb_trsport_only  = 0.0
    sclrm_trsport_only(:, 1) = 0.0
  endif

  if ( nqni > 0 ) then
    if (do_ice_nucl_wpdf .eqv. .false.)  &
      call error_mesg ('clubb_driver_mod', 'nqni > 0, but do_ice_nucl_wpdf is false', FATAL)
    Icedrop_act_clubb        = 0.0
    Ncrystal_max             = 0.0
    ice_before_clubb         = 0.0
    ice_evap_clubb           = 0.0
    qnidt_clubb_trsport_only = 0.0
    sclrm_trsport_only(:, 2) = 0.0
  endif
 
  ! initialize large-scale condensation (liquid + ice), which is the difference of
  dcond_ls_liquid  = 0.0
  dcond_ls_ice     = 0.0


  omega                = omega_avg
  env_qv_scale         = 1.0
  env_condensate_scale = 1.0
  rdt_orig             = rdt

  if ( do_conv_flag_clubb > 0 ) then

    if ( present(mc_full) ) then
      omega = omega + mc_full*grav
    endif
    if ( present( conv_frac_clubb ) ) then
      omega = omega / ( 1.0 -  conv_frac_clubb )
    endif

    if( do_conv_flag_clubb >1 ) then

      if ( present( convective_humidity_ratio_clubb ) ) then
        env_qv_scale =  convective_humidity_ratio_clubb
        where ( convective_humidity_ratio_clubb .lt. 0. )
          env_qv_scale = 1.0
        end where
      endif

      if ( present( conv_frac_clubb ) ) then
        env_condensate_scale = 1.0 -  conv_frac_clubb
      endif

    endif

  endif
  
! ---> h1g, 2012-08-23
  if ( do_clubb_conservation_checks )   tmp2d_check = 0.0
! <--- h1g, 2012-08-23

! ---> h1g, 2012-08-28
      if ( .not. l_host_applies_sfc_fluxes )  then
         if ( .not. ( present(tdt_shf) ) )  &
             call error_mesg ('clubb_driver_mod','l_host_applies_sfc_fluxes is false, but tdt_shf is absent', FATAL)
         if ( .not. ( present(qdt_lhf) ) )  &
             call error_mesg ('clubb_driver_mod','l_host_applies_sfc_fluxes is false, but qdt_lhf is absent', FATAL)
      endif
! <--- h1g, 2012-08-28


  do iy = 1, jdim
    do ix = 1, idim
 
    ! Re-adjust model height (from host model)
    call host2clubb_half(zhalf(ix,iy,:), momentum_heights)
    call host2clubb_full(zfull(ix,iy,:), thermodynamic_heights)
    call setup_grid_heights(.true., 2, avg_deltaz, momentum_heights(1), momentum_heights, thermodynamic_heights)

    ! Load higher order terms and high-res results
    call clubb_3d_2_1d(ix, iy, rdiag)

    ! Pressure on thermodynamic points [Pa]
    call host2clubb_full(pfull(ix,iy,:), p_in_Pa)
    p_in_Pa( 1 ) = phalf(ix, iy,  kdim+1) ! surface pressure  [Pa]

    ! Exner on thermodynamic points
    exner = ( p_in_Pa/p0 )**kappa
    exner(1) = exner(2)
 
    ! Pressure on momentum points[Pa]
    call host2clubb_half(phalf(ix,iy,:), p_in_Pa_zm)

    ! Exner on momentum points
    exner_zm = ( p_in_Pa_zm/p0 )**kappa

    ! Calculate tempearture on thermodynamic levels
    call host2clubb_full(t(ix,iy,:), temp_clubb(:))
    temp_clubb(1) = temp_clubb(2)
          
    ! Calculate liquid potential tempearture on thermodynamic levels
! ---> h1g, 2012-06-14
  if( do_liquid_only_in_clubb ) then
    tmp_host(:) = t(ix,iy,:)  &
                  - Lv/Cp*(r(ix,iy,:,nql) )  &
                          /env_condensate_scale(ix,iy,:)
  else
    tmp_host(:) = t(ix,iy,:)  &
                  - Lv/Cp*(r(ix,iy,:,nql) + r(ix,iy,:,nqi))  &
                          /env_condensate_scale(ix,iy,:)
  endif
! <--- h1g, 2012-06-14
    call host2clubb_full(tmp_host, thlm)
    thlm = thlm/exner
    thlm(1) = thlm(2)
 
    ! Calculate total water content on thermodynamic levels [kg/kg]
! ---> h1g, 2012-06-14
  if( do_liquid_only_in_clubb ) then
    tmp_host(:) = r(ix,iy,:,nsphum)/env_qv_scale(ix,iy,:)       &
                + r(ix,iy,:,nql)/env_condensate_scale(ix,iy,:)
  else
    tmp_host(:) = r(ix,iy,:,nsphum)/env_qv_scale(ix,iy,:)       &
                + r(ix,iy,:,nql)/env_condensate_scale(ix,iy,:)  &
                + r(ix,iy,:,nqi)/env_condensate_scale(ix,iy,:)
  endif
! <--- h1g, 2012-06-14
    call host2clubb_full(tmp_host, rtm)
    rtm(1) = rtm(2)

    ! Calculate liquid water content on thermodynamic levels [kg/kg], 
    ! add liquid and ice together as liquid as inputs for clubb
! ---> h1g, 2012-06-14
  if( do_liquid_only_in_clubb ) then
    tmp_host(:) = r(ix,iy,:,nql)/env_condensate_scale(ix,iy,:)
  else
    tmp_host(:) = r(ix,iy,:,nql)/env_condensate_scale(ix,iy,:) &
                + r(ix,iy,:,nqi)/env_condensate_scale(ix,iy,:)
  endif
! <--- h1g, 2012-06-14
    call host2clubb_full(tmp_host,rcm)
    rcm(1) = rcm(2)
 
    ! Calculate cloud fraction thermodynamic levels
    tmp_host(:) = r(ix,iy,:,nqa)
    call host2clubb_full(tmp_host,cloud_frac)
    cloud_frac(1) = cloud_frac(2)
 
    ! Calculate thv_ds_zt, thv_ds_zm for clubb
    thv_ds_zt = temp_clubb/exner * (1.0+ep2*(rtm-rcm))**kappa
    do iz_clubb = 1,kdim
      thv_ds_zm(iz_clubb) = zt2zm(temp_clubb/exner, iz_clubb) &
                       * ( 1.0 + ep2 * max( zt2zm( rtm - rcm, iz_clubb), &
                                            zero_threshold ) )**kappa
    enddo
    iz_clubb = kdim+1
    thv_ds_zm(iz_clubb) = 2.0*thv_ds_zm(iz_clubb-1)-thv_ds_zm(iz_clubb-2)

    ! Calculate droplet number concentration
    if (nqn > 0) then
      ! domain average droplet number concentration
      tmp_host(:) = r(ix,iy,:,nqn)/env_condensate_scale(ix,iy,:)
      if (use_sclr_HOC) then
        call host2clubb_full(tmp_host,sclrm(:,1))
        do iz_clubb = 1,kdim+1
          if (cloud_frac(iz_clubb) <= cloud_frac_min) sclrm(iz_clubb,1) = 0.0
        enddo
      else
        call host2clubb_full(tmp_host,edsclrm(:,1))
        do iz_clubb = 1,kdim+1
          if (cloud_frac(iz_clubb) <= cloud_frac_min) edsclrm( iz_clubb, 1) = 0.0
        enddo
      endif
    endif

    ! Calculate ice number concentration
    if (nqni > 0) then
      ! domain average ice number concentration
      tmp_host(:)= r(ix,iy,:,nqni)/env_condensate_scale(ix,iy,:) 
      if (use_sclr_HOC) then
        call host2clubb_full(tmp_host,sclrm(:,2))
        do iz_clubb = 1,kdim+1
          if (cloud_frac(iz_clubb) <= cloud_frac_min) sclrm(iz_clubb,2) = 0.0
        enddo
      else
        call host2clubb_full(tmp_host,edsclrm(:,2))
        do iz_clubb = 1,kdim+1
          if (cloud_frac(iz_clubb) <= cloud_frac_min) edsclrm(iz_clubb,2) = 0.0
        enddo
      endif !  use_sclr_HOC
    endif !  nqni > 0

    ! In order to have better conservation properties, air density is calculated 
    ! using hydro-static approximation, rather than gas state equation
    ! at thermo-dynamic levels
    do iz_clubb = 2,kdim+1
      rho(iz_clubb) = - (p_in_Pa_zm (iz_clubb) - p_in_Pa_zm (iz_clubb-1)) &
                      * gr%invrs_dzt(iz_clubb)/grav
    enddo
    rho(1) = rho(2)  ! surface air density  [kg/m3]

    ! at momentum levels
    do iz_clubb = 2, kdim
      rho_zm( iz_clubb ) = - (p_in_Pa(iz_clubb+1) - p_in_Pa(iz_clubb)) &
                           * gr%invrs_dzm(iz_clubb)/grav
    enddo
    rho_zm(1)      = rho(1)
    rho_zm(kdim+1) = rho(kdim+1)

    rho_ds_zt = rho
    rho_ds_zm = rho_zm
    invrs_rho_ds_zt = 1.0 / rho_ds_zt
    invrs_rho_ds_zm = 1.0 / rho_ds_zm

    ! horizontal velocity on thermodynamic points
    call host2clubb_full(u(ix,iy,:), um)
    um(1) = um(2)
    call host2clubb_full(v(ix,iy,:), vm)
    vm(1) = vm(2)

    ! vertical velocity on thermodynamic points
    call host2clubb_full( omega(ix,iy,:), wm_zt)
    ! Boundary conditions on subsidence (thermodynamic grid)
    wm_zt    = -wm_zt/grav/rho
    wm_zt(1) = 0.0        ! Below surface

    ! vertical velocity on momentum points
    wm_zm = zt2zm(wm_zt)

    ! Boundary conditions on subsidence (mom. grid)
    wm_zm(1)      = 0.0               ! At surface
    wm_zm(kdim+1) = 0.0               ! Model top

    ! Apply host tendencies at once
    if (.not.spread_host_tdcy) then
      call add_host_tdcy(                                                    &
             ix, iy, dtmain,                                                 & ! in
             udt, vdt, tdt, rdt, env_qv_scale, env_condensate_scale,         & ! in
             u_star, b_star, q_star,                                         & ! in
             um, vm, thlm, rtm, cloud_frac, rcm, edsclrm, sclrm,             & ! inout
             upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc )                       ! out
    end if

    ! initiallize eddy diffusivity coefficients to zero
    diff_t_1d = 0.0
    qcvar_clubb_1d = qcvar_clubb_max
    qcvar_cf_clubb_1d = 0.0 
    cf_avg    = 0.0
    rcm_avg   = 0.0
    rcp2_avg  = 0.0
    qcvar_clubb_avg = 0.0
    qcvar_clubb_tmp = qcvar_clubb_max


    err_code = 0
    do isub_clubb=1,nsub_clubb  

      Time_clubb = Time + set_time( int( isub_clubb*clubb_dt ) )

      ! Apply host tendencies by spreading them over the time step
      if (spread_host_tdcy) then
        call add_host_tdcy(                                                    &
               ix, iy, clubb_dt,                                               & ! in
               udt, vdt, tdt, rdt, env_qv_scale, env_condensate_scale,         & ! in
               u_star, b_star, q_star,                                         & ! in
               um, vm, thlm, rtm, cloud_frac, rcm, edsclrm, sclrm,             & ! inout
               upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc )                       ! out
      end if

      if (nqn > 0) then
        ! droplet evaporation
        rcm_before_clubb        = rcm
        cloud_frac_before_clubb = cloud_frac
        if (use_sclr_HOC) then
          drp_before_clubb(:) = max( 0.0, sclrm(:,1) )
        else
          drp_before_clubb(:) = max( 0.0, edsclrm(:,1) )
        endif
      endif ! nqn > 0

      if (nqni > 0) then
        ! ice particle evaporation
        if (use_sclr_HOC) then
          ice_before_clubb(:) = max( 0.0, sclrm(:,2) )
        else
          ice_before_clubb(:) = max( 0.0, edsclrm(:,2) )
        endif
      endif ! nqni > 0

      ! Activate clubb internal stats
      if ( do_stats .and. ix.eq.istats .and. iy.eq.jstats ) then
        call get_time( Time_clubb - Time_stats_init, itime_elapsed )
        time_elapsed = itime_elapsed
        l_stats = .true.
        call stats_begin_timestep( time_elapsed )
      else
        l_stats = .false.
        l_stats_samp = .false.
      end if

      thlm_forcing = 0.0
      rtm_forcing  = 0.0
      um_forcing   = 0.0
      vm_forcing   = 0.0
      wprtp_forcing   = 0.0 
      wpthlp_forcing  = 0.0 
      rtp2_forcing    = 0.0 
      thlp2_forcing   = 0.0 
      rtpthlp_forcing = 0.0 

      wpsclrp_sfc(1:sclr_dim)       = 0.0
      wpedsclrp_sfc(1:sclr_dim)     = 0.0
      sclrm_forcing(:,1:sclr_dim)   = 0.0
      edsclrm_forcing(:,1:sclr_dim) = 0.0

! ---> h1g, 2012-06-07, add  rfrzm: Total ice-phase water mixing ratio (consistent with clubb(r5845))
      rfrzm = 0.0
! <--- h1g, 2012-06-07


! ---> h1g, 2012-08-23
      if ( .not. l_host_applies_sfc_fluxes )  then
           wpthlp_sfc = tdt_shf(ix, iy) * gr%dzt(1)
           wprtp_sfc  = qdt_lhf(ix, iy) * gr%dzt(1)
           if ( do_clubb_conservation_checks ) &
             tmp2d_check( ix, iy ) = tmp2d_check( ix, iy ) - ( wpthlp_sfc*Cp )*rho(2)
      endif
! <--- h1g, 2012-08-23


      call mpp_clock_begin( clubb_core_clock )

      !     if( abs(lon(ix,iy)*RAD_TO_DEG - 289.82) <= 0.1 .and. abs(lat(ix,iy)*RAD_TO_DEG - 51.87) <= 0.1 .and. isub_clubb>= 9) then
      !   write(*,'(a, 3i5,  9f16.8, i5)') 'before ', ix, iy, isub_clubb,  &
      !        lon(ix,iy)*RAD_TO_DEG,  lat(ix,iy)*RAD_TO_DEG, b_star( ix, iy ),  wpthlp_sfc,  q_star( ix, iy ),  wprtp_sfc, u_star( ix, iy ),  upwp_sfc, vpwp_sfc,  err_code
      !   print*, 'wm_zt = ',  wm_zt
      !   print*, '  p_in_Pa  =' , p_in_Pa
      !   print*, '  rho = ', rho
      !   print*, '  thv_ds_zt  = ',  thv_ds_zt
      !   print*, '   um = ', um
      !   print*, '   vm = ', vm
      !   print*, '   upwp =', upwp
      !   print*, '   vpwp =', vpwp
      !   print*, '    thlm = ',  thlm
      !   print*, '   rtm = ', rtm


      !    endif


      call advance_clubb_core(  &
             .true., clubb_dt, 0.0, momentum_heights(1),  &          ! Intent(in)
             thlm_forcing, rtm_forcing, um_forcing, vm_forcing,  &   ! Intent(in)
             sclrm_forcing, edsclrm_forcing, wprtp_forcing, &        ! Intent(in) 
             wpthlp_forcing, rtp2_forcing, thlp2_forcing, &          ! Intent(in) 
             rtpthlp_forcing, wm_zm, wm_zt, &                        ! Intent(in) 
             wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,  &           ! Intent(in)
             wpsclrp_sfc, wpedsclrp_sfc,  &                          ! Intent(in)
             p_in_Pa, rho_zm, rho, exner,  &                         ! Intent(in)
             rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,  &               ! Intent(in)
             invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt,  &               ! Intent(in)
             rfrzm, &                                                ! Intent(in), Total ice-phase water mixing ratio [kg/kg], currently 0
             um, vm, upwp, vpwp, up2, vp2,  &                        ! Intent(inout)
             thlm, rtm, wprtp, wpthlp,  &                            ! Intent(inout)
             wp2, wp3, rtp2, thlp2, rtpthlp,  &                      ! Intent(inout)
             sclrm,  &                                               ! Intent(inout)
             rcp2_zt, sclrm_trsport_only,  &                         ! Intent(inout)
             sclrp2, sclrprtp, sclrpthlp,  &                         ! Intent(inout)
             wpsclrp, edsclrm, err_code,  &                          ! Intent(inout)
             RH_crit, do_liquid_only_in_clubb,   &                   ! Intent(out)
             rcm, wprcp, cloud_frac,  &                              ! Intent(out)
             rcm_in_layer, cloud_cover,  &                           ! Intent(out)
             khzm, khzt, &                                           ! Intent(out)
             pdf_params )                                            ! Intent(out)

      !    if( abs(lon(ix,iy)*RAD_TO_DEG - 289.82) <= 0.1 .and. abs(lat(ix,iy)*RAD_TO_DEG - 51.87) <= 0.1 ) &
      !   write(*,'(a, 3i5,  6f16.8, i5)') 'after ',  ix, iy, isub_clubb,  &
      !        lon(ix,iy)*RAD_TO_DEG,  lat(ix,iy)*RAD_TO_DEG,   wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,  err_code


      call mpp_clock_end( clubb_core_clock )
      if ( fatal_error(err_code) ) then
        call error_mesg ('clubb_driver_mod', 'fatal error in advance_clubb_core', FATAL)
      endif

      ! store eddy diffusivity coefficient
      diff_t_1d = diff_t_1d + khzt
      

      ! option to use cloud_cover and rcm_in_layer instead of cloud_frac and rcm
      if (use_cloud_cover) then
        cloud_frac = cloud_cover
        rcm = rcm_in_layer
      endif

      ! store cloud fraction (cloud_frac), cloud water content (rcm), grid box cloud water variance (rcp2)
      cf_avg    = cf_avg + cloud_frac
      rcm_avg   = rcm_avg + rcm
      rcp2_avg  = rcp2_avg + rcp2_zt


      ! compute qcvar
      do iz_clubb = 2, kdim+1
         if ( use_qcvar_in_cloud ) then
           if ( rcm(iz_clubb)*rcm(iz_clubb)/qcvar_clubb_max < &
                (cloud_frac(iz_clubb)*rcp2_zt(iz_clubb)-(1.0-cloud_frac(iz_clubb))*rcm(iz_clubb)*rcm(iz_clubb)) ) then

                qcvar_clubb_tmp(iz_clubb) = rcm(iz_clubb)*rcm(iz_clubb) &
                  /(cloud_frac(iz_clubb)*rcp2_zt(iz_clubb)-(1.0-cloud_frac(iz_clubb))*rcm(iz_clubb)*rcm(iz_clubb))

                qcvar_clubb_tmp(iz_clubb) = qcvar_clubb_tmp(iz_clubb) * qcvar_factor_const
                qcvar_clubb_tmp(iz_clubb) = max(qcvar_clubb_tmp(iz_clubb), qcvar_clubb_min)
           else

             if( cloud_frac(iz_clubb)*rcp2_zt(iz_clubb)-(1.0-cloud_frac(iz_clubb))*rcm(iz_clubb)*rcm(iz_clubb) <= 0.0 ) then
                qcvar_clubb_tmp(iz_clubb) = qcvar_clubb_min
                if( do_scale_negative_qcvar ) then
                   if( cloud_frac(iz_clubb)*rcp2_zt(iz_clubb)-(1.0-cloud_frac(iz_clubb))*rcm(iz_clubb)*rcm(iz_clubb) < 0.0 & 
                       .and. rcm(iz_clubb) >= rcm_min ) then
                     qcvar_clubb_tmp(iz_clubb) = qcvar_clubb_tmp(iz_clubb)                            &
                                                /(1.0 - qcvar_factor*rcm(iz_clubb)*rcm(iz_clubb)      & 
  /(cloud_frac(iz_clubb)*rcp2_zt(iz_clubb)-(1.0-cloud_frac(iz_clubb))*rcm(iz_clubb)*rcm(iz_clubb)) )

                     qcvar_clubb_tmp(iz_clubb) =  max(qcvar_clubb_tmp(iz_clubb), qcvar_clubb_min/scale_negative_qcvar_factor)
                   endif
                endif
             endif

           endif 
         else
           if ( rcm(iz_clubb)*rcm(iz_clubb)/qcvar_clubb_max < rcp2_zt(iz_clubb) ) then
              qcvar_clubb_tmp(iz_clubb) = rcm(iz_clubb)*rcm(iz_clubb)/rcp2_zt(iz_clubb)
              qcvar_clubb_tmp(iz_clubb) = qcvar_clubb_tmp(iz_clubb) * qcvar_factor_const
              qcvar_clubb_tmp(iz_clubb) = max(qcvar_clubb_tmp(iz_clubb), qcvar_clubb_min)
           endif
         endif         
      enddo
      qcvar_clubb_avg = qcvar_clubb_avg + qcvar_clubb_tmp



      ! end of timestep for clubb internal stats
      if ( do_stats .and. ix.eq.istats .and. iy.eq.jstats ) then
        call stats_end_timestep()
      end if

      ! convert liquid potential temperature to air temperature, and calculate temperature tendency
      thlm(1) =  thlm(2)
      temp_clubb = thlm2T_in_K(thlm, exner, rcm)

!--> Note: can we move the drop and ice crystal evaporation code into a separate
!          subroutine?

     ! droplet number evaporation due to macro-physics
     if(nqn>0)  then 
      if(use_sclr_HOC) then
       if(do_drp_evap_cf) then 
          do iz_clubb=2, kdim+1
             iz_host = kdim + 2 - iz_clubb
             if( ( cloud_frac_before_clubb(iz_clubb) - cloud_frac(iz_clubb) ) > 0.0 .and. cloud_frac(iz_clubb)>0.0) then
! droplet evaporation
                 drp_evap_clubb( ix, iy,  iz_host) = &
                    drp_evap_clubb( ix, iy,  iz_host)   &
                   - min( sclrm( iz_clubb, 1 ), drp_before_clubb( iz_clubb )    &
                      * ( cloud_frac_before_clubb(iz_clubb) -cloud_frac(iz_clubb) ) &
                      /cloud_frac_before_clubb( iz_clubb ) )

                 sclrm( iz_clubb, 1 ) = sclrm( iz_clubb, 1 ) &
                                    - min( sclrm( iz_clubb, 1 ), drp_before_clubb( iz_clubb )    &
                                       * ( cloud_frac_before_clubb(iz_clubb) -cloud_frac(iz_clubb) ) &
                                      /cloud_frac_before_clubb( iz_clubb ) )

             endif

             if( rcm(iz_clubb) <= 0.0 .or. cloud_frac(iz_clubb) <=cloud_frac_min ) then
! in clear areas, droplet concentration should be 0
                 drp_evap_clubb( ix, iy,  iz_host) = &
                   drp_evap_clubb( ix, iy,  iz_host)   &
                  - sclrm( iz_clubb, 1 )
 
                 sclrm( iz_clubb, 1 ) = 0.0

             endif
          enddo !enddo iz_clubb
       else !NOT do_drp_evap_cf
          do iz_clubb=2, kdim+1
             iz_host = kdim + 2 - iz_clubb
             if( ( rcm_before_clubb(iz_clubb) - rcm(iz_clubb) ) > 0.0 .and. rcm(iz_clubb)>0.0) then
! droplet evaporation
               drp_evap_clubb( ix, iy,  iz_host) = &
                 drp_evap_clubb( ix, iy,  iz_host)   &
                 - min( sclrm( iz_clubb, 1 ), drp_before_clubb( iz_clubb )    &
                   * ( rcm_before_clubb(iz_clubb) -rcm(iz_clubb) ) &
                   /rcm_before_clubb( iz_clubb ) )

               sclrm( iz_clubb, 1 ) = sclrm( iz_clubb, 1 ) &
                                    - min( sclrm( iz_clubb, 1 ), drp_before_clubb( iz_clubb )    &
                                      * ( rcm_before_clubb(iz_clubb) -rcm(iz_clubb) ) &
                                      /rcm_before_clubb( iz_clubb ) )

             endif

             if( rcm(iz_clubb) <= 0.0 .or. cloud_frac(iz_clubb) <=cloud_frac_min ) then
! in clear areas, droplet concentration should be 0
                drp_evap_clubb( ix, iy,  iz_host) = &
                  drp_evap_clubb( ix, iy,  iz_host)   &
                 - sclrm( iz_clubb, 1 )
 
                sclrm( iz_clubb, 1 ) = 0.0

             endif
          enddo !enddo iz_clubb
       endif !endif do_drp_evap_cf

     else ! not use scalar HOC

       if(do_drp_evap_cf) then 
          do iz_clubb=2, kdim+1
            iz_host = kdim + 2 - iz_clubb
            if( ( cloud_frac_before_clubb(iz_clubb) - cloud_frac(iz_clubb) ) > 0.0 .and. cloud_frac(iz_clubb)>0.0) then
! droplet evaporation
                 drp_evap_clubb( ix, iy,  iz_host) = &
                    drp_evap_clubb( ix, iy,  iz_host)   &
                   - min( edsclrm( iz_clubb, 1 ), drp_before_clubb( iz_clubb )    &
                      * ( cloud_frac_before_clubb(iz_clubb) -cloud_frac(iz_clubb) ) &
                         /cloud_frac_before_clubb( iz_clubb ) )

                  edsclrm( iz_clubb, 1 ) = edsclrm( iz_clubb, 1 ) &
                                            - min( edsclrm( iz_clubb, 1 ), drp_before_clubb( iz_clubb )    &
                                              * ( cloud_frac_before_clubb(iz_clubb) -cloud_frac(iz_clubb) ) &
                                             /cloud_frac_before_clubb( iz_clubb ) )

            endif

            if( rcm(iz_clubb) <= 0.0 .or. cloud_frac(iz_clubb) <=cloud_frac_min ) then
! in clear areas, droplet concentration should be 0
               drp_evap_clubb( ix, iy,  iz_host) = &
                   drp_evap_clubb( ix, iy,  iz_host)   &
                  - edsclrm( iz_clubb, 1 )
 
               edsclrm( iz_clubb, 1 ) = 0.0

             endif
          enddo !enddo iz_clubb
        else !NOT do_drp_evap_cf
          do iz_clubb=2, kdim+1
             iz_host = kdim + 2 - iz_clubb
             if( ( rcm_before_clubb(iz_clubb) - rcm(iz_clubb) ) > 0.0 .and. rcm(iz_clubb)>0.0) then
! droplet evaporation
                  drp_evap_clubb( ix, iy,  iz_host) = &
                    drp_evap_clubb( ix, iy,  iz_host)   &
                   - min( edsclrm( iz_clubb, 1 ), drp_before_clubb( iz_clubb )    &
                      * ( rcm_before_clubb(iz_clubb) -rcm(iz_clubb) ) &
                      /rcm_before_clubb( iz_clubb ) )

                   edsclrm( iz_clubb, 1 ) = edsclrm( iz_clubb, 1 ) &
                                                 - min( edsclrm( iz_clubb, 1 ), drp_before_clubb( iz_clubb )    &
                                                 * ( rcm_before_clubb(iz_clubb) -rcm(iz_clubb) ) &
                                                 /rcm_before_clubb( iz_clubb ) )

             endif

             if( rcm(iz_clubb) <= 0.0 .or. cloud_frac(iz_clubb) <=cloud_frac_min ) then
! in clear areas, droplet concentration should be 0
                 drp_evap_clubb( ix, iy,  iz_host) = &
                   drp_evap_clubb( ix, iy,  iz_host)   &
                  - edsclrm( iz_clubb, 1 )
 
                  edsclrm( iz_clubb, 1 ) = 0.0

             endif
          enddo
        endif !endif do_drp_evap_cf
      endif !use_sclr_HOC
     endif !nqn > 0


! for ice crystal evaporation due to macro-physics
     if(nqni>0)  then 
      if(use_sclr_HOC) then
       if(do_drp_evap_cf) then 
          do iz_clubb=2, kdim+1
             iz_host = kdim + 2 - iz_clubb
             if( ( cloud_frac_before_clubb(iz_clubb) - cloud_frac(iz_clubb) ) > 0.0 .and. cloud_frac(iz_clubb)>0.0) then
! ice crystal evaporation
                 ice_evap_clubb( ix, iy,  iz_host) = &
                    ice_evap_clubb( ix, iy,  iz_host)   &
                   - min( sclrm( iz_clubb, 2 ), ice_before_clubb( iz_clubb )    &
                      * ( cloud_frac_before_clubb(iz_clubb) -cloud_frac(iz_clubb) ) &
                      /cloud_frac_before_clubb( iz_clubb ) )

                 sclrm( iz_clubb, 2 ) = sclrm( iz_clubb, 2 ) &
                                    - min( sclrm( iz_clubb, 2 ), ice_before_clubb( iz_clubb )    &
                                       * ( cloud_frac_before_clubb(iz_clubb) -cloud_frac(iz_clubb) ) &
                                      /cloud_frac_before_clubb( iz_clubb ) )

             endif

             if( rcm(iz_clubb) <= 0.0 .or. cloud_frac(iz_clubb) <=cloud_frac_min ) then
! in clear areas, ice crystal concentration should be 0
                 ice_evap_clubb( ix, iy,  iz_host) = &
                   ice_evap_clubb( ix, iy,  iz_host)   &
                  - sclrm( iz_clubb, 2 )

                 sclrm( iz_clubb, 2 ) = 0.0

             endif
          enddo !enddo iz_clubb
       else !NOT do_drp_evap_cf
          do iz_clubb=2, kdim+1
             iz_host = kdim + 2 - iz_clubb
             if( ( rcm_before_clubb(iz_clubb) - rcm(iz_clubb) ) > 0.0 .and. rcm(iz_clubb)>0.0) then
! ice crystal evaporation
                 ice_evap_clubb( ix, iy,  iz_host) = &
                    ice_evap_clubb( ix, iy,  iz_host)   &
                 - min( sclrm( iz_clubb, 2 ), ice_before_clubb( iz_clubb )    &
                  * ( rcm_before_clubb(iz_clubb) -rcm(iz_clubb) ) &
                   /rcm_before_clubb( iz_clubb ) )

                 sclrm( iz_clubb, 2 ) = sclrm( iz_clubb, 2 ) &
                                    - min( sclrm( iz_clubb, 2 ), ice_before_clubb( iz_clubb )    &
                                       * ( rcm_before_clubb(iz_clubb) -rcm(iz_clubb) ) &
                                        /rcm_before_clubb( iz_clubb ) )
             endif

             if( rcm(iz_clubb) <= 0.0 .or. cloud_frac(iz_clubb) <=cloud_frac_min ) then
! in clear areas, ice crystal concentration should be 0
                 ice_evap_clubb( ix, iy,  iz_host) = &
                   ice_evap_clubb( ix, iy,  iz_host)   &
                  - sclrm( iz_clubb, 2 )

                 sclrm( iz_clubb, 2 ) = 0.0
             endif
          enddo !enddo iz_clubb
       endif !endif do_drp_evap_cf

     else ! not use scalar HOC

       if(do_drp_evap_cf) then 
          do iz_clubb=2, kdim+1
            iz_host = kdim + 2 - iz_clubb
            if( ( cloud_frac_before_clubb(iz_clubb) - cloud_frac(iz_clubb) ) > 0.0 .and. cloud_frac(iz_clubb)>0.0) then
! ice crystal evaporation
                 ice_evap_clubb( ix, iy,  iz_host) = &
                    ice_evap_clubb( ix, iy,  iz_host)   &
                   - min( edsclrm( iz_clubb, 2 ), ice_before_clubb( iz_clubb )    &
                      * ( cloud_frac_before_clubb(iz_clubb) -cloud_frac(iz_clubb) ) &
                      /cloud_frac_before_clubb( iz_clubb ) )

                 edsclrm( iz_clubb, 2 ) = edsclrm( iz_clubb, 2 ) &
                                    - min( edsclrm( iz_clubb, 2 ), ice_before_clubb( iz_clubb )    &
                                       * ( cloud_frac_before_clubb(iz_clubb) -cloud_frac(iz_clubb) ) &
                                      /cloud_frac_before_clubb( iz_clubb ) )
            endif

            if( rcm(iz_clubb) <= 0.0 .or. cloud_frac(iz_clubb) <=cloud_frac_min ) then
! in clear areas, ice crystal concentration should be 0
               ice_evap_clubb( ix, iy,  iz_host) = &
                   ice_evap_clubb( ix, iy,  iz_host)   &
                  - edsclrm( iz_clubb, 2 )

               edsclrm( iz_clubb, 2 ) = 0.0
             endif
          enddo !enddo iz_clubb
        else !NOT do_drp_evap_cf
          do iz_clubb=2, kdim+1
             iz_host = kdim + 2 - iz_clubb
             if( ( rcm_before_clubb(iz_clubb) - rcm(iz_clubb) ) > 0.0 .and. rcm(iz_clubb)>0.0) then
! ice crystal evaporation
                 ice_evap_clubb( ix, iy,  iz_host) = &
                    ice_evap_clubb( ix, iy,  iz_host)   &
                                   - min( edsclrm( iz_clubb, 2 ), ice_before_clubb( iz_clubb )    &
                                      * ( rcm_before_clubb(iz_clubb) -rcm(iz_clubb) ) &
                                      /rcm_before_clubb( iz_clubb ) )

                 edsclrm( iz_clubb, 2 ) = edsclrm( iz_clubb, 2 ) &
                                   - min( edsclrm( iz_clubb, 2 ), ice_before_clubb( iz_clubb )    &
                                      * ( rcm_before_clubb(iz_clubb) -rcm(iz_clubb) ) &
                                      /rcm_before_clubb( iz_clubb ) )

             endif

             if( rcm(iz_clubb) <= 0.0 .or. cloud_frac(iz_clubb) <=cloud_frac_min ) then
! in clear areas, ice crystal concentration should be 0
                 ice_evap_clubb( ix, iy,  iz_host) = &
                   ice_evap_clubb( ix, iy,  iz_host)   &
                  - edsclrm( iz_clubb, 2 )

                 edsclrm( iz_clubb, 2 ) = 0.0
             endif
          enddo
        endif !endif do_drp_evap_cf
      endif !use_sclr_HOC
     endif !nqni > 0

!<-- Note: end of drop and ice crystal evaporation code

      end do ! isub_clubb=1,nsub_clubb, end of clubb calculation
    !  write(*,'(a, 2i5,  6f16.8, i5)') 'time_loop ',  ix, iy,  &
    !          lon(ix,iy)*RAD_TO_DEG,  lat(ix,iy)*RAD_TO_DEG

      ! Option to check for unphysical temperature values
      if (icheck_temp .gt. 0 ) then
        do k=1,kdim+1
          if (temp_clubb(k).lt.-150.0+273.15 .or. temp_clubb(k).gt.90+273.15) then
            write(*,'(a,f6.2,3i4,2f16.8)') 'clubb_driver: bad temperature ',   &
             temp_clubb(k),ix,iy,k,  &
             lon(ix,iy),lat(ix,iy)
          endif
        enddo
      endif
 
      ! **********************************************************
      ! calculate u component tendencies
      call tndy_clubb (dtmain, um, u( ix, iy, : ),     & ! Intent(in) 
                       udt( ix, iy, : ),               & ! Intent(inout) 
                       uten_clubb ( ix, iy, :) )         ! Intent(out) 

      ! calculate v component tendencies
      call tndy_clubb (dtmain, vm, v( ix, iy, : ),  & ! Intent(in) 
                       vdt( ix, iy, : ),            & ! Intent(inout) 
                       vten_clubb ( ix, iy, :) )      ! Intent(out) 

      ! calculate temperature tendencies
      tten_clubb ( ix, iy, :) = 0.0       
      call tndy_clubb (dtmain, temp_clubb, t( ix, iy, : ),  & ! Intent(in)
                       tdt( ix, iy, : ),                    & ! Intent(inout)
                       tten_clubb ( ix, iy, :) )              ! Intent(out)

      ! tendency for droplet evaporation
      if ( nqn > 0 ) then
        drp_evap_clubb = drp_evap_clubb/dtmain
        sclrm_trsport_only(:,1) = sclrm_trsport_only(:,1)/dtmain

        ! calculate droplet number tendency due to pure transport in clubb
        ! convert from clubb coordinate to host SCM coordinate
        do iz_clubb=2, kdim+1
          iz_host = kdim + 2 - iz_clubb
          qndt_clubb_trsport_only( ix, iy,  iz_host) = sclrm_trsport_only( iz_clubb, 1)
        enddo
      endif ! nqn > 0

      ! tendency for ice crystal evaporation
      if ( nqni > 0 ) then
        ice_evap_clubb   = ice_evap_clubb/dtmain
        sclrm_trsport_only(:,2) = sclrm_trsport_only(:,2)/dtmain

        ! calculate ice number tendency due to pure transport in clubb
        ! convert from clubb coordinate to host SCM coordinate
        do iz_clubb=2, kdim+1
          iz_host = kdim + 2 - iz_clubb
          qnidt_clubb_trsport_only( ix, iy,  iz_host) = sclrm_trsport_only( iz_clubb, 2)
        enddo
      endif ! nqni > 0

      ! calculate tracer tendencies
      do it = 1, ntp

        if( it == nsphum )  trs_clubb =  rtm - rcm
        if( it == nql )     trs_clubb =  rcm
        if( it == nqa )     trs_clubb =  cloud_frac

        ! passive tracer evaporation during transport 
        if ( nqn > 0 ) then
          if ( it == nqn )   then
            if ( use_sclr_HOC ) then
              trs_clubb( : ) = sclrm( : , 1 )
            else
              trs_clubb( : ) = edsclrm( : , 1 )
            endif ! use_sclr_HOC
            rdt( ix , iy , : , it) &
            = rdt( ix , iy , : , it) + drp_evap_clubb( ix , iy , : )
          endif ! it == nqn
        endif ! nqn > 0

        if ( nqni > 0 ) then
          if ( it == nqni )   then
            if ( use_sclr_HOC ) then
              trs_clubb( : ) = sclrm( : , 2 )
            else
              trs_clubb( : ) = edsclrm( : , 2 )
            endif ! use_sclr_HOC
            rdt( ix , iy , : , it) &
            = rdt( ix , iy , : , it) + ice_evap_clubb( ix , iy , : )
          endif ! it == nqni
        endif ! nqni > 0

        ! water vapor (nsphum), cloud fraction (nqa), droplet number concentration (nqn), 
        ! ice number concentration (nqni) tendendy in clubb is handled here
        rdt_clubb ( ix, iy, : , it) = 0.0
        if ( it == nqa ) then
          call tndy_clubb (dtmain, trs_clubb, r( ix, iy, : , it),    & ! Intent(in)
                           rdt( ix, iy, : , it),                     & ! Intent(inout)
                           rdt_clubb ( ix, iy, : , it) )               ! Intent(out)
        endif

        if ( it==nsphum .or. it==nqn .or. it==nqni ) then
          call tndy_clubb (dtmain, trs_clubb, r( ix, iy, : , it),    & ! Intent(in)
                           rdt( ix, iy, : , it),                     & ! Intent(inout)
                           rdt_clubb ( ix, iy, : , it) )               ! Intent(out)
          ! ---> h1g, scaled by enviromental fraction (1.0 - convective_fraction)
          rdt_clubb ( ix, iy, : , it) = rdt_clubb ( ix, iy, : , it)     &
                                        *env_condensate_scale(ix, iy, :)

          rdt( ix, iy, : , it)        = rdt_orig( ix, iy, : , it)   &
                                        +rdt_clubb ( ix, iy, : , it)
          ! <--- h1g, scaled by enviromental fraction (1.0 - convective_fraction)
        endif

        if ( it==nql ) then
! ---> h1g, 2012-06-14
         if(  do_liquid_only_in_clubb ) then
          tmp_host = rdt( ix, iy, : , nql)
          call tndy_clubb (dtmain, trs_clubb, & 
                           r( ix, iy, : , nql),                         & ! Intent(in)
                           tmp_host,                                    & ! Intent(inout)
                           rdt_clubb ( ix, iy, : , nql) )                 ! Intent(out)
          rdt_clubb ( ix, iy, : , nql) = rdt_clubb ( ix, iy, : , nql)    &
                                         *env_condensate_scale(ix, iy, :)
          rdt( ix, iy, : , nql) = rdt_orig( ix, iy, : , nql)             &
                                  +rdt_clubb ( ix, iy, : , nql)
         else
          tmp_host = rdt( ix, iy, : , nql) + rdt( ix, iy, : , nqi)
          call tndy_clubb (dtmain, trs_clubb, & 
                           r( ix, iy, : , nql) + r( ix, iy, : , nqi),   & ! Intent(in)
                           tmp_host,                                    & ! Intent(inout)
                           rdt_clubb ( ix, iy, : , nql) )                 ! Intent(out)
          rdt( ix, iy, : , nql) = tmp_host - rdt( ix, iy, : , nqi)
          rdt_clubb ( ix, iy, : , nql) = rdt_clubb ( ix, iy, : , nql)    &
                                         *env_condensate_scale(ix, iy, :)
          rdt( ix, iy, : , nql) = rdt_orig( ix, iy, : , nql)             &
                                  +rdt_clubb ( ix, iy, : , nql)
         endif
! <--- h1g, 2012-06-14
        endif  ! it==nql

      enddo ! do it = 1, ntp

      ! compute average eddy diffusivity cefficients and copy to host grid
      diff_t_1d = diff_t_1d / nsub_clubb
      call clubb2host_full ( diff_t_1d,                                  & ! Intent(in)
                             diff_t_clubb(ix, iy, :))   ! Intent(out)
      ! compute qcvar
      cf_avg   = cf_avg   / nsub_clubb
      rcm_avg  = rcm_avg  / nsub_clubb
      rcp2_avg = rcp2_avg / nsub_clubb

      qcvar_clubb_avg = qcvar_clubb_avg / nsub_clubb

      if ( .not. use_avg4qcvar ) then
        cf_avg   = cloud_frac
        rcm_avg  = rcm
        rcp2_avg = rcp2_zt
      endif

      do iz_clubb = 2, kdim+1
         if ( use_qcvar_in_cloud ) then
           if ( rcm_avg(iz_clubb)*rcm_avg(iz_clubb)/qcvar_clubb_max < &
                (cf_avg(iz_clubb)*rcp2_avg(iz_clubb)-(1.0-cf_avg(iz_clubb))*rcm_avg(iz_clubb)*rcm_avg(iz_clubb)) ) then

                qcvar_clubb_1d(iz_clubb) = rcm_avg(iz_clubb)*rcm_avg(iz_clubb) &
                  /(cf_avg(iz_clubb)*rcp2_avg(iz_clubb)-(1.0-cf_avg(iz_clubb))*rcm_avg(iz_clubb)*rcm_avg(iz_clubb))
                qcvar_clubb_1d(iz_clubb) = qcvar_clubb_1d(iz_clubb) * qcvar_factor_const
                qcvar_clubb_1d(iz_clubb) = max(qcvar_clubb_1d(iz_clubb), qcvar_clubb_min)
           else
             if( cf_avg(iz_clubb)*rcp2_avg(iz_clubb)-(1.0-cf_avg(iz_clubb))*rcm_avg(iz_clubb)*rcm_avg(iz_clubb) <= 0.0 ) then
                qcvar_clubb_1d(iz_clubb) = qcvar_clubb_min
                if( do_scale_negative_qcvar ) then
                   if( cf_avg(iz_clubb)*rcp2_avg(iz_clubb)-(1.0-cf_avg(iz_clubb))*rcm_avg(iz_clubb)*rcm_avg(iz_clubb) < 0.0 & 
                       .and. rcm(iz_clubb) >= rcm_min ) then

                     qcvar_clubb_1d(iz_clubb) = qcvar_clubb_1d(iz_clubb)                   &
             * exp( qcvar_factor*rcm_avg(iz_clubb)*rcm_avg(iz_clubb)                              &
     /(cf_avg(iz_clubb)*rcp2_avg(iz_clubb)-(1.0-cf_avg(iz_clubb))*rcm_avg(iz_clubb)*rcm_avg(iz_clubb)) )

                     qcvar_clubb_1d(iz_clubb) =  max(qcvar_clubb_1d(iz_clubb), qcvar_clubb_min*exp(-scale_negative_qcvar_factor))


  !                    write(*,'(a,12(f7.3,1x))') 'Neg. qcvar:',                                            &
  !                        lon(ix,iy)*RAD_TO_DEG,  lat(ix,iy)*RAD_TO_DEG,                                    &
  !                                             thermodynamic_heights(iz_clubb)*0.001,                       &
  !                      rcm_avg(iz_clubb)*1.e5,  cf_avg(iz_clubb)*100.0,   rcp2_avg(iz_clubb)*1.e10,        &
  !                        rcm_avg(iz_clubb)*rcm_avg(iz_clubb)                                               &
  !   /(cf_avg(iz_clubb)*rcp2_avg(iz_clubb)-(1.0-cf_avg(iz_clubb))*rcm_avg(iz_clubb)*rcm_avg(iz_clubb)),     &
  !         qcvar_clubb_min * exp( qcvar_factor*rcm_avg(iz_clubb)*rcm_avg(iz_clubb)                          &
  !   /(cf_avg(iz_clubb)*rcp2_avg(iz_clubb)-(1.0-cf_avg(iz_clubb))*rcm_avg(iz_clubb)*rcm_avg(iz_clubb)) ),   &
  !           qcvar_clubb_min/exp(-scale_negative_qcvar_factor),                                                   &
  !                   qcvar_clubb_1d(iz_clubb)

                   endif
                endif
             endif
           endif
         else
           if ( rcm_avg(iz_clubb)*rcm_avg(iz_clubb)/qcvar_clubb_max < rcp2_avg(iz_clubb) ) then
              qcvar_clubb_1d(iz_clubb) = rcm_avg(iz_clubb)*rcm_avg(iz_clubb)/rcp2_avg(iz_clubb)
              qcvar_clubb_1d(iz_clubb) = qcvar_clubb_1d(iz_clubb) * qcvar_factor_const
              qcvar_clubb_1d(iz_clubb) = max(qcvar_clubb_1d(iz_clubb), qcvar_clubb_min)
           endif
         endif
         qcvar_cf_clubb_1d(iz_clubb) = qcvar_clubb_1d(iz_clubb) * cf_avg(iz_clubb)
      enddo

      if ( use_avg_qcvar ) then
         qcvar_clubb_1d    = qcvar_clubb_avg
         qcvar_cf_clubb_1d = qcvar_clubb_1d * cf_avg
      endif 

      if( present( qcvar_clubb ) ) &
      call clubb2host_full ( qcvar_clubb_1d,        & ! Intent(in)
                             qcvar_clubb(ix, iy, :))   ! Intent(out) 

      call clubb2host_full ( rcm,                  & ! Intent(in)
                             rcm_clubb(ix, iy, :))   ! Intent(out) 

      call clubb2host_full ( cloud_frac,           & ! Intent(in)
                             cf_clubb(ix, iy, :))   ! Intent(out)

      call clubb2host_full ( cf_avg,               & ! Intent(in)
                             cf_avg_clubb_3d(ix, iy, :))   ! Intent(out)

      call clubb2host_full ( qcvar_cf_clubb_1d,        & ! Intent(in)
                             qcvar_cf_clubb_3d(ix, iy, :))   ! Intent(out)

      ! dump higher-order terms
      call clubb2host_full ( rcp2_zt,              & ! Intent(in)
                             rcp2_clubb(ix, iy, :))  ! Intent(out) 

      call clubb2host_full ( wp3,                                        & ! Intent(in)
                             wp3_clubb(ix, iy, :))   ! Intent(out)

      call clubb2host_half ( wp2,                                        & ! Intent(in)
                             wp2_clubb(ix, iy, :))   ! Intent(out)

      call clubb2host_half ( upwp,                                        & ! Intent(in)
                             upwp_clubb(ix, iy, :))   ! Intent(out)

      call clubb2host_half ( vpwp,                                        & ! Intent(in)
                             vpwp_clubb(ix, iy, :))   ! Intent(out)

      call clubb2host_half ( up2,                                         & ! Intent(in)
                             up2_clubb(ix, iy, :))    ! Intent(out)
  
      call clubb2host_half ( vp2,                                         & ! Intent(in)
                             vp2_clubb(ix, iy, :))    ! Intent(out)

      call clubb2host_half ( wprtp,                                       & ! Intent(in)
                             wprtp_clubb(ix, iy, :))  ! Intent(out)
 
      call clubb2host_half ( wpthlp,                                      & ! Intent(in)
                             wpthlp_clubb(ix, iy,:))  ! Intent(out)
 
      call clubb2host_half ( rtp2,                                        & ! Intent(in)
                             rtp2_clubb(ix, iy, :))   ! Intent(out)
       
      call clubb2host_half ( thlp2,                                       & ! Intent(in)
                             thlp2_clubb(ix, iy, :))  ! Intent(out)
 
      call clubb2host_half ( rtpthlp,                                      & ! Intent(in)
                             rtpthlp_clubb(ix, iy,:))  ! Intent(out)

      call clubb2host_full ( wm_zt,                                      & ! Intent(in)
                             wm_clubb(ix, iy,:))     ! Intent(out)

      ! calculate some vertically averaged terms, such as wp2_vert_avg_clubb
      tmp0D = 0.0
      do iz_clubb = 2, kdim
         wp2_vert_avg_clubb(ix, iy)  = wp2_vert_avg_clubb(ix,iy) + wp2(iz_clubb)*rho_zm(iz_clubb)/gr%invrs_dzm(iz_clubb)
         up2_vert_avg_clubb(ix, iy)  = up2_vert_avg_clubb(ix,iy) + up2(iz_clubb)*rho_zm(iz_clubb)/gr%invrs_dzm(iz_clubb)
         vp2_vert_avg_clubb(ix, iy)  = vp2_vert_avg_clubb(ix,iy) + vp2(iz_clubb)*rho_zm(iz_clubb)/gr%invrs_dzm(iz_clubb)
         rtp2_vert_avg_clubb(ix, iy) = rtp2_vert_avg_clubb(ix,iy) + rtp2(iz_clubb)*rho_zm(iz_clubb)/gr%invrs_dzm(iz_clubb)
         thlp2_vert_avg_clubb(ix, iy)= thlp2_vert_avg_clubb(ix,iy) + thlp2(iz_clubb)*rho_zm(iz_clubb)/gr%invrs_dzm(iz_clubb)

         tmp0D = tmp0D + rho_zm(iz_clubb)/gr%invrs_dzm(iz_clubb)
      enddo
      wp2_vert_avg_clubb(ix, iy)   = wp2_vert_avg_clubb(ix,iy)/tmp0D
      up2_vert_avg_clubb(ix, iy)   = up2_vert_avg_clubb(ix,iy)/tmp0D
      vp2_vert_avg_clubb(ix, iy)   = vp2_vert_avg_clubb(ix,iy)/tmp0D
      rtp2_vert_avg_clubb(ix, iy)  = rtp2_vert_avg_clubb(ix,iy)/tmp0D
      thlp2_vert_avg_clubb(ix, iy) = thlp2_vert_avg_clubb(ix,iy)/tmp0D

      ! Aerosol activation code (?)
      if ( nqn > 0 )  then
        ! convert aerosol mass from host from clubb levels
        do i_totmass = 1, n_totmass
          call host2clubb_full( totalmass1(ix, iy,:, i_totmass), &    !intent (in)
                              aeromass_clubb (:, i_totmass) )         !intent (out)
        enddo  ! i_totmass
       
        do i_imass = 1, n_imass
          call host2clubb_full( imass1(ix, iy,:, i_imass), &    !intent (in)
                              aeroimass_clubb (:, i_imass) )         !intent (out)
        enddo  

        if (do_aeromass_clubb_const) then
          aeromass_clubb (:, :) = aeromass_clubb_const
        end if
 
        if ( do_ice_nucl_wpdf ) THEN
          ! have both warm and ice nucleation 
          call  warm_ice_nucl_mns_clubb( aeromass_clubb,     &    !  Intent(in)
                                         aeroimass_clubb,    &    !  Intent(in)
                                         temp_clubb,         &    !  Intent(in)
                                         Ndrop_max,          &    !  Intent(out)
                                         Ncrystal_max,       &    !  Intent(out)
                                         RH_crit  )               !  Intent(out)
        else
          Ncrystal_max = 0.0
          RH_crit = 1.0
  
          if ( do_quadrature_gauss )   &
            call aer_act_clubb_quadrature_Gauss(Time_next, aeromass_clubb, temp_clubb, &   ! Intent(in)
                                                Ndrop_max)                                 ! Intent(out)
  
          if ( do_diffK_gauss )   &
            call aer_act_clubb_diffK_Gauss( aeromass_clubb, temp_clubb,                &   ! Intent(in)
                                            Ndrop_max )                                    ! Intent(out)
  
          if ( do_BL_gauss )   &
            call  aer_act_clubb_BL_Gauss( aeromass_clubb,        &   ! Intent(in)
                                          Ndrop_max )                ! Intent(out)

        endif ! do_ice_nucl_wpdf

        ! convert the unit from #/cm3 to #/kg air 
        Ndrop_max = Ndrop_max * 1.0e6 / rho
   
        ! update tracer concentrations: droplet # concentration
        do iz_clubb = 2,kdim+1

          ! h1g, 2011-04-20,  no liquid drop nucleation if T < -40 C
          if ( temp_clubb(iz_clubb) <= 233.15 )  Ndrop_max( iz_clubb ) = 0.0  ! if T<-40C, no liquid drop nucleation

          if ( rcm(iz_clubb) > 0.0 .and. cloud_frac(iz_clubb) > cloud_frac_min ) then
            if ( use_sclr_HOC ) then
              sclrm( iz_clubb, 1 ) = max(Ndrop_max( iz_clubb ), sclrm( iz_clubb, 1 ) )
            else
              edsclrm( iz_clubb, 1 ) = max(Ndrop_max( iz_clubb ), edsclrm( iz_clubb, 1 ) )
            endif
          endif

        enddo

        do iz_clubb = 2,kdim+1
          iz_host = kdim + 2 - iz_clubb
          if ( rcm(iz_clubb) > 0.0 .and. cloud_frac(iz_clubb) > cloud_frac_min ) then
            ! Ndrop_act_clubb: in-cloud averaged activated droplet number concentration (in the unit of #/kg)
            Ndrop_act_clubb(ix, iy, iz_host) = Ndrop_max(iz_clubb)/cloud_frac(iz_clubb)
          endif
        enddo

        if ( use_sclr_HOC ) then
          trs_clubb(:) =  sclrm(:,1)
        else
          trs_clubb(:) =  edsclrm(:,1)
        endif

        call tndy_clubb (dtmain, trs_clubb, r(ix, iy,:, nqn),  & ! Intent(in)
                         rdt( ix, iy, : , nqn),                & ! Intent(inout)
                         aer_ccn_act_clubb( ix, iy,:))           ! Intent(out)

        ! calculate drop concentration net flux
        drp_flux_clubb(  ix, iy, 1 ) = 0.0
        do iz_host = 2, kdim+1
          !-->cjg: bug fix?
          ! iz_clubb = (kdim+1 - iz_host) + 1 + 1
          iz_clubb = kdim + 2 - iz_host
          !<--cjg: bug fix

          if ( iz_clubb .gt. kdim ) then
            drp_flux_clubb(  ix, iy,  iz_host ) = 0.0
          else
            sum_sclrm_trsport = 0.0
            !--> cjg: bizarre loop
            do isum_sclrm = iz_clubb+1,iz_clubb+1
              sum_sclrm_trsport = sum_sclrm_trsport + sclrm_trsport_only( isum_sclrm , 1)
            enddo ! isum_sclrm
 
            drp_flux_clubb( ix, iy, iz_host ) =           &
                            drp_flux_clubb(  ix, iy, iz_host-1 )   &
                            + sum_sclrm_trsport                                               &
                            *( zhalf( ix, iy, iz_host-1) - zhalf( ix, iy, iz_host) )
          endif ! iz_clubb  .gt.  kdim
        enddo ! iz_host

        if ( do_ice_nucl_wpdf ) THEN
          ! convert the unit from #/cm3 to #/kg air for ice crystal number concentration
          Ncrystal_max =  Ncrystal_max / rho
    
          do iz_clubb = 2,kdim+1
            if ( rcm(iz_clubb) > 0.0 .and. cloud_frac(iz_clubb) > cloud_frac_min ) then
              if ( use_sclr_HOC ) then
                sclrm( iz_clubb, 2 ) = max(Ncrystal_max( iz_clubb ), sclrm( iz_clubb, 2 ) )
              else
                edsclrm( iz_clubb, 2 ) = max(Ncrystal_max( iz_clubb ), edsclrm( iz_clubb, 2 ) )
              endif
            endif
          enddo

          !   **********  in current version   ********** 
          do iz_clubb = 2,kdim+1
            iz_host = kdim + 2 - iz_clubb
            if (rcm(iz_clubb) > 0.0 .and. cloud_frac(iz_clubb) > cloud_frac_min) then
              ! Icedrop_act_clubb: in-cloud  averaged nucleated ice-crystal number concentration (in the unit of #/kg)  
              Icedrop_act_clubb( ix, iy, iz_host ) = Ncrystal_max( iz_clubb )/cloud_frac(iz_clubb) 
            endif
          enddo

          if ( use_sclr_HOC ) then
            trs_clubb( : ) =  sclrm( : , 2 )
          else
            trs_clubb( : ) =  edsclrm( : , 2 )
          endif

          aer_ice_act_clubb( ix, iy, : )  = 0.0
          
          !call tndy_clubb (dtmain, trs_clubb, r( ix, iy, : , nqni),    & ! Intent(in)
          !                 rdt( ix, iy, : , nqni),                     & ! Intent(inout)
          !                 aer_ice_act_clubb( ix, iy, : ) )              ! Intent(out)

          ! calculate ice crystal concentration flux
          Icedrop_flux_clubb(  ix, iy, 1 ) = 0.0
          do  iz_host = 2, kdim+1
            !-->cjg: bug fix?
            ! iz_clubb = (kdim+1 - iz_host) + 1 + 1
            iz_clubb = kdim + 2 - iz_host
            !<--cjg: bug fix

            if ( iz_clubb .gt. kdim ) then
              Icedrop_flux_clubb(  ix, iy,  iz_host ) = 0.0
            else
              sum_sclrm_trsport = 0.0
              !--> cjg: bizarre loop
              do isum_sclrm = iz_clubb+1,iz_clubb+1
                sum_sclrm_trsport = sum_sclrm_trsport + sclrm_trsport_only( isum_sclrm , 2)
              enddo ! isum_sclrm
 
              Icedrop_flux_clubb(  ix, iy, iz_host ) =           &
                            Icedrop_flux_clubb(  ix, iy, iz_host-1 )   &
                           + sum_sclrm_trsport                                                     &
                           *( zhalf( ix, iy, iz_host-1) - zhalf( ix, iy, iz_host) )
            endif ! iz_clubb .gt. kdim
          enddo ! iz_host
        endif ! do_ice_nucl_wpdf

      endif  ! nqn > 0
      ! end passive tracers

      call clubb_1D_2_3D( ix, iy, rdiag )

    enddo ! enddo of ix
  enddo ! enddo of iy

  if (do_alt_cloud > 0) then

    call alt_cloud( do_alt_cloud,dtmain,pfull,  &
                    t,q,r(:,:,:,nqa),r(:,:,:,nql),r(:,:,:,nqi),  &
                    tdt,qdt,rdt(:,:,:,nqa),rdt(:,:,:,nql),rdt(:,:,:,nqi) )

  endif

  ! ----- diagnostics -----
  if (id_udt_clubb > 0) then
    used = send_data( id_udt_clubb, uten_clubb,  &
                      Time_next, is, js, 1 )
  endif

  if (id_vdt_clubb > 0) then
    used = send_data( id_vdt_clubb, vten_clubb,  &
                      Time_next, is, js, 1 )
  endif

  do it = 1, ntp
    if ( id_tracer_clubb(it)>0 ) then
      used = send_data( id_tracer_clubb(it), rdt_clubb ( : , : , : , it),  &
                        Time_next, is, js, 1 )
    endif 
  enddo

  if ( do_clubb_conservation_checks )  then
    do iz_host = 1, kdim
      do iy = 1, jdim
        do ix = 1, idim
          pmass_3d( ix, iy,  iz_host) = (phalf( ix, iy,  iz_host+1) &
                                        - phalf( ix, iy, iz_host) )/grav
        enddo
      enddo
    enddo


! ---> h1g, 2012-08-23
      if ( .not. l_host_applies_sfc_fluxes )  then
         tmp2d_check(:, :) = tmp2d_check(:, :) / nsub_clubb
      endif
! <--- h1g, 2012-08-23

    do iz_clubb= 2, kdim+1
       iz_host = kdim + 2 - iz_clubb
       do iy = 1, jdim
          do ix = 1, idim
            tmp2d_check(  ix, iy ) =   tmp2d_check(  ix, iy ) &
                      +   pmass_3d(ix, iy,  iz_host) &
                         *(   tten_clubb(  ix, iy,  iz_host )*Cp &
                             - rdt_clubb(  ix, iy,  iz_host , nql)*Lv )
          enddo
       enddo
    enddo

! ---> h1g    Very slightly adjust tendencies to force exact   ***
!                     enthalpy  conservation, 2010-06-15
    do iy = 1, jdim
      do ix = 1, idim
       ! temperature fix is only applied for the clubb vertical domain
         tten_clubb_fix =  tmp2d_check( ix, iy ) & 
                          / ( phalf(ix, iy, kdim+1) - phalf(ix, iy, 1)) &
                          * grav / Cp

         tmp2d_check(ix, iy) = 0.0

       ! recalculate column enthalpy tendency
         do  iz_clubb= 2, kdim+1
             iz_host = kdim + 2 - iz_clubb
             tten_clubb(ix, iy, iz_host) = tten_clubb(ix, iy, iz_host) - tten_clubb_fix

            ! add temperature fix to the temperature tendency  "tdt"
             tdt(ix, iy, iz_host) = tdt(ix, iy, iz_host) - tten_clubb_fix

             tmp2d_check(ix, iy) = tmp2d_check(ix, iy) &
                                  +pmass_3d(ix, iy,  iz_host) &
                                  *( tten_clubb( ix, iy, iz_host)*Cp &
                                   - rdt_clubb( ix, iy, iz_host, nql)*Lv)
         enddo
      enddo
    enddo
! <--- h1g,  2010-06-15 
    if (id_enth_clubb_col > 0 ) then
      used = send_data(id_enth_clubb_col, tmp2d_check, Time_next, is, js)
    endif

    tmp2d_check = 0.0
    if ( do_clubb_diss_heat ) then
     do iz_clubb= 2, kdim+1
      iz_host = kdim + 2 - iz_clubb
      do iy = 1, jdim
        do ix = 1, idim
         uv_dissipation =  (u(ix, iy, iz_host) + (udt(ix, iy, iz_host) - 0.5*uten_clubb(ix, iy, iz_host))*dtmain ) &
                          *uten_clubb(ix, iy, iz_host) &
                         + (v(ix, iy, iz_host) + (vdt(ix, iy, iz_host) - 0.5*vten_clubb(ix, iy, iz_host))*dtmain ) &
                          *vten_clubb(ix, iy, iz_host)

         tten_clubb(ix, iy, iz_host) = tten_clubb(ix, iy, iz_host) - uv_dissipation/Cp
         tdt       (ix, iy, iz_host) = tdt       (ix, iy, iz_host) - uv_dissipation/Cp

         tmp2d_check(ix, iy) = tmp2d_check(ix, iy) - pmass_3d(ix, iy, iz_host) * uv_dissipation
        enddo
      enddo
     enddo
    endif  ! do_clubb_diss_heat
    if ( id_diss_heat_clubb > 0 ) then
         used = send_data (id_diss_heat_clubb, tmp2d_check, Time_next, is, js)
    endif  ! id_diss_heat_clubb > 0

    if ( id_wat_clubb_col > 0 ) then
      tmp2d_check = 0.0
      do iz_host = 1, kdim
        do iy = 1, jdim
          do ix = 1, idim
            tmp2d_check(  ix, iy ) =   tmp2d_check(  ix, iy ) &
                      +   pmass_3d(ix, iy,  iz_host) &
                         *(   rdt_clubb(  ix, iy,  iz_host , nsphum) &
                            + rdt_clubb(  ix, iy,  iz_host , nql) )
          enddo
        enddo
      enddo
      used = send_data ( id_wat_clubb_col, tmp2d_check, Time_next, is, js )
    endif  ! id_wat_clubb_col > 0

    do it = 1, ntp
      tmp2d_check = 0.0
      do iz_host = 1, kdim
        do iy = 1, jdim
          do ix = 1, idim
            tmp2d_check(ix, iy) = tmp2d_check(ix, iy) + pmass_3d(ix, iy, iz_host) * rdt_clubb(ix, iy, iz_host, it)
          enddo
        enddo
      enddo

      if ( id_tracer_clubb_col(it)>0 ) used = send_data (id_tracer_clubb_col(it), tmp2d_check, Time_next, is, js)
    enddo ! it

  endif ! do_clubb_conservation_checks

  if (id_tdt_clubb > 0) then
    used = send_data ( id_tdt_clubb,   tten_clubb, Time_next, is, js, 1 )
  endif

  if ( nqn > 0 )   then
    if ( id_drp_evap_clubb>0 ) then
      used = send_data(id_drp_evap_clubb, drp_evap_clubb, Time_next, is, js, 1)
    endif

    if ( id_aer_ccn_act_clubb > 0) then
      used = send_data(id_aer_ccn_act_clubb, aer_ccn_act_clubb, Time_next, is, js, 1)
    endif

    if ( id_Ndrop_act_clubb > 0) then
      used = send_data(id_Ndrop_act_clubb, Ndrop_act_clubb, Time_next, is, js, 1 )
    endif

    ! calculate droplet number turbulence flux (m/cm3/s)
    if ( id_drp_flux_clubb > 0 ) then
      used = send_data(id_drp_flux_clubb, drp_flux_clubb, Time_next, is, js, 1 )
    endif !  id_drp_flux_clubb > 0

    if ( id_qndt_clubb_trsport_only > 0 ) then
      used = send_data(id_qndt_clubb_trsport_only, qndt_clubb_trsport_only, Time_next, is, js, 1)
    endif
  endif   !  nqn > 0
  
  if ( nqni > 0 )   then
    if ( id_ice_evap_clubb > 0) then
      used = send_data(id_ice_evap_clubb, ice_evap_clubb, Time_next, is, js, 1)
    endif

    if ( id_aer_ice_act_clubb>0 ) then
      used = send_data(id_aer_ice_act_clubb, aer_ice_act_clubb, Time_next, is, js, 1)
    endif

    if ( id_Icedrop_act_clubb>0 ) then
      used = send_data(id_Icedrop_act_clubb, Icedrop_act_clubb, Time_next, is, js, 1)
    endif

    ! calculate droplet number turbulence flux (m/cm3/s)
    if ( id_icedrop_flux_clubb > 0 ) then
      used = send_data(id_icedrop_flux_clubb, icedrop_flux_clubb, Time_next, is, js, 1)
    endif !  id_drp_flux_clubb > 0

    if ( id_qnidt_clubb_trsport_only > 0 ) then
      used = send_data ( id_qnidt_clubb_trsport_only, qnidt_clubb_trsport_only, Time_next, is, js, 1)
    endif
  endif   !  nqni > 0

  ! output higher order terms and fluxes
  ! at full levels
  if ( id_wp3_clubb > 0 ) then
    used = send_data(id_wp3_clubb, wp3_clubb, Time_next, is, js, 1)
  endif

  if ( id_wm_clubb > 0 ) then
    used = send_data(id_wm_clubb, wm_clubb, Time_next, is, js, 1)
  endif

  if ( id_omega_clubb > 0 ) then
    used = send_data(id_omega_clubb, omega, Time_next, is, js, 1)
  endif

  if ( id_qcvar_clubb > 0 ) then
    used = send_data(id_qcvar_clubb, qcvar_clubb, Time_next, is, js, 1)
  endif

  if ( id_qcvar_cf_clubb > 0 ) then
    used = send_data(id_qcvar_cf_clubb, qcvar_cf_clubb_3d, Time_next, is, js, 1)
  endif

  if ( id_rcm_clubb >0 ) then
    used = send_data(id_rcm_clubb, rcm_clubb, Time_next, is, js, 1)
  endif

  if ( id_rcp2_clubb >0 ) then
    used = send_data(id_rcp2_clubb, rcp2_clubb, Time_next, is, js, 1)
  endif

  if ( id_cf_clubb >0 ) then
    used = send_data(id_cf_clubb, cf_clubb, Time_next, is, js, 1)
  endif

  if ( id_cf_avg_clubb >0 ) then
    used = send_data(id_cf_avg_clubb, cf_avg_clubb_3d, Time_next, is, js, 1)
  endif


  ! at half levels
  if ( id_wp2_clubb > 0 ) then
    used = send_data ( id_wp2_clubb, wp2_clubb, Time_next, is, js, 1)
  endif

  if ( id_upwp_clubb > 0 ) then
    used = send_data ( id_upwp_clubb, upwp_clubb, Time_next, is, js, 1)
  endif

  if ( id_vpwp_clubb > 0 ) then
    used = send_data ( id_vpwp_clubb, vpwp_clubb, Time_next, is, js, 1)
  endif

  if ( id_up2_clubb > 0 ) then
    used = send_data ( id_up2_clubb, up2_clubb, Time_next, is, js, 1)
  endif

  if ( id_vp2_clubb > 0 ) then
    used = send_data ( id_vp2_clubb, vp2_clubb, Time_next, is, js, 1)
  endif

  if ( id_wprtp_clubb > 0 ) then
    used = send_data ( id_wprtp_clubb, wprtp_clubb, Time_next, is, js, 1)
  endif

  if ( id_wpthlp_clubb > 0 ) then
    used = send_data ( id_wpthlp_clubb, wpthlp_clubb, Time_next, is, js, 1)
  endif

  if ( id_rtp2_clubb > 0 ) then
    used = send_data ( id_rtp2_clubb, rtp2_clubb, Time_next, is, js, 1)
  endif

  if ( id_thlp2_clubb > 0 ) then
    used = send_data ( id_thlp2_clubb, thlp2_clubb, Time_next, is, js, 1)
  endif

  if ( id_rtpthlp_clubb > 0 ) then
    used = send_data ( id_rtpthlp_clubb, rtpthlp_clubb, Time_next, is, js, 1)
  endif

  ! column vertically interated averages
  if ( id_wp2_vert_avg_clubb > 0 ) then
    used = send_data ( id_wp2_vert_avg_clubb, wp2_vert_avg_clubb, Time_next, is, js)
  endif

  if ( id_up2_vert_avg_clubb > 0 ) then
    used = send_data ( id_up2_vert_avg_clubb, up2_vert_avg_clubb, Time_next, is, js)
  endif

  if ( id_vp2_vert_avg_clubb > 0 ) then
    used = send_data ( id_vp2_vert_avg_clubb, vp2_vert_avg_clubb, Time_next, is, js)
  endif

  if ( id_rtp2_vert_avg_clubb > 0 ) then
    used = send_data ( id_rtp2_vert_avg_clubb, rtp2_vert_avg_clubb, Time_next, is, js)
  endif

  if ( id_thlp2_vert_avg_clubb > 0 ) then
    used = send_data ( id_thlp2_vert_avg_clubb, thlp2_vert_avg_clubb, Time_next, is, js)
  endif

  return
  end subroutine clubb

  !=====================================================================
  !=====================================================================

  subroutine clubb_init(id, jd, kd, lon, lat, axes, Time, phalf, &
                              l_host_applies_sfc_fluxes_in , &
                                   do_ice_nucl_wpdf_in, do_liq_num_in)

  !---------------------------------------------------------------------
  ! Initialize and set-up clubb
  !
  !  id, jd, kd                         input
  !    subdomain dimensions
  !
  !  lon, lat                           input
  !    longitude, latitude in radians
  !
  !  axes                               input
  !    axis indices, (/x,y,pf,ph/)
  !    (returned from diag axis manager)
  !
  !  Time                               input
  !    current time (time_type)
  !
  !  phalf                              input
  !    pressure at half levels in pascals
  !    [real, dimension(nlon,nlat,nlev+1)]
  !
  !---------------------------------------------------------------------

  ! ----- Calling arguments -----

  integer, intent(in)                  :: id, jd, kd
  real, dimension(:,:), intent(in)     :: lon, lat
  integer, dimension(4), intent(in)    :: axes
  type(time_type), intent(in)          :: Time
  real, dimension(:,:,:), intent(in)   :: phalf
  logical,   intent(in)                :: l_host_applies_sfc_fluxes_in
  logical,   intent(in)                :: do_ice_nucl_wpdf_in   
  logical,   intent(in)                :: do_liq_num_in         
 
  ! ----- Local variables -----

  integer, dimension(3) :: half = (/1,2,4/)
  integer :: ierr, io, unit            !open namelist error

  integer :: i, j, n

  integer :: isclr_dim
  character(len=64)            :: fname='INPUT/clubb.res.nc'
  character(len=64)            :: sclr_prefix,  sclr_txt
  integer :: year, month, day, hour, minute, second
  real(kind=time_precision) :: time_current
  type(time_type)           :: Time_init

  ! ----- Begin code -----

  if (module_is_initialized) then
    return
  else
    module_is_initialized = .true.
  endif


  clubb_core_clock = mpp_clock_id( '   clubb_driver: core',    &
                                   grain=CLOCK_MODULE_DRIVER )

  ! ----- Read namelist -----

  if ( file_exist( 'input.nml' ) ) then

    unit = open_namelist_file ( )
    ierr = 1        
    do while( ierr /= 0 )
      read ( unit,  nml = clubb_setting_nml, iostat = io, end = 10 ) 
      ierr = check_nml_error (io, 'clubb_setting_nml')
    end do
10  call close_file( unit )

    unit = open_namelist_file ( )
    ierr = 1
    do while( ierr /= 0 )
     read ( unit,  nml = clubb_stats_setting_nml, iostat = io, end = 20 ) 
     ierr = check_nml_error (io, 'clubb_stats_setting_nml')
    end do
20  call close_file( unit )

  end if

!   Save needed module variables
  do_liq_num = do_liq_num_in
  do_ice_nucl_wpdf = do_ice_nucl_wpdf_in
  l_host_applies_sfc_fluxes = l_host_applies_sfc_fluxes_in


  ! Get tracer indices
  nsphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )
  nql    = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
  nqi    = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
  nqa    = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
  nqn    = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
  nqni   = get_tracer_index ( MODEL_ATMOS, 'ice_num' )

  nupwp    = get_tracer_index ( MODEL_ATMOS, 'upwp' )
  nvpwp    = get_tracer_index ( MODEL_ATMOS, 'vpwp' )
  nup2     = get_tracer_index ( MODEL_ATMOS, 'up2' )
  nvp2     = get_tracer_index ( MODEL_ATMOS, 'vp2' )
  nwp2     = get_tracer_index ( MODEL_ATMOS, 'wp2' )
  nwprtp   = get_tracer_index ( MODEL_ATMOS, 'wprtp' )
  nwpthlp  = get_tracer_index ( MODEL_ATMOS, 'wpthlp' )
  nrtp2    = get_tracer_index ( MODEL_ATMOS, 'rtp2' )
  nthlp2   = get_tracer_index ( MODEL_ATMOS, 'thlp2' )
  nrtpthlp = get_tracer_index ( MODEL_ATMOS, 'rtpthlp' )
  nwp3     = get_tracer_index ( MODEL_ATMOS, 'wp3' )
  nrhcrit1 = get_tracer_index ( MODEL_ATMOS, 'rhcrit1' )
  nrhcrit2 = get_tracer_index ( MODEL_ATMOS, 'rhcrit2' )

  if( do_liq_num .and. nqn<= 0 )  &
    call error_mesg ('clubb_driver_mod','do_liq_num is true, but nqn<= 0', FATAL)
 
  ierr = 0
  if ( nupwp <= 0 ) ierr = 1
  if ( nvpwp <= 0 ) ierr = 1
  if ( nup2  <= 0 ) ierr = 1
  if ( nvp2  <= 0 ) ierr = 1
  if ( nwp2  <= 0 ) ierr = 1
  if ( nwprtp <= 0 ) ierr = 1
  if ( nwpthlp <= 0 ) ierr = 1
  if ( nrtp2 <= 0 ) ierr = 1
  if ( nthlp2 <= 0 ) ierr = 1
  if ( nrtpthlp <= 0 ) ierr = 1
  if ( nwp3 <= 0 ) ierr = 1
  if ( nrhcrit1 <= 0 ) ierr = 1
  if ( nrhcrit2 <= 0 ) ierr = 1
  if ( ierr > 0 ) then
    call error_mesg ('clubb_driver_mod',                                      &
                     'one or more CLUBB variables not found in field_table',  &
                     FATAL)
  end if

  ! ----- write version number and namelist -----
  call write_version_number( version, tagname )
  if(mpp_pe() == mpp_root_pe() ) write(stdlog(),nml=clubb_setting_nml)

  ! set-up clubb
  call clubb_setup(id, jd, phalf)

  ! initialize aerosol activation modules
  call aer_ccn_act_init
  if ( do_ice_nucl_wpdf ) call ice_nucl_wpdf_init
  if (var_limit_ice .EQ. -999. ) THEN
      do_var_limit_ice = .false.
  ELSE
      do_var_limit_ice = .true.
  END IF

  !set-up diagostics
  call get_number_tracers(MODEL_ATMOS, num_tracers=nt, num_prog=ntp)

  allocate( id_tracer_clubb(nt) )
  allocate( id_tracer_clubb_col(nt) )

  id_udt_clubb = register_diag_field ( mod_name,    &
    'udt_clubb', axes(1:3), Time,                  &
    'u-component tendency from clubb (boundary layer cloud)', &
    'm/s/s', missing_value=missing_value)

  id_vdt_clubb = register_diag_field ( mod_name,    &
    'vdt_clubb', axes(1:3), Time,                   &
    'v-component tendency from clubb (boundary layer cloud)', &
    'm/s/s', missing_value=missing_value)

  id_tdt_clubb = register_diag_field ( mod_name,    &
    'tdt_clubb', axes(1:3), Time,                   &
    'temperature tendency from clubb (boundary layer cloud)', &
    'K/s', missing_value=missing_value)
 
  do it = 1, ntp
    call get_tracer_names(MODEL_ATMOS, it, name = tracer_name,  &
                          units = tracer_units)
        
    diaglname = trim(tracer_name)//  &
               ' tendency from clubb (boundary layer cloud)'
    id_tracer_clubb(it) =    &
        register_diag_field ( mod_name, &
        TRIM(tracer_name)//'dt_clubb',  &
        axes(1:3), Time, trim(diaglname), &
        TRIM(tracer_units)//'/s',  &
        missing_value=missing_value)
  enddo

  do it = 1, ntp
    call get_tracer_names (MODEL_ATMOS, it, name = tracer_name,  &
                           units = tracer_units)

    diaglname = trim(tracer_name)//' column tendency from clubb'  

    id_tracer_clubb_col(it) =    &
                          register_diag_field ( mod_name, &
                          TRIM(tracer_name)//'_clubb_col',  &
                          axes(1:2), Time, trim(diaglname), &
                          TRIM(tracer_units)//'/s',  &
                          missing_value=missing_value) 
  enddo

! --->h1g, 2012-05-29
  id_diss_heat_clubb = register_diag_field ( mod_name, &
    'diss_heat_clubb', axes(1:2), Time, &
    'Column dissipative heat tendency from clubb','W/m2' )
! <---h1g, 2012-05-29

  id_enth_clubb_col = register_diag_field ( mod_name, &
    'enth_clubb_col', axes(1:2), Time, &
    'Column enthalpy tendency from clubb','W/m2' )
 
  id_wat_clubb_col = register_diag_field ( mod_name, &
    'wat_clubb_col', axes(1:2), Time, &
    'Column total water tendency from clubb','kg/m2/s' )

  id_drp_evap_clubb = register_diag_field ( mod_name,    &
    'drp_evapdt_clubb', axes(1:3), Time,                  &
    'droplet evaporation rate from clubb (boundary layer cloud)', &
    '#/kg/s', missing_value=missing_value)

  id_aer_ccn_act_clubb = register_diag_field ( mod_name,    &
    'aer_ccn_actdt_clubb', axes(1:3), Time,                  &
    'droplet activation rate from clubb (boundary layer cloud)', &
    '#/kg/s', missing_value=missing_value)

  id_Ndrop_act_clubb = register_diag_field ( mod_name,    &
    'Ndrop_act_clubb', axes(1:3), Time,                  &
    'Maximum in-cloud droplet activation number concentration in clubb (boundary layer cloud)', &
    '#/kg', missing_value=missing_value)

  id_drp_flux_clubb = register_diag_field ( mod_name,    &
    'drp_flux_clubb', axes(half), Time,                  &
    'droplet number turbulent flux from clubb (boundary layer cloud)', &
    'm/kg/s', missing_value=missing_value)

  id_qndt_clubb_trsport_only = register_diag_field ( mod_name,    &
    'qndt_clubb_trsport_only', axes(1:3), Time,                  &
    'tendency of droplet number from pure transport in clubb (boundary layer cloud)', &
     'm/kg/s', missing_value=missing_value)

  id_ice_evap_clubb = register_diag_field ( mod_name,    &
    'ice_evapdt_clubb', axes(1:3), Time,                  &
    'ice crystal evaporation rate from clubb (boundary layer cloud)', &
    '#/kg/s', missing_value=missing_value)

  id_aer_ice_act_clubb = register_diag_field ( mod_name,    &
    'aer_ice_actdt_clubb', axes(1:3), Time,                  &
    'ice nucleation rate from clubb (boundary layer cloud)', &
    '#/kg/s', missing_value=missing_value)

  id_Icedrop_act_clubb = register_diag_field ( mod_name,    &
    'Icedrop_act_clubb', axes(1:3), Time,                  &
    'Maximum ice nucleation number concentration in clubb (boundary layer cloud)', &
    '#/kg', missing_value=missing_value)

  id_icedrop_flux_clubb = register_diag_field ( mod_name,    &
    'icedrop_flux_clubb', axes(half), Time,                  &
    'ice number turbulent flux from clubb (boundary layer cloud)', &
    'm/kg/s', missing_value=missing_value)

  id_qnidt_clubb_trsport_only = register_diag_field ( mod_name,    &
    'qnidt_clubb_trsport_only', axes(1:3), Time,                  &
    'tendency of ice number from pure transport in clubb (boundary layer cloud)', &
    'm/kg/s', missing_value=missing_value)

  id_sulfate = register_diag_field ( mod_name, 'sulfate', &
    axes(1:3), Time, 'sulfate mass conentration',     &
    'ug so4/m3', missing_value=missing_value )
 
  id_seasalt_sub = register_diag_field ( mod_name, 'seasalt_sub', &
    axes(1:3), Time, 'sub-micron sea salt mass conentration',     &
    'ug/m3', missing_value=missing_value )
 
  id_seasalt_sup = register_diag_field ( mod_name, 'seasalt_sup', &
    axes(1:3), Time, 'super-micron sea salt mass conentration',     &
    'ug/m3', missing_value=missing_value )
 
  id_om = register_diag_field ( mod_name, 'OM', &
    axes(1:3), Time, 'OM mass conentration',     &
    'ug/m3', missing_value=missing_value )

  ! register higher order terms and fluxes
  ! at full levels
  id_wp3_clubb = register_diag_field ( mod_name,    &
    'wp3_clubb', axes(1:3), Time,                  &
    'w third order moment from clubb (boundary layer cloud)', &
    'm^3/s^3', missing_value=missing_value)

  id_wm_clubb = register_diag_field ( mod_name,    &
    'wm_clubb', axes(1:3), Time,                  &
    'vertical velocity W adjusted by (donner) convective mass flux', &
    'm/s', missing_value=missing_value)

  id_omega_clubb = register_diag_field ( mod_name,    &
    'omega_clubb', axes(1:3), Time,                  &
    'vertical velocity Omega adjusted by (donner) convective mass flux', &
    'Pa/s', missing_value=missing_value)

  id_qcvar_clubb = register_diag_field ( mod_name,    &
    'qcvar_clubb', axes(1:3), Time,                  &
    'inverse of the relative variance of in-cloud water content', &
    'none', missing_value=missing_value)

  id_qcvar_cf_clubb = register_diag_field ( mod_name,    &
    'qcvar_cf_clubb', axes(1:3), Time,                  &
    'inverse of the relative variance of in-cloud water content, and weighted by cloud fraction', &
    'none', missing_value=missing_value)


  id_rcm_clubb = register_diag_field ( mod_name,    &
    'rcm_clubb', axes(1:3), Time,                  &
    'grid box mean cloud water content from CLUBB', &
    'kg/kg', missing_value=missing_value)

  id_rcp2_clubb = register_diag_field ( mod_name,    &
    'rcp2_clubb', axes(1:3), Time,                  &
    'grid box cloud water variance from CLUBB', &
    'kg^2/kg^2', missing_value=missing_value)

  id_cf_clubb = register_diag_field ( mod_name,    &
    'cf_clubb', axes(1:3), Time,                  &
    'cloud fraction from CLUBB', &
    'none', missing_value=missing_value)


  id_cf_avg_clubb = register_diag_field ( mod_name,    &
    'cf_avg_clubb', axes(1:3), Time,                  &
    'time-averaged cloud fraction from CLUBB', &
    'none', missing_value=missing_value)


  ! at half levels
  id_wp2_clubb = register_diag_field ( mod_name,    &
    'wp2_clubb', axes(half), Time,                  &
    'w second order moment from clubb (boundary layer cloud) at half levels', &
    'm^2/s^2', missing_value=missing_value)

  id_upwp_clubb = register_diag_field ( mod_name,    &
    'upwp_clubb', axes(half), Time,                  &
    'u and w co-variance from clubb (boundary layer cloud) at half levels', &
    'm^2/s^2', missing_value=missing_value)

  id_vpwp_clubb = register_diag_field ( mod_name,    &
    'vpwp_clubb', axes(half), Time,                  &
    'v and w co-variance from clubb (boundary layer cloud) at half levels', &
    'm^2/s^2', missing_value=missing_value)

  id_up2_clubb = register_diag_field ( mod_name,    &
    'up2_clubb', axes(half), Time,                  &
    'u second order moment from clubb (boundary layer cloud) at half levels', &
    'm^2/s^2', missing_value=missing_value)
 
  id_vp2_clubb = register_diag_field ( mod_name,    &
    'vp2_clubb', axes(half), Time,                  &
    'v second order moment from clubb (boundary layer cloud) at half levels', &
    'm^2/s^2', missing_value=missing_value)

  id_wprtp_clubb = register_diag_field ( mod_name,    &
    'wprtp_clubb', axes(half), Time,                  &
    'w and rtm co-variance from clubb (boundary layer cloud) at half levels', &
    'm/s', missing_value=missing_value)

  id_wpthlp_clubb = register_diag_field ( mod_name,    &
    'wpthlp_clubb', axes(half), Time,                  &
    'w and thlm co-variance from clubb (boundary layer cloud) at half levels', &
    'm/s K', missing_value=missing_value)

  id_rtp2_clubb = register_diag_field ( mod_name,    &
    'rtp2_clubb', axes(half), Time,                  &
    'rtm  second order moment from clubb (boundary layer cloud) at half levels', &
    ' ', missing_value=missing_value)

  id_thlp2_clubb = register_diag_field ( mod_name,    &
    'thlp2_clubb', axes(half), Time,                  &
    'thlm second order moment from clubb (boundary layer cloud) at half levels', &
    'K^2', missing_value=missing_value)

  id_rtpthlp_clubb = register_diag_field ( mod_name,    &
    'rtpthlp_clubb', axes(half), Time,                  &
    'rtm and thlm co-variance from clubb (boundary layer cloud) at half levels', &
    'K^2', missing_value=missing_value)


  id_wp2_vert_avg_clubb = register_diag_field ( mod_name,    &
    'wp2_vert_avg_clubb', axes(1:2), Time,                  &
    'column vertically averaged wp2 from clubb', &
    'm^2/s^2', missing_value=missing_value)

  id_up2_vert_avg_clubb = register_diag_field ( mod_name,    &
    'up2_vert_avg_clubb', axes(1:2), Time,                  &
    'column vertically averaged up2 from clubb', &
    'm^2/s^2', missing_value=missing_value)

  id_vp2_vert_avg_clubb = register_diag_field ( mod_name,    &
    'vp2_vert_avg_clubb', axes(1:2), Time,                  &
    'column vertically averaged vp2 from clubb', &
    'm^2/s^2', missing_value=missing_value)

  id_rtp2_vert_avg_clubb = register_diag_field ( mod_name,    &
    'rtp2_vert_avg_clubb', axes(1:2), Time,                  &
    'column vertically averaged rtp2 from clubb', &
    'kg^2/kg^2', missing_value=missing_value)

  id_thlp2_vert_avg_clubb = register_diag_field ( mod_name,    &
    'thlp2_vert_avg_clubb', axes(1:2), Time,                  &
    'column vertically averaged thlp2 from clubb', &
    'K^2', missing_value=missing_value)


  ! Initialize clubb internal stats
  if (do_stats) then

    ! For GCM identify column of interest, for SCM set to 1,1

#ifdef SCM
    istats = 1
    jstats = 1
    fstats = 1
#else
    do i=1,id
      do j=1,jd
        do n=1,size(lon_stats)
          if ( abs(lon(i,j)*RAD_TO_DEG-lon_stats(n)).lt.0.5 .and.  &
               abs(lat(i,j)*RAD_TO_DEG-lat_stats(n)).lt.0.5 ) then
            istats = i
            jstats = j
            fstats = n
          endif
        enddo
      enddo
    enddo
#endif

    ! If stats is active on current PE
    if ( istats.ge.0 .and. jstats.ge.0 ) then

      write(*,'(a,2i4,2f12.4,a)') 'clubb_driver: stats active for ',   &
      istats, jstats,                                                  &
      RAD_TO_DEG*lon(istats,jstats), RAD_TO_DEG*lat(istats,jstats),    &
      trim(fname_prefix(fstats))

      unit_stats = open_namelist_file()
      close(unit_stats)

      Time_stats_init = Time
      call get_date( Time_stats_init, year, month, day, hour, minute, second )
      time_current = 3600.0*hour + 60.0*minute + second

      call stats_init( unit_stats, fname_prefix(fstats), "",                       &
                       do_stats, stats_fmt, stats_tsamp, stats_tout, 'input.nml',  &
                       gr%nz, gr%zt, gr%zm,                                        &
                       gr%nz, gr%zt, gr%nz, gr%zm,                                 &
                       day, month, year,                                           & 
                       RAD_TO_DEG*lat(istats:istats,jstats),                       &
                       RAD_TO_DEG*lon(istats:istats,jstats),                       &
                       time_current, clubb_dt)

    endif

  endif  ! do_stats

  return
  end subroutine clubb_init
  !=====================================================================

  !=====================================================================
  subroutine clubb_setup(id, jd, phalf)

  !---------------------------------------------------------------------
  !  id, jd                             input
  !    subdomain dimensions
  !
  !  phalf                              input
  !    pressure at half levels in pascals
  !    [real, dimension(nlon,nlat,nlev+1)]
  !
  !---------------------------------------------------------------------

  ! ----- Calling arguments -----

  integer, intent(in)                 :: id, jd
  real, dimension(:,:,:), intent(in)  ::  phalf

  ! ----- Local variables -----

  real, parameter :: Hscale = 7500.0
  integer :: ntot_num_trs      ! total # of number concentrations
  integer :: err_code          ! valid run?

  integer :: kdim
  real, dimension(nparams) :: params  ! Array of the model constants, initialize in "clubb_setup"

  integer, parameter        :: sclr_max = 10     ! Maximum number of passive scalars (arbitrary)
  real, dimension(sclr_max) :: sclr_tol = 1.e-2  ! Thresholds on the passive scalars [units vary]

  ! ----- Begin code -----

  kdim = size(phalf,3)-1

  ! Define model constant parameters
  ntot_num_trs = 0
  if (nqn > 0)  ntot_num_trs = ntot_num_trs + 1
  if (nqni > 0) ntot_num_trs = ntot_num_trs + 1

  ! allocate momentum and thermodynamic heights
  allocate( momentum_heights(1:kdim+1) )
  allocate( thermodynamic_heights(1:kdim+1) )

  ! using scale height (Hscale=7500) to get a generic momentum height
  do iz_clubb =1, kdim
    momentum_heights(iz_clubb)= Hscale*log(phalf(id, jd,kdim+1)/phalf(id, jd,kdim+2-iz_clubb))
  enddo
  iz_clubb=kdim+1
 momentum_heights(iz_clubb)=2.0*momentum_heights(iz_clubb-1)-momentum_heights(iz_clubb-2)

  ! using linear interpolation to get a generic thermodynamic height
  do iz_clubb =2, kdim+1
    thermodynamic_heights(iz_clubb)=0.5*( momentum_heights(iz_clubb) &
                                         +momentum_heights(iz_clubb-1) )
  enddo
  iz_clubb=1
  thermodynamic_heights(iz_clubb)=  2.0*thermodynamic_heights(iz_clubb+1) &
                                       -thermodynamic_heights(iz_clubb+2)

  ! Read in model parameter values
  call read_parameters(1,"input.nml", params )   ! intent (out)

  call set_clubb_debug_level( debug_level ) ! Intent(in)

  call  setup_clubb_core                                    &
        ( kdim+1, 300.0, 0.0,                               & ! In: nzmax, T0_in, ts_nudge_in
          0, ntot_num_trs,                                  & ! In: hydromet_dim_in, sclr_dim_in
          sclr_tol, ntot_num_trs, params,                   & ! In: sclr_tol_in, edsclr_dim_in, params
          l_host_applies_sfc_fluxes,                        & ! In: l_host_applies_sfc_fluxes
          .false., saturation_formula_in,                   & ! In: l_uv_nudge, saturation_formula
          .true.,                                           & ! In: I_sat_sphum
          .true., 2, avg_deltaz,                            & ! In: l_implemented, grid_type, deltaz
          momentum_heights(1), momentum_heights(kdim+1),    & ! In: zm_init, zm_top
          momentum_heights, thermodynamic_heights,          & ! In: momentum_heights, thermodynamic_heights
          host_dx, host_dy, momentum_heights(1),            & ! In: host_dx, host_dy, sfc_elevation
          cloud_frac_min_in ,                               & ! In: cloud_frac_min
          err_code )                                          ! Out

  return
  end subroutine clubb_setup
  !=====================================================================

  !=====================================================================
  subroutine clubb_end

  use clubb_core, only: cleanup_clubb_core  ! Procedure

  implicit none
  character(len=64)     :: fname='RESTART/clubb.res.nc'
  character(len=64)     :: sclr_prefix,  sclr_txt
  integer               ::  isclr_dim
   
  ! ----- verify that the module is initialized -----
  if ( .not. module_is_initialized) then
    call error_mesg('clubb_driver_mod','module has not been initialized', FATAL)
    return
  else
    module_is_initialized = .false.
  endif

  ! If stats is active on current PE
  if ( istats.ge.0 .and. jstats.ge.0 ) then
    l_stats = .true.
    call stats_finalize()
  endif

  call cleanup_clubb_core( .true. )

  ! De-allocate the array for the passive scalar tolerances

  deallocate( momentum_heights )
  deallocate( thermodynamic_heights )
  deallocate (id_tracer_clubb)
  deallocate (id_tracer_clubb_col)

  call aer_ccn_act_end
  if ( do_ice_nucl_wpdf ) call ice_nucl_wpdf_end

  ! ----- mark the module as uninitialized -----
  module_is_initialized = .false.

  return
  end subroutine clubb_end
  !=====================================================================

  !=====================================================================
  subroutine host2clubb_full (var_host, var_clubb)
  real, intent(in) ,  dimension(:)     ::  var_host
  real, intent(out) , dimension(:)     ::  var_clubb

  ! Fill clubb grid ghost points by extrapolation
  var_clubb(1) = 2.0*var_host(kdim) - var_host(kdim-1)
  ! Copy interior clubb grid points
  do iz_clubb=2,kdim+1
    var_clubb(iz_clubb) = var_host(kdim+2-iz_clubb)
  enddo

  return
  end subroutine host2clubb_full
  !=====================================================================

  !=====================================================================
  subroutine host2clubb_half (var_host, var_clubb)
  real, intent(in),  dimension(:)     ::  var_host
  real, intent(out), dimension(:)     ::  var_clubb

  do iz_clubb=1, kdim+1
     iz_host   = kdim + 2 - iz_clubb
     var_clubb( iz_clubb) = var_host( iz_host )
  enddo

  return
  end subroutine host2clubb_half
  !=====================================================================

  !=====================================================================
  subroutine clubb2host_full (var_clubb, var_host)
  real, intent(in),  dimension(:)     ::  var_clubb
  real, intent(out), dimension(:)     ::  var_host

  do iz_clubb=2, kdim+1
     iz_host = kdim + 2 - iz_clubb
     var_host( iz_host ) = var_clubb( iz_clubb)      
  enddo

  return
  end subroutine clubb2host_full 
  !=====================================================================

  !=====================================================================
  subroutine clubb2host_half (var_clubb, var_host)
  real, intent(in),  dimension(:)     ::  var_clubb
  real, intent(out), dimension(:)     ::  var_host
    
  do iz_clubb=1, kdim+1
     iz_host = kdim + 2 - iz_clubb
     var_host( iz_host ) = var_clubb( iz_clubb)
  enddo

  return
  end subroutine clubb2host_half
  !=====================================================================



!#######################################################################
subroutine clubb_3D_2_1D( ix, iy, r )
implicit none
integer,  intent(in) ::  ix, iy
real,     intent(in), dimension(:,:,:,ntp+1:) ::  r

!=======================================================================
do iz_clubb=2, kdim+1
  ! iz_host = iz_clubb - 1
   iz_host = kdim + 2 - iz_clubb
   wp2( iz_clubb )    = r( ix, iy, iz_host, nwp2 )
   wp3( iz_clubb )    = r( ix, iy, iz_host, nwp3 )
   upwp( iz_clubb )   = r( ix, iy, iz_host, nupwp )
   vpwp( iz_clubb )   = r( ix, iy, iz_host, nvpwp )
   wprtp(iz_clubb )   = r( ix, iy, iz_host, nwprtp )
   wpthlp(iz_clubb )  = r( ix, iy, iz_host, nwpthlp )
   rtp2( iz_clubb )   = r( ix, iy, iz_host, nrtp2 )
   thlp2(iz_clubb )   = r( ix, iy, iz_host, nthlp2 )
   rtpthlp(iz_clubb ) = r( ix, iy, iz_host, nrtpthlp )
   up2( iz_clubb )    = r( ix, iy, iz_host, nup2 )
   vp2( iz_clubb )    = r( ix, iy, iz_host, nvp2 )

   if ( sclr_dim > 0) then
      RH_crit( iz_clubb, 1, 1 ) = r( ix, iy, iz_host, nrhcrit1 )
      RH_crit( iz_clubb, 1, 2 ) = r( ix, iy, iz_host, nrhcrit2 )
   endif
enddo

! fill ghost points
wp2(1) = wp2(2)
wp3(1) = wp3(2)
upwp(1) = upwp(2)
vpwp(1) = vpwp(2)
wprtp(1) = wprtp(2)
wpthlp(1) = wpthlp(2)
rtp2(1) = rtp2(2)
thlp2(1) = thlp2(2)
rtpthlp(1) = rtpthlp(2)
up2(1) = up2(2)
vp2(1) = vp2(2)
if ( sclr_dim > 0) then
     RH_crit(1, 1, 1:2) = RH_crit(2, 1, 1:2)
endif

return
end subroutine clubb_3D_2_1D
!#######################################################################    


!#######################################################################
subroutine clubb_1D_2_3D( ix, iy, r )
implicit none
integer,  intent(in) ::  ix, iy
real, intent(inout), dimension(:,:,:,ntp+1:) :: r

!=======================================================================
do iz_clubb=2, kdim+1
  ! iz_host = iz_clubb - 1
   iz_host = kdim + 2 - iz_clubb
   r(ix,iy,iz_host,nwp2    ) =     wp2(iz_clubb)
   r(ix,iy,iz_host,nwp3    ) =     wp3(iz_clubb)
   r(ix,iy,iz_host,nupwp   ) =    upwp(iz_clubb)
   r(ix,iy,iz_host,nvpwp   ) =    vpwp(iz_clubb)
   r(ix,iy,iz_host,nwprtp  ) =   wprtp(iz_clubb)
   r(ix,iy,iz_host,nwpthlp ) =  wpthlp(iz_clubb)
   r(ix,iy,iz_host,nrtp2   ) =    rtp2(iz_clubb)
   r(ix,iy,iz_host,nthlp2  ) =   thlp2(iz_clubb)
   r(ix,iy,iz_host,nrtpthlp) = rtpthlp(iz_clubb)
   r(ix,iy,iz_host,nup2    ) =     up2(iz_clubb)
   r(ix,iy,iz_host,nvp2    ) =     vp2(iz_clubb)
   if ( sclr_dim > 0) then
     r(ix,iy,iz_host,nrhcrit1) = RH_crit(iz_clubb,1,1)
     r(ix,iy,iz_host,nrhcrit2) = RH_crit(iz_clubb,1,2)
   endif
enddo

return
end subroutine clubb_1D_2_3D
!#######################################################################    



!#######################################################################
subroutine tndy_clubb ( dtmain, var_clubb,     &
                        var_host, vardt_host,  &
                        vardt_clubb )
!=======================================================================
real :: dtmain
real, intent(in),    dimension(:)    ::   var_clubb
real, intent(in),    dimension(:)    ::   var_host
real, intent(inout), dimension(:)    ::   vardt_host
real, intent(inout), dimension(:)    ::   vardt_clubb
 
do iz_clubb=2, kdim+1
   iz_host = kdim + 2 - iz_clubb
   vardt_clubb( iz_host ) =                                             &
            ( var_clubb( iz_clubb ) - var_host( iz_host ) )/dtmain   &
           - vardt_host( iz_host)
 
   vardt_host( iz_host) = ( var_clubb( iz_clubb ) - var_host( iz_host ) )/dtmain
enddo
  
return
end subroutine tndy_clubb
!#######################################################################    



recursive function erff(x) result(y)
! Error function from Numerical Recipes.
! erf(x) = 1 - erfc(x)

real dumerfc, x
real t, z, y

z = abs(x)
t = 1.0 / ( 1.0 + 0.5 * z )

dumerfc =     t * exp(-z * z - 1.26551223 + t *      &
          ( 1.00002368 + t * ( 0.37409196 + t *    &
          ( 0.09678418 + t * (-0.18628806 + t *    &
          ( 0.27886807 + t * (-1.13520398 + t *    &
          ( 1.48851587 + t * (-0.82215223 + t * 0.17087277 )))))))))

if ( x.lt.0.0 ) dumerfc = 2.0 - dumerfc
 
y = 1.0 - dumerfc

end function erff
!#######################################################################    



!#######################################################################    
 subroutine aer_act_clubb_BL_Gauss( aeromass_clubb, & ! Intent(in)
                                    Ndrop_max )       ! Intent(out)
! using Boucher an Lohmann empirical relationship between Nd and sulface mass concentration
! Nd = 10^(2.21+0.41 log(mso4))
!=======================================================================
use  aer_ccn_act_k_mod,   only: aer_ccn_act_k, aer_ccn_act_wpdf_k 

implicit none
real, intent(inout),  dimension(:, :)    ::    aeromass_clubb
real, intent(inout),  dimension(:)       ::    Ndrop_max
real ::   drop
  
!=======================================================================
!  if      aeromass_clubb = 1 ug/m3,  drop = 162 /cm3
!  if      aeromass_clubb = 1 ug/m3,  drop = 314 /cm3
do iz_clubb = 2, kdim+1
   drop = 2.21 + 0.41 * log10(  aeromass_clubb ( iz_clubb, 1 )*1.e12  )
   drop = exp( log(10.0) * drop )
   Ndrop_max ( iz_clubb ) = drop  &
          * ( pdf_params( iz_clubb)%mixt_frac * pdf_params( iz_clubb)%cloud_frac1  &
             +( 1.0-pdf_params( iz_clubb)%mixt_frac ) * pdf_params( iz_clubb)%cloud_frac2 )     
enddo  
return
end subroutine  aer_act_clubb_BL_Gauss
!#######################################################################    


subroutine aer_act_clubb_diffK_Gauss( aeromass_clubb, temp_clubb_act, &   ! Intent(in)
                                                      Ndrop_max )         ! Intent(out)
use aer_ccn_act_k_mod,   only: aer_ccn_act_k, aer_ccn_act_wpdf_k 

implicit none
real,  intent(in),      dimension(:)     ::    temp_clubb_act
real, intent(inout),    dimension(:, :)  ::    aeromass_clubb
real, intent(out),    dimension(:)       ::    Ndrop_max

real                ::  drop
integer             ::  Tym, ier
character(len=256)  ::  ermesg
!=======================================================================
Tym = size(aeromass_clubb, 2)

do iz_clubb = 2, kdim+1
   call aer_ccn_act_wpdf_k(temp_clubb_act( iz_clubb), p_in_Pa( iz_clubb), &! intent (in)
                           wm_zt( iz_clubb),          Var_w,              &! intent (in)
                           aeromass_clubb( iz_clubb, : ), Tym,            &! intent (in)
                           drop,   ier,   ermesg )                         ! intent (out)

    Ndrop_max ( iz_clubb ) =   drop &
            * ( pdf_params(iz_clubb)%mixt_frac * pdf_params(iz_clubb)%cloud_frac1  &
               +( 1.0-pdf_params(iz_clubb)%mixt_frac)*pdf_params(iz_clubb)%cloud_frac2 )
enddo  
return
end subroutine  aer_act_clubb_diffK_Gauss
!###################################################################    






subroutine aer_act_clubb_quadrature_Gauss(Time_next, aeromass_clubb, temp_clubb_act, & ! Intent(in)
                                          Ndrop_max   )                                ! Intent(out)
!=======================================================================
use aer_ccn_act_k_mod,   only: aer_ccn_act_k, aer_ccn_act_wpdf_k 

implicit none
! ---> h1g, 2010-08-24, dumping Nact
type(time_type),         intent(in)      ::    Time_next
! <--- h1g, 2010-08-24

real,  intent(in),    dimension(:)       ::    temp_clubb_act
real, intent(inout),  dimension(:, :)    ::    aeromass_clubb
real, intent(out),    dimension(:)       ::    Ndrop_max

real                ::  drop
integer             ::  Tym, ier
character(len=256)  ::  ermesg

real               :: P1_updraft,   P2_updraft  ! probability of updraft
real, parameter    :: P_updraft_eps = 1.e-16    ! updraft probability threshold
real, parameter    :: wp2_eps = 0.0001          ! w variance threshold
 
!=======================================================================
Tym = size(aeromass_clubb, 2)

do iz_clubb = 2, kdim+1

   if( pdf_params( iz_clubb)%varnce_w1 > wp2_eps) then
       P1_updraft = 0.5+0.5*erff( pdf_params( iz_clubb)%w1/sqrt(2.0*pdf_params( iz_clubb)%varnce_w1) )
       P1_updraft = P1_updraft * pdf_params( iz_clubb)%mixt_frac * pdf_params( iz_clubb)%cloud_frac1
   else
      if( pdf_params( iz_clubb)%w1 > 0.0) &
          P1_updraft = pdf_params( iz_clubb)%mixt_frac * pdf_params( iz_clubb)%cloud_frac1
   endif


   if( pdf_params( iz_clubb)%varnce_w2 > wp2_eps) then
       P2_updraft = 0.5+0.5*erff( pdf_params( iz_clubb)%w2/sqrt(2.0*pdf_params( iz_clubb)%varnce_w2) )
       P2_updraft = P2_updraft * ( 1.0-pdf_params( iz_clubb )%mixt_frac ) * pdf_params( iz_clubb)%cloud_frac2
   else
       if( pdf_params( iz_clubb)%w2 > 0.0) &
           P2_updraft = ( 1.0-pdf_params( iz_clubb )%mixt_frac ) * pdf_params( iz_clubb)%cloud_frac2
   endif

   if( P1_updraft + P2_updraft   > P_updraft_eps  ) then
       P1_updraft =  P1_updraft / (  P1_updraft + P2_updraft  )
       P2_updraft =  P2_updraft / (  P1_updraft + P2_updraft  )
   else
       P1_updraft = 0.0
       P2_updraft = 0.0
   endif
 
   call aer_ccn_act_wpdf_k( temp_clubb_act( iz_clubb),  p_in_Pa( iz_clubb),              &! intent (in)
                            pdf_params( iz_clubb)%w1,   pdf_params( iz_clubb)%varnce_w1, &! intent (in)
                            aeromass_clubb( iz_clubb, : ), Tym,                          &! intent (in)
                            drop,   ier,   ermesg )           ! intent (out)
    
    Ndrop_max ( iz_clubb ) = drop * P1_updraft

    call aer_ccn_act_wpdf_k( temp_clubb_act( iz_clubb),  p_in_Pa( iz_clubb),               &! intent (in)
                             pdf_params( iz_clubb)%w2,   pdf_params( iz_clubb)%varnce_w2,  &! intent (in)
                             aeromass_clubb( iz_clubb, : ), Tym,                           &! intent (in)
                             drop,   ier,   ermesg)           ! intent (out)
! in-cloud activated droplet concentration
    Ndrop_max ( iz_clubb ) = Ndrop_max( iz_clubb ) + drop * P2_updraft

! get the layer-averaged activated droplet concentration (/cm3)
    Ndrop_max ( iz_clubb ) = Ndrop_max ( iz_clubb ) *  &
                 (      pdf_params( iz_clubb)%mixt_frac  * pdf_params( iz_clubb)%cloud_frac1 + &
                   (1.- pdf_params( iz_clubb)%mixt_frac) * pdf_params( iz_clubb)%cloud_frac2 )

enddo
return
end subroutine aer_act_clubb_quadrature_Gauss
!#######################################################################    


! this subroutine is to droplet and ice crystal activation

subroutine warm_ice_nucl_mns_clubb( aeromass_clubb,      &    !  Intent(in)
                                    aeroimass_clubb,     &    !  Intent(in)
                                    temp_clubb_act,      &    !  Intent(in)  
                                    Ndrop_max,           &    !  Intent(out)  
                                    Ncrystal_max,        &    !  Intent(out)
                                    RH_crit_clubb )           !  Intent(out)

use  aer_ccn_act_k_mod,   only: aer_ccn_act_k,  aer_ccn_act_wpdf_k 
use  sat_vapor_pres_mod,  only: compute_qs
USE  polysvp_mod,         only: polysvp_l,  polysvp_i
USE  ice_nucl_mod,        only: ice_nucl_wpdf


implicit none          
real,  intent(in),    dimension(:)       ::    temp_clubb_act
real,  intent(inout), dimension(:, :)    ::    aeromass_clubb
real,  intent(inout), dimension(:, :)    ::    aeroimass_clubb
     
real, intent(out),    dimension(:)       ::    Ndrop_max     ! domain warm cloud droplet number concentration
real, intent(out),    dimension(:)       ::    Ncrystal_max  ! domain ice crystal number concentration
!real, intent(inout),  dimension( gr%nz, 1:min(1,sclr_dim), 2 )  :: RH_crit_clubb       ! critical relative humidity for nucleation
real, intent(inout),  dimension(:, :, :) :: RH_crit_clubb       ! critical relative humidity for nucleation
! -------------------------------------------------------------------------------------------------------------
!local variables for droplet nucleation
real               :: drop
integer            :: Tym, ier
character(len=256) :: ermesg

real               :: P1_updraft,   P2_updraft ! probability of updraft
real, parameter    :: P_updraft_eps = 1.e-16   ! updraft probability threshold
real, parameter    :: wp2_eps = 0.0001         ! w variance threshold

! -------------------------------------------------------------------------------------------------------------
!local variables for ice nucleation
real            ::      rh_crit_1d
real            ::      rh_crit_min_1d

real            ::      ni_sulf, ni_dust, ni_bc

real, dimension(  n_totmass  ) :: totalmass !r1 , mass_ratio
real, dimension( n_imass ) :: imass
real, parameter :: d378 = 0.378

real  ::   qs
real  ::   eslt,  qs_d, qvsl,  esit,  qvsi,  qvt,  u_i,  u_l
real  ::   wp2i
real  ::   Ncrystal1, Ncrystal2
real  ::   t1_combined,  t2_combined,  t3_combined
  
real, dimension( kdim+1 ) ::  hom
!===========================================================

Tym = size(aeromass_clubb, 2)

do iz_clubb = 2, kdim+1

  P1_updraft = 0.0
  P2_updraft = 0.0
  if( pdf_params( iz_clubb)%varnce_w1 > wp2_eps) then
    P1_updraft = 0.5+0.5*erff( pdf_params( iz_clubb)%w1/sqrt(2.0*pdf_params( iz_clubb)%varnce_w1) )
    P1_updraft = P1_updraft * pdf_params( iz_clubb)%mixt_frac * pdf_params( iz_clubb)%cloud_frac1
  else
    if( pdf_params( iz_clubb)%w1 > 0.0) &
                 P1_updraft = pdf_params( iz_clubb)%mixt_frac * pdf_params( iz_clubb)%cloud_frac1
  endif


  if( pdf_params( iz_clubb)%varnce_w2 > wp2_eps) then
    P2_updraft = 0.5+0.5*erff( pdf_params( iz_clubb)%w2/sqrt(2.0*pdf_params( iz_clubb)%varnce_w2) )
    P2_updraft = P2_updraft * ( 1.0-pdf_params( iz_clubb )%mixt_frac ) * pdf_params( iz_clubb)%cloud_frac2
  else
    if( pdf_params( iz_clubb)%w2 > 0.0) &
                 P2_updraft = ( 1.0-pdf_params( iz_clubb )%mixt_frac ) * pdf_params( iz_clubb)%cloud_frac2
  endif

  if( P1_updraft + P2_updraft > P_updraft_eps ) then
    P1_updraft =  P1_updraft / (  P1_updraft + P2_updraft  )
    P2_updraft =  P2_updraft / (  P1_updraft + P2_updraft  )
  endif

! for the first Gaussian distribution portion of the warm cloud nucleation 
  call aer_ccn_act_wpdf_k ( temp_clubb_act( iz_clubb), p_in_Pa( iz_clubb), &  !  intent(in)
                            pdf_params( iz_clubb)%w1,                      &  !  intent(in)
                            pdf_params( iz_clubb)%varnce_w1,               &  !  intent(in)
                            aeromass_clubb( iz_clubb, : ), Tym,            &  !  intent(in)
                            drop,   ier,   ermesg )                           !  intent(out)

  Ndrop_max ( iz_clubb ) = drop * P1_updraft


! for the first Gaussian distribution portion of the ice cloud nucleation 
  imass(:)     = aeroimass_clubb(iz_clubb,:)
  totalmass(:) = aeromass_clubb(iz_clubb,:)

  eslt         =  polysvp_l( temp_clubb_act(iz_clubb) )  !satuarated vapor pressure with respect to liquid
  qs_d         =  p_in_Pa(iz_clubb) - d378*eslt
  qs_d         =  max(qs_d,eslt)
  qvsl         = 0.622 *eslt / qs_d

  esit         = polysvp_i( temp_clubb_act(iz_clubb) )  !satuarated vapor pressure with respect to ice
  qs_d         = p_in_Pa(iz_clubb) - d378*esit
  qs_d         = max(qs_d,esit)
  qvsi         = 0.622 *esit / qs_d

!cms 4/21/2009 changed nothing
  t1_combined = 273.16
  t2_combined = 268.16
  t3_combined = 238.16 
  qvt = qvsl
  if ( rh_act_opt .EQ. 1 ) THEN
     qvt = rtm(iz_clubb) - rcm(iz_clubb)
  else
!cms 4/23/2009  changed 
!environmental qv
     if ( cloud_frac(iz_clubb) .LT. cf_thresh_nucl ) THEN
         ! call compute_qs( temp_clubb_act(iz_clubb),  p_in_Pa(iz_clubb), qs)
          if( temp_clubb_act(iz_clubb) > t1_combined ) then
             qs = qvsl
          elseif( temp_clubb_act(iz_clubb) > t2_combined ) then
             qs =  qvsl * (temp_clubb_act(iz_clubb) - t2_combined)/(t1_combined - t2_combined) &
                  +qvsi * (t1_combined - temp_clubb_act(iz_clubb))/(t1_combined - t2_combined)
             qs = qs* RH_factor
          elseif( temp_clubb_act(iz_clubb) > t3_combined ) then
             qs =   qvsi &
                  + qvsi * (RH_crit_clubb(  iz_clubb, 1, 1 ) - 1.0 ) &
                    * (t2_combined - temp_clubb_act(iz_clubb))/(t2_combined - t3_combined)
             qs = qs* RH_factor
          else
             qs = qvsi * RH_crit_clubb(  iz_clubb, 1, 1 )
             qs = qs* RH_factor
          endif
          qvt =  qs
     else
          qvt =  rtm(iz_clubb) - rcm(iz_clubb)
     endif
  endif

  u_i =  qvt/qvsi
  u_l =  qvt/qvsl

  rh_crit_1d     = 1.
  rh_crit_min_1d = 1.

  hom     = 0.0
  ni_sulf = 0.
  ni_dust = 0.
  ni_bc = 0.

  IF (do_var_limit_ice) THEN
      wp2i = MAX (pdf_params( iz_clubb )%varnce_w1, var_limit_ice**2)
  ELSE
      wp2i = pdf_params( iz_clubb )%varnce_w1
  END IF
  call ice_nucl_wpdf( temp_clubb_act(iz_clubb), u_i, u_l,  &! intent (in) 
                      pdf_params( iz_clubb )%w1,         &! intent (in) 
                      wp2i,                              &! intent (in) 
                      thermodynamic_heights( iz_clubb ), &! intent (in)
                      totalmass, imass,                  &! intent (in)
                      n_totmass, n_imass,                &! intent (in)
                      Ncrystal1,                         &! intent (out)
                      drop,                              &! intent (in)
                      hom(iz_clubb),                     &! intent (inout)
                      rh_crit_1d, rh_crit_min_1d,        &! intent (out)
                      ni_sulf, ni_dust, ni_bc )           ! intent (out)

  if( sclr_dim > 0) then
    if( temp_clubb_act( iz_clubb ) .LT. 250. ) THEN
        RH_crit_clubb(  iz_clubb, 1, 1 ) = MAX(rh_crit_1d, 1.)
    else
        RH_crit_clubb(  iz_clubb, 1, 1 ) = 1.
    endif
  endif



! for the second Gaussian distribution portion of the warm cloud nucleation 
  call aer_ccn_act_wpdf_k ( temp_clubb_act( iz_clubb), p_in_Pa( iz_clubb), &   ! intent (in)
                            pdf_params( iz_clubb)%w2,                      &   ! intent (in) 
                            pdf_params( iz_clubb)%varnce_w2,               &   ! intent (in)
                            aeromass_clubb( iz_clubb, : ), Tym,            &   ! intent (in)
                            drop, ier, ermesg )                                ! intent (out)

! in-cloud activated droplet concentration
  Ndrop_max ( iz_clubb ) = Ndrop_max ( iz_clubb ) + drop *  P2_updraft

! get the layer-averaged activated droplet concentration (/cm3)
  Ndrop_max ( iz_clubb ) = Ndrop_max ( iz_clubb ) *  &
                   (    pdf_params( iz_clubb)%mixt_frac * pdf_params( iz_clubb)%cloud_frac1 + &
                   (1.- pdf_params( iz_clubb)%mixt_frac)* pdf_params( iz_clubb)%cloud_frac2 )


! for the second Gaussian distribution portion of the ice cloud nucleation
  if ( rh_act_opt .EQ. 1 ) THEN
     qvt = rtm(iz_clubb) - rcm(iz_clubb)
  else
     if ( cloud_frac(iz_clubb) .LT. cf_thresh_nucl ) THEN
         ! call compute_qs( temp_clubb_act(iz_clubb),  p_in_Pa(iz_clubb), qs)
          if( temp_clubb_act(iz_clubb) > t1_combined ) then
             qs = qvsl
          elseif( temp_clubb_act(iz_clubb) > t2_combined ) then
             qs =  qvsl * (temp_clubb_act(iz_clubb) - t2_combined)/(t1_combined - t2_combined) &
                  +qvsi * (t1_combined - temp_clubb_act(iz_clubb))/(t1_combined - t2_combined)
             qs = qs* RH_factor
          elseif( temp_clubb_act(iz_clubb) > t3_combined ) then
             qs =   qvsi &
                  + qvsi * (RH_crit_clubb(  iz_clubb, 1, 2 ) - 1.0 ) &
                    * (t2_combined - temp_clubb_act(iz_clubb))/(t2_combined - t3_combined)
             qs = qs* RH_factor
          else
             qs = qvsi * RH_crit_clubb(  iz_clubb, 1, 2 )
             qs = qs* RH_factor
          endif
          qvt =  qs
     else
          qvt =  rtm(iz_clubb) - rcm(iz_clubb)
     endif
  endif

  u_i =  qvt/qvsi
  u_l =  qvt/qvsl

  rh_crit_1d     = 1.
  rh_crit_min_1d = 1.
  IF (do_var_limit_ice) THEN
      wp2i = MAX (pdf_params( iz_clubb )%varnce_w2, var_limit_ice**2)
  ELSE
      wp2i = pdf_params( iz_clubb )%varnce_w2
  END IF
  call ice_nucl_wpdf( temp_clubb_act(iz_clubb), u_i, u_l,  &! intent (in) 
                      pdf_params( iz_clubb )%w2,         &! intent (in) 
                      wp2i,                              &! intent (in) 
                      thermodynamic_heights( iz_clubb ), &! intent (in)
                      totalmass, imass,                  &! intent (in)
                      n_totmass, n_imass,                &! intent (in)
                      Ncrystal2,                         &! intent (out)
                      drop,                              &! intent (in)
                      hom(iz_clubb),                     &! intent (inout)
                      rh_crit_1d, rh_crit_min_1d,        &! intent (out)
                      ni_sulf, ni_dust, ni_bc )           ! intent (out)

  if( sclr_dim > 0) then
    if( temp_clubb_act( iz_clubb ) .LT. 250. ) THEN
        RH_crit_clubb(  iz_clubb, 1, 2 ) = MAX(rh_crit_1d, 1.)
    else
        RH_crit_clubb(  iz_clubb, 1, 2 ) = 1.
    endif
  endif

  Ncrystal_max(iz_clubb) = Ncrystal1 * pdf_params( iz_clubb)%mixt_frac * pdf_params( iz_clubb) %cloud_frac1&
                         + Ncrystal2 * (1.0 - pdf_params( iz_clubb)%mixt_frac) * pdf_params( iz_clubb)%cloud_frac2
enddo

return
end subroutine  warm_ice_nucl_mns_clubb
!#######################################################################    





  ! ----------------------------------------------------------------------------
  ! This subroutine adds model tendencies to the clubb prognostic variables
  
  subroutine add_host_tdcy(                                                    &
               ix, iy, clubb_dt,                         &
               udt, vdt, tdt, rdt, env_qv_scale, env_condensate_scale,         &
               u_star, b_star, q_star,                                         &
               um, vm, thlm, rtm, cloud_frac, rcm, edsclrm, sclrm,             &
               upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc )

  implicit none

  ! Calling arguments
  integer, intent(in)                  :: ix, iy
  real, intent(in)                     :: clubb_dt   ! time-step to apply tendencies
  real, intent(in), dimension(:,:,:)   :: udt, vdt, tdt
  real, intent(in), dimension(:,:,:,:) :: rdt
  real, intent(in), dimension(:,:,:)   :: env_qv_scale, env_condensate_scale
  real, intent(in), dimension(:,:)     :: u_star, b_star, q_star

  real, intent(inout), dimension(:)    :: um, vm, thlm, rtm, cloud_frac, rcm
  real, intent(inout), dimension(:,:)  :: edsclrm, sclrm
  real, intent(out)                    :: upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc

  ! Local variables
  real, dimension(size(thlm,1)) :: tmp_clubb

  ! Cloud fraction
  tmp_clubb = 0.0
  call host2clubb_full(rdt(ix, iy, :, nqa),  &  ! intent (in)
                       tmp_clubb)                                     ! intent (out)
  cloud_frac    = cloud_frac + clubb_dt*tmp_clubb
  cloud_frac(1) = cloud_frac(2)
  where ( cloud_frac(:) <=cloud_frac_min )
    cloud_frac(:) = 0.0
  end where

  ! Water vapor
  tmp_clubb = 0.0
  call host2clubb_full ( rdt(ix, iy, : , nsphum)    &
                         /env_qv_scale(ix, iy, :),  &  !intent (in)
                         tmp_clubb )                                         !intent (out)
  rtm = rtm + clubb_dt * tmp_clubb

  ! Cloud water
  tmp_clubb = 0.0
  call host2clubb_full ( rdt( ix, iy, : , nql)              &
                         /env_condensate_scale(ix, iy, :),  &  !intent (in)
                         tmp_clubb )                                                 !intent (out) 
  rcm = rcm + clubb_dt * tmp_clubb
  rtm = rtm + clubb_dt * tmp_clubb

  ! Cloud ice
  tmp_clubb = 0.0
  call host2clubb_full ( rdt(ix, iy, : , nqi)               &
                         /env_condensate_scale(ix, iy, :),  &  !intent (in)
                         tmp_clubb )                                                 !intent (out)
  ! regard ice water as liquid water
! --->h1g, 2012-06-14
if( .not. do_liquid_only_in_clubb ) then
  rcm = rcm + clubb_dt * tmp_clubb
  rtm = rtm + clubb_dt * tmp_clubb
endif
! <---h1g, 2012-06-14
  rcm(1) = rcm(2)
  rtm(1) = rtm(2)

  ! Droplet number concentration
  if (nqn > 0) then
    tmp_clubb = 0.0
    call host2clubb_full( rdt(ix, iy, : , nqn)               &
                          /env_condensate_scale(ix, iy, :),  & !intent (in)
                          tmp_clubb )                                                !intent (out)
    if (use_sclr_HOC) then
      sclrm( :, 1 ) = sclrm( :, 1 ) + clubb_dt * tmp_clubb(:)
      sclrm( 1, 1 ) = sclrm( 2, 1 )
      do  iz_clubb = 1, kdim+1
        if( cloud_frac(iz_clubb) <= cloud_frac_min ) sclrm( iz_clubb , 1) = 0.0
      enddo
    else
      edsclrm( :, 1) = edsclrm( : , 1 ) + clubb_dt * tmp_clubb(:)
      edsclrm( 1, 1) = edsclrm( 2, 1 )
      do  iz_clubb = 1, kdim+1
        if( cloud_frac(iz_clubb) <= cloud_frac_min ) edsclrm( iz_clubb , 1) = 0.0
      enddo              
    endif
  endif

  ! Ice number concentration
  if (nqni > 0) then
    tmp_clubb = 0.0
    call host2clubb_full( rdt( ix, iy, : , nqni)             &
                          /env_condensate_scale(ix, iy, :),  & !intent (in)
                          tmp_clubb)                                                 !intent (out)
    if (use_sclr_HOC) then
      sclrm( : , 2 ) =  sclrm( : , 2 ) + clubb_dt * tmp_clubb( : )
      sclrm( 1 , 2 ) = sclrm( 2 , 2 )
      do  iz_clubb = 1, kdim+1
        if ( cloud_frac(iz_clubb) <= cloud_frac_min ) sclrm( iz_clubb , 2 ) = 0.0
      enddo
    else
      edsclrm( : , 2) = edsclrm( : , 2) + clubb_dt * tmp_clubb( :)
      edsclrm( 1 , 2) = edsclrm( 2 , 2)
      do  iz_clubb = 1, kdim+1
        if( cloud_frac(iz_clubb) <= cloud_frac_min ) edsclrm( iz_clubb , 2) = 0.0
      enddo
    endif
  endif

  ! Temperature
  call host2clubb_full( tdt(ix, iy, :),  & !intent (in)
                        tmp_clubb )                              !intent (out)
  temp_clubb = temp_clubb + clubb_dt * tmp_clubb

  ! Compute thlm
  thlm = (temp_clubb - Lv*rcm/Cp)/exner
  thlm(1) = thlm(2)

  ! Calculate surface u (upwp_sfc) and v (vpwp_sfc) momentum fluxes from u_star (input)
  upwp_sfc = -um(2) * u_star(ix, iy)**2  &
             / (max(gust_const, sqrt( um(2)*um(2) + vm(2)*vm(2)) ))
  vpwp_sfc = -vm(2) * u_star(ix, iy)**2  &
             / (max(gust_const, sqrt( um(2)*um(2) + vm(2)*vm(2)) ))

  ! Calculate surface sensible (wpthlp_sfc) and latend (wprtp_sfc) fluxes from u_star, b_star, and q_star (inputs)
  wpthlp_sfc = b_star( ix, iy ) * u_star( ix, iy ) *292./grav
  wprtp_sfc  = q_star( ix, iy ) * u_star( ix, iy )

  ! u wind [m/s]
  tmp_clubb = 0.0
  call host2clubb_full( udt( ix, iy, : ),  & !intent (in)
                        tmp_clubb)                                 !intent (out)
  um = um + clubb_dt * tmp_clubb
  um( 1 ) = um( 2 )

  ! v wind [m/s]
  tmp_clubb = 0.0
  call host2clubb_full( vdt( ix, iy, : ),  & !intent (in)
                        tmp_clubb )                                !intent (out)
  vm = vm + clubb_dt * tmp_clubb
  vm( 1 ) = vm( 2 )

  return
  end subroutine add_host_tdcy

#else
! NULL routines return error if called but not compiled for clubb
  subroutine clubb(is, ie, js, je, lon, lat,                  &
                   Time_next,                                 &
                   dtmain,                                    &
                   phalf, pfull, zhalf, zfull, omega_avg,     &
                   t, q, r, u, v,                             &
                   u_star, b_star, q_star,                    &
                   tdt, qdt, rdt, rdiag, udt, vdt,                   &
                   dcond_ls_liquid, dcond_ls_ice,             &
                   Ndrop_act_clubb, Icedrop_act_clubb,        &
                   ndust, rbar_dust,                          &
                   diff_t_clubb,                              &
                   qcvar_clubb,                               &
                   tdt_shf,  qdt_lhf ,                        &                   
                   Aerosol, mask,                             &
                   mc_full,                                   &
                   conv_frac_clubb,                           &
                   convective_humidity_ratio_clubb)

  !---------------------------------------------------------------------
  !
  !  is,ie,js,je                        input
  !    starting/ending subdomain i,j indices of data in the physics
  !    window being integrated
  !
  !  lon, lat                           input
  !    longitude, latitude in radians
  !
  !  Time_next                          input
  !    next time, used for diagnostics
  !
  !  dtmain                             input
  !    time step in seconds
  !
  !  phalf                              input
  !    pressure at half levels in pascals
  !
  !  pfull                              input
  !    pressure at full levels in pascals
  !
  !  zhalf                              input
  !    height at half levels in meters
  !
  !  zfull                              input
  !    height at full levels in meters
  !
  !  omega_avg                          input
  !    grid average omega at full levels in pascals per second
  !
  !  t, q                               input
  !    temperature (t) [deg k] and specific humidity
  !    of water vapor (q) [kg/kg] at full model levels,
  !
  !  r                                  input/output
  !    tracer fields at full model levels, unit varies,
  !
  !  u, v                               input
  !    zonal and meridional wind [m/s] at full model levels
  !  
  !  u_star                             input
  !    friction velocity [m/s]
  !
  !  b_star                             input
  !    buoyancy scale    [m/s^2/K]
  !
  !  q_star                             input
  !    moisture scale    dimensionless
  !
  !  tdt, qdt                           input/output
  !    temperature (tdt) [deg k/sec] and specific humidity of water 
  !    vapor (qdt) tendency [1/sec]
  !
  !  rdt                                input/output
  !    tracer tendencies , unit varies, 
  !
  !  udt, vdt                           input/output
  !    zonal and meridional wind tendencies [m/s/s]
  !
  !  dcond_ls_liquid                    output
  !    large-scale liquid condensation rate
  !
  !  dcond_ls_ice                       output
  !    large-scale ice condensation rate
  !
  !  Ndrop_act_clubb                    output
  !    in-cloud (NOT DOMAIN) averaged activated droplet number 
  !    concentration (#/kg)
  !
  !  Icedrop_act_clubb                  output
  !    in-cloud (NOT DOMAIN) averaged activated ice-crystal number 
  !    concentration (#/kg)
  !
  !  diff_t_clubb                       output
  !    tracer eddy diffusivity from clubb
  !
  !  Aerosol                            input, optional
  !    different aerosol species mass concentration
  !
  !  mask                               input, optional
  !    real array indicating the point is above the surface if equal 
  !    to 1.0 and indicating the point is below the surface if equal 
  !    to 0.
  !
  !  mc_full                            input, optional
  !    (total) net convective mass flux [kg/m2/s]
  !    used to adjust environmental vertical velocity 
  !    due to convection:
  !      omega = omega_avg + mc_full*grav
  !
  !  conv_frac_clubb                    input, optional
  !    total convective cloud fraction
  !
  !  convective_humidity_ratio_clubb    input, optional
  !    ratio of the grid average specific humidity 
  !    to environmental specific humidity 
  !
  !---------------------------------------------------------------------

  ! ----- Calling arguments -----

  integer, intent(in)                           ::  is, ie, js, je
  real, intent(in), dimension(:,:)              ::  lon, lat
  type(time_type), intent(in)                   ::  Time_next
  real, intent(in)                              ::  dtmain
  real, intent(in), dimension(:,:,:)            ::  phalf, pfull, zhalf, zfull, omega_avg
  real, intent(in), dimension(:,:,:)            ::  t, q, u, v
  real, intent(in), dimension(:,:,:,:)          ::  r
  real, intent(in), dimension(:,:)              ::  u_star, b_star, q_star
  real, intent(inout), dimension(:,:,:)         ::  tdt, qdt, udt, vdt
  real, intent(inout), dimension(:,:,:,:)       ::  rdt
  real, intent(inout), dimension(:,:,:,:)       ::  rdiag
  real, intent(out), dimension(:,:,:)           ::  dcond_ls_liquid
  real, intent(out), dimension(:,:,:)           ::  dcond_ls_ice
  real, intent(out), dimension(:,:,:)           ::  Ndrop_act_clubb
  real, intent(out), dimension(:,:,:)           ::  Icedrop_act_clubb
  real, intent(out), dimension(:,:,:)           ::  ndust, rbar_dust
  real, intent(out), dimension(:,:,:)           ::  diff_t_clubb
  real, intent(out), optional, dimension(:,:,:) ::  qcvar_clubb   
  type(aerosol_type), intent(in), optional      ::  Aerosol
  real, intent(in), optional, dimension(:,:,:)  ::  mask
  real, intent(in), optional, dimension(:,:,:)  ::  mc_full
  real, intent(in), optional, dimension(:,:,:)  ::  conv_frac_clubb
  real, intent(in), optional, dimension(:,:,:)  ::  convective_humidity_ratio_clubb
  real, intent(in), optional, dimension(:,:)    ::  tdt_shf,  qdt_lhf

  call error_mesg ('clubb_driver_mod', 'Not compiled with -DCLUBB', FATAL)
  end subroutine clubb

  !=====================================================================
  !=====================================================================

  subroutine clubb_init(id, jd, kd, lon, lat, axes, Time, phalf )

  !---------------------------------------------------------------------
  ! Initialize and set-up clubb
  !
  !  id, jd, kd                         input
  !    subdomain dimensions
  !
  !  lon, lat                           input
  !    longitude, latitude in radians
  !
  !  axes                               input
  !    axis indices, (/x,y,pf,ph/)
  !    (returned from diag axis manager)
  !
  !  Time                               input
  !    current time (time_type)
  !
  !  phalf                              input
  !    pressure at half levels in pascals
  !    [real, dimension(nlon,nlat,nlev+1)]
  !
  !---------------------------------------------------------------------

  ! ----- Calling arguments -----

  integer, intent(in)                  :: id, jd, kd
  real, dimension(:,:), intent(in)     :: lon, lat
  integer, dimension(4), intent(in)    :: axes
  type(time_type), intent(in)          :: Time
  real, dimension(:,:,:), intent(in)   :: phalf
 
  call error_mesg ('clubb_driver_mod', 'Not compiled with -DCLUBB', FATAL)
  end subroutine clubb_init
  !=====================================================================

  !=====================================================================
  subroutine clubb_setup(id, jd, phalf)

  !---------------------------------------------------------------------
  !  id, jd                             input
  !    subdomain dimensions
  !
  !  phalf                              input
  !    pressure at half levels in pascals
  !    [real, dimension(nlon,nlat,nlev+1)]
  !
  !---------------------------------------------------------------------

  ! ----- Calling arguments -----

  integer, intent(in)                 :: id, jd
  real, dimension(:,:,:), intent(in)  ::  phalf
  call error_mesg ('clubb_driver_mod', 'Not compiled with -DCLUBB', FATAL)
  end subroutine clubb_setup
  !=====================================================================

  !=====================================================================
  subroutine clubb_end
  call error_mesg ('clubb_driver_mod', 'Not compiled with -DCLUBB', FATAL)
  end subroutine clubb_end
  !=====================================================================
#endif
end module  clubb_driver_mod
