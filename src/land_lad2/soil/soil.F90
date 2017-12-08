! ============================================================================
! soil model module
! ============================================================================
module soil_mod

#include "../shared/debug.inc"

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only: error_mesg, file_exist, check_nml_error, &
     stdlog, close_file, mpp_pe, mpp_root_pe, FATAL, WARNING, NOTE
use time_manager_mod,   only: time_type_to_real
use diag_manager_mod,   only: diag_axis_init
use constants_mod,      only: tfreeze, hlv, hlf, dens_h2o
use tracer_manager_mod, only: NO_TRACER

use land_constants_mod, only : NBANDS, BAND_VIS, BAND_NIR, seconds_per_year
use soil_tile_mod, only : GW_LM2, GW_LINEAR, GW_HILL_AR5, GW_HILL, GW_TILED, &
     soil_tile_type, soil_pars_type, read_soil_data_namelist, &
     soil_data_thermodynamics, &
     soil_data_hydraulic_properties, soil_data_psi_for_rh, &
     soil_data_gw_hydraulics, soil_data_gw_hydraulics_ar5, &
     soil_data_vwc_for_init_only, &
     soil_data_init_derive_subsurf_pars, &
     soil_data_init_derive_subsurf_pars_ar5, soil_data_init_derive_subsurf_pars_tiled, &
     max_lev, psi_wilt, cpw, clw, csw, g_iso, g_vol, g_geo, use_brdf, g_RT, aspect, &
     gw_scale_length, gw_scale_relief, gw_scale_soil_depth, &
     slope_exp, gw_scale_perm, k0_macro_x, retro_a0n1, &
     soil_type_file, &
     soil_tile_stock_pe, initval, comp, soil_theta, soil_ice_porosity

use soil_carbon_mod, only: poolTotalCarbon, soilMaxCohorts, &
     update_pool, add_litter, add_carbon_to_cohorts, &
     carbon_leaching_with_litter,transfer_pool_fraction, n_c_types, &
     soil_carbon_option, SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE, &
     C_CEL, C_LIG, C_MIC, A_function, debug_pool, adjust_pool_ncohorts

use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_enum_type, &
     first_elmt, prev_elmt, loop_over_tiles
use land_utils_mod, only : put_to_tiles_r0d_fptr, put_to_tiles_r1d_fptr
use land_tile_diag_mod, only : diag_buff_type, &
     register_tiled_static_field, register_tiled_diag_field, &
     send_tile_data, send_tile_data_r0d_fptr, send_tile_data_r1d_fptr, &
     send_tile_data_i0d_fptr, &
     add_tiled_diag_field_alias, add_tiled_static_field_alias, &
     set_default_diag_filter, cmor_name, cmor_mrsos_depth
use land_data_mod, only : lnd, lnd_sg, log_version
use land_io_mod, only : read_field
use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_tile_data, add_int_tile_data, get_tile_data, get_int_tile_data, &
     add_restart_axis, field_exists
use vegn_data_mod, only: K1, K2
use vegn_tile_mod, only : vegn_tile_type, vegn_uptake_profile, vegn_hydraulic_properties
use land_debug_mod, only : is_watch_point, get_current_point, set_current_point, &
     check_conservation, check_var_range, is_watch_cell
use uptake_mod, only : UPTAKE_LINEAR, UPTAKE_DARCY2D, UPTAKE_DARCY2D_LIN, &
     uptake_init, darcy2d_uptake, darcy2d_uptake_solver

use hillslope_mod, only : do_hillslope_model, max_num_topo_hlsps, &
     num_vertclusters, hlsp_coldfracs, use_geohydrodata, & !pond, &
     horiz_wt_depth_to_init, calculate_wt_init, simple_inundation
use land_io_mod, only : &
     init_cover_field
use soil_tile_mod, only : n_dim_soil_types, soil_to_use, &
     soil_index_constant, input_cover_types
use hillslope_hydrology_mod, only: hlsp_hydro_lev_init, hlsp_hydrology_2, &
     stiff_explicit_gwupdate
use river_mod, only : river_tracer_index

! Test tridiagonal solution for advection
use land_numerics_mod, only : tridiag
implicit none
private

! ==== public interfaces =====================================================
public :: read_soil_namelist
public :: soil_cover_cold_start
public :: retrieve_soil_tags
public :: soil_init
public :: soil_end
public :: save_soil_restart

public :: soil_sfc_water
public :: soil_evap_limits
public :: soil_step_1
public :: soil_step_2
public :: soil_step_3
public :: Dsdt
public :: add_root_litter
public :: add_root_exudates
public :: redistribute_peat_carbon
! =====end of public interfaces ==============================================



! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'soil'
#include "../shared/version_variable.inc"

! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
logical :: lm2                  = .false.
logical :: use_E_min            = .false.     ! prohibit condensation
logical :: use_E_max            = .true.      ! theoretical effiltration capacity flag
real    :: init_temp            = 288.        ! cold-start soil T
real    :: init_w               = 150.        ! cold-start w(l)/dz(l)
real    :: init_wtdep           = -1.         ! positive value activates hydrostatic IC,
                                              ! overriding init_w
real    :: init_groundwater     =   0.        ! cold-start gw storage
real    :: lrunf_ie_min         = -1.0e-4     ! trigger for clip and runoff
real    :: lrunf_ie_tol         =  1.e-12
character(len=16) :: albedo_to_use = ''       ! or 'albedo-map' or 'brdf-maps'
character(len=24) :: uptake_to_use = 'linear' ! or 'darcy2d', or 'darcy2d-linearized'
logical :: uptake_oneway        = .false.     ! if true, roots cannot loose water to soil
logical :: uptake_from_sat      = .true.      ! if false, the uptake from saturated soil is prohibited
logical :: allow_negative_rie   = .false.
logical :: baseflow_where_frozen = .false.
logical :: write_when_flagged   = .false.
logical :: corrected_lm2_gw     = .true.
logical :: use_fringe           = .false.
logical :: push_down_sfc_excess = .true.
logical :: lrunf_from_div       = .true.
logical :: use_tall_plants      = .false.
logical :: cold_infilt          = .false.
logical :: bottom_up_cold_infilt= .false.
logical :: use_depth_to_wt_4    = .false.
logical :: always_use_bsw       = .false.
logical :: harmonic_mean_K      = .false.
logical :: verbose              = .false.
logical :: no_min_Dpsi          = .false.
logical :: div_bug_fix          = .false.
logical :: require_pos_hcap     = .false.
logical :: use_richards_clean   = .false.
real    :: active_layer_drainage_acceleration = 0.
real    :: hlf_factor           = 1.
real    :: gw_flux_max          = 1.e10
real    :: aquifer_heat_cap     = 0.         ! in equivalent liquid water amount, kg/m2
logical :: use_tridiag_foradvec = .false.    ! use tridiagonal solution for advection
                        ! information for soil carbon acceleration
logical :: write_soil_carbon_restart = .FALSE. ! indicates whether to write
                        ! information for soil carbon acceleration
logical :: horiz_init_wt        = .false.   ! initialize horizontal water table, if gw_option == GW_TILED
logical :: use_coldstart_wtt_data = .false. ! read additional data for soil initialization
character(len=256)  :: coldstart_datafile = 'INPUT/soil_wtt.nc'
logical :: allow_neg_rnu        = .false.   ! Refill from stream if wl < 0 with warning, i.e. during spinup.
logical :: allow_neg_wl         = .false.   ! Warn rather than abort if wl < 0, even if .not. allow_neg_rnu
logical :: fix_neg_subsurface_wl       = .false.
logical :: prohibit_negative_water_div = .false. ! if TRUE, div_bf abd dif_if are
  ! set to zero in case water content of *any* layer is negative
real    :: zeta_bar_override    = -1.
real    :: cold_depth           = 0.
real    :: Wl_min               = -1.e20
real    :: bwood_macinf         = -1.
integer :: max_iter_trans = 100 ! max number of iterations for psi_crown_min
integer :: layer_for_gw_switch = 1000000 ! to accelerate permafrost gw shutoff
real    :: eps_trans     = 1.e-7 ! convergence crit for psi_crown_min
logical :: supercooled_rnu = .true. ! Excess ice converted to supercooled water for runoff.
real    :: wet_depth = 0.6 ! [m] water table depth threshold for diagnosing wetland fraction
real    :: thetathresh = -0.01 ! [-] threshold for negative soil liquid water saturation
  ! before warning or abort
real    :: negrnuthresh = -0.1 ! [mm/s] threshold for negative lrunf_nu
  ! before warning or abort

real :: max_soil_C_density = 50.0   !(kgC/m3) -- for redistribution of peat
real :: max_litter_thickness = 0.05 ! m of litter layer thickness before it gets redistributed

namelist /soil_nml/ lm2, use_E_min, use_E_max,           &
                    init_temp,      &
                    init_w,   init_wtdep,    &
                    init_groundwater, lrunf_ie_min, lrunf_ie_tol, &
                    cpw, clw, csw, &
                    albedo_to_use, &
                    uptake_to_use, uptake_oneway, uptake_from_sat, &
                    allow_negative_rie, &
                    baseflow_where_frozen, &
                    write_when_flagged, &
                    corrected_lm2_gw, &
                    use_fringe, &
                    use_tall_plants, cold_infilt, bottom_up_cold_infilt, &
                    use_depth_to_wt_4, always_use_bsw, &
                    harmonic_mean_K, verbose, no_min_Dpsi, div_bug_fix, &
                    require_pos_hcap, use_richards_clean, &
                    push_down_sfc_excess, lrunf_from_div, &
                    active_layer_drainage_acceleration, hlf_factor, &
                    gw_flux_max, aquifer_heat_cap, use_tridiag_foradvec, &
                    horiz_init_wt, use_coldstart_wtt_data, coldstart_datafile, &
                    allow_neg_rnu, allow_neg_wl, fix_neg_subsurface_wl, prohibit_negative_water_div, &
                    zeta_bar_override, &
                    cold_depth, Wl_min, &
                    bwood_macinf, &
                    max_iter_trans, layer_for_gw_switch, eps_trans, &
                    supercooled_rnu, wet_depth, thetathresh, negrnuthresh, &
                    write_soil_carbon_restart, &
                    max_soil_C_density, max_litter_thickness
!---- end of namelist --------------------------------------------------------

logical         :: module_is_initialized =.FALSE.
real            :: delta_time ! fast (physical) time step, s
real            :: dt_fast_yr ! fast (physical) time step, yr (year is defined as 365 days)
logical         :: use_single_geo
integer         :: num_l              ! # of water layers
real            :: dz    (max_lev)    ! thicknesses of layers
real            :: zfull (max_lev)
real            :: zhalf (max_lev+1)
real            :: Eg_min

integer         :: uptake_option = -1
integer         :: gw_option = -1

integer :: n_river_tracers = 0
integer :: i_river_DOC     = NO_TRACER

! ---- diagnostic field IDs
integer :: id_fast_soil_C, id_slow_soil_C, id_protected_C, id_fsc, id_ssc,&
    id_leaflitter_deadmic, id_leaflitter_livemic, id_leaflitter_fast_C, id_leaflitter_slow_C, id_nleaflittercohorts, &
    id_finewoodlitter_deadmic, id_finewoodlitter_livemic, id_finewoodlitter_fast_C, id_finewoodlitter_slow_C, id_nfinewoodlittercohorts, &
    id_coarsewoodlitter_deadmic, id_coarsewoodlitter_livemic, id_coarsewoodlitter_fast_C, id_coarsewoodlitter_slow_C, id_ncoarsewoodlittercohorts, &
    id_deadmic, id_livemic, id_nsoilcohorts, id_Qmax, id_protectedC,  id_deadmic_total, id_livemic_total,&
    id_total_soil_C,id_dissolved_total,id_coarseWoodlitter_total_C,id_fineWoodlitter_total_C,id_leaflitter_total_C,&
    id_total_carbon_layered,&
    id_lwc, id_swc, id_psi, id_temp, &
    id_ie, id_sn, id_bf, id_if, id_al, id_nu, id_sc, &
    id_hie, id_hsn, id_hbf, id_hif, id_hal, id_hnu, id_hsc, &
    id_heat_cap, id_thermal_cond, id_type, id_tau_gw, id_slope_l, &
    id_slope_Z, id_zeta_bar, id_e_depth, id_vwc_sat, id_vwc_fc, &
    id_vwc_wilt, id_K_sat, id_K_gw, id_w_fc, id_alpha, &
    id_refl_dry_dif, id_refl_dry_dir, id_refl_sat_dif, id_refl_sat_dir, &
    id_f_iso_dry, id_f_vol_dry, id_f_geo_dry, &
    id_f_iso_sat, id_f_vol_sat, id_f_geo_sat, &
    id_evap, id_uptk_n_iter, id_uptk, id_psi_x0, id_uptk_residual, &
    id_sws_n_iter, id_psi_x0_sws, &
    id_excess, id_deficit, id_deficit_2, id_deficit_3, id_deficit_4, id_zeta, id_tau, &
    id_psi_bot, id_sat_frac, id_stor_frac, id_sat_depth, id_sat_dept2, &
    id_cf_1, id_cf_3, id_wt_1, id_wt_2, id_wt_2a, id_wt_2b, id_wt_3, id_wt2_3, id_wt_4, &
    id_div_bf, id_div_if, id_div_al, &
    id_z_cap, id_active_layer, id_surface_water, id_inun_frac, id_rsn_frac, id_flow, id_reflux, &
    id_fast_C_leaching,id_slow_C_leaching,id_livemic_C_leaching,id_deadmic_C_leaching,&
    id_protected_C_leaching,id_protected_total,&
    id_fast_dissolved_C,id_slow_dissolved_C,id_deadmic_dissolved_C,&
    id_leaflitter_fast_dissolved_C,id_leaflitter_slow_dissolved_C,id_leaflitter_deadmic_dissolved_C,&
    id_finewoodlitter_fast_dissolved_C,id_finewoodlitter_slow_dissolved_C,id_finewoodlitter_deadmic_dissolved_C,&
    id_coarsewoodlitter_fast_dissolved_C,id_coarsewoodlitter_slow_dissolved_C,id_coarsewoodlitter_deadmic_dissolved_C,&
    id_rsoil_leaflitter_deadmic, id_rsoil_leaflitter_fast, id_rsoil_leaflitter_slow, &
    id_rsoil_finewoodlitter_deadmic, id_rsoil_finewoodlitter_fast, id_rsoil_finewoodlitter_slow, &
    id_rsoil_coarsewoodlitter_deadmic, id_rsoil_coarsewoodlitter_fast, id_rsoil_coarsewoodlitter_slow, &
    id_rsoil_fast, id_rsoil_slow,id_rsoil_deadmic,id_asoil,id_rsoil,&
    id_leaflitter_dissolved_fast,id_leaflitter_dissolved_slow,id_leaflitter_dissolved_deadmic,&
    id_leaflitter_deposited_fast,id_leaflitter_deposited_slow,id_leaflitter_deposited_deadmic,&
        id_finewoodlitter_dissolved_fast,id_finewoodlitter_dissolved_slow,id_finewoodlitter_dissolved_deadmic,&
    id_finewoodlitter_deposited_fast,id_finewoodlitter_deposited_slow,id_finewoodlitter_deposited_deadmic,&
        id_coarsewoodlitter_dissolved_fast,id_coarsewoodlitter_dissolved_slow,id_coarsewoodlitter_dissolved_deadmic,&
    id_coarsewoodlitter_deposited_fast,id_coarsewoodlitter_deposited_slow,id_coarsewoodlitter_deposited_deadmic,&
    id_dissolved_fast,id_dissolved_slow,id_dissolved_deadmic,&
    id_deposited_fast,id_deposited_slow,id_deposited_deadmic, &
    id_leaflitter_fast_C_leaching,id_leaflitter_slow_C_leaching,id_leaflitter_deadmic_C_leaching,&
    id_finewoodlitter_fast_C_leaching,id_finewoodlitter_slow_C_leaching,id_finewoodlitter_deadmic_C_leaching,&
    id_coarsewoodlitter_fast_C_leaching,id_coarsewoodlitter_slow_C_leaching,id_coarsewoodlitter_deadmic_C_leaching,&
    id_fast_DOC_div_loss,id_slow_DOC_div_loss,id_deadmic_DOC_div_loss, &
    id_slomtot, id_wet_frac, id_macro_infilt, &
    id_surf_DOC_loss, id_total_C_leaching, id_total_DOC_div_loss

! diag IDs of CMOR variables
integer :: id_mrlsl, id_mrsfl, id_mrsll, id_mrsol, id_mrso, id_mrsos, id_mrlso, id_mrfso, &
    id_mrs1mLut, id_mrro, id_mrros, id_csoil, id_rh, &
    id_csoilfast, id_csoilmedium, id_csoilslow

! test tridiagonal solver for advection
integer :: id_st_diff

integer, allocatable, dimension(:,:), private :: soil_tags ! module copy of soil tags for cold start
real, allocatable :: mrsos_weight(:) ! weights for mrsos averaging
real, allocatable :: mrs1m_weight(:) ! weights for mrs1m averaging
! ==== end of module variables ===============================================

contains

! ============================================================================
subroutine read_soil_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: l

  call read_soil_data_namelist(num_l,dz,use_single_geo,gw_option)

  call log_version(version, module_name, &
  __FILE__)
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=soil_nml, iostat=io)
  ierr = check_nml_error(io, 'soil_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=soil_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'soil_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=soil_nml)
  endif

  ! ---- set up vertical discretization
  zhalf(1) = 0
  do l = 1, num_l;
     zhalf(l+1) = zhalf(l) + dz(l)
     zfull(l) = 0.5*(zhalf(l+1) + zhalf(l))
  enddo

  ! ---- convert symbolic names of options into numeric IDs to speed up
  ! the selection during run-time
  if (trim(uptake_to_use)=='linear') then
     uptake_option = UPTAKE_LINEAR
  else if (trim(uptake_to_use)=='darcy2d') then
     uptake_option = UPTAKE_DARCY2D
  else if (trim(uptake_to_use)=='darcy2d-linearized') then
     uptake_option = UPTAKE_DARCY2D_LIN
  else
     call error_mesg('soil_init',&
          'soil uptake option uptake_to_use="'//&
          trim(uptake_to_use)//'" is invalid, use "linear", "darcy2d" or "darcy2d-linearized"',&
          FATAL)
  endif

  if (use_E_min) then
     Eg_min = 0.
  else
     Eg_min = -HUGE(Eg_min)
  endif

  ! Configuration checking
  if (gw_option == GW_TILED .and. .not. use_tridiag_foradvec) then
     call error_mesg(module_name, 'read_soil_namelist: over-riding "use_tridiag_foradvec" value set '// &
                     'in namelist or module default. Tiled groundwater model / full hillslope model requires '// &
                     '"use_tridiag_foradvec" == .true.', NOTE)
     use_tridiag_foradvec = .true.
  end if

end subroutine read_soil_namelist


! ============================================================================
! initialize soil model
subroutine soil_init (predefined_tiles, id_ug,id_band,id_zfull)
  logical,intent(in)  :: predefined_tiles !<If true, this subroutine assumes that 
        ! the distribution of basic soil parameters is set from external data set, 
        ! and only initializes derived soil properties
  integer,intent(in)  :: id_ug    !<Unstructured axis id.
  integer,intent(in)  :: id_band  !<ID of spectral band axis
  integer,intent(out) :: id_zfull !<ID of vertical soil axis

  ! ---- local vars
  type(land_tile_enum_type)     :: ce   ! tile list enumerator
  type(land_tile_type), pointer :: tile ! pointer to current tile
  ! input data buffers for respective variables:
  real, allocatable :: gw_param(:), gw_param2(:), gw_param3(:), albedo(:,:)
  real, allocatable :: f_iso(:,:), f_vol(:,:), f_geo(:,:), refl_dif(:,:)

  real :: local_wt_depth ! [m] water table depth for tile (+ for below surface)
  real, allocatable :: ref_soil_t(:) ! reference soil temperature (based on 5 m or surface air temperature)
                                     ! for cold-start initialization
  real, allocatable :: wetmask(:)    ! input mask for zones with high water table
  logical :: drypoint                ! This point is predicted to have a falling water table.

  integer :: i, k, ll ! indices
  real :: psi(num_l), mwc(num_l)

  type(land_restart_type) :: restart, restart1
  logical :: restart_exists
  character(*), parameter :: restart_file_name = 'INPUT/soil.res.nc'

  module_is_initialized = .TRUE.
  delta_time = time_type_to_real(lnd%dt_fast)
  dt_fast_yr = delta_time/seconds_per_year

  call uptake_init(num_l,dz,zfull)
  call hlsp_hydro_lev_init(num_l,dz,zfull)

  ! initialize river tracer indices
  i_river_DOC  = river_tracer_index('doc')

  ! -------- initialize soil model diagnostic fields
  call soil_diag_init(id_ug,id_band,id_zfull)

  if (predefined_tiles) then
     ! Initialize the soil derived parameters
     ce = first_elmt(land_tile_map)
     do while(loop_over_tiles(ce,tile))
        if (.not.associated(tile%soil)) cycle
        call soil_data_init_derive_subsurf_pars_tiled(tile%soil, use_geohydrodata)
     end do
  else if (use_single_geo) then
     if (gw_option == GW_TILED) then ! and use_single_geo
        ! Error checking
        if (.not. use_geohydrodata) then
           call error_mesg(module_name, 'soil_init: incompatible namelist options selected. gw_option =='// &
                           ' tiled, use_geohydrodata == .false., and use_single_geo == .true.', FATAL)
        else
           call error_mesg(module_name, 'soil_init: Warning: using tiled hillslope groundwater model '// &
                           'with single global values for soil hydrological properties (i.e. "use_single_geo).', &
                           WARNING)
        endif
     endif
  else 
     ! NOT predefined soil parameters and NOT use_single_geo:
     ! read spatially distributed fields for groundwater parameters, as requested
     select case (gw_option)
     case (GW_LINEAR,GW_LM2)
        allocate(gw_param(lnd%ls:lnd%le))
        call read_field( 'INPUT/groundwater_residence.nc','tau', lnd%lon, lnd%lat, &
             gw_param, interp='bilinear' )
        call put_to_tiles_r0d_fptr( gw_param, land_tile_map, soil_tau_groundwater_ptr )
        deallocate(gw_param)
     case (GW_HILL, GW_HILL_AR5)
        allocate(gw_param (lnd%ls:lnd%le))
        allocate(gw_param2(lnd%ls:lnd%le))
        allocate(gw_param3(lnd%ls:lnd%le))
        call read_field( 'INPUT/geohydrology.nc','hillslope_length',  lnd%lon, lnd%lat, &
          gw_param, interp='bilinear' )
        call put_to_tiles_r0d_fptr( gw_param*gw_scale_length, land_tile_map, soil_hillslope_length_ptr )
        call read_field( 'INPUT/geohydrology.nc','slope', lnd%lon, lnd%lat, &
          gw_param2, interp='bilinear' )
        gw_param = gw_param*gw_param2
        call put_to_tiles_r0d_fptr( gw_param*gw_scale_relief, land_tile_map, soil_hillslope_relief_ptr )

        if (retro_a0n1 .or. gw_option.eq.GW_HILL_AR5) then
            gw_param = 0.
            call put_to_tiles_r0d_fptr( gw_param, land_tile_map, soil_hillslope_a_ptr )
            gw_param = 1.
            call put_to_tiles_r0d_fptr( gw_param, land_tile_map, soil_hillslope_n_ptr )
!            call read_field( 'INPUT/geohydrology.nc','hillslope_zeta_bar', &
!              lnd_sg%lon, lnd_sg%lat, gw_param, interp='bilinear' )
            gw_param = 0.5
            call put_to_tiles_r0d_fptr( gw_param, land_tile_map, soil_hillslope_zeta_bar_ptr )
        else
            call read_field( 'INPUT/geohydrology.nc','hillslope_a', &
              lnd%lon, lnd%lat, gw_param, interp='bilinear' )
            call put_to_tiles_r0d_fptr( gw_param, land_tile_map, soil_hillslope_a_ptr )
            call read_field( 'INPUT/geohydrology.nc','hillslope_n', &
              lnd%lon, lnd%lat, gw_param2, interp='bilinear' )
            call put_to_tiles_r0d_fptr( gw_param2, land_tile_map, soil_hillslope_n_ptr )
            gw_param3 = (1./(gw_param2+1.)+gw_param/(gw_param2+2.))/(1.+gw_param/2.)
            call put_to_tiles_r0d_fptr( gw_param3, land_tile_map, soil_hillslope_zeta_bar_ptr )
        endif

        call read_field( 'INPUT/geohydrology.nc','soil_e_depth', &
          lnd%lon, lnd%lat, gw_param, interp='bilinear' )
        if (slope_exp.gt.0.01) then
            call put_to_tiles_r0d_fptr( gw_param*gw_scale_soil_depth*(0.08/gw_param2)**slope_exp, &
                                                  land_tile_map, soil_soil_e_depth_ptr )
        else
            call put_to_tiles_r0d_fptr( gw_param*gw_scale_soil_depth, land_tile_map, soil_soil_e_depth_ptr )
        endif
        if (gw_option /= GW_HILL_AR5) then
            call read_field( 'INPUT/geohydrology.nc','perm', lnd%lon, lnd%lat, &
                 gw_param, interp='bilinear' )
            call put_to_tiles_r0d_fptr(9.8e9*gw_scale_perm*gw_param, land_tile_map, &
                                            soil_k_sat_gw_ptr )
        endif
        deallocate(gw_param, gw_param2, gw_param3)
        ce = first_elmt(land_tile_map)
        do while(loop_over_tiles(ce,tile))
            if (.not.associated(tile%soil)) cycle
            select case (gw_option)
            case (GW_HILL)
                call soil_data_init_derive_subsurf_pars(tile%soil)
            case (GW_HILL_AR5)
                call soil_data_init_derive_subsurf_pars_ar5(tile%soil)
            end select
        enddo
     case (GW_TILED)
        if (use_geohydrodata) then
           allocate(gw_param (lnd%ls:lnd%le), gw_param2(lnd%ls:lnd%le))
           call read_field( 'INPUT/geohydrology.nc','hillslope_length',  lnd%lon, lnd%lat, &
             gw_param, interp='bilinear' )
           call put_to_tiles_r0d_fptr( gw_param*gw_scale_length, land_tile_map, soil_hillslope_length_ptr )
           call read_field( 'INPUT/geohydrology.nc','slope', lnd%lon, lnd%lat, &
             gw_param2, interp='bilinear' )
           gw_param = gw_param*gw_param2
           call put_to_tiles_r0d_fptr( gw_param*gw_scale_relief, land_tile_map, soil_hillslope_relief_ptr )
           call read_field( 'INPUT/geohydrology.nc','hillslope_zeta_bar', &
             lnd%lon, lnd%lat, gw_param, interp='bilinear' )
           if (zeta_bar_override.gt.0.) gw_param=zeta_bar_override
           call put_to_tiles_r0d_fptr( gw_param, land_tile_map, soil_hillslope_zeta_bar_ptr )
           call read_field( 'INPUT/geohydrology.nc','soil_e_depth', &
             lnd%lon, lnd%lat, gw_param, interp='bilinear' )

           if (slope_exp.gt.0.01) then
           ! ZMS It's probably inconsistent to leave in this if statement.
               call error_mesg(module_name, 'soil_init: "slope_exp" > 0.0 requested even though '// &
                               'running with tiled hillslope model.  This may be inconsistent.', WARNING)
               call put_to_tiles_r0d_fptr( gw_param*gw_scale_soil_depth*(0.08/gw_param2)**slope_exp, &
                                                     land_tile_map, soil_soil_e_depth_ptr )
           else
               call put_to_tiles_r0d_fptr( gw_param*gw_scale_soil_depth, land_tile_map, soil_soil_e_depth_ptr )
           endif
           call read_field( 'INPUT/geohydrology.nc','perm', lnd%lon, lnd%lat, &
                  gw_param, interp='bilinear' )
           call put_to_tiles_r0d_fptr(9.8e9*gw_scale_perm*gw_param, land_tile_map, &
                                          soil_k_sat_gw_ptr )
           deallocate(gw_param, gw_param2)
        end if
        ce = first_elmt(land_tile_map)
        do while(loop_over_tiles(ce,tile))
            if (.not.associated(tile%soil)) cycle
            call soil_data_init_derive_subsurf_pars_tiled(tile%soil, use_geohydrodata)
        end do
     end select ! gw_option
  endif

  ! -------- set dry soil albedo values, if requested
  if (trim(albedo_to_use)=='albedo-map') then
     allocate(albedo(lnd%ls:lnd%le,NBANDS))
     call read_field( 'INPUT/soil_albedo.nc','SOIL_ALBEDO_VIS',&
          lnd%lon, lnd%lat, albedo(:,BAND_VIS),'bilinear')
     call read_field( 'INPUT/soil_albedo.nc','SOIL_ALBEDO_NIR',&
          lnd%lon, lnd%lat, albedo(:,BAND_NIR),'bilinear')
     call put_to_tiles_r1d_fptr( albedo, land_tile_map, soil_refl_dry_dir_ptr )
     call put_to_tiles_r1d_fptr( albedo, land_tile_map, soil_refl_dry_dif_ptr )
     ! for now, put the same value into the saturated soil albedo, so that
     ! the albedo does not depend on soil wetness
     call put_to_tiles_r1d_fptr( albedo, land_tile_map, soil_refl_sat_dir_ptr )
     call put_to_tiles_r1d_fptr( albedo, land_tile_map, soil_refl_sat_dif_ptr )
     deallocate(albedo)
  else if (trim(albedo_to_use)=='brdf-maps') then
     use_brdf = .true.
     allocate(   f_iso(lnd%ls:lnd%le,NBANDS))
     allocate(   f_vol(lnd%ls:lnd%le,NBANDS))
     allocate(   f_geo(lnd%ls:lnd%le,NBANDS))
     allocate(refl_dif(lnd%ls:lnd%le,NBANDS))
     call read_field( 'INPUT/soil_brdf.nc','f_iso_vis',&
          lnd%lon, lnd%lat, f_iso(:,BAND_VIS),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_vol_vis',&
          lnd%lon, lnd%lat, f_vol(:,BAND_VIS),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_geo_vis',&
          lnd%lon, lnd%lat, f_geo(:,BAND_VIS),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_iso_nir',&
          lnd%lon, lnd%lat, f_iso(:,BAND_NIR),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_vol_nir',&
          lnd%lon, lnd%lat, f_vol(:,BAND_NIR),'bilinear')
     call read_field( 'INPUT/soil_brdf.nc','f_geo_nir',&
          lnd%lon, lnd%lat, f_geo(:,BAND_NIR),'bilinear')
     refl_dif = g_iso*f_iso + g_vol*f_vol + g_geo*f_geo
     call put_to_tiles_r1d_fptr( f_iso,    land_tile_map, soil_f_iso_dry_ptr )
     call put_to_tiles_r1d_fptr( f_vol,    land_tile_map, soil_f_vol_dry_ptr )
     call put_to_tiles_r1d_fptr( f_geo,    land_tile_map, soil_f_geo_dry_ptr )
     call put_to_tiles_r1d_fptr( refl_dif, land_tile_map, soil_refl_dry_dif_ptr )
     ! for now, put the same value into the saturated soil albedo, so that
     ! the albedo does not depend on soil wetness
     call put_to_tiles_r1d_fptr( f_iso,    land_tile_map, soil_f_iso_sat_ptr )
     call put_to_tiles_r1d_fptr( f_vol,    land_tile_map, soil_f_vol_sat_ptr )
     call put_to_tiles_r1d_fptr( f_geo,    land_tile_map, soil_f_geo_sat_ptr )
     call put_to_tiles_r1d_fptr( refl_dif, land_tile_map, soil_refl_sat_dif_ptr )
     deallocate(f_iso, f_vol, f_geo, refl_dif)
  else if (trim(albedo_to_use)=='') then
     ! do nothing, that is leave soil albedo parameters as defined based on the data table
  else
     call error_mesg('soil_init',&
          'option albedo_to_use="'// trim(albedo_to_use)//&
          '" is invalid, use "albedo-map", "brdf-maps", or empty line ("")',&
          FATAL)
  endif

  ! Call calculate_wt_init outside tile loop so that it is done once per hillslope
  if (init_wtdep .gt. 0. .and. gw_option == GW_TILED) then
     call calculate_wt_init(init_wtdep)
  end if

  if (use_coldstart_wtt_data) then
     allocate(ref_soil_t(lnd%ls:lnd%le), wetmask(lnd%ls:lnd%le))
     call read_field( coldstart_datafile, 'REFSOILT', &
             lnd%lon, lnd%lat, ref_soil_t, interp='bilinear' )
     call read_field( coldstart_datafile, 'WETMASK', &
             lnd%lon, lnd%lat, wetmask, interp='bilinear' )
  end if

  ! -------- initialize soil state --------
  ce = first_elmt(land_tile_map, ls=lnd%ls) ! Use global indices here because element indices
                                            ! needed.
  do while(loop_over_tiles(ce,tile,k=k,l=ll))
     if (.not.associated(tile%soil)) cycle
     call set_current_point(ll,k)
     if (init_wtdep .gt. 0.) then
        if (.not. use_coldstart_wtt_data) then
           if (horiz_init_wt .and. gw_option == GW_TILED) then
              call horiz_wt_depth_to_init(tile%soil, prev_elmt(ce), local_wt_depth)
              ! Note: if restart_exists, then this function returns dummy local_wt_depth == 0.
              ! prev_elmt(ce) passed because indices will be needed.
              psi = zfull(1:num_l) - local_wt_depth
           else
              psi = zfull(1:num_l) - init_wtdep
           end if
        else
           if (wetmask(ll) > 0.5) then ! wet point
              drypoint = .false.
           else
              drypoint = .true.
           end if
           if (gw_option == GW_TILED) then
              call horiz_wt_depth_to_init(tile%soil, prev_elmt(ce), local_wt_depth, dry=drypoint)
              psi = zfull(1:num_l) - local_wt_depth
           else if (drypoint) then
              psi = zfull(1:num_l) - tile%soil%pars%hillslope_relief*tile%soil%pars%hillslope_zeta_bar
           else
              psi = zfull(1:num_l) - init_wtdep
           end if
        end if
        call soil_data_vwc_for_init_only(tile%soil, psi, mwc)
        mwc = mwc * dens_h2o
     else if (init_w .ge. 0.) then
        mwc = init_w
     else ! negative init_w is to be intrepreted as prescribed saturation
        mwc = -init_w*tile%soil%pars%vwc_sat*dens_h2o
     endif
     if (.not. use_coldstart_wtt_data) then
        if (init_temp.ge.tile%soil%pars%tfreeze) then
           tile%soil%wl = mwc*dz(1:num_l)
           tile%soil%ws = 0
        else
           tile%soil%wl = 0
           tile%soil%ws = mwc*dz(1:num_l)
        endif
        tile%soil%T             = init_temp
        tile%soil%groundwater   = init_groundwater
        tile%soil%groundwater_T = init_temp
        tile%soil%uptake_T           = init_temp
     else
        call init_soil_twc(tile%soil, ref_soil_t(ll), mwc)
     end if
  end do

  if (use_coldstart_wtt_data) then
     deallocate(ref_soil_t, wetmask)
  end if

  call open_land_restart(restart,restart_file_name,restart_exists)
  if (restart_exists) then
     call error_mesg('soil_init', 'reading NetCDF restart "'//trim(restart_file_name)//'"', NOTE)
     call get_tile_data(restart, 'temp', 'zfull', soil_T_ptr  )
     call get_tile_data(restart, 'wl', 'zfull', soil_wl_ptr )
     call get_tile_data(restart, 'ws', 'zfull', soil_ws_ptr )
     call get_tile_data(restart, 'groundwater', 'zfull', soil_groundwater_ptr )
     call get_tile_data(restart, 'groundwater_T', 'zfull', soil_groundwater_T_ptr)
     if(field_exists(restart, 'uptake_T')) &
          call get_tile_data(restart, 'uptake_T', soil_uptake_T_ptr)

     if (soil_carbon_option==SOILC_CENTURY.or.soil_carbon_option==SOILC_CENTURY_BY_LAYER) then
        if(field_exists(restart, 'fsc')) then
           call get_tile_data(restart,'fsc','zfull',soil_fast_soil_C_ptr)
           call get_tile_data(restart,'ssc','zfull',soil_slow_soil_C_ptr)
        else
           ! try to read fsc and ssc from vegetation restart
           call open_land_restart(restart1,'INPUT/vegn2.res.nc',restart_exists)
           if (restart_exists) then
              ! read old (scalar) fsc and ssc into the first element of the fast_soil_C
              ! and slow_soil_C arrays
              call get_tile_data(restart1,'fsc',soil_fast_soil_C_ptr,1)
              call get_tile_data(restart1,'ssc',soil_slow_soil_C_ptr,1)
           endif
           call free_land_restart(restart1)
        endif
     endif

     if (field_exists(restart,'fast_soil_C')) then
        ! we are dealing with CORPSE restart
        ce = first_elmt(land_tile_map)
        do while(loop_over_tiles(ce,tile))
            if (.not.associated(tile%soil)) cycle
            call adjust_pool_ncohorts(tile%soil%leafLitter)
            call adjust_pool_ncohorts(tile%soil%fineWoodLitter)
            call adjust_pool_ncohorts(tile%soil%coarseWoodLitter)
            do i = 1,num_l
               call adjust_pool_ncohorts(tile%soil%soil_C(i))
            enddo
        end do

        call get_tile_data(restart, 'fast_soil_C', 'zfull','soilCCohort', sc_soil_C_ptr, C_CEL)
        call get_tile_data(restart, 'slow_soil_C', 'zfull','soilCCohort', sc_soil_C_ptr, C_LIG)
        call get_tile_data(restart, 'deadMic',     'zfull','soilCCohort', sc_soil_C_ptr, C_MIC)
        call get_tile_data(restart, 'fastProtectedC', 'zfull','soilCCohort', sc_protected_C_ptr, C_CEL)
        call get_tile_data(restart, 'slowProtectedC', 'zfull','soilCCohort', sc_protected_C_ptr, C_LIG)
        call get_tile_data(restart, 'deadMicrobeProtectedC', 'zfull','soilCCohort', sc_protected_C_ptr, C_MIC)

        call get_tile_data(restart,'liveMic', 'zfull','soilCCohort',soilc_livingMicrobeC_ptr)
        call get_tile_data(restart,'CO2', 'zfull','soilCCohort',soilc_CO2_ptr)
        call get_tile_data(restart,'Rtot', 'zfull','soilCCohort',soilc_Rtot_ptr)
        call get_tile_data(restart,'originalCohortC', 'zfull','soilCCohort', soilc_originalLitterC_ptr)

        call get_tile_data(restart, 'soil_DOC_fast', 'zfull',soil_fast_DOC_ptr)
        call get_tile_data(restart, 'soil_DOC_slow', 'zfull',soil_slow_DOC_ptr)
        call get_tile_data(restart, 'soil_DOC_deadmic', 'zfull',soil_deadmicrobe_DOC_ptr)

        if(field_exists(restart, 'fast_DOC_leached')) then
           call get_tile_data(restart,'fast_DOC_leached', soil_fast_DOC_leached_ptr)
           call get_tile_data(restart,'slow_DOC_leached', soil_slow_DOC_leached_ptr)
           call get_tile_data(restart,'deadmic_DOC_leached', soil_deadmic_DOC_leached_ptr)
        endif

        call get_tile_data(restart, 'leaf_litter_fast_C', 'soilCCohort', soilc_leafLitter_litterC_ptr, C_CEL)
        call get_tile_data(restart, 'leaf_litter_slow_C', 'soilCCohort', soilc_leafLitter_litterC_ptr, C_LIG)
        call get_tile_data(restart, 'leaf_litter_deadMic_C', 'soilCCohort', soilc_leafLitter_litterC_ptr, C_MIC)
        call get_tile_data(restart, 'leaf_litter_liveMic_C', 'soilCCohort', soilc_leafLitter_livingMicrobeC_ptr)
        call get_tile_data(restart, 'leaf_litter_CO2', 'soilCCohort', soilc_leafLitter_CO2_ptr)
        call get_tile_data(restart, 'leaf_litter_Rtot', 'soilCCohort', soilc_leafLitter_Rtot_ptr)
        call get_tile_data(restart, 'leaf_litter_originalCohortC', 'soilCCohort',soilc_leafLitter_originalLitterC_ptr)
        call get_tile_data(restart, 'leaf_litter_fastProtectedC', 'soilCCohort',soilc_leafLitter_protectedC_ptr, C_CEL)
        call get_tile_data(restart, 'leaf_litter_slowProtectedC', 'soilCCohort',soilc_leafLitter_protectedC_ptr, C_LIG)
        call get_tile_data(restart, 'leaf_litter_deadMicrobeProtectedC', 'soilCCohort', soilc_leafLitter_protectedC_ptr, C_MIC)

        call get_tile_data(restart, 'leaf_litter_DOC_fast',soilc_leaflitter_DOC_ptr, C_CEL)
        call get_tile_data(restart, 'leaf_litter_DOC_slow',soilc_leaflitter_DOC_ptr, C_LIG)
        call get_tile_data(restart, 'leaf_litter_DOC_deadmic',soilc_leaflitter_DOC_ptr, C_MIC)

        call get_tile_data(restart, 'fineWood_litter_fast_C', 'soilCCohort', soilc_fineWoodLitter_litterC_ptr, C_CEL)
        call get_tile_data(restart, 'fineWood_litter_slow_C', 'soilCCohort', soilc_fineWoodLitter_litterC_ptr, C_LIG)
        call get_tile_data(restart, 'fineWood_litter_deadMic_C', 'soilCCohort' , soilc_fineWoodLitter_litterC_ptr, C_MIC)
        call get_tile_data(restart, 'fineWood_litter_liveMic_C', 'soilCCohort' , soilc_fineWoodLitter_livingMicrobeC_ptr)
        call get_tile_data(restart, 'fineWood_litter_CO2', 'soilCCohort', soilc_fineWoodLitter_CO2_ptr)
        call get_tile_data(restart, 'fineWood_litter_Rtot', 'soilCCohort', soilc_fineWoodLitter_Rtot_ptr)
        call get_tile_data(restart, 'fineWood_litter_originalCohortC', 'soilCCohort',soilc_fineWoodLitter_originalLitterC_ptr)
        call get_tile_data(restart, 'fineWood_litter_fastProtectedC', 'soilCCohort',soilc_fineWoodLitter_protectedC_ptr, C_CEL)
        call get_tile_data(restart, 'fineWood_litter_slowProtectedC', 'soilCCohort',soilc_fineWoodLitter_protectedC_ptr, C_LIG)
        call get_tile_data(restart, 'fineWood_litter_deadMicrobeProtectedC', 'soilCCohort',soilc_fineWoodLitter_protectedC_ptr, C_MIC)

        call get_tile_data(restart, 'fineWood_litter_DOC_fast',soilc_fineWoodlitter_DOC_ptr, C_CEL)
        call get_tile_data(restart, 'fineWood_litter_DOC_slow',soilc_fineWoodlitter_DOC_ptr, C_LIG)
        call get_tile_data(restart, 'fineWood_litter_DOC_deadmic',soilc_fineWoodlitter_DOC_ptr, C_MIC)

        call get_tile_data(restart, 'coarseWood_litter_fast_C', 'soilCCohort', soilc_coarseWoodLitter_litterC_ptr, C_CEL)
        call get_tile_data(restart, 'coarseWood_litter_slow_C', 'soilCCohort', soilc_coarseWoodLitter_litterC_ptr, C_LIG)
        call get_tile_data(restart, 'coarseWood_litter_deadMic_C', 'soilCCohort', soilc_coarseWoodLitter_litterC_ptr, C_MIC)
        call get_tile_data(restart, 'coarseWood_litter_liveMic_C', 'soilCCohort', soilc_coarseWoodLitter_livingMicrobeC_ptr)
        call get_tile_data(restart, 'coarseWood_litter_CO2', 'soilCCohort', soilc_coarseWoodLitter_CO2_ptr)
        call get_tile_data(restart, 'coarseWood_litter_Rtot', 'soilCCohort', soilc_coarseWoodLitter_Rtot_ptr)
        call get_tile_data(restart, 'coarseWood_litter_originalCohortC', 'soilCCohort', soilc_coarseWoodLitter_originalLitterC_ptr)
        call get_tile_data(restart, 'coarseWood_litter_fastProtectedC', 'soilCCohort', soilc_coarseWoodLitter_protectedC_ptr, C_CEL)
        call get_tile_data(restart, 'coarseWood_litter_slowProtectedC', 'soilCCohort', soilc_coarseWoodLitter_protectedC_ptr, C_LIG)
        call get_tile_data(restart, 'coarseWood_litter_deadMicrobeProtectedC', 'soilCCohort', soilc_coarseWoodLitter_protectedC_ptr, C_MIC)

        call get_tile_data(restart, 'coarseWood_litter_DOC_fast',soilc_coarseWoodlitter_DOC_ptr, C_CEL)
        call get_tile_data(restart, 'coarseWood_litter_DOC_slow',soilc_coarseWoodlitter_DOC_ptr, C_LIG)
        call get_tile_data(restart, 'coarseWood_litter_DOC_deadmic',soilc_coarseWoodlitter_DOC_ptr, C_MIC)

        if(field_exists(restart, 'is_peat')) then
           call get_int_tile_data(restart, 'is_peat','zfull', soil_is_peat_ptr)
        endif
     endif
  else
     call error_mesg('soil_init', 'cold-starting soil', NOTE)
  endif
  call free_land_restart(restart)

  ! read soil carbon restart, if present
  call open_land_restart(restart,'INPUT/soil_carbon.res.nc',restart_exists)
  if (restart_exists) then
     call error_mesg('veg_data_init','reading soil_carbon restart',NOTE)
     call get_tile_data(restart,'asoil_in','zfull',soil_asoil_in_ptr)
     call get_tile_data(restart,'fsc_in','zfull',soil_fsc_in_ptr)
     call get_tile_data(restart,'ssc_in','zfull',soil_ssc_in_ptr)
     if (soil_carbon_option==SOILC_CORPSE .and. field_exists(restart, 'deadmic_in')) then
        call get_tile_data(restart,'deadmic_in','zfull',soil_deadmic_in_ptr)
        call get_tile_data(restart,'fast_protected_in','zfull',soil_fast_protected_in_ptr)
        call get_tile_data(restart,'slow_protected_in','zfull', soil_slow_protected_in_ptr)
        call get_tile_data(restart,'deadmic_protected_in','zfull', soil_deadmic_protected_in_ptr)
        call get_tile_data(restart,'leaflitter_fsc_in', soil_leaflitter_fsc_in_ptr)
        call get_tile_data(restart,'leaflitter_ssc_in', soil_leaflitter_ssc_in_ptr)
        call get_tile_data(restart,'leaflitter_deadmic_in', soil_leaflitter_deadmic_in_ptr)
        call get_tile_data(restart,'finewoodlitter_fsc_in', soil_finewoodlitter_fsc_in_ptr)
        call get_tile_data(restart,'finewoodlitter_ssc_in', soil_finewoodlitter_ssc_in_ptr)
        call get_tile_data(restart,'finewoodlitter_deadmic_in', soil_finewoodlitter_deadmic_in_ptr)
        call get_tile_data(restart,'coarsewoodlitter_fsc_in', soil_coarsewoodlitter_fsc_in_ptr)
        call get_tile_data(restart,'coarsewoodlitter_ssc_in', soil_coarsewoodlitter_ssc_in_ptr)
        call get_tile_data(restart,'coarsewoodlitter_deadmic_in', soil_coarsewoodlitter_deadmic_in_ptr)
        call get_tile_data(restart,'fast_turnover_accumulated','zfull', soil_fast_turnover_accumulated_ptr)
        call get_tile_data(restart,'slow_turnover_accumulated','zfull', soil_slow_turnover_accumulated_ptr)
        call get_tile_data(restart,'deadmic_turnover_accumulated','zfull', soil_deadmic_turnover_accumulated_ptr)
        call get_tile_data(restart,'fast_protected_turnover_accumulated','zfull', soil_fast_protected_turnover_accumulated_ptr)
        call get_tile_data(restart,'slow_protected_turnover_accumulated','zfull', soil_slow_protected_turnover_accumulated_ptr)
        call get_tile_data(restart,'deadmic_protected_turnover_accumulated','zfull', soil_deadmic_protected_turnover_accumulated_ptr)
        call get_tile_data(restart,'leaflitter_fast_turnover_accumulated', soil_leaflitter_fast_turnover_accumulated_ptr)
        call get_tile_data(restart,'leaflitter_slow_turnover_accumulated', soil_leaflitter_slow_turnover_accumulated_ptr)
        call get_tile_data(restart,'leaflitter_deadmic_turnover_accumulated', soil_leaflitter_deadmic_turnover_accumulated_ptr)
        call get_tile_data(restart,'finewoodlitter_fast_turnover_accumulated', soil_finewoodlitter_fast_turnover_accumulated_ptr)
        call get_tile_data(restart,'finewoodlitter_slow_turnover_accumulated', soil_finewoodlitter_slow_turnover_accumulated_ptr)
        call get_tile_data(restart,'finewoodlitter_deadmic_turnover_accumulated', soil_finewoodlitter_deadmic_turnover_accumulated_ptr)
        call get_tile_data(restart,'coarsewoodlitter_fast_turnover_accumulated', soil_coarsewoodlitter_fast_turnover_accumulated_ptr)
        call get_tile_data(restart,'coarsewoodlitter_slow_turnover_accumulated', soil_coarsewoodlitter_slow_turnover_accumulated_ptr)
        call get_tile_data(restart,'coarsewoodlitter_deadmic_turnover_accumulated', soil_coarsewoodlitter_deadmic_turnover_accumulated_ptr)
     endif
  endif
  call free_land_restart(restart)

  ! ---- static diagnostic section
  call send_tile_data_r0d_fptr(id_tau_gw,       soil_tau_groundwater_ptr)
  call send_tile_data_r0d_fptr(id_slope_l,      soil_hillslope_length_ptr)
  call send_tile_data_r0d_fptr(id_slope_Z,      soil_hillslope_relief_ptr)
  call send_tile_data_r0d_fptr(id_zeta_bar,     soil_hillslope_zeta_bar_ptr)
  call send_tile_data_r0d_fptr(id_e_depth,      soil_soil_e_depth_ptr)
  call send_tile_data_r0d_fptr(id_zeta,         soil_zeta_ptr)
  call send_tile_data_r0d_fptr(id_tau,          soil_tau_ptr)
  call send_tile_data_r0d_fptr(id_vwc_wilt,     soil_vwc_wilt_ptr)
  call send_tile_data_r0d_fptr(id_vwc_fc,       soil_vwc_fc_ptr)
  call send_tile_data_r0d_fptr(id_vwc_sat,      soil_vwc_sat_ptr)
  call send_tile_data_r0d_fptr(id_K_sat,        soil_k_sat_ref_ptr)
  call send_tile_data_r0d_fptr(id_K_gw,         soil_k_sat_gw_ptr)
  call send_tile_data_r1d_fptr(id_w_fc,         soil_w_fc_ptr)
  call send_tile_data_r1d_fptr(id_alpha,        soil_alpha_ptr)
  call send_tile_data_r1d_fptr(id_refl_dry_dir, soil_refl_dry_dir_ptr)
  call send_tile_data_r1d_fptr(id_refl_dry_dif, soil_refl_dry_dif_ptr)
  call send_tile_data_r1d_fptr(id_refl_sat_dir, soil_refl_sat_dir_ptr)
  call send_tile_data_r1d_fptr(id_refl_sat_dif, soil_refl_sat_dif_ptr)
  call send_tile_data_r1d_fptr(id_f_iso_dry, soil_f_iso_dry_ptr)
  call send_tile_data_r1d_fptr(id_f_vol_dry, soil_f_vol_dry_ptr)
  call send_tile_data_r1d_fptr(id_f_geo_dry, soil_f_geo_dry_ptr)
  call send_tile_data_r1d_fptr(id_f_iso_sat, soil_f_iso_sat_ptr)
  call send_tile_data_r1d_fptr(id_f_vol_sat, soil_f_vol_sat_ptr)
  call send_tile_data_r1d_fptr(id_f_geo_sat, soil_f_geo_sat_ptr)
  call send_tile_data_i0d_fptr(id_type,         soil_tag_ptr)
end subroutine soil_init


! ============================================================================
subroutine soil_diag_init(id_ug,id_band,id_zfull)
  integer,intent(in)  :: id_ug    !<Unstructured axis id.
  integer,intent(in)  :: id_band  ! ID of spectral band axis
  integer,intent(out) :: id_zfull ! ID of vertical soil axis
!----------

  ! ---- local vars
  integer :: axes(2)
  integer :: id_zhalf
  integer :: l ! layer index
  character(80) :: str ! buffer for forming long names (e.g. see mrsos)

  ! define vertical axis and its edges
  id_zhalf = diag_axis_init ( &
       'zhalf_soil', zhalf(1:num_l+1), 'meters', 'z', 'half level',  -1, set_name='soil' )
  id_zfull = diag_axis_init ( &
       'zfull_soil', zfull(1:num_l),   'meters', 'z', 'full level',  -1, set_name='soil', &
       edges=id_zhalf )

  ! define array of axis indices
  axes = (/id_ug,id_zfull/)

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('soil')

  ! define diagnostic fields
  id_fast_soil_C = register_tiled_diag_field ( module_name, 'fast_soil_C', axes,  &
       lnd%time, 'fast soil carbon', 'kg C/m3', missing_value=-100.0 )
  id_slow_soil_C = register_tiled_diag_field ( module_name, 'slow_soil_C', axes,  &
       lnd%time, 'slow soil carbon', 'kg C/m3', missing_value=-100.0 )
  id_deadmic = register_tiled_diag_field ( module_name, 'dead_microbe_C', axes,  &
       lnd%time, 'dead microbe soil carbon per layer', 'kg C/m3', missing_value=-100.0 )
  id_protectedC = register_tiled_diag_field ( module_name, 'protected_soil_C', axes,  &
       lnd%time, 'protected soil carbon per layer', 'kg C/m3', missing_value=-100.0 )
  id_livemic = register_tiled_diag_field ( module_name, 'livemic',  &
       (/id_ug,id_zfull/), lnd%time, 'live microbe soil carbon', 'kg C/m3', &
       missing_value=-100.0 )
  id_fast_dissolved_C = register_tiled_diag_field ( module_name, 'fast_dissolved_C', axes,  &
       lnd%time, 'fast dissolved carbon per layer', 'kg C/m3', missing_value=-100.0 )
  id_slow_dissolved_C = register_tiled_diag_field ( module_name, 'slow_dissolved_C', axes,  &
       lnd%time, 'slow dissolved carbon per layer', 'kg C/m3', missing_value=-100.0 )
  id_deadmic_dissolved_C = register_tiled_diag_field ( module_name, 'dead_microbe_dissolved_C', axes,  &
       lnd%time, 'dead microbe dissolved carbon per layer', 'kg C/m3', missing_value=-100.0 )
  id_total_carbon_layered = register_tiled_diag_field ( module_name, 'total_soil_carbon_layered', axes,  &
       lnd%time, 'total soil carbon per layer', 'kg C/m3', missing_value=-100.0 )
  id_fast_DOC_div_loss = register_tiled_diag_field ( module_name, 'fast_DOC_div_loss', (/id_ug/),  &
       lnd%time, 'total fast DOC divergence loss', 'kg C/m2', missing_value=-100.0 )
  id_slow_DOC_div_loss = register_tiled_diag_field ( module_name, 'slow_DOC_div_loss', (/id_ug/),  &
       lnd%time, 'total slow DOC divergence loss', 'kg C/m2', missing_value=-100.0 )
  id_deadmic_DOC_div_loss = register_tiled_diag_field ( module_name, 'deadmic_DOC_div_loss', (/id_ug/),  &
       lnd%time, 'total dead microbe DOC divergence loss', 'kg C/m2', missing_value=-100.0 )
  id_total_DOC_div_loss = register_tiled_diag_field ( module_name, 'total_DOC_div', axes(1:1), &
       lnd%time, 'total rate of DOC divergence loss', 'kg C/m^2/s', missing_value=initval)
  id_rsoil = register_tiled_diag_field ( module_name, 'rsoil', (/id_ug/), &
       lnd%time, 'soil respiration', 'kg C/(m2 year)', missing_value=-100.0 )
  id_rsoil_fast = register_tiled_diag_field ( module_name, 'rsoil_fast',  axes, &
       lnd%time, 'fast soil carbon respiration', 'kg C/(m3 year)', missing_value=-100.0 )
  id_rsoil_slow = register_tiled_diag_field ( module_name, 'rsoil_slow',  axes, &
       lnd%time, 'slow soil carbon respiration', 'kg C/(m3 year)', missing_value=-100.0 )
  id_rsoil_deadmic = register_tiled_diag_field ( module_name, 'rsoil_deadmic',  axes, &
       lnd%time, 'dead microbe soil carbon respiration', 'kg C/(m3 year)', missing_value=-100.0 )
  id_rsoil_leaflitter_fast = register_tiled_diag_field ( module_name, 'rsoil_leaflitter_fast',  &
       (/id_ug/), lnd%time, 'surface leaf litter fast C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_leaflitter_slow = register_tiled_diag_field ( module_name, 'rsoil_leaflitter_slow',  &
       (/id_ug/), lnd%time, 'surface leaf litter slow C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_leaflitter_deadmic = register_tiled_diag_field ( module_name, 'rsoil_leaflitter_deadmic',  &
       (/id_ug/), lnd%time, 'surface leaf litter dead microbe C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_coarsewoodlitter_fast = register_tiled_diag_field ( module_name, 'rsoil_coarsewoodlitter_fast',  &
       (/id_ug/), lnd%time, 'surface coarse wood litter fast C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_coarsewoodlitter_slow = register_tiled_diag_field ( module_name, 'rsoil_coarsewoodlitter_slow',  &
       (/id_ug/), lnd%time, 'surface coarse wood litter slow C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_coarsewoodlitter_deadmic = register_tiled_diag_field ( module_name, 'rsoil_coarsewoodlitter_deadmic',  &
       (/id_ug/), lnd%time, 'surface coarse wood litter dead microbe C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_finewoodlitter_fast = register_tiled_diag_field ( module_name, 'rsoil_finewoodlitter_fast',  &
       (/id_ug/), lnd%time, 'surface fine wood litter fast C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_finewoodlitter_slow = register_tiled_diag_field ( module_name, 'rsoil_finewoodlitter_slow',  &
       (/id_ug/), lnd%time, 'surface fine wood litter slow C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_finewoodlitter_deadmic = register_tiled_diag_field ( module_name, 'rsoil_finewoodlitter_deadmic',  &
       (/id_ug/), lnd%time, 'surface fine wood litter dead microbe C respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_dissolved_fast = register_tiled_diag_field ( module_name, 'fast_dissolve_rate',  &
       axes, lnd%time, 'fast soil carbon dissolving rate', 'kg C/(m3 year)', &
       missing_value=-100.0 )
  id_dissolved_slow = register_tiled_diag_field ( module_name, 'slow_dissolve_rate',  &
       axes, lnd%time, 'slow soil carbon dissolving rate', 'kg C/(m3 year)', &
       missing_value=-100.0 )
  id_dissolved_deadmic = register_tiled_diag_field ( module_name, 'deadmic_dissolve_rate',  &
       axes, lnd%time, 'dead microbe soil carbon dissolving rate', 'kg C/(m3 year)', &
       missing_value=-100.0 )
  id_deposited_fast = register_tiled_diag_field ( module_name, 'fast_deposition_rate',  &
       axes, lnd%time, 'fast soil carbon deposition from DOC rate', 'kg C/(m3 year)', &
       missing_value=-100.0 )
  id_deposited_slow = register_tiled_diag_field ( module_name, 'slow_deposition_rate',  &
       axes, lnd%time, 'slow soil carbon deposition from DOC rate', 'kg C/(m3 year)', &
       missing_value=-100.0 )
  id_deposited_deadmic = register_tiled_diag_field ( module_name, 'deadmic_deposition_rate',  &
       axes, lnd%time, 'dead microbe soil carbon deposition from DOC rate', 'kg C/(m3 year)', &
       missing_value=-100.0 )
  id_leaflitter_dissolved_fast = register_tiled_diag_field ( module_name, 'leaflitter_fast_dissolve_rate',  &
       axes(1:1), lnd%time, 'fast leaf litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_leaflitter_dissolved_slow = register_tiled_diag_field ( module_name, 'leaflitter_slow_dissolve_rate',  &
       axes(1:1), lnd%time, 'slow leaf litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_leaflitter_dissolved_deadmic = register_tiled_diag_field ( module_name, 'leaflitter_deadmic_dissolve_rate',  &
       axes(1:1), lnd%time, 'dead microbe leaf litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_leaflitter_deposited_fast = register_tiled_diag_field ( module_name, 'leaflitter_fast_deposition_rate',  &
       axes(1:1), lnd%time, 'fast leaf litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_leaflitter_deposited_slow = register_tiled_diag_field ( module_name, 'leaflitter_slow_deposition_rate',  &
       axes(1:1), lnd%time, 'slow leaf litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_leaflitter_deposited_deadmic = register_tiled_diag_field ( module_name, 'leaflitter_deadmic_deposition_rate',  &
       axes(1:1), lnd%time, 'dead microbe leaf litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_finewoodlitter_dissolved_fast = register_tiled_diag_field ( module_name, 'finewoodlitter_fast_dissolve_rate',  &
       axes(1:1), lnd%time, 'fast fine wood litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_finewoodlitter_dissolved_slow = register_tiled_diag_field ( module_name, 'finewoodlitter_slow_dissolve_rate',  &
       axes(1:1), lnd%time, 'slow fine wood litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_finewoodlitter_dissolved_deadmic = register_tiled_diag_field ( module_name, 'finewoodlitter_deadmic_dissolve_rate',  &
       axes(1:1), lnd%time, 'dead microbe fine wood litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_finewoodlitter_deposited_fast = register_tiled_diag_field ( module_name, 'finewoodlitter_fast_deposition_rate',  &
       axes(1:1), lnd%time, 'fast fine wood litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_finewoodlitter_deposited_slow = register_tiled_diag_field ( module_name, 'finewoodlitter_slow_deposition_rate',  &
       axes(1:1), lnd%time, 'slow fine wood litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_finewoodlitter_deposited_deadmic = register_tiled_diag_field ( module_name, 'finewoodlitter_deadmic_deposition_rate',  &
       axes(1:1), lnd%time, 'dead microbe fine wood litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
    id_coarsewoodlitter_dissolved_fast = register_tiled_diag_field ( module_name, 'coarsewoodlitter_fast_dissolve_rate',  &
       axes(1:1), lnd%time, 'fast coarse wood litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_coarsewoodlitter_dissolved_slow = register_tiled_diag_field ( module_name, 'coarsewoodlitter_slow_dissolve_rate',  &
       axes(1:1), lnd%time, 'slow coarse wood litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_coarsewoodlitter_dissolved_deadmic = register_tiled_diag_field ( module_name, 'coarsewoodlitter_deadmic_dissolve_rate',  &
       axes(1:1), lnd%time, 'dead microbe coarse wood litter carbon dissolving rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_coarsewoodlitter_deposited_fast = register_tiled_diag_field ( module_name, 'coarsewoodlitter_fast_deposition_rate',  &
       axes(1:1), lnd%time, 'fast coarse wood litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_coarsewoodlitter_deposited_slow = register_tiled_diag_field ( module_name, 'coarsewoodlitter_slow_deposition_rate',  &
       axes(1:1), lnd%time, 'slow coarse wood litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_coarsewoodlitter_deposited_deadmic = register_tiled_diag_field ( module_name, 'coarsewoodlitter_deadmic_deposition_rate',  &
       axes(1:1), lnd%time, 'dead microbe coarse wood litter carbon deposition from DOC rate', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_Qmax = register_tiled_diag_field ( module_name, 'Qmax', axes(1:1),  &
       lnd%time, 'Maximum clay sorptive capacity', 'kg C/m3', missing_value=-100.0 )
  id_leaflitter_fast_C = register_tiled_diag_field ( module_name, 'fast_leaflitter_C', axes(1:1),  &
       lnd%time, 'fast leaf litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_leaflitter_slow_C = register_tiled_diag_field ( module_name, 'slow_leaflitter_C', axes(1:1),  &
       lnd%time, 'slow leaf litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_leaflitter_deadmic = register_tiled_diag_field ( module_name, 'leaflitter_dead_microbe_C', axes(1:1),  &
       lnd%time, 'dead microbe leaf litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_leaflitter_livemic = register_tiled_diag_field ( module_name, 'leaflitter_live_microbe_C', axes(1:1),  &
       lnd%time, 'live microbe leaf litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_coarsewoodlitter_fast_C = register_tiled_diag_field ( module_name, 'fast_coarsewoodlitter_C', axes(1:1),  &
       lnd%time, 'fast coarse wood litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_coarsewoodlitter_slow_C = register_tiled_diag_field ( module_name, 'slow_coarsewoodlitter_C', axes(1:1),  &
       lnd%time, 'slow coarse wood litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_coarsewoodlitter_deadmic = register_tiled_diag_field ( module_name, 'coarsewoodlitter_dead_microbe_C', axes(1:1),  &
       lnd%time, 'dead microbe coarse wood litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_coarsewoodlitter_livemic = register_tiled_diag_field ( module_name, 'coarsewoodlitter_live_microbe_C', axes(1:1),  &
       lnd%time, 'live microbe coarse wood litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_finewoodlitter_fast_C = register_tiled_diag_field ( module_name, 'fast_finewoodlitter_C', axes(1:1),  &
       lnd%time, 'fast fine wood litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_finewoodlitter_slow_C = register_tiled_diag_field ( module_name, 'slow_finewoodlitter_C', axes(1:1),  &
       lnd%time, 'slow fine wood litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_finewoodlitter_deadmic = register_tiled_diag_field ( module_name, 'finewoodlitter_dead_microbe_C', axes(1:1),  &
       lnd%time, 'dead microbe fine wood litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_finewoodlitter_livemic = register_tiled_diag_field ( module_name, 'finewoodlitter_live_microbe_C', axes(1:1),  &
       lnd%time, 'live microbe fine wood litter carbon', 'kg C/m2', missing_value=-100.0 )
  id_leaflitter_fast_dissolved_C = register_tiled_diag_field ( module_name, 'fast_leaflitter_dissolved_C', axes(1:1),  &
       lnd%time, 'fast leaf litter dissolved carbon', 'kg C/m2', missing_value=-100.0 )
  id_leaflitter_slow_dissolved_C = register_tiled_diag_field ( module_name, 'slow_leaflitter_dissolved_C', axes(1:1),  &
       lnd%time, 'slow leaf litter dissolved carbon', 'kg C/m2', missing_value=-100.0 )
  id_leaflitter_deadmic_dissolved_C = register_tiled_diag_field ( module_name, 'leaflitter_dead_microbe_dissolved_C', axes(1:1),  &
       lnd%time, 'dead microbe leaf litter dissolved carbon', 'kg C/m2', missing_value=-100.0 )
  id_finewoodlitter_fast_dissolved_C = register_tiled_diag_field ( module_name, 'fast_finewoodlitter_dissolved_C', axes(1:1),  &
       lnd%time, 'fast fine wood litter dissolved carbon', 'kg C/m2', missing_value=-100.0 )
  id_finewoodlitter_slow_dissolved_C = register_tiled_diag_field ( module_name, 'slow_finewoodlitter_dissolved_C', axes(1:1),  &
       lnd%time, 'slow fine wood litter dissolved carbon', 'kg C/m2', missing_value=-100.0 )
  id_finewoodlitter_deadmic_dissolved_C = register_tiled_diag_field ( module_name, 'finewoodlitter_dead_microbe_dissolved_C', axes(1:1),  &
       lnd%time, 'dead microbe fine wood litter dissolved carbon', 'kg C/m2', missing_value=-100.0 )

  id_coarsewoodlitter_fast_dissolved_C = register_tiled_diag_field ( module_name, 'fast_coarsewoodlitter_dissolved_C', axes(1:1),  &
       lnd%time, 'fast coarse woodlitter dissolved carbon', 'kg C/m2', missing_value=-100.0 )
  id_coarsewoodlitter_slow_dissolved_C = register_tiled_diag_field ( module_name, 'slow_coarsewoodlitter_dissolved_C', axes(1:1),  &
       lnd%time, 'slow coarse wood litter dissolved carbon', 'kg C/m2', missing_value=-100.0 )
  id_coarsewoodlitter_deadmic_dissolved_C = register_tiled_diag_field ( module_name, 'coarsewoodlitter_dead_microbe_dissolved_C', axes(1:1),  &
       lnd%time, 'dead microbe coarse wood litter dissolved carbon', 'kg C/m2', missing_value=-100.0 )
  id_livemic = register_tiled_diag_field ( module_name, 'live_microbe_C', axes,  &
       lnd%time, 'Total live microbe soil carbon', 'kg C/m3', missing_value=-100.0 )
  id_nsoilcohorts = register_tiled_diag_field ( module_name, 'n_soil_cohorts', axes,  &
       lnd%time, 'number of soil cohorts', missing_value=-100.0 )
  id_nleaflittercohorts = register_tiled_diag_field ( module_name, 'n_leaflitter_cohorts', axes(1:1),  &
       lnd%time, 'number of leaf litter cohorts', missing_value=-100.0 )
  id_nfinewoodlittercohorts = register_tiled_diag_field ( module_name, 'n_finewoodlitter_cohorts', axes(1:1),  &
       lnd%time, 'number of fine wood litter cohorts', missing_value=-100.0 )
  id_ncoarsewoodlittercohorts = register_tiled_diag_field ( module_name, 'n_coarsewoodlitter_cohorts', axes(1:1),  &
       lnd%time, 'number of coarse wood litter cohorts', missing_value=-100.0 )
  id_deadmic_total = register_tiled_diag_field ( module_name, 'deadmic_total', axes(1:1),  &
       lnd%time, 'total dead microbe soil carbon', 'kg C/m2', missing_value=-100.0 )
  id_livemic_total = register_tiled_diag_field ( module_name, 'livemic_total', axes(1:1),  &
       lnd%time, 'total live microbe soil carbon', 'kg C/m2', missing_value=-100.0 )
  id_protected_total = register_tiled_diag_field ( module_name, 'protected_total', axes(1:1),  &
       lnd%time, 'total protected soil carbon', 'kg C/m2', missing_value=-100.0 )
  id_dissolved_total = register_tiled_diag_field ( module_name, 'dissolved_total', axes(1:1),  &
       lnd%time, 'total dissolved soil carbon', 'kg C/m2', missing_value=-100.0 )
  id_total_soil_C = register_tiled_diag_field ( module_name, 'total_soil_C', axes(1:1),  &
       lnd%time, 'total soil carbon', 'kg C/m2', missing_value=-100.0 )
  id_fast_C_leaching = register_tiled_diag_field ( module_name, 'fast_C_leaching', axes, &
       lnd%time, 'net layer fast soil C leaching',  'kg/(m2 s)', missing_value=-100.0)
  id_slow_C_leaching = register_tiled_diag_field ( module_name, 'slow_C_leaching', axes, &
       lnd%time, 'net layer slow soil C leaching',  'kg/(m2 s)', missing_value=-100.0)
  id_deadmic_C_leaching = register_tiled_diag_field ( module_name, 'deadmic_C_leaching', axes, &
       lnd%time, 'net layer dead microbe soil C leaching',  'kg/(m2 s)', missing_value=-100.0)
  id_total_C_leaching = register_tiled_diag_field ( module_name, 'total_C_leaching', axes, &
       lnd%time, 'net layer total vertical soil C leaching', 'kg/(m2 s)', missing_value=initval)
  id_livemic_C_leaching = register_tiled_diag_field ( module_name, 'livemic_C_leaching', axes, &
       lnd%time, 'net layer live microbe C leaching',  'kg/(m2 s)', missing_value=-100.0)
  !id_protected_C_leaching = register_tiled_diag_field ( module_name, 'protected_C_leaching', axes, &
  !     lnd%time, 'net layer protected soil C leaching',  'kg/(m2 s)', missing_value=-100.0)
  id_leaflitter_fast_C_leaching = register_tiled_diag_field ( module_name, 'fast_leaflitter_C_leaching', axes(1:1), &
        lnd%time, 'Leaf litter fast C leaching','kg/(m2 s)', missing_value=-100.0)
  id_leaflitter_slow_C_leaching = register_tiled_diag_field ( module_name, 'slow_leaflitter_C_leaching', axes(1:1), &
        lnd%time, 'Leaf litter slow C leaching','kg/(m2 s)', missing_value=-100.0)
  id_leaflitter_deadmic_C_leaching = register_tiled_diag_field ( module_name, 'deadmic_leaflitter_C_leaching', axes(1:1), &
        lnd%time, 'Leaf litter dead microbe C leaching','kg/(m2 s)', missing_value=-100.0)
  id_coarsewoodlitter_fast_C_leaching = register_tiled_diag_field ( module_name, 'fast_coarsewoodlitter_C_leaching', axes(1:1), &
        lnd%time, 'Coarse wood litter fast C leaching','kg/(m2 s)', missing_value=-100.0)
  id_coarsewoodlitter_slow_C_leaching = register_tiled_diag_field ( module_name, 'slow_coarsewoodlitter_C_leaching', axes(1:1), &
        lnd%time, 'Coarse wood litter slow C leaching','kg/(m2 s)', missing_value=-100.0)
  id_coarsewoodlitter_deadmic_C_leaching = register_tiled_diag_field ( module_name, 'deadmic_coarsewoodlitter_C_leaching', axes(1:1), &
        lnd%time, 'Coarse wood litter dead microbe C leaching','kg/(m2 s)', missing_value=-100.0)
  id_finewoodlitter_fast_C_leaching = register_tiled_diag_field ( module_name, 'fast_finewoodlitter_C_leaching', axes(1:1), &
        lnd%time, 'Fine wood litter fast C leaching','kg/(m2 s)', missing_value=-100.0)
  id_finewoodlitter_slow_C_leaching = register_tiled_diag_field ( module_name, 'slow_finewoodlitter_C_leaching', axes(1:1), &
        lnd%time, 'Fine wood litter slow C leaching','kg/(m2 s)', missing_value=-100.0)
  id_finewoodlitter_deadmic_C_leaching = register_tiled_diag_field ( module_name, 'deadmic_finewoodlitter_C_leaching', axes(1:1), &
        lnd%time, 'Fine wood litter dead microbe C leaching','kg/(m2 s)', missing_value=-100.0)
  ! ZMS
  id_slomtot = register_tiled_diag_field ( module_name, 'total_lit_SOM_C', axes(1:1), &
       lnd%time, 'vertical sum of all litter and soil carbon pools', 'kg C/m^2', missing_value=-100.0)
  id_surf_DOC_loss = register_tiled_diag_field ( module_name, 'surf_DOC_loss', axes(1:1), &
       lnd%time, 'loss of top layer DOC to surface runoff due to efflux', 'kg C/m^2/s', &
       missing_value=initval)
  id_fsc = register_tiled_diag_field ( module_name, 'fsc', axes(1:1),  &
       lnd%time, 'total fast soil carbon', 'kg C/m2', missing_value=-100.0 )
  id_ssc = register_tiled_diag_field ( module_name, 'ssc', axes(1:1),  &
       lnd%time, 'total slow soil carbon', 'kg C/m2', missing_value=-100.0 )

!FIGURE OUT WHAT IS WRONG WITH THESE DIAGNOSITCS!
  id_lwc = register_tiled_diag_field ( module_name, 'soil_liq', axes,  &
       lnd%time, 'bulk density of liquid water', 'kg/m3', missing_value=-100.0 )
  id_swc  = register_tiled_diag_field ( module_name, 'soil_ice',  axes,  &
       lnd%time, 'bulk density of solid water', 'kg/m3',  missing_value=-100.0 )

  id_psi = register_tiled_diag_field ( module_name, 'soil_psi', axes,  &
       lnd%time, 'soil-water matric head', 'm', missing_value=-100.0 )
  id_temp  = register_tiled_diag_field ( module_name, 'soil_T',  axes,       &
       lnd%time, 'temperature',            'degK',  missing_value=-100.0 )
  id_ie  = register_tiled_diag_field ( module_name, 'soil_rie',  axes(1:1),  &
       lnd%time, 'inf exc runf',            'kg/(m2 s)',  missing_value=-100.0 )
  id_sn  = register_tiled_diag_field ( module_name, 'soil_rsn',  axes(1:1),  &
       lnd%time, 'satn runf',            'kg/(m2 s)',  missing_value=-100.0 )
  id_bf  = register_tiled_diag_field ( module_name, 'soil_rbf',  axes(1:1),  &
       lnd%time, 'baseflow',            'kg/(m2 s)',  missing_value=-100.0 )
  id_if  = register_tiled_diag_field ( module_name, 'soil_rif',  axes(1:1),  &
       lnd%time, 'interflow',            'kg/(m2 s)',  missing_value=-100.0 )
  id_al  = register_tiled_diag_field ( module_name, 'soil_ral',  axes(1:1),  &
       lnd%time, 'active layer flow',    'kg/(m2 s)',  missing_value=-100.0 )
  id_nu  = register_tiled_diag_field ( module_name, 'soil_rnu',  axes(1:1),  &
       lnd%time, 'numerical runoff',    'kg/(m2 s)',  missing_value=-100.0 )
  id_sc  = register_tiled_diag_field ( module_name, 'soil_rsc',  axes(1:1),  &
       lnd%time, 'lm2 groundwater runoff',    'kg/(m2 s)',  missing_value=-100.0 )
  id_hie  = register_tiled_diag_field ( module_name, 'soil_hie',  axes(1:1), &
       lnd%time, 'heat ie runf',            'W/m2',  missing_value=-100.0 )
  id_hsn  = register_tiled_diag_field ( module_name, 'soil_hsn',  axes(1:1), &
       lnd%time, 'heat sn runf',            'W/m2',  missing_value=-100.0 )
  id_hbf  = register_tiled_diag_field ( module_name, 'soil_hbf',  axes(1:1), &
       lnd%time, 'heat bf runf',            'W/m2',  missing_value=-100.0 )
  id_hif  = register_tiled_diag_field ( module_name, 'soil_hif',  axes(1:1), &
       lnd%time, 'heat if runf',            'W/m2',  missing_value=-100.0 )
  id_hal  = register_tiled_diag_field ( module_name, 'soil_hal',  axes(1:1), &
       lnd%time, 'heat al runf',            'W/m2',  missing_value=-100.0 )
  id_hnu  = register_tiled_diag_field ( module_name, 'soil_hnu',  axes(1:1), &
       lnd%time, 'heat nu runoff',          'W/m2',  missing_value=-100.0 )
  id_hsc  = register_tiled_diag_field ( module_name, 'soil_hsc',  axes(1:1), &
       lnd%time, 'heat sc runoff',          'W/m2',  missing_value=-100.0 )
  id_evap  = register_tiled_diag_field ( module_name, 'soil_evap',  axes(1:1), &
       lnd%time, 'soil evap',            'kg/(m2 s)',  missing_value=-100.0 )
  id_excess  = register_tiled_diag_field ( module_name, 'sfc_excess',  axes(1:1),  &
       lnd%time, 'sfc excess pushed down',    'kg/(m2 s)',  missing_value=-100.0 )

  id_uptk_n_iter  = register_tiled_diag_field ( module_name, 'uptake_n_iter',  axes(1:1), &
       lnd%time, 'number of iterations for soil uptake',  missing_value=-100.0 )
  id_uptk = register_tiled_diag_field ( module_name, 'soil_uptk', axes, &
       lnd%time, 'uptake of water by roots', 'kg/(m2 s)',  missing_value=-100.0 )
  id_psi_x0 = register_tiled_diag_field ( module_name, 'soil_psix0', axes(1:1), &
       lnd%time, 'xylem potential at z=0', 'm',  missing_value=-100.0 )
  id_sws_n_iter  = register_tiled_diag_field ( module_name, 'sws_n_iter',  axes(1:1), &
       lnd%time, 'number of iterations for soil water supply',  missing_value=-100.0 )
  id_psi_x0_sws = register_tiled_diag_field ( module_name, 'soil_psix0_sws', axes(1:1), &
       lnd%time, 'xylem potential at z=0 for max transpiration', 'm',  missing_value=-100.0 )
  id_deficit = register_tiled_diag_field ( module_name, 'soil_def', axes(1:1), &
       lnd%time, 'groundwater storage deficit', '-',  missing_value=-100.0 )
  id_deficit_2 = register_tiled_diag_field ( module_name, 'soil_def2', axes(1:1), &
       lnd%time, 'groundwater storage deficit2', '-',  missing_value=-100.0 )
  id_deficit_3 = register_tiled_diag_field ( module_name, 'soil_def3', axes(1:1), &
       lnd%time, 'groundwater storage deficit3', '-',  missing_value=-100.0 )
  id_deficit_4 = register_tiled_diag_field ( module_name, 'soil_def4', axes(1:1), &
       lnd%time, 'groundwater storage deficit4', '-',  missing_value=-100.0 )
  id_psi_bot = register_tiled_diag_field ( module_name, 'soil_psi_n', axes(1:1), &
       lnd%time, 'psi at bottom of soil column', 'm',  missing_value=-100.0 )
  id_sat_frac = register_tiled_diag_field ( module_name, 'soil_fsat', axes(1:1), &
       lnd%time, 'fraction of soil area saturated at surface', '-',  missing_value=-100.0 )
  id_stor_frac = register_tiled_diag_field ( module_name, 'soil_fgw', axes(1:1), &
       lnd%time, 'groundwater storage frac above base elev', '-',  missing_value=-100.0 )
  id_sat_depth = register_tiled_diag_field ( module_name, 'soil_wtdep', axes(1:1), &
       lnd%time, 'depth below sfc to saturated soil', 'm',  missing_value=-100.0 )
  id_sat_dept2 = register_tiled_diag_field ( module_name, 'soil_wtdp2', axes(1:1), &
       lnd%time, 'alt depth below sfc to saturated soil', 'm',  missing_value=-100.0 )
  id_z_cap = register_tiled_diag_field ( module_name, 'soil_zcap', axes(1:1), &
       lnd%time, 'depth below sfc to capillary fringe', 'm',  missing_value=-100.0 )

  id_div_bf = register_tiled_diag_field ( module_name, 'soil_dvbf', axes, &
       lnd%time, 'baseflow by layer', 'kg/(m2 s)',  missing_value=-100.0 )
  id_div_if = register_tiled_diag_field ( module_name, 'soil_dvif', axes, &
       lnd%time, 'interflow by layer', 'kg/(m2 s)',  missing_value=-100.0 )
  id_div_al = register_tiled_diag_field ( module_name, 'soil_dval', axes, &
       lnd%time, 'active-layer flow by layer', 'kg/(m2 s)',  missing_value=-100.0 )

  id_cf_1 = register_tiled_diag_field ( module_name, 'soil_cf_1', axes(1:1), &
       lnd%time, 'soil_cf_1', 'm',  missing_value=-100.0 )
  id_cf_3 = register_tiled_diag_field ( module_name, 'soil_cf_3', axes(1:1), &
       lnd%time, 'soil_cf_3', 'm',  missing_value=-100.0 )
  id_wt_1 = register_tiled_diag_field ( module_name, 'soil_wt_1', axes(1:1), &
       lnd%time, 'soil_wt_1', 'm',  missing_value=-100.0 )
  id_wt_2 = register_tiled_diag_field ( module_name, 'soil_wt_2', axes(1:1), &
       lnd%time, 'soil_wt_2', 'm',  missing_value=-100.0 )
  id_wt_2a = register_tiled_diag_field ( module_name, 'soil_wt_2a', axes(1:1), &
       lnd%time, 'Water Table Depth from Surface to Saturation', 'm',  missing_value=-100.0 )
  id_wt_2b = register_tiled_diag_field ( module_name, 'soil_wt_2b', axes(1:1), &
       lnd%time, 'Water Table Depth from Surface to Liquid Saturation', 'm',  missing_value=-100.0 )
  id_wt_3 = register_tiled_diag_field ( module_name, 'soil_wt_3', axes(1:1), &
       lnd%time, 'soil_wt_3', 'm',  missing_value=-100.0 )
  id_wt2_3 = register_tiled_diag_field ( module_name, 'soil_wt2_3', axes(1:1), &
       lnd%time, 'soil_wt2_3', 'm',  missing_value=-100.0 )
  id_wt_4 = register_tiled_diag_field ( module_name, 'soil_wt_4', axes(1:1), &
       lnd%time, 'Interpolated psi = 0 from Bottom Up', 'm',  missing_value=-100.0 )

  id_active_layer = register_tiled_diag_field ( module_name, 'soil_alt', axes(1:1), &
       lnd%time, 'active-layer thickness', 'm',  missing_value=-100.0 )
  id_heat_cap = register_tiled_diag_field ( module_name, 'soil_heat_cap',  &
       axes, lnd%time, 'heat capacity of dry soil','J/(m3 K)', missing_value=-100.0 )
  id_thermal_cond =  register_tiled_diag_field ( module_name, 'soil_tcon', &
       axes, lnd%time, 'soil thermal conductivity', 'W/(m K)',  missing_value=-100.0 )

  id_surface_water = register_tiled_diag_field (module_name, 'surface_water', &
       axes(1:1), lnd%time, 'surface water storage', 'm', missing_value=-100.0 )
  id_inun_frac = register_tiled_diag_field (module_name, 'inun_fraction', &
       axes(1:1), lnd%time, 'inundated area fraction', '-', missing_value=-100.0 )
  if (gw_option == GW_TILED) then
     id_wet_frac = register_tiled_diag_field (module_name, 'wet_fraction', &
       axes(1:1), lnd%time, 'diagnostic wetland fraction', '-', missing_value=-100.0 )
  end if
  if (gw_option == GW_TILED .and. simple_inundation) then
      id_rsn_frac = register_tiled_diag_field (module_name, 'surface_runoff_frac', &
         axes(1:1), lnd%time, 'effective fraction of throughfall converted to sat-excess surface runoff', '-', missing_value=-100.0 )
  end if
  id_flow = register_tiled_diag_field (module_name, 'flow', axes, &
       lnd%time, 'vertical soil water flow at interface above (+ downward)', 'mm/s', missing_value=initval )
  id_reflux = register_tiled_diag_field (module_name, 'reflux', axes(1:1), &
       lnd%time, 'upwards flow of soil water at surface; zero if flow into surface', 'mm/s', missing_value=-100.0 )
  id_macro_infilt = register_tiled_diag_field (module_name, 'macro_inf', axes(1:1), &
       lnd%time, 'infiltration (decrease to IE runoff) at soil surface due to vertical macroporosity', 'mm/s', missing_value=-100.0 )

  id_type = register_tiled_static_field ( module_name, 'soil_type',  &
       axes(1:1), 'soil type', missing_value=-1.0 )
  id_tau_gw = register_tiled_static_field ( module_name, 'tau_gw',  &
       axes(1:1), 'groundwater residence time', 's', missing_value=-100.0 )
  id_slope_l = register_tiled_static_field ( module_name, 'slope_l',  &
       axes(1:1), 'hillslope length', 'm', missing_value=-100.0 )
  id_slope_Z = register_tiled_static_field ( module_name, 'soil_rlief',  &
       axes(1:1), 'hillslope relief', 'm', missing_value=-100.0 )
  id_zeta_bar = register_tiled_static_field ( module_name, 'zeta_bar',  &
       axes(1:1), 'hillslope zeta bar', '-', missing_value=-100.0 )
  id_e_depth = register_tiled_static_field ( module_name, 'soil_depth',  &
       axes(1:1), 'soil hydraulic e-folding depth', 'm', missing_value=-100.0 )
  id_zeta = register_tiled_static_field ( module_name, 'soil_zeta',      &
       axes(1:1), 'soil depth/topo relief', '-',  missing_value=-100.0 )
  id_tau = register_tiled_static_field ( module_name, 'soil_tau',        &
       axes(1:1), 'gw transmissivity/soil transmissivity', '-',  missing_value=-100.0 )
  id_vwc_wilt = register_tiled_static_field ( module_name, 'soil_wilt',  &
       axes(1:1), 'wilting water content', '-', missing_value=-100.0 )
  id_vwc_fc = register_tiled_static_field ( module_name, 'soil_fc',  &
       axes(1:1), 'field capacity', '-', missing_value=-100.0 )
  id_vwc_sat = register_tiled_static_field ( module_name, 'soil_sat',  &
       axes(1:1), 'soil porosity', '-', missing_value=-100.0 )
  id_K_sat = register_tiled_static_field ( module_name, 'soil_Ksat',  &
       axes(1:1), 'soil sat. hydraulic conductivity', 'kg /(m2 s)', missing_value=-100.0 )
  id_K_gw  = register_tiled_static_field ( module_name, 'soil_K_gw',  &
       axes(1:1), 'deep hydraulic conductivity', 'kg /(m2 s)', missing_value=-100.0 )
!----------
  id_w_fc = register_tiled_static_field ( module_name, 'w_fc',  &
       axes, 'soil field capacity', missing_value=-1.0 )
  id_alpha = register_tiled_static_field ( module_name, 'soil_alpha',  &
       axes, 'soil microscopic length scale', missing_value=-1.0 )
  id_refl_dry_dir = register_tiled_static_field ( module_name, 'refl_dry_dir',  &
       (/id_ug,id_band/), 'reflectance of dry soil for direct light', &
       missing_value=-1.0 )
  id_refl_dry_dif = register_tiled_static_field ( module_name, 'refl_dry_dif',  &
       (/id_ug,id_band/), 'reflectance of dry soil for diffuse light', &
       missing_value=-1.0 )
  id_refl_sat_dir = register_tiled_static_field ( module_name, 'refl_sat_dir',  &
       (/id_ug,id_band/), 'reflectance of saturated soil for direct light', &
       missing_value=-1.0 )
  id_refl_sat_dif = register_tiled_static_field ( module_name, 'refl_sat_dif',  &
       (/id_ug,id_band/), 'reflectance of saturated soil for diffuse light', &
       missing_value=-1.0 )
  id_f_iso_dry = register_tiled_static_field ( module_name, 'f_iso_dry',  &
       (/id_ug,id_band/), 'isotropic brdf weight, dry soil', &
       missing_value=-1.0 )
  id_f_vol_dry = register_tiled_static_field ( module_name, 'f_vol_dry',  &
       (/id_ug,id_band/), 'volumetric brdf weight, dry soil', &
       missing_value=-1.0 )
  id_f_geo_dry = register_tiled_static_field ( module_name, 'f_geo_dry',  &
       (/id_ug,id_band/), 'geometric brdf weight, dry soil', &
       missing_value=-1.0 )
  id_f_iso_sat = register_tiled_static_field ( module_name, 'f_iso_sat',  &
       (/id_ug,id_band/), 'isotropic brdf weight, saturated soil', &
       missing_value=-1.0 )
  id_f_vol_sat = register_tiled_static_field ( module_name, 'f_vol_sat',  &
       (/id_ug,id_band/), 'volumetric brdf weight, saturated soil', &
       missing_value=-1.0 )
  id_f_geo_sat = register_tiled_static_field ( module_name, 'f_geo_sat',  &
       (/id_ug,id_band/), 'geometric brdf weight, saturated soil', &
       missing_value=-1.0 )

  id_asoil = register_tiled_diag_field ( module_name, 'asoil', &
       (/id_ug/), lnd%time, 'aerobic activity modifier', &
       missing_value=-100.0 )

  ! the following fields are for compatibility with older diag tables only
  call add_tiled_static_field_alias ( id_slope_Z, module_name, 'slope_Z',  &
       axes(1:1), 'hillslope relief (obsolete, use "soil_rlief" instead)',&
       'm', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_e_depth, module_name, 'e_depth',  &
       axes(1:1), 'soil e-folding depth (obsolete, use "soil_depth" instead)', &
       'm', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_vwc_wilt, module_name, 'vwc_wilt',  &
       axes(1:1), 'wilting water content (obsolete, use "soil_wilt" instead)', &
       '-', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_vwc_fc, module_name, 'vwc_fc',  &
       axes(1:1), 'field capacity (obsolete, use "soil_fc" instead)', &
       '-', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_vwc_sat, module_name, 'vwc_sat',  &
       axes(1:1), 'soil porosity (obsolete, use "soil_sat")', &
       '-', missing_value=-100.0 )
  call add_tiled_static_field_alias ( id_K_sat, module_name, 'K_sat',  &
       axes(1:1), 'soil sat. hydraulic conductivity (obsolte, use "soil_Ksat" instead)', &
       'kg /(m2 s)', missing_value=-100.0 )

#ifdef ZMSDEBUG_TRIDIAGTEST
  ! For testing tridiagonal solution for advection
  id_st_diff = register_tiled_diag_field ( module_name, 'soil_T_diff', axes,  &
      lnd%time, 'soil Temperature difference after advection with tridiagonal solution', &
      'K', missing_value=-100.0 )
#endif

  ! CMOR variables
  ! tsl (soil temperature) should be reported as missing for the non-soil grid cells;
  ! also averaging over non-soil tiles does not make sense and therefore is not done.
  ! The only difference with soil_T is metadata (units and standard_name).
  call add_tiled_diag_field_alias ( id_temp, cmor_name, 'tsl', axes(1:2),  &
       lnd%time, 'Temperature of Soil', 'K', missing_value=-100.0, &
       standard_name='soil_temperature', fill_missing=.FALSE.)
  ! set up weights for mrsos averaging
  allocate(mrsos_weight(num_l), mrs1m_weight(num_l))
  do l = 1,num_l
     mrsos_weight(l) = min(1.0,max(0.0,(cmor_mrsos_depth-zhalf(l))/dz(l)))
     mrs1m_weight(l) = min(1.0,max(0.0,(1.0             -zhalf(l))/dz(l)))
  enddo
  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('land')
  id_mrlsl = register_tiled_diag_field ( cmor_name, 'mrlsl', axes,  &
       lnd%time, 'Water Content of Soil Layer', 'kg m-2', missing_value=-100.0, &
       standard_name='moisture_content_of_soil_layer', fill_missing=.TRUE.)
  id_mrsfl = register_tiled_diag_field ( cmor_name, 'mrsfl', axes,  &
       lnd%time, 'Frozen Water Content of Soil Layer', 'kg m-2', missing_value=-100.0, &
       standard_name='frozen_moisture_content_of_soil_layer', fill_missing=.TRUE.)
  id_mrsll = register_tiled_diag_field ( cmor_name, 'mrsll', axes,  &
       lnd%time, 'Liquid Water Content of Soil Layer', 'kg m-2', missing_value=-100.0, &
       standard_name='liquid_moisture_content_of_soil_layer', fill_missing=.TRUE.)
  id_mrsol = register_tiled_diag_field ( cmor_name, 'mrsol', axes,  &
       lnd%time, 'Total Water Content of Soil Layer', 'kg m-2', missing_value=-100.0, &
       standard_name='total_moisture_content_of_soil_layer', fill_missing=.TRUE.)
  id_mrso  = register_tiled_diag_field ( cmor_name, 'mrso', axes(1:1),  &
       lnd%time, 'Total Soil Moisture Content', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_moisture_content', fill_missing=.TRUE.)
  call add_tiled_diag_field_alias ( id_mrso, cmor_name, 'mrsoLut', axes(1:1),  &
       lnd%time, 'Total Soil Moisture Content', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_moisture_content_lut', fill_missing=.FALSE.)
  id_mrsos  = register_tiled_diag_field ( cmor_name, 'mrsos', axes(1:1),  &
       lnd%time, 'Moisture in Upper Portion of Soil Column', &
       'kg m-2', missing_value=-100.0, standard_name='moisture_content_of_soil_layer', &
       fill_missing=.TRUE.)
  call  add_tiled_diag_field_alias ( id_mrsos, cmor_name, 'mrsosLut', axes(1:1),  &
       lnd%time, 'Moisture in Upper Portion of Soil Column of Land Use Tile', &
       'kg m-2', missing_value=-100.0, standard_name='moisture_content_of_soil_layer', &
       fill_missing=.FALSE.)
  id_mrfso = register_tiled_diag_field ( cmor_name, 'mrfso', axes(1:1),  &
       lnd%time, 'Soil Frozen Water Content', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_frozen_water_content', fill_missing=.TRUE.)
  id_mrlso = register_tiled_diag_field ( cmor_name, 'mrlso', axes(1:1),  &
       lnd%time, 'Soil Liquid Water Content', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_frozen_water_content', fill_missing=.TRUE.)
  id_mrros = register_tiled_diag_field ( cmor_name, 'mrros',  axes(1:1),  &
       lnd%time, 'Surface Runoff', 'kg m-2 s-1',  missing_value=-100.0, &
       standard_name='surface_runoff_flux', fill_missing=.TRUE.)
  id_mrro = register_tiled_diag_field ( cmor_name, 'mrro',  axes(1:1),  &
       lnd%time, 'Total Runoff', 'kg m-2 s-1',  missing_value=-100.0, &
       standard_name='runoff_flux', fill_missing=.TRUE.)
  call add_tiled_diag_field_alias ( id_mrro, cmor_name, 'mrroLut',  axes(1:1),  &
       lnd%time, 'Total Runoff From Land Use Tile', 'kg m-2 s-1',  missing_value=-100.0, &
       standard_name='runoff_flux_lut', fill_missing=.FALSE.)

  id_csoil = register_tiled_diag_field ( cmor_name, 'cSoil', axes(1:1),  &
       lnd%time, 'Carbon in Soil Pool', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_carbon_content', fill_missing=.TRUE.)
  call add_tiled_diag_field_alias ( id_csoil, cmor_name, 'cSoilLut', axes(1:1),  &
       lnd%time, 'Carbon  In Soil Pool On Land Use Tiles', 'kg m-2', missing_value=-100.0, &
       standard_name='soil_carbon_content_lut', fill_missing=.FALSE.)

  id_csoilfast = register_tiled_diag_field ( cmor_name, 'cSoilFast', axes(1:1),  &
       lnd%time, 'Carbon in Fast Soil Pool', 'kg C m-2', missing_value=-100.0, &
       standard_name='carbon_in_fast_soil_pool', fill_missing=.TRUE.)
  id_csoilmedium = register_tiled_diag_field ( cmor_name, 'cSoilMedium', axes(1:1),  &
       lnd%time, 'Carbon in Medium Soil Pool', 'kg C m-2', missing_value=-100.0, &
       standard_name='carbon_in_medium_soil_pool', fill_missing=.TRUE.)
  id_csoilslow = register_tiled_diag_field ( cmor_name, 'cSoilSlow', axes(1:1),  &
       lnd%time, 'Carbon in Slow Soil Pool', 'kg C m-2', missing_value=-100.0, &
       standard_name='carbon_in_fast_soil_pool', fill_missing=.TRUE.)
  id_rh = register_tiled_diag_field ( cmor_name, 'rh', (/id_ug/), &
       lnd%time, 'Heterotrophic Respiration', 'kg m-2 s-1', missing_value=-1.0, &
       standard_name='heterotrophic_respiration_carbon_flux', fill_missing=.TRUE.)
  call add_tiled_diag_field_alias ( id_rh, cmor_name, 'rhLut', axes(1:1),  &
       lnd%time, 'Soil Heterotrophic Respiration On Land Use Tile', 'kg m-2 s-1', &
       standard_name='heterotrophic_respiration_carbon_flux', fill_missing=.FALSE., &
       missing_value=-100.0)
  id_mrs1mLut = register_tiled_diag_field ( cmor_name, 'mrs1mLut', axes(1:1), &
       lnd%time, 'Moisture in Top 1 Meter of Land Use Tile Soil Column', 'kg m-2', &
       missing_value=-100.0, standard_name='moisture_content_of_soil_layer', &
       fill_missing=.FALSE.)

end subroutine soil_diag_init


! ============================================================================
subroutine soil_end ()

  module_is_initialized =.FALSE.

end subroutine soil_end


! ============================================================================
subroutine save_soil_restart (tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars ----------------------------------------------------------
  character(267) :: filename
  type(land_restart_type) :: restart ! restart file i/o object
  type(land_tile_enum_type)     :: ce   ! tile list enumerator
  type(land_tile_type), pointer :: tile ! pointer to current tile
  integer :: i

  call error_mesg('soil_end','writing NetCDF restart',NOTE)
! must set domain so that io_domain is available
! Note that filename is updated for tile & rank numbers during file creation
  filename = trim(timestamp)//'soil.res.nc'
  call init_land_restart(restart, filename, soil_tile_exists, tile_dim_length)
  call add_restart_axis(restart,'zfull',zfull(1:num_l),'Z','m','full level',sense=-1)
  if (soil_carbon_option==SOILC_CORPSE) then
     call add_restart_axis(restart,'soilCCohort',(/(float(i),i=1,soilMaxCohorts)/),'CC')
  endif

  ! write out fields
  call add_tile_data(restart,'temp'         , 'zfull', soil_T_ptr,  'soil temperature','degrees_K')
  call add_tile_data(restart,'wl'           , 'zfull', soil_wl_ptr, 'liquid water content','kg/m2')
  call add_tile_data(restart,'ws'           , 'zfull', soil_ws_ptr, 'solid water content','kg/m2')
  call add_tile_data(restart,'groundwater'  , 'zfull', soil_groundwater_ptr, units='kg/m2' )
  call add_tile_data(restart,'groundwater_T', 'zfull', soil_groundwater_T_ptr, 'groundwater temperature','degrees_K' )
  call add_tile_data(restart,'uptake_T', soil_uptake_T_ptr, 'temperature of transpiring water', 'degrees_K')
  select case(soil_carbon_option)
  case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
     call add_tile_data(restart,'fsc', 'zfull', soil_fast_soil_C_ptr ,'fast soil carbon', 'kg C/m2')
     call add_tile_data(restart,'ssc', 'zfull', soil_slow_soil_C_ptr ,'slow soil carbon', 'kg C/m2')
  case (SOILC_CORPSE)
     ce = first_elmt(land_tile_map)
     do while(loop_over_tiles(ce,tile))
        if (.not.associated(tile%soil)) cycle
        call adjust_pool_ncohorts(tile%soil%leafLitter)
        call adjust_pool_ncohorts(tile%soil%fineWoodLitter)
        call adjust_pool_ncohorts(tile%soil%coarseWoodLitter)
        do i = 1,num_l
           call adjust_pool_ncohorts(tile%soil%soil_C(i))
        enddo
     end do

     call add_tile_data(restart,'fast_soil_C','zfull','soilCCohort', &
                        sc_soil_C_ptr,C_CEL,'Fast soil carbon','kg/m2')
     call add_tile_data(restart,'slow_soil_C', 'zfull','soilCCohort', &
                        sc_soil_C_ptr,C_LIG,'Slow soil carbon','kg/m2')
     call add_tile_data(restart,'deadMic', 'zfull','soilCCohort', &
                        sc_soil_C_ptr,C_MIC,'Dead microbe carbon','kg/m2')
     call add_tile_data(restart,'fastProtectedC', 'zfull','soilCCohort', &
                        sc_protected_C_ptr,C_CEL,'Protected fast carbon','kg/m2')
     call add_tile_data(restart,'slowProtectedC', 'zfull','soilCCohort', &
                        sc_protected_C_ptr,C_LIG,'Protected slow carbon','kg/m2')
     call add_tile_data(restart,'deadMicrobeProtectedC', 'zfull','soilCCohort', &
                        sc_protected_C_ptr,C_MIC,'Protected dead microbe carbon','kg/m2')

     call add_tile_data(restart,'liveMic', 'zfull','soilCCohort', &
                        soilc_livingMicrobeC_ptr,'Living microbial carbon','kg/m2')
     call add_tile_data(restart,'CO2', 'zfull','soilCCohort', &
                        soilc_CO2_ptr,'Cohort CO2 generated','kg/m2')
     call add_tile_data(restart,'Rtot', 'zfull','soilCCohort', &
                        soilc_Rtot_ptr,'Total degradation','kg/m2')
     call add_tile_data(restart,'originalCohortC', 'zfull','soilCCohort', &
                        soilc_originalLitterC_ptr,'Cohort original carbon','g/m2')

     call add_tile_data(restart,'soil_DOC_fast', 'zfull', soil_fast_DOC_ptr ,'Dissolved fast carbon','kg/m2')
     call add_tile_data(restart,'soil_DOC_slow', 'zfull', soil_slow_DOC_ptr ,'Dissolved slow carbon','kg/m2')
     call add_tile_data(restart,'soil_DOC_deadmic', 'zfull', &
                        soil_deadmicrobe_DOC_ptr ,'Dissolved dead microbe carbon','kg/m2')

     call add_tile_data(restart,'fast_DOC_leached',     soil_fast_DOC_leached_ptr, &
                        'Cumulative fast DOC leached out of the column', 'kg/m2')
     call add_tile_data(restart,'slow_DOC_leached',     soil_slow_DOC_leached_ptr, &
                        'Cumulative slow DOC leached out of the column', 'kg/m2')
     call add_tile_data(restart,'deadmic_DOC_leached',     soil_deadmic_DOC_leached_ptr, &
                        'Cumulative dead microbe DOC leached out of the column', 'kg/m2')

     call add_tile_data(restart,'leaf_litter_fast_C', 'soilCCohort', &
                        soilc_leafLitter_litterC_ptr,C_CEL,'Leaf litter fast C','kg/m2')
     call add_tile_data(restart,'leaf_litter_slow_C', 'soilCCohort', &
                        soilc_leafLitter_litterC_ptr,C_LIG,'Leaf litter slow C','kg/m2')
     call add_tile_data(restart,'leaf_litter_deadMic_C', 'soilCCohort', &
                        soilc_leafLitter_litterC_ptr,C_MIC,'Leaf litter dead microbe C','kg/m2')
     call add_tile_data(restart,'leaf_litter_liveMic_C','soilCCohort', &
                        soilc_leafLitter_livingMicrobeC_ptr,'Leaf litter live microbe C','kg/m2')
     call add_tile_data(restart,'leaf_litter_CO2','soilCCohort', &
                        soilc_leafLitter_CO2_ptr,'Leaf litter CO2 generated','kg/m2')
     call add_tile_data(restart,'leaf_litter_Rtot','soilCCohort', &
                        soilc_leafLitter_Rtot_ptr,'Leaf litter total degradation','kg/m2')
     call add_tile_data(restart,'leaf_litter_originalCohortC','soilCCohort', &
                        soilc_leafLitter_originalLitterC_ptr,'Leaf litter cohort original carbon','kg/m2')
     call add_tile_data(restart,'leaf_litter_fastProtectedC', 'soilCCohort', &
                        soilc_leafLitter_protectedC_ptr,C_CEL,'Leaf litter fast protected C','kg/m2')
     call add_tile_data(restart,'leaf_litter_slowProtectedC', 'soilCCohort', &
                        soilc_leafLitter_protectedC_ptr,C_LIG,'Leaf litter slow protected C','kg/m2')
     call add_tile_data(restart,'leaf_litter_deadMicrobeProtectedC','soilCCohort', &
                        soilc_leafLitter_protectedC_ptr,C_MIC,'Leaf litter dead microbe protected C','kg/m2')

     call add_tile_data(restart,'leaf_litter_DOC_fast',soilc_leafLitter_DOC_ptr,C_CEL,'Dissolved leaf litter fast carbon','kg/m2')
     call add_tile_data(restart,'leaf_litter_DOC_slow',soilc_leafLitter_DOC_ptr,C_LIG,'Dissolved leaf litter slow carbon','kg/m2')
     call add_tile_data(restart,'leaf_litter_DOC_deadmic',soilc_leafLitter_DOC_ptr,C_MIC,'Dissolved leaf litter dead microbe carbon','kg/m2')

     call add_tile_data(restart,'fineWood_litter_fast_C', 'soilCCohort', &
                        soilc_fineWoodLitter_litterC_ptr,C_CEL,'Fine wood litter fast C','kg/m2')
     call add_tile_data(restart,'fineWood_litter_slow_C', 'soilCCohort', &
                        soilc_fineWoodLitter_litterC_ptr,C_LIG,'Fine wood litter slow C','kg/m2')
     call add_tile_data(restart,'fineWood_litter_deadMic_C', 'soilCCohort', &
                        soilc_fineWoodLitter_litterC_ptr,C_MIC,'Fine wood litter dead microbe C','kg/m2')
     call add_tile_data(restart,'fineWood_litter_liveMic_C', 'soilCCohort', &
                        soilc_fineWoodLitter_livingMicrobeC_ptr,'Fine wood litter live microbe C','kg/m2')
     call add_tile_data(restart,'fineWood_litter_CO2','soilCCohort', &
                        soilc_fineWoodLitter_CO2_ptr,'Fine wood litter CO2 generated','kg/m2')
     call add_tile_data(restart,'fineWood_litter_Rtot','soilCCohort', &
                        soilc_fineWoodLitter_Rtot_ptr,'Fine wood litter total degradation','kg/m2')
     call add_tile_data(restart,'fineWood_litter_originalCohortC','soilCCohort', &
                        soilc_fineWoodLitter_originalLitterC_ptr,'Fine wood litter cohort original carbon','kg/m2')
     call add_tile_data(restart,'fineWood_litter_fastProtectedC', 'soilCCohort', &
                        soilc_fineWoodLitter_protectedC_ptr,C_CEL,'Fine wood litter fast protected C','kg/m2')
     call add_tile_data(restart,'fineWood_litter_slowProtectedC', 'soilCCohort', &
                        soilc_fineWoodLitter_protectedC_ptr,C_LIG,'Fine wood litter slow protected C','kg/m2')
     call add_tile_data(restart,'fineWood_litter_deadMicrobeProtectedC', 'soilCCohort', &
                        soilc_fineWoodLitter_protectedC_ptr,C_MIC,'Fine wood litter dead microbe protected C','kg/m2')

     call add_tile_data(restart,'fineWood_litter_DOC_fast',soilc_fineWoodLitter_DOC_ptr,C_CEL,&
                        'Dissolved fine wood litter fast carbon','kg/m2')
     call add_tile_data(restart,'fineWood_litter_DOC_slow',soilc_fineWoodLitter_DOC_ptr,C_LIG,&
                        'Dissolved fine wood litter slow carbon','kg/m2')
     call add_tile_data(restart,'fineWood_litter_DOC_deadmic',soilc_fineWoodLitter_DOC_ptr,C_MIC,&
                        'Dissolved fine wood litter dead microbe carbon','kg/m2')

     call add_tile_data(restart,'coarseWood_litter_fast_C', 'soilCCohort', &
                        soilc_coarseWoodLitter_litterC_ptr,C_CEL,'Coarse wood litter fast C','kg/m2')
     call add_tile_data(restart,'coarseWood_litter_slow_C', 'soilCCohort', &
                        soilc_coarseWoodLitter_litterC_ptr,C_LIG,'Coarse wood litter slow C','kg/m2')
     call add_tile_data(restart,'coarseWood_litter_deadMic_C', 'soilCCohort', &
                        soilc_coarseWoodLitter_litterC_ptr,C_MIC,'Coarse wood litter dead microbe C','kg/m2')
     call add_tile_data(restart,'coarseWood_litter_liveMic_C', 'soilCCohort', &
                        soilc_coarseWoodLitter_livingMicrobeC_ptr,'Coarse wood litter live microbe C','kg/m2')
     call add_tile_data(restart,'coarseWood_litter_CO2','soilCCohort', &
                        soilc_coarseWoodLitter_CO2_ptr,'Coarse wood litter CO2 generated','kg/m2')
     call add_tile_data(restart,'coarseWood_litter_Rtot','soilCCohort', &
                        soilc_coarseWoodLitter_Rtot_ptr,'Coarse wood litter total degradation','kg/m2')
     call add_tile_data(restart,'coarseWood_litter_originalCohortC','soilCCohort', &
                        soilc_coarseWoodLitter_originalLitterC_ptr,'Coarse wood litter cohort original carbon','kg/m2')
     call add_tile_data(restart,'coarseWood_litter_fastProtectedC', 'soilCCohort', &
                        soilc_coarseWoodLitter_protectedC_ptr,C_CEL,'Coarse wood litter fast protected C','kg/m2')
     call add_tile_data(restart,'coarseWood_litter_slowProtectedC', 'soilCCohort', &
                        soilc_coarseWoodLitter_protectedC_ptr,C_LIG,'Coarse wood litter slow protected C','kg/m2')
     call add_tile_data(restart,'coarseWood_litter_deadMicrobeProtectedC', 'soilCCohort', &
                        soilc_coarseWoodLitter_protectedC_ptr,C_MIC,'Coarse wood litter dead microbe protected C','kg/m2')

     call add_tile_data(restart,'coarseWood_litter_DOC_fast',soilc_coarseWoodLitter_DOC_ptr,C_CEL, &
                        'Dissolved coarse wood litter fast carbon','kg/m2')
     call add_tile_data(restart,'coarseWood_litter_DOC_slow',soilc_coarseWoodLitter_DOC_ptr,C_LIG, &
                        'Dissolved coarse wood litter slow carbon','kg/m2')
     call add_tile_data(restart,'coarseWood_litter_DOC_deadmic',soilc_coarseWoodLitter_DOC_ptr,C_MIC, &
                        'Dissolved coarse wood litter dead microbe carbon','kg/m2')

     call add_int_tile_data(restart,'is_peat', 'zfull', soil_is_peat_ptr ,'Is layer peat?','Boolean')
  case default
     call error_mesg('save_soil_restart','soil_carbon_option is invalid. This should never happen. Contact developer', FATAL)
  end select
  call save_land_restart(restart)
  call free_land_restart(restart)

  if (write_soil_carbon_restart) then
     filename = trim(timestamp)//'soil_carbon.res.nc'
     call init_land_restart(restart, filename, soil_tile_exists, tile_dim_length)
     call add_restart_axis(restart,'zfull',zfull(1:num_l),'Z','m','full level',sense=-1)

     call add_tile_data(restart,'asoil_in', 'zfull', soil_asoil_in_ptr ,'aerobic activity modifier', 'unitless')
     call add_tile_data(restart,'fsc_in', 'zfull', soil_fsc_in_ptr ,'fast soil carbon input', 'kg C/m2')
     call add_tile_data(restart,'ssc_in', 'zfull', soil_ssc_in_ptr ,'slow soil carbon input', 'kg C/m2')
     if (soil_carbon_option == SOILC_CORPSE) then
        call add_tile_data(restart,'deadmic_in', 'zfull', &
                           soil_deadmic_in_ptr ,'dead microbe soil carbon input', 'kg C/m2')
        call add_tile_data(restart,'fast_protected_in', 'zfull', &
                           soil_fast_protected_in_ptr ,'protected fast soil carbon input', 'kg C/m2')
        call add_tile_data(restart,'slow_protected_in', 'zfull', &
                           soil_slow_protected_in_ptr ,'protected slow soil carbon input', 'kg C/m2')
        call add_tile_data(restart,'deadmic_protected_in', 'zfull', &
                           soil_deadmic_protected_in_ptr ,'protected dead microbe soil carbon input', 'kg C/m2')
        call add_tile_data(restart,'fast_turnover_accumulated', 'zfull', &
                           soil_fast_turnover_accumulated_ptr ,'fast soil carbon turnover', 'year-1')
        call add_tile_data(restart,'slow_turnover_accumulated', 'zfull', &
                           soil_slow_turnover_accumulated_ptr ,'slow soil carbon turnover', 'year-1')
        call add_tile_data(restart,'deadmic_turnover_accumulated', 'zfull', &
                           soil_deadmic_turnover_accumulated_ptr ,'dead microbe soil carbon turnover', 'year-1')
        call add_tile_data(restart,'fast_protected_turnover_accumulated', 'zfull', &
                           soil_fast_protected_turnover_accumulated_ptr ,'fast protected soil carbon turnover', 'year-1')
        call add_tile_data(restart,'slow_protected_turnover_accumulated', 'zfull', &
                           soil_slow_protected_turnover_accumulated_ptr ,'slow protected soil carbon turnover', 'year-1')
        call add_tile_data(restart,'deadmic_protected_turnover_accumulated', 'zfull', &
             soil_deadmic_protected_turnover_accumulated_ptr ,'dead microbe protected soil carbon turnover', 'year-1')

        call add_tile_data(restart,'leaflitter_fast_turnover_accumulated',&
             soil_leaflitter_fast_turnover_accumulated_ptr,'fast leaf litter carbon turnover', 'year-1')
        call add_tile_data(restart,'leaflitter_slow_turnover_accumulated', &
             soil_leaflitter_slow_turnover_accumulated_ptr,'slow leaf litter carbon turnover', 'year-1')
        call add_tile_data(restart,'leaflitter_deadmic_turnover_accumulated', &
             soil_leaflitter_deadmic_turnover_accumulated_ptr,'dead microbe leaf litter carbon turnover', 'year-1')
        call add_tile_data(restart,'leaflitter_fsc_in',soil_leaflitter_fsc_in_ptr,'fast leaf litter carbon input', 'kg C/m2')
        call add_tile_data(restart,'leaflitter_ssc_in',soil_leaflitter_ssc_in_ptr,'slow leaf litter carbon input', 'kg C/m2')
        call add_tile_data(restart,'leaflitter_deadmic_in',&
             soil_leaflitter_deadmic_in_ptr,'dead microbe leaf litter carbon input', 'kg C/m2')

        call add_tile_data(restart,'finewoodlitter_fast_turnover_accumulated',&
             soil_finewoodlitter_fast_turnover_accumulated_ptr,'fast fine wood litter carbon turnover', 'year-1')
        call add_tile_data(restart,'finewoodlitter_slow_turnover_accumulated',&
             soil_finewoodlitter_slow_turnover_accumulated_ptr,'slow fine wood litter carbon turnover', 'year-1')
        call add_tile_data(restart,'finewoodlitter_deadmic_turnover_accumulated',&
             soil_finewoodlitter_deadmic_turnover_accumulated_ptr,'dead microbe fine wood litter carbon turnover', 'year-1')
        call add_tile_data(restart,'finewoodlitter_fsc_in',&
             soil_finewoodlitter_fsc_in_ptr,'fast fine wood litter carbon input', 'kg C/m2')
        call add_tile_data(restart,'finewoodlitter_ssc_in',&
             soil_finewoodlitter_ssc_in_ptr,'slow fine wood litter carbon input', 'kg C/m2')
        call add_tile_data(restart,'finewoodlitter_deadmic_in',&
             soil_finewoodlitter_deadmic_in_ptr,'dead microbe fine wood litter carbon input', 'kg C/m2')

        call add_tile_data(restart,'coarsewoodlitter_fast_turnover_accumulated',&
             soil_coarsewoodlitter_fast_turnover_accumulated_ptr,'fast coarse wood litter carbon turnover', 'year-1')
        call add_tile_data(restart,'coarsewoodlitter_slow_turnover_accumulated',&
             soil_coarsewoodlitter_slow_turnover_accumulated_ptr,'slow coarse wood litter carbon turnover', 'year-1')
        call add_tile_data(restart,'coarsewoodlitter_deadmic_turnover_accumulated', &
             soil_coarsewoodlitter_deadmic_turnover_accumulated_ptr,'dead microbe coarse wood litter carbon turnover', 'year-1')
        call add_tile_data(restart,'coarsewoodlitter_fsc_in', &
             soil_coarsewoodlitter_fsc_in_ptr,'fast coarse wood litter carbon input', 'kg C/m2')
        call add_tile_data(restart,'coarsewoodlitter_ssc_in', &
             soil_coarsewoodlitter_ssc_in_ptr,'slow coarse wood litter carbon input', 'kg C/m2')
        call add_tile_data(restart,'coarsewoodlitter_deadmic_in',&
             soil_coarsewoodlitter_deadmic_in_ptr,'dead microbe coarse wood litter carbon input', 'kg C/m2')
     endif
     call save_land_restart(restart)
     call free_land_restart(restart)
  endif
end subroutine save_soil_restart

! ============================================================================
! compute beta function
! after Manabe (1969), but distributed vertically.
subroutine soil_data_beta ( soil, vegn, diag, soil_beta, soil_water_supply, &
                            soil_uptake_T, soil_rh, soil_rh_psi )
  type(soil_tile_type), intent(inout) :: soil
  type(vegn_tile_type), intent(in)    :: vegn
  type(diag_buff_type), intent(inout) :: diag
  real, intent(out) :: soil_beta
  real, intent(out) :: soil_water_supply ! max rate of water supply to roots, kg/(m2 s)
  real, intent(out) :: soil_uptake_T ! an estimate of temperature of the water
             ! taken up by transpiration. In case of 'linear' uptake it is an exact
             ! value; in case of 'darcy*' treatments the actual uptake profile
             ! is calculated only in step 2, so the value returned is an estimate
  real, intent(out) :: soil_rh
  real, intent(out) :: soil_rh_psi

  ! ---- local vars
  integer :: l, iter
  real, dimension(num_l) :: &
       uptake_frac_max, & ! root distribution
       vegn_uptake_term, &
       vlc, vsc, & ! volumetric fractions of water and ice in the layer
       VRL, & ! vertical distribution of volumetric root length, m/m3
       u, du  ! uptake and its derivative (the latter is not used)
  real :: psi_for_rh
  real :: K_r, r_r ! root properties
  real :: psi_crown_min, grav_head, plant_height, xylem_resist, xylem_area_frac, sws, dum4
  real :: DsapDpsi, psi_left, psi_right, psi_mid, f_left, f_right, f_last, f_mid
  real :: z  !  soil depth

  call vegn_uptake_profile (vegn, dz(1:num_l), uptake_frac_max, vegn_uptake_term )

  vlc=0;vsc=0
  do l = 1, num_l
    vlc(l) = max(0., soil%wl(l) / (dens_h2o*dz(l)))
    vsc(l) = max(0., soil%ws(l) / (dens_h2o*dz(l)))
  enddo

  soil%uptake_frac = 0
  do l = 1, num_l
     soil%uptake_frac(l) = uptake_frac_max(l) &
          * max(0.0, min(1.0,(vlc(l)-soil%w_wilt(l))/&
               (0.75*(soil%w_fc(l)-soil%w_wilt(l)))))
  enddo
  soil_beta = sum(soil%uptake_frac)
  do l = 1, num_l
     if (soil_beta /= 0) then
          soil%uptake_frac(l) = soil%uptake_frac(l) / soil_beta
     else
          soil%uptake_frac(l) = uptake_frac_max(l)
     endif
  enddo
  if (lm2) soil%uptake_frac = uptake_frac_max

  ! calculate relative humidity at soil surface
  call soil_data_psi_for_rh ( soil, vlc, vsc, soil%psi, psi_for_rh )
  soil_rh = exp(psi_for_rh*g_RT)
  soil_rh_psi = g_RT*soil_rh

  ! calculate total water supply
  select case (uptake_option)
  case(UPTAKE_LINEAR)
     soil_water_supply = 0
     z = 0
     do l = 1, num_l
        soil_water_supply = soil_water_supply + &
          vegn_uptake_term(l)*max(0.0,soil%wl(l)/dz(l)-soil%w_wilt(l)*dens_h2o)
        z = z + dz(l)
     enddo
     soil_water_supply = z * soil_water_supply
     soil_water_supply = soil_water_supply/delta_time
     soil_uptake_T = sum(soil%uptake_frac*soil%T)
  case(UPTAKE_DARCY2D, UPTAKE_DARCY2D_LIN)
     call vegn_hydraulic_properties (vegn, dz(1:num_l), always_use_bsw, &
                                     VRL, K_r, r_r, &
                                     plant_height, xylem_area_frac, xylem_resist, &
                                     dum4 )
     grav_head = 0
     iter = 0
     psi_mid = 0
     if (use_tall_plants) grav_head = plant_height
     if (xylem_area_frac.gt.0. .and. xylem_resist.gt.0. .and. plant_height.gt.0. ) then
        DsapDpsi = DENS_H2O*xylem_area_frac/(plant_height*xylem_resist)
        psi_left = psi_wilt + grav_head
        call darcy2d_uptake ( soil, psi_left, VRL, &
               K_r, r_r, uptake_option, uptake_oneway, uptake_from_sat, u, du )
        f_left = max(0.0,sum(u)) - (psi_left-psi_wilt-grav_head)*DsapDpsi
        psi_right = 0.
        call darcy2d_uptake ( soil, psi_right, VRL, &
               K_r, r_r, uptake_option, uptake_oneway, uptake_from_sat, u, du )
        f_right = max(0.0,sum(u)) - (psi_right-psi_wilt-grav_head)*DsapDpsi
        f_last = f_left
        do iter = 1, max_iter_trans
           if (is_watch_point()) then
               write(*,*)'##### sws iteration iter=',iter
               __DEBUG5__(psi_left,f_left,psi_right,f_right,f_last)
           endif
           psi_mid = 0.5*(psi_left+psi_right)
           call darcy2d_uptake ( soil, psi_mid, VRL, &
                  K_r, r_r, uptake_option, uptake_oneway, uptake_from_sat, u, du )
           f_mid = max(0.0,sum(u)) - (psi_mid-psi_wilt-grav_head)*DsapDpsi
           if (abs(f_mid-f_last).lt.eps_trans) exit
           f_last = f_mid
           if (f_mid.lt.0.) then
               psi_right = psi_mid
               f_right = f_mid
           else
               psi_left = psi_mid
               f_left = f_mid
           endif
        enddo
     else
        call darcy2d_uptake ( soil, psi_wilt+grav_head, VRL, &
              K_r, r_r, uptake_option, uptake_oneway, uptake_from_sat, u, du )
     endif
     soil_water_supply = max(0.0,sum(u))
     soil_uptake_T = soil%uptake_T
     call send_tile_data(id_sws_n_iter, real(iter), diag)
     call send_tile_data(id_psi_x0_sws, psi_mid, diag)
  end select
end subroutine soil_data_beta

! ============================================================================
subroutine soil_sfc_water(soil, grnd_liq, grnd_ice, grnd_subl, grnd_tf)
  type(soil_tile_type), intent(in) :: soil
  real, intent(out) :: &
     grnd_liq, grnd_ice, & ! surface liquid and ice, respectively, kg/m2
     grnd_subl, &          ! fraction of vapor flux that sublimates
     grnd_tf               ! freezing temperature

  grnd_liq = 0
  grnd_ice = 0
  grnd_subl = soil_subl_frac(soil)
  ! set soil freezing temperature
  grnd_tf = soil%pars%tfreeze
end subroutine soil_sfc_water

! ============================================================================
real function soil_subl_frac(soil)
  type(soil_tile_type), intent(in) :: soil

  real :: grnd_liq, grnd_ice ! surface liquid and ice, respectively, kg/m2
  grnd_liq  = max(soil%wl(1), 0.)
  grnd_ice  = max(soil%ws(1), 0.)
  soil_subl_frac = 0.0
  if (grnd_liq + grnd_ice > 0) &
      soil_subl_frac = grnd_ice / (grnd_liq + grnd_ice)
end function soil_subl_frac

! ============================================================================
subroutine soil_evap_limits(soil, soil_E_min, soil_E_max)
  type(soil_tile_type), intent(in) :: soil
  real, intent(out) :: soil_E_min, soil_E_max

  real :: vlc
  soil_E_min = Eg_min
  if (use_E_max) then
     vlc = max(0.0, soil%wl(1) / (dens_h2o * dz(1)))
     soil_E_max = (soil%pars%k_sat_ref*soil%alpha(1)**2) &
               * (-soil%pars%psi_sat_ref/soil%alpha(1)) &
               * ((4.+soil%pars%chb)*vlc/ &
                ((3.+soil%pars%chb)*soil%pars%vwc_sat))**(3.+soil%pars%chb) &
                / ((1.+3./soil%pars%chb)*dz(1))
  else
     soil_E_max =  HUGE(soil_E_max)
  endif
end subroutine soil_evap_limits

! ============================================================================
! update soil properties explicitly for time step.
! MAY WISH TO INTRODUCE 'UNIVERSAL' SENSITIVITIES FOR SIMPLICITY.
! T-DEPENDENCE OF HYDRAULIC PROPERTIES COULD BE DONE LESS FREQUENTLY.
! integrate soil-heat conduction equation upward from bottom of soil
! to surface, delivering linearization of surface ground heat flux.
subroutine soil_step_1 ( soil, vegn, diag, &
                         soil_uptake_T, soil_beta, soil_water_supply, &
                         soil_rh, soil_rh_psi, &
                         soil_G0, soil_DGDT )
  type(soil_tile_type), intent(inout) :: soil
  type(vegn_tile_type), intent(in)    :: vegn
  type(diag_buff_type), intent(inout) :: diag
  real, intent(out) :: &
       soil_uptake_T, & ! estimate of the temperature of the water taken up by transpiration
       soil_beta, &
       soil_water_supply, & ! supply of water to vegetation per unit total active root biomass, kg/m2
       soil_rh,   & ! soil surface relative humidity
       soil_rh_psi,& ! derivative of soil_rh w.r.t. soil surface matric head
       soil_G0, soil_DGDT ! linearization of ground heat flux
  ! ---- local vars
  real :: bbb, denom, dt_e
  real, dimension(num_l) :: aaa, ccc, thermal_cond, heat_capacity, vlc, vsc
  integer :: l

  if(is_watch_point()) then
     write(*,*) 'soil%tag', soil%tag
     write(*,*) 'soil%pars%k_sat_ref', soil%pars%k_sat_ref
     write(*,*) 'soil%pars%psi_sat_ref', soil%pars%psi_sat_ref
     write(*,*) 'soil%pars%chb', soil%pars%chb
     write(*,*) 'soil%pars%w_sa', soil%pars%vwc_sat
  endif
! ----------------------------------------------------------------------------
! in preparation for implicit energy balance, determine various measures
! of water availability, so that vapor fluxes will not exceed mass limits
! ----------------------------------------------------------------------------

  call soil_data_beta ( soil, vegn, diag, soil_beta, soil_water_supply, soil_uptake_T, &
                        soil_rh, soil_rh_psi )

  do l = 1, num_l
     vlc(l) = max(0.0, soil%wl(l) / (dens_h2o * dz(l)))
     vsc(l) = max(0.0, soil%ws(l) / (dens_h2o * dz(l)))
  enddo
  call soil_data_thermodynamics ( soil, vlc, vsc, thermal_cond )

  do l = 1, num_l
     heat_capacity(l) = soil%heat_capacity_dry(l) *dz(l) &
          + clw*soil%wl(l) + csw*soil%ws(l)
  enddo

  if(num_l > 1) then
     do l = 1, num_l-1
        dt_e = 2 / ( dz(l+1)/thermal_cond(l+1) &
                     + dz(l)/thermal_cond(l)   )
        aaa(l+1) = - dt_e * delta_time / heat_capacity(l+1)
        ccc(l)   = - dt_e * delta_time / heat_capacity(l)
     enddo

     bbb = 1.0 - aaa(num_l)
     denom = bbb
     dt_e = aaa(num_l)*(soil%T(num_l) - soil%T(num_l-1)) &
               + soil%geothermal_heat_flux * delta_time / heat_capacity(num_l)
     soil%e(num_l-1) = -aaa(num_l)/denom
     soil%f(num_l-1) = dt_e/denom
     do l = num_l-1, 2, -1
        bbb = 1.0 - aaa(l) - ccc(l)
        denom = bbb + ccc(l)*soil%e(l)
        dt_e = - ( ccc(l)*(soil%T(l+1) - soil%T(l)  ) &
                  -aaa(l)*(soil%T(l)   - soil%T(l-1)) )
        soil%e(l-1) = -aaa(l)/denom
        soil%f(l-1) = (dt_e - ccc(l)*soil%f(l))/denom
     end do
     denom = delta_time/(heat_capacity(1) )
     soil_G0   = ccc(1)*(soil%T(2)- soil%T(1) + soil%f(1)) / denom
     soil_DGDT = (1 - ccc(1)*(1-soil%e(1))) / denom
  else  ! one-level case
     denom = delta_time/heat_capacity(1)
     soil_G0    = 0.
     soil_DGDT  = 1. / denom
  end if

  if(is_watch_point()) then
     write(*,*) '#### soil_step_1 checkpoint 1 ####'
     write(*,*) 'mask    ', .true.
     write(*,*) 'uptake_T', soil_uptake_T
     write(*,*) 'beta    ', soil_beta
     write(*,*) 'rh      ', soil_rh
     write(*,*) 'G0      ', soil_G0
     write(*,*) 'DGDT    ', soil_DGDT
     __DEBUG1__(soil_water_supply)
  endif

  call send_tile_data(id_thermal_cond, thermal_cond, diag)

end subroutine soil_step_1


! ============================================================================
! apply boundary flows to soil water and move soil water vertically.
  subroutine soil_step_2 ( soil, vegn, diag, snow_lprec, snow_hlprec,  &
                           vegn_uptk, &
                           subs_DT, subs_M_imp, subs_evap, &
                           use_tfreeze_in_grnd_latent, &
                           soil_levap, soil_fevap, soil_melt, &
                           soil_lrunf, soil_hlrunf, soil_Ttop, soil_Ctop, &
                           soil_frunf, soil_hfrunf, soil_tr_runf)
  type(soil_tile_type), intent(inout) :: soil
  type(vegn_tile_type), intent(in)    :: vegn
  type(diag_buff_type), intent(inout) :: diag
  real, intent(in) :: &
       snow_lprec, & ! ?? solid / liquid throughfall infiltrating the snow [mm/s]
       snow_hlprec, & ! ?? heat associated with snow_lprec [W/m^2]
       vegn_uptk, &  ! vegetation soil water uptake flux [mm/s]
       subs_DT,       & ! ?? soil surface layer temperature tendency [K]
       subs_M_imp,       &! rate of phase change of non-evaporated soil water ?? [mm/s]
       subs_evap         ! ?? solution for soil surface evaporation [mm/s]
  logical, intent(in) :: use_tfreeze_in_grnd_latent
  real, intent(out) :: &
       soil_levap, & ! ?? liquid soil surface evaporation [mm/s]
       soil_fevap, & ! ?? solid soil surface sublimation [mm/s]
       soil_melt, &  ! ?? net melt of non-evaporated soil ice over timestep [mm]
       soil_lrunf, & ! ?? total liquid runoff from soil [mm/s]
       soil_hlrunf, & ! ?? heat associated with runoff from soil [W/m^2]
       soil_Ttop, & ! ?? soil surface layer temperature [K]
       soil_Ctop, & ! ?? soil surface layer heat capacity [J/m^2.K]
       soil_frunf, & ! ?? frozen runoff from soil [mm/s]
       soil_hfrunf, & ! ?? heat associated with frozen runoff from soil [W/m^2]
       soil_tr_runf(:) ! dissolved organic carbon runoff from soil [kgC/m^2/s]

  ! ---- local vars ----------------------------------------------------------
  real, dimension(num_l)   :: del_t, &! ?? temperature tendency [K]
       psi, &! soil moisture potential wrt local elevation [m]
       DThDP, & ! ?? deriv. of vol. liq. cont. wrt psi [1/m]
       K_z, K_x, &! soil horiz. and vert. hydraulic conductivity, respectively [mm/s]
       DKDP, & ! ?? deriv. of hyd. cond. wrt psi [kg/m^3]
       vlc, & ! volumetric liquid water content [-]
       vsc, & ! volumetric solid water content [-]
       dW_l, & ! tendency of soil water mass [mm]
       DPsi
  real, dimension(num_l+1) :: flow, & ! downwards flow at layer interface above [mm/timestep]
       infilt
  real, dimension(num_l  ) :: div, & ! total divergence of soil water [mm/s]
       div_it    ! divergence of water due to inter-tile flow (incl. to stream)
  ! set in hlsp_hydrology_1 [mm/s]
  real, dimension(num_l  ) :: hdiv_it, &! divergence of heat due to inter-tile water flow [W/m^2]
       div_bf, & ! baseflow [mm/s]
       div_if, & ! interlow [mm/s]
       div_al, & ! div from active layer [mm/s]
       dq, div_active, &
       air_depth, macro_frac, extra

  real  :: &
       lprec_eff, & ! infiltrating throughfall (less saturated runoff) [mm/s], and
       hlprec_eff, & !  associated heat [W/m^2]
       tflow, & ! assumed temperature of liquid flow entering top layer [K]
       hcap, & ! heat capacity [J/m^2.K]
       dheat, &
       melt_per_deg, melt, macro_inf, &
       extra_cum, ddW, sum_air, denom, h1, h2, &
       flow_macro, & ! infiltration-excess runoff penetrating the soil due to macroporosity [mm/s]
       lrunf_sn, & ! runoff from saturated ground [mm/s]
       lrunf_ie,lrunf_bf,lrunf_if,lrunf_al, & ! ??, baseflow (incl inter-tile), interflow, active-layer runoff [mm/s]
       lrunf_nu, & ! runoff due to numerical saturation excess [mm/s]
       lrunf_sc, & ! sub-column runoff [mm/s]
       frunf, & ! frozen runoff (due to sat excess) [mm/s]
       hfrunf, & ! heat from frozen runoff (due to sat excess) [W/m^2]
       soil_subl, & ! fraction of water vapor that leaves surface as sublimation
       d_GW, &
       hlrunf_sn,hlrunf_ie,hlrunf_bf,hlrunf_if,hlrunf_al,hlrunf_nu,hlrunf_sc, & ! runoff heat fluxes [W/m^2]
       c0, c1, c2, Dpsi_min, Dpsi_max, &
       sat_area_frac, sat_thick, sum_trans, &
       gw_flux, depth_to_wt2_3, &
       depth_to_gw_flow, depth_to_gw_flow_3, &
       active_layer_thickness, d_psi, d_psi_s, psi_star, &
       depth_to_cf_1, depth_to_cf_3, &
       depth_to_wt_1, depth_to_wt_2, depth_to_wt_2a, depth_to_wt_2b, depth_to_wt_3, depth_to_wt_4, &
       storage_2, deficit_2, deficit_3, deficit_4, deficit, dum1, dum2, dum3
  logical :: stiff
  logical :: hlsp_stiff ! were all the horizontal conductivities zero at the call to hlsp_hydrology_1?
  real :: zimh, ziph, dTr_g(num_l), dTr_s(num_l)
  integer :: n_iter, l, l_max_active_layer
  real :: &
       VRL(num_l), & ! volumetric root length, m/m3
       K_r, & ! root membrame permeability, kg/(m3 s)
       r_r, & ! root radius, m
       bwood, & ! heartwood biomass kg C/m2
       uptake(num_l),   & ! uptake by roots per layer, kg/(m2 s)
       uptake_pos,      & ! sum of the positive uptake, kg/(m2 s)
       uptake_T_new, & ! updated average temperature of uptaken water, deg K
       uptake_T_corr,& ! correction for uptake temperature, deg K
       Tu,           & ! temperature of water taken up from (or added to) a layer, deg K
       psi_x0          ! water potential inside roots (in xylem) at zero depth, m
  ! For testing tridiagonal solution
  real, dimension(num_l)   :: t_soil_tridiag ! soil temperature based on generic tridiagonal solution [K]
  real, dimension(num_l)   :: t_diff ! difference from original advection subroutine [K]

  real :: DOC_leached(n_c_types,num_l), div_DOC_loss(n_c_types,num_l),  &     ! C leaching
         leaflitter_DOC_loss(n_c_types),woodlitter_DOC_loss(n_c_types)        ! Surface litter C leaching loss

  real :: surface_water ! diagnostic surface water storage [m]
  real :: inundated_frac ! diagnostic inundated area fraction [-]
  real :: wet_frac ! diagnostic wetland area fraction
  real :: sat_runf_frac ! [-] effective saturated fraction used for calculating surface runoff, used
  ! when simple inundation scheme is applied
  real :: flow_s(num_l) ! flow downwards at layer interfaces [mm/s] (flow/delta_time excluding bottom interface)
  real :: reflux        ! [mm/s] upwards flow at soil surface in excess of rejected throughfall
#ifdef ZMSDEBUG
  ! Checking diagnostics from Richards that are used in advection solution.
  real :: w1, w2
#endif
  real, dimension(num_l)   :: wl_before ! water content before call to Richards [mm]
  real, parameter :: wthresh = 1.e-9   ! [mm] tolerance for roundoff error for water balance
  real :: wsum1, wsum2   ! total water stocks for balance check [mm, or kg/m^2]
  real :: sliq, sice     ! for call to soil_tile_stock_pe
  character(len=512) :: mesg ! for error message
  integer :: ipt, jpt, kpt, fpt ! for debug
  real :: surf_DOC_loss(n_c_types)! [kg C/m^2] DOC loss from top soil layer to surface runoff due
                                  ! to efflux
  real :: total_C_leaching(num_l) ! [kg C/m^2/s] net total vertical DOC leaching by layer
  real :: total_DOC_div           ! [kg C/m^2/s] net total DOC divergence loss rate
  ! --------------------------------------------------------------------------
  div_active(:) = 0.0

  !.........................................................................
  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 1 #####'
     write(*,*) 'mask    ', .true.
     __DEBUG1__(subs_evap)
     __DEBUG1__(snow_lprec)
     __DEBUG1__(vegn_uptk)
     __DEBUG1__(subs_M_imp)
     call dpri('theta_s ',soil%pars%vwc_sat); write(*,*)
     do l = 1, num_l
        write(*,'(a,i2.2)',advance='NO') 'level=', l
        call dpri(' T =', soil%T(l))
        call dpri(' Th=', (soil%ws(l)+soil%wl(l))/(dens_h2o*dz(l)))
        call dpri(' wl=', soil%wl(l))
        call dpri(' ws=', soil%ws(l))
        call dpri(' gw=', soil%groundwater(l))
        write(*,*)
     enddo
  endif
  !.........................................................................

!  if (do_component_balchecks) then
  ! Sum total water mass at beginning of soil_step_2
  call soil_tile_stock_pe (soil, sliq, sice )
  wsum1 = sliq + sice
!  end if

  ! ---- record fluxes -----------------------------------------------------
  soil_subl = soil_subl_frac(soil)
  soil_levap  = subs_evap*(1-soil_subl)
  soil_fevap  = subs_evap*   soil_subl
  soil_melt   = subs_M_imp / delta_time

  ! ---- load surface temp change and perform back substitution ------------
  del_t(1) = subs_DT
  soil%T(1) = soil%T(1) + del_t(1)
  if ( num_l > 1) then
     do l = 1, num_l-1
        del_t(l+1) = soil%e(l) * del_t(l) + soil%f(l)
        soil%T(l+1) = soil%T(l+1) + del_t(l+1)
     enddo
  endif

  !.........................................................................
  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 2 #####'
     do l = 1, num_l
        write(*,'(a,i2.2)',advance='NO') 'level=',l
        call dpri('T=', soil%T(l))
        call dpri('del_t=', del_t(l))
        call dpri('e=', soil%e(l))
        call dpri('f=', soil%f(l))
        write(*,*)
     enddo
  endif
  !.........................................................................

  ! ---- extract evap from soil, adjusting T, and do implicit melt ---------
  IF (LM2) THEN  ! (extract surface E--is there any?--uniformly from bucket)
     do l = 1, num_l
        soil%wl(l) = soil%wl(l) &
                      - soil%uptake_frac(l)*soil_levap*delta_time
     enddo
  ELSE
     soil%wl(1) = soil%wl(1) - soil_levap*delta_time
     soil%ws(1) = soil%ws(1) - soil_fevap*delta_time
  ENDIF
  hcap = soil%heat_capacity_dry(1)*dz(1) &
                       + clw*soil%wl(1) + csw*soil%ws(1)
  ! T adjustment for nonlinear terms (del_T)*(del_W)
  dheat = delta_time*(clw*soil_levap+csw*soil_fevap)*del_T(1)
  ! take out extra heat not claimed in advance for evaporation
  if (use_tfreeze_in_grnd_latent) dheat = dheat &
          - delta_time*((cpw-clw)*soil_levap+(cpw-csw)*soil_fevap) &
                                 *(soil%T(1)-del_T(1)-tfreeze)
  soil%T(1)  = soil%T(1)  + dheat/hcap
  soil%wl(1) = soil%wl(1) + subs_M_imp
  soil%ws(1) = soil%ws(1) - subs_M_imp
  soil%T(1)  = tfreeze + (hcap*(soil%T(1)-tfreeze) ) &
                              / ( hcap + (clw-csw)*subs_M_imp )

  ! ---- calculate actual uptake and update its T --------------------------
  psi_x0 = 0
  select case(uptake_option)
  case ( UPTAKE_LINEAR )
     uptake_T_corr = 0
     n_iter = 0
     uptake = soil%uptake_frac*vegn_uptk
     soil%psi_x0 = 0.
     bwood = 0.
  case ( UPTAKE_DARCY2D, UPTAKE_DARCY2D_LIN )
     ! for Darcy-flow uptake, find the root water potential to satify actual
     ! transpiration by the vegetation
     ! **** introduce optional arguments for dum1 etc ? ********
     call vegn_hydraulic_properties (vegn, dz(1:num_l), always_use_bsw, &
                               VRL, K_r, r_r, dum1, dum2, dum3, bwood)

     call darcy2d_uptake_solver (soil, vegn_uptk, VRL, K_r, r_r, &
             uptake_option, uptake_oneway, uptake_from_sat, uptake, psi_x0, n_iter)
     soil%psi_x0 = psi_x0

     uptake_pos = sum(uptake(:),mask=uptake(:)>0)
     if (uptake_pos > 0) then
        ! calculate actual temperature of uptake
        uptake_T_new  = sum(uptake*soil%T,mask=uptake>0)/uptake_pos
        ! and temperature correction
        uptake_T_corr = soil%uptake_T - uptake_T_new
        if(is_watch_point()) then
           __DEBUG3__(soil%uptake_T, uptake_T_new, uptake_T_corr)
        endif
        ! save new uptake for the next time step to serve as an estimate of uptake
        ! temperature
        soil%uptake_T = uptake_T_new
     else
        uptake_T_corr = 0.0
        ! and do not change the soil%uptake_T
     endif
  case default
     call error_mesg('soil_step_2', 'invalid soil uptake option', FATAL)
  end select

  !.........................................................................
  if (is_watch_point())then
     write(*,*) ' ##### soil_step_2 checkpoint 2.1 #####'
     __DEBUG2__(vegn_uptk,sum(uptake))
     do l = 1,num_l
        write(*,'(a,i2.2)',advance='NO')'level=',l
        call dpri('uptake=',uptake(l))
        call dpri('dwl=',-uptake(l)*delta_time)
        call dpri('wl=',soil%wl(l))
        call dpri('new wl=',soil%wl(l) - uptake(l)*delta_time)
        write(*,*)
     enddo
  endif
  !.........................................................................

  call send_tile_data(id_uptk_n_iter, real(n_iter), diag)
  call send_tile_data(id_uptk, uptake, diag)
  call send_tile_data(id_psi_x0, psi_x0, diag)

  ! ---- perform the uptake ------------------------------------------------
  do l = 1, num_l
     ! calculate the temperature of water that is taken from the layer (or added
     ! to the layer), including energy balance correction
     if (uptake(l) > 0) then
        Tu = soil%T(l) + uptake_T_corr
     else
        Tu = soil%uptake_T + uptake_T_corr
     endif
     hcap = soil%heat_capacity_dry(l)*dz(l) &
           + clw*soil%wl(l) + csw*soil%ws(l)
     soil%T(l) = soil%T(l) - &
           uptake(l)*delta_time*clw*( Tu-soil%T(l) ) / &
           ( hcap - uptake(l)*delta_time*clw )
     soil%wl(l) = soil%wl(l) - uptake(l)*delta_time
  enddo

  !.........................................................................
  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 3 #####'
     do l = 1, num_l
        write(*,'(x,a,x,i2.2)',advance='NO')' level=', l
        call dpri(' T =',soil%T(l))
        call dpri(' Th=', (soil%ws(l)+soil%wl(l))/(dens_h2o*dz(l)))
        call dpri(' wl=', soil%wl(l))
        call dpri(' ws=', soil%ws(l))
        write(*,*)
     enddo
  endif
  !.........................................................................

  ! ---- push down any excess surface water, with heat ---------------------
  IF (PUSH_DOWN_SFC_EXCESS) THEN
     CALL SOIL_PUSH_DOWN_EXCESS ( soil, diag, lrunf_nu, hlrunf_nu, frunf, hfrunf)
  ELSE
     lrunf_nu=0; hlrunf_nu=0; frunf=0; hfrunf = 0
  ENDIF

  ! ---- fetch soil hydraulic properties -----------------------------------
  do l = 1, num_l
    vlc(l) = max(0., soil%wl(l) / (dens_h2o*dz(l)))
    vsc(l) = max(0., soil%ws(l) / (dens_h2o*dz(l)))
  enddo
  call soil_data_hydraulic_properties (soil, vlc, vsc, &
                   psi, DThDP, K_z, K_x, DKDP, Dpsi_min, Dpsi_max )

  ! ---- compute various measures of water table depth ---------------------
  sat_thick = 0.
  do l=num_l,1,-1
     if(vsc(l)+vlc(l).le.soil%pars%vwc_sat) exit
     sat_thick = sat_thick + dz(l)
  enddo
  depth_to_cf_1 = zhalf(num_l+1) - sat_thick
  depth_to_wt_1 = depth_to_cf_1 - soil%pars%psi_sat_ref/soil%alpha(max(l,1))

  depth_to_wt_2 = zfull(num_l)-psi(num_l)

  depth_to_wt_2a = 0.
  depth_to_wt_3  = 0.0
  do l=1,num_l
     if (soil%wl(l)+soil%ws(l) .lt. &
                      soil%pars%vwc_sat*dens_h2o*dz(l)) then
        depth_to_wt_2a = depth_to_wt_2a + dz(l)
        if (gw_option /= GW_TILED .and. l.eq.num_l) depth_to_wt_2a = -1.
        ! ZMS changed definition here so time-average will remain deep water table
        ! and for use in calculating inundated fraction in hlsp_hydrology_2
     else
        exit
     endif
  enddo
  depth_to_wt_2b = 0. ! This version requires an entirely unfrozen saturated layer
  do l=1,num_l
     if (soil%wl(l) .lt. soil%pars%vwc_sat*dens_h2o*dz(l)) then
        depth_to_wt_2b = depth_to_wt_2b + dz(l)
     else
        exit
     end if
  end do

  storage_2 = 1 - depth_to_wt_2  &
           /(soil%pars%hillslope_zeta_bar*soil%pars%hillslope_relief)
  storage_2 = min( max( 0., storage_2 ) , 1.)
  deficit_2 = 1 - storage_2

  if (vsc(num_l).gt.0.) then   ! permafrost
     if (gw_option /= GW_TILED) then
        depth_to_wt2_3 = 0.
        depth_to_cf_3 = 0.
        depth_to_wt_3 = 0.
     else
        depth_to_wt2_3 = zhalf(l+1)
        depth_to_cf_3 = zhalf(l+1)
        depth_to_wt_3 = zhalf(l+1)
     end if
  else                       ! liquid water at depth
     depth_to_cf_3 = 0.
     if (use_fringe) then
        do l = num_l, 1, -1
           if ( l.eq.num_l .and. psi(l).le.soil%pars%psi_sat_ref/soil%alpha(l) ) then
              depth_to_cf_3 = zfull(l) + soil%pars%psi_sat_ref/soil%alpha(l) - psi(l)
              exit
           else if (psi(l).le.soil%pars%psi_sat_ref/soil%alpha(l) ) then
              d_psi = psi(l+1) - psi(l)
              d_psi_s = (soil%pars%psi_sat_ref/soil%alpha(l+1)) &
                       -(soil%pars%psi_sat_ref/soil%alpha(l))
              psi_star = (psi(l)*d_psi_s - &
                          d_psi*(soil%pars%psi_sat_ref/soil%alpha(l)))&
                          / (d_psi_s - d_psi)
              depth_to_cf_3 = zfull(l) + (zfull(l+1)-zfull(l)) &
                                    * (psi_star-psi(l)) / d_psi
              exit
           else if (l.eq.1) then
              d_psi = psi(l+1) - psi(l)
              d_psi_s = (soil%pars%psi_sat_ref/soil%alpha(l+1)) &
                       -(soil%pars%psi_sat_ref/soil%alpha(l))
              psi_star = (psi(l)*d_psi_s - &
                          d_psi*(soil%pars%psi_sat_ref/soil%alpha(l)))&
                          / (d_psi_s - d_psi)
              depth_to_cf_3 = zfull(l) + (zfull(l+1)-zfull(l)) &
                                    * (psi_star-psi(l)) / d_psi
              depth_to_cf_3 = max(0.,depth_to_cf_3)
           endif
        enddo
     endif
     depth_to_wt_3 = max(0., zhalf(num_l+1)-(psi(num_l)+dz(num_l)/2.))
     depth_to_wt2_3 = depth_to_wt_3
  endif

  if (use_fringe) then
     depth_to_gw_flow_3 = depth_to_cf_3
  else
     depth_to_gw_flow_3 = depth_to_wt_3
  endif
  deficit_3 = depth_to_gw_flow_3  &
           /(soil%pars%hillslope_zeta_bar*soil%pars%hillslope_relief)

  if (vsc(min(num_l,layer_for_gw_switch)).gt.0.) then   ! permafrost
     if (gw_option /= GW_TILED) then
        depth_to_wt_4 = 0.
     else
        depth_to_wt_4 = zhalf(num_l+1)
     end if
  else ! liquid water at depth
     depth_to_wt_4 = 0.
     do l = num_l, 1, -1
        if ( l.eq.num_l .and. psi(l).le.0. ) then
           depth_to_wt_4 = zfull(l) - psi(l)
           exit
        else if (psi(l).le.0. ) then
           d_psi = psi(l+1) - psi(l)
           depth_to_wt_4 = zfull(l) - (zfull(l+1)-zfull(l)) &
                                 * (psi(l)) / d_psi
           exit
        else if (l.eq.1) then
           d_psi = psi(l+1) - psi(l)
           depth_to_wt_4 = zfull(l) - (zfull(l+1)-zfull(l)) &
                                 * (psi(l)) / d_psi
           depth_to_wt_4 = max(0.,depth_to_wt_4)
        endif
     enddo
  endif
  deficit_4 = depth_to_wt_4  &
           /(soil%pars%hillslope_zeta_bar*soil%pars%hillslope_relief)

! ---- get saturated area and column flow divergences --------------------
  SELECT CASE(gw_option)

  CASE(GW_LM2)

   div_bf=0; div_if=0; div_al=0; sat_area_frac = 0; div_it = 0.; hdiv_it = 0.

  CASE(GW_LINEAR)

      IF (CORRECTED_LM2_GW) THEN
         do l = 1, num_l
            if (vlc(l) .ge. soil%pars%vwc_sat .and. vsc(l).le.0.) &
               div_bf(l) = 0.15*dens_h2o*dz(l)/soil%pars%tau_groundwater
         enddo
      ELSE
         do l = 1, num_l
            if ((vsc(l)+vlc(l)) .ge. soil%pars%vwc_sat) &
                div_bf(l) = 0.15*dens_h2o*dz(l)*(vlc(l)/(vsc(l)+vlc(l)))  &
                                    /soil%pars%tau_groundwater
         enddo
      ENDIF
      div_if = 0
      div_al = 0; div_it = 0.; hdiv_it = 0.
      sat_thick = zhalf(num_l+1) - depth_to_cf_1
      sat_area_frac = min((sat_thick/zhalf(num_l+1))**soil%pars%rsa_exp,1.)

  CASE(GW_HILL_AR5)

      call soil_data_gw_hydraulics_ar5(soil, storage_2, &
                                               gw_flux, sat_area_frac)
      dq = 0.
      sat_thick = 0.
      do l=num_l,1,-1
         if(psi(l).le.0.) exit
         if (vsc(l).le.0.) dq(l) = dz(l)
         sat_thick = sat_thick + dz(l)
      enddo
      div_bf = 0.
      if (sat_thick.gt.0.) div_bf = (dq/sat_thick)*gw_flux

      div_active = 0.
      l_max_active_layer = 0
      do l=1,num_l
         if(vsc(l).gt.0.) exit
         l_max_active_layer = l
      enddo
      if (l_max_active_layer.lt.num_l .and. l_max_active_layer.gt.0) then
          do l = 1, l_max_active_layer
             if(vlc(l).gt.0) &
             div_active(l) = K_x(l) * soil%pars%hillslope_relief*dz(l) &
              / (soil%pars%hillslope_length*soil%pars%hillslope_length)
          enddo
      endif

      div_al = 0
      where (div_bf.eq.0.) div_al = div_active*active_layer_drainage_acceleration
   div_if = 0; div_it = 0.; hdiv_it = 0.

  CASE(GW_HILL)
     if (any(soil%wl+soil%ws.lt.0.) .and. prohibit_negative_water_div) then
        sat_area_frac = 0.
        div_bf = 0.
        div_if = 0.
     else if (vsc(min(num_l,layer_for_gw_switch)).gt.0.) then   ! permafrost
        sat_area_frac = 0.
        div_bf = 0.
        div_if = 0.
     else                       ! liquid water at depth
        if (use_depth_to_wt_4) then
           depth_to_gw_flow = depth_to_wt_4
           deficit = deficit_4
        else
           depth_to_gw_flow = depth_to_gw_flow_3
           deficit = deficit_3
        endif
        call soil_data_gw_hydraulics(soil, deficit, gw_flux, sat_area_frac)
        gw_flux = min(gw_flux, gw_flux_max)
        dTr_g = 0.
        dTr_s = 0.
        dTr_g(num_l) = 1.
        l = num_l
        ziph = sum(dz(1:num_l))
        zimh = ziph - dz(num_l)
        if (depth_to_gw_flow .lt. zimh) then
           dTR_g(l) = dz(l)
           dTr_s(l) = (exp(-zimh/soil%pars%soil_e_depth))
           do l = num_l-1, 1, -1
              if (vsc(l).gt.0.) exit
              ziph = zimh
              zimh = ziph - dz(l)
              if (depth_to_gw_flow .lt. zimh) then
                 dTR_g(l) = dz(l)
                 dTr_s(l) = exp(-zimh/soil%pars%soil_e_depth)-exp(-ziph/soil%pars%soil_e_depth)
              else if (depth_to_gw_flow .lt. ziph) then
                 dTR_g(l) =(ziph-depth_to_gw_flow)
                 dTr_s(l) = exp(-depth_to_gw_flow/soil%pars%soil_e_depth)-exp(-ziph/soil%pars%soil_e_depth)
              else
                 exit
              endif
           enddo
        endif
        sum_trans = sum(dTr_g)
        if (sum_trans.ne.0.) then
           dTR_g = dTR_g / sum_trans
           dTR_g = dTR_g * soil%pars%k_sat_gw*aspect*soil%pars%hillslope_length
        endif
        dTR_s = dTR_s * (soil%pars%k_sat_sfc+k0_macro_x)*soil%pars%soil_e_depth
        sum_trans = sum(dTR_g) + sum(dTr_s)
        if (sum_trans.ne.0.) then
           div_bf = gw_flux * dTR_g /sum_trans
           div_if = gw_flux * dTR_s /sum_trans
        else
           div_bf = 0.
           div_if = 0.
        endif
     endif

     div_al = 0; div_it = 0.; hdiv_it = 0.
     l_max_active_layer = 0   ! "active layer" either over permafrost or perched
     do l=1,num_l
        if(vsc(l).gt.0.) exit
        l_max_active_layer = l
     enddo
     if (l_max_active_layer.lt.num_l .and. l_max_active_layer.gt.0) then
        do l = 1, l_max_active_layer
           div_al(l) = active_layer_drainage_acceleration &
                 * K_x(l) * soil%pars%hillslope_relief*dz(l) &
                 / (soil%pars%hillslope_length*soil%pars%hillslope_length)
        enddo
     endif

  CASE (GW_TILED)
     call hlsp_hydrology_2(soil, psi, vlc, vsc, div_it, hdiv_it, &
        sat_area_frac, inundated_frac, storage_2, depth_to_wt_2b, surface_water, sat_runf_frac)
     if (depth_to_wt_2b >= 0. .and. depth_to_wt_2b <= wet_depth) then
        wet_frac = 1.
     else
        wet_frac = 0.
     end if
     div_if=0.
     div_al=0.
     div_bf=0.

  END SELECT

  div = div_bf + div_if + div_al + div_it ! div includes inter-tile flow
  lrunf_bf = sum(div_bf + div_it) ! baseflow runoff includes inter-tile flow
  lrunf_if = sum(div_if)
  lrunf_al = sum(div_al)

  if (snow_lprec.ne.0.) then
     if (gw_option == GW_TILED .and. simple_inundation) then
        ! use sat_runf_frac calculated in hlsp_hydrology_2. Includes microtopography and finite runoff rate.
        lrunf_sn = sat_runf_frac * snow_lprec
     else
        lrunf_sn = sat_area_frac * snow_lprec
        ! if gw_option == GW_TILED, sat_area_frac = 0 or 1 depending on saturation status of top layer
     end if
     hlrunf_sn = lrunf_sn*snow_hlprec/snow_lprec
  else
     lrunf_sn = 0.
     hlrunf_sn = 0.
  endif
  hlrunf_ie=0
  lprec_eff = snow_lprec - lrunf_sn
  hlprec_eff = snow_hlprec - hlrunf_sn

  if(is_watch_point()) then
     do l = 1, num_l
        write(*,'(a,1x,i2.2,100(2x,g23.16))')'div_ac,div_bf,div_if,div_al,div', &
                       l,div_active(l),div_bf(l),div_if(l),div_al(l),div(l)
     enddo
     do l = 1, num_l
        write(*,'(a,1x,i2.2,100(2x,g23.16))')'vsc,psi,dz',l,vsc(l),psi(l),dz(l)
     enddo
     write(*,*)'lrunf_bf',lrunf_bf
     write(*,*)'tau_gw',soil%pars%tau_groundwater
     write(*,*)'dens_h2o',dens_h2o
  endif

  wl_before(1:num_l) = soil%wl(1:num_l)

  ! ---- soil-water flow ----------------------------------------------------
  IF (LM2) THEN
     flow(1) = 0
     do l = 1, num_l
        infilt(l) = soil%uptake_frac(l)*lprec_eff *delta_time
        flow(l+1) = max(0., soil%wl(l) + flow(l) &
              + infilt(l) - soil%w_fc(l)*dz(l)*dens_h2o)
        dW_l(l) = flow(l) - flow(l+1) + infilt(l)
        soil%wl(l) = soil%wl(l) + dW_l(l)
     enddo
     do l = 1, num_l
        flow(l) = flow(l) + infilt(l)
     enddo
     dW_l=0
     dpsi=0
     c0 = delta_time/soil%pars%tau_groundwater
     c1 = exp(-c0)
     c2 = (1-c1)/c0
     l = 1
     d_GW = c1 * soil%groundwater(l) + c2 * flow(num_l+1) &
                           - soil%groundwater(l)
     soil%groundwater(l) = soil%groundwater(l) + d_GW
     lrunf_sc  = (1-c1)*soil%groundwater(l)/delta_time &
                         + (1-c2)*flow(num_l+1)/delta_time
     lrunf_ie=0
  ELSE
     lrunf_sc = 0
     d_GW = 0
     stiff = all(DThDP.eq.0)
     IF(stiff) THEN
        ! Note: in order to maintain energy and water conservation with between tile-flows
        ! in case of the tiled groundwater / full hillslope model, we require stiff ==>
        ! soil%T(l) < tfreeze for all l.
        ! ZMS: since this could change during timestep, impose some explicit updates in case this
        ! is now the case but the flows are nonzero.
        hlsp_stiff = all(div_it==0.) .and. all(hdiv_it==0.)
        if (is_watch_cell()) then
           write(*,*)'Stiff point in watch-cell at hidx_j=', soil%hidx_j
           __DEBUG2__(div_it, hdiv_it)
           if (hlsp_stiff) write(*,*) 'hlsp_stiff'
        end if
        if (.not. hlsp_stiff) then
           call stiff_explicit_gwupdate(soil, div_it, hdiv_it, div, lrunf_bf)
        end if
        flow = 0.
        dW_l = 0.
        if (.not. hlsp_stiff) then
           div  = 0.
           lrunf_bf = 0.
        end if
        div_bf=0.; div_if=0; div_al=0
        lrunf_if = 0; lrunf_al = 0
        lrunf_ie = lprec_eff
        hlrunf_ie = hlprec_eff
        psi=zfull(1:num_l)
        dpsi=0.
     ELSE
        if (USE_RICHARDS_CLEAN) then
           call RICHARDS_clean(soil, psi, DThDP, K_z, DKDP, div, &
                  lprec_eff, Dpsi_min, Dpsi_max, delta_time, &
                  dPsi, dW_l, flow, lrunf_ie)
        else
           call RICHARDS(soil, psi, DThDP, K_z, DKDP, div, &
                  lprec_eff, Dpsi_min, Dpsi_max, delta_time, &
                  dPsi, dW_l, flow, lrunf_ie)
        endif
        do l = 1, num_l
           soil%wl(l) = soil%wl(l) + dW_l(l)
        enddo
#ifdef ZMSDEBUG
        do l = 1, num_l
           w1 = wl_before(l)
           w2 = soil%wl(l) - dW_l(l)
           call check_conservation('soil_step_2: Richards Eqn. Diagnostics, wl_before', 'Water', w1, w2, &
                wthresh, WARNING)
           w1 = dW_l(l)
           if (l < num_l) then
              w2 = flow(l) - flow(l+1) - div(l)*delta_time
           else
              w2 = flow(l) - div(l)*delta_time
           end if
           call check_conservation('soil_step_2: Richards Eqn. Diagnostics, dW_l', 'Water', w1, w2, &
                wthresh, WARNING)
        end do
#endif
      ENDIF
  ENDIF

  ! Check for negative wl
  do l = 1, num_l
     if (soil%wl(l) < 0.) then
        if (soil%wl(l) / (dens_h2o*dz(l)*soil%pars%vwc_sat) < thetathresh) then
           call get_current_point(ipt, jpt, kpt, fpt)
           write(mesg,*) 'soil%wl(l) < 0! l,i,j,k,face:', l, ipt, jpt, kpt, fpt, '. degree of saturation = ', &
                soil%wl(l) / (dens_h2o*dz(l)*soil%pars%vwc_sat), '. If ".not. allow_neg_wl", '// &
                'model will abort.'
           if (.not. allow_neg_wl) then
              call error_mesg(module_name, mesg, FATAL)
           else
              if (verbose) call error_mesg(module_name, mesg, WARNING)
           end if
        end if
     end if
  end do

  ! Check for negative hcap
  if (require_pos_hcap) then
     do l = 1, num_l
        if (soil%wl(l) < 0.) then
            ! Make sure bulk heat capacity stays above zero
            hcap = soil%heat_capacity_dry(l)*dz(l) &
                   + clw*soil%wl(l) + csw*soil%ws(l)
            if (hcap .le. 0.) then
               call get_current_point(ipt, jpt, kpt, fpt)
               write(mesg,*) 'soil%wl(l) < 0! l,i,j,k,face:', l, ipt, jpt, kpt, fpt, '. degree of saturation = ', &
                      soil%wl(l) / (dens_h2o*dz(l)*soil%pars%vwc_sat), '. This makes hcap = ', hcap, &
                      ', which is < 0! Model Aborting!'
               call error_mesg(module_name, mesg, FATAL)
            end if
        end if
     end do
  endif
  ! Check total
  if (lrunf_nu < negrnuthresh) then
     call get_current_point(ipt, jpt, kpt, fpt)
     write(mesg,*) 'soil%wl(l) < 0 at one or more l at i,j,k,face:', ipt, jpt, kpt, fpt, '.', &
          ' Total lrunf_nu required to set to zero = ', lrunf_nu, ' mm/s.'
     call error_mesg(module_name, mesg, WARNING)
  end if


  ! ---- heat advection by water flow ---------------------------------------
  if  (snow_lprec.ne.0.) then
    tflow = tfreeze + snow_hlprec/(clw*snow_lprec)
  else
    tflow = tfreeze
  endif

  if(is_watch_point()) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.4 ***** '
     write(*,*) ' tfreeze', tfreeze
     write(*,*) '  tflow ', tflow
     write(*,*) ' snow_hlprec', snow_hlprec
  endif

#ifndef ZMSDEBUG_TRIDIAGTEST
  if (use_tridiag_foradvec) then
     call advection_tri(soil, flow, dW_l, tflow, d_GW, div, delta_time, t_soil_tridiag, hdiv_it)
  else
     call advection(soil, flow, dW_l, tflow, d_GW, div, delta_time)
  end if
#else
  call advection_tri(soil, flow, dW_l, tflow, d_GW, div, delta_time, t_soil_tridiag, hdiv_it)
  call advection(soil, flow, dW_l, tflow, d_GW, div, delta_time)
  ! Calculate difference between two solutions
  t_diff(1:num_l) = t_soil_tridiag(1:num_l) - soil%T(1:num_l)
  call send_tile_data(id_st_diff, t_diff, diag)
#endif

  if (lprec_eff.ne.0. .and. flow(1).ge.0. ) then
     hlrunf_ie = lrunf_ie*hlprec_eff/lprec_eff
  else if (flow(1).lt.0. ) then
     hlrunf_ie = hlprec_eff - (flow(1)/delta_time)*clw &
                         *(soil%T(1)-tfreeze)
  else
     hlrunf_ie = 0.
  endif

  ! Initialize for use in output below
  macro_inf = 0.
  extra_cum = 0.
  ! ---- allow infiltration-excess runoff to enter soil via macroporosity
  IF (BOTTOM_UP_COLD_INFILT.and.bwood.gt.bwood_macinf.and.lrunf_ie.gt.0.) then
     do l = 1, num_l
        air_depth(l) = soil%pars%vwc_sat*dz(l)-(soil%wl(l)+soil%ws(l))/dens_h2o
        air_depth(l) = min(air_depth(l), soil%pars%vwc_sat*dz(l))
        air_depth(l) = max(air_depth(l), 0.)
        if (zfull(l).gt.COLD_DEPTH) air_depth(l) = 0.
     enddo
     sum_air = sum(air_depth)
     macro_inf = min(lrunf_ie*delta_time,sum_air*dens_h2o)
     tflow = tfreeze + hlrunf_ie/(clw*lrunf_ie)
     ! compute the fill amounts from bottom up
     extra_cum = macro_inf
     dW_l = 0.
     do l = num_l, 1, -1
        if (extra_cum.gt.0.) then
           dW_l(l) = min(extra_cum, air_depth(l)*dens_h2o)
           extra_cum = extra_cum - dW_l(l)
        else
           exit
        endif
     enddo
     ! place the water and heat
     do l= 1, num_l
        h1 = soil%heat_capacity_dry(l)*dz(l) &
             + csw*soil%ws(l) + clw*soil%wl(l)
        h2 = clw*dW_l(l)
        soil%T(l) = (h1 * soil%T(l) &
                      + h2 * tflow )  / (h1+h2)
        soil%wl(l) = soil%wl(l) + dW_l(l)
     enddo
     ! extra_cum should now be 0, but included below for safety
     lrunf_ie = lrunf_ie - (macro_inf-extra_cum)/delta_time
     hlrunf_ie = hlrunf_ie - clw*(tflow-tfreeze)*(macro_inf-extra_cum)/delta_time
  ELSE if (cold_infilt.and.lrunf_ie.gt.0.) then
     do l = 1, num_l
        air_depth(l) = soil%pars%vwc_sat*dz(l)-(soil%wl(l)+soil%ws(l))/dens_h2o
        air_depth(l) = min(air_depth(l), soil%pars%vwc_sat*dz(l))
        air_depth(l) = max(air_depth(l), 0.)
        macro_frac(l) = exp(-zhalf(l)  /soil%pars%soil_e_depth) &
                       -exp(-zhalf(l+1)/soil%pars%soil_e_depth)
     enddo
     sum_air = sum(air_depth)
     macro_inf = min(lrunf_ie*delta_time,sum_air*dens_h2o)
     tflow = tfreeze + hlrunf_ie/(clw*lrunf_ie)
     if (sum_air.gt.0.) then
        denom = sum(air_depth*macro_frac/dz(1:num_l))
        ! tentatively distribute macropore infiltration in proportion
        ! to available pore volume and macropore density
        do l= 1, num_l
           dW_l(l) = macro_inf*(air_depth(l)*macro_frac(l)/dz(l))/denom
        enddo
        extra = 0.
        ! find excess infiltration by layer
        do l = 1, num_l
           extra(l) = max(dW_l(l) - air_depth(l)*dens_h2o, 0.)
           dW_l(l) = dW_l(l) - extra(l)
        enddo
        ! sweep any excess downward
        extra_cum = 0.
        do l = 1, num_l
           extra_cum = extra_cum + extra(l)
           if (air_depth(l)*dens_h2o.gt.dW_l(l)) then
              ddW = min(extra_cum, air_depth(l)*dens_h2o-dW_l(l))
              extra_cum = extra_cum - ddw
              dW_l(l) = dW_l(l) + ddW
           endif
        enddo
        ! sweep any remaining excess upward
        do l = num_l, 1, -1
           extra_cum = extra_cum + extra(l)
           if (air_depth(l)*dens_h2o.gt.dW_l(l)) then
              ddW = min(extra_cum, air_depth(l)*dens_h2o-dW_l(l))
              extra_cum = extra_cum - ddw
              dW_l(l) = dW_l(l) + ddW
           endif
        enddo
        ! place the water and heat
        do l= 1, num_l
           h1 = soil%heat_capacity_dry(l)*dz(l) &
                + csw*soil%ws(l) + clw*soil%wl(l)
           h2 = clw*dW_l(l)
           soil%T(l) = (h1 * soil%T(l) &
                               + h2 * tflow )  / (h1+h2)
           soil%wl(l) = soil%wl(l) + dW_l(l)
        enddo
     endif
     ! extra_cum should now be 0, but included below for safety
     lrunf_ie = lrunf_ie - (macro_inf-extra_cum)/delta_time
     hlrunf_ie = hlrunf_ie - clw*(tflow-tfreeze)*(macro_inf-extra_cum)/delta_time
  ENDIF

  flow_macro = (macro_inf-extra_cum)/delta_time

  hlrunf_bf = clw*sum(div_bf*(soil%T-tfreeze)) + sum(hdiv_it)
  ! div_bf is 0. for GW_TILED, else hdiv_it == zero
  hlrunf_if = clw*sum(div_if*(soil%T-tfreeze))
  hlrunf_al = clw*sum(div_al*(soil%T-tfreeze))
  hlrunf_sc = clw*lrunf_sc  *(soil%groundwater_T(1)-tfreeze)
  if (lrunf_from_div) then
     soil_lrunf  =  lrunf_sn +  lrunf_ie +  sum(div) +  lrunf_nu +  lrunf_sc
     if (gw_option /= GW_TILED) then
        soil_hlrunf = hlrunf_sn + hlrunf_ie +  clw*sum(div*(soil%T-tfreeze)) &
                                                      + hlrunf_nu + hlrunf_sc
     else
        soil_hlrunf = hlrunf_sn + hlrunf_ie +  sum(hdiv_it) &
             + hlrunf_nu + hlrunf_sc
     end if
  else
     soil_lrunf  =  lrunf_sn +  lrunf_ie +  lrunf_bf +  lrunf_if &
                             +  lrunf_al +  lrunf_nu +  lrunf_sc
     soil_hlrunf = hlrunf_sn + hlrunf_ie + hlrunf_bf + hlrunf_if &
                             + hlrunf_al + hlrunf_nu + hlrunf_sc
  endif
  soil_frunf = frunf
  soil_hfrunf = hfrunf

  if (is_watch_cell()) then
     write(*,*)'Soil runoff from point in watch_cell'
     __DEBUG3__(soil%hidx_j, soil_lrunf, lrunf_bf)
     __DEBUG5__(lrunf_sn, lrunf_ie, lrunf_if, lrunf_al, lrunf_sc)
  end if


  do l = 1, num_l
     ! ---- compute explicit melt/freeze --------------------------------------
     hcap = soil%heat_capacity_dry(l)*dz(l) &
              + clw*soil%wl(l) + csw*soil%ws(l)
     melt_per_deg = hcap/(hlf_factor*hlf)
     if       (soil%ws(l)>0 .and. soil%T(l)>soil%pars%tfreeze) then
       melt =  min(soil%ws(l), (soil%T(l)-soil%pars%tfreeze)*melt_per_deg)
     else if (soil%wl(l)>0 .and. soil%T(l)<soil%pars%tfreeze) then
       melt = -min(soil%wl(l), (soil%pars%tfreeze-soil%T(l))*melt_per_deg)
     else
       melt = 0
     endif
     soil%wl(l) = soil%wl(l) + melt
     soil%ws(l) = soil%ws(l) - melt
     soil%T(l) = tfreeze &
        + (hcap*(soil%T(l)-tfreeze) - hlf_factor*hlf*melt) &
                             / ( hcap + (clw-csw)*melt )
     soil_melt = soil_melt + melt / delta_time
  enddo

  if(is_watch_point()) then
     write(*,*) ' ##### soil_step_2 checkpoint 5 #####'
     do l = 1, num_l
        write(*,'(a,i2.2,100(2x,a,g23.16))') ' level=', l,&
             ' T =', soil%T(l),&
             ' Th=', (soil%ws(l) +soil%wl(l))/(dens_h2o*dz(l)),&
             ' wl=', soil%wl(l),&
             ' ws=', soil%ws(l),&
             ' gw=', soil%groundwater(l)
     enddo
     if (soil_carbon_option==SOILC_CORPSE) &
         call debug_pool(soil%leafLitter, 'leafLitter')
  endif

  active_layer_thickness = 0.
  do l = 1, num_l
     if (soil%ws(l).gt.0.) then
        active_layer_thickness = active_layer_thickness &
          + dz(l)*max(0.,soil%wl(l))/(soil%wl(l)+soil%ws(l))
        exit
     endif
     active_layer_thickness = active_layer_thickness + dz(l)
  enddo

  soil_Ttop = soil%T(1)
  soil_Ctop = soil%heat_capacity_dry(1)*dz(1) &
    + clw*soil%wl(1) + csw*soil%ws(1)

  soil%psi=psi+dPsi
!  if (do_component_balchecks) then
     ! Sum total water mass at end of soil_step_2
     call soil_tile_stock_pe (soil, sliq, sice )
     wsum2 = sliq + sice

     ! Subtract influxes and add outfluxes
     wsum2 = wsum2 - delta_time *  snow_lprec &
          + delta_time * ( subs_evap + vegn_uptk + soil_lrunf + soil_frunf)

     call check_conservation('soil_mod: soil_step_2', 'Water', wsum1, wsum2, wthresh, FATAL)
! endif


   if (is_watch_point()) then
      write(*,*)'##### soil_step_2 checkpoint 6 #####'
      __DEBUG1__(flow)
      __DEBUG1__(div)
      __DEBUG1__(wl_before)
      __DEBUG1__(gw_option)
      if (soil_carbon_option==SOILC_CORPSE) then
         call debug_pool(soil%leafLitter,       'leafLitter')
         call debug_pool(soil%fineWoodLitter,   'fineWoodLitter')
         call debug_pool(soil%coarseWoodLitter, 'coarseWoodLitter')
         do l = 1, num_l
            call debug_pool(soil%soil_C(l), 'soil_C(l)')
         enddo
      endif
      do l = 1, size(soil%div_hlsp_DOC,2)
         __DEBUG1__(soil%div_hlsp_DOC(:,l))
      enddo
   endif

!New version that combines the two leaching steps and should do a better job of moving DOC from litter layer
!For now, we are assuming that only leaf litter gets leached
!ZMS Edited to allow for tiled fluxes. Also pass in water content before Richards.
   if (gw_option == GW_TILED) then
      call carbon_leaching_with_litter(soil%soil_C(:),soil%leafLitter,soil%coarsewoodLitter,flow, &
            max(0.0,flow(1)),div,dz(1:num_l), wl_before, &
            delta_time,DOC_leached,leaflitter_DOC_loss,woodlitter_DOC_loss,div_DOC_loss, .TRUE., &
            soil%div_hlsp_DOC, surf_DOC_loss)
   else
      call carbon_leaching_with_litter(soil%soil_C(:),soil%leafLitter,soil%coarsewoodLitter,flow, &
            max(0.0,flow(1)),div,dz(1:num_l), wl_before, &
            delta_time,DOC_leached,leaflitter_DOC_loss,woodlitter_DOC_loss,div_DOC_loss, .FALSE.)
      surf_DOC_loss(:) = 0.
   end if

   soil%fast_DOC_leached=soil%fast_DOC_leached+sum(div_DOC_loss(1,:)) + surf_DOC_loss(1)
   soil%slow_DOC_leached=soil%slow_DOC_leached+sum(div_DOC_loss(2,:)) + surf_DOC_loss(2)
   soil%deadmic_DOC_leached=soil%deadmic_DOC_leached+sum(div_DOC_loss(3,:)) + surf_DOC_loss(3)
   ! Diagnostic. Later pass this back to land_model for transfer to rivers.
   total_DOC_div = sum(surf_DOC_loss(:))
   do l=1,num_l
      total_DOC_div = total_DOC_div + sum(div_DOC_loss(:,l))
   end do
   total_DOC_div = total_DOC_div/delta_time
   if (i_river_DOC/=NO_TRACER) &
       soil_tr_runf(i_river_DOC) = total_DOC_div

   if (is_watch_point()) then
      if (soil_carbon_option==SOILC_CORPSE) then
         write(*,*)'##### soil_step_2 checkpoint 7 #####'
         call debug_pool(soil%leafLitter,       'leafLitter')
         call debug_pool(soil%fineWoodLitter,   'fineWoodLitter')
         call debug_pool(soil%coarseWoodLitter, 'coarseWoodLitter')
         __DEBUG3__(leaflitter_DOC_loss,woodlitter_DOC_loss,total_DOC_div)
         do l = 1, num_l
            call debug_pool(soil%soil_C(l), 'soil_C(l)')
         enddo
      endif
   endif

! ----------------------------------------------------------------------------
! given solution for surface energy balance, write diagnostic output.
!

  ! ---- diagnostic section
  call send_tile_data(id_temp, soil%T, diag)
  if (id_lwc > 0) call send_tile_data(id_lwc,  soil%wl/dz(1:num_l), diag)
  if (id_swc > 0) call send_tile_data(id_swc,  soil%ws/dz(1:num_l), diag)
  if (id_psi > 0) call send_tile_data(id_psi,  psi+dPsi, diag)

  ! CMOR variables
  if (id_mrlsl > 0) call send_tile_data(id_mrlsl, soil%wl+soil%ws, diag)
  if (id_mrsfl > 0) call send_tile_data(id_mrsfl, soil%ws, diag)
  if (id_mrsll > 0) call send_tile_data(id_mrsll, soil%wl, diag)
  if (id_mrsol > 0) call send_tile_data(id_mrsol, soil%wl+soil%ws, diag)
  if (id_mrso > 0)  call send_tile_data(id_mrso,  sum(soil%wl+soil%ws), diag)
  if (id_mrsos > 0) call send_tile_data(id_mrsos, sum((soil%wl+soil%ws)*mrsos_weight), diag)
  if (id_mrs1mLut > 0) call send_tile_data(id_mrs1mLut, sum((soil%wl+soil%ws)*mrs1m_weight), diag)
  if (id_mrfso > 0) call send_tile_data(id_mrfso, sum(soil%ws), diag)
  if (id_mrlso > 0) call send_tile_data(id_mrlso, sum(soil%wl), diag)
  call send_tile_data(id_mrros, lrunf_ie+lrunf_sn, diag)
  call send_tile_data(id_mrro,  lrunf_ie+lrunf_sn+lrunf_bf+lrunf_nu, diag)

  ! ZMS uncomment for back-compatibility with diag tables
  if (gw_option == GW_TILED) then
     call send_tile_data(id_deficit, deficit, diag)
     call send_tile_data(id_sat_depth, depth_to_wt_3, diag)
     call send_tile_data(id_sat_dept2, depth_to_wt2_3, diag)
     call send_tile_data(id_z_cap, depth_to_cf_3, diag)
     if (depth_to_wt_2a .ge. -0.5) &
          call send_tile_data(id_sat_depth, depth_to_wt_2a, diag)

  end if
  !f1p: save sat_area_frac for use in tracer deposition calculations
  soil%sat_area_frac = sat_area_frac

  call send_tile_data(id_cf_1, depth_to_cf_1, diag)
  call send_tile_data(id_cf_3, depth_to_cf_3, diag)
  call send_tile_data(id_wt_1, depth_to_wt_1, diag)
  call send_tile_data(id_wt_2, depth_to_wt_2, diag)
  call send_tile_data(id_wt_2a, depth_to_wt_2a, diag)
  call send_tile_data(id_wt_2b, depth_to_wt_2b, diag)
  call send_tile_data(id_wt_3, depth_to_wt_3, diag)
  call send_tile_data(id_wt2_3, depth_to_wt2_3, diag)
  call send_tile_data(id_wt_4, depth_to_wt_4, diag)
  call send_tile_data(id_deficit_2, deficit_2, diag)
  call send_tile_data(id_deficit_3, deficit_3, diag)
  call send_tile_data(id_deficit_4, deficit_4, diag)
  call send_tile_data(id_sat_frac, sat_area_frac, diag)
  call send_tile_data(id_div_bf, div_bf, diag)
  call send_tile_data(id_div_if, div_if, diag)
  call send_tile_data(id_div_al, div_al, diag)

  call send_tile_data(id_ie,   lrunf_ie, diag)
  call send_tile_data(id_sn,   lrunf_sn, diag)
  call send_tile_data(id_bf,   lrunf_bf, diag)
  call send_tile_data(id_if,   lrunf_if, diag)
  call send_tile_data(id_al,   lrunf_al, diag)
  call send_tile_data(id_nu,   lrunf_nu, diag)
  call send_tile_data(id_sc,   lrunf_sc, diag)
  call send_tile_data(id_hie,  hlrunf_ie, diag)
  call send_tile_data(id_hsn,  hlrunf_sn, diag)
  call send_tile_data(id_hbf,  hlrunf_bf, diag)
  call send_tile_data(id_hif,  hlrunf_if, diag)
  call send_tile_data(id_hal,  hlrunf_al, diag)
  call send_tile_data(id_hnu,  hlrunf_nu, diag)
  call send_tile_data(id_hsc,  hlrunf_sc, diag)
  if (id_evap > 0) call send_tile_data(id_evap,  soil_levap+soil_fevap, diag)

   if (id_leaflitter_fast_C_leaching > 0) call send_tile_data(id_leaflitter_fast_C_leaching,leaflitter_DOC_loss(1)/delta_time,diag)
   if (id_leaflitter_slow_C_leaching > 0) call send_tile_data(id_leaflitter_slow_C_leaching,leaflitter_DOC_loss(2)/delta_time,diag)
   if (id_leaflitter_deadmic_C_leaching > 0) call send_tile_data(id_leaflitter_deadmic_C_leaching,leaflitter_DOC_loss(3)/delta_time,diag)

   if (id_coarsewoodlitter_fast_C_leaching > 0) call send_tile_data(id_coarsewoodlitter_fast_C_leaching,woodlitter_DOC_loss(1)/delta_time,diag)
   if (id_coarsewoodlitter_slow_C_leaching > 0) call send_tile_data(id_coarsewoodlitter_slow_C_leaching,woodlitter_DOC_loss(2)/delta_time,diag)
   if (id_coarsewoodlitter_deadmic_C_leaching > 0) call send_tile_data(id_coarsewoodlitter_deadmic_C_leaching,woodlitter_DOC_loss(3)/delta_time,diag)


  call send_tile_data(id_heat_cap, soil%heat_capacity_dry, diag)
  call send_tile_data(id_active_layer, active_layer_thickness, diag)
  if (gw_option == GW_TILED) then
     call send_tile_data(id_surface_water, surface_water, diag)
     call send_tile_data(id_inun_frac, inundated_frac, diag)
     call send_tile_data(id_wet_frac, wet_frac, diag)
     if (simple_inundation) call send_tile_data(id_rsn_frac, sat_runf_frac, diag)
  end if
  if (id_flow > 0) then
     flow_s(:) = flow(1:num_l) / delta_time
     call send_tile_data(id_flow, flow_s, diag)
  end if
  if (id_reflux > 0) then
     reflux = max(lrunf_ie - lprec_eff, 0.)
     call send_tile_data(id_reflux, reflux, diag)
  end if
  call send_tile_data(id_macro_infilt, flow_macro, diag)
  call send_tile_data(id_fast_C_leaching, DOC_leached(1,:)/delta_time,diag)
  call send_tile_data(id_slow_C_leaching, DOC_leached(2,:)/delta_time,diag)
  call send_tile_data(id_deadmic_C_leaching, DOC_leached(3,:)/delta_time,diag)
  do l=1,num_l
     total_C_leaching(l) = sum(DOC_leached(:,l))/delta_time
  end do
  call send_tile_data(id_total_C_leaching, total_C_leaching, diag)
  if (gw_option == GW_TILED) then
     call send_tile_data(id_surf_DOC_loss, sum(surf_DOC_loss(:))/delta_time,diag)
  end if
  call send_tile_data(id_total_DOC_div_loss, total_DOC_div, diag)

  if (.not. LM2) call send_tile_data(id_psi_bot, soil%psi(num_l), diag)

end subroutine soil_step_2

! ============================================================================
subroutine soil_step_3(soil, diag)
  type(soil_tile_type), intent(in) :: soil
  type(diag_buff_type), intent(inout) :: diag

  real :: sum_fsc, sum_ssc, sum_deadmic, sum_livemic, sum_protectedC !, slomtot
  real :: fast_C(num_l), slow_C(num_l), deadMicrobeC(num_l), liveMicrobeC(num_l), protectedC(num_l)
  real :: litter_fast_C, litter_slow_C, litter_deadmic, litter_livemic
  integer :: layer, ncohorts(num_l), litter_ncohorts
  real, dimension(num_l) :: fast_dissolved,slow_dissolved,deadmic_dissolved
  real :: total_fast, total_slow, total_deadmic, total_livemic, total_protected, total_dissolved, total_carbon
  real :: total_carbon_layered(num_l)

  select case (soil_carbon_option)
  case (SOILC_CENTURY,SOILC_CENTURY_BY_LAYER)
     if (id_fsc>0) call send_tile_data(id_fsc, sum(soil%fast_soil_C(:)), diag)
     if (id_ssc>0) call send_tile_data(id_ssc, sum(soil%slow_soil_C(:)), diag)
     if (id_csoil>0) call send_tile_data(id_csoil, sum(soil%fast_soil_C(:))+sum(soil%slow_soil_C(:)), diag)
     call send_tile_data(id_fast_soil_C, soil%fast_soil_C(:)/dz(1:num_l), diag)
     call send_tile_data(id_slow_soil_C, soil%slow_soil_C(:)/dz(1:num_l), diag)

     ! --- CMOR vars
     if (id_csoilfast>0) call send_tile_data(id_csoilfast, sum(soil%fast_soil_C(:)), diag)
     if (id_csoilmedium>0) call send_tile_data(id_csoilmedium, sum(soil%slow_soil_C(:)), diag)
     call send_tile_data(id_csoilslow, 0.0, diag)
     ! --- end of CMOR vars

  case (SOILC_CORPSE)
     total_carbon_layered=0.0
     total_fast=0.0
     total_slow=0.0
     total_deadmic=0.0
     total_livemic=0.0
     total_protected=0.0
     total_dissolved=0.0
     total_carbon=0.0

     DO layer=1,num_l
       call poolTotalCarbon(soil%soil_C(layer),fast_C(layer),slow_C(layer),&
       deadMicrobeC(layer),liveMicrobeC(layer),protectedC(layer),&
       fast_dissolved(layer),slow_dissolved(layer),deadmic_dissolved(layer),ncohorts(layer),total_carbon_layered(layer))
     ENDDO

     total_fast=sum(fast_C)
     total_slow=sum(slow_C)
     total_deadmic=sum(deadMicrobeC)
     total_livemic=sum(liveMicrobeC)
     total_protected=sum(protectedC)
     total_dissolved=sum(fast_dissolved+slow_dissolved+deadmic_dissolved)

     if (id_fast_soil_C > 0) call send_tile_data(id_fast_soil_C, fast_C/dz(1:num_l), diag)
     if (id_slow_soil_C > 0) call send_tile_data(id_slow_soil_C, slow_C/dz(1:num_l), diag)
     if (id_deadmic > 0) call send_tile_data(id_deadmic, deadMicrobeC/dz(1:num_l), diag)
     if (id_livemic > 0) call send_tile_data(id_livemic, liveMicrobeC/dz(1:num_l), diag)
     if (id_protectedC > 0) call send_tile_data(id_protectedC, protectedC/dz(1:num_l), diag)
     if (id_nsoilcohorts > 0) call send_tile_data(id_nsoilcohorts, real(ncohorts), diag)
     if (id_fast_dissolved_C > 0) call send_tile_data(id_fast_dissolved_C, fast_dissolved/dz(1:num_l), diag)
     if (id_slow_dissolved_C > 0) call send_tile_data(id_slow_dissolved_C, slow_dissolved/dz(1:num_l), diag)
     if (id_deadmic_dissolved_C > 0) call send_tile_data(id_deadmic_dissolved_C, deadmic_dissolved/dz(1:num_l), diag)
     if (id_total_carbon_layered > 0) call send_tile_data(id_total_carbon_layered, total_carbon_layered/dz(1:num_l),diag)

     if (id_fast_DOC_div_loss > 0) call send_tile_data(id_fast_DOC_div_loss, soil%fast_DOC_leached,diag)
     if (id_slow_DOC_div_loss > 0) call send_tile_data(id_slow_DOC_div_loss, soil%slow_DOC_leached,diag)
     if (id_deadmic_DOC_div_loss > 0) call send_tile_data(id_deadmic_DOC_div_loss, soil%deadmic_DOC_leached,diag)

     call poolTotalCarbon(soil%leafLitter,litter_fast_C,litter_slow_C,litter_deadmic,litter_livemic,ncohorts=litter_ncohorts)
     total_fast=total_fast+litter_fast_C
     total_slow=total_slow+litter_slow_C
     total_deadmic=total_deadmic+litter_deadmic
     total_livemic=total_livemic+litter_livemic
     total_dissolved=total_dissolved+sum(soil%leafLitter%dissolved_carbon(:))


     if (id_nleaflittercohorts > 0) call send_tile_data(id_nleaflittercohorts, real(litter_ncohorts), diag)
     if (id_leaflitter_fast_C > 0) call send_tile_data(id_leaflitter_fast_C, litter_fast_C, diag)
     if (id_leaflitter_slow_C > 0) call send_tile_data(id_leaflitter_slow_C, litter_slow_C, diag)
     if (id_leaflitter_deadmic > 0) call send_tile_data(id_leaflitter_deadmic, litter_deadmic, diag)
     if (id_leaflitter_livemic > 0) call send_tile_data(id_leaflitter_livemic, litter_livemic, diag)
     if (id_leaflitter_fast_dissolved_C > 0) call send_tile_data(id_leaflitter_fast_dissolved_C, soil%leafLitter%dissolved_carbon(1), diag)
     if (id_leaflitter_slow_dissolved_C > 0) call send_tile_data(id_leaflitter_slow_dissolved_C, soil%leafLitter%dissolved_carbon(2), diag)
     if (id_leaflitter_deadmic_dissolved_C > 0) call send_tile_data(id_leaflitter_deadmic_dissolved_C, soil%leafLitter%dissolved_carbon(3), diag)
     if (id_leaflitter_total_C > 0) call send_tile_data(id_leaflitter_total_C, &
           litter_fast_C+litter_slow_C+litter_deadmic+litter_livemic+sum(soil%leafLitter%dissolved_carbon(:)),diag)

     call poolTotalCarbon(soil%fineWoodLitter,litter_fast_C,litter_slow_C,litter_deadmic,litter_livemic,ncohorts=litter_ncohorts)
     total_fast=total_fast+litter_fast_C
     total_slow=total_slow+litter_slow_C
     total_deadmic=total_deadmic+litter_deadmic
     total_livemic=total_livemic+litter_livemic
     total_dissolved=total_dissolved+sum(soil%fineWoodLitter%dissolved_carbon(:))

     if (id_nfineWoodlittercohorts > 0) call send_tile_data(id_nfineWoodlittercohorts, real(litter_ncohorts), diag)
     if (id_fineWoodlitter_fast_C > 0) call send_tile_data(id_fineWoodlitter_fast_C, litter_fast_C, diag)
     if (id_fineWoodlitter_slow_C > 0) call send_tile_data(id_fineWoodlitter_slow_C, litter_slow_C, diag)
     if (id_fineWoodlitter_deadmic > 0) call send_tile_data(id_fineWoodlitter_deadmic, litter_deadmic, diag)
     if (id_fineWoodlitter_livemic > 0) call send_tile_data(id_fineWoodlitter_livemic, litter_livemic, diag)
     if (id_fineWoodlitter_fast_dissolved_C > 0) call send_tile_data(id_fineWoodlitter_fast_dissolved_C, soil%fineWoodLitter%dissolved_carbon(1), diag)
     if (id_fineWoodlitter_slow_dissolved_C > 0) call send_tile_data(id_fineWoodlitter_slow_dissolved_C, soil%fineWoodLitter%dissolved_carbon(2), diag)
     if (id_fineWoodlitter_deadmic_dissolved_C > 0) call send_tile_data(id_fineWoodlitter_deadmic_dissolved_C, soil%fineWoodLitter%dissolved_carbon(3), diag)
     if (id_fineWoodlitter_total_C > 0) call send_tile_data(id_fineWoodlitter_total_C, &
           litter_fast_C+litter_slow_C+litter_deadmic+litter_livemic+sum(soil%fineWoodLitter%dissolved_carbon(:)),diag)


     call poolTotalCarbon(soil%coarseWoodLitter,litter_fast_C,litter_slow_C,litter_deadmic,litter_livemic,ncohorts=litter_ncohorts)
     total_fast=total_fast+litter_fast_C
     total_slow=total_slow+litter_slow_C
     total_deadmic=total_deadmic+litter_deadmic
     total_livemic=total_livemic+litter_livemic
     total_dissolved=total_dissolved+sum(soil%leafLitter%dissolved_carbon(:))

     if (id_ncoarseWoodlittercohorts > 0) call send_tile_data(id_ncoarseWoodlittercohorts, real(litter_ncohorts), diag)
     if (id_coarseWoodlitter_fast_C > 0) call send_tile_data(id_coarseWoodlitter_fast_C, litter_fast_C, diag)
     if (id_coarseWoodlitter_slow_C > 0) call send_tile_data(id_coarseWoodlitter_slow_C, litter_slow_C, diag)
     if (id_coarseWoodlitter_deadmic > 0) call send_tile_data(id_coarseWoodlitter_deadmic, litter_deadmic, diag)
     if (id_coarseWoodlitter_livemic > 0) call send_tile_data(id_coarseWoodlitter_livemic, litter_livemic, diag)
     if (id_coarseWoodlitter_fast_dissolved_C > 0) call send_tile_data(id_coarseWoodlitter_fast_dissolved_C, soil%coarseWoodLitter%dissolved_carbon(1), diag)
     if (id_coarseWoodlitter_slow_dissolved_C > 0) call send_tile_data(id_coarseWoodlitter_slow_dissolved_C, soil%coarseWoodLitter%dissolved_carbon(2), diag)
     if (id_coarseWoodlitter_deadmic_dissolved_C > 0) call send_tile_data(id_coarseWoodlitter_deadmic_dissolved_C, soil%coarseWoodLitter%dissolved_carbon(3), diag)
     if (id_coarseWoodlitter_total_C > 0) call send_tile_data(id_coarseWoodlitter_total_C, &
           litter_fast_C+litter_slow_C+litter_deadmic+litter_livemic+sum(soil%coarseWoodLitter%dissolved_carbon(:)),diag)

     sum_fsc = total_fast
     sum_ssc = total_slow
     sum_deadmic = total_deadmic
     sum_livemic = total_livemic
     sum_protectedC = total_protected
     total_carbon=total_fast+total_slow+total_deadmic+total_livemic+total_dissolved+total_protected
     call send_tile_data(id_fsc, sum_fsc, diag)
     call send_tile_data(id_ssc, sum_ssc, diag)
     ! --- CMOR vars
     call send_tile_data(id_csoilfast,   sum_fsc, diag)
     call send_tile_data(id_csoilmedium, sum_ssc, diag)
     call send_tile_data(id_csoilslow,   0.0,     diag)
     ! --- end of CMOR vars
     if (id_csoil > 0) call send_tile_data(id_csoil, sum_ssc+sum_fsc, diag)
     if (id_deadmic_total > 0) call send_tile_data(id_deadmic_total, sum_deadmic, diag)
     if (id_livemic_total > 0) call send_tile_data(id_livemic_total, sum_livemic, diag)
     if (id_slomtot > 0) then
   !     slomtot = sum_fsc + sum_ssc + sum_deadmic + sum_livemic + sum_protectedC + &
   !               litter_fast_C + litter_slow_C + litter_deadmic + litter_livemic
        call send_tile_data(id_slomtot, total_carbon, diag)
     end if

     if (id_protected_total > 0) call send_tile_data(id_protected_total, sum_protectedC, diag)
     if (id_dissolved_total > 0) call send_tile_data(id_dissolved_total, total_dissolved, diag)
     if (id_total_soil_C > 0) call send_tile_data(id_total_soil_C, total_carbon, diag)
  end select

end subroutine soil_step_3


! ============================================================================
subroutine Dsdt(vegn, soil, diag, soilt, theta)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  type(diag_buff_type), intent(inout) :: diag
  real                , intent(in)    :: soilt ! average soil temperature, deg K
  real                , intent(in)    :: theta ! average soil moisture

  select case (soil_carbon_option)
  case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
     call Dsdt_CENTURY(vegn, soil, diag, soilt, theta)
  case (SOILC_CORPSE)
     call Dsdt_CORPSE(vegn, soil, diag)
  case default
     call error_mesg('Dsdt','soil_carbon_option is invalid. This should never happen. Contact developer', FATAL)
  end select
end subroutine Dsdt


! ============================================================================
subroutine Dsdt_CORPSE(vegn, soil, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  type(diag_buff_type), intent(inout) :: diag

  real :: leaflitter_fast_C_loss_rate, leaflitter_slow_C_loss_rate, leaflitter_deadmic_C_loss_rate
  real :: finewoodlitter_fast_C_loss_rate, finewoodlitter_slow_C_loss_rate, finewoodlitter_deadmic_C_loss_rate
  real :: coarsewoodlitter_fast_C_loss_rate, coarsewoodlitter_slow_C_loss_rate, coarsewoodlitter_deadmic_C_loss_rate
  real :: fast_C_loss_rate(size(soil%soil_C))
  real :: slow_C_loss_rate(size(soil%soil_C))
  real :: dead_microbe_C_loss_rate(size(soil%soil_C))
  real, dimension(size(soil%soil_C)) :: decomp_T,decomp_theta,ice_porosity
  real :: A          (size(soil%soil_C)) ! decomp rate reduction due to moisture and temperature

  integer :: badCohort   ! For soil carbon pool carbon balance and invalid number check
  integer :: k
  real :: CO2prod,protected_produced(3,size(soil%soil_C)),protected_turnover_rate(3,size(soil%soil_C))
  real :: leaflitter_protected_produced(3),leaflitter_protected_turnover_rate(3)
  real :: finewoodlitter_protected_produced(3),finewoodlitter_protected_turnover_rate(3)
  real :: coarsewoodlitter_protected_produced(3),coarsewoodlitter_protected_turnover_rate(3)
  real :: leaflitter_C_dissolved(3),leaflitter_C_deposited(3),C_dissolved(3,num_l),C_deposited(3,num_l)
  real :: finewoodlitter_C_dissolved(3),finewoodlitter_C_deposited(3)
  real :: coarsewoodlitter_C_dissolved(3),coarsewoodlitter_C_deposited(3)
  real :: deadmic_produced(size(soil%soil_C)), leaflitter_deadmic_produced, finewoodlitter_deadmic_produced, coarsewoodlitter_deadmic_produced
  real :: total_fast_C, total_slow_C, total_deadmic_C,total_livemic_C,temp_fast,temp_slow,temp_livemic,temp_deadmic,temp_protected
  real :: temp_protected_fast, temp_protected_slow, temp_protected_deadmic
  integer :: point_i,point_j,point_k,point_face

  A(:) = A_function(soil%T(:), soil_theta(soil))
  decomp_T = soil%T(:)
  decomp_theta = soil_theta(soil)
  ice_porosity = soil_ice_porosity(soil)

  vegn%rh=0.0
  total_fast_C=0.0
  total_slow_C=0.0
  total_deadmic_C=0.0
  total_livemic_C=0.0

  !  First surface litter is decomposed

  call update_pool(soil%leafLitter,decomp_T(1),decomp_theta(1),1.0-(decomp_theta(1)+ice_porosity(1)),&
            soil%wl(1),soil%ws(1),dt_fast_yr,dz(1),&
            leaflitter_fast_C_loss_rate,leaflitter_slow_C_loss_rate, leaflitter_deadmic_C_loss_rate, CO2prod, &
            leaflitter_deadmic_produced, leaflitter_protected_produced, leaflitter_protected_turnover_rate, leaflitter_C_dissolved, leaflitter_C_deposited, badCohort)
  IF (badCohort.ne.0) THEN
        call get_current_point(point_i,point_j,point_k,point_face)
        WRITE (*,*), 'Found bad cohort in leaf litter.  Point i,j,k,face:',point_i,point_j,point_k,point_face
        WRITE (*,*), 'T=',decomp_T(1),'theta=',decomp_theta(1),'dt=',dt_fast_yr
        call error_mesg('Dsdt','Found bad cohort in leaf litter',FATAL)
  ENDIF

  call update_pool(soil%fineWoodLitter,decomp_T(1),decomp_theta(1),1.0-(decomp_theta(1)+ice_porosity(1)),&
            soil%wl(1),soil%ws(1),dt_fast_yr,dz(1),&
            fineWoodlitter_fast_C_loss_rate,fineWoodlitter_slow_C_loss_rate, fineWoodlitter_deadmic_C_loss_rate, CO2prod, &
            fineWoodlitter_deadmic_produced, fineWoodlitter_protected_produced, fineWoodlitter_protected_turnover_rate, fineWoodlitter_C_dissolved, fineWoodlitter_C_deposited, badCohort)
  IF (badCohort.ne.0) THEN
        call get_current_point(point_i,point_j,point_k,point_face)
        WRITE (*,*), 'Found bad cohort in fineWood litter.  Point i,j,k,face:',point_i,point_j,point_k,point_face
        WRITE (*,*), 'T=',decomp_T(1),'theta=',decomp_theta(1),'dt=',dt_fast_yr
        call error_mesg('Dsdt','Found bad cohort in fineWood litter',FATAL)
  ENDIF

  call update_pool(soil%coarseWoodLitter,decomp_T(1),decomp_theta(1),1.0-(decomp_theta(1)+ice_porosity(1)),&
            soil%wl(1),soil%ws(1),dt_fast_yr,dz(1),&
            coarseWoodlitter_fast_C_loss_rate,coarseWoodlitter_slow_C_loss_rate, coarseWoodlitter_deadmic_C_loss_rate, CO2prod, &
            coarseWoodlitter_deadmic_produced, coarseWoodlitter_protected_produced, coarseWoodlitter_protected_turnover_rate, coarseWoodlitter_C_dissolved, coarseWoodlitter_C_deposited, badCohort)
  IF (badCohort.ne.0) THEN
        call get_current_point(point_i,point_j,point_k,point_face)
        WRITE (*,*), 'Found bad cohort in coarseWood litter.  Point i,j,k,face:',point_i,point_j,point_k,point_face
        WRITE (*,*), 'T=',decomp_T(1),'theta=',decomp_theta(1),'dt=',dt_fast_yr
        call error_mesg('Dsdt','Found bad cohort in coarseWood litter',FATAL)
  ENDIF

  ! loss of C to atmosphere
  vegn%rh=vegn%rh + CO2prod/dt_fast_yr

  call poolTotalCarbon(soil%leafLitter,fastC=temp_fast,slowC=temp_slow,deadMicrobeC=temp_deadmic,liveMicrobeC=temp_livemic)
  total_fast_C=total_fast_C+temp_fast
  total_slow_C=total_slow_C+temp_slow
  total_deadmic_C=total_deadmic_C+temp_deadmic
  total_livemic_C=total_livemic_C+temp_livemic

  !Accumulate turnover rates for determining steady state pools
  if(temp_fast>0)soil%leaflitter_fast_turnover_accumulated=soil%leaflitter_fast_turnover_accumulated+leaflitter_fast_C_loss_rate/temp_fast
  if(temp_slow>0)soil%leaflitter_slow_turnover_accumulated=soil%leaflitter_slow_turnover_accumulated+leaflitter_slow_C_loss_rate/temp_slow
  if(temp_deadmic>0)soil%leaflitter_deadmic_turnover_accumulated=soil%leaflitter_deadmic_turnover_accumulated+leaflitter_deadmic_C_loss_rate/temp_deadmic
  soil%leaflitter_deadmic_in=soil%leaflitter_deadmic_in+leaflitter_deadmic_produced
  soil%leaflitter_fsc_in=soil%leaflitter_fsc_in+leaflitter_C_deposited(1)-leaflitter_C_dissolved(1)
  soil%leaflitter_ssc_in=soil%leaflitter_ssc_in+leaflitter_C_deposited(2)-leaflitter_C_dissolved(2)
  soil%leaflitter_deadmic_in=soil%leaflitter_deadmic_in+leaflitter_C_deposited(3)-leaflitter_C_dissolved(3)

  call poolTotalCarbon(soil%finewoodLitter,fastC=temp_fast,slowC=temp_slow,deadMicrobeC=temp_deadmic,liveMicrobeC=temp_livemic)
  total_fast_C=total_fast_C+temp_fast
  total_slow_C=total_slow_C+temp_slow
  total_deadmic_C=total_deadmic_C+temp_deadmic
  total_livemic_C=total_livemic_C+temp_livemic

  !Accumulate turnover rates for determining steady state pools
  if(temp_fast>0)soil%finewoodlitter_fast_turnover_accumulated=soil%finewoodlitter_fast_turnover_accumulated+finewoodlitter_fast_C_loss_rate/temp_fast
  if(temp_slow>0)soil%finewoodlitter_slow_turnover_accumulated=soil%finewoodlitter_slow_turnover_accumulated+finewoodlitter_slow_C_loss_rate/temp_slow
  if(temp_deadmic>0)soil%finewoodlitter_deadmic_turnover_accumulated=soil%finewoodlitter_deadmic_turnover_accumulated+finewoodlitter_deadmic_C_loss_rate/temp_deadmic
  soil%finewoodlitter_deadmic_in=soil%finewoodlitter_deadmic_in+finewoodlitter_deadmic_produced
  soil%finewoodlitter_fsc_in=soil%finewoodlitter_fsc_in+finewoodlitter_C_deposited(1)-finewoodlitter_C_dissolved(1)
  soil%finewoodlitter_ssc_in=soil%finewoodlitter_ssc_in+finewoodlitter_C_deposited(2)-finewoodlitter_C_dissolved(2)
  soil%finewoodlitter_deadmic_in=soil%finewoodlitter_deadmic_in+finewoodlitter_C_deposited(3)-finewoodlitter_C_dissolved(3)


  call poolTotalCarbon(soil%coarsewoodLitter,fastC=temp_fast,slowC=temp_slow,deadMicrobeC=temp_deadmic,liveMicrobeC=temp_livemic)
  total_fast_C=total_fast_C+temp_fast
  total_slow_C=total_slow_C+temp_slow
  total_deadmic_C=total_deadmic_C+temp_deadmic
  total_livemic_C=total_livemic_C+temp_livemic

  !Accumulate turnover rates for determining steady state pools
  if(temp_fast>0)soil%coarsewoodlitter_fast_turnover_accumulated=soil%coarsewoodlitter_fast_turnover_accumulated+coarsewoodlitter_fast_C_loss_rate/temp_fast
  if(temp_slow>0)soil%coarsewoodlitter_slow_turnover_accumulated=soil%coarsewoodlitter_slow_turnover_accumulated+coarsewoodlitter_slow_C_loss_rate/temp_slow
  if(temp_deadmic>0)soil%coarsewoodlitter_deadmic_turnover_accumulated=soil%coarsewoodlitter_deadmic_turnover_accumulated+coarsewoodlitter_deadmic_C_loss_rate/temp_deadmic
  soil%coarsewoodlitter_deadmic_in=soil%coarsewoodlitter_deadmic_in+coarsewoodlitter_deadmic_produced
  soil%coarsewoodlitter_fsc_in=soil%coarsewoodlitter_fsc_in+coarsewoodlitter_C_deposited(1)-coarsewoodlitter_C_dissolved(1)
  soil%coarsewoodlitter_ssc_in=soil%coarsewoodlitter_ssc_in+coarsewoodlitter_C_deposited(2)-coarsewoodlitter_C_dissolved(2)
  soil%coarsewoodlitter_deadmic_in=soil%coarsewoodlitter_deadmic_in+coarsewoodlitter_C_deposited(3)-coarsewoodlitter_C_dissolved(3)


  ! Next we have to go through layers and decompose the soil carbon pools
  do k=1,size(soil%soil_C)
    call update_pool(soil%soil_C(k),decomp_T(k),decomp_theta(k),1.0-(decomp_theta(k)+ice_porosity(k)),&
        soil%wl(k),soil%ws(k),dt_fast_yr,dz(k),&
    fast_C_loss_rate(k), slow_C_loss_rate(k), dead_microbe_C_loss_rate(k),CO2prod,&
    deadmic_produced(k),protected_produced(:,k),protected_turnover_rate(:,k),C_dissolved(:,k),C_deposited(:,k),badCohort)
    IF (badCohort.ne.0) THEN
        call get_current_point(point_i,point_j,point_k,point_face)
        WRITE (*,*), 'Found bad cohort in layer',k,'Point i,j,k,face:',point_i,point_j,point_k,point_face
        WRITE (*,*), 'T=',decomp_T(k),'theta=',decomp_theta(k),'dt=',dt_fast_yr
        call error_mesg('Dsdt','Found bad cohort',FATAL)
    ENDIF

    vegn%rh=vegn%rh + CO2prod/dt_fast_yr
    call poolTotalCarbon(soil%soil_C(k),fastC=temp_fast,slowC=temp_slow,deadMicrobeC=temp_deadmic,&
            liveMicrobeC=temp_livemic,protectedC=temp_protected,&
            fast_protectedC=temp_protected_fast,slow_protectedC=temp_protected_slow,deadmic_protectedC=temp_protected_deadmic)
    total_fast_C=total_fast_C+temp_fast
    total_slow_C=total_slow_C+temp_slow
    total_deadmic_C=total_deadmic_C+temp_deadmic
    total_livemic_C=total_livemic_C+temp_livemic

    !Accumulate turnover rates for determining steady state pools
    if(temp_fast>0)soil%fast_turnover_accumulated(k)=soil%fast_turnover_accumulated(k)+fast_C_loss_rate(k)/temp_fast
    if(temp_slow>0)soil%slow_turnover_accumulated(k)=soil%slow_turnover_accumulated(k)+slow_C_loss_rate(k)/temp_slow
    if(temp_deadmic>0)soil%deadmic_turnover_accumulated(k)=soil%deadmic_turnover_accumulated(k)+dead_microbe_C_loss_rate(k)/temp_deadmic
    soil%fast_protected_in(k)=soil%fast_protected_in(k)+protected_produced(1,k)
    soil%slow_protected_in(k)=soil%slow_protected_in(k)+protected_produced(2,k)
    soil%deadmic_protected_in(k)=soil%deadmic_protected_in(k)+protected_produced(3,k)
    if(temp_protected_fast>0) soil%fast_protected_turnover_accumulated(k)=soil%fast_protected_turnover_accumulated(k)+protected_turnover_rate(1,k)/temp_protected_fast
    if(temp_protected_slow>0) soil%slow_protected_turnover_accumulated(k)=soil%slow_protected_turnover_accumulated(k)+protected_turnover_rate(2,k)/temp_protected_slow
    if(temp_protected_deadmic>0) soil%deadmic_protected_turnover_accumulated(k)=soil%deadmic_protected_turnover_accumulated(k)+protected_turnover_rate(3,k)/temp_protected_deadmic
    soil%deadmic_in(k)=soil%deadmic_in(k)+deadmic_produced(k)
    soil%fsc_in(k)=soil%fsc_in(k)+C_deposited(1,k)-C_dissolved(1,k)
    soil%ssc_in(k)=soil%ssc_in(k)+C_deposited(2,k)-C_dissolved(2,k)
    soil%deadmic_in(k)=soil%deadmic_in(k)+C_deposited(3,k)-C_dissolved(3,k)
  enddo


  ! for budget check
  vegn%fsc_out = vegn%fsc_out + (sum(fast_C_loss_rate(:)) + leaflitter_fast_C_loss_rate + finewoodlitter_fast_C_loss_rate + coarsewoodlitter_fast_C_loss_rate)*dt_fast_yr
  vegn%ssc_out = vegn%ssc_out + (sum(slow_C_loss_rate(:)) + leaflitter_slow_C_loss_rate + finewoodlitter_slow_C_loss_rate + coarsewoodlitter_slow_C_loss_rate)*dt_fast_yr;
  vegn%deadmic_out = vegn%deadmic_out + (sum(dead_microbe_C_loss_rate(:)) + leaflitter_deadmic_C_loss_rate + coarsewoodlitter_deadmic_C_loss_rate + finewoodlitter_deadmic_C_loss_rate)*dt_fast_yr


  ! accumulate decomposition rate reduction for the soil carbon restart output
  soil%asoil_in(:) = soil%asoil_in(:) + A(:)



  ! TODO: arithmetic averaging of A doesn't seem correct; we need to invent something better,
  !       e.g. weight it with the carbon loss, or something like that

  ! ---- diagnostic section
  if (id_rsoil_fast>0)  call send_tile_data(id_rsoil_fast, fast_C_loss_rate(:)/dz(1:num_l), diag)
  if (id_rsoil_slow>0)  call send_tile_data(id_rsoil_slow, slow_C_loss_rate(:)/dz(1:num_l), diag)
  if (id_rsoil_deadmic>0) call send_tile_data(id_rsoil_deadmic, dead_microbe_C_loss_rate(:)/dz(1:num_l), diag)
  if (id_rsoil_leaflitter_fast>0) call send_tile_data(id_rsoil_leaflitter_fast, leaflitter_fast_C_loss_rate, diag)
  if (id_rsoil_leaflitter_slow>0) call send_tile_data(id_rsoil_leaflitter_slow, leaflitter_slow_C_loss_rate, diag)
  if (id_rsoil_leaflitter_deadmic>0) call send_tile_data(id_rsoil_leaflitter_deadmic, leaflitter_deadmic_C_loss_rate, diag)
  if (id_rsoil_finewoodlitter_fast>0) call send_tile_data(id_rsoil_finewoodlitter_fast, finewoodlitter_fast_C_loss_rate, diag)
  if (id_rsoil_finewoodlitter_slow>0) call send_tile_data(id_rsoil_finewoodlitter_slow, finewoodlitter_slow_C_loss_rate, diag)
  if (id_rsoil_finewoodlitter_deadmic>0) call send_tile_data(id_rsoil_finewoodlitter_deadmic, finewoodlitter_deadmic_C_loss_rate, diag)
  if (id_rsoil_coarsewoodlitter_fast>0) call send_tile_data(id_rsoil_coarsewoodlitter_fast, coarsewoodlitter_fast_C_loss_rate, diag)
  if (id_rsoil_coarsewoodlitter_slow>0) call send_tile_data(id_rsoil_coarsewoodlitter_slow, coarsewoodlitter_slow_C_loss_rate, diag)
  if (id_rsoil_coarsewoodlitter_deadmic>0) call send_tile_data(id_rsoil_coarsewoodlitter_deadmic, coarsewoodlitter_deadmic_C_loss_rate, diag)
  call send_tile_data(id_rsoil, vegn%rh, diag)
  call send_tile_data(id_rh, vegn%rh/seconds_per_year, diag)
  ! TODO: arithmetic averaging of A doesn't seem correct; we need to invent something better,
  !       e.g. weight it with the carbon loss, or something like that
  if (id_asoil>0) call send_tile_data(id_asoil, sum(A(:))/size(A(:)), diag)


  if (id_leaflitter_dissolved_fast>0) call send_tile_data(id_leaflitter_dissolved_fast,leaflitter_C_dissolved(1)/dt_fast_yr,diag)
  if (id_leaflitter_dissolved_slow>0) call send_tile_data(id_leaflitter_dissolved_slow,leaflitter_C_dissolved(2)/dt_fast_yr,diag)
  if (id_leaflitter_dissolved_deadmic>0) call send_tile_data(id_leaflitter_dissolved_deadmic,leaflitter_C_dissolved(3)/dt_fast_yr,diag)
  if (id_finewoodlitter_dissolved_fast>0) call send_tile_data(id_finewoodlitter_dissolved_fast,finewoodlitter_C_dissolved(1)/dt_fast_yr,diag)
  if (id_finewoodlitter_dissolved_slow>0) call send_tile_data(id_finewoodlitter_dissolved_slow,finewoodlitter_C_dissolved(2)/dt_fast_yr,diag)
  if (id_finewoodlitter_dissolved_deadmic>0) call send_tile_data(id_finewoodlitter_dissolved_deadmic,finewoodlitter_C_dissolved(3)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_dissolved_fast>0) call send_tile_data(id_coarsewoodlitter_dissolved_fast,coarsewoodlitter_C_dissolved(1)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_dissolved_slow>0) call send_tile_data(id_coarsewoodlitter_dissolved_slow,coarsewoodlitter_C_dissolved(2)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_dissolved_deadmic>0) call send_tile_data(id_coarsewoodlitter_dissolved_deadmic,coarsewoodlitter_C_dissolved(3)/dt_fast_yr,diag)
  if (id_dissolved_fast>0) call send_tile_data(id_dissolved_fast,C_dissolved(1,:)/dt_fast_yr/dz(1:num_l),diag)
  if (id_dissolved_slow>0) call send_tile_data(id_dissolved_slow,C_dissolved(2,:)/dt_fast_yr/dz(1:num_l),diag)
  if (id_dissolved_deadmic>0) call send_tile_data(id_dissolved_deadmic,C_dissolved(3,:)/dt_fast_yr/dz(1:num_l),diag)

  if (id_leaflitter_deposited_fast>0) call send_tile_data(id_leaflitter_deposited_fast,leaflitter_C_deposited(1)/dt_fast_yr,diag)
  if (id_leaflitter_deposited_slow>0) call send_tile_data(id_leaflitter_deposited_slow,leaflitter_C_deposited(2)/dt_fast_yr,diag)
  if (id_leaflitter_deposited_deadmic>0) call send_tile_data(id_leaflitter_deposited_deadmic,leaflitter_C_deposited(3)/dt_fast_yr,diag)
  if (id_finewoodlitter_deposited_fast>0) call send_tile_data(id_finewoodlitter_deposited_fast,finewoodlitter_C_deposited(1)/dt_fast_yr,diag)
  if (id_finewoodlitter_deposited_slow>0) call send_tile_data(id_finewoodlitter_deposited_slow,finewoodlitter_C_deposited(2)/dt_fast_yr,diag)
  if (id_finewoodlitter_deposited_deadmic>0) call send_tile_data(id_finewoodlitter_deposited_deadmic,finewoodlitter_C_deposited(3)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_deposited_fast>0) call send_tile_data(id_coarsewoodlitter_deposited_fast,coarsewoodlitter_C_deposited(1)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_deposited_slow>0) call send_tile_data(id_coarsewoodlitter_deposited_slow,coarsewoodlitter_C_deposited(2)/dt_fast_yr,diag)
  if (id_coarsewoodlitter_deposited_deadmic>0) call send_tile_data(id_coarsewoodlitter_deposited_deadmic,coarsewoodlitter_C_deposited(3)/dt_fast_yr,diag)

  if (id_deposited_fast>0) call send_tile_data(id_deposited_fast,C_deposited(1,:)/dt_fast_yr/dz(1:num_l),diag)
  if (id_deposited_slow>0) call send_tile_data(id_deposited_slow,C_deposited(2,:)/dt_fast_yr/dz(1:num_l),diag)
  if (id_deposited_deadmic>0) call send_tile_data(id_deposited_deadmic,C_deposited(3,:)/dt_fast_yr/dz(1:num_l),diag)
end subroutine Dsdt_CORPSE


! ============================================================================
subroutine Dsdt_CENTURY(vegn, soil, diag, soilt, theta)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  type(diag_buff_type), intent(inout) :: diag
  real                , intent(in)    :: soilt ! average soil temperature, deg K
  real                , intent(in)    :: theta ! average soil moisture

  real :: fast_C_loss(size(soil%fast_soil_C))
  real :: slow_C_loss(size(soil%slow_soil_C))
  real :: A          (size(soil%slow_soil_C)) ! decomp rate reduction due to moisture and temperature

  select case (soil_carbon_option)
  case(SOILC_CENTURY)
      A(:) = A_function(soilt, theta)
  case(SOILC_CENTURY_BY_LAYER)
      A(:) = A_function(soil%T, soil_theta(soil))
  case default
    call error_mesg('Dsdt_CENTURY','The value of soil_carbon_option is invalid. This should never happen. See developer.',FATAL)
  end select

  fast_C_loss = soil%fast_soil_C(:)*A*K1*dt_fast_yr;
  slow_C_loss = soil%slow_soil_C(:)*A*K2*dt_fast_yr;

  soil%fast_soil_C = soil%fast_soil_C - fast_C_loss;
  soil%slow_soil_C = soil%slow_soil_C - slow_C_loss;

  ! for budget check
  vegn%fsc_out = vegn%fsc_out + sum(fast_C_loss(:));
  vegn%ssc_out = vegn%ssc_out + sum(slow_C_loss(:));

  ! loss of C to atmosphere and leaching
  vegn%rh = sum(fast_C_loss(:)+slow_C_loss(:))/dt_fast_yr;

  ! accumulate decomposition rate reduction for the soil carbon restart output
  soil%asoil_in(:) = soil%asoil_in(:) + A(:)

  ! ---- diagnostic section
  if (id_rsoil_fast>0)  call send_tile_data(id_rsoil_fast, fast_C_loss(:)/(dz(1:num_l)*dt_fast_yr), diag)
  if (id_rsoil_slow>0)  call send_tile_data(id_rsoil_slow, slow_C_loss(:)/(dz(1:num_l)*dt_fast_yr), diag)
  call send_tile_data(id_rsoil, vegn%rh, diag)
  call send_tile_data(id_rh, vegn%rh/seconds_per_year, diag)
  ! TODO: arithmetic averaging of A doesn't seem correct; we need to invent something better,
  !       e.g. weight it with the carbon loss, or something like that
  if (id_asoil>0) call send_tile_data(id_asoil, sum(A(:))/size(A(:)), diag)

end subroutine Dsdt_CENTURY


! ============================================================================
subroutine soil_push_down_excess ( soil, diag, lrunf_nu, hlrunf_nu, frunf, hfrunf)
  type(soil_tile_type), intent(inout) :: soil
  type(diag_buff_type), intent(inout) :: diag
  real, intent(out) :: lrunf_nu, hlrunf_nu
  real, intent(out) :: frunf, hfrunf ! frozen runoff (mm/s); frozen runoff heat (W/m^2)

  ! ---- local vars ----------------------------------------------------------
  real      :: &
     liq_frac, excess_wat, excess_liq, excess_ice, excess_t, &
     h1, h2, summax, space_avail, liq_placed, ice_placed
  integer :: k,l

  liq_frac=0;excess_wat=0;excess_liq=0;excess_ice=0;h1=0;h2=0
  l = 1
  summax = max(0.,soil%wl(l))+max(0.,soil%ws(l))
  if (summax > 0) then
     liq_frac = max(0.,soil%wl(l)) / summax
  else
     liq_frac = 1
  endif
  excess_wat = max(0., soil%wl(l) + soil%ws(l) &
       - dens_h2o*dz(l)*soil%vwc_max(l) )
  excess_liq = excess_wat*liq_frac
  excess_ice = excess_wat-excess_liq
  excess_t   = soil%T(l)
  soil%wl(l) = soil%wl(l) - excess_liq
  soil%ws(l) = soil%ws(l) - excess_ice
  call send_tile_data(id_excess, excess_wat/delta_time, diag)

!  if(is_watch_cell()) then
  if(is_watch_point()) then
     write(*,*) ' ##### push_down_excess input #####'
      call get_current_point(k=k)
      write(*,*) 'For watch_cell, tile=',k
     __DEBUG3__(l,summax,liq_frac)
     __DEBUG3__(soil%vwc_max(l),excess_liq,excess_ice)
     __DEBUG2__(dens_h2o,dz(l))
  endif

  do l = 2, num_l
     if (excess_liq+excess_ice>0) then
        space_avail = dens_h2o*dz(l)*soil%vwc_max(l) &
             - (soil%wl(l) + soil%ws(l))
        liq_placed = max(min(space_avail, excess_liq), 0.)
        ice_placed = max(min(space_avail-liq_placed, excess_ice), 0.)
        h1 = (soil%heat_capacity_dry(l)*dz(l) &
             + csw*soil%ws(l) + clw*soil%wl(l))
        h2 = liq_placed*clw+ice_placed*csw
        soil%T(l) = (h1 * soil%T(l) &
             + h2 * excess_T )  / (h1+h2)
        soil%wl(l) = soil%wl(l) + liq_placed
        soil%ws(l) = soil%ws(l) + ice_placed
        excess_liq = excess_liq - liq_placed
        excess_ice = excess_ice - ice_placed
     endif
  enddo

! to avoid adding frozen runoff to soil interface, melt all remaining
! excess ice, even if it results in supercooled liquid runoff
  if (supercooled_rnu) then
     lrunf_nu  = (excess_liq+excess_ice) / delta_time
     hlrunf_nu = (  excess_liq*clw*(excess_T-tfreeze)  &
                  + excess_ice*csw*(excess_T-tfreeze)  &
                  - hlf_factor*hlf*excess_ice                   ) / delta_time
     frunf = 0.
     hfrunf = 0.
  else
!! ZMS: Change this to propagate frozen runoff to avoid lake crashes from supercooled water.
     lrunf_nu = excess_liq / delta_time
     hlrunf_nu = excess_liq*clw*(excess_T-tfreeze) / delta_time
     frunf = excess_ice / delta_time
     hfrunf = excess_ice*csw*(excess_T-tfreeze) / delta_time
  end if

!  if(is_watch_cell()) then
  if(is_watch_point()) then
     write(*,*) ' ##### push_down_excess output #####'
     __DEBUG2__(lrunf_nu,hlrunf_nu)
     write(*,*) 'For watch_cell'
     __DEBUG3__(soil%hidx_k, frunf, hfrunf)
     do l = 1, num_l
        write(*,'(x,a,x,i2.2)',advance='NO')' level=', l
        call dpri(' T =',soil%T(l))
        call dpri(' Th=', (soil%ws(l)+soil%wl(l))/(dens_h2o*dz(l)))
        call dpri(' wl=', soil%wl(l))
        call dpri(' ws=', soil%ws(l))
        write(*,*)
     enddo
  endif
end subroutine soil_push_down_excess

! ============================================================================
  subroutine RICHARDS_clean(soil, psi, DThDP, hyd_cond, DKDP, div, &
                lprec_eff, Dpsi_min, Dpsi_max, dt_richards, &
                 dPsi, dW_l, flow, lrunf_ie)
! ZMS How to set prevent lrunf_ie when surface soil water < pond?

  type(soil_tile_type), intent(inout)   :: soil
  real, intent(in),  dimension(num_l)   :: psi, DThDP, hyd_cond, DKDP, div
  real, intent(in)                      :: lprec_eff, Dpsi_min, Dpsi_max
  real, intent(in)                      :: dt_richards
  real, intent(out), dimension(num_l)   :: dPsi, dW_l
  real, intent(out), dimension(num_l+1) :: flow
  real, intent(out)                     :: lrunf_ie
  ! ---- local vars ----------------------------------------------------------
  integer l, ipt, jpt, kpt, fpt, l_dest
  real, dimension(num_l-1) :: del_z, K, DKDPm, DKDPp, grad, eee, fff
  real aaa, bbb, ccc, ddd, xxx, dpsi_alt, w_shortage, adj
  logical flag

  flag = .false.
  flow(1) = dt_richards*lprec_eff
  do l = 1, num_l-1
     del_z(l) = zfull(l+1)-zfull(l)
     if (harmonic_mean_K) then
         K(l) = (dz(l)+dz(l+1))/(dz(l)/hyd_cond(l)+dz(l+1)/hyd_cond(l+1))
     else
         K(l) = 0.5*(hyd_cond(l)+hyd_cond(l+1))
     endif
     DKDPm(l) = 0. ! 0.5*DKDP(l)
     DKDPp(l) = 0. ! 0.5*DKDP(l+1)
!        K(l) = hyd_cond(l)
!        DKDPm(l) = DKDP(l)
!        DKDPp(l) = 0
      grad(l)  = (psi(l+1)-psi(l))/del_z(l) - 1
  enddo

  if (is_watch_point()) then
     write(*,*) '##### soil_step_2 checkpoint 3.1 #####'
     do l = 1, num_l
        write(*,'(x,a,x,i2.2,x,a,100(x,g23.16))') 'level=', l, 'DThDP,hyd_cond,psi,DKDP', &
             DThDP(l),&
             hyd_cond(l),&
             psi(l),&
             DKDP(l)
     enddo
     do l = 1, num_l-1
        write(*,'(a,i2.2,1x,a,100(2x,g23.16))') 'interface=', l, 'K,DKDPm,DKDPp,grad,del_z', &
             K(l),&
             DKDPm(l),&
             DKDPp(l),&
             grad(l)
     enddo
  endif


  l = num_l
  xxx = dens_h2o*dz(l)*DThDP(l)/dt_richards
  aaa =     - ( K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
  bbb = xxx - (-K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1) )
  ddd = - K(l-1) *grad(l-1) - div(l)
  eee(l-1) = -aaa/bbb
  fff(l-1) =  ddd/bbb

  if(is_watch_point()) then
       write(*,'(a,i2.2,100(2x,g23.16))') 'l,a,b, ,d', l,aaa, bbb,ddd
  endif

  do l = num_l-1, 2, -1
    xxx = dens_h2o*dz(l)*DThDP(l)/dt_richards
    aaa = - ( K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
    bbb = xxx-( -K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1)&
                -K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
    ccc =   - (  K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
    ddd =       K(l)*grad(l) - K(l-1)*grad(l-1) &
                          - div(l)
    eee(l-1) =                    -aaa/(bbb+ccc*eee(l))
    fff(l-1) =  (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
    if(is_watch_point()) then
       write(*,'(a,i2.2,100(2x,g23.16))') 'l,a,b,c,d', l,aaa, bbb,ccc,ddd
    endif
  enddo

  l = 1
  xxx = dens_h2o*dz(l)*DThDP(l)/dt_richards
  bbb = xxx - ( -K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
  ccc =     - (  K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
  ddd =          flow(1)/dt_richards +    K(l)     *grad(l) &
                          - div(l)

  if (Dpsi_min.ge.Dpsi_max) call error_mesg(module_name, '=== Dpsi_min.ge.Dpsi_max', FATAL)

  IF (bbb+ccc*eee(l) .EQ. 0.) THEN
      call get_current_point(ipt,jpt,kpt,fpt)
      write(*,*) '===richards b+ce=0 ===','at point ',ipt,jpt,kpt,fpt
        write(*,*) 'bbb+ccc*eee(l) .EQ. 0.'
        write(*,*) 'bbb', bbb
        write(*,*) 'ccc', ccc
        write(*,*) 'ddd', ddd
        write(*,*) 'eee(l)', eee(l)
        write(*,*) 'fff(l)', fff(l)
        write(*,*) 'dPsi(l)', dPsi(l)
        write(*,*) 'dPsi(l)', dPsi(l)
        write(*,*) 'Dpsi_min', Dpsi_min
        write(*,*) 'Dpsi_max', Dpsi_max
    call error_mesg(module_name, 'b+ce=0 in soil-water equations', FATAL)
  ENDIF

     dPsi(l) = (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
     if(is_watch_point()) then
        write(*,*) 'bbb+ccc*eee(l) .NE. 0.'
        write(*,*) 'bbb', bbb
        write(*,*) 'ccc', ccc
        write(*,*) 'ddd', ddd
        write(*,*) 'eee(l)', eee(l)
        write(*,*) 'fff(l)', fff(l)
        write(*,*) 'dPsi(l)', dPsi(l)
        write(*,*) 'dPsi(l)', dPsi(l)
        write(*,*) 'Dpsi_min', Dpsi_min
        write(*,*) 'Dpsi_max', Dpsi_max
     endif
     if (dPsi(l).lt.Dpsi_max) then
        lrunf_ie = 0.
     else
        dPsi(l) = Dpsi_max
        if (div_bug_fix) then
            flow(l) = (dPsi(l)*(bbb+ccc*eee(l))+ccc*fff(l)+div(l) &
                     - K(l)*grad(l))*dt_richards
        else
            flow(l) = (dPsi(l)*(bbb+ccc*eee(l))+ccc*fff(l) &
                     - K(l)*grad(l))*dt_richards
        endif
        lrunf_ie = lprec_eff - flow(l)/dt_richards
     endif

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,'(a,i2.2,100(2x,g23.16))') 'l,  b,c,d', l, bbb,ccc,ddd
     write(*,*) ' ##### soil_step_2 checkpoint 3.2 #####'
     write(*,*) 'ie:', lrunf_ie
     do l = 1, num_l-1
        write(*,'(a,i2.2,100(2x,g23.16))') 'l,eee(l),fff(l)',l,eee(l),fff(l)
     enddo
     write(*,*) 'DThDP(1)', DThDP(1)
     write(*,*) 'K(1)', K(1)
     write(*,*) 'grad(1)', grad(1)
     write(*,*) 'ddd(1)', ddd
     write(*,*) 'ccc(1)', ccc
     write(*,*) 'bbb(1)', bbb
     write(*,*) 'dPsi(1)', dPsi(1)
     write(*,*) 'Psi(1)', Psi(1)
     write(*,*) 'div(1)', div(1)
  endif

  do l = 1, num_l-1
     dPsi(l+1) = eee(l)*dPsi(l) + fff(l)
     flow(l+1) = dt_richards*( &
         -K(l)*(grad(l)&
         +(DPsi(l+1)-DPsi(l))/ del_z(l)) &
         -grad(l)*(DKDPp(l)*Dpsi(l+1)+ &
                         DKDPm(l)*Dpsi(l) )  )
     dW_l(l) = flow(l) - flow(l+1) - div(l)*dt_richards
  enddo
  flow(num_l+1) = 0.
  dW_l(num_l) = flow(num_l) - flow(num_l+1) &
                          - div(num_l)*dt_richards

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,*) ' ##### soil_step_2 checkpoint 3.21 #####'
     do l = 1, num_l
        write(*,'(i2.2,100(2x,a,g23.16))') l,&
             ' dW_l=', dW_l(l),&
             ' flow=', flow(l),&
             ' div=', div(l)
     enddo
  endif

! Check/adjust for negative water content at surface.
  if (dPsi(1).lt.Dpsi_min) then
     if (verbose) then
         call get_current_point(ipt,jpt,kpt,fpt)
         write(*,*) '=== warning: dPsi=',dPsi(1),'<min=',dPsi_min,'at',ipt,jpt,kpt,fpt
       endif
     w_shortage = -(soil%wl(1)+dW_l(1))
     l_dest = 1
     call move_up(dW_l, flow, w_shortage, num_l, l_dest)
  endif

! Adjust for negative water content in subsurface.
  if (fix_neg_subsurface_wl) then
    do l=2, num_l
      if ((soil%wl(l)+dW_l(l))/(dens_h2o*dz(l)*soil%pars%vwc_sat) < thetathresh) then
        call get_current_point(ipt,jpt,kpt,fpt)
        write(*,*) '=== warning: fixing neg wl=',soil%wl(l)+dW_l(l),'at',l,ipt,jpt,kpt,fpt
        l_dest = l
        w_shortage = -(soil%wl(l_dest)+dW_l(l_dest))
        call move_up(dW_l, flow, w_shortage, num_l, l_dest)
      endif
    enddo
  endif

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,*) ' ##### soil_step_2 checkpoint 3.22 #####'
     do l = 1, num_l
        write(*,'(i2.2,100(2x,a,g23.16))') l,&
             ' dW_l=', dW_l(l),&
             ' flow=', flow(l),&
             ' div=', div(l)
     enddo
  endif

! Check for negative runoff:
  IF (lrunf_ie < lrunf_ie_min) THEN
     call get_current_point(ipt,jpt,kpt,fpt)
     write(*,*) 'note: at point ',ipt,jpt,kpt,fpt,'lrunf_ie=',lrunf_ie,' < lrunf_ie_min=',lrunf_ie_min
     call error_mesg(module_name, 'lrunf_ie < lrunf_ie_min', FATAL)
  ENDIF

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.3 ***** '
     write(*,*) 'psi_sat',soil%pars%psi_sat_ref
     write(*,*) 'Dpsi_max',Dpsi_max
     do l = 1, num_l
        write(*,'(i2.2,100(2x,a,g23.16))') l, &
             'Th=', (soil%ws(l) +soil%wl(l)+dW_l(l))/(dens_h2o*dz(l)), &
             'wl=', soil%wl(l)+dW_l(l), &
             'ws=', soil%ws(l), &
             'dW_l=', dW_l(l), &
             'dPsi=', dPsi(l), &
             'flow=', flow(l)
     enddo
  endif

end subroutine richards_clean

! ============================================================================
  subroutine move_up(dW_l, flow, w_shortage, num_l, l_dest)
  real, intent(inout), dimension(num_l)   :: dW_l
  real, intent(inout), dimension(num_l+1) :: flow
  real, intent(in)                        ::  w_shortage
  integer, intent(in)                     ::  num_l, l_dest
  ! ---- local vars ----------------------------------------------------------
  integer l, l_source
  real dW_l_source, w_to_move_up
     l_source = l_dest
     dW_l_source = -1.e20
     do l = l_dest+1, num_l
        if (dW_l(l).gt.dW_l_source) then
           l_source = l
           dW_l_source = dW_l(l)
        endif
     enddo
     w_to_move_up = min(dW_l_source, w_shortage)
     w_to_move_up = max(w_to_move_up, 0.)
     write(*,*) 'l_dest,l_source=',l_dest,l_source
     write(*,*) 'dW_l_source=',dW_l_source
     write(*,*) 'w_shortage',w_shortage
     write(*,*) 'w_to_move_up=',w_to_move_up
     if (l_source.gt.l_dest) then
        dW_l(l_dest)   = dW_l(l_dest)   + w_to_move_up
        dW_l(l_source) = dW_l(l_source) - w_to_move_up
        do l = l_dest+1, l_source
           flow(l) = flow(l) - w_to_move_up
        enddo
     endif
  end subroutine move_up
! ============================================================================
  subroutine RICHARDS(soil, psi, DThDP, hyd_cond, DKDP, div, &
                lprec_eff, Dpsi_min, Dpsi_max, dt_richards, &
                 dPsi, dW_l, flow, lrunf_ie)
! ZMS How to set prevent lrunf_ie when surface soil water < pond?

  type(soil_tile_type), intent(inout)   :: soil
  real, intent(in),  dimension(num_l)   :: psi, DThDP, hyd_cond, DKDP, div
  real, intent(in)                      :: lprec_eff, Dpsi_min, Dpsi_max
  real, intent(in)                      :: dt_richards
  real, intent(out), dimension(num_l)   :: dPsi, dW_l
  real, intent(out), dimension(num_l+1) :: flow
  real, intent(out)                     :: lrunf_ie
  ! ---- local vars ----------------------------------------------------------
  integer l, ipt, jpt, kpt, fpt, l_dest
  real, dimension(num_l-1) :: del_z, K, DKDPm, DKDPp, grad, eee, fff
  real aaa, bbb, ccc, ddd, xxx, dpsi_alt, w_shortage, adj
  logical flag

  flag = .false.
  flow(1) = dt_richards*lprec_eff
  do l = 1, num_l-1
     del_z(l) = zfull(l+1)-zfull(l)
     if (harmonic_mean_K) then
         K(l) = (dz(l)+dz(l+1))/(dz(l)/hyd_cond(l)+dz(l+1)/hyd_cond(l+1))
     else
         K(l) = 0.5*(hyd_cond(l)+hyd_cond(l+1))
     endif
     DKDPm(l) = 0. ! 0.5*DKDP(l)
     DKDPp(l) = 0. ! 0.5*DKDP(l+1)
!        K(l) = hyd_cond(l)
!        DKDPm(l) = DKDP(l)
!        DKDPp(l) = 0
      grad(l)  = (psi(l+1)-psi(l))/del_z(l) - 1
  enddo

  if (is_watch_point()) then
     write(*,*) '##### soil_step_2 checkpoint 3.1 #####'
     do l = 1, num_l
        write(*,'(x,a,x,i2.2,x,a,100(x,g23.16))') 'level=', l, 'DThDP,hyd_cond,psi,DKDP', &
             DThDP(l),&
             hyd_cond(l),&
             psi(l),&
             DKDP(l)
     enddo
     do l = 1, num_l-1
        write(*,'(a,i2.2,1x,a,100(2x,g23.16))') 'interface=', l, 'K,DKDPm,DKDPp,grad,del_z', &
             K(l),&
             DKDPm(l),&
             DKDPp(l),&
             grad(l)
     enddo
  endif


  l = num_l
  xxx = dens_h2o*dz(l)*DThDP(l)/dt_richards
  aaa =     - ( K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
  bbb = xxx - (-K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1) )
  ddd = - K(l-1) *grad(l-1) - div(l)
  eee(l-1) = -aaa/bbb
  fff(l-1) =  ddd/bbb

  if(is_watch_point()) then
       write(*,'(a,i2.2,100(2x,g23.16))') 'l,a,b, ,d', l,aaa, bbb,ddd
  endif

  do l = num_l-1, 2, -1
    xxx = dens_h2o*dz(l)*DThDP(l)/dt_richards
    aaa = - ( K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
    bbb = xxx-( -K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1)&
                -K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
    ccc =   - (  K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
    ddd =       K(l)*grad(l) - K(l-1)*grad(l-1) &
                          - div(l)
    eee(l-1) =                    -aaa/(bbb+ccc*eee(l))
    fff(l-1) =  (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
    if(is_watch_point()) then
       write(*,'(a,i2.2,100(2x,g23.16))') 'l,a,b,c,d', l,aaa, bbb,ccc,ddd
    endif
  enddo

  l = 1
  xxx = dens_h2o*dz(l)*DThDP(l)/dt_richards
  bbb = xxx - ( -K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
  ccc =     - (  K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
  ddd =          flow(1)/dt_richards +    K(l)     *grad(l) &
                          - div(l)

  if (Dpsi_min.ge.Dpsi_max) call error_mesg(module_name, '=== Dpsi_min.ge.Dpsi_max', FATAL)

  IF (bbb+ccc*eee(l) .NE. 0.) THEN
     dPsi(l) = (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
     if(is_watch_point()) then
        write(*,*) 'bbb+ccc*eee(l) .NE. 0.'
        write(*,*) 'bbb', bbb
        write(*,*) 'ccc', ccc
        write(*,*) 'ddd', ddd
        write(*,*) 'eee(l)', eee(l)
        write(*,*) 'fff(l)', fff(l)
        write(*,*) 'dPsi(l)', dPsi(l)
        write(*,*) 'dPsi(l)', dPsi(l)
        write(*,*) 'Dpsi_min', Dpsi_min
        write(*,*) 'Dpsi_max', Dpsi_max
     endif
     if (verbose.and.dPsi(l).lt.Dpsi_min) then
         call get_current_point(ipt,jpt,kpt,fpt)
         write(*,*) '=== warning: dPsi=',dPsi(l),'<min=',dPsi_min,'at',ipt,jpt,kpt,fpt
     endif
     if ((dPsi(l).gt.Dpsi_min.or.no_min_Dpsi) .and. dPsi(l).lt.Dpsi_max) then
        lrunf_ie = 0.
     else
        dPsi(l) = min (dPsi(l), Dpsi_max)
        if (dPsi(l).lt.Dpsi_min.and.(.not.no_min_Dpsi)) then
           flag = .true.
           dPsi(l) = Dpsi_min
        endif
        if (div_bug_fix) then
           flow(l) = (dPsi(l)*(bbb+ccc*eee(l))+ccc*fff(l)+div(l) &
                   - K(l)*grad(l))*dt_richards
        else
           flow(l) = (dPsi(l)*(bbb+ccc*eee(l))+ccc*fff(l) &
                   - K(l)*grad(l))*dt_richards
        endif
        lrunf_ie = lprec_eff - flow(l)/dt_richards
        if (lrunf_ie.lt.-lrunf_ie_tol) then
           if (verbose) then
              call get_current_point(ipt,jpt,kpt,fpt)
              write(*,*) '=== warning: rie= ',lrunf_ie,'<0 at',ipt,jpt,kpt,fpt
           endif
           if (.not.allow_negative_rie) then
              flag = .true.
              dpsi_alt = (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
              if (verbose) write(*,*) '=== rie= ',lrunf_ie,'<0 reset to 0 and dPsi=',dPsi(l),' reset to ',dpsi_alt,'at ',ipt,jpt,kpt,fpt
              dPsi(l) = dpsi_alt
              lrunf_ie = 0.
              flow(l) = lprec_eff*dt_richards
           endif
        endif
     endif
  ELSE
     call error_mesg(module_name, 'b+ce=0 in soil-water equations', FATAL)
     if (verbose) then
       call get_current_point(ipt,jpt,kpt,fpt)
       write(*,*) '===richards b+ce=0 ===','at point ',ipt,jpt,kpt,fpt
       endif
     if(is_watch_point()) then
        write(*,*) 'bbb+ccc*eee(l) .EQ. 0.'
        write(*,*) 'bbb', bbb
        write(*,*) 'ccc', ccc
        write(*,*) 'ddd', ddd
        write(*,*) 'eee(l)', eee(l)
        write(*,*) 'fff(l)', fff(l)
        write(*,*) 'dPsi(l)', dPsi(l)
        write(*,*) 'dPsi(l)', dPsi(l)
        write(*,*) 'Dpsi_min', Dpsi_min
        write(*,*) 'Dpsi_max', Dpsi_max
     endif
     dPsi(l) = Dpsi_max
     flow(l) = (dPsi(l)*(bbb+ccc*eee(l))+ccc*fff(l) &
                   - K(l)*grad(l))*dt_richards
     lrunf_ie = lprec_eff - flow(l)/dt_richards
     if (lrunf_ie.lt.-lrunf_ie_tol) then
        if (verbose) then
           call get_current_point(ipt,jpt,kpt,fpt)
           write(*,*) '===richards b+ce=0 AND lrunf_ie.lt.-lrunf_ie_tol at point ',ipt,jpt,kpt,fpt
           write(*,*) '===richards b+ce=0 AND lrunf_ie.lt.-lrunf_ie_tol at point ',ipt,jpt,kpt,fpt,'rie=',lrunf_ie
        endif
        if(.not.allow_negative_rie) then
           flag = .true.
           dpsi_alt = 0.
           if (verbose) then
             write(*,*) '===richards b+ce=0 AND lrunf_ie.lt.-lrunf_ie_tol at point ',ipt,jpt,kpt,fpt,'rie=',lrunf_ie,' reset to 0'
             write(*,*) '===richards b+ce=0 AND lrunf_ie.lt.-lrunf_ie_tol at point ',ipt,jpt,kpt,fpt,'dPsi=',dPsi(l),' reset to',dpsi_alt
             write(*,*) '===richards b+ce=0 AND lrunf_ie.lt.-lrunf_ie_tol at point ',ipt,jpt,kpt,fpt,'dPsi_min/max=',dPsi_min,dPsi_max
           endif
           dPsi(l) = dpsi_alt
           lrunf_ie = 0.
           flow(l) = lprec_eff*dt_richards
        endif
     endif
  ENDIF

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,'(a,i2.2,100(2x,g23.16))') 'l,  b,c,d', l, bbb,ccc,ddd
     write(*,*) ' ##### soil_step_2 checkpoint 3.2 #####'
     write(*,*) 'ie:', lrunf_ie
     do l = 1, num_l-1
        write(*,'(a,i2.2,100(2x,g23.16))') 'l,eee(l),fff(l)',l,eee(l),fff(l)
     enddo
     write(*,*) 'DThDP(1)', DThDP(1)
     write(*,*) 'K(1)', K(1)
     write(*,*) 'grad(1)', grad(1)
     write(*,*) 'ddd(1)', ddd
     write(*,*) 'ccc(1)', ccc
     write(*,*) 'bbb(1)', bbb
     write(*,*) 'dPsi(1)', dPsi(1)
     write(*,*) 'Psi(1)', Psi(1)
     write(*,*) 'div(1)', div(1)
  endif

  do l = 1, num_l-1
     dPsi(l+1) = eee(l)*dPsi(l) + fff(l)
     flow(l+1) = dt_richards*( &
         -K(l)*(grad(l)&
         +(DPsi(l+1)-DPsi(l))/ del_z(l)) &
         -grad(l)*(DKDPp(l)*Dpsi(l+1)+ &
                         DKDPm(l)*Dpsi(l) )  )
     dW_l(l) = flow(l) - flow(l+1) - div(l)*dt_richards
  enddo
  flow(num_l+1) = 0.
  dW_l(num_l) = flow(num_l) - flow(num_l+1) &
                          - div(num_l)*dt_richards

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,*) ' ##### soil_step_2 checkpoint 3.21 #####'
     do l = 1, num_l
        write(*,'(i2.2,100(2x,a,g23.16))') l,&
             ' dW_l=', dW_l(l),&
             ' flow=', flow(l),&
             ' div=', div(l)
     enddo
  endif

  if (flag) then
     w_shortage=-(soil%wl(1)+dW_l(1))
     l_dest = 1
     call move_up(dW_l, flow, w_shortage, num_l, l_dest)
  endif

! Adjust for negative water content in subsurface.
  if (fix_neg_subsurface_wl) then
    do l=2, num_l
      if ((soil%wl(l)+dW_l(l))/(dens_h2o*dz(l)*soil%pars%vwc_sat) < thetathresh) then
        call get_current_point(ipt,jpt,kpt,fpt)
        write(*,*) '=== warning: fixing neg wl=',soil%wl(l)+dW_l(l),'at',l,ipt,jpt,kpt,fpt
        l_dest = l
        w_shortage = -(soil%wl(l_dest)+dW_l(l_dest))
        call move_up(dW_l, flow, w_shortage, num_l, l_dest)
      endif
    enddo
  endif


  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,*) ' ##### soil_step_2 checkpoint 3.22 #####'
     do l = 1, num_l
        write(*,'(i2.2,100(2x,a,g23.16))') l,&
             ' dW_l=', dW_l(l),&
             ' flow=', flow(l),&
             ' div=', div(l)
     enddo
  endif

! In rare situations where lrunf_ie is large and negative, clip any liquid supersaturation
! layer by layer and recompute lrunf_ie (this is not good, since it ignores 'comp'):
  IF (lrunf_ie < lrunf_ie_min) THEN
     call get_current_point(ipt,jpt,kpt,fpt)
     write(*,*) 'note: at point ',ipt,jpt,kpt,fpt,' clip triggered by lrunf_ie=',lrunf_ie
     call error_mesg(module_name, 'lrunf_ie < lrunf_ie_min', FATAL)
     do l = num_l, 1, -1
        adj = max(dW_l(l)+soil%ws(l)+soil%wl(l) &
             - soil%pars%vwc_sat*dz(l)*dens_h2o, 0. )

        if(is_watch_point()) then
           write(*,*) '3.22 l=', l,&
                ' soil%wl=',soil%wl(l),  &
                ' soil%ws=',soil%ws(l) , &
                ' soil%pars%vwc_sat=', soil%pars%vwc_sat, &
                ' dz=', dz(l), &
                ' adj=', adj
        endif

        adj = min(adj, max(0.,soil%wl(l)))

        if(is_watch_point()) then
           write(*,*) '3.23 l=', l, ' adj=', adj
        endif

        dW_l(l) = dW_l(l) - adj
        flow(l) = flow(l+1) + dW_l(l) + div(l)*dt_richards
     enddo
     lrunf_ie = lprec_eff - flow(1)/dt_richards

  ENDIF

  if(is_watch_point().or.(flag.and.write_when_flagged)) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.3 ***** '
     write(*,*) 'psi_sat',soil%pars%psi_sat_ref
     write(*,*) 'Dpsi_max',Dpsi_max
     do l = 1, num_l
        write(*,'(i2.2,100(2x,a,g23.16))') l, &
             'Th=', (soil%ws(l) +soil%wl(l)+dW_l(l))/(dens_h2o*dz(l)), &
             'wl=', soil%wl(l)+dW_l(l), &
             'ws=', soil%ws(l), &
             'dW_l=', dW_l(l), &
             'dPsi=', dPsi(l), &
             'flow=', flow(l)
     enddo
  endif

end subroutine richards

! ============================================================================
  subroutine advection(soil, flow, dW_l, tflow, d_GW, div, delta_time)
  type(soil_tile_type), intent(inout) :: soil
  real, intent(in), dimension(:) :: flow  ! water tendency downwards into layer [mm]
  real, intent(in), dimension(:) :: dW_l  ! net water tendency in layer [mm]
  real, intent(in), dimension(:) :: div   ! horizontal water flux divergence [mm/s]
  real, intent(in)               :: delta_time ! land model timestep [s]
  real, intent(in) :: &
     tflow, & ! temperature of surface downwards flow [K]
     d_GW
  ! ---- local vars ----------------------------------------------------------
  real, dimension(num_l)   :: u_minus, u_plus, del_t
  real, dimension(num_l-1) :: eee, fff
  real hcap, aaa, bbb, ccc
  integer l
  ! For energy conservation
  real :: esum1, esum2 ! [W/m^2] heat content of soil before and after solution
  real, parameter :: ethresh = 1.e-4 ! [W/m^2] Allowable error in energy solution for roundoff

!  if (do_component_balchecks .and. .not. LM2) then

     ! Sum energy content in soil before solution
     esum1 = clw*max(flow(1), 0.)*(tflow-tfreeze) ! initialize to incoming surface energy tendency
     do l = 1, num_l
        esum1 = esum1 + (soil%heat_capacity_dry(l)*dz(l) + csw*soil%ws(l) + &
                           clw*(soil%wl(l) - dW_l(l)) ) * (soil%T(l)-tfreeze)
                           ! use water content before Richards eq. update to be consistent with
                           ! heat advection solution
     end do

!  end if

! Upstream weighting of advection. Preserving u_plus here for now.
  u_minus = 1.
  where (flow(1:num_l).lt.0.) u_minus = 0.
  do l = 1, num_l-1
     u_plus(l) = 1. - u_minus(l+1)
  enddo
  hcap = (soil%heat_capacity_dry(num_l)*dz(num_l) &
                              + csw*soil%ws(num_l))/clw
  aaa = -flow(num_l) * u_minus(num_l)
  bbb =  hcap + soil%wl(num_l) - dW_l(num_l) - aaa
  eee(num_l-1) = -aaa/bbb
  fff(num_l-1) = aaa*(soil%T(num_l)-soil%T(num_l-1)) / bbb

  do l = num_l-1, 2, -1
    hcap = (soil%heat_capacity_dry(l)*dz(l) &
                              + csw*soil%ws(l))/clw
    aaa = -flow(l)   * u_minus(l)
    ccc =  flow(l+1) * u_plus (l)
    bbb =  hcap + soil%wl(l) - dW_l(l) - aaa - ccc
    eee(l-1) = -aaa / ( bbb +ccc*eee(l) )
    fff(l-1) = (   aaa*(soil%T(l)-soil%T(l-1))    &
                       + ccc*(soil%T(l)-soil%T(l+1))    &
                       - ccc*fff(l) ) / ( bbb +ccc*eee(l) )
  enddo

  hcap = (soil%heat_capacity_dry(1)*dz(1) + csw*soil%ws(1))/clw
  aaa = -flow(1) * u_minus(1)
  ccc =  flow(2) * u_plus (1)
  bbb =  hcap + soil%wl(1) - dW_l(1) - aaa - ccc

  del_t(1) =  (  aaa*(soil%T(1)-tflow          ) &
                     + ccc*(soil%T(1)-soil%T(2)) &
                     - ccc*fff(1) ) / (bbb+ccc*eee(1))
  soil%T(1) = soil%T(1) + del_t(1)

  if(is_watch_point()) then
     write(*,*) ' ***** soil_step_2 checkpoint 3.4.1 ***** '
     write(*,*) 'hcap', hcap
     write(*,*) 'aaa', aaa
     write(*,*) 'bbb', bbb
     write(*,*) 'ccc', ccc
     write(*,*) 'del_t(1)', del_t(1)
     write(*,*) ' T(1)', soil%T(1)
  endif

  do l = 1, num_l-1
     del_t(l+1) = eee(l)*del_t(l) + fff(l)
     soil%T(l+1) = soil%T(l+1) + del_t(l+1)
  enddo

!  if (do_component_balchecks .and. .not. LM2) then
     ! Sum energy content in soil after solution
     esum2 = -clw * min(flow(1), 0.) * (soil%T(1) - tfreeze)
     if (flow(1) < 0. .and. is_watch_point()) then
        write(*,*) 'Apparent reflux in advection'
        write(*,*) esum2, ' J/m^2 energy refluxed.'
     end if
     do l = 1, num_l
        esum2 = esum2 + (soil%heat_capacity_dry(l)*dz(l) + csw*soil%ws(l) + clw*soil%wl(l)) * &
                      (soil%T(l)-tfreeze) ! no phase change here
        ! Add energy lost to stream via divergence
        esum2 = esum2 + clw*div(l) * delta_time * (soil%T(l)-tfreeze)
     end do

     call check_conservation(module_name // ': Advection subroutine.', 'Energy', esum1, esum2, ethresh, &
                             FATAL)

!  end if

  ! (lumped=lm2 groundwater stored in l=1 prog variable, liquid only)
  if (soil%groundwater(1).ne. 0.) soil%groundwater_T(1) =    &
       + ((aquifer_heat_cap+soil%groundwater(1)-d_GW)  &
                                 *soil%groundwater_T(1) &
        + flow(num_l+1)*soil%T(num_l)) &
         /((aquifer_heat_cap+soil%groundwater(1)-d_GW) + flow(num_l+1))

end subroutine advection


! ============================================================================
! Thermal solution for advection of heat by soil water flow,
! using a call to the generic tridiagonal solver.
! Tested in diagnostic mode; required for gw_option == GW_TILED.
subroutine advection_tri(soil, flow, dW_l, tflow, d_GW, div, delta_time, t_soil_tridiag, hdiv_it)
   type(soil_tile_type), intent(inout) :: soil
   real, intent(in), dimension(:) :: flow  ! water tendency downwards into layer [mm]
   real, intent(in), dimension(:) :: dW_l  ! net water tendency in layer [mm]
   real, intent(in), dimension(:) :: div   ! horizontal water flux divergence [mm/s]
                                           ! used with simple or integrated groundwater model
   real, intent(in)               :: delta_time ! land model timestep [s]
   real, intent(out), dimension(:) :: t_soil_tridiag ! soil temperature solution [K]
   real, intent(in) :: &
      tflow, & ! temperature of surface downwards flow [K]
      d_GW     ! groundwater tendency (for LM2) [mm ??]
   real, intent(in), dimension(:) :: hdiv_it ! divergence of heat due to inter-tile water flow [W/m^2]
   ! ---- local vars ----------------------------------------------------------
   real, dimension(num_l)   :: u_minus, &! coefficient for flow into layer from above [-]
                              u_plus, & ! coefficient for flow into layer from below [-] (=1-u_minus for layer below)
                              del_t, &  ! layer T tendency [K]
                              aaa, bbb, ccc, ddd  ! Tridiagonal matrix coefficients on del_t
   real :: hcapr ! heat capacity ratio of layer, c_solid / c_l * Delta z [m]
   integer :: l ! layer index
   ! For energy conservation
   real :: esum1, esum2 ! [W/m^2] heat content of soil before and after solution
   real, parameter :: ethresh = 1.e-4 ! [W/m^2] Allowable error in energy solution for roundoff

!   if (do_component_balchecks) then
      esum1 = clw*max(flow(1), 0.)*(tflow-tfreeze) ! initialize to incoming surface energy tendency
      do l = 1, num_l
         esum1 = esum1 + (soil%heat_capacity_dry(l)*dz(l) + csw*soil%ws(l) + &
                           clw*(soil%wl(l) - dW_l(l)) ) * (soil%T(l)-tfreeze)
                           ! use water content before Richards eq. update to be consistent with
                           ! heat advection solution
         if (is_watch_point()) write(*,*)'l, esum1:', l, esum1
      end do
!   end if

   ! Upstream weighting of advection. Preserving u_plus here for now.
   u_minus = 1.0; u_plus = 0.0
   where (flow(1:num_l).lt.0.) u_minus = 0.
   do l = 1, num_l-1
      u_plus(l) = 1. - u_minus(l+1)
   enddo

   ! Set Tridiagonal coefficients
   ! Solving equations aaa(l)*delta_t(l-1) + bbb(l)*delta_t(l) + ccc(l)*delta_t(l+1) = ddd(l)
   do l = 1, num_l
      hcapr = (soil%heat_capacity_dry(l)*dz(l) &
                              + csw*soil%ws(l))/clw
      if (l == 1) then
         aaa(l) = 0. ! T_0 is prescribed to tflow
      else
         aaa(l) = - u_minus(l) * flow(l)
      end if

      bbb(l) = hcapr + soil%wl(l) - (1.-u_minus(l)) * flow(l) + &
            (1.-u_plus(l)) * flow(l+1) + delta_time * div(l)

      if (l < num_l) then
         ccc(l) = u_plus(l) * flow(l+1)
      else
         ccc(l) = 0. ! zero flux bottom boundary; flow(l+1) = 0.
      end if

      if (l == 1) then
         ddd(l) = u_minus(l) * flow(l) * tflow  &
                  + soil%T(l) * ( (1.-u_minus(l)) * flow(l)  &
                     - (1-u_plus(l)) * flow(l+1) - dW_l(l) - delta_time * div(l) ) &
                  - u_plus(l) * flow(l+1) * soil%T(l+1)
      else if (l < num_l) then
         ddd(l) = u_minus(l) * flow(l) * soil%T(l-1)  &
                  + soil%T(l) * ( (1.-u_minus(l)) * flow(l)  &
                     - (1-u_plus(l)) * flow(l+1) - dW_l(l) - delta_time * div(l) ) &
                  - u_plus(l) * flow(l+1) * soil%T(l+1)
      else ! l == num_l
         ddd(l) = u_minus(l) * flow(l) * soil%T(l-1)  &
                  + soil%T(l) * ( (1.-u_minus(l)) * flow(l)  &
                      - dW_l(l) - delta_time * div(l) )
      end if

      ! update coefficients if full tiled hillslope model is used
      if (gw_option == GW_TILED) then
         ! Remove div terms. Energy flux will be prescribed explicitly based on inter-tile flows
         ! calculated in hlsp_hydrology_1.
         bbb(l) = bbb(l) - delta_time * div(l)
         ddd(l) = ddd(l) + soil%T(l) * delta_time * div(l)
         ! Add explicit energy sink term to ddd(l)
         ddd(l) = ddd(l) - delta_time / clw * hdiv_it(l) - delta_time * div(l) * tfreeze
         ! kg/m^2 K          s        /(J/kg/K) * J/m^2/s      s     *  kg/m^2/s   *K
      end if
   end do

   ! Update temperature
   call tridiag(aaa, bbb, ccc, ddd, del_t)
   t_soil_tridiag(1:num_l) = soil%T(1:num_l) + del_t(1:num_l)

   if (use_tridiag_foradvec) then
      soil%T(1:num_l) = t_soil_tridiag(1:num_l)
   end if

   if(is_watch_point()) then
      write(*,*) ' ***** soil_step_2 checkpoint 3.4.1 ***** '
      write(*,*) 'hcapr(num_l)', hcapr
      write(*,*) 'aaa', aaa(:)
      write(*,*) 'bbb', bbb(:)
      write(*,*) 'ccc', ccc(:)
      write(*,*) 'ddd', ddd(:)
      write(*,*) ' T', soil%T(:)
      write(*,*) 'hdiv_it', hdiv_it(:)
   endif

!   if (do_component_balchecks) then
      ! Sum energy content in soil after solution
      esum2 = -clw * min(flow(1), 0.) * (t_soil_tridiag(1) - tfreeze) ! if upwards flow at surface
                                      ! will be energy lost to sat runoff ?
      if (flow(1) < 0. .and. is_watch_point()) then
         write(*,*) 'Apparent reflux in advection'
         write(*,*) esum2, ' J/m^2 energy refluxed.'
      end if
      do l = 1, num_l
         esum2 = esum2 + (soil%heat_capacity_dry(l)*dz(l) + csw*soil%ws(l) + clw*soil%wl(l)) * &
                        (t_soil_tridiag(l)-tfreeze) ! no phase change here
         ! Add energy lost via divergence
         if (gw_option /= GW_TILED) then
            esum2 = esum2 + clw*div(l) * delta_time * (t_soil_tridiag(l)-tfreeze)
         else
            esum2 = esum2 + hdiv_it(l) * delta_time
         end if
         if (is_watch_point()) write(*,*)'l, esum2:', l, esum2
      end do

      call check_conservation(module_name // ': Advection_tri subroutine.', 'Energy', esum1, esum2, ethresh, &
                              FATAL)

!   end if


   ! (lumped=lm2 groundwater stored in l=1 prog variable, liquid only)
   if (soil%groundwater(1).ne. 0.) soil%groundwater_T(1) =    &
         + ((aquifer_heat_cap+soil%groundwater(1)-d_GW)  &
                            *soil%groundwater_T(1) &
         + flow(num_l+1)*soil%T(num_l)) &
         /((aquifer_heat_cap+soil%groundwater(1)-d_GW) + flow(num_l+1))

end subroutine advection_tri


! ============================================================================
! Spread new root C through profile, using vertical root profile from vegn_uptake_profile
subroutine add_root_litter(soil,vegn,newlitterC)
    type(soil_tile_type), intent(inout)  :: soil
    type(vegn_tile_type), intent(in)     :: vegn
    real,intent(in) :: newlitterC(:)

    real,dimension(num_l) :: uptake_frac_max, vegn_uptake_term
    integer :: nn

    call vegn_uptake_profile (vegn, dz(1:num_l), uptake_frac_max, vegn_uptake_term )
    if(abs(sum(uptake_frac_max(1:num_l)) - 1.0) > 1e-10 ) then
        print *,'1 - sum(vegn_uptake_frac_max)',1.0-sum(uptake_frac_max(1:num_l))
        call error_mesg('add_root_litter','total of vegn_uptake_frac_max not 1',FATAL)
    endif
    do nn=1,num_l
        call add_litter(soil%soil_C(nn),newLitterC*uptake_frac_max(nn))
        soil%fsc_in(nn)=soil%fsc_in(nn)+newLitterC(1)*uptake_frac_max(nn)
        soil%ssc_in(nn)=soil%ssc_in(nn)+newLitterC(2)*uptake_frac_max(nn)
    enddo
end subroutine add_root_litter


! ============================================================================
! Spread root exudate C through profile, using vertical root profile from vegn_uptake_profile
! Differs from add_root_litter -- C is distributed through existing cohorts, not deposited as new cohort
subroutine add_root_exudates(soil,vegn,exudateC)
  type(soil_tile_type), intent(inout)  :: soil
  type(vegn_tile_type), intent(in)     :: vegn
  real,intent(in) :: exudateC

  real,dimension(num_l) :: uptake_frac_max, vegn_uptake_term
  integer :: nn
  real, parameter :: tol = 1e-10 ! tolerance of profile integral test.

  select case (soil_carbon_option)
  case (SOILC_CENTURY)
     soil%fast_soil_C(1) = soil%fast_soil_C(1) + exudateC

  case (SOILC_CENTURY_BY_LAYER)
     call vegn_uptake_profile (vegn, dz(1:num_l), uptake_frac_max, vegn_uptake_term )
     call check_var_range(sum(uptake_frac_max)-1.0, -tol,+tol, 'add_root_exudates', 'sum(vegn_uptake_frac_max)-1', FATAL)

     do nn=1,num_l
        soil%fast_soil_C(nn) = soil%fast_soil_C(nn) + exudateC*uptake_frac_max(nn)
     enddo

  case (SOILC_CORPSE)
     call vegn_uptake_profile (vegn, dz(1:num_l), uptake_frac_max, vegn_uptake_term )
     call check_var_range(sum(uptake_frac_max)-1.0, -tol,+tol, 'add_root_exudates', 'sum(vegn_uptake_frac_max)-1', FATAL)

     if (is_watch_point()) then
        write (*,*) '##### add_root_exudates #####'
        __DEBUG1__(exudateC)
        __DEBUG1__(uptake_frac_max)
     endif

     do nn=1,num_l
        if (is_watch_point()) then
           call debug_pool(soil%soil_C(nn),'soil_C(nn) before')
        endif
        call add_carbon_to_cohorts(soil%soil_C(nn),litterC=(/exudateC*uptake_frac_max(nn),0.0,0.0/))
        soil%fsc_in(nn)=soil%fsc_in(nn)+exudateC*uptake_frac_max(nn)
        if (is_watch_point()) then
           call debug_pool(soil%soil_C(nn),'soil_C(nn) after ')
        endif
     enddo
  end select
end subroutine add_root_exudates


! ============================================================================
subroutine redistribute_peat_carbon(soil)
    type(soil_tile_type), intent(inout) :: soil

    integer :: nn
    real :: layer_total_C,layer_total_C_2,layer_max_C,layer_extra_C,fraction_to_remove
    real :: total_C_before,total_C_after
    real :: leaflitter_total_C, woodlitter_total_C

    !For conservation check.
    total_C_before=0.0
    do nn=1,num_l
    call poolTotalCarbon(soil%soil_C(num_l),layer_total_C)
    total_C_before=total_C_before+layer_total_C
    enddo

    call poolTotalCarbon(soil%leaflitter,totalCarbon=leaflitter_total_C)
    call poolTotalCarbon(soil%coarseWoodLitter,totalCarbon=woodlitter_total_C)
    layer_total_C=leaflitter_total_C+woodlitter_total_C

    layer_max_C=max_litter_thickness*max_soil_C_density
    layer_extra_C = layer_total_C-layer_max_C
    if(layer_extra_C>0) then
        fraction_to_remove=1.0-layer_max_C/layer_total_C
        call transfer_pool_fraction(soil%leaflitter,soil%soil_C(1),fraction_to_remove)
        call transfer_pool_fraction(soil%coarsewoodlitter,soil%soil_C(1),fraction_to_remove)
    endif

    !Move carbon down if it exceeds layer_max_C
    do nn=1,num_l-1
        call poolTotalCarbon(soil%soil_C(nn),totalCarbon=layer_total_C)
        layer_max_C=dz(nn)*max_soil_C_density
        layer_extra_C=layer_total_C-layer_max_C
        if (layer_extra_C>0) then
            fraction_to_remove=1.0-layer_max_C/layer_total_C
            call transfer_pool_fraction(soil%soil_C(nn),soil%soil_C(nn+1),fraction_to_remove)
            soil%is_peat(nn)=1
        endif

        if (layer_extra_C < 0 .and. (soil%is_peat(nn).ne.0) .and. (soil%is_peat(nn+1).ne.0)) then
             call poolTotalCarbon(soil%soil_C(nn+1),totalCarbon=layer_total_C_2)
             fraction_to_remove = -layer_extra_C/layer_total_C_2
             if (fraction_to_remove > 0.5) then
                soil%is_peat(nn+1)=0
             else
                call transfer_pool_fraction(soil%soil_C(nn+1),soil%soil_C(nn),fraction_to_remove)
             endif
        endif
    enddo

    total_C_after=0.0
    do nn=1,num_l
    call poolTotalCarbon(soil%soil_C(num_l),layer_total_C)
    total_C_after=total_C_after+layer_total_C
    enddo

    if (abs(total_C_before-total_C_after)>1e-10) then
            print *,'Carbon before:',total_C_before
            print *,'Carbon after:',total_C_after
            call error_mesg('redistribute_peat_carbon','Carbon not conserved after downward move',FATAL)
    endif
end subroutine redistribute_peat_carbon

! ============================================================================
! tile existence detector: returns a logical value indicating whether component
! model tile exists or not
logical function soil_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   soil_tile_exists = associated(tile%soil)
end function soil_tile_exists


!! ZMS Functions moved here from soil_tile_mod to avoid circular dependencies
!! with hillslope_mod and hillslope_tile_mod.
!============================================================================
function soil_cover_cold_start(land_mask, lonb, latb) result(soil_frac)
! creates and initializes a field of fractional soil coverage
  logical, intent(in) :: land_mask(:)    ! land mask
  real,    intent(in) :: lonb(:,:), latb(:,:) ! boundaries of the grid cells
  real,    pointer :: soil_frac (:,:) ! output: map of soil fractional coverage
  integer :: k

  if (.not. do_hillslope_model) then
     allocate( soil_frac(size(land_mask(:)),n_dim_soil_types))
     allocate( soil_tags(size(land_mask(:)),n_dim_soil_types))
     do k = 1, n_dim_soil_types
        soil_tags(:, k) = k
     end do
  else
     allocate( soil_frac(size(land_mask(:)), &
               n_dim_soil_types * num_vertclusters * max_num_topo_hlsps))
     allocate( soil_tags(size(land_mask(:)), &
               n_dim_soil_types * num_vertclusters * max_num_topo_hlsps))
     do k = 1,size(soil_tags,2)
        if (mod(k, n_dim_soil_types) == 0) then
           soil_tags(:,k) = n_dim_soil_types
        else
           soil_tags(:,k) = mod(k, n_dim_soil_types)
        end if
     end do
  end if

  call init_cover_field(soil_to_use, trim(soil_type_file), 'cover','frac', &
       lonb, latb, soil_index_constant, input_cover_types, soil_frac)

  if (do_hillslope_model) then
     call hlsp_coldfracs(soil_frac, n_dim_soil_types)
  end if

end function soil_cover_cold_start


! ============================================================================
subroutine retrieve_soil_tags(soiltags)
! passes back field of soil cover types from soil_cover_cold_start
  integer, pointer :: soiltags(:,:)

  if (size(soiltags,1) /= size(soil_tags,1) .or. size(soiltags,2) /= size(soil_tags,2)) &
     call error_mesg(module_name,'Wrong dimension size in "soil_tile_mod:retrieve_soil_tags.',FATAL)

  soiltags(:,:) = soil_tags(:,:)

  deallocate(soil_tags)
end subroutine retrieve_soil_tags

! ============================================================================
! Calculates soil temperature profile for soil tile, if use_coldstart_wtt_data
subroutine init_soil_twc(soil, ref_soil_t, mwc)
   type(soil_tile_type), intent(inout) :: soil
   real, intent(in)   :: ref_soil_t ! local reference temperature from input data [K]
   real, intent(in), dimension(num_l)  :: mwc  ! total water content to init [kg/m^3]
   integer :: l ! soil level
   real, dimension(num_l) :: vlc, vsc ! volumetric liquid and solid water contents [-]
   real, dimension(num_l)  :: thermal_cond ! soil thermal conductivity [W/m/K]
   real  :: tres       ! thermal resistance [K / (W/m^2)]

   ! First tentatively initialize soil water / ice content, assuming isothermal profile.
   if (ref_soil_t.ge.soil%pars%tfreeze) then
      soil%wl(1:num_l) = mwc(1:num_l)*dz(1:num_l)
      soil%ws(1:num_l) = 0.
   else
      soil%wl(1:num_l) = 0.
      soil%ws(1:num_l) = mwc(1:num_l)*dz(1:num_l)
   endif
   ! Initialize T and groundwater
   soil%T             = ref_soil_t
   soil%groundwater   = init_groundwater
   soil%groundwater_T = ref_soil_t
   soil%uptake_T           = ref_soil_t
   ! Initialize thermal conductivity
   vlc(:) = soil%wl(1:num_l)/ (dz(1:num_l)*dens_h2o*soil%pars%vwc_sat)
   vsc(:) = soil%ws(1:num_l)/ (dz(1:num_l)*dens_h2o*soil%pars%vwc_sat)
   call soil_data_thermodynamics ( soil, vlc, vsc, thermal_cond)

   ! Walk down through soil and maintain geothermal heat flux profile.
   do l=2,num_l
      tres = 0.5*(dz(l-1)/thermal_cond(l-1) + dz(l)/thermal_cond(l)) ! [K / (W/m^2)] =  [m / (W/m/K)]
      soil%T(l) = soil%T(l-1) + soil%geothermal_heat_flux * tres

      ! Check for switch to above-freezing
      if (soil%T(l) >= soil%pars%tfreeze .and. soil%T(l-1) < soil%pars%tfreeze) then
         ! Reset wl & ws, and recalculate thermal conductivity, from here down.
         soil%wl(l:num_l) = mwc(l:num_l)*dz(l:num_l)
         soil%ws(l:num_l) = 0.
         vlc(l:num_l) = soil%wl(l:num_l)/ (dz(l:num_l)*dens_h2o*soil%pars%vwc_sat)
         vsc(l:num_l) = soil%ws(l:num_l)/ (dz(l:num_l)*dens_h2o*soil%pars%vwc_sat)
         call soil_data_thermodynamics ( soil, vlc, vsc, thermal_cond)
      end if
   end do

   ! Debug
   ! current point set above call in soil_init
   if (is_watch_point()) then
      write(*,*)'Initializing soil temperature and water for watch_point.'
      __DEBUG1__(ref_soil_t)
      do l=1,num_l
         __DEBUG4__(l, vlc(l), vsc(l), soil%T(l))
      end do
   end if


end subroutine init_soil_twc


! ============================================================================
! cohort accessor functions: given a pointer to cohort, return a pointer to a
! specific member of the cohort structure
#define DEFINE_SOIL_ACCESSOR_0D(xtype,x) subroutine soil_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%x;endif;\
end subroutine
#define DEFINE_SOIL_ACCESSOR_1D(xtype,x) subroutine soil_ ## x ## _ptr(t,i,p);\
type(land_tile_type),pointer::t;integer,intent(in)::i;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%x(i);endif;\
end subroutine
#define DEFINE_SOIL_COMPONENT_ACCESSOR_0D(xtype,component,x) subroutine soil_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%component%x;endif;\
end subroutine
#define DEFINE_SOIL_COMPONENT_ACCESSOR_1D(xtype,component,x) subroutine soil_ ## x ## _ptr(t,i,p);\
type(land_tile_type),pointer::t;integer,intent(in)::i;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%component%x(i);endif;\
end subroutine

DEFINE_SOIL_ACCESSOR_1D(real,w_fc)
DEFINE_SOIL_ACCESSOR_1D(real,alpha)
DEFINE_SOIL_ACCESSOR_0D(real,uptake_T)
DEFINE_SOIL_ACCESSOR_0D(integer,tag)
DEFINE_SOIL_ACCESSOR_1D(real,fast_soil_C)
DEFINE_SOIL_ACCESSOR_1D(real,slow_soil_C)
DEFINE_SOIL_ACCESSOR_1D(real,T)
DEFINE_SOIL_ACCESSOR_1D(real,wl)
DEFINE_SOIL_ACCESSOR_1D(real,ws)
DEFINE_SOIL_ACCESSOR_1D(real,groundwater)
DEFINE_SOIL_ACCESSOR_1D(real,groundwater_T)
DEFINE_SOIL_ACCESSOR_1D(integer,is_peat)

DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,tau_groundwater)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_length)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_relief)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_a)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_n)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,hillslope_zeta_bar)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,soil_e_depth)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,zeta)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,tau)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,k_sat_gw)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,vwc_wilt)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,vwc_fc)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,vwc_sat)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,k_sat_ref)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars,Qmax)

DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_dry_dir)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_dry_dif)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_sat_dir)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,refl_sat_dif)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_iso_dry)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_vol_dry)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_geo_dry)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_iso_sat)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_vol_sat)
DEFINE_SOIL_COMPONENT_ACCESSOR_1D(real,pars,f_geo_sat)

DEFINE_SOIL_ACCESSOR_0D(real,fast_DOC_leached)
DEFINE_SOIL_ACCESSOR_0D(real,slow_DOC_leached)
DEFINE_SOIL_ACCESSOR_0D(real,deadmic_DOC_leached)

DEFINE_SOIL_ACCESSOR_1D(real,asoil_in)
DEFINE_SOIL_ACCESSOR_1D(real,fsc_in)
DEFINE_SOIL_ACCESSOR_1D(real,ssc_in)
DEFINE_SOIL_ACCESSOR_1D(real,deadmic_in)
DEFINE_SOIL_ACCESSOR_1D(real,fast_protected_in)
DEFINE_SOIL_ACCESSOR_1D(real,slow_protected_in)
DEFINE_SOIL_ACCESSOR_1D(real,deadmic_protected_in)
DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_fsc_in)
DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_ssc_in)
DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_deadmic_in)
DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_fsc_in)
DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_ssc_in)
DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_deadmic_in)
DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_fsc_in)
DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_ssc_in)
DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_deadmic_in)
DEFINE_SOIL_ACCESSOR_1D(real,fast_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_1D(real,slow_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_1D(real,deadmic_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_1D(real,fast_protected_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_1D(real,slow_protected_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_1D(real,deadmic_protected_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_fast_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_slow_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,leaflitter_deadmic_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_fast_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_slow_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,finewoodlitter_deadmic_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_fast_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_slow_turnover_accumulated)
DEFINE_SOIL_ACCESSOR_0D(real,coarsewoodlitter_deadmic_turnover_accumulated)

! stuff below is for CORPSE

#define DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR1(xtype,x) subroutine soilc_ ## x ## _ptr(t,i,j,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;integer,intent(in)::i,j;p=>NULL();if(associated(t))then;\
if(associated(t%soil))p=>t%soil%soil_C(i)%litterCohorts(j)%x;endif;\
endsubroutine
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR1(real,livingMicrobeC)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR1(real,Rtot)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR1(real,CO2)
DEFINE_SOIL_LAYER_COHORT_COMPONENT_ACCESSOR1(real,originalLitterC)

subroutine sc_soil_C_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_C(i)%litterCohorts(j)%litterC(k)
  endif
end subroutine

subroutine sc_protected_C_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%soil_C(i)%litterCohorts(j)%protectedC(k)
  endif
end subroutine

#define DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR0(xtype,pool,x) subroutine soilc_ ## pool ## _ ## x ## _ptr(t,i,p);\
type(land_tile_type),pointer::t;integer,intent(in)::i;xtype,pointer::p;p=>NULL();if(associated(t))then;\
if(associated(t%soil))p=>t%soil%pool%litterCohorts(i)%x;endif;\
endsubroutine

DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR0(real,leafLitter,livingMicrobeC)
DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR0(real,leafLitter,CO2)
DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR0(real,leafLitter,Rtot)
DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR0(real,leafLitter,originalLitterC)

DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR0(real,fineWoodLitter,livingMicrobeC)
DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR0(real,fineWoodLitter,CO2)
DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR0(real,fineWoodLitter,Rtot)
DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR0(real,fineWoodLitter,originalLitterC)

DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR0(real,coarseWoodLitter,livingMicrobeC)
DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR0(real,coarseWoodLitter,CO2)
DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR0(real,coarseWoodLitter,Rtot)
DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR0(real,coarseWoodLitter,originalLitterC)

#define DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR1(xtype,pool,x) subroutine soilc_ ## pool ## _ ## x ## _ptr(t,i,j,p);\
type(land_tile_type),pointer::t;integer,intent(in)::i,j;xtype,pointer::p;p=>NULL();if(associated(t))then;\
if(associated(t%soil))p=>t%soil%pool%litterCohorts(i)%x(j);endif;\
endsubroutine

DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR1(real,leafLitter,litterC)
DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR1(real,leafLitter,protectedC)

DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR1(real,fineWoodLitter,litterC)
DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR1(real,fineWoodLitter,protectedC)

DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR1(real,coarseWoodLitter,litterC)
DEFINE_SOIL_C_POOL_COHORT_COMPONENT_ACCESSOR1(real,coarseWoodLitter,protectedC)

#define DEFINE_SOIL_C_POOL_DOC_ACCESSOR(xtype,pool) subroutine soilc_ ## pool ## _DOC_ptr(t,i,p);\
type(land_tile_type),pointer::t;integer,intent(in)::i;real,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%pool%dissolved_carbon(i);endif;\
endsubroutine

DEFINE_SOIL_C_POOL_DOC_ACCESSOR(real,leafLitter)
DEFINE_SOIL_C_POOL_DOC_ACCESSOR(real,fineWoodLitter)
DEFINE_SOIL_C_POOL_DOC_ACCESSOR(real,coarseWoodLitter)

subroutine soil_fast_DOC_ptr(t,i,p)
   type(land_tile_type), pointer    :: t
   integer,              intent(in) :: i
   real,                 pointer    :: p

   p=>NULL()
   if(associated(t))then
      if(associated(t%soil))p=>t%soil%soil_C(i)%dissolved_carbon(1)
   endif
end subroutine

subroutine soil_slow_DOC_ptr(t,i,p)
   type(land_tile_type), pointer    :: t
   integer,              intent(in) :: i
   real,                 pointer    :: p

   p=>NULL()
   if(associated(t))then
      if(associated(t%soil))p=>t%soil%soil_C(i)%dissolved_carbon(2)
   endif
end subroutine

subroutine soil_deadMicrobe_DOC_ptr(t,i,p)
   type(land_tile_type), pointer    :: t
   integer,              intent(in) :: i
   real,                 pointer    :: p

   p=>NULL()
   if(associated(t))then
      if(associated(t%soil))p=>t%soil%soil_C(i)%dissolved_carbon(3)
   endif
end subroutine

end module soil_mod
