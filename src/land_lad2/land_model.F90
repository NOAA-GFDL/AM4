! ============================================================================
! top-level core of the Land Dynamics (LaD) model code
! ============================================================================
module land_model_mod

#include "shared/debug.inc"

use time_manager_mod, only : time_type, get_time, increment_time, time_type_to_real, &
     operator(+)
use mpp_domains_mod, only : domain2d, mpp_get_ntile_count, mpp_pass_SG_to_UG, mpp_pass_ug_to_sg, &
                            mpp_get_UG_domain_tile_pe_inf, mpp_get_UG_domain_ntiles

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use mpp_mod, only : mpp_max, mpp_sum, mpp_chksum, MPP_FILL_INT, MPP_FILL_DOUBLE
use fms_io_mod, only : restart_file_type, free_restart_type
use fms_mod, only : error_mesg, FATAL, WARNING, NOTE, mpp_pe, &
     mpp_root_pe, file_exist, check_nml_error, close_file, &
     stdlog, stderr, mpp_clock_id, mpp_clock_begin, mpp_clock_end, string, &
     stdout, CLOCK_FLAG_DEFAULT, CLOCK_COMPONENT, CLOCK_ROUTINE
use field_manager_mod, only : MODEL_LAND, MODEL_ATMOS
use data_override_mod, only : data_override_ug
use diag_manager_mod, only : diag_axis_init, register_static_field, &
     register_diag_field, send_data, diag_field_add_attribute
use constants_mod, only : radius, hlf, hlv, hls, tfreeze, pi, rdgas, rvgas, cp_air, &
     stefan
use astronomy_mod, only : astronomy_init, diurnal_solar
use tracer_manager_mod, only : NO_TRACER, get_tracer_index, get_tracer_names
use sphum_mod, only : qscomp

use land_constants_mod, only : NBANDS, BAND_VIS, BAND_NIR, d608, mol_air, mol_C, mol_co2
use land_tracers_mod, only : land_tracers_init, land_tracers_end, ntcana, isphum, ico2
use land_tracer_driver_mod, only: land_tracer_driver_init, land_tracer_driver_end, &
     update_cana_tracers
use glacier_mod, only : read_glac_namelist, glac_init, glac_end, &
     glac_step_1, glac_step_2, save_glac_restart, glac_sfc_water
use lake_mod, only : read_lake_namelist, lake_init, lake_end, &
     lake_sfc_water, lake_step_1, lake_step_2, save_lake_restart, &
     lake_rh_feedback, LAKE_RH_BETA
use soil_mod, only : read_soil_namelist, soil_init, soil_end, &
     soil_sfc_water, soil_evap_limits, soil_step_1, soil_step_2, soil_step_3, save_soil_restart
use soil_carbon_mod, only : read_soil_carbon_namelist, n_c_types
use lake_mod, only : lake_init_predefined
use snow_mod, only : read_snow_namelist, snow_init, snow_end, &
     snow_get_depth_area, snow_step_1, snow_step_2, &
     save_snow_restart, sweep_tiny_snow, snow_sfc_water
use vegetation_mod, only : read_vegn_namelist, vegn_init, vegn_end, vegn_get_cover, &
     vegn_radiation, vegn_properties, vegn_step_1, vegn_step_2, vegn_step_3, &
     update_vegn_slow, save_vegn_restart
use cana_tile_mod, only : canopy_air_mass, canopy_air_mass_for_tracers, cana_tile_heat
use canopy_air_mod, only : read_cana_namelist, cana_init, cana_end, cana_state,&
     cana_step_2, cana_roughness, &
     save_cana_restart
use river_mod, only : river_init, river_end, update_river, river_stock_pe, &
     save_river_restart, river_tracers_init, num_river_tracers, river_tracer_index
use topo_rough_mod, only : topo_rough_init, topo_rough_end, update_topo_rough
use soil_tile_mod, only : soil_tile_stock_pe, soil_tile_heat, soil_roughness, &
     soil_radiation
use soil_mod, only      : soil_cover_cold_start, retrieve_soil_tags ! moved here
                          ! to eliminate circular dependencies with hillslope mods
use vegn_tile_mod, only : vegn_cover_cold_start, vegn_data_rs_min, &
     update_derived_vegn_data, vegn_tile_stock_pe, vegn_tile_heat
use lake_tile_mod, only : lake_cover_cold_start, lake_tile_stock_pe, &
     lake_tile_heat, lake_radiation, lake_roughness
use glac_tile_mod, only : glac_cover_cold_start, glac_tile_stock_pe, &
     glac_tile_heat, glac_radiation, glac_roughness
use snow_tile_mod, only : snow_tile_stock_pe, snow_tile_heat, snow_roughness, &
     snow_radiation, snow_active
use land_numerics_mod, only : ludcmp, lubksb, &
     horiz_remap_type, horiz_remap_new, horiz_remap, horiz_remap_del, &
     horiz_remap_print
use land_io_mod, only : read_land_io_namelist, input_buf_size, new_land_io
use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_list_type, &
     land_tile_enum_type, new_land_tile, insert, nitems, &
     first_elmt, get_tile_tags, land_tile_carbon, land_tile_heat, &
     get_tile_water, init_tile_map, free_tile_map, max_n_tiles, &
     tile_exists_func, loop_over_tiles
use land_data_mod, only : land_data_type, atmos_land_boundary_type, &
     land_state_type, land_state_type_sg, land_data_init, land_data_end, lnd, lnd_sg, &
     log_version
use nf_utils_mod,  only : nfu_inq_var, nfu_inq_dim, nfu_get_var
use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_tile_data, add_int_tile_data, get_tile_data, &
     field_exists, print_netcdf_error
use land_tile_diag_mod, only : tile_diag_init, tile_diag_end, &
     set_default_diag_filter, get_area_id, &
     register_tiled_diag_field, register_tiled_area_fields, &
     add_tiled_diag_field_alias, &
     send_tile_data, dump_tile_diag_fields, &
     OP_AVERAGE, OP_SUM, cmor_name, send_global_land_diag
use land_debug_mod, only : land_debug_init, land_debug_end, set_current_point, &
     is_watch_point, get_watch_point, check_temp_range, current_face, &
     get_current_point, check_conservation, water_cons_tol, carbon_cons_tol, &
     is_watch_cell, is_watch_time, do_check_conservation, check_var_range, log_date
use static_vegn_mod, only : write_static_vegn
use land_transitions_mod, only : &
     land_transitions_init, land_transitions_end, land_transitions, &
     save_land_transitions_restart
use stock_constants_mod, only: ISTOCK_WATER, ISTOCK_HEAT, ISTOCK_SALT
use hillslope_mod, only: retrieve_hlsp_indices, save_hlsp_restart, hlsp_end, &
     read_hlsp_namelist, hlsp_init, hlsp_config_check
use hillslope_hydrology_mod, only: hlsp_hydrology_1, hlsp_hydro_init
use hillslope_mod, only: hlsp_init_predefined
use vegn_data_mod, only : LU_CROP, LU_PAST, LU_NTRL, LU_SCND, LU_URBN, &
    SP_C4GRASS, SP_C3GRASS, SP_TEMPDEC, SP_TROPICAL, SP_EVERGR
use predefined_tiles_mod, only: land_cover_cold_start_0d_predefined_tiles,&
                                open_database_predefined_tiles,&
                                close_database_predefined_tiles

use fms_io_mod,      only: fms_io_unstructured_read
use mpp_domains_mod, only: domainUG
use mpp_domains_mod, only: mpp_get_UG_compute_domain
use mpp_domains_mod, only: mpp_get_UG_domain_grid_index
use diag_axis_mod,   only: diag_axis_add_attribute

implicit none
private

! ==== public interfaces =====================================================
public land_model_init          ! initialize the land model
public land_model_end           ! finish land model calculations
public land_model_restart       ! saves the land model restart(s)
public update_land_model_fast   ! time-step integration
public update_land_model_slow   ! time-step integration
public atmos_land_boundary_type ! data from coupler to land
public land_data_type           ! data from land to coupler
public land_data_type_chksum    ! routine to print checksums for land_data_type
public atm_lnd_bnd_type_chksum  ! routine to print checksums for atmos_land_boundary_type

public :: Lnd_stock_pe          ! return stocks of conservative quantities

! re-export land diagnostic subroutines for tiled diag in flux exchange
public set_default_diag_filter, register_tiled_diag_field, send_tile_data
public send_global_land_diag
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'land'
#include "shared/version_variable.inc"

! ==== module variables ======================================================

! ---- namelist --------------------------------------------------------------
logical :: use_old_conservation_equations  = .false.
logical :: lm2                             = .false.
logical :: give_stock_details              = .false.
logical :: use_tfreeze_in_grnd_latent      = .false.
logical :: use_atmos_T_for_precip_T        = .false.
logical :: use_atmos_T_for_evap_T          = .false.
real    :: cpw = 1952.  ! specific heat of water vapor at constant pressure
real    :: clw = 4218.  ! specific heat of water (liquid)
real    :: csw = 2106.  ! specific heat of water (ice)
real    :: min_sum_lake_frac = 1.e-8
real    :: min_frac = 0.0 ! minimum fraction of soil, lake, and glacier that is not discarded on cold start
real    :: gfrac_tol         = 1.e-6
real    :: discharge_tol = -1.e20
real    :: con_fac_large = 1.e6
real    :: con_fac_small = 1.e-6
real    :: tau_snow_T_adj = -1.0 ! time scale of snow temperature adjustment
              ! for the snow-free surface (s); negative means no adjustment
logical :: prohibit_negative_canopy_water = .FALSE. ! if true, then in case of negative canopy
              ! water the evaporation is fixed and the equations are re-solved.
              ! Default retrieves old behavior.
character(16) :: nearest_point_search = 'global' ! specifies where to look for
              ! nearest points for missing data, "global" or "face"
logical :: print_remapping = .FALSE. ! if true, full land cover remapping
              ! information is printed on the cold start
logical :: predefined_tiles= .FALSE. ! If true, the tiles for each grid cell for
              ! each grid cell are read from an external file
integer :: layout(2) = (/0,0/)
integer :: io_layout(2) = (/0,0/)
integer :: npes_io_group = 0
  ! mask_table contains information for masking domain ( n_mask, layout and mask_list).
  !   A text file to specify n_mask, layout and mask_list to reduce number of processor
  !   usage by masking out some domain regions which contain all ocean points.
  !   The default file name of mask_table is "INPUT/land_mask_table". Please note that
  !   the file name must begin with "INPUT/". The first
  !   line of mask_table will be number of region to be masked out. The second line
  !   of the mask_table will be the layout of the model. User need to set land_model_nml
  !   variable layout to be the same as the second line of the mask table.
  !   The following n_mask line will be the position of the processor and face number
  !   to be masked out (First entry is i-direction postion, second entry is
  !   j-direction position, third entry is face number ).
  !   The mask_table could be created by tools check_mask.
  !   For example the mask_table will be as following if n_mask=4, layout=4,6 and
  !   the processor (1,2) and (3,6) will be masked out.
  !     2
  !     4,6
  !     1,2,1
  !     3,6,2
  !     4,5,5
  !     2,4,6

  character(len=128) :: mask_table = "INPUT/land_mask_table"

logical :: use_coast_rough  = .false. ! if true, override roughness on coasts with prescribed value
real    :: coast_rough_mom  = 1.0e-4  ! prescribed coastal roughness for momentum
real    :: coast_rough_heat = 1.35e-5 ! prescribed coastal roughness for heat and tracers
real    :: max_coast_frac   = 1.0     ! threshold defining which point is coastal
logical :: use_coast_topo_rough = .true. ! if false, the topographic roughness scaling
                                      ! is not used over coastal points
real    :: precip_warning_tol = -1.0e-18 ! if liquid or solid precip (input
           ! from atmos) is below this value, a warning is printed

namelist /land_model_nml/ use_old_conservation_equations, &
                          lm2, give_stock_details, &
                          use_tfreeze_in_grnd_latent, &
                          use_atmos_T_for_precip_T, &
                          use_atmos_T_for_evap_T, &
                          cpw, clw, csw, min_sum_lake_frac, min_frac, &
                          gfrac_tol, discharge_tol, &
                          con_fac_large, con_fac_small, &
                          tau_snow_T_adj, prohibit_negative_canopy_water, &
                          nearest_point_search, print_remapping, &
                          use_coast_rough, coast_rough_mom, coast_rough_heat, &
                          max_coast_frac, use_coast_topo_rough, &
                          layout, io_layout, mask_table, &
                          precip_warning_tol, npes_io_group, predefined_tiles
! ---- end of namelist -------------------------------------------------------

logical  :: module_is_initialized = .FALSE.
logical  :: stock_warning_issued  = .FALSE.
logical  :: update_cana_co2 ! if false, cana_co2 is not updated during the model run.
character(len=256) :: grid_spec_file="INPUT/grid_spec.nc"
real     :: delta_time ! duration of main land time step (s)

! ---- indices of river tracers
integer  :: n_river_tracers
integer  :: i_river_ice, i_river_heat, i_river_DOC

! ---- diag field IDs --------------------------------------------------------
integer :: &
 ! COLUMN        VEGN        SNOW      GLAC/LAKE/SOIL  CANOPY-AIR  RIVER
  id_VWS,                                               id_VWSc,           &
  id_LWS,      id_LWSv,     id_LWSs,     id_LWSg,                          &
  id_FWS,      id_FWSv,     id_FWSs,     id_FWSg,                          &
  id_HS,       id_HSv,      id_HSs,      id_HSg,        id_HSc,            &
  id_precip,                                                               &
  id_hprec,                                                                &
  id_lprec,    id_lprecv,   id_lprecs,   id_lprecg,                        &
  id_hlprec,   id_hlprecv,  id_hlprecs,  id_hlprecg,                       &
  id_fprec,    id_fprecv,   id_fprecs,                                     &
  id_hfprec,   id_hfprecv,  id_hfprecs,                                    &
  id_evap,                                                                 &
  id_hevap,                                                                &
  id_levap,    id_levapv,   id_levaps,   id_levapg,                        &
  id_hlevap,   id_hlevapv,  id_hlevaps,  id_hlevapg,                       &
  id_fevap,    id_fevapv,   id_fevaps,   id_fevapg,                        &
  id_hfevap,   id_hfevapv,  id_hfevaps,  id_hfevapg,                       &
  id_runf,                                                                 &
  id_hrunf,                                                                &
  id_lrunf,                 id_lrunfs,   id_lrunfg,                        &
  id_hlrunf,                id_hlrunfs,  id_hlrunfg,                       &
  id_frunf,                 id_frunfs,                                     &
  id_hfrunf,                id_hfrunfs,                                    &
  id_melt,     id_meltv,    id_melts,    id_meltg,                         &
  id_fsw,      id_fswv,     id_fsws,     id_fswg,                          &
  id_flw,      id_flwv,     id_flws,     id_flwg,                          &
  id_sens,     id_sensv,    id_senss,    id_sensg,                         &
!
  id_e_res_1,  id_e_res_2,  id_cd_m,     id_cd_t,                          &
  id_cellarea, id_landfrac,                                                &
  id_geolon_t, id_geolat_t,                                                &
  id_frac,     id_area,     id_ntiles,                                     &
  id_z0m,      id_z0s,      id_con_g_h,                                    &
  id_transp,                id_wroff,    id_sroff,                         &
  id_htransp,  id_huptake,  id_hroff,    id_gsnow,    id_gequil,           &
  id_grnd_flux,                                                            &
  id_soil_water_supply,     id_levapg_max,                                 &
  id_water,    id_snow,                                                    &
  id_Trad,     id_Tca,      id_qca,      id_qco2_dvmr,                     &
  id_swdn_dir, id_swdn_dif, id_swup_dir, id_swup_dif, id_lwdn,             &
  id_fco2,                                                                 &
  id_vegn_cover,    id_cosz,                                               &
  id_albedo_dir,    id_albedo_dif,                                         &
  id_vegn_refl_dir, id_vegn_refl_dif, id_vegn_refl_lw,                     &
  id_vegn_tran_dir, id_vegn_tran_dif, id_vegn_tran_lw,                     &
  id_vegn_sctr_dir,                                                        &
  id_subs_refl_dir, id_subs_refl_dif, id_subs_emis, id_grnd_T,             &
  id_water_cons,    id_carbon_cons, id_DOCrunf
! diagnostic ids for canopy air tracers (moist mass ratio)
integer, allocatable :: id_cana_tr(:)
! diag IDs of CMOR variables
integer :: id_sftlf, id_sftgif
integer :: id_pcp, id_prra, id_prveg, id_tran, id_evspsblveg, id_evspsblsoi, id_nbp, &
           id_snw, id_snd, id_snc, id_lwsnl, id_snm, id_sweLut, id_cLand, &
           id_hflsLut, id_rlusLut, id_rsusLut, id_tslsiLut
integer :: id_cropFrac, id_cropFracC3, id_cropFracC4, id_pastureFrac, id_residualFrac, &
           id_grassFrac, id_grassFracC3, id_grassFracC4, &
           id_treeFrac, id_c3pftFrac, id_c4pftFrac, id_nwdFracLut, &
           id_fracLut_psl, id_fracLut_crp, id_fracLut_pst, id_fracLut_urb

! init_value is used to fill most of the allocated boundary condition arrays.
! It is supposed to be double-precision signaling NaN, to trigger a trap when
! the program is compiled with trapping non-initialized values.
! See http://ftp.uniovi.es/~antonio/uned/ieee754/IEEE-754references.html
! real, parameter :: init_value = Z'FFF0000000000001'
real, parameter :: init_value = 0.0

! ---- global clock IDs
integer :: landClock, landFastClock, landSlowClock

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),__FILE__,__LINE__)

contains


! ============================================================================
subroutine land_model_init &
     (cplr2land, land2cplr, time_init, time, dt_fast, dt_slow)
! initialize land model using grid description file as an input. This routine
! reads land grid boundaries and area of land from a grid description file

! NOTES: theoretically, the grid description file can specify any regular
! rectangular grid for land, not just lon/lat grid. Therefore the variables
! "xbl" and "ybl" in NetCDF grid spec file are not necessarily lon and lat
! boundaries of the grid.
!   However, at this time the module land_properties assumes that grid _is_
! lon/lat and therefore the entire module also have to assume that the land
! grid is lon/lat.
!   lon/lat grid is also assumed for the diagnostics, but this is probably not
! so critical.
  type(atmos_land_boundary_type), intent(inout) :: cplr2land ! boundary data
  type(land_data_type)          , intent(inout) :: land2cplr ! boundary data
  type(time_type), intent(in) :: time_init ! initial time of simulation (?)
  type(time_type), intent(in) :: time      ! current time
  type(time_type), intent(in) :: dt_fast   ! fast time step
  type(time_type), intent(in) :: dt_slow   ! slow time step

  ! ---- local vars ----------------------------------------------------------
  integer :: ncid, varid
  integer :: unit, ierr, io
  integer :: id_band, id_zfull ! IDs of land diagnostic axes
  integer :: id_ug !<Unstructured axis id.
  logical :: used                        ! return value of send_data diagnostics routine
  integer :: i,j,k,l
  type(land_tile_type), pointer :: tile
  type(land_tile_enum_type) :: ce
  integer :: ico2_atm ! index of CO2 tracer in the atmos, or NO_TRACER

  type(land_restart_type) :: restart
  character(*), parameter :: restart_file_name='INPUT/land.res.nc'
  logical :: restart_exists

  ! IDs of local clocks
  integer :: landInitClock

  module_is_initialized = .TRUE.

  ! [1] print out version number
  call log_version (version, module_name, &
  __FILE__)

  ! initialize land model clocks
  landClock      = mpp_clock_id('Land'               ,CLOCK_FLAG_DEFAULT,CLOCK_COMPONENT)
  landFastClock  = mpp_clock_id('Update-Land-Fast'   ,CLOCK_FLAG_DEFAULT,CLOCK_ROUTINE)
  landSlowClock  = mpp_clock_id('Update-Land-Slow'   ,CLOCK_FLAG_DEFAULT,CLOCK_ROUTINE)
  landInitClock  = mpp_clock_id('Land init'          ,CLOCK_FLAG_DEFAULT,CLOCK_ROUTINE)

  call mpp_clock_begin(landInitClock)

  ! [2] read land model namelist
#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=land_model_nml, iostat=io)
     ierr = check_nml_error(io, 'land_model_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=land_model_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'land_model_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=land_model_nml)
     call close_file (unit)
  endif

  ! initialize astronomy, in case it is not initialized, e.g. when using atmos_null
  call astronomy_init()

  ! initialize land state data, including grid geometry and processor decomposition
  call land_data_init(layout, io_layout, time, dt_fast, dt_slow, mask_table,npes_io_group)

  ! initialize land debug output
  call land_debug_init()

  ! initialize tile-specific diagnostics internals
  call tile_diag_init()

  ! initialize some tracer indices
  call land_tracers_init()
  call river_tracers_init()

  ! read sub-model namelists: then need to be read before component initialization
  ! because they can affect the way cover and tiling is initialized on cold start.
  ! Also, some of them register diagnostic sub-sampling selectors, so they better
  ! be after land_tile_diag_init
  call read_land_io_namelist()
  call read_soil_namelist()
  call read_hlsp_namelist() ! Must be called after read_soil_namelist
  call read_vegn_namelist()
  call read_soil_carbon_namelist()
  call read_lake_namelist()
  call read_glac_namelist()
  call read_snow_namelist()
  call read_cana_namelist()

  delta_time  = time_type_to_real(lnd%dt_fast) ! store in a module variable for convenience
  call init_tile_map()

  ! initialize subgrid tile distribution
  call open_land_restart(restart,restart_file_name,restart_exists)
  if(restart_exists) then
     ! read map of tiles -- retrieve information from
     call land_cover_warm_start(restart)
     ! initialize land model data
     if (field_exists(restart, 'lwup'   )) call get_tile_data(restart,'lwup',   land_lwup_ptr)
     if (field_exists(restart, 'e_res_1')) call get_tile_data(restart,'e_res_1',land_e_res_1_ptr)
     if (field_exists(restart, 'e_res_2')) call get_tile_data(restart,'e_res_2',land_e_res_2_ptr)
  else
     ! initialize map of tiles -- construct it by combining tiles
     ! from component models
     call error_mesg('land_model_init', 'cold-starting land cover map', NOTE)
     if (predefined_tiles) then
        call land_cover_cold_start_predefined()
     else
        call land_cover_cold_start()
     endif
  endif
  call free_land_restart(restart)

  ! initialize land model diagnostics -- must be before *_data_init so that
  ! *_data_init can write static fields if necessary
  call land_diag_init( lnd%coord_glonb, lnd%coord_glatb, lnd%coord_glon, lnd%coord_glat, &
      time, lnd%domain, id_band, id_ug )

  ! set the land diagnostic axes ids for the flux exchange
  land2cplr%axes = (/id_ug/)
  ! send some static diagnostic fields to output
  if ( id_cellarea > 0 ) used = send_data ( id_cellarea, lnd%cellarea,     lnd%time )
  if ( id_landfrac > 0 ) used = send_data ( id_landfrac, lnd%landfrac,     lnd%time )
  if ( id_geolon_t > 0 ) used = send_data ( id_geolon_t, lnd%lon*180.0/PI, lnd%time )
  if ( id_geolat_t > 0 ) used = send_data ( id_geolat_t, lnd%lat*180.0/PI, lnd%time )

  ! CMOR variables
  if ( id_sftgif > 0 ) call send_cellfrac_data(id_sftgif,is_glacier)
  if ( id_sftlf > 0 )  used = send_data(id_sftlf,lnd%landfrac*100, lnd%time)

  ! initialize individual sub-models
  if (predefined_tiles) then
     call hlsp_init_predefined(id_ug) ! Must be called before soil_init
  else
     call hlsp_init(id_ug) ! Must be called before soil_init
  endif
  call soil_init(predefined_tiles,id_ug,id_band,id_zfull)

  call hlsp_hydro_init(id_ug,id_zfull) ! Must be called after soil_init
  call vegn_init(id_ug,id_band)
  if (predefined_tiles) then
     call lake_init_predefined(id_ug)
  else
     call lake_init(id_ug)
  endif
  call glac_init(id_ug)
  call snow_init()
  call cana_init()
  call topo_rough_init(lnd%time, lnd_sg%lonb, lnd_sg%latb, lnd_sg%domain, lnd%domain, id_ug)
  call river_init( lnd_sg%lon, lnd_sg%lat, &
                   lnd%time, lnd%dt_fast, lnd_sg%domain,     &
                   lnd%domain, lnd_sg%landfrac,              &
                   discharge_tol, clw, csw )

  ! initialize river tracer indices
  n_river_tracers = num_river_tracers()
  i_river_ice  = river_tracer_index('ice')
  i_river_heat = river_tracer_index('het')
  i_river_DOC  = river_tracer_index('doc')
  if (i_river_ice  == NO_TRACER) call error_mesg ('land_model_init', &
                                      'required river tracer for ice not found', FATAL)
  if (i_river_heat == NO_TRACER) call error_mesg ('land_model_init', &
                                      'required river tracer for heat not found', FATAL)

  call land_transitions_init(id_ug, id_cellarea)

  call hlsp_config_check () ! Needs to be done after land_transitions_init and vegn_init

  call land_tracer_driver_init(id_ug)

  ! initialize boundary data
  ! allocate storage for the boundary data
  call realloc_land2cplr ( land2cplr )
  call realloc_cplr2land ( cplr2land )

  ! [8.2] set the land mask to FALSE everywhere -- update_land_bc_fast
  ! will set it to true where necessary
  land2cplr%mask = .FALSE.
  land2cplr%tile_size = 0.0
  ! [8.3] get the current state of the land boundary for the coupler
  ce = first_elmt(land_tile_map, ls=lnd%ls )
  do while(loop_over_tiles(ce,tile, l,k))
     call set_current_point(l,k)
     call update_land_bc_fast (tile, l,k, land2cplr, is_init=.true.)
  enddo

  ! [8.4] update topographic roughness scaling
  call update_land_bc_slow( land2cplr )

  ! mask error checking
  do l=lnd%ls,lnd%le
     if(lnd%landfrac(l)>0.neqv.ANY(land2cplr%mask(l,:))) then
        call error_mesg('land_model_init','land masks from grid spec and from land restart do not match',FATAL)
     endif
  enddo

  ! [9] check the properties of co2 exchange with the atmosphere and set appropriate
  ! flags. Since co2 tracer is always present in the land tracer table, we use
  ! the presence of the atmos co2 tracers as an indication that the exchange of
  ! CO2 with the atmos is set up.
  ico2_atm = get_tracer_index(MODEL_ATMOS, 'co2')
  if (canopy_air_mass_for_tracers==0.and.ico2_atm==NO_TRACER) then
     call error_mesg('land_model_init', &
          'canopy_air_mass_for_tracers is set to zero, and CO2 exchange with the atmosphere is not set up: '// &
          'canopy air CO2 concentration will not be updated',NOTE)
     update_cana_co2 = .FALSE.
  else
     update_cana_co2 = .TRUE.
  end if

  call mpp_clock_end(landInitClock)


end subroutine land_model_init


! ============================================================================
subroutine land_model_end (cplr2land, land2cplr)
  type(atmos_land_boundary_type), intent(inout) :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr

  module_is_initialized = .FALSE.

  call error_mesg('land_model_end','writing NetCDF restart',NOTE)
  call land_model_restart()

  ! we still want to call the *_end procedures for component models, even
  ! if the number of tiles in this domain is zero, in case they are doing
  ! something else besides saving the restart, of if they want to save
  ! restart anyway
  call land_tracer_driver_end()
  call land_transitions_end()
  call glac_end ()
  call lake_end ()
  call soil_end ()
  call hlsp_end ()
  call snow_end ()
  call vegn_end ()
  call cana_end ()
  call topo_rough_end()
  call river_end()

  ! deallocate storage allocated for tracers
  deallocate(id_cana_tr)

  call dealloc_land2cplr(land2cplr, dealloc_discharges=.TRUE.)
  call dealloc_cplr2land(cplr2land)

  call land_tracers_end()

  call tile_diag_end()

  ! deallocate tiles
  call free_tile_map()
  call land_data_end()

  ! finish up the land debugging diagnostics
  call land_debug_end()

end subroutine land_model_end


! ============================================================================
subroutine land_model_restart(timestamp)
  character(*), intent(in), optional :: timestamp ! timestamp to add to the file name

  ! ---- local vars
  character(267) :: filename
  type(land_restart_type) :: restart ! restart file i/o object
  integer :: tile_dim_length ! length of tile dimension in output files
                             ! global max of number of tiles per gridcell
  integer :: l
  character(256) :: timestamp_

  ! [1] count all land tiles and determine the length of tile dimension
  ! sufficient for the current domain
  tile_dim_length = 0
  do l = lnd%ls, lnd%le
     tile_dim_length = max(tile_dim_length,nitems(land_tile_map(l)))
  enddo

  ! [2] calculate the tile dimension length by taking the max across all domains
  call mpp_max(tile_dim_length)
  if (tile_dim_length==0) then
     call error_mesg('land_model_restart',&
       'No land points exist in global domain (either entire map, or cubic sphere face), therefore no land restarts will be saved',&
       WARNING)
     return
  endif

  ! [3] create tile output file
  timestamp_=''
  if (present(timestamp)) then
     if(trim(timestamp)/='') timestamp_=trim(timestamp)//'.'
  endif
  ! Note that filename is updated for tile & rank numbers during file creation
  filename = trim(timestamp_)//'land.res.nc'
  call init_land_restart(restart, filename, land_tile_exists, tile_dim_length)

  ! [4] write data fields
  ! write fractions and tile tags
  call add_tile_data(restart,'frac',land_frac_ptr,'fractional area of tile')
  call add_int_tile_data(restart,'glac',glac_tag_ptr,'tag of glacier tiles')
  call add_int_tile_data(restart,'lake',lake_tag_ptr,'tag of lake tiles')
  call add_int_tile_data(restart,'soil',soil_tag_ptr,'tag of soil tiles')
  call add_int_tile_data(restart,'vegn',vegn_tag_ptr,'tag of vegetation tiles')
  ! write the upward long-wave flux
  call add_tile_data(restart,'lwup',land_lwup_ptr,'upward long-wave flux')
  ! write energy residuals
  call add_tile_data(restart,'e_res_1',land_e_res_1_ptr,&
       'energy residual in canopy air energy balance equation', 'W/m2')
  call add_tile_data(restart,'e_res_2',land_e_res_2_ptr,&
       'energy residual in canopy energy balance equation', 'W/m2')

  ! [5] close file
  call save_land_restart(restart)
  call free_land_restart(restart)

  ! [6] save component models' restarts
  call save_land_transitions_restart(timestamp_)
  call save_glac_restart(tile_dim_length,timestamp_)
  call save_lake_restart(tile_dim_length,timestamp_)
  call save_soil_restart(tile_dim_length,timestamp_)
  call save_hlsp_restart(tile_dim_length,timestamp_)
  call save_snow_restart(tile_dim_length,timestamp_)
  call save_vegn_restart(tile_dim_length,timestamp_)
  call save_cana_restart(tile_dim_length,timestamp_)
  call save_river_restart(timestamp_)

end subroutine land_model_restart

! ============================================================================
subroutine land_cover_cold_start()
  real, dimension(:,:), pointer :: &
       glac, soil, lake, vegn ! arrays of fractions for respective sub-models
  integer, pointer, dimension(:,:) :: soiltags ! array of soil type tags
  integer, pointer, dimension(:,:) :: hlsp_pos ! hillslope position index
  integer, pointer, dimension(:,:) :: hlsp_par ! hillslope parent index
  real   , pointer, dimension(:,:) :: rbuffer ! real buffer for remap
  logical, dimension(lnd%le-lnd%ls+1) :: &
       land_mask, valid_data, invalid_data
  integer :: iwatch,jwatch,kwatch,face,lwatch
  integer :: l,p,lll
  integer :: tile_root_pe,tile_npes
  integer :: ps,pe ! boundaries of PE list for remapping
  type(horiz_remap_type) :: map

  ! calculate the global land mask
  land_mask = lnd%area > 0

  ! get the global maps of fractional covers for each of the sub-models
  glac=>glac_cover_cold_start(land_mask,lnd_sg%lonb,lnd_sg%latb)
  lake=>lake_cover_cold_start(land_mask,lnd_sg%lonb,lnd_sg%latb,lnd_sg%domain)
  soil=>soil_cover_cold_start(land_mask,lnd_sg%lonb,lnd_sg%latb)
  vegn=>vegn_cover_cold_start(land_mask,lnd_sg%lonb,lnd_sg%latb)

  ! Because of hillslope model, soil tiles may not be returned in order of soil type.
  allocate(soiltags(size(soil,1), size(soil,2)))
  call retrieve_soil_tags(soiltags)
  ! Tiles will be constructed with hillslope data.
  allocate(hlsp_pos(size(soil,1), size(soil,2)))
  allocate(hlsp_par(size(soil,1), size(soil,2)))
  call retrieve_hlsp_indices(hlsp_pos, hlsp_par)

  ! remove any input lake fraction in coastal cells
  where (lnd%landfrac.lt. 1.-gfrac_tol) lake(:,1) = 0.
  ! NOTE that the lake area in the coastal cells can be set to non-zero
  ! again by the "ground fraction reconciliation code" below. Strictly
  ! speaking the above line of code should be replaced with the section
  ! commented out with "!-zero" below, but we preserve the old way to avoid
  ! backward incompatibility with older runs. This needs updating in the
  ! future when the decision about what to do with lakes in coastal cells is
  ! made.

  ! reconcile ground fractions with the land mask within compute domain
  valid_data = land_mask.and.(sum(glac,2)+sum(lake,2)+sum(soil,2)>0)
  invalid_data = land_mask.and..not.valid_data

  call get_watch_point(iwatch,jwatch,kwatch,face,lwatch)
  if (face==lnd%face.and.(lnd%ls<=lwatch.and.lwatch<=lnd%le) ) then
     write(*,*)'###### land_cover_cold_start: input data #####'
     write(*,'(99(a,i4.2,x))')'iwatch=',iwatch,'jwatch=',jwatch,'face=',face
     write(*,'(99(a,g23.16,x))')'lon=',lnd%lon(lwatch)*180/PI,'lat=',lnd%lat(lwatch)*180/PI
     ! calculate local compute domain indices; we assume glac,lake,soil,vegn all
     ! have the same lbounds
     l = lwatch-lnd%ls+lbound(glac,1)
     __DEBUG1__(lnd%ls)
     write(*,'(a,99(a,i4.2,x))')'local indices:','l=',l
     __DEBUG3__(lnd%landfrac(lwatch),land_mask(l),valid_data(l))
     __DEBUG1__(glac(l,:))
     __DEBUG1__(lake(l,:))
     __DEBUG1__(soil(l,:))
     __DEBUG1__(vegn(l,:))
  endif

  if (trim(nearest_point_search)=='global') then
     ps=0 ; pe=size(lnd%pelist)-1
  else if (trim(nearest_point_search)=='face') then
     call mpp_get_UG_domain_tile_pe_inf(lnd%domain, tile_root_pe, npes=tile_npes)
     ps = -1
     do p = 0, size(lnd%pelist(:))-1
        if(lnd%pelist(p) == tile_root_pe) then
           ps = p
           exit
        endif
     enddo
     if( ps == -1) call error_mesg('land_cover_cold_start',&
           'tile_root_pe is not in the lnd%pelist', FATAL)
     pe = ps + tile_npes-1
  else
     call error_mesg('land_cover_cold_start',&
          'option nearest_point_search="'//trim(nearest_point_search)//&
          '" is illegal, use "global" or "face"',&
          FATAL)
  endif
  call horiz_remap_new(invalid_data,valid_data,lnd%lon,lnd%lat,lnd%domain,&
          lnd%pelist(ps:pe),map)
  if (print_remapping) call horiz_remap_print(map,'land cover remap:')
  call horiz_remap(map,lnd%domain,glac)
  call horiz_remap(map,lnd%domain,lake)
  call horiz_remap(map,lnd%domain,soil)
  allocate(rbuffer(size(soil,1), size(soil,2)))
  ! ZMS This is awkward: perhaps horiz_remap should handle integers?
  rbuffer(:,:) = real(soiltags(:,:))
  call horiz_remap(map,lnd%domain,rbuffer)
  soiltags(:,:) = nint(rbuffer(:,:))
  rbuffer(:,:) = real(hlsp_pos(:,:))
  call horiz_remap(map,lnd%domain,rbuffer)
  hlsp_pos(:,:) = nint(rbuffer(:,:))
  rbuffer(:,:) = real(hlsp_par(:,:))
  call horiz_remap(map,lnd%domain,rbuffer)
  hlsp_par(:,:) = nint(rbuffer(:,:))
  call horiz_remap_del(map)

!-zero  ! remove any input lake fraction in coastal cells
!-zero  do j = lnd_sg%js,lnd_sg%je
!-zero  do i = lnd_sg%is,lnd_sg%ie
!-zero     call set_current_point(i,j,1)
!-zero     if (lnd_sg%landfrac(i,j) < 1-gfrac_tol) then
!-zero        lake(i,j,:) = 0.0
!-zero        if(is_watch_point())then
!-zero           write(*,*)'###### land_cover_cold_start: lake fraction is set to zero #####'
!-zero        endif
!-zero     endif
!-zero  enddo
!-zero  enddo

  ! reconcile vegetation fractions with the land mask within compute domain
  valid_data = sum(vegn,2) > 0
  invalid_data = .FALSE.
  do l = 1,size(land_mask(:))
     if(.not.land_mask(l)) cycle ! skip ocean points
     if(valid_data(l)) cycle ! do not need to do anything with valid points
     if(sum(glac(l,:))+sum(lake(l,:))>=1) &
          cycle                ! skip points fully covered by glaciers or lakes
     invalid_data(l)=.TRUE.
  enddo
  call horiz_remap_new(invalid_data,valid_data,lnd%lon,lnd%lat,lnd%domain,&
       lnd%pelist(ps:pe),map)
  if (print_remapping) call horiz_remap_print(map,'vegetation cover remap:')
  call horiz_remap(map,lnd%domain,vegn)
  call horiz_remap_del(map)

  ! create tiles
  do l = 1,size(land_mask(:))
     lll = l+lnd%ls-1
     if(.not.land_mask(l)) cycle ! skip ocean points
     call set_current_point(lll,1)
     call land_cover_cold_start_0d &
          (land_tile_map(lll),glac(l,:),lake(l,:),soil(l,:),soiltags(l,:),&
               hlsp_pos(l,:), hlsp_par(l,:), vegn(l,:))
     if(nitems(land_tile_map(lll))==0) then
        call error_mesg('land_cover_cold_start',&
             'No tiles were created for a valid land point at i='&
             //trim(string(lnd%i_index(lll)))//' j='//trim(string(lnd%j_index(lll)))//' face='//trim(string(lnd%face)), FATAL)
     endif
  enddo

  deallocate(glac,lake,soil,soiltags,hlsp_pos,hlsp_par,vegn,rbuffer)
end subroutine land_cover_cold_start

! ============================================================================
subroutine land_cover_cold_start_predefined()
  integer :: l,h5id

  ! Open access to model input database
  call open_database_predefined_tiles(h5id)

  do l = lnd%ls, lnd%le
    if(.not.lnd%area(l)>0) cycle ! skip ocean points
    call set_current_point(l,1)
    call land_cover_cold_start_0d_predefined_tiles(land_tile_map(l),lnd,l,h5id)
  enddo

  ! Close access to model input database
  call close_database_predefined_tiles(h5id)
end subroutine land_cover_cold_start_predefined

! ============================================================================
subroutine land_cover_cold_start_0d (set,glac0,lake0,soil0,soiltags0,&
                                     hlsp_pos0,hlsp_par0,vegn0)
  type(land_tile_list_type), intent(inout) :: set
  real, dimension(:)       , intent(in) :: &
       glac0,lake0,soil0,vegn0 ! fractions of area
  integer, dimension(:)    , intent(in) :: &
       soiltags0, hlsp_pos0, hlsp_par0 ! soil and hillslope tags

  ! ---- local vars
  real :: glac(size(glac0(:))), lake(size(lake0(:))), &
          soil(size(soil0(:))), vegn(size(vegn0(:)))
  type(land_tile_type), pointer :: tile
  integer :: i,j,k
  real :: factor ! normalizing factor for the tile areas
  real :: frac
  type(land_tile_enum_type) :: first_non_vegn ! position of first non-vegetated tile in the list
  type(land_tile_enum_type) :: ce
  real :: athresh = 1.e-10 ! area threshold allowable for deviation from 1
  real :: area ! area sum for gridcell

  glac = glac0; lake = lake0; soil = soil0; vegn = vegn0
  if (sum(glac)>1) &
       glac=glac/sum(glac)
  if (sum(lake)+sum(glac)>1)&
       lake = lake*(1-sum(glac))/sum(lake)
  if (sum(lake)<min_sum_lake_frac) lake=0
  if (sum(soil)+sum(glac)+sum(lake)>1)&
       soil = soil*(1-sum(lake)-sum(glac))/sum(soil)
  ! make sure that the sum of the fractions of the soil, lake, and glaciers are
  ! either one or zero
  factor = sum(soil)+sum(glac)+sum(lake)
  if(factor>0)then
     glac = glac/factor
     lake = lake/factor
     soil = soil/factor
  endif

  ! remove soil/glac/lake fractions that are too small
  if (min_frac>0) then
     where (glac<min_frac) glac = 0
     where (lake<min_frac) lake = 0
     where (soil<min_frac) soil = 0
     ! do the renormalization again
     factor = sum(soil)+sum(glac)+sum(lake)
     if(factor>0)then
	glac = glac/factor
	lake = lake/factor
	soil = soil/factor
     endif
  endif

  if(is_watch_point()) then
     write(*,*)'#### land_cover_cold_start_0d input data ####'
     __DEBUG1__(glac0)
     __DEBUG1__(lake0)
     __DEBUG1__(soil0)
     __DEBUG1__(vegn0)
     __DEBUG1__(factor)
     write(*,*)'#### land_cover_cold_start_0d renormalized fractions ####'
     __DEBUG1__(glac)
     __DEBUG1__(lake)
     __DEBUG1__(soil)
     __DEBUG1__(vegn)
  endif

  do i = 1,size(glac)
     if (glac(i)>0) then
        tile => new_land_tile(frac=glac(i),glac=i)
        call insert(tile,set)
        if(is_watch_point()) then
           write(*,*)'created glac tile: frac=',glac(i),' tag=',i
        endif
     endif
  enddo
  do i = 1,size(lake)
     if (lake(i)>0) then
        tile => new_land_tile(frac=lake(i),lake=i)
        call insert(tile,set)
        if(is_watch_point()) then
           write(*,*)'created lake tile: frac=',lake(i),' tag=',i
        endif
     endif
  enddo

  factor = sum(soil)*sum(vegn)
  if (factor/=0) factor = 1/factor
  factor = factor*(1-sum(glac)-sum(lake))
  ! vegetation tiles, if any, are inserted in front of non-vegetated tiles;
  ! this really does not matter except for the static vegetation override
  ! case with the data saved by lm3v -- there the vegetation tiles are
  ! in front, so it works more consistently where lad2 has more than
  ! one tile (e.g. glac/soil or lake/soil), if lad2 vegetation tiles are
  ! also in front of the list.
  first_non_vegn=first_elmt(set)
  do i = 1,size(soil)
  do j = 1,size(vegn)
     frac = soil(i)*vegn(j)*factor
     if(frac>0) then
        tile  => new_land_tile(frac=frac,soil=soiltags0(i),vegn=j,&
                               htag_j=hlsp_pos0(i),htag_k=hlsp_par0(i))
        call insert(tile,first_non_vegn)
        if(is_watch_point()) then
           write(*,*)'created soil tile: frac=', frac, ' soil tag=',soiltags0(i), ' veg tag=',j
        endif
     endif
  enddo
  enddo

  ! Check that fractions are set correctly.
  ce = first_elmt(set) ; area = 0.
  do while (loop_over_tiles(ce, tile))
     area = area + tile%frac
  end do

  if (abs(area - 1.) > athresh) then
     call error_mesg(module_name, 'Area fractions do not add up to 1 for tile list in land_cover_cold_start_0d!', &
                     FATAL)
  end if

end subroutine land_cover_cold_start_0d

! ============================================================================
subroutine land_cover_warm_start(restart)
  type(land_restart_type), intent(in) :: restart
  if (new_land_io) then
     call land_cover_warm_start_new(restart)
  else
     call land_cover_warm_start_orig(restart)
  endif
end subroutine land_cover_warm_start

! ============================================================================
! reads the land restart file and restores the tiling structure from this file
subroutine land_cover_warm_start_new (restart)
  type(land_restart_type), intent(in) :: restart

  ! ---- local vars
  integer, allocatable :: glac(:), lake(:), soil(:), vegn(:) ! tile tags
  real,    allocatable :: frac(:) ! fraction of land covered by tile
  integer :: ntiles    ! total number of land tiles in the input file
  integer :: k,it,l,g, npts
  type(land_tile_type), pointer :: tile

  ntiles = size(restart%tidx)
  allocate(glac(ntiles), lake(ntiles), soil(ntiles), vegn(ntiles), frac(ntiles))

  call fms_io_unstructured_read(restart%basename, "frac", frac, lnd%domain, timelevel=1)
  call fms_io_unstructured_read(restart%basename, "glac", glac, lnd%domain, timelevel=1)
  call fms_io_unstructured_read(restart%basename, "lake", lake, lnd%domain, timelevel=1)
  call fms_io_unstructured_read(restart%basename, "soil", soil, lnd%domain, timelevel=1)
  call fms_io_unstructured_read(restart%basename, "vegn", vegn, lnd%domain, timelevel=1)

  npts = lnd%nlon*lnd%nlat
  ! create tiles
  do it = 1,ntiles
     k = restart%tidx(it)
     if (k<0) cycle ! skip negative indices
     g = modulo(k,npts)+1
     if (g<lnd%gs.or.g>lnd%ge) cycle ! skip points outside of domain
     l = lnd%l_index(g)
     ! the size of the tile set at the point (i,j) must be equal to k
     tile=>new_land_tile(frac=frac(it),&
              glac=glac(it),lake=lake(it),soil=soil(it),vegn=vegn(it))
     call insert(tile,land_tile_map(l))
  enddo
  deallocate(glac, lake, soil, vegn, frac)
end subroutine land_cover_warm_start_new


! ============================================================================
! reads the land restart file and restores the tiling structure from this file
subroutine land_cover_warm_start_orig (restart)
  type(land_restart_type), intent(in) :: restart

  ! ---- local vars
  integer, allocatable :: idx(:) ! compressed tile index
  integer, allocatable :: glac(:), lake(:), soil(:), snow(:), cana(:), vegn(:) ! tile tags
  real,    allocatable :: frac(:) ! fraction of land covered by tile
  integer :: ncid ! unit number of the input file
  integer :: ntiles    ! total number of land tiles in the input file
  integer :: bufsize   ! size of the input buffer
  integer :: dimids(1) ! id of tile dimension
  character(NF_MAX_NAME) :: tile_dim_name ! name of the tile dimension and respective variable
  integer :: k,it,npts,g,l
  type(land_tile_type), pointer :: tile;
  integer :: start, count ! slab for reading
  ! netcdf variable IDs
  integer :: id_idx, id_frac, id_glac, id_lake, id_soil, id_vegn

  __NF_ASRT__(nf_open(restart%filename,NF_NOWRITE,ncid))
  ! allocate the input data
  __NF_ASRT__(nfu_inq_var(ncid,'frac',id=id_frac,varsize=ntiles,dimids=dimids))
   ! allocate input buffers for compression index and the variable
  bufsize=min(input_buf_size,ntiles)
  allocate(idx (bufsize), glac(bufsize), lake(bufsize), soil(bufsize), &
           snow(bufsize), cana(bufsize), vegn(bufsize), frac(bufsize)  )
  ! get the name of the fist (and only) dimension of the variable 'frac' -- this
  ! is supposed to be the compressed dimension, and associated variable will
  ! hold the compressed indices
  __NF_ASRT__(nfu_inq_dim(ncid,dimids(1),name=tile_dim_name))
  __NF_ASRT__(nfu_inq_var(ncid,tile_dim_name,id=id_idx))
  ! get the IDs of the variables to read
  __NF_ASRT__(nfu_inq_var(ncid,'glac',id=id_glac))
  __NF_ASRT__(nfu_inq_var(ncid,'lake',id=id_lake))
  __NF_ASRT__(nfu_inq_var(ncid,'soil',id=id_soil))
  __NF_ASRT__(nfu_inq_var(ncid,'vegn',id=id_vegn))

  npts = lnd%nlon*lnd%nlat
  do start = 1,ntiles,bufsize
    count = min(bufsize,ntiles-start+1)
    ! read the compressed tile indices
    __NF_ASRT__(nf_get_vara_int(ncid,id_idx,(/start/),(/count/),idx))
    ! read input data -- fractions and tags
    __NF_ASRT__(nf_get_vara_double(ncid,id_frac,(/start,1/),(/count,1/),frac))
    __NF_ASRT__(nf_get_vara_int(ncid,id_glac,(/start,1/),(/count,1/),glac))
    __NF_ASRT__(nf_get_vara_int(ncid,id_lake,(/start,1/),(/count,1/),lake))
    __NF_ASRT__(nf_get_vara_int(ncid,id_soil,(/start,1/),(/count,1/),soil))
    __NF_ASRT__(nf_get_vara_int(ncid,id_vegn,(/start,1/),(/count,1/),vegn))
    ! create tiles
    do it = 1,count
       k = idx(it)
       if (k<0) cycle ! skip negative indices
       g = modulo(k,npts)+1
       if (g<lnd%gs.or.g>lnd%ge) cycle ! skip points outside of domain
       l = lnd%l_index(g)
       ! the size of the tile set at the point (i,j) must be equal to k
       tile=>new_land_tile(frac=frac(it),&
                glac=glac(it),lake=lake(it),soil=soil(it),vegn=vegn(it))
       call insert(tile,land_tile_map(l))
    enddo
  enddo
  __NF_ASRT__(nf_close(ncid))
  deallocate(idx, glac, lake, soil, snow, cana, vegn, frac)
end subroutine


! ============================================================================
subroutine update_land_model_fast ( cplr2land, land2cplr )
  type(atmos_land_boundary_type), intent(in)    :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr

  ! ---- local vars
  real :: &
     ISa_dn_dir(NBANDS), & ! downward direct sw radiation at the top of the canopy
     ISa_dn_dif(NBANDS), & ! downward diffuse sw radiation at the top of the canopy
     cana_q

  ! variables for stock calculations
  real :: &
     cana_VMASS, cana_HEAT,             &
     vegn_LMASS, vegn_FMASS, vegn_HEAT, &
     snow_LMASS, snow_FMASS, snow_HEAT, &
     subs_LMASS, subs_FMASS, subs_HEAT, &
     glac_LMASS, glac_FMASS, glac_HEAT, &
     lake_LMASS, lake_FMASS, lake_HEAT, &
     soil_LMASS, soil_FMASS, soil_HEAT

  real, dimension(lnd_sg%is:lnd_sg%ie,lnd_sg%js:lnd_sg%je) :: runoff_sg
  real, dimension(lnd_sg%is:lnd_sg%ie,lnd_sg%js:lnd_sg%je,n_river_tracers) :: runoff_c_sg

  real, dimension(lnd%ls:lnd%le) :: &
       runoff            ! total (liquid+snow) runoff accumulated over tiles in cell
  real, dimension(lnd%ls:lnd%le,n_river_tracers) :: &
       runoff_c          ! runoff of tracers accumulated over tiles in cell (including ice and heat)
  logical :: used          ! return value of send_data diagnostics routine
  real, allocatable :: runoff_1d(:),runoff_snow_1d(:),runoff_heat_1d(:)
  integer :: i,j,k,l     ! lon, lat, and tile indices
  integer :: i_species ! river tracer iterator
  integer :: i1        ! index used to iterate over grid cells efficiently
  integer :: is,ie,js,je ! horizontal bounds of the override buffer
  type(land_tile_enum_type) :: ce ! tile enumerator
  type(land_tile_type), pointer :: tile ! pointer to current tile

  ! variables for data override
  real, allocatable :: phot_co2_data(:)    ! buffer for data
  logical           :: phot_co2_overridden ! flag indicating successful override
  integer           :: iwatch,jwatch,kwatch,face


  ! start clocks
  call mpp_clock_begin(landClock)
  call mpp_clock_begin(landFastClock)

  ! to avoid output of static vegetation after the transitions worked and
  ! changed the tiling structure, static vegetation output is done here.
  call write_static_vegn()

  ! override data at the beginning of the time step
  is= lnd_sg%is ; ie = lnd_sg%ie
  js= lnd_sg%js ; je = lnd_sg%je
  allocate(phot_co2_data(lnd%ls:lnd%le))
  phot_co2_data = 0
  call data_override_ug('LND','phot_co2',phot_co2_data,lnd%time, &
       override=phot_co2_overridden)

  ! clear the runoff values, for accumulation over the tiles
  runoff = 0 ; runoff_c = 0

  ! Calculate groundwater and associated heat fluxes between tiles within each gridcell.
  call hlsp_hydrology_1(n_c_types)

  ! main tile loop
!$OMP parallel do default(none) shared(lnd,land_tile_map,cplr2land,land2cplr,phot_co2_overridden, &
!$OMP                                  phot_co2_data,runoff,runoff_c,id_area,id_z0m,id_z0s,       &
!$OMP                                  id_Trad,id_Tca,id_qca,isphum,id_cd_m,id_cd_t) &
!$OMP                                  private(i1,i,j,k,ce,tile,ISa_dn_dir,ISa_dn_dif)
  do l = lnd%ls, lnd%le
     i = lnd%i_index(l)
     j = lnd%j_index(l)
!     __DEBUG4__(is,js,i-is+lnd_sg%is,j-js+lnd_sg%js)
     ce = first_elmt(land_tile_map(l))
     do while (loop_over_tiles(ce,tile,k=k))
        ! set this point coordinates as current for debug output
        call set_current_point(i,j,k,l)

        ISa_dn_dir(BAND_VIS) = cplr2land%sw_flux_down_vis_dir(l,k)
        ISa_dn_dir(BAND_NIR) = cplr2land%sw_flux_down_total_dir(l,k)&
                              -cplr2land%sw_flux_down_vis_dir(l,k)
        ISa_dn_dif(BAND_VIS) = cplr2land%sw_flux_down_vis_dif(l,k)
        ISa_dn_dif(BAND_NIR) = cplr2land%sw_flux_down_total_dif(l,k)&
                              -cplr2land%sw_flux_down_vis_dif(l,k)

        call update_land_model_fast_0d(tile, l, k, land2cplr, &
           cplr2land%lprec(l,k),  cplr2land%fprec(l,k), cplr2land%tprec(l,k), &
           cplr2land%t_flux(l,k), cplr2land%dhdt(l,k), &
           cplr2land%tr_flux(l,k,:), cplr2land%dfdtr(l,k,:), &
           ISa_dn_dir, ISa_dn_dif, cplr2land%lwdn_flux(l,k), &
           cplr2land%ustar(l,k), cplr2land%p_surf(l,k), cplr2land%drag_q(l,k), &
           phot_co2_overridden, phot_co2_data(l),&
           runoff(l), runoff_c(l,:) &
        )
        ! some of the diagnostic variables are sent from here, purely for coding
        ! convenience: the compute domain-level 2d and 3d vars are generally not
        ! available inside update_land_model_fast_0d, so the diagnostics for those
        ! was left here.
        call send_tile_data(id_area, tile%frac*lnd%area(l),        tile%diag)
        call send_tile_data(id_z0m,  land2cplr%rough_mom(l,k),     tile%diag)
        call send_tile_data(id_z0s,  land2cplr%rough_heat(l,k),    tile%diag)
        call send_tile_data(id_Trad, land2cplr%t_surf(l,k),        tile%diag)
        call send_tile_data(id_Tca,  land2cplr%t_ca(l,k),          tile%diag)
        call send_tile_data(id_qca,  land2cplr%tr(l,k,isphum),     tile%diag)
        call send_tile_data(id_cd_m, cplr2land%cd_m(l,k),          tile%diag)
        call send_tile_data(id_cd_t, cplr2land%cd_t(l,k),          tile%diag)
    enddo
  enddo

  !--- pass runoff from unstructured grid to structured grid.
  runoff_sg = 0 ; runoff_c_sg = 0
  call mpp_pass_UG_to_SG(lnd%domain, runoff,   runoff_sg  )
  call mpp_pass_UG_to_SG(lnd%domain, runoff_c, runoff_c_sg)

  call get_watch_point(iwatch,jwatch,kwatch,face)
  if (face==lnd%face.and.(lnd_sg%is<=iwatch.and.iwatch<=lnd_sg%ie).and.&
                         (lnd_sg%js<=jwatch.and.jwatch<=lnd_sg%je)) then
     __DEBUG1__(runoff_sg(iwatch,jwatch))
     __DEBUG1__(runoff_c_sg(iwatch,jwatch,:))
  endif

  !--- update river state
  call update_river(runoff_sg, runoff_c_sg, land2cplr)

  ce = first_elmt(land_tile_map, ls=lbound(cplr2land%t_flux,1) )
  do while(loop_over_tiles(ce,tile,l,k))
     cana_VMASS = 0. ;                   cana_HEAT = 0.
     vegn_LMASS = 0. ; vegn_FMASS = 0. ; vegn_HEAT = 0.
     snow_LMASS = 0. ; snow_FMASS = 0. ; snow_HEAT = 0.
     subs_LMASS = 0. ; subs_FMASS = 0. ; subs_HEAT = 0.
     glac_LMASS = 0. ; glac_FMASS = 0. ; glac_HEAT = 0.
     lake_LMASS = 0. ; lake_FMASS = 0. ; lake_HEAT = 0.
     soil_LMASS = 0. ; soil_FMASS = 0. ; soil_HEAT = 0.
     if (associated(tile%cana)) then
         call cana_state( tile%cana, cana_q=cana_q )
         cana_VMASS = canopy_air_mass*cana_q
         cana_HEAT  = cana_tile_heat(tile%cana)
       endif
     if (associated(tile%vegn)) then
         call vegn_tile_stock_pe(tile%vegn, vegn_LMASS, vegn_FMASS)
         vegn_HEAT = vegn_tile_heat(tile%vegn)
       endif
     if(associated(tile%snow)) then
         call snow_tile_stock_pe(tile%snow, snow_LMASS, snow_FMASS)
         snow_HEAT = snow_tile_heat(tile%snow)
     endif
     if (associated(tile%glac)) then
        call glac_tile_stock_pe(tile%glac, subs_LMASS, subs_FMASS)
        subs_HEAT  = glac_tile_heat(tile%glac)
        glac_LMASS = subs_LMASS
        glac_FMASS = subs_FMASS
        glac_HEAT  = subs_HEAT
     else if (associated(tile%lake)) then
        call lake_tile_stock_pe(tile%lake, subs_LMASS, subs_FMASS)
        subs_HEAT  = lake_tile_heat(tile%lake)
        lake_LMASS = subs_LMASS
        lake_FMASS = subs_FMASS
        lake_HEAT  = subs_HEAT
     else if (associated(tile%soil)) then
        call soil_tile_stock_pe(tile%soil, subs_LMASS, subs_FMASS)
        subs_HEAT  = soil_tile_heat(tile%soil)
        soil_LMASS = subs_LMASS
        soil_FMASS = subs_FMASS
        soil_HEAT  = subs_HEAT
     endif

     call send_tile_data(id_VWS,  cana_VMASS, tile%diag)
     call send_tile_data(id_VWSc, cana_VMASS, tile%diag)
     call send_tile_data(id_LWS,  vegn_LMASS+snow_LMASS+subs_LMASS, tile%diag)
     call send_tile_data(id_LWSv, vegn_LMASS, tile%diag)
     call send_tile_data(id_LWSs, snow_LMASS, tile%diag)
     call send_tile_data(id_LWSg, subs_LMASS, tile%diag)
     call send_tile_data(id_FWS,  vegn_FMASS+snow_FMASS+subs_FMASS, tile%diag)
     call send_tile_data(id_FWSv, vegn_FMASS, tile%diag)
     call send_tile_data(id_FWSs, snow_FMASS, tile%diag)
     call send_tile_data(id_FWSg, subs_FMASS, tile%diag)
     call send_tile_data(id_HS,  vegn_HEAT+snow_HEAT+subs_HEAT+cana_HEAT, tile%diag)
     call send_tile_data(id_HSv, vegn_HEAT, tile%diag)
     call send_tile_data(id_HSs, snow_HEAT, tile%diag)
     call send_tile_data(id_HSg, subs_HEAT, tile%diag)
     call send_tile_data(id_HSc, cana_HEAT, tile%diag)
     call send_tile_data(id_water, subs_LMASS+subs_FMASS, tile%diag)
     call send_tile_data(id_snow,  snow_LMASS+snow_FMASS, tile%diag)

     ! CMOR variables
     call send_tile_data(id_snw, snow_FMASS, tile%diag)
     call send_tile_data(id_lwsnl, snow_LMASS, tile%diag)
     ! factor 1000.0 kg/m3 is the liquid water density; it converts mass of water into depth
     call send_tile_data(id_sweLut, max(snow_FMASS+snow_LMASS,0.0)/1000.0, tile%diag)
  enddo

  ! advance land model time
  lnd%time = lnd%time + lnd%dt_fast

  ! send the accumulated diagnostics to the output
  call dump_tile_diag_fields(land_tile_map,lnd%time)

  ! send CMOR cell fraction fields
  call send_cellfrac_data(id_cropFrac,     is_crop)
  call send_cellfrac_data(id_cropFracC3,   is_crop_C3)
  call send_cellfrac_data(id_cropFracC4,   is_crop_C4)
  call send_cellfrac_data(id_pastureFrac,  is_pasture)
  call send_cellfrac_data(id_residualFrac, is_residual)
  call send_cellfrac_data(id_treeFrac,     is_tree)
  call send_cellfrac_data(id_grassFrac,    is_ntrlgrass)
  call send_cellfrac_data(id_grassFracC3,  is_ntrlgrass_C3)
  call send_cellfrac_data(id_grassFracC4,  is_ntrlgrass_C4)
  call send_cellfrac_data(id_c3pftFrac,    is_C3)
  call send_cellfrac_data(id_c4pftFrac,    is_C4)
  ! LUMIP land use fractions
  call send_cellfrac_data(id_fracLut_psl,  is_psl,     scale=1.0)
  call send_cellfrac_data(id_fracLut_crp,  is_crop,    scale=1.0)
  call send_cellfrac_data(id_fracLut_pst,  is_pasture, scale=1.0)
  call send_cellfrac_data(id_fracLut_urb,  is_urban,   scale=1.0)
  ! deallocate override buffer
  deallocate(phot_co2_data)

  call mpp_clock_end(landFastClock)
  call mpp_clock_end(landClock)
end subroutine update_land_model_fast


! ============================================================================
subroutine update_land_model_fast_0d(tile, l, k, land2cplr, &
   precip_l, precip_s, atmos_T, &
   Ha0, DHaDTc, tr_flux, dfdtr, &
   ISa_dn_dir, ISa_dn_dif, ILa_dn, &
   ustar, p_surf, drag_q, &
   phot_co2_overridden, phot_co2_data, &
   runoff, runoff_c &
   )
  type (land_tile_type), pointer :: tile
  type(land_data_type), intent(inout) :: land2cplr
  integer, intent(in) :: l,k ! coordinates
  real, intent(in) :: &
       precip_l, precip_s, & ! liquid and solid precipitation, kg/(m2 s)
       atmos_T, &        ! incoming precipitation temperature (despite its name), deg K
       Ha0, DHaDTc, &    ! sensible heat flux from the canopy air to the atmosphere
       tr_flux(:), dfdtr(:), &  ! tracer flux from canopy air to the atmosphere
       ISa_dn_dir(NBANDS), & ! downward direct sw radiation at the top of the canopy
       ISa_dn_dif(NBANDS), & ! downward diffuse sw radiation at the top of the canopy
       ILa_dn,             & ! downward lw radiation at the top of the canopy
       ustar,              & ! friction velocity, m/s
       p_surf,             & ! surface pressure, Pa
       drag_q,             & !
       phot_co2_data         ! data input for the CO2 for photosynthesis

  logical, intent(in):: phot_co2_overridden
  real, intent(inout) :: &
       runoff, &   ! total runoff of H2O
       runoff_c(:) ! runoff of tracers (including ice/snow and heat)


  ! ---- local constants
  ! indices of variables and equations for implicit time stepping solution :
  integer, parameter :: iqc=1, iTc=2, iTv=3, iwl=4, iwf=5

  ! ---- local vars
  real :: A(5,5),B0(5),B1(5),B2(5) ! implicit equation matrix and right-hand side vectors
  real :: A00(5,5),B10(5),B00(5) ! copy of the above, only for debugging
  integer :: indx(5) ! permutation vector
  ! linearization coefficients of various fluxes between components of land
  ! surface scheme
  real :: &
       Ea0,   DEaDqc, &  ! water vapor flux from canopy air to the atmosphere
       fco2_0,Dfco2Dq,&  ! co2 flux from canopy air to the atmosphere
       G0,    DGDTg,  &  ! ground heat flux
       Hv0,   DHvDTv,   DHvDTc, & ! sens heat flux from vegetation
       Et0,   DEtDTv,   DEtDqc,   DEtDwl,   DEtDwf,  & ! transpiration
       Eli0,  DEliDTv,  DEliDqc,  DEliDwl,  DEliDwf, & ! evaporation of intercepted water
       Esi0,  DEsiDTv,  DEsiDqc,  DEsiDwl,  DEsiDwf, & ! sublimation of intercepted snow
       Hg0,   DHgDTg,   DHgDTc, & ! linearization of the sensible heat flux from ground
       Eg0,   DEgDTg,   DEgDqc, DEgDpsig, & ! linearization of evaporation from ground
       flwv0,  DflwvDTg,  DflwvDTv,& ! linearization of net LW radiation to the canopy
       flwg0,  DflwgDTg,  DflwgDTv,& ! linearization of net LW radiation to the canopy
       vegn_drip_l, vegn_drip_s, & ! drip rate of water and snow, respectively, kg/(m2 s)
       qsat,  DqsatDTg, & ! saturated specifuc humidity at the ground and its derivative w.r.t. temperature
       rho, & ! density of canopy air
       vegn_lai

  ! increments of respective variables over time step, results of the implicit
  ! time step:
  real :: delta_qc, delta_Tc, delta_Tv, delta_wl, delta_ws, delta_Tg, delta_psig, delta_co2
  real :: flwg ! updated value of long-wave ground energy balance
  real :: sum0, sum1

  real :: &
       grnd_T, gT, & ! ground temperature and its value used for sensible heat advection
       vegn_T, vT, & ! vegetation (canopy) temperature
       cana_T, cT, & ! canopy air temperature
       evap_T, eT, & ! temperature assigned to vapor going between land and atmosphere
       soil_uptake_T, & ! average temperature of water taken up by the vegetation
       vegn_Wl,  vegn_Ws, & ! water and snow mass of the canopy
       vegn_ifrac, & ! intercepted fraction of liquid or frozen precipitation
       vegn_hcap,      & ! vegetation heat capacity, including intercepted water and snow
       vegn_fco2, & ! co2 flux from the vegetation, kg CO2/(m2 s)
       hlv_Tv, hlv_Tu, & ! latent heat of vaporization at vegn and uptake temperatures, respectively
       hls_Tv, &         ! latent heat of sublimation at vegn temperature
       grnd_q,         & ! explicit specific humidity at ground surface
       grnd_rh,        & ! explicit relative humidity at ground surface
       grnd_rh_psi,    & ! psi derivative of relative humidity at ground surface
       grnd_flux, &
       soil_beta, &
       RSv(NBANDS), & ! net short-wave radiation balance of the canopy, W/m2
       con_g_h, con_g_v, & ! turbulent cond. between ground and canopy air, for heat and vapor respectively
       con_v_h, con_v_v, & ! turbulent cond. between canopy and canopy air, for heat and vapor respectively
       stomatal_cond, & ! integral stomatal conductance of canopy (that is, multiplied by LAI), for water vapor, m/s
       snow_area, &
       cana_q, & ! specific humidity of canopy air
       cana_co2, & ! co2 moist mixing ratio in canopy air, kg CO2/kg wet air
       cana_co2_mol, & ! co2 dry mixing ratio in canopy air, mol CO2/mol dry air
       fswg, evapg, sensg, &
       subs_G, subs_G2, Mg_imp, snow_G_Z, snow_G_TZ, &
       snow_avrg_T, delta_T_snow,  & ! vertically-average snow temperature and it's change due to s
       vegn_ovfl_l,  vegn_ovfl_s,  & ! overflow of liquid and solid water from the canopy
       vegn_ovfl_Hl, vegn_ovfl_Hs, & ! heat flux from canopy due to overflow
       delta_fprec, & ! correction of below-canopy solid precip in case it's average T > tfreeze

       hprec,              & ! sensible heat flux carried by precipitation
       hevap,              & ! sensible heat flux carried by total evapotranspiration
       land_evap,          & ! total vapor flux from land to atmosphere
       land_sens,          & ! turbulent sensible heat flux from land to atmosphere
       vegn_flw,vegn_sens,snow_sens,snow_levap,snow_fevap,snow_melt,&
       snow_lprec, snow_hlprec,snow_lrunf,vegn_levap,vegn_fevap,vegn_uptk,&
       vegn_fsw, vegn_melt,vegn_lprec,vegn_fprec,vegn_hlprec,vegn_hfprec,&
       precip_T,pT,snow_fsw,snow_flw,snow_frunf,snow_hlrunf,&
       snow_hfrunf,subs_fsw,subs_flw,subs_sens,&
       subs_DT, subs_M_imp, subs_evap, snow_Tbot, snow_Cbot, snow_C, subs_levap,&
       subs_fevap,subs_melt,subs_lrunf,subs_hlrunf, subs_frunf, subs_hfrunf,&
       subs_Ttop,subs_Ctop, new_T, &
       subs_tr_runf(n_river_tracers) ! runoff of tracers from soil
  real :: soil_water_supply ! supply of water to roots, per unit active root biomass, kg/m2
  real :: snow_T, snow_rh
  integer :: ii, jj, i, j ! indices for debug output
  integer :: ierr
  integer :: tr, n ! tracer indices
  logical :: conserve_glacier_mass, snow, redo_leaf_water
  integer :: canopy_water_step
  real :: subs_z0m, subs_z0s, snow_z0m, snow_z0s, grnd_z0s
  ! variables for conservation checks
  real :: lmass0, fmass0, heat0, cmass0, v0
  real :: lmass1, fmass1, heat1, cmass1
  character(*), parameter :: tag = 'update_land_model_fast_0d'
  real :: lswept, fswept, hlswept, hfswept ! amounts of liquid and frozen snow, and corresponding
                                           ! heat swept with tiny snow


  i = lnd%i_index(l)
  j = lnd%j_index(l)

  if(is_watch_point()) then
     write(*,*)
     call log_date('#### update_land_model_fast_0d begins:',lnd%time)
  endif
  ! sanity checks of some input values
  call check_var_range(precip_l, precip_warning_tol, HUGE(1.0), 'land model input', 'precip_l', WARNING)
  call check_var_range(precip_s, precip_warning_tol, HUGE(1.0), 'land model input', 'precip_s', WARNING)
  call check_temp_range(atmos_T,                   'land model input', 'atmos_T')
  call check_var_range(ISa_dn_dir, 0.0, 1360.0,    'land model input', 'sw.down.dir', WARNING)
  call check_var_range(ISa_dn_dif, 0.0, 1360.0,    'land model input', 'sw.down.dif', WARNING)
  call check_var_range(ILa_dn,     0.0, 1360.0,    'land model input', 'lw.down',     WARNING)
  call check_var_range(ustar,      0.0, HUGE(1.0), 'land model input', 'ustar',       WARNING)
  call check_var_range(drag_q,     0.0, HUGE(1.0), 'land model input', 'drag_q',      WARNING)
  call check_var_range(p_surf,     0.0, HUGE(1.0), 'land model input', 'p_surf',      WARNING)
  ! not checking fluxes and their derivatives, since they can be either positive
  ! or negative, and it is hard to determine valid ranges for them.

  Ea0    = tr_flux(isphum) ; DEaDqc  = dfdtr(isphum)
  fco2_0 = tr_flux(ico2)   ; Dfco2Dq = dfdtr(ico2)

  if(do_check_conservation) then
     ! + conservation check, part 1: calculate the pre-transition totals
     call get_tile_water(tile,lmass0,fmass0)
     cmass0 = land_tile_carbon(tile)
     ! - end of conservation check, part 1
  endif

  ! if requested (in snow_nml), sweep tiny snow before calling step_1 subroutines to
  ! avoid numerical issues.
  call sweep_tiny_snow(tile%snow, lswept, fswept, hlswept, hfswept)

  soil_uptake_T = tfreeze ! just to avoid using un-initialized values
  soil_water_supply = 0.0
  if (associated(tile%glac)) then
     call glac_step_1 ( tile%glac, &
          grnd_rh, snow_G_Z, snow_G_TZ, conserve_glacier_mass  )
     grnd_rh_psi = 0
  else if (associated(tile%lake)) then
     call lake_step_1 ( ustar, p_surf, &
          lnd%lat(l), tile%lake, &
          grnd_rh, snow_G_Z, snow_G_TZ)
     grnd_rh_psi = 0
  else if (associated(tile%soil)) then
     call soil_step_1 ( tile%soil, tile%vegn, tile%diag, &
          soil_uptake_T, soil_beta, soil_water_supply, &
          grnd_rh, grnd_rh_psi, snow_G_Z, snow_G_TZ)
  else
     call get_current_point(face=ii)
     call error_mesg('update_land_model_fast','none of the surface tiles exist at ('//&
          trim(string(i))//','//trim(string(j))//','//trim(string(k))//&
          ', face='//trim(string(ii))//')',FATAL)
  endif

  if(do_check_conservation) then
     ! + heat conservation check, part 1; land_tile_heat has to be called after
     !   soil_step_1, because soil dry heat capacity is initialized there
     heat0  = land_tile_heat(tile)
     ! - end of conservation check, part 1
  endif

  call snow_step_1 ( tile%snow, snow_G_Z, snow_G_TZ, snow_rh, snow_area, G0, DGDTg )
  snow = snow_active(tile%snow)
  if (snow) then
     grnd_rh   = snow_rh
     grnd_rh_psi = 0
  endif
  snow_T = tile%snow%T(1)

  call cana_state(tile%cana, cana_T, cana_q, cana_co2)

  if (associated(tile%vegn)) then
     ! Calculate net short-wave radiation input to the vegetation
     RSv    = tile%Sv_dir*ISa_dn_dir + tile%Sv_dif*ISa_dn_dif
     ! calculate roughness of the surface under vegetation
     call soil_roughness(tile%soil, subs_z0s, subs_z0m)
     call snow_roughness(tile%snow, snow_z0s, snow_z0m)
     grnd_z0s = exp( (1-snow_area)*log(subs_z0s) + snow_area*log(snow_z0s))

     ! cana_co2 is moist mass mixing ratio [kg CO2/kg wet air], convert it to dry
     ! volumetric mixing ratio [mol CO2/mol dry air]
     cana_co2_mol = cana_co2*mol_air/mol_CO2/(1-cana_q)
     if (phot_co2_overridden) cana_co2_mol = phot_co2_data
     call vegn_step_1 ( tile%vegn, tile%soil, tile%diag, &
        p_surf, &
        ustar, &
        drag_q, &
        ISa_dn_dir+ISa_dn_dif, RSv, precip_l, precip_s, &
        tile%land_d, tile%land_z0s, tile%land_z0m, grnd_z0s, &
        soil_beta, soil_water_supply,&
        cana_T, cana_q, cana_co2_mol, &
        ! output
        con_g_h, con_g_v, &
        con_v_h, con_v_v, &
        stomatal_cond,  &
        vegn_T, vegn_Wl, vegn_Ws, & ! temperature, water and snow mass on the canopy
        vegn_ifrac, vegn_lai, &
        vegn_drip_l, vegn_drip_s,&
        vegn_hcap, & ! total vegetation heat capacity (including intercepted water/snow)
        Hv0,   DHvDTv,   DHvDTc,            &
        Et0,   DEtDTv,   DEtDqc,   DEtDwl,   DEtDwf,  &
        Eli0,  DEliDTv,  DEliDqc,  DEliDwl,  DEliDwf, &
        Esi0,  DEsiDTv,  DEsiDqc,  DEsiDwl,  DEsiDwf  )
        if (LM2) then
           con_g_h = con_g_h * con_fac_large
           if (snow) then
              con_g_v = con_g_v * con_fac_large
           else
              con_g_v = con_g_v * con_fac_small
           endif
        endif
  else
     RSv    = 0
     con_g_h = con_fac_large ; con_g_v = con_fac_large
     if(associated(tile%glac).and.conserve_glacier_mass.and..not.snow) &
          con_g_v = con_fac_small
     con_v_h = 0.0 ; con_v_v = 0.0; stomatal_cond = 0.0
     vegn_T  = cana_T ; vegn_Wl = 0 ; vegn_Ws = 0
     vegn_ifrac  = 0 ; vegn_lai    = 0
     vegn_drip_l = 0 ; vegn_drip_s = 0
     vegn_hcap = 1.0
     Hv0 =0;  DHvDTv =0;  DHvDTc=0;
     Et0 =0;  DEtDTv =0;  DEtDqc=0;   DEtDwl=0;   DEtDwf=0
     Eli0=0;  DEliDTv=0;  DEliDqc=0;  DEliDwl=0;  DEliDwf=0
     Esi0=0;  DEsiDTv=0;  DEsiDqc=0;  DEsiDwl=0;  DEsiDwf=0
  endif
  ! calculate net shortwave for ground and canopy
  fswg     = SUM(tile%Sg_dir*ISa_dn_dir + tile%Sg_dif*ISa_dn_dif)
  vegn_fsw = SUM(RSv)

  grnd_T = land_grnd_T(tile)

  ! + cana_step_1
  rho      =  p_surf/(rdgas*tile%cana%T*(1+d608*tile%cana%tr(isphum)))
  Hg0      =  rho*cp_air*con_g_h*(grnd_T - tile%cana%T)
  DHgDTg   =  rho*cp_air*con_g_h
  DHgDTc   = -rho*cp_air*con_g_h

  if (associated(tile%lake).and.lake_rh_feedback == LAKE_RH_BETA) then
     ! adjust the conductance so that the water vapor flux to the atmopshere is
     ! E = beta*rho*CD*|v|*(qsat - qatm), with beta=grnd_rh
     grnd_rh  = min(grnd_rh, 1-1.0e-6) ! to protect from infinite conductance
     con_g_v  = grnd_rh/(1-grnd_rh)*DEaDqc/rho
     grnd_rh  = 1.0
  endif
  call check_temp_range(grnd_T,'update_land_model_fast_0d','grnd_T')
  call qscomp(grnd_T, p_surf, qsat, DqsatDTg)
  grnd_q   = grnd_rh * qsat
  Eg0      =  rho*con_g_v*(grnd_q  - tile%cana%tr(isphum))
  DEgDTg   =  rho*con_g_v*DqsatDTg*grnd_rh
  DEgDqc   = -rho*con_g_v
  DEgDpsig =  rho*con_g_v*qsat*grnd_rh_psi
  ! - cana_step_1

! [X.X] using long-wave optical properties, calculate the explicit long-wave
!       radiative balances and their derivatives w.r.t. temperatures
  call land_lw_balance(ILa_dn, vegn_T, grnd_T, &
        tile%vegn_tran_lw, tile%vegn_refl_lw, tile%surf_refl_lw, &
        flwv0, flwg0, DflwvDTv, DflwvDTg, DflwgDTv, DflwgDTg )

  if (use_atmos_T_for_precip_T) then
    precip_T = atmos_T
  else
    precip_T = cana_T
  endif
  if (use_atmos_T_for_evap_T) then
    evap_T = atmos_T
  else
    evap_T = cana_T
  endif
  if (use_old_conservation_equations) then
    hlv_Tv = hlv       - (cpw-clw)*tfreeze + cpw*vegn_T
    hls_Tv = hlv + hlf - (cpw-csw)*tfreeze + cpw*vegn_T
    hlv_Tu = hlv       - (cpw-clw)*tfreeze + cpw*vegn_T - clw*soil_uptake_T
    pT = precip_T
    cT = cana_T
    eT = evap_T
    gT = grnd_T
    vT = vegn_T
  else
    hlv_Tv = hlv    + cpw*(vegn_T-tfreeze)
    hls_Tv = hlf    + hlv_Tv
    hlv_Tu = hlv_Tv - clw*(soil_uptake_T-tfreeze)
    pT = precip_T-tfreeze
    cT = cana_T-tfreeze
    eT = evap_T-tfreeze
    gT = grnd_T-tfreeze
    vT = vegn_T-tfreeze
  endif

  do canopy_water_step = 1,2
     if(is_watch_point()) then
        write(*,*)'#### input data for the matrix ####'
        __DEBUG1__(delta_time)
        __DEBUG4__(vegn_T,vT,vegn_Wl,vegn_Ws)
        __DEBUG3__(grnd_T,gT,grnd_rh)
        __DEBUG3__(cana_T,cT,cana_q)
        __DEBUG2__(evap_T,eT)
        __DEBUG4__(precip_l, vegn_drip_l, pT, precip_T)
        __DEBUG2__(precip_s, vegn_drip_s)
        __DEBUG2__(vegn_ifrac, vegn_lai)
        __DEBUG1__(ILa_dn)
        __DEBUG2__(ISa_dn_dir(1),ISa_dn_dir(2))
        __DEBUG2__(ISa_dn_dif(1),ISa_dn_dif(2))
        __DEBUG2__(fswg, vegn_fsw)
        __DEBUG1__(vegn_hcap)
        __DEBUG3__(hlv_Tv, hlv_Tu, hls_Tv)
        __DEBUG2__(G0, DGDTg)
        __DEBUG2__(Ha0, DHaDTc)
        __DEBUG2__(Ea0, DEaDqc)
        __DEBUG3__(Hv0, DHvDTv, DHvDTc)
        __DEBUG5__(Et0,  DEtDTv,  DEtDqc,  DEtDwl,  DEtDwf)
        __DEBUG5__(Eli0, DEliDTv, DEliDqc, DEliDwl, DEliDwf)
        __DEBUG5__(Esi0, DEsiDTv, DEsiDqc, DEsiDwl, DEsiDwf)
        __DEBUG3__(Hg0, DHgDTg, DHgDTc)
        __DEBUG3__(Eg0, DEgDTg, DEgDqc)
        __DEBUG3__(flwv0, DflwvDTg, DflwvDTv)
        __DEBUG3__(flwg0, DflwgDTg, DflwgDTv)
        __DEBUG2__(tile%e_res_1,tile%e_res_2)
     endif

! [X.1] form the system of equations for implicit scheme, such that A*X = B1*delta_Tg+B2*delta_psig+B0
! [X.1.1] equation of canopy air mass balance
     A(iqc,iqc) = canopy_air_mass/delta_time-DEtDqc-DEliDqc-DEsiDqc-DEgDqc+DEaDqc
     A(iqc,iTc) = 0
     A(iqc,iTv) = -DEtDTv-DEliDTv-DEsiDTv
     A(iqc,iwl) = -DEtDwl-DEliDwl-DEsiDwl
     A(iqc,iwf) = -DEtDwf-DEliDwf-DEsiDwf
     B0(iqc)  = Esi0+Eli0+Et0+Eg0-Ea0
     B1(iqc)  = DEgDTg
     B2(iqc)  = DEgDpsig
! [X.1.2] equation of canopy air energy balance
#ifdef USE_DRY_CANA_MASS
     A(iTc,iqc) = canopy_air_mass*cpw*cT/delta_time &
#else
     A(iTc,iqc) = canopy_air_mass*(cpw-cp_air)*cT/delta_time &
#endif
       - cpw*vT*(DEtDqc+DEliDqc+DEsiDqc) - cpw*gT*DEgDqc + cpw*eT*DEaDqc
#ifdef USE_DRY_CANA_MASS
     A(iTc,iTc) = canopy_air_mass*cp_air/delta_time-DHvDTc-DHgDTc+DHaDTc
#else
     A(iTc,iTc) = canopy_air_mass*(cp_air+cana_q*(cpw-cp_air))/delta_time-DHvDTc-DHgDTc+DHaDTc
#endif
     A(iTc,iTv) = -DHvDTv-cpw*vT*(DEtDTv+DEliDTv+DEsiDTv)
     A(iTc,iwl) =        -cpw*vT*(DEtDwl+DEliDwl+DEsiDwl)
     A(iTc,iwf) =        -cpw*vT*(DEtDwf+DEliDwf+DEsiDwf)
     B0(iTc)  = Hv0 + Hg0 - Ha0 + cpw*(vT*(Et0+Eli0+Esi0)+gT*Eg0-eT*Ea0) - tile%e_res_1
     B1(iTc)  = DHgDTg + cpw*gT*DEgDTg
     B2(iTc)  =          cpw*gT*DEgDpsig
! [X.1.3] equation of canopy energy balance
     A(iTv,iqc) = hlv_Tu*DEtDqc + hlv_Tv*DEliDqc + hls_Tv*DEsiDqc
     A(iTv,iTc) = DHvDTc
     A(iTv,iTv) = vegn_hcap/delta_time-DflwvDTv + DHvDTv + &
          hlv_Tu*DEtDTv + hlv_Tv*DEliDTv + hls_Tv*DEsiDTv + clw*vegn_drip_l + csw*vegn_drip_s
     A(iTv,iwl) = clw*vT/delta_time + hlv_Tu*DEtDwl + hlv_Tv*DEliDwl + hls_Tv*DEsiDwl
     A(iTv,iwf) = csw*vT/delta_time + hlv_Tu*DEtDwf + hlv_Tv*DEliDwf + hls_Tv*DEsiDwf
     B0(iTv)  = vegn_fsw + flwv0 - Hv0 - hlv_Tu*Et0 - Hlv_Tv*Eli0 - hls_Tv*Esi0 &
          + clw*precip_l*vegn_ifrac*pT + csw*precip_s*vegn_ifrac*pT &
          - clw*vegn_drip_l*vT - csw*vegn_drip_s*vT - tile%e_res_2
     B1(iTv)  = DflwvDTg
     B2(iTv)  = 0
! [X.1.4] equation of intercepted liquid water mass balance
     A(iwl,iqc) = DEliDqc
     A(iwl,iTc) = 0
     A(iwl,iTv) = DEliDTv
     A(iwl,iwl) = 1.0/delta_time + DEliDwl
     A(iwl,iwf) = DEliDwf
     B0(iwl)  = -Eli0 + precip_l*vegn_ifrac - vegn_drip_l
     B1(iwl)  = 0
     B2(iwl)  = 0
! [X.1.5] equation of intercepted frozen water mass balance
     A(iwf,iqc) = DEsiDqc
     A(iwf,iTc) = 0
     A(iwf,iTv) = DEsiDTv
     A(iwf,iwl) = DEsiDwl
     A(iwf,iwf) = 1.0/delta_time + DEsiDwf
     B0(iwf)  = -Esi0 + precip_s*vegn_ifrac - vegn_drip_s
     B1(iwf)  = 0
     B2(iwf)  = 0
! [X.1.6] if LAI becomes zero (and, therefore, all fluxes from vegetation and their
! derivatives must be zero too) we get a degenerate case. Still, the drip may be non-zero
! because some water may remain from before leaf drop, and non-zero energy residual can be
! carried over from the previous time step.
! To prevent temperature from going haywire in those cases, we simply replace the equations
! of canopy energy and mass balance with the following:
! vegn_T + delta_Tv = cana_T + delta_Tc
! delta_Wl = -vegn_drip_l*delta_time
! delta_Ws = -vegn_drip_s*delta_time
! the residual vegn_Wl and vegn_Ws, if any, are taken care of by the overflow calculations
     if(vegn_hcap==0) then
        ! vegn_T + delta_Tv = cana_T + delta_Tc
        A(iTv,:)   = 0
        A(iTv,iTc) = -1
        A(iTv,iTv) = +1
        B0(iTv) = cana_T - vegn_T
        B1(iTv) = 0
        ! delta_Wl = -vegn_drip_l*delta_time
        A(iwl,:)   = 0
        A(iwl,iwl) = 1
        B0(iwl) = -vegn_drip_l*delta_time
        B1(iwl) = 0
        ! delta_Ws = -vegn_drip_s*delta_time
        A(iwf,:)   = 0
        A(iwf,iwf) = 1
        B0(iwf) = -vegn_drip_s*delta_time
        B1(iwf) = 0
     endif

     if(is_watch_point()) then
        write(*,*)'#### A, B0, B1, B2 ####'
        do ii = 1, size(A,1)
           write(*,'(99g23.16)')(A(ii,jj),jj=1,size(A,2)),B0(ii),B1(ii),B2(ii)
        enddo
     endif

     A00 = A
     B00 = B0
     B10 = B1

! [X.2] solve the system for free terms and delta_Tg and delta_psig terms, getting
!       linear equation for delta_Tg and delta_psig
     call ludcmp(A,indx, ierr)
     if (ierr/=0)&
          write(*,*) 'Matrix is singular',i,j,k
     call lubksb(A,indx,B0)
     call lubksb(A,indx,B1)
     call lubksb(A,indx,B2)

     if(is_watch_point()) then
        write(*,*)'#### solution: B0, B1, B2 ####'
        do ii = 1, size(A,1)
           __DEBUG3__(B0(ii),B1(ii),B2(ii))
        enddo
!!$        write(*,*)'#### solution check ####'
!!$        do ii = 1, size(A,1)
!!$           sum0 = 0; sum1 = 0;
!!$           do jj = 1, size(A,2)
!!$              sum0 = sum0 + A00(ii,jj)*B0(jj)
!!$              sum1 = sum1 + A00(ii,jj)*B1(jj)
!!$           enddo
!!$           write(*,'(99g)')sum0-B00(ii),sum1-B10(ii)
!!$        enddo
     endif
     ! the result of this solution is a set of expressions for delta_xx in terms
     ! of delta_Tg and delta_psig:
     ! delta_xx(i) = B0(i) + B1(i)*delta_Tg + B2(i)*delta_psig. Note that A, B0, B1 and B2
     ! are destroyed in the process: A is replaced with LU-decomposition, and
     ! B0, B1, B2 are replaced with solutions
     !
     ! to find delta_Tg and delta_psig, solve (non-linear) equation of energy balance at
     ! the surface.
     call land_surface_energy_balance( tile, &
          fswg, &
          flwg0 + b0(iTv)*DflwgDTv, DflwgDTg + b1(iTv)*DflwgDTv, b2(iTv)*DflwgDTv, &
          Hg0   + b0(iTc)*DHgDTc,   DHgDTg   + b1(iTc)*DHgDTc,   b2(iTc)*DHgDTc,   &
          Eg0   + b0(iqc)*DEgDqc,   DEgDTg   + b1(iqc)*DEgDqc,   DEgDpsig + b2(iqc)*DEgDqc,   &
          G0,                       DGDTg, &
          ! output
          delta_Tg, delta_psig, Mg_imp )

! [X.5] calculate final value of other tendencies
     delta_qc = B0(iqc) + B1(iqc)*delta_Tg + B2(iqc)*delta_psig
     delta_Tc = B0(iTc) + B1(iTc)*delta_Tg + B2(iTc)*delta_psig
     delta_Tv = B0(iTv) + B1(iTv)*delta_Tg + B2(iTv)*delta_psig
     delta_wl = B0(iwl) + B1(iwl)*delta_Tg + B2(iwl)*delta_psig
     delta_ws = B0(iwf) + B1(iwf)*delta_Tg + B2(iwf)*delta_psig

! [X.6] calculate updated values of energy balance components used in further
!       calculations
     flwg       = flwg0 + DflwgDTg*delta_Tg + DflwgDTv*delta_Tv
     evapg      = Eg0   + DEgDTg*delta_Tg   + DEgDpsig*delta_psig + DEgDqc*delta_qc
     sensg      = Hg0   + DHgDTg*delta_Tg   + DHgDTc*delta_Tc
     grnd_flux  = G0    + DGDTg*delta_Tg
     vegn_sens  = Hv0   + DHvDTv*delta_Tv   + DHvDTc*delta_Tc
     vegn_levap = Eli0  + DEliDTv*delta_Tv  + DEliDqc*delta_qc + DEliDwl*delta_wl + DEliDwf*delta_ws
     vegn_fevap = Esi0  + DEsiDTv*delta_Tv  + DEsiDqc*delta_qc + DEsiDwl*delta_wl + DEsiDwf*delta_ws
     vegn_uptk  = Et0   + DEtDTv*delta_Tv   + DEtDqc*delta_qc  + DEtDwl*delta_wl  + DEtDwf*delta_ws
     vegn_flw   = flwv0 + DflwvDTv*delta_Tv + DflwvDTg*delta_Tg
     land_evap  = Ea0   + DEaDqc*delta_qc
     land_sens  = Ha0   + DHaDTc*delta_Tc
! [X.7] calculate energy residuals due to cross-product of time tendencies
#ifdef USE_DRY_CANA_MASS
     tile%e_res_1 = canopy_air_mass*cpw*delta_qc*delta_Tc/delta_time
#else
     tile%e_res_1 = canopy_air_mass*(cpw-cp_air)*delta_qc*delta_Tc/delta_time
#endif
     tile%e_res_2 = delta_Tv*(clw*delta_Wl+csw*delta_Ws)/delta_time
! calculate the final value upward long-wave radiation flux from the land, to be
! returned to the flux exchange.
     tile%lwup = ILa_dn - vegn_flw - flwg

     if(is_watch_point())then
        write(*,*)'#### ground balance'
        __DEBUG2__(fswg,flwg)
        __DEBUG2__(sensg,evapg*hlv)
        ! note that evapg*hlv is only approximate latent heat flux
        __DEBUG1__(grnd_flux)
        __DEBUG1__(Mg_imp)
        write(*,*)'#### implicit time steps'
        __DEBUG3__(delta_Tg, grnd_T,  grnd_T+delta_Tg )
        __DEBUG1__(delta_psig                         )
        __DEBUG3__(delta_qc, cana_q,  cana_q+delta_qc )
        __DEBUG3__(delta_Tc, cana_T,  cana_T+delta_Tc )
        __DEBUG3__(delta_Tv, vegn_T,  vegn_T+delta_Tv )
        __DEBUG3__(delta_wl, vegn_Wl, vegn_Wl+delta_wl)
        __DEBUG3__(delta_ws, vegn_Ws, vegn_Ws+delta_ws)
        __DEBUG2__(tile%e_res_1, tile%e_res_2)
        write(*,*)'#### resulting fluxes'
        __DEBUG4__(flwg, evapg, sensg, grnd_flux)
        __DEBUG3__(vegn_levap,vegn_fevap,vegn_uptk)
        __DEBUG2__(vegn_sens,vegn_flw)
        __DEBUG1__(Ea0+DEaDqc*delta_qc)
        __DEBUG2__(tile%cana%tr(isphum),cana_q)
     endif

     if (.not.prohibit_negative_canopy_water) exit ! do no corrections
     redo_leaf_water = .FALSE.
     if (vegn_Wl+delta_wl<0) then
        redo_leaf_water = .TRUE.
        Eli0 = vegn_Wl/delta_time + precip_l*vegn_ifrac - vegn_drip_l
        DEliDTv = 0.0;  DEliDqc = 0.0
        DEliDwl = 0.0;  DEliDwf = 0.0
     endif
     if (vegn_Ws+delta_ws<0) then
        redo_leaf_water = .TRUE.
        Esi0 = vegn_Ws/delta_time + precip_s*vegn_ifrac - vegn_drip_s
        DEsiDTv = 0.0;  DEsiDqc = 0.0
        DEsiDwl = 0.0;  DEsiDwf = 0.0
     endif
     if (.not.redo_leaf_water) exit ! from loop
  enddo ! canopy_water_step

  call cana_step_2 ( tile%cana, delta_Tc, delta_qc )

  if(associated(tile%vegn)) then
     call vegn_step_2 ( tile%vegn, tile%diag, &
          delta_Tv, delta_wl, delta_ws, &
          vegn_melt,  &
          vegn_ovfl_l,   vegn_ovfl_s, &
          vegn_ovfl_Hl, vegn_ovfl_Hs )
     ! calculate total amount of liquid and solid precipitation below the canopy
     vegn_lprec  = (1-vegn_ifrac)*precip_l + vegn_drip_l + vegn_ovfl_l
     vegn_fprec  = (1-vegn_ifrac)*precip_s + vegn_drip_s + vegn_ovfl_s
     ! calculate heat carried by liquid and solid precipitation below the canopy
     vegn_hlprec = clw*((1-vegn_ifrac)*precip_l*(precip_T-tfreeze) &
                      + vegn_drip_l*(vegn_T+delta_Tv-tfreeze)) &
                      + vegn_ovfl_Hl
     vegn_hfprec = csw*((1-vegn_ifrac)*precip_s*(precip_T-tfreeze) &
                      + vegn_drip_s*(vegn_T+delta_Tv-tfreeze)) &
                      + vegn_ovfl_Hs
     ! make sure the temperature of the snow falling below canopy is below freezing
     ! this correction was introduced in an attempt to fix the problem with fictitious
     ! heat accumulating in near-zero-mass snow; however it does not seem to make a
     ! difference.
     if(vegn_hfprec>0)then
        ! solid precipitation from vegetation carries positive energy -- we cannot have
        ! that, because that would bring snow T above tfreeze, so convert excess to
        ! liquid
        delta_fprec = min(vegn_fprec,vegn_hfprec/hlf)
        vegn_fprec = vegn_fprec - delta_fprec
        vegn_lprec = vegn_lprec + delta_fprec
        vegn_hfprec = vegn_hfprec - hlf*delta_fprec
        ! we do not need to correct the vegn_hlprec since the temperature of additional
        ! liquid precip is tfreeze, and therefore its contribution to vegn_hlprec is
        ! exactly zero
     endif
     ! possibly we need to correct for the opposite situation: negative energy carried
     ! by liquid precipitation.
  else
     vegn_lprec  = precip_l
     vegn_fprec  = precip_s
     vegn_hlprec = precip_l*clw*(precip_T-tfreeze)
     vegn_hfprec = precip_s*csw*(precip_T-tfreeze)
     ! the fields below are only used in diagnostics
     vegn_melt   = 0
     vegn_fsw    = 0
  endif

  call snow_step_2 ( tile%snow, &
       vegn_lprec, vegn_fprec, vegn_hlprec, vegn_hfprec, &
       delta_Tg, Mg_imp, evapg, fswg, flwg, sensg, &
       use_tfreeze_in_grnd_latent, &
       ! output:
       subs_DT, subs_M_imp, subs_evap, subs_fsw, subs_flw, subs_sens, &
       snow_fsw, snow_flw, snow_sens, &
       snow_levap, snow_fevap, snow_melt, &
       snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, &
       snow_hlrunf, snow_hfrunf, snow_Tbot, snow_Cbot, snow_C, snow_avrg_T )
       snow_lrunf  = snow_lrunf  + lswept/delta_time
       snow_frunf  = snow_frunf  + fswept/delta_time
       snow_hlrunf = snow_hlrunf + hlswept/delta_time
       snow_hfrunf = snow_hfrunf + hfswept/delta_time
  if(is_watch_point()) then
     write(*,*) 'subs_M_imp', subs_M_imp
  endif

  if (snow) then
     subs_G = snow_G_Z+snow_G_TZ*subs_DT
  else
     subs_G = 0
  endif

  if (associated(tile%glac)) then
     call glac_step_2 ( tile%glac, tile%diag, snow_lprec, snow_hlprec, &
          subs_DT, subs_M_imp, subs_evap, &
          subs_levap, subs_fevap, &
          subs_melt, subs_lrunf, subs_hlrunf, subs_Ttop, subs_Ctop )
     subs_frunf = 0.
     subs_hfrunf = 0.
     subs_tr_runf = 0.0
  else if (associated(tile%lake)) then
     call lake_step_2 ( tile%lake, tile%diag, snow_lprec, snow_hlprec, &
          subs_DT, subs_M_imp, subs_evap, &
          use_tfreeze_in_grnd_latent, subs_levap, subs_fevap, &
          subs_melt, subs_Ttop, subs_Ctop )
     subs_lrunf = 0.
     subs_hlrunf = 0.
     subs_frunf = 0.
     subs_hfrunf = 0.
     subs_tr_runf = 0.0
  else if (associated(tile%soil)) then
     call soil_step_2 ( tile%soil, tile%vegn, tile%diag, snow_lprec, snow_hlprec, &
          vegn_uptk, subs_DT, subs_M_imp, subs_evap, &
          use_tfreeze_in_grnd_latent, &
          ! output:
          subs_levap, subs_fevap, &
          subs_melt, subs_lrunf, subs_hlrunf, subs_Ttop, subs_Ctop, &
          subs_frunf, subs_hfrunf, subs_tr_runf)
  endif

! TEMP FIX: MAIN PROG SHOULD NOT TOUCH CONTENTS OF PROG VARS. ******
! ALSO, DIAGNOSTICS IN COMPONENT MODULES SHOULD _FOLLOW_ THIS ADJUSTMENT******
  if (LM2) then
     tile%snow%T = subs_Ttop
     subs_G2 = 0.
  else
     if (sum(tile%snow%ws(:))>0)then
        new_T = (subs_Ctop*subs_Ttop +snow_Cbot*snow_Tbot) &
                        / (subs_Ctop+snow_Cbot)
        tile%snow%T(size(tile%snow%T)) = new_T
        if(associated(tile%glac)) tile%glac%T(1) = new_T
        if(associated(tile%lake)) tile%lake%T(1) = new_T
        if(associated(tile%soil)) tile%soil%T(1) = new_T
        subs_G2 = subs_Ctop*(new_T-subs_Ttop)/delta_time
     else
        if(tau_snow_T_adj>=0) then
           delta_T_snow = subs_Ctop*(subs_Ttop-snow_avrg_T)/&
                (subs_Ctop*tau_snow_T_adj/delta_time+subs_Ctop+snow_C)
           tile%snow%T(:) = snow_avrg_T + delta_T_snow

           new_T = subs_Ttop-snow_C/subs_Ctop*delta_T_snow
           if(associated(tile%glac)) tile%glac%T(1) = new_T
           if(associated(tile%lake)) tile%lake%T(1) = new_T
           if(associated(tile%soil)) tile%soil%T(1) = new_T
           subs_G2 = subs_Ctop*(new_T-subs_Ttop)/delta_time
        else
           subs_G2 = 0.
        endif
     endif
  endif

  vegn_fco2 = 0
  if (associated(tile%vegn)) then
     ! do the calculations that require updated land surface prognostic variables
     call vegn_step_3 (tile%vegn, tile%soil, tile%cana%T, precip_l+precip_s, &
          vegn_fco2, tile%diag)
     ! if vegn is present, then soil must be too
     call soil_step_3(tile%soil, tile%diag)
  endif
  ! update co2 concentration in the canopy air. It would be more consistent to do that
  ! in the same place and fashion as the rest of prognostic variables: that is, have the
  ! vegn_step_1 (and perhaps other *_step_1 procedures) calculate fluxes and their
  ! derivatives, then solve the linear equation(s), and finally have cana_step_2 update
  ! the concentration.
  if(update_cana_co2) then
     delta_co2 = (vegn_fco2 - fco2_0)/(canopy_air_mass_for_tracers/delta_time+Dfco2Dq)
     tile%cana%tr(ico2) = tile%cana%tr(ico2) + delta_co2
  else
     delta_co2 = 0
  endif
  if(is_watch_point())then
     __DEBUG1__(tile%cana%tr(ico2))
     __DEBUG3__(fco2_0,Dfco2Dq,vegn_fco2)
  endif

  call update_cana_tracers(tile, tr_flux, dfdtr, &
           precip_l, precip_s, p_surf, ustar, con_g_v, con_v_v, stomatal_cond )

  call update_land_bc_fast (tile, l, k, land2cplr)

  ! accumulate runoff variables over the tiles
  runoff      = runoff      + (snow_frunf  + subs_lrunf  + snow_lrunf + subs_frunf)*tile%frac
  do tr = 1,n_river_tracers
     if (tr==i_river_heat) then
        runoff_c(tr) = runoff_c(tr) + (snow_hfrunf + subs_hlrunf + snow_hlrunf + subs_hfrunf)*tile%frac
     else if (tr==i_river_ice) then
        runoff_c(tr) = runoff_c(tr) + (snow_frunf + subs_frunf)*tile%frac
     else
        runoff_c(tr) = runoff_c(tr) + subs_tr_runf(tr) * tile%frac
     endif
  enddo
  hprec = (clw*precip_l+csw*precip_s)*(precip_T-tfreeze)
  hevap = cpw*land_evap*(evap_T-tfreeze)

  if (is_watch_cell()) then
     write(*,*)'Accumulated runoff for watch_cell'
     __DEBUG3__(k, runoff, runoff_c)
  end if

  if(do_check_conservation) then
     ! + conservation check, part 2: calculate totals in final state, and compare
     ! with previous totals
     call get_tile_water(tile,lmass1,fmass1)
     call check_conservation (tag,'water', &
         lmass0+fmass0+(precip_l+precip_s-land_evap-(snow_frunf+subs_lrunf+snow_lrunf))*delta_time, &
         lmass1+fmass1, water_cons_tol)
     v0=lmass0+fmass0+(precip_l+precip_s-land_evap-(snow_frunf+subs_lrunf+snow_lrunf))*delta_time
     call send_tile_data(id_water_cons, (lmass1+fmass1-v0)/delta_time, tile%diag)

     cmass1 = land_tile_carbon(tile)
     call check_conservation (tag,'carbon', &
        cmass0-(fco2_0+Dfco2Dq*delta_co2)*mol_C/mol_CO2*delta_time, &
        cmass1, carbon_cons_tol)
     v0 = cmass0-(fco2_0+Dfco2Dq*delta_co2)*mol_C/mol_CO2*delta_time
     call send_tile_data(id_carbon_cons, (cmass1-v0)/delta_time, tile%diag)

     heat1  = land_tile_heat(tile)
     ! latent heat is missing below, and it's not trivial to add, because there are
     ! multiple components with their own vaporization heat
   !  call check_conservation (tag,'heat content', &
   !      heat0+(hprec-land_sens-hevap &
   !           +sum(ISa_dn_dir*(1-tile%land_refl_dir)+ISa_dn_dif*(1-tile%land_refl_dif)) &
   !           +ILa_dn-tile%lwup &
   !           -(snow_hfrunf + subs_hlrunf + snow_hlrunf) &
   !           )*delta_time, &
   !      heat1, 1e-16)
     ! - end of conservation check, part 2
  endif

  ! ---- diagnostic section ----------------------------------------------
  call send_tile_data(id_frac,    tile%frac,                          tile%diag)
  call send_tile_data(id_ntiles,  1.0,                                tile%diag)
  call send_tile_data(id_precip,  precip_l+precip_s,                  tile%diag)
  call send_tile_data(id_hprec,   hprec,                              tile%diag)
  call send_tile_data(id_lprec,   precip_l,                           tile%diag)
  call send_tile_data(id_lprecv,  precip_l-vegn_lprec,                tile%diag)
  call send_tile_data(id_lprecs,  vegn_lprec-snow_lprec,              tile%diag)
  call send_tile_data(id_lprecg,  snow_lprec,                         tile%diag)
  call send_tile_data(id_hlprec,  clw*precip_l*(precip_T-tfreeze),    tile%diag)
  call send_tile_data(id_hlprecv, clw*precip_l*(precip_T-tfreeze)-vegn_hlprec, &
                                                                      tile%diag)
  call send_tile_data(id_hlprecs, vegn_hlprec-snow_hlprec,            tile%diag)
  call send_tile_data(id_hlprecg, snow_hlprec,                        tile%diag)
  call send_tile_data(id_fprec,   precip_s,                           tile%diag)
  call send_tile_data(id_fprecv,  precip_s-vegn_fprec,                tile%diag)
  call send_tile_data(id_fprecs,  vegn_fprec,                         tile%diag)
  call send_tile_data(id_hfprec,  csw*precip_s*(precip_T-tfreeze),    tile%diag)
  call send_tile_data(id_hfprecv, csw*precip_s*(precip_T-tfreeze)-vegn_hfprec, &
                                                                      tile%diag)
  call send_tile_data(id_hfprecs, vegn_hfprec,                        tile%diag)
  call send_tile_data(id_evap,    land_evap,                          tile%diag)
  call send_tile_data(id_hevap,   hevap,                              tile%diag)
  call send_tile_data(id_levap,   vegn_levap+snow_levap+subs_levap+vegn_uptk, &
                                                                      tile%diag)
  call send_tile_data(id_levapv,  vegn_levap,                         tile%diag)
  call send_tile_data(id_levaps,  snow_levap,                         tile%diag)
  call send_tile_data(id_levapg,  subs_levap,                         tile%diag)
  call send_tile_data(id_hlevap,  cpw*vegn_levap*(vegn_T-tfreeze) &
                                    +cpw*snow_levap*(snow_T-tfreeze) &
                                    +cpw*subs_levap*(grnd_T-tfreeze), tile%diag)
  call send_tile_data(id_hlevapv, cpw*vegn_levap*(vegn_T-tfreeze),    tile%diag)
  call send_tile_data(id_hlevaps, cpw*snow_levap*(snow_T-tfreeze),    tile%diag)
  call send_tile_data(id_hlevapg, cpw*subs_levap*(grnd_T-tfreeze),    tile%diag)
  call send_tile_data(id_fevap,   vegn_fevap+snow_fevap+subs_fevap,   tile%diag)
  call send_tile_data(id_fevapv,  vegn_fevap,                         tile%diag)
  call send_tile_data(id_fevaps,  snow_fevap,                         tile%diag)
  call send_tile_data(id_fevapg,  subs_fevap,                         tile%diag)
  call send_tile_data(id_hfevap,  cpw*vegn_fevap*(vegn_T-tfreeze) &
                                    +cpw*snow_fevap*(snow_T-tfreeze) &
                                    +cpw*subs_fevap*(grnd_T-tfreeze), tile%diag)
  call send_tile_data(id_hfevapv, cpw*vegn_fevap*(vegn_T-tfreeze),    tile%diag)
  call send_tile_data(id_hfevaps, cpw*snow_fevap*(snow_T-tfreeze),    tile%diag)
  call send_tile_data(id_hfevapg, cpw*subs_fevap*(grnd_T-tfreeze),    tile%diag)
  call send_tile_data(id_runf,    snow_lrunf+snow_frunf+subs_lrunf,   tile%diag)
  call send_tile_data(id_hrunf,   snow_hlrunf+snow_hfrunf+subs_hlrunf,tile%diag)
  call send_tile_data(id_lrunf,   snow_lrunf+subs_lrunf,              tile%diag)
  call send_tile_data(id_lrunfs,  snow_lrunf,                         tile%diag)
  call send_tile_data(id_lrunfg,  subs_lrunf,                         tile%diag)
  call send_tile_data(id_hlrunf,  snow_hlrunf+subs_hlrunf,            tile%diag)
  call send_tile_data(id_hlrunfs, snow_hlrunf,                        tile%diag)
  call send_tile_data(id_hlrunfg, subs_hlrunf,                        tile%diag)
  call send_tile_data(id_frunf,   snow_frunf + subs_frunf,            tile%diag)
  call send_tile_data(id_frunfs,  snow_frunf,                         tile%diag)
  call send_tile_data(id_hfrunf,  snow_hfrunf + subs_hfrunf,          tile%diag)
  call send_tile_data(id_hfrunfs, snow_hfrunf,                        tile%diag)
  ! TODO: generalize diagnostic for runoff of tracers
  if (i_river_DOC /= NO_TRACER) &
        call send_tile_data(id_DOCrunf, subs_tr_runf(i_river_DOC),    tile%diag)
  call send_tile_data(id_melt,    vegn_melt+snow_melt+subs_melt,      tile%diag)
  call send_tile_data(id_meltv,   vegn_melt,                          tile%diag)
  call send_tile_data(id_melts,   snow_melt,                          tile%diag)
  call send_tile_data(id_meltg,   subs_melt,                          tile%diag)
  call send_tile_data(id_fsw,     vegn_fsw+snow_fsw+subs_fsw,         tile%diag)
  call send_tile_data(id_fswv,    vegn_fsw,                           tile%diag)
  call send_tile_data(id_fsws,    snow_fsw,                           tile%diag)
  call send_tile_data(id_fswg,    subs_fsw,                           tile%diag)
  call send_tile_data(id_flw,     vegn_flw+snow_flw+subs_flw,         tile%diag)
  call send_tile_data(id_flwv,    vegn_flw,                           tile%diag)
  call send_tile_data(id_flws,    snow_flw,                           tile%diag)
  call send_tile_data(id_flwg,    subs_flw,                           tile%diag)
  call send_tile_data(id_sens,    land_sens,                          tile%diag)
  call send_tile_data(id_sensv,   vegn_sens,                          tile%diag)
  call send_tile_data(id_senss,   snow_sens,                          tile%diag)
  call send_tile_data(id_sensg,   subs_sens,                          tile%diag)
  call send_tile_data(id_e_res_1, tile%e_res_1,                       tile%diag)
  call send_tile_data(id_e_res_2, tile%e_res_2,                       tile%diag)
  call send_tile_data(id_con_g_h, con_g_h,                            tile%diag)
  call send_tile_data(id_transp,  vegn_uptk,                          tile%diag)
  call send_tile_data(id_wroff,   snow_lrunf+subs_lrunf,              tile%diag)
  call send_tile_data(id_sroff,   snow_frunf,                         tile%diag)
  call send_tile_data(id_htransp, cpw*vegn_uptk*(vegn_T-tfreeze),     tile%diag)
  call send_tile_data(id_huptake, clw*vegn_uptk*(soil_uptake_T-tfreeze), &
                                                                      tile%diag)
  call send_tile_data(id_hroff,   snow_hlrunf+subs_hlrunf+snow_hfrunf, &
                                                                      tile%diag)
  call send_tile_data(id_gsnow,   subs_G,                             tile%diag)
  call send_tile_data(id_gequil,  subs_G2,                            tile%diag)
  call send_tile_data(id_grnd_flux, grnd_flux,                        tile%diag)
  call send_tile_data(id_soil_water_supply, soil_water_supply,        tile%diag)
!  if(grnd_E_max.lt.0.5*HUGE(grnd_E_Max)) &
!      call send_tile_data(id_levapg_max, grnd_E_max,                  tile%diag)
  do tr = 1,ntcana
     call send_tile_data(id_cana_tr(tr), tile%cana%tr(tr),            tile%diag)
  enddo
  call send_tile_data(id_qco2_dvmr,&
       tile%cana%tr(ico2)*mol_air/mol_co2/(1-tile%cana%tr(isphum)),   tile%diag)
  call send_tile_data(id_fco2,    vegn_fco2*mol_C/mol_co2,            tile%diag)
  call send_tile_data(id_swdn_dir, ISa_dn_dir,                        tile%diag)
  call send_tile_data(id_swdn_dif, ISa_dn_dif,                        tile%diag)
  call send_tile_data(id_swup_dir, ISa_dn_dir*tile%land_refl_dir,     tile%diag)
  call send_tile_data(id_swup_dif, ISa_dn_dif*tile%land_refl_dif,     tile%diag)
  call send_tile_data(id_lwdn,     ILa_dn,                            tile%diag)
  call send_tile_data(id_subs_emis,1-tile%surf_refl_lw,               tile%diag)

  ! CMOR variables
  call send_tile_data(id_pcp,    precip_l+precip_s,                   tile%diag)
  call send_tile_data(id_prra,   precip_l,                            tile%diag)
  call send_tile_data(id_prveg, (precip_l+precip_s)*vegn_ifrac,       tile%diag)
  ! note that the expression below gives only approximate value of latent heat
  ! flux, since in the model specific heat of vaporization is not equal to hlv;
  ! it depends on temperature and phase state.
  call send_tile_data(id_hflsLut, land_evap*hlv,                      tile%diag)
  call send_tile_data(id_rlusLut, tile%lwup,                          tile%diag)
  if (id_rsusLut>0) call send_tile_data(id_rsusLut, &
    sum(ISa_dn_dir*tile%land_refl_dir+ISa_dn_dif*tile%land_refl_dif), tile%diag)
  call send_tile_data(id_tran,  vegn_uptk,                            tile%diag)
  ! evspsblsoi is evaporation from *soil*, so we send zero from glaciers and lakes;
  ! the result is averaged over the entire land surface, as required by CMIP. evspsblveg
  ! does not need this distinction because it is already zero over glaciers and lakes.
  if (associated(tile%soil)) then
     call send_tile_data(id_evspsblsoi, subs_levap+subs_fevap,        tile%diag)
  else
     call send_tile_data(id_evspsblsoi, 0.0,                          tile%diag)
  endif
  call send_tile_data(id_evspsblveg,  vegn_levap+vegn_fevap,          tile%diag)
  call send_tile_data(id_nbp,    -vegn_fco2*mol_C/mol_co2,            tile%diag)
  call send_tile_data(id_snm, snow_melt,                              tile%diag)
  if (id_cLand > 0) &
      call send_tile_data(id_cLand, land_tile_carbon(tile),           tile%diag)
  if (id_tslsiLut>0) &
      call send_tile_data(id_tslsiLut, (tile%lwup/stefan)**0.25,      tile%diag)
  if(is_tree(tile)) then
     call send_tile_data(id_nwdFracLut,    0.0,                       tile%diag)
  else
     call send_tile_data(id_nwdFracLut,    1.0,                       tile%diag)
  endif

end subroutine update_land_model_fast_0d


! ============================================================================
subroutine update_land_model_slow ( cplr2land, land2cplr )
  type(atmos_land_boundary_type), intent(inout) :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr

  ! ---- local vars
  integer :: i, j, k, l
  type(land_tile_type), pointer :: tile
  type(land_tile_enum_type) :: ce

  call mpp_clock_begin(landClock)
  call mpp_clock_begin(landSlowClock)

  call land_transitions( lnd%time )
  call update_vegn_slow( )
  ! send the accumulated diagnostics to the output
  call dump_tile_diag_fields(land_tile_map, lnd%time)

  ! land_transitions may have changed the number of tiles per grid cell: reallocate
  ! boundary conditions, if necessary
  call realloc_cplr2land( cplr2land )
  call realloc_land2cplr( land2cplr )

  ! set the land mask to FALSE everywhere -- update_land_bc_fast will set it to
  ! true where necessary
  land2cplr%mask = .FALSE.
  land2cplr%tile_size = 0.0

  ! get the current state of the land boundary for the coupler
  ce = first_elmt(land_tile_map, ls=lnd%ls)
  do while(loop_over_tiles(ce,tile,l,k))
     call set_current_point(l,k)
     call update_land_bc_fast (tile, l,k, land2cplr, is_init=.true.)
  enddo

  call update_land_bc_slow( land2cplr )

  call mpp_clock_end(landClock)
  call mpp_clock_end(landSlowClock)
end subroutine update_land_model_slow


! ============================================================================
! solve for surface temperature. ensure that melt does not exceed available
! snow or soil ice (which would create energy-balance problems). also ensure
! that melt and temperature are consistent and that evaporation from liquid
! soil water does not exceed exfiltration rate limit, if one is supplied.
! because the possible combinations of active constraints has multiplied
! greatly, we do not allow phase change for any surface (i.e., soil) at which
! we might apply a constraint on Eg.
subroutine land_surface_energy_balance ( tile, &
     ! surface parameters
     ! components of the ground energy balance linearization. Note that those
     ! are full derivatives, which include the response of the other surface
     ! scheme parameters to the change in ground temperature.
     fswg,            & ! net short-wave
     flwg0, DflwgDTg, DflwgDpsig, & ! net long-wave
     Hg0,   DHgDTg,   DHgDpsig,   & ! sensible heat
     Eg0,   DEgDTg,   DEgDpsig,   & ! latent heat
     G0,    DGDTg,    & ! sub-surface heat
     delta_Tg,        & ! surface temperature change for the time step
     delta_psig,      &
     Mg_imp          )  ! implicit melt, kg/m2

  type(land_tile_type), intent(in) :: tile
  real, intent(in) :: &
     fswg,            & ! net short-wave
     flwg0, DflwgDTg, DflwgDpsig, & ! net long-wave
     Hg0,   DHgDTg,   DHgDpsig,   & ! sensible heat
     Eg0,   DEgDTg,   DEgDpsig,   & ! latent heat
     G0,    DGDTg                   ! sub-surface heat
  real, intent(out) :: &
     delta_Tg,        & ! change in surface temperature
     delta_psig,      & ! change in surface soil-water matric head
     Mg_imp             ! mass of surface ice melted (or water frozen) during the
                        ! time step, kg/m2

  real :: grnd_B     ! surface energy balance
  real :: grnd_DBDTg ! full derivative of grnd_B w.r.t. surface temperature
  real :: grnd_DBDpsig ! full derivative of grnd_B w.r.t. surface soil-water matric head
  real :: grnd_E_force, Eg_trial, Eg_check, determinant
  real :: &
     grnd_T,          & ! ground temperature
     grnd_E_min,      & ! Eg floor of 0 if condensation is prohibited
     grnd_E_max,      & ! exfiltration rate limit, kg/(m2 s)
     grnd_liq, grnd_ice, & ! amount of water available for freeze or melt on the surface, kg/m2
     grnd_latent,     & ! specific heat of vaporization for the ground
     grnd_Tf,         & ! ground freezing temperature
     grnd_subl

  grnd_T = land_grnd_T(tile)
  call land_sfc_water(tile, grnd_liq, grnd_ice, grnd_subl, grnd_Tf)
  if (use_tfreeze_in_grnd_latent) then
    grnd_latent = hlv + hlf*grnd_subl
  else
    grnd_latent = hlv + (cpw-clw)*(grnd_T-tfreeze) &
               + (hlf + (clw-csw)*(grnd_T-tfreeze)) * grnd_subl
  endif
  call grnd_evap_limits(tile, grnd_E_min, grnd_E_max)

  grnd_B      = fswg + flwg0      - Hg0      - grnd_latent*Eg0      - G0
  grnd_DBDTg  =        DflwgDTg   - DHgDTg   - grnd_latent*DEgDTg   - DGDTg
  grnd_DBDpsig =       DflwgDpsig - DHgDpsig - grnd_latent*DEgDpsig

  if(is_watch_point())then
     write(*,*)'#### ground balance input'
     __DEBUG1__(grnd_T)
     __DEBUG2__(grnd_liq, grnd_ice)
     __DEBUG1__(grnd_latent)
     __DEBUG1__(grnd_Tf)
     __DEBUG1__(grnd_E_min)
     __DEBUG1__(grnd_E_max)
     __DEBUG1__(fswg)
     __DEBUG3__(flwg0, DflwgDTg,DflwgDpsig)
     __DEBUG3__(Hg0,   DHgDTg,  DHgDpsig  )
     __DEBUG3__(Eg0,   DEgDTg,  DEgDpsig  )
     __DEBUG2__(G0,    DGDTg)
     write(*,*)'#### end of ground balance input'
     __DEBUG3__(grnd_B, grnd_DBDTg, grnd_DBDpsig)
  endif

  ! determine the ground temperature change under the assumptions that
  ! (1) no phase change occurs at surface (always true for soil now), and
  ! (2) delta_psig is zero (always true for snow, lake, glacier now)
  delta_Tg = - grnd_B/grnd_DBDTg
  delta_psig = 0
  ! calculate phase change on the ground, if necessary
  if     (grnd_ice>0.and.grnd_T+delta_Tg>grnd_Tf) then ! melt > 0
     Mg_imp =  min(grnd_ice,  grnd_DBDTg*(grnd_Tf-grnd_T-delta_Tg)*delta_time/hlf)
  elseif (grnd_liq>0.and.grnd_T+delta_Tg<grnd_Tf) then ! melt < 0
     Mg_imp = -min(grnd_liq, -grnd_DBDTg*(grnd_Tf-grnd_T-delta_Tg)*delta_time/hlf)
  else
     Mg_imp = 0
  endif
  ! adjust temperature change for the phase change
  delta_Tg = -(grnd_B - Mg_imp*hlf/delta_time)/grnd_DBDTg
  Eg_trial = Eg0 + DEgDTg*delta_Tg

  if(is_watch_point())then
     write(*,*)'#### ground balance solution with psig constant:'
     __DEBUG2__(grnd_B, grnd_DBDTg)
     __DEBUG3__(Mg_imp, delta_Tg, delta_psig)
     __DEBUG3__(grnd_E_min, Eg_trial, grnd_E_max)
  endif

  ! Solution above (assuming no change in psig) is acceptable if
  ! it does not imply unacceptable value of Eg. this is always
  ! true for lake, glacier, snow. If implied Eg is outside prescribed
  ! bounds (a possibility only for soil), then Eg is set at bound (grnd_E_force)
  ! and the required Tg and psig are found. To accomplish this,
  ! we solve the system
  ! grnd_B + grnd_DBDTg*delta_Tg + grnd_DBDpsig*delta_psig = 0
  ! Eg0    +     DEgDTg*delta_Tg +     DEgDpsig*delta_psig = grnd_E_force
  ! There is no need to revisit the solution for phase change, which
  ! is only done explicitly for soil.

  if (Eg_trial.lt.grnd_E_min .or. Eg_trial.gt.grnd_E_max) then
      grnd_E_force = max(grnd_E_min, Eg_trial)
      grnd_E_force = min(grnd_E_max, grnd_E_force)
      determinant = grnd_DBDTg*DEgDpsig-grnd_DBDpsig*DEgDTg
      delta_Tg   = - (DEgDpsig*grnd_B + grnd_DBDpsig*(grnd_E_force-Eg0))/determinant
      delta_psig =   (DEgDTg  *grnd_B + grnd_DBDTg  *(grnd_E_force-Eg0))/determinant
      Eg_check = Eg0 + DEgDTg*delta_Tg + DEgDpsig*delta_psig
      Mg_imp = 0
         if(is_watch_point())then
            write(*,*)'#### trial solution violated Eg limit, new soln:'
            __DEBUG2__(grnd_B,grnd_DBDTg)
            __DEBUG3__(Mg_imp,delta_Tg,delta_psig)
            __DEBUG3__(grnd_E_min, Eg_check, grnd_E_max)
         endif
    endif
end subroutine land_surface_energy_balance

! ============================================================================
subroutine land_lw_balance(lwdn_atm, vegn_T, surf_T, &
  vegn_tran_lw, vegn_refl_lw, surf_refl_lw, &
  flwv, flwg, DflwvDTv, DflwvDTg, DflwgDTv, DflwgDTg )

  real, intent(in) :: lwdn_atm ! downward long-wave radiation on top of the canopy, W/m2
  real, intent(in) :: vegn_T ! canopy temperatures, deg K
  real, intent(in) :: surf_T    ! ground surface temperature, deg K
  real, intent(in) :: vegn_tran_lw ! transmittance of cohort canopies to long-wave radiation
  real, intent(in) :: vegn_refl_lw ! reflectance of cohort canopies to long-wave radiation
  real, intent(in) :: surf_refl_lw ! reflectance of ground surface to long-wave

  real, intent(out) :: flwv ! long-wave balances of canopy, W/m2
  real, intent(out) :: flwg ! ground surface long-wave balance, W/m2
  real, intent(out) :: DflwvDTv ! derivative of canopy balance w.r.t. canopy temperatures, W/(m2 K)
  real, intent(out) :: DflwvDTg ! derivative of canopy balance w.r.t. ground surface temperature, W/(m2 K)
  real, intent(out) :: DflwgDTv ! derivative of ground surface balance w.r.t. canopy temperature, W/(m2 K)
  real, intent(out) :: DflwgDTg ! derivative of ground surface balance w.r.t. ground surface temperature, W/(m2 K)

  real :: vegn_emis_lw, surf_emis_lw ! emissivities of ground and surface
  real :: vegn_emsn,    surf_emsn    ! emission by vegetation and surface, respectively
  real :: denom ! denominator in the LW radiative balance calculations

  vegn_emis_lw = 1 - vegn_refl_lw - vegn_tran_lw
  surf_emis_lw = 1 - surf_refl_lw

  denom = 1-vegn_refl_lw*surf_refl_lw

  vegn_emsn = vegn_emis_lw * stefan * vegn_T**4
  surf_emsn = surf_emis_lw * stefan * surf_T**4

  flwv = lwdn_atm * vegn_emis_lw*(1+vegn_tran_lw*surf_refl_lw/denom) &
       + vegn_emsn * (surf_refl_lw*vegn_emis_lw/denom-2) &
       + surf_emsn * vegn_emis_lw/denom
  DflwvDTg = vegn_emis_lw/denom                 * surf_emis_lw * stefan * 4 * surf_T**3
  DflwvDTv = (surf_refl_lw*vegn_emis_lw/denom-2)* vegn_emis_lw * stefan * 4 * vegn_T**3

  flwg = (lwdn_atm*vegn_tran_lw + vegn_emsn)*(1-surf_refl_lw)/denom &
       - surf_emsn*(1-vegn_refl_lw)/denom
  DflwgDTg = -(1-vegn_refl_lw)/denom * surf_emis_lw * stefan * 4 * surf_T**3
  DflwgDTv =  (1-surf_refl_lw)/denom * vegn_emis_lw * stefan * 4 * vegn_T**3
  if (is_watch_point()) then
     __DEBUG2__(vegn_emis_lw,surf_emis_lw)
     __DEBUG2__(vegn_emsn,surf_emsn)
  endif
end subroutine land_lw_balance

! ============================================================================
! set up constants for linearization of radiative transfer, using information
! provided by soil, snow and vegetation modules.
subroutine land_sw_radiation (lm2, &
     subs_refl_dir, subs_refl_dif, subs_refl_lw, &
     snow_refl_dir, snow_refl_dif, snow_refl_lw, &
     snow_area, &
     vegn_refl_dif, vegn_tran_dif, &
     vegn_refl_dir, vegn_tran_dir, vegn_tran_dir_dir,  &
     ! output
     Sg_dir, Sg_dif, Sv_dir, Sv_dif, &
     land_albedo_dir, land_albedo_dif )

  logical, intent(in) :: lm2
  real, intent(in) :: &
       subs_refl_dir(NBANDS), subs_refl_dif(NBANDS), subs_refl_lw,  & ! sub-snow reflectances for direct, diffuse, and LW radiation respectively
       snow_refl_dir(NBANDS), snow_refl_dif(NBANDS), snow_refl_lw,  & ! snow reflectances for direct, diffuse, and LW radiation respectively
       snow_area, &
       vegn_refl_dif(NBANDS), vegn_tran_dif(NBANDS), & ! vegn reflectance & transmittance for diffuse light
       vegn_refl_dir(NBANDS), vegn_tran_dir(NBANDS), & ! vegn reflectance & transmittance for direct light
       vegn_tran_dir_dir(NBANDS)

  real, intent(out) :: &
     Sg_dir(NBANDS), Sg_dif(NBANDS), & ! fraction of downward short-wave absorbed by ground and snow
     Sv_dir(NBANDS), Sv_dif(NBANDS), & ! fraction of downward short-wave absorbed by vegetation
     land_albedo_dir(NBANDS), land_albedo_dif(NBANDS)

  ! ---- local vars
  real :: &
       grnd_refl_dir(NBANDS), & ! SW reflectances of ground surface (by spectral band)
       grnd_refl_dif(NBANDS), & ! SW reflectances of ground surface (by spectral band)
       grnd_refl_lw             ! LW reflectance of ground surface
  real :: &
       subs_up_from_dir(NBANDS), subs_up_from_dif(NBANDS), &
       subs_dn_dir_from_dir(NBANDS),  subs_dn_dif_from_dif(NBANDS), subs_dn_dif_from_dir(NBANDS)

  grnd_refl_dir = subs_refl_dir + (snow_refl_dir - subs_refl_dir) * snow_area
  grnd_refl_dif = subs_refl_dif + (snow_refl_dif - subs_refl_dif) * snow_area
  grnd_refl_lw  = subs_refl_lw  + (snow_refl_lw  - subs_refl_lw ) * snow_area

  ! ---- shortwave -----------------------------------------------------------
  ! allocation to canopy and ground, based on solution for single
  ! vegetation layer of limited cover. Both ground and vegetation are gray.
  Sv_dir = 0; Sv_dif = 0
  if (lm2) then
     subs_dn_dir_from_dir = vegn_tran_dir_dir
     subs_dn_dif_from_dir = vegn_tran_dir
     subs_dn_dif_from_dif = vegn_tran_dif
     subs_up_from_dir = grnd_refl_dir*subs_dn_dir_from_dir + &
                        grnd_refl_dif*subs_dn_dif_from_dir
     subs_up_from_dif = grnd_refl_dif*subs_dn_dif_from_dif
     land_albedo_dir  = subs_up_from_dir+vegn_refl_dir
     land_albedo_dif  = subs_up_from_dif+vegn_refl_dif
  else
     subs_dn_dir_from_dir = vegn_tran_dir_dir
     subs_dn_dif_from_dir = (vegn_tran_dir + vegn_refl_dif*grnd_refl_dir*vegn_tran_dir_dir)&
                          / (1 - grnd_refl_dif*vegn_refl_dif)
     subs_dn_dif_from_dif = vegn_tran_dif &
                          / (1 - grnd_refl_dif*vegn_refl_dif)
     subs_up_from_dir = grnd_refl_dir*subs_dn_dir_from_dir + &
                        grnd_refl_dif*subs_dn_dif_from_dir
     subs_up_from_dif = grnd_refl_dif*subs_dn_dif_from_dif
     land_albedo_dir  = subs_up_from_dir*vegn_tran_dif + vegn_refl_dir
     land_albedo_dif  = subs_up_from_dif*vegn_tran_dif + vegn_refl_dif
  endif

  Sg_dir = subs_dn_dir_from_dir + subs_dn_dif_from_dir - subs_up_from_dir
  Sg_dif = subs_dn_dif_from_dif - subs_up_from_dif
  Sv_dir = 1 - Sg_dir - land_albedo_dir
  Sv_dif = 1 - Sg_dif - land_albedo_dif
end subroutine land_sw_radiation

! ============================================================================
subroutine update_land_bc_fast (tile, l ,k, land2cplr, is_init)
  type(land_tile_type), intent(inout) :: tile
  integer             , intent(in) :: l,k
  type(land_data_type), intent(inout) :: land2cplr
  logical, optional :: is_init

  ! ---- local vars
  real ::  subs_z0m, subs_z0s, &
           snow_z0s, snow_z0m, &
         snow_area, snow_depth

  real :: subs_refl_dir(NBANDS), subs_refl_dif(NBANDS) ! direct and diffuse albedos
  real :: subs_refl_lw ! reflectance for thermal radiation
  real :: snow_refl_dir(NBANDS), snow_refl_dif(NBANDS) ! direct and diffuse albedos of snow
  real :: snow_refl_lw ! snow reflectance for thermal radiation
  real :: snow_emis ! snow emissivity
  real :: grnd_emis ! ground emissivity
  ! NOTE :  grnd_emis is used only to satisfy xxxx_radiation interfaces; its value is ignored, but
  ! 1-refl is used instead. snow_emis is used in the the vegn_radiation, but it should not be since
  ! properties of intercepted snowpack are, in general, different from the snow on the ground
  real :: snow_area_rad ! "snow area for radiation calculations" -- introduced
                        ! to reproduce lm2 behavior
  real :: vegn_refl_lw, vegn_tran_lw ! reflectance and transmittance of vegetation for thermal radiation
  real :: vegn_refl_dif(NBANDS), vegn_tran_dif(NBANDS) ! reflectance and transmittance of vegetation for diffuse light
  real :: vegn_refl_dir(NBANDS), vegn_tran_dir(NBANDS) ! reflectance and transmittance of vegetation for direct light
  real :: vegn_tran_dir_dir(NBANDS) ! (?)
  real :: &
         vegn_Tv,     &
         vegn_cover,  &
         vegn_height, vegn_lai, vegn_sai, vegn_d_leaf, cana_co2
  logical :: do_update

  real :: cosz    ! cosine of solar zenith angle
  real :: fracday ! daytime fraction of time interval
  real :: rrsun   ! earth-sun distance (r) relative to semi-major axis
                  ! of orbital ellipse (a) : (a/r)**2
  integer :: face ! for debugging
  integer :: tr   ! tracer index
  integer :: i, j
  vegn_Tv = 0

  i = lnd%i_index(l)
  j = lnd%j_index(l)
  do_update = .not.present(is_init)

  if (is_watch_point()) then
     write(*,*)
     call log_date('#### update_land_bc_fast begins:',lnd%time)
  endif

  ! on initialization the albedos are calculated for the current time step ( that is, interval
  ! lnd_sg%time, lnd_sg%time+lnd_sg%dt_fast); in the course of the run this subroutine is called
  ! at the end of time step (but before time is advanced) to calculate the radiative properties
  ! for the _next_ time step
  if (do_update) then
     call diurnal_solar(lnd%lat(l), lnd%lon(l), lnd%time+lnd%dt_fast, &
          cosz, fracday, rrsun, lnd%dt_fast)
  else
     call diurnal_solar(lnd%lat(l), lnd%lon(l), lnd%time, &
          cosz, fracday, rrsun, lnd%dt_fast)
  endif

  if (associated(tile%glac)) then
     call glac_radiation(tile%glac, cosz, subs_refl_dir, subs_refl_dif, subs_refl_lw, grnd_emis)
     call glac_roughness(tile%glac, subs_z0s, subs_z0m )
  else if (associated(tile%lake)) then
     call lake_radiation(tile%lake, cosz, subs_refl_dir, subs_refl_dif, subs_refl_lw, grnd_emis)
     call lake_roughness(tile%lake, subs_z0s, subs_z0m )
  else if (associated(tile%soil)) then
     call soil_radiation(tile%soil, cosz, subs_refl_dir, subs_refl_dif, subs_refl_lw, grnd_emis)
     call soil_roughness(tile%soil, subs_z0s, subs_z0m )
  else
     call get_current_point(face=face)
     call error_mesg('update_land_model_fast','none of the surface tiles exist at ('//&
             trim(string(i))//','//trim(string(j))//','//trim(string(k))//&
             ', face='//trim(string(face))//')',FATAL)
  endif

  call snow_radiation ( tile%snow%T(1), cosz, snow_refl_dir, snow_refl_dif, snow_refl_lw, snow_emis)
  call snow_get_depth_area ( tile%snow, snow_depth, snow_area )
  call snow_roughness ( tile%snow, snow_z0s, snow_z0m )

  if (associated(tile%vegn)) then
     call update_derived_vegn_data(tile%vegn)
     ! USE OF SNOWPACK RAD PROPERTIES FOR INTERCEPTED SNOW IS ERRONEOUS,
     ! NEEDS TO BE CHANGED. TEMPORARY.
     call vegn_radiation ( tile%vegn, cosz, snow_depth, snow_refl_dif, snow_emis, &
                   vegn_refl_dif, vegn_tran_dif, &
                   vegn_refl_dir, vegn_tran_dir, vegn_tran_dir_dir, &
                   vegn_refl_lw, vegn_tran_lw)
     ! (later see if we can remove vegn_cover from c-a-radiation...) TEMPORARY
     call vegn_get_cover ( tile%vegn, snow_depth, vegn_cover)
     ! vegn_properties returns integral properties of the canopy, relevant for the
     ! calculations of the land roughness and displacement
     call vegn_properties ( tile%vegn, vegn_cover, vegn_height, vegn_lai, vegn_sai, vegn_d_leaf)
  else
     ! set radiative properties for null vegetation
     vegn_refl_dif     = 0
     vegn_tran_dif     = 1
     vegn_refl_dir     = 0
     vegn_tran_dir     = 0
     vegn_tran_dir_dir = 1
     vegn_refl_lw      = 0
     vegn_tran_lw      = 1
     ! set cover for null vegetation
     vegn_cover        = 0
     ! set other parameters for null vegetation
     vegn_height       = 0
     vegn_lai          = 0
     vegn_sai          = 0
     vegn_d_leaf       = 0
  endif

  ! store the values of long-wave optical properties to be used in the update_land_model_fast
  tile%surf_refl_lw = subs_refl_lw  + (snow_refl_lw  - subs_refl_lw ) * snow_area
  tile%vegn_refl_lw = vegn_refl_lw
  tile%vegn_tran_lw = vegn_tran_lw


  if(is_watch_point()) then
     write(*,*) '#### update_land_bc_fast ### checkpoint 1 ####'
     __DEBUG3__(cosz, fracday, rrsun)
     __DEBUG2__(vegn_lai,vegn_sai)
     __DEBUG1__(subs_refl_dif)
     __DEBUG1__(subs_refl_dir)
     __DEBUG1__(vegn_refl_dif)
     __DEBUG1__(vegn_tran_dif)
     __DEBUG1__(vegn_refl_dir)
     __DEBUG1__(vegn_tran_dir)
     __DEBUG1__(vegn_tran_dir_dir)
     __DEBUG2__(vegn_refl_lw,vegn_tran_lw)
     write(*,*) '#### update_land_bc_fast ### end of checkpoint 1 ####'
  endif

  snow_area_rad = snow_area
  if (lm2) then
     if(associated(tile%glac)                    ) snow_area_rad = 0
     if(associated(tile%soil).and.vegn_cover>0.01) snow_area_rad = 1
  endif

  call land_sw_radiation( lm2, &
       subs_refl_dir, subs_refl_dif, subs_refl_lw, &
       snow_refl_dir, snow_refl_dif, snow_refl_lw, &
       snow_area_rad,  &
       vegn_refl_dif, vegn_tran_dif, &
       vegn_refl_dir, vegn_tran_dir, vegn_tran_dir_dir, &
       ! output
       tile%Sg_dir, tile%Sg_dif, tile%Sv_dir, tile%Sv_dif, &
       tile%land_refl_dir, tile%land_refl_dif )

  call cana_roughness( lm2, &
     subs_z0m, subs_z0s, &
     snow_z0m, snow_z0s, snow_area, &
     vegn_cover,  vegn_height, vegn_lai, vegn_sai, &
     tile%land_d, tile%land_z0m, tile%land_z0s, &
     associated(tile%lake).or.associated(tile%glac))

  if(is_watch_point()) then
     write(*,*) '#### update_land_bc_fast ### checkpoint 2 ####'
     write(*,*) 'Sg_dir', tile%Sg_dir
     write(*,*) 'Sg_dif', tile%Sg_dif
     write(*,*) 'Sv_dir', tile%Sv_dir
     write(*,*) 'Sv_dif', tile%Sv_dif
     write(*,*) 'land_albedo_dir', tile%land_refl_dir
     write(*,*) 'land_albedo_dif', tile%land_refl_dif
     write(*,*) 'land_z0m', tile%land_z0m
     write(*,*) '#### update_land_bc_fast ### end of checkpoint 2 ####'
  endif

  land2cplr%t_surf         (l,k) = tfreeze
  land2cplr%t_ca           (l,k) = tfreeze
  land2cplr%tr             (l,k, isphum) = 0.0
  land2cplr%albedo         (l,k) = 0.0
  land2cplr%albedo_vis_dir (l,k) = 0.0
  land2cplr%albedo_nir_dir (l,k) = 0.0
  land2cplr%albedo_vis_dif (l,k) = 0.0
  land2cplr%albedo_nir_dif (l,k) = 0.0
  land2cplr%rough_mom      (l,k) = 0.1
  land2cplr%rough_heat     (l,k) = 0.1

  ! Calculate radiative surface temperature. lwup cannot be calculated here
  ! based on the available temperatures because it's a result of the implicit
  ! time step: lwup = lwup0 + DlwupDTg*delta_Tg + ..., so we have to carry it
  ! from the update_land_fast
  ! Consequence: since update_landbc_fast is called once at the very beginning of
  ! every run (before update_land_fast is called) lwup from the previous step
  ! must be stored in the in the restart for reproducibility
  land2cplr%t_surf(l,k) = ( tile%lwup/stefan ) ** 0.25

  ! set the boundary conditions for the flux exchange
  land2cplr%mask           (l,k) = .TRUE.
  land2cplr%tile_size      (l,k) = tile%frac

  land2cplr%t_ca(l,k) = tile%cana%T
  do tr = 1,ntcana
      land2cplr%tr(l,k,tr) = tile%cana%tr(tr)
  enddo

  land2cplr%albedo_vis_dir (l,k) = tile%land_refl_dir(BAND_VIS)
  land2cplr%albedo_nir_dir (l,k) = tile%land_refl_dir(BAND_NIR)
  land2cplr%albedo_vis_dif (l,k) = tile%land_refl_dif(BAND_VIS)
  land2cplr%albedo_nir_dif (l,k) = tile%land_refl_dif(BAND_NIR)
  land2cplr%albedo         (l,k) = SUM(tile%land_refl_dir + tile%land_refl_dif)/4 ! incorrect, replace with proper weighting later
  land2cplr%rough_mom      (l,k) = tile%land_z0m
  land2cplr%rough_heat     (l,k) = tile%land_z0s

  if (use_coast_rough .and. lnd%landfrac(l)<max_coast_frac) then
     land2cplr%rough_mom   (l,k) = coast_rough_mom
     land2cplr%rough_heat  (l,k) = coast_rough_heat
  endif

  if(is_watch_point()) then
     write(*,*)'#### update_land_bc_fast ### output ####'
     write(*,*)'land2cplr%mask',land2cplr%mask(l,k)
     write(*,*)'land2cplr%tile_size',land2cplr%tile_size(l,k)
     write(*,*)'land2cplr%t_surf',land2cplr%t_surf(l,k)
     write(*,*)'land2cplr%t_ca',land2cplr%t_ca(l,k)
     write(*,*)'land2cplr%albedo',land2cplr%albedo(l,k)
     write(*,*)'land2cplr%rough_mom',land2cplr%rough_mom(l,k)
     write(*,*)'land2cplr%rough_heat',land2cplr%rough_heat(l,k)
     write(*,*)'land2cplr%tr',land2cplr%tr(l,k,:)
     write(*,*)'#### update_land_bc_fast ### end of output ####'
  endif

  ! ---- diagnostic section
  call send_tile_data(id_vegn_cover, vegn_cover, tile%diag)
  call send_tile_data(id_cosz, cosz, tile%diag)
  call send_tile_data(id_albedo_dir, tile%land_refl_dir, tile%diag)
  call send_tile_data(id_albedo_dif, tile%land_refl_dif, tile%diag)
  call send_tile_data(id_vegn_refl_dir, vegn_refl_dir,     tile%diag)
  call send_tile_data(id_vegn_refl_dif, vegn_refl_dif, tile%diag)
  call send_tile_data(id_vegn_refl_lw,  vegn_refl_lw, tile%diag)
  call send_tile_data(id_vegn_tran_dir, vegn_tran_dir_dir, tile%diag)
  call send_tile_data(id_vegn_tran_dif, vegn_tran_dif, tile%diag)
  call send_tile_data(id_vegn_tran_lw,  vegn_tran_lw, tile%diag)
  call send_tile_data(id_vegn_sctr_dir, vegn_tran_dir,     tile%diag)
  call send_tile_data(id_subs_refl_dir, subs_refl_dir, tile%diag)
  call send_tile_data(id_subs_refl_dif, subs_refl_dif, tile%diag)
  call send_tile_data(id_grnd_T,  land_grnd_T(tile),   tile%diag)
  ! CMOR variables
  call send_tile_data(id_snd, max(snow_depth,0.0),     tile%diag)
  call send_tile_data(id_snc, snow_area*100,           tile%diag)
  ! --- debug section
  call check_temp_range(land2cplr%t_ca(l,k),'update_land_bc_fast','T_ca')

end subroutine update_land_bc_fast


! ============================================================================
subroutine update_land_bc_slow (land2cplr)
  type(land_data_type), intent(inout) :: land2cplr

  ! ---- local vars
  integer :: i,j,k,face,l ! coordinates of the watch point, for debug printout

  call update_topo_rough(land2cplr%rough_scale)
  where (land2cplr%mask) &
       land2cplr%rough_scale = max(land2cplr%rough_mom,land2cplr%rough_scale)
  if (.not.use_coast_topo_rough) then
     do k = 1, size(land2cplr%rough_mom,2)
        where(lnd%landfrac<max_coast_frac) &
            land2cplr%rough_scale(:,k) = land2cplr%rough_mom(:,k)
     enddo
  endif
  call get_watch_point(i,j,k,face,l)
  if ( lnd%face==face.and.                   &
       lnd%ls<=l.and.l<=lnd%le.and.          &
       k<=size(land2cplr%rough_scale,2).and. &
       is_watch_time()) then
     write(*,*)'#### update_land_bc_slow ### output ####'
     write(*,*)'land2cplr%rough_scale',land2cplr%rough_scale(l,k)
     write(*,*)'#### update_land_bc_slow ### end of output ####'
  endif

end subroutine update_land_bc_slow

! ============================================================================
subroutine land_sfc_water(tile, grnd_liq, grnd_ice, grnd_subl, grnd_tf)
  type(land_tile_type), intent(in) :: tile
  real, intent(out) :: &
     grnd_liq, grnd_ice, & ! surface liquid and ice, respectively, kg/m2
     grnd_subl, &          ! fraction of vapor flux that sublimates
     grnd_tf               ! freezing temperature

  if (associated(tile%glac)) call glac_sfc_water(tile%glac, grnd_liq, grnd_ice, grnd_subl, grnd_tf)
  if (associated(tile%lake)) call lake_sfc_water(tile%lake, grnd_liq, grnd_ice, grnd_subl, grnd_tf)
  if (associated(tile%soil)) call soil_sfc_water(tile%soil, grnd_liq, grnd_ice, grnd_subl, grnd_tf)

  if (snow_active(tile%snow)) call snow_sfc_water(tile%snow, grnd_liq, grnd_ice, grnd_subl, grnd_tf)
end subroutine land_sfc_water

! ============================================================================
subroutine grnd_evap_limits(tile, grnd_E_min, grnd_E_max)
  type(land_tile_type), intent(in) :: tile
  real, intent(out) :: grnd_E_min, grnd_E_max

  grnd_E_min = -HUGE(grnd_E_min);
  grnd_E_max =  HUGE(grnd_E_max);
  if (associated(tile%soil).and..not.snow_active(tile%snow)) &
      call soil_evap_limits(tile%soil, grnd_E_min, grnd_E_max)
end subroutine grnd_evap_limits

! ============================================================================
real function land_grnd_T(tile)
  type(land_tile_type), intent(in) :: tile

  if (associated(tile%glac)) land_grnd_T = tile%glac%T(1)
  if (associated(tile%lake)) land_grnd_T = tile%lake%T(1)
  if (associated(tile%soil)) land_grnd_T = tile%soil%T(1)

  if (snow_active(tile%snow)) land_grnd_T = tile%snow%T(1)
end function land_grnd_T

! ============================================================================
subroutine Lnd_stock_pe(bnd,index,value)

type(land_data_type), intent(in)  :: bnd
integer             , intent(in)  :: index
real                , intent(out) :: value ! Domain water (Kg) or heat (Joules)

integer :: i,j,n
type(land_tile_enum_type)     :: ce
type(land_tile_type), pointer :: tile
character(len=128) :: message
integer :: ls,le,l
real :: area_factor, river_value
! *_twd are tile water densities (kg water per m2 of tile)
real twd_gas_cana,               twd_liq_glac, twd_sol_glac, &
     twd_liq_lake, twd_sol_lake, twd_liq_soil, twd_sol_soil, &
     twd_liq_snow, twd_sol_snow, twd_liq_vegn, twd_sol_vegn
! *_gcwd are grid-cell water densities (kg water per m2 of land in grid cell)
real gcwd_cana, gcwd_glac, gcwd_lake, gcwd_soil, gcwd_snow, gcwd_vegn
! v_* are global masses of water
real v_cana, v_glac, v_lake, v_soil, v_snow, v_vegn
real cana_q, a_globe

value = 0.0
v_cana = 0.
v_glac = 0.
v_lake = 0.
v_soil = 0.
v_snow = 0.
v_vegn = 0.
if(.not.bnd%pe) return
ls = lnd%ls
le = lnd%le
! The following is a dirty getaround
if(lnd%cellarea(ls) < 1.0) then
  area_factor = 4*pi*radius**2 ! lnd_sg%area is fraction of globe
else
  area_factor = 1.0 ! lnd_sg%area is actual area (m**2)
endif

select case(index)
case(ISTOCK_WATER)
  do l = ls, le
    ce = first_elmt(land_tile_map(l))
    gcwd_cana = 0.0; gcwd_glac = 0.0; gcwd_lake = 0.0
    gcwd_soil = 0.0; gcwd_snow = 0.0; gcwd_vegn = 0.0
    do while(loop_over_tiles(ce,tile))
      twd_gas_cana = 0.0
      twd_liq_glac = 0.0 ; twd_sol_glac = 0.0
      twd_liq_lake = 0.0 ; twd_sol_lake = 0.0
      twd_liq_soil = 0.0 ; twd_sol_soil = 0.0
      twd_liq_snow = 0.0 ; twd_sol_snow = 0.0
      twd_liq_vegn = 0.0 ; twd_sol_vegn = 0.0
      if(associated(tile%cana)) then
        call cana_state ( tile%cana, cana_q=cana_q )
        twd_gas_cana = canopy_air_mass*cana_q
        endif
      if(associated(tile%glac)) &
        call glac_tile_stock_pe(tile%glac, twd_liq_glac, twd_sol_glac)
      if(associated(tile%lake)) &
        call lake_tile_stock_pe(tile%lake, twd_liq_lake, twd_sol_lake)
      if(associated(tile%soil)) &
        call soil_tile_stock_pe(tile%soil, twd_liq_soil, twd_sol_soil)
      if(associated(tile%snow)) &
        call snow_tile_stock_pe(tile%snow, twd_liq_snow, twd_sol_snow)
      if(associated(tile%vegn)) &
        call vegn_tile_stock_pe(tile%vegn, twd_liq_vegn, twd_sol_vegn)
      gcwd_cana = gcwd_cana +  twd_gas_cana                 * tile%frac
      gcwd_glac = gcwd_glac + (twd_liq_glac + twd_sol_glac) * tile%frac
      gcwd_lake = gcwd_lake + (twd_liq_lake + twd_sol_lake) * tile%frac
      gcwd_soil = gcwd_soil + (twd_liq_soil + twd_sol_soil) * tile%frac
      gcwd_snow = gcwd_snow + (twd_liq_snow + twd_sol_snow) * tile%frac
      gcwd_vegn = gcwd_vegn + (twd_liq_vegn + twd_sol_vegn) * tile%frac
    enddo
    v_cana = v_cana + gcwd_cana * lnd%area(l)*area_factor
    v_glac = v_glac + gcwd_glac * lnd%area(l)*area_factor
    v_lake = v_lake + gcwd_lake * lnd%area(l)*area_factor
    v_soil = v_soil + gcwd_soil * lnd%area(l)*area_factor
    v_snow = v_snow + gcwd_snow * lnd%area(l)*area_factor
    v_vegn = v_vegn + gcwd_vegn * lnd%area(l)*area_factor
  enddo
  value  = v_cana + v_glac + v_lake + v_soil + v_snow + v_vegn
a_globe = 4. * pi * radius**2
case(ISTOCK_HEAT)
  if(.not.stock_warning_issued) then
    call error_mesg('Lnd_stock_pe','Heat stock not yet implemented',NOTE)
    stock_warning_issued = .true.
  endif
! do j = js, je
! do i = is, ie
!   ce = first_elmt(land_tile_map(i,j))
!   te = tail_elmt (land_tile_map(i,j))
!   grid_cell_heat_density = 0.0
!   do while(ce /= te)
!     tile => current_tile(ce)
!     tile_heat_density = 0.0
!     if(associated(tile%soil)) then
!       do n=1, size(tile%soil%temp)
!       tile_heat_density = tile_heat_density + (tile%soil%T(n)-tfreeze)* &
!                  (tile%soil%heat_capacity_dry(n)*dz(n) + &
!                   clw*tile%soil%wl(n)             + &
!                   csw*tile%soil%ws(n))
!       enddo
!       tile_heat_density = tile_heat_density + clw*soil%groundwater(1)*(soil%groundwater_T(1)-tfreeze) ! Why is this outside n loop?
!     endif
!     grid_cell_heat_density = grid_cell_heat_density + tile_heat_density * tile%frac
!     ce=next_elmt(ce)
!   enddo
!   grid_cell_heat = grid_cell_heat_density * lnd_sg%area(i,j)*area_factor
!   value = value + grid_cell_heat
! enddo
! enddo
case(ISTOCK_SALT)
  if(.not.stock_warning_issued) then
    call error_mesg('Lnd_stock_pe','Salt stock not yet implemented',NOTE)
    stock_warning_issued = .true.
  endif
case default
  write(message,'(i2,a,i2,a,i2,a,i2,a)') &
  index,' is an invalid stock index. Must be ISTOCK_WATER or ISTOCK_HEAT or ISTOCK_SALT (',ISTOCK_WATER,' or ',ISTOCK_HEAT,' or ', ISTOCK_SALT,')'
  call error_mesg('Lnd_stock_pe',message,FATAL)
end select

call river_stock_pe(index, river_value)
value = value + river_value

if (index.eq.ISTOCK_WATER.and.give_stock_details) then
    call mpp_sum(river_value, pelist=lnd%pelist)
    call mpp_sum(v_cana, pelist=lnd%pelist)
    call mpp_sum(v_glac, pelist=lnd%pelist)
    call mpp_sum(v_lake, pelist=lnd%pelist)
    call mpp_sum(v_soil, pelist=lnd%pelist)
    call mpp_sum(v_snow, pelist=lnd%pelist)
    call mpp_sum(v_vegn, pelist=lnd%pelist)
    write (message,'(a,f10.5)') 'total land storage:',v_cana/a_globe+v_glac/a_globe+ &
        v_lake/a_globe+v_soil/a_globe+v_snow/a_globe+v_vegn/a_globe+river_value/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '...canopy air:',v_cana/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '......glacier:',v_glac/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '.........lake:',v_lake/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '.........soil:',v_soil/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '.........snow:',v_snow/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '.........vegn:',v_vegn/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
    write (message,'(a,7f10.5)') '.......rivers:',river_value/a_globe
    call error_mesg('Lnd_stock_pe',message,NOTE)
endif

end subroutine Lnd_stock_pe


! ============================================================================
! initialize horizontal axes for land grid so that all sub-modules can use them,
! instead of creating their own
subroutine land_diag_init(clonb, clatb, clon, clat, time, domain, id_band, id_ug)
 !Inputs/outputs
  real,dimension(:),intent(in) :: clonb   !<longitudes of grid cells vertices
  real,dimension(:),intent(in) :: clatb   !<latitudes of grid cells vertices
  real,dimension(:),intent(in) :: clon    !<Longitude of grid cell centers.
  real,dimension(:),intent(in) :: clat    !<Latitude of grid cell centers
  type(time_type),intent(in)   :: time    !<Initial time for diagnostic fields.
  type(domainUG), intent(in)   :: domain  !<
  integer,intent(out)          :: id_band !<"band" axis id.
  integer,intent(out)          :: id_ug   !<Unstructured axis id.

  ! ---- local vars ----------------------------------------------------------
  integer :: nlon, nlat       ! sizes of respective axes
  integer             :: axes(1)        ! Array of axes for 1-D unstructured fields.
  integer             :: ug_dim_size    ! Size of the unstructured axis
  integer,allocatable :: ug_dim_data(:) ! Unstructured axis data.
  integer             :: id_lon, id_lonb
  integer             :: id_lat, id_latb
  integer :: i
  character(32) :: name       ! tracer name

  ! Register the unstructured axis for the unstructured domain.
  call mpp_get_UG_compute_domain(domain, size=ug_dim_size)
  if (.not. allocated(ug_dim_data)) then
      allocate(ug_dim_data(ug_dim_size))
  endif
  call mpp_get_UG_domain_grid_index(domain, ug_dim_data)
  !--- grid_index needs to be starting from 0.
  ug_dim_data = ug_dim_data - 1
  id_ug = diag_axis_init("grid_index",  real(ug_dim_data), "none", "U", long_name="grid indices", &
                         set_name="land", DomainU=domain, aux="geolon_t geolat_t")
  if (allocated(ug_dim_data)) then
      deallocate(ug_dim_data)
  endif

 ! Register horizontal axes that are required by the post-processing so that the output 
 ! files can be "decompressed": converted from unstructured back to lon-lat or cubic sphere.
 ! The "grid_xt" and "grid_yt" axes should run from 1 to the total number of x- and 
 ! y-points on cubic sphere face. It is assumed that all faces tiles contain the same 
 ! number of x- and y-points.
  nlon = size(clon)
  nlat = size(clat)
  if(mpp_get_UG_domain_ntiles(lnd%domain)==1) then
     ! grid has just one tile, so we assume that the grid is regular lat-lon
     ! define geographic axes and its edges
     id_lonb = diag_axis_init ('lonb', clonb, 'degrees_E', 'X', 'longitude edges', set_name='land')
     id_lon  = diag_axis_init ('lon',  clon,  'degrees_E', 'X', 'longitude', set_name='land',  edges=id_lonb)
     id_latb = diag_axis_init ('latb', clatb, 'degrees_N', 'Y', 'latitude edges', set_name='land')
     id_lat  = diag_axis_init ('lat',  clat,  'degrees_N', 'Y', 'latitude', set_name='land', edges=id_latb)
     ! add "compress" attribute to the unstructured grid axis
     call diag_axis_add_attribute(id_ug, "compress", "lat lon")
  else
     id_lon = diag_axis_init ( 'grid_xt', (/(real(i),i=1,nlon)/), 'degrees_E', 'X', &
          'T-cell longitude', set_name='land' )
     id_lat = diag_axis_init ( 'grid_yt', (/(real(i),i=1,nlat)/), 'degrees_N', 'Y', &
          'T-cell latitude', set_name='land' )
     ! add "compress" attribute to the unstructured grid axis
     call diag_axis_add_attribute(id_ug, "compress", "grid_yt grid_xt")
  endif

  id_band = diag_axis_init ('band',  (/1.0,2.0/), 'unitless', 'Z', 'spectral band', set_name='land' )

  ! Set up an array of axes ids, for convenience.
  axes(1) = id_ug

  ! register auxiliary coordinate variables

  id_geolon_t = register_static_field ( module_name, 'geolon_t', axes, &
       'longitude of grid cell centers', 'degrees_E', missing_value = -1.0e+20 )
  id_geolat_t = register_static_field ( module_name, 'geolat_t', axes, &
       'latitude of grid cell centers', 'degrees_N', missing_value = -1.0e+20 )

  ! register static diagnostic fields

  id_cellarea = register_static_field ( module_name, 'cell_area', axes, &
       'total area in grid cell', 'm2', missing_value=-1.0 )
  id_landfrac = register_static_field ( module_name, 'land_frac', axes, &
       'fraction of land in grid cell','unitless', missing_value=-1.0, area=id_cellarea)
  call diag_field_add_attribute(id_landfrac,'ocean_fillvalue',0)

  ! register areas and fractions for the rest of the diagnostic fields
  call register_tiled_area_fields(module_name, axes, time, id_area, id_frac)

  ! set the default filter (for area and subsampling) for consequent calls to
  ! register_tiled_diag_field
  call set_default_diag_filter('land')

  ! register regular (dynamic) diagnostic fields

  id_ntiles = register_tiled_diag_field(module_name,'ntiles', axes,  &
       time, 'number of tiles', 'unitless', missing_value=-1.0, op=OP_SUM)


  id_VWS = register_tiled_diag_field ( module_name, 'VWS', axes, time, &
             'vapor storage on land', 'kg/m2', missing_value=-1.0e+20 )
  id_VWSc    = register_tiled_diag_field ( module_name, 'VWSc', axes, time, &
             'vapor mass in canopy air', 'kg/m2', missing_value=-1.0e+20 )
  id_LWS = register_tiled_diag_field ( module_name, 'LWS', axes, time, &
             'liquid storage on land', 'kg/m2', missing_value=-1.0e+20 )
  id_LWSv    = register_tiled_diag_field ( module_name, 'LWSv', axes, time, &
             'liquid interception storage', 'kg/m2', missing_value=-1.0e+20 )
  id_LWSs    = register_tiled_diag_field ( module_name, 'LWSs', axes, time, &
             'liquid storage in snowpack', 'kg/m2', missing_value=-1.0e+20 )
  id_LWSg    = register_tiled_diag_field ( module_name, 'LWSg', axes, time, &
             'liquid ground storage', 'kg/m2', missing_value=-1.0e+20 )
  id_FWS = register_tiled_diag_field ( module_name, 'FWS', axes, time, &
             'frozen storage on land', 'kg/m2', missing_value=-1.0e+20 )
  id_FWSv    = register_tiled_diag_field ( module_name, 'FWSv', axes, time, &
             'frozen interception storage', 'kg/m2', missing_value=-1.0e+20 )
  id_FWSs    = register_tiled_diag_field ( module_name, 'FWSs', axes, time, &
             'frozen storage in snowpack', 'kg/m2', missing_value=-1.0e+20 )
  id_FWSg    = register_tiled_diag_field ( module_name, 'FWSg', axes, time, &
             'frozen ground storage', 'kg/m2', missing_value=-1.0e+20 )
  id_HS = register_tiled_diag_field ( module_name, 'HS', axes, time, &
             'land heat storage', 'J/m2', missing_value=-1.0e+20 )
  id_HSv     = register_tiled_diag_field ( module_name, 'HSv', axes, time, &
             'interception heat storage', 'J/m2', missing_value=-1.0e+20 )
  id_HSs     = register_tiled_diag_field ( module_name, 'HSs', axes, time, &
             'heat storage in snowpack', 'J/m2', missing_value=-1.0e+20 )
  id_HSg     = register_tiled_diag_field ( module_name, 'HSg', axes, time, &
             'ground heat storage', 'J/m2', missing_value=-1.0e+20 )
  id_HSc     = register_tiled_diag_field ( module_name, 'HSc', axes, time, &
             'canopy-air heat storage', 'J/m2', missing_value=-1.0e+20 )

  id_precip = register_tiled_diag_field ( module_name, 'precip', axes, time, &
             'precipitation rate', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_hprec = register_tiled_diag_field ( module_name, 'hprec', axes, time, &
             'sensible heat of precipitation', 'W/m2', missing_value=-1.0e+20 )
  id_lprec = register_tiled_diag_field ( module_name, 'lprec_l', axes, time, &
             'rainfall to land', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_lprecv = register_tiled_diag_field ( module_name, 'lprecv', axes, time, &
             'net rainfall to vegetation', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_lprecs = register_tiled_diag_field ( module_name, 'lprecs', axes, time, &
             'rainfall to snow, minus drainage', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_lprecg = register_tiled_diag_field ( module_name, 'lprecg', axes, time, &
             'effective rainfall to ground sfc', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_hlprec  = register_tiled_diag_field ( module_name, 'hlprec', axes, time, &
             'total liq precipitation heat', 'W/m2', missing_value=-1.0e+20)
  id_hlprecv = register_tiled_diag_field ( module_name, 'hlprecv', axes, time, &
             'net liq heat to vegetation', 'W/m2', missing_value=-1.0e+20)
  id_hlprecs = register_tiled_diag_field ( module_name, 'hlprecs', axes, time, &
             'net liq heat to snow', 'W/m2', missing_value=-1.0e+20)
  id_hlprecg = register_tiled_diag_field ( module_name, 'hlprecg', axes, time, &
             'net liq heat to ground sfc', 'W/m2', missing_value=-1.0e+20)
  id_fprec = register_tiled_diag_field ( module_name, 'fprec_l', axes, time, &
             'snowfall to land', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_fprecv = register_tiled_diag_field ( module_name, 'fprecv', axes, time, &
             'net snowfall to vegetation', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_fprecs = register_tiled_diag_field ( module_name, 'fprecs', axes, time, &
             'effective snowfall to snowpack', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_hfprec = register_tiled_diag_field ( module_name, 'hfprec', axes, time, &
             'sens heat of snowfall', 'W/m2', missing_value=-1.0e+20)
  id_hfprecv = register_tiled_diag_field ( module_name, 'hfprecv', axes, time, &
             'net sol heat to vegetation', 'W/m2', missing_value=-1.0e+20)
  id_hfprecs = register_tiled_diag_field ( module_name, 'hfprecs', axes, time, &
             'net sol heat to snow', 'W/m2', missing_value=-1.0e+20)
  id_evap = register_tiled_diag_field ( module_name, 'evap', axes, time, &
             'vapor flux up from land', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_hevap = register_tiled_diag_field ( module_name, 'hevap', axes, time, &
             'sensible heat of evap', 'W/m2', missing_value=-1.0e+20 )
  id_levap   = register_tiled_diag_field ( module_name, 'levap', axes, time, &
             'vapor flux from all liq (inc Tr)', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_levapv = register_tiled_diag_field ( module_name, 'levapv', axes, time, &
             'vapor flux leaving intercepted liquid', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_levaps = register_tiled_diag_field ( module_name, 'levaps', axes, time, &
             'vapor flux leaving snow liquid', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_levapg = register_tiled_diag_field ( module_name, 'levapg', axes, time, &
             'vapor flux from ground liquid', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_hlevap = register_tiled_diag_field ( module_name, 'hlevap', axes, time, &
             'vapor flux heat from liq source', 'W/m2', missing_value=-1.0e+20)
  id_hlevapv = register_tiled_diag_field ( module_name, 'hlevapv', axes, time, &
             'vapor heat from liq interc', 'W/m2', missing_value=-1.0e+20)
  id_hlevaps = register_tiled_diag_field ( module_name, 'hlevaps', axes, time, &
             'vapor heat from snow liq', 'W/m2', missing_value=-1.0e+20)
  id_hlevapg = register_tiled_diag_field ( module_name, 'hlevapg', axes, time, &
             'vapor heat from ground liq', 'W/m2', missing_value=-1.0e+20)
  id_fevap   = register_tiled_diag_field ( module_name, 'fevap', axes, time, &
             'vapor flux from all ice', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_fevapv = register_tiled_diag_field ( module_name, 'fevapv', axes, time, &
             'vapor flux leaving vegn ice', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_fevaps = register_tiled_diag_field ( module_name, 'fevaps', axes, time, &
             'vapor flux leaving snow ice', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_fevapg = register_tiled_diag_field ( module_name, 'fevapg', axes, time, &
             'vapor flux leaving ground ice', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_hfevap = register_tiled_diag_field ( module_name, 'hfevap', axes, time, &
             'vapor flux heat from solid source', 'W/m2', missing_value=-1.0e+20)
  id_hfevapv = register_tiled_diag_field ( module_name, 'hfevapv', axes, time, &
             'vapor heat from sol interc', 'W/m2', missing_value=-1.0e+20)
  id_hfevaps = register_tiled_diag_field ( module_name, 'hfevaps', axes, time, &
             'vapor heat from snow sol', 'W/m2', missing_value=-1.0e+20)
  id_hfevapg = register_tiled_diag_field ( module_name, 'hfevapg', axes, time, &
             'vapor heat from ground sol', 'W/m2', missing_value=-1.0e+20)
  id_runf   = register_tiled_diag_field ( module_name, 'runf', axes, time, &
             'total runoff', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_hrunf   = register_tiled_diag_field ( module_name, 'hrunf', axes, time, &
             'sensible heat of total runoff', 'W/m2', missing_value=-1.0e+20 )
  id_lrunf   = register_tiled_diag_field ( module_name, 'lrunf', axes, time, &
             'total rate of liq runoff', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_lrunfs  = register_tiled_diag_field ( module_name, 'lrunfs', axes, time, &
             'rate of liq runoff via calving', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_lrunfg  = register_tiled_diag_field ( module_name, 'lrunfg', axes, time, &
             'rate of liq runoff, ground', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_hlrunf  = register_tiled_diag_field ( module_name, 'hlrunf', axes, time, &
             'heat of liq runoff', 'W/m2', missing_value=-1.0e+20 )
  id_hlrunfs  = register_tiled_diag_field ( module_name, 'hlrunfs', axes, time, &
             'heat of liq runoff from snow pack', 'W/m2', missing_value=-1.0e+20 )
  id_hlrunfg  = register_tiled_diag_field ( module_name, 'hlrunfg', axes, time, &
             'heat of liq surface runoff', 'W/m2', missing_value=-1.0e+20 )
  id_frunf   = register_tiled_diag_field ( module_name, 'frunf', axes, time, &
             'total rate of solid runoff', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_frunfs  = register_tiled_diag_field ( module_name, 'frunfs', axes, time, &
             'rate of solid calving', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_hfrunf  = register_tiled_diag_field ( module_name, 'hfrunf', axes, time, &
             'heat of total ice runoff', 'W/m2', missing_value=-1.0e+20 )
  id_hfrunfs  = register_tiled_diag_field ( module_name, 'hfrunfs', axes, time, &
             'heat of sol snow runoff', 'W/m2', missing_value=-1.0e+20 )
  id_DOCrunf = register_tiled_diag_field ( module_name, 'lrunf_DOC', axes, time, &
             'total rate of DOC runoff', 'kgC/(m2 s)', missing_value=-1.0e+20 )
  id_melt    = register_tiled_diag_field ( module_name, 'melt', axes, time, &
             'total rate of melt', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_meltv   = register_tiled_diag_field ( module_name, 'meltv', axes, time, &
             'rate of melt, interception', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_melts   = register_tiled_diag_field ( module_name, 'melts', axes, time, &
             'rate of snow melt', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_meltg   = register_tiled_diag_field ( module_name, 'meltg', axes, time, &
             'rate of substrate thaw', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_fsw     = register_tiled_diag_field ( module_name, 'fsw', axes, time, &
             'net sw rad to land', 'W/m2', missing_value=-1.0e+20 )
  id_fswv    = register_tiled_diag_field ( module_name, 'fswv', axes, time, &
             'net sw rad to vegetation', 'W/m2', missing_value=-1.0e+20 )
  id_fsws    = register_tiled_diag_field ( module_name, 'fsws', axes, time, &
             'net sw rad to snow', 'W/m2', missing_value=-1.0e+20 )
  id_fswg    = register_tiled_diag_field ( module_name, 'fswg', axes, time, &
             'net sw rad to ground', 'W/m2', missing_value=-1.0e+20 )
  id_flw     = register_tiled_diag_field ( module_name, 'flw', axes, time, &
             'net lw rad to land', 'W/m2', missing_value=-1.0e+20 )
  id_flwv    = register_tiled_diag_field ( module_name, 'flwv', axes, time, &
             'net lw rad to vegetation', 'W/m2', missing_value=-1.0e+20 )
  id_flws    = register_tiled_diag_field ( module_name, 'flws', axes, time, &
             'net lw rad to snow', 'W/m2', missing_value=-1.0e+20 )
  id_flwg    = register_tiled_diag_field ( module_name, 'flwg', axes, time, &
             'net lw rad to ground', 'W/m2', missing_value=-1.0e+20 )
  id_cd_m    = register_tiled_diag_field ( module_name, 'cd_m', axes, time, &
       'drag coefficient for momentum', missing_value=-1e20)
  id_cd_t    = register_tiled_diag_field ( module_name, 'cd_t', axes, time, &
       'drag coefficient for heat and tracers', missing_value=-1e20)
  id_sens    = register_tiled_diag_field ( module_name, 'sens', axes, time, &
             'sens heat flux from land', 'W/m2', missing_value=-1.0e+20 )
  id_sensv   = register_tiled_diag_field ( module_name, 'sensv', axes, time, &
             'sens heat flux from vegn', 'W/m2', missing_value=-1.0e+20 )
  id_senss   = register_tiled_diag_field ( module_name, 'senss', axes, time, &
             'sens heat flux from snow', 'W/m2', missing_value=-1.0e+20 )
  id_sensg   = register_tiled_diag_field ( module_name, 'sensg', axes, time, &
             'sens heat flux from ground', 'W/m2', missing_value=-1.0e+20 )
  id_e_res_1 = register_tiled_diag_field ( module_name, 'e_res_1', axes, time, &
       'canopy air energy residual due to nonlinearities', 'W/m2', missing_value=-1e20)
  id_e_res_2 = register_tiled_diag_field ( module_name, 'e_res_2', axes, time, &
       'canopy energy residual due to nonlinearities', 'W/m2', missing_value=-1e20)
  id_z0m     = register_tiled_diag_field ( module_name, 'z0m', axes, time, &
             'momentum roughness of land', 'm', missing_value=-1.0e+20 )
  id_z0s     = register_tiled_diag_field ( module_name, 'z0s', axes, time, &
             'scalar roughness of land', 'm', missing_value=-1.0e+20 )
  id_con_g_h = register_tiled_diag_field ( module_name, 'con_g_h', axes, time, &
       'conductance for sensible heat between ground surface and canopy air', &
       'm/s', missing_value=-1.0 )
  id_transp  = register_tiled_diag_field ( module_name, 'transp', axes, time, &
             'Transpiration', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_wroff = register_tiled_diag_field ( module_name, 'wroff', axes, time, &
             'rate of liquid runoff to rivers', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_sroff = register_tiled_diag_field ( module_name, 'sroff', axes, time, &
             'rate of solid runoff to rivers', 'kg/(m2 s)', missing_value=-1.0e+20 )
  id_htransp = register_tiled_diag_field ( module_name, 'htransp', axes, time, &
             'heat of transpired vapior', 'W/m2', missing_value=-1.0e+20 )
  id_huptake = register_tiled_diag_field ( module_name, 'huptk', axes, time, &
             'heat of soil water uptake', 'W/m2', missing_value=-1.0e+20 )
  id_hroff = register_tiled_diag_field ( module_name, 'hroff', axes, time, &
             'sensible heat of runoff', 'W/m2', missing_value=-1.0e+20 )
  id_gsnow   = register_tiled_diag_field ( module_name, 'gsnow', axes, time, &
             'sens heat into ground from snow', 'W/m2', missing_value=-1.0e+20 )
  call add_tiled_diag_field_alias ( id_gsnow, module_name, 'gflux', axes, time, &
             'obsolete, please use "gsnow" instead', 'W/m2', missing_value=-1.0e+20 )
  id_gequil   = register_tiled_diag_field ( module_name, 'gequil', axes, time, &
             'snow-subs equilibration flux', 'W/m2', missing_value=-1.0e+20 )
  id_grnd_flux = register_tiled_diag_field ( module_name, 'grnd_flux', axes, time, &
             'sensible heat into ground from surface', 'W/m2', missing_value=-1.0e+20 )
  id_soil_water_supply = register_tiled_diag_field ( module_name, 'soil_water_supply', axes, time, &
       'maximum rate of soil water supply to vegetation', 'kg/(m2 s)', missing_value=-1e20)
  id_levapg_max = register_tiled_diag_field ( module_name, 'Eg_max', axes, time, &
             'soil_water limit on vapor flux from ground liquid', 'kg/(m2 s)', missing_value=-1.0e+20)
  id_water = register_tiled_diag_field ( module_name, 'water', axes, time, &
             'column-integrated soil water', 'kg/m2', missing_value=-1.0e+20 )
  id_snow = register_tiled_diag_field ( module_name, 'snow', axes, time, &
             'column-integrated snow water', 'kg/m2', missing_value=-1.0e+20 )
  id_Trad    = register_tiled_diag_field ( module_name, 'Trad', axes, time, &
             'radiative sfc temperature', 'degK', missing_value=-1.0e+20 )
  id_Tca     = register_tiled_diag_field ( module_name, 'Tca', axes, time, &
             'canopy-air temperature', 'degK', missing_value=-1.0e+20 )
  id_qca     = register_tiled_diag_field ( module_name, 'qca', axes, time, &
             'canopy-air specific humidity', 'kg/kg', missing_value=-1.0 )
  allocate(id_cana_tr(ntcana))
  do i = 1,ntcana
     call get_tracer_names(MODEL_LAND,i,name=name)
     id_cana_tr(i) = register_tiled_diag_field ( module_name, 'cana_'//trim(name), axes, time, &
             'canopy air '//trim(name)//' moist mass mixing ratio', 'kg/kg', missing_value=-1.0 )
  enddo
  id_qco2_dvmr = register_tiled_diag_field ( module_name, 'qco2_dvmr', axes, time, &
             'canopy-air CO2 dry volumetric mixing ratio', 'mol CO2/mol air', missing_value=-1.0 )
  id_fco2    = register_tiled_diag_field ( module_name, 'fco2', axes, time, &
             'flux of CO2 to canopy air', 'kg C/(m2 s)', missing_value=-1.0 )
   id_swdn_dir = register_tiled_diag_field ( module_name, 'swdn_dir', (/id_ug,id_band/), time, &
       'downward direct short-wave radiation flux to the land surface', 'W/m2', missing_value=-999.0)
   id_swdn_dif = register_tiled_diag_field ( module_name, 'swdn_dif', (/id_ug,id_band/), time, &
       'downward diffuse short-wave radiation flux to the land surface', 'W/m2', missing_value=-999.0)
   id_swup_dir = register_tiled_diag_field ( module_name, 'swup_dir', (/id_ug,id_band/), time, &
       'direct short-wave radiation flux reflected by the land surface', 'W/m2', missing_value=-999.0)
   id_swup_dif = register_tiled_diag_field ( module_name, 'swup_dif', (/id_ug,id_band/), time, &
       'diffuse short-wave radiation flux reflected by the land surface', 'W/m2', missing_value=-999.0)
  id_lwdn = register_tiled_diag_field ( module_name, 'lwdn', axes, time, &
       'downward long-wave radiation flux to the land surface', 'W/m2', missing_value=-999.0)
  id_vegn_cover = register_tiled_diag_field ( module_name, 'vegn_cover', axes, time, &
             'fraction covered by vegetation', missing_value=-1.0 )
  id_cosz = register_tiled_diag_field ( module_name, 'coszen', axes, time, &
       'cosine of zenith angle', missing_value=-2.0 )
  id_albedo_dir = register_tiled_diag_field ( module_name, 'albedo_dir', &
       (/id_ug,id_band/), time, &
       'land surface albedo for direct light', missing_value=-1.0 )
  id_albedo_dif = register_tiled_diag_field ( module_name, 'albedo_dif', &
       (/id_ug,id_band/), time, &
       'land surface albedo for diffuse light', missing_value=-1.0 )
  id_vegn_refl_dir = register_tiled_diag_field(module_name, 'vegn_refl_dir', &
       (/id_ug,id_band/), time, &
       'black-background canopy reflectivity for direct light',missing_value=-1.0)
  id_vegn_refl_dif = register_tiled_diag_field(module_name, 'vegn_refl_dif', &
       (/id_ug,id_band/), time, &
       'black-background canopy reflectivity for diffuse light',missing_value=-1.0)
  id_vegn_refl_lw = register_tiled_diag_field ( module_name, 'vegn_refl_lw', axes, time, &
       'canopy reflectivity for thermal radiation', missing_value=-1.0)
  id_vegn_tran_dir = register_tiled_diag_field(module_name, 'vegn_tran_dir', &
       (/id_ug,id_band/), time, &
       'part of direct light that passes through canopy unscattered',missing_value=-1.0)
  id_vegn_tran_dif = register_tiled_diag_field(module_name, 'vegn_tran_dif', &
       (/id_ug,id_band/), time, &
       'black-background canopy transmittance for diffuse light',missing_value=-1.0)
  id_vegn_tran_lw = register_tiled_diag_field ( module_name, 'vegn_tran_lw', axes, time, &
       'canopy transmittance for thermal radiation', missing_value=-1.0)
  id_vegn_sctr_dir = register_tiled_diag_field(module_name, 'vegn_sctr_dir', &
       (/id_ug,id_band/), time, &
       'part of direct light scattered downward by canopy',missing_value=-1.0)
  id_subs_refl_dir = register_tiled_diag_field(module_name, 'subs_refl_dir', &
       (/id_ug,id_band/), time, &
       'substrate reflectivity for direct light',missing_value=-1.0)
  id_subs_refl_dif = register_tiled_diag_field(module_name, 'subs_refl_dif', &
       (/id_ug,id_band/), time, &
       'substrate reflectivity for diffuse light',missing_value=-1.0)
  id_subs_emis = register_tiled_diag_field(module_name, 'subs_emis', axes, time, &
       'substrate emissivity for long-wave radiation',missing_value=-1.0)
  id_grnd_T = register_tiled_diag_field ( module_name, 'Tgrnd', axes, time, &
       'ground surface temperature', 'degK', missing_value=-1.0 )

  id_water_cons = register_tiled_diag_field ( module_name, 'water_cons', axes, time, &
       'water non-conservation in update_land_model_fast_0d', 'kg/(m2 s)', missing_value=-1.0 )
  id_carbon_cons = register_tiled_diag_field ( module_name, 'carbon_cons', axes, time, &
       'carbon non-conservation in update_land_model_fast_0d', 'kgC/(m2 s)', missing_value=-1.0 )

  ! CMOR variables
  id_pcp = register_tiled_diag_field ( cmor_name, 'pcp', axes, time, &
             'Total Precipitation', 'kg m-2 s-1', missing_value=-1.0e+20, &
             standard_name='total_precipitation_flux', fill_missing=.TRUE.)
  id_prra = register_tiled_diag_field ( cmor_name, 'prra', axes, time, &
             'Rainfall Rate', 'kg m-2 s-1', missing_value=-1.0e+20, &
             standard_name='rainfall_flux', fill_missing=.TRUE.)
  id_prveg = register_tiled_diag_field ( cmor_name, 'prveg', axes, time, &
             'Precipitation onto Canopy', 'kg m-2 s-1', missing_value=-1.0e+20, &
             standard_name='precipitation_flux_onto_canopy', fill_missing=.TRUE.)
  id_evspsblveg = register_tiled_diag_field ( cmor_name, 'evspsblveg', axes, time, &
             'Evaporation from Canopy', 'kg m-2 s-1', missing_value=-1.0e+20, &
             standard_name='water_evaporation_flux_from_canopy', fill_missing=.TRUE.)
  id_evspsblsoi = register_tiled_diag_field ( cmor_name, 'evspsblsoi', axes, time, &
             'Water Evaporation from Soil', 'kg m-2 s-1', missing_value=-1.0e+20, &
             standard_name='water_evaporation_flux_from_soil', fill_missing=.TRUE.)
  id_hflsLut = register_tiled_diag_field ( cmor_name, 'hflsLut', axes, time, &
             'Latent Heat Flux on Land Use Tile', 'W m-2', missing_value=-1.0e+20, &
             standard_name='surface_upward_latent_heat_flux', fill_missing=.FALSE.)
  id_nbp = register_tiled_diag_field ( cmor_name, 'nbp', axes, time, &
             'Carbon Mass Flux out of Atmosphere due to Net Biospheric Production on Land', &
             'kg C m-2 s-1', missing_value=-1.0, &
             standard_name='surface_net_downward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_all_land_processes', &
             fill_missing=.TRUE.)
  id_rlusLut = register_tiled_diag_field ( cmor_name, 'rlusLut', axes, time, &
             'Surface Upwelling Longwave on Land Use Tile', 'W m-2', &
             missing_value=-1.0e+20, standard_name='surface_upwelling_longwave_flux_in_air')
  id_rsusLut = register_tiled_diag_field ( cmor_name, 'rsusLut', axes, time, &
             'Surface Upwelling Shortwave on Land Use Tile', 'W m-2', &
             missing_value=-1.0e+20, standard_name='surface_upwelling_shortwave_flux_in_air')
  id_tslsiLut = register_tiled_diag_field ( cmor_name, 'tslsiLut', axes, time, &
             'Surface Skin Temperature on Land Use Tile', 'K', &
             missing_value=-1.0, standard_name='surface_temperature')
  id_tran  = register_tiled_diag_field ( cmor_name, 'tran', axes, time, &
             'Transpiration', 'kg m-2 s-1', missing_value=-1.0e+20, &
             standard_name='transpiration_flux', fill_missing=.TRUE.)
  id_snw = register_tiled_diag_field ( cmor_name, 'snw', axes, time, &
             'Surface Snow Amount','kg m-2', standard_name='surface_snow_amount', &
             missing_value=-1.0e+20, fill_missing=.TRUE.)
  id_snd = register_tiled_diag_field ( cmor_name, 'snd', axes, time, &
             'Snow Depth','m', standard_name='surface_snow_thickness', &
             missing_value=-1.0e+20, fill_missing=.TRUE.)
  id_snc = register_tiled_diag_field ( cmor_name, 'snc', axes, time, &
             'Snow Area Fraction','%', standard_name='surface_snow_area_fraction', &
             missing_value=-1.0e+20, fill_missing=.TRUE.)
  id_lwsnl = register_tiled_diag_field ( cmor_name, 'lwsnl', axes, time, &
             'Liquid Water Content of Snow Layer','kg m-2', &
             standard_name='liquid_water_content_of_snow_layer', &
             missing_value=-1.0e+20, fill_missing=.TRUE.)
  id_snm = register_tiled_diag_field ( cmor_name, 'snm', axes, time, &
             'Surface Snow Melt','kg m-2 s-1', standard_name='surface_snow_melt_flux', &
             missing_value=-1.0e+20, fill_missing=.TRUE.)
  call add_tiled_diag_field_alias(id_sens, cmor_name, 'hfssLut', axes, time, &
      'Sensible Heat Flux on Land Use Tile', 'W m-2', missing_value=-1.0e+20, &
      standard_name='surface_upward_sensible_heat_flux')

  id_sweLut = register_tiled_diag_field ( cmor_name, 'sweLut', axes, time, &
             'Snow Water Equivalent on Land Use Tile','m', standard_name='snow_water_equivalent_lut', &
             missing_value=-1.0e+20, fill_missing=.FALSE. )

  id_cLand = register_tiled_diag_field ( cmor_name, 'cLand', axes, time, &
             'Total carbon in all terrestrial carbon pools', 'kg C m-2', &
             standard_name='total_land_carbon', &
             missing_value=-1.0, fill_missing=.TRUE. )
  ! add alias for compatibility with older diag tables
  call add_tiled_diag_field_alias(id_cLand, module_name, 'Ctot', axes, time, &
     'total land carbon', 'kg C/m2', missing_value=-1.0)

  id_nwdFracLut = register_tiled_diag_field ( cmor_name, 'nwdFracLut', axes, time, &
             'Fraction of Land Use Tile Tile That is Non-Woody Vegetation', 'fraction', &
             standard_name='under_review', missing_value=-1.0)

  id_sftlf = register_static_field ( cmor_name, 'sftlf', axes, &
             'Land Area Fraction','%', standard_name='land_area_fraction', &
             area=id_cellarea)
  call diag_field_add_attribute(id_sftlf,'cell_methods','area: mean')
  id_sftgif = register_static_field ( cmor_name, 'sftgif', axes, &
             'Fraction of Grid Cell Covered with Glacier','%', standard_name='land_ice_area_fraction', &
             area=id_cellarea)
  call diag_field_add_attribute(id_sftgif,'cell_methods','area: mean')
  id_cropFrac = register_diag_field ( cmor_name, 'cropFrac', axes, time, &
             'Crop Fraction','%', standard_name='crop_fraction', &
             area=id_cellarea)
  call diag_field_add_attribute(id_cropFrac,'cell_methods','area: mean')
  id_cropFracC3 = register_diag_field ( cmor_name, 'cropFracC3', axes, time, &
             'C3 Crop Fraction','%', standard_name='crop_fraction_c3', &
             area=id_cellarea)
  call diag_field_add_attribute(id_cropFracC3,'cell_methods','area: mean')
  id_cropFracC4 = register_diag_field ( cmor_name, 'cropFracC4', axes, time, &
             'C4 Crop Fraction','%', standard_name='crop_fraction_c4', &
             area=id_cellarea)
  call diag_field_add_attribute(id_cropFracC4,'cell_methods','area: mean')
  id_pastureFrac = register_diag_field ( cmor_name, 'pastureFrac', axes, time, &
             'Anthropogenic Pasture Fraction','%', standard_name='anthropogenic_pasture_fraction', &
             area=id_cellarea)
  call diag_field_add_attribute(id_pastureFrac,'cell_methods','area: mean')
  id_residualFrac = register_static_field ( cmor_name, 'residualFrac', axes, &
             'Fraction of Grid Cell that is Land but Neither Vegetation-Covered nor Bare Soil','%', &
             standard_name='fraction_of_land_which_is_non_vegetation_and_non_bare_soil', &
             area=id_cellarea)
  call diag_field_add_attribute(id_cropFrac,'cell_methods','area: mean')
  id_treeFrac = register_diag_field ( cmor_name, 'treeFrac', axes, time, &
             'Tree Cover Fraction','%', standard_name='tree_cover_fraction', area=id_cellarea)
  call diag_field_add_attribute(id_treeFrac,'cell_methods','area: mean')
  id_grassFrac = register_diag_field ( cmor_name, 'grassFrac', axes, time, &
             'Natural Grass Fraction','%', standard_name='natural_grass_fraction', &
             area=id_cellarea)
  call diag_field_add_attribute(id_grassFrac,'cell_methods','area: mean')
  id_grassFracC3 = register_diag_field ( cmor_name, 'grassFracC3', axes, time, &
             'C3 Natural Grass Fraction','%', standard_name='natural_grass_fraction_c3', &
             area=id_cellarea)
  call diag_field_add_attribute(id_grassFracC3,'cell_methods','area: mean')
  id_grassFracC4 = register_diag_field ( cmor_name, 'grassFracC4', axes, time, &
             'C4 Natural Grass Fraction','%', standard_name='natural_grass_fraction_c4', &
             area=id_cellarea)
  call diag_field_add_attribute(id_grassFracC4,'cell_methods','area: mean')
  id_c3pftFrac = register_diag_field ( cmor_name, 'c3PftFrac', axes, time, &
             'Total C3 PFT Cover Fraction','%', standard_name='total_c3_pft_cover_fraction', &
             area=id_cellarea)
  call diag_field_add_attribute(id_c3pftFrac,'cell_methods','area: mean')
  id_c4pftFrac = register_diag_field ( cmor_name, 'c4PftFrac', axes, time, &
             'Total C4 PFT Cover Fraction','%', standard_name='total_c4_pft_cover_fraction', &
             area=id_cellarea)
  call diag_field_add_attribute(id_c4pftFrac,'cell_methods','area: mean')
  ! LUMIP land fractions
  id_fracLut_psl = register_diag_field ( cmor_name, 'fracLut_psl', axes, time, &
             'Fraction of Grid Cell for Each Land Use Tile','fraction', &
             standard_name='under_review', area=id_cellarea)
  id_fracLut_crp = register_diag_field ( cmor_name, 'fracLut_crop', axes, time, &
             'Fraction of Grid Cell for Each Land Use Tile','fraction', &
             standard_name='under_review', area=id_cellarea)
  id_fracLut_pst = register_diag_field ( cmor_name, 'fracLut_past', axes, time, &
             'Fraction of Grid Cell for Each Land Use Tile','fraction', &
             standard_name='under_review', area=id_cellarea)
  id_fracLut_urb = register_diag_field ( cmor_name, 'fracLut_urbn', axes, time, &
             'Fraction of Grid Cell for Each Land Use Tile','fraction', &
             standard_name='under_review', area=id_cellarea)
  call diag_field_add_attribute(id_fracLut_psl,'cell_methods','area: mean')
  call diag_field_add_attribute(id_fracLut_crp,'cell_methods','area: mean')
  call diag_field_add_attribute(id_fracLut_pst,'cell_methods','area: mean')
  call diag_field_add_attribute(id_fracLut_urb,'cell_methods','area: mean')
end subroutine land_diag_init

! ==============================================================================
subroutine send_cellfrac_data(id, f, scale)
  integer, intent(in) :: id ! id of the diagnostic field
  procedure(tile_exists_func) :: f ! existence detector function
  real, intent(in), optional  :: scale ! scaling factor, for unit conversions

  ! ---- local vars
  integer :: l,k
  logical :: used
  type(land_tile_type), pointer :: tile
  type(land_tile_enum_type) :: ce
  real :: frac(lnd%ls:lnd%le)
  real :: scale_

  if (.not.id>0) return ! do nothing if the field was not registered
  scale_ = 100.0 ! by fractions are in percent
  if (present(scale)) scale_ = scale

  frac(:) = 0.0
  ce = first_elmt(land_tile_map, ls=lnd%ls)
  do while (loop_over_tiles(ce, tile, l,k))
     ! accumulate fractions
     if (f(tile)) then
        frac(l) = frac(l)+tile%frac*scale_*lnd%landfrac(l)
     endif
  enddo
  used = send_data(id, frac, lnd%time)
end subroutine send_cellfrac_data

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_crop(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  if (.not.associated(tile%vegn)) return
  answer = (tile%vegn%landuse == LU_CROP)
end function is_crop

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_crop_C3(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  if (.not.associated(tile%vegn)) return
  answer = (tile%vegn%landuse == LU_CROP).and. &
           (tile%vegn%cohorts(1)%species /= SP_C4GRASS)
end function is_crop_C3

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_crop_C4(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  if (.not.associated(tile%vegn)) return
  answer = (tile%vegn%landuse == LU_CROP).and. &
           (tile%vegn%cohorts(1)%species == SP_C4GRASS)
end function is_crop_C4

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_pasture(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  if (.not.associated(tile%vegn)) return
  answer = (tile%vegn%landuse == LU_PAST)
end function is_pasture

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_urban(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  if (.not.associated(tile%vegn)) return
  answer = (tile%vegn%landuse == LU_URBN)
end function is_urban

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_psl(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  if (.not.associated(tile%vegn)) return
  answer = (tile%vegn%landuse == LU_NTRL.or.tile%vegn%landuse == LU_SCND)
end function is_psl

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_tree(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  if (.not.associated(tile%vegn)) return
  answer = (tile%vegn%cohorts(1)%species == SP_TEMPDEC).or. &
           (tile%vegn%cohorts(1)%species == SP_TROPICAL).or. &
           (tile%vegn%cohorts(1)%species == SP_EVERGR)
end function is_tree

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_ntrlgrass(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  if (.not.associated(tile%vegn)) return
  answer = ((tile%vegn%cohorts(1)%species == SP_C3GRASS).or. &
            (tile%vegn%cohorts(1)%species == SP_C4GRASS) &
           ).and.( &
            (tile%vegn%landuse==LU_NTRL).or. &
            (tile%vegn%landuse==LU_SCND) &
           )
end function is_ntrlgrass

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_ntrlgrass_C3(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  if (.not.associated(tile%vegn)) return
  answer = (tile%vegn%cohorts(1)%species == SP_C3GRASS).and.( &
            (tile%vegn%landuse==LU_NTRL).or. &
            (tile%vegn%landuse==LU_SCND) &
           )
end function is_ntrlgrass_C3

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_ntrlgrass_C4(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  if (.not.associated(tile%vegn)) return
  answer = (tile%vegn%cohorts(1)%species == SP_C4GRASS).and.( &
            (tile%vegn%landuse==LU_NTRL).or. &
            (tile%vegn%landuse==LU_SCND) &
           )
end function is_ntrlgrass_C4

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_C4(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  if (.not.associated(tile%vegn)) return
  answer = (tile%vegn%cohorts(1)%species == SP_C4GRASS)
end function is_C4

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_C3(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  if (.not.associated(tile%vegn)) return
  answer = (tile%vegn%cohorts(1)%species /= SP_C4GRASS)
end function is_C3

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_residual(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  answer = (associated(tile%glac).or.associated(tile%lake))
end function is_residual

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function is_glacier(tile) result(answer); logical :: answer
  type(land_tile_type), pointer :: tile

  answer = .FALSE.
  if (.not.associated(tile)) return
  answer = associated(tile%glac)
end function is_glacier

subroutine realloc_land2cplr ( bnd )
  type(land_data_type), intent(inout) :: bnd     ! data to allocate

  ! ---- local vars
  integer :: n_tiles

  call dealloc_land2cplr(bnd, dealloc_discharges=.FALSE.)

  bnd%domain = lnd_sg%domain
  bnd%ug_domain = lnd%domain
  n_tiles = max_n_tiles()


  ! allocate data according to the domain boundaries
  allocate( bnd%mask(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%tile_size(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%t_surf(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%t_ca(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%tr(lnd%ls:lnd%le,n_tiles,ntcana) )
  allocate( bnd%albedo(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%albedo_vis_dir(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%albedo_nir_dir(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%albedo_vis_dif(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%albedo_nir_dif(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%rough_mom(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%rough_heat(lnd%ls:lnd%le,n_tiles) )
  allocate( bnd%rough_scale(lnd%ls:lnd%le,n_tiles) )

  bnd%mask              = .FALSE.
  bnd%tile_size         = init_value
  bnd%t_surf            = init_value
  bnd%t_ca              = init_value
  bnd%tr                = init_value
  bnd%albedo            = init_value
  bnd%albedo_vis_dir    = init_value
  bnd%albedo_nir_dir    = init_value
  bnd%albedo_vis_dif    = init_value
  bnd%albedo_nir_dif    = init_value
  bnd%rough_mom         = init_value
  bnd%rough_heat        = init_value
  bnd%rough_scale       = init_value

  ! in contrast to the rest of the land boundary condition fields, discharges
  ! are specified per grid cell, not per tile; therefore they should not be
  ! re-allocated when the number of tiles changes. In fact, they must not be
  ! changed at all here because their values are assigned in update_land_model_fast,
  ! not in update_land_bc_*, and therefore would be lost if re-allocated.
  if (.not.associated(bnd%discharge)) then
     allocate( bnd%discharge          (lnd_sg%is:lnd_sg%ie, lnd_sg%js:lnd_sg%je) )
     allocate( bnd%discharge_heat     (lnd_sg%is:lnd_sg%ie, lnd_sg%js:lnd_sg%je) )
     allocate( bnd%discharge_snow     (lnd_sg%is:lnd_sg%ie, lnd_sg%js:lnd_sg%je) )
     allocate( bnd%discharge_snow_heat(lnd_sg%is:lnd_sg%ie, lnd_sg%js:lnd_sg%je) )

     ! discharge and dischargee_snow must be, in contrast to the rest of the boundary
     ! values, filled with zeroes. The reason is because not all of the usable elements
     ! are updated by the land model (only coastal points are).
     bnd%discharge           = 0.0
     bnd%discharge_heat      = 0.0
     bnd%discharge_snow      = 0.0
     bnd%discharge_snow_heat = 0.0
  endif

end subroutine realloc_land2cplr

! ============================================================================
! deallocates boundary data memory
! NOTE that the discharges should be deallocated only at the final clean-up
! stage; during the model run they should be preserved unchanged even when
! other fields are reallocated.
#define __DEALLOC__(x) if (associated(x)) deallocate(x)

subroutine dealloc_land2cplr ( bnd, dealloc_discharges )
  type(land_data_type), intent(inout) :: bnd  ! data to de-allocate
  logical, intent(in) :: dealloc_discharges

  __DEALLOC__( bnd%tile_size )
  __DEALLOC__( bnd%tile_size )
  __DEALLOC__( bnd%t_surf )
  __DEALLOC__( bnd%t_ca )
  __DEALLOC__( bnd%tr )
  __DEALLOC__( bnd%albedo )
  __DEALLOC__( bnd%albedo_vis_dir )
  __DEALLOC__( bnd%albedo_nir_dir )
  __DEALLOC__( bnd%albedo_vis_dif )
  __DEALLOC__( bnd%albedo_nir_dif )
  __DEALLOC__( bnd%rough_mom )
  __DEALLOC__( bnd%rough_heat )
  __DEALLOC__( bnd%rough_scale )
  __DEALLOC__( bnd%mask )

  if (dealloc_discharges) then
     __DEALLOC__( bnd%discharge           )
     __DEALLOC__( bnd%discharge_heat      )
     __DEALLOC__( bnd%discharge_snow      )
     __DEALLOC__( bnd%discharge_snow_heat )
  end if

end subroutine dealloc_land2cplr
! ============================================================================
! allocates boundary data for land domain and current number of tiles;
! initializes data for data override.
! NOTE: previously the body of the procedure was in the flux_exchange_init,
! currently it is called from land_model_init
subroutine realloc_cplr2land( bnd )
  type(atmos_land_boundary_type), intent(inout) :: bnd

  ! ---- local vars
  integer :: kd

  call dealloc_cplr2land(bnd)

  ! allocate data according to the domain boundaries
  kd = max_n_tiles()

  allocate( bnd%t_flux(lnd%ls:lnd%le,kd) )
  allocate( bnd%lw_flux(lnd%ls:lnd%le,kd) )
  allocate( bnd%sw_flux(lnd%ls:lnd%le,kd) )
  allocate( bnd%lprec(lnd%ls:lnd%le,kd) )
  allocate( bnd%fprec(lnd%ls:lnd%le,kd) )
  allocate( bnd%tprec(lnd%ls:lnd%le,kd) )
  allocate( bnd%dhdt(lnd%ls:lnd%le,kd) )
  allocate( bnd%dhdq(lnd%ls:lnd%le,kd) )
  allocate( bnd%drdt(lnd%ls:lnd%le,kd) )
  allocate( bnd%p_surf(lnd%ls:lnd%le,kd) )
  allocate( bnd%tr_flux(lnd%ls:lnd%le,kd,ntcana) )
  allocate( bnd%dfdtr(lnd%ls:lnd%le,kd,ntcana) )

  allocate( bnd%lwdn_flux(lnd%ls:lnd%le,kd) )
  allocate( bnd%swdn_flux(lnd%ls:lnd%le,kd) )
  allocate( bnd%sw_flux_down_vis_dir(lnd%ls:lnd%le,kd) )
  allocate( bnd%sw_flux_down_total_dir(lnd%ls:lnd%le,kd) )
  allocate( bnd%sw_flux_down_vis_dif(lnd%ls:lnd%le,kd) )
  allocate( bnd%sw_flux_down_total_dif(lnd%ls:lnd%le,kd) )
  allocate( bnd%cd_t(lnd%ls:lnd%le,kd) )
  allocate( bnd%cd_m(lnd%ls:lnd%le,kd) )
  allocate( bnd%bstar(lnd%ls:lnd%le,kd) )
  allocate( bnd%ustar(lnd%ls:lnd%le,kd) )
  allocate( bnd%wind(lnd%ls:lnd%le,kd) )
  allocate( bnd%z_bot(lnd%ls:lnd%le,kd) )

  allocate( bnd%drag_q(lnd%ls:lnd%le,kd) )

  bnd%t_flux                 = init_value
  bnd%lw_flux                = init_value
  bnd%sw_flux                = init_value
  bnd%lprec                  = init_value
  bnd%fprec                  = init_value
  bnd%tprec                  = init_value
  bnd%dhdt                   = init_value
  bnd%dhdq                   = init_value
  bnd%drdt                   = init_value
  bnd%p_surf                 = init_value
  bnd%tr_flux                = init_value
  bnd%dfdtr                  = init_value

  bnd%lwdn_flux              = init_value
  bnd%swdn_flux              = init_value
  bnd%sw_flux_down_vis_dir   = init_value
  bnd%sw_flux_down_total_dir = init_value
  bnd%sw_flux_down_vis_dif   = init_value
  bnd%sw_flux_down_total_dif = init_value
  bnd%cd_t                   = init_value
  bnd%cd_m                   = init_value
  bnd%bstar                  = init_value
  bnd%ustar                  = init_value
  bnd%wind                   = init_value
  bnd%z_bot                  = init_value

  bnd%drag_q                 = init_value

end subroutine realloc_cplr2land

! ============================================================================
subroutine dealloc_cplr2land( bnd )
  type(atmos_land_boundary_type), intent(inout) :: bnd

  __DEALLOC__( bnd%t_flux )
  __DEALLOC__( bnd%lw_flux )
  __DEALLOC__( bnd%sw_flux )
  __DEALLOC__( bnd%lprec )
  __DEALLOC__( bnd%fprec )
  __DEALLOC__( bnd%tprec )
  __DEALLOC__( bnd%dhdt )
  __DEALLOC__( bnd%dhdq )
  __DEALLOC__( bnd%drdt )
  __DEALLOC__( bnd%p_surf )
  __DEALLOC__( bnd%lwdn_flux )
  __DEALLOC__( bnd%swdn_flux )
  __DEALLOC__( bnd%sw_flux_down_vis_dir )
  __DEALLOC__( bnd%sw_flux_down_total_dir )
  __DEALLOC__( bnd%sw_flux_down_vis_dif )
  __DEALLOC__( bnd%sw_flux_down_total_dif )
  __DEALLOC__( bnd%cd_t )
  __DEALLOC__( bnd%cd_m )
  __DEALLOC__( bnd%bstar )
  __DEALLOC__( bnd%ustar )
  __DEALLOC__( bnd%wind )
  __DEALLOC__( bnd%z_bot )
  __DEALLOC__( bnd%tr_flux )
  __DEALLOC__( bnd%dfdtr )
  __DEALLOC__( bnd%drag_q )

end subroutine dealloc_cplr2land
#undef __DEALLOC__

! ===========================================================================
!  Prints checksums of the various fields in the atmos_land_boundary_type.
subroutine atm_lnd_bnd_type_chksum(id, timestep, albt)
    character(len=*), intent(in) :: id  ! Label to differentiate where this
                      ! routine is being called from.
    integer         , intent(in) :: timestep ! An integer to indicate which
                      ! timestep this routine is being called for.
    type(atmos_land_boundary_type), intent(in) :: albt
    integer ::   n, outunit

    outunit = stdout()

    write(outunit,*) 'BEGIN CHECKSUM(atmos_land_boundary_type):: ', id, timestep
    write(outunit,100) 'albt%t_flux                ', mpp_chksum( albt%t_flux)
    write(outunit,100) 'albt%lw_flux               ', mpp_chksum( albt%lw_flux)
    write(outunit,100) 'albt%lwdn_flux             ', mpp_chksum( albt%lwdn_flux)
    write(outunit,100) 'albt%sw_flux               ', mpp_chksum( albt%sw_flux)
    write(outunit,100) 'albt%swdn_flux               ', mpp_chksum( albt%swdn_flux)
    write(outunit,100) 'albt%lprec                 ', mpp_chksum( albt%lprec)
    write(outunit,100) 'albt%fprec                 ', mpp_chksum( albt%fprec)
    write(outunit,100) 'albt%tprec                 ', mpp_chksum( albt%tprec)
    write(outunit,100) 'albt%sw_flux_down_vis_dir  ', mpp_chksum( albt%sw_flux_down_vis_dir)
    write(outunit,100) 'albt%sw_flux_down_total_dir', mpp_chksum( albt%sw_flux_down_total_dir)
    write(outunit,100) 'albt%sw_flux_down_vis_dif  ', mpp_chksum( albt%sw_flux_down_vis_dif)
    write(outunit,100) 'albt%sw_flux_down_total_dif', mpp_chksum( albt%sw_flux_down_total_dif)
    write(outunit,100) 'albt%dhdt                  ', mpp_chksum( albt%dhdt)
    write(outunit,100) 'albt%dhdq                  ', mpp_chksum( albt%dhdq)
    write(outunit,100) 'albt%drdt                  ', mpp_chksum( albt%drdt)
    write(outunit,100) 'albt%cd_m                  ', mpp_chksum( albt%cd_m)
    write(outunit,100) 'albt%cd_t                  ', mpp_chksum( albt%cd_t)
    write(outunit,100) 'albt%ustar                 ', mpp_chksum( albt%ustar)
    write(outunit,100) 'albt%bstar                 ', mpp_chksum( albt%bstar)
    write(outunit,100) 'albt%wind                  ', mpp_chksum( albt%wind)
    write(outunit,100) 'albt%z_bot                 ', mpp_chksum( albt%z_bot)
    write(outunit,100) 'albt%drag_q                ', mpp_chksum( albt%drag_q)
    write(outunit,100) 'albt%p_surf                ', mpp_chksum( albt%p_surf)
    do n = 1,size(albt%tr_flux,3)
    write(outunit,100) 'albt%tr_flux               ', mpp_chksum( albt%tr_flux(:,:,n))
    enddo
    do n = 1,size(albt%dfdtr,3)
    write(outunit,100) 'albt%dfdtr                 ', mpp_chksum( albt%dfdtr(:,:,n))
    enddo

100 FORMAT("CHECKSUM::",A32," = ",Z20)

end subroutine atm_lnd_bnd_type_chksum


! ===========================================================================
!  Prints checksums of the various fields in the land_data_type.
subroutine land_data_type_chksum(id, timestep, land)
    character(len=*), intent(in) :: id ! Label to differentiate where this
        ! routine in being called from.
    integer         , intent(in) :: timestep ! An integer to indicate which
        ! timestep this routine is being called for.
    type(land_data_type), intent(in) :: land
    integer ::   n, outunit, k, j

    outunit = stdout()

    write(outunit,*) 'BEGIN CHECKSUM(land_data_type):: ', id, timestep

    write(outunit,100) 'land%tile_size         ',mpp_chksum(land%tile_size)
    write(outunit,100) 'land%t_surf            ',mpp_chksum(land%t_surf)
    write(outunit,100) 'land%t_ca              ',mpp_chksum(land%t_ca)
    write(outunit,100) 'land%albedo            ',mpp_chksum(land%albedo)
    write(outunit,100) 'land%albedo_vis_dir    ',mpp_chksum(land%albedo_vis_dir)
    write(outunit,100) 'land%albedo_nir_dir    ',mpp_chksum(land%albedo_nir_dir)
    write(outunit,100) 'land%albedo_vis_dif    ',mpp_chksum(land%albedo_vis_dif)
    write(outunit,100) 'land%albedo_nir_dif    ',mpp_chksum(land%albedo_nir_dif)
    write(outunit,100) 'land%rough_mom         ',mpp_chksum(land%rough_mom)
    write(outunit,100) 'land%rough_heat        ',mpp_chksum(land%rough_heat)
    write(outunit,100) 'land%rough_scale       ',mpp_chksum(land%rough_scale)

    do n = 1, size(land%tr,3)
    write(outunit,100) 'land%tr                ',mpp_chksum(land%tr(:,:,n))
    enddo
    write(outunit,100) 'land%discharge         ',mpp_chksum(land%discharge)
    write(outunit,100) 'land%discharge_snow    ',mpp_chksum(land%discharge_snow)
    write(outunit,100) 'land%discharge_heat    ',mpp_chksum(land%discharge_heat)


100 FORMAT("CHECKSUM::",A32," = ",Z20)
end subroutine land_data_type_chksum

! the code below defines the accessor routines that are used to access fields of the
! tile data structure in collective operations, like restart i/o. Fore example, a statement
! DEFINE_LAND_ACCESSOR_0D(real,lwup)
! defines a function "land_lwup_ptr" that, given a tile, returns a pointer to the field
! called "lwup" in this tile. The procedure implementing a collective operation would
! enumerate all the tiles within the domain, call accessor routine for each of them, and
! get or set the value pointed to by the accessor routine.
#define DEFINE_LAND_ACCESSOR_0D(xtype,x) subroutine land_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))p=>t%x;\
end subroutine

DEFINE_LAND_ACCESSOR_0D(real,frac)
DEFINE_LAND_ACCESSOR_0D(real,lwup)
DEFINE_LAND_ACCESSOR_0D(real,e_res_1)
DEFINE_LAND_ACCESSOR_0D(real,e_res_2)

! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function land_tile_exists(tile)
  type(land_tile_type), pointer :: tile
  land_tile_exists = associated(tile)
end function land_tile_exists

#define DEFINE_TAG_ACCESSOR(x) subroutine  x ## _tag_ptr(t,p);\
type(land_tile_type),pointer::t;integer,pointer::p;p=>NULL();if(associated(t))\
then;if (associated(t%x)) p=>t%x%tag;endif;\
end subroutine

DEFINE_TAG_ACCESSOR(glac)
DEFINE_TAG_ACCESSOR(lake)
DEFINE_TAG_ACCESSOR(soil)
DEFINE_TAG_ACCESSOR(vegn)

end module land_model_mod

