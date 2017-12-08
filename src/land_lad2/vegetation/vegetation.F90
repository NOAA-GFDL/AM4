module vegetation_mod

#include "../shared/debug.inc"

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only: error_mesg, NOTE,FATAL, file_exist, close_file, &
      check_nml_error, stdlog
use mpp_mod, only: mpp_sum, mpp_max, mpp_pe, mpp_root_pe, mpp_sync, stdout
use time_manager_mod, only: time_type, time_type_to_real, get_date, operator(-)
use constants_mod,    only: tfreeze, rdgas, rvgas, hlv, hlf, cp_air, PI
use sphum_mod, only: qscomp

use vegn_tile_mod, only: vegn_tile_type, &
     vegn_seed_demand, vegn_seed_supply, vegn_add_bliving, &
     cpw, clw, csw
use soil_tile_mod, only: soil_tile_type, soil_ave_temp, &
                         soil_ave_theta0, soil_ave_theta1, soil_psi_stress
use land_constants_mod, only : NBANDS, BAND_VIS, d608, mol_C, mol_CO2, mol_air, &
     seconds_per_year
use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_enum_type, &
     first_elmt, loop_over_tiles, land_tile_heat, land_tile_carbon, get_tile_water
use land_tile_diag_mod, only : register_tiled_static_field, register_tiled_diag_field, &
     add_tiled_diag_field_alias, set_default_diag_filter, send_tile_data, diag_buff_type, &
     OP_STD, OP_VAR, cmor_name
use land_data_mod, only : lnd, log_version
use land_io_mod, only : read_field

use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_restart_axis, add_tile_data, add_int_tile_data, add_scalar_data, &
     get_tile_data, get_int_tile_data, field_exists

use vegn_data_mod, only : SP_C4GRASS, LEAF_ON, LU_NTRL, read_vegn_data_namelist, &
     tau_drip_l, tau_drip_s, T_transp_min, cold_month_threshold, soil_carbon_depth_scale, &
     fsc_pool_spending_time, ssc_pool_spending_time, harvest_spending_time, &
     N_HARV_POOLS, HARV_POOL_NAMES, HARV_POOL_PAST, HARV_POOL_CROP, HARV_POOL_CLEARED, &
     HARV_POOL_WOOD_FAST, HARV_POOL_WOOD_MED, HARV_POOL_WOOD_SLOW, agf_bs
use vegn_cohort_mod, only : vegn_cohort_type, &
     vegn_data_heat_capacity, vegn_data_intrcptn_cap, update_species,&
     get_vegn_wet_frac, vegn_data_cover
use canopy_air_mod, only : cana_turbulence

use cohort_io_mod, only :  read_create_cohorts, create_cohort_dimension, &
     add_cohort_data, add_int_cohort_data, get_cohort_data, get_int_cohort_data
use land_debug_mod, only : is_watch_point, set_current_point, check_temp_range, &
     carbon_cons_tol, water_cons_tol, check_conservation, do_check_conservation
use vegn_radiation_mod, only : vegn_radiation_init, vegn_radiation
use vegn_photosynthesis_mod, only : vegn_photosynthesis_init, vegn_photosynthesis
use static_vegn_mod, only : read_static_vegn_namelist, static_vegn_init, static_vegn_end, &
     read_static_vegn
use vegn_dynamics_mod, only : vegn_dynamics_init, vegn_carbon_int, vegn_growth, &
     vegn_daily_npp, vegn_phenology, vegn_biogeography
use vegn_disturbance_mod, only : vegn_disturbance_init, vegn_nat_mortality, &
     vegn_disturbance, update_fuel
use vegn_harvesting_mod, only : &
     vegn_harvesting_init, vegn_harvesting_end, vegn_harvesting
use soil_carbon_mod, only : add_litter, poolTotalCarbon, cull_cohorts, &
     soil_carbon_option, SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE
use soil_mod, only : add_root_litter, redistribute_peat_carbon

use fms_io_mod, only: fms_io_unstructured_read

implicit none
private

! ==== public interfaces =====================================================
public :: read_vegn_namelist
public :: vegn_init
public :: vegn_end
public :: save_vegn_restart

public :: vegn_get_cover
public :: vegn_radiation
public :: vegn_properties

public :: vegn_step_1
public :: vegn_step_2
public :: vegn_step_3

public :: update_vegn_slow
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'vegn'
#include "../shared/version_variable.inc"

! values for internal selector of CO2 option used for photosynthesis
integer, parameter :: VEGN_PHOT_CO2_PRESCRIBED  = 1
integer, parameter :: VEGN_PHOT_CO2_INTERACTIVE = 2


! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
logical :: lm2               = .false.
real    :: init_Wl           = 0
real    :: init_Ws           = 0
real    :: init_Tv           = 288.
real    :: init_cohort_bl    = 0.05 ! initial biomass of leaves, kg C/m2
real    :: init_cohort_blv   = 0.0  ! initial biomass of labile store, kg C/m2
real    :: init_cohort_br    = 0.05 ! initial biomass of fine roots, kg C/m2
real    :: init_cohort_bsw   = 0.05 ! initial biomass of sapwood, kg C/m2
real    :: init_cohort_bwood = 0.05 ! initial biomass of heartwood, kg C/m2
real    :: init_cohort_cmc   = 0.0  ! initial intercepted water
character(32) :: rad_to_use = 'big-leaf' ! or 'two-stream'
character(32) :: snow_rad_to_use = 'ignore' ! or 'paint-leaves'
character(32) :: photosynthesis_to_use = 'simple' ! or 'leuning'
character(32) :: co2_to_use_for_photosynthesis = 'prescribed' ! or 'interactive'
   ! specifies what co2 concentration to use for photosynthesis calculations:
   ! 'prescribed'  : a prescribed value is used, equal to co2_for_photosynthesis
   !      specified below.
   ! 'interactive' : concentration of co2 in canopy air is used
real    :: co2_for_photosynthesis = 350.0e-6 ! concentration of co2 for photosynthesis
   ! calculations, mol/mol. Ignored if co2_to_use_for_photosynthesis is not 'prescribed'
logical :: do_cohort_dynamics   = .TRUE. ! if true, do vegetation growth
logical :: do_patch_disturbance = .TRUE. !
logical :: do_phenology         = .TRUE.
logical :: xwilt_available      = .TRUE.
logical :: do_biogeography      = .TRUE.
logical :: do_seed_transport    = .TRUE.
real    :: min_Wl=-1.0, min_Ws=-1.0 ! threshold values for condensation numerics, kg/m2:
   ! if water or snow on canopy fall below these values, the derivatives of
   ! condensation are set to zero, thereby prohibiting switching from condensation to
   ! evaporation in one time step.
real    :: tau_smooth_ncm = 0.0 ! Time scale for ncm smoothing
   ! (low-pass filtering), years. 0.0 retrieves previous behavior (no smoothing)
real :: rav_lit_0         = 0.0 ! constant litter resistance to vapor
real :: rav_lit_vi        = 0.0 ! litter resistance to vapor per LAI+SAI
real :: rav_lit_fsc       = 0.0 ! litter resistance to vapor per fsc
real :: rav_lit_ssc       = 0.0 ! litter resistance to vapor per ssc
real :: rav_lit_deadmic   = 0.0 ! litter resistance to vapor per dead microbe C
real :: rav_lit_bwood     = 0.0 ! litter resistance to vapor per bwood

logical :: do_peat_redistribution = .FALSE.

logical :: biodata_bug    = .FALSE. ! if true, initialization of t_ann, t_cold, p_ann,
                                    ! and p_cold from biodata is not done
namelist /vegn_nml/ &
    lm2, init_Wl, init_Ws, init_Tv, cpw, clw, csw, &
    init_cohort_bl, init_cohort_blv, init_cohort_br, init_cohort_bsw, &
    init_cohort_bwood, init_cohort_cmc, &
    rad_to_use, snow_rad_to_use, photosynthesis_to_use, &
    co2_to_use_for_photosynthesis, co2_for_photosynthesis, &
    do_cohort_dynamics, do_patch_disturbance, do_phenology, &
    xwilt_available, &
    do_biogeography, do_seed_transport, &
    min_Wl, min_Ws, tau_smooth_ncm, &
    rav_lit_0, rav_lit_vi, rav_lit_fsc, rav_lit_ssc, rav_lit_deadmic, rav_lit_bwood,&
    do_peat_redistribution, &
    biodata_bug

!---- end of namelist --------------------------------------------------------

logical         :: module_is_initialized =.FALSE.
type(time_type) :: time ! *** NOT YET USED
real            :: delta_time      ! fast time step
real            :: dt_fast_yr      ! fast time step in years
integer         :: vegn_phot_co2_option = -1 ! internal selector of co2 option
                                   ! used for photosynthesis
! diagnostic field ids
integer :: id_vegn_type, id_temp, id_wl, id_ws, id_height, &
   id_lai, id_sai, id_lai_var, id_lai_std, id_leaf_size, &
   id_root_density, id_root_zeta, id_rs_min, id_leaf_refl, id_leaf_tran,&
   id_leaf_emis, id_snow_crit, id_stomatal, id_an_op, id_an_cl, &
   id_bl, id_blv, id_br, id_bsw, id_bwood, id_btot, id_species, id_status, &
   id_con_v_h, id_con_v_v, id_fuel, id_harv_pool(N_HARV_POOLS), &
   id_harv_rate(N_HARV_POOLS), id_t_harv_pool, id_t_harv_rate, &
   id_csmoke_pool, id_csmoke_rate, id_fsc_in, id_fsc_out, id_ssc_in, &
   id_ssc_out, id_deadmic_in, id_deadmic_out, id_veg_in, id_veg_out, &
   id_fsc_pool_ag, id_fsc_rate_ag, id_fsc_pool_bg, id_fsc_rate_bg,&
   id_ssc_pool_ag, id_ssc_rate_ag, id_ssc_pool_bg, id_ssc_rate_bg,&
   id_leaflitter_buffer_ag, id_coarsewoodlitter_buffer_ag,id_leaflitter_buffer_rate_ag, id_coarsewoodlitter_buffer_rate_ag,& ! id_coarsewoodlitter_buffer_rate_ag is 34 characters long (pjp)
   id_t_ann, id_t_cold, id_p_ann, id_ncm, &
   id_lambda, id_afire, id_atfall, id_closs, id_cgain, id_wdgain, id_leaf_age, &
   id_phot_co2, id_theph, id_psiph, id_evap_demand
! CMOR variables
integer :: id_cProduct, id_cAnt, &
   id_fFire, id_fGrazing, id_fHarvest, id_fLuc, &
   id_cLeaf, id_cWood, id_cRoot, id_cMisc
! ==== end of module variables ===============================================

contains

! ============================================================================
subroutine read_vegn_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  logical :: use_static_veg ! if true, switch off vegetation dynamics

  call read_vegn_data_namelist()
  call read_static_vegn_namelist(use_static_veg)

  call log_version(version, module_name, &
  __FILE__)
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=vegn_nml, iostat=io)
    ierr = check_nml_error(io, 'vegn_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=vegn_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'vegn_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif

  unit=stdlog()

  ! switch off vegetation dynamics if static vegetation is set
  if (use_static_veg) then
     call error_mesg('vegn_init', &
          'use_static_veg=.TRUE., switching off vegetation dynamics', NOTE)
     write(unit,*)'use_static_veg=.TRUE., switching off vegetation dynamics'
     do_cohort_dynamics   = .FALSE.
     do_patch_disturbance = .FALSE.
     do_phenology         = .FALSE.
     do_biogeography      = .FALSE.
     do_seed_transport    = .FALSE.
  endif

  if (mpp_pe() == mpp_root_pe()) then
     write(unit, nml=vegn_nml)
  endif

  ! convert symbolic names of photosynthesis CO2 options into numeric IDs to
  ! speed up selection during run-time
  if (trim(co2_to_use_for_photosynthesis)=='prescribed') then
     vegn_phot_co2_option = VEGN_PHOT_CO2_PRESCRIBED
  else if (trim(co2_to_use_for_photosynthesis)=='interactive') then
     vegn_phot_co2_option = VEGN_PHOT_CO2_INTERACTIVE
  else
     call error_mesg('vegn_init',&
          'vegetation photosynthesis option co2_to_use_for_photosynthesis="'//&
          trim(co2_to_use_for_photosynthesis)//'" is invalid, use "prescribed" or "interactive"',&
          FATAL)
  endif

  ! ---- initialize vegetation radiation options
  call vegn_radiation_init(rad_to_use, snow_rad_to_use)

  ! ---- initialize vegetation photosynthesis options
  call vegn_photosynthesis_init(photosynthesis_to_use)

end subroutine read_vegn_namelist


! ============================================================================
! initialize vegetation
subroutine vegn_init(id_ug,id_band)
  integer,intent(in) :: id_ug   !<Unstructured axis id.
  integer,intent(in) :: id_band ! ID of spectral band axis

  ! ---- local vars
  integer :: unit         ! unit for various i/o
  type(land_tile_enum_type)     :: ce    ! tile list enumerator
  type(land_tile_type), pointer :: tile  ! pointer to current tile
  type(vegn_cohort_type), pointer :: cohort! pointer to initial cohort for cold-start
  integer :: n_accum
  integer :: nmn_acm
  type(land_restart_type) :: restart1, restart2
  logical :: restart_1_exists, restart_2_exists
  real, allocatable :: t_ann(:),t_cold(:),p_ann(:),ncm(:) ! buffers for biodata reading
  logical :: did_read_biodata
  integer :: i,j,l ! indices of current tile

  module_is_initialized = .TRUE.

  ! ---- make module copy of time and calculate time step ------------------
  delta_time = time_type_to_real(lnd%dt_fast)
  dt_fast_yr = delta_time/seconds_per_year


  ! ---- initialize vegn state ---------------------------------------------
  n_accum = 0
  nmn_acm = 0
  call open_land_restart(restart1,'INPUT/vegn1.res.nc',restart_1_exists)
  call open_land_restart(restart2,'INPUT/vegn2.res.nc',restart_2_exists)
  if (restart_1_exists) then
     call error_mesg('vegn_init',&
          'reading NetCDF restarts "INPUT/vegn1.res.nc" and "INPUT/vegn2.res.nc"',&
          NOTE)

     ! read the cohort index and generate appropriate number of cohorts
     ! for each vegetation tile
     call read_create_cohorts(restart1)
     call get_cohort_data(restart1, 'tv', cohort_tv_ptr)
     call get_cohort_data(restart1, 'wl', cohort_wl_ptr)
     call get_cohort_data(restart1, 'ws', cohort_ws_ptr)

     ! read global variables
     call fms_io_unstructured_read(restart2%basename, &
                                   "n_accum", &
                                   n_accum, &
                                   lnd%domain)
     call fms_io_unstructured_read(restart2%basename, &
                                   "nmn_acm", &
                                   nmn_acm, &
                                   lnd%domain)

     call get_int_cohort_data(restart2, 'species', cohort_species_ptr)
     call get_cohort_data(restart2, 'hite', cohort_height_ptr)
     call get_cohort_data(restart2, 'bl', cohort_bl_ptr)
     call get_cohort_data(restart2, 'blv', cohort_blv_ptr)
     call get_cohort_data(restart2, 'br', cohort_br_ptr)
     call get_cohort_data(restart2, 'bsw', cohort_bsw_ptr)
     call get_cohort_data(restart2, 'bwood', cohort_bwood_ptr)
     call get_cohort_data(restart2, 'bliving', cohort_bliving_ptr)
     call get_int_cohort_data(restart2, 'status', cohort_status_ptr)
     if(field_exists(restart2,'leaf_age')) &
          call get_cohort_data(restart2,'leaf_age',cohort_leaf_age_ptr)
     call get_cohort_data(restart2, 'npp_prev_day', cohort_npp_previous_day_ptr)

     if(field_exists(restart2,'landuse')) &
          call get_int_tile_data(restart2,'landuse',vegn_landuse_ptr)
     call get_tile_data(restart2,'age',vegn_age_ptr)
     if(field_exists(restart2,'fsc_pool_ag')) then
       call get_tile_data(restart2,'fsc_pool_ag',vegn_fsc_pool_ag_ptr)
       call get_tile_data(restart2,'fsc_rate_ag',vegn_fsc_rate_ag_ptr)
       call get_tile_data(restart2,'fsc_pool_bg',vegn_fsc_pool_bg_ptr)
       call get_tile_data(restart2,'fsc_rate_bg',vegn_fsc_rate_bg_ptr)
       call get_tile_data(restart2,'ssc_pool_ag',vegn_ssc_pool_ag_ptr)
       call get_tile_data(restart2,'ssc_rate_ag',vegn_ssc_rate_ag_ptr)
       call get_tile_data(restart2,'ssc_pool_bg',vegn_ssc_pool_bg_ptr)
       call get_tile_data(restart2,'ssc_rate_bg',vegn_ssc_rate_bg_ptr)
       call get_tile_data(restart2,'leaflitter_buffer_ag',vegn_leaflitter_buffer_ag_ptr)
       call get_tile_data(restart2,'coarsewoodlitter_buffer_ag',vegn_coarsewoodlitter_buffer_ag_ptr)
       call get_tile_data(restart2,'leaflitter_buffer_rate_ag',vegn_leaflitter_buffer_rate_ag_ptr)
       call get_tile_data(restart2,'coarsewoodlitter_buffer_rate_ag',vegn_coarsewoodlitter_buffer_rate_ag_ptr)
     else
       call get_tile_data(restart2,'fsc_pool',vegn_fsc_pool_bg_ptr)
       call get_tile_data(restart2,'fsc_rate',vegn_fsc_rate_bg_ptr)
       call get_tile_data(restart2,'ssc_pool',vegn_ssc_pool_bg_ptr)
       call get_tile_data(restart2,'ssc_rate',vegn_ssc_rate_bg_ptr)
     endif
     ! monthly-mean values
     call get_tile_data(restart2,'tc_av', vegn_tc_av_ptr)
     if(field_exists(restart2,'theta_av_phen')) then
        call get_tile_data(restart2,'theta_av_phen', vegn_theta_av_phen_ptr)
        call get_tile_data(restart2,'theta_av_fire', vegn_theta_av_fire_ptr)
        call get_tile_data(restart2,'psist_av', vegn_psist_av_ptr)
     else
        call get_tile_data(restart2,'theta_av', vegn_theta_av_phen_ptr)
        call get_tile_data(restart2,'theta_av', vegn_theta_av_fire_ptr)
        ! psist_av remains at initial value (equal to 0)
     endif
     call get_tile_data(restart2,'tsoil_av', vegn_tsoil_av_ptr)
     call get_tile_data(restart2,'precip_av', vegn_precip_av_ptr)
     call get_tile_data(restart2,'lambda', vegn_lambda_ptr)
     call get_tile_data(restart2,'fuel', vegn_fuel_ptr)
     ! annual-mean values
     call get_tile_data(restart2,'t_ann', vegn_t_ann_ptr)
     call get_tile_data(restart2,'t_cold', vegn_t_cold_ptr)
     call get_tile_data(restart2,'p_ann', vegn_p_ann_ptr)
     call get_tile_data(restart2,'ncm', vegn_ncm_ptr)
     ! accumulated values for annual averaging
     call get_tile_data(restart2,'t_ann_acm', vegn_t_ann_acm_ptr)
     call get_tile_data(restart2,'t_cold_acm', vegn_t_cold_acm_ptr)
     call get_tile_data(restart2,'p_ann_acm', vegn_p_ann_acm_ptr)
     call get_tile_data(restart2,'ncm_acm', vegn_ncm_acm_ptr)
     ! burned carbon pool and rate
     if(field_exists(restart2,'csmoke_pool')) &
          call get_tile_data(restart2,'csmoke_pool',vegn_csmoke_pool_ptr)
     if(field_exists(restart2,'csmoke_rate')) &
          call get_tile_data(restart2,'csmoke_rate',vegn_csmoke_rate_ptr)
     ! harvesting pools and rates
     do i = 1, N_HARV_POOLS
        if (field_exists(restart2,trim(HARV_POOL_NAMES(i))//'_harv_pool')) &
             call get_tile_data(restart2,trim(HARV_POOL_NAMES(i))//'_harv_pool',vegn_harv_pool_ptr,i)
        if (field_exists(restart2,trim(HARV_POOL_NAMES(i))//'_harv_rate')) &
             call get_tile_data(restart2,trim(HARV_POOL_NAMES(i))//'_harv_rate',vegn_harv_rate_ptr,i)
     enddo
  else
     call error_mesg('vegn_init',&
          'cold-starting vegetation',&
          NOTE)
  endif
  call free_land_restart(restart1)
  call free_land_restart(restart2)

  ! read climatological fields for initialization of species distribution
  if (file_exist('INPUT/biodata.nc'))then
     allocate(&
          t_ann (lnd%ls:lnd%le),&
          t_cold(lnd%ls:lnd%le),&
          p_ann (lnd%ls:lnd%le),&
          ncm   (lnd%ls:lnd%le) )
     call read_field( 'INPUT/biodata.nc','T_ANN', lnd%lon, lnd%lat, t_ann,  interp='nearest')
     call read_field( 'INPUT/biodata.nc','T_COLD',lnd%lon, lnd%lat, t_cold, interp='nearest')
     call read_field( 'INPUT/biodata.nc','P_ANN', lnd%lon, lnd%lat, p_ann,  interp='nearest')
     call read_field( 'INPUT/biodata.nc','NCM',   lnd%lon, lnd%lat, ncm,    interp='nearest')
     did_read_biodata = .TRUE.
     call error_mesg('vegn_init','did read INPUT/biodata.nc',NOTE)
  else
     did_read_biodata = .FALSE.
     call error_mesg('vegn_init','did NOT read INPUT/biodata.nc',NOTE)
  endif
  ! Go through all tiles and initialize the cohorts that have not been initialized yet --
  ! this allows to read partial restarts. Also initialize accumulation counters to zero
  ! or the values from the restarts.
  ce = first_elmt(land_tile_map, ls=lnd%ls)
  do while(loop_over_tiles(ce,tile,l))
     if (.not.associated(tile%vegn)) cycle

     tile%vegn%n_accum = n_accum
     tile%vegn%nmn_acm = nmn_acm

     if (tile%vegn%n_cohorts>0) cycle ! skip initialized tiles

     ! create and initialize cohorts for this vegetation tile
     ! for now, just create a new cohort with default values of biomasses
     tile%vegn%n_cohorts = 1
     allocate(tile%vegn%cohorts(tile%vegn%n_cohorts))
     cohort => tile%vegn%cohorts(1)
     cohort%Wl      = init_Wl
     cohort%Ws      = init_Ws
     cohort%Tv      = init_Tv

     cohort%bl      = init_cohort_bl
     cohort%blv     = init_cohort_blv
     cohort%br      = init_cohort_br
     cohort%bsw     = init_cohort_bsw
     cohort%bwood   = init_cohort_bwood
     cohort%bliving = cohort%bl+cohort%br+cohort%blv+cohort%bsw
     cohort%npp_previous_day = 0.0
     cohort%status  = LEAF_ON
     cohort%leaf_age = 0.0
     if(did_read_biodata.and.do_biogeography) then
        call update_species(cohort,t_ann(l),t_cold(l),p_ann(l),ncm(l),LU_NTRL)
        if (.not.biodata_bug) then
           tile%vegn%t_ann  = t_ann (l)
           tile%vegn%t_cold = t_cold(l)
           tile%vegn%p_ann  = p_ann (l)
           tile%vegn%ncm    = ncm   (l)
        endif
     else
        cohort%species = tile%vegn%tag
     endif
  enddo

  ! initialize carbon integrator
  call vegn_dynamics_init(id_ug,lnd%time,delta_time)

  ! initialize static vegetation
  call static_vegn_init ()
  call read_static_vegn (lnd%time)

  ! initialize harvesting options
  call vegn_harvesting_init()

  ! initialize distrurbances 
  call vegn_disturbance_init()
  
  ! initialize vegetation diagnostic fields
  call vegn_diag_init(id_ug,id_band,lnd%time)

  ! ---- diagnostic section
  ce = first_elmt(land_tile_map, ls=lnd%ls)
  do while(loop_over_tiles(ce, tile))
     if (.not.associated(tile%vegn)) cycle ! skip non-vegetation tiles
     call send_tile_data(id_vegn_type,  real(tile%vegn%tag), tile%diag)
  enddo

  if (allocated(t_ann))  deallocate(t_ann)
  if (allocated(t_cold)) deallocate(t_cold)
  if (allocated(p_ann))  deallocate(p_ann)
  if (allocated(ncm))    deallocate(ncm)

end subroutine vegn_init

! ============================================================================
subroutine vegn_diag_init(id_ug,id_band,time)
  integer        ,intent(in) :: id_ug   !<Unstructured axis id.
  integer        ,intent(in) :: id_band ! ID of spectral band axis
  type(time_type),intent(in) :: time    ! initial time for diagnostic fields

  ! ---- local vars
  integer :: i

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('soil')

  id_vegn_type = register_tiled_static_field ( module_name, 'vegn_type',  &
       (/id_ug/), 'vegetation type', missing_value=-1.0 )

  id_temp = register_tiled_diag_field ( module_name, 'temp',  &
       (/id_ug/), time, 'canopy temperature', 'degK', missing_value=-1.0 )
  id_wl = register_tiled_diag_field ( module_name, 'wl',  &
       (/id_ug/), time, 'canopy liquid water content', 'kg/m2', missing_value=-1.0 )
  id_ws = register_tiled_diag_field ( module_name, 'ws',  &
       (/id_ug/), time, 'canopy solid water content', 'kg/m2', missing_value=-1.0 )

  id_height = register_tiled_diag_field ( module_name, 'height',  &
       (/id_ug/), time, 'vegetation height', 'm', missing_value=-1.0 )
  id_lai    = register_tiled_diag_field ( module_name, 'lai',  &
       (/id_ug/), time, 'leaf area index', 'm2/m2', missing_value=-1.0 )
  id_lai_var = register_tiled_diag_field ( module_name, 'lai_var',  &
       (/id_ug/), time, 'variance of leaf area index across tiles in grid cell', 'm4/m4', &
       missing_value=-1.0 , op=OP_VAR)
  id_lai_std = register_tiled_diag_field ( module_name, 'lai_std',  &
       (/id_ug/), time, 'standard deviation of leaf area index across tiles in grid cell', 'm2/m2', &
       missing_value=-1.0, op=OP_STD)
  id_sai    = register_tiled_diag_field ( module_name, 'sai',  &
       (/id_ug/), time, 'stem area index', 'm2/m2', missing_value=-1.0 )
  id_leaf_size = register_tiled_diag_field ( module_name, 'leaf_size',  &
       (/id_ug/), time, missing_value=-1.0 )
  id_root_density = register_tiled_diag_field ( module_name, 'root_density',  &
       (/id_ug/), time, 'total biomass below ground', 'kg/m2', missing_value=-1.0 )
  id_root_zeta = register_tiled_diag_field ( module_name, 'root_zeta',  &
       (/id_ug/), time, 'e-folding depth of root biomass', 'm',missing_value=-1.0 )
  id_rs_min = register_tiled_diag_field ( module_name, 'rs_min',  &
       (/id_ug/), time, missing_value=-1.0 )
  id_leaf_refl = register_tiled_diag_field ( module_name, 'leaf_refl',  &
       (/id_ug,id_band/), time, 'reflectivity of leaf', missing_value=-1.0 )
  id_leaf_tran = register_tiled_diag_field ( module_name, 'leaf_tran',  &
       (/id_ug,id_band/), time, 'transmittance of leaf', missing_value=-1.0 )
  id_leaf_emis = register_tiled_diag_field ( module_name, 'leaf_emis',  &
       (/id_ug/), time, 'leaf emissivity', missing_value=-1.0 )
  id_snow_crit = register_tiled_diag_field ( module_name, 'snow_crit',  &
       (/id_ug/), time, missing_value=-1.0 )
  id_stomatal = register_tiled_diag_field ( module_name, 'stomatal_cond',  &
       (/id_ug/), time, 'vegetation stomatal conductance', missing_value=-1.0 )
  id_evap_demand = register_tiled_diag_field ( module_name, 'evap_demand',  &
       (/id_ug/), time, 'plant evaporative water demand',&
       'kg/(m2 s)', missing_value=-1e20 )
  id_an_op = register_tiled_diag_field ( module_name, 'an_op',  &
       (/id_ug/), time, 'net photosynthesis with open stomata', &
       '(mol CO2)(m2 of leaf)^-1 year^-1', missing_value=-1e20 )
  id_an_cl = register_tiled_diag_field ( module_name, 'an_cl',  &
       (/id_ug/), time, 'net photosynthesis with closed stomata', &
       '(mol CO2)(m2 of leaf)^-1 year^-1', missing_value=-1e20 )

  id_bl = register_tiled_diag_field ( module_name, 'bl',  &
       (/id_ug/), time, 'biomass of leaves', 'kg C/m2', missing_value=-1.0 )
  id_blv = register_tiled_diag_field ( module_name, 'blv',  &
       (/id_ug/), time, 'biomass in labile store', 'kg C/m2', missing_value=-1.0 )
  id_br = register_tiled_diag_field ( module_name, 'br',  &
       (/id_ug/), time, 'biomass of fine roots', 'kg C/m2', missing_value=-1.0 )
  id_bsw = register_tiled_diag_field ( module_name, 'bsw',  &
       (/id_ug/), time, 'biomass of sapwood', 'kg C/m2', missing_value=-1.0 )
  id_bwood = register_tiled_diag_field ( module_name, 'bwood',  &
       (/id_ug/), time, 'biomass of heartwood', 'kg C/m2', missing_value=-1.0 )
  id_btot = register_tiled_diag_field ( module_name, 'btot',  &
       (/id_ug/), time, 'total biomass', 'kg C/m2', missing_value=-1.0 )
  id_fuel = register_tiled_diag_field ( module_name, 'fuel',  &
       (/id_ug/), time, 'mass of fuel', 'kg C/m2', missing_value=-1.0 )
  id_lambda = register_tiled_diag_field (module_name, 'lambda',(/id_ug/), &
       time, 'drought', 'months', missing_value=-100.0)

  id_species = register_tiled_diag_field ( module_name, 'species',  &
       (/id_ug/), time, 'vegetation species number', missing_value=-1.0 )
  id_status = register_tiled_diag_field ( module_name, 'status',  &
       (/id_ug/), time, 'status of leaves', missing_value=-1.0 )
  id_theph = register_tiled_diag_field ( module_name, 'theph',  &
       (/id_ug/), time, 'theta for phenology', missing_value=-1.0 )
  id_psiph = register_tiled_diag_field ( module_name, 'psiph',  &
       (/id_ug/), time, 'psi stress for phenology', missing_value=-1.0 )
  id_leaf_age = register_tiled_diag_field ( module_name, 'leaf_age',  &
       (/id_ug/), time, 'age of leaves since bud burst', 'days', missing_value=-1.0 )!ens

  id_con_v_h = register_tiled_diag_field ( module_name, 'con_v_h', (/id_ug/), &
       time, 'conductance for sensible heat between canopy and canopy air', &
       'm/s', missing_value=-1.0 )
  id_con_v_v = register_tiled_diag_field ( module_name, 'con_v_v', (/id_ug/), &
       time, 'conductance for water vapor between canopy and canopy air', &
       'm/s', missing_value=-1.0 )

  id_cgain = register_tiled_diag_field ( module_name, 'cgain', (/id_ug/), &
       time, 'carbon gain', 'kg C/m2', missing_value=-100.0 )
  id_closs = register_tiled_diag_field ( module_name, 'closs', (/id_ug/), &
       time, 'carbon loss', 'kg C/m2', missing_value=-100.0 )
  id_wdgain = register_tiled_diag_field ( module_name, 'wdgain', (/id_ug/), &
       time, 'wood biomass gain', 'kg C/m2', missing_value=-100.0 )

  id_t_ann  = register_tiled_diag_field ( module_name, 't_ann', (/id_ug/), &
       time, 'annual mean temperature', 'degK', missing_value=-999.0 )
  id_t_cold  = register_tiled_diag_field ( module_name, 't_cold', (/id_ug/), &
       time, 'average temperature of the coldest month', 'degK', missing_value=-999.0 )
  id_p_ann  = register_tiled_diag_field ( module_name, 'p_ann', (/id_ug/), &
       time, 'annual mean precipitation', 'kg/(m2 s)', missing_value=-999.0 )
  id_ncm = register_tiled_diag_field ( module_name, 'ncm', (/id_ug/), &
       time, 'number of cold months', 'dimensionless', missing_value=-999.0 )

  id_t_harv_pool = register_tiled_diag_field( module_name, 'harv_pool', (/id_ug/), &
       time, 'total harvested carbon', 'kg C/m2', missing_value=-999.0)
  id_t_harv_rate = register_tiled_diag_field( module_name, 'harv_rate', (/id_ug/), &
       time, 'total rate of release of harvested carbon to the atmosphere', &
       'kg C/(m2 year)', missing_value=-999.0)
  do i = 1,N_HARV_POOLS
     id_harv_pool(i) = register_tiled_diag_field( module_name, &
          trim(HARV_POOL_NAMES(i))//'_harv_pool', (/id_ug/), time, &
          'harvested carbon', 'kg C/m2', missing_value=-999.0)
     id_harv_rate(i) = register_tiled_diag_field( module_name, &
          trim(HARV_POOL_NAMES(i))//'_harv_rate', (/id_ug/), time, &
          'rate of release of harvested carbon to the atmosphere', 'kg C/(m2 year)', &
          missing_value=-999.0)
  enddo

  id_fsc_pool_ag = register_tiled_diag_field (module_name, 'fsc_pool_ag', (/id_ug/), &
       time, 'intermediate pool of above-ground fast soil carbon', 'kg C/m2', missing_value=-999.0)
  id_fsc_rate_ag = register_tiled_diag_field (module_name, 'fsc_rate_ag', (/id_ug/), &
       time, 'rate of conversion of above-ground fsc_pool to the fast soil_carbon', 'kg C/(m2 yr)', &
       missing_value=-999.0)
  id_ssc_pool_ag = register_tiled_diag_field (module_name, 'ssc_pool_ag', (/id_ug/), &
       time, 'intermediate pool of above-ground slow soil carbon', 'kg C/m2', missing_value=-999.0)
  id_ssc_rate_ag = register_tiled_diag_field (module_name, 'ssc_rate_ag', (/id_ug/), &
       time, 'rate of conversion of above-ground ssc_pool to the fast soil_carbon', 'kg C/(m2 yr)', &
       missing_value=-999.0)
  id_leaflitter_buffer_ag = register_tiled_diag_field (module_name, 'leaflitter_buffer_ag', (/id_ug/), &
       time, 'intermediate pool of leaf litter carbon', 'kg C/m2', missing_value=-999.0)
  id_leaflitter_buffer_rate_ag = register_tiled_diag_field (module_name, 'leaflitter_buffer_rate_ag', (/id_ug/), &
       time, 'rate of conversion of above-ground leaf litter buffer to the fast soil_carbon', 'kg C/(m2 yr)', &
       missing_value=-999.0)
  id_coarsewoodlitter_buffer_ag = register_tiled_diag_field (module_name, 'coarsewoodlitter_buffer_ag', (/id_ug/), &
       time, 'intermediate pool of coarsewood litter carbon', 'kg C/m2', missing_value=-999.0)
  id_coarsewoodlitter_buffer_rate_ag = register_tiled_diag_field (module_name, 'coarsewoodlitter_buffer_rate_ag', (/id_ug/), &
       time, 'rate of conversion of above-ground coarsewood litter buffer to the fast soil_carbon', 'kg C/(m2 yr)', &
       missing_value=-999.0)
  id_fsc_pool_bg = register_tiled_diag_field (module_name, 'fsc_pool_bg', (/id_ug/), &
       time, 'intermediate pool of below-ground fast soil carbon', 'kg C/m2', missing_value=-999.0)
  id_fsc_rate_bg = register_tiled_diag_field (module_name, 'fsc_rate_bg', (/id_ug/), &
       time, 'rate of conversion of below-ground fsc_pool to the fast soil_carbon', 'kg C/(m2 yr)', &
       missing_value=-999.0)
  id_ssc_pool_bg = register_tiled_diag_field (module_name, 'ssc_pool_bg', (/id_ug/), &
       time, 'intermediate pool of below-ground slow soil carbon', 'kg C/m2', missing_value=-999.0)
  id_ssc_rate_bg = register_tiled_diag_field (module_name, 'ssc_rate_bg', (/id_ug/), &
       time, 'rate of conversion of below-ground ssc_pool to the fast soil_carbon', 'kg C/(m2 yr)', &
       missing_value=-999.0)

  id_csmoke_pool = register_tiled_diag_field ( module_name, 'csmoke', (/id_ug/), &
       time, 'carbon lost through fire', 'kg C/m2', missing_value=-999.0)
  id_csmoke_rate = register_tiled_diag_field ( module_name, 'csmoke_rate', (/id_ug/), &
       time, 'rate of release of carbon lost through fire to the atmosphere', &
       'kg C/(m2 yr)', missing_value=-999.0)

  id_ssc_in = register_tiled_diag_field ( module_name, 'ssc_in',  (/id_ug/), &
     time,  'soil slow carbon in', 'kg C/m2', missing_value=-999.0 )
  id_ssc_out = register_tiled_diag_field ( module_name, 'ssc_out',  (/id_ug/), &
     time,  'soil slow carbon out', 'kg C/m2', missing_value=-999.0 )
  id_deadmic_out = register_tiled_diag_field ( module_name, 'deadmic_out',  (/id_ug/), &
     time,  'daed microbe carbon out', 'kg C/m2', missing_value=-999.0 )
  id_fsc_in = register_tiled_diag_field ( module_name, 'fsc_in',  (/id_ug/), &
     time,  'soil fast carbon in', 'kg C/m2', missing_value=-999.0 )
  id_fsc_out = register_tiled_diag_field ( module_name, 'fsc_out',  (/id_ug/), &
     time,  'soil fast carbon out', 'kg C/m2', missing_value=-999.0 )
  id_veg_in = register_tiled_diag_field ( module_name, 'veg_in',  (/id_ug/), &
     time,  'vegetation carbon in', 'kg C/m2', missing_value=-999.0 )
  id_veg_out = register_tiled_diag_field ( module_name, 'veg_out',  (/id_ug/), &
     time,  'vegetation carbon out', 'kg C/m2', missing_value=-999.0 )

  id_afire = register_tiled_diag_field (module_name, 'afire', (/id_ug/), &
       time, 'area been fired', missing_value=-100.0)
  id_atfall = register_tiled_diag_field (module_name, 'atfall',(/id_ug/), &
       time, 'area been disturbed', missing_value=-100.0)

  id_phot_co2 = register_tiled_diag_field (module_name, 'qco2_phot',(/id_ug/), &
       time, 'CO2 mixing ratio for photosynthesis calculations', 'mol CO2/mol dry air', &
       missing_value=-1.0)

  ! CMOR variables
  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('land')
  call add_tiled_diag_field_alias(id_lai, cmor_name, 'lai', (/id_ug/), &
       time, 'Leaf Area Index', '1', missing_value = -1.0, &
       standard_name = 'leaf_area_index', fill_missing = .TRUE.)
  call add_tiled_diag_field_alias(id_lai, cmor_name, 'laiLut', (/id_ug/), &
       time, 'leaf area index on land use tile', '1', missing_value = -1.0, &
       standard_name = 'leaf_area_index_lut', fill_missing = .FALSE.)
  call add_tiled_diag_field_alias ( id_btot, cmor_name, 'cVeg', (/id_ug/), &
       time, 'Carbon Mass in Vegetation', 'kg C m-2', missing_value=-1.0, &
       standard_name='vegetation_carbon_content', fill_missing=.TRUE.)
  call add_tiled_diag_field_alias (id_btot, cmor_name, 'cVegLut', (/id_ug/), &
       time, 'Carbon Mass in Vegetation', 'kg C m-2', missing_value=-1.0, &
       standard_name='vegetation_carbon_content', fill_missing=.FALSE.)
  id_cProduct = register_tiled_diag_field( cmor_name, 'cProduct', (/id_ug/), &
       time, 'Carbon in Products of Land Use Change', 'kg C m-2', missing_value=-999.0, &
       standard_name='carbon_in_producs_of_luc', fill_missing=.TRUE.)
  id_cAnt = register_tiled_diag_field( cmor_name, 'cAnt', (/id_ug/), &
       time, 'Carbon in Anthropogenic Pool', 'kg C m-2', missing_value=-999.0, &
       fill_missing=.TRUE.) ! standard_name not known at this time
  call add_tiled_diag_field_alias(id_cAnt, cmor_name, 'cAntLut', (/id_ug/), &
       time, 'Carbon in Anthropogenic Pools Associated with Land Use Tiles', 'kg C m-2', &
       missing_value=-999.0, fill_missing = .TRUE.) ! standard_name not known at this time
  id_fFire = register_tiled_diag_field ( cmor_name, 'fFire', (/id_ug/), &
       time, 'CO2 Emission from Fire', 'kg C m-2 s-1', missing_value=-1.0, &
       standard_name='co2_emission_from_fire', fill_missing=.TRUE.)
  id_fGrazing = register_tiled_diag_field( cmor_name, 'fGrazing', (/id_ug/), &
       time, 'CO2 flux to Atmosphere from Grazing', 'kg C m-2 s-1', missing_value=-1.0, &
       standard_name='co2_flux_to_atmosphere_from_grazing', fill_missing=.TRUE.)
  id_fHarvest = register_tiled_diag_field( cmor_name, 'fHarvest', (/id_ug/), &
       time, 'CO2 flux to Atmosphere from Crop Harvesting', 'kg C m-2 s-1', missing_value=-1.0, &
       standard_name='co2_flux_to_atmosphere_from_crop_harvesting', fill_missing=.TRUE.)
  id_fLuc = register_tiled_diag_field( cmor_name, 'fLuc', (/id_ug/), &
       time, 'CO2 flux to Atmosphere from Land Use Change', 'kg C m-2 s-1', missing_value=-1.0, &
       standard_name='co2_flux_to_atmosphere_from_land_use_change', fill_missing=.TRUE.)
  id_cLeaf = register_tiled_diag_field ( cmor_name, 'cLeaf',  (/id_ug/), &
       time, 'Carbon in Leaves', 'kg C m-2', missing_value=-1.0, &
       standard_name='carbon_in_leaves', fill_missing=.TRUE.)
  id_cWood = register_tiled_diag_field ( cmor_name, 'cWood',  (/id_ug/), &
       time, 'Carbon in Wood', 'kg C m-2', missing_value=-1.0, &
       standard_name='carbon_in_wood', fill_missing=.TRUE.)
  id_cRoot = register_tiled_diag_field ( cmor_name, 'cRoot',  (/id_ug/), &
       time, 'Carbon in Roots', 'kg C m-2', missing_value=-1.0, &
       standard_name='carbon_in_roots', fill_missing=.TRUE.)
  id_cMisc = register_tiled_diag_field ( cmor_name, 'cMisc',  (/id_ug/), &
       time, 'Carbon in Other Living Compartments', 'kg C m-2', missing_value=-1.0, &
       standard_name='carbon_in_other_living_compartments', fill_missing=.TRUE.)
end subroutine


! ============================================================================
! write restart file and release memory
subroutine vegn_end ()

  module_is_initialized =.FALSE.

  ! finalize harvesting
  call vegn_harvesting_end ()

  ! finalize static vegetation, if necessary
  call static_vegn_end ()
end subroutine vegn_end


! ============================================================================
subroutine save_vegn_restart(tile_dim_length,timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars
  integer ::  i
  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile
  integer :: n_accum, nmn_acm

  character(267) :: filename
  type(land_restart_type) :: restart1, restart2 ! restart file i/o object

  call error_mesg('vegn_end','writing NetCDF restart',NOTE)

  ! create output file, including internal structure necessary for tile output
  filename = trim(timestamp)//'vegn1.res.nc'
  call init_land_restart(restart1, filename, vegn_tile_exists, tile_dim_length)

  ! create compressed dimension for vegetation cohorts -- must be called even
  ! if restart has not been created, because it calls mpp_max and that should
  ! be called on all PEs to work
  call create_cohort_dimension(restart1)

  call add_cohort_data(restart1,'tv',cohort_tv_ptr,'vegetation temperature','degrees_K')
  call add_cohort_data(restart1,'wl',cohort_wl_ptr,'vegetation liquid water content','kg/m2')
  call add_cohort_data(restart1,'ws',cohort_ws_ptr,'vegetation solid water content','kg/m2')
  call save_land_restart(restart1)
  call free_land_restart(restart1)


  filename = trim(timestamp)//'vegn2.res.nc'
  call init_land_restart(restart2, filename, vegn_tile_exists, tile_dim_length)
  ! create compressed dimension for vegetation cohorts -- see note above
  call create_cohort_dimension(restart2)

  ! store global variables
  ! find first tile and get n_accum and nmn_acm from it
  n_accum = 0; nmn_acm = 0
  ce = first_elmt(land_tile_map)
  do while ( loop_over_tiles(ce,tile))
     if(associated(tile%vegn)) then
        n_accum = tile%vegn%n_accum
        nmn_acm = tile%vegn%nmn_acm
     endif
  enddo
  ! n_accum and nmn_acm are currently the same for all tiles; we only call mpp_max
  ! to handle the situation when there are no tiles in the current domain
  call mpp_max(n_accum); call mpp_max(nmn_acm)
  call add_scalar_data(restart2,'n_accum',n_accum,'number of accumulated steps')

  call add_scalar_data(restart2,'nmn_acm',nmn_acm,'number of accumulated months')
  call add_int_cohort_data(restart2,'species', cohort_species_ptr, 'vegetation species')
  call add_cohort_data(restart2,'hite', cohort_height_ptr, 'vegetation height','m')
  call add_cohort_data(restart2,'bl', cohort_bl_ptr, 'biomass of leaves per individual','kg C/m2')
  call add_cohort_data(restart2,'blv', cohort_blv_ptr, 'biomass of virtual leaves (labile store) per individual','kg C/m2')
  call add_cohort_data(restart2,'br', cohort_br_ptr, 'biomass of fine roots per individual','kg C/m2')
  call add_cohort_data(restart2,'bsw', cohort_bsw_ptr, 'biomass of sapwood per individual','kg C/m2')
  call add_cohort_data(restart2,'bwood', cohort_bwood_ptr, 'biomass of heartwood per individual','kg C/m2')
  call add_cohort_data(restart2,'bliving', cohort_bliving_ptr, 'total living biomass per individual','kg C/m2')
  call add_int_cohort_data(restart2,'status', cohort_status_ptr, 'leaf status')
  call add_cohort_data(restart2,'leaf_age',cohort_leaf_age_ptr, 'age of leaves since bud burst', 'days')

  call add_cohort_data(restart2,'npp_prev_day', cohort_npp_previous_day_ptr, 'previous day NPP','kg C/(m2 year)')

  call add_int_tile_data(restart2,'landuse',vegn_landuse_ptr,'vegetation land use type')
  call add_tile_data(restart2,'age',vegn_age_ptr,'vegetation age', 'yr')
  call add_tile_data(restart2,'fsc_pool_ag',vegn_fsc_pool_ag_ptr, &
       'intermediate pool for AG fast soil carbon input', 'kg C/m2')
  call add_tile_data(restart2,'fsc_rate_ag',vegn_fsc_rate_ag_ptr, &
       'conversion rate of AG fsc_pool to fast soil carbon', 'kg C/(m2 yr)')
  call add_tile_data(restart2,'ssc_pool_ag',vegn_ssc_pool_ag_ptr, &
       'intermediate pool for AG slow soil carbon input', 'kg C/m2')
  call add_tile_data(restart2,'ssc_rate_ag',vegn_ssc_rate_ag_ptr, &
       'conversion rate of AG ssc_pool to slow soil carbon', 'kg C/(m2 yr)')
  call add_tile_data(restart2,'fsc_pool_bg',vegn_fsc_pool_bg_ptr, &
       'intermediate pool for BG fast soil carbon input', 'kg C/m2')
  call add_tile_data(restart2,'fsc_rate_bg',vegn_fsc_rate_bg_ptr, &
       'conversion rate of BG fsc_pool to fast soil carbon', 'kg C/(m2 yr)')
  call add_tile_data(restart2,'ssc_pool_bg',vegn_ssc_pool_bg_ptr, &
       'intermediate pool for BG slow soil carbon input', 'kg C/m2')
  call add_tile_data(restart2,'ssc_rate_bg',vegn_ssc_rate_bg_ptr, &
       'conversion rate of BG ssc_pool to slow soil carbon', 'kg C/(m2 yr)')

  call add_tile_data(restart2,'leaflitter_buffer_ag',vegn_leaflitter_buffer_ag_ptr, &
       'intermediate pool for AG leaf litter carbon input', 'kg C/m2')
  call add_tile_data(restart2,'leaflitter_buffer_rate_ag',vegn_leaflitter_buffer_rate_ag_ptr, &
       'conversion rate of AG leaf litter to litter carbon pool', 'kg C/(m2 yr)')
  call add_tile_data(restart2,'coarsewoodlitter_buffer_ag',vegn_coarsewoodlitter_buffer_ag_ptr, &
       'intermediate pool for AG coarsewood litter carbon input', 'kg C/m2')
  call add_tile_data(restart2,'coarsewoodlitter_buffer_rate_ag',vegn_coarsewoodlitter_buffer_rate_ag_ptr, &
       'conversion rate of AG coarsewood litter to litter carbon pool', 'kg C/(m2 yr)')

  ! monthly-mean values
  call add_tile_data(restart2,'tc_av', vegn_tc_av_ptr,'average canopy air temperature','degK')
  call add_tile_data(restart2,'theta_av_phen', vegn_theta_av_phen_ptr,'average soil moisture for phenology')
  call add_tile_data(restart2,'theta_av_fire', vegn_theta_av_fire_ptr,'average soil moisture for fire')
  call add_tile_data(restart2,'psist_av', vegn_psist_av_ptr,'average soil-water-stress index')
  call add_tile_data(restart2,'tsoil_av', vegn_tsoil_av_ptr,'average bulk soil temperature for soil carbon','degK')
  call add_tile_data(restart2,'precip_av', vegn_precip_av_ptr,'average total precipitation','kg/(m2 s)')
  call add_tile_data(restart2,'lambda', vegn_lambda_ptr,'dryness parameter')
  call add_tile_data(restart2,'fuel', vegn_fuel_ptr,'fuel density','kg C/m2')
  ! annual-mean values
  call add_tile_data(restart2,'t_ann', vegn_t_ann_ptr,'average annual canopy air temperature','degK')
  call add_tile_data(restart2,'t_cold', vegn_t_cold_ptr,'average canopy air temperature of coldest month','degK')
  call add_tile_data(restart2,'p_ann', vegn_p_ann_ptr,'average annual precipitation','kg/(m2 s)')
  call add_tile_data(restart2,'ncm', vegn_ncm_ptr,'number of cold months')
  ! accumulated values for annual averaging
  call add_tile_data(restart2,'t_ann_acm', vegn_t_ann_acm_ptr,'accumulated annual canopy air temperature','degK')
  call add_tile_data(restart2,'t_cold_acm', vegn_t_cold_acm_ptr,'accumulated temperature of coldest month','degK')
  call add_tile_data(restart2,'p_ann_acm', vegn_p_ann_acm_ptr,'accumulated precipitation','kg/(m2 s)')
  call add_tile_data(restart2,'ncm_acm', vegn_ncm_acm_ptr,'accumulated number of cold months')

  ! burned carbon pool and rate
  call add_tile_data(restart2,'csmoke_pool',vegn_csmoke_pool_ptr,'carbon lost through fires', 'kg C/m2')
  call add_tile_data(restart2,'csmoke_rate',vegn_csmoke_rate_ptr,'rate of release of carbon lost through fires to the atmosphere', 'kg C/(m2 yr)')

  ! harvesting pools and rates
  do i = 1, N_HARV_POOLS
     call add_tile_data(restart2, trim(HARV_POOL_NAMES(i))//'_harv_pool', &
          vegn_harv_pool_ptr, i, 'harvested carbon','kg C/m2')
     call add_tile_data(restart2, trim(HARV_POOL_NAMES(i))//'_harv_rate', &
          vegn_harv_rate_ptr, i, 'rate of release of harvested carbon to the atmosphere','kg C/(m2 yr)')
  enddo

  call save_land_restart(restart2)
  call free_land_restart(restart2)
end subroutine save_vegn_restart

! ============================================================================
subroutine vegn_get_cover(vegn, snow_depth, vegn_cover)
  type(vegn_tile_type), intent(inout)  :: vegn ! it is only inout because vegn%data%cover
                                    ! changes cohort; can it be avoided?
  real,                 intent(in)  :: snow_depth
  real,                 intent(out) :: vegn_cover

  real :: vegn_cover_snow_factor

  call vegn_data_cover(vegn%cohorts(1), snow_depth, vegn_cover, vegn_cover_snow_factor)

end subroutine vegn_get_cover


! ============================================================================
subroutine vegn_properties ( vegn, vegn_cover, vegn_height, vegn_lai, vegn_sai, vegn_d_leaf)
  type(vegn_tile_type), intent(in) :: vegn
  real, intent(out) :: &
       vegn_cover, vegn_height, vegn_lai, vegn_sai, vegn_d_leaf

  vegn_cover  = vegn%cohorts(1)%cover
  vegn_lai    = vegn%cohorts(1)%lai
  vegn_sai    = vegn%cohorts(1)%sai
  vegn_height = vegn%cohorts(1)%height
  vegn_d_leaf = vegn%cohorts(1)%leaf_size

end subroutine vegn_properties


! ============================================================================
subroutine vegn_step_1 ( vegn, soil, diag, &
        p_surf, ustar, drag_q, &
        SWdn, RSv, precip_l, precip_s, &
        land_d, land_z0s, land_z0m, grnd_z0s, &
        soil_beta, soil_water_supply, &
        cana_T, cana_q, cana_co2_mol, &
        ! output
        con_g_h, con_g_v, & ! aerodynamic conductance between canopy air and canopy, for heat and vapor flux
        con_v_h, con_v_v, & ! aerodynamic conductance between canopy and canopy air, for heat and vapor
        stomatal_cond, &    ! integral stomatal conductance of canopy (that is, multiplied by LAI), for water vapor, m/s
        vegn_T,vegn_Wl,  vegn_Ws,           & ! temperature, water and snow mass of the canopy
        vegn_ifrac,                         & ! intercepted fraction of liquid and frozen precipitation
        vegn_lai,                           & ! leaf area index
        drip_l, drip_s,                     & ! water and snow drip rate from precipitation, kg/(m2 s)
        vegn_hcap,                          & ! vegetation heat capacity
        Hv0,   DHvDTv,   DHvDTc,            & ! sens heat flux
        Et0,   DEtDTv,   DEtDqc,   DEtDwl,   DEtDwf,  & ! transpiration
        Eli0,  DEliDTv,  DEliDqc,  DEliDwl,  DEliDwf, & ! evaporation of intercepted water
        Efi0,  DEfiDTv,  DEfiDqc,  DEfiDwl,  DEfiDwf  ) ! sublimation of intercepted snow
  type(vegn_tile_type), intent(inout) :: vegn ! vegetation data
  type(soil_tile_type), intent(inout) :: soil ! soil data
  ! TODO: possibly move calculation of soil-related stuff from calling subroutine to here
  !       now since we have soil tiled passed to us
  type(diag_buff_type), intent(inout) :: diag ! diagnostic buffer
  real, intent(in) :: &
       p_surf,    & ! surface pressure, N/m2
       ustar,     & ! friction velocity, m/s
       drag_q,    & ! bulk drag coefficient for specific humidity
       SWdn(NBANDS), & ! downward SW radiation at the top of the canopy, W/m2
       RSv (NBANDS), & ! net SW radiation balance of the canopy, W/m2
       precip_l, precip_s, & ! liquid and solid precipitation rates, kg/(m2 s)
       land_d, land_z0s, land_z0m, & ! land displacement height and roughness, m
       grnd_z0s, & ! roughness of ground surface (including snow effect)
       soil_beta, & ! relative water availability
       soil_water_supply, & ! max rate of water supply to the roots, kg/(m2 s)
       cana_T,    & ! temperature of canopy air, deg K
       cana_q,    & ! specific humidity of canopy air, kg/kg
       cana_co2_mol ! co2 mixing ratio in the canopy air, mol CO2/mol dry air
  ! output -- coefficients of linearized expressions for fluxes
  real, intent(out) ::   &
       vegn_T,vegn_Wl,  vegn_Ws,& ! temperature, water and snow mass of the canopy
       vegn_ifrac, & ! intercepted fraction of liquid and frozen precipitation
       vegn_lai, & ! vegetation leaf area index
       drip_l, drip_s, & ! water and snow drip rate from precipitation, kg/(m2 s)
       vegn_hcap, & ! total vegetation heat capacity, including intercepted water and snow
       con_g_h, con_g_v, & ! aerodynamic conductance between ground and canopy air, for heat and vapor
       con_v_h, con_v_v, & ! aerodynamic conductance between canopy and canopy air, for heat and vapor
       stomatal_cond, &    ! integral stomatal conductance of canopy (that is, multiplied by LAI), for water vapor, m/s
       Hv0,   DHvDTv,   DHvDTc, & ! sens heat flux
       Et0,   DEtDTv,   DEtDqc,   DEtDwl,   DEtDwf,  & ! transpiration
       Eli0,  DEliDTv,  DEliDqc,  DEliDwl,  DEliDwf, & ! evaporation of intercepted water
       Efi0,  DEfiDTv,  DEfiDqc,  DEfiDwl,  DEfiDwf    ! sublimation of intercepted snow

  ! ---- local vars
  real :: &
       ft,DftDwl,DftDwf, & ! fraction of canopy not covered by intercepted water/snow, and its
                    ! derivatives w.r.t. intercepted water masses
       fw,DfwDwl,DfwDwf, & ! fraction of canopy covered by intercepted water, and its
                    ! derivatives w.r.t. intercepted water masses
       fs,DfsDwl,DfsDwf, & ! fraction of canopy covered by intercepted snow, and its
                    ! derivatives w.r.t. intercepted water masses
       rav_lit,   & ! additional resistance of litter to vapor transport
       total_cond, &! overall conductance from inside stomata to canopy air
       qvsat,     & ! sat. specific humidity at the leaf T
       DqvsatDTv, & ! derivative of qvsat w.r.t. leaf T
       rho,       & ! density of canopy air
       phot_co2,  & ! co2 mixing ratio for photosynthesis, mol CO2/mol dry air
       evap_demand, & ! evaporative water demand, kg/(m2 s)
       photosynt, & ! photosynthesis
       photoresp    ! photo-respiration
  real :: litter_fast_C, litter_slow_C, litter_deadmic_C ! For rav_lit calculations
  type(vegn_cohort_type), pointer :: cohort

  ! get the pointer to the first (and, currently, the only) cohort
  cohort => vegn%cohorts(1)

  if(is_watch_point()) then
     write(*,*)'#### vegn_step_1 input ####'
     __DEBUG3__(p_surf, ustar, drag_q)
     __DEBUG1__(SWdn)
     __DEBUG1__(RSv)
     __DEBUG2__(precip_l, precip_s)
     __DEBUG4__(land_d, land_z0s, land_z0m, grnd_z0s)
     __DEBUG2__(soil_beta, soil_water_supply)
     __DEBUG3__(cana_T, cana_q, cana_co2_mol)
     write(*,*)'#### end of vegn_step_1 input ####'
     __DEBUG3__(cohort%height, cohort%lai, cohort%sai)
     __DEBUG2__(cohort%cover,cohort%leaf_size)
     __DEBUG1__(cohort%Tv)
  endif

  ! check the range of input temperature
  call check_temp_range(cohort%Tv,'vegn_step_1','cohort%Tv')

  ! calculate the fractions of intercepted precipitation
  vegn_ifrac = cohort%cover

  ! get the lai
  vegn_lai = cohort%lai

  ! calculate the aerodynamic conductance coefficients
  call cana_turbulence(ustar, &
     cohort%cover, cohort%height, cohort%lai, cohort%sai, cohort%leaf_size, &
     land_d, land_z0m, land_z0s, grnd_z0s, &
     con_v_h, con_v_v, con_g_h, con_g_v)

  ! take into account additional resistance of litter to the water vapor flux.
  ! not a good parameterization, but just using for sensitivity analyses now.
  ! ignores differing biomass and litter turnover rates.
  select case(soil_carbon_option)
  case(SOILC_CENTURY,SOILC_CENTURY_BY_LAYER)
     litter_fast_C    = soil%fast_soil_C(1)
     litter_slow_C    = soil%slow_soil_C(1)
     litter_deadmic_C = 0.0
  case(SOILC_CORPSE)
     call poolTotalCarbon(soil%leafLitter,fastC=litter_fast_C,slowC=litter_slow_C,deadMicrobeC=litter_deadmic_C)
  case default
     call error_mesg('vegn_step_1','The value of soil_carbon_option is invalid. This should never happen. Contact developer.',FATAL)
  end select
  rav_lit = rav_lit_0 + rav_lit_vi * (cohort%lai+cohort%sai) &
                      + rav_lit_fsc * litter_fast_C &
                      + rav_lit_ssc * litter_slow_C &
                      + rav_lit_deadmic * litter_deadmic_C &
                      + rav_lit_bwood * cohort%bwood
  con_g_v = con_g_v/(1.0+rav_lit*con_g_v)

  ! calculate the vegetation photosynthesis and associated stomatal conductance
  if (vegn_phot_co2_option == VEGN_PHOT_CO2_INTERACTIVE) then
     phot_co2 = cana_co2_mol
  else
     phot_co2 = co2_for_photosynthesis
  endif
  call vegn_photosynthesis ( vegn, &
     SWdn(BAND_VIS), RSv(BAND_VIS), cana_q, phot_co2, p_surf, drag_q, &
     soil_beta, soil_water_supply, &
     evap_demand, stomatal_cond, photosynt, photoresp )

  call get_vegn_wet_frac ( cohort, fw, DfwDwl, DfwDwf, fs, DfsDwl, DfsDwf )
  ! transpiring fraction and its derivatives
  ft     = 1 - fw - fs
  DftDwl = - DfwDwl - DfsDwl
  DftDwf = - DfwDwf - DfsDwf
  call qscomp(cohort%Tv, p_surf, qvsat, DqvsatDTv)

  rho = p_surf/(rdgas*cana_T *(1+d608*cana_q))

  ! get the vegetation temperature
  vegn_T  =  cohort%Tv
  ! get the amount of intercepted water and snow
  vegn_Wl =  cohort%Wl
  vegn_Ws =  cohort%Ws
  ! calculate the drip rates
  drip_l  = max(vegn_Wl,0.0)/tau_drip_l
  drip_s  = max(vegn_Ws,0.0)/tau_drip_s
  ! correct the drip rates so that the amount of water and snow accumulated over time step
  ! is no larger then the canopy water-holding capacity
  drip_l = max((vegn_Wl+precip_l*delta_time*vegn_ifrac-cohort%Wl_max)/delta_time,drip_l)
  drip_s = max((vegn_Ws+precip_s*delta_time*vegn_ifrac-cohort%Ws_max)/delta_time,drip_s)

  ! calculate the total heat capacity
  call vegn_data_heat_capacity (cohort, vegn_hcap)
  vegn_hcap = vegn_hcap + clw*cohort%Wl + csw*cohort%Ws
  ! calculate the coefficient of sensible heat flux linearization
  Hv0     =  2*rho*cp_air*con_v_h*(cohort%Tv - cana_T)
  DHvDTv  =  2*rho*cp_air*con_v_h
  DHvDTc  = -2*rho*cp_air*con_v_h
  ! calculate the coefficients of the transpiration linearization
  if(con_v_v==0.and.stomatal_cond==0) then
     total_cond = 0.0
  else
     total_cond = stomatal_cond*con_v_v/(stomatal_cond+con_v_v)
  endif

  if(qvsat>cana_q)then
     ! flux is directed from the surface: transpiration is possible, and the
     ! evaporation of intercepted water depends on the fraction of wet/snow
     ! covered canopy.

     ! prohibit transpiration if leaf temperature below some predefined minimum
     ! typically (268K, but check namelist)
     if(cohort%Tv < T_transp_min) total_cond = 0
     ! calculate the transpiration linearization coefficients
     Et0     =  rho*total_cond*ft*(qvsat - cana_q)
     DEtDTv  =  rho*total_cond*ft*DqvsatDTv
     DEtDqc  = -rho*total_cond*ft
     DEtDwl  =  rho*total_cond*DftDwl*(qvsat - cana_q)
     DEtDwf  =  rho*total_cond*DftDwf*(qvsat - cana_q)
     ! calculate the coefficients of the intercepted liquid evaporation linearization
     Eli0    =  rho*con_v_v*fw*(qvsat - cana_q)
     DEliDTv =  rho*con_v_v*fw*DqvsatDTv
     DEliDqc = -rho*con_v_v*fw
     DEliDwl =  rho*con_v_v*DfwDwl*(qvsat-cana_q)
     DEliDwf =  rho*con_v_v*DfwDwf*(qvsat-cana_q)
     ! calculate the coefficients of the intercepted snow evaporation linearization
     Efi0    =  rho*con_v_v*fs*(qvsat - cana_q)
     DEfiDTv =  rho*con_v_v*fs*DqvsatDTv
     DEfiDqc = -rho*con_v_v*fs
     DEfiDwl =  rho*con_v_v*DfsDwl*(qvsat-cana_q)
     DEfiDwf =  rho*con_v_v*DfsDwf*(qvsat-cana_q)
  else
     ! Flux is directed TOWARD the surface: no transpiration (assuming plants do not
     ! take water through stomata), and condensation does not depend on the fraction
     ! of wet canopy -- dew formation occurs on the entire surface

     ! prohibit transpiration:
     Et0     = 0
     DEtDTv  = 0; DEtDwl = 0; DEtDwf = 0;
     DEtDqc  = 0
     ! calculate dew or frost formation rates, depending on the temperature
     Eli0    = 0; Efi0    = 0
     DEliDTv = 0; DEfiDTv = 0
     DEliDqc = 0; DEfiDqc = 0
     DEliDwl = 0; DEfiDwl = 0
     DEliDwf = 0; DEfiDwf = 0
     ! calculate the coefficients of the intercepted liquid condensation linearization
     if(vegn_T >= tfreeze) then
        Eli0    =  rho*con_v_v*(qvsat - cana_q)
        DEliDTv =  rho*con_v_v*DqvsatDTv
        DEliDqc = -rho*con_v_v
     else
        ! calculate the coefficients of the intercepted snow condensation linearization
        Efi0    =  rho*con_v_v*(qvsat - cana_q)
        DEfiDTv =  rho*con_v_v*DqvsatDTv
        DEfiDqc = -rho*con_v_v
     endif
     ! prohibit switching from condensation to evaporation if the water content
     ! is below certain threshold
     if (vegn_Wl < min_Wl) then
        Eli0 = 0 ; DEliDTv = 0 ; DEliDqc = 0 ; DEliDwl = 0 ; DEliDwf = 0
     endif
     if (vegn_Ws < min_Ws) then
        Efi0 = 0 ; DEfiDTv = 0 ; DEfiDqc = 0 ; DEfiDwl = 0 ; DEfiDwf = 0
     endif

  endif
  ! ---- diagnostic section
  call send_tile_data(id_evap_demand, evap_demand, diag)
  call send_tile_data(id_stomatal, stomatal_cond, diag)
  call send_tile_data(id_an_op, cohort%An_op, diag)
  call send_tile_data(id_an_cl, cohort%An_cl, diag)
  call send_tile_data(id_con_v_h, con_v_h, diag)
  call send_tile_data(id_con_v_v, con_v_v, diag)
  call send_tile_data(id_phot_co2, phot_co2, diag)

end subroutine vegn_step_1


! ============================================================================
! Given the surface solution, substitute it back into the vegetation equations
! to determine new vegetation state.
subroutine vegn_step_2 ( vegn, diag, &
     delta_Tv, delta_wl, delta_wf, &
     vegn_melt, &
     vegn_ovfl_l,  vegn_ovfl_s,  & ! overflow of liquid and solid water from the canopy, kg/(m2 s)
     vegn_ovfl_Hl, vegn_ovfl_Hs  ) ! heat flux carried from canopy by overflow, W/(m2 s)

  ! ---- arguments
  type(vegn_tile_type) , intent(inout) :: vegn
  type(diag_buff_type) , intent(inout) :: diag
  real, intent(in) :: &
       delta_Tv, & ! change in vegetation temperature, degK
       delta_wl, & ! change in intercepted liquid water mass, kg/m2
       delta_wf    ! change in intercepted frozen water mass, kg/m2
  real, intent(out) :: &
       vegn_melt, &
       vegn_ovfl_l,   vegn_ovfl_s,   & ! overflow of liquid and solid water from the canopy
       vegn_ovfl_Hl, vegn_ovfl_Hs      ! heat flux from canopy due to overflow

  ! ---- local variables
  real :: &
     vegn_Wl_max, &  ! max. possible amount of liquid water in the canopy
     vegn_Ws_max, &  ! max. possible amount of solid water in the canopy
     mcv, &
     cap0, melt_per_deg, &
     Wl, Ws  ! positively defined amounts of water and snow on canopy
  type(vegn_cohort_type), pointer :: cohort

  ! get the pointer to the first (and, currently, the only) cohort
  cohort => vegn%cohorts(1)

  if (is_watch_point()) then
     write(*,*)'#### vegn_step_2 input ####'
     __DEBUG3__(delta_Tv, delta_wl, delta_wf)
     __DEBUG1__(cohort%Tv)
  endif

  ! update vegetation state
  cohort%Tv = cohort%Tv + delta_Tv
  cohort%Wl = cohort%Wl + delta_wl
  cohort%Ws = cohort%Ws + delta_wf

  call vegn_data_intrcptn_cap(cohort, vegn_Wl_max, vegn_Ws_max)
  call vegn_data_heat_capacity(cohort, mcv)


  ! ---- update for evaporation and interception -----------------------------
  cap0 = mcv + clw*cohort%Wl + csw*cohort%Ws

  if(is_watch_point()) then
     write (*,*)'#### vegn_step_2 #### 1'
     __DEBUG1__(cap0)
     __DEBUG1__(cohort%Tv)
     __DEBUG2__(cohort%Wl, cohort%Ws)
  endif
  ! melt on the vegetation should probably be prohibited altogether, since
  ! the amount of melt or freeze calculated this way is severely underestimated
  ! (depending on the overall vegetation heat capacity) which leads to extended
  ! periods when the canopy temperature is fixed at freezing point.
  if (lm2) then
     vegn_melt = 0
  else
     ! ---- freeze/melt of intercepted water
     ! heat capacity of leaf + intercepted water/snow _can_ go below zero if the
     ! total water content goes below zero as a result of implicit time step.
     ! If it does, we just prohibit melt, setting it to zero.
     if(cap0 > 0)then
        melt_per_deg = cap0 / hlf
        if (cohort%Ws>0 .and. cohort%Tv>tfreeze) then
           vegn_melt =  min(cohort%Ws, (cohort%Tv-tfreeze)*melt_per_deg)
        else if (cohort%Wl>0 .and. cohort%Tv<tfreeze) then
           vegn_melt = -min(cohort%Wl, (tfreeze-cohort%Tv)*melt_per_deg)
        else
           vegn_melt = 0
        endif
        cohort%Ws = cohort%Ws - vegn_melt
        cohort%Wl = cohort%Wl + vegn_melt
        if (vegn_melt/=0) &
             cohort%Tv = tfreeze + (cap0*(cohort%Tv-tfreeze) - hlf*vegn_melt) &
             / ( cap0 + (clw-csw)*vegn_melt )
        vegn_melt = vegn_melt / delta_time
     else
        vegn_melt = 0
     endif
  endif

  if(is_watch_point()) then
     write (*,*)'#### vegn_step_2 #### 1'
     __DEBUG1__(cap0)
     __DEBUG1__(cohort%Tv)
     __DEBUG3__(vegn_melt, cohort%Wl, cohort%Ws)
  endif

  ! ---- update for overflow -------------------------------------------------
  Wl = max(cohort%Wl,0.0); Ws = max(cohort%Ws,0.0)
  vegn_ovfl_l = max (0.,Wl-vegn_Wl_max)/delta_time
  vegn_ovfl_s = max (0.,Ws-vegn_Ws_max)/delta_time
  vegn_ovfl_Hl = clw*vegn_ovfl_l*(cohort%Tv-tfreeze)
  vegn_ovfl_Hs = csw*vegn_ovfl_s*(cohort%Tv-tfreeze)

  cohort%Wl = cohort%Wl - vegn_ovfl_l*delta_time
  cohort%Ws = cohort%Ws - vegn_ovfl_s*delta_time

  if(is_watch_point()) then
     write(*,*)'#### vegn_step_2 output #####'
     __DEBUG3__(vegn_melt, vegn_ovfl_l, vegn_ovfl_s)
     __DEBUG2__(vegn_ovfl_Hl,vegn_ovfl_Hs)
  endif

  ! ---- diagnostic section
  call send_tile_data(id_temp,   cohort%Tv, diag)
  call send_tile_data(id_wl,     cohort%Wl, diag)
  call send_tile_data(id_ws,     cohort%Ws, diag)

  call send_tile_data(id_height, cohort%height, diag)
  call send_tile_data(id_lai, cohort%lai, diag)
  call send_tile_data(id_lai_var, cohort%lai, diag)
  call send_tile_data(id_lai_std, cohort%lai, diag)
  call send_tile_data(id_sai, cohort%sai, diag)
  call send_tile_data(id_leaf_size, cohort%leaf_size, diag)
  call send_tile_data(id_root_density, cohort%root_density, diag)
  call send_tile_data(id_root_zeta, cohort%root_zeta, diag)
  call send_tile_data(id_rs_min, cohort%rs_min, diag)
  call send_tile_data(id_leaf_refl, cohort%leaf_refl, diag)
  call send_tile_data(id_leaf_tran, cohort%leaf_tran, diag)
  call send_tile_data(id_leaf_emis, cohort%leaf_emis, diag)
  call send_tile_data(id_snow_crit, cohort%snow_crit, diag)
end subroutine vegn_step_2


! ============================================================================
! do the vegetation calculations that require updated (end-of-timestep) values
! of prognostic land variables
subroutine vegn_step_3(vegn, soil, cana_T, precip, vegn_fco2, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  type(soil_tile_type), intent(inout) :: soil
  real, intent(in) :: cana_T ! canopy temperature, deg K
  real, intent(in) :: precip ! total (rain+snow) precipitation, kg/(m2 s)
  real, intent(out) :: vegn_fco2 ! co2 flux from vegetation, kg CO2/(m2 s)
  type(diag_buff_type), intent(inout) :: diag

  ! ---- local vars
  real :: tsoil ! average temperature of soil for soil carbon decomposition, deg K
  real :: theta ! average soil wetness, unitless
  real :: psist ! psi stress index
  real :: depth_ave! depth for averaging soil moisture based on Jackson function for root distribution
  real :: percentile = 0.95

  tsoil = soil_ave_temp (soil,soil_carbon_depth_scale)
  ! depth for 95% of root according to Jackson distribution
  depth_ave = -log(1.-percentile)*vegn%cohorts(1)%root_zeta

  theta = soil_ave_theta1(soil, depth_ave)

  if(is_watch_point()) then
     write(*,*)'#### vegn_step_3 drought input ####'
     __DEBUG3__(depth_ave, tsoil, theta)
  endif

  call vegn_carbon_int(vegn, soil, tsoil, theta, diag)
  ! decrease, if necessary, csmoke spending rate so that csmoke pool
  ! is never depleted below zero
  vegn%csmoke_rate = max( 0.0, &
       min( vegn%csmoke_rate, &
            vegn%csmoke_pool/dt_fast_yr)&
       )
  ! update smoke pool -- stored amount of carbon lost to fire
  vegn%csmoke_pool = vegn%csmoke_pool - &
       vegn%csmoke_rate*dt_fast_yr
  ! decrease harvested rates so that pools are not depleted below zero
  vegn%harv_rate(:) = max( 0.0, &
                           min(vegn%harv_rate(:), vegn%harv_pool(:)/dt_fast_yr) &
                         )
  ! update harvested pools -- amounts of stored harvested carbon by category
  vegn%harv_pool(:) = vegn%harv_pool(:) - &
       vegn%harv_rate(:)*dt_fast_yr
  ! --- calculate total co2 flux from vegetation
  vegn_fco2 = -vegn%nep + vegn%csmoke_rate + sum(vegn%harv_rate(:))
  ! --- convert it to kg CO2/(m2 s)
  vegn_fco2 = vegn_fco2*mol_CO2/(mol_C*seconds_per_year)

  ! --- accumulate values for climatological averages
  vegn%tc_av     = vegn%tc_av + cana_T
  vegn%tsoil_av  = vegn%tsoil_av + tsoil
  vegn%precip_av = vegn%precip_av + precip
  if (xwilt_available) then
     theta = soil_ave_theta1(soil,depth_ave)
  else
     theta = soil_ave_theta0(soil,vegn%cohorts(1)%root_zeta)
  endif
  vegn%theta_av_phen = vegn%theta_av_phen + theta
  vegn%theta_av_fire = vegn%theta_av_fire + soil_ave_theta1(soil,depth_ave)
  psist = soil_psi_stress(soil,vegn%cohorts(1)%root_zeta)
  vegn%psist_av  = vegn%psist_av + psist

  vegn%n_accum   = vegn%n_accum+1

  call send_tile_data(id_theph, theta, diag)
  call send_tile_data(id_psiph, psist, diag)

end subroutine vegn_step_3


! ============================================================================
! update slow components of the vegetation model
subroutine update_vegn_slow( )

  ! ---- local vars ----------------------------------------------------------
  integer :: second, minute, hour, day0, day1, month0, month1, year0, year1
  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile
  integer :: i,j,k,l ! current point indices
  integer :: ii ! pool iterator
  integer :: n ! number of cohorts
  real    :: weight_ncm ! low-pass filter value for the number of cold months
  character(64) :: timestamp

  ! variables for conservation checks
  real :: lmass0, fmass0, heat0, cmass0
  real :: lmass1, fmass1, heat1, cmass1
  character(64) :: tag

  ! get components of calendar dates for this and previous time step
  call get_date(lnd%time,             year0,month0,day0,hour,minute,second)
  call get_date(lnd%time-lnd%dt_slow, year1,month1,day1,hour,minute,second)

  if(month0 /= month1) then
     ! heartbeat
     write(timestamp,'("Current date is ",i4.4,"-",i2.2,"-",i2.2)')year0,month0,day0
     call error_mesg('update_vegn_slow',trim(timestamp),NOTE)
  endif

  ce = first_elmt(land_tile_map, lnd%ls)
  do while (loop_over_tiles(ce,tile, l,k))
     call set_current_point(l,k) ! this is for debug output only
     if(.not.associated(tile%vegn)) cycle ! skip the rest of the loop body

     if(do_check_conservation) then
        ! + conservation check, part 1: calculate the pre-transition totals
        call get_tile_water(tile,lmass0,fmass0)
        heat0  = land_tile_heat  (tile)
        cmass0 = land_tile_carbon(tile)
        ! - end of conservation check, part 1
     endif

     if (day1 /= day0) then
        call vegn_daily_npp(tile%vegn)
     endif

     ! monthly averaging
     if (month1 /= month0) then
        ! compute averages from accumulated monthly values
        tile%vegn%tc_av     = tile%vegn%tc_av     / tile%vegn%n_accum
        tile%vegn%tsoil_av  = tile%vegn%tsoil_av  / tile%vegn%n_accum
        tile%vegn%theta_av_phen  = tile%vegn%theta_av_phen  / tile%vegn%n_accum
        tile%vegn%theta_av_fire  = tile%vegn%theta_av_fire  / tile%vegn%n_accum
	tile%vegn%psist_av  = tile%vegn%psist_av  / tile%vegn%n_accum
        tile%vegn%precip_av = tile%vegn%precip_av / tile%vegn%n_accum
        ! accumulate annual values
        tile%vegn%p_ann_acm = tile%vegn%p_ann_acm+tile%vegn%precip_av
        tile%vegn%t_ann_acm = tile%vegn%t_ann_acm+tile%vegn%tc_av
        if ( tile%vegn%tc_av < cold_month_threshold ) &
             tile%vegn%ncm_acm = tile%vegn%ncm_acm+1
        tile%vegn%t_cold_acm = min(tile%vegn%t_cold_acm, tile%vegn%tc_av)

        tile%vegn%nmn_acm = tile%vegn%nmn_acm+1 ! increase the number of accumulated months
     endif

     ! annual averaging
     if (year1 /= year0) then
        ! The ncm smoothing is coded as a low-pass exponential filter. See, for example
        ! http://en.wikipedia.org/wiki/Low-pass_filter
        weight_ncm = 1/(1+tau_smooth_ncm)
        if(tile%vegn%nmn_acm /= 0) then
           ! calculate annual averages from accumulated values
           tile%vegn%p_ann  = tile%vegn%p_ann_acm/tile%vegn%nmn_acm
           tile%vegn%t_ann  = tile%vegn%t_ann_acm/tile%vegn%nmn_acm
           tile%vegn%t_cold = tile%vegn%t_cold_acm
           tile%vegn%ncm    = weight_ncm*tile%vegn%ncm_acm + (1-weight_ncm)*tile%vegn%ncm
           ! reset accumulated values
           tile%vegn%ncm_acm    = 0
           tile%vegn%p_ann_acm  = 0
           tile%vegn%t_ann_acm  = 0
           tile%vegn%t_cold_acm = HUGE(tile%vegn%t_cold_acm)
        endif
!!$        call calc_miami_npp(tile%vegn)
        tile%vegn%nmn_acm = 0
     endif

     if (year1 /= year0 .and. do_biogeography) then
        call vegn_biogeography(tile%vegn)
     endif

     if (year1 /= year0 .and. do_peat_redistribution) then
        call redistribute_peat_carbon(tile%soil)
     endif

     if (month1 /= month0.and.do_patch_disturbance) then
        call update_fuel(tile%vegn,tile%soil%w_wilt(1)/tile%soil%pars%vwc_sat)
        ! assume that all layers are the same soil type and wilting is vertically homogeneous
     endif

     if (day1 /= day0 .and. do_cohort_dynamics) then
        n = tile%vegn%n_cohorts
        call send_tile_data(id_cgain,sum(tile%vegn%cohorts(1:n)%carbon_gain),tile%diag)
        call send_tile_data(id_closs,sum(tile%vegn%cohorts(1:n)%carbon_loss),tile%diag)
        call send_tile_data(id_wdgain,sum(tile%vegn%cohorts(1:n)%bwood_gain),tile%diag)
        call vegn_growth(tile%vegn)
        call vegn_nat_mortality(tile%vegn,tile%soil,86400.0)
     endif

     if  (month1 /= month0 .and. do_phenology) then
        call vegn_phenology (tile%vegn, tile%soil)
        ! assume that all layers are the same soil type and wilting is vertically homogeneous
     endif

     if (year1 /= year0 .and. do_patch_disturbance) then
        call vegn_disturbance(tile%vegn, tile%soil, seconds_per_year)
     endif

     call vegn_harvesting(tile%vegn, tile%soil, year0/=year1, month0/=month1, day0/=day1)

     if (year1 /= year0) then
        ! update rates of carbon transfer from intermediate pools to soil and litter
        tile%vegn%fsc_rate_ag = tile%vegn%fsc_pool_ag/fsc_pool_spending_time
        tile%vegn%ssc_rate_ag = tile%vegn%ssc_pool_ag/ssc_pool_spending_time
        tile%vegn%fsc_rate_bg = tile%vegn%fsc_pool_bg/fsc_pool_spending_time
        tile%vegn%ssc_rate_bg = tile%vegn%ssc_pool_bg/ssc_pool_spending_time

        tile%vegn%leaflitter_buffer_rate_ag = tile%vegn%leaflitter_buffer_ag/fsc_pool_spending_time
        tile%vegn%coarsewoodlitter_buffer_rate_ag = tile%vegn%coarsewoodlitter_buffer_ag/ssc_pool_spending_time
        ! update rates of harvest spending
        where(harvest_spending_time(:)>0)
           tile%vegn%harv_rate(:) = &
                tile%vegn%harv_pool(:)/harvest_spending_time(:)
        elsewhere
           tile%vegn%harv_rate(:) = 0.0
        end where
     endif

     if (do_check_conservation) then
        ! + conservation check, part 2: calculate totals in final state, and compare
        ! with previous totals
        tag = 'update_vegn_slow'
        call get_tile_water(tile,lmass1,fmass1)
        heat1  = land_tile_heat  (tile)
        cmass1 = land_tile_carbon(tile)
        call check_conservation (tag,'liquid water', lmass0, lmass1, water_cons_tol)
        call check_conservation (tag,'frozen water', fmass0, fmass1, water_cons_tol)
        call check_conservation (tag,'carbon'      , cmass0, cmass1, carbon_cons_tol)
        ! call check_conservation (tag,'heat content', heat0 , heat1 , 1e-16)
        ! - end of conservation check, part 2
     endif

     ! ---- diagnostic section
     call send_tile_data(id_t_ann,   tile%vegn%t_ann,   tile%diag)
     call send_tile_data(id_t_cold,  tile%vegn%t_cold,  tile%diag)
     call send_tile_data(id_lambda,  tile%vegn%lambda,  tile%diag)
     call send_tile_data(id_p_ann,   tile%vegn%p_ann,   tile%diag)
     call send_tile_data(id_ncm,     real(tile%vegn%ncm), tile%diag)
     call send_tile_data(id_afire,   tile%vegn%disturbance_rate(1), tile%diag)
     call send_tile_data(id_atfall,  tile%vegn%disturbance_rate(0), tile%diag)

     do ii = 1,N_HARV_POOLS
        call send_tile_data(id_harv_pool(ii),tile%vegn%harv_pool(ii),tile%diag)
        call send_tile_data(id_harv_rate(ii),tile%vegn%harv_rate(ii),tile%diag)
     enddo
     if (id_t_harv_pool>0) call send_tile_data(id_t_harv_pool,sum(tile%vegn%harv_pool(:)),tile%diag)
     if (id_t_harv_rate>0) call send_tile_data(id_t_harv_rate,sum(tile%vegn%harv_rate(:)),tile%diag)
     call send_tile_data(id_csmoke_pool,tile%vegn%csmoke_pool,tile%diag)
     call send_tile_data(id_csmoke_rate,tile%vegn%csmoke_rate,tile%diag)

     call send_tile_data(id_fsc_pool_ag,tile%vegn%fsc_pool_ag,tile%diag)
     call send_tile_data(id_fsc_rate_ag,tile%vegn%fsc_rate_ag,tile%diag)
     call send_tile_data(id_ssc_pool_ag,tile%vegn%ssc_pool_ag,tile%diag)
     call send_tile_data(id_ssc_rate_ag,tile%vegn%ssc_rate_ag,tile%diag)
     call send_tile_data(id_fsc_pool_bg,tile%vegn%fsc_pool_ag,tile%diag)
     call send_tile_data(id_fsc_rate_bg,tile%vegn%fsc_rate_ag,tile%diag)
     call send_tile_data(id_ssc_pool_bg,tile%vegn%ssc_pool_ag,tile%diag)
     call send_tile_data(id_ssc_rate_bg,tile%vegn%ssc_rate_ag,tile%diag)

     call send_tile_data(id_leaflitter_buffer_ag,tile%vegn%leaflitter_buffer_ag,tile%diag)
     call send_tile_data(id_leaflitter_buffer_rate_ag,tile%vegn%leaflitter_buffer_rate_ag,tile%diag)
     call send_tile_data(id_coarsewoodlitter_buffer_ag,tile%vegn%coarsewoodlitter_buffer_ag,tile%diag)
     call send_tile_data(id_coarsewoodlitter_buffer_rate_ag,tile%vegn%coarsewoodlitter_buffer_rate_ag,tile%diag)

     n=tile%vegn%n_cohorts
     call send_tile_data(id_bl,      sum(tile%vegn%cohorts(1:n)%bl),     tile%diag)
     call send_tile_data(id_blv,     sum(tile%vegn%cohorts(1:n)%blv),    tile%diag)
     call send_tile_data(id_br,      sum(tile%vegn%cohorts(1:n)%br),     tile%diag)
     call send_tile_data(id_bsw,     sum(tile%vegn%cohorts(1:n)%bsw),    tile%diag)
     call send_tile_data(id_bwood,   sum(tile%vegn%cohorts(1:n)%bwood),  tile%diag)
     if (id_btot>0) call send_tile_data(id_btot, &
                sum(tile%vegn%cohorts(1:n)%bl    &
                   +tile%vegn%cohorts(1:n)%blv   &
                   +tile%vegn%cohorts(1:n)%br    &
                   +tile%vegn%cohorts(1:n)%bsw   &
                   +tile%vegn%cohorts(1:n)%bwood ), tile%diag)

     call send_tile_data(id_fuel,    tile%vegn%fuel, tile%diag)
     call send_tile_data(id_species, real(tile%vegn%cohorts(1)%species), tile%diag)
     call send_tile_data(id_status,  real(tile%vegn%cohorts(1)%status),  tile%diag)
     call send_tile_data(id_leaf_age,real(tile%vegn%cohorts(1)%leaf_age),  tile%diag)!ens

     ! carbon budget tracking
     call send_tile_data(id_fsc_in,  sum(tile%soil%fsc_in(:)),  tile%diag)
     call send_tile_data(id_fsc_out, tile%vegn%fsc_out, tile%diag)
     call send_tile_data(id_ssc_in,  sum(tile%soil%ssc_in(:)),  tile%diag)
     call send_tile_data(id_ssc_out, tile%vegn%ssc_out, tile%diag)
     call send_tile_data(id_deadmic_out, tile%vegn%deadmic_out, tile%diag)
     call send_tile_data(id_veg_in,  tile%vegn%veg_in,  tile%diag)
     call send_tile_data(id_veg_out, tile%vegn%veg_out, tile%diag)

     ! CMOR variables
     if (id_cProduct>0) then
        cmass1 = 0.0
        do i = 1, N_HARV_POOLS
           if (i/=HARV_POOL_CLEARED) cmass1 = cmass1 + tile%vegn%harv_pool(i)
        enddo
        call send_tile_data(id_cProduct, cmass1, tile%diag)
     endif
     if (id_cAnt>0) then
        call send_tile_data(id_cAnt, sum(tile%vegn%harv_pool(:)), tile%diag)
     endif
     call send_tile_data(id_fFire, tile%vegn%csmoke_rate/seconds_per_year, tile%diag)
     call send_tile_data(id_fGrazing, tile%vegn%harv_rate(HARV_POOL_PAST)/seconds_per_year, tile%diag)
     call send_tile_data(id_fHarvest, tile%vegn%harv_rate(HARV_POOL_CROP)/seconds_per_year, tile%diag)
     call send_tile_data(id_fLuc, &
         (tile%vegn%harv_rate(HARV_POOL_CLEARED) &
         +tile%vegn%harv_rate(HARV_POOL_WOOD_FAST) &
         +tile%vegn%harv_rate(HARV_POOL_WOOD_MED) &
         +tile%vegn%harv_rate(HARV_POOL_WOOD_SLOW) &
         )/seconds_per_year, tile%diag)
     if (id_cLeaf>0) call send_tile_data(id_cLeaf, sum(tile%vegn%cohorts(1:n)%bl), tile%diag)
     if (id_cWood>0) call send_tile_data(id_cWood, &
         (sum(tile%vegn%cohorts(1:n)%bwood)+sum(tile%vegn%cohorts(1:n)%bsw))*agf_bs, tile%diag)
     if (id_cRoot>0) call send_tile_data(id_cRoot, &
         (sum(tile%vegn%cohorts(1:n)%bwood)+sum(tile%vegn%cohorts(1:n)%bsw))*(1-agf_bs) &
            +sum(tile%vegn%cohorts(1:n)%br), &
         tile%diag)
     if (id_cMisc>0) call send_tile_data(id_cMisc, sum(tile%vegn%cohorts(1:n)%blv), tile%diag)

     ! ---- end of diagnostic section

     ! reset averages and number of steps to 0 before the start of new month
     if (month1 /= month0) then
        tile%vegn%n_accum  = 0
        tile%vegn%tc_av    = 0.
        tile%vegn%tsoil_av = 0.
        tile%vegn%theta_av_phen = 0.
        tile%vegn%theta_av_fire = 0.
        tile%vegn%psist_av = 0.
        tile%vegn%precip_av= 0.
     endif

     !reset fuel and drought months before the start of new year
     if (year1 /= year0) then
        tile%vegn%lambda     = 0
        tile%vegn%fuel       = 0
     endif

     if(soil_carbon_option==SOILC_CORPSE) then
        !Knock soil carbon cohorts down to their maximum number
        call cull_cohorts(tile%soil%leafLitter)
        call cull_cohorts(tile%soil%fineWoodLitter)
        call cull_cohorts(tile%soil%coarseWoodLitter)
        do ii=1,size(tile%soil%soil_C)
              call cull_cohorts(tile%soil%soil_C(ii))
        enddo
     endif
  enddo

  ! seed transport
  if (year1 /= year0 .and. do_seed_transport) then
     call vegn_seed_transport()
  endif

  ! override with static vegetation
  if(day1/=day0) &
       call  read_static_vegn(lnd%time)
end subroutine update_vegn_slow


! ============================================================================
subroutine vegn_seed_transport()

  ! local vars
  type(land_tile_enum_type) :: ce
  type(land_tile_type), pointer :: tile
  integer :: l ! current point indices
  real :: total_seed_supply
  real :: total_seed_demand
  real :: f_supply ! fraction of the supply that gets spent
  real :: f_demand ! fraction of the demand that gets satisfied

  total_seed_supply = 0.0; total_seed_demand = 0.0
  ce = first_elmt(land_tile_map, lnd%ls)
  do while (loop_over_tiles(ce,tile,l))
     if(.not.associated(tile%vegn)) cycle ! skip the rest of the loop body

     total_seed_supply = total_seed_supply + vegn_seed_supply(tile%vegn)*tile%frac*lnd%area(l)
     total_seed_demand = total_seed_demand + vegn_seed_demand(tile%vegn)*tile%frac*lnd%area(l)
  enddo
  ! sum totals globally
  call mpp_sum(total_seed_demand, pelist=lnd%pelist)
  call mpp_sum(total_seed_supply, pelist=lnd%pelist)
  ! if either demand or supply are zeros we don't need (or can't) transport anything
  if (total_seed_demand==0.or.total_seed_supply==0)then
     return
  end if

  ! calculate the fraction of the supply that is going to be used
  f_supply = MIN(total_seed_demand/total_seed_supply, 1.0)
  ! calculate the fraction of the demand that is going to be satisfied
  f_demand = MIN(total_seed_supply/total_seed_demand, 1.0)
  ! note that either f_supply or f_demand is 1; the mass conservation law in the
  ! following calculations is satisfied since
  ! f_demand*total_seed_demand - f_supply*total_seed_supply == 0

  ! redistribute part (or possibly all) of the supply to satisfy part (or possibly all)
  ! of the demand
  ce = first_elmt(land_tile_map)
  do while (loop_over_tiles(ce,tile))
     if(associated(tile%vegn)) call vegn_add_bliving(tile%vegn, &
          f_demand*vegn_seed_demand(tile%vegn)-f_supply*vegn_seed_supply(tile%vegn))
  enddo
end subroutine vegn_seed_transport


! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function vegn_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   vegn_tile_exists = associated(tile%vegn)
end function vegn_tile_exists


! ============================================================================
! cohort accessor functions: given a pointer to cohort, return a pointer to a
! specific member of the cohort structure
#define DEFINE_VEGN_ACCESSOR_0D(xtype,x) subroutine vegn_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%vegn))p=>t%vegn%x;endif;\
end subroutine

#define DEFINE_VEGN_ACCESSOR_1D(xtype,x) subroutine vegn_ ## x ## _ptr(t,i,p);\
type(land_tile_type),pointer::t;integer,intent(in)::i;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%vegn))p=>t%vegn%x(i);endif;\
end subroutine

#define DEFINE_COHORT_ACCESSOR(xtype,x) subroutine cohort_ ## x ## _ptr(c,p);\
type(vegn_cohort_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%x;\
end subroutine

#define DEFINE_COHORT_COMPONENT_ACCESSOR(xtype,component,x) subroutine cohort_ ## x ## _ptr(c,p);\
type(vegn_cohort_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%component%x;\
end subroutine

DEFINE_VEGN_ACCESSOR_0D(integer,landuse)
DEFINE_VEGN_ACCESSOR_0D(real,age)
DEFINE_VEGN_ACCESSOR_0D(real,fsc_pool_ag)
DEFINE_VEGN_ACCESSOR_0D(real,fsc_rate_ag)
DEFINE_VEGN_ACCESSOR_0D(real,ssc_pool_ag)
DEFINE_VEGN_ACCESSOR_0D(real,ssc_rate_ag)
DEFINE_VEGN_ACCESSOR_0D(real,fsc_pool_bg)
DEFINE_VEGN_ACCESSOR_0D(real,fsc_rate_bg)
DEFINE_VEGN_ACCESSOR_0D(real,ssc_pool_bg)
DEFINE_VEGN_ACCESSOR_0D(real,ssc_rate_bg)

DEFINE_VEGN_ACCESSOR_0D(real,leaflitter_buffer_ag)
DEFINE_VEGN_ACCESSOR_0D(real,leaflitter_buffer_rate_ag)
DEFINE_VEGN_ACCESSOR_0D(real,coarsewoodlitter_buffer_ag)
DEFINE_VEGN_ACCESSOR_0D(real,coarsewoodlitter_buffer_rate_ag)

DEFINE_VEGN_ACCESSOR_0D(real,tc_av)
DEFINE_VEGN_ACCESSOR_0D(real,theta_av_phen)
DEFINE_VEGN_ACCESSOR_0D(real,theta_av_fire)
DEFINE_VEGN_ACCESSOR_0D(real,psist_av)
DEFINE_VEGN_ACCESSOR_0D(real,tsoil_av)
DEFINE_VEGN_ACCESSOR_0D(real,precip_av)
DEFINE_VEGN_ACCESSOR_0D(real,fuel)
DEFINE_VEGN_ACCESSOR_0D(real,lambda)
DEFINE_VEGN_ACCESSOR_0D(real,t_ann)
DEFINE_VEGN_ACCESSOR_0D(real,p_ann)
DEFINE_VEGN_ACCESSOR_0D(real,t_cold)
DEFINE_VEGN_ACCESSOR_0D(real,ncm)
DEFINE_VEGN_ACCESSOR_0D(real,t_ann_acm)
DEFINE_VEGN_ACCESSOR_0D(real,p_ann_acm)
DEFINE_VEGN_ACCESSOR_0D(real,t_cold_acm)
DEFINE_VEGN_ACCESSOR_0D(real,ncm_acm)
DEFINE_VEGN_ACCESSOR_0D(real,csmoke_pool)
DEFINE_VEGN_ACCESSOR_0D(real,csmoke_rate)

DEFINE_VEGN_ACCESSOR_1D(real,harv_pool)
DEFINE_VEGN_ACCESSOR_1D(real,harv_rate)

DEFINE_COHORT_ACCESSOR(integer,species)
DEFINE_COHORT_ACCESSOR(real,bl)
DEFINE_COHORT_ACCESSOR(real,br)
DEFINE_COHORT_ACCESSOR(real,blv)
DEFINE_COHORT_ACCESSOR(real,bsw)
DEFINE_COHORT_ACCESSOR(real,bwood)
DEFINE_COHORT_ACCESSOR(real,bliving)
DEFINE_COHORT_ACCESSOR(integer,status)
DEFINE_COHORT_ACCESSOR(real,leaf_age)
DEFINE_COHORT_ACCESSOR(real,npp_previous_day)

DEFINE_COHORT_ACCESSOR(real,tv)
DEFINE_COHORT_ACCESSOR(real,wl)
DEFINE_COHORT_ACCESSOR(real,ws)

DEFINE_COHORT_ACCESSOR(real,height)

end module vegetation_mod
