module moist_processes_mod

!-----------------------------------------------------------------------
!
!         interface module for moisture processes
!         ---------------------------------------
!        1) sets up needed derived-type variables related to
!                             condensation / convection parameterizations
!        2) calls convection_driver to process model convection
!        3) calls lscloud_driver to process large-scale clouds
!        4) outputs combined convective / large scale diagnostics
!
!-----------------------------------------------------------------------

! fms modules
use sat_vapor_pres_mod,    only: compute_qs
use time_manager_mod,      only: time_type
use diag_manager_mod,      only: register_diag_field, send_data, &
                                 get_diag_field_id, DIAG_FIELD_NOT_FOUND
use diag_axis_mod,         only: get_axis_num
use diag_data_mod,         only: CMOR_MISSING_VALUE

use mpp_mod,               only: input_nml_file
use fms_mod,               only: error_mesg, FATAL, NOTE,        &
                                 file_exist, check_nml_error,    &
                                 open_namelist_file, close_file, &
                                 write_version_number, stdout,   &
                                 mpp_pe, mpp_root_pe, stdlog,    &
                                 mpp_clock_id, mpp_clock_begin,  &
                                 mpp_clock_end, CLOCK_MODULE,    &
                                 MPP_CLOCK_SYNC, read_data, write_data
use field_manager_mod,     only: MODEL_ATMOS
use tracer_manager_mod,    only: get_tracer_index,&
                                 get_tracer_names, &
                                 NO_TRACER
use constants_mod,         only: CP_AIR, GRAV, HLV, HLS, HLF, &
                                 TFREEZE, WTMAIR, SECONDS_PER_DAY
! atmos_param modules
use physics_types_mod,    only : physics_control_type,    &
                                 physics_tendency_block_type, &
                                 precip_flux_type, &
                                 phys_mp_exch_type, &
                                 physics_input_block_type
use physics_radiation_exch_mod,       &
                          only : clouds_from_moist_block_type, &
                                 exchange_control_type
use lscloud_driver_mod,   only : lscloud_driver_init, lscloud_driver, &
                                 lscloud_driver_time_vary,  &
                                 lscloud_driver_endts, &
                                 lscloud_driver_end
use convection_driver_mod,only : convection_driver_init, &
                                 convection_driver, &
                                 convection_driver_time_vary,  &
                                 cape_cin_diagnostics, &
                                 convection_driver_endts, &
                                 convection_driver_restart, &
                                 convection_driver_end, &
                                 id_pr_g, id_prc_g, id_prsn_g
use diag_integral_mod,    only : diag_integral_field_init, &
                                 sum_diag_integral_field
use atmos_global_diag_mod, only: register_global_diag_field, &
                                 buffer_global_diag, &
                                 send_global_diag
use vert_diff_driver_mod, only : surf_diff_type
use aerosol_types_mod,    only : aerosol_type
use moist_proc_utils_mod, only : tempavg, column_diag, rh_calc,  &
                                 MP_input_type, MP_nml_type,  &
                                 mp_tendency_type, mp_removal_type, &
                                 mp_removal_control_type, &
                                 mp_conv2ls_type, mp_output_type

! atmos_shared modules
use atmos_tracer_utilities_mod, only : get_cmip_param, get_chem_param
use atmos_dust_mod,       only : atmos_dust_init, dust_tracers,   &
                                 n_dust_tracers, do_dust,   &
                                 atmos_dust_wetdep_flux_set
use atmos_sea_salt_mod,   only : atmos_sea_salt_init, seasalt_tracers,  &
                                 n_seasalt_tracers,do_seasalt
use atmos_cmip_diag_mod,   only: register_cmip_diag_field_2d, &
                                 register_cmip_diag_field_3d, &
                                 send_cmip_data_3d, &
                                 cmip_diag_id_type, &
                                 query_cmip_diag_id


implicit none
private

!-----------------------------------------------------------------------
!-------------------- public data/interfaces ---------------------------

public   moist_processes_init, moist_processes_time_vary, moist_processes,&
         moist_processes_restart, &
         moist_processes_endts, moist_processes_end,   &
         set_cosp_precip_sources, define_cosp_precip_fluxes

!-----------------------------------------------------------------------
!-------------------- private data -------------------------------------

private combined_MP_diagnostics, MP_alloc, MP_dealloc, create_Nml_mp, &
        diag_field_init, height_adjust

!--------------------- version number ----------------------------------

   character(len=128) ::  version = '$Id$'
   character(len=128) :: tagname = '$Name$'
   character(len=5), private :: mod_name = 'moist'
   character(len=7), private :: mod_name_tr = 'tracers'
   logical            :: moist_allocated = .false.
   logical            :: module_is_initialized = .false.
!-------------------- namelist data (private) --------------------------

!---------------- namelist variable definitions ------------------------
!
!   do_unified_clouds =
!              switch to turn on/off a unified (LS + conv) cloud
!                scheme (not yet available)
!                [logical, default: do_unified_clouds=false ]
!   do_lsc   = switch to turn on/off large scale condensation
!                [logical, default: do_lsc=false ]
!   do_mca   = switch to turn on/off moist convective adjustment;
!                [logical, default: do_mca=false ]
!   do_ras   = switch to turn on/off relaxed arakawa shubert
!                [logical, default: do_ras=false ]
!   do_uw_conv = switch to turn on/off Univ of Wash shallow convect scheme
!                [logical, default: do_uw_conv=false ]
!   do_donner_deep = switch to turn on/off donner deep convection scheme
!                [logical, default: do_donner_deep=false ]
!   do_dryadj = switch to turn on/off dry adjustment scheme
!                [logical, default: do_dryadj=false ]
!   do_bm    = switch to turn on/off betts-miller scheme
!                [logical, default: do_bm=false ]
!   do_bmmass  = switch to turn on/off betts-miller massflux scheme
!                [logical, default: do_bmmass=false ]
!   do_bmomp  = switch to turn on/off olivier's version of the betts-miller
!                scheme (with separated boundary layer)
!                [logical, default: do_bmomp=false ]
!   do_simple = switch to turn on alternative definition of specific
!                humidity. When true, specific humidity =
!                (rdgas/rvgas)*esat/pressure
!   do_rh_clouds = switch to turn on/off simple relative humidity cloud
!                scheme
!                [logical, default: do_rh_clouds=false ]
!   pdepth   = boundary layer depth in pascals for determining mean
!                temperature tfreeze (used for snowfall determination)
!                [real, default =150.e2 Pa]
!   limit_conv_cloud_frac =
!                [logical, default: limit_conv_cloud_frac=false]
!   include_donmca_in_cosp =
!                [logical, default: include_donmca_in_cosp = true]
!  <DATA NAME="use_online_aerosol" TYPE="logical"  DEFAULT=".true.">
!   the online aerosol fields should be used for nucleation source
!   (rather than climo values) ?
!  <DATA NAME="use_sub_seasalt" TYPE="logical"  DEFAULT=".false.">
!   only sub-micron seasalt particles should be used for nucleation ?
!  </DATA>
!  <DATA NAME="sea_salt_scale" UNITS="" TYPE="real" DEFAULT="0.1">
!   scaling factor used when offline seasalt aerosol is used for
!   nucleation
!  </DATA>
!  <DATA NAME="om_to_oc" UNITS="" TYPE="real"  DEFAULT="1.67">
!   scaling factor used to convert offline organic matter to organic
!   carbon for nucleation
!  </DATA>

logical :: do_unified_clouds = .false.
logical :: do_lsc = .false.
logical :: do_mca=.false.
logical :: do_ras=.false.
logical :: do_uw_conv=.false.
logical :: do_donner_deep=.false.
logical :: do_dryadj=.false.
logical :: do_bm=.false.
logical :: do_bmmass =.false.
logical :: do_bmomp  =.false.
logical :: do_simple =.false.
logical :: do_rh_clouds=.false.
real    :: pdepth = 150.e2
logical :: limit_conv_cloud_frac = .false.
logical :: include_donmca_in_cosp = .true.
logical :: use_online_aerosol = .true.
logical :: use_sub_seasalt = .false.
real    :: sea_salt_scale = 0.1
real    :: om_to_oc = 1.67
logical :: do_height_adjust = .false.
   logical :: do_diag_clouds=.false.

namelist /moist_processes_nml/ do_unified_clouds, do_lsc, do_mca, do_ras,   &
                  do_uw_conv, do_donner_deep, do_dryadj, do_bm,             &
                  do_bmmass, do_bmomp, do_simple, do_rh_clouds,             &
                  do_diag_clouds,                                           &
                  pdepth, limit_conv_cloud_frac, include_donmca_in_cosp,    &
                  use_online_aerosol, use_sub_seasalt, sea_salt_scale,      &
                  om_to_oc, do_height_adjust


!-------------------- diagnostics fields -------------------------------

integer ::   &
           id_snow_tot, id_tot_cld_amt, id_precip  , id_WVP,  &
           id_AWP, id_gust_conv, id_tot_cloud_area,  id_tot_liq_amt,  &
           id_tot_ice_amt, id_tot_h2o, id_tot_vapor, &
           id_lsc_cloud_area,  id_lsc_liq_amt,  id_lsc_ice_amt,  &
           id_conv_cloud_area, id_conv_liq_amt, id_conv_ice_amt, &
           id_LWP_all_clouds,  id_IWP_all_clouds, id_WP_all_clouds, &
           id_rh,  id_qs, id_rh_cmip, id_enth_moist_col,   &
           id_wat_moist_col

integer :: id_prsn, id_pr, id_prw, id_prrc, id_prra, id_prsnc, &
      id_clt,  id_clwvi, id_clivi

integer :: id_max_enthalpy_imbal, id_max_water_imbal
integer :: id_wetdep_om, id_wetdep_SOA, id_wetdep_bc, &
           id_wetdep_so4, id_wetdep_so2, id_wetdep_DMS, &
           id_wetdep_NH4NO3, id_wetdep_seasalt, id_wetdep_dust

type(cmip_diag_id_type) :: ID_cl, ID_clw, ID_cli, ID_hur

integer, dimension(:), allocatable :: id_wetdep
integer, dimension(:), allocatable :: id_wetdep_uw, id_wetdep_donner, &
                                      id_wetdepc_donner, id_wetdepm_donner  !f1p
integer, dimension(:), allocatable :: id_wetdep_kg_m2_s
real, dimension(:), allocatable    :: conv_wetdep, conv_wetdep_kg_m2_s, nb_N_ox, nb_N_red, nb_N

real :: missing_value = -999.

! cmip names, long_names, standard names for wetdep diag fields
integer :: id_wetpoa_cmip, id_wetsoa_cmip, id_wetbc_cmip, id_wetdust_cmip, &
           id_wetss_cmip, id_wetso4_cmip, id_wetso2_cmip, id_wetdms_cmip, id_wetnh4_cmip
character(len=8), dimension(9) :: cmip_names = (/"poa ","soa ","bc  ","dust","ss  ","so4 ","so2 ","dms ","nh4 "/)
character(len=64), dimension(9) :: cmip_longnames = &
                                  (/"Dry Aerosol Primary Organic Matter  ", &
                                    "Dry Aerosol Secondary Organic Matter", &
                                    "Black Carbon Aerosol Mass           ", &
                                    "Dust                                ", &
                                    "Seasalt                             ", &
                                    "SO4                                 ", &
                                    "SO2                                 ", &
                                    "DMS                                 ", &
                                    "NH4+NH3                             "/)
character(len=64), dimension(9) :: cmip_stdnames = &
                                  (/"primary_particulate_organic_matter_dry_aerosol  ", &
                                    "secondary_particulate_organic_matter_dry_aerosol", &
                                    "elemental_carbon_dry_aerosol                    ", &
                                    "dust_dry_aerosol                                ", &
                                    "seasalt_dry_aerosol                             ", &
                                    "sulfate_dry_aerosol                             ", &
                                    "sulfur_dioxide                                  ", &
                                    "dimethyl_sulfide                                ", &
                                    "ammonium_dry_aerosol                            "/)

!-------------------- individual scheme tracers ------------------------

integer :: nbcphobic =0
integer :: nbcphilic =0
integer :: nomphobic =0
integer :: nomphilic =0
integer :: nDMS      =0
integer :: nSO2      =0
integer :: nSO4      =0
integer :: nSOA      =0
integer :: nNH4NO3   =0
integer :: nNH4      =0
integer :: nH2O2     =0


!------------------- other global variables and parameters -------------

real    :: strat_precip_in_cosp = 0.
real    :: donner_precip_in_cosp = 0.
real    :: uw_precip_in_cosp = 0.


!------------------ allocatable moist processes variables --------------

real, allocatable, dimension(:,:)   :: max_enthalpy_imbal,   &
                                       max_water_imbal

real, allocatable, dimension(:,:)   :: prec_intgl

!-----------------------------------------------------------------------

type(MP_nml_type)      :: Nml_mp
!ype(MP_removal_type)  :: Removal_mp
type(MP_removal_control_type)  :: Removal_mp_control

logical :: wetdep_diagnostics_desired = .false.

!-----------------------------------------------------------------------
!   variables extracted from control variables during _init
!-----------------------------------------------------------------------
logical :: donner_meso_is_largescale
logical :: doing_prog_clouds
integer :: do_clubb
logical :: do_cosp
logical :: use_tau
integer :: nsphum, nql, nqi, nqa, nqn, nqni, nqr, nqs, nqg
integer :: num_prog_tracers



                             contains



!########################################################################

subroutine moist_processes_init ( id, jd, kd, lonb, latb, &
                                 lon, lat, phalf, pref,  axes, Time, &
                                 Physics_control, Exch_ctrl)

!-----------------------------------------------------------------------
integer,                  intent(in)    :: id, jd, kd, axes(4)
real, dimension(:,:),     intent(in)    :: lonb, latb
real,dimension(:,:),      intent(in)    :: lon,  lat    ! h1g
real,dimension(:,:,:),    intent(in)    :: phalf        ! h1g
real, dimension(:),       intent(in)    :: pref
type(time_type),          intent(in)    :: Time
type (physics_control_type), intent(in) :: Physics_control
type (exchange_control_type), intent(inout) :: Exch_ctrl
!-----------------------------------------------------------------------
!
!      input
!     --------
!
!      id, jd        number of horizontal grid points in the global
!                    fields along the x and y axis, repectively.
!                      [integer]
!
!      kd            number of vertical points in a column of atmosphere
!      lonb
!      latb
!      lon
!      lat
!      phalf
!      pref
!      axes
!      Time
!      Physics_control
!      Exch_ctrl
!-----------------------------------------------------------------------

!-------------------------
!   local variables:

      integer :: unit, io, ierr, logunit
      integer :: n, k

!-----------------------------------------------------------------------

      if ( module_is_initialized ) return

!-----------------------------------------------------------------------
!   define some variables needed in moist processes that come from the
!   exchange_control_type variable. these variables are needed in both
!   radiation and physics modules.
!-----------------------------------------------------------------------
      do_clubb = Exch_ctrl%do_clubb
      do_cosp = Exch_ctrl%do_cosp
      doing_prog_clouds = Exch_ctrl%doing_prog_clouds
      donner_meso_is_largescale = Exch_ctrl%donner_meso_is_largescale

!---------------------------------------------------------------------
!   define  variables needed here that are needed in other physics routines
!   and are passed in via a physics_control_type variable.
!---------------------------------------------------------------------
      nsphum = Physics_control%nsphum
      nql    = Physics_control%nql
      nqi    = Physics_control%nqi
      nqa    = Physics_control%nqa
      nqn    = Physics_control%nqn
      nqni   = Physics_control%nqni
      nqr    = Physics_control%nqr
      nqs    = Physics_control%nqs
      nqg    = Physics_control%nqg

      num_prog_tracers = Physics_control%num_prog_tracers
      use_tau = Physics_control%use_tau

!-----------------------------------------------------------------------
!   process the moist_processes_nml.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=moist_processes_nml, iostat=io)
        ierr = check_nml_error(io,'moist_processes_nml')
#else

        unit = open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
          read  (unit, nml=moist_processes_nml, iostat=io, end=10)
          ierr = check_nml_error(io,'moist_processes_nml')
        enddo
  10    call close_file (unit)
#endif

!--------- write version and namelist to standard log ------------

        call write_version_number ( version, tagname )
        logunit = stdlog()
        if ( mpp_pe() == mpp_root_pe() ) &
                         write ( logunit, nml=moist_processes_nml )
      endif

!------------------------------------------------------------------------
!   define the number of cloud schemes active in the model. Define logicals
!   indicating the status of each available cloud scheme in the current
!   model configuration. (prognostic LS cloud scheme already defined
!   in physics_driver based on presence or absence of cloud tracers.)
!------------------------------------------------------------------------
      Exch_ctrl%ncld = 0

      if (doing_prog_clouds) then
        Exch_ctrl%ncld = Exch_ctrl%ncld + 1
      endif
      if (do_donner_deep) then
        Exch_ctrl%ncld = Exch_ctrl%ncld + 2
        Exch_ctrl%doing_donner =  .true.
      endif
      if (do_uw_conv) then
        Exch_ctrl%ncld = Exch_ctrl%ncld + 1
        Exch_ctrl%doing_uw_conv =  .true.
      endif

!----------------------------------------------------------------------
!   create an mp_nml_type variable (Nml_mp) so that moist_processes_nml
!   variables may be made available to other related modules as needed,
!   obviating the need for the occurrence of the same variable in mutiple
!   namelists.
!----------------------------------------------------------------------
      call create_Nml_mp

!-----------------------------------------------------------------------
!   consistency checks for thes namelist variables
!-----------------------------------------------------------------------
      if ( (do_rh_clouds) .and. doing_prog_clouds ) &
        call error_mesg ('moist_processes_init', &
       'rh_clouds cannot be active when prognostic clouds are', FATAL)

      if (do_donner_deep .and. do_rh_clouds) &
           call error_mesg ('moist_processes_init',  &
            'Cannot currently activate donner_deep_mod with rh_clouds', &
                                                                   FATAL)

!-------------------------------------------------------------------------
!   initialize quantities for global precip field
!-------------------------------------------------------------------------
      call diag_integral_field_init ('prec', 'f6.3')
      allocate (prec_intgl(id,jd))


!----------------------------------------------------------------------
!    define indices for the various potentially-available tracers.
!----------------------------------------------------------------------
      nbcphobic = get_tracer_index(MODEL_ATMOS,'bcphob')
      nbcphilic = get_tracer_index(MODEL_ATMOS,'bcphil')
      nomphobic = get_tracer_index(MODEL_ATMOS,'omphob')
      nomphilic = get_tracer_index(MODEL_ATMOS,'omphil')
      call atmos_dust_init (lonb, latb, axes, Time )
      call atmos_sea_salt_init (lonb, latb, axes, Time )

      nDMS      = get_tracer_index(MODEL_ATMOS,'simpleDMS')
      if (nDMS == NO_TRACER) then
        nDMS    = get_tracer_index(MODEL_ATMOS,'dms')
      endif

      nH2O2      = get_tracer_index(MODEL_ATMOS,'simpleH2O2')
      if ( nH2O2 == NO_TRACER ) THEN
        nH2O2     = get_tracer_index(MODEL_ATMOS,'H2O2')
      end if

      nSO2      = get_tracer_index(MODEL_ATMOS,'simpleSO2')
      if (nSO2 == NO_TRACER) then
        nSO2    = get_tracer_index(MODEL_ATMOS,'so2')
      endif

      nSO4      = get_tracer_index(MODEL_ATMOS,'simpleSO4')
      if (nSO4 == NO_TRACER) then
        nSO4    = get_tracer_index(MODEL_ATMOS,'so4')
      endif

      nSOA      = get_tracer_index(MODEL_ATMOS,'SOA')
      nNH4NO3   = get_tracer_index(MODEL_ATMOS,'nh4no3')
      nNH4      = get_tracer_index(MODEL_ATMOS,'nh4')


!---------------------------------------------------------------------
!    allocate and initialize arrays to hold maximum enthalpy and water
!    imbalances in each column over the course of the current job
!    submission.
!---------------------------------------------------------------------
      allocate (max_enthalpy_imbal (id, jd))
      allocate (max_water_imbal (id, jd))
      max_enthalpy_imbal = 0.
      max_water_imbal = 0.

!------------------------------------------------------------------------
!    allocate and initialize an Mp_removal_control_type variable which will
!    indicate for each tracer whether it is to be transported by the
!    various available convective schemes. Also included is a counter of
!    the number of tracers being affected by each available convective
!    scheme.
!------------------------------------------------------------------------
!     allocate (Removal_mp%control%tracers_in_donner(num_prog_tracers))
!     allocate (Removal_mp%control%tracers_in_ras(num_prog_tracers))
!     allocate (Removal_mp%control%tracers_in_uw(num_prog_tracers))
!     allocate (Removal_mp%control%tracers_in_mca(num_prog_tracers))
!     Removal_mp%control%tracers_in_donner = .false.
!     Removal_mp%control%tracers_in_ras    = .false.
!     Removal_mp%control%tracers_in_uw     = .false.
!     Removal_mp%control%tracers_in_mca    = .false.
!     Removal_mp%control%num_mca_tracers   = 0
!     Removal_mp%control%num_ras_tracers   = 0
!     Removal_mp%control%num_donner_tracers   = 0
!     Removal_mp%control%num_uw_tracers   = 0
      allocate (Removal_mp_control%tracers_in_donner(num_prog_tracers))
      allocate (Removal_mp_control%tracers_in_ras(num_prog_tracers))
      allocate (Removal_mp_control%tracers_in_uw(num_prog_tracers))
      allocate (Removal_mp_control%tracers_in_mca(num_prog_tracers))
      Removal_mp_control%tracers_in_donner = .false.
      Removal_mp_control%tracers_in_ras    = .false.
      Removal_mp_control%tracers_in_uw     = .false.
      Removal_mp_control%tracers_in_mca    = .false.
      Removal_mp_control%num_mca_tracers   = 0
      Removal_mp_control%num_ras_tracers   = 0
      Removal_mp_control%num_donner_tracers   = 0
      Removal_mp_control%num_uw_tracers   = 0

!-----------------------------------------------------------------------
!    call convection_driver_init to initialize the convection scheme(s).
!-----------------------------------------------------------------------
      call convection_driver_init (id, jd, kd, axes, Time, &
                          Physics_control, Exch_ctrl, Nml_mp,   &
!                         Removal_mp%control, lonb, latb, pref)
                          Removal_mp_control, lonb, latb, pref)

!-----------------------------------------------------------------------
!    call lscloud_driver_init to initialize the large-scale cloud scheme.
!-----------------------------------------------------------------------
      call lscloud_driver_init (id,jd,kd, axes, Time, Exch_ctrl, Nml_mp, &
                                 Physics_control, lon, lat, phalf, pref)

!-----------------------------------------------------------------------
!   initialize quantities for diagnostics output
!-----------------------------------------------------------------------
      call diag_field_init ( axes, Time )

!-----------------------------------------------------------------------
!   mark the module as initialized.
!-----------------------------------------------------------------------
      module_is_initialized = .true.

!-----------------------------------------------------------------------


end subroutine moist_processes_init

!#####################################################################

subroutine moist_processes_time_vary (dt)

real, intent(in) :: dt

!-----------------------------------------------------------------------

      call convection_driver_time_vary (dt)
      call lscloud_driver_time_vary (dt)

end subroutine moist_processes_time_vary



!#######################################################################

subroutine moist_processes ( is, ie, js, je, npz, Time, dt, land, ustar,  &
                             bstar, qstar, area, lon, lat,   &
                             Physics_input_block, Moist_clouds_block,   &
                             Physics_tendency_block, Phys_mp_exch,  &
                             Surf_diff, Removal_mp, shflx, lhflx,  &
                             lprec, fprec, gust_cv,  Aerosol)

!-----------------------------------------------------------------------
!
!    in:  is,ie      starting and ending i indices for window
!
!         js,je      starting and ending j indices for window
!
!         npz        number of model layers
!
!         Time       time used for diagnostics [time_type]
!
!         dt         time step (from t(n-1) to t(n+1) if leapfrog)
!                    in seconds   [real]
!
!         land       fraction of surface covered by land
!                      [real, dimension(nlon,nlat)]
!
!         ustar
!         bstar
!         qstar
!         area       grid box area (in m2)
!                      [real, dimension(nlon,nlat)]
!
!         lon        longitude in radians           ! h1g
!                      [real, dimension(nlon,nlat)] ! h1g
!
!         lat        latitude in radians
!                      [real, dimension(nlon,nlat)]
!
!         Physics_input_block
!
! inout:  Moist_clouds_block
!         Physics_tendency_block
!         Phys_mp_exch
!
!
!   out:  lprec      liquid precipitiaton rate (rain) in kg/m2/s
!                      [real, dimension(nlon,nlat)]
!
!         fprec      frozen precipitation rate (snow) in kg/m2/s
!                      [real, dimension(nlon,nlat)]
!
!         gust_cv    gustiness from convection  in m/s
!                      [real, dimension(nlon,nlat)]
!
!       optional
!  -----------------
!
!    in:  Aerosol
!
!-----------------------------------------------------------------------
integer,         intent(in)             :: is,ie,js,je, npz
type(time_type), intent(in)             :: Time
real, intent(in)                        :: dt
real, intent(in) , dimension(:,:)       :: land,        ustar, bstar, qstar
real, intent(in) , dimension(:,:)       :: area, lon, lat
type(physics_input_block_type),    &
                             intent(in) :: Physics_input_block

type(clouds_from_moist_block_type),   &
                          intent(inout) :: Moist_clouds_block
type(physics_tendency_block_type),    &
                          intent(inout) :: Physics_tendency_block
type(phys_mp_exch_type),  intent(inout) :: Phys_mp_exch
type(surf_diff_type),     intent(in)    :: Surf_diff
type(mp_removal_type),    intent(inout) :: Removal_mp
real, dimension(:,:),     intent(in)    :: shflx, lhflx

real, intent(out), dimension(:,:)       :: lprec, fprec, gust_cv
type(aerosol_type),intent(in), optional :: Aerosol

!-----------------------------------------------------------------------
!   local variables:

      type(MP_input_type)    :: Input_mp
      type(MP_output_type)   :: Output_mp
      type(MP_tendency_type) :: Tend_mp
      type(MP_conv2ls_type)  :: C2ls_mp

      real, dimension(ie-is+1, je-js+1, npz)  :: tdt_init
      real, dimension(ie-is+1, je-js+1, npz)  :: tdt_dif, qdt_dif  !miz
      real, dimension(ie-is+1, je-js+1, npz, num_prog_tracers) :: qdt_init
      integer                :: i, j, k, tr
      integer                :: istrat

!---------------------------------------------------------------------
!    verify that the module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('moist_processes_mod',  &
                 'moist_processes_init has not been called.', FATAL)
      endif

!------------------------------------------------------------------------
!   define index to large-scale cloud data in moist_cloud_block_type
!   variable.
!------------------------------------------------------------------------
      istrat = Moist_clouds_block%index_strat

!------------------------------------------------------------------------
!   save values of tdt and qdt upon entry to moist_processes so that
!   the total tendency due to moist_processes may be calculated.
!------------------------------------------------------------------------
      tdt_init = Physics_tendency_block%t_dt
      qdt_init = Physics_tendency_block%q_dt
      tdt_dif  = Physics_tendency_block%t_dt    !miz
      qdt_dif  = Physics_tendency_block%q_dt(:,:,:,nsphum) +   &    !miz
                                Physics_tendency_block%q_dt(:,:,:,nql)  + &
                                   Physics_tendency_block%q_dt(:,:,:,nqi)

!-------------------------------------------------------------------------
!    call MP_alloc to allocate and initialize  (or associate) elements of
!    the derived type variables used in this subroutine.
!-------------------------------------------------------------------------
      call MP_alloc (Physics_input_block, Physics_tendency_block, &
                     Phys_mp_exch, dt, area, lon, lat, land, ustar,  &
                     bstar, qstar, Input_mp, Tend_mp, C2ls_mp, Output_mp,&
                     Removal_mp)

!----------------------------------------------------------------------
!    call routines to process the model clouds. If a unified cloud scheme
!    exists, call its driver; otherwise call separate drivers for the
!    convective and large-scale cloud schemes.
!----------------------------------------------------------------------
      if (Nml_mp%do_unified_clouds) then
        call error_mesg ('moist_processes/moist_processes.F90',   &
                               'unified clouds not yet available', FATAL)
      else
        call convection_driver    &
                    (is, ie, js, je,  Time, dt, Input_mp,   &
                     Tend_mp, C2ls_mp, Output_mp, Removal_mp, &
                     Removal_mp_control,  &
                     Surf_diff, Phys_mp_exch, shflx, lhflx, tdt_dif,  &
                     qdt_dif, Moist_clouds_block, Aerosol=Aerosol)

        call lscloud_driver    &
                    (is, ie, js, je, Time, dt, Input_mp, &
                     Physics_tendency_block%qdiag, Tend_mp, C2ls_mp, &
                     Output_mp, Removal_mp,    &
                     Moist_clouds_block%cloud_data(istrat), &
                     Aerosol = Aerosol)
      endif

!------------------------------------------------------------------------
!   call combined_MP_diagnostics to generate diagnostics that contain
!   combined convective and large-scale cloud information and that define
!   conditions after moist processes.
!------------------------------------------------------------------------
      call combined_MP_diagnostics   &
              (is, ie, js, je, Time, tdt_init, qdt_init, Input_mp,   &
               Moist_clouds_block, Output_mp, Removal_mp)

!------------------------------------------------------------------------
!    define needed output arguments. redefine r to be the value after
!    modification in moist_processes (only by clubb ? -- need to check why
!    this is) and pass it back to physics_driver.
!------------------------------------------------------------------------
      Phys_mp_exch%convect = Output_mp%convect
      lprec   = Output_mp%lprec
      fprec   = Output_mp%fprec
      gust_cv = Output_mp%gust_cv
      if (do_clubb > 0) then
        Phys_mp_exch%diff_t_clubb = Output_mp%diff_t_clubb
      endif
      Phys_mp_exch%diff_cu_mo = Output_mp%diff_cu_mo

!---------------------------------------------------------------------
!    deallocate the derived type variables resident in moist_processes_mod.
!---------------------------------------------------------------------
      call MP_dealloc (Input_mp, Tend_mp, C2ls_mp, Output_mp, Removal_mp)

!-----------------------------------------------------------------------


end subroutine moist_processes

!#####################################################################

subroutine moist_processes_endts

logical :: used

!----------------------------------------------------------------------
!    call convection_driver_endts to complete calcs on this timestep.
!----------------------------------------------------------------------
      call convection_driver_endts

!----------------------------------------------------------------------
!    call lscloud_driver_endts to complete calcs on this timestep.
!----------------------------------------------------------------------
      call lscloud_driver_endts

!----------------------------------------------------------------------
!    output the global precip integral. reset the sum to 0. for use on
!    next step.
!----------------------------------------------------------------------
      call sum_diag_integral_field ('prec', prec_intgl)
      prec_intgl = 0.0
      if (id_pr_g   > 0) used = send_global_diag (id_pr_g)
      if (id_prc_g  > 0) used = send_global_diag (id_prc_g)
      if (id_prsn_g > 0) used = send_global_diag (id_prsn_g)

end subroutine moist_processes_endts


!#######################################################################

subroutine moist_processes_end ( )


      if( .not.module_is_initialized ) return

!--------------------------------------------------------------------
!    call lscloud_driver_end to complete operations in that module.
!--------------------------------------------------------------------
      call lscloud_driver_end

!--------------------------------------------------------------------
!    call convection_driver_end to complete operations in that module.
!--------------------------------------------------------------------
      call convection_driver_end

!--------------------------------------------------------------------
!    deallocate module variables.
!--------------------------------------------------------------------
      deallocate (max_water_imbal)
      deallocate (max_enthalpy_imbal)

      deallocate (Removal_mp_control%tracers_in_donner )   !---> h1g, 2017-02-02
      deallocate (Removal_mp_control%tracers_in_ras )      !---> h1g, 2017-02-02
      deallocate (Removal_mp_control%tracers_in_uw  )      !---> h1g, 2017-02-02
      deallocate (Removal_mp_control%tracers_in_mca )      !---> h1g, 2017-02-02

      deallocate (prec_intgl)                              !---> h1g, 2017-02-02

      if (allocated(id_wetdep_kg_m2_s))   deallocate(id_wetdep_kg_m2_s)
      if (allocated(conv_wetdep_kg_m2_s)) deallocate(conv_wetdep_kg_m2_s)
      if (allocated(conv_wetdep))         deallocate(conv_wetdep)
!--------------------------------------------------------------------

      module_is_initialized = .false.

!-----------------------------------------------------------------------


end subroutine moist_processes_end

!#######################################################################

subroutine set_cosp_precip_sources (cosp_precip_sources)

character(len=16),        intent(in) :: cosp_precip_sources

!------------------------------------------------------------------------
!    define which cloud-producing schemes will have their clouds seen by
!    the COSP simulator.
!------------------------------------------------------------------------
      if (trim(cosp_precip_sources)  == 'stratdeepuw') then
        strat_precip_in_cosp = 1.
        donner_precip_in_cosp = 1.
        uw_precip_in_cosp = 1.
      else if (trim(cosp_precip_sources)  == 'stratdeep') then
        strat_precip_in_cosp = 1.
        donner_precip_in_cosp = 1.
      else if (trim(cosp_precip_sources)  == 'stratuw') then
        strat_precip_in_cosp = 1.
        uw_precip_in_cosp = 1.
      else if (trim(cosp_precip_sources)  == 'deepuw') then
        donner_precip_in_cosp = 1.
        uw_precip_in_cosp = 1.
      else if (trim(cosp_precip_sources)  == 'strat') then
        strat_precip_in_cosp = 1.
      else if (trim(cosp_precip_sources)  == 'deep') then
        donner_precip_in_cosp = 1.
      else if (trim(cosp_precip_sources)  == 'uw') then
        uw_precip_in_cosp = 1.
      else if (trim(cosp_precip_sources)  == 'noprecip') then
!     COSP will be run without any precip input (default settings)
      else
        call error_mesg ('moist_processes_mod:set_cosp_precip_sources', &
        'cosp_precip_sources does not match any currently allowed string',&
                                                                 FATAL)
      endif

end subroutine set_cosp_precip_sources

!#######################################################################


subroutine define_cosp_precip_fluxes (is, js, Precip_flux, Removal_mp)

integer,                   intent(in)    :: is, js
type(precip_flux_type),    intent(inout) :: Precip_flux
type(mp_removal_type),     intent(inout) :: Removal_mp

      integer :: i, j ,k
      integer :: ii, jj

!----------------------------------------------------------------------
!    define the grid-box precip flux as the average of the interface
!    fluxes.
!----------------------------------------------------------------------
      do k=1, size(Removal_mp%frz_meso,3)
        do j=1, size(Removal_mp%frz_meso,2)
          do i=1, size(Removal_mp%frz_meso,1)
            ii = i+is-1
            jj = j+js-1
            if (             donner_meso_is_largescale) then
              Precip_flux%fl_lsrain(ii,jj,k) = 0.5* &
                   ((strat_precip_in_cosp*  &
                                  (Removal_mp%rain3d(i,j,k) +    &
                                        Removal_mp%rain3d(i,j,k+1)) + &
                     donner_precip_in_cosp*  &
                                 (Removal_mp%liq_mesoh(i,j,k) +    &
                                          Removal_mp%liq_mesoh(i,j,k+1))))
              Precip_flux%fl_lssnow(ii,jj,k) = 0.5*   &
                   ((strat_precip_in_cosp*   &
                                   (Removal_mp%snowclr3d(i,j,k) +    &
                                         Removal_mp%snowclr3d(i,j,k+1)) + &
                     donner_precip_in_cosp*  &
                                    (Removal_mp%frz_mesoh(i,j,k) +  &
                                         Removal_mp%frz_mesoh(i,j,k+1))))
              Precip_flux%fl_ccrain(ii,jj,k) =  0.5*  &
                   ((donner_precip_in_cosp* &
                                    (Removal_mp%liq_cellh(i,j,k) +  &
                                        Removal_mp%liq_cellh(i,j,k+1)) +  &
                     uw_precip_in_cosp*   &
                                (Removal_mp%liq_precflxh(i,j,k) +  &
                                       Removal_mp%liq_precflxh(i,j,k+1))))
              Precip_flux%fl_ccsnow(ii,jj,k) =  0.5*  &
                   ((donner_precip_in_cosp*  &
                                 (Removal_mp%frz_cellh(i,j,k) +  &
                                     Removal_mp%frz_cellh(i,j,k+1))  +  &
                     uw_precip_in_cosp*  &
                                 (Removal_mp%ice_precflxh(i,j,k) +   &
                                      Removal_mp%ice_precflxh(i,j,k+1))))
            else
              Precip_flux%fl_lsrain(ii,jj,k) =  0.5*   &
                     strat_precip_in_cosp*   &
                                   (Removal_mp%rain3d(i,j,k) +     &
                                               Removal_mp%rain3d(i,j,k+1))
              Precip_flux%fl_lssnow(ii,jj,k) =  0.5*  &
                     strat_precip_in_cosp*    &
                                  (Removal_mp%snowclr3d(i,j,k) +    &
                                            Removal_mp%snowclr3d(i,j,k+1))
              Precip_flux%fl_ccrain(ii,jj,k) =    0.5* &
                   ((donner_precip_in_cosp*  &
                              (Removal_mp%liq_mesoh(i,j,k) +    &
                                   Removal_mp%liq_mesoh(i,j,k+1) +  &
                                       Removal_mp%liq_cellh(i,j,k) +    &
                                        Removal_mp%liq_cellh(i,j,k+1)) +  &
                     uw_precip_in_cosp*   &
                             (Removal_mp%liq_precflxh(i,j,k) +    &
                                      Removal_mp%liq_precflxh(i,j,k+1))))
               Precip_flux%fl_ccsnow(ii,jj,k) =  0.5*  &
                   ((donner_precip_in_cosp*  &
                              (Removal_mp%frz_mesoh(i,j,k) +    &
                                  Removal_mp%frz_mesoh(i,j,k+1) +  &
                                     Removal_mp%frz_cellh(i,j,k) +  &
                                        Removal_mp%frz_cellh(i,j,k+1)) +  &
                     uw_precip_in_cosp*  &
                               (Removal_mp%ice_precflxh(i,j,k) +    &
                                      Removal_mp%ice_precflxh(i,j,k+1))))
            endif

            if (include_donmca_in_cosp .and. &
                              donner_precip_in_cosp .eq. 1.0) then
               Precip_flux%fl_donmca_snow(ii,jj,k) =   0.5*  &
                            (Removal_mp%mca_frzh(i,j,k) +   &
                                           Removal_mp%mca_frzh(i,j,k+1))
               Precip_flux%fl_donmca_rain(ii,jj,k) =   0.5*  &
                            (Removal_mp%mca_liqh(i,j,k) +   &
                                           Removal_mp%mca_liqh(i,j,k+1))
            endif
          end do
        end do
      end do

!-------------------------------------------------------------------------
!    deallocate the components of Removal_mp% that are used here and no
!    longer needed.
!-------------------------------------------------------------------------
      deallocate (Removal_mp%ice_precflxh)
      deallocate (Removal_mp%liq_precflxh)
      deallocate (Removal_mp%frz_mesoh   )
      deallocate (Removal_mp%liq_mesoh   )
      deallocate (Removal_mp%frz_cellh   )
      deallocate (Removal_mp%liq_cellh   )
      deallocate (Removal_mp%mca_frzh    )
      deallocate (Removal_mp%mca_liqh    )
      deallocate (Removal_mp%rain3d      )
      deallocate (Removal_mp%snowclr3d      )

!---------------------------------------------------------------------

end subroutine define_cosp_precip_fluxes


!#######################################################################

subroutine combined_MP_diagnostics    &
        (is, ie, js, je, Time, tdt_init, qdt_init,    &
         Input_mp, Moist_clouds_block, Output_mp, Removal_mp)

integer,                   intent(in)    :: is, ie, js, je
type(time_type),           intent(in)    :: Time
real, dimension(:,:,:,:),  intent(in   ) :: qdt_init
real, dimension(:,:,:  ),  intent(in   ) :: tdt_init
type(mp_input_type),       intent(in)    :: Input_mp

type(clouds_from_moist_block_type),          &
                           intent(inout) :: Moist_clouds_block
type(mp_output_type),      intent(inout) :: Output_mp
type(mp_removal_type),     intent(inout) :: Removal_mp


!----------------------------------------------------------------------
!   local variables:

      logical, dimension   &
            (size(Output_mp%rdt,1),size(Output_mp%rdt,2)) :: ltemp
      real, dimension  &
            (size(Output_mp%rdt,1),size(Output_mp%rdt,2) ) ::  &
                                     precip, temp_2d, tca2
      real, dimension   &
            (size(Output_mp%rdt,1),size(Output_mp%rdt,2))  :: &
                              total_wetdep_dust, total_wetdep_seasalt
      real, dimension   &
            (size(Output_mp%rdt,1),size(Output_mp%rdt,2),  &
                                              size(Output_mp%rdt,4)) ::  &
                       total_wetdep, total_wetdep_uw, total_wetdep_donner,&
                       total_wetdepc_donner, total_wetdepm_donner
      real, dimension  &
            (size(Output_mp%rdt,1),size(Output_mp%rdt,2),   &
                                             size(Output_mp%rdt,3) ) ::  &
                      total_conv_cloud, conv_cld_frac, tot_conv_liq,   &
                      temp_3d1, temp_3d2, temp_3d3, total_cloud_area, &
                      tot_conv_ice, RH, qsat

      integer :: m, mm, n, k, kx, nbin_dust, nbin_seasalt
      logical :: used
      real    :: temp

      integer :: i_lsc, i_meso, i_cell, i_shallow

      i_lsc = Moist_clouds_block%index_strat
      i_cell = Moist_clouds_block%index_donner_cell
      i_meso = Moist_clouds_block%index_donner_meso
      i_shallow = Moist_clouds_block%index_uw_conv

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                  GENERAL MOISTURE DIAGNOSTICS
!
!    output diagnostics reflecting the combination of convective and
!    large-scale parameterizations.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


      kx = size(Output_mp%rdt,3)

!-----------------------------------------------------------------------
!    if any wet deposition diagnostics are desired, execute the following.
!-----------------------------------------------------------------------
      if (wetdep_diagnostics_desired) then

!--------------------------------------------------------------------
!    obtain wet deposition removal for each tracer species n from:
!     a) all precip sources combined;
!     b) from uw convection;
!     c) from donner convection.
!    the wet deposition removal of each species by ls precip is available
!    in ls_wetdep. The arrays holding the wet depo removal by the various
!    convective schemes only have entries for those tracers designated as
!    being affected by the particular scheme,
!--------------------------------------------------------------------
!       total_wetdep = Removal_mp%ls_wetdep
        total_wetdep = 0.
        total_wetdep_dust(:,:) = 0.
        total_wetdep_seasalt(:,:) = 0.

        m = 1
        mm = 1
        do n=1,size(Output_mp%rdt,4)
          if (Removal_mp_control%tracers_in_donner(n) .and.   &
                                                    do_donner_deep) then
            total_wetdep(:,:,n) = total_wetdep(:,:,n) +  &
                                            Removal_mp%donner_wetdep(:,:,m)
            total_wetdep_donner(:,:,n) = Removal_mp%donner_wetdep(:,:,m)
            total_wetdepc_donner(:,:,n) = Removal_mp%donner_wetdepc(:,:,m)
            total_wetdepm_donner(:,:,n) = Removal_mp%donner_wetdepm(:,:,m)
            m = m + 1
          else
            total_wetdep_donner(:,:,n) = 0.
          endif
          if (Removal_mp_control%tracers_in_uw(n) .and. do_uw_conv) then
            total_wetdep(:,:,n) = total_wetdep(:,:,n) +  &
                                               Removal_mp%uw_wetdep(:,:,mm)
            total_wetdep_uw    (:,:,n) = Removal_mp%uw_wetdep(:,:,mm)
            mm = mm + 1
          else
            total_wetdep_uw    (:,:,n) = 0.
          endif
          total_wetdep(:,:,n) = conv_wetdep(n)*total_wetdep(:,:,n)
        end do

!------------------------------------------------------------------------
!    Add in wet dep from large-scale clouds, which already has the proper
!    conversion.
!------------------------------------------------------------------------
        total_wetdep =  total_wetdep + Removal_mp%ls_wetdep

!-----------------------------------------------------------------------
!    process wet deposition removal for each tracer, if requested.
!-----------------------------------------------------------------------
        do n=1, size(Output_mp%rdt,4)
          if (id_wetdep(n) > 0) then
            used = send_data (id_wetdep(n), total_wetdep(:,:,n),  &
                                                           Time, is, js)
          endif

          if (id_wetdep_kg_m2_s(n) > 0) then
            used = send_data (id_wetdep_kg_m2_s(n), total_wetdep(:,:,n)*conv_wetdep_kg_m2_s(n),  &
                                                                                    Time, is, js)
          endif
        end do

        do n=1, size(Output_mp%rdt,4)
          if (id_wetdep_donner(n) > 0) then
             used = send_data (id_wetdep_donner(n),   &
                      total_wetdep_donner(:,:,n)*conv_wetdep(n),  &
                                                             Time, is, js)
          endif
        end do

        do n=1, size(Output_mp%rdt,4)
          if (id_wetdepm_donner(n) > 0) then
            used = send_data (id_wetdepm_donner(n),    &
                       total_wetdepm_donner(:,:,n)*conv_wetdep(n),  &
                                                            Time, is, js)
           endif
        end do

        do n=1, size(Output_mp%rdt,4)
          if (id_wetdepc_donner(n) > 0) then
            used = send_data (id_wetdepc_donner(n),   &
                      total_wetdepc_donner(:,:,n)*conv_wetdep(n),  &
                                                             Time, is, js)
          endif
        end do

        do n=1, size(Output_mp%rdt,4)
          if (id_wetdep_uw(n) > 0) then
            used = send_data (id_wetdep_uw(n),   &
                      total_wetdep_uw(:,:,n)*conv_wetdep(n),  &
                                                             Time, is, js)
          endif
        end do

!-------------------------------------------------------------------------
!   process diagnostics for specific tracer species / families.
!-------------------------------------------------------------------------
        if (id_wetdep_om > 0) then
          used = send_data (id_wetdep_om,  &
               total_wetdep       (:,:,nomphilic) + &
               total_wetdep       (:,:,nomphobic) , &
                                                Time, is,js)
        endif
     if (id_wetpoa_cmip > 0) then
       used = send_data (id_wetpoa_cmip,  &
               total_wetdep(:,:,nomphilic) + total_wetdep(:,:,nomphobic), Time, is,js)
     endif

     if (id_wetdep_SOA > 0) then
       used = send_data (id_wetdep_SOA, total_wetdep(:,:,nSOA) , Time, is,js)
     endif
     if (id_wetsoa_cmip > 0) then
       used = send_data (id_wetsoa_cmip, total_wetdep(:,:,nSOA) , Time, is,js)
     endif

     if (id_wetdep_bc > 0) then
       used = send_data (id_wetdep_bc,  &
               total_wetdep       (:,:,nbcphilic) + &
               total_wetdep       (:,:,nbcphobic) , &
                                                Time, is,js)
     endif
     if (id_wetbc_cmip > 0) then
       used = send_data (id_wetbc_cmip,  &
               total_wetdep(:,:,nbcphilic) + total_wetdep(:,:,nbcphobic), Time, is,js)
     endif

     if (id_wetdep_so4 > 0 .or. id_wetso4_cmip > 0) then
       temp_2d = 0.0
       if( do_donner_deep ) temp_2d = temp_2d + (96.0/WTMAIR)*total_wetdep_donner(:,:,nso4)
       if( do_uw_conv  )    temp_2d = temp_2d + (96.0/WTMAIR)*total_wetdep_uw    (:,:,nso4)
       if( doing_prog_clouds )       temp_2d = temp_2d + 0.096*Removal_mp%ls_wetdep(:,:,nso4)
       if (id_wetdep_so4  > 0) used = send_data (id_wetdep_so4,  temp_2d, Time, is,js)
       if (id_wetso4_cmip > 0) used = send_data (id_wetso4_cmip, temp_2d, Time, is,js)
     endif

     if (id_wetdep_so2 > 0 .or. id_wetso2_cmip > 0) then
       temp_2d = 0.0
       if( do_donner_deep ) temp_2d = temp_2d + (64.0/WTMAIR)*total_wetdep_donner(:,:,nso2)
       if( do_uw_conv  )    temp_2d = temp_2d + (64.0/WTMAIR)*total_wetdep_uw    (:,:,nso2)
       if( doing_prog_clouds )       temp_2d = temp_2d + 0.064*Removal_mp%ls_wetdep(:,:,nso2)
       if (id_wetdep_so2  > 0) used = send_data (id_wetdep_so2,  temp_2d, Time, is,js)
       if (id_wetso2_cmip > 0) used = send_data (id_wetso2_cmip, temp_2d, Time, is,js)
     endif

     if (id_wetdep_DMS > 0 .or. id_wetdms_cmip > 0) then
       temp_2d = 0.0
       if( do_donner_deep ) temp_2d = temp_2d + (62.0/WTMAIR)*total_wetdep_donner(:,:,nDMS)
       if( do_uw_conv  )    temp_2d = temp_2d + (62.0/WTMAIR)*total_wetdep_uw    (:,:,nDMS)
       if( doing_prog_clouds )       temp_2d = temp_2d + 0.062*Removal_mp%ls_wetdep(:,:,nDMS)
       if (id_wetdep_DMS  > 0) used = send_data (id_wetdep_DMS,  temp_2d, Time, is,js)
       if (id_wetdms_cmip > 0) used = send_data (id_wetdms_cmip, temp_2d, Time, is,js)
     endif

     if (id_wetdep_NH4NO3 > 0 .or. id_wetnh4_cmip > 0) then
       temp_2d = 0.0
       if( do_donner_deep ) temp_2d = temp_2d + (18.0/WTMAIR)*(total_wetdep_donner(:,:,nnH4NO3) + &
                                                               total_wetdep_donner(:,:,nNH4) )
       if( do_uw_conv  )    temp_2d = temp_2d + (18.0/WTMAIR)*(total_wetdep_uw(:,:,nNH4NO3) + &
                                                               total_wetdep_uw(:,:,nNH4) )
       if( doing_prog_clouds )       temp_2d = temp_2d + 0.018*(Removal_mp%ls_wetdep(:,:,nNH4NO3) + Removal_mp%ls_wetdep(:,:,nNH4))
       if (id_wetdep_NH4NO3 > 0) used = send_data (id_wetdep_NH4NO3, temp_2d, Time, is,js)
       if (id_wetnh4_cmip   > 0) used = send_data (id_wetnh4_cmip,   temp_2d, Time, is,js)
     endif

     if (id_wetdep_seasalt > 0 .or. id_wetss_cmip > 0) then
       do n=1, n_seasalt_tracers
         nbin_seasalt=seasalt_tracers(n)%tr
         total_wetdep_seasalt(:,:)=total_wetdep_seasalt(:,:)+total_wetdep(:,:,nbin_seasalt)
       enddo
       if (id_wetdep_seasalt > 0) used = send_data (id_wetdep_seasalt, total_wetdep_seasalt, Time, is,js)
       if (id_wetss_cmip     > 0) used = send_data (id_wetss_cmip,     total_wetdep_seasalt, Time, is,js)
     endif

     if (id_wetdep_dust > 0 .or. id_wetdust_cmip > 0) then
     do n=1, n_dust_tracers
        nbin_dust=dust_tracers(n)%tr
        total_wetdep_dust(:,:)=total_wetdep_dust(:,:)+total_wetdep(:,:,nbin_dust)
     enddo
     call atmos_dust_wetdep_flux_set(total_wetdep_dust, is,ie,js,je)
       if (id_wetdep_dust  > 0) used = send_data (id_wetdep_dust,  total_wetdep_dust, Time, is,js)
       if (id_wetdust_cmip > 0) used = send_data (id_wetdust_cmip, total_wetdep_dust, Time, is,js)
     endif

      endif ! (wetdep_diagnostics_desired)

!---------------------------------------------------------------------
!    total precipitation (all sources):
!---------------------------------------------------------------------
      precip = Output_mp%fprec + Output_mp%lprec
      if (id_precip > 0) then
        used = send_data (id_precip, precip, Time, is, js)
      endif
      if (id_pr > 0) then
        used = send_data (id_pr, precip, Time, is, js)
      endif

!---------------------------------------------------------------------
!    snowfall rate due to all sources:
!---------------------------------------------------------------------
      used = send_data (id_snow_tot, Output_mp%fprec, Time, is, js)
      used = send_data (id_prsn, Output_mp%fprec, Time, is, js)
   if (id_prsn_g > 0)  call buffer_global_diag (id_prsn_g, Output_mp%fprec, Time, is, js)

!---------------------------------------------------------------------
!    rainfall rate due to all sources:
!---------------------------------------------------------------------
    used = send_data (id_prra, Output_mp%lprec, Time, is, js)
!---------------------------------------------------------------------
!    column integrated enthalpy and total water tendencies due to
!    moist processes and their imbalances:
!---------------------------------------------------------------------

!------------------------------------------------------------------------
! enthalpy

      if (id_enth_moist_col > 0 .or. id_max_enthalpy_imbal > 0) then
        temp_3d1 = Output_mp%tdt - tdt_init
        temp_3d2 = Output_mp%rdt(:,:,:,             nql) -   &
                                        qdt_init(:,:,:,             nql)
        temp_3d3 = Output_mp%rdt(:,:,:,             nqi) -   &
                                        qdt_init(:,:,:,             nqi)
        temp_2d(:,:) = -HLV*precip -HLF*Output_mp%fprec
        call column_diag(id_enth_moist_col, is, js, Time, temp_3d1,   &
                         CP_AIR, temp_3d2, -HLV, temp_3d3, -HLS,   &
                                                  Input_mp%pmass, temp_2d)

!-----------------------------------------------------------------------
!    this diagnostic captures the extreme value seen during the model
!    segment currently being run.
!-----------------------------------------------------------------------
        if (id_max_enthalpy_imbal > 0) then
          max_enthalpy_imbal(is:ie,js:je) =   &
                      max( abs(temp_2d), max_enthalpy_imbal(is:ie,js:je) )
          used = send_data(id_max_enthalpy_imbal,   &
                             max_enthalpy_imbal(is:ie,js:je), Time, is, js)
        endif
      endif

!------------------------------------------------------------------------
!  total water

      if (id_wat_moist_col > 0 .or. id_max_water_imbal > 0) then
        temp_3d1 = Output_mp%rdt(:,:,:,1) - qdt_init(:,:,:,1)
        temp_3d2 = Output_mp%rdt(:,:,:,             nql) -   &
                                         qdt_init(:,:,:,             nql)
        temp_3d3 = Output_mp%rdt(:,:,:,             nqi) -   &
                                          qdt_init(:,:,:,             nqi)
        temp_2d(:,:) = precip
        call column_diag(id_wat_moist_col, is, js, Time, temp_3d1, 1.0, &
                     temp_3d2, 1.0, temp_3d3, 1.0, Input_mp%pmass, temp_2d)

!-----------------------------------------------------------------------
!    this diagnostic captures the extreme value seen during the model
!    segment currently being run.
!-----------------------------------------------------------------------
        if (id_max_water_imbal > 0) then
          max_water_imbal(is:ie,js:je) =   &
                        max( abs(temp_2d), max_water_imbal(is:ie,js:je) )
          used = send_data(id_max_water_imbal,    &
                               max_water_imbal(is:ie,js:je), Time, is, js)
        endif
      endif

!----------------------------------------------------------------------
!    define total_cloud_area (ls plus convective).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!    add convective cloud area contributions to total_cloud_area. define
!    total convective cloud area (conv_cld_frac) and grid-box mean
!    convective cloud condensate mass (tot_conv_cloud, tot_conv_liq,
!    tot_conv_ice).
!----------------------------------------------------------------------
      total_cloud_area = 0.
      total_conv_cloud =  0.
      conv_cld_frac = 0.
      tot_conv_liq = 0.
      tot_conv_ice = 0.
      if (i_lsc > 0) then
        total_cloud_area = total_cloud_area +   &
                           Moist_clouds_block%cloud_data(i_lsc)%cloud_area
      endif
      if (i_cell > 0) then
        total_cloud_area = total_cloud_area +   &
                         Moist_clouds_block%cloud_data(i_cell)%cloud_area
        total_conv_cloud =  total_conv_cloud + &
             Moist_clouds_block%cloud_data(i_cell)%cloud_area*   &
                     Moist_clouds_block%cloud_data(i_cell)%ice_amt  +  &
             Moist_clouds_block%cloud_data(i_cell)%cloud_area*   &
                     Moist_clouds_block%cloud_data(i_cell)%liquid_amt
        conv_cld_frac = conv_cld_frac +   &
               Moist_clouds_block%cloud_data(i_cell)%cloud_area
        tot_conv_liq =  tot_conv_liq +  &
          Moist_clouds_block%cloud_data(i_cell)%cloud_area*    &
                     Moist_clouds_block%cloud_data(i_cell)%liquid_amt
        tot_conv_ice =  tot_conv_ice +  &
          Moist_clouds_block%cloud_data(i_cell)%cloud_area*    &
                     Moist_clouds_block%cloud_data(i_cell)%ice_amt
      endif
      if (i_meso > 0) then
        total_cloud_area = total_cloud_area +   &
                          Moist_clouds_block%cloud_data(i_meso)%cloud_area
        total_conv_cloud =  total_conv_cloud + &
             Moist_clouds_block%cloud_data(i_meso)%cloud_area*   &
                     Moist_clouds_block%cloud_data(i_meso)%ice_amt  +  &
             Moist_clouds_block%cloud_data(i_meso)%cloud_area*   &
                     Moist_clouds_block%cloud_data(i_meso)%liquid_amt
        conv_cld_frac = conv_cld_frac +   &
               Moist_clouds_block%cloud_data(i_meso)%cloud_area
        tot_conv_liq =  tot_conv_liq +  &
          Moist_clouds_block%cloud_data(i_meso)%cloud_area*    &
                     Moist_clouds_block%cloud_data(i_meso)%liquid_amt
        tot_conv_ice =  tot_conv_ice +  &
          Moist_clouds_block%cloud_data(i_meso)%cloud_area*    &
                     Moist_clouds_block%cloud_data(i_meso)%ice_amt
      endif
      if (i_shallow > 0) then
        total_cloud_area = total_cloud_area +   &
                       Moist_clouds_block%cloud_data(i_shallow)%cloud_area
        total_conv_cloud =  total_conv_cloud + &
             Moist_clouds_block%cloud_data(i_shallow)%cloud_area*   &
                     Moist_clouds_block%cloud_data(i_shallow)%ice_amt  +  &
             Moist_clouds_block%cloud_data(i_shallow)%cloud_area*   &
                     Moist_clouds_block%cloud_data(i_shallow)%liquid_amt
        conv_cld_frac = conv_cld_frac +   &
               Moist_clouds_block%cloud_data(i_shallow)%cloud_area
        tot_conv_liq =  tot_conv_liq +  &
          Moist_clouds_block%cloud_data(i_shallow)%cloud_area*    &
                     Moist_clouds_block%cloud_data(i_shallow)%liquid_amt
        tot_conv_ice =  tot_conv_ice +  &
          Moist_clouds_block%cloud_data(i_shallow)%cloud_area*    &
                     Moist_clouds_block%cloud_data(i_shallow)%ice_amt
      endif



!------------------------------------------------------------------------
!    generate ls cloud diagnostics, normalized by total cloud, as requested
!    in CMIP5.
!------------------------------------------------------------------------
      if (i_lsc > 0) then
        if (id_lsc_cloud_area > 0 .and.   &
                               doing_prog_clouds) then
          used = send_data (id_lsc_cloud_area,   &
                      100.*Moist_clouds_block%cloud_data(i_lsc)%cloud_area,  &
                                                         Time, is, js, 1)
        end if
        if (id_lsc_liq_amt > 0 .and. doing_prog_clouds ) then
          used = send_data (id_lsc_liq_amt,  &
            Moist_clouds_block%cloud_data(i_lsc)%liquid_amt/    &
                                               (1.0 + total_conv_cloud), &
                                                         Time, is, js, 1)
        endif

        if ( id_lsc_ice_amt  > 0 .and.   &
                                             doing_prog_clouds ) then
          used = send_data (id_lsc_ice_amt,   &
            Moist_clouds_block%cloud_data(i_lsc)%ice_amt/     &
                                               (1.0 + total_conv_cloud), &
                                                         Time, is, js, 1)
        endif
      endif  ! (i_lsc > 0)

!-------------------------------------------------------------------------
!   generate the 3d total and convective cloud fraction diagnostics.
!-------------------------------------------------------------------------
    if ( id_tot_cloud_area > 0 ) &
      used = send_data (id_tot_cloud_area, 100.*total_cloud_area,  &
                                                         Time, is, js, 1)

    if (query_cmip_diag_id(ID_cl)) then
      used = send_cmip_data_3d (ID_cl, 100.*total_cloud_area,  &
                                          Time, is, js, 1)!, rmask=mask)
    endif

    if ( id_conv_cloud_area > 0 ) &
    used = send_data (id_conv_cloud_area, 100.*conv_cld_frac, &
                                           Time, is, js, 1)!, rmask=mask)

!---------------------------------------------------------------------
!    define the total 2D cloud area.
!---------------------------------------------------------------------
      if (id_tot_cld_amt > 0 ) then
        tca2 = 1.0
        do k=1,kx
          tca2(:,:) = tca2(:,:)*(1.0 - total_cloud_area(:,:,k))
        end do
        tca2 = 100.*(1. - tca2)
        used = send_data (id_tot_cld_amt, tca2, Time, is, js)
      endif

      if (id_clt > 0 ) then
        tca2 = 1.0
        do k=1,kx
          tca2(:,:) = tca2(:,:)*(1.0 - total_cloud_area(:,:,k))
        end do
        tca2 = 100.*(1. - tca2) ! cmip6 = Cloud Cover Percentage
        used = send_data (id_clt, tca2, Time, is, js)
      endif

      IF (i_lsc > 0) then
!---------------------------------------------------------------------
!    define the total and convective liquid and liquid water path.
!---------------------------------------------------------------------
        if (id_tot_liq_amt > 0 ) &
          used = send_data (id_tot_liq_amt, &
          (Moist_clouds_block%cloud_data(i_lsc)%liquid_amt + tot_conv_liq)/ &
                                                (1.0 + total_conv_cloud), &
                                                           Time, is, js, 1)
        if (query_cmip_diag_id(ID_clw)) then
          used = send_cmip_data_3d (ID_clw, &
    !                (lsc_liquid + tot_conv_liq )/(1.0 + total_conv_cloud), &
                 (Moist_clouds_block%cloud_data(i_lsc)%liquid_amt + tot_conv_liq + Moist_clouds_block%cloud_data(i_lsc)%rain)/(1.0 + total_conv_cloud), &
                                           Time, is, js, 1)!, rmask=mask)
        endif

        if (id_conv_liq_amt > 0 ) &
          used = send_data (id_conv_liq_amt, &
                    tot_conv_liq /(1.0 + total_conv_cloud), &
                                                          Time, is, js, 1)

        if (id_LWP_all_clouds > 0 ) &
          call column_diag (id_LWP_all_clouds, is, js, Time, &
          Moist_clouds_block%cloud_data(i_lsc)%liquid_amt+tot_conv_liq+   &
                         Moist_clouds_block%cloud_data(i_lsc)%rain,    &
                                                      1.0, Input_mp%pmass)

!---------------------------------------------------------------------
!    define the total and convective ice and ice water path.
!---------------------------------------------------------------------
        if (id_tot_ice_amt > 0 ) &
          used = send_data (id_tot_ice_amt, &
          (Moist_clouds_block%cloud_data(i_lsc)%ice_amt + tot_conv_ice)/   &
                              (1.0 + total_conv_cloud), &
                                                         Time, is, js, 1)

        if (query_cmip_diag_id(ID_cli)) then
           used = send_cmip_data_3d (ID_cli, &
!                    (lsc_ice + tot_conv_ice)/(1.0 + total_conv_cloud), &
         (Moist_clouds_block%cloud_data(i_lsc)%ice_amt + tot_conv_ice+ Moist_clouds_block%cloud_data(i_lsc)%snow)/(1.0 + total_conv_cloud), &
                                            Time, is, js, 1)!, rmask=mask)
        endif

        if (id_conv_ice_amt > 0 ) &
          used = send_data (id_conv_ice_amt, &
                  tot_conv_ice/(1.0 + total_conv_cloud), &
                                                         Time, is, js, 1)

        if (id_IWP_all_clouds > 0 ) &
          call column_diag (id_IWP_all_clouds, is, js, Time, &
         Moist_clouds_block%cloud_data(i_lsc)%ice_amt+tot_conv_ice+    &
         Moist_clouds_block%cloud_data(i_lsc)%snow, 1.0, Input_mp%pmass)

        if (id_clivi > 0 ) &
          call column_diag (id_clivi, is, js, Time, &
           Moist_clouds_block%cloud_data(i_lsc)%ice_amt+tot_conv_ice+    &
           Moist_clouds_block%cloud_data(i_lsc)%snow, 1.0, Input_mp%pmass)

!---------------------------------------------------------------------
!    define the total water substance and condensate water path.
!---------------------------------------------------------------------
        used = send_data (id_tot_h2o  , &
              (Input_mp%qin(:,:,:) +   &
               Moist_clouds_block%cloud_data(i_lsc)%ice_amt +   &
                         tot_conv_ice +    &
               Moist_clouds_block%cloud_data(i_lsc)%liquid_amt +  &
                         tot_conv_liq)/(1.0 + total_conv_cloud), &
                                                          Time, is, js, 1)

        if (id_WP_all_clouds > 0 ) &
          call column_diag(id_WP_all_clouds, is, js, Time, &
               Moist_clouds_block%cloud_data(i_lsc)%ice_amt +   &
                            tot_conv_ice +    &
               Moist_clouds_block%cloud_data(i_lsc)%liquid_amt +  &
                            tot_conv_liq + &
               Moist_clouds_block%cloud_data(i_lsc)%rain +   &
               Moist_clouds_block%cloud_data(i_lsc)%snow, &
                                                       1.0, Input_mp%pmass)

        if (id_clwvi > 0 ) &
          call column_diag(id_clwvi, is, js, Time, &
             Moist_clouds_block%cloud_data(i_lsc)%ice_amt +   &
                          tot_conv_ice +    &
             Moist_clouds_block%cloud_data(i_lsc)%liquid_amt +  &
                          tot_conv_liq + &
             Moist_clouds_block%cloud_data(i_lsc)%rain +   &
             Moist_clouds_block%cloud_data(i_lsc)%snow, &
                                                       1.0, Input_mp%pmass)

      ELSE  ! (i_lsc > 0)
!---------------------------------------------------------------------
!    define the total and convective liquid and liquid water path.
!---------------------------------------------------------------------
        if (id_tot_liq_amt > 0 ) &
          used = send_data (id_tot_liq_amt, &
            tot_conv_liq/(1.0 + total_conv_cloud), Time, is, js, 1)

        if (query_cmip_diag_id(ID_clw)) then
          used = send_cmip_data_3d (ID_clw, &
!          (lsc_liquid + tot_conv_liq )/(1.0 + total_conv_cloud), &
           (Moist_clouds_block%cloud_data(i_lsc)%liquid_amt + tot_conv_liq + Moist_clouds_block%cloud_data(i_lsc)%rain)/(1.0 + total_conv_cloud), &
                                           Time, is, js, 1)!, rmask=mask)
    endif

        if (id_conv_liq_amt > 0 ) &
          used = send_data (id_conv_liq_amt, &
                tot_conv_liq /(1.0 + total_conv_cloud), Time, is, js, 1)

        if (id_LWP_all_clouds > 0 ) &
          call column_diag (id_LWP_all_clouds, is, js, Time, &
                           tot_conv_liq, 1.0, Input_mp%pmass)

!---------------------------------------------------------------------
!    define the total and convective ice and ice water path.
!---------------------------------------------------------------------
        if (id_tot_ice_amt > 0 ) &
          used = send_data (id_tot_ice_amt, &
            tot_conv_ice/(1.0 + total_conv_cloud), Time, is, js, 1)

    if (query_cmip_diag_id(ID_cli)) then
       used = send_cmip_data_3d (ID_cli, &
!                (lsc_ice + tot_conv_ice)/(1.0 + total_conv_cloud), &
     (Moist_clouds_block%cloud_data(i_lsc)%ice_amt + tot_conv_ice+ Moist_clouds_block%cloud_data(i_lsc)%snow)/(1.0 + total_conv_cloud), &
                                            Time, is, js, 1)!, rmask=mask)
    endif

        if (id_conv_ice_amt > 0 ) &
          used = send_data (id_conv_ice_amt, &
             tot_conv_ice/(1.0 + total_conv_cloud), Time, is, js, 1)

        if (id_IWP_all_clouds > 0 ) &
          call column_diag (id_IWP_all_clouds, is, js, Time, &
                            tot_conv_ice, 1.0, Input_mp%pmass)

        if (id_clivi > 0 ) &
          call column_diag (id_clivi, is, js, Time, &
                          tot_conv_ice, 1.0, Input_mp%pmass)

!---------------------------------------------------------------------
!    define the total water substance and condensate water path.
!---------------------------------------------------------------------
        used = send_data (id_tot_h2o  , &
              (Input_mp%qin(:,:,:) + tot_conv_ice +    &
                  tot_conv_liq)/(1.0 + total_conv_cloud), Time, is, js, 1)

        if (id_WP_all_clouds > 0 ) &
          call column_diag(id_WP_all_clouds, is, js, Time, &
                      tot_conv_ice + tot_conv_liq, 1.0, Input_mp%pmass)

        if (id_clwvi > 0 ) &
          call column_diag(id_clwvi, is, js, Time, &
                      tot_conv_ice + tot_conv_liq, 1.0, Input_mp%pmass)

      ENDIF ! (i_lsc > 0)

!---------------------------------------------------------------------
!    define the water vapor path and total vapor.
!---------------------------------------------------------------------
      if (id_WVP > 0) &
           call column_diag(id_WVP, is, js, Time, Input_mp%qin, 1.0,  &
                                                         Input_mp%pmass)
      if (id_prw > 0) &
           call column_diag(id_prw, is, js, Time, Input_mp%qin, 1.0,  &
                                                         Input_mp%pmass)

      used = send_data (id_tot_vapor, Input_mp%qin, Time, is, js, 1)

!---------------------------------------------------------------------
!    column integrated cloud mass:
!---------------------------------------------------------------------
      if (id_AWP > 0 .and. doing_prog_clouds) &
        call column_diag(id_AWP, is, js, Time,    &
              Input_mp%tracer(:,:,:,             nqa), 1.0, Input_mp%pmass)

!---------------------------------------------------------------------
!    output the global integral of precipitation in units of mm/day.
!---------------------------------------------------------------------
      prec_intgl(is:ie,js:je) = precip(:,:)*SECONDS_PER_DAY
      if (id_pr_g > 0)  call buffer_global_diag (id_pr_g,  precip(:,:), Time, is, js)
!---------------------------------------------------------------------
!    relative humidity:
!---------------------------------------------------------------------
      if (id_rh > 0) then
        if (.not. (       do_rh_clouds                           )) then
          call rh_calc (Input_mp%pfull, Input_mp%tin(:,:,:),  &
                   Input_mp%qin(:,:,:), RH(:,:,:), do_simple )
           used = send_data (id_rh, rh*100., Time, is, js, 1)
        endif
      endif

!---------------------------------------------------------------------
!    relative humidity (CMIP formulation):
!---------------------------------------------------------------------
      if (id_rh_cmip > 0) then
        if (.not. (       do_rh_clouds                           )) then
          call rh_calc (Input_mp%pfull, input_mp%tin, Input_mp%qin, RH, &
                                   .false.,      do_cmip=.true.)
          used = send_data (id_rh_cmip, rh*100., Time, is, js, 1)
        endif
      endif

    if (query_cmip_diag_id(ID_hur)) then
      if (.not. (do_rh_clouds .or. do_diag_clouds)) then
        call rh_calc (Input_mp%pfull, input_mp%tin, Input_mp%qin, RH, &
                                        .false., do_cmip=.true.)
      endif
      used = send_cmip_data_3d (ID_hur, rh*100., Time, is, js, 1)!, rmask=mask)
    endif

!---------------------------------------------------------------------
!    saturation specific humidity:
!---------------------------------------------------------------------
      if (id_qs > 0) then
        call compute_qs (Input_mp%tin, Input_mp%pfull, qsat, &
                                                      q = Input_mp%qin)
        used = send_data (id_qs, qsat(:,:,:), Time, is, js, 1)
      endif

!------------------------------------------------------------------------
!   call routine which calculates CAPE and CIN using model fields after
!   adjustment by moist_processes.
!------------------------------------------------------------------------
      call cape_cin_diagnostics (is,ie,js,je, Input_mp, Time)

!---------------------------------------------------------------------


end subroutine combined_MP_diagnostics



!#######################################################################

subroutine MP_alloc (Physics_input_block, Physics_tendency_block, &
                     Phys_mp_exch, dt, area, lon, lat, land, ustar, &
                     bstar, qstar, Input_mp,Tend_mp, C2ls_mp, Output_mp, &
                     Removal_mp)

type(physics_input_block_type),                          &
                          intent(in)    :: Physics_input_block
type(physics_tendency_block_type),                       &
                          intent(in)    :: Physics_tendency_block
type(phys_mp_exch_type),  intent(in)    :: Phys_mp_exch
real,                     intent(in)    :: dt
real, dimension(:,:),     intent(in)    :: area, lon, lat
real, dimension(:,:),     intent(in)    :: land, ustar, bstar, qstar
type(MP_input_type),      intent(inout) :: Input_mp
type(MP_output_type),     intent(inout) :: Output_mp
type(MP_tendency_type),   intent(inout) :: Tend_mp
type(MP_conv2ls_type),    intent(inout) :: C2ls_mp
type(MP_removal_type),    intent(inout) :: Removal_mp


!------------------------------------------------------------------------
!   local variables:

      real, dimension(size(land,1),size(land,2) )  :: tsnow
      real, dimension(size(land,1),size(land,2),    &
                  size(Physics_input_block%t,3) )  :: d_zfull
      real, dimension(size(land,1),size(land,2),    &
                  size(Physics_input_block%t,3)+1 )  :: d_zhalf
      integer :: ix, jx, kx, nt
      integer :: k, tr

!-------------------------------------------------------------------------
!    define input array sizes.
!-------------------------------------------------------------------------
      ix = size(Physics_input_block%t,1)
      jx = size(Physics_input_block%t,2)
      kx = size(Physics_input_block%t,3)
      nt = size(Physics_tendency_block%q_dt,4)

!------------------------------------------------------------------------
!    allocate and initialize (or associate where possible) an mp_input_type
!    variable which will contain atmospheric field inputs needed in
!    moist_processes.
!------------------------------------------------------------------------
      Input_mp%phalf => Physics_input_block%p_half
      Input_mp%pfull => Physics_input_block%p_full
      Input_mp%zhalf => Physics_input_block%z_half
      Input_mp%zfull => Physics_input_block%z_full
      allocate (Input_mp%tin   (ix,jx,kx  ))
      allocate (Input_mp%qin   (ix,jx,kx  ))
      allocate (Input_mp%uin   (ix,jx,kx  ))
      allocate (Input_mp%vin   (ix,jx,kx  ))
      Input_mp%t => Physics_input_block%t
      Input_mp%q => Physics_input_block%q(:,:,:,1)
      Input_mp%u => Physics_input_block%u
      Input_mp%v => Physics_input_block%v
      Input_mp%w => Physics_input_block%w
      if (associated(Physics_input_block%um)) then
      Input_mp%tm => Physics_input_block%tm
      Input_mp%qm => Physics_input_block%qm (:,:,:,1)
      Input_mp%um => Physics_input_block%um
      Input_mp%vm => Physics_input_block%vm
      Input_mp%rm => Physics_input_block%qm
      else
      Input_mp%tm => Physics_input_block%t
      Input_mp%qm => Physics_input_block%q (:,:,:,1)
      Input_mp%um => Physics_input_block%u
      Input_mp%vm => Physics_input_block%v
      Input_mp%rm => Physics_input_block%q
      endif

      Input_mp%omega => Physics_input_block%omega

      Input_mp%r => Physics_input_block%q


      Input_mp%radturbten => Phys_mp_exch%radturbten
      Input_mp%diff_t => Phys_mp_exch%diff_t
      allocate (Input_mp%tracer(ix,jx,kx, size(Physics_input_block%q,4) ))
      allocate (Input_mp%area  (ix,jx  ))   ; Input_mp%area   = area
      allocate (Input_mp%lon   (ix,jx  ))   ; Input_mp%lon    = lon
      allocate (Input_mp%lat   (ix,jx  ))   ; Input_mp%lat    = lat
      allocate (Input_mp%land  (ix,jx  ))   ; Input_mp%land   = land
      Input_mp%cush  => Phys_mp_exch%cush
      Input_mp%cbmf  => Phys_mp_exch%cbmf
      Input_mp%pblht => Phys_mp_exch%pbltop
      allocate (Input_mp%ustar (ix,jx  ))   ; Input_mp%ustar  = ustar
      allocate (Input_mp%bstar (ix,jx  ))   ; Input_mp%bstar  = bstar
      allocate (Input_mp%qstar (ix,jx  ))   ; Input_mp%qstar  = qstar
      Input_mp%tdt_shf    => Phys_mp_exch%tdt_shf
      Input_mp%qdt_lhf    => Phys_mp_exch%qdt_lhf
      allocate (Input_mp%coldT(ix,jx  ))    ; Input_mp%coldT = .false.
      allocate (Input_mp%pmass  (ix,jx,kx))

!---------------------------------------------------------------------
!    define input fields to be used, either the tau time level fields,
!    or the tau time level values updated with the time tendencies thus
!    far calculated on the current step. control is through variable
!    use_tau which was obtained from Physics%control during initialization.
!---------------------------------------------------------------------
      if (use_tau) then
        Input_mp%tin = Physics_input_block%t
        Input_mp%qin= Physics_input_block%q(:,:,:,1)
        Input_mp%uin = Physics_input_block%u
        Input_mp%vin = Physics_input_block%v
        do tr=1,size(Physics_input_block%q,4)
          Input_mp%tracer(:,:,:,tr) = Physics_input_block%q(:,:,:,tr)
        end do
      else
      if (associated(Physics_input_block%um)) then
        Input_mp%tin = Physics_input_block%tm +    &
                                        Physics_tendency_block%t_dt*dt
        Input_mp%qin = Physics_input_block%qm(:,:,:,1) +    &
                                 Physics_tendency_block%q_dt(:,:,:,1)*dt
        Input_mp%uin = Physics_input_block%um +   &
                                 Physics_tendency_block%u_dt*dt
        Input_mp%vin = Physics_input_block%vm +   &
                                 Physics_tendency_block%v_dt*dt
        do tr=1,size(Physics_tendency_block%q_dt,4)
          Input_mp%tracer(:,:,:,tr) = Physics_input_block%qm(:,:,:,tr) + &
                               Physics_tendency_block%q_dt(:,:,:,tr)*dt
        end do
        do tr=size(Physics_tendency_block%q_dt,4) +1,   &
                                              size(Physics_input_block%q,4)
          Input_mp%tracer(:,:,:,tr) = Physics_input_block%q(:,:,:,tr)
        end do

        if (do_height_adjust) then
          call height_adjust (Physics_input_block%tm,    &
                              Physics_input_block%qm(:,:,:,1),  &
                              Physics_input_block%qm,           &
                              Input_mp%tin, Input_mp%qin, Input_mp%tracer,&
                              Input_mp%phalf, Input_mp%pfull, &
                              Physics_input_block%z_half,    &
                              Physics_input_block%z_full,    &
                              d_zhalf, d_zfull)
          Input_mp%zhalf = Physics_input_block%z_half + d_zhalf
          Input_mp%zfull = Physics_input_block%z_full + d_zfull
        else
          Input_mp%zhalf = Physics_input_block%z_half
          Input_mp%zfull = Physics_input_block%z_full
        endif

       else ! associated
        Input_mp%tin = Physics_input_block%t  +    &
                                       Physics_tendency_block%t_dt*dt
        Input_mp%qin = Physics_input_block%q (:,:,:,1) +   &
                                 Physics_tendency_block%q_dt(:,:,:,1)*dt
        Input_mp%uin = Physics_input_block%u  +   &
                                          Physics_tendency_block%u_dt*dt
        Input_mp%vin = Physics_input_block%v  +     &
                                           Physics_tendency_block%v_dt*dt
        do tr=1,size(Physics_tendency_block%q_dt,4)
          Input_mp%tracer(:,:,:,tr) = Physics_input_block%q (:,:,:,tr) +  &
                                   Physics_tendency_block%q_dt(:,:,:,tr)*dt
        end do
        do tr=size(Physics_tendency_block%q_dt,4) +1,   &
                                             size(Physics_input_block%q,4)
          Input_mp%tracer(:,:,:,tr) = Physics_input_block%q(:,:,:,tr)
        end do

        if (do_height_adjust) then
          call height_adjust (Physics_input_block%t ,    &
                              Physics_input_block%q (:,:,:,1),  &
                              Physics_input_block%q ,           &
                              Input_mp%tin, Input_mp%qin, Input_mp%tracer,&
                              Input_mp%phalf, Input_mp%pfull, &
                              Physics_input_block%z_half,    &
                              Physics_input_block%z_full,    &
                              d_zhalf, d_zfull)
          Input_mp%zhalf = Physics_input_block%z_half + d_zhalf
          Input_mp%zfull = Physics_input_block%z_full + d_zfull
        else
          Input_mp%zhalf = Physics_input_block%z_half
          Input_mp%zfull = Physics_input_block%z_full
        endif

       endif ! associated
      endif  !use_tau

!----------------------------------------------------------------------
!    compute the mean temperature in the lower atmosphere (the lowest
!    pdepth Pa), to be used to determine whether rain or snow reaches
!    the surface. define a logical variable coldT indicating whether
!    snow or rain falls in the column.
!    ????    SHOULD TIN BE USED RATHER THAN t ??
!----------------------------------------------------------------------
      call tempavg (Nml_mp%pdepth, Input_mp%phalf, Input_mp%t, tsnow)
      where (tsnow(:,:) <= TFREEZE)
        Input_mp%coldT(:,:) = .true.
      endwhere

!----------------------------------------------------------------------
!    compute the mass in each model layer.
!----------------------------------------------------------------------
      do k=1,kx
        Input_mp%pmass(:,:,k) =    &
                    (Input_mp%phalf(:,:,k+1) - Input_mp%phalf(:,:,k))/GRAV
      end do

!------------------------------------------------------------------------
!    allocate and initialize an mp_tend_type variable which will
!    contain tendencies from moist_processes,
!------------------------------------------------------------------------
      allocate (Tend_mp%ttnd     (ix,jx,kx))   ; Tend_mp%ttnd = 0.
      allocate (Tend_mp%qtnd     (ix,jx,kx))   ; Tend_mp%qtnd = 0.
      allocate (Tend_mp%ttnd_conv(ix,jx,kx))   ; Tend_mp%ttnd_conv = 0.
      allocate (Tend_mp%qtnd_conv(ix,jx,kx))   ; Tend_mp%qtnd_conv = 0.
      allocate (Tend_mp%qtnd_wet (ix,jx,kx))   ; Tend_mp%qtnd_wet = 0.
      allocate (Tend_mp%wetdeptnd(ix,jx,kx))   ; Tend_mp%wetdeptnd  = 0.
      allocate (Tend_mp%qldt_conv(ix,jx,kx))  ; Tend_mp%qldt_conv   = 0.
      allocate (Tend_mp%qidt_conv(ix,jx,kx))   ; Tend_mp%qidt_conv   = 0.
      allocate (Tend_mp%qadt_conv(ix,jx,kx))   ; Tend_mp%qadt_conv   = 0.
      allocate (Tend_mp%qndt_conv(ix,jx,kx))   ; Tend_mp%qndt_conv   = 0.
      allocate (Tend_mp%qnidt_conv(ix,jx,kx))   ; Tend_mp%qnidt_conv   = 0.
      allocate (Tend_mp%q_tnd     (ix,jx,kx,nt))   ; Tend_mp%q_tnd     = 0.

!------------------------------------------------------------------------
!    allocate and initialize an mp_c2ls_type variable which will
!    hold quantities which need to be passed form the convection driver to
!    the large-scale driver.
!------------------------------------------------------------------------
      allocate (C2ls_mp%donner_humidity_area (ix,jx,kx))
                                          C2ls_mp%donner_humidity_area = 0.
      allocate (C2ls_mp%donner_humidity_factor (ix,jx,kx))
                                        C2ls_mp%donner_humidity_factor = 0.
      allocate (C2ls_mp%convective_humidity_area (ix,jx,kx))
                                      C2ls_mp%convective_humidity_area = 0.
      allocate (C2ls_mp%convective_humidity_ratio (ix,jx,kx))
                                     C2ls_mp%convective_humidity_ratio = 0.
      allocate (C2ls_mp%conv_frac_clubb    (ix,jx,kx))
                                          C2ls_mp%conv_frac_clubb      = 0.
      allocate (C2ls_mp%convective_humidity_ratio_clubb (ix,jx,kx))
                              C2ls_mp%convective_humidity_ratio_clubb = 0.
      allocate (C2ls_mp%wet_data (ix,jx,kx,nt))   ; C2ls_mp%wet_data = 0.
      allocate (C2ls_mp%cloud_wet (ix,jx,kx))   ; C2ls_mp%cloud_wet = 0.
      allocate (C2ls_mp%cloud_frac (ix,jx,kx))   ; C2ls_mp%cloud_frac = 0.
      allocate (C2ls_mp%mc_full (ix,jx,kx))   ; C2ls_mp%mc_full    = 0.
      allocate (C2ls_mp%mc_half(ix,jx,kx+1))   ; C2ls_mp%mc_half    = 0.

!------------------------------------------------------------------------
!    allocate and initialize an Mp_output_type variable which will
!    hold quantities which need to be passed from moist_processes back
!    to physics_driver.
!------------------------------------------------------------------------
      Output_mp%tdt  => Physics_tendency_block%t_dt
      Output_mp%udt  => Physics_tendency_block%u_dt
      Output_mp%vdt  => Physics_tendency_block%v_dt
      Output_mp%rdt  => Physics_tendency_block%q_dt
      Output_mp%convect  => Phys_mp_exch%convect
      Output_mp%convect = .false.
      allocate ( Output_mp%lprec  (ix,jx))    ; Output_mp%lprec   = 0.
      allocate ( Output_mp%fprec  (ix,jx))    ; Output_mp%fprec   = 0.
      allocate ( Output_mp%gust_cv(ix,jx))    ; Output_mp%gust_cv = 0.
      Output_mp%diff_t_clubb => Phys_mp_exch%diff_t_clubb
                                              Output_mp%diff_t_clubb =0.
      Output_mp%diff_cu_mo  => Phys_mp_exch%diff_cu_mo
                                             Output_mp%diff_cu_mo  = 0.

!------------------------------------------------------------------------
!    allocate and initialize an mp_removal_type variable which will
!    hold quantities related to the water and tracer removed from the
!    atmosphere by the convective and large-scale cloud processes.
!------------------------------------------------------------------------

      allocate ( Removal_mp%ice_precflx(ix,jx,kx))
                                              Removal_mp%ice_precflx= 0.
      allocate ( Removal_mp%liq_precflx(ix,jx,kx))
                                              Removal_mp%liq_precflx= 0.
      allocate ( Removal_mp%ice_precflxh(ix,jx,kx+1))
                                              Removal_mp%ice_precflxh= 0.
      allocate ( Removal_mp%liq_precflxh(ix,jx,kx+1))
                                              Removal_mp%liq_precflxh= 0.
      allocate ( Removal_mp%frz_meso(ix,jx,kx)) ; Removal_mp%frz_meso= 0.
      allocate ( Removal_mp%liq_meso(ix,jx,kx)) ; Removal_mp%liq_meso= 0.
      allocate ( Removal_mp%frz_cell(ix,jx,kx)) ; Removal_mp%frz_cell= 0.
      allocate ( Removal_mp%liq_cell(ix,jx,kx)) ; Removal_mp%liq_cell= 0.
      allocate ( Removal_mp%frz_mesoh(ix,jx,kx+1))
                                                 Removal_mp%frz_mesoh= 0.
      allocate ( Removal_mp%liq_mesoh(ix,jx,kx+1))
                                                 Removal_mp%liq_mesoh= 0.
      allocate ( Removal_mp%frz_cellh(ix,jx,kx+1))
                                                 Removal_mp%frz_cellh= 0.
      allocate ( Removal_mp%liq_cellh(ix,jx,kx+1))
                                                 Removal_mp%liq_cellh= 0.
      allocate ( Removal_mp%mca_frz  (ix,jx,kx))
                                                 Removal_mp%mca_frz  = 0.
      allocate ( Removal_mp%mca_liq  (ix,jx,kx))
                                                 Removal_mp%mca_liq  = 0.
      allocate ( Removal_mp%mca_frzh (ix,jx,kx+1))
                                                  Removal_mp%mca_frzh = 0.
      allocate ( Removal_mp%mca_liqh (ix,jx,kx+1))
                                                 Removal_mp%mca_liqh = 0.
      allocate ( Removal_mp%rain3d   (ix,jx,kx+1))
                                                 Removal_mp%rain3d   = 0.
      allocate ( Removal_mp%snow3d   (ix,jx,kx+1))
                                                 Removal_mp%snow3d   = 0.
      allocate ( Removal_mp%snowclr3d   (ix,jx,kx+1))
                                               Removal_mp%snowclr3d   = 0.
      allocate ( Removal_mp%uw_wetdep (ix,jx,  &
                                       Removal_mp_control%num_uw_tracers))
                                               Removal_mp%uw_wetdep   = 0.
      allocate ( Removal_mp%donner_wetdep  &
                                  (ix,jx,    &
                                   Removal_mp_control%num_donner_tracers))
                                          Removal_mp%donner_wetdep   = 0.
      allocate ( Removal_mp%donner_wetdepm  &
                                  (ix,jx,    &
                                   Removal_mp_control%num_donner_tracers))
                                         Removal_mp%donner_wetdepm   = 0.
      allocate ( Removal_mp%donner_wetdepc  &
                                  (ix,jx,    &
                                   Removal_mp_control%num_donner_tracers))
                                         Removal_mp%donner_wetdepc   = 0.
      allocate ( Removal_mp%ls_wetdep   (ix,jx,nt))
                                               Removal_mp%ls_wetdep   = 0.

!-------------------------------------------------------------------------

end subroutine MP_alloc


!########################################################################

subroutine MP_dealloc (Input_mp, Tend_mp, C2ls_mp, Output_mp, Removal_mp)

type(MP_input_type),    intent(inout) :: Input_mp
type(MP_tendency_type), intent(inout) :: Tend_mp
type(MP_conv2ls_type),  intent(inout) :: C2ls_mp
type(MP_output_type),   intent(inout) :: Output_mp
type(MP_removal_type),  intent(inout) :: Removal_mp

!------------------------------------------------------------------------
!    deallocate the components of the derived type variables defined in
!    moist_processes.
!------------------------------------------------------------------------

      Input_mp%phalf => null()
      Input_mp%pfull => null()
      Input_mp%zhalf => null()
      Input_mp%zfull => null()

      deallocate (Input_mp%tin   )
      deallocate (Input_mp%qin   )
      deallocate (Input_mp%uin   )
      deallocate (Input_mp%vin   )

      Input_mp%t => null()
      Input_mp%q => null()
      Input_mp%u => null()
      Input_mp%v => null()
      Input_mp%w => null()
      Input_mp%tm => null()
      Input_mp%qm => null()
      Input_mp%um => null()
      Input_mp%vm => null()

      Input_mp%r => null()
      Input_mp%rm => null()
      Input_mp%omega => null()

      deallocate (Input_mp%area  )
      deallocate (Input_mp%lon   )
      deallocate (Input_mp%lat   )
      deallocate (Input_mp%tracer)
      deallocate (Input_mp%land  )
      deallocate (Input_mp%ustar )
      deallocate (Input_mp%bstar )
      deallocate (Input_mp%qstar )
      deallocate (Input_mp%pmass  )
      deallocate (Input_mp%coldT  )

      Input_mp%radturbten => null()
      Input_mp%diff_t => null()
      Input_mp%cush  => null()
      Input_mp%cbmf  => null()
      Input_mp%pblht => null()
      Input_mp%tdt_shf => null()
      Input_mp%qdt_lhf => null()

      deallocate (Tend_mp%ttnd)
      deallocate (Tend_mp%qtnd)
      deallocate (Tend_mp%ttnd_conv)
      deallocate (Tend_mp%qtnd_conv)
      deallocate (Tend_mp%qtnd_wet)
      deallocate (Tend_mp%wetdeptnd )
      deallocate (Tend_mp%qldt_conv )
      deallocate (Tend_mp%qidt_conv )
      deallocate (Tend_mp%qadt_conv )
      deallocate (Tend_mp%qndt_conv )
      deallocate (Tend_mp%qnidt_conv )
      deallocate (Tend_mp%q_tnd      )

      deallocate (C2ls_mp%donner_humidity_area)
      deallocate (C2ls_mp%donner_humidity_factor)
      deallocate (C2ls_mp%convective_humidity_area)
      deallocate (C2ls_mp%convective_humidity_ratio)
      deallocate (C2ls_mp%conv_frac_clubb   )
      deallocate (C2ls_mp%convective_humidity_ratio_clubb)
      deallocate (C2ls_mp%wet_data              )
      deallocate (C2ls_mp%cloud_wet             )
      deallocate (C2ls_mp%cloud_frac            )
      deallocate (C2ls_mp%mc_full               )
      deallocate (C2ls_mp%mc_half               )

      deallocate (Removal_mp%ice_precflx)
      deallocate (Removal_mp%liq_precflx)

      deallocate (Removal_mp%frz_meso    )
      deallocate (Removal_mp%liq_meso    )
      deallocate (Removal_mp%frz_cell    )
      deallocate (Removal_mp%liq_cell    )

      deallocate (Removal_mp%mca_frz     )
      deallocate (Removal_mp%mca_liq     )

      deallocate (Removal_mp%uw_wetdep      )
      deallocate (Removal_mp%donner_wetdep      )
      deallocate (Removal_mp%donner_wetdepm     )
      deallocate (Removal_mp%donner_wetdepc     )
      deallocate (Removal_mp%ls_wetdep      )

      Output_mp%tdt => null()
      Output_mp%udt => null()
      Output_mp%vdt => null()
      Output_mp%rdt  => null()
      deallocate (Output_mp%lprec  )
      deallocate (Output_mp%fprec  )
      deallocate (Output_mp%gust_cv)

      Output_mp%convect => null()
      Output_mp%diff_t_clubb => null()
      Output_mp%diff_cu_mo    => null()

!--------------------------------------------------------------------

end subroutine MP_dealloc



!########################################################################

subroutine create_Nml_mp


      Nml_mp%do_mca =  do_mca
      Nml_mp%do_lsc =  do_lsc
      Nml_mp%do_ras =  do_ras
      Nml_mp%do_uw_conv  = do_uw_conv
      Nml_mp%limit_conv_cloud_frac = limit_conv_cloud_frac
      Nml_mp%do_dryadj = do_dryadj
      Nml_mp%pdepth = pdepth
      Nml_mp%include_donmca_in_cosp  = include_donmca_in_cosp
      Nml_mp%do_simple = do_simple
      Nml_mp%do_rh_clouds = do_rh_clouds
      Nml_mp%do_donner_deep = do_donner_deep
      Nml_mp%do_bm =  do_bm
      Nml_mp%do_bmmass = do_bmmass
      Nml_mp%do_bmomp = do_bmomp
      Nml_mp%do_unified_clouds = do_unified_clouds
      Nml_mp%use_online_aerosol = use_online_aerosol
      Nml_mp%use_sub_seasalt = use_sub_seasalt
      Nml_mp%sea_salt_scale = sea_salt_scale
      Nml_mp%om_to_oc = om_to_oc

!--------------------------------------------------------------------

end subroutine create_Nml_mp


!#######################################################################

subroutine diag_field_init ( axes, Time )

integer,         intent(in) :: axes(4)
type(time_type), intent(in) :: Time
integer                     :: id, ic
integer                     :: id_wetdep_cmip

!-----------------------------------------------------------------------
!   local variables:

      character(len=32)     :: tracer_units, tracer_name
      character(len=128)    :: diaglname, diaglname_uw, diaglname_donner
      integer, dimension(3) :: half = (/1,2,4/)
      integer               :: n, nn, outunit

      character(len=256) :: cmip_name, cmip_longname, cmip_longname2
      logical :: cmip_is_aerosol
      real    :: tracer_mw

!------------ register the diagnostic fields of this module -------------

      id_snow_tot  = register_diag_field ( mod_name, &
        'snow_tot ', axes(1:2), Time, &
        'Frozen precip rate from all sources',       'kg(h2o)/m2/s', &
                              interp_method = "conserve_order1" )

      id_prra  = register_cmip_diag_field_2d ( mod_name, 'prra', Time, &
                                        'Rainfall Rate', 'kg m-2 s-1', &
                                      standard_name = 'rainfall_flux', &
                                     interp_method = "conserve_order1" )

      id_prsn  = register_cmip_diag_field_2d ( mod_name, 'prsn', Time, &
                                        'Snowfall Flux', 'kg m-2 s-1', &
                                      standard_name = 'snowfall_flux', &
                                     interp_method = "conserve_order1" )


      id_max_enthalpy_imbal    = register_diag_field    &
        (mod_name, 'max_enth_imbal', axes(1:2), Time,  &
        'max enthalpy  imbalance from moist_processes  ', 'W/m2',   &
                                     missing_value=missing_value)

      id_max_water_imbal    = register_diag_field    &
        (mod_name, 'max_water_imbal', axes(1:2), Time,   &
        'max water  imbalance from moist_processes  ', 'kg(h2o)/m2/s',  &
                                             missing_value=missing_value)

      id_enth_moist_col = register_diag_field ( mod_name, &
        'enth_moist_col', axes(1:2), Time, &
        'Column enthalpy imbalance from moist processes','W/m2' )

      id_wat_moist_col = register_diag_field ( mod_name, &
        'wat_moist_col', axes(1:2), Time, &
        'Column total water imbalance from moist processes','kg/m2/s' )

      id_precip = register_diag_field ( mod_name, &
        'precip', axes(1:2), Time, &
        'Total precipitation rate',                     'kg/m2/s', &
                                      interp_method = "conserve_order1" )

      id_pr = register_cmip_diag_field_2d ( mod_name, 'pr', Time, &
                                  'Precipitation',  'kg m-2 s-1', &
                              standard_name='precipitation_flux', &
                                interp_method = "conserve_order1" )

      id_WVP = register_diag_field ( mod_name, &
        'WVP', axes(1:2), Time, &
        'Column integrated water vapor',                'kg/m2'  )

      id_prw = register_cmip_diag_field_2d ( mod_name, 'prw', Time, &
                                      'Water Vapor Path', 'kg m-2', &
                   standard_name = 'atmosphere_water_vapor_content' )

!-----------------------------------------------------------------------
!    the following diagnostics only available with prognostic large-scale
!    cloud scheme:
!-----------------------------------------------------------------------
      if ( doing_prog_clouds ) then
        id_AWP = register_diag_field ( mod_name, &
         'AWP', axes(1:2), Time, &
         'Column integrated cloud mass ',                'kg/m2'   )

        id_tot_cld_amt = register_diag_field    &
              (mod_name, 'cld_amt_2d', axes(1:2), Time, &
                'total cloud amount', 'percent')

        id_clt = register_cmip_diag_field_2d (mod_name, 'clt', Time, &
                                      'Total Cloud Cover Percentage', '%', &
                                standard_name= 'cloud_area_fraction' )

        id_tot_cloud_area = register_diag_field ( mod_name, &
          'tot_cloud_area', axes(1:3), Time, &
          'Cloud area -- all clouds', 'percent',    &
                                           missing_value=missing_value )

        ID_cl = register_cmip_diag_field_3d ( mod_name, 'cl', Time, &
                                  'Percentage Cloud Cover', '%',    &
            standard_name='cloud_area_fraction_in_atmosphere_layer' )

        id_tot_h2o     = register_diag_field ( mod_name, &
          'tot_h2o', axes(1:3), Time, &
          'total h2o -- all phases', 'kg/kg', missing_value=missing_value)

        id_tot_vapor     = register_diag_field ( mod_name, &
          'tot_vapor', axes(1:3), Time, &
          'total vapor', 'kg/kg', missing_value=missing_value)

        id_tot_liq_amt = register_diag_field ( mod_name, &
          'tot_liq_amt', axes(1:3), Time, &
          'Liquid amount -- all clouds', 'kg/kg',    &
                                              missing_value=missing_value)

        ID_clw = register_cmip_diag_field_3d ( mod_name, 'clw', Time, &
          'Mass Fraction of Cloud Liquid Water', 'kg kg-1',   &
          standard_name='mass_fraction_of_cloud_liquid_water_in_air' )

        id_tot_ice_amt = register_diag_field ( mod_name, &
          'tot_ice_amt', axes(1:3), Time, &
          'Ice amount -- all clouds', 'kg/kg',  &
                                            missing_value=missing_value )

        ID_cli = register_cmip_diag_field_3d ( mod_name, 'cli', Time, &
          'Mass Fraction of Cloud Ice', 'kg kg-1',   &
          standard_name='mass_fraction_of_cloud_ice_in_air' )

        id_lsc_cloud_area = register_diag_field ( mod_name, &
          'lsc_cloud_area', axes(1:3), Time, &
          'Large-scale cloud area', 'percent',    &
                                            missing_value=missing_value )

        id_lsc_liq_amt = register_diag_field ( mod_name, &
          'lsc_liq_amt', axes(1:3), Time, &
           'Large-scale cloud liquid amount', 'kg/kg',   &
                                            missing_value=missing_value )

        id_lsc_ice_amt = register_diag_field ( mod_name, &
          'lsc_ice_amt', axes(1:3), Time, &
          'Large-scale cloud ice amount', 'kg/kg',    &
                                            missing_value=missing_value )

        id_conv_cloud_area = register_diag_field ( mod_name, &
          'conv_cloud_area', axes(1:3), Time, &
          'Convective cloud area', 'percent', missing_value=missing_value )

        id_conv_liq_amt = register_diag_field ( mod_name, &
          'conv_liq_amt', axes(1:3), Time, &
          'Convective cloud liquid amount', 'kg/kg',    &
                                             missing_value=missing_value )

        id_conv_ice_amt = register_diag_field ( mod_name, &
          'conv_ice_amt', axes(1:3), Time, &
          'Convective cloud ice amount', 'kg/kg',    &
                                             missing_value=missing_value)

        id_WP_all_clouds = register_diag_field ( mod_name, &
          'WP_all_clouds', axes(1:2), Time, &
          'Total  water path -- all clouds + ls precip',        'kg/m2'   )

        id_clwvi = register_cmip_diag_field_2d ( mod_name, 'clwvi', Time, &
                                'Condensed Water Path', 'kg m-2', &
         standard_name='atmosphere_cloud_condensed_water_content' )

        id_LWP_all_clouds = register_diag_field ( mod_name, &
          'LWP_all_clouds', axes(1:2), Time, &
          'Liquid water path -- all clouds + ls rain',          'kg/m2'   )

        id_IWP_all_clouds = register_diag_field ( mod_name, &
          'IWP_all_clouds', axes(1:2), Time, &
          'Ice water path -- all clouds + ls snow',           'kg/m2'   )

        id_clivi = register_cmip_diag_field_2d ( mod_name, 'clivi', Time, &
                                      'Ice Water Path', 'kg m-2', &
                     standard_name='atmosphere_mass_content_of_cloud_ice' )

      endif

      id_rh = register_diag_field ( mod_name, &
        'rh', axes(1:3), Time, &
         'relative humidity',                            'percent',  &
                                             missing_value=missing_value )

      id_rh_cmip = register_diag_field ( mod_name, &
        'rh_cmip', axes(1:3), Time, &
        'relative humidity',                            'percent',  &
                             missing_value=missing_value               )

      ID_hur = register_cmip_diag_field_3d ( mod_name, 'hur', Time, &
        'Relative Humidity',                            '%',  &
         standard_name='relative_humidity' )

      id_qs = register_diag_field ( mod_name, &
        'qs', axes(1:3), Time, &
             'saturation specific humidity',                 'kg/kg',    &
                               missing_value=missing_value               )

!---------------------------------------------------------------------
!    register the diagnostics associated with convective tracer
!    transport.
!---------------------------------------------------------------------
      id_wetdep_om = register_diag_field ( mod_name, &
        'om_wet_dep',  axes(1:2), Time,  &
        'total om wet deposition', 'kg/m2/s',  missing_value=missing_value)
      if (id_wetdep_om > 0) wetdep_diagnostics_desired = .true.

      id_wetdep_SOA = register_diag_field ( mod_name, &
        'SOA_wet_dep',  axes(1:2), Time,  &
        'total SOA wet deposition', 'kg/m2/s', missing_value=missing_value)
      if (id_wetdep_SOA > 0) wetdep_diagnostics_desired = .true.

      id_wetdep_bc = register_diag_field ( mod_name, &
        'bc_wet_dep',  axes(1:2), Time,  &
        'total bc wet deposition', 'kg/m2/s',  missing_value=missing_value)
      if (id_wetdep_bc > 0) wetdep_diagnostics_desired = .true.

      id_wetdep_so4 = register_diag_field ( mod_name, &
        'so4_wet_dep',  axes(1:2), Time,  &
        'total so4 wet deposition', 'kg/m2/s', missing_value=missing_value)
      if (id_wetdep_so4 > 0) wetdep_diagnostics_desired = .true.

      id_wetdep_so2 = register_diag_field ( mod_name, &
        'so2_wet_dep',  axes(1:2), Time,  &
        'total so2 wet deposition', 'kg/m2/s', missing_value=missing_value)
      if (id_wetdep_so2 > 0) wetdep_diagnostics_desired = .true.

      id_wetdep_DMS = register_diag_field ( mod_name, &
        'DMS_wet_dep',  axes(1:2), Time,  &
        'total DMS wet deposition', 'kg/m2/s', missing_value=missing_value)
      if (id_wetdep_DMS > 0) wetdep_diagnostics_desired = .true.

      if (nNH4NO3 .ne. NO_TRACER .and. nNH4 .ne. NO_TRACER) then
        id_wetdep_NH4NO3 =  register_diag_field ( mod_name, &
          'totNH4_wet_dep',  axes(1:2), Time,  &
          'total NH4 + NH3 wet deposition', 'kg/m2/s',  &
                                               missing_value=missing_value)
      else
        id_wetdep_NH4NO3 = 0
      endif
      if (id_wetdep_NH4NO3 > 0) wetdep_diagnostics_desired = .true.

      id_wetdep_seasalt   =  register_diag_field ( mod_name, &
        'ssalt_wet_dep',  axes(1:2), Time,  &
        'total seasalt wet deposition', 'kg/m2/s',  &
                                               missing_value=missing_value)
      if (id_wetdep_seasalt > 0) wetdep_diagnostics_desired = .true.

      id_wetdep_dust   =  register_diag_field ( mod_name, &
        'dust_wet_dep',  axes(1:2), Time,  &
        'total dust wet deposition', 'kg/m2/s',   &
                                              missing_value=missing_value)
      if (id_wetdep_dust > 0) wetdep_diagnostics_desired = .true.


     !-------- cmip wet deposition fields  ---------
      do ic = 1, size(cmip_names,1)
        if (TRIM(cmip_names(ic)) .eq. 'nh4' .and. (nNH4NO3 .eq. NO_TRACER .or. nNH4 .eq. NO_TRACER)) then
          id_wetnh4_cmip = 0; cycle  ! skip when tracers are not in field table
        endif

        id_wetdep_cmip = register_cmip_diag_field_2d ( mod_name, 'wet'//TRIM(cmip_names(ic)), Time,  &
                                  'Wet Deposition Rate of '//TRIM(cmip_longnames(ic)), 'kg m-2 s-1', &
                    standard_name='tendency_of_atmosphere_mass_content_of_'//TRIM(cmip_stdnames(ic))//'_particles_due_to_wet_deposition' )
        if (TRIM(cmip_names(ic)) .eq. 'poa'  ) id_wetpoa_cmip  = id_wetdep_cmip
        if (TRIM(cmip_names(ic)) .eq. 'soa'  ) id_wetsoa_cmip  = id_wetdep_cmip
        if (TRIM(cmip_names(ic)) .eq. 'bc'   ) id_wetbc_cmip   = id_wetdep_cmip
        if (TRIM(cmip_names(ic)) .eq. 'dust' ) id_wetdust_cmip = id_wetdep_cmip
        if (TRIM(cmip_names(ic)) .eq. 'ss'   ) id_wetss_cmip   = id_wetdep_cmip
        if (TRIM(cmip_names(ic)) .eq. 'so4'  ) id_wetso4_cmip  = id_wetdep_cmip
        if (TRIM(cmip_names(ic)) .eq. 'so2'  ) id_wetso2_cmip  = id_wetdep_cmip
        if (TRIM(cmip_names(ic)) .eq. 'dms'  ) id_wetdms_cmip  = id_wetdep_cmip
        if (TRIM(cmip_names(ic)) .eq. 'nh4'  ) id_wetnh4_cmip  = id_wetdep_cmip
        if (id_wetdep_cmip > 0) wetdep_diagnostics_desired = .true.
      enddo
     !-----------------------------------------------

!------------------------------------------------------------------------
!    define wet dep diagnostic for each tracer. units and conversion factor
!    will differ dependent on tracer type (vmr vs. mmr).
!------------------------------------------------------------------------
      allocate (id_wetdep          (num_prog_tracers))
      allocate (id_wetdep_uw       (num_prog_tracers))
      allocate (id_wetdep_donner   (num_prog_tracers))
      allocate (id_wetdepc_donner  (num_prog_tracers))
      allocate (id_wetdepm_donner  (num_prog_tracers))
      allocate (id_wetdep_kg_m2_s  (num_prog_tracers))
      allocate (conv_wetdep        (num_prog_tracers))
      allocate (conv_wetdep_kg_m2_s(num_prog_tracers))
      allocate (nb_N_red           (num_prog_tracers))
      allocate (nb_N_ox            (num_prog_tracers))
      allocate (nb_N               (num_prog_tracers))
      id_wetdep         = -1
      id_wetdep_uw      = -1
      id_wetdep_donner  = -1
      id_wetdepc_donner = -1
      id_wetdepm_donner = -1
      id_wetdep_kg_m2_s = -1

      outunit = stdout()
      do n = 1, num_prog_tracers
        call get_tracer_names (MODEL_ATMOS, n, name = tracer_name, units=tracer_units)
        call  get_cmip_param (n, cmip_name=cmip_name, cmip_longname=cmip_longname, cmip_longname2=cmip_longname2)
        call  get_chem_param (n, mw=tracer_mw, is_aerosol=cmip_is_aerosol, nb_N=nb_N(n), nb_N_Ox=nb_N_Ox(n), nb_N_red=nb_N_red(n))

        write(outunit,'(a,g14.6)') trim(tracer_name)//', tracer_mw=',tracer_mw

        if (cmip_is_aerosol) then
            id_wetdep_kg_m2_s(n) = register_cmip_diag_field_2d ( mod_name, &
                               trim(tracer_name)//'_wetdep_kg_m2_s', Time, &
                               'Wet Deposition Rate of '//TRIM(cmip_longname2), 'kg m-2 s-1', &
                   standard_name='tendency_of_atmosphere_mass_content_of_'//TRIM(cmip_name)//'_dry_aerosol_due_to_wet_deposition')
        else
            id_wetdep_kg_m2_s(n) = register_cmip_diag_field_2d ( mod_name, &
                               trim(tracer_name)//'_wetdep_kg_m2_s', Time, &
                               'Wet Deposition Rate of '//TRIM(cmip_longname2), 'kg m-2 s-1', &
                           standard_name='tendency_of_atmosphere_mass_content_of_'//TRIM(cmip_name)//'_due_to_wet_deposition')
        end if
        if (id_wetdep_kg_m2_s(n) > 0) wetdep_diagnostics_desired = .true.

        if (id_wetdep_kg_m2_s(n) > 0) then
          if (tracer_mw < 0.0) then
            call error_mesg ('moist_processes', 'mw needs to be defined for tracer: '//trim(tracer_name), FATAL)
         !else
         !  write(outunit,'(a,g14.6)') trim(tracer_name)//', tracer_mw=',tracer_mw
          end if
        end if

        diaglname = trim(tracer_name)//  &
                        ' wet deposition from all precip'

        diaglname_uw = trim(tracer_name)//  &
                        ' wet deposition from uw'

        diaglname_donner = trim(tracer_name)//  &
                        ' wet deposition from donner'

        if ( tracer_units .eq. "vmr" ) then
          id_wetdep(n) = register_diag_field ( mod_name, &
              TRIM(tracer_name)//'_wet_depo',  &
              axes(1:2), Time, trim(diaglname), &
              'mole/m2/s',  missing_value=missing_value)

          id_wetdep_uw(n) = &
                       register_diag_field ( mod_name, &
                       TRIM(tracer_name)//'_wet_depo_uw',  &
                       axes(1:2), Time, trim(diaglname_uw), &
                       'mole/m2/s',  &
                       missing_value=missing_value)

          id_wetdep_donner(n) = &
                       register_diag_field ( mod_name, &
                       TRIM(tracer_name)//'_wet_depo_donner',  &
                       axes(1:2), Time, trim(diaglname_donner), &
                       'mole/m2/s',  &
                       missing_value=missing_value)

          id_wetdepm_donner(n) = &
                       register_diag_field ( mod_name, &
                       TRIM(tracer_name)//'_wet_depo_m_donner',  &
                       axes(1:2), Time, trim(diaglname_donner)//" meso", &
                       'mole/m2/s',  &
                       missing_value=missing_value)

          id_wetdepc_donner(n) = &
                       register_diag_field ( mod_name, &
                       TRIM(tracer_name)//'_wet_depo_c_donner',  &
                       axes(1:2), Time, trim(diaglname_donner)//" deep ",&
                       'mole/m2/s',  &
                       missing_value=missing_value)

          conv_wetdep(n) = 1d3/WTMAIR
          conv_wetdep_kg_m2_s(n) = tracer_mw*1e-3  ! std units are mol/m2/s


        elseif ( tracer_units.eq. "mmr" ) then
          id_wetdep(n) = &
                          register_diag_field ( mod_name, &
                          TRIM(tracer_name)//'_wet_depo',  &
                          axes(1:2), Time, trim(diaglname), &
                          'kg/m2/s',  &
                          missing_value=missing_value)

          id_wetdep_uw(n) = &
                       register_diag_field ( mod_name, &
                       TRIM(tracer_name)//'_wet_depo_uw',  &
                       axes(1:2), Time, trim(diaglname_uw), &
                       'kg/m2/s',  &
                       missing_value=missing_value)

          id_wetdep_donner(n) = &
                       register_diag_field ( mod_name, &
                       TRIM(tracer_name)//'_wet_depo_donner',  &
                       axes(1:2), Time, trim(diaglname_donner), &
                       'kg/m2/s',  &
                       missing_value=missing_value)

          id_wetdepm_donner(n) = &
                       register_diag_field ( mod_name, &
                       TRIM(tracer_name)//'_wet_depo_m_donner',  &
                       axes(1:2), Time, trim(diaglname_donner)//" meso", &
                       'kg/m2/s',  &
                       missing_value=missing_value)

          id_wetdepc_donner(n) = &
                      register_diag_field ( mod_name, &
                      TRIM(tracer_name)//'_wet_depo_c_donner',  &
                      axes(1:2), Time, trim(diaglname_donner)//" deep " , &
                      'kg/m2/s',  &
                      missing_value=missing_value)

        conv_wetdep(n) = 1.
        conv_wetdep_kg_m2_s(n) = 1. ! no conversion needed

        else
          write(outunit,'(a)') 'unsupported tracer: '//trim(tracer_name)//', units='//trim(tracer_units)
          conv_wetdep(n) = 0.
          conv_wetdep_kg_m2_s(n) = 0.
        end if

        if (id_wetdep(n) > 0) wetdep_diagnostics_desired = .true.

      end do

!---------------------------------------------------------------------

end subroutine diag_field_init


!#######################################################################

!#######################################################################
! <SUBROUTINE NAME="moist_processes_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments:
!   timestamp (optional, intent(in)) : A character string that represents the model time,
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix.
! </DESCRIPTION>
!
subroutine moist_processes_restart(timestamp)
  character(len=*), intent(in), optional :: timestamp

! if (doing_prog_clouds)       call strat_cloud_restart(timestamp)
! if (do_diag_clouds) call diag_cloud_restart(timestamp)
!  if (do_donner_deep) call donner_deep_restart(timestamp)
  call convection_driver_restart (timestamp)

end subroutine moist_processes_restart
! </SUBROUTINE> NAME="moist_processes_restart"


!#######################################################################

subroutine height_adjust(t, qv, r, tn, qvn, rn, &
                         phalf, pfull, zhalf, zfull, d_zhalf, d_zfull)

  real, intent(in),    dimension(:,:,:)   :: t, qv, tn, qvn
  real, intent(in),    dimension(:,:,:,:) :: r, rn
  real, intent(in),    dimension(:,:,:)   :: phalf, pfull, zhalf, zfull
  real, intent(out),   dimension(:,:,:)   :: d_zhalf, d_zfull


  integer :: i, j, k, ix, jx, kx, nql, nqi, nqa
  real, dimension(size(t,1), size(t,2), size(t,3))   :: tv, tvn, dz, dz_n
  real, dimension(size(t,1), size(t,2), size(t,3))   :: zfull_n
  real, dimension(size(t,1), size(t,2), size(t,3)+1) :: zhalf_n

      tv=0.; tvn=0.; dz=0.; dz_n=0; zhalf_n=0.; zfull_n=0.;
      ix = size(t,1)
      jx = size(t,2)
      kx = size(t,3)
      nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
      nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )

      tv  = t(:,:,:) *(1 + 0.608*qv (:,:,:) - r (:,:,:,nql)-r (:,:,:,nqi))
      tvn = tn(:,:,:)*(1 + 0.608*qvn(:,:,:) - rn(:,:,:,nql)-rn(:,:,:,nqi))

      zhalf_n(:,:,kx+1)=zhalf(:,:,kx+1);
      do k=kx,1,-1
         dz     (:,:,k)=zhalf(:,:,k)-zhalf(:,:,k+1)
         dz_n   (:,:,k)=dz(:,:,k)*tvn(:,:,k)/tv(:,:,k)
         zhalf_n(:,:,k)=zhalf_n(:,:,k+1)+dz_n(:,:,k)
         zfull_n(:,:,k)=(zhalf_n(:,:,k)+zhalf_n(:,:,k+1))*0.5
      enddo

      d_zhalf(:,:,:)=zhalf_n(:,:,:)-zhalf(:,:,:)
      d_zfull(:,:,:)=zfull_n(:,:,:)-zfull(:,:,:)

!      do k=1,kx
!         dlp (:,:,k)=log(phalf(:,:,k+1))-log(phalf(:,:,k))
!        tmp1(:,:,k)=RDGAS/GRAV*tv (:,:,k)*dlp(:,:,k)
!        tmp2(:,:,k)=RDGAS/GRAV*tvn(:,:,k)*dlp(:,:,k)
!        tv    (:,:,k)=zhalf_n(:,:,k)-zhalf(:,:,k)
!        tvn   (:,:,k)=zfull_n(:,:,k)-zfull(:,:,k)
!        tmp   (:,:,k)=(zhalf(:,:,k+1)+zhalf(:,:,k))*0.5-zfull(:,:,k)
!      enddo

end subroutine height_adjust


!#########################################################################

                        end module moist_processes_mod
