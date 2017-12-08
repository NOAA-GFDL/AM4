                    module lscloud_driver_mod

!-----------------------------------------------------------------------
!
!         interface module for large-scale moisture processes
!         ---------------------------------------
!             calls legacy strat_cloud code if desired
!             determines aerosol available for condensation
!             enforces realizability conditions
!             calls bulk large-scale condensation routine if requested
!             calls cloud macrophysics driver
!             calls cloud microphysics driver
!             computes tracer removal in ls clouds by wet deposition  
!             outputs large-scale cloud diagnostics and budget terms
!
!-----------------------------------------------------------------------

! fms modules
use time_manager_mod,      only: time_type
use diag_manager_mod,      only: register_diag_field, send_data
use mpp_mod,               only: input_nml_file
use fms_mod,               only: error_mesg, FATAL, NOTE,        &
                                 file_exist, check_nml_error,    &
                                 open_namelist_file, close_file, &
                                 write_version_number,           &
                                 open_file, &
                                 stdout, open_ieee32_file, &
                                 mpp_pe, mpp_root_pe, stdlog,    &
                                 mpp_clock_id, mpp_clock_begin,  &
                                 mpp_clock_end, CLOCK_MODULE,    &
                                 CLOCK_MODULE_DRIVER, &
                                 MPP_CLOCK_SYNC
use field_manager_mod,     only: MODEL_ATMOS
use tracer_manager_mod,    only: get_tracer_names, get_tracer_index, &
                                 NO_TRACER
use constants_mod,         only: CP_AIR, GRAV, HLV, HLS, HLF, &
                                 RDGAS, RVGAS
use physics_types_mod,     only: physics_control_type

! lscloud_driver modules
use lscale_cond_mod,       only: lscale_cond_init, lscale_cond_end, &
                                 lscale_cond
use lscloud_debug_mod,     only: lscloud_debug_init, lscloud_debug, &
                                 write_debug_output, output_refusals,  &
                                 lscloud_debug_setup
use lscloud_types_mod,     only: diag_id_type, diag_pt_type, &
                                 lscloud_types_init, &
                                 lscloud_nml_type, & 
                                 lsc_constants_type, lscloud_debug_type,&
                                 atmos_state_type, &
                                 particles_type, cloud_state_type, &
                                 cloud_processes_type, precip_state_type
use polysvp_mod,           only: compute_qs_a, polysvp_init, &
                                 polysvp_end
use lscloud_netcdf_mod,    only: lscloud_netcdf, lscloud_netcdf_init, &
                                 lscloud_netcdf_end

! other atmos_param modules
use rh_clouds_mod,         only: rh_clouds_init, rh_clouds_sum, &
                                 rh_clouds_end
use physics_radiation_exch_mod,        &
                           only: cloud_scheme_data_type,   &
                                 exchange_control_type
use ls_cloud_macrophysics_mod,      &
                           only: ls_cloud_macrophysics_init,  &
                                 ls_cloud_macrophysics_time_vary, &
                                 ls_cloud_macrophysics, &
                                 ls_cloud_macrophysics_end 
use ls_cloud_microphysics_mod,      &
                           only: ls_cloud_microphysics_init, &
                                 ls_cloud_microphysics_time_vary, &
                                 ls_cloud_microphysics,   &
                                 ls_cloud_microphysics_end 
use aerosol_cloud_mod,     only: aerosol_cloud_init, &
                                 determine_available_aerosol, &
                                 aerosol_cloud_end 
use strat_cloud_mod,       only: strat_cloud_init, strat_cloud_end, &
                                 strat_cloud
use moist_proc_utils_mod,  only: column_diag, rh_calc, & 
                                 mp_nml_type, mp_input_type, &
                                 mp_removal_type, mp_tendency_type,  &
                                 mp_conv2ls_type, mp_output_type,   &
                                 mp_lsdiag_type, mp_lsdiag_control_type

! atmos_shared modules
use aerosol_types_mod,     only: aerosol_type
use atmos_tracer_utilities_mod,       &
                           only: wet_deposition
use diag_axis_mod,         only: get_axis_num
use atmos_global_diag_mod, only: register_global_diag_field, &
                                 buffer_global_diag, &
                                 send_global_diag
use atmos_cmip_diag_mod,   only: register_cmip_diag_field_2d, &
                                 register_cmip_diag_field_3d, &
                                 send_cmip_data_3d, &
                                 cmip_diag_id_type, &
                                 query_cmip_diag_id
use diag_data_mod,         only: CMOR_MISSING_VALUE
implicit none
private

!-----------------------------------------------------------------------
!-------------------- public data/interfaces ---------------------------

public   lscloud_driver_init, lscloud_driver_time_vary,      &
         lscloud_driver, lscloud_driver_endts, lscloud_driver_end
  
private  diag_field_init, lscloud_alloc,                  &
         impose_realizability, impose_realizability_clubb, &
         detailed_diagnostics, update_fields_and_tendencies,&
         compute_ls_wetdep, basic_diagnostics, lscloud_dealloc

!-----------------------------------------------------------------------
!-------------------- private data -------------------------------------

!--------------------- version number ----------------------------------
character(len=128) :: &
version = '$Id: $'
character(len=128) :: tagname = '$Name: $'

!------------------------------------------------------------------------
!---namelist-------------------------------------------------------------


!-----------------------------------------------------------------------
!    do_legacy_strat_cloud 
!                   activate the older version of the strat_cloud module
!                   rather than the latest ? (default = .false)

!  <DATA NAME="Dmin" UNITS="dimensionless" TYPE="real"  DEFAULT="1.0e-07">
!   minimum permissible dissipation in analytic integration of qa, ql, qi
!   equations. This constant only affects the method by which the
!   prognostic equations are integrated.
!   NOTE: Dmin will be MACHINE DEPENDENT and occur when
!   a) 1. -exp(-Dmin) = 0. instead of Dmin in the limit of very small Dmin.
!   AND
!   b) 1. - exp(-D) < D for all D > Dmin
!  </DATA>
!  <DATA NAME="cfact" UNITS="" TYPE="real" DEFAULT="1.0">
!   factor in bergeron process calculation when the effect of dust
!   is not considered; values > 1 enhance the conversion to ice.
!  </DATA>
!  <DATA NAME="microphys_scheme" TYPE="character" DEFAULT="rotstayn_klein">
!   the microphysics scheme being used (currently either
!   "morrison_gettelman" or "rotstayn_klein" or "mg_ncar" or "ncar")
!  </DATA>
!  <DATA NAME="macrophys_scheme" TYPE="character" DEFAULT="tiedtke">
!   the macrophysics scheme being used (currently either
!   "tiedtke" or "           ")
!  </DATA>
!  <DATA NAME="aerosol_activation_scheme" TYPE="character" DEFAULT="dqa">
!   the aerosol activation scheme being used (currently either
!   "dqa" or "total")
!  </DATA>
!  <DATA NAME="do_dust_berg" TYPE="logical" DEFAULT=".false.">
!   sub-micron dust particles are used as ice nuclei for bergeron
!   process evaluation?
!  </DATA>
!  <DATA NAME="super_ice_opt" TYPE="integer" DEFAULT="0">
!   flag to indicate how to treat supersaturation; 0 => don't allow.
!  </DATA>
!  <DATA NAME="do_ice_nucl_wpdf" TYPE="logical" DEFAULT=".true.">
!   should we use activate ice nuclei by assuming a pdf ?
!  </DATA>
!  <DATA NAME="do_pdf_clouds" TYPE="logical" DEFAULT=".false.">
!   the statistical cloud scheme should be run?
!  </DATA>
!  <DATA NAME="betaP" UNITS="" TYPE="integer" DEFAULT="5">
!   p-parameter to the beta distribution - used when do_pdf_clouds is true
!  </DATA>
!  <DATA NAME="qthalfwidth" UNITS="" TYPE="real" DEFAULT="0.1">
!   half-width to the qt PDF - used when do_pdf_clouds is true and
!   diagnostic variance
!   The fraction of qtbar (mean total water in the grid box) that the
!   maximum and minimum of the distribution differ from qtbar. That is,
!   total water at the sub-grid scale may take on values anywhere between
!   (1.-qthalfwidth)*qtbar and (1.+qthalfwidth)*qtbar
!  </DATA>
!  <DATA NAME="nsublevels" UNITS="" TYPE="integer" DEFAULT="1">
!   number of sublevels to vertical sub-grid cloud structure - used when
!   do_pdf_cloud is true
!  </DATA>
!  <DATA NAME="kmap" UNITS="" TYPE="integer" DIM="" DEFAULT="1">
!   PPM partial remap integer - used when do_pdf_cloud is true and if
!   vertical subgrid structure is used
!  </DATA>
!  <DATA NAME="kord" UNITS="" TYPE="integer"  DEFAULT="7">
!   PPM method number - used when do_pdf_cloud is true and if vertical
!   subgrid structure is used
!  </DATA>
!  <DATA NAME="pdf_org" TYPE="logical" DEFAULT=".true.">
!   define pdf based on input water fields (rather than updated ones) ?
!  </DATA>
!-----------------------------------------------------------------------



logical            :: do_legacy_strat_cloud = .false.
real               :: Dmin = 1.0e-7
real               :: cfact = 1.0
character(len=64)  :: microphys_scheme = 'rotstayn_klein'
character(len=64)  :: macrophys_scheme = 'tiedtke'
character(len=64)  :: aerosol_activation_scheme = 'dqa'
logical            :: do_dust_berg = .false.
integer            :: super_ice_opt = 0
logical            :: do_ice_nucl_wpdf = .true.
logical            :: do_pdf_clouds = .false.
integer            :: betaP = 5
real               :: qthalfwidth = 0.1
integer            :: nsublevels = 1
integer            :: kmap = 1
integer            :: kord = 7
logical            :: pdf_org = .true.
logical :: use_cf_metadata = .false.


namelist / lscloud_driver_nml / do_legacy_strat_cloud, Dmin, cfact, &
                                microphys_scheme, macrophys_scheme, &
                                aerosol_activation_scheme, &
                                do_dust_berg, &
                                super_ice_opt, do_ice_nucl_wpdf, &
                                do_pdf_clouds, betaP, qthalfwidth, &
                                nsublevels, kmap, kord, pdf_org, &
                                use_cf_metadata



!------------------------------------------------------------------------
!     ------ constants used in the module -------
!-----------------------------------------------------------------------
real, parameter :: d608 = (RVGAS - RDGAS)/RDGAS


!-------------------- clock definitions --------------------------------

integer :: lscloud_driver_clock, stratcloud_clock, ls_macrophysics_clock,&
           ls_microphysics_clock, lscalecond_clock, &
           polysvp_clock, lscloud_debug_clock, aerosol_cloud_clock, &
           lscloud_driver_init_clock, lscloud_driver_term_clock, &
           lscloud_netcdf_clock

integer :: lscloud_alloc_clock, realiz_clock, realiz_clubb_clock, &
           detail_diag_clock, update_fields_clock, ls_wetdep_clock, &
           basic_diags_clock, dealloc_clock, adjust_cond_clock, &
           adjust_part_num_clock
            

!-------------------- diagnostics fields -------------------------------

integer :: id_LWP, id_IWP, id_tdt_ls, id_qdt_ls, id_prec_ls, id_snow_ls, &
           id_qadt_ls, id_qldt_ls, id_qndt_ls, id_qidt_ls, id_qnidt_ls, &
           id_qa_ls_col, id_ql_ls_col, id_qn_ls_col, id_qi_ls_col, &
           id_qni_ls_col, id_q_ls_col, id_t_ls_col, &
           id_enth_ls_col, id_wat_ls_col, &
           id_lsc_precip, id_lsc_freq,                           &
           id_lscale_rain3d, id_lscale_snow3d, id_lscale_precip3d
 
integer :: id_qvout, id_qaout, id_qlout, id_qiout
integer :: id_qnout, id_qniout
integer :: id_f_snow_berg, id_f_snow_berg_cond, id_f_snow_berg_wtd

integer, dimension(:), allocatable :: id_wet_deposition    

real :: missing_value = -999.
character(len=5), private :: mod_name = 'moist'

!------------------- other global variables and parameters -------------

integer, parameter :: n_totmass  = 4
integer, parameter :: n_imass = 12

!type(mp_lsdiag_type)       :: Lsdiag_mp
type(mp_lsdiag_control_type) :: Lsdiag_mp_control
type(lsc_constants_type)   :: Constants_lsc
type(lscloud_nml_type)     :: Nml_lsc

real     :: N_land, N_ocean, qcvar
logical  :: do_liq_num
real     :: qmin
logical  :: do_ice_num
integer  :: do_clubb
logical  :: limit_conv_cloud_frac
logical  :: do_rh_clouds
logical  :: module_is_initialized = .false.
logical  :: doing_prog_clouds
logical  :: do_lsc
logical  :: do_simple

real     :: dtcloud, inv_dtcloud
logical  :: do_predicted_ice_number
integer  :: nsphum, nql, nqi, nqa, nqn, nqni, nqr, nqs, nqg
integer  :: nso2, nso4
integer  :: num_prog_tracers
logical  :: debug

type(cmip_diag_id_type) :: ID_tntscp, ID_tnhusscp


                             contains


!#######################################################################

subroutine lscloud_driver_init (id, jd, kd, axes, Time, &
                                Exch_ctrl, Nml_mp, Physics_control, &
                                lon, lat, phalf, pref)

integer,                 intent(in)     :: id, jd, kd
integer,                 intent(in)     :: axes(4)
type(time_type),         intent(in)     :: Time
type(exchange_control_type), intent(in) :: Exch_ctrl
type(mp_nml_type),       intent(inout)  :: Nml_mp
type(physics_control_type), intent(in)  :: Physics_control
real,dimension(:,:),     intent(in)     :: lon,  lat    ! h1g
real,dimension(:,:,:),   intent(in)     :: phalf        ! h1g
real, dimension(:),      intent(in)     :: pref


! --- internal variables ---
      integer            :: unit, io, ierr, logunit
      character(len=128) :: errstring     ! Output status: non-blank for 
                                          ! error return
      integer            :: lscalecond_init_clock, strat_init_clock, &
                            debug_init_clock, polysvp_init_clock, &
                            lscloud_types_init_clock,   &
                            aerosol_cloud_init_clock,  &
                            macrophysics_init_clock, &
                            microphysics_init_clock, &
                            lscloud_netcdf_init_clock, &
                            diag_field_init_clock
      integer            :: k

!------------------------------------------------------------------------
!   if module is already initialized, return.
!------------------------------------------------------------------------
      if (module_is_initialized) return

!------------------------------------------------------------------------
!   initialize clock for this module's initialization and then activate.
!------------------------------------------------------------------------
      lscloud_driver_init_clock = mpp_clock_id(   &
                   '  Lscloud_driver: Initialization' , &
                                               grain=CLOCK_MODULE_DRIVER )
      call mpp_clock_begin (lscloud_driver_init_clock)

!------------------------------------------------------------------------
!    extract variables needed later from derived type inputs.
!------------------------------------------------------------------------
      qmin = Exch_ctrl%qmin
      N_land = Exch_ctrl%N_land
      N_ocean = Exch_ctrl%N_ocean
      qcvar = Exch_ctrl%qcvar
      do_liq_num = Exch_ctrl%do_liq_num
      do_ice_num = Exch_ctrl%do_ice_num
      do_clubb   = Exch_ctrl%do_clubb  

      limit_conv_cloud_frac = Nml_mp%limit_conv_cloud_frac
      do_rh_clouds = Nml_mp%do_rh_clouds
      do_lsc = Nml_mp%do_lsc
      do_simple = Nml_mp%do_simple
      doing_prog_clouds = Exch_ctrl%doing_prog_clouds

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

!-------------------------------------------------------------------------
!    initialize module clocks.
!-------------------------------------------------------------------------
      lscalecond_init_clock = mpp_clock_id(   &
                   'Lscloud_driver: lscale_cond:Initialization' , &
                                               grain=CLOCK_MODULE_DRIVER )
      strat_init_clock = mpp_clock_id(   &
                   'Lscloud_driver: strat_cloud:Initialization' , &
                                               grain=CLOCK_MODULE_DRIVER )
      debug_init_clock = mpp_clock_id(   &
                   'Lscloud_driver: lscloud_debug:Initialization' , &
                                               grain=CLOCK_MODULE_DRIVER )
      polysvp_init_clock = mpp_clock_id(   &
                   'Lscloud_driver: polysvp:Initialization' , &
                                               grain=CLOCK_MODULE_DRIVER )
      lscloud_types_init_clock = mpp_clock_id(   &
                   'Lscloud_driver: lscloud_types:Initialization' , &
                                               grain=CLOCK_MODULE_DRIVER )
      aerosol_cloud_init_clock = mpp_clock_id(   &
                   'Lscloud_driver: aerosol_cloud:Initialization' , &
                                               grain=CLOCK_MODULE_DRIVER )
      macrophysics_init_clock = mpp_clock_id(   &
                   'Lscloud_driver: ls_cloud_macrophysics:Initialization',&
                                               grain=CLOCK_MODULE_DRIVER )
      microphysics_init_clock = mpp_clock_id(   &
                   'Lscloud_driver: ls_cloud_microphysics:Initialization',&
                                               grain=CLOCK_MODULE_DRIVER )
      lscloud_netcdf_init_clock = mpp_clock_id(   &
                   'Lscloud_driver: lscloud_netcdf:Initialization' , &
                                               grain=CLOCK_MODULE_DRIVER )
      lscloud_driver_clock = mpp_clock_id( '   Lscloud_driver: '   ,&
                                                grain=CLOCK_MODULE_DRIVER )
      lscalecond_clock = mpp_clock_id(   &
                   '   Lscloud_driver: lscale_cond' , &
                                              grain=CLOCK_MODULE_DRIVER )
      stratcloud_clock = mpp_clock_id  &
                              ( '   Lscloud_driver: strat_cloud' ,&
                                                grain=CLOCK_MODULE_DRIVER )
      polysvp_clock = mpp_clock_id  &
                              ( '   Lscloud_driver: polysvp' ,&
                                                grain=CLOCK_MODULE_DRIVER )
      lscloud_debug_clock = mpp_clock_id  &
                              ( '   Lscloud_driver: lscloud_debug' ,&
                                                grain=CLOCK_MODULE_DRIVER )
      aerosol_cloud_clock = mpp_clock_id  &
                              ( '   Lscloud_driver: aerosol_cloud' ,&
                                                grain=CLOCK_MODULE_DRIVER )
      ls_macrophysics_clock = mpp_clock_id   &
                               ( 'Lscloud_driver: ls_cloud_macrophysics' ,&
                                               grain=CLOCK_MODULE_DRIVER )
      ls_microphysics_clock = mpp_clock_id   &
                               ( 'Lscloud_driver: ls_cloud_microphysics' ,&
                                               grain=CLOCK_MODULE_DRIVER )
      lscloud_netcdf_clock = mpp_clock_id   &
                               ( 'Lscloud_driver: lscloud_netcdf' ,&
                                               grain=CLOCK_MODULE_DRIVER )
      lscloud_driver_term_clock = mpp_clock_id   &
                               ( '   Lscloud_driver: Termination'   ,&
                                                grain=CLOCK_MODULE_DRIVER )

!------------------------------------------------------------------------
!    define clocks for the private routines of this module.
!------------------------------------------------------------------------
      lscloud_alloc_clock = mpp_clock_id(   &
                   'Lscloud_driver: lscloud_alloc  ', &
                                               grain=CLOCK_MODULE_DRIVER )
      realiz_clock = mpp_clock_id(   &
                   'Lscloud_driver: impose_real    ', &
                                               grain=CLOCK_MODULE_DRIVER )
      realiz_clubb_clock = mpp_clock_id(   &
                   'Lscloud_driver: impo_real_clubb', &
                                               grain=CLOCK_MODULE_DRIVER )
      detail_diag_clock =  mpp_clock_id(   &
                   'Lscloud_driver: detail_diag    ', &
                                               grain=CLOCK_MODULE_DRIVER )
      update_fields_clock = mpp_clock_id(   &
                   'Lscloud_driver: update_fields  ', &
                                               grain=CLOCK_MODULE_DRIVER )
      ls_wetdep_clock = mpp_clock_id(   &
                   'Lscloud_driver: ls_wetdep      ', &
                                               grain=CLOCK_MODULE_DRIVER )
      basic_diags_clock = mpp_clock_id(   &
                   'Lscloud_driver: basic_diags    ', &
                                               grain=CLOCK_MODULE_DRIVER )
      dealloc_clock = mpp_clock_id(   &
                   'Lscloud_driver: lscloud_dealloc', &
                                               grain=CLOCK_MODULE_DRIVER )
      adjust_cond_clock = mpp_clock_id(   &
                   'Lscloud_driver: adjust_condensate', &
                                               grain=CLOCK_MODULE_DRIVER )
      adjust_part_num_clock = mpp_clock_id(   &
                   'Lscloud_driver: adjust_part_num', &
                                               grain=CLOCK_MODULE_DRIVER )

!------------------------------------------------------------------------
!    if prognostic clouds are not active, most of this routine is not
!    executed.
!------------------------------------------------------------------------
      if (doing_prog_clouds) then

!----------------------------------------------------------------------- 
!    process namelist.
!-----------------------------------------------------------------------  
#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=lscloud_driver_nml, iostat=io)
        ierr = check_nml_error(io,'lscloud_driver_nml')
#else
        if ( file_exist('input.nml')) then
          unit = open_namelist_file ()
          ierr=1; do while (ierr /= 0)
          read  (unit, nml=lscloud_driver_nml, iostat=io, end=10)
          ierr = check_nml_error(io,'lscloud_driver_nml')
          enddo
10        call close_file (unit)
        endif
#endif

!----------------------------------------------------------------------- 
!    write version and namelist to stdlog.
!-----------------------------------------------------------------------  
        call write_version_number(Version, Tagname)
        logunit = stdlog()
        if ( mpp_pe() == mpp_root_pe() )  &
                    write (logunit, nml=lscloud_driver_nml)

!-----------------------------------------------------------------------
!    obtain indices for tracers no2 and no4.
!-----------------------------------------------------------------------
        nso2      = get_tracer_index(MODEL_ATMOS,'simpleSO2')
        if (nso2 == NO_TRACER) then
          nso2      = get_tracer_index(MODEL_ATMOS,'so2')
        endif
 
        nso4      = get_tracer_index(MODEL_ATMOS,'simpleSO4')
        if (nso4 == NO_TRACER) then
          nso4      = get_tracer_index(MODEL_ATMOS,'so4')
        endif

!------------------------------------------------------------------------
!    define logicals defining macrophysics scheme which is active.
!------------------------------------------------------------------------
        if (trim(macrophys_scheme) == 'tiedtke') then
          Constants_lsc%tiedtke_macrophysics = .true.
        else
          Constants_lsc%tiedtke_macrophysics = .false.
        endif

!-----------------------------------------------------------------------
!    prevent legacy strat cloud from being activated with clubb.
!-----------------------------------------------------------------------
        if (do_clubb > 0) then
          Constants_lsc%tiedtke_macrophysics = .false.
          do_legacy_strat_cloud = .false.
        endif
 
!-----------------------------------------------------------------------
!   perform realizability checks -- make sure one and only one large-scale
!   cloud scheme has been activated, and that no untested combinations
!   are being attempted.
!-----------------------------------------------------------------------
        if ( do_clubb > 0 .and. do_lsc) then
          call error_mesg ('lscloud_driver_mod', &
                 ' cannot do large-scale condensation when&
                                & CLUBB is active', FATAL)
        endif

        if ( Constants_lsc%tiedtke_macrophysics .and. do_lsc) then
          call error_mesg ('lscloud_driver_mod', &
                 ' cannot do large-scale condensation when&
                                & tiedtke_macrophysics is active', FATAL)
        endif

        if ( .not. Constants_lsc%tiedtke_macrophysics .and.   &
                                              do_legacy_strat_cloud) then
          call error_mesg ('lscloud_driver_mod', &
                 ' do_legacy_strat_cloud cannot be true when&
                                & tiedtke_macrophysics is false', FATAL)
        endif

        if (trim(microphys_scheme) == 'lin' .and. do_clubb > 0) then
          call error_mesg ('lscloud_driver_mod', &
                        'cannot run lin microphysics with CLUBB', FATAL)
        endif

        if (trim(microphys_scheme) == 'lin' .and. do_lsc) then
          call error_mesg ('lscloud_driver_mod', &
                        'cannot run lin microphysics with large-&
                                   &scale condensation', FATAL)
        endif

        if (do_clubb > 0 .and. .not. do_liq_num) then
          call error_mesg ('lscloud_driver_mod', &
              'can only execute clubb with prognostic droplet number', &
                                                                    FATAL)
        endif ! do_clubb

        if (trim(microphys_scheme) == 'lin') then

          if (Constants_lsc%tiedtke_macrophysics) then
            call error_mesg ('lscloud_driver_mod', &
              'cannot have both tiedtke_macrophysics and lin active &
                      &at once: setting tiedtke_macrophysics   &
                       &and do_legacy_strat_cloud = F', NOTE)
            Constants_lsc%tiedtke_macrophysics = .false.
            do_legacy_strat_cloud = .false.

            logunit = stdlog()
            if ( mpp_pe() == mpp_root_pe() ) &
              write ( logunit, '(a)')    &
               'variables tiedtke_macrophysics and   &
                     &do_legacy_strat_cloud set to .false,   &
                             &since lin microphysics has been activated.'
          endif ! tiedtke_macrophysics
        endif  ! trim(microphys_scheme) ==  'lin'

!-----------------------------------------------------------------------
!    check for acceptable namelist values:
!    qthalfwidth must be greater than 0.001
!    nsublevels must be greater than 0
!-----------------------------------------------------------------------
        if (qthalfwidth .lt. 1.e-03) then
          call error_mesg ( 'lscloud_driver_mod', &
                       'qthalfwidth must be greater than 0.001', FATAL)
        endif
        if (nsublevels .lt. 1) then
          call error_mesg ( 'lscloud_driver_mod', &
                            'nsublevels must be greater than 0', FATAL)
        endif

!-----------------------------------------------------------------------
!    put lscloud_driver_nml variables into a derived type variable for
!    passing to other modules as needed.
!-----------------------------------------------------------------------
        Nml_lsc%do_legacy_strat_cloud = do_legacy_strat_cloud
        Nml_lsc%Dmin = Dmin 
        Nml_lsc%cfact = cfact
        Nml_lsc%super_ice_opt = super_ice_opt
        Nml_lsc%do_ice_nucl_wpdf = do_ice_nucl_wpdf
        Nml_lsc%do_dust_berg = do_dust_berg
        Nml_lsc%do_pdf_clouds = do_pdf_clouds
        Nml_lsc%betaP = betaP
        Nml_lsc%qthalfwidth = qthalfwidth
        Nml_lsc%nsublevels = nsublevels
        Nml_lsc%kmap = kmap
        Nml_lsc%kord = kord
        Nml_lsc%pdf_org = pdf_org

!------------------------------------------------------------------------
!    define logicals defining active microphysics scheme and its type.
!------------------------------------------------------------------------
        if (trim(microphys_scheme) =='rotstayn_klein') then
          Constants_lsc%do_rk_microphys = .true.
          Constants_lsc%do_mg_microphys = .false.
          Constants_lsc%do_mg_ncar_microphys = .false.
          Constants_lsc%do_ncar_microphys = .false.
          do_predicted_ice_number = .false.
          Constants_lsc%do_lin_cld_microphys = .false.
        else if (trim(microphys_scheme) == 'morrison_gettelman') then
          Constants_lsc%do_rk_microphys = .false.
          Constants_lsc%do_mg_microphys = .true.
          Constants_lsc%do_mg_ncar_microphys = .false.
          Constants_lsc%do_ncar_microphys = .false.
          do_predicted_ice_number = .true.
          Constants_lsc%do_lin_cld_microphys = .false.
        else if (trim(microphys_scheme) == 'mg_ncar') then
          Constants_lsc%do_rk_microphys = .false.
          Constants_lsc%do_mg_microphys = .false.
          Constants_lsc%do_mg_ncar_microphys = .true.
          Constants_lsc%do_ncar_microphys = .false.
          do_predicted_ice_number = .true.
          Constants_lsc%do_lin_cld_microphys = .false.
        else if (trim(microphys_scheme) == 'ncar') then
          Constants_lsc%do_rk_microphys = .false.
          Constants_lsc%do_mg_microphys = .false.
          Constants_lsc%do_mg_ncar_microphys = .false.
          Constants_lsc%do_ncar_microphys = .true.
          do_predicted_ice_number = .true.
          Constants_lsc%do_lin_cld_microphys = .false.
        else if (trim(microphys_scheme) == 'lin') then
          Constants_lsc%do_rk_microphys = .false.
          Constants_lsc%do_mg_microphys = .false.
          Constants_lsc%do_mg_ncar_microphys = .false.
          Constants_lsc%do_ncar_microphys = .false.
          Constants_lsc%do_lin_cld_microphys = .true.
! this version of lin could not be active with prog drop number (and thus
!   with predicted ice number)
          do_predicted_ice_number = .false.
        else
          call error_mesg ('lscloud_driver_init', &
                'invalid expression supplied for nml variable &
                                           &microphys_scheme', FATAL)
        endif

!------------------------------------------------------------------------
!    define logicals defining aerosol activation scheme which is active.
!------------------------------------------------------------------------
        if (do_clubb > 0) then
          Constants_lsc%dqa_activation = .false.
          Constants_lsc%total_activation = .false.
        else
          if (trim(aerosol_activation_scheme) == 'dqa') then
            Constants_lsc%dqa_activation = .true.
            Constants_lsc%total_activation = .false.
          else if (trim(aerosol_activation_scheme) == 'total') then
            Constants_lsc%dqa_activation = .false.
            Constants_lsc%total_activation = .true.
          else
            call error_mesg ('lscloud_driver_init',   &
           'invalid value for aerosol_activation_scheme specified', FATAL)
          endif
        endif

!------------------------------------------------------------------------
!    call strat_cloud_init to do initialization there. Call is only needed
!    when executing legacy strat_cloud code.
!------------------------------------------------------------------------
        if (do_legacy_strat_cloud) then
          call mpp_clock_begin ( strat_init_clock)
          call strat_cloud_init (Nml_mp, Nml_lsc, Exch_ctrl,   &
                                                       Physics_control)
          call mpp_clock_end ( strat_init_clock)
        endif

!------------------------------------------------------------------------
!    call lscloud_debug_init to initialize the debug module.
!------------------------------------------------------------------------
        call mpp_clock_begin (debug_init_clock)
        call lscloud_debug_init (debug)
        call mpp_clock_end (debug_init_clock)

!------------------------------------------------------------------------
!    initialize peripheral routines used with large-scale clouds.
!------------------------------------------------------------------------
        if (Constants_lsc%tiedtke_macrophysics) then
          call mpp_clock_begin (polysvp_init_clock)
          call polysvp_init (Nml_lsc)
          call mpp_clock_end (polysvp_init_clock)
        endif

        call mpp_clock_begin (lscloud_types_init_clock)
        call lscloud_types_init
        call mpp_clock_end (lscloud_types_init_clock)

!------------------------------------------------------------------------
!    initialize aerosol module unless doing legacy strat_cloud. in that 
!    case, this functionality is included within the legacy strat_cloud
!    module.
!------------------------------------------------------------------------
        call mpp_clock_begin (aerosol_cloud_init_clock)
        call aerosol_cloud_init (Constants_lsc, Nml_lsc, Nml_mp, Exch_ctrl)
        call mpp_clock_end (aerosol_cloud_init_clock)

!----------------------------------------------------------------------
!    initialize the cloud macrophysics module.
!----------------------------------------------------------------------
        call mpp_clock_begin ( macrophysics_init_clock)
        call ls_cloud_macrophysics_init   &
                      (id, jd, kd, lon, lat, axes, Time, phalf, Nml_mp, &
                        Constants_lsc, Physics_control, Nml_lsc, Exch_ctrl)
        call mpp_clock_end ( macrophysics_init_clock)

!----------------------------------------------------------------------
!    initialize the cloud microphysics module.
!----------------------------------------------------------------------
        call mpp_clock_begin ( microphysics_init_clock)
        call ls_cloud_microphysics_init   &
                     (Nml_mp, Constants_lsc, Physics_control, id, jd,   &
                      kd, Time,axes, pref, Nml_lsc, Exch_ctrl)
        call mpp_clock_end ( microphysics_init_clock)

!------------------------------------------------------------------------
!    define various control variables for the non-prognostic clouds case.
!------------------------------------------------------------------------
      else ! (doing_prog_clouds)
        Constants_lsc%tiedtke_macrophysics = .false.
        Constants_lsc%do_rk_microphys = .false.
        Constants_lsc%do_mg_microphys = .false.
        Constants_lsc%do_mg_ncar_microphys = .false.
        Constants_lsc%do_ncar_microphys = .false.
        do_predicted_ice_number = .false.
        Constants_lsc%do_lin_cld_microphys = .false.
        Constants_lsc%dqa_activation = .false.
        Constants_lsc%total_activation = .false.

!------------------------------------------------------------------------ 
!    if large-scale bulk condensation is active, initialize it here,   
!    along with the cloud scheme to be used to obtain a cloud field for
!    the radiation code to work with.  
!------------------------------------------------------------------------
        if (do_lsc) then
          call mpp_clock_begin (lscalecond_init_clock)
          call lscale_cond_init ()

!----------------------------------------------------------------------
!    initialize the rh_clouds module, if needed.
!----------------------------------------------------------------------
          if (Nml_mp%do_rh_clouds) then
            call rh_clouds_init (id, jd, kd)
          endif
          call mpp_clock_end   (lscalecond_init_clock)
        endif
      endif ! (doing_prog_clouds)
                                               
!-----------------------------------------------------------------------
!    call lscloud_netcdf_init to set up the netcdf diagnostic output
!    requested via the diag_table. 
!-----------------------------------------------------------------------
      call mpp_clock_begin ( lscloud_netcdf_init_clock)
      call lscloud_netcdf_init (axes, Time, Lsdiag_mp_control%diag_id,  &
                                          Lsdiag_mp_control%diag_pt,   &
                                          Lsdiag_mp_control%n_diag_4d, &
                                          Lsdiag_mp_control%n_diag_4d_kp1)
      call mpp_clock_end ( lscloud_netcdf_init_clock)

!------------------------------------------------------------------------
!    initiate clock to time routine which initializes netcdf diagnostics.
!------------------------------------------------------------------------
      diag_field_init_clock = mpp_clock_id(   &
                   'Lscloud_driver: diag_field_init', &
                                               grain=CLOCK_MODULE_DRIVER )
      call mpp_clock_begin (diag_field_init_clock)

!-------------------------------------------------------------------------
!    call diag_field_init to register desired diagnostics.
!-------------------------------------------------------------------------
      call diag_field_init (axes, Time)
      call mpp_clock_end (diag_field_init_clock)

!-----------------------------------------------------------------------
!    mark the module as imnitialized and turn off the module timer.
!-----------------------------------------------------------------------
      module_is_initialized = .true.
      call mpp_clock_end (lscloud_driver_init_clock)

!------------------------------------------------------------------------

end subroutine lscloud_driver_init



!########################################################################

subroutine lscloud_driver_time_vary (dt)

real,                    intent(in) :: dt

!------------------------------------------------------------------------
!    if legacy strat_cloud code is being executed, there is no separate 
!    time-varying routine to set time-dependent, spatially constant
!    variables, so this routine is skipped.
!------------------------------------------------------------------------
      if (doing_prog_clouds .and. (.not. do_legacy_strat_cloud) ) then

!-----------------------------------------------------------------------
!    set the current time step (in case it varies with time).
!-----------------------------------------------------------------------
        dtcloud = dt
        inv_dtcloud = 1.0/dt

!-----------------------------------------------------------------------
!    call routines to set any time-dependent, spatially independent
!    variables used by the macrophysics and microphysics modules.
!    here the physics timestep is needed.
!-----------------------------------------------------------------------
        call mpp_clock_begin (ls_macrophysics_clock)
        call ls_cloud_macrophysics_time_vary (dtcloud)
        call mpp_clock_end (ls_macrophysics_clock)

        call mpp_clock_begin (ls_microphysics_clock)
        call ls_cloud_microphysics_time_vary (dtcloud)
        call mpp_clock_end (ls_microphysics_clock)
      endif
 
!-----------------------------------------------------------------------


end subroutine lscloud_driver_time_vary



!#########################################################################

subroutine lscloud_driver  (is, ie, js, je, Time, dt, Input_mp, &
                            rdiag, Tend_mp, C2ls_mp, Output_mp, &
                            Removal_mp, Cld_props, Aerosol)

!-----------------------------------------------------------------------
!
!    in:  is,ie      starting and ending i indices for window
!
!         js,je      starting and ending j indices for window
!
!         Time       time used for diagnostics [time_type]
!
!         dt         time step (from t(n-1) to t(n+1) if leapfrog)
!                    in seconds   [real]
!
!         Input_mp
!
!         rdiag
!  
!         Tend_mp
!
!         C2ls_mp
! 
!         Output_mp
!         
!         Removal_mp 
!
!         Cld_props
!
!         Aerosol
!
!
!-----------------------------------------------------------------------
integer,                     intent(in)           :: is,ie,js,je
type(time_type),             intent(in)           :: Time
real,                        intent(in)           :: dt
type(mp_input_type),         intent(inout)        :: Input_mp
type(mp_tendency_type),      intent(inout)        :: Tend_mp
type(mp_conv2ls_type),       intent(inout)        :: C2ls_mp
type(mp_output_type),        intent(inout)        :: Output_mp
type(mp_removal_type),       intent(inout)        :: Removal_mp
type(cloud_scheme_data_type),intent(inout)        :: Cld_props
real, dimension(:,:,:,size(Output_mp%rdt,4)+1:),      &
                             intent(inout)        :: rdiag
type(aerosol_type),          intent(in), optional :: Aerosol

!-----------------------------------------------------------------------
!   derived type variables used to pass fields between modules:
!------------------------------------------------------------------------
      type(atmos_state_type)     :: Atmos_state
      type(cloud_state_type)     :: Cloud_state
      type(particles_type)       :: Particles 
      type(precip_state_type)    :: Precip_state
      type(cloud_processes_type) :: Cloud_processes
      type(mp_lsdiag_type)       :: Lsdiag_mp
  
!------------------------------------------------------------------------
!   local variables:
!------------------------------------------------------------------------
      integer :: ix, jx, kx, nt, tr

!---------------------------------------------------------------------
!    begin the timing of the lscloud_driver subroutine.
!---------------------------------------------------------------------
      call mpp_clock_begin (lscloud_driver_clock)

!-----------------------------------------------------------------------
!    define input array sizes.
!-----------------------------------------------------------------------
      ix = size(Input_mp%t,1) 
      jx = size(Input_mp%t,2) 
      kx = size(Input_mp%t,3) 
      nt = size(Output_mp%rdt,4)

!---------------------------------------------------------------------
!    verify that the module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('lscloud_driver_mod',  &
                 'lscloud_driver_init has not been called.', FATAL)
      endif

!------------------------------------------------------------------------
!    allocate and initialize arrays to hold the diagnostic variables.
!------------------------------------------------------------------------
      allocate (Lsdiag_mp%diag_3d(ix, jx, 0:Lsdiag_mp_control%n_diag_4d))
      allocate    &
            (Lsdiag_mp%diag_4d(ix, jx, kx, 0:Lsdiag_mp_control%n_diag_4d))
      allocate (Lsdiag_mp%diag_4d_kp1  &
                        (ix, jx, kx+1, 0:Lsdiag_mp_control%n_diag_4d_kp1))
      Lsdiag_mp%diag_3d(:,:,0:) = 0.
      Lsdiag_mp%diag_4d(:,:,:,0:) = 0.
      Lsdiag_mp%diag_4d_kp1(:,:,:,0:) = 0.

!---------------------------------------------------------------------
!   reset the tendencies; on input they contain the convective scheme 
!   tendencies.
!---------------------------------------------------------------------
      Tend_mp%ttnd   = 0.
      Tend_mp%qtnd   = 0.
      Tend_mp%q_tnd  = 0.

!---------------------------------------------------------------------
!    allocate the derived type arrays used in the large-scale cloud
!    integration. initiate the clock to time this routine.
!-----------------------------------------------------------------------
      call mpp_clock_begin (lscloud_alloc_clock)
      call lscloud_alloc     &
           (ix, jx, kx, Input_mp, C2ls_mp, Particles, Cloud_processes, &
                                    Precip_state, Atmos_state, Cloud_state)
      call mpp_clock_end (lscloud_alloc_clock)

!------------------------------------------------------------------------
!    integrate the large-scale cloud equations.
!------------------------------------------------------------------------
      if (doing_prog_clouds) then


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!              SCHEMES USING PROGNOSTIC CLOUD VARIABLES
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


        if (do_legacy_strat_cloud) then    

!---------------------------------------------------------------------
!    start clock to time the legacy strat_cloud code.
!---------------------------------------------------------------------
          call mpp_clock_begin (stratcloud_clock)
          call strat_cloud    &
                (Lsdiag_mp, Lsdiag_mp_control, Time, is, ie, js, je, dt,  &
                 C2ls_mp, Atmos_state, Cloud_state, Input_mp, Tend_mp,   &
                 Cloud_processes, Removal_mp, Particles, Precip_state,   &
                 Aerosol = Aerosol)

!---------------------------------------------------------------------
!    stop clock for the legacy strat_cloud code.
!---------------------------------------------------------------------
          call mpp_clock_end (stratcloud_clock)

!----------------------------------------------------------------------
!    if legacy code not being executed, then macrophysics and microphysics
!    are being treated in separate modules. For the macrophysics, current
!    choices are tiedtke (tiedtke_macrophysics = .T.), or clubb 
!    macrophysics (do_clubb ==2).
!----------------------------------------------------------------------
        else

!------------------------------------------------------------------------
!    if tiedtke clouds are active, calculate needed svp arrays, initialize 
!    some debug option fields, make sure the cloud tracer arrays are 
!    consistent and of acceptable magnitude, and define particle
!    number variables in units needed for later use.
!------------------------------------------------------------------------
          if (Constants_lsc%tiedtke_macrophysics) then

!-----------------------------------------------------------------------
!    call compute_qs_a to calculate saturation vapor pressure under
!    existing model assumptions and cloud treatment (differs dependent on
!    whether pdf clouds active, whether supersaturation is being allowed).
!    also calculate relative humidity in grid box outside of cloud (so
!    grid box mean rh is preserved).
!-----------------------------------------------------------------------
            call mpp_clock_begin (polysvp_clock)
            call compute_qs_a (   &
                  ix, jx, kx, Input_mp%tin, Input_mp%qin, Input_mp%pfull, &
                  Cloud_state%ql_in, Cloud_state%qi_in, &
                  C2ls_mp%convective_humidity_ratio, Atmos_state)
            call mpp_clock_end (polysvp_clock)

!-----------------------------------------------------------------------
!    call lscloud_debug_setup to output header info when the debug 
!    option is active.
!-----------------------------------------------------------------------
            call mpp_clock_begin (lscloud_debug_clock)
            call lscloud_debug_setup (is, ie, js, je, Input_mp,    &
                                                               Atmos_state)
            call mpp_clock_end (lscloud_debug_clock)

!-----------------------------------------------------------------------
!    call impose_realizability to account for the fact that other processes
!    may have created negative tracer or extremely small values of tracer 
!    fields. start clock to time this subroutine.
!----------------------------------------------------------------------
            call mpp_clock_begin (realiz_clock)
            call impose_realizability (    &
                   Atmos_state, Cloud_state, Tend_mp%qtnd, Tend_mp%ttnd, &
                   Lsdiag_mp%diag_4d, Lsdiag_mp_control%diag_id,   &
                                               Lsdiag_mp_control%diag_pt)
            call mpp_clock_end (realiz_clock)

!------------------------------------------------------------------------
!    define particle numbers (in units of number / cm**3) after 
!    any realizability adjustment.
!------------------------------------------------------------------------
            if (do_liq_num) then
              Particles%N3d = Cloud_state%qn_upd*Atmos_state%airdens*1.e-6
              Particles%N3di = Cloud_state%qni_upd*   &
                                                 Atmos_state%airdens*1.e-6
            endif
          endif

!------------------------------------------------------------------------
!    if either tiedtke or clubb macrophysics are active, call 
!    determine_available_aerosol to determine the available condensation 
!    nuclei.
!------------------------------------------------------------------------
          if (Constants_lsc%tiedtke_macrophysics .or. do_clubb == 2) then
            call mpp_clock_begin (aerosol_cloud_clock)
            call determine_available_aerosol (    &
                   ix, jx, kx, Lsdiag_mp, Lsdiag_mp_control,  &
                   Atmos_state, Particles, Aerosol)
            call mpp_clock_end (aerosol_cloud_clock)
          endif

!---------------------------------------------------------------------
!    code is now ready to calculate macrophysics tendencies. start the 
!    clock for the ls cloud macrophysics code and call 
!    ls_cloud_macrophysics to handle grid-scale condensation effects.
!    either tiedtke or clubb  macrophysics will be employed, dependent on
!    the nml specification provided.
!---------------------------------------------------------------------
          call mpp_clock_begin (ls_macrophysics_clock)
          call ls_cloud_macrophysics (    &
                   is, ie, js, je, Time, dt, rdiag, Input_mp, Output_mp, &
                   Tend_mp, C2ls_mp, Lsdiag_mp_control, Lsdiag_mp, &
                   Atmos_state, Cloud_state, &
                   Particles, Precip_state, Cloud_processes, Aerosol)    

!---------------------------------------------------------------------
!    end the timing of the ls cloud macrophysics code section.
!---------------------------------------------------------------------
          call mpp_clock_end (ls_macrophysics_clock)

!------------------------------------------------------------------------
!    if  clubb is activated, the input fields are updated and subroutine 
!    impose_realizability_clubb is called to validate the consistency and 
!    magnitude of the cloud tracers before microphysics is calculated.
!------------------------------------------------------------------------
          if (do_clubb == 2) then
            Input_mp%tin = Input_mp%t + Output_mp%tdt*dt
            Input_mp%qin = Input_mp%q + Output_mp%rdt(:,:,:,1)*dt
            Input_mp%uin = Input_mp%u + Output_mp%udt*dt
            Input_mp%vin = Input_mp%v + Output_mp%vdt*dt
            do tr=1,size(Output_mp%rdt,4)
              Input_mp%tracer(:,:,:,tr) = Input_mp%r(:,:,:,tr) +   &
                                               Output_mp%rdt(:,:,:,tr)*dt
            end do
            do tr=size(Output_mp%rdt,4) +1, size(Input_mp%r,4)
              Input_mp%tracer(:,:,:,tr) = Input_mp%r(:,:,:,tr)
            end do

!----------------------------------------------------------------------
!    start clock to time realizability subroutine.
!----------------------------------------------------------------------
            call mpp_clock_begin (realiz_clubb_clock)
            call impose_realizability_clubb (   &
                 Input_mp, Tend_mp, Cloud_state,   &
                 C2ls_mp%convective_humidity_area, dtcloud, Lsdiag_mp, &
                                                       Lsdiag_mp_control)
            call mpp_clock_end (realiz_clubb_clock)

!------------------------------------------------------------------------
!    define particle numbers (in units of number / cm**3) after 
!    realizability adjustment.
!------------------------------------------------------------------------
            if (do_liq_num) then
              Particles%N3d = Cloud_state%qn_upd*Atmos_state%airdens*1.e-6
              Particles%N3di =    &
                             Cloud_state%qni_upd*Atmos_state%airdens*1.e-6
            endif
          endif

!-----------------------------------------------------------------------
!    when prognostic clouds are active, call ls_cloud_microphysics to 
!    compute the cloud microphysical effects using the requested 
!    microphysics scheme:
!    a) rotstayn-klein, b) MG (from Marc Salzmann, c) mg-ncar, a
!    newer version of the MG code with those changes made by Marc,
!    d) ncar, the latest available ncar microphysics version ,
!    e) lin_cld_microphysics.
!-----------------------------------------------------------------------
          if (doing_prog_clouds)  then

!---------------------------------------------------------------------
!    begin the timing of the ls cloud microphysics code section.
!---------------------------------------------------------------------
            call mpp_clock_begin (ls_microphysics_clock)
            call ls_cloud_microphysics    &
                (is, ie, js, je, Time, dt,  &
                 Input_mp, Output_mp, C2ls_mp, Tend_mp, Lsdiag_mp, &
                 Lsdiag_mp_control, &
                 Atmos_state, Cloud_state, Particles, Precip_state,  &
                                     Cloud_processes, Removal_mp, Aerosol)

!---------------------------------------------------------------------
!    end the timing of the ls cloud microphysics code section.
!---------------------------------------------------------------------
            call mpp_clock_end (ls_microphysics_clock)
          endif ! doing_prog_clouds

!-----------------------------------------------------------------------
!    call detailed_diagnostics to generate and output cloud budget
!    diagnostics and detailed netcdf and/or debug variable fields.
!    note that budgets may not be complete for CLUBB, but should balance
!    in non-CLUBB cases.
!-----------------------------------------------------------------------
          if (Constants_lsc%tiedtke_macrophysics .or. do_clubb > 0) then
            call mpp_clock_begin (detail_diag_clock)
            call detailed_diagnostics (      &
                is, ie, js, je, Time, Lsdiag_mp_control%n_diag_4d,      &
                Lsdiag_mp%diag_4d, Lsdiag_mp%diag_4d_kp1,&
                Lsdiag_mp%diag_3d, Lsdiag_mp_control%diag_pt,   &
                Lsdiag_mp_control%diag_id, &
                C2ls_mp, Input_mp, Atmos_state, Cloud_state, Particles, &
                Precip_state, Cloud_processes, Tend_mp, Removal_mp) 
            call mpp_clock_end (detail_diag_clock)
          endif

        endif  ! (do_legacy_strat_cloud)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!              SCHEMES USING NON-PROGNOSTIC CLOUDS
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!    if a non-prognostic cloud scheme is active, then call lscale_cond 
!    to calculate the temperature and specific humidity tendencies 
!    related to the latent heat release associated with the removal of
!    any large-scale "supersaturation".
!-----------------------------------------------------------------------
      else if (do_lsc) then   ! (doing_prog_clouds)
        call mpp_clock_begin (lscalecond_clock)
        call lscale_cond (Input_mp%tin, Input_mp%qin, Input_mp%pfull,   &
                        Input_mp%phalf, Input_mp%coldT,  &
                         Precip_state%surfrain, Precip_state%surfsnow,  &
                         Tend_mp%ttnd, Tend_mp%qtnd)
        call mpp_clock_end (lscalecond_clock)
      endif ! (doing_prog_clouds)

!------------------------------------------------------------------------
!    call update_fields_and_tendencies to apply computed cloud tendencies.
!------------------------------------------------------------------------
      call mpp_clock_begin (update_fields_clock)
      call update_fields_and_tendencies    &
                (is, ie, js, je, Time, 1.0/dt, Input_mp, Cld_props,  &
                                       Precip_state, Tend_mp, Output_mp)
      call mpp_clock_end (update_fields_clock)

!------------------------------------------------------------------------
!    deallocate the arrays associated with the large-scale 
!    cloud diagnostics.
!------------------------------------------------------------------------
      deallocate ( Lsdiag_mp%diag_4d )
      deallocate ( Lsdiag_mp%diag_4d_kp1 )
      deallocate ( Lsdiag_mp%diag_3d )
    

!------------------------------------------------------------------------
!    call compute_ls_wetdep to calculate wet deposition removal by lscloud
!    and to output related diagnostics.
!------------------------------------------------------------------------
      call mpp_clock_begin (ls_wetdep_clock)
      call compute_ls_wetdep (is, js, Time, Tend_mp, C2ls_mp, Input_mp,  &
                              Removal_mp, Cloud_processes%f_snow_berg, &
                              dt, Output_mp, Precip_state)
      call mpp_clock_end (ls_wetdep_clock)

!---------------------------------------------------------------------
!    output basic diagnostics defining the  model condensate and the
!    effect of the large-scale cloud/condensation scheme on atmospheric 
!    fields.
!---------------------------------------------------------------------
      call mpp_clock_begin (basic_diags_clock)
      call basic_diagnostics (is, js, Time, Tend_mp, Input_mp, &
                                                Precip_state, Removal_mp) 
      call mpp_clock_end (basic_diags_clock)

!-------------------------------------------------------------------------
!    call lscloud_dealloc to deallocate the local derived type variable 
!    components.
!-------------------------------------------------------------------------
      call mpp_clock_begin (dealloc_clock)
      call lscloud_dealloc (Atmos_state, Particles, Cloud_state,  &
                            Precip_state, Cloud_processes)
      call mpp_clock_end (dealloc_clock)

!---------------------------------------------------------------------
!    end the timing of the large-scale condensation code section.
!---------------------------------------------------------------------
      call mpp_clock_end (lscloud_driver_clock)
 
!-----------------------------------------------------------------------


end subroutine lscloud_driver 


!########################################################################

subroutine lscloud_driver_endts



      return



end subroutine lscloud_driver_endts


!########################################################################

subroutine lscloud_driver_end 


!------------------------------------------------------------------------

      integer :: lscalecond_term_clock, strat_term_clock,   &
                 polysvp_term_clock, aerosol_cloud_term_clock,  &
                 macrophysics_term_clock,  microphysics_term_clock, &
                 lscloud_netcdf_term_clock

!----------------------------------------------------------------------
!    activate lscloud_driver termination clock.
!----------------------------------------------------------------------
      call mpp_clock_begin (lscloud_driver_term_clock)

!------------------------------------------------------------------------
!    define clocks for termination.
!------------------------------------------------------------------------
      lscalecond_term_clock = mpp_clock_id(   &
                 '   Lscloud_driver: lscale_cond:Termination' , &
                                               grain=CLOCK_MODULE_DRIVER )
      strat_term_clock = mpp_clock_id(   &
                 '   Lscloud_driver: strat_cloud:Termination' , &
                                               grain=CLOCK_MODULE_DRIVER )
      polysvp_term_clock = mpp_clock_id(   &
                 '   Lscloud_driver: polysvp:Termination' , &
                                               grain=CLOCK_MODULE_DRIVER )
      aerosol_cloud_term_clock = mpp_clock_id(   &
                 '   Lscloud_driver: aerosol_cloud:Termination' , &
                                               grain=CLOCK_MODULE_DRIVER )
      macrophysics_term_clock = mpp_clock_id(   &
                 '   Lscloud_driver: ls_cloud_macrophysics:Termination' , &
                                               grain=CLOCK_MODULE_DRIVER )
      microphysics_term_clock = mpp_clock_id(   &
                 '   Lscloud_driver: ls_cloud_microphysics:Termination' , &
                                               grain=CLOCK_MODULE_DRIVER )
      lscloud_netcdf_term_clock = mpp_clock_id(   &
                 '   Lscloud_driver: lscloud_netcdf:Termination' , &
                                               grain=CLOCK_MODULE_DRIVER )

!-----------------------------------------------------------------------
!    close the modules which were opened by this module. time their
!    termination.
!-----------------------------------------------------------------------
      if (doing_prog_clouds) then 
        if (Nml_lsc%do_legacy_strat_cloud) then
          call mpp_clock_begin ( strat_term_clock)
          call strat_cloud_end 
          call mpp_clock_end ( strat_term_clock)
        else
          if (Constants_lsc%tiedtke_macrophysics) then
            call mpp_clock_begin ( polysvp_term_clock)
            call polysvp_end
            call mpp_clock_end ( polysvp_term_clock)
          endif
          if (Constants_lsc%tiedtke_macrophysics .or. do_clubb ==2) then
            call mpp_clock_begin ( aerosol_cloud_term_clock)
            call aerosol_cloud_end                 
            call mpp_clock_end ( aerosol_cloud_term_clock)
          endif
          call mpp_clock_begin ( macrophysics_term_clock)
          call ls_cloud_macrophysics_end                                
          call mpp_clock_end ( macrophysics_term_clock)
          call mpp_clock_begin ( microphysics_term_clock)
          call ls_cloud_microphysics_end 
          call mpp_clock_end ( microphysics_term_clock)
          if (Constants_lsc%tiedtke_macrophysics .or. do_clubb ==2) then
            call mpp_clock_begin ( lscloud_netcdf_term_clock)
            call lscloud_netcdf_end
            call mpp_clock_end ( lscloud_netcdf_term_clock)
          endif
        endif
      else  ! (doing_prog_clouds)
        if (do_lsc) then
          call mpp_clock_begin ( lscalecond_term_clock )
          call lscale_cond_end()
          if (do_rh_clouds) then
            call rh_clouds_end
          endif 
          call mpp_clock_end   ( lscalecond_term_clock )
        endif
      endif  ! (doing_prog_clouds)

!-----------------------------------------------------------------------
!    mark the module as uninitialized and terminate the clock for this
!    subroutine.
!-----------------------------------------------------------------------
      module_is_initialized = .false.
      call mpp_clock_end (lscloud_driver_term_clock)

!-----------------------------------------------------------------------

end subroutine lscloud_driver_end


!######################################################################



!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!
!   PRIVATE SUBROUTINES
!
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||




!#######################################################################

subroutine diag_field_init (axes, Time)

integer,                 intent(in) :: axes(4)
type(time_type),         intent(in) :: Time

!------------------------------------------------------------------------
!   local variables:

      character(len=32)     :: tracer_units, tracer_name
      character(len=128)    :: diaglname
      integer, dimension(3) :: half = (/1,2,4/)
      integer               :: n, nn

!---------------------------------------------------------------------
!    register diagnostic fields.
!---------------------------------------------------------------------

      id_lsc_precip = register_diag_field ( mod_name, &
        'lsc_precip', axes(1:2), Time, &
        'Precipitation rate from lsc ',       'kg/m2/s' )

      id_lsc_freq = register_diag_field ( mod_name, &
        'lsc_freq', axes(1:2), Time, &
        'frequency of precip from lsc ',       'number' , &
         missing_value = missing_value                       )

      id_lscale_rain3d= register_diag_field ( mod_name, &
        'lscale_rain3d', axes(half), Time, &
        'Rain fall rate from lscale  -3D ',   'kg(h2o)/m2/s' )

      id_lscale_snow3d= register_diag_field ( mod_name, &
        'lscale_snow3d', axes(half), Time, &
        'Snow fall rate from lscale -3D',       'kg(h2o)/m2/s' )
   
      id_lscale_precip3d= register_diag_field ( mod_name, &
        'lscale_precip3d', axes(1:3), Time, &
        'LS Precip falling out of gridbox',       'kg(h2o)/m2/s' , &
        mask_variant = .true., missing_value = missing_value)

      if (do_lsc ) then
        id_tdt_ls = register_diag_field ( mod_name, &
          'tdt_ls', axes(1:3), Time, &
          'Temperature tendency from large-scale cond',   'deg_K/s',  &
          missing_value=missing_value               )

        id_qdt_ls = register_diag_field ( mod_name, &
          'qdt_ls', axes(1:3), Time, &
          'Spec humidity tendency from large-scale cond', 'kg/kg/s',  &
          missing_value=missing_value               )

        id_prec_ls = register_diag_field ( mod_name, &
          'prec_ls', axes(1:2), Time, &
          'Precipitation rate from large-scale cond',     'kg/m2/s', &
          interp_method = "conserve_order1" )

        id_snow_ls = register_diag_field ( mod_name, &
          'snow_ls', axes(1:2), Time, &
          'Frozen precip rate from large-scale cond',     'kg/m2/s', &
          interp_method = "conserve_order1" )

        id_q_ls_col = register_diag_field ( mod_name, &
          'q_ls_col', axes(1:2), Time, &
          'Water vapor path tendency from large-scale cond','kg/m2/s' )
   
        id_t_ls_col = register_diag_field ( mod_name, &
          't_ls_col', axes(1:2), Time, &
          'Column static energy tendency from large-scale cond','W/m2' )
      endif

      ! register cmip diagnostics for large-scale clouds/precip
      if ( do_lsc .or. doing_prog_clouds ) then
        ID_tntscp = register_cmip_diag_field_3d ( mod_name, 'tntscp', Time, &
           'Tendency of Air Temperature Due to Stratiform Clouds and Precipitation', 'K s-1', &
           standard_name='tendency_of_air_temperature_due_to_stratiform_clouds_and_precipitation' )

        ID_tnhusscp = register_cmip_diag_field_3d ( mod_name, 'tnhusscp', Time, &
           'Tendency of Specific Humidity Due to Stratiform Clouds and Precipitation', 's-1', &
           standard_name='tendency_of_specific_humidity_due_to_stratiform_clouds_and_precipitation' )
      endif

      if (doing_prog_clouds ) then
       if (use_cf_metadata) then
         id_LWP = register_cmip_diag_field_2d ( mod_name, 'lwp', Time, &
                                        'Liquid Water Path', 'kg m-2', &
                 standard_name='atmosphere_mass_content_of_cloud_liquid_water' )

        id_IWP = register_cmip_diag_field_2d ( mod_name, 'iwp', Time, &
                                          'Ice Water Path', 'kg m-2', &
                         standard_name='atmosphere_mass_content_of_cloud_ice' )

       else
        id_LWP = register_diag_field ( mod_name, &
          'LWP', axes(1:2), Time, &
          'Liquid water path',                            'kg/m2'   )

        id_IWP = register_diag_field ( mod_name, &
          'IWP', axes(1:2), Time, &
          'Ice water path',                               'kg/m2'   )
       endif
        id_tdt_ls = register_diag_field ( mod_name, &
          'tdt_ls', axes(1:3), Time, &
          'Temperature tendency from strat cloud',        'deg_K/s',  &
          missing_value=missing_value               )

        id_qdt_ls = register_diag_field ( mod_name, &
          'qdt_ls', axes(1:3), Time, &
          'Spec humidity tendency from strat cloud',      'kg/kg/s',  &
          missing_value=missing_value               )

        id_prec_ls = register_diag_field ( mod_name, &
          'prec_ls', axes(1:2), Time, &
          'Precipitation rate from strat cloud',          'kg/m2/s' )

        id_snow_ls = register_diag_field ( mod_name, &
          'snow_ls', axes(1:2), Time, &
          'Frozen precip rate from strat cloud',          'kg/m2/s' )

        id_q_ls_col = register_diag_field ( mod_name, &
          'q_ls_col', axes(1:2), Time, &
          'Water vapor path tendency from strat cloud',   'kg/m2/s' )
   
        id_t_ls_col = register_diag_field ( mod_name, &
          't_ls_col', axes(1:2), Time, &
          'Column static energy tendency from strat cloud','W/m2' )
   
        id_qldt_ls = register_diag_field ( mod_name, &
          'qldt_ls', axes(1:3), Time, &
          'Liquid water tendency from strat cloud',       'kg/kg/s',  &
          missing_value=missing_value               )

        if (do_liq_num) then
          id_qndt_ls = register_diag_field ( mod_name, &
            'qndt_ls', axes(1:3), Time, &
            'Drop number tendency from strat cloud',        '#/kg/s',  &
            missing_value=missing_value               )
        endif

        id_qidt_ls = register_diag_field ( mod_name, &
          'qidt_ls', axes(1:3), Time, &
          'Ice water tendency from strat cloud',          'kg/kg/s',  &
          missing_value=missing_value               )

        if (do_ice_num) then
          id_qnidt_ls = register_diag_field ( mod_name, &
            'qnidt_ls', axes(1:3), Time, &
            'Ice number tendency from strat cloud',          '#/kg/s',  &
            missing_value=missing_value               )
        endif

        id_qadt_ls = register_diag_field ( mod_name, &
          'qadt_ls', axes(1:3), Time, &
          'Cloud fraction tendency from strat cloud',     '1/sec',    &
          missing_value=missing_value               )

        id_ql_ls_col = register_diag_field ( mod_name, &
          'ql_ls_col', axes(1:2), Time, &
          'Liquid water path tendency from strat cloud',   'kg/m2/s' )
   
        if (do_liq_num) then
          id_qn_ls_col = register_diag_field ( mod_name, &
            'qn_ls_col', axes(1:2), Time, &
            'Column drop number tendency from strat cloud',  '#/m2/s' )
        endif

        if (do_ice_num) then
          id_qni_ls_col = register_diag_field ( mod_name, &
          'qni_ls_col', axes(1:2), Time, &
          'Column ice particle number tendency from strat cloud', '#/m2/s')
        endif

        id_qi_ls_col = register_diag_field ( mod_name, &
          'qi_ls_col', axes(1:2), Time, &
          'Ice water path tendency from strat cloud',      'kg/m2/s' )
   
        id_qa_ls_col = register_diag_field ( mod_name, &
          'qa_ls_col', axes(1:2), Time, &
          'Cloud mass tendency from strat cloud',          'kg/m2/s' )
      
        id_enth_ls_col = register_diag_field ( mod_name, &
          'enth_ls_col', axes(1:2), Time, &
          'Column enthalpy tendency from strat cloud','W/m2' )
 
        id_wat_ls_col = register_diag_field ( mod_name, &
          'wat_ls_col', axes(1:2), Time, &
          'Column total water tendency from strat cloud','kg/m2/s' )

      endif  !(doing_prog_clouds)

      id_qvout = register_diag_field ( mod_name, &
        'qvout', axes(1:3), Time, 'qv after strat_cloud', 'kg/kg', &
        missing_value=missing_value               )

      id_qaout = register_diag_field ( mod_name, &
        'qaout', axes(1:3), Time, 'qa after strat_cloud', 'none', &
        missing_value=missing_value               )

      id_qlout = register_diag_field ( mod_name, &
        'qlout', axes(1:3), Time, 'ql after strat_cloud', 'kg/kg', &
        missing_value=missing_value               )

      id_qiout = register_diag_field ( mod_name, &
        'qiout', axes(1:3), Time, 'qi after strat_cloud', 'kg/kg', &
        missing_value=missing_value               )

      if (do_liq_num) then
        id_qnout = register_diag_field ( mod_name, &
          'qnout', axes(1:3), Time, 'qn after strat_cloud', '#/kg', &
          missing_value=missing_value               )
      endif

      if (do_ice_num) then
        id_qniout = register_diag_field ( mod_name, &
          'qniout', axes(1:3), Time, 'qni after strat_cloud', '#/kg', &
          missing_value=missing_value               )
      endif

      id_f_snow_berg   =  register_diag_field ( mod_name, &
        'f_snow_berg', axes(1:3), Time,  &
        'fraction of snow/ice produced having IFN', 'fraction',  &
        missing_value=missing_value)

      id_f_snow_berg_cond   =  register_diag_field ( mod_name, &
        'f_snow_berg_cond',  axes(1:3), Time,  &
        'conditional fraction of snow/ice produced having IFN',  &
        'fraction',  mask_variant = .true., missing_value=missing_value)

      id_f_snow_berg_wtd   =  register_diag_field ( mod_name, &
        'f_snow_berg_wtd',  axes(1:3), Time,  &
        'product of snow/ice produced having IFN and ls precip falling&
        & out of gridbox', 'kg(h2o)/m2/s', mask_variant = .true.,   &
        missing_value=missing_value)

!---------------------------------------------------------------------
!    register the diagnostics associated with wet deposition of each 
!    tracer species. 
!---------------------------------------------------------------------
      allocate (id_wet_deposition(num_prog_tracers))
      id_wet_deposition = -1

      do n = 1, num_prog_tracers
        call get_tracer_names (MODEL_ATMOS, n, name = tracer_name,  &
                                                    units = tracer_units)
        id_wet_deposition(n) = register_diag_field( mod_name, &
          trim(tracer_name)//'_wetdep', axes(1:3), Time, &
          trim(tracer_name)//' tendency from wet deposition',  &
          TRIM(tracer_units)//'/sec', missing_value=missing_value )
      end do

!------------------------------------------------------------------------


end subroutine diag_field_init



!########################################################################

subroutine lscloud_alloc (    &
         idim, jdim, kdim, Input_mp, C2ls_mp, Particles, Cloud_processes, &
                                    Precip_state, Atmos_state, Cloud_state)

integer,                    intent(in)    :: idim,jdim,kdim
type(mp_input_type),        intent(in)    :: Input_mp
type(mp_conv2ls_type),      intent(in)    :: C2ls_mp
type(particles_type),       intent(inout) :: Particles
type(cloud_processes_type), intent(inout) :: Cloud_processes
type(precip_state_type),    intent(inout) :: Precip_state   
type(atmos_state_type),     intent(inout) :: Atmos_state
type(cloud_state_type),     intent(inout) :: Cloud_state


!----------------------------------------------------------------------
!   local variables:

      integer :: i,j,k
      real, dimension(idim,jdim,kdim) :: airdens_aerosol, T_aerosol

!-----------------------------------------------------------------------
!    allocate and initialize the components of the particles_type 
!    variable Particles.
!-----------------------------------------------------------------------
      allocate (Particles%concen_dust_sub   (idim, jdim, kdim) )
      allocate (Particles%drop1             (idim, jdim, kdim) )
      allocate (Particles%drop2             (idim, jdim, kdim) )
      allocate (Particles%crystal1          (idim, jdim, kdim) )
      allocate (Particles%Ndrop_act_CLUBB   (idim, jdim, kdim) )
      allocate (Particles%Icedrop_act_CLUBB (idim, jdim, kdim) )
      allocate (Particles%rbar_dust         (idim, jdim, kdim) )
      allocate (Particles%ndust             (idim, jdim, kdim) )
      allocate (Particles%hom               (idim, jdim, kdim) )
      allocate (Particles%N3D               (idim, jdim, kdim) )
      allocate (Particles%N3Di              (idim, jdim, kdim) )
      allocate (Particles%totalmass1        (idim, jdim, kdim, n_totmass) )
      allocate (Particles%imass1            (idim, jdim, kdim, n_imass) )
 
      Particles%concen_dust_sub   = 0.
      Particles%drop1    = 0.
      Particles%drop2    = 0.
      Particles%crystal1    = 0.
      Particles%Ndrop_act_CLUBB  = 0.
      Particles%Icedrop_act_CLUBB  = 0.
      Particles%rbar_dust   = 0.
      Particles%ndust   = 0.
      Particles%hom   = 0.
      if (do_liq_num) then
        Particles%N3d   = 0.
      else
        do k=1,kdim
          Particles%N3d(:,:,k)   = N_land*Input_mp%land(:,:) + &
                                          N_ocean*(1.0-Input_mp%land(:,:))
        end do
      endif
      Particles%N3di  = 0.
      Particles%imass1 = 0.
      Particles%totalmass1 = 0.

!-----------------------------------------------------------------------
!    allocate and initialize the components of the cloud_processes_type
!    variable Cloud_processes.
!-----------------------------------------------------------------------
      allocate  (Cloud_processes%da_ls            (idim, jdim, kdim) )
      allocate  (Cloud_processes%D_eros           (idim, jdim, kdim) )
      allocate  (Cloud_processes%qvg              (idim, jdim, kdim) )
      allocate  (Cloud_processes%dcond_ls         (idim, jdim, kdim) )
      allocate  (Cloud_processes%dcond_ls_liquid  (idim, jdim, kdim) )
      allocate  (Cloud_processes%dcond_ls_ice     (idim, jdim, kdim) )
      allocate  (Cloud_processes%dcond_ls_tot     (idim, jdim, kdim) )
      allocate  (Cloud_processes%delta_cf         (idim, jdim, kdim) )
      allocate  (Cloud_processes%f_snow_berg      (idim, jdim, kdim) )
 
      Cloud_processes%da_ls          = 0.
      Cloud_processes%D_eros         = 0.
      Cloud_processes%qvg            = 0.
      Cloud_processes%dcond_ls       = 0.
      Cloud_processes%dcond_ls_liquid  = 0.
      Cloud_processes%dcond_ls_ice   = 0.
      Cloud_processes%dcond_ls_tot   = 0.
      Cloud_processes%delta_cf       = 0.
      Cloud_processes%f_snow_berg    = 0.

!-----------------------------------------------------------------------
!    allocate and initialize the components of the precip_state_type 
!    variable Precip_state.
!-----------------------------------------------------------------------
      allocate (Precip_state%lsc_snow      (idim, jdim, kdim) )
      allocate (Precip_state%lsc_rain      (idim, jdim, kdim) )
      allocate (Precip_state%lsc_snow_size (idim, jdim, kdim) )
      allocate (Precip_state%lsc_rain_size (idim, jdim, kdim) )
      allocate (Precip_state%qsout3d_mg    (idim, jdim, kdim) )
      allocate (Precip_state%qrout3d_mg    (idim, jdim, kdim) )
      allocate (Precip_state%precip        (idim, jdim) )
      allocate (Precip_state%surfrain      (idim, jdim) )
      allocate (Precip_state%surfsnow      (idim, jdim) )

      Precip_state%lsc_snow      = 0.
      Precip_state%lsc_rain      = 0.
      Precip_state%lsc_snow_size = 0.
      Precip_state%lsc_rain_size = 0.
      Precip_state%qsout3d_mg    = 0.
      Precip_state%qrout3d_mg    = 0.
      Precip_state%precip        = 0.
      Precip_state%precip        = 0.
      Precip_state%surfrain      = 0.
      Precip_state%surfsnow      = 0.


!-----------------------------------------------------------------------
!    allocate and initialize the components of the atmos_state_type 
!    variable Atmos_state.
!-----------------------------------------------------------------------
      allocate (Atmos_state%airdens        (idim, jdim, kdim) )
      allocate (Atmos_state%tn             (idim, jdim, kdim) )
      allocate (Atmos_state%qvn            (idim, jdim, kdim) )
      allocate (Atmos_state%qs             (idim, jdim, kdim) )
      allocate (Atmos_state%dqsdT          (idim, jdim, kdim) )
      allocate (Atmos_state%qsi            (idim, jdim, kdim) )
      allocate (Atmos_state%qsl            (idim, jdim, kdim) )
      allocate (Atmos_state%rh_crit        (idim, jdim, kdim) )
      allocate (Atmos_state%rh_crit_min    (idim, jdim, kdim) )
      allocate (Atmos_state%gamma          (idim, jdim, kdim) )
      allocate (Atmos_state%esat0          (idim, jdim, kdim) )
      allocate (Atmos_state%U_ca           (idim, jdim, kdim) )
      allocate (Atmos_state%delp           (idim, jdim, kdim) )
      allocate (Atmos_state%U01            (idim, jdim, kdim) )
      allocate (Atmos_state%pthickness     (idim, jdim, kdim) )

      if (do_clubb > 0) then
        T_aerosol = Input_mp%t
      else
        T_aerosol = Input_mp%tin
      endif

!-----------------------------------------------------------------------
!    calculate air density.   
!-----------------------------------------------------------------------

      if (nql == NO_TRACER .or. nqi == NO_TRACER)  then
        where (C2ls_mp%convective_humidity_ratio .gt. 0.) 
          Atmos_state%airdens = Input_mp%pfull/(RDGAS*Input_mp%tin*  &
            (1. + (d608*Input_mp%qin/C2ls_mp%convective_humidity_ratio) ))
        elsewhere
          Atmos_state%airdens = Input_mp%pfull/(RDGAS*Input_mp%tin)
        end where
        if (do_clubb > 0) then
          airdens_aerosol = Input_mp%pfull/(RDGAS*Input_mp%t)
        else
          airdens_aerosol = Atmos_state%airdens
        endif
      else  ! (NO_TRACER)
        where (C2ls_mp%convective_humidity_ratio .gt. 0.) 
          Atmos_state%airdens = Input_mp%pfull/(RDGAS*Input_mp%tin*  &
           (1. + (d608*Input_mp%qin/C2ls_mp%convective_humidity_ratio) -  &
                Input_mp%tracer(:,:,:,nql) - Input_mp%tracer(:,:,:,nqi))) 
        elsewhere
          Atmos_state%airdens = Input_mp%pfull/(RDGAS*Input_mp%tin*  &
           (1.  - Input_mp%tracer(:,:,:,nql) - Input_mp%tracer(:,:,:,nqi)))
        end where
        if (do_clubb > 0) then
          airdens_aerosol = Input_mp%pfull/(RDGAS*Input_mp%t*  &
                (1.  - Input_mp%r(:,:,:,nql) - Input_mp%r(:,:,:,nqi)))
        else
          airdens_aerosol = Atmos_state%airdens
        endif
      endif ! (NO_TRACER)

!------------------------------------------------------------------------
!    define layer thickness.
!------------------------------------------------------------------------
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            if (Input_mp%phalf(i,j,k) < 1.0) then
              Atmos_state%pthickness(i,j,k) =  &
                       (Input_mp%phalf(i,j,k+1) - Input_mp%phalf(i,j,k))/&
                                               GRAV/airdens_aerosol(i,j,k)
            else
              Atmos_state%pthickness(i,j,k) =   &
                 log(Input_mp%phalf(i,j,k+1)/Input_mp%phalf(i,j,k))*  &
                                      8.314*T_aerosol(i,j,k)/(9.8*0.02888) 
            end if
          end do
        end do
      end do
      Atmos_state%tn            = Input_mp%tin
      Atmos_state%qvn           = 0.
      Atmos_state%qs            = 0.
      Atmos_state%dqsdT         = 0.
      Atmos_state%qsi           = 0.
      Atmos_state%qsl           = 0.
      Atmos_state%rh_crit       = 1.
      Atmos_state%rh_crit_min   = 1.
      Atmos_state%gamma         = 0.
      Atmos_state%esat0         = 0.
      Atmos_state%U_ca          = 0.
      do k = 1, kdim
        Atmos_state%delp(:,:,k) =   &
                          Input_mp%phalf(:,:,k+1) - Input_mp%phalf(:,:,k)
      enddo     
      Atmos_state%U01           = 0.

!-----------------------------------------------------------------------
!    allocate and initialize the components of the cloud_state_type 
!    variable Cloud_state.
!-----------------------------------------------------------------------
      allocate (Cloud_state%ql_upd    (idim, jdim, kdim) )
      allocate (Cloud_state%qi_upd    (idim, jdim, kdim) )
      allocate (Cloud_state%qa_upd    (idim, jdim, kdim) )
      allocate (Cloud_state%qn_upd    (idim, jdim, kdim) )
      allocate (Cloud_state%qni_upd   (idim, jdim, kdim) )
  
      allocate (Cloud_state%ql_mean   (idim, jdim, kdim) )
      allocate (Cloud_state%qi_mean   (idim, jdim, kdim) )
      allocate (Cloud_state%qa_mean   (idim, jdim, kdim) )
      allocate (Cloud_state%qn_mean   (idim, jdim, kdim) )
      allocate (Cloud_state%qni_mean  (idim, jdim, kdim) )
 
      allocate (Cloud_state%ql_in     (idim, jdim, kdim) )
      allocate (Cloud_state%qi_in     (idim, jdim, kdim) )
      allocate (Cloud_state%qa_in     (idim, jdim, kdim) )
      allocate (Cloud_state%qn_in     (idim, jdim, kdim) )
      allocate (Cloud_state%qni_in    (idim, jdim, kdim) )

      allocate (Cloud_state%SL_out    (idim, jdim, kdim) )
      allocate (Cloud_state%SI_out    (idim, jdim, kdim) )
      allocate (Cloud_state%SA_out    (idim, jdim, kdim) )
      allocate (Cloud_state%SN_out    (idim, jdim, kdim) )
      allocate (Cloud_state%SNi_out   (idim, jdim, kdim) )
      allocate (Cloud_state%qcvar_clubb   (idim, jdim, kdim) )
      allocate (Cloud_state%relvarn       (idim, jdim, kdim) )

      allocate (Cloud_state%qa_upd_0     (idim, jdim, kdim) )
      allocate (Cloud_state%SA_0         (idim, jdim, kdim) )

      Cloud_state%ql_upd    = 0.
      Cloud_state%qi_upd    = 0.
      Cloud_state%qa_upd    = 0.
      Cloud_state%qn_upd    = 0.
      Cloud_state%qni_upd   = 0.

      Cloud_state%ql_mean    = 0.      
      Cloud_state%qi_mean    = 0.      
      Cloud_state%qa_mean    = 0.      
      Cloud_state%qn_mean    = 0.      
      Cloud_state%qni_mean   = 0.      

      if (nql == NO_TRACER) then
        Cloud_state%ql_in = 0.    
      else
        Cloud_state%ql_in = Input_mp%tracer(:,:,:,nql)
      endif
      if (nqi == NO_TRACER) then
        Cloud_state%qi_in = 0.    
      else
        Cloud_state%qi_in = Input_mp%tracer(:,:,:,nqi)
      endif
      if (nqa == NO_TRACER) then
        Cloud_state%qa_in = 0.    
      else
        Cloud_state%qa_in = Input_mp%tracer(:,:,:,nqa)
      endif
      if (nqn == NO_TRACER) then
        Cloud_state%qn_in = 0.    
      else
        Cloud_state%qn_in = Input_mp%tracer(:,:,:,nqn)
      end if

      if (nqni == NO_TRACER) then
        Cloud_state%qni_in = 0.    
      else
        Cloud_state%qni_in = Input_mp%tracer(:,:,:,nqni)
      end if

      Cloud_state%SL_out  = 0.
      Cloud_state%SI_out  = 0.
      Cloud_state%SA_out  = 0.
      Cloud_state%SN_out  = 0.
      Cloud_state%SNi_out = 0.
    
      Cloud_state%qcvar_clubb  = 0.
      Cloud_state%relvarn      = 0.
      Cloud_state%qa_upd_0 = 0.
      Cloud_state%SA_0        = 0.


!------------------------------------------------------------------------


end subroutine lscloud_alloc 



!#######################################################################

subroutine impose_realizability (   &
                           Atmos_state, Cloud_state, SQ_out, ST_out,  &
                                                diag_4d, diag_id, diag_pt)

!-----------------------------------------------------------------------
!    account for the fact that other processes may have created negative 
!    tracer or extremely small values of tracer fields. the general reason
!    for the extremely small values of the tracer fields is due to vertical
!    diffusion, advection of condensate or cumulus induced subsidence (also
!    a form of advection) of condensate.
!    in this step any values of the prognostic variables which are less 
!    than qmin are reset to zero, while conserving total moisture.
!    note that this is done slightly different for the Tiedtke cloud 
!    fraction than it is for pdf clouds. In the former, the filling 
!    requires that cloud liquid, cloud ice, and cloud fraction are greater
!    than qmin. For PDF clouds, cloud fraction need not be considered 
!    since it is diagnosed later from the PDF cloud field.
!-----------------------------------------------------------------------

!-------------------------------------------------------------------------
type( atmos_state_type),    intent(inout) :: Atmos_state
type( cloud_state_type),    intent(inout) :: Cloud_state
real, dimension(:,:,:),     intent(inout) :: ST_out, SQ_out
real, dimension(:,:,:,0:), intent(inout) :: diag_4d
type(diag_id_type),         intent(in)    :: diag_id
type(diag_pt_type),         intent(in)    :: diag_pt
!-------------------------------------------------------------------------

!------------------------------------------------------------------------
!   local variables:

      logical, dimension(size(ST_out,1), size(ST_out,2),   &
                           size(ST_out,3)) :: ql_too_small, qi_too_small 

      integer      :: idim,jdim,kdim
      integer      :: i,j,k


!----------------------------------------------------------------------
!    define dimensions.
!----------------------------------------------------------------------
      idim = size(ST_out,1)
      jdim = size(ST_out,2)
      kdim = size(ST_out,3)

!-----------------------------------------------------------------------
!    for the non-pdf scheme,  assure that cloud fraction is greater than 
!    qmin.  if it is not, set it to 0.0 and save the tendency and updated
!    value. save a diagnostic if requested.
!-----------------------------------------------------------------------
      if (.not. do_pdf_clouds) then
        where (Cloud_state%qa_in .le. qmin)
          Cloud_state%SA_out = Cloud_state%SA_out - Cloud_state%qa_in
          Cloud_state%qa_upd = 0.
        elsewhere
          Cloud_state%qa_upd = Cloud_state%qa_in     
        end where
        if (diag_id%qadt_fill + diag_id%qa_fill_col > 0 ) then
          where (Cloud_state%qa_in .le.        qmin)
            diag_4d(:,:,:,diag_pt%qadt_fill) =  -Cloud_state%qa_in*   &
                                                               inv_dtcloud
          endwhere
        end if

!------------------------------------------------------------------------ 
!    define the max cloud area permitted (U01), while maintaining the 
!    grid-box-mean RH, under the assumption that the cloudy air is satur-
!    ated and the temperature inside and outside of the cloud are ~ the 
!    same. save a diagnostic indicating the amount the cloud area was 
!    reduced due to this requirenment,
!------------------------------------------------------------------------ 
        Atmos_state%U01 = min(Atmos_state%U01, 1.)
        if (diag_id%qadt_rhred + diag_id%qa_rhred_col >0) then
          where (Cloud_state%qa_upd .gt. Atmos_state%U01)
            diag_4d(:,:,:,diag_pt%qadt_rhred) = -(Cloud_state%qa_upd -  &
                                              Atmos_state%U01)*inv_dtcloud
          endwhere
        end if
        
        where (Cloud_state%qa_upd .gt. Atmos_state%U01)
          Cloud_state%SA_out = Cloud_state%SA_out + Atmos_state%U01 -   &
                                                      Cloud_state%qa_upd
          Cloud_state%qa_upd = Atmos_state%U01      
        end where
      endif

!-------------------------------------------------------------------------
!    define the conditions under which liquid and ice filling must be done,
!    for both the pdf and non-pdf cases.
!    for the non-pdf scheme, the filling requires that cloud liquid, cloud
!    ice, and cloud fraction are greater than qmin.
!    for pdf clouds, cloud fraction need not be considered since it is 
!    diagnosed later from the PDF cloud field.
!-------------------------------------------------------------------------
      if (.not. do_pdf_clouds) then
        if (do_liq_num) then
          ql_too_small = (Cloud_state%ql_in .le. qmin .or.   &
                          Cloud_state%qa_in .le. qmin .or.   &
                          Cloud_state%qn_in .le. qmin)
        else
          ql_too_small = (Cloud_state%ql_in .le. qmin .or.   &
                          Cloud_state%qa_in .le. qmin)
        endif
        if (do_predicted_ice_number) then
          qi_too_small = (Cloud_state%qi_in .le.  qmin .or.   &
                          Cloud_state%qa_in .le. qmin .or.   &
                          Cloud_state%qni_in .le. qmin)
        else
          qi_too_small = (Cloud_state%qi_in .le.  qmin .or.   &
                          Cloud_state%qa_in .le. qmin)
        endif
      else
        if (do_liq_num) then
          ql_too_small = (Cloud_state%ql_in .le.  qmin  .or.   &
                          Cloud_state%qn_in .le.  qmin)
        else
          ql_too_small = (Cloud_state%ql_in .le.  qmin)
        endif
        if (do_predicted_ice_number) then
          qi_too_small = (Cloud_state%qi_in .le.  qmin .or.   &
                          Cloud_state%qni_in .le. qmin)
        else
          qi_too_small = (Cloud_state%qi_in .le.  qmin )
        endif
      endif
  
!------------------------------------------------------------------------
!    call subroutine adjust_condensate to conservatively fill ql if needed.
!------------------------------------------------------------------------
      call adjust_condensate (ql_too_small, Cloud_state%SL_out, &
            SQ_out, ST_out, Cloud_state%ql_in, HLV, Cloud_state%ql_upd)

!------------------------------------------------------------------------
!    save diagnostics defining the cloud liquid filling amount. 
!------------------------------------------------------------------------
      if ( diag_id%qldt_fill  + diag_id%ql_fill_col > 0 ) then
        where (ql_too_small )
          diag_4d(:,:,:,diag_pt%qldt_fill) =   &
                            -1.*Cloud_state%ql_in*inv_dtcloud
        endwhere
      end if
      if ( diag_id%qdt_liquid_init > 0 ) then
        where (ql_too_small )
          diag_4d(:,:,:,diag_pt%qdt_liquid_init) =   &
                            Cloud_state%ql_in*inv_dtcloud
        endwhere
      end if

!---------------------------------------------------------------------
!   define diagnostics for the initial large-scale liquid and ice cloud 
!   fractions, after they have been adjusted to lie within valid limits.
!---------------------------------------------------------------------
      if (diag_id%cf_liq_init   > 0)  then
        where (ql_too_small)
          diag_4d(:,:,:,diag_pt%cf_liq_init) = 0.
        elsewhere
          diag_4d(:,:,:,diag_pt%cf_liq_init) = MIN(Cloud_state%qa_in, 1.0)
        end where
      endif 

      if (diag_id%cf_ice_init > 0) then
        where (qi_too_small)
          diag_4d(:,:,:,diag_pt%cf_ice_init) = 0.
        elsewhere
          diag_4d(:,:,:,diag_pt%cf_ice_init) = MIN(Cloud_state%qa_in, 1.0)
        end where
      endif

!------------------------------------------------------------------------
!    adjust cloud droplet numbers as needed when those fields are being 
!    predicted. if droplet number is not being predicted, values were 
!    set at allocation.  
!------------------------------------------------------------------------
      if (do_liq_num) then
        call adjust_particle_number (ql_too_small, Cloud_state%SN_out, &
                cloud_state%qn_in, Cloud_state%qn_upd)

!------------------------------------------------------------------------
!    save a diagnostic defining the cloud droplet number filling amount. 
!------------------------------------------------------------------------
        if ( diag_id%qndt_fill  + diag_id%qn_fill_col + &
                    diag_id%qldt_fill + diag_id%ql_fill_col > 0 ) then
          where (ql_too_small )
            diag_4d(:,:,:,diag_pt%qndt_fill) =    &
                           -1.*Cloud_state%qn_in*inv_dtcloud
          endwhere
        end if
      endif

!------------------------------------------------------------------------
!    call subroutine adjust_condensate to conservatively fill qi if needed.
!------------------------------------------------------------------------
      call adjust_condensate (qi_too_small, Cloud_state%SI_out, &
            SQ_out, ST_out, CLoud_state%qi_in, HLS, Cloud_state%qi_upd)

!------------------------------------------------------------------------
!    save diagnostics defining the cloud ice filling amount. 
!------------------------------------------------------------------------
      if (diag_id%qidt_fill  + diag_id%qi_fill_col > 0 ) then
        where (qi_too_small )
          diag_4d(:,:,:,diag_pt%qidt_fill) =     &
                          -1.*Cloud_state%qi_in*inv_dtcloud
        endwhere
      end if
      if (diag_id%qdt_ice_init > 0 ) then
        where (qi_too_small )
          diag_4d(:,:,:,diag_pt%qdt_ice_init) =     &
                             Cloud_state%qi_in*inv_dtcloud
        endwhere
      end if

!------------------------------------------------------------------------
!    adjust ice particle numbers as needed when those fields
!    are being predicted. 
!------------------------------------------------------------------------
      if (do_ice_num) then
        call mpp_clock_begin (lscloud_debug_clock)
        call write_debug_output (" SNi 00 ", Cloud_state%SNi_out, &
                                                  scalar_in = inv_dtcloud)
        call mpp_clock_end (lscloud_debug_clock)
        call adjust_particle_number (qi_too_small, Cloud_state%SNi_out, &
                                  Cloud_state%qni_in, Cloud_state%qni_upd)

!------------------------------------------------------------------------
!    save a diagnostic defining the ice crystal number filling amount. 
!------------------------------------------------------------------------
        if (diag_id%qnidt_fill  + diag_id%qni_fill_col > 0 ) then
          where (qi_too_small) 
            diag_4d(:,:,:,diag_pt%qnidt_fill) =  &
                          -1.*Cloud_state%qni_in*inv_dtcloud
          endwhere
        end if

!------------------------------------------------------------------------
!    output debug information if desired. 
!------------------------------------------------------------------------
        call mpp_clock_begin (lscloud_debug_clock)
        call write_debug_output (" SNi 01 ", Cloud_state%SNi_out, &
                                    scalar_in = inv_dtcloud)
        call write_debug_output ("     qnidt_fill ", &
                                       diag_4d(:,:,:,diag_pt%qnidt_fill))
        call mpp_clock_end (lscloud_debug_clock)
      endif

!-----------------------------------------------------------------------


end subroutine impose_realizability 



!########################################################################

subroutine impose_realizability_clubb   &
                            (Input_mp, Tend_mp, Cloud_state, ahuco3d,    &
                                    dtcloud, Lsdiag_mp, Lsdiag_mp_control)

real,                          intent(in)    :: dtcloud
type(mp_tendency_type),        intent(inout) :: Tend_mp
type(mp_input_type),           intent(in   ) :: Input_mp
type(mp_lsdiag_type),          intent(inout) :: Lsdiag_mp
type(mp_lsdiag_control_type),  intent(inout) :: Lsdiag_mp_control
type(cloud_state_type),        intent(inout) :: Cloud_state
real, dimension(:,:,:),        intent(in   ) :: ahuco3d            

!-----------------------------------------------------------------------
!   local variables:

      logical, dimension(size(ahuco3d,1), size(ahuco3d,2),   &
                           size(ahuco3d,3)) :: ql_too_small, qi_too_small 
      integer :: idim, jdim, kdim
      integer :: i, j, k

!------------------------------------------------------------------------
!    define local variables.
!------------------------------------------------------------------------
      idim = size(Input_mp%tracer,1)
      jdim = size(Input_mp%tracer,2)
      kdim = size(Input_mp%tracer,3)

!-----------------------------------------------------------------------
!    account for the fact that other processes may have created negative 
!    tracer or extremely small values of tracer fields. in this step any 
!    values of the prognostic variables which are less than qmin are 
!    reset to zero, while conserving total moisture.
!----------------------------------------------------------------------
      where (Input_mp%tracer(:,:,:,nqa) .le. qmin)
        Tend_mp%q_tnd(:,:,:,nqa) = Tend_mp%q_tnd(:,:,:,nqa) -   &
                                              Input_mp%tracer(:,:,:,nqa)
        Cloud_state%qa_upd = 0.
      elsewhere
        Cloud_state%qa_upd = Input_mp%tracer(:,:,:,nqa)
      end where

!------------------------------------------------------------------------
!    total cloud fraction should be no greater than 1.0, i.e, sum of 
!    large-scale cloud fraction (qa_upd) and convective cloud fraction 
!    (ahuco) qa_upd+ahuco <= 1.0
!------------------------------------------------------------------------
      where (Cloud_state%qa_upd .gt. (1.0 - ahuco3d) )
        Tend_mp%q_tnd(:,:,:,nqa) = Tend_mp%q_tnd(:,:,:,nqa) +   &
                                    (1.0 - ahuco3d) - Cloud_state%qa_upd
        Cloud_state%qa_upd = (1.0 - ahuco3d)
      end where

!-----------------------------------------------------------------------
!    define conditions under which condensate is removed before calling
!    microphysics.
!-----------------------------------------------------------------------
      do k = 1,kdim
        do j = 1,jdim
          do i = 1,idim
            ql_too_small(i,j,k) =   &
                     Input_mp%tracer(i,j,k,nql) .le. qmin .or.   &
                     Input_mp%tracer(i,j,k,nqa) .le. qmin .or.   &
                     Input_mp%tracer(i,j,k,nqn)*1.e-6 .le. qmin 
            qi_too_small(i,j,k) =   &
                     Input_mp%tracer(i,j,k,nqi).le.qmin .or.   &
                     Input_mp%tracer(i,j,k,nqa).le.qmin .or.   &
                     Input_mp%tracer(i,j,k,nqni)*1.e-3.le.qmin 
          end do
        end do
      end do

!------------------------------------------------------------------------
!    call subroutine adjust_condensate to conservatively fill ql if needed.
!------------------------------------------------------------------------
      call adjust_condensate (ql_too_small, Tend_mp%q_tnd(:,:,:,nql), &
            Tend_mp%qtnd, Tend_mp%ttnd, Input_mp%tracer(:,:,:,nql), &
            HLV, Cloud_state%ql_upd)

!------------------------------------------------------------------------
!    adjust cloud droplet numbers as needed when those fields are being 
!    predicted. if droplet number is not being predicted, values were 
!    set at allocation.  
!------------------------------------------------------------------------
      call adjust_particle_number (ql_too_small,   &
                   Tend_mp%q_tnd(:,:,:,nqn), Input_mp%tracer(:,:,:,nqn), &
                                                      Cloud_state%qn_upd)

!------------------------------------------------------------------------
!    save diagnostics defining the cloud liquid and cloud particle
!    number filling amount. 
!------------------------------------------------------------------------
      do k = 1,kdim
        do j = 1,jdim
          do i = 1,idim
            if (ql_too_small(i,j,k)) then
              if (Lsdiag_mp_control%diag_id%qdt_liquid_init > 0)  &
                   Lsdiag_mp%diag_4d(i,j,k,  &
                           Lsdiag_mp_control%diag_pt%qdt_liquid_init)  =  &
                                       Input_mp%tracer(i,j,k,nql)/dtcloud
              if (Lsdiag_mp_control%diag_id%qndt_fill  +   &
                  Lsdiag_mp_control%diag_id%qn_fill_col + &
                  Lsdiag_mp_control%diag_id%qldt_fill +   &
                  Lsdiag_mp_control%diag_id%ql_fill_col > 0 )    &
                    Lsdiag_mp%diag_4d(i,j,k,  &
                             Lsdiag_mp_control%diag_pt%qndt_fill) = &
                                      -Input_mp%tracer(i,j,k,nqn)/dtcloud 

            endif
          end do
        end do
      end do

!------------------------------------------------------------------------
!    call subroutine adjust_condensate to conservatively fill qi if needed.
!------------------------------------------------------------------------
      call adjust_condensate (qi_too_small, Tend_mp%q_tnd(:,:,:,nqi), &
            Tend_mp%qtnd, Tend_mp%ttnd, Input_mp%tracer(:,:,:,nqi), &
            HLS, Cloud_state%qi_upd)

!------------------------------------------------------------------------
!    save diagnostics defining the cloud ice  filling amount. 
!------------------------------------------------------------------------
      do k = 1,kdim
        do j = 1,jdim
          do i = 1,idim
            if (qi_too_small(i,j,k)) then
              if (Lsdiag_mp_control%diag_id%qdt_ice_init > 0)   &
                     Lsdiag_mp%diag_4d(i,j,k,  &
                            Lsdiag_mp_control%diag_pt%qdt_ice_init ) =   &
                                        Input_mp%tracer(i,j,k,nqi)/dtcloud
            endif
          end do
        end do
      end do

!------------------------------------------------------------------------
!    adjust ice particle numbers as needed when those fields are being 
!    predicted. 
!------------------------------------------------------------------------
      if (do_ice_num) then
        call adjust_particle_number (qi_too_small,   &
                   Tend_mp%q_tnd(:,:,:,nqni), Input_mp%tracer(:,:,:,nqni),&
                                                      Cloud_state%qni_upd)
        do k = 1,kdim
          do j = 1,jdim
            do i = 1,idim
              if (qi_too_small(i,j,k)) then

!------------------------------------------------------------------------
!    save a diagnostic defining the ice crystal number filling amount. 
!------------------------------------------------------------------------

                if (Lsdiag_mp_control%diag_id%qnidt_fill  +    &
                    Lsdiag_mp_control%diag_id%qni_fill_col > 0 )   &
                   Lsdiag_mp%diag_4d(i,j,k,   &
                             Lsdiag_mp_control%diag_pt%qnidt_fill) =   &
                                      -Input_mp%tracer(i,j,k,nqni)/dtcloud
               endif
             end do
           end do
         end do
       endif


!-----------------------------------------------------------------------



end subroutine impose_realizability_clubb



!########################################################################

subroutine detailed_diagnostics (     &
                   is, ie, js, je, Time, n_diag_4d, diag_4d, diag_4d_kp1, &
                   diag_3d, diag_pt, diag_id, C2ls_mp, Input_mp,   &
                   Atmos_state, Cloud_state, Particles, Precip_state,   &
                   Cloud_processes, Tend_mp, Removal_mp)

!-----------------------------------------------------------------------
!     subroutine detailed_diagnostics computes tendencies and imbalances 
!     for the equations related to the large-scale clouds (cloud area,
!     cloud liquid, cloud ice, cloud droplet number, cloud ice particle
!     number, water vapor, temperature, precip, rain and snow), outputs 
!     any debug variables requested to either stdout or an ascii data file
!     (option to write to data file is being removed), computes column 
!     integrated diagnostics, and calls lscloud_netcdf to output the 
!     relevant netcdf diagnostics. 

!     note that budgets for CLUBB will not balance (more work to be done),
!     but should balance in non-CLUBB cases.

!     note that for the budget imbalance terms to be valid ALL terms in
!     the particular budget equation must be present in the diag_table.
!RSH:  DO THIS OR NOT ???
!     a diag_table with all such terms is available with code checkout at
!     /atmos_param/strat_cloud/diag_table_sc.
!-----------------------------------------------------------------------

!------------------------------------------------------------------------
type(mp_tendency_type),     intent(inout) :: Tend_mp
type(mp_removal_type),      intent(inout) :: Removal_mp
type(time_type),            intent (in)   :: Time
integer,                    intent (in)   :: is,ie,js,je
integer,                    intent(in)    :: n_diag_4d
real, dimension(:,:,:,0:),  intent(inout) :: diag_4d, diag_4d_kp1
real, dimension(:,:,0:),    intent(inout) :: diag_3d
type(diag_pt_type),         intent(inout) :: diag_pt
type(diag_id_type),         intent(inout) :: diag_id
type(atmos_state_type),     intent(inout) :: Atmos_state
type(mp_input_type),        intent(inout) :: Input_mp   
type(mp_conv2ls_type),      intent(inout) :: C2ls_mp   
type(cloud_state_type),     intent(inout) :: Cloud_state
type(particles_type),       intent(inout) :: Particles
type(precip_state_type),    intent(inout) :: Precip_state
type(cloud_processes_type), intent(inout) :: Cloud_processes


!------------------------------------------------------------------------
!---local variables------------------------------------------------------

!------------------------------------------------------------------------
!  variables used in calculation of particle number diagnostics:
      real, dimension(size(Input_mp%tin,1),size(Input_mp%tin,2))  ::    &
              N3D_col, N3Di_col, N3D_col250, gb_N3D_col, gb_N3Di_col
      real, dimension(size(Input_mp%tin,1),size(Input_mp%tin,2), &
                                     size(Input_mp%tin,3))  ::    &
              dum, qa_new, qn_new, ql_new, qi_new, qni_new

!------------------------------------------------------------------------
      integer :: idim, jdim, kdim
      integer :: i, j ,k, nn
      logical :: used

!-------------------------------------------------------------------------
!    define spatial dimensions.                
!-------------------------------------------------------------------------
      idim = SIZE(Input_mp%tin, 1)
      jdim = SIZE(Input_mp%tin, 2)
      kdim = SIZE(Input_mp%tin, 3) 

!-----------------------------------------------------------------------
!    define diagnostics for the time tendencies of the model prognostic
!    variables.
!-----------------------------------------------------------------------
      if (diag_id%SA3d + diag_id%SA2d > 0) then
        diag_4d(:,:,:,diag_pt%SA3d) = Tend_mp%q_tnd(:,:,:,nqa)*inv_dtcloud
      endif
      if (diag_id%ST3d + diag_id%ST2d > 0) then
        diag_4d(:,:,:,diag_pt%ST3d) = Tend_mp%ttnd(:,:,:)*inv_dtcloud
      endif
      if (diag_id%SQ3d + diag_id%SQ2d > 0) then
        diag_4d(:,:,:,diag_pt%SQ3d) = Tend_mp%qtnd(:,:,:)*inv_dtcloud
      endif
      if (diag_id%SL3d + diag_id%SL2d > 0) then
        diag_4d(:,:,:,diag_pt%SL3d) =  Tend_mp%q_tnd(:,:,:,nql)*inv_dtcloud
      endif
      if (diag_id%SI3d + diag_id%SI2d > 0) then
        diag_4d(:,:,:,diag_pt%SI3d) = Tend_mp%q_tnd(:,:,:,nqi)*inv_dtcloud
      endif
      if (diag_id%SN3d + diag_id%SN2d > 0) then
        diag_4d(:,:,:,diag_pt%SN3d) = Cloud_state%SN_out*inv_dtcloud
      endif
      if (diag_id%SNI3d + diag_id%SNI2d > 0) then
        diag_4d(:,:,:,diag_pt%SNI3d) = Cloud_state%SNI_out*inv_dtcloud
      endif
        
!-----------------------------------------------------------------------
!    define diagnostics reflecting any budget imbalances between the
!    total tendency to be applied on the time step and the individual terms
!    making up that tendency, ie, is the actual change fully accounted for
!    by the current budget terms? all terms on the rt side of these
!    expressions must be in the diag_table for the left side term to be
!    meaningful.
!-----------------------------------------------------------------------
      if (diag_id%SA_imb + diag_id%SA_imb_col > 0) then
        diag_4d(:,:,:,diag_pt%SA_imb) =    &
             diag_4d(:,:,:,diag_pt%SA3d) -   (       &
                diag_4d(:,:,:,diag_pt%qadt_lsform)   &
             +  diag_4d(:,:,:,diag_pt%qadt_lsdiss)   &
             +  diag_4d(:,:,:,diag_pt%qadt_rhred)    &
             +  diag_4d(:,:,:,diag_pt%qadt_eros)     &
             +  diag_4d(:,:,:,diag_pt%qadt_fill)     &
             +  diag_4d(:,:,:,diag_pt%qadt_super)    &
             +  diag_4d(:,:,:,diag_pt%qadt_destr)    &
             +  diag_4d(:,:,:,diag_pt%qadt_limits)   &
             +  diag_4d(:,:,:,diag_pt%qadt_ahuco)    &
                                                          )
      endif
      if (diag_id%SL_imb + diag_id%SL_imb_col > 0) then
        diag_4d(:,:,:,diag_pt%SL_imb) =  &
             diag_4d(:,:,:,diag_pt%SL3d) -   (            &
                diag_4d(:,:,:,diag_pt%qldt_cond )         &
              + diag_4d(:,:,:,diag_pt%qldt_evap  )        &
              + diag_4d(:,:,:,diag_pt%qldt_eros  )        &
              + diag_4d(:,:,:,diag_pt%qldt_berg)          &
              + diag_4d(:,:,:,diag_pt%qldt_freez )        &
              + diag_4d(:,:,:,diag_pt%liq_adj    )        &
              + diag_4d(:,:,:,diag_pt%qldt_rime  )        &
              + diag_4d(:,:,:,diag_pt%qldt_accr)          &
              + diag_4d(:,:,:,diag_pt%qldt_auto)          &
              + diag_4d(:,:,:,diag_pt%qldt_fill  )        &
              + diag_4d(:,:,:,diag_pt%qldt_destr )        &
              + diag_4d(:,:,:,diag_pt%qldt_freez2)        &
              + diag_4d(:,:,:,diag_pt%qldt_sedi  )        &
              + diag_4d(:,:,:,diag_pt%qldt_accrs)         &
              + diag_4d(:,:,:,diag_pt%qldt_bergs)         &
              + diag_4d(:,:,:,diag_pt%qldt_HM_splinter)   &
              - diag_4d(:,:,:,diag_pt%qidt_melt2 )        &
              - diag_4d(:,:,:,diag_pt%qidt_accrs)         &
              - diag_4d(:,:,:,diag_pt%qdt_cleanup_liquid) &    
                                                             )
      endif
      if (diag_id%SI_imb + diag_id%SI_imb_col > 0) then
        diag_4d(:,:,:,diag_pt%SI_imb) =     &
             diag_4d(:,:,:,diag_pt%SI3d) -   (          &
              - diag_4d(:,:,:,diag_pt%qldt_berg)        &
              - diag_4d(:,:,:,diag_pt%qldt_freez )      &
              - diag_4d(:,:,:,diag_pt%qldt_rime  )      &
              - diag_4d(:,:,:,diag_pt%qldt_freez2)      &
              - diag_4d(:,:,:,diag_pt%qldt_HM_splinter) &
              + diag_4d(:,:,:,diag_pt%qidt_dep  )       &
              + diag_4d(:,:,:,diag_pt%qidt_subl  )      &
              + diag_4d(:,:,:,diag_pt%qidt_fall  )      &
              + diag_4d(:,:,:,diag_pt%qidt_eros  )      &
              + diag_4d(:,:,:,diag_pt%qidt_melt  )      &
              + diag_4d(:,:,:,diag_pt%qidt_melt2 )      &
              + diag_4d(:,:,:,diag_pt%qidt_fill  )      &
              + diag_4d(:,:,:,diag_pt%qidt_destr )      &
              + diag_4d(:,:,:,diag_pt%qidt_qvdep )      &
              + diag_4d(:,:,:,diag_pt%qidt_auto)        &
              + diag_4d(:,:,:,diag_pt%qidt_accr)        &
              + diag_4d(:,:,:,diag_pt%qidt_accrs)       &
              + diag_4d(:,:,:,diag_pt%ice_adj    )      &
              - diag_4d(:,:,:,diag_pt%qdt_cleanup_ice)  &
                                                          )
      endif
      if (diag_id%SN_imb + diag_id%SN_imb_col > 0) then
        diag_4d(:,:,:,diag_pt%SN_imb) =  &
             diag_4d(:,:,:,diag_pt%SN3d) -   ( &
                diag_4d(:,:,:,diag_pt%qndt_cond  )      &
              + diag_4d(:,:,:,diag_pt%qndt_evap  )      &
              + diag_4d(:,:,:,diag_pt%qndt_fill  )      &
              + diag_4d(:,:,:,diag_pt%qndt_berg  )      &
              + diag_4d(:,:,:,diag_pt%qndt_destr )      &
              + diag_4d(:,:,:,diag_pt%qndt_super )      &
              + diag_4d(:,:,:,diag_pt%qndt_freez )      &
              + diag_4d(:,:,:,diag_pt%qndt_sacws )      &
              + diag_4d(:,:,:,diag_pt%qndt_sacws_o)     &
              + diag_4d(:,:,:,diag_pt%qndt_eros  )      &
              + diag_4d(:,:,:,diag_pt%qndt_pra   )      &
              + diag_4d(:,:,:,diag_pt%qndt_auto  )      &
              + diag_4d(:,:,:,diag_pt%qndt_nucclim)     &
              + diag_4d(:,:,:,diag_pt%qndt_sedi )       &
              + diag_4d(:,:,:,diag_pt%qndt_melt)        &
              + diag_4d(:,:,:,diag_pt%qndt_ihom)        &
              + diag_4d(:,:,:,diag_pt%qndt_size_adj)    &
              + diag_4d(:,:,:,diag_pt%qndt_fill2)       &
              + diag_4d(:,:,:,diag_pt%qndt_contact_frz) &
              + diag_4d(:,:,:,diag_pt%qndt_cleanup)     &
              + diag_4d(:,:,:,diag_pt%qndt_cleanup2)    &
                                                            )
      endif
      if (diag_id%SNi_imb + diag_id%SNi_imb_col > 0) then
        diag_4d(:,:,:,diag_pt%SNi_imb) =     &
             diag_4d(:,:,:,diag_pt%SNi3d) -   ( &
                diag_4d(:,:,:,diag_pt%qnidt_fill )     &
              + diag_4d(:,:,:,diag_pt%qnidt_nnuccd)    &
              + diag_4d(:,:,:,diag_pt%qnidt_nsubi)     &
              + diag_4d(:,:,:,diag_pt%qnidt_nerosi)    &
              + diag_4d(:,:,:,diag_pt%qnidt_nprci)     &
              + diag_4d(:,:,:,diag_pt%qnidt_nprai)     &
              + diag_4d(:,:,:,diag_pt%qnidt_nucclim1)  &
              + diag_4d(:,:,:,diag_pt%qnidt_nucclim2)  &
              + diag_4d(:,:,:,diag_pt%qnidt_sedi  )    &
              + diag_4d(:,:,:,diag_pt%qnidt_melt  )    &
              + diag_4d(:,:,:,diag_pt%qnidt_size_adj ) &
              + diag_4d(:,:,:,diag_pt%qnidt_fill2  )   &
              + diag_4d(:,:,:,diag_pt%qnidt_super )    &
              + diag_4d(:,:,:,diag_pt%qnidt_ihom )     &
              + diag_4d(:,:,:,diag_pt%qnidt_destr )    &
              + diag_4d(:,:,:,diag_pt%qnidt_cleanup)   &
              + diag_4d(:,:,:,diag_pt%qnidt_cleanup2)  &
              + diag_4d(:,:,:,diag_pt%qnidt_nsacwi)    &
                                                           )
      endif
      if (diag_id%SQ_imb + diag_id%SQ_imb_col > 0) then
        diag_4d(:,:,:,diag_pt%SQ_imb) =     &
             diag_4d(:,:,:,diag_pt%SQ3d) -   (                &
              - diag_4d(:,:,:,diag_pt%qldt_cond  )            &
              - diag_4d(:,:,:,diag_pt%qldt_evap  )            &
              - diag_4d(:,:,:,diag_pt%qldt_eros  )            &
              - diag_4d(:,:,:,diag_pt%liq_adj    )            &
              - diag_4d(:,:,:,diag_pt%qldt_fill  )            &
              - diag_4d(:,:,:,diag_pt%qldt_destr )            &
              - diag_4d(:,:,:,diag_pt%qidt_dep  )             &
              - diag_4d(:,:,:,diag_pt%qidt_subl  )            &
              - diag_4d(:,:,:,diag_pt%qidt_eros  )            &
              - diag_4d(:,:,:,diag_pt%qidt_fill  )            &
              - diag_4d(:,:,:,diag_pt%qidt_destr )            &
              - diag_4d(:,:,:,diag_pt%qidt_qvdep )            &
              - diag_4d(:,:,:,diag_pt%ice_adj    )            &
              + diag_4d(:,:,:,diag_pt%rain_evap  )            &
              + diag_4d(:,:,:,diag_pt%qdt_sedi_ice2vapor)     & 
              + diag_4d(:,:,:,diag_pt%qdt_sedi_liquid2vapor)  &  
              + diag_4d(:,:,:,diag_pt%qdt_cleanup_ice)        &
              + diag_4d(:,:,:,diag_pt%qdt_cleanup_liquid)     &    
              + diag_4d(:,:,:,diag_pt%qdt_snow_sublim  )      &
              + diag_4d(:,:,:,diag_pt%qdt_snow2vapor    )     &
                                                                )
      endif
      if (diag_id%ST_imb + diag_id%ST_imb_col > 0) then
        diag_4d(:,:,:,diag_pt%ST_imb) =     &
             diag_4d(:,:,:,diag_pt%ST3d) -    (              &
              - HLF*diag_4d(:,:,:,diag_pt%qldt_berg)         &
              - HLF*diag_4d(:,:,:,diag_pt%qldt_freez )       &
              - HLF*diag_4d(:,:,:,diag_pt%qldt_rime  )       &
              - HLF*diag_4d(:,:,:,diag_pt%qldt_freez2)       &
              - HLF*diag_4d(:,:,:,diag_pt%qldt_accrs)        &
              - HLF*diag_4d(:,:,:,diag_pt%qldt_bergs)        &
              - HLF*diag_4d(:,:,:,diag_pt%qldt_HM_splinter)  &
              + HLF*diag_4d(:,:,:,diag_pt%qidt_melt  )       &
              + HLF*diag_4d(:,:,:,diag_pt%qidt_melt2 )       &
              + HLF*diag_4d(:,:,:,diag_pt%qidt_accrs)        &
              + HLF*diag_4d(:,:,:,diag_pt%rain_freeze)       &
              - HLF*diag_4d(:,:,:,diag_pt%srfrain_accrs )    &
              - HLF*diag_4d(:,:,:,diag_pt%srfrain_freez )    &
              - HLF*diag_4d(:,:,:,diag_pt%snow_melt)         &
              
              + HLV*diag_4d(:,:,:,diag_pt%qldt_cond  )           &
              + HLV*diag_4d(:,:,:,diag_pt%qldt_evap  )           &
              + HLV*diag_4d(:,:,:,diag_pt%qldt_eros  )           &
              + HLV*diag_4d(:,:,:,diag_pt%liq_adj    )           &
              + HLV*diag_4d(:,:,:,diag_pt%qldt_fill  )           &
              + HLV*diag_4d(:,:,:,diag_pt%qldt_destr )           &
              - HLV*diag_4d(:,:,:,diag_pt%rain_evap  )           &
              - HLV*diag_4d(:,:,:,diag_pt%qdt_sedi_liquid2vapor) &  
              - HLV*diag_4d(:,:,:,diag_pt%qdt_cleanup_liquid)    &      

              + HLS*diag_4d(:,:,:,diag_pt%qidt_dep  )          &
              + HLS*diag_4d(:,:,:,diag_pt%qidt_subl  )         &
              + HLS*diag_4d(:,:,:,diag_pt%qidt_eros  )         &
              + HLS*diag_4d(:,:,:,diag_pt%qidt_fill  )         &
              + HLS*diag_4d(:,:,:,diag_pt%qidt_destr )         &
              + HLS*diag_4d(:,:,:,diag_pt%qidt_qvdep )         &
              + HLS*diag_4d(:,:,:,diag_pt%ice_adj    )         &
              - HLS*diag_4d(:,:,:,diag_pt%qdt_sedi_ice2vapor)  & 
              - HLS*diag_4d(:,:,:,diag_pt%qdt_cleanup_ice)     &
              - HLS*diag_4d(:,:,:,diag_pt%qdt_snow_sublim  )   &
              - HLV*diag_4d(:,:,:,diag_pt%qdt_snow2vapor    )  &
                                                              )/CP_AIR 
      endif 

!------------------------------------------------------------------------
!    define diagnostics for the 3D rain and snow fields.
!------------------------------------------------------------------------
      if (diag_id%rain3d > 0) then 
        diag_4d_kp1(:,:,:,diag_pt%rain3d) = Removal_mp%rain3d(:,:,:) 
      endif 
      if (diag_id%snow3d > 0) then 
        diag_4d_kp1(:,:,:,diag_pt%snow3d) = Removal_mp%snow3d(:,:,:) 
      endif 
      if (diag_id%qrout + diag_id%qrout_col > 0) then
        diag_4d(:,:,:,diag_pt%qrout) = Precip_state%lsc_rain 
      end if
      if (diag_id%qsout  + diag_id%qsout_col > 0) then
        diag_4d(:,:,:,diag_pt%qsout) = Precip_state%lsc_snow 
      end if

!------------------------------------------------------------------------
!    define droplet and ice particle number diagnostics.
!------------------------------------------------------------------------
      if (diag_id%droplets > 0) then
        diag_4d(:,:,:,diag_pt%droplets) = Particles%N3D(:,:,:)
      end if
      if (diag_id%nice > 0) then
        diag_4d(:,:,:,diag_pt%nice) = Particles%N3Di(:,:,:)
      end if
      if (diag_id%droplets_wtd > 0) then
        diag_4d(:,:,:,diag_pt%droplets_wtd) =   &
                            Particles%N3D(:,:,:)*Cloud_state%ql_in(:,:,:)
      end if
      if (diag_id%ql_wt > 0) then
        diag_4d(:,:,:,diag_pt%ql_wt) = Cloud_state%ql_in(:,:,:)
      end if

!-------------------------------------------------------------------------
!    call lscloud_debug to output data to  a debug file, if requested.
!-------------------------------------------------------------------------
      call mpp_clock_begin (lscloud_debug_clock)
      call lscloud_debug (Tend_mp%ttnd,Tend_mp%qtnd, Cloud_state,  &
             Removal_mp,  Precip_state, Input_mp, C2ls_mp, Atmos_State)
      call mpp_clock_end (lscloud_debug_clock)

!------------------------------------------------------------------------
!    if it has been requested, output the number of columns in which 
!    microphysics was skipped due to negative h2o in column (only
!    relevant with "mg" microphysics.
!------------------------------------------------------------------------
      call mpp_clock_begin (lscloud_debug_clock)
      call output_refusals 
      call mpp_clock_end (lscloud_debug_clock)

!------------------------------------------------------------------------
!    generate column integrated diagnostics.
!------------------------------------------------------------------------
      do nn=1, n_diag_4d
        do k=kdim,1, -1
          diag_3d(:,:,nn) = diag_3d(:,:,nn) +   &
                                diag_4d(:,:,k,nn)*Input_mp%pmass(:,:,k)
        enddo
      enddo

!------------------------------------------------------------------------
!    define additional column diagnostics that are not simply the  
!    pressure-depth-weighted sum of the column values:
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!       1) cloud ice and cloud liquid fallout for NCAR microphysics: 
!------------------------------------------------------------------------
      if (Constants_lsc%do_mg_microphys .or.   &
                               Constants_lsc%do_ncar_microphys .or.   &
                                   Constants_lsc%do_mg_ncar_microphys) then
        if (diag_id%cld_liq_imb + diag_id%cld_liq_imb_col > 0) then
          diag_3d(:,:,diag_pt%cld_liq_imb) =   - ( &
                diag_3d(:,:,diag_pt%qldt_sedi )  +   &
                diag_3d(:,:,diag_pt%qdt_sedi_liquid2vapor)  ) 
        endif
        if (diag_id%cld_ice_imb + diag_id%cld_ice_imb_col > 0) then
          diag_3d(:,:,diag_pt%cld_ice_imb) =   - ( &
                diag_3d(:,:,diag_pt%qidt_fall )  +   &
                diag_3d(:,:,diag_pt%qdt_sedi_ice2vapor)  ) 
        endif
      endif

!------------------------------------------------------------------------
!       2) negative values of rain and snow at the surface
!------------------------------------------------------------------------
      IF ( diag_id%neg_rain > 0   ) &
          diag_3d(:,:,diag_pt%neg_rain) =   &
                          diag_4d(:,:,1,diag_pt%neg_rain) 
      IF ( diag_id%neg_snow > 0   ) &
          diag_3d(:,:,diag_pt%neg_snow) =   &
                          diag_4d(:,:,1,diag_pt%neg_snow) 

!------------------------------------------------------------------------
!       3) imbalance between rainfall and its source/sink terms
!------------------------------------------------------------------------
      if (diag_id%rain_imb > 0) then
        diag_3d(:,:,diag_pt%rain_imb) =    &
             Precip_state%surfrain(:,:)*inv_dtcloud -   &
                diag_3d(:,:,diag_pt%cld_liq_imb)   +  (  &
                diag_3d(:,:,diag_pt%qldt_accr)      &
              + diag_3d(:,:,diag_pt%qldt_auto )     &
              + diag_3d(:,:,diag_pt%qidt_melt  )    &
              + diag_3d(:,:,diag_pt%rain_evap)      &
              + diag_3d(:,:,diag_pt%rain_freeze)    &
              - diag_3d(:,:,diag_pt%srfrain_accrs)  &
              - diag_3d(:,:,diag_pt%srfrain_freez)  &
              + diag_3d(:,:,diag_pt%neg_rain)       &
              - diag_3d(:,:,diag_pt%snow_melt )     &
              + diag_3d(:,:,diag_pt%qdt_snow2vapor) &
                                                         )
      endif

!------------------------------------------------------------------------
!       4) imbalance between snowfall and its source/sink terms
!------------------------------------------------------------------------
      if (diag_id%snow_imb > 0) then
        if (Constants_lsc%do_rk_microphys) then
        diag_3d(:,:,diag_pt%snow_imb) =    &
             Precip_state%surfsnow(:,:)*inv_dtcloud -    &
               diag_3d(:,:,diag_pt%cld_ice_imb)  + (  &
                diag_3d(:,:,diag_pt%qldt_accrs)      &
              + diag_3d(:,:,diag_pt%qldt_bergs)      &
              + diag_3d(:,:,diag_pt%qidt_fall)       &
              + diag_3d(:,:,diag_pt%qidt_auto  )     &
              + diag_3d(:,:,diag_pt%qidt_accr  )     &
              - diag_3d(:,:,diag_pt%rain_freeze)     &
              + diag_3d(:,:,diag_pt%srfrain_accrs)   &
              + diag_3d(:,:,diag_pt%srfrain_freez)   &
              + diag_3d(:,:,diag_pt%snow_melt )      &
              + diag_3d(:,:,diag_pt%neg_snow)        &
              + diag_3d(:,:,diag_pt%qdt_snow_sublim) &
                                                         )
        else
        diag_3d(:,:,diag_pt%snow_imb) =    &
             Precip_state%surfsnow(:,:)*inv_dtcloud -    &
               diag_3d(:,:,diag_pt%cld_ice_imb)  + (  &
                diag_3d(:,:,diag_pt%qldt_accrs)      &
              + diag_3d(:,:,diag_pt%qldt_bergs)      &
              + diag_3d(:,:,diag_pt%qidt_auto  )     &
              + diag_3d(:,:,diag_pt%qidt_accr  )     &
              - diag_3d(:,:,diag_pt%rain_freeze)     &
              + diag_3d(:,:,diag_pt%srfrain_accrs)   &
              + diag_3d(:,:,diag_pt%srfrain_freez)   &
              + diag_3d(:,:,diag_pt%snow_melt )      &
              + diag_3d(:,:,diag_pt%neg_snow)        &
              + diag_3d(:,:,diag_pt%qdt_snow_sublim) &
                                                         )
        endif
      endif

!------------------------------------------------------------------------
!       5)  changes to sfc precip required for h2o balance in a column
!------------------------------------------------------------------------
      if (diag_id%rain_mass_conv > 0   ) &
          diag_3d(:,:,diag_pt%rain_mass_conv) =   &
                          diag_4d(:,:,1,diag_pt%rain_mass_conv) 
      if (diag_id%snow_mass_conv > 0   ) &
          diag_3d(:,:,diag_pt%snow_mass_conv) =   &
                          diag_4d(:,:,1,diag_pt%snow_mass_conv) 

!------------------------------------------------------------------------
!       6) in-cloud droplet column burden
!          fixed as per email from yim to rsh of 11/3/2011
!------------------------------------------------------------------------
      if (diag_id%droplets_col > 0 .or.   &
          diag_id%gb_droplets_col > 0  .or. &
          diag_id%droplets_col250 > 0 ) then
        if (do_liq_num) then
          N3D_col(:,:) = 0.
          N3D_col250(:,:) = 0.
          gb_N3D_col(:,:) = 0.
          do k =1,kdim
            do j=1,jdim
              do i=1,idim
                qa_new(i,j,k) = Cloud_state%qa_in(i,j,k) +    &
                                                 Tend_mp%q_tnd(i,j,k,nqa)
                ql_new(i,j,k) = Cloud_state%ql_in(i,j,k) +    &
                                                  Tend_mp%q_tnd(i,j,k,nql)
                qn_new(i,j,k) = Cloud_state%qn_in(i,j,k) +    &
                                                   Tend_mp%q_tnd(i,j,k,nqn)
! may change answers for clubb diagnostics:
                if (ql_new(i,j,k) > qmin .and. &
                    qa_new(i,j,k) > qmin .and. &
                    qn_new(i,j,k) > qmin ) then      
                  dum(i,j,k) = qn_new(i,j,k)*Input_mp%pmass(i,j,k)*1.e-4

!------------------------------------------------------------------------
!    count only columns with qa > 5% for N3d_col and N3d_col250. 
!------------------------------------------------------------------------
                  if (qa_new(i,j,k) > 0.05) then 
                    N3D_col(i,j) = N3D_col(i,j) +   &
                                         dum(i,j,k)/min(qa_new(i,j,k), 1.)
                    if (Input_mp%tin(i,j,k) +   &
                                 Tend_mp%ttnd(i,j,k) .ge. 250.) then
                      N3D_col250(i,j)  = N3D_col250(i,j) +   &
                                         dum(i,j,k)/min(qa_new(i,j,k), 1.)
                    endif
                  endif

!------------------------------------------------------------------------
!    count all columns with qa > qmin for gb_N3d_col. 
!------------------------------------------------------------------------
                  gb_N3D_col(i,j) = gb_N3D_col(i,j) + dum(i,j,k)
                endif
              end do
            end do
          end do
          diag_3d(:,:,diag_pt%droplets_col) = N3D_col
          diag_3d(:,:,diag_pt%droplets_col250) = N3D_col250
          diag_3d(:,:,diag_pt%gb_droplets_col) = gb_N3D_col
        endif
      endif

!------------------------------------------------------------------------
!       7) yim:  in-cloud ice crystal column burden
!------------------------------------------------------------------------
      if (diag_id%nice_col > 0 .or. diag_id%gb_nice_col > 0) then
        if (do_ice_num) then
          N3Di_col(:,:) = 0.
          gb_N3Di_col(:,:) = 0.
          do k=1,kdim
            do j=1,jdim
              do i=1,idim
                qa_new(i,j,k) = Cloud_state%qa_in(i,j,k) +   &
                                 Tend_mp%q_tnd(i,j,k,nqa)
                qi_new(i,j,k) = Cloud_state%qi_in(i,j,k) +   &
                                 Tend_mp%q_tnd(i,j,k,nqi)
                qni_new(i,j,k) = Cloud_state%qni_in(i,j,k) +  &
                                 Tend_mp%q_tnd(i,j,k,nqni)
! may change answers for clubb diagnostics:
                if (qi_new(i,j,k) > qmin .and. &
                    qa_new(i,j,k) > qmin .and. &
                    qni_new(i,j,k)  > qmin ) then
                  dum(i,j,k) =  qni_new(i,j,k)*Atmos_state%airdens(i,j,k)*&
                                        Input_mp%pmass(i,j,k)*1.e-6

!------------------------------------------------------------------------
!   count only columns with qa > 5%  for N3Di_col.
!------------------------------------------------------------------------
                  if (qa_new(i,j,k) > 0.05) then 
                    N3Di_col(i,j) = N3Di_col(i,j) +   &
                                        dum(i,j,k) /min(qa_new(i,j,k), 1.)
                  endif

!------------------------------------------------------------------------
!   count all columns with qa > qmin  for gb_N3Di_col.
!------------------------------------------------------------------------
                  gb_N3Di_col(i,j) = gb_N3Di_col(i,j) + dum(i,j,k)
                endif
              end do
            end do
          end do
          diag_3d(:,:,diag_pt%nice_col) = N3Di_col
          diag_3d(:,:,diag_pt%gb_nice_col) = gb_N3Di_col
        endif
      endif

!-------------------------------------------------------------------------
!    call lscloud_netcdf to output the requested netcdf diagnostics.
!-------------------------------------------------------------------------
      call mpp_clock_begin (lscloud_netcdf_clock)
      call lscloud_netcdf (diag_id, diag_pt, diag_4d, diag_4d_kp1,   &
                                            diag_3d, Time, is, js, kdim )
      call mpp_clock_end (lscloud_netcdf_clock)

!-------------------------------------------------------------------------
!    output diagnostics related to the fraction of ice with ice-forming
!    nuclei.
!-------------------------------------------------------------------------
      used = send_data (id_f_snow_berg, Cloud_processes%f_snow_berg,  &
                                                          Time, is, js, 1)
      used = send_data (id_f_snow_berg_cond,   &
                          Cloud_processes%f_snow_berg, Time, is, js, 1,   &
                               mask = Cloud_processes%f_snow_berg /= 0.0)
      used = send_data (id_f_snow_berg_wtd,   &
               Cloud_processes%f_snow_berg(:,:,:)*  &
               (Removal_mp%rain3d(:,:,2:) + Removal_mp%snow3d(:,:,2:)), &
                Time, is, js, 1, mask = Cloud_processes%f_snow_berg /= 0.0)

!---------------------------------------------------------------------
!     precipitating condensate from prognostic clouds (scaling factor
!     for f_snow_berg_wtd)
!---------------------------------------------------------------------
      used = send_data(id_lscale_precip3d,   &
              (Removal_mp%snow3d(:,:,2:) + Removal_mp%rain3d(:,:,2:)), &
                Time, is, js, 1, mask = Cloud_processes%f_snow_berg /= 0.0)

!-----------------------------------------------------------------------


end subroutine detailed_diagnostics


!#######################################################################

subroutine update_fields_and_tendencies (  &
                   is, ie, js, je, Time, dtinv, Input_mp, Cld_props,  &
                                         Precip_state, Tend_mp, Output_mp)

!------------------------------------------------------------------------
!    subroutine update_fields_and_tendencies updates the model variables
!    (T, q, ql, qi, qa, qn, qni) with the computed tendencies from the
!    large-scale prognostic cloud scheme. it also outputs these updated
!    fields if desired. additionally, time increments due to prognostic
!    large-scale clouds are converted to time tendencies, and added to
!    those arrays accumulating the total tendency terms from 
!    moist_processes, for return to the calling routine. also, the output
!    from moist_processes that will be needed by the radiation package is
!    collected and bundled in the Cld_props derived type for return to
!    physics_driver_mod.
!------------------------------------------------------------------------

integer,                     intent(in)     :: is, ie, js, je
type(time_type),             intent(in)     :: Time
type(mp_input_type),         intent(inout)  :: Input_mp
type(cloud_scheme_data_type), intent(inout) :: Cld_props
type(precip_state_type),     intent(inout)  :: Precip_state
type(mp_tendency_type),      intent(inout)  :: Tend_mp
type(mp_output_type),        intent(inout)  :: Output_mp
real,                        intent(in)     :: dtinv

!----------------------------------------------------------------------
!   local variables:

      real, dimension(size(Input_mp%tin,1), size(Input_mp%tin,2), &
                         size(Input_mp%tin,3))  :: rh
      logical :: used


!----------------------------------------------------------------------
!    update the cloud liquid, ice, area, droplet number and crystal number
!    (if relevant) and the temperature and specific humidity fields.
!----------------------------------------------------------------------
      if (doing_prog_clouds) then
        Input_mp%tracer(:,:,:,nql) = Input_mp%tracer(:,:,:,nql) +   &
                                                 Tend_mp%q_tnd(:,:,:,nql)  
        Input_mp%tracer(:,:,:,nqi) = Input_mp%tracer(:,:,:,nqi) +   &
                                                 Tend_mp%q_tnd(:,:,:,nqi)
        Input_mp%tracer(:,:,:,nqa) = Input_mp%tracer(:,:,:,nqa) +   &
                                                 Tend_mp%q_tnd(:,:,:,nqa)
        if (do_liq_num)   &
          Input_mp%tracer(:,:,:,nqn) =  Input_mp%tracer(:,:,:,nqn) +   &
                                                 Tend_mp%q_tnd(:,:,:,nqn)
        if (do_ice_num)   &
          Input_mp%tracer(:,:,:,nqni) = Input_mp%tracer(:,:,:,nqni) +  &
                                                 Tend_mp%q_tnd(:,:,:,nqni)
      endif

      Input_mp%tin = Input_mp%tin + Tend_mp%ttnd 
      Input_mp%qin = Input_mp%qin + Tend_mp%qtnd
        
!-----------------------------------------------------------------------
!    output the basic variable fields after the above large-scale cloud 
!    tendencies have been added.
!-----------------------------------------------------------------------
      used = send_data (id_qvout, Input_mp%qin, Time, is, js, 1)
      if (doing_prog_clouds) then
        used = send_data (id_qaout, Input_mp%tracer(:,:,:,nqa),   &
                                                        Time, is, js, 1)
        used = send_data (id_qlout, Input_mp%tracer(:,:,:,nql),   &
                                                        Time, is, js, 1)
        used = send_data (id_qiout, Input_mp%tracer(:,:,:,nqi),   &
                                                        Time, is, js, 1)

        if (do_liq_num) &
          used = send_data (id_qnout, Input_mp%tracer(:,:,:,nqn), &
                                                        Time, is, js, 1)

        if (do_ice_num) &
          used = send_data (id_qniout, Input_mp%tracer(:,:,:,nqni),   &
                                                        Time, is, js, 1)
      endif

!----------------------------------------------------------------------
!    convert time increments to time tendencies.
!----------------------------------------------------------------------
      Tend_mp%ttnd = Tend_mp%ttnd*dtinv 
      Tend_mp%qtnd = Tend_mp%qtnd*dtinv
      Precip_state%surfrain = Precip_state%surfrain*dtinv 
      Precip_state%surfsnow = Precip_state%surfsnow*dtinv
      if (doing_prog_clouds) then
        Tend_mp%q_tnd(:,:,:,nql) = Tend_mp%q_tnd(:,:,:,nql)*dtinv
        Tend_mp%q_tnd(:,:,:,nqi) = Tend_mp%q_tnd(:,:,:,nqi)*dtinv
        Tend_mp%q_tnd(:,:,:,nqa) = Tend_mp%q_tnd(:,:,:,nqa)*dtinv
        if (do_liq_num) Tend_mp%q_tnd(:,:,:,nqn) =   &
                                          Tend_mp%q_tnd(:,:,:,nqn)*dtinv
        if (do_ice_num) Tend_mp%q_tnd(:,:,:,nqni) =   &
                                          Tend_mp%q_tnd(:,:,:,nqni)*dtinv
      endif
   
!----------------------------------------------------------------------
!    update the total tendency terms (temperature, vapor specific 
!    humidity, cloud liquid, cloud ice, cloud area, droplet number,
!    ice particle number, liquid precip, frozen precip) with the 
!    contributions from the prognostic large-scale cloud  scheme.
!----------------------------------------------------------------------
      Output_mp%tdt = Output_mp%tdt + Tend_mp%ttnd 
      Output_mp%rdt(:,:,:,1) = Output_mp%rdt(:,:,:,1) + Tend_mp%qtnd
      if (doing_prog_clouds) then
        Output_mp%rdt(:,:,:,nql) = Output_mp%rdt(:,:,:,nql) +   &
                                                 Tend_mp%q_tnd(:,:,:,nql)
        Output_mp%rdt(:,:,:,nqi) = Output_mp%rdt(:,:,:,nqi) +   &
                                                 Tend_mp%q_tnd(:,:,:,nqi)
        Output_mp%rdt(:,:,:,nqa) = Output_mp%rdt(:,:,:,nqa) +   &
                                                 Tend_mp%q_tnd(:,:,:,nqa)
        if (do_liq_num)   &
                Output_mp%rdt(:,:,:,nqn) = Output_mp%rdt(:,:,:,nqn) +   &
                                                 Tend_mp%q_tnd(:,:,:,nqn)
        if (do_ice_num)   &
                Output_mp%rdt(:,:,:,nqni) = Output_mp%rdt(:,:,:,nqni) + &
                                                 Tend_mp%q_tnd(:,:,:,nqni)
      endif
      Output_mp%lprec = Output_mp%lprec + Precip_state%surfrain
      Output_mp%fprec = Output_mp%fprec + Precip_state%surfsnow

!----------------------------------------------------------------------
!    save the large-scale cloud fields needed for use by the
!    radiation package.
!----------------------------------------------------------------------
      if (doing_prog_clouds) then
        Cld_props%cloud_area = Input_mp%tracer(:,:,:,nqa)
        Cld_props%liquid_amt = Input_mp%tracer(:,:,:,nql)
        Cld_props%ice_amt    = Input_mp%tracer(:,:,:,nqi)
        if (do_liq_num)   &
                Cld_props%droplet_number = Input_mp%tracer(:,:,:,nqn)
        if (do_ice_num)   &
                Cld_props%ice_number = Input_mp%tracer(:,:,:,nqni)
        Cld_props%snow      = Precip_state%lsc_snow
        Cld_props%rain      = Precip_state%lsc_rain
        Cld_props%snow_size = Precip_state%lsc_snow_size
        Cld_props%rain_size = Precip_state%lsc_rain_size
      endif

!--------------------------------------------------------------------
!    if rh_clouds is active, call rh_calc to determine the grid box
!    relative humidity. call rh_clouds_sum to pass this field to 
!    rh_clouds_mod so it may be used to determine the grid boxes which
!    will contain clouds for the radiation package.
!---------------------------------------------------------------------
      if (do_lsc .and. do_rh_clouds) then            
        call rh_calc (Input_mp%pfull, Input_mp%tin, Input_mp%qin, rh, &
                        do_simple)
        call rh_clouds_sum (is, js, rh)
      endif

!----------------------------------------------------------------------



end subroutine update_fields_and_tendencies


!#########################################################################

subroutine compute_ls_wetdep (   &
                  is, js, Time, Tend_mp, C2ls_mp, Input_mp, Removal_mp,   &
                                 f_snow_berg, dt, Output_mp, Precip_state)

!---------------------------------------------------------------------
!    calculate the wet deposition associated with the large scale 
!    condensation. 
!---------------------------------------------------------------------

integer,                 intent(in)    :: is, js
type(time_type),         intent(in)    :: Time
type(mp_tendency_type),  intent(inout) :: Tend_mp
type(mp_conv2ls_type),   intent(inout) :: C2ls_mp
type(mp_input_type),     intent(in)    :: Input_mp
type(mp_removal_type),   intent(inout) :: Removal_mp
real, dimension(:,:,:),  intent(in)    :: f_snow_berg
real,                    intent(in)    :: dt
type(mp_output_type),    intent(inout) :: Output_mp
type(precip_state_type), intent(inout) :: Precip_state

!------------------------------------------------------------------------
!   local variables:

      real, dimension(size(Input_mp%t,1),size(Input_mp%t,2),  &
                                          size(Input_mp%t,3))  ::   &
                                                    so2_so4_evap, so2_so4
      integer :: n
      integer :: kx
      logical :: used

      kx = size(Input_mp%t,3)

!---------------------------------------------------------------------
!    define the water tendency (vapor, condensate) from the large-scale 
!    cloud / condensation scheme.
!---------------------------------------------------------------------
      Tend_mp%qtnd_wet = Tend_mp%qtnd
      if (doing_prog_clouds) then
        Tend_mp%qtnd_wet = Tend_mp%qtnd_wet + Tend_mp%q_tnd(:,:,:,nql) +  &
                                              Tend_mp%q_tnd(:,:,:,nqi)

!-----------------------------------------------------------------------
!    sum up the precipitation formed over timestep.
!-----------------------------------------------------------------------
        if (Constants_lsc%do_lin_cld_microphys) then
          C2ls_mp%cloud_wet = Input_mp%tracer(:,:,:,nqr) +   &
                              Input_mp%tracer(:,:,:,nqs) +   &
                              Input_mp%tracer(:,:,:,nqg)
        else
          C2ls_mp%cloud_wet(:,:,:) =   &
            Removal_mp%rain3d(:,:,2:kx+1) - Removal_mp%rain3d(:,:,1:kx) + &
            Removal_mp%snow3d(:,:,2:kx+1) - Removal_mp%snow3d(:,:,1:kx)

!------------------------------------------------------------------------
!    convert from kg/m2/s to kg/kg
!------------------------------------------------------------------------
          C2ls_mp%cloud_wet(:,:,:) = C2ls_mp%cloud_wet(:,:,:)*dt/   &
                      Input_mp%pmass(:,:,:) 
        endif

!-----------------------------------------------------------------------
!    add the cloud amount at end of timestep.
!-----------------------------------------------------------------------
        C2ls_mp%cloud_wet(:,:,:) = C2ls_mp%cloud_wet(:,:,:) +   &
                                   Input_mp%tracer(:,:,:,nql) +  &
                                   Input_mp%tracer(:,:,:,nqi)

!-----------------------------------------------------------------------
!    define the large-scale cloud fraction, making sure it is between
!    0 and 1.
!-----------------------------------------------------------------------
        C2ls_mp%cloud_frac(:,:,:) =    &
                   max( min( Input_mp%tracer(:,:,:,nqa), 1. ), 0. )

!-----------------------------------------------------------------------
!    define values for the non-prognostic cloud case.
!-----------------------------------------------------------------------
      else
        C2ls_mp%cloud_wet(:,:,:) = 0.5e-3
        C2ls_mp%cloud_frac(:,:,:) = 1.
      end if

!----------------------------------------------------------------------
!    initialize field to hold so2--> so4 tendency returned from
!    subroutine wet_deposition.
!----------------------------------------------------------------------
      so2_so4_evap = 0.

!----------------------------------------------------------------------
!    loop over each tracer, skipping the cloud tracers, calling the
!    wet deposition routine to compute the amount of that tracer removed 
!    from the column in each layer (Tend_mp%wetdeptnd(i,j,k)) and in the
!    column (Removal_mp%ls_wetdep(i,j,n)).
!----------------------------------------------------------------------
      do n=1,size(Output_mp%rdt,4)
        if ( n /=              nsphum .and. &
             n /=              nql .and.  &
             n /=              nqi .and.  &
             n /=              nqa .and.  &
             n /=              nqn .and.  &
             n /=              nqni       &
                                       ) then
          Tend_mp%wetdeptnd(:,:,:) = 0.0
          call wet_deposition (        &
              n, Input_mp%t, Input_mp%pfull, Input_mp%phalf,   &
              Input_mp%zfull, Input_mp%zhalf, Precip_state%surfrain,  &
              Precip_state%surfsnow, Tend_mp%qtnd_wet, C2ls_mp%cloud_wet, &
              C2ls_mp%cloud_frac, f_snow_berg, Removal_mp%rain3d,  &
              Removal_mp%snow3d, Input_mp%tracer(:,:,:,n),  &
              Tend_mp%wetdeptnd, Time, 'lscale', is, js, dt,  &
              sum_wdep_out=Removal_mp%ls_wetdep(:,:,n),  &
              so2_so4_out = so2_so4(:,:,:) )

!-----------------------------------------------------------------------
!    add the wet deposition tendency for the tracer to the accumulated
!    total tracer tendency at each level.
!-----------------------------------------------------------------------
          Output_mp%rdt (:,:,:,n) =  Output_mp%rdt(:,:,:,n) -   &
                                                  Tend_mp%wetdeptnd(:,:,:)

!-----------------------------------------------------------------------
!    add the large-scale wet deposition tendency for the tracer to the 
!    previously-obtained convective wet deposition tendency 
!    (C2ls_mp%wet_data).
!-----------------------------------------------------------------------
          C2ls_mp%wet_data(:,:,:,n) = C2ls_mp%wet_data(:,:,:,n) +   &
                                                  Tend_mp%wetdeptnd(:,:,:)

          if ( n .eq. nso2 ) then
            so2_so4_evap = so2_so4
          end if

!------------------------------------------------------------------------
!    output the total wet deposition diagnostic for this tracer, 
!    if desired.
!------------------------------------------------------------------------
          used = send_data( id_wet_deposition(n),   &
                                   C2ls_mp%wet_data(:,:,:,n), Time,   &
                                                       is_in=is,js_in=js )
        end if
      end do

!-----------------------------------------------------------------------
!    correct so2 and so4 tendency. so2 is converted to so4.
!-----------------------------------------------------------------------

      if ( nso2 .ne.  NO_TRACER .and. nso4 .ne. NO_TRACER) then
        Output_mp%rdt(:,:,:,nso2) = Output_mp%rdt(:,:,:,nso2) -   &
                                                             so2_so4_evap
        Output_mp%rdt(:,:,:,nso4) = Output_mp%rdt(:,:,:,nso4) +   &
                                                             so2_so4_evap
      end if

!-----------------------------------------------------------------------


end subroutine compute_ls_wetdep


!#######################################################################

subroutine basic_diagnostics (is, js, Time, Tend_mp, Input_mp,    &
                                                 Precip_state, Removal_mp)

integer,                 intent(in)    :: is, js
type(time_type),         intent(in)    :: Time
type(mp_tendency_type),  intent(inout) :: Tend_mp
type(mp_input_type),     intent(inout) :: Input_mp
type(precip_state_type), intent(inout) :: Precip_state
type(mp_removal_type),   intent(inout) :: Removal_mp


!------------------------------------------------------------------------
!   local variables:

      real, dimension(size(Input_mp%t,1),size(Input_mp%t,2)) :: temp_2d

      logical :: used

!---------------------------------------------------------------------
!    output large-scale LWP and IWP.
!---------------------------------------------------------------------
      if (id_LWP > 0 .and. doing_prog_clouds) &
         call column_diag    &
             (id_LWP, is, js, Time, Input_mp%tracer(:,:,:,nql), 1.0,   &
                                                          Input_mp%pmass)
      if (id_IWP > 0 .and. doing_prog_clouds) &
          call column_diag   &
             (id_IWP, is, js, Time, Input_mp%tracer(:,:,:,nqi), 1.0,   &
                                                           Input_mp%pmass)

!---------------------------------------------------------------------
!    temperature change due to large-scale clouds:
!---------------------------------------------------------------------
      used = send_data (id_tdt_ls, Tend_mp%ttnd(:,:,:), Time, is, js, 1)
      used = send_cmip_data_3d (ID_tntscp, Tend_mp%ttnd, Time, is, js, 1)!, rmask=mask)

!---------------------------------------------------------------------
!    dry static energy tendency due to large-scale clouds:
!---------------------------------------------------------------------
      if (id_t_ls_col > 0) &
           call column_diag (id_t_ls_col, is, js, Time,   &
                              Tend_mp%ttnd(:,:,:), CP_AIR, Input_mp%pmass)

!---------------------------------------------------------------------
!    water vapor path tendency due to large-scale clouds:
!---------------------------------------------------------------------
      if (id_q_ls_col > 0) &
           call column_diag(id_q_ls_col, is, js, Time,    &
                                Tend_mp%qtnd(:,:,:), 1.0, Input_mp%pmass) 
 
!---------------------------------------------------------------------
!    specific humidity change due to large-scale clouds:
!---------------------------------------------------------------------
      used = send_data (id_qdt_ls, Tend_mp%qtnd(:,:,:), Time, is, js, 1)
      used = send_cmip_data_3d (ID_tnhusscp, Tend_mp%qtnd, Time, is, js, 1)!, rmask=mask)
!---------------------------------------------------------------------
!    surface precip due to large-scale clouds:
!---------------------------------------------------------------------
      used = send_data (id_lsc_precip,   &
                         Precip_state%surfrain + Precip_state%surfsnow, &
                                                            Time, is, js)
        
!---------------------------------------------------------------------
!    precip frequency due to large-scale clouds:
!---------------------------------------------------------------------
      if (id_lsc_freq > 0) then
        where (Precip_state%surfrain > 0. .or. Precip_state%surfsnow > 0.0)
          temp_2d = 1.
        elsewhere
          temp_2d = 0.
        end where
        used = send_data (id_lsc_freq, temp_2d, Time, is, js)
      endif

!---------------------------------------------------------------------
!    total precipitation rate due to large-scale clouds:
!---------------------------------------------------------------------
      used = send_data (id_prec_ls,    &
               Precip_state%surfrain + Precip_state%surfsnow, Time, is, js)

!---------------------------------------------------------------------
!    snowfall rate due to large-scale clouds:
!---------------------------------------------------------------------
      used = send_data (id_snow_ls, Precip_state%surfsnow, Time, is, js)

!---------------------------------------------------------------------
!    define diagnostics specific to a prognostic cloud formulation:
!---------------------------------------------------------------------
      if (doing_prog_clouds) then

!---------------------------------------------------------------------
!    cloud liquid, ice, area and ice crystal number tendencies due to 
!    prognostic cloud parameterization:
!---------------------------------------------------------------------
        used = send_data (id_qldt_ls,  &
                    Tend_mp%q_tnd(:,:,:,nql), Time, is, js, 1)
        if (do_liq_num)   &
          used = send_data (id_qndt_ls,   &
                   Tend_mp%q_tnd(:,:,:,nqn ), Time, is, js, 1)
        used = send_data (id_qidt_ls,   &
                   Tend_mp%q_tnd(:,:,:,nqi), Time, is, js, 1)
        used = send_data (id_qadt_ls,   &
                   Tend_mp%q_tnd(:,:,:,nqa), Time, is, js, 1)
        if (do_ice_num)    &
          used = send_data (id_qnidt_ls,   &
                   Tend_mp%q_tnd(:,:,:,nqni), Time, is, js, 1)

!---------------------------------------------------------------------
!    cloud liquid and ice water path  and particle number tendencies due 
!    to prognostic cloud parameterization:
!---------------------------------------------------------------------
        if (id_ql_ls_col > 0) &
          call column_diag(id_ql_ls_col, is, js, Time,   &
               Tend_mp%q_tnd(:,:,:,nql), 1.0, Input_mp%pmass) 
        if (id_qi_ls_col > 0) &
          call column_diag(id_qi_ls_col, is, js, Time,   &
               Tend_mp%q_tnd(:,:,:,nqi), 1.0, Input_mp%pmass) 
        if (do_liq_num .and. id_qn_ls_col > 0) &
          call column_diag(id_qn_ls_col, is, js, Time,   &
                Tend_mp%q_tnd(:,:,:,nqn), 1.0, Input_mp%pmass)
        if (do_ice_num .and. id_qni_ls_col > 0) &
          call column_diag(id_qni_ls_col, is, js, Time,   &
               Tend_mp%q_tnd(:,:,:,nqni), 1.0, Input_mp%pmass)
      
!---------------------------------------------------------------------
!    column integrated enthalpy and total water tendencies due to 
!    prognostic cloud parameterization:
!---------------------------------------------------------------------
        if (id_enth_ls_col > 0) then
          temp_2d = -HLV*Precip_state%surfrain -HLS*Precip_state%surfsnow
          call column_diag(id_enth_ls_col, is, js, Time,  &
                        Tend_mp%ttnd(:,:,:), CP_AIR,   &
                        Tend_mp%q_tnd(:,:,:,nql), -HLV,  &
                        Tend_mp%q_tnd(:,:,:,nqi), -HLS,   &
                                                  Input_mp%pmass, temp_2d) 
        endif
 
        if (id_wat_ls_col > 0) then
          temp_2d = Precip_state%surfrain+Precip_state%surfsnow
          call column_diag(id_wat_ls_col, is, js, Time,   &
                     Tend_mp%qtnd(:,:,:), 1.0, &
                     Tend_mp%q_tnd(:,:,:,nql), 1.0,  &
                     Tend_mp%q_tnd(:,:,:,nqi), 1.0,   &
                                                Input_mp%pmass, temp_2d) 
        endif

!---------------------------------------------------------------------
!    stratiform cloud volume tendency due to prognostic cloud 
!    parameterization:
!---------------------------------------------------------------------
        if (id_qa_ls_col > 0) &
          call column_diag(id_qa_ls_col, is, js, Time,   &
               Tend_mp%q_tnd(:,:,:,nqa), 1.0, Input_mp%pmass)

!---------------------------------------------------------------------
!    large scale surface rainfall from prognostic clouds
!---------------------------------------------------------------------
        used = send_data(id_lscale_rain3d, Removal_mp%rain3d,   &
                                                        Time, is, js, 1)

!---------------------------------------------------------------------
!    large scale surface snowfall from prognostic clouds
!---------------------------------------------------------------------
        used = send_data(id_lscale_snow3d, Removal_mp%snow3d,   &
                                                         Time, is, js, 1)

      endif

!-----------------------------------------------------------------------


end subroutine basic_diagnostics 



!########################################################################

subroutine lscloud_dealloc (Atmos_state, Particles, Cloud_State, &
                                            Precip_state, Cloud_Processes)
      
!-----------------------------------------------------------------------
type(atmos_state_type),     intent(inout) :: Atmos_state
type(particles_type),       intent(inout) :: Particles   
type(cloud_state_type),     intent(inout) :: Cloud_state
type(precip_state_type),    intent(inout) :: Precip_state
type(cloud_processes_type), intent(inout) :: Cloud_processes


!-----------------------------------------------------------------------
!    deallocate the components of the atmos_state_type variable 
!    Atmos_State.
!-----------------------------------------------------------------------
      deallocate (Atmos_state%airdens     )
      deallocate (Atmos_state%tn          )
      deallocate (Atmos_state%qvn         )
      deallocate (Atmos_state%qs          )
      deallocate (Atmos_state%dqsdT       )
      deallocate (Atmos_state%qsi         )
      deallocate (Atmos_state%qsl         )
      deallocate (Atmos_state%rh_crit     )
      deallocate (Atmos_state%rh_crit_min )
      deallocate (Atmos_state%gamma       )
      deallocate (Atmos_state%esat0       )
      deallocate (Atmos_state%U_ca        )
      deallocate (Atmos_state%delp        )
      deallocate (Atmos_state%U01         )
      deallocate (Atmos_state%pthickness  )
 
!-----------------------------------------------------------------------
!    deallocate the components of the particles_type variable Particles.
!-----------------------------------------------------------------------
      deallocate (Particles%concen_dust_sub )
      deallocate (Particles%drop1           )
      deallocate (Particles%drop2           )
      deallocate (Particles%crystal1        )
      deallocate (Particles%Ndrop_act_CLUBB )
      deallocate (Particles%Icedrop_act_CLUBB )
      deallocate (Particles%rbar_dust       )
      deallocate (Particles%ndust           )
      deallocate (Particles%hom             )
      deallocate (Particles%N3D             )
      deallocate (Particles%N3Di            )
      deallocate (Particles%totalmass1      )
      deallocate (Particles%imass1      )

!-----------------------------------------------------------------------
!    deallocate the components of the cloud_state_type variable 
!    Cloud_state.
!-----------------------------------------------------------------------
      deallocate (Cloud_state%ql_upd  )
      deallocate (Cloud_state%qi_upd  )
      deallocate (Cloud_state%qa_upd  )
      deallocate (Cloud_state%qn_upd  )
      deallocate (Cloud_state%qni_upd )
      deallocate (Cloud_state%ql_mean )
      deallocate (Cloud_state%qi_mean )
      deallocate (Cloud_state%qa_mean )
      deallocate (Cloud_state%qn_mean )
      deallocate (Cloud_state%qni_mean)
      deallocate (Cloud_state%ql_in   )
      deallocate (Cloud_state%qi_in   )
      deallocate (Cloud_state%qa_in   )
      deallocate (Cloud_state%qn_in   )
      deallocate (Cloud_state%qni_in  )
      deallocate (Cloud_state%SL_out  )
      deallocate (Cloud_state%SI_out  )
      deallocate (Cloud_state%SA_out  )
      deallocate (Cloud_state%SN_out  )
      deallocate (Cloud_state%SNi_out )
      deallocate (Cloud_state%qcvar_clubb )
      deallocate (Cloud_state%relvarn     )
      deallocate (Cloud_state%qa_upd_0)
      deallocate (Cloud_state%SA_0    )

!-----------------------------------------------------------------------
!    deallocate the components of the precip_state_type variable 
!    Precip_state.
!-----------------------------------------------------------------------
      deallocate (Precip_state%lsc_snow      )
      deallocate (Precip_state%lsc_rain      )
      deallocate (Precip_state%lsc_snow_size )
      deallocate (Precip_state%lsc_rain_size )
      deallocate (Precip_state%qsout3d_mg    )
      deallocate (Precip_state%qrout3d_mg    )
      deallocate (Precip_state%surfrain      )
      deallocate (Precip_state%surfsnow      )
      deallocate (Precip_state%precip        )

!-----------------------------------------------------------------------
!    deallocate the components of the cloud_processes_type variable 
!    Cloud_processes.
!-----------------------------------------------------------------------
      deallocate (Cloud_processes%da_ls       )
      deallocate (Cloud_processes%D_eros      )
      deallocate (Cloud_processes%qvg         )
      deallocate (Cloud_processes%dcond_ls    )
      deallocate (Cloud_processes%dcond_ls_liquid)
      deallocate (Cloud_processes%dcond_ls_ice)
      deallocate (Cloud_processes%dcond_ls_tot)
      deallocate (Cloud_processes%delta_cf    )
      deallocate (Cloud_processes%f_snow_berg )

!-------------------------------------------------------------------------

end subroutine lscloud_dealloc


!#######################################################################

subroutine adjust_condensate (   &
                mask, delta_cond, delta_q, delta_T, cond_in, lh, cond_out)

!------------------------------------------------------------------------
!    subroutine adjust_condensate eliminates unacceptably small values of
!    condensate in a way which assures mass and enthalpy conservation.
!    tendency terms and output condensate fields are adjusted. 
!------------------------------------------------------------------------

logical, dimension(:,:,:),  intent(in)    :: mask
real, dimension(:,:,:),     intent(in)    :: cond_in
real, dimension(:,:,:),     intent(inout) :: delta_cond, delta_q, delta_T
real, dimension(:,:,:),     intent(out)   :: cond_out
real,                       intent(in)    :: lh

!---------------------------------------------------------------------
!   local variables:

      integer :: idim, jdim, kdim    
      integer :: i, j, k

      call mpp_clock_begin (adjust_cond_clock)

      idim = size(mask,1)
      jdim = size(mask,2)
      kdim = size(mask,3)

!---------------------------------------------------------------------
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            if (mask(i,j,k)) then
              delta_cond(i,j,k) = delta_cond(i,j,k) -   cond_in(i,j,k)
              delta_q(i,j,k) = delta_q(i,j,k) + cond_in(i,j,k)
              delta_T(i,j,k) = delta_T(i,j,k) - lh*cond_in(i,j,k)/CP_AIR
              cond_out(i,j,k) = 0.
            else
              cond_out(i,j,k) = cond_in(i,j,k)
            endif
          end do
        end do
      end do

      call mpp_clock_end (adjust_cond_clock)

!------------------------------------------------------------------------


end subroutine adjust_condensate 


!#########################################################################

!----------------------------------------------------------------------

subroutine adjust_particle_number (mask, delta_particles,   &
                                             particles_in, particles_out)

logical, dimension(:,:,:),  intent(in)    :: mask
real, dimension(:,:,:),     intent(in)    :: particles_in
real, dimension(:,:,:),     intent(inout) :: delta_particles
real, dimension(:,:,:),     intent(out)   :: particles_out

!---------------------------------------------------------------------
!   local variables:

      integer :: idim, jdim, kdim    
      integer :: i, j, k

      call mpp_clock_begin (adjust_part_num_clock)

      idim = size(mask,1)
      jdim = size(mask,2)
      kdim = size(mask,3)

!---------------------------------------------------------------------
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              if (mask(i,j,k)) then
                delta_particles(i,j,k) = delta_particles(i,j,k) -  &
                                                 particles_in(i,j,k)
                particles_out(i,j,k) = 0.
              else
                particles_out(i,j,k) = particles_in(i,j,k)
              endif
            end do
          end do
        end do

      call mpp_clock_end (adjust_part_num_clock)

!------------------------------------------------------------------------


end subroutine adjust_particle_number


!########################################################################




                   end module lscloud_driver_mod


