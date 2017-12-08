
                    module convection_driver_mod

!-----------------------------------------------------------------------
!
!         interface module for convective processes
!         ---------------------------------------
!             dry convective adjustment
!             moist convective adjustment
!             relaxed arakawa-schubert
!             donner deep convection
!             betts-miller convective adjustment
!             uw shallow convection
!             cumulus momentum transport
!
!-----------------------------------------------------------------------

! fms modules
use sat_vapor_pres_mod,    only: compute_qs
use time_manager_mod,      only: time_type, get_time
use diag_manager_mod,      only: register_diag_field, send_data, &
                                 get_diag_field_id, DIAG_FIELD_NOT_FOUND
use diag_data_mod,         only: CMOR_MISSING_VALUE
use mpp_mod,               only: input_nml_file
use fms_mod,               only: error_mesg, FATAL, WARNING,NOTE,&
                                 file_exist, check_nml_error,    &
                                 open_namelist_file, close_file, &
                                 write_version_number,           &
                                 mpp_pe, mpp_root_pe, stdlog,    &
                                 mpp_clock_id, mpp_clock_begin,  &
                                 mpp_clock_end, CLOCK_MODULE,    &
                                 CLOCK_MODULE_DRIVER, &
                                 MPP_CLOCK_SYNC, read_data, write_data
use field_manager_mod,     only: MODEL_ATMOS
use tracer_manager_mod,    only: get_tracer_index,&
                                 get_tracer_names, &
                                 query_method, &
                                 NO_TRACER
use constants_mod,         only: CP_AIR, HLV, HLS, HLF, RDGAS, RVGAS, &
                                 SECONDS_PER_DAY, KAPPA
! atmos_param modules
use physics_types_mod,     only: physics_control_type, phys_mp_exch_type
use vert_diff_driver_mod,  only: surf_diff_type
use physics_radiation_exch_mod,        &
                           only: clouds_from_moist_block_type, &
                                 exchange_control_type
use betts_miller_mod,      only: betts_miller, betts_miller_init
use bm_massflux_mod,       only: bm_massflux, bm_massflux_init
use bm_omp_mod,            only: bm_omp, bm_omp_init
use donner_deep_mod,       only: donner_deep_init,               &
                                 donner_deep_time_vary,  &
                                 donner_deep_endts,         &
                                 donner_deep, donner_deep_end,   &
                                 donner_deep_restart
use moist_conv_mod,        only: moist_conv, moist_conv_init
use uw_conv_mod,           only: uw_conv_end, uw_conv_init
use ras_mod,               only: ras_end, ras_init
use dry_adj_mod,           only: dry_adj, dry_adj_init
use detr_ice_num_mod,      only: detr_ice_num, detr_ice_num_init,   &
                                 detr_ice_num_end
use rh_clouds_mod,         only: rh_clouds_init, rh_clouds_end, &
                                 rh_clouds_sum
use cu_mo_trans_mod,       only: cu_mo_trans_init, cu_mo_trans,   &
                                 cu_mo_trans_end
use moz_hook_mod,          only: moz_hook
use aerosol_types_mod,     only: aerosol_type
use moist_proc_utils_mod,  only: capecalcnew, column_diag, rh_calc, &
                                 mp_nml_type, mp_input_type,   &
                                 mp_tendency_type, mp_removal_type, &
                                 mp_removal_control_type, &
                                 mp_conv2ls_type, mp_output_type
use moistproc_kernels_mod, only: moistproc_mca, moistproc_ras,        &
                                 moistproc_cmt, moistproc_uw_conv,    &
                                 moistproc_scale_uw, moistproc_scale_donner
! atmos_shared modules
use atmos_tracer_utilities_mod,         &
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
implicit none
private

!-----------------------------------------------------------------------
!-------------------- public data/interfaces ---------------------------

   public  convection_driver_init, convection_driver,   &
           convection_driver_time_vary, convection_driver_endts,  &
           convection_driver_end, convection_driver_restart,   &
           cape_cin_diagnostics
  
   private diag_field_init, prevent_neg_precip_fluxes,  &
           compute_convective_area

!-----------------------------------------------------------------------
!-------------------- private data -------------------------------------

!--------------------- version number ----------------------------------
character(len=128) ::  version = '$Id: $'
character(len=128) :: tagname = '$Name: $'
!-------------------- namelist data (private) --------------------------

!---------------- namelist variable definitions -----------------------


!
!   do_limit_donner = limit Donner deeo tendencies to prevent the
!                formation of grid points with negative water vapor,
!                liquid or ice.
!
!   do_limit_uw = limit UW shallow tendencies to prevent the formation
!                of grid points with negative total water specific 
!                humidities. This situation can occur because both
!                shallow and deep convection operate on the same
!                soundings without knowledge of what the other is doing
!
!   do_unified_convective_closure = use cloud base mass flux calculated
!                in uw_conv module as value for donner deep parameter-
!                ization; adjust cbmf available for uw shallow appropr-
!                iately. only available when uw shallow and donner deep
!                are the active convective schemes

!   cmt_mass_flux_source = parameterization(s) being used to supply the 
!                mass flux profiles seen by the cumulus momentum transport
!                module; currently either 'ras', 'donner', 'uw', 
!                'donner_and_ras', 'donner_and_uw', 'ras_and_uw', 
!                'donner_and_ras_and_uw' or 'all'
!   do_gust_cv = switch to use convective gustiness (default = false)
!   do_gust_cv_new = new switch to use convective gustiness    
!                                                    (default = false)
!   gustmax    = maximum convective gustiness (m/s)
!   gustconst  = precip rate which defines precip rate which begins to
!                matter for convective gustiness (kg/m2/sec)
!                default = 1 cm/day = 10 mm/da
integer, public    :: id_pr_g, id_prc_g, id_prsn_g, id_prsnc, id_prrc
logical            :: do_limit_donner =.true. 
logical            :: do_limit_uw = .true.    
logical            :: do_unified_convective_closure = .false.
character(len=64)  :: cmt_mass_flux_source = 'all'
logical            :: do_gust_cv = .false.
logical            :: do_gust_cv_new = .false.
real               :: gustmax = 3.         ! maximum gustiness wind (m/s)
real               :: gustconst = 10./SECONDS_PER_DAY 
logical            :: do_donner_before_uw = .true.
logical            :: use_updated_profiles_for_uw = .false.
logical            :: only_one_conv_scheme_per_column = .false.
logical            :: do_cmt=.true.
logical            :: force_donner_moist_conserv = .false.
logical            :: do_donner_conservation_checks = .true.
logical            :: do_donner_mca=.false.
logical            :: using_fms = .true.
logical            :: detrain_liq_num = .false.
logical            :: detrain_ice_num = .false.
real               :: conv_frac_max = 0.99
logical            :: remain_detrain_bug = .false.
logical            :: keep_icenum_detrain_bug = .false.

namelist /convection_driver_nml/    &
              do_limit_donner, do_limit_uw, &
              do_unified_convective_closure, &
              cmt_mass_flux_source, do_gust_cv,  do_gust_cv_new, &
              gustmax, gustconst, &
              do_donner_before_uw,  use_updated_profiles_for_uw,  &
              only_one_conv_scheme_per_column,  do_cmt, &
              force_donner_moist_conserv, do_donner_conservation_checks, &
              do_donner_mca, using_fms, detrain_liq_num, detrain_ice_num, &
              conv_frac_max, remain_detrain_bug, keep_icenum_detrain_bug


!-------------------- clock definitions --------------------------------

integer :: convection_clock, donner_clock, mca_clock, ras_clock, &
           donner_mca_clock, bm_clock, cmt_clock, shallowcu_clock

!-------------------- diagnostics variables ----------------------------

integer :: id_cell_cld_frac,  id_meso_cld_frac, id_donner_humidity_area
integer :: id_tdt_conv, id_qdt_conv, id_prec_conv, id_snow_conv, &
           id_conv_freq, id_gust_conv, id_mc_full, id_mc_half,   &
           id_tdt_dadj, id_mc_donner, id_mc_donner_half, &
           id_mc_conv_up, id_conv_cld_base, id_conv_cld_top, &
           id_tdt_deep_donner, id_qdt_deep_donner, &
           id_qadt_deep_donner, id_qldt_deep_donner, &
           id_qidt_deep_donner, &
           id_qndt_deep_donner,  id_qnidt_deep_donner, &
           id_tdt_mca_donner, id_qdt_mca_donner, &
           id_prec_deep_donner, id_precret_deep_donner,id_prec_mca_donner,&
           id_tdt_uw, id_qdt_uw, &
           id_qadt_uw, id_qldt_uw, id_qidt_uw, id_qndt_uw, id_qnidt_uw, &
           id_prec1_deep_donner, id_snow_deep_donner, id_snow_mca_donner, &
           id_qadt_conv, id_qldt_conv, id_qndt_conv, id_qidt_conv, &
           id_qnidt_conv, &
           id_qa_conv_col, id_ql_conv_col, id_qn_conv_col,  &
           id_qni_conv_col, id_qi_conv_col, &
           id_bmflag, id_klzbs, id_invtaubmt, id_invtaubmq, &
           id_massflux, id_tref, id_qref, &
           id_tp, id_rp, id_lcl, id_lfc, id_lzb, &
           id_q_conv_col, id_t_conv_col,              &
           id_enth_conv_col, id_wat_conv_col, &
           id_enth_donner_col, id_wat_donner_col, &
           id_enth_donner_col2,  &
           id_enth_donner_col3,  &
           id_enth_donner_col4,  &
           id_enth_donner_col5,  &
           id_enth_donner_col6,  &
           id_enth_donner_col7,  &
           id_enth_mca_donner_col, id_wat_mca_donner_col, &
           id_enth_uw_col, id_wat_uw_col, &
           id_scale_donner, id_scale_uw, &
           id_ras_precip, id_ras_freq, id_don_precip, id_don_freq, &
           id_uw_precip, id_uw_snow, id_uw_freq, id_prod_no, &
           id_m_cdet_donner, id_m_cellup, id_conv_rain3d, id_conv_snow3d
integer :: id_cape, id_cin
integer :: id_vaporint, id_condensint, id_precipint, id_diffint
integer :: id_vertmotion
integer :: id_max_enthalpy_imbal_don, id_max_water_imbal_don
integer :: id_enthint, id_lprcp, id_lcondensint, id_enthdiffint
 
integer :: id_prc, id_ci, id_ccb, id_cct

integer, dimension(:), allocatable ::    &
                                      id_tracerdt_conv,  &
                                      id_tracerdt_conv_col, &
                                      id_conv_tracer,  &
                                      id_conv_tracer_col, &
                                      id_tracerdt_mcadon, &
                                      id_tracerdt_mcadon_col

real                      :: missing_value = -999.
character(len=5), private :: mod_name = 'moist'
character(len=8), private :: mod_name_tr = 'moist_tr'  ! fixes name conflict with
                                                       ! tracers and moist



!------------------- other global variables  ------------------------

real, allocatable, dimension(:,:)   ::  &
                          max_enthalpy_imbal_don, max_water_imbal_don

logical :: do_liq_num
logical :: do_ice_num
logical :: do_cosp
integer :: do_clubb
integer :: do_lsc
logical :: do_mca 
logical :: do_ras 
logical :: do_uw_conv
logical :: do_donner_deep
logical :: do_dryadj
logical :: limit_conv_cloud_frac
logical :: include_donmca_in_cosp
logical :: do_rh_clouds
logical :: do_bm
logical :: do_bmmass
logical :: do_bmomp
logical :: do_simple

logical :: doing_prog_clouds
logical :: doing_diffusive
integer :: nsphum, nql, nqi, nqa, nqn, nqni, nqr,nqs, nqg
integer :: num_prog_tracers
logical, dimension(:), allocatable :: cloud_tracer
logical :: cmt_uses_donner, cmt_uses_ras, cmt_uses_uw
logical :: module_is_initialized = .false.

type(cmip_diag_id_type) :: ID_tntc, ID_tnhusc, ID_mc


                             contains


!#######################################################################

subroutine convection_driver_init (id, jd, kd, axes, Time,   &
                    Physics_control, Exch_ctrl, Nml_mp, Control,   &
                     lonb, latb, pref )
 
integer,                       intent(in)    :: id, jd, kd
integer,                       intent(in)    :: axes(4)
type(time_type),               intent(in)    :: Time
type(physics_control_type),    intent(in)    :: Physics_control
type(exchange_control_type),   intent(in)    :: Exch_ctrl
type(mp_nml_type),             intent(in)    :: Nml_mp
type(mp_removal_control_type), intent(inout) :: Control   
real, dimension(:,:),          intent(in)    :: lonb, latb
real, dimension(:),            intent(in)    :: pref
 
      integer :: secs, days  
      integer :: n, unit, io, ierr, logunit
      character(len=80)  :: scheme

!!-----------------------------------------------------------------------

      if ( module_is_initialized ) return

!-----------------------------------------------------------------------
!   process the moist_processes_nml.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=convection_driver_nml, iostat=io)
        ierr = check_nml_error(io,'convection_driver_nml')
#else
        unit = open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
          read  (unit, nml=convection_driver_nml, iostat=io, end=10)
          ierr = check_nml_error(io,'convection_driver_nml')
        enddo
  10    call close_file (unit)
#endif

!--------- write version and namelist to standard log ------------

        call write_version_number ( version, tagname )
        logunit = stdlog()
        if ( mpp_pe() == mpp_root_pe() ) &
                        write ( logunit, nml=convection_driver_nml )
      endif

!----------------------------------------------------------------------
!    define needed local module variables from the derived types passed in.
!----------------------------------------------------------------------
      do_liq_num = Exch_ctrl%do_liq_num
      do_ice_num = Exch_ctrl%do_ice_num
      do_clubb = Exch_ctrl%do_clubb
      do_lsc = Nml_mp%do_lsc
      do_cosp = Exch_ctrl%do_cosp
      do_mca = Nml_mp%do_mca
      do_ras = Nml_mp%do_ras
      do_uw_conv = Nml_mp%do_uw_conv
      do_donner_deep = Nml_mp%do_donner_deep
      limit_conv_cloud_frac = Nml_mp%limit_conv_cloud_frac
      do_dryadj = Nml_mp%do_dryadj
      include_donmca_in_cosp = Nml_mp%include_donmca_in_cosp
      do_rh_clouds = Nml_mp%do_rh_clouds
      do_bm = Nml_mp%do_bm
      do_bmmass = Nml_mp%do_bmmass
      do_bmomp  = Nml_mp%do_bmomp 
      do_simple = Nml_mp%do_simple
      doing_prog_clouds = Exch_ctrl%doing_prog_clouds
      nsphum = Physics_control%nsphum
      nql = Physics_control%nql
      nqi = Physics_control%nqi
      nqa = Physics_control%nqa
      nqn = Physics_control%nqn
      nqni = Physics_control%nqni
      nqr = Physics_control%nqr
      nqs = Physics_control%nqs
      nqg = Physics_control%nqg
      num_prog_tracers = Physics_control%num_prog_tracers

      allocate (cloud_tracer(size(Physics_control%cloud_tracer)))
      cloud_tracer = Physics_control%cloud_tracer

!-----------------------------------------------------------------------
!   check to make sure only compatible convection schemes have
!   been activated.
!-----------------------------------------------------------------------
      if (do_mca .and. do_ras ) call error_mesg   &
         ('convection_driver_init',  &
               'both do_mca and do_ras cannot be specified', FATAL)

      if (do_mca .and. do_bm ) call error_mesg   &
         ('convection_driver_init',  &
                    'both do_mca and do_bm cannot be specified', FATAL)
      if (do_ras .and. do_bm ) call error_mesg   &
         ('convection_driver_init',  &
                     'both do_bm and do_ras cannot be specified', FATAL)
      if (do_bm .and. do_bmmass ) call error_mesg   &
         ('convection_driver_init',  &
                   'both do_bm and do_bmmass cannot be specified', FATAL)
      if (do_bm .and. do_bmomp ) call error_mesg   &
         ('convection_driver_init',  &
                   'both do_bm and do_bmomp cannot be specified', FATAL)
      if (do_bmomp .and. do_bmmass ) call error_mesg   &
         ('convection_driver_init',  &
             'both do_bmomp and do_bmmass cannot be specified', FATAL)
      if (do_bmmass .and. do_mca ) call error_mesg   &
         ('convection_driver_init',  &
             'both do_bmmass and do_mca cannot be specified', FATAL)
      if (do_bmmass .and. do_ras ) call error_mesg   &
         ('convection_driver_init',  &
                  'both do_bmmass and do_ras cannot be specified', FATAL)
      if (do_bmomp .and. do_mca ) call error_mesg   &
         ('convection_driver_init',  &
                 'both do_bmomp and do_mca cannot be specified', FATAL)
      if (do_bmomp .and. do_ras ) call error_mesg   &
         ('convection_driver_init',  &
                 'both do_bmomp and do_ras cannot be specified', FATAL)

      if (do_mca .and. do_donner_deep) call error_mesg &
         ('convection_driver_init',  &
             'do_donner_deep and do_mca cannot both be active', FATAL)

!----------------------------------------------------------------------
!   check for inconsistent / invalid settings involving nml variables.
!----------------------------------------------------------------------
      if (do_donner_deep) then 
        if (include_donmca_in_cosp .and. (.not. do_donner_mca) ) &
          call error_mesg ('convection_driver_init', &
             'trying to include donmca in COSP when donmca is inactive', &
                                                                  FATAL)
        if (do_cosp .and. .not. (do_donner_conservation_checks)) then
          do_donner_conservation_checks = .true.
          call error_mesg ('convection_driver_init', &
              'setting do_donner_conservation_checks to true so that &
                 &needed fields for COSP are produced.', NOTE)
        endif
      endif

      if (force_donner_moist_conserv .and.    &
                                 .not. do_donner_conservation_checks) then
        call error_mesg ('convection_driver_init', &
              'when force_donner_moist_conserv is .true., &
                &do_donner_conservation_checks must be .true.', FATAL)
      endif
 
      if (use_updated_profiles_for_uw .and.   &     
                            .not. (do_donner_before_uw) ) then
        call error_mesg ('convection_driver_init', &
         'use_updated_profiles_for_uw is only meaningful when &
                              &do_donner_before_uw is true', FATAL)
      endif

      if (only_one_conv_scheme_per_column .and.   &
                .not. (do_donner_before_uw) ) then
        call error_mesg ('convection_driver_init', &
          'only_one_conv_scheme_per_column is only meaningful when &
                             &do_donner_before_uw is true', FATAL)
      endif
 
      if (limit_conv_cloud_frac .and.  (.not. do_donner_before_uw)) then
        call error_mesg ('convection_driver_init', &
            'when limit_conv_cloud_frac is .true., &
                             &do_donner_before_uw must be .true.', FATAL)
      endif

      if (do_unified_convective_closure) then
        call error_mesg ('convection_driver_init', &
         'do_unified_convective_closure is currently not allowed', FATAL)
      endif

      if (do_cmt) then
        if ( .not. do_ras .and. .not. do_donner_deep  .and. &
                                          .not. do_uw_conv) then
          call error_mesg ( 'convection_driver_init', &
                'do_cmt is active but no cumulus schemes activated', &
                                                              FATAL)
        endif
      endif

!----------------------------------------------------------------------
!    for each tracer, determine if it is to be transported by convect-
!    ion, and the convection schemes that are to transport it. set a 
!    logical flag to .true. for each tracer that is to be transported by
!    each scheme and increment the count of tracers to be transported
!    by that scheme.
!----------------------------------------------------------------------
      do n=1, num_prog_tracers
        if (query_method ('convection', MODEL_ATMOS, n, scheme)) then
          select case (scheme)
            case ("none")
            case ("donner")
               Control%num_donner_tracers = Control%num_donner_tracers + 1
               Control%tracers_in_donner(n) = .true.
            case ("mca")
               Control%num_mca_tracers = Control%num_mca_tracers + 1
               Control%tracers_in_mca(n) = .true.
            case ("ras")
               Control%num_ras_tracers = Control%num_ras_tracers + 1
               Control%tracers_in_ras(n) = .true.
            case ("uw")
               Control%num_uw_tracers = Control%num_uw_tracers + 1
               Control%tracers_in_uw(n) = .true.
            case ("donner_and_ras")
               Control%num_donner_tracers = Control%num_donner_tracers + 1
               Control%tracers_in_donner(n) = .true.
               Control%num_ras_tracers = Control%num_ras_tracers + 1
               Control%tracers_in_ras(n) = .true.
            case ("donner_and_mca")
               Control%num_donner_tracers = Control%num_donner_tracers + 1
               Control%tracers_in_donner(n) = .true.
               Control%num_mca_tracers = Control%num_mca_tracers + 1

               Control%tracers_in_mca(n) = .true.
            case ("mca_and_ras")
               Control%num_mca_tracers = Control%num_mca_tracers + 1
               Control%tracers_in_mca(n) = .true.
               Control%num_ras_tracers = Control%num_ras_tracers + 1
               Control%tracers_in_ras(n) = .true.
            case ("all")
               Control%num_donner_tracers = Control%num_donner_tracers + 1
               Control%tracers_in_donner(n) = .true.
               Control%num_mca_tracers = Control%num_mca_tracers + 1
               Control%tracers_in_mca(n) = .true.
               Control%num_ras_tracers = Control%num_ras_tracers + 1
               Control%tracers_in_ras(n) = .true.
               Control%num_uw_tracers = Control%num_uw_tracers + 1
               Control%tracers_in_uw(n) = .true.
            case ("all_nodonner")
               Control%num_mca_tracers = Control%num_mca_tracers + 1
               Control%tracers_in_mca(n) = .true.
               Control%num_ras_tracers = Control%num_ras_tracers + 1
               Control%tracers_in_ras(n) = .true.
               Control%num_uw_tracers = Control%num_uw_tracers + 1
               Control%tracers_in_uw(n) = .true.
            case default  ! corresponds to "none"
          end select
        endif
      end do

!-----------------------------------------------------------------------
!   initialize the convective schemes active in this experiment.
!-----------------------------------------------------------------------
      if (do_dryadj) call dry_adj_init ()
      if (do_cmt)    call cu_mo_trans_init (axes, Time, doing_diffusive)
      if (do_bm)     call betts_miller_init () 
      if (do_bmmass) call bm_massflux_init()
      if (do_bmomp)  call bm_omp_init () 

!-------------------------------------------------------------------------
!   define logicals indicating which convective schemes are to be seen
!   by the cumulus momentum transport scheme.
!-------------------------------------------------------------------------
      if (do_cmt) then
        if (trim(cmt_mass_flux_source) == 'ras') then
          cmt_uses_ras = .true.
          cmt_uses_donner = .false.
          cmt_uses_uw = .false.
          if (.not. do_ras) then
            call error_mesg ('moist_processes_mod', &
              'if cmt_uses_ras = T, then do_ras must be T', FATAL)
          endif

        else if (trim(cmt_mass_flux_source) == 'donner') then
          cmt_uses_ras = .false.
          cmt_uses_donner = .true.
          cmt_uses_uw = .false.
          if (.not. do_donner_deep)  then
            call error_mesg ('convection_driver_init', &
              'if cmt_uses_donner = T, then do_donner_deep must be T', &
                                                                    FATAL)
          endif

        else if (trim(cmt_mass_flux_source) == 'uw') then
          cmt_uses_ras = .false.
          cmt_uses_donner = .false.
          cmt_uses_uw = .true.
          if (.not. do_uw_conv)  then
            call error_mesg ('convection_driver_init', &
                'if cmt_uses_uw = T, then do_uw_conv must be T', FATAL)
          endif

        else if (trim(cmt_mass_flux_source) == 'donner_and_ras') then
          cmt_uses_ras = .true.
          if (.not. do_ras) then
            call error_mesg ('convection_driver_init', &
              'if cmt_uses_ras = T, then do_ras must be T', FATAL)
          endif
          cmt_uses_donner = .true.
          if (.not.  do_donner_deep)  then
            call error_mesg ('convection_driver_init', &
              'if cmt_uses_donner = T, then do_donner_deep must be T', &
                                                                    FATAL)
          endif
          cmt_uses_uw = .false.

        else if (trim(cmt_mass_flux_source) == 'donner_and_uw') then
          cmt_uses_uw = .true.
          if (.not. do_uw_conv)  then
            call error_mesg ('convection_driver_init', &
                'if cmt_uses_uw = T, then do_uw_conv must be T', FATAL)
          endif
          cmt_uses_donner = .true.
          if (.not.  do_donner_deep)  then
            call error_mesg ('convection_driver_init', &
              'if cmt_uses_donner = T, then do_donner_deep must be T', &
                                                                    FATAL)
          endif
          cmt_uses_ras = .false.

        else if (trim(cmt_mass_flux_source) == 'ras_and_uw') then
          cmt_uses_ras = .true.
          if (.not. do_ras) then
            call error_mesg ('convection_driver_init', &
              'if cmt_uses_ras = T, then do_ras must be T', FATAL)
          endif
          cmt_uses_uw = .true.
          if (.not. do_uw_conv)  then
            call error_mesg ('convection_driver_init', &
                'if cmt_uses_uw = T, then do_uw_conv must be T', FATAL)
          endif
          cmt_uses_donner = .false.

        else if   &
              (trim(cmt_mass_flux_source) == 'donner_and_ras_and_uw') then
          cmt_uses_ras = .true.
          if (.not. do_ras) then
            call error_mesg ('convection_driver_init', &
              'if cmt_uses_ras = T, then do_ras must be T', FATAL)
          endif
          cmt_uses_donner = .true.
          if (.not.  do_donner_deep)  then
            call error_mesg ('convection_driver_init', &
              'if cmt_uses_donner = T, then do_donner_deep must be T', &
                                                                    FATAL)
          endif
          cmt_uses_uw = .true.
          if (.not. do_uw_conv)  then
            call error_mesg ('convection_driver_init', &
                'if cmt_uses_uw = T, then do_uw_conv must be T', FATAL)
          endif
        else if (trim(cmt_mass_flux_source) == 'all') then
          if (do_ras) then
            cmt_uses_ras = .true.
          else
            cmt_uses_ras = .false.
          endif
          if (do_donner_deep)  then
            cmt_uses_donner = .true.
          else
            cmt_uses_donner = .false.
          endif
          if (do_uw_conv)  then
            cmt_uses_uw = .true.
          else
            cmt_uses_uw = .false.
          endif
        else
          call error_mesg ('convection_driver_init', &
             'invalid specification of cmt_mass_flux_source', FATAL)
        endif

        if (cmt_uses_uw .and. .not. doing_diffusive) then
          call error_mesg ('convection_driver_init', &
             'currently cannot do non-local cmt with uw as mass &
                                                &flux_source', FATAL)
        endif
      endif  !(do_cmt)

!--------------------------------------------------------------------
!    continue the initialization of the convection scheme modules.
!--------------------------------------------------------------------
      if (do_donner_deep) then
        call get_time (Time, secs, days)
        call donner_deep_init (lonb, latb, pref, axes, secs, days,  &
                               Control%tracers_in_donner,  &
                               do_donner_conservation_checks, &
                               do_unified_convective_closure, &
                               using_fms)
        if (do_donner_conservation_checks) then
          allocate (max_enthalpy_imbal_don (id, jd))
          allocate (max_water_imbal_don (id, jd))
          max_enthalpy_imbal_don = 0.
          max_water_imbal_don = 0.
        endif
      endif ! (do_donner_deep)
 
      if (do_ras)  then
        call ras_init (doing_prog_clouds, do_liq_num, axes, Time, Nml_mp, &
                                                    Control%tracers_in_ras)
      endif

      if (do_uw_conv) call uw_conv_init (doing_prog_clouds, axes, Time,   &
                                         kd, Nml_mp, Control%tracers_in_uw)

      if (do_mca .or. do_donner_mca)  then
        call  moist_conv_init (axes, Time, Control%tracers_in_mca)
      endif

!-----------------------------------------------------------------------
!   initialize clocks.
!-----------------------------------------------------------------------
      convection_clock = mpp_clock_id( '   Physics_up: Moist Proc: Conv' ,&
                                             grain=CLOCK_MODULE_DRIVER )
      donner_clock     = mpp_clock_id( '   Moist Processes: Donner_deep' ,&
                                             grain=CLOCK_MODULE_DRIVER )
      mca_clock        = mpp_clock_id( '   Moist Processes: MCA'         ,&
                                             grain=CLOCK_MODULE_DRIVER )
      donner_mca_clock = mpp_clock_id( '   Moist Processes: Donner_MCA'  ,&
                                             grain=CLOCK_MODULE_DRIVER )
      ras_clock        = mpp_clock_id( '   Moist Processes: RAS'         ,&
                                             grain=CLOCK_MODULE_DRIVER )
      shallowcu_clock  = mpp_clock_id( '   Moist Processes: Shallow_cu'  ,&
                                             grain=CLOCK_MODULE_DRIVER )
      cmt_clock        = mpp_clock_id( '   Moist Processes: CMT'         ,&
                                             grain=CLOCK_MODULE_DRIVER )
      bm_clock         = mpp_clock_id( '   Moist Processes: Betts-Miller',&
                                             grain=CLOCK_MODULE_DRIVER )
 
!------------------------------------------------------------------------
!   call diag_field_init to register the netcdf diagnostic fields.
!------------------------------------------------------------------------
      call diag_field_init (axes, Time, Control)

!------------------------------------------------------------------------
!   initialize the ice number detrain module.
!------------------------------------------------------------------------
      call detr_ice_num_init

!---------------------------------------------------------------------
      module_is_initialized = .true.


end subroutine convection_driver_init



!########################################################################

subroutine convection_driver   &
                   (is, ie, js, je, Time, dt, Input_mp, Tend_mp, C2ls_mp, &
                    Output_mp, Removal_mp,  Removal_mp_control,   &
                    Surf_diff, Phys_mp_exch, &
                    shflx, lhflx, tdt_dif, qdt_dif, Moist_clouds_block,  &
                    Aerosol, mask, kbot)

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
!
!       optional
!  -----------------
! 
!
!
!-----------------------------------------------------------------------
integer,                intent(in)    :: is,ie,js,je
type(time_type),        intent(in)    :: Time
real,                   intent(in)    :: dt
type(mp_input_type),    intent(inout) :: Input_mp
type(mp_output_type),   intent(inout) :: Output_mp
type(mp_tendency_type), intent(inout) :: Tend_mp
type(mp_removal_type),  intent(inout) :: Removal_mp
type(mp_removal_control_type),  intent(inout) :: Removal_mp_control
type(surf_diff_type),   intent(in)    :: Surf_diff
type(phys_mp_exch_type),intent(inout) :: Phys_mp_exch
real, dimension(:,:),   intent(in)    :: shflx, lhflx
real, dimension(:,:,:), intent(in)    :: tdt_dif, qdt_dif
type(mp_conv2ls_type),  intent(inout) :: C2ls_mp
type(clouds_from_moist_block_type), &
                        intent(inout) :: Moist_clouds_block

type(aerosol_type),      intent(in), optional :: Aerosol
real,  dimension(:,:,:), intent(in), optional :: mask
integer, dimension(:,:), intent(in), optional :: kbot


!-----------------------------------------------------------------------
!    local variables:


   real, dimension(size(Input_mp%t,1),size(Input_mp%t,2)) ::      &
                                cape, cin,  precip, total_precip,   &
                                lheat_precip, precip_returned,    &
                                precip_adjustment, vert_motion, &
                                rain, snow, rain_don, snow_don, rain_ras, &
                                snow_ras, rain_donmca, snow_donmca, &
                                bmflag, klzbs, invtaubmt, invtaubmq, &
                                scale,  sumneg, freq_count,  enthint, &
                                lcondensint, enthdiffint, vaporint, &
                                condensint, precipint, diffint, &
                                sfc_sh_flux, sfc_vapor_flux,   &
                                adjust_frac, temp_2d, &
                                rain_uw, snow_uw, scale_donner

   logical, dimension(size(Input_mp%t,1),size(Input_mp%t,2)) ::    &
                                conv_calc_completed, ltemp

   integer, dimension(size(Output_mp%rdt,1),size(Output_mp%rdt,2)) :: &
                                cldtop, cldbot, klzb, klcl, klfc, &
                                maxTe_launch_level

   real, dimension(size(Input_mp%t,1),size(Input_mp%t,2),  &
                                         size(Input_mp%phalf,3)) ::    &
                                rain3d, snow3d, snowclr3d, cmf, &
                                mc, m_cellup, mc_cmt, mc_donner_half

   real, dimension(size(Input_mp%t,1),size(Input_mp%t,2),   &
                                         size(Input_mp%t,3)) ::  &
                                f_snow_berg, prod_no, ttnd_adjustment, &
                                available_cf_for_uw, total_cloud_area,  &
                                rin, rh, temp_3d1, temp_3d2, &
                                tin_orig, qin_orig, & 
                                t_ref, q_ref, ttnd_don, qtnd_don, &
                                delta_temp, delta_q, delta_vapor,   &
                                delta_qn, delta_qni, &
                                tp, rp, thetae, &
                                liquid_precip, frozen_precip, &
                                det0, det_cmt, mc_donner, mc_donner_up, &
                                m_cdet_donner, massflux, nllin, nilin, &
                                ttnd_uw, qtnd_uw, utnd_uw, vtnd_uw,   &
                                qltnd_uw, qitnd_uw, qatnd_uw, qntnd_uw, &
                                qnitnd_uw, delta_ql, delta_qi, delta_qa, &
                                qlin, qiin, qain

   real, dimension(size(Input_mp%t,1),size(Input_mp%t,2),  &
                                         size(Input_mp%t,3),   &
                                         size(Output_mp%rdt,4)   ) ::  &
                                rdt_init

   real, dimension(size(Input_mp%t,1),size(Input_mp%t,2),   &
                                         size(Input_mp%t,3),    &
                                         size(Input_mp%r,4)) ::   &
                                tracer_orig

   real, dimension(size(Input_mp%t,1),size(Input_mp%t,2),    &
                                         size(Input_mp%t,3),    &
                               Removal_mp_control%num_donner_tracers) ::  &
                                qtr, donner_tracer

   real, dimension(size(Input_mp%t,1),size(Input_mp%t,2),   &
                               Removal_mp_control%num_donner_tracers) ::  &
                                tr_flux  

   real, dimension(size(Input_mp%t,1),size(Input_mp%t,2),    &
                                         size(Input_mp%t,3),   &
                               Removal_mp_control%num_uw_tracers) ::  &
                                qtruw

   real, dimension(size(Input_mp%t,1),size(Input_mp%t,2),   &
                                         size(Output_mp%rdt,4)   ) ::   &
                                total_wetdep, total_wetdep_uw, &
                                total_wetdep_donner

   character(len=128) :: warn_mesg

   integer         :: secs, days
   integer         :: n, nn, i, j, k, ix, jx, kx, nt, tr
   integer         :: m, mm
   logical         :: used, avgbl
   real            :: dtinv
   real, parameter :: boltz = 1.38044e-16
   real            :: qnew
   real            :: current_total_sec
   integer         :: current_sec, current_days
   integer         :: i_cell, i_meso, i_shallow
   real            :: temp


!------------------------------------------------------------------------
!   local variables:
!
!     sfc_sh_flux      sensible heat flux across the surface
!                      [ watts / m**2 ]
!     sfc_vapor_flux   water vapor flux across the surface
!                      [ kg(h2o) / (m**2 sec) ]
!     tr_flux          tracer fux across the surface
!                      [ kg(tracer) / (m**2 sec) ]
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    verify that the module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('convection_driver_mod',  &
                 'convection_driver_init has not been called.', FATAL)
      endif

!-----------------------------------------------------------------------
!    define array dimensions.
!-----------------------------------------------------------------------
      ix = size(Input_mp%t,1) 
      jx = size(Input_mp%t,2) 
      kx = size(Input_mp%t,3) 
      nt = size(Output_mp%rdt,4)

!-----------------------------------------------------------------------
!    define index into cloud data for the different convective types.
!-----------------------------------------------------------------------
      i_shallow = Moist_clouds_block%index_uw_conv
      i_meso    = Moist_clouds_block%index_donner_meso
      i_cell    = Moist_clouds_block%index_donner_cell

!-----------------------------------------------------------------------
!    save the value of rdt upon entry to the routine do the total
!    convective contribution to the tendency may later be defined.
!    save input tracer fields.
!-----------------------------------------------------------------------
      rdt_init = Output_mp%rdt
      tracer_orig = Input_mp%tracer

!--------------------------------------------------------------------
!    initialize some variables. define the inverse of the time step.
!--------------------------------------------------------------------
      conv_calc_completed = .false.
      available_cf_for_uw = 1.0
      dtinv = 1.0/dt

!--------------------------------------------------------------------
!    initialize various arrays which are used in this subroutine.
!--------------------------------------------------------------------
      t_ref = 0.
      q_ref = 0.
      ttnd_don = 0.
      qtnd_don = 0.
      cmf = 0.
      delta_temp = 0.
      delta_q    = 0.
      delta_vapor = 0.
      delta_qn   = 0.
      delta_qni  = 0.
      liquid_precip = 0.
      frozen_precip = 0.
      det0 = 0.
      det_cmt = 0.
      mc_donner = 0.
      mc_donner_half = 0.
      mc_donner_up = 0.
      m_cdet_donner = 0.
      massflux = 0.
      qtr = 0.
     
      rain_uw = 0.
      snow_uw = 0.
      ttnd_uw = 0.
      qtnd_uw = 0.
      utnd_uw = 0.
      vtnd_uw = 0.
      qltnd_uw = 0.
      qitnd_uw = 0.
      qatnd_uw = 0.
      qntnd_uw = 0.
      qnitnd_uw = 0.
      qtruw = 0.
      delta_ql = 0.
      delta_qi = 0.
      delta_qa = 0.
 
      rain_donmca  = 0.0
      snow_donmca  = 0.0

      precip       = 0.0 
      rain3d       = 0.0
      snow3d       = 0.0

!----------------------------------------------------------------------
!    output any requested convectively-transported tracer fields 
!    and / or their column sums before convective transport.
!----------------------------------------------------------------------
      do n=1, num_prog_tracers
        used = send_data (id_conv_tracer(n), Input_mp%tracer(:,:,:,n),  &
                                                        Time, is, js, 1)
        if (id_conv_tracer_col(n) > 0)  &
          call column_diag(id_conv_tracer_col(n), is, js, Time, &
                           Input_mp%tracer(:,:,:,n), 1.0, Input_mp%pmass) 
      end do

!---------------------------------------------------------------------
!    begin the convective parameterizations clock.
!---------------------------------------------------------------------
      call mpp_clock_begin (convection_clock)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                   DRY CONVECTION PARAMETERIZATION
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!---------------------------------------------------------------------
!    if dry adjustment is desired call subroutine dry_adj to obtain
!    the temperature tendencies which must be applied to adjust each
!    column to a non-superadiabatic lapse rate. 
!---------------------------------------------------------------------
      if (do_dryadj) then
        call dry_adj (Input_mp%tin, Input_mp%pfull, Input_mp%phalf,   &
                                                 delta_temp, mask)

!-------------------------------------------------------------------
!    add the temperature change due to dry adjustment to the current
!    temperature. convert the temperature change to a heating rate and
!    add that to the temperature tendency array accumulating the ten-
!    dencies due to all physics processes.
!-------------------------------------------------------------------
        Input_mp%tin  = Input_mp%tin + delta_temp
        Tend_mp%ttnd  = delta_temp*dtinv
        Output_mp%tdt = Output_mp%tdt + Tend_mp%ttnd

!---------------------------------------------------------------------
!    output the temperature tendency from dry adjustment, if desired.
!---------------------------------------------------------------------
        used = send_data (id_tdt_dadj, Tend_mp%ttnd, Time, is, js, 1 )

!---------------------------------------------------------------------
!    add the temperature time tendency from dry adjustment to the array
!    accumulating the total temperature time tendency from convection.
!---------------------------------------------------------------------
        Tend_mp%ttnd_conv = Tend_mp%ttnd_conv + Tend_mp%ttnd
      endif


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                  MOIST CONVECTION PARAMETERIZATIONS
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                0. UW SHALLOW CONVECTION PARAMETERIZATION
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!---------------------------------------------------------------------
!    if it is desired to do shallow convection before donner deep 
!    convection, turn on the shallowcu clock and make the call to
!    uw_conv.
!---------------------------------------------------------------------
      if (.not. do_donner_before_uw) then
        call mpp_clock_begin (shallowcu_clock)
        if (do_uw_conv) then

!---------------------------------------------------------------------
!    be sure the shallow cloud field arguments associated with the uw_conv
!    parameterization are present.
!---------------------------------------------------------------------
          if (i_shallow /= 0) then
          else
            call error_mesg ('convection_driver_mod', 'convection_driver: &
                       &improper arguments for uw_shallow clouds', FATAL)
          endif

!---------------------------------------------------------------------
!    call the uw_conv wrapper routine.
!---------------------------------------------------------------------
          call moistproc_uw_conv(      &
               Time, is, ie, js, je, dt, Input_mp%tin, Input_mp%qin, &
               Input_mp%uin, Input_mp%vin, Input_mp%tracer,    &
               Input_mp%pfull, Input_mp%phalf, Input_mp%zfull,   &
               Input_mp%zhalf, Input_mp%omega, Input_mp%pblht,        &
               Input_mp%ustar, Input_mp%bstar, Input_mp%qstar,   &
               shflx, lhflx,         &
               Input_mp%land, Input_mp%coldT, Aerosol,       &
               Surf_diff%tdt_rad, Surf_diff%tdt_dyn, Surf_diff%qdt_dyn, &
               Surf_diff%dgz_dyn, Surf_diff%ddp_dyn, tdt_dif, qdt_dif, &
               Phys_mp_exch%hmint, Input_mp%lat, Input_mp%lon,   &
               Input_mp%cush, Input_mp%cbmf, Phys_mp_exch%cgust,   &
               Phys_mp_exch%tke, Phys_mp_exch%pblhto, &
               Phys_mp_exch%rkmo, Phys_mp_exch%taudpo, &
               Phys_mp_exch%exist_shconv, Phys_mp_exch%exist_dpconv, &
               Phys_mp_exch%pblht_prev, Phys_mp_exch%hlsrc_prev, &
               Phys_mp_exch%qtsrc_prev, Phys_mp_exch%cape_prev, &
               Phys_mp_exch%cin_prev, Phys_mp_exch%tke_prev, &
               cmf, conv_calc_completed,   &
               available_cf_for_uw, Output_mp%tdt, Output_mp%rdt(:,:,:,1),&
               Output_mp%udt, Output_mp%vdt, Output_mp%rdt,    &
               Tend_mp%ttnd_conv, Tend_mp%qtnd_conv, Output_mp%lprec,   &
               Output_mp%fprec, precip, Removal_mp%liq_precflx,  &
               Removal_mp%ice_precflx, rain_uw, snow_uw, ttnd_uw,    &
               qtnd_uw, utnd_uw, vtnd_uw, qtruw, qltnd_uw, qitnd_uw,   &
               qatnd_uw, qntnd_uw, qnitnd_uw, doing_prog_clouds,   &
               do_limit_uw, do_liq_num, num_prog_tracers,  &
               Removal_mp_control%tracers_in_uw,    &
               Removal_mp_control%num_uw_tracers,   &
               Moist_clouds_block%cloud_data(i_shallow),  &
               Removal_mp%uw_wetdep, do_ice_num, detrain_ice_num)      
        endif  !(       do_uw_conv)
        call mpp_clock_end   (shallowcu_clock)
      else
!----------------------------------------------------------------------
!    save the original input values of t and q for use when uw is called 
!    in case the option to not use updated profiles in uw is activated. 
!    (otherwise these fields will be modified by any other calls to 
!    convective processes made before the uw_conv call.)
!----------------------------------------------------------------------
        tin_orig = Input_mp%tin
        qin_orig = Input_mp%qin
      endif  ! (.not do_donner_before_uw)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                A. DONNER DEEP CONVECTION PARAMETERIZATION
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!---------------------------------------------------------------------
!    if donner_deep convection is activated, execute the following code.
!    activate the donner clock.
!---------------------------------------------------------------------
      if (do_donner_deep) then
        call mpp_clock_begin (donner_clock)

!---------------------------------------------------------------------
!    be sure the donner cloud field arguments are present.
!---------------------------------------------------------------------
        if (i_cell /= 0 .and. i_meso /= 0 ) then  
        else
          call error_mesg ('convection_driver_mod',   &
               'input args for donner clouds not correct', FATAL)
        endif

!--------------------------------------------------------------------
!    if prognostic clouds are present, define the cloud liquid and 
!    cloud ice specific humidities and cloud area associated with them,
!    so that they may be input to donner_deep_mod. if not using prognostic
!    clouds,  define these arrays to be zero. 
!--------------------------------------------------------------------
        if (doing_prog_clouds) then
          qlin = Input_mp%tracer(:,:,:,nql)
          qiin = Input_mp%tracer(:,:,:,nqi)
          qain = Input_mp%tracer(:,:,:,nqa)
          if (do_liq_num ) nllin =  Input_mp%tracer(:,:,:,nqn)
          if (do_ice_num ) nilin =  Input_mp%tracer(:,:,:,nqni)
        endif

!--------------------------------------------------------------------
!    convert vapor specific humidity to vapor mixing ratio so it may
!    be input to donner_deep_mod.
!--------------------------------------------------------------------
        rin = Input_mp%qin/(1.0 - Input_mp%qin)

!---------------------------------------------------------------------
!    if any tracers are to be transported by donner convection, 
!    check each active tracer to find those to be transported and fill 
!    the donner_tracer array with these fields.
!---------------------------------------------------------------------
        if (Removal_mp_control%num_donner_tracers > 0) then
          nn = 1
          do n=1,num_prog_tracers
            if (Removal_mp_control%tracers_in_donner(n)) then
              donner_tracer(:,:,:,nn) = Input_mp%tracer(:,:,:,n)
              nn = nn + 1
            endif
          end do
        else
          donner_tracer(:,:,:,:) = 0.0  
        endif

!---------------------------------------------------------------------
!  NOTE 1: sfc_sh_flux, sfc_vapor_flux, tr_flux are the surface fluxes
!          that will have been obtained from the flux exchange module
!          and passed on to moist_processes and then to donner_deep.
!          FOR NOW, these values are defined herein, and given
!          values of 0.0. Could / should be passed in appropriately 
!          from the Surf_diff variable.
!---------------------------------------------------------------------
        sfc_sh_flux    = 0.0
        sfc_vapor_flux = 0.0
        if (Removal_mp_control%num_donner_tracers > 0) then
          nn = 1
          do n=1, num_prog_tracers
            if (Removal_mp_control%tracers_in_donner(n)) then
              tr_flux(:,:,nn) = 0.0                                 
              nn = nn + 1
            endif
          end do
        else
          tr_flux = 0.
        endif

!-----------------------------------------------------------------------
!    define boundary layer kinetic energy to pass to donner deep routine.
!-----------------------------------------------------------------------
        temp_2d = Input_mp%pblht
        temp_2d = min(max(temp_2d, 0.0),5000.)
        temp_2d = Input_mp%ustar**3. +   &
                             0.6*Input_mp%ustar*Input_mp%bstar*temp_2d
        where (temp_2d .gt. 0.)
          temp_2d = temp_2d**(2./3.)
        end where
        temp_2d = MAX (1.e-6, temp_2d)

!----------------------------------------------------------------------
!    define model time in days and secs from base time.
!----------------------------------------------------------------------
        call get_time (Time, secs, days)

!---------------------------------------------------------------------
!    call donner_deep to compute the effects of deep convection on the 
!    temperature, vapor mixing ratio, tracers, cloud liquid, cloud ice
!    cloud area and precipitation fields.
!---------------------------------------------------------------------
        if (doing_prog_clouds) then
          call donner_deep     &
               (is, ie, js, je, dt, Input_mp%tin, rin, Input_mp%pfull,  &
                Input_mp%phalf, Input_mp%zfull, Input_mp%zhalf,    &
                Input_mp%omega, Input_mp%pblht, temp_2d, Input_mp%qstar, &
                Input_mp%cush, Input_mp%coldT, Input_mp%land,     &
                sfc_sh_flux, sfc_vapor_flux, tr_flux, donner_tracer,   &
                secs, days, Input_mp%cbmf,            &
                Moist_clouds_block%cloud_data(i_cell)%cloud_area, &
                Moist_clouds_block%cloud_data(i_cell)%liquid_amt, &
                Moist_clouds_block%cloud_data(i_cell)%liquid_size, &
                Moist_clouds_block%cloud_data(i_cell)%ice_amt    , &
                Moist_clouds_block%cloud_data(i_cell)%ice_size   , &
                Moist_clouds_block%cloud_data(i_cell)%droplet_number, &
                Moist_clouds_block%cloud_data(i_meso)%cloud_area, &
                Moist_clouds_block%cloud_data(i_meso)%liquid_amt, &
                Moist_clouds_block%cloud_data(i_meso)%liquid_size, &
                Moist_clouds_block%cloud_data(i_meso)%ice_amt    , &
                Moist_clouds_block%cloud_data(i_meso)%ice_size   , &
                Moist_clouds_block%cloud_data(i_meso)%droplet_number, &
                Moist_clouds_block%cloud_data(i_meso)%nsum_out, &
                maxTe_launch_level,   &
                precip_returned, delta_temp, delta_vapor,   &
                m_cdet_donner, m_cellup, mc_donner, mc_donner_up,    &
                mc_donner_half, C2ls_mp%donner_humidity_area,    &
                C2ls_mp%donner_humidity_factor, qtr,  &
                Removal_mp%donner_wetdep, lheat_precip, vert_motion,    &
                total_precip, liquid_precip, frozen_precip, &
                Removal_mp%frz_meso,  Removal_mp%liq_meso, &
                Removal_mp%frz_cell, Removal_mp%liq_cell, &
                qlin, qiin, qain, delta_ql,                 &!optional
                delta_qi, delta_qa)                          !optional  
        else
          call donner_deep    &
               (is, ie, js, je, dt, Input_mp%tin, rin, Input_mp%pfull,   &
                Input_mp%phalf, Input_mp%zfull, Input_mp%zhalf,   &
                Input_mp%omega, Input_mp%pblht, temp_2d, Input_mp%qstar, &
                Input_mp%cush, Input_mp%coldT, Input_mp%land,   &
                sfc_sh_flux, sfc_vapor_flux, tr_flux, donner_tracer,   &
                secs, days,  Input_mp%cbmf,           &
                Moist_clouds_block%cloud_data(i_cell)%cloud_area, &
                Moist_clouds_block%cloud_data(i_cell)%liquid_amt, &
                Moist_clouds_block%cloud_data(i_cell)%liquid_size, &
                Moist_clouds_block%cloud_data(i_cell)%ice_amt    , &
                Moist_clouds_block%cloud_data(i_cell)%ice_size   , &
                Moist_clouds_block%cloud_data(i_cell)%droplet_number, &
                Moist_clouds_block%cloud_data(i_meso)%cloud_area, &
                Moist_clouds_block%cloud_data(i_meso)%liquid_amt, &
                Moist_clouds_block%cloud_data(i_meso)%liquid_size, &
                Moist_clouds_block%cloud_data(i_meso)%ice_amt    , &
                Moist_clouds_block%cloud_data(i_meso)%ice_size   , &
                Moist_clouds_block%cloud_data(i_meso)%droplet_number, &
                Moist_clouds_block%cloud_data(i_meso)%nsum_out, &
                maxTe_launch_level,   &
                precip_returned, delta_temp, delta_vapor,   &
                m_cdet_donner, m_cellup, mc_donner, mc_donner_up,    &
                mc_donner_half, C2ls_mp%donner_humidity_area,    &
                C2ls_mp%donner_humidity_factor, qtr,    &
                Removal_mp%donner_wetdep, lheat_precip, vert_motion,    &
                total_precip, liquid_precip, frozen_precip, &
                Removal_mp%frz_meso,  Removal_mp%liq_meso,  &
                Removal_mp%frz_cell, Removal_mp%liq_cell)
        endif

!---------------------------------------------------------------------
!    update the current timestep tracer changes with the contributions 
!    just obtained from donner transport.
!---------------------------------------------------------------------
        if (Removal_mp_control%num_donner_tracers > 0) then
          nn = 1
          do n=1, num_prog_tracers
            if (Removal_mp_control%tracers_in_donner(n)) then
              Output_mp%rdt(:,:,:,n) = Output_mp%rdt(:,:,:,n) +   &
                                                             qtr(:,:,:,nn)
              nn = nn + 1
            endif
          end do
        endif

!--------------------------------------------------------------------
!    obtain updated vapor specific humidity (qnew) resulting from deep 
!    convection. define the vapor specific humidity change due to deep 
!    convection (qtnd).
!--------------------------------------------------------------------
        do k=1,kx
          do j=1,jx
            do i=1,ix
              if (delta_vapor(i,j,k) /= 0.0) then
!was qnew... now temp
                temp = (rin(i,j,k) + delta_vapor(i,j,k))/   &
                       (1.0 + (rin(i,j,k) + delta_vapor(i,j,k)))
                delta_q(i,j,k) = temp - Input_mp%qin(i,j,k)
              else
                delta_q(i,j,k) = 0.
              endif
            enddo
          enddo
        end do

!------------------------------------------------------------------------
!    if conservation checks on the water and enthalpy changes produced
!    within the donner deep convectuion scheme have been requested, 
!    compute the needed vertical integrals here.
!------------------------------------------------------------------------
        if (do_donner_conservation_checks) then
          vaporint = 0.
          lcondensint = 0.
          condensint = 0.
          diffint = 0.
          enthint = 0.
          enthdiffint = 0.
    
          do k=1,kx
            vaporint = vaporint + Input_mp%pmass (:,:,k)*delta_q(:,:,k)
            enthint = enthint + CP_AIR*     &
                                   Input_mp%pmass(:,:,k)*delta_temp(:,:,k)
            condensint = condensint + Input_mp%pmass(:,:,k) *  &
                                       (delta_ql(:,:,k) + delta_qi(:,:,k))
            lcondensint = lcondensint + Input_mp%pmass(:,:,k) *  &
                               (HLV*delta_ql(:,:,k) + HLS*delta_qi(:,:,k))
          end do
          precipint = total_precip/seconds_per_day
          diffint = (vaporint + condensint)*dtinv  + precipint
          enthdiffint = (enthint - lcondensint)*dtinv -    &
                                      lheat_precip/seconds_per_day -     &
                                            vert_motion/seconds_per_day 

!------------------------------------------------------------------------
!    update the variable collecting the maximum imbalance over the entire
!    model run, if the present imbalance value is larger than the 
!    previously recorded.
!------------------------------------------------------------------------
          do j=1,size(enthdiffint,2)
            do i=1,size(enthdiffint,1)
              max_enthalpy_imbal_don(i+is-1,j+js-1) =    &
                         max( abs(enthdiffint(i,j)), &
                                    max_enthalpy_imbal_don(i+is-1,j+js-1) )
              max_water_imbal_don(i+is-1,j+js-1) =     &
                         max( abs(diffint(i,j)), &
                                    max_water_imbal_don(i+is-1,j+js-1) )
            end do
          end do

!------------------------------------------------------------------------
!    output diagnostics related to water and enthalpy conservation.
!------------------------------------------------------------------------
          used = send_data(id_max_enthalpy_imbal_don,    &
                       max_enthalpy_imbal_don(is:ie,js:je), Time, is, js)
          used = send_data(id_max_water_imbal_don,     &
                          max_water_imbal_don(is:ie,js:je), Time, is, js)
          used = send_data(id_vaporint, vaporint*dtinv, Time, is, js)
          used = send_data(id_condensint, condensint*dtinv, Time, is, js)
          used = send_data(id_vertmotion, vert_motion/seconds_per_day,  &
                                                             Time, is, js)
          used = send_data(id_precipint, precipint, Time, is, js)
          used = send_data(id_diffint, diffint, Time, is, js)
          used = send_data(id_enthint, enthint*dtinv, Time, is, js)
          used = send_data(id_lcondensint, lcondensint*dtinv, Time, is, js)
          used = send_data(id_lprcp, lheat_precip/seconds_per_day,   &
                                                              Time, is, js)
          used = send_data(id_enthdiffint, enthdiffint, Time, is, js)
        endif

!---------------------------------------------------------------------
!    scale Donner tendencies to prevent the formation of negative
!    total water specific humidities
!---------------------------------------------------------------------
        if (doing_prog_clouds .and. do_limit_donner) then
          call moistproc_scale_donner (   &
                is, ie, js, je, dt, Input_mp%qin, delta_temp, delta_q, &
                precip_returned, total_precip, lheat_precip,    &
                liquid_precip, frozen_precip, Input_mp%pmass,   &
                num_prog_tracers,   &
                Removal_mp_control%tracers_in_donner, delta_ql, delta_qi, &
                delta_qa, qlin, qiin, qtr, scale_donner)
        else
          scale_donner = 1.0
        end if 
        used = send_data (id_scale_donner, scale_donner, Time, is, js )

!---------------------------------------------------------------------
!    recalculate the precip using the delta specific humidity tenden-
!    cies. define precip_adjustment as the change in precipitation 
!    resulting from the recalculation.
!---------------------------------------------------------------------
        if (force_donner_moist_conserv) then

!---------------------------------------------------------------------
!    calculate the adjustment to the precipitation that would be needed 
!    in order to balance the change in water content in the column.
!    if this is smaller than 1.e-10, ignore the imbalance. Also calculate
!    the fraction of the returned precip that this corresponds to.
!---------------------------------------------------------------------
          temp_2d = 0.
          do k=1,kx
            temp_2d (:,:) = temp_2d (:,:) + (-delta_q(:,:,k) -  &
                                delta_ql(:,:,k) -delta_qi(:,:,k))*  &
                                                    Input_mp%pmass(:,:,k)
          end do
          precip_adjustment = (temp_2d - precip_returned)

          do j=1,jx
            do i=1,ix
              if (ABS(precip_adjustment(i,j)) < 1.0e-10) then
                precip_adjustment (i,j) = 0.0
              endif
              if ( precip_adjustment(i,j) < 0.0 .and. &
                 (precip_adjustment(i,j)+precip_returned(i,j)) < 0.0 ) then
! If precip_returned is greater than the "change in water content" balance,
!and
! there is not enough water available to beg/borrow/steal from, we need to zero
! out the various tendencies. i.e. donner_deep will not contribute to changing 
! the column 
!           write (warn_mesg,'(2i4,2e12.4)') i,j,precip_adjustment(i,j), precip_returned(i,j)
!           call error_mesg ('moist_processes_mod', 'moist_processes: &
!                 &Change in water content does not balance precip &
!                 &from donner_deep routine.'//trim(warn_mesg), WARNING)
                scale(i,j) = 0.0
                delta_vapor(i,j,:) = 0.0
                delta_q(i,j,:) = 0.0
                delta_qi(i,j,:) = 0.0
                delta_ql(i,j,:) = 0.0
                delta_qa(i,j,:) = 0.0
                total_precip(i,j) = 0.0
                precip_returned(i,j) = 0.0
                liquid_precip(i,j,:) = 0.0
                frozen_precip(i,j,:) = 0.0
                lheat_precip(i,j) = 0.0
              endif
            end do
          end do
          do j=1,jx
            do i=1,ix
              if (precip_returned(i,j) > 0.0) then
                adjust_frac(i,j) =     &
                               precip_adjustment(i,j)/precip_returned(i,j)
              else
                adjust_frac(i,j) = 0.
              endif
            end do
          end do

!----------------------------------------------------------------------
!    now adjust the temperature to balance the precip adjustment
!    and so conserve enthalpy in the column, and  define the new values
!    of liquid and frozen precipitation after adjustment.
!--------------------------------------------------------------------- 
          do k=1,kx
            ttnd_adjustment(:,:,k) = &
                      ((HLV*liquid_precip(:,:,k)*adjust_frac(:,:) + &
                        HLS*frozen_precip(:,:,k)*adjust_frac(:,:))  &
                       *dt/seconds_per_day)/CP_AIR
            liquid_precip(:,:,k) = liquid_precip(:,:,k)*  &
                                                  (1.0+adjust_frac(:,:))
            frozen_precip(:,:,k) = frozen_precip(:,:,k)*   &
                                                  (1.0+adjust_frac(:,:))
          end do
        else ! (force_donner_moist_conserv)

!------------------------------------------------------------------------
!    define the adjustments to be 0.0 if conservation is not being forced.
!------------------------------------------------------------------------
          precip_adjustment = 0.0
          adjust_frac       = 0.0
          ttnd_adjustment = 0.
        endif  ! (force_donner_moist_conserv)

!-------------------------------------------------------------------------
!    define the column rainfall and snowfall from the donner scheme.
!-------------------------------------------------------------------------
        rain_don = 0.0
        snow_don = 0.0
        do k=1,kx
          rain_don = rain_don + liquid_precip(:,:,k)*  &
                                    Input_mp%pmass(:,:,k)/seconds_per_day
          snow_don = snow_don + frozen_precip(:,:,k)*  &
                                    Input_mp%pmass(:,:,k)/seconds_per_day
        end do

!----------------------------------------------------------------------
!   modify the 3d precip fluxes used by COSP to account for the 
!   conservation adjustment.
!----------------------------------------------------------------------
        if (do_cosp) then
          do k=1, kx
            do j=1,jx  
              do i=1,ix  
                Removal_mp%frz_meso(i,j,k) = Removal_mp%frz_meso(i,j,k)*  &
                                 Input_mp%pmass(i,j,k)*scale(i,j)* &
                                     (1.0+adjust_frac(i,j))/SECONDS_PER_DAY
                Removal_mp%liq_meso(i,j,k) = Removal_mp%liq_meso(i,j,k)*  &
                                 Input_mp%pmass(i,j,k)*scale(i,j)* &
                                     (1.0+adjust_frac(i,j))/SECONDS_PER_DAY
                Removal_mp%frz_cell(i,j,k) = Removal_mp%frz_cell(i,j,k)*  &
                                 Input_mp%pmass(i,j,k)*scale(i,j)* &
                                     (1.0+adjust_frac(i,j))/SECONDS_PER_DAY
                Removal_mp%liq_cell(i,j,k) = Removal_mp%liq_cell(i,j,k)*  &
                                 Input_mp%pmass(i,j,k)*scale(i,j)* &
                                     (1.0+adjust_frac(i,j))/SECONDS_PER_DAY
              end do
            end do
          end do
        endif
    
!-------------------------------------------------------------------------
!    if the option to allow only one convective scheme per column is
!    active, mark those columns which underwent donner convection.
!-------------------------------------------------------------------------
        if (only_one_conv_scheme_per_column) then
          conv_calc_completed = (rain_don + snow_don) > 0.0
        endif

!---------------------------------------------------------------------
!    convert the deltas in temperature, vapor specific humidity and 
!    precipitation resulting from deep convection to time tendencies 
!    of these quantities.
!---------------------------------------------------------------------
        ttnd_don = delta_temp*dtinv 
        ttnd_don = ttnd_don + ttnd_adjustment*dtinv
        qtnd_don = delta_q*dtinv

!--------------------------------------------------------------------
!    add the tendencies of temperature and specific humidity resulting
!    from the deep convection component of the donner parameterization
!    to the total convective tendency. 
!--------------------------------------------------------------------
        Tend_mp%ttnd_conv = Tend_mp%ttnd_conv + ttnd_don
        Tend_mp%qtnd_conv = Tend_mp%qtnd_conv + qtnd_don

!--------------------------------------------------------------------
!    add the tendencies of temperature and specific humidity resulting
!    from the deep convection component of the donner parameterization
!    to the total tendencies from all physics processes.
!--------------------------------------------------------------------
        Output_mp%tdt = Output_mp%tdt + ttnd_don 
        Output_mp%rdt(:,:,:,1) = Output_mp%rdt(:,:,:,1) + qtnd_don

!--------------------------------------------------------------------
!    add the liquid (rain) and frozen (snow) precipitation generated by
!    deep convection on this step to the arrays accumulating precip-
!    itation from all sources (lprec, fprec).
!--------------------------------------------------------------------

!if (minval(rain_don) < 0.0 ) then
!           write (warn_mesg,'(2i4,e12.4)') minloc(rain_don), minval(rain_don)
!           call error_mesg ('moist_processes_mod', 'moist_processes: &
!                 &Donner_deep rain is negative.'//trim(warn_mesg), WARNING)
!endif
!if (minval(snow_don) < 0.0 )  then
!           write (warn_mesg,'(2i4,e12.4)') minloc(snow_don), minval(snow_don)
!           call error_mesg ('moist_processes_mod', 'moist_processes: &
!                 &Donner_deep snow is negative.'//trim(warn_mesg), WARNING)
!endif

        Output_mp%lprec  = Output_mp%lprec + rain_don
        Output_mp%fprec  = Output_mp%fprec + snow_don

!--------------------------------------------------------------------
!    output diagnostics for the time tendencies of temperature, vapor 
!    specific humidity and large scale cloud fields, and various precip 
!    and mass flux diagnostics due to donner deep convection.
!--------------------------------------------------------------------
        used = send_data (id_tdt_deep_donner, ttnd_don, Time, is, js, 1)
        used = send_data (id_qdt_deep_donner, qtnd_don, Time, is, js, 1)
        used = send_data (id_qadt_deep_donner, delta_qa*dtinv,   &
                                                        Time, is, js, 1)
        used = send_data (id_qldt_deep_donner, delta_ql*dtinv,    &
                                                        Time, is, js, 1)
        used = send_data (id_qidt_deep_donner, delta_qi*dtinv,    &
                                                        Time, is, js, 1)

        used = send_data (id_mc_donner, mc_donner, Time, is, js, 1)
        used = send_data (id_mc_donner_half, mc_donner_half,    &
                                                        Time, is, js, 1 )
        used = send_data (id_m_cdet_donner, m_cdet_donner,    &
                                                        Time,  is, js, 1 )
        used = send_data (id_m_cellup, m_cellup, Time, is, js, 1 )
        used = send_data (id_snow_deep_donner, snow_don, Time, is, js)
        used = send_data (id_prec_deep_donner, rain_don + snow_don,   &
                                                        Time, is, js )
        used = send_data (id_prec1_deep_donner, precip_adjustment,  &
                              Time, is, js, mask = precip_returned > 0.0)
        used = send_data (id_precret_deep_donner, precip_returned,  &
                              Time, is, js)       

!------------------------------------------------------------------------
!    if donner conservation checks have been done, output various
!    diagnostics describing the results. 
!------------------------------------------------------------------------
        if (do_donner_conservation_checks) then
          used = send_data (id_enth_donner_col2, -hlv*rain_don,    &
                                                            Time, is, js)
          used = send_data (id_enth_donner_col3, -hls*snow_don,    &
                                                            Time, is, js)
          if (id_enth_donner_col4 > 0)   &
                     call column_diag(id_enth_donner_col4, is, js, Time, &
                                   ttnd_don(:,:,:), CP_AIR, Input_mp%pmass)
          if (id_enth_donner_col5 > 0)    &
                     call column_diag(id_enth_donner_col5, is, js, Time, &
                             delta_ql(:,:,:), -HLV*dtinv,   &
                             delta_qi(:,:,:), -HLS*dtinv, Input_mp%pmass)
          if (id_enth_donner_col6 > 0)     &
                     call column_diag(id_enth_donner_col6, is, js, Time, &
                                 ttnd_adjustment, CP_AIR, Input_mp%pmass)
          used = send_data (id_enth_donner_col7, adjust_frac, Time, is, js)
       
!------------------------------------------------------------------------
!    compute and output column enthalpy change due to donner deep 
!    convection.
!------------------------------------------------------------------------
          temp_2d = 0.
          do k=1,kx
            temp_2d(:,:) = temp_2d(:,:)   + &
                          (-HLV*liquid_precip(:,:,k)/seconds_per_day -  &
                            hls*frozen_precip(:,:,k)/seconds_per_day  + &
                            CP_AIR*ttnd_don(:,:,k)  -  &
                          (HLV*delta_ql(:,:,k)*dtinv +   &
                           HLS*delta_qi(:,:,k)*dtinv))*   &
                                                     Input_mp%pmass(:,:,k)
          end do
          used = send_data (id_enth_donner_col, temp_2d, Time, is, js)

!------------------------------------------------------------------------
!    compute and output column water change due to donner deep convection.
!------------------------------------------------------------------------
          if (id_wat_donner_col > 0) then
            temp_2d = rain_don + snow_don
            call column_diag(id_wat_donner_col, is, js, Time, qtnd_don, &
                             1.0, delta_ql, dtinv, delta_qi, dtinv, &
                             Input_mp%pmass, temp_2d)
          endif
        endif ! (donner_conservation_checks)

!------------------------------------------------------------------------
!    output additional diagnostics related to the clouds associated with
!    donner convection.
!------------------------------------------------------------------------
        used = send_data (id_cell_cld_frac,   &
                     Moist_clouds_block%cloud_data(i_cell)%cloud_area, &
                                                         Time, is, js, 1 )
        used = send_data (id_meso_cld_frac,   &
                     Moist_clouds_block%cloud_data(i_meso)%cloud_area, &
                                                         Time, is, js, 1)
        used = send_data (id_donner_humidity_area,    &
                     C2ls_mp%donner_humidity_area(:,:,:), Time, is, js, 1 )

!1/19/16 MOD RSH
!---------------------------------------------------------------------
!    update the values of temperature and vapor specific humidity to
!    include the effects of deep convection. if mca was included in the
!    donner deep scheme, then this update has already been done.
!---------------------------------------------------------------------
        if (keep_icenum_detrain_bug) then 
        else
           Input_mp%tin(:,:,:) = Input_mp%tin(:,:,:) + delta_temp(:,:,:)
           Input_mp%qin(:,:,:) = Input_mp%qin(:,:,:) + delta_q(:,:,:)
        endif

!-----------------------------------------------------------------------
!    turn off the donner clock.
!-----------------------------------------------------------------------
        call mpp_clock_end (donner_clock)

!----------------------------------------------------------------------
!    this section calculates the moist convective adjustment associated
!    with the donner convection scheme. It is activated / deactivated
!    by moist_processes_nml variable do_donner_mca.
!----------------------------------------------------------------------
        if (do_donner_mca) then

!----------------------------------------------------------------------
!    if donner mca is active, turn on its clock.
!----------------------------------------------------------------------
          call mpp_clock_begin (donner_mca_clock)

!--------------------------------------------------------------------
!    call subroutine moist_conv to handle any shallow convection present 
!    in the grid. this call is made without the optional lsc variables so 
!    that no convective detrainment (and corresponding change in 
!    large-scale cloud amount and area) occurs, consistent with this call 
!    being intended to handle only shallow convection. The temp and vapor
!    fields are updated with any changes from deep convection before the 
!    routine is called.
!--------------------------------------------------------------------
          if (keep_icenum_detrain_bug) then
            Input_mp%tin = Input_mp%tin + delta_temp
            Input_mp%qin = Input_mp%qin + delta_q
          endif
          call moist_conv (    &
                  Input_mp%tin, Input_mp%qin, Input_mp%pfull,  &
                  Input_mp%phalf, Input_mp%coldT, ttnd_don, qtnd_don, &
                  rain_donmca, snow_donmca, dtinv, Time, is, js,     &
                  donner_tracer, qtr, Lbot=kbot, mask=mask)

!-----------------------------------------------------------------------
!    if the effects of the mca component of donner are to be seen by
!    COSP, define the associated precip fluxes.
!-----------------------------------------------------------------------
          if (do_cosp .and. include_donmca_in_cosp) then
            do j=1,jx 
              do i=1,ix 
                if (Input_mp%coldT(i,j)) then
                  do k=1,kx
                    Removal_mp%mca_frz(i,j,k) =    &
                               -1.0*qtnd_don(i,j,k)*Input_mp%pmass(i,j,k)
                    Removal_mp%mca_liq(i,j,k) = 0.
                  end do
                else
                  do k=1,kx
                    Removal_mp%mca_frz(i,j,k) = 0.
                    Removal_mp%mca_liq(i,j,k) =     &
                               -1.0*qtnd_don(i,j,k)*Input_mp%pmass(i,j,k)
                  end do
                endif
              end do
            end do
          else
            Removal_mp%mca_frz = 0.
            Removal_mp%mca_liq = 0.
          endif

!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    just obtained from moist convective adjustment. currently there
!    is no tracer transport by this process, so qtr will be 0.0 for all
!    tracers.
!---------------------------------------------------------------------
          nn = 1
          do n=1, num_prog_tracers
            if (Removal_mp_control%tracers_in_donner(n)) then
              Output_mp%rdt(:,:,:,n) = Output_mp%rdt(:,:,:,n) +   &
                                                           qtr(:,:,:,nn)
              nn = nn + 1
            endif
          end do

!--------------------------------------------------------------------
!    add the heating and moistening rates from the mca portion of
!    donner convection to the arrays accumulating total convective
!    changes and those accumulating total physics changes.
!--------------------------------------------------------------------
          Tend_mp%ttnd_conv = Tend_mp%ttnd_conv + ttnd_don
          Tend_mp%qtnd_conv = Tend_mp%qtnd_conv + qtnd_don
          Output_mp%tdt = Output_mp%tdt + ttnd_don
          Output_mp%rdt(:,:,:,1) = Output_mp%rdt(:,:,:,1) + qtnd_don

!--------------------------------------------------------------------
!    add the liquid (rain) and frozen (snow) precipitation generated by
!    the moist convective adjustment pass of the donner parameterization
!    on this step to the arrays accumulating precipitation from all 
!    sources (lprec, fprec).
!--------------------------------------------------------------------
          Output_mp%lprec  = Output_mp%lprec + rain_donmca
          Output_mp%fprec  = Output_mp%fprec + snow_donmca

!--------------------------------------------------------------------
!    output the time tendencies of temperature, vapor specific humid-
!    ity, precipitation and snow due to the moist convective 
!    adjustment pass of the donner parameterization.
!--------------------------------------------------------------------
          used = send_data (id_tdt_mca_donner, ttnd_don, Time, is, js, 1)
          used = send_data (id_qdt_mca_donner, qtnd_don, Time, is, js, 1)
          used = send_data (id_prec_mca_donner, rain_donmca+snow_donmca, &
                                                             Time, is, js)
          used = send_data (id_snow_mca_donner, snow_donmca, Time, is, js)

!------------------------------------------------------------------------
!    output the column imbalances of enthalpy and water resulting from the
!    mca component of donner convection. 
!------------------------------------------------------------------------
          if (id_enth_mca_donner_col > 0) then
            temp_2d = -HLV*rain_donmca -HLS*snow_donmca
            call column_diag(id_enth_mca_donner_col, is, js, Time,   &
                               ttnd_don, CP_AIR, Input_mp%pmass, temp_2d)
          endif

          if (id_wat_mca_donner_col > 0) then
            temp_2d = rain_donmca + snow_donmca
            call column_diag(id_wat_mca_donner_col, is, js, Time,   &
                               qtnd_don, 1.0,  Input_mp%pmass, temp_2d)
          endif

!--------------------------------------------------------------------
!    output the time tendencies of tracer and of column tracer 
!    due to the moist convective adjustment pass of the donner 
!    parameterization. currently moist convective adjustment does not
!    affect the tracer fields, so these fields are always 0.0.
!--------------------------------------------------------------------
          do n = 1, Removal_mp_control%num_donner_tracers
            if ( id_tracerdt_mcadon(n) > 0 ) &
              used = send_data(id_tracerdt_mcadon(n), qtr(:,:,:,n),   &
                                                         Time, is, js, 1 )
            if (id_tracerdt_mcadon_col(n) > 0 )  &
              call column_diag(id_tracerdt_mcadon_col(n), is, js, Time, &
                                      qtr(:,:,:,n), 1.0, Input_mp%pmass)
          enddo

!-----------------------------------------------------------------------
!    turn off the donner mca clock.
!-----------------------------------------------------------------------
          call mpp_clock_end (donner_mca_clock)
        endif !(do_donner_mca) 

!-----------------------------------------------------------------------
!    if a realizability constraint is to be placed on total cloud fraction,
!    define the area available for clouds from other schemes (usually
!    uw shallow) after the donner cloud area has been accounted for.
!    Note also that if the entire area at any level is taken up by donner 
!    clouds, then uw clouds will not be allowed in the box 
!    ( set conv_calc_completed = T).
!-----------------------------------------------------------------------
        if (limit_conv_cloud_frac) then
          ltemp = ANY(C2ls_mp%donner_humidity_area(:,:,:) >= 0.999,   &
                                                                  dim = 3)
          where (ltemp(:,:)) conv_calc_completed(:,:) = .true.
          available_cf_for_uw = MAX(0.999 -    &
                                 C2ls_mp%donner_humidity_area(:,:,:), 0.0)
        endif

!-----------------------------------------------------------------------
!    update the largescale cloud fields and their total tendencies from
!    physics  with the tendencies resulting from the donner deep 
!    convection scheme.
!-----------------------------------------------------------------------
        if (doing_prog_clouds) then
          Input_mp%tracer(:,:,:,nql) = qlin + delta_ql
          Input_mp%tracer(:,:,:,nqi) = qiin + delta_qi
          Input_mp%tracer(:,:,:,nqa) = qain + delta_qa
          Output_mp%rdt(:,:,:,nql) = Output_mp%rdt(:,:,:,nql) +  &
                                                   delta_ql*dtinv
          Output_mp%rdt(:,:,:,nqi) = Output_mp%rdt(:,:,:,nqi) +   &
                                                   delta_qi*dtinv
          Output_mp%rdt(:,:,:,nqa) = Output_mp%rdt(:,:,:,nqa) +  &
                                                   delta_qa*dtinv


!------------------------------------------------------------------------
!    calculate the amount of ice particles detrained from the donner
!    convective clouds. Modify the ice particle number and ice particle
!    number tendency from physics to account for this detrainment.  output
!    a diagnostic if desired.
!------------------------------------------------------------------------
          if (do_ice_num .and. detrain_ice_num) then
            CALL detr_ice_num (Input_mp%tin, delta_qi, delta_qni)   
            Input_mp%tracer(:,:,:,nqni) =  nilin  + delta_qni 
            Output_mp%rdt(:,:,:,nqni) = Output_mp%rdt(:,:,:,nqni) +    &
                                                    delta_qni*dtinv
            used = send_data (id_qnidt_deep_donner, delta_qni*dtinv,   &
                                                          Time, is, js, 1)
          endif  

!-------------------------------------------------------------------------
!    detrain liquid droplets if desired. the original code had a bug which
!    may be preserved for test purposes with the remain_detrain_bug nml
!    variable. assume 10 micron mean volume radius for detrained droplets. 
!    Modify the particle number and particle number tendency from physics 
!    to account for this detrainment. output a diagnostic if desired.
!-------------------------------------------------------------------------
          if (do_liq_num .and. detrain_liq_num) then
            if (remain_detrain_bug ) then
              delta_qn =  delta_ql/1000.*3./(4.*3.14*10.e-15)
            else
              delta_qn =  delta_ql/1000.*3./(4.*3.14e-15)
            endif !(        remain_detrain_bug )
            Input_mp%tracer(:,:,:,nqn) =  nllin + delta_qn 
            Output_mp%rdt(:,:,:,nqn) = Output_mp%rdt(:,:,:,nqn) +   &
                                                  delta_qn*dtinv
            used = send_data (id_qndt_deep_donner,   &
                                        delta_qn*dtinv, Time, is, js, 1)
          endif
        endif  ! doing_prog_clouds

!---------------------------------------------------------------------
!    update the values of temperature and vapor specific humidity to
!    include the effects of deep convection. if mca was included in the
!    donner deep scheme, then this update has already been done.
!---------------------------------------------------------------------
        if (keep_icenum_detrain_bug) then
          if (.not. do_donner_mca) then
            Input_mp%tin(:,:,:) = Input_mp%tin(:,:,:) + delta_temp(:,:,:)
            Input_mp%qin(:,:,:) = Input_mp%qin(:,:,:) + delta_q(:,:,:)
          endif
        endif
      endif !(do_donner_deep)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                0. UW SHALLOW CONVECTION PARAMETERIZATION
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!------------------------------------------------------------------------
!    if uw shallow convection is to be done after donner deep, do it here.
!    start shallowcu clock.
!------------------------------------------------------------------------
      if (do_donner_before_uw) then
        if (do_uw_conv) then
          call mpp_clock_begin (shallowcu_clock)

!---------------------------------------------------------------------
!    be sure all optional arguments associated with the uw_conv param-
!    eterization are present.
!---------------------------------------------------------------------
          if (i_shallow /= 0) then
          else
            call error_mesg ('convection_driver_mod',  &
                        'improper arguments for uw_shallow clouds', FATAL)
          endif

!---------------------------------------------------------------------
!    update tracer fields with tendencies due to donner convection and 
!    wet deposition by donner deep precipitation if these updated fields
!    are what is to be seen by uw convection.
!---------------------------------------------------------------------
          if (use_updated_profiles_for_uw) then 
            do n=1,nt  
              if (.not. cloud_tracer(n)) then
                  Input_mp%tracer(:,:,:,n) = tracer_orig(:,:,:,n) +   &
                        (Output_mp%rdt(:,:,:,n) - rdt_init(:,:,:,n)) *dt
              endif
            end do

!---------------------------------------------------------------------
!    call the uw_conv wrapper routine.
!---------------------------------------------------------------------
            call moistproc_uw_conv  &
                   (Time, is, ie, js, je, dt, Input_mp%tin, Input_mp%qin, &
                    Input_mp%uin, Input_mp%vin, Input_mp%tracer,    &
                    Input_mp%pfull, Input_mp%phalf, Input_mp%zfull,   &
                    Input_mp%zhalf, Input_mp%omega, Input_mp%pblht,   &
                    Input_mp%ustar, Input_mp%bstar, Input_mp%qstar,   &
                    shflx, lhflx, Input_mp%land, Input_mp%coldT, Aerosol, &
                    Surf_diff%tdt_rad, Surf_diff%tdt_dyn,   &
                    Surf_diff%qdt_dyn, Surf_diff%dgz_dyn,   &
                    Surf_diff%ddp_dyn, tdt_dif, qdt_dif, &
                    Phys_mp_exch%hmint, Input_mp%lat, Input_mp%lon,   &
                    Input_mp%cush, Input_mp%cbmf, Phys_mp_exch%cgust,   &
                    Phys_mp_exch%tke, Phys_mp_exch%pblhto, &
                    Phys_mp_exch%rkmo, Phys_mp_exch%taudpo, &
                    Phys_mp_exch%exist_shconv, Phys_mp_exch%exist_dpconv, &
                    Phys_mp_exch%pblht_prev, Phys_mp_exch%hlsrc_prev, &
                    Phys_mp_exch%qtsrc_prev, Phys_mp_exch%cape_prev, &
                    Phys_mp_exch%cin_prev, Phys_mp_exch%tke_prev, &
                    cmf, conv_calc_completed, available_cf_for_uw,   &
                    Output_mp%tdt, Output_mp%rdt(:,:,:,1), Output_mp%udt, &
                    Output_mp%vdt, Output_mp%rdt, Tend_mp%ttnd_conv,  &
                    Tend_mp%qtnd_conv, Output_mp%lprec, Output_mp%fprec, &
                    precip, Removal_mp%liq_precflx,   &
                    Removal_mp%ice_precflx, rain_uw, snow_uw, ttnd_uw,   &
                    qtnd_uw, utnd_uw, vtnd_uw, qtruw, qltnd_uw, qitnd_uw, &
                    qatnd_uw, qntnd_uw, qnitnd_uw, doing_prog_clouds,  &
                    do_limit_uw, do_liq_num, num_prog_tracers,  &
                    Removal_mp_control%tracers_in_uw,   &
                    Removal_mp_control%num_uw_tracers,   &
                    Moist_clouds_block%cloud_data(i_shallow),  &
                    Removal_mp%uw_wetdep, do_ice_num, detrain_ice_num)

!---------------------------------------------------------------------
!    if not using updated profiles for uw shallow, call the uw_conv 
!    wrapper routine with the original input fields of t,q, and tracer.
!---------------------------------------------------------------------
          else ! ( use_updated_profiles_for_uw)
            call moistproc_uw_conv   &
                   (Time, is, ie, js, je, dt, tin_orig, qin_orig, &
                    Input_mp%uin, Input_mp%vin, tracer_orig,    &
                    Input_mp%pfull, Input_mp%phalf, Input_mp%zfull,  &
                    Input_mp%zhalf, Input_mp%omega, Input_mp%pblht,  &
                    Input_mp%ustar, Input_mp%bstar, Input_mp%qstar,  &
                    shflx, lhflx, Input_mp%land, Input_mp%coldT, Aerosol, &
                    Surf_diff%tdt_rad, Surf_diff%tdt_dyn,   &
                    Surf_diff%qdt_dyn, Surf_diff%dgz_dyn,   &
                    Surf_diff%ddp_dyn, tdt_dif, qdt_dif, &
                    Phys_mp_exch%hmint, Input_mp%lat, Input_mp%lon,   &
                    Input_mp%cush, Input_mp%cbmf, Phys_mp_exch%cgust,   &
                    Phys_mp_exch%tke, Phys_mp_exch%pblhto, &
                    Phys_mp_exch%rkmo, Phys_mp_exch%taudpo, &
                    Phys_mp_exch%exist_shconv, Phys_mp_exch%exist_dpconv, &
                    Phys_mp_exch%pblht_prev, Phys_mp_exch%hlsrc_prev, &
                    Phys_mp_exch%qtsrc_prev, Phys_mp_exch%cape_prev, &
                    Phys_mp_exch%cin_prev, Phys_mp_exch%tke_prev, &
                    cmf, conv_calc_completed, available_cf_for_uw,   &
                    Output_mp%tdt, Output_mp%rdt(:,:,:,1), Output_mp%udt, &
                    Output_mp%vdt, Output_mp%rdt, Tend_mp%ttnd_conv,   &
                    Tend_mp%qtnd_conv, Output_mp%lprec, Output_mp%fprec, &
                    precip, Removal_mp%liq_precflx,   &
                    Removal_mp%ice_precflx, rain_uw, snow_uw, ttnd_uw,   &
                    qtnd_uw, utnd_uw, vtnd_uw, qtruw, qltnd_uw, qitnd_uw, &
                    qatnd_uw, qntnd_uw, qnitnd_uw, doing_prog_clouds,   &
                    do_limit_uw, do_liq_num, num_prog_tracers,  &
                    Removal_mp_control%tracers_in_uw,    &
                    Removal_mp_control%num_uw_tracers,   &
                    Moist_clouds_block%cloud_data(i_shallow),  &
                    Removal_mp%uw_wetdep, do_ice_num, detrain_ice_num)
          endif ! (use_updated_profiles_for_uw)
          call mpp_clock_end (shallowcu_clock)
        endif !(do_uw_conv)
      endif !(do_donner_before_uw)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                B. MOIST CONVECTIVE ADJUSTMENT             
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!----------------------------------------------------------------------
!    if moist convective adjustment is active,activate its clock and call 
!    its wrapper routine.
!----------------------------------------------------------------------
      if (do_mca) then
        call mpp_clock_begin (mca_clock)
        call moistproc_mca     &
             (Time, is, js, Input_mp%tin, Input_mp%qin, Input_mp%tracer, &
              Input_mp%pfull, Input_mp%phalf, Input_mp%coldT, dtinv, &
              Output_mp%tdt, Output_mp%rdt(:,:,:,1), Output_mp%rdt,  &
              Tend_mp%q_tnd, Tend_mp%ttnd_conv, Tend_mp%qtnd_conv,      &
              Output_mp%lprec, Output_mp%fprec, doing_prog_clouds,  &
              num_prog_tracers, Removal_mp_control%tracers_in_mca,   &
              Removal_mp_control%num_mca_tracers )
        call mpp_clock_end (mca_clock)
      endif 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!           X. BETTS-MILLER CONVECTION SCHEME 
!			
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!----------------------------------------------------------------------
!    if one of the betts-miller convection schemes is active, activate 
!    its clock and call the appropriate subroutine.
!----------------------------------------------------------------------
      if ( any((/do_bm, do_bmmass, do_bmomp/)) ) then
        call mpp_clock_begin (bm_clock)

! betts-miller cumulus param scheme
        if (do_bm) then 
          call betts_miller     &
               (dt, Input_mp%tin, Input_mp%qin, Input_mp%pfull, &
                Input_mp%phalf, Input_mp%coldT, rain, snow, Tend_mp%ttnd, &
                Tend_mp%qtnd, q_ref, bmflag, klzbs, cape, cin, t_ref, &
                invtaubmt, invtaubmq, mask=mask)
        endif

! betts-miller-style massflux cumulus param scheme
        if (do_bmmass) then  
          call bm_massflux    &
               (dt, Input_mp%tin, Input_mp%qin, Input_mp%pfull,   &
                Input_mp%phalf, Input_mp%coldT, rain, snow, Tend_mp%ttnd, &
                Tend_mp%qtnd, q_ref, bmflag, klzbs, t_ref, massflux,  &
                mask=mask)
        endif

! olivier's betts-miller cumulus param scheme
        if (do_bmomp) then 
          call bm_omp    &
               (dt, Input_mp%tin, Input_mp%qin, Input_mp%pfull,  &
                Input_mp%phalf, Input_mp%coldT, rain, snow, Tend_mp%ttnd, &
                Tend_mp%qtnd, q_ref, bmflag, klzbs, t_ref, mask=mask)
        endif

!----------------------------------------------------------------------
!    update input values and compute tendency.
!----------------------------------------------------------------------
        Input_mp%tin = Input_mp%tin + Tend_mp%ttnd
        Input_mp%qin = Input_mp%qin + Tend_mp%qtnd
        Tend_mp%ttnd = Tend_mp%ttnd*dtinv
        Tend_mp%qtnd = Tend_mp%qtnd*dtinv
        rain = rain*dtinv
        snow = snow*dtinv

!----------------------------------------------------------------------
!    add tendencies and generated precipitation to arrays accumulating 
!    tendencies from physics.
!----------------------------------------------------------------------
        Output_mp%tdt = Output_mp%tdt + Tend_mp%ttnd
        Output_mp%rdt(:,:,:,1) = Output_mp%rdt(:,:,:,1) + Tend_mp%qtnd
        Output_mp%lprec = Output_mp%lprec + rain
        Output_mp%fprec = Output_mp%fprec + snow
        precip = precip + rain + snow
                                                                         
!-------------------------------------------------------------------------
!     compute rh clouds if they are active with betts-miller. first 
!     calculate the relative humidity, then pass it to rh_clouds_mod to be
!     stored till needed.
!-------------------------------------------------------------------------
        if (do_rh_clouds) then
          call rh_calc   &
             (Input_mp%pfull, Input_mp%tin, Input_mp%qin, RH, do_simple)
          call rh_clouds_sum (is, js, RH) 
        end if

!-----------------------------------------------------------------------
!     save desired betts-miller diagnostics.
!-----------------------------------------------------------------------
        used = send_data (id_tref, t_ref, Time, is, js, 1 )
        used = send_data (id_qref, q_ref, Time, is, js, 1 )
        used = send_data (id_bmflag, bmflag, Time, is, js)
        used = send_data (id_klzbs, klzbs, Time, is, js)
        used = send_data (id_invtaubmt, invtaubmt, Time, is, js)
        used = send_data (id_invtaubmq, invtaubmq, Time, is, js)
        used = send_data (id_massflux, massflux, Time, is, js, 1)
        call mpp_clock_end (bm_clock)
      endif ! if ( any((/do_bm, do_bmmass, do_bmomp/)) )


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!           C. RELAXED ARAKAWA-SCHUBERT PARAMETERIZATION
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!    execute relaxed arakawa/schubert cumulus parameterization scheme,
!    if desired.
!-----------------------------------------------------------------------
      if (do_ras) then
        call mpp_clock_begin (ras_clock)

!------------------------------------------------------------------------
!    call the wrapper routine for the ras parameterization.
!------------------------------------------------------------------------
        call moistproc_ras   &
            (Time, is, js, dt, Input_mp%coldT, Input_mp%tin, Input_mp%qin,&
             Input_mp%uin, Input_mp%vin, Input_mp%tracer, Input_mp%pfull, &
             Input_mp%phalf, Input_mp%zhalf, Output_mp%tdt,   &
             Output_mp%rdt(:,:,:,1), Output_mp%udt, Output_mp%vdt,  &
             Output_mp%rdt, Tend_mp%q_tnd, Tend_mp%ttnd, Tend_mp%qtnd, &
             Tend_mp%ttnd_conv, Tend_mp%qtnd_conv, mc, det0,  &
             Output_mp%lprec, Output_mp%fprec, rain_ras, snow_ras,  &
             rain3d, snow3d, Aerosol, doing_prog_clouds, do_liq_num,  &
             num_prog_tracers, Removal_mp_control%tracers_in_ras,   &
             Removal_mp_control%num_ras_tracers,             &
             do_ice_num, detrain_ice_num)
        call mpp_clock_end (ras_clock)
      else

!---------------------------------------------------------------------
!    if ras_mod is not activated, set the ras mass flux and precip fields
!    to 0.
!---------------------------------------------------------------------
        mc   = 0.0
        rain_ras = 0.0
        snow_ras = 0.0
      endif  ! (       do_ras)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!      END OF INDIVIDUAL CONVECTIVE SCHEMES
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!---------------------------------------------------------------------
!    now that all potential convection schemes have been processed, 
!    calculate cumulus momentum transport, if desired. 
!---------------------------------------------------------------------
      if (do_cmt) then

!----------------------------------------------------------------------
!    activate the cmt clock. initialize the output field (not needed, was
!    done when allocated).
!----------------------------------------------------------------------
        call mpp_clock_begin (cmt_clock)

!----------------------------------------------------------------------
!    if doing nonlocal cmt, call cu_mo_trans for each convective scheme
!    separately.
!----------------------------------------------------------------------
        if (.not. doing_diffusive) then
          if (cmt_uses_ras) then
            call moistproc_cmt    &
               ( Time, is, js, Input_mp%tin, Input_mp%uin, Input_mp%vin, &
                 Input_mp%tracer, Input_mp%pfull, Input_mp%phalf, &
                 Input_mp%zfull, Input_mp%zhalf, Input_mp%pmass,  &
                 Output_mp%tdt, Output_mp%udt, Output_mp%vdt,   &
                 Output_mp%rdt, Tend_mp%ttnd_conv, dt, mc, det0,  &
                 Output_mp%diff_cu_mo, num_prog_tracers)
          endif !(cmt_uses_ras)
 
          if (cmt_uses_donner) then
            call moistproc_cmt    &
               ( Time, is, js, Input_mp%tin, Input_mp%uin, Input_mp%vin, &
                 Input_mp%tracer, Input_mp%pfull, Input_mp%phalf, &
                 Input_mp%zfull, Input_mp%zhalf, Input_mp%pmass,   &
                 Output_mp%tdt, Output_mp%udt, Output_mp%vdt,   &
                 Output_mp%rdt, Tend_mp%ttnd_conv, dt, m_cellup,  &
                 M_cdet_donner, Output_mp%diff_cu_mo,  num_prog_tracers)
          endif
 
          if (cmt_uses_uw) then

!----------------------------------------------------------------------
!    CURRENTLY no detrained mass flux is provided from uw_conv; should only
!    use with 'diffusive' cmt scheme, not the non-local. (attempt to
!    use non-local will cause FATAL in _init routine.)
!----------------------------------------------------------------------
          endif

        else ! (we are doing_diffusive)

!-----------------------------------------------------------------------
!    if using diffusive cmt, call cu_mo_trans once with combined mass
!    fluxes from all desired convective schemes.
!-----------------------------------------------------------------------
          mc_cmt = 0.
          if (cmt_uses_ras) then
            mc_cmt = mc_cmt + mc
          endif
          if (cmt_uses_donner) then
            mc_cmt = mc_cmt + m_cellup 
          endif
          if (cmt_uses_uw) then
            do k=2,kx
              mc_cmt(:,:,k) = mc_cmt(:,:,k) + cmf(:,:,k-1)
            end do
          endif

          call moistproc_cmt      &
               ( Time, is, js, Input_mp%tin, Input_mp%uin, Input_mp%vin, &
                 Input_mp%tracer, Input_mp%pfull, Input_mp%phalf, &
                 Input_mp%zfull, Input_mp%zhalf, Input_mp%pmass,   &
                 Output_mp%tdt, Output_mp%udt, Output_mp%vdt,   &
                 Output_mp%rdt, Tend_mp%ttnd_conv, dt, mc_cmt, det_cmt,  &
                 Output_mp%diff_cu_mo, num_prog_tracers)
        endif ! (.not. doing_diffusive)

        call mpp_clock_end (cmt_clock)
      endif  !(do_cmt)

!------------------------------------------------------------------------
!    initialize fields needed for call to wet deposition routine.
!------------------------------------------------------------------------
      f_snow_berg = 0.
      C2ls_mp%wet_data = 0.0
      C2ls_mp%cloud_frac = 0.1
      C2ls_mp%cloud_wet = 1.e-3
      Tend_mp%qtnd_wet(:,:,:) = Tend_mp%qtnd(:,:,:)
      if (doing_prog_clouds) then
        Tend_mp%qtnd_wet(:,:,:) = Tend_mp%qtnd_wet(:,:,:) +   &
                                  Tend_mp%q_tnd(:,:,:,nql) +    &
                                  Tend_mp%q_tnd(:,:,:,nqi)
      end if

!---------------------------------------------------------------------
!    for each tracer for which wet deposition has been requested, call
!    subrouitine wet_deposition to calculate the tracer tendency due to 
!    wet deposition (wetdeptnd) caused by the convectively generated 
!    precipitation (rain, snow). 
!---------------------------------------------------------------------
      do n=1, nt
        if (.not. cloud_tracer(n)) then
          Tend_mp%wetdeptnd(:,:,:) = 0.0
          call wet_deposition   &
              (n, Input_mp%t, Input_mp%pfull, Input_mp%phalf,   &
               Input_mp%zfull, Input_mp%zhalf, rain_ras, snow_ras, &
               Tend_mp%qtnd_wet, C2ls_mp%cloud_wet, C2ls_mp%cloud_frac,  &
               f_snow_berg, rain3d, snow3d, Input_mp%tracer(:,:,:,n),   &
               Tend_mp%wetdeptnd, Time, 'convect', is, js, dt )

!-----------------------------------------------------------------------
!    add this tendency to the tracer tendency due to all physics (rdt). 
!    save it also in an array which will be combined with any wet 
!    deposition resulting from large-scale precip producing the total wet 
!    deposition for the tracer (wet_data).
!---------------------------------------------------------------------
          Output_mp%rdt (:,:,:,n) = Output_mp%rdt(:,:,:,n) -   &
                                                 Tend_mp%wetdeptnd(:,:,:)
          C2ls_mp%wet_data(:,:,:,n) = Tend_mp%wetdeptnd(:,:,:)
        endif
      end do

!------------------------------------------------------------------------
!    define total convective mass flux from all sources, at both full
!    levels and at half levels.
!------------------------------------------------------------------------
      C2ls_mp%mc_full(:,:,1)=0.; 
      C2ls_mp%mc_half(:,:,1)=0.; 
      do k=2,kx   
        C2ls_mp%mc_full(:,:,k) = 0.5*(mc(:,:,k) + mc(:,:,k+1)) +   &
                                 0.5*(cmf(:,:,k)+cmf(:,:,k-1)) +   &
                                 mc_donner(:,:,k)
      end do
      do k=2,kx+1   
        C2ls_mp%mc_half(:,:,k) = mc(:,:,k) +    &
                                 cmf(:,:,k-1)+   &
                                 mc_donner_half(:,:,k)
      end do

!---------------------------------------------------------------------
!    output diagnostics:
!    total cumulus mass flux on full levels,
!    total cumulus mass flux on half levels,
!    total cumulus mass flux on half levels (CMOR standard).
!---------------------------------------------------------------------
      used = send_data (id_mc_full, C2ls_mp%mc_full, Time, is, js, 1)
      used = send_data (id_mc_half, C2ls_mp%mc_half, Time, is, js, 1)
      used = send_cmip_data_3d (ID_mc, C2ls_mp%mc_half, Time, is, js, 1)!, rmask=mask)
!      used = send_data (id_mc     , C2ls_mp%mc_half, Time, is, js, 1)

!---------------------------------------------------------------------
!    total convective updraft mass flux (uw + donner cell up + 
!    donner meso up)
!---------------------------------------------------------------------
      if (id_mc_conv_up > 0 ) &
        used = send_data (id_mc_conv_up, cmf(:,:,:) + &
                                  mc_donner_up(:,:,:), Time, is, js, 1 )

!------------------------------------------------------------------------ 
!    define convective cloud base and cloud top. these are needed if 
!    diagnostics defining the temporal and spatial location of convection 
!    are desired or if tracer "no" is present, so that the nox tendency 
!    due to lightning may be calculated.
!------------------------------------------------------------------------ 
      if (get_tracer_index(MODEL_ATMOS,'no') .ne. NO_TRACER &
           .or. id_conv_freq > 0 &
           .or. id_ci > 0 &
           .or. id_conv_cld_base > 0 &
           .or. id_ccb           > 0 &
           .or. id_cct           > 0 &
           .or. id_conv_cld_top > 0 ) then

        cldbot = 0
        cldtop = 0
        do j = 1,jx
          do i = 1,ix
            do k = 1,kx
              if (C2ls_mp%mc_full(i,j,k) /= 0.0 ) then
                cldtop(i,j) = k
                exit
              endif
            enddo
            do k = size(Input_mp%r,3),1,-1
              if (C2ls_mp%mc_full(i,j,k) /= 0.0 ) then
                cldbot(i,j) = k
                exit
              endif
            enddo
          enddo
        enddo
      end if

!------------------------------------------------------------------------
!    output diagnostics related to convective cloud base and cloud top.
!    both FMS-standard and CMOR-standard output variables are currently
!    present.
!------------------------------------------------------------------------
      if ( id_conv_cld_base > 0 ) then
        temp_2d = missing_value
        do j = 1,jx
          do i = 1,ix
            if ( cldbot(i,j) > 0 ) temp_2d(i,j) =    &
                                         Input_mp%pfull(i,j,cldbot(i,j))
          end do
        end do
        used = send_data(id_conv_cld_base, temp_2d, Time, is_in=is,   &
                                             js_in=js,  mask = cldbot > 0)
      end if

      if ( id_ccb > 0 ) then
        temp_2d = CMOR_MISSING_VALUE
        do j = 1,jx  
          do i = 1,ix
            if ( cldbot(i,j) > 0 ) temp_2d(i,j) =    &
                                         Input_mp%pfull(i,j,cldbot(i,j))
          end do
        end do
        used = send_data(id_ccb, temp_2d, Time, is_in=is,   &
                                            js_in=js,  mask = cldbot > 0)
      end if

      if ( id_conv_cld_top > 0 ) then
        temp_2d = missing_value
        do j = 1,jx
          do i = 1,ix
            if ( cldtop(i,j) > 0 ) temp_2d(i,j) =   &
                                         Input_mp%pfull(i,j,cldtop(i,j))
          end do
        end do
        used = send_data(id_conv_cld_top, temp_2d, Time, is_in=is, &
                                           js_in=js,  mask = cldtop > 0)
      end if

      if ( id_cct > 0 ) then
        temp_2d = CMOR_MISSING_VALUE
        do j = 1,jx
          do i = 1,ix
            if ( cldtop(i,j) > 0 ) temp_2d(i,j) =    &
                                         Input_mp%pfull(i,j,cldtop(i,j))
          end do
        end do
        used = send_data(id_cct, temp_2d, Time, is_in=is, &
                                           js_in=js,  mask = cldtop > 0)
      end if

!-----------------------------------------------------------------------
!    calculate NOx tendency from lightning and add it to the tendency
!    field.
!-----------------------------------------------------------------------
      if ( get_tracer_index(MODEL_ATMOS,'no') .ne. NO_TRACER ) then
        call moz_hook       &
              (cldtop, cldbot, Input_mp%land, Input_mp%zfull,   &
               Input_mp%zhalf, Input_mp%t, prod_no, Input_mp%area,   &
               Input_mp%lat, Time, is, js)
!  conc_air
        Output_mp%rdt(:,:,:,get_tracer_index(MODEL_ATMOS,'no')) =  &
              Output_mp%rdt(:,:,:,get_tracer_index(MODEL_ATMOS,'no')) + &
                  prod_no* ((boltz * Input_mp%t) / (10. * Input_mp%pfull)) 
        used = send_data(id_prod_no, prod_no, Time, is_in=is, js_in=js)
      endif

!-----------------------------------------------------------------------
!    define the total precipitation rate (precip).
!-----------------------------------------------------------------------
      precip = Output_mp%lprec + Output_mp%fprec

!---------------------------------------------------------------------
!    apply changes resulting from uw_conv. scale if necessary to maintain
!    realizability. output a diagnostic for the scale.
!---------------------------------------------------------------------
      if (do_uw_conv) then
        if (do_limit_uw) then
          call moistproc_scale_uw    &
              (is, ie, js, je, dt, Input_mp%qin, Input_mp%tracer,   &
               Output_mp%tdt, Output_mp%rdt(:,:,:,1), Output_mp%udt,   &
               Output_mp%vdt, Output_mp%rdt, Tend_mp%ttnd_conv,   &
               Tend_mp%qtnd_conv, Output_mp%lprec, Output_mp%fprec,   &
               precip, qtruw, rain_uw, snow_uw, ttnd_uw, qtnd_uw,   &
               utnd_uw, vtnd_uw, qltnd_uw, qitnd_uw, qatnd_uw, qntnd_uw, &
               qnitnd_uw, doing_prog_clouds, do_liq_num, num_prog_tracers,&
               Removal_mp_control%tracers_in_uw, scale, do_ice_num)
        else  
          scale = 1.0
        endif 
        used = send_data (id_scale_uw, scale, Time, is, js )

!-------------------------------------------------------------------------
!    update input fields with changes from uw_conv
!-------------------------------------------------------------------------
        Input_mp%tin = Input_mp%tin + ttnd_uw*dt
        Input_mp%qin = Input_mp%qin + qtnd_uw*dt
        Input_mp%uin = Input_mp%uin + utnd_uw*dt
        Input_mp%vin = Input_mp%vin + vtnd_uw*dt
        Input_mp%tracer(:,:,:,nql) = Input_mp%tracer(:,:,:,nql) +    &
                                                       qltnd_uw*dt
        Input_mp%tracer(:,:,:,nqi) = Input_mp%tracer(:,:,:,nqi) +     &
                                                       qitnd_uw*dt
        Input_mp%tracer(:,:,:,nqa) = Input_mp%tracer(:,:,:,nqa) +    &
                                                       qatnd_uw*dt
        if (do_liq_num) then
          Input_mp%tracer(:,:,:,nqn) = Input_mp%tracer(:,:,:,nqn) +   &
                                                       qntnd_uw*dt
        endif
        if (do_ice_num) then
          Input_mp%tracer(:,:,:,nqni) = Input_mp%tracer(:,:,:,nqni) +   &
                                                       qnitnd_uw*dt
        endif
      endif !(uw_conv)

!-----------------------------------------------------------------------
! Note: the following two blocks of code were moved from before to after 
! the call to moistproc_scale_uw in order to account for UW convection 
! contributions to 'precip'. This will change answers if UW is active
! and do_gust_cv = .true. or entrain_nml:convect_shutoff = .true. (cjg)
!
!    calculate convective gustiness, if desired.
!-----------------------------------------------------------------------
      if (do_gust_cv) then
        where((precip) > 0.0)
          Output_mp%gust_cv = gustmax*sqrt( precip/(gustconst + precip) )
        end where
      end if

      if (do_gust_cv_new) then
        Output_mp%gust_cv = sqrt(Phys_mp_exch%cgust)
      end if

!---------------------------------------------------------------------
!    save a field indicating whether or not convection has occurred
!    within the column.
!---------------------------------------------------------------------
      where (precip > 0.) Output_mp%convect = .true.

!---------------------------------------------------------------------
!    update tracer fields with tendencies due to convection and wet 
!    deposition by convective precipitation.
!---------------------------------------------------------------------
      do n=1,nt                       
        if (.not. cloud_tracer(n)) then
          Input_mp%tracer(:,:,:,n) = tracer_orig(:,:,:,n) +   &
                          (Output_mp%rdt(:,:,:,n) - rdt_init(:,:,:,n)) *dt
        endif
      end do

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                   CONVECTION DIAGNOSTICS      
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!    precip from the various convection schemes.
!-----------------------------------------------------------------------
      used = send_data (id_ras_precip, rain_ras + snow_ras, Time, is, js)
      used = send_data (id_don_precip, rain_don + snow_don +  &
                                       rain_donmca + snow_donmca, &
                                                             Time, is, js)
      used = send_data (id_uw_precip, rain_uw + snow_uw, Time, is, js)
      used = send_data (id_uw_snow, snow_uw, Time, is, js)

!-----------------------------------------------------------------------
!    uw_conv diagnostics: prognostic variable tendencies and enthalpy
!    and water column tendencies.
!-----------------------------------------------------------------------
      if (do_uw_conv) then
        used = send_data (id_tdt_uw, ttnd_uw, Time, is, js, 1)
        used = send_data (id_qdt_uw, qtnd_uw, Time, is, js, 1)
        used = send_data (id_qadt_uw, qatnd_uw, Time, is, js, 1)
        used = send_data (id_qldt_uw, qltnd_uw, Time, is, js, 1)
        used = send_data (id_qidt_uw, qitnd_uw, Time, is, js, 1)
        if (do_liq_num) then
          used = send_data (id_qndt_uw, qntnd_uw, Time, is, js, 1)
        endif
        if (do_ice_num) then
          used = send_data (id_qnidt_uw, qnitnd_uw, Time, is, js, 1)
        end if
      endif

      if (id_enth_uw_col > 0) then
        temp_2d = -HLV*rain_uw -HLS*snow_uw
        call column_diag (id_enth_uw_col, is, js, Time, ttnd_uw, CP_AIR, &
               qltnd_uw, -HLV, qitnd_uw, -HLS, Input_mp%pmass, temp_2d)
      endif

      if (id_wat_uw_col > 0) then
        temp_2d = rain_uw + snow_uw
        call column_diag(id_wat_uw_col, is, js, Time, qtnd_uw, 1.0,  &
               qltnd_uw, 1.0, qitnd_uw, 1.0, Input_mp%pmass, temp_2d)
      endif
        
!----------------------------------------------------------------------
!    convection scheme frequency diagnostics.
!----------------------------------------------------------------------
      if (id_ras_freq > 0) then
        ltemp = rain_ras > 0. .or. snow_ras > 0.0
        where (ltemp) 
          temp_2d = 1.
        elsewhere
          temp_2d = 0.
        end where
        used = send_data (id_ras_freq, temp_2d,Time, is, js)
      endif

      if (id_don_freq > 0) then
        ltemp = rain_don > 0. .or. snow_don > 0.0 .or. &
                rain_donmca > 0. .or. snow_donmca > 0.0
        where (ltemp) 
          temp_2d = 1.
        elsewhere
          temp_2d = 0.
        end where
        used = send_data (id_don_freq, temp_2d, Time, is, js)
      endif

      if (id_uw_freq > 0) then
        ltemp = rain_uw > 0. .or. snow_uw > 0.0
        where (ltemp) 
          temp_2d = 1.
        elsewhere
          temp_2d = 0.
        end where
        used = send_data (id_uw_freq, temp_2d, Time, is, js)
      endif

!---------------------------------------------------------------------
!    temperature change due to dry and moist convection:
!---------------------------------------------------------------------
      used = send_data (id_tdt_conv, Tend_mp%ttnd_conv, Time, is, js, 1)
      used = send_cmip_data_3d (ID_tntc, Tend_mp%ttnd_conv, Time, is, js, 1)!, rmask=mask)

!---------------------------------------------------------------------
!    vapor specific humidity change due to convection:
!---------------------------------------------------------------------
      used = send_data (id_qdt_conv, Tend_mp%qtnd_conv, Time, is, js, 1)
      used = send_cmip_data_3d (ID_tnhusc, Tend_mp%qtnd_conv, Time, is, js, 1)!, rmask=mask)

!---------------------------------------------------------------------
!    total precipitation due to convection (both FMS and CMOR standards):
!---------------------------------------------------------------------
      used = send_data (id_prec_conv, precip, Time, is, js)
      used = send_data (id_prc, precip, Time, is, js)
      if (id_prc_g > 0) call buffer_global_diag (id_prc_g, precip(:,:), Time, is, js)

!---------------------------------------------------------------------
!    frozen precipitation (snow) due to convection:
!---------------------------------------------------------------------
      used = send_data (id_snow_conv, Output_mp%fprec, Time, is, js)
      used = send_data (id_prsnc, Output_mp%fprec, Time, is, js)
!---------------------------------------------------------------------
!    liquid precipitation (rain) due to convection:
!---------------------------------------------------------------------
      used = send_data (id_prrc, Output_mp%lprec, Time, is, js)

!---------------------------------------------------------------------
!    convective frequency (both FMS and CMOR standards).
!---------------------------------------------------------------------
      if (id_conv_freq > 0) then
        ltemp = precip > 0. .or. cldtop > 0
        where (ltemp)
          freq_count = 1.
        elsewhere
          freq_count = 0.
        end where
        used = send_data (id_conv_freq, freq_count, Time, is, js )
      endif

      if (id_ci > 0) then
        ltemp = precip > 0. .or. cldtop > 0
        where (ltemp)
          freq_count = 1.
        elsewhere
          freq_count = 0.
        end where
        used = send_data (id_ci, freq_count, Time, is, js )
      endif

!---------------------------------------------------------------------
!    diagnostic for 3D convective precip
!---------------------------------------------------------------------
      used = send_data ( id_conv_rain3d, rain3d, Time, is, js, 1 )
      used = send_data ( id_conv_snow3d, snow3d, Time, is, js, 1 )

!---------------------------------------------------------------------
!    surface wind gustiness due to convection:
!---------------------------------------------------------------------
      used = send_data (id_gust_conv, Output_mp%gust_cv, Time, is, js)

!---------------------------------------------------------------------
!    water vapor path tendency due to convection:
!---------------------------------------------------------------------
      if (id_q_conv_col > 0)   &
           call column_diag (id_q_conv_col, is, js, Time, &
                              Tend_mp%qtnd_conv, 1.0, Input_mp%pmass)
  
!---------------------------------------------------------------------
!    dry static energy tendency due to dry and moist convection:
!---------------------------------------------------------------------
      if (id_t_conv_col > 0)   &
            call column_diag (id_t_conv_col, is, js, Time, &
                               Tend_mp%ttnd_conv, CP_AIR, Input_mp%pmass)
   
!---------------------------------------------------------------------
!    define the total prognostic cloud liquid, ice, drop number, 
!    ice number and area tendencies due to convection.
!---------------------------------------------------------------------
      if (doing_prog_clouds) then
        Tend_mp%qldt_conv = Output_mp%rdt(:,:,:,nql) - rdt_init(:,:,:,nql)
        Tend_mp%qidt_conv = Output_mp%rdt(:,:,:,nqi) - rdt_init(:,:,:,nqi)
        if (do_liq_num) Tend_mp%qndt_conv =    &
                           Output_mp%rdt(:,:,:,nqn) - rdt_init(:,:,:,nqn)
        if (do_ice_num) Tend_mp%qnidt_conv =    &
                         Output_mp%rdt(:,:,:,nqni) - rdt_init(:,:,:,nqni)
        Tend_mp%qadt_conv = Output_mp%rdt(:,:,:,nqa) - rdt_init(:,:,:,nqa)

!---------------------------------------------------------------------
!    output diagnostics for cloud liquid tendency and liquid water path 
!    tendency due to convection.
!---------------------------------------------------------------------
        if (id_qldt_conv > 0 .or. id_ql_conv_col > 0) then
          used = send_data (id_qldt_conv, Tend_mp%qldt_conv,    &
                                                       Time, is, js, 1)
          if (id_ql_conv_col > 0)    &
                call column_diag (id_ql_conv_col, is, js, Time,   &
                                  Tend_mp%qldt_conv, 1.0, Input_mp%pmass)

        endif

!---------------------------------------------------------------------
!    output diagnostics for cloud drop number tendency and cloud drop 
!    number path tendency due to convection.
!---------------------------------------------------------------------
        if (do_liq_num) then
          if (id_qndt_conv > 0 .or. id_qn_conv_col > 0) then
            used = send_data (id_qndt_conv, Tend_mp%qndt_conv,    &
                                                       Time, is, js, 1)
            if (id_qn_conv_col > 0)     &
                call column_diag (id_qn_conv_col, is, js, Time,   &
                                  Tend_mp%qndt_conv, 1.0, Input_mp%pmass)
          endif
        endif

!---------------------------------------------------------------------
!    output diagnostics for cloud ice tendency and cloud ice water path 
!    tendency due to convection.
!---------------------------------------------------------------------
        if (id_qidt_conv > 0 .or. id_qi_conv_col > 0) then
          used = send_data (id_qidt_conv, Tend_mp%qidt_conv, Time,   &
                                                              is, js, 1)
          if (id_qi_conv_col > 0)    &
                call column_diag (id_qi_conv_col, is, js, Time,    &
                                  Tend_mp%qidt_conv, 1.0, Input_mp%pmass)
        endif        


!---------------------------------------------------------------------
!    output diagnostics for cloud ice number tendency and cloud ice number
!    path tendency due to convection.
!---------------------------------------------------------------------
        if (do_ice_num) then
          if (id_qnidt_conv > 0 .or. id_qni_conv_col > 0) then
            used = send_data (id_qnidt_conv, Tend_mp%qnidt_conv,    &
                                                         Time, is, js, 1)
            if (id_qni_conv_col > 0)   &
                 call column_diag (id_qni_conv_col, is, js, Time,   &
                                Tend_mp%qnidt_conv, 1.0, Input_mp%pmass)
          endif
        endif

!---------------------------------------------------------------------
!    output diagnostics for cloud area tendency and column integrated 
!    cloud mass tendency due to convection.
!---------------------------------------------------------------------
        if (id_qadt_conv > 0 .or.  id_qa_conv_col > 0 ) then
          used = send_data (id_qadt_conv, Tend_mp%qadt_conv,    &
                                                        Time, is, js, 1)
          if (id_qa_conv_col > 0)     &
               call column_diag (id_qa_conv_col, is, js, Time,    &
                                Tend_mp%qadt_conv, 1.0, Input_mp%pmass)
        endif
      endif !(doing_prog_clouds)
         
!---------------------------------------------------------------------
!    compute the column integrated enthalpy and total water tendencies 
!    due to convection parameterizations, if those diagnostics are desired.
!---------------------------------------------------------------------
      if (id_enth_conv_col > 0 .or. id_wat_conv_col > 0) then
        temp_3d1 = Output_mp%rdt(:,:,:,nql) - rdt_init(:,:,:,nql)
        temp_3d2 = Output_mp%rdt(:,:,:,nqi) - rdt_init(:,:,:,nqi)

        if (id_enth_conv_col > 0) then
          temp_2d = -HLV*precip -HLF*Output_mp%fprec
          call column_diag    &
             (id_enth_conv_col, is, js, Time, Tend_mp%ttnd_conv, CP_AIR, &
                  temp_3d1, -HLV, temp_3d2, -HLS, Input_mp%pmass, temp_2d)
        endif

        if (id_wat_conv_col > 0) then
          temp_2d = precip
          call column_diag   &
             (id_wat_conv_col, is, js, Time, Tend_mp%qtnd_conv, 1.0,   &
                   temp_3d1, 1.0, temp_3d2, 1.0, Input_mp%pmass, temp_2d)
        endif
      endif

!---------------------------------------------------------------------
!    compute the tracer tendencies due to convection for any tracers that
!    are to be transported by any convective parameterization.
!---------------------------------------------------------------------
!      do n=1,size(Output_mp%rdt,4)
      do n=1,nt
        if (Removal_mp_control%tracers_in_donner(n) .or.   &
            Removal_mp_control%tracers_in_ras(n) .or.  &
            Removal_mp_control%tracers_in_mca(n) .or.   &
            Removal_mp_control%tracers_in_uw(n))    then
 
!---------------------------------------------------------------------
!    output diagnostics for tracer tendency and column integrated 
!    tracer tendency due to convection.
!---------------------------------------------------------------------
          if (id_tracerdt_conv(n) > 0 .or.    &
                                 id_tracerdt_conv_col(n) > 0) then
            temp_3d1 = Output_mp%rdt(:,:,:,n) - rdt_init(:,:,:,n)
            used = send_data (id_tracerdt_conv(n), temp_3d1,    &
                                                       Time, is, js, 1 )

            if (id_tracerdt_conv_col(n) > 0) &
              call column_diag    &
                  (id_tracerdt_conv_col(n), is, js, Time, temp_3d1,   &
                                                      1.0, Input_mp%pmass)
          endif        
        endif
      end do

!----------------------------------------------------------------------
!    define the precip fluxes needed for input to the COSP simulator 
!    package.
!---------------------------------------------------------------------
      if (do_cosp) then

!---------------------------------------------------------------------
!    define precip fluxes from convective schemes at each layer 
!    interface.  (index 1 is model lid)
!---------------------------------------------------------------------
        do k=2, size(Input_mp%t,3)+1
          Removal_mp%liq_mesoh(:,:,k) = Removal_mp%liq_mesoh (:,:,k-1) + &
                                        Removal_mp%liq_meso (:,:,k-1)
          Removal_mp%frz_mesoh(:,:,k) = Removal_mp%frz_mesoh (:,:,k-1) + &
                                        Removal_mp%frz_meso (:,:,k-1)
          Removal_mp%liq_cellh(:,:,k) = Removal_mp%liq_cellh (:,:,k-1) + &
                                        Removal_mp%liq_cell (:,:,k-1)
          Removal_mp%frz_cellh(:,:,k) = Removal_mp%frz_cellh (:,:,k-1) + &
                                        Removal_mp%frz_cell (:,:,k-1)
          Removal_mp%ice_precflxh(:,:,k) =                  &
                                     Removal_mp%ice_precflxh(:,:,k-1) +  &
                                          Removal_mp%ice_precflx(:,:,k-1)
          Removal_mp%liq_precflxh(:,:,k) =        &
                              Removal_mp%liq_precflxh(:,:,k-1) +   &
                                          Removal_mp%liq_precflx(:,:,k-1)
          if (include_donmca_in_cosp) then
            Removal_mp%mca_liqh(:,:,k) = Removal_mp%mca_liqh (:,:,k-1) + &
                                              Removal_mp%mca_liq(:,:,k-1)
            Removal_mp%mca_frzh(:,:,k) = Removal_mp%mca_frzh (:,:,k-1) + &
                                              Removal_mp%mca_frz(:,:,k-1)
          endif
        end do

!--------------------------------------------------------------------
!    adjust precip fluxes to account for any negative values produced.
!    precip contribution is determined as the negative of the moisture
!    tendency, so at top of clouds a positive moisture tendency some-
!    times results in a negative precipitation contribution. 
!----------------------------------------------------------------------
        call prevent_neg_precip_fluxes (Removal_mp%liq_mesoh)
        call prevent_neg_precip_fluxes (Removal_mp%frz_mesoh)
        call prevent_neg_precip_fluxes (Removal_mp%liq_cellh)
        call prevent_neg_precip_fluxes (Removal_mp%frz_cellh)
        call prevent_neg_precip_fluxes (Removal_mp%ice_precflxh)
        call prevent_neg_precip_fluxes (Removal_mp%liq_precflxh)
        if (include_donmca_in_cosp) then
          call prevent_neg_precip_fluxes (Removal_mp%mca_liqh)
          call prevent_neg_precip_fluxes (Removal_mp%mca_frzh)
        endif
      endif ! do_cosp

!-----------------------------------------------------------------------
!    compute the grid box area taken up by convective clouds.  If CLUBB is
!    active, then a second call is made using slightly different inputs. 
!-----------------------------------------------------------------------
      if (do_clubb == 2) then

!---> h1g, 2017-01-31
        if( do_uw_conv .and. do_donner_deep ) then
          call compute_convective_area     &
                 (Input_mp%t, Input_mp%pfull, Input_mp%q, do_uw_conv,  &
                  do_donner_deep, C2ls_mp%donner_humidity_area,   &
                  C2ls_mp%donner_humidity_factor, conv_frac_max,   &
                  C2ls_mp%convective_humidity_ratio_clubb, &
                  C2ls_mp%conv_frac_clubb,   &
                  shallow_cloud_area=  &
                       Moist_clouds_block%cloud_data(i_shallow)%cloud_area, &
                  cell_cld_frac =    &
                       Moist_clouds_block%cloud_data(i_cell)%cloud_area)
        else if (do_donner_deep) then
          call compute_convective_area     &
                 (Input_mp%t, Input_mp%pfull, Input_mp%q, do_uw_conv,  &
                  do_donner_deep, C2ls_mp%donner_humidity_area,   &
                  C2ls_mp%donner_humidity_factor, conv_frac_max,   &
                  C2ls_mp%convective_humidity_ratio_clubb, &
                  C2ls_mp%conv_frac_clubb,   &
                  cell_cld_frac =    &
                       Moist_clouds_block%cloud_data(i_cell)%cloud_area)

        else if (do_uw_conv) then
          call compute_convective_area     &
                 (Input_mp%t, Input_mp%pfull, Input_mp%q, do_uw_conv,  &
                  do_donner_deep, C2ls_mp%donner_humidity_area,   &
                  C2ls_mp%donner_humidity_factor, conv_frac_max,   &
                  C2ls_mp%convective_humidity_ratio_clubb, &
                  C2ls_mp%conv_frac_clubb,   &
                  shallow_cloud_area=  &
                       Moist_clouds_block%cloud_data(i_shallow)%cloud_area )

        else
           call compute_convective_area     &
                 (Input_mp%t, Input_mp%pfull, Input_mp%q, do_uw_conv,  &
                  do_donner_deep, C2ls_mp%donner_humidity_area,   &
                  C2ls_mp%donner_humidity_factor, conv_frac_max,   &
                  C2ls_mp%convective_humidity_ratio_clubb, &
                  C2ls_mp%conv_frac_clubb)

        endif
!<--- h1g, 2017-01-31

      endif
 
      if (.not. do_lsc) then

!---> h1g, 2017-01-31
        if( do_uw_conv .and. do_donner_deep ) then
          call compute_convective_area     &
                  (Input_mp%tin, Input_mp%pfull, Input_mp%qin, do_uw_conv, &
                   do_donner_deep, C2ls_mp%donner_humidity_area,          &
                   C2ls_mp%donner_humidity_factor, 1.0,    &
                   C2ls_mp%convective_humidity_ratio, &
                   C2ls_mp%convective_humidity_area,   &
                   shallow_cloud_area=   &
                      Moist_clouds_block%cloud_data(i_shallow)%cloud_area, &
                   cell_cld_frac=    &
                      Moist_clouds_block%cloud_data(i_cell)%cloud_area)
        else if (do_donner_deep) then
          call compute_convective_area     &
                  (Input_mp%tin, Input_mp%pfull, Input_mp%qin, do_uw_conv, &
                   do_donner_deep, C2ls_mp%donner_humidity_area,          &
                   C2ls_mp%donner_humidity_factor, 1.0,    &
                   C2ls_mp%convective_humidity_ratio, &
                   C2ls_mp%convective_humidity_area,   &
                   cell_cld_frac=    &
                      Moist_clouds_block%cloud_data(i_cell)%cloud_area)
        else if (do_uw_conv) then
          call compute_convective_area     &
                  (Input_mp%tin, Input_mp%pfull, Input_mp%qin, do_uw_conv, &
                   do_donner_deep, C2ls_mp%donner_humidity_area,          &
                   C2ls_mp%donner_humidity_factor, 1.0,    &
                   C2ls_mp%convective_humidity_ratio, &
                   C2ls_mp%convective_humidity_area,   &
                   shallow_cloud_area=   &
                      Moist_clouds_block%cloud_data(i_shallow)%cloud_area)
        else
          call compute_convective_area     &
                  (Input_mp%tin, Input_mp%pfull, Input_mp%qin, do_uw_conv, &
                   do_donner_deep, C2ls_mp%donner_humidity_area,          &
                   C2ls_mp%donner_humidity_factor, 1.0,    &
                   C2ls_mp%convective_humidity_ratio, &
                   C2ls_mp%convective_humidity_area)
        endif
!<--- h1g, 2017-01-31

      endif
              
!---------------------------------------------------------------------
!    end the timing of the convection code section.
!---------------------------------------------------------------------
      call mpp_clock_end (convection_clock)

!------------------------------------------------------------------------



end subroutine convection_driver



!######################################################################

subroutine cape_cin_diagnostics (is, ie, js, je, Input_mp, Time)

integer,             intent(in) :: is,ie,js,je
type(mp_input_type), intent(in) :: Input_mp
type(time_type),     intent(in) :: Time

      real, dimension(size(Input_mp%tin,1), size(Input_mp%tin,2),  &
                                            size(Input_mp%tin,3)) ::  &
                                                            rin, rp, tp
      real, dimension(size(Input_mp%tin,1), size(Input_mp%tin,2)) ::  &
                                                                 cape, cin
      integer, dimension(size(Input_mp%tin,1), size(Input_mp%tin,2)) ::  &
                                                klcl, klfc, klzb

      logical :: avgbl, used
      integer :: i, j, ix, jx, kx
 
!------------------------------------------------------
!    compute and write out CAPE and CIN.
!------------------------------------------------------
 
      if ( id_cape > 0 .or. id_cin > 0) then
        kx = size(Input_mp%tin,3)
        ix = size(Input_mp%tin,1)
        jx = size(Input_mp%tin,2)

!----------------------------------------------
!    calculate mixing ratio.
!----------------------------------------------
        rin = Input_mp%qin/(1.0 - Input_mp%qin) 

!-----------------------------------------------------------------------
!    call routine to calculate cape and cin.
!-----------------------------------------------------------------------
        avgbl = .false.
        do j = 1,jx 
          do i = 1,ix 
            call capecalcnew   &
                ( kx, Input_mp%pfull(i,j,:), Input_mp%phalf(i,j,:),   &
                  CP_AIR, RDGAS, RVGAS, HLV, KAPPA, Input_mp%tin(i,j,:), &
                  rin(i,j,:), avgbl, cape(i,j), cin(i,j), tp(i,j,:), &
                  rp(i,j,:), klcl(i,j), klfc(i,j), klzb(i,j))
          end do
        end do

!-------------------------------------------------------------------------
!    output diagnostics.
!-------------------------------------------------------------------------
        if (id_cape > 0) used = send_data ( id_cape, cape, Time, is, js )
        if ( id_cin > 0 ) used = send_data ( id_cin, cin, Time, is, js )
        if ( id_tp  > 0 ) used = send_data ( id_tp,  tp, Time, is, js )
        if ( id_rp  > 0 ) used = send_data ( id_rp,  rp, Time, is, js )
        if ( id_lcl > 0 ) used = send_data ( id_lcl, 1.0*klcl, Time,  &
                                                                 is, js )
        if ( id_lfc > 0 ) used = send_data ( id_lfc, 1.0*klfc, Time,  &
                                                                 is, js )
        if ( id_lzb > 0 ) used = send_data ( id_lzb, 1.0*klzb, Time,  &
                                                                 is, js )
      end if

!-----------------------------------------------------------------------



end subroutine cape_cin_diagnostics



!#######################################################################

subroutine diag_field_init ( axes, Time, Control)

integer,         intent(in) :: axes(4)
type(time_type), intent(in) :: Time
type(mp_removal_control_type), intent(in) :: Control

      character(len=32)     :: tracer_units, tracer_name
      character(len=128)    :: diaglname
      integer, dimension(3) :: half = (/1,2,4/)
      integer               :: n, nn

!------------ initialize diagnostic fields in this module -------------

!----- initialize global integrals for netCDF output -----
   id_pr_g = register_global_diag_field ('pr', Time, 'Precipitation', &
                     'kg m-2 s-1', standard_name='precipitation_flux', buffer=.true. )
   id_prc_g = register_global_diag_field ('prc', Time, 'Convective Precipitation', &
                     'kg m-2 s-1', standard_name='convective_precipitation_flux', buffer=.true. )
   id_prsn_g = register_global_diag_field ('prsn', Time, 'Snowfall Flux', 'kg m-2 s-1', &
                                           standard_name='snowfall_flux', buffer=.true. )
!-------------------------------------------------------------------------
!    diagnostics related to total convective tendencies of temperature,
!    vapor and precipitation.
!-------------------------------------------------------------------------
      id_tdt_conv = register_diag_field ( mod_name, &
                   'tdt_conv', axes(1:3), Time, &
                   'Temperature tendency from convection ',  'deg_K/s',  &
                   missing_value=missing_value               )



      ID_tntc = register_cmip_diag_field_3d ( mod_name, 'tntc', Time, &
                  'Tendency of Air Temperature Due to Convection ', 'K s-1', &
                  standard_name='tendency_of_air_temperature_due_to_convection' )

      id_qdt_conv = register_diag_field ( mod_name, &
                  'qdt_conv', axes(1:3), Time, &
                  'Spec humidity tendency from convection ',  'kg/kg/s',  &
                  missing_value=missing_value               )

      ID_tnhusc = register_cmip_diag_field_3d ( mod_name, 'tnhusc', Time, &
                  'Tendency of Specific Humidity Due to Convection ', 's-1', &
                  standard_name='tendency_of_specific_humidity_due_to_convection' )

      id_q_conv_col = register_diag_field ( mod_name, &
                      'q_conv_col', axes(1:2), Time, &
                      'Water vapor path tendency from convection ',  &
                      'kg/m2/s' )
   
      id_t_conv_col = register_diag_field ( mod_name, &
                      't_conv_col', axes(1:2), Time, &
                      'Column static energy tendency from convection ', &
                      'W/m2' )
   
      id_enth_conv_col = register_diag_field ( mod_name, &
                         'enth_conv_col', axes(1:2), Time, &
                         'Column enthalpy tendency from convection',  &
                         'W/m2' )
 
      id_wat_conv_col = register_diag_field ( mod_name, &
                        'wat_conv_col', axes(1:2), Time, &
                        'Column total water tendency from convection', &
                        'kg(h2o)/m2/s' )

      id_prec_conv = register_diag_field ( mod_name, &
                     'prec_conv', axes(1:2), Time, &
                     'Precipitation rate from convection ',    &
                     'kg(h2o)/m2/s', interp_method = "conserve_order1" )

      id_prc = register_cmip_diag_field_2d ( mod_name, 'prc', Time, &
                        'Convective Precipitation',   'kg m-2 s-1', &
                   standard_name = 'convective_precipitation_flux', &
                                  interp_method = "conserve_order1" ) 

      id_prrc = register_cmip_diag_field_2d ( mod_name, 'prrc', Time, &
                            'Convective Rainfall Rate', 'kg m-2 s-1', &
                            standard_name='convective_rainfall_flux', &
                                      interp_method="conserve_order1" )

      id_snow_conv = register_diag_field ( mod_name, &
                     'snow_conv', axes(1:2), Time, &
                     'Frozen precip rate from convection ',  &
                     'kg(h2o)/m2/s', interp_method = "conserve_order1" )

      id_conv_freq = register_diag_field ( mod_name, &
                     'conv_freq', axes(1:2), Time, &
                     'frequency of convection ',       '1', &
                     missing_value = missing_value                       )

      id_prsnc = register_cmip_diag_field_2d ( mod_name, 'prsnc', Time, &
                              'Convective Snowfall Flux', 'kg m-2 s-1', &
                              standard_name='convective_snowfall_flux', &
                                        interp_method="conserve_order1" )

      id_ci = register_cmip_diag_field_2d ( mod_name, 'ci', Time, &
                    'Fraction of Time Convection Occurs in Cell',  '1.0', &
                         standard_name='convection_time_fraction' )

      id_gust_conv = register_diag_field ( mod_name, &
                     'gust_conv', axes(1:2), Time, &
                     'Gustiness resulting from convection ',       'm/s' )

      id_conv_rain3d= register_diag_field ( mod_name, &
                      'conv_rain3d', axes(half), Time, &
                      'Rain fall rate from convection -3D ',    &
                      'kg(h2o)/m2/s' )

      id_conv_snow3d= register_diag_field ( mod_name, &
                      'conv_snow3d', axes(half), Time, &
                      'Snow fall rate from convection -3D',   &
                      'kg(h2o)/m2/s' )

!----------------------------------------------------------------------
!    tendencies of cloud tracers resulting from convection.
!----------------------------------------------------------------------
      if (doing_prog_clouds ) then

        id_qldt_conv = register_diag_field ( mod_name, &
                       'qldt_conv', axes(1:3), Time, &
                       'Liquid water tendency from convection',  &
                       'kg/kg/s', missing_value=missing_value  )

        if (do_liq_num) then
          id_qndt_conv = register_diag_field ( mod_name, &
                         'qndt_conv', axes(1:3), Time, &
                         'Liquid drop tendency from convection', '#/kg/s',&
                         missing_value=missing_value               )
        endif

        id_qidt_conv = register_diag_field ( mod_name, &
                       'qidt_conv', axes(1:3), Time, &
                       'Ice water tendency from convection', 'kg/kg/s',  &
                        missing_value=missing_value               )

        id_qadt_conv = register_diag_field ( mod_name, &
                       'qadt_conv', axes(1:3), Time, &
                       'Cloud fraction tendency from convection', '1/sec',&
                        missing_value=missing_value               )

        id_ql_conv_col = register_diag_field ( mod_name, &
                         'ql_conv_col', axes(1:2), Time, &
                         'Liquid water path tendency from convection',  &
                         'kg/m2/s' )
   
        if (do_liq_num) then
          id_qn_conv_col = register_diag_field ( mod_name, &
                           'qn_conv_col', axes(1:2), Time, &
                           'Liquid drp tendency from convection',  &
                           'kg/m2/s' )
        endif
 
        id_qi_conv_col = register_diag_field ( mod_name, &
                         'qi_conv_col', axes(1:2), Time, &
                         'Ice water path tendency from convection',  &
                         'kg/m2/s' )
   
        id_qa_conv_col = register_diag_field ( mod_name, &
                         'qa_conv_col', axes(1:2), Time, &
                         'Cloud mass tendency from convection', 'kg/m2/s' )
      
        if (do_ice_num) then
          id_qnidt_conv = register_diag_field ( mod_name, &
                          'qnidt_conv', axes(1:3), Time, &
                          'Ice number tendency from convection', '#/kg/s',&
                          missing_value=missing_value               )

          id_qni_conv_col = register_diag_field ( mod_name, &
                            'qni_conv_col', axes(1:2), Time, &
                            'Ice number tendency from convection',   &
                            'kg/m2/s' )
        endif
      endif ! (doing_prog_clouds)

!-----------------------------------------------------------------------
!    diagnostics for cloud base and cloud top.
!-----------------------------------------------------------------------
      id_conv_cld_base = register_diag_field ( mod_name, &
                         'conv_cld_base', axes(1:2), Time, &
                         'pressure at convective cloud base',   'Pa', &
                         mask_variant = .true., &
                         missing_value=missing_value               )

      id_ccb = register_cmip_diag_field_2d ( mod_name, 'ccb', Time, &
                     'Air Pressure at Convective Cloud Base', 'Pa', &
                     standard_name = 'air_pressure_at_convective_cloud_base', &
                     mask_variant = .true. )

      id_conv_cld_top = register_diag_field ( mod_name, &
                        'conv_cld_top', axes(1:2), Time, &
                        'pressure at convective cloud top',   'Pa', &
                        mask_variant = .true., &
                        missing_value=missing_value               )

      id_cct = register_cmip_diag_field_2d ( mod_name, 'cct', Time, &
                      'Air Pressure at Convective Cloud Top', 'Pa', &
                      standard_name = 'air_pressure_at_convective_cloud_top', &
                      mask_variant = .true. )

!-----------------------------------------------------------------------
!    convective mass flux diagnostics.
!-----------------------------------------------------------------------
      id_mc_full = register_diag_field ( mod_name, &
                   'mc_full', axes(1:3), Time, &
                   'Net Mass Flux from convection',   'kg/m2/s', &
                   missing_value=missing_value               )

      id_mc_half = register_diag_field ( mod_name, &
                   'mc_half', axes(half), Time, &
                   'Net Mass Flux from convection on half levs',   &
                   'kg/m2/s', missing_value=missing_value               )

      ID_mc = register_cmip_diag_field_3d ( mod_name, 'mc', Time, &
            'Convective Mass Flux',   'kg m-2 s-1', &
            standard_name='atmosphere_net_upward_convective_mass_flux', &
            axis="half" )
   
      id_mc_conv_up = register_diag_field ( mod_name, &
                      'mc_conv_up', axes(1:3), Time, &
                      'Upward Mass Flux from convection',   'kg/m2/s', &
                      missing_value=missing_value               )

!---------------------------------------------------------------------
!    register diagnostics for lightning NOx.
!---------------------------------------------------------------------
      if (get_tracer_index(MODEL_ATMOS,'no') > 0) &
        id_prod_no = register_diag_field ( 'tracers', &
                     'hook_no', axes(1:3), Time, &
                     'hook_no',   'molec/cm3/s')

!-------------------------------------------------------------------------
!    register diagnostics specific to the Betts-Miller experiments.
!-------------------------------------------------------------------------
      if ( any((/do_bm, do_bmmass, do_bmomp/)) ) then
        id_qref = register_diag_field ( mod_name, &
                  'qref', axes(1:3), Time, &
                  'Adjustment reference specific humidity profile', &
                  'kg/kg',  missing_value=missing_value               )

        id_tref = register_diag_field ( mod_name, &
                  'tref', axes(1:3), Time, &
                  'Adjustment reference temperature profile', &
                  'K',  missing_value=missing_value                   )

        id_bmflag = register_diag_field (mod_name, &
                    'bmflag', axes(1:2), Time, &
                    'Betts-Miller flag', &
                    'no units', missing_value=missing_value            )

        id_klzbs  = register_diag_field  (mod_name, &
                    'klzbs', axes(1:2), Time, &
                    'klzb', &
                    'no units', missing_value=missing_value            )

      endif

      id_cape = register_diag_field ( mod_name, &
                'cape', axes(1:2), Time, &
                'Convectively available potential energy',      'J/Kg')
      
      id_cin = register_diag_field ( mod_name, &
               'cin', axes(1:2), Time, &                       
               'Convective inhibition',                        'J/Kg')

      id_tp = register_diag_field ( mod_name, &
              'tp', axes(1:3), Time, &
              'Temperature of lifted parcel',                    'K')
      id_rp = register_diag_field ( mod_name, &
              'rp', axes(1:3), Time, &
              'Humidity of lifted parcel',                   'kg/kg')
      id_lcl = register_diag_field ( mod_name, &
               'klcl', axes(1:2), Time, &
               'Index of LCL',                        'none')
      id_lfc = register_diag_field ( mod_name, &
               'klfc', axes(1:2), Time, &
               'Index of LFC',                        'none')
      id_lzb = register_diag_field ( mod_name, &
               'klzb', axes(1:2), Time, &
               'Index of LZB',                        'none')

      if (do_bm ) then
        id_invtaubmt  = register_diag_field  (mod_name, &
                        'invtaubmt', axes(1:2), Time, &
                        'Inverse temperature relaxation time', &
                        '1/s', missing_value=missing_value            )

        id_invtaubmq = register_diag_field  (mod_name, &
                       'invtaubmq', axes(1:2), Time, &
                       'Inverse humidity relaxation time', &
                       '1/s', missing_value=missing_value            )
      end if  ! if ( do_bm )

      if (do_bmmass) then
        id_massflux = register_diag_field (mod_name, &
                      'massflux', axes(1:3), Time, &
                      'Massflux implied by temperature adjustment', &
                      'm/s', missing_value=missing_value                 )
      end if  ! if ( do_bmmass )

!------------------------------------------------------------------------
!    register diagnostics specific to the ras parameterization.
!------------------------------------------------------------------------
!RSH activate this if when convection code redone:
!     if (do_ras) then
        id_ras_precip = register_diag_field ( mod_name, &
                       'ras_precip', axes(1:2), Time, &
                       'Precipitation rate from ras ',       'kg/m2/s' )

        id_ras_freq = register_diag_field ( mod_name, &
                      'ras_freq', axes(1:2), Time, &
                      'frequency of precip from ras ',       'number' , &
                      missing_value = missing_value                       )
!     endif

!------------------------------------------------------------------------
!    register diagnostics specific to the donner parameterization.
!------------------------------------------------------------------------
!RSH activate this if when convection code redone:
! following if activated 10/8/16:
      if (do_donner_deep) then
        id_don_precip = register_diag_field ( mod_name, &
                        'don_precip', axes(1:2), Time, &
                        'Precipitation rate from donner ',      'kg/m2/s' )

        id_don_freq = register_diag_field ( mod_name, &
                      'don_freq', axes(1:2), Time, &
                      'frequency of precip from donner ',       'number', &
                      missing_value = missing_value                       )

        id_enth_donner_col2 = register_diag_field ( mod_name, &
                              'enth_donner_col2', axes(1:2), Time, &
                              'column enthalpy tendency from Donner liq&
                              & precip','W/m2' )
 
        id_enth_donner_col3 = register_diag_field ( mod_name, &
                              'enth_donner_col3', axes(1:2), Time, &
                              'Column enthalpy tendency from Donner &
                              &frzn precip','W/m2' )
 
        id_enth_donner_col4 = register_diag_field ( mod_name, &
                              'enth_donner_col4', axes(1:2), Time, &
                              'Atmospheric column enthalpy tendency from&
                              & Donner convection', 'W/m2' )
 
        id_enth_donner_col5 = register_diag_field ( mod_name, &
                              'enth_donner_col5', axes(1:2), Time, &
                              'Column enthalpy tendency due to condensate&
                              & xfer from Donner to lsc','W/m2' )

        id_enth_donner_col6 = register_diag_field ( mod_name, &
                              'enth_donner_col6', axes(1:2), Time, &
                              'Column enthalpy tendency from donner &
                              &moisture  conservation  adjustment','W/m2' )
 
        id_enth_donner_col7 = register_diag_field ( mod_name, &
                              'enth_donner_col7', axes(1:2), Time, &
                              'Precip adjustment needed to balance donner&
                              & moisture  adjustment','kg(h2o)/m2/s' )

        id_enth_donner_col = register_diag_field ( mod_name, &
                             'enth_donner_col', axes(1:2), Time, &
                             'Column enthalpy imbalance from Donner &
                             &convection','W/m2' )

        id_wat_donner_col = register_diag_field ( mod_name, &
                            'wat_donner_col', axes(1:2), Time, &
                            'Column total water tendency from Donner&
                            & convection','kg(h2o)/m2/s' )
  
        id_enth_mca_donner_col = register_diag_field ( mod_name, &
                                 'enth_mca_donner_col', axes(1:2), Time, &
                                 'Column enthalpy imbalance from Donner&
                                 & MCA convection','W/m2' )

        id_wat_mca_donner_col = register_diag_field ( mod_name, &
                                'wat_mca_donner_col', axes(1:2), Time, &
                                'Column total water imbalance from Donner&
                                & MCA convection', 'kg(h2o)/m2/s' )

        id_scale_donner = register_diag_field ( mod_name, &
                         'scale_donner', axes(1:2), Time, &
!RSH: FIx in final version:
!                        'Scaling factor applied to donner convection&
                         'Scaling factor applied to UW convection&
                         & tendencies','1' )

        id_tdt_deep_donner= register_diag_field ( mod_name, &
                            'tdt_deep_donner', axes(1:3), Time, &
                            ' heating rate - deep portion', 'deg K/s', &
                            missing_value=missing_value               )

        id_qdt_deep_donner = register_diag_field ( mod_name, &
                             'qdt_deep_donner', axes(1:3), Time, &
                             ' moistening rate - deep portion', 'kg/kg/s',&
                             missing_value=missing_value               )

        id_qadt_deep_donner = register_diag_field ( mod_name, &
                              'qadt_deep_donner', axes(1:3), Time, &
                              ' cloud amount tendency - deep portion',  &
                              '1/s', missing_value=missing_value      )

        id_qldt_deep_donner = register_diag_field ( mod_name, &
                              'qldt_deep_donner', axes(1:3), Time, &
                              ' cloud liquid tendency - deep portion',  &
                              'kg/kg/s', missing_value=missing_value    )

        id_qidt_deep_donner = register_diag_field ( mod_name, &
                              'qidt_deep_donner', axes(1:3), Time, &
                              ' ice water tendency - deep portion',  &
                              'kg/kg/s', missing_value=missing_value    )
        if (do_liq_num) &
          id_qndt_deep_donner = register_diag_field ( mod_name, &
                                'qndt_deep_donner', axes(1:3), Time, &
                                'deep convection cloud drop tendency',  &
                                '#/kg/s', missing_value=missing_value    )

        if (do_ice_num) &
          id_qnidt_deep_donner = register_diag_field ( mod_name, &
                                 'qnidt_deep_donner', axes(1:3), Time, &
                                 ' ice number tendency - deep portion', &
                                 '#/kg/s', missing_value=missing_value   )

        id_tdt_mca_donner = register_diag_field ( mod_name, &
                            'tdt_mca_donner', axes(1:3), Time, &
                            ' heating rate - mca  portion', 'deg K/s', &
                            missing_value=missing_value               )

        id_qdt_mca_donner = register_diag_field ( mod_name, &
                            'qdt_mca_donner', axes(1:3), Time, &
                            ' moistening rate - mca  portion', 'kg/kg/s', &
                            missing_value=missing_value               )

        id_prec_deep_donner = register_diag_field ( mod_name, &
                              'prc_deep_donner', axes(1:2), Time, &
                              ' total precip rate - deep portion',  &
                              'kg/m2/s', missing_value=missing_value, &
                              interp_method = "conserve_order1"        )

        id_precret_deep_donner = register_diag_field ( mod_name, &
                                 'prc_ret_deep_donner', axes(1:2), Time, &
                                 ' precip_returned - per timestep',   &
                                 'kg/m2/timestep', &
                                 missing_value=missing_value, &
                                 interp_method = "conserve_order1"     )

        id_prec1_deep_donner = register_diag_field ( mod_name, &
                               'prc1_deep_donner', axes(1:2), Time, &
                               ' change in precip for conservation&
                               & in donner', 'kg/m2/s ', &
                               missing_value=missing_value,  &
                               mask_variant = .true., &
                               interp_method = "conserve_order1"  )

        id_prec_mca_donner = register_diag_field ( mod_name, &
                             'prc_mca_donner', axes(1:2), Time, &
                             ' total precip rate - mca  portion',  &
                             'kg/m2/s', missing_value=missing_value, &
                             interp_method = "conserve_order1"        )

        id_snow_deep_donner = register_diag_field ( mod_name, &
                              'snow_deep_donner', axes(1:2), Time, &
                              ' frozen precip rate - deep portion',   &
                              'kg/m2/s', missing_value=missing_value, &
                              interp_method = "conserve_order1"        )

        id_snow_mca_donner = register_diag_field ( mod_name, &
                             'snow_mca_donner', axes(1:2), Time, &
                             ' frozen precip rate -  mca portion',  &
                             'kg/m2/s', missing_value=missing_value, &
                             interp_method = "conserve_order1"        )

        id_mc_donner = register_diag_field ( mod_name, &
                       'mc_donner', axes(1:3), Time, &
                       'Net Mass Flux from donner',   'kg/m2/s', &
                       missing_value=missing_value               )

        id_mc_donner_half = register_diag_field ( mod_name, &
                            'mc_donner_half', axes(half), Time, &
                            'Net Mass Flux from donner at half levs',  &
                            'kg/m2/s', missing_value=missing_value      )

        id_m_cdet_donner = register_diag_field ( mod_name, &
                           'm_cdet_donner', axes(1:3), Time, &
                           'Detrained Cell Mass Flux from donner',  &
                           'kg/m2/s', missing_value=missing_value      )

        id_m_cellup = register_diag_field ( mod_name, &
                     'm_cellup', axes(half), Time, &
                     'Upward Cell Mass Flux from donner',   'kg/m2/s', &
                     missing_value=missing_value               )

        id_cell_cld_frac = register_diag_field ( mod_name, &
                           'cell_cld_frac', axes(1:3), Time, & 
                           'cell cloud fraction from donner',   '', &
                           missing_value=missing_value               )

        id_meso_cld_frac = register_diag_field ( mod_name, &
                           'meso_cld_frac', axes(1:3), Time, & 
                           'meso-scale cloud fraction from donner',   '', &
                           missing_value=missing_value               )

        id_donner_humidity_area = register_diag_field ( mod_name, &
                                  'donner_humidity_area', axes(1:3), Time,&
                                  'donner humidity area',  '', &
                                   missing_value=missing_value          )

        if (do_donner_conservation_checks) then

          id_enthint = register_diag_field    &
                       (mod_name, 'enthint_don', axes(1:2), Time,  &
                       'atmospheric column enthalpy change from donner', &
                       'W/m2', missing_value=missing_value)

          id_lcondensint = register_diag_field    &
                           (mod_name, 'lcondensint_don', axes(1:2), Time, &
                           'enthalpy transferred by condensate from &
                           &donner to lscale', 'W/m2',  &
                           missing_value=missing_value)

          id_lprcp = register_diag_field    &
                     (mod_name, 'lprcpint_don', axes(1:2),   &
                     Time, 'enthalpy removed by donner precip', 'W/m2',   &
                     missing_value=missing_value)

          id_vertmotion = register_diag_field    &
                          (mod_name, 'vertmotion_don', axes(1:2), Time,  &
                          'enthalpy change due to cell and meso motion &
                          &in donner', 'W/m2',  &
                          missing_value=missing_value)

          id_enthdiffint = register_diag_field    &
                           (mod_name, 'enthdiffint_don', axes(1:2),   &
                           Time, 'enthalpy  imbalance due to donner',  &
                           'W/m2', missing_value=missing_value)

          id_vaporint = register_diag_field    &
                        (mod_name, 'vaporint_don', axes(1:2),   &
                        Time, 'column water vapor change', 'kg(h2o)/m2/s',&
                        missing_value=missing_value)

          id_max_enthalpy_imbal_don = register_diag_field    &
                                      (mod_name, 'max_enth_imbal_don',  &
                                      axes(1:2), Time,   &
                                      'max enthalpy  imbalance from&
                                      & donner', 'W/m**2',  &
                                      missing_value=missing_value)

          id_max_water_imbal_don = register_diag_field    &
                                   (mod_name, 'max_water_imbal_don', &
                                   axes(1:2), Time, 'max water imbalance&
                                   & from donner', 'kg(h2o)/m2/s', &
                                   missing_value=missing_value)

          id_condensint = register_diag_field    &
                          (mod_name, 'condensint_don', axes(1:2), Time,  &
                          'column condensate exported from donner&
                          & to lscale', 'kg(h2o)/m2/s',  &
                          missing_value=missing_value )

          id_precipint = register_diag_field    &
                         (mod_name, 'precipint_don', axes(1:2),   &
                         Time, 'column precip from donner',  &
                         'kg(h2o)/m2/s', missing_value=missing_value)

          id_diffint= register_diag_field    &
                      (mod_name, 'diffint_don', axes(1:2),   &
                      Time, 'water imbalance due to donner', &
                      'kg(h2o)/m2/s', missing_value=missing_value)

        endif
      endif

!------------------------------------------------------------------------
!    register diagnostics specific to the uw parameterization.
!------------------------------------------------------------------------
!RSH activate this if when convection code redone:
!     if (do_uw_conv) then
        id_uw_precip = register_diag_field ( mod_name, &
                       'uw_precip', axes(1:2), Time, &
                       'Precipitation rate from uw shallow',  'kg/m2/s', &
                       interp_method = "conserve_order1" )

        id_uw_snow = register_diag_field ( mod_name, &
                     'uw_snow', axes(1:2), Time, &
                     'Snow rate from uw shallow',       'kg/m2/s' , &
                     interp_method = "conserve_order1" )

        id_uw_freq = register_diag_field ( mod_name, &
                     'uw_freq', axes(1:2), Time, &
                     'frequency of precip from uw shallow ',  'number' , &
                     missing_value = missing_value                       )

        id_enth_uw_col = register_diag_field ( mod_name, &
                         'enth_uw_col', axes(1:2), Time, &
                         'Column enthalpy tendency from UW convection', &
                         'W/m2' )
 
        id_wat_uw_col = register_diag_field ( mod_name, &
                        'wat_uw_col', axes(1:2), Time, &
                        'Column total water tendency from UW convection',&
                        'kg(h2o)/m2/s' )

        id_scale_uw = register_diag_field ( mod_name, &
                      'scale_uw', axes(1:2), Time, &
                      'Scaling factor applied to UW convection&
                      & tendencies','1' )
          
        id_tdt_uw = register_diag_field ( mod_name, &
                    'tdt_uw', axes(1:3), Time, &
                    'UW convection heating rate', 'deg K/s', &
                    missing_value=missing_value               )

        id_qdt_uw = register_diag_field ( mod_name, &
                    'qdt_uw', axes(1:3), Time, &
                    'UW convection moistening rate', 'kg/kg/s', &
                    missing_value=missing_value               )

        id_qadt_uw = register_diag_field ( mod_name, &
                     'qadt_uw', axes(1:3), Time, &
                     'UW convection cloud amount tendency', '1/s', &
                     missing_value=missing_value               )

        id_qldt_uw = register_diag_field ( mod_name, &
                     'qldt_uw', axes(1:3), Time, &
                     'UW convection cloud liquid tendency', 'kg/kg/s', &
                     missing_value=missing_value               )

        id_qidt_uw = register_diag_field ( mod_name, &
                     'qidt_uw', axes(1:3), Time, &
                     'UW convection ice water tendency', 'kg/kg/s', &
                     missing_value=missing_value               )

        if (do_liq_num) &
          id_qndt_uw = register_diag_field ( mod_name, &
                       'qndt_uw', axes(1:3), Time, &
                       'UW convection cloud drop tendency', '#/kg/s', &
                       missing_value=missing_value               )

        if (do_ice_num) &
          id_qnidt_uw = register_diag_field ( mod_name, &
                        'qnidt_uw', axes(1:3), Time, &
                        'UW convection ice number tendency', '#/kg/s', &
                        missing_value=missing_value               )

!     endif

!----------------------------------------------------------------------
!    dry adjustment diagnostic.
!----------------------------------------------------------------------
      id_tdt_dadj = register_diag_field ( mod_name, &
                    'tdt_dadj', axes(1:3), Time, &
                    'Temperature tendency from dry conv adj', 'deg_K/s',  &
                    missing_value=missing_value               )

!---------------------------------------------------------------------
!    allocate and initialize arrays to hold the diagnostic ids for each 
!    active tracer. diagnostics for tendency due to convection, 
!    column tendency due to convection, the tracer amount and tracer 
!    column amount are available.
!---------------------------------------------------------------------
      allocate (id_tracerdt_conv    (num_prog_tracers))
      allocate (id_tracerdt_conv_col(num_prog_tracers))
      allocate (id_conv_tracer           (num_prog_tracers))
      allocate (id_conv_tracer_col(num_prog_tracers))

      id_tracerdt_conv = NO_TRACER
      id_tracerdt_conv_col = NO_TRACER
      id_conv_tracer = NO_TRACER
      id_conv_tracer_col = NO_TRACER
 
!------------------------------------------------------------------------
!    define the diagnostics names that are requested and register the
!    diagnostics for those tracers that were specified to be affected
!    by a convection scheme.
!------------------------------------------------------------------------
      do n = 1,num_prog_tracers
        call get_tracer_names (MODEL_ATMOS, n, name = tracer_name,  &
                                                  units = tracer_units)
        if (Control%tracers_in_donner(n) .or. &
            Control%tracers_in_ras(n)      .or.  &
            Control%tracers_in_mca(n)      .or.  &
            Control%tracers_in_uw(n)) then
          diaglname = trim(tracer_name)//  &
                        ' total tendency from moist convection'
          id_tracerdt_conv(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'dt_conv',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

          diaglname = trim(tracer_name)//  &
                       ' total path tendency from moist convection'
          id_tracerdt_conv_col(n) =  &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'dt_conv_col', &
                         axes(1:2), Time, trim(diaglname), &
                         TRIM(tracer_units)//'*(kg/m2)/s',   &
                         missing_value=missing_value)
        endif
 
!----------------------------------------------------------------------
!    output the distribution and column values of any tracer for which
!    they are requested, even if not transported by convection.
!----------------------------------------------------------------------
        diaglname = trim(tracer_name)
        id_conv_tracer(n) =    &
                        register_diag_field ( mod_name_tr, &
                        TRIM(tracer_name),  &
                        axes(1:3), Time, trim(diaglname), &
                        TRIM(tracer_units)      ,  &
                        missing_value=missing_value)
        diaglname =  ' column integrated' // trim(tracer_name)
        id_conv_tracer_col(n) =  &
                        register_diag_field ( mod_name_tr, &
                        TRIM(tracer_name)//'_col', &
                        axes(1:2), Time, trim(diaglname), &
                        TRIM(tracer_units)      ,   &
                        missing_value=missing_value)
      end do

!------------------------------------------------------------------
!    register the diagnostics which will report the tendencies due to
!    mca component of donner convection.
!------------------------------------------------------------------
      if (do_donner_deep) then
        allocate (id_tracerdt_mcadon  (Control%num_donner_tracers))
        allocate (id_tracerdt_mcadon_col(Control%num_donner_tracers))
 
        nn = 1
        do n = 1,num_prog_tracers
          call get_tracer_names (MODEL_ATMOS, n, name = tracer_name,  &
                                                    units = tracer_units)
          if (Control%tracers_in_donner(n) ) then
            diaglname = trim(tracer_name)//  &
                       ' tendency from donner-mca'
            id_tracerdt_mcadon(nn) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_donmca',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                        missing_value=missing_value)

            diaglname = trim(tracer_name)//  &
                       ' total path tendency from donner-mca'
            id_tracerdt_mcadon_col(nn) =  &
                        register_diag_field ( mod_name, &
                        TRIM(tracer_name)//'_donmca_col', &
                        axes(1:2), Time, trim(diaglname), &
                          TRIM(tracer_units)//'*(kg/m2)/s',   &
                        missing_value=missing_value)
            nn = nn + 1
          endif
        end do
      endif


!---------------------------------------------------------------------


end subroutine diag_field_init



!######################################################################

subroutine convection_driver_time_vary (dt)

real, intent(in) :: dt

!--------------------------------------------------------------------
      if (do_donner_deep) then
        call donner_deep_time_vary (dt)
      endif

!--------------------------------------------------------------------


end subroutine convection_driver_time_vary



!######################################################################

subroutine convection_driver_endts 

!-----------------------------------------------------------------------

      if (do_donner_deep) then
        call donner_deep_endts
      endif 

!-----------------------------------------------------------------------


end subroutine convection_driver_endts


!######################################################################

subroutine convection_driver_end 

!--------------------------------------------------------------------- 
!    call the destructor routines for the active convection modules.
!--------------------------------------------------------------------- 
      if (do_donner_deep) call donner_deep_end 
      if (do_ras        ) call ras_end 
      if (do_uw_conv    ) call uw_conv_end 
      if (do_cmt        ) call cu_mo_trans_end 
      call detr_ice_num_end 

!----------------------------------------------------------------------
!    deallocate module variables.
!----------------------------------------------------------------------
      if (do_donner_deep .and. do_donner_conservation_checks) then 
        deallocate (max_water_imbal_don) 
        deallocate (max_enthalpy_imbal_don)
      endif

      deallocate (cloud_tracer)

      deallocate(id_tracerdt_conv)        ! h1g, 2017-02-02
      deallocate (id_tracerdt_conv_col)   ! h1g, 2017-02-02
      deallocate (id_conv_tracer)         ! h1g, 2017-02-02
      deallocate (id_conv_tracer_col)     ! h1g, 2017-02-02
      if (do_donner_deep) then            ! h1g, 2017-02-02
        deallocate ( id_tracerdt_mcadon )        ! h1g, 2017-02-02
        deallocate ( id_tracerdt_mcadon_col )    ! h1g, 2017-02-02
      endif


!--------------------------------------------------------------------

   
 end subroutine convection_driver_end



!######################################################################

subroutine convection_driver_restart (timestamp)

character(len=*), intent(in), optional :: timestamp

      if (do_donner_deep) call donner_deep_restart(timestamp)


end subroutine convection_driver_restart


!######################################################################

subroutine prevent_neg_precip_fluxes (fluxh)

real, dimension(:,:,:), intent(inout) :: fluxh

      integer :: i,j,k
      real, dimension(size(fluxh,1), size(fluxh,2)) :: sumneg

!-----------------------------------------------------------------------
!    move down each column looking for negative precip fluxes at each
!    level. if found, the negative flux is eliminated by reducing the 
!    incoming precip flux from above.
!-----------------------------------------------------------------------
      sumneg(:,:) = 0.
      do k=2, size(fluxh,3)
        do j=1,size(fluxh,2)
          do i=1,size(fluxh,1)
            if (fluxh(i,j,k) > 0.0) then
              if (fluxh(i,j,k) > ABS(sumneg(i,j))) then
                fluxh(i,j,k) = fluxh(i,j,k) + sumneg(i,j)
                sumneg(i,j) = 0.
              else
                sumneg(i,j) = sumneg(i,j) + fluxh(i,j,k)
                fluxh(i,j,k) = 0.
              endif
            else
              sumneg(i,j) = sumneg(i,j) + fluxh(i,j,k)
              fluxh(i,j,k) = 0.
            endif
          end do
        end do
      end do

!----------------------------------------------------------------------


end subroutine prevent_neg_precip_fluxes



!#######################################################################

subroutine compute_convective_area     &
          (t, pfull, q, do_uw_conv, do_donner_deep, donner_humidity_area,&
           donner_humidity_factor, max_cnv_frac, humidity_ratio,   &
           convective_area, shallow_cloud_area, cell_cld_frac)

!-------------------------------------------------------------------------
!    subroutine compute_convective_area defines the grid box area affected
!    by the convective clouds and the ratio of the grid-box relative
!    humidity to the humidity in the environment of the convective
!    clouds. 
!-------------------------------------------------------------------------

real, dimension(:,:,:), intent(in)          :: t, q, pfull,   &
                                               donner_humidity_area, &
                                               donner_humidity_factor
logical,                intent(in)          :: do_uw_conv, do_donner_deep
real,                   intent(in)          :: max_cnv_frac
real, dimension(:,:,:), intent(out)         :: humidity_ratio,  &
                                               convective_area
real, dimension(:,:,:), intent(in),optional :: shallow_cloud_area, &
                                               cell_cld_frac

!------------------------------------------------------------------------
!   local variables:
!

      real, dimension(size(t,1), size(t,2), size(t,3)) :: qs, qrf, &
                                                          env_fraction, &
                                                          env_qv
      integer :: i,j,k
      integer :: ix, jx, kx

!-----------------------------------------------------------------------
!    define array dimensions.
!-----------------------------------------------------------------------
      ix = size(t,1)
      jx = size(t,2)
      kx = size(t,3)

!----------------------------------------------------------------------
!    define a realizable grid box specific humidity (qrf) and the 
!    saturation specific humidity (qs).
!------------------------------------------------------------------
      qrf = MAX(q, 0.0)
      call compute_qs (t, pfull, qs)

!----------------------------------------------------------------------
!    define the grid box area whose humidity is affected by the 
!    convective clouds (convective_area). define the environmental 
!    rh (env_qv) which is the value needed in the non-convective cloud
!    area in order to have the computed grid-box relative humidity. the
!    convective cloud area is assumed saturated for the uw clouds, in the
!    donner cell clouds and in the region of donner meso updraft, but is
!    assumed subsaturated in the donner meso downdraft layer above cloud 
!    base, with the degree of saturation given by the 
!    donner_humidity_factor (mesoscale area times assumed RH).
!-------------------------------------------------------------------
      if (do_uw_conv .and. do_donner_deep) then
        convective_area = donner_humidity_area + shallow_cloud_area
        env_qv = qrf - qs*(cell_cld_frac + donner_humidity_factor +  &
                                                shallow_cloud_area)
      else if (do_donner_deep) then
        convective_area = donner_humidity_area
        env_qv = qrf - qs*(cell_cld_frac + donner_humidity_factor)
      else if (do_uw_conv) then
        convective_area = shallow_cloud_area
        env_qv = qrf - shallow_cloud_area*qs
      else
        convective_area = 0.0
        env_qv = qrf
      endif

!-----------------------------------------------------------------------
!    define the convective area, limited by nml parameter max_cnv_frac.
!    the convective environment fraction is the remainder of the box.
!-----------------------------------------------------------------------
      do k=1, kx
        do j=1,jx   
          do i=1,ix  
            convective_area(i,j,k) = min (convective_area(i,j,k), &
                                                             max_cnv_frac)
            env_fraction(i,j,k) = 1.0 - convective_area(i,j,k)

!---------------------------------------------------------------------
!    define the ratio of the grid-box relative humidity to the humidity
!    in the environment of the convective clouds. this can be done only
!    if the grid box contains vapor (qrf > 0.0) and there is some vapor
!    outside of the convective clouds (env_qv > 0.).
!----------------------------------------------------------------------
            if (qrf(i,j,k) /= 0.0 .and. env_qv(i,j,k) > 0.0) then
 
!--------------------------------------------------------------------
!    there must also be grid box area not filled with convective clouds.
!--------------------------------------------------------------------  
              if (env_fraction(i,j,k) > 0.0) then
                humidity_ratio(i,j,k) =    &
                   MAX (qrf(i,j,k)*env_fraction(i,j,k)/env_qv(i,j,k), 1.0)
 
!---------------------------------------------------------------------
!    if the grid box is filled with convective clouds, set humidity ratio 
!    to a flag value.
!----------------------------------------------------------------------
              else
                humidity_ratio(i,j,k) = -10.0
              endif

!--------------------------------------------------------------------
!    if there either is no vapor in the gridbox or the vapor has been 
!    taken up by the convective clouds, set the humidity_ratio to 1.0.
!---------------------------------------------------------------------
            else
              humidity_ratio(i,j,k) = 1.0
            endif
          end do
        end do
      end do

!------------------------------------------------------------------------


end subroutine compute_convective_area


!#######################################################################

end module convection_driver_mod



