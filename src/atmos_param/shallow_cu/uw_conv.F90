#include <fms_platform.h>
MODULE UW_CONV_MOD

  use           mpp_mod, only : mpp_pe, mpp_root_pe, stdlog
  use      Constants_Mod, ONLY: tfreeze,HLv,HLf,HLs,CP_AIR,GRAV,Kappa,rdgas,rvgas
  use   Diag_Manager_Mod, ONLY: register_diag_field, send_data
  use   Time_Manager_Mod, ONLY: time_type, get_time
  use           mpp_mod, only : input_nml_file
  use           fms_mod, only : write_version_number, open_namelist_file, check_nml_error,&
                                FILE_EXIST, ERROR_MESG,  &
                                lowercase, &
                                CLOSE_FILE, FATAL, NOTE
  use  field_manager_mod, only: MODEL_ATMOS
  use  tracer_manager_mod, only: get_tracer_names, query_method, &
                                 get_tracer_index, NO_TRACER
  use atmos_cmip_diag_mod, only: register_cmip_diag_field_2d
  use  sat_vapor_pres_mod,only : sat_vapor_pres_init
  use atmos_tracer_utilities_mod, only : get_wetdep_param
  use moist_proc_utils_mod, only : mp_nml_type

  use  aerosol_types_mod, only : aerosol_type

  use  aer_ccn_act_mod, only :   aer_ccn_act_init
  use  conv_utilities_mod,only :   uw_params_init
  use  conv_utilities_k_mod,only : sd_init_k, sd_copy_k, sd_end_k,  &
                                   ac_init_k, ac_clear_k, ac_end_k, &
                                   pack_sd_k, adi_cloud_k, extend_sd_k,&
                                   exn_init_k, exn_end_k, findt_init_k,&
                                   findt_end_k, &
                                   check_tracer_realizability, &
                                   qt_parcel_k, qt_parcel_deep_k, &
                                   adicloud, sounding, uw_params

  use  conv_plumes_k_mod,only    : cp_init_k, cp_end_k, cp_clear_k, &
                                   ct_init_k, ct_end_k, ct_clear_k, &
                                   cumulus_tend_k, cumulus_plume_k, &
                                   cplume, ctend, cpnlist, cwetdep_type

  use  conv_closures_mod,only    : cclosure_bretherton,   &
                                   cclosure_relaxcbmf, &
                                   cclosure_relaxwfn,  &
                                   cclosure_implicit, cclosure

  use  deep_conv_mod,only        : deepc, cpn_copy, dpconv0, dpconv1, dpconv2, dpconv3
  use random_numbers_mod,    only: randomNumberStream, getRandomNumbers, &
                                   initializeRandomNumberStream, constructSeed

!---------------------------------------------------------------------
  implicit none
  private
!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

!---------------------------------------------------------------------
!-------  interfaces --------

  public  :: uw_conv, uw_conv_init, uw_conv_end

  real, parameter :: aday = 1.
  real, parameter :: mv = -999.
  logical         :: module_is_initialized = .false.

  character(len=7) :: mod_name = 'uw_conv'

  !namelist parameters for UW convection scheme
  integer :: iclosure = 0      ! 0: Bretherton UWShCu orginal / -CIN/TKE based
                               ! 1: Emanuel-Rayment: quasiequilibrium PBL
  real    :: rkm_sh1  = 10.0
  real    :: rkm_sh   = 3.0    ! fractional lateral mixing rate for shallow
  real    :: cldhgt_max   = 50.e3
  real    :: cldhgt_max_shallow = 0.
  real    :: landfact_m   = 0.5
  integer :: idpchoice = 0
  logical :: do_deep = .false.
  logical :: do_coldT = .true.
  logical :: do_lands = .false.
  logical :: do_peff_land = .false.
  logical :: do_uwcmt = .false.
  logical :: do_fast  = .false.
  logical :: do_ice   = .true.
  logical :: do_ppen  = .true.
  logical :: do_forcedlifting = .false.
  logical :: do_minmse = .false.
  real    :: atopevap = 0.
  logical :: apply_tendency = .true.
  logical :: prevent_unreasonable = .true.
  logical :: reproduce_old_version = .false.
  real    :: aerol = 1.e-12
  real    :: tkemin   = 1.e-6
  real    :: wmin_ratio = 0.05
  logical :: do_auto_aero = .false.
  logical :: do_rescale   = .false.
  logical :: do_rescale_t = .false.
  logical :: do_debug     = .false.
  real    :: cush_ref     = 0.
  real    :: plev_cin     = 60000.
  real    :: plev_omg     = 60000.
  real    :: pblht0 = 500.
  real    :: lofactor0 = 1.
  integer :: lochoice  = 0
  real    :: wrel_min = 1.
  real    :: N0 = 100.e6
  logical :: do_umf_pbl = .false.
  logical :: do_qctflx_zero = .false.
  logical :: do_hlflx_zero  = .true.
  logical :: do_new_qnact   = .false.
  logical :: do_2nd_act     = .true.
  logical :: do_downdraft   = .false.
  logical :: do_conv_micro_N = .false.
  logical :: do_varying_rpen = .false.
  logical :: do_new_convcld   = .true.
  logical :: do_subcloud_flx = .false.
  logical :: do_new_subflx   = .false.
  logical :: do_detran_zero = .false.
  logical :: do_prog_gust = .false.
  logical :: do_gust_qt = .false.
  logical :: use_new_let = .false.
  logical :: use_lcl_only =.false.
  logical :: do_new_pevap =.false.
  logical :: stop_at_let  =.false.
  logical :: use_pblhttke_avg =.false.
  logical :: use_hlqtsrc_avg =.false.
  logical :: use_capecin_avg =.false.
  logical :: zero_out_conv_area = .false.
  logical :: do_mse_budget = .false.
  logical :: do_eis_limit  = .false.
  logical :: do_eis_limitn = .false.
  logical :: do_lts_limit  = .false.
  logical :: do_lts_limitn = .false.
  integer :: src_choice = 0
  integer :: gqt_choice = 0
  integer :: rpen_choice = 0
  real    :: eis_max    = 10.
  real    :: eis_min    = 0.
  real    :: plev_for   = 50000.
  real    :: plev_umf   = 70000.
  real    :: shallow_umf_thresh = 0.001
  logical :: do_plev_umf = .false.
  real    :: duration = 10800
  real    :: tau_gust = 7200
  real    :: geff   = 1.
  real    :: cgust0 = 1.
  real    :: cgust_max = 10.
  real    :: sigma0 = 0.5
  real    :: tmax0  = 363.15

  integer :: tracer_check_type = -999 !legacy
  !< select realizability checks to be applied to tracers
  !! -999 (default): apply min/max checks (with scaling of tendencies), no filling
  !!              1: omit min/max checks, apply filling (using sjl_fillz), no scaling
  !!              2: omit min/max checks, no filling (apply scaling to avoid negatives)

  logical :: use_turb_tke = .false.  !h1g, 2015-08-11

  NAMELIST / uw_conv_nml / iclosure, rkm_sh1, rkm_sh, cldhgt_max, plev_cin, plev_omg, eis_max,eis_min,do_peff_land, &
       do_deep, idpchoice, do_coldT, do_lands, do_uwcmt, do_varying_rpen, do_new_convcld, do_mse_budget, rpen_choice, &
       do_fast, do_ice, do_ppen, do_forcedlifting, do_gust_qt, use_new_let, do_hlflx_zero, do_new_qnact, do_2nd_act, N0, &
       atopevap, apply_tendency, prevent_unreasonable, do_minmse, aerol, tkemin, cldhgt_max_shallow,do_umf_pbl,&
       wmin_ratio,  landfact_m, pblht0, lofactor0, lochoice, do_downdraft, &
       do_auto_aero, do_rescale, do_rescale_t, wrel_min, do_debug, tmax0, do_conv_micro_N, &
       cush_ref, do_prog_gust, tau_gust, geff, cgust0, cgust_max, sigma0,  do_qctflx_zero, do_detran_zero, &
       duration, do_subcloud_flx, do_new_subflx, src_choice, gqt_choice,   &
       zero_out_conv_area, tracer_check_type, use_turb_tke, use_lcl_only, do_new_pevap, plev_for, stop_at_let, &
       use_pblhttke_avg, use_hlqtsrc_avg, use_capecin_avg, reproduce_old_version, do_plev_umf, plev_umf, shallow_umf_thresh, &
       do_eis_limit, do_eis_limitn, do_lts_limit, do_lts_limitn

  !namelist parameters for UW convective plume
  real    :: rle      = 0.10   ! for critical stopping distance for entrainment
  real    :: rpen     = 5.0    ! for entrainment efficiency
  real    :: rmaxfrac = 0.15   ! maximum allowable updraft fraction
  real    :: wmin     = 0.5    ! minimum vertical velocity for computing updraft fraction
  real    :: wmax     = 50     ! maximum allowable vertical velocity
  real    :: rbuoy    = 1.0    ! for nonhydrostatic pressure effects on updraft
  real    :: rdrag    = 1.0
  real    :: frac_drs = 0.0    !
  real    :: frac_dr0 = 1.0    !
  real    :: bigc     = 0.0    ! for momentum transfer default set to 0 assuming as tracers
  real    :: auto_th0 = 0.5e-3 ! threshold for precipitation
  real    :: auto_rate= 1.e-3
  real    :: tcrit    = -60.0  ! critical temperature
  real    :: deltaqc0 = 0.5e-3
  logical :: do_pdfpcp= .false.
  logical :: do_pmadjt= .false.
  logical :: do_emmax = .false.
  logical :: do_pnqv  = .false.
  logical :: do_tten_max = .false.
  logical :: include_emf_s = .false.
  real    :: rad_crit = 14.0   ! critical droplet radius
  real    :: emfrac_max = 1.0
  integer :: mixing_assumption = 0
  integer :: mp_choice = 1
  integer :: de_choice = 0
  real    :: Nl_land   = 300.e6
  real    :: Nl_ocean  = 100.e6
  real    :: qi_thresh = 1.e-4
  real    :: r_thresh  = 10.e-6
  logical :: do_pevap = .false.
  real    :: cfrac     = 0.1
  real    :: hcevap    = 0.8
  real    :: hcevappbl = 1.0
  real    :: pblfac    = 0.0
  logical :: do_new_pblfac = .false.
  logical :: do_weffect = .false.
  logical :: do_limit_wmax =.false.
  logical :: do_limit_fdr  =.false.
  real    :: weffect    = 0.5
  real    :: peff_l     = 1.0
  real    :: peff_i     = 1.0
  real    :: scaleh0    = 1000.
  real    :: beta       = 1000.
  real    :: tten_max   = 1000.
  real    :: rkm_max    = 20.
  real    :: rkm_min    = 0.0001
  real    :: nom_ratio  = 1.

  NAMELIST / uw_plume_nml / rle, rpen, rmaxfrac, wmin, wmax, rbuoy, rdrag, frac_drs, frac_dr0, bigc, do_limit_wmax, &
       auto_th0, auto_rate, tcrit, deltaqc0, do_pdfpcp, do_pmadjt, do_emmax, do_pnqv, do_tten_max, rad_crit, emfrac_max, &
       mixing_assumption, mp_choice, de_choice, Nl_land, Nl_ocean, qi_thresh, r_thresh, do_pevap, cfrac, hcevap, rkm_max, &
       hcevappbl, pblfac, include_emf_s, do_weffect, do_new_pblfac, weffect, peff_l, peff_i, tten_max, scaleh0, beta, rkm_min, &
       nom_ratio, do_limit_fdr
  !namelist parameters for UW convective closure
  integer :: igauss   = 1      ! options for cloudbase massflux closure
                               ! 1: cin/gaussian closure, using TKE to compute CIN.
                               ! 2: cin/gaussian closure, using W* to compute CIN.
                               ! 0: cin/tke mapse-style closure;
  real    :: rkfre    = 0.05   ! vertical velocity variance as fraction of tke
  real    :: tau_sh   = 7200.  !
  real    :: wcrit_min= 0.
  real    :: mass_fact= 0.25
  logical :: do_old_cbmfmax = .true.

  NAMELIST / uw_closure_nml / igauss, rkfre, tau_sh, wcrit_min, mass_fact, do_old_cbmfmax


!========Option for deep convection=======================================
  real    :: cbmf0         = 0.0001
  real    :: rkm_dp1       = 10.
  real    :: rkm_dp2       = 1.
  real    :: crh_th_ocean  = 0.5
  real    :: crh_th_land   = 0.5
  real    :: crh_max       = 1.0001
  real    :: cape_th       = 10.
  real    :: cin_th        = 5.
  real    :: cwfn_th       = 0.
  real    :: tau_dp        = 7200.
  real    :: rpen_d        = 5.0
  integer :: mixing_assumption_d = 0
  integer :: norder      = 1
  logical :: do_ppen_d   = .true.
  logical :: do_pevap_d  = .true.
  real    :: cfrac_d     = 0.05
  real    :: hcevap_d    = 0.8
  real    :: pblfac_d    = 0.0
  real    :: hcevappbl_d = 1.0
  real    :: dcapedm_th  = 0
  real    :: dcwfndm_th  = 0
  real    :: frac_limit_d = 0.25
  real    :: lofactor_d   = 1.0
  real    :: auto_th0_d   = 1.0e-3
  real    :: tcrit_d      = -120
  real    :: peff_l_d     = 3.8e-5
  real    :: peff_i_d     = 3.8e-5
  integer :: src_choice_d = 0
  logical :: do_forcedlifting_d = .false.
  logical :: do_lod_rkm   = .false.
  logical :: do_lod_cfrac = .false.
  logical :: do_lod_tcrit = .false.
  logical :: do_lod_cape  = .false.
  logical :: do_lod_tau   = .false.
  logical :: do_lod_cush  = .false.
  logical :: do_stochastic_rkm = .false.
  logical :: do_cgust_dp  = .false.
  logical :: include_emf_d = .true.
  real    :: gustmax      = 3.  ! maximum gustiness wind (m/s)
  real    :: cpool_gust   = 10. ! constant for cool pool effect W/m2, default=10 W/m2
  logical :: do_forced_conv = .false.
  real    :: frac_rkm_pert = 1.
  real    :: tau_dp_fact = 10
  real    :: cin_fact    = 1
  real    :: wcrit_min_gust = 0.2
  integer :: cgust_choice = 0
  NAMELIST / deep_conv_nml / cbmf0, rkm_dp1, rkm_dp2, do_forced_conv, include_emf_d, &
                 crh_th_ocean, crh_th_land, do_forcedlifting_d, frac_limit_d, wcrit_min_gust, cin_fact,&
                 cape_th, cin_th, cwfn_th, tau_dp, rpen_d, mixing_assumption_d, norder, dcwfndm_th, &
                 do_ppen_d, do_pevap_d, cfrac_d, hcevap_d, pblfac_d, hcevappbl_d, lofactor_d, dcapedm_th, &
                 auto_th0_d, tcrit_d, do_lod_rkm, do_lod_cfrac, do_lod_tcrit, do_lod_cape, &
                 peff_l_d, peff_i_d, do_lod_tau, do_lod_cush, cgust_choice, tau_dp_fact, crh_max, &
                 do_stochastic_rkm, frac_rkm_pert, do_cgust_dp, gustmax, cpool_gust, src_choice_d
!========Option for deep convection=======================================

!===option for idealized forcing===================
  logical :: do_imposing_forcing = .false.
  real    :: tdt_rate = 0.0
  real    :: qdt_rate = 0.0
  real    :: pres_min = 0.0
  real    :: pres_max = 0.0
  integer :: klevel = 10
  logical :: use_klevel   = .true.
  logical :: do_imposing_rad_cooling = .false.
  real    :: cooling_rate = -1.5 !K/day
  real    :: t_thresh = 207.5    !K
  real    :: t_strato = 200.0    !K
  real    :: tau_rad  = 5.0      !day
  NAMELIST / idealized_forcing_nml / do_imposing_forcing, tdt_rate, qdt_rate, pres_min, pres_max, &
           klevel, use_klevel, do_imposing_rad_cooling, cooling_rate, t_thresh, t_strato, tau_rad
!===option for idealized forcing====================

!------------------------------------------------------------------------

  integer :: nqv, nql, nqi, nqa ,nqn
  logical :: do_qn = .false.    ! use droplet number tracer field ?

  integer :: id_tdt_uwc, id_qdt_uwc, id_udt_uwc, id_vdt_uwc, id_prec_uwc, id_snow_uwc, &
             id_tdt_uws, id_qdt_uws, id_udt_uws, id_vdt_uws, id_prec_uws, id_snow_uws, &
             id_pct_uwc, id_pcb_uwc, id_pct_uws, id_pcb_uws, id_pct_uwd, id_pcb_uwd,   &
             id_cqa_uwc, id_cql_uwc, id_cqi_uwc, id_cqn_uwc, id_cltc_uwc,              &
             id_cqa_uws, id_cql_uws, id_cqi_uws, id_cqn_uws,                           &
       id_cin_uwc, id_cbmf_uwc, id_tke_uwc, id_tkep_uwc, id_plcl_uwc, id_zlcl_uwc, id_zinv_uwc,  &
       id_cush_uws,  id_plfc_uwc, id_enth_uwc,  &
       id_qldt_uwc, id_qidt_uwc, id_qadt_uwc, id_qndt_uwc, id_qtdt_uwc, id_cmf_uwc, &
       id_qldt_uws, id_qidt_uws, id_qadt_uws, id_qndt_uws, id_qtdt_uws, id_cmf_uws, id_wuo_uws,  &
       id_fwu_uws, id_fqa_uws, id_fql_uws, id_fqi_uws, id_fqn_uws, &
       id_fer_uws,  id_fdr_uws, id_fdrs_uws,      &
       id_hlflx_uwc, id_qtflx_uwc, id_nqtflx_uwc, &
       id_cape_uwc, id_dcin_uwc, id_dcape_uwc, id_crh_uwc, id_fcrh_uws, id_fcrh_uwd, id_pblht_uwc, &
       id_ocode_uwc, id_plnb_uwc, id_wrel_uwc, id_ufrc_uwc, id_qtmp_uwc,id_gust_uwc, &
       id_tdt_pevap_uwc, id_qdt_pevap_uwc, id_xpsrc_uwc, id_xhlsrc_uwc, id_xqtsrc_uwc,&
       id_qldet_uws, id_qidet_uws, id_qadet_uws, id_qndet_uws, id_dting_uwc, &
       id_cfq_uws, id_feq_uws, id_feq_uwc, id_hmo_uwc, id_hms_uwc, id_abu_uwc, id_peo_uwc, &
       id_tten_rad_uwc, id_tdt_forc_uwc, id_qdt_forc_uwc, id_tdt_diss_uwc, &
       id_dbuodp_uws, id_dbuodp_uwd, id_dthvdp_uwc,                        &
       id_dpint0, id_hfint0, id_hfintn0, id_dpint, id_hfint, id_hfintn,    &
       id_dgz_conv_int, id_tdt_conv_int, id_qdt_conv_int, id_hdt_conv_int,               &
       id_pflx_uwc, id_lhflx_uwc, id_shflx_uwc,                            &
       id_tdt_rad_uwc, id_tdt_dyn_uwc, id_tdt_dif_uwc, id_qvdt_dyn_uwc, id_qvdt_dif_uwc, &
       id_dgz_dyn_uwc, id_ddp_dyn_int, id_hdt_tot_int, id_hdt_ver_int,                   &
       id_hdt_vadv_int, id_hdt_hadv_int, id_hdt_sum_int,                                 &
       id_tdt_rad_int, id_tdt_dyn_int, id_tdt_dif_int, id_dgz_dyn_int, id_hdp_dyn_int,   &
       id_qvdt_dyn_int, id_qvdt_dif_int, id_hdt_forc_int, id_dgz_phy_int, id_dgz_phy_int0,  &
       id_hdt_vadv_int0, id_hdt_hadv_int0, id_hdt_sum_int0, id_ddp_dyn_int0,                &
       id_tdt_rad_int0, id_tdt_dyn_int0, id_tdt_dif_int0, id_dgz_dyn_int0, id_hdp_dyn_int0, &
       id_qvdt_dyn_int0, id_qvdt_dif_int0, id_hdt_forc_int0,                                &
       id_pblht_avg, id_hlsrc_avg, id_qtsrc_avg, id_cape_avg, id_cin_avg, id_omg_avg,       &
       id_cpool_uwc, id_lts_uwc, id_eis_uwc, id_nbuo_uws, id_buo_uws,      &
       id_qtflx_up_uwc, id_qtflx_dn_uwc, id_omega_up_uwc, id_omega_dn_uwc, &
       id_omgmc_up_uwc, id_rkm_uws, id_frkm_uws, id_scale_uwc, id_scaletr_uwc


  integer, allocatable :: id_tracerdt_uwc(:), id_tracerdt_uwc_col(:), &
                          id_tracerdtwet_uwc(:), id_tracerdtwet_uwc_col(:), &
                          id_tracerdt_uwc_nc(:), id_tracerdt_uwc_col_nc(:), id_rn(:)
  integer, allocatable :: id_trevp_uwc(:), id_trevp_uwd(:)

!========Option for deep convection=======================================
  integer :: id_tdt_uwd, id_qdt_uwd, id_qtdt_uwd, id_prec_uwd, id_snow_uwd,   &
       id_qldt_uwd, id_qidt_uwd, id_qndt_uwd, id_qadt_uwd,                    &
       id_qldet_uwd, id_qidet_uwd, id_qndet_uwd, id_qadet_uwd,                &
       id_cmf_uwd, id_wuo_uwd, id_fwu_uwd, id_fqa_uwd, id_fql_uwd, id_fqi_uwd, id_fqn_uwd, &
       id_fer_uwd, id_cbmf_uwd, id_enth_uwd, &
       id_fdr_uwd, id_fdrs_uwd, id_cqa_uwd, id_cql_uwd, id_cqi_uwd, id_cqn_uwd, &
       id_hlflx_uwd, id_qtflx_uwd, id_nqtflx_uwd, id_dcin_uwd, &
       id_dcapedm_uwd, id_dcwfndm_uwd, id_ocode_uwd, id_cush_uwd,      &
       id_tdt_pevap_uwd, id_qdt_pevap_uwd, id_rkm_uwd, id_frkm_uwd, id_buo_uwd,  &
       id_nbuo_uwd, id_cfq_uwd, id_feq_uwd,               &
       id_rand_uwd, id_pwfn_uwd, id_cwfn_uwd, id_dcwfndt_dpc, &
       id_dcwfndt_fre, id_dcwfndt_pbl, id_cwfn3d_uwd, id_cape3d_uwd
!========Option for deep convection=======================================

  type(cwetdep_type), dimension(:), allocatable :: wetdep
  type(uw_params),  save  :: Uw_p
  character(len=32), dimension(:), allocatable   :: tracername
  character(len=32), dimension(:), allocatable   :: tracer_units

  logical :: use_online_aerosol
  logical :: use_sub_seasalt
  real    :: sea_salt_scale
  real    :: om_to_oc

contains

!#####################################################################
!#####################################################################

  SUBROUTINE UW_CONV_INIT(do_strat, axes, Time, kd, Nml_mp, tracers_in_uw)
    logical,         intent(in) :: do_strat
    integer,         intent(in) :: axes(4), kd
    type(time_type), intent(in) :: Time
    type(mp_nml_type), intent(in) :: Nml_mp
    logical,         intent(in) :: tracers_in_uw(:)

!---------------------------------------------------------------------
!  intent(in) variables:
!
!      tracers_in_uw
!                   logical array indicating which of the activated
!                   tracers are to be transported by UW convection
!
!-------------------------------------------------------------------

    integer   :: unit, io

    integer   :: ntracers, n, nn, ierr, logunit
    logical   :: flag
    character(len=200) :: text_in_scheme, control
    real :: frac_junk, frac_junk2

    integer, dimension(3) :: full = (/1,2,3/), half = (/1,2,4/)

    ntracers = count(tracers_in_uw)

    call uw_params_init   (Uw_p)

!   Initialize lookup tables needed for findt and exn
!   sat_vapor_pres needs to be initialized if not already done
    call sat_vapor_pres_init
    call exn_init_k (Uw_p)
    call findt_init_k (Uw_p)

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=uw_closure_nml, iostat=io)
      ierr = check_nml_error(io,'uw_closure_nml')
      read (input_nml_file, nml=uw_conv_nml, iostat=io)
      ierr = check_nml_error(io,'uw_conv_nml')
      read (input_nml_file, nml=uw_plume_nml, iostat=io)
      ierr = check_nml_error(io,'uw_plume_nml')
      read (input_nml_file, nml=deep_conv_nml, iostat=io)
      ierr = check_nml_error(io,'deep_conv_nml')
      read (input_nml_file, nml=idealized_forcing_nml, iostat=io)
      ierr = check_nml_error(io,'idealized_forcing_nml')
#else
    if( FILE_EXIST( 'input.nml' ) ) then
       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = uw_closure_nml, iostat = io, end = 10 )
          ierr = check_nml_error(io,'uw_closure_nml')
       end do
10     call close_file ( unit )

       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = uw_conv_nml, iostat = io, end = 20 )
          ierr = check_nml_error(io,'uw_conv_nml')
       end do
20     call close_file ( unit )

       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = uw_plume_nml, iostat = io, end = 30 )
          ierr = check_nml_error(io,'uw_plume_nml')
       end do
30     call close_file ( unit )

!========Option for deep convection=======================================
       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = deep_conv_nml, iostat = io, end = 40 )
          ierr = check_nml_error(io,'deep_conv_nml')
       end do
40     call close_file ( unit )
!========Option for deep convection=======================================
!========Option for idealized forcing=====================================
       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = idealized_forcing_nml, iostat = io, end = 40 )
          ierr = check_nml_error(io,'idealized_forcing_nml')
       end do
50     call close_file ( unit )
!========Option for idealized forcing=====================================
    end if
#endif

    use_online_aerosol = Nml_mp%use_online_aerosol
    use_sub_seasalt = Nml_mp%use_sub_seasalt
    sea_salt_scale = Nml_mp%sea_salt_scale
    om_to_oc = Nml_mp%om_to_oc

    call write_version_number (version, tagname)
    logunit = stdlog()
    WRITE( logunit, nml = uw_closure_nml )
    WRITE( logunit, nml = uw_conv_nml )
    WRITE( logunit, nml = uw_plume_nml )
    WRITE( logunit, nml = deep_conv_nml )
    WRITE( logunit, nml = idealized_forcing_nml )

    if ( use_online_aerosol ) call aer_ccn_act_init

    nqv = get_tracer_index ( MODEL_ATMOS, 'sphum' )
    nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
    nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
    nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
    nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
    if (nqn /= NO_TRACER) do_qn = .true.
    if (ntracers > 0) then
      allocate ( tracername   (ntracers) )
      allocate ( tracer_units (ntracers) )
      allocate ( wetdep       (ntracers) )
      nn = 1
      do n=1,size(tracers_in_uw(:))
         if (tracers_in_uw(n)) then
             call get_tracer_names (MODEL_ATMOS, n,  &
                                    name = tracername(nn), &
                                    units = tracer_units(nn))
             flag = query_method( 'wet_deposition', MODEL_ATMOS, n, &
                                  text_in_scheme, control )
             call get_wetdep_param( text_in_scheme, control, &
                                    wetdep(nn)%scheme, &
                                    wetdep(nn)%Henry_constant, &
                                    wetdep(nn)%Henry_variable, &
                                    frac_junk, frac_junk2, &
                                    wetdep(nn)%alpha_r, &
                                    wetdep(nn)%alpha_s, &
                                    wetdep(nn)%Lwetdep, &
                                    wetdep(nn)%Lgas, &
                                    wetdep(nn)%Laerosol, &
                                    wetdep(nn)%Lice, &
                                    frac_in_cloud_uw=wetdep(nn)%frac_in_cloud )
             wetdep(nn)%scheme = lowercase( wetdep(nn)%scheme )
             nn = nn + 1
          endif
       end do
    endif

    id_xpsrc_uwc  = register_diag_field (mod_name,'xpsrc_uwc', axes(1:2), Time, &
         'xpsrc', 'hPa' )
    id_xhlsrc_uwc = register_diag_field (mod_name,'xhlsrc_uwc', axes(1:2), Time, &
         'xhlsrc', 'J/kg' )
    id_xqtsrc_uwc = register_diag_field (mod_name,'xqtsrc_uwc', axes(1:2), Time, &
         'xqtsrc', 'kg/kg' )

    id_tdt_pevap_uwc = register_diag_field ( mod_name, 'tdt_pevap_uwc', axes(1:3), Time, &
         'Temperature tendency due to pevap from uw_conv', 'K/s', missing_value=mv )
    id_qdt_pevap_uwc = register_diag_field ( mod_name, 'qdt_pevap_uwc', axes(1:3), Time, &
         'Spec. humidity tendency due to pevap from uw_conv', 'kg/kg/s', missing_value=mv)

    id_tdt_uwc = register_diag_field ( mod_name, 'tdt_uwc', axes(1:3), Time, &
         'Temperature tendency from uw_conv', 'K/s', missing_value=mv )
    id_qdt_uwc = register_diag_field ( mod_name, 'qdt_uwc', axes(1:3), Time, &
         'Spec. humidity tendency from uw_conv', 'kg/kg/s', missing_value=mv)
    id_udt_uwc = register_diag_field ( mod_name, 'udt_uwc', axes(1:3), Time, &
         'U tendency from uw_conv', 'm/s2', missing_value=mv )
    id_vdt_uwc = register_diag_field ( mod_name, 'vdt_uwc', axes(1:3), Time, &
         'V tendency from uw_conv', 'm/s2', missing_value=mv)
    id_cmf_uwc = register_diag_field ( mod_name, 'cmf_uwc', axes(half), Time, &
         'Total convective mass flux from uw_conv', 'kg/m2/s', missing_value=mv)
    id_cmf_uws = register_diag_field ( mod_name, 'cmf_uws', axes(half), Time, &
         'Convective mass flux from shallow plume uw_conv', 'kg/m2/s', missing_value=mv)
    id_cfq_uws = register_diag_field ( mod_name, 'cfq_uws', axes(half), Time,   &
         'Convective frequency for shallow plume', 'none', missing_value=mv)
    id_peo_uwc = register_diag_field ( mod_name, 'peo_uwc', axes(1:3), Time,   &
         'Convective precipitation efficiency', 'none', missing_value=mv)
    id_hmo_uwc = register_diag_field ( mod_name, 'hmo_uwc', axes(1:3), Time,   &
         'moist static energy', 'J/kg', missing_value=mv)
    id_hms_uwc = register_diag_field ( mod_name, 'hms_uwc', axes(1:3), Time,   &
         'moist static energy', 'J/kg', missing_value=mv)
    id_abu_uwc = register_diag_field ( mod_name, 'abu_uwc', axes(1:3), Time,   &
         'adiabatic buoyancy', 'K', missing_value=mv)
    id_dthvdp_uwc = register_diag_field ( mod_name, 'dthvdp_uwc', axes(1:3), Time,   &
         'dthvdp_uwc', 'K/Pa', missing_value=mv)
    id_dbuodp_uws = register_diag_field ( mod_name, 'dbuodp_uws', axes(1:3), Time,   &
         'dbuodp_uws', 'm/s2', missing_value=mv)
    id_dbuodp_uwd = register_diag_field ( mod_name, 'dbuodp_uwd', axes(1:3), Time,   &
         'dbuodp_uwd', 'm/s2', missing_value=mv)
     id_buo_uws= register_diag_field ( mod_name, 'buo_uws', axes(1:3), Time,    &
            'shallow plume buoyancy', 'K', missing_value=mv)
    id_wuo_uws = register_diag_field ( mod_name, 'wuo_uws', axes(half), Time,   &
         'Updraft velocity from shallow plume', 'm/s', missing_value=mv)
    id_fwu_uws = register_diag_field ( mod_name, 'fwu_uws', axes(half), Time,   &
         'Frequency weighted updraft velocity from shallow plume', 'm/s', missing_value=mv)
    id_fqa_uws = register_diag_field ( mod_name, 'fqa_uws', axes(half), Time,   &
         'Frequency weighted updraft fraction from shallow plume', 'none', missing_value=mv)
    id_fql_uws = register_diag_field ( mod_name, 'fql_uws', axes(half), Time,   &
         'Frequency weighted updraft liquid from shallow plume', 'kg/kg', missing_value=mv)
    id_fqi_uws = register_diag_field ( mod_name, 'fqi_uws', axes(half), Time,   &
         'Frequency weighted updraft ice from shallow plume', 'kg/kg', missing_value=mv)
    id_fqn_uws = register_diag_field ( mod_name, 'fqn_uws', axes(half), Time,   &
         'Frequency weighted updraft liquid drop no from shallow plume', '#/kg', missing_value=mv)

    id_fer_uws = register_diag_field ( mod_name, 'fer_uws', axes(1:3), Time, &
         'fractional entrainment rate from shallow plume', '1/Pa', missing_value=mv)
    id_fdr_uws = register_diag_field ( mod_name, 'fdr_uws', axes(1:3), Time, &
         'fractional detrainment rate from shallow plume', '1/Pa', missing_value=mv)
    id_fdrs_uws = register_diag_field (mod_name,'fdrs_uws', axes(1:3), Time, &
         'fractional detrainment rate for saturated air from shallow plume', '1/Pa', missing_value=mv)

    id_cqa_uwc = register_diag_field ( mod_name, 'cqa_uwc', axes(1:3), Time, &
         'convective cloud area fraction from uw_conv', 'none', missing_value=mv)
    id_cql_uwc = register_diag_field ( mod_name, 'cql_uwc', axes(1:3), Time, &
         'mass fraction of convective cloud liquid water from uw_conv', 'kg/kg', missing_value=mv)
    id_cqi_uwc = register_diag_field ( mod_name, 'cqi_uwc', axes(1:3), Time, &
         'mass fraction of convective cloud ice water from uw_conv', 'kg/kg', missing_value=mv)
    id_cqn_uwc = register_diag_field ( mod_name, 'cqn_uwc', axes(1:3), Time, &
         'Updraft liquid drop number from uw_conv', '/kg', missing_value=mv)

    id_cqa_uws = register_diag_field ( mod_name, 'cqa_uws', axes(half), Time, &
         'Updraft fractional area from shallow plume', 'none', missing_value=mv)
    id_cql_uws = register_diag_field ( mod_name, 'cql_uws', axes(half), Time, &
         'Updraft liquid water mixing ratio from shallow plume', 'kg/kg', missing_value=mv)
    id_cqi_uws = register_diag_field ( mod_name, 'cqi_uws', axes(half), Time, &
         'Updraft ice water mixing ratio from shallow plume', 'kg/kg', missing_value=mv)
    id_cqn_uws = register_diag_field ( mod_name, 'cqn_uws', axes(half), Time, &
         'Updraft liquid drop number from shallow plume', '/kg', missing_value=mv)

    id_cltc_uwc = register_cmip_diag_field_2d (mod_name, 'cltc_uw', Time,  &
                                'Convective Cloud Cover Percentage', '%',  &
                             standard_name='convective_cloud_area_fraction')

    id_hlflx_uwc=register_diag_field (mod_name,'hlflx_uwc',axes(half),Time, &
         'liquid water static energy flux from uw_conv', 'W/m2', missing_value=mv)
    id_qtflx_uwc = register_diag_field (mod_name,'qtflx_uwc',axes(half),Time, &
         'total water flux from uw_conv', 'kg/m2/s', missing_value=mv)

  if(do_mse_budget) then

    id_dpint0 = register_diag_field (mod_name,'dpint0', axes(1:2), Time, &
         'column integrated mass', 'kg',                                 &
         interp_method = "conserve_order1" )
    id_hfint0 = register_diag_field (mod_name,'hfint0', axes(1:2), Time, &
         'column integrated MSE', 'J/m2',                                &
         interp_method = "conserve_order1" )
    id_hfintn0 = register_diag_field (mod_name,'hfintn0', axes(1:2), Time, &
         'column integrated MSE normalized by total column mass', 'J/kg',  &
         interp_method = "conserve_order1" )
    id_hdt_tot_int = register_diag_field (mod_name,'hdt_tot_int', axes(1:2), Time,          &
         'total vertically integrated MSE tendency based on consecutive time step', 'W/m2', &
         interp_method = "conserve_order1" )
    id_hdt_ver_int = register_diag_field (mod_name,'hdt_ver_int', axes(1:2), Time, &
         'total vertically integrated MSE budget for verfication purpose', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_hdt_sum_int0 = register_diag_field (mod_name,'hdt_sum_int0', axes(1:2), Time, &
         'sum of the column integrated MSE tendencies',   'W/m2',              &
         interp_method = "conserve_order1" )
    id_hdt_vadv_int0 = register_diag_field (mod_name,'hdt_vadv_int0', axes(1:2), Time, &
         'column integrated MSE tendency due to vertical advection',   'W/m2',   &
         interp_method = "conserve_order1" )
    id_hdt_hadv_int0 = register_diag_field (mod_name,'hdt_hadv_int0', axes(1:2), Time, &
         'column integrated MSE tendency due to horizontal advection', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_tdt_rad_int0 = register_diag_field (mod_name,'tdt_rad_int0', axes(1:2), Time, &
         'column integrated MSE tendency due to radiation', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_tdt_dyn_int0 = register_diag_field (mod_name,'tdt_dyn_int0', axes(1:2), Time, &
         'column integrated MSE tendency due to dynamical tendency for T', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_tdt_dif_int0 = register_diag_field (mod_name,'tdt_dif_int0', axes(1:2), Time, &
         'column integrated MSE tendency due to diffusion tendency for T', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_qvdt_dyn_int0 = register_diag_field (mod_name,'qvdt_dyn_int0', axes(1:2), Time, &
         'column integrated MSE tendency due to dynamical tendency for vapor', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_qvdt_dif_int0 = register_diag_field (mod_name,'qvdt_dif_int0', axes(1:2), Time, &
         'column integrated MSE tendency due to diffusion tendency for vapor', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_dgz_dyn_int0 = register_diag_field (mod_name,'dgz_dyn_int0', axes(1:2), Time, &
         'column integrated MSE tendency due to dynamical tendency for geopotential height', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_dgz_phy_int0 = register_diag_field (mod_name,'dgz_phy_int0', axes(1:2), Time, &
         'column integrated MSE tendency due to height adjustment before moist convection', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_hdp_dyn_int0 = register_diag_field (mod_name,'hdp_dyn_int0', axes(1:2), Time, &
         'column integrated MSE tendency due to mass divergence from the dynamics', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_ddp_dyn_int0 = register_diag_field (mod_name,'ddp_dyn_int0', axes(1:2), Time, &
         'column integrated mass tendency due to dynamics', 'Pa/s',   &
         interp_method = "conserve_order1" )

    id_dpint = register_diag_field (mod_name,'dpint', axes(1:3), Time, &
         'vertically integrated mass', 'kg',                           &
         interp_method = "conserve_order1" )
    id_hfint = register_diag_field (mod_name,'hfint', axes(1:3), Time, &
         'vertically integrated MSE', 'J/m2',                       &
         interp_method = "conserve_order1" )
    id_hfintn = register_diag_field (mod_name,'hfintn', axes(1:3), Time, &
         'vertically integrated MSE normalized by total column mass', 'J/kg',      &
         interp_method = "conserve_order1" )
    id_dgz_conv_int = register_diag_field (mod_name,'dgz_conv_int', axes(1:3), Time, &
         'vertically integrated geopotential height tendency due to convection', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_tdt_conv_int = register_diag_field (mod_name,'tdt_conv_int', axes(1:3), Time, &
         'vertically integrated MSE tendency due to convective tendency in T', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_qdt_conv_int = register_diag_field (mod_name,'qdt_conv_int', axes(1:3), Time, &
         'vertically integrated MSE tendency due to convective tendency in qv and qi', 'W/m2', &
         interp_method = "conserve_order1" )
    id_hdt_conv_int = register_diag_field (mod_name,'hdt_conv_int', axes(1:3), Time, &
         'vertically integrated MSE tendency due to convection', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_hdt_forc_int = register_diag_field (mod_name,'hdt_forc_int', axes(1:3), Time, &
         'vertically integrated MSE tendency due to all forcings before convection', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_hdt_vadv_int = register_diag_field (mod_name,'hdt_vadv_int', axes(1:3), Time, &
         'vertically integrated MSE tendency due to vertical advection',   'W/m2',   &
         interp_method = "conserve_order1" )
    id_hdt_hadv_int = register_diag_field (mod_name,'hdt_hadv_int', axes(1:3), Time, &
         'vertically integrated MSE tendency due to horizontal advection', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_hdt_sum_int = register_diag_field (mod_name,'hdt_sum_int', axes(1:3), Time, &
         'sum of the vertically integrated MSE tendencies',   'W/m2',              &
         interp_method = "conserve_order1" )
    id_tdt_rad_int = register_diag_field (mod_name,'tdt_rad_int', axes(1:3), Time, &
         'vertically integrated MSE tendency due to radiation', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_tdt_dyn_int = register_diag_field (mod_name,'tdt_dyn_int', axes(1:3), Time, &
         'vertically integrated MSE tendency due to dynamical tendency for T', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_tdt_dif_int = register_diag_field (mod_name,'tdt_dif_int', axes(1:3), Time, &
         'vertically integrated MSE tendency due to diffusion tendency for T', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_qvdt_dyn_int = register_diag_field (mod_name,'qvdt_dyn_int', axes(1:3), Time, &
         'vertically integrated MSE tendency due to dynamical tendency for vapor', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_qvdt_dif_int = register_diag_field (mod_name,'qvdt_dif_int', axes(1:3), Time, &
         'vertically integrated MSE tendency due to diffusion tendency for vapor', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_dgz_dyn_int = register_diag_field (mod_name,'dgz_dyn_int', axes(1:3), Time, &
         'vertically integrated MSE tendency due to dynamical tendency for geopotential height', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_dgz_phy_int = register_diag_field (mod_name,'dgz_phy_int', axes(1:3), Time, &
         'vertically integrated MSE tendency due to height adjustment before moist convection', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_hdp_dyn_int = register_diag_field (mod_name,'hdp_dyn_int', axes(1:3), Time, &
         'vertically integrated MSE tendency due to mass divergence from the dynamics', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_ddp_dyn_int = register_diag_field (mod_name,'ddp_dyn_int', axes(1:3), Time, &
         'vertically integrated mass tendency due to dynamics', 'Pa/s',   &
         interp_method = "conserve_order1" )

    id_tdt_rad_uwc = register_diag_field (mod_name,'tdt_rad_uwc',axes(1:3),Time, &
         'radiative cooling rate', 'K/s', missing_value=mv)
    id_tdt_dyn_uwc = register_diag_field (mod_name,'tdt_dyn_uwc',axes(1:3),Time, &
         'dynamics T tendency', 'K/s', missing_value=mv)
    id_tdt_dif_uwc = register_diag_field (mod_name,'tdt_dif_uwc',axes(1:3),Time, &
         'diffusion T tendency', 'K/s', missing_value=mv)
    id_qvdt_dyn_uwc = register_diag_field (mod_name,'qvdt_dyn_uwc',axes(1:3),Time, &
         'dynamics qv tendency', 'kg/kg/s', missing_value=mv)
    id_qvdt_dif_uwc = register_diag_field (mod_name,'qvdt_dif_uwc',axes(1:3),Time, &
         'diffusion qv tendency', 'kg/kg/s', missing_value=mv)
    id_dgz_dyn_uwc = register_diag_field (mod_name,'dgz_dyn_uwc',axes(1:3),Time, &
         'geopotential height tendency due to dynamics', 'm2/s2/s', missing_value=mv)

    id_lhflx_uwc = register_diag_field (mod_name,'lhflx_uwc', axes(1:2), Time, &
         'surface latent heat flux from uw_conv', 'W/m2',                      &
         interp_method = "conserve_order1" )
    id_shflx_uwc = register_diag_field (mod_name,'shflx_uwc', axes(1:2), Time, &
         'surface buoyancy flux from uw_conv', 'W/m2',                       &
         interp_method = "conserve_order1" )

    id_nqtflx_uwc = register_diag_field (mod_name,'nqtflx_uwc',axes(half),Time, &
         'net total water flux from uw_conv', 'kg/m2/s', missing_value=mv)
    id_qtflx_up_uwc = register_diag_field (mod_name,'qtflx_up_uwc',axes(half),Time, &
         'Total water flux from resolved upward flow', 'kg/m2/s', missing_value=mv)
    id_qtflx_dn_uwc = register_diag_field (mod_name,'qtflx_dn_uwc',axes(half),Time, &
         'Total water flux from resolved downward flow', 'kg/m2/s', missing_value=mv)
    id_omgmc_up_uwc = register_diag_field (mod_name,'omgmc_up_uwc',axes(half),Time, &
         'Total upward mass flux', 'kg/m2/s', missing_value=mv)
    id_omega_up_uwc = register_diag_field (mod_name,'omega_up_uwc',axes(half),Time, &
         'Total mass flux from resolved upward flow', 'kg/m2/s', missing_value=mv)
    id_omega_dn_uwc = register_diag_field (mod_name,'omega_dn_uwc',axes(half),Time, &
         'Total mass flux from resolved downward flow', 'kg/m2/s', missing_value=mv)
    id_pflx_uwc = register_diag_field (mod_name,'pflx_uwc',axes(half),Time, &
         'vertical distribution of precipitation flux', 'kg/m2/s', missing_value=mv)

  endif !do_mse_budget

    id_nbuo_uws = register_diag_field (mod_name,'nbuo_uws', axes(1:2), Time, &
         'negative buoyancy for penetrative plume', 'K', interp_method = "conserve_order1" )

    id_lts_uwc = register_diag_field (mod_name,'lts_uwc', axes(1:2), Time, &
         'low tropospheric stability', 'K', interp_method = "conserve_order1" )
    id_eis_uwc = register_diag_field (mod_name,'eis_uwc', axes(1:2), Time, &
         'estimated inversion strength', 'K', interp_method = "conserve_order1" )

    id_hlsrc_avg = register_diag_field (mod_name,'hlsrc_avg', axes(1:2), Time, &
         'hlsrc_avg', 'J/m2', interp_method = "conserve_order1" )
    id_qtsrc_avg = register_diag_field (mod_name,'qtsrc_avg', axes(1:2), Time, &
         'qtsrc_avg', 'kg/kg', interp_method = "conserve_order1" )
    id_cape_avg = register_diag_field (mod_name,'cape_avg', axes(1:2), Time, &
         'cape_avg', 'J/kg', interp_method = "conserve_order1" )
    id_cin_avg = register_diag_field (mod_name,'cin_avg', axes(1:2), Time, &
         'cin_avg', 'J/kg', interp_method = "conserve_order1" )
    id_pblht_avg = register_diag_field (mod_name,'pblht_avg', axes(1:2), Time, &
         'pblht_avg', 'm', interp_method = "conserve_order1" )
    id_omg_avg = register_diag_field (mod_name,'omg_avg', axes(1:2), Time, &
         'omg_avg', 'Pa/s', interp_method = "conserve_order1" )

    id_cpool_uwc = register_diag_field (mod_name,'cpool_uwc', axes(1:2), Time, &
         'cold pool', 'm2/s3', &
         interp_method = "conserve_order1" )

    id_prec_uwc = register_diag_field (mod_name,'prec_uwc', axes(1:2), Time, &
         'Precipitation rate from uw_conv', 'kg/m2/sec',                     &
         interp_method = "conserve_order1" )
    id_snow_uwc = register_diag_field (mod_name,'snow_uwc', axes(1:2), Time, &
         'Frozen precip. rate from uw_conv', 'kg/m2/sec',                       &
         interp_method = "conserve_order1" )
    id_prec_uws = register_diag_field (mod_name,'prec_uws', axes(1:2), Time, &
         'Precipitation rate from shallow plume', 'kg/m2/sec',                     &
         interp_method = "conserve_order1" )
    id_snow_uws = register_diag_field (mod_name,'snow_uws', axes(1:2), Time, &
         'Frozen precip. rate from shallow plume', 'kg/m2/sec',                       &
         interp_method = "conserve_order1" )

    id_pct_uwc = register_diag_field ( mod_name, 'pct_uwc', axes(1:2), Time, &
         'Cloud-top pressure from uw_conv', 'hPa' )
    id_pcb_uwc = register_diag_field ( mod_name, 'pcb_uwc', axes(1:2), Time, &
         'Cloud-base pressure from uw_conv', 'hPa' )
    id_pct_uws = register_diag_field ( mod_name, 'pct_uws', axes(1:2), Time, &
         'Cloud-top pressure from shallow plume', 'hPa' )
    id_pcb_uws = register_diag_field ( mod_name, 'pcb_uws', axes(1:2), Time, &
         'Cloud-base pressure from shallow plume', 'hPa' )

    id_cin_uwc = register_diag_field ( mod_name, 'cin_uwc', axes(1:2), Time, &
         'CIN from uw_conv', 'J/kg' )
    id_cape_uwc= register_diag_field ( mod_name,'cape_uwc', axes(1:2), Time, &
         'CAPE from uw_conv', 'J/kg' )
    id_gust_uwc= register_diag_field ( mod_name,'gust_uwc', axes(1:2), Time, &
         'gustiness from uw_conv', 'm2/s2' )
    id_crh_uwc= register_diag_field ( mod_name,'crh_uwc', axes(1:2), Time, &
         'Column RH from uw_conv', '%' )
    id_fcrh_uws= register_diag_field ( mod_name,'fcrh_uws', axes(1:2), Time, &
         'shallow frequency weighted column RH from uw_conv', '%' )
    id_fcrh_uwd= register_diag_field ( mod_name,'fcrh_uwd', axes(1:2), Time, &
         'deep frequency weighted column RH from uw_conv', '%' )
    id_pblht_uwc= register_diag_field ( mod_name,'pblht_uwc', axes(1:2), Time, &
         'PBL height from uw_conv', 'm' )
    id_cbmf_uwc = register_diag_field (mod_name,'cbmf_uwc', axes(1:2), Time, &
         'Cloud-base mass flux from uw_conv', 'kg/m2/s' )
    id_wrel_uwc = register_diag_field (mod_name,'wrel_uwc', axes(1:2), Time, &
         'Release level vertical velocity from uw_conv', 'm/s' )
    id_ufrc_uwc = register_diag_field (mod_name,'ufrc_uwc', axes(1:2), Time, &
         'Release level updraft fraction from uw_conv', 'none' )
    id_tke_uwc = register_diag_field ( mod_name, 'tke_uwc', axes(1:2), Time, &
         'PBL mean TKE from uw_conv', 'm2/s2' )
    id_tkep_uwc = register_diag_field ( mod_name, 'tkep_uwc', axes(1:2), Time, &
         'prognostic estimate of PBL mean TKE from uw_conv', 'm2/s2' )
    id_plcl_uwc = register_diag_field (mod_name,'plcl_uwc', axes(1:2), Time, &
         'LCL pressure from uw_conv', 'hPa' )
    id_zlcl_uwc = register_diag_field (mod_name,'zlcl_uwc', axes(1:2), Time, &
         'LCL depth from uw_conv', 'm' )
    id_plfc_uwc = register_diag_field (mod_name,'plfc_uwc', axes(1:2), Time, &
         'LFC pressure from uw_conv', 'hPa' )
    id_plnb_uwc = register_diag_field (mod_name,'plnb_uwc', axes(1:2), Time, &
         'LNB pressure from uw_conv', 'hPa' )
    id_zinv_uwc = register_diag_field (mod_name,'zinv_uwc', axes(1:2), Time, &
         'Inversion pressure from uw_conv', 'm' )
    id_cush_uws = register_diag_field (mod_name,'cush_uws', axes(1:2), Time, &
         'Convective scale height from uw_conv', 'm' )
    id_dcin_uwc = register_diag_field (mod_name, 'dcin_uwc', axes(1:2), Time, &
         'dCIN/cbmf from uw_conv', 'J/kg/(kg/m2/s)' )
    id_dcape_uwc= register_diag_field (mod_name, 'dcape_uwc', axes(1:2), Time, &
         'dCAPE/cbmf from uw_conv', 'J/kg/(kg/m2/s)' )
    id_enth_uwc = register_diag_field (mod_name,'enth_uwc', axes(1:2), Time, &
         'Column-integrated enthalpy tendency from uw_conv', 'W/m2' )
    id_qtmp_uwc = register_diag_field (mod_name,'qtmp_uwc', axes(1:2), Time, &
         'Column-integrated water tendency from uw_conv', 'kg/m2/s' )
    id_dting_uwc = register_diag_field (mod_name,'dting_uwc', axes(1:2), Time, &
         'Column-integrated heating rate from uw_conv', 'W/m2' )
    id_ocode_uwc = register_diag_field (mod_name,'ocode_uwc', axes(1:2), Time, &
         'Out code from uw_conv', 'none' )
    id_feq_uwc = register_diag_field (mod_name, 'feq_uwc',   axes(1:2), Time,   &
         'fraction of time convection occurs from uw_conv', 'none', missing_value=mv)
    id_feq_uws = register_cmip_diag_field_2d (mod_name, 'sci_uw', Time,  &
                    'Fraction of Time Shallow Convection Occurs', '1.0', &
                         standard_name='shallow_convection_time_fraction')
    id_rkm_uws = register_diag_field (mod_name, 'rkm_uws', axes(1:2), Time, &
            'lateral mixing rate parameter for shallow_plume', 'none' )
    id_frkm_uws = register_diag_field (mod_name, 'frkm_uws', axes(1:2), Time, &
            'frequency weighted lateral mixing rate parameter for shallow_plume', 'none' )
    id_scale_uwc = register_diag_field (mod_name, 'scale_uwc', axes(1:2), Time, &
            'scale_uwc in shallow_conv', '' )
    id_scaletr_uwc = register_diag_field (mod_name, 'scaletr_uwc', axes(1:2), Time, &
            'scaletr_uwc in shallow_conv', '' )
    if ( do_strat ) then
       id_qldt_uwc= register_diag_field (mod_name,'qldt_uwc',axes(1:3),Time, &
            'liquid water tendency from uw_conv', 'kg/kg/s', missing_value=mv)
       id_qidt_uwc= register_diag_field (mod_name,'qidt_uwc',axes(1:3),Time, &
            'ice water tendency from uw_conv', 'kg/kg/s', missing_value=mv)
       id_qadt_uwc= register_diag_field (mod_name,'qadt_uwc',axes(1:3),Time, &
            'cloud amount tendency from uw_conv', '1/s', missing_value=mv )
       id_qndt_uwc= register_diag_field (mod_name,'qndt_uwc',axes(1:3),Time, &
            'cloud droplet number fraction tendency from uw_conv', '#/kg/s', missing_value=mv )
       id_qtdt_uwc= register_diag_field (mod_name,'qtdt_uwc',axes(1:3),Time, &
            'total water tendency from uw_conv', 'kg/kg/s', missing_value=mv)
       id_qldt_uws= register_diag_field (mod_name,'qldt_uws',axes(1:3),Time, &
            'liquid water tendency from shallow plume', 'kg/kg/s', missing_value=mv)
       id_qidt_uws= register_diag_field (mod_name,'qidt_uws',axes(1:3),Time, &
            'ice water tendency from shallow plume', 'kg/kg/s', missing_value=mv)
       id_qadt_uws= register_diag_field (mod_name,'qadt_uws',axes(1:3),Time, &
            'cloud amount tendency from shallow plume', '1/s', missing_value=mv )
       id_qndt_uws= register_diag_field (mod_name,'qndt_uws',axes(1:3),Time, &
            'cloud droplet number fraction tendency from shallow plume', '#/kg/s', missing_value=mv )
       id_qtdt_uws= register_diag_field (mod_name,'qtdt_uws',axes(1:3),Time, &
            'total water tendency from shallow plume', 'kg/kg/s', missing_value=mv)
       id_qadet_uws = register_diag_field (mod_name,'qadet_uws',axes(1:3),Time, &
            'qa detrainment from shallow plume', '1/s', missing_value=mv)
       id_qldet_uws = register_diag_field (mod_name,'qldet_uws',axes(1:3),Time, &
            'ql detrainment from shallow plume', 'kg/kg/s', missing_value=mv)
       id_qidet_uws = register_diag_field (mod_name,'qidet_uws',axes(1:3),Time, &
            'qi detrainment from shallow plume', 'kg/kg/s', missing_value=mv)
       id_qndet_uws = register_diag_field (mod_name,'qndet_uws',axes(1:3),Time, &
            'qn detrainment from shallow plume', '#/kg/s', missing_value=mv)
    end if

    if (do_imposing_rad_cooling) then
       id_tten_rad_uwc = register_diag_field ( mod_name, 'tten_rad_uwc', axes(1:3), Time, &
         'Idealized radiative temperature tendency from uw_conv', 'K/s', missing_value=mv )
    end if
    if (do_imposing_forcing) then
       id_tdt_forc_uwc = register_diag_field ( mod_name, 'tdt_forc_uwc', axes(1:3), Time, &
         'Idealized temperature forcing from uw_conv', 'K/s', missing_value=mv )
       id_qdt_forc_uwc = register_diag_field ( mod_name, 'qdt_forc_uwc', axes(1:3), Time, &
         'Idealized humidity forcing from uw_conv', 'K/s', missing_value=mv )
    end if
    if (do_uwcmt) then
        id_tdt_diss_uwc = register_diag_field ( mod_name, 'tdt_diss_uwc', axes(1:3), Time, &
         'Temperature tendency due to dissipation from uw_conv', 'K/s', missing_value=mv )
    end if
!========Option for deep convection=======================================
    if (do_deep) then
       id_pct_uwd = register_diag_field (mod_name, 'pct_uwd', axes(1:2), Time, &
         'Cloud-top pressure from deep plume', 'hPa' )
       id_pcb_uwd = register_diag_field (mod_name, 'pcb_uwd', axes(1:2), Time, &
         'Cloud-base pressure from deep plume', 'hPa' )
       id_feq_uwd = register_diag_field (mod_name, 'feq_uwd', axes(1:2), Time, &
         'fraction of time deep plume occurs', 'none', missing_value=mv)

       id_tdt_pevap_uwd = register_diag_field ( mod_name, 'tdt_pevap_uwd', axes(1:3), Time, &
            'Temperature tendency due to pevap from deep plume', 'K/s', missing_value=mv )
       id_qdt_pevap_uwd = register_diag_field ( mod_name, 'qdt_pevap_uwd', axes(1:3), Time, &
            'Spec. humidity tendency due to pevap from deep plume', 'kg/kg/s', missing_value=mv)

       id_tdt_uwd = register_diag_field ( mod_name, 'tdt_uwd', axes(1:3), Time, &
            'Temperature tendency from deep plume', 'K/s', missing_value=mv )
       id_qdt_uwd = register_diag_field ( mod_name, 'qdt_uwd', axes(1:3), Time, &
            'Spec. humidity tendency from deep plume', 'kg/kg/s', missing_value=mv)
       id_qtdt_uwd= register_diag_field ( mod_name, 'qtdt_uwd', axes(1:3), Time, &
            'Total water spec. humidity tendency from deep plume', 'kg/kg/s', missing_value=mv)
       id_cmf_uwd = register_diag_field ( mod_name, 'cmf_uwd', axes(half), Time, &
            'Convective mass flux from deep plume', 'kg/m2/s', missing_value=mv)
       id_cfq_uwd = register_diag_field ( mod_name, 'cfq_uwd', axes(half), Time,   &
            'Convective frequency for deep plume', 'none', missing_value=mv)
       id_wuo_uwd = register_diag_field ( mod_name, 'wuo_uwd', axes(half), Time,   &
            'Updraft velocity from deep plume', 'm/s', missing_value=mv)
       id_fwu_uwd = register_diag_field ( mod_name, 'fwu_uwd', axes(half), Time,   &
            'Frequency weighted updraft velocity from deep plume', 'm/s', missing_value=mv)
       id_fqa_uwd = register_diag_field ( mod_name, 'fqa_uwd', axes(half), Time,   &
            'Frequency weighted updraft fraction from deep plume', 'none', missing_value=mv)
       id_fql_uwd = register_diag_field ( mod_name, 'fql_uwd', axes(half), Time,   &
            'Frequency weighted updraft liquid from deep plume', 'kg/kg', missing_value=mv)
       id_fqi_uwd = register_diag_field ( mod_name, 'fqi_uwd', axes(half), Time,   &
            'Frequency weighted updraft ice from deep plume', 'kg/kg', missing_value=mv)
       id_fqn_uwd = register_diag_field ( mod_name, 'fqn_uwd', axes(half), Time,   &
           'Frequency weighted updraft liquid drop no from deep plume', '#/kg', missing_value=mv)

       id_buo_uwd= register_diag_field ( mod_name, 'buo_uwd', axes(1:3), Time,   &
            'deep plume buoyancy', 'K', missing_value=mv)
       id_fer_uwd = register_diag_field ( mod_name, 'fer_uwd', axes(1:3), Time, &
            'fractional entrainment rate from deep plume', '1/Pa', missing_value=mv)
       id_fdr_uwd = register_diag_field ( mod_name, 'fdr_uwd', axes(1:3), Time, &
            'fractional detrainment rate from deep plume', '1/Pa', missing_value=mv)
       id_fdrs_uwd = register_diag_field (mod_name,'fdrs_uwd', axes(1:3), Time, &
            'fractional detrainment rate for saturated air from deep plume', '1/Pa', missing_value=mv)
       id_cqa_uwd = register_diag_field ( mod_name, 'cqa_uwd', axes(half), Time, &
            'Updraft fraction from deep plume', 'none', missing_value=mv)
       id_cql_uwd = register_diag_field ( mod_name, 'cql_uwd', axes(half), Time, &
            'Updraft liquid water mixing ratio from deep plume', 'kg/kg', missing_value=mv)
       id_cqi_uwd = register_diag_field ( mod_name, 'cqi_uwd', axes(half), Time, &
            'Updraft ice water mixing ratio from deep plume', 'kg/kg', missing_value=mv)
       id_cqn_uwd = register_diag_field ( mod_name, 'cqn_uwd', axes(half), Time, &
            'Updraft liquid drop number from deep plume', '/kg', missing_value=mv)
       id_hlflx_uwd=register_diag_field (mod_name,'hlflx_uwd',axes(1:3),Time, &
            'liquid water static energy flux from deep plume', 'W/m2', missing_value=mv)
       id_qtflx_uwd = register_diag_field (mod_name,'qtflx_uwd',axes(1:3),Time, &
            'total water flux from deep plume', 'W/m2', missing_value=mv)
       id_nqtflx_uwd = register_diag_field (mod_name,'nqtflx_uwd',axes(1:3),Time, &
            'net total water flux from deep plume', 'W/m2', missing_value=mv)
       id_prec_uwd = register_diag_field (mod_name,'prec_uwd', axes(1:2), Time, &
            'Precipitation rate from deep plume', 'kg/m2/sec' )
       id_snow_uwd = register_diag_field (mod_name,'snow_uwd', axes(1:2), Time, &
            'Frozen precip. rate from deep plume', 'kg/m2/sec' )
       id_cbmf_uwd = register_diag_field (mod_name,'cbmf_uwd', axes(1:2), Time, &
            'Cloud-base mass flux from deep plume', 'kg/m2/s' )
       id_cwfn_uwd = register_diag_field (mod_name,'cwfn_uwd', axes(1:2), Time, &
            'Cloud work function from deep plume', 'J/kg' )
       id_dcapedm_uwd= register_diag_field (mod_name, 'dcapedm_uwd', axes(1:2), Time, &
            'dCAPE/cbmf from deep plume', 'J/kg/(kg/m2/s)' )
       id_dcwfndm_uwd= register_diag_field (mod_name, 'dcwfndm_uwd', axes(1:2), Time, &
            'dCWFN/cbmf from deep plume', '(J/kg)/(kg/m2/s)' )
       id_cush_uwd = register_diag_field (mod_name, 'cush_uwd',  axes(1:2), Time, &
            'convective depth from deep plume', 'm' )
       id_enth_uwd = register_diag_field (mod_name,'enth_uwd', axes(1:2), Time, &
            'Column-integrated enthalpy tendency from deep plume', 'K/s' )
       id_ocode_uwd = register_diag_field (mod_name,'ocode_uwd', axes(1:2), Time, &
            'Out code from deep plume', 'none' )
       id_rkm_uwd = register_diag_field (mod_name,'rkm_uwd', axes(1:2), Time, &
            'lateral mixing rate parameter for deep plume', 'none' )
       id_frkm_uwd = register_diag_field (mod_name,'frkm_uwd', axes(1:2), Time, &
            'frequency weighted lateral mixing rate parameter for deep plume', 'none' )
       id_rand_uwd = register_diag_field (mod_name,'rand_uwd', axes(1:2), Time, &
         'rand_uwd', 'none' )
       id_nbuo_uwd = register_diag_field (mod_name,'nbuo_uwd', axes(1:2), Time, &
            'negative buoyancy for penetrative plume', 'K', interp_method = "conserve_order1" )
       if ( do_strat ) then
          id_qldt_uwd= register_diag_field (mod_name,'qldt_uwd',axes(1:3),Time, &
               'liquid water tendency from deep plume', 'kg/kg/s', missing_value=mv)
          id_qidt_uwd= register_diag_field (mod_name,'qidt_uwd',axes(1:3),Time, &
               'ice water tendency from deep plume', 'kg/kg/s', missing_value=mv)
          id_qadt_uwd= register_diag_field (mod_name,'qadt_uwd',axes(1:3),Time, &
               'cloud fraction tendency from deep plume', '1/s', missing_value=mv )
          id_qldet_uwd = register_diag_field (mod_name,'qldet_uwd',axes(1:3),Time, &
               'ql detrainment from deep plume', 'kg/kg/s', missing_value=mv)
          id_qidet_uwd = register_diag_field (mod_name,'qidet_uwd',axes(1:3),Time, &
               'qi detrainment from deep plume', 'kg/kg/s', missing_value=mv)
          id_qadet_uwd = register_diag_field (mod_name,'qadet_uwd',axes(1:3),Time, &
               'qa detrainment from deep plume', '1/s', missing_value=mv)
          id_qndet_uwd = register_diag_field (mod_name,'qndet_uwd',axes(1:3),Time, &
               'qn detrainment from deep plume', '#/kg/s', missing_value=mv)
      end if
    end if
!========Option for deep convection=======================================


    if ( ntracers>0 ) then
      allocate(id_tracerdt_uwc(ntracers), id_tracerdt_uwc_col(ntracers) )
      allocate(id_tracerdt_uwc_nc(ntracers), id_tracerdt_uwc_col_nc(ntracers))
      allocate(id_rn(ntracers))
      allocate(id_tracerdtwet_uwc(ntracers), id_tracerdtwet_uwc_col(ntracers))
       allocate(id_trevp_uwc(ntracers))
       allocate(id_trevp_uwd(ntracers))
      do nn = 1,ntracers
         id_rn(nn) = &
            register_diag_field (mod_name, trim(tracername(nn))//'_rscale', &
                                    axes(1:3), Time, &
                                  trim(tracername(nn)) //' correction', &
                                  'none', missing_value=mv)

         id_tracerdt_uwc(nn) = &
            register_diag_field (mod_name, trim(tracername(nn))//'dt_uwc', &
                                    axes(1:3), Time, &
                                  trim(tracername(nn)) //' tendency from uw_conv', &
                                  trim(tracer_units(nn))//'/s', missing_value=mv)
            id_tracerdt_uwc_col(nn) = &
              register_diag_field (mod_name, trim(tracername(nn))//'dt_uwc_col', &
                                     axes(1:2), Time, &
                                   trim(tracername(nn)) //' column tendency from uw_conv', &
                                   trim(tracer_units(nn))//'*(kg/m2)/s', missing_value=mv)
         id_tracerdt_uwc_nc(nn) = &
            register_diag_field (mod_name, trim(tracername(nn))//'dt_uwc_nc', &
                                    axes(1:3), Time, &
                                  trim(tracername(nn)) //' tendency from uw_conv before correction', &
                                  trim(tracer_units(nn))//'/s', missing_value=mv)
            id_tracerdt_uwc_col_nc(nn) = &
              register_diag_field (mod_name, trim(tracername(nn))//'dt_uwc_col_nc', &
                                     axes(1:2), Time, &
                                   trim(tracername(nn)) //' column tendency from uw_conv before correction', &
                                   trim(tracer_units(nn))//'*(kg/m2)/s', missing_value=mv)
           id_tracerdtwet_uwc(nn) = &
              register_diag_field (mod_name, trim(tracername(nn))//'dt_uwc_wet', &
                                    axes(1:3), Time, &
                                   trim(tracername(nn)) //' tendency from uw_conv wetdep', &
                                   trim(tracer_units(nn))//'/s', missing_value=mv)
            id_tracerdtwet_uwc_col(nn) = &
              register_diag_field (mod_name, trim(tracername(nn))//'dt_uwc_wet_col', &
                                   axes(1:2), Time, &
                                   trim(tracername(nn)) //' column tendency from uw_conv wetdep', &
                                   trim(tracer_units(nn))//'*(kg/m2)/s', missing_value=mv)
           id_trevp_uwc(nn) = &
            register_diag_field (mod_name, trim(tracername(nn))//'dt_evp_uwc', &
                                    axes(1:3), Time, &
                                  trim(tracername(nn)) //' tendency from uw_conv precip_revap', &
                                  trim(tracer_units(nn))//'/s', missing_value=mv)
           id_trevp_uwd(nn) = &
            register_diag_field (mod_name, trim(tracername(nn))//'dt_evp_uwd', &
                                    axes(1:3), Time, &
                                  trim(tracername(nn)) //' tendency from uw_conv precip_revap', &
                                  trim(tracer_units(nn))//'/s', missing_value=mv)
        end do
     end if

    select case (tracer_check_type)
       case(1)
          call error_mesg ('uw_conv', &
             'tracer checks: no min/max check, yes filling (WARNING! non-conservative)', NOTE)
       case(2)
          call error_mesg ('uw_conv', 'tracer checks: no min/max check, no filling', NOTE)
       case DEFAULT
          call error_mesg ('uw_conv', 'tracer checks: DEFAULT = min/max check, no filling', NOTE)
    end select

    module_is_initialized = .true.


  end SUBROUTINE UW_CONV_INIT

!#####################################################################
!#####################################################################

  subroutine uw_conv_end
    call exn_end_k
    call findt_end_k
    module_is_initialized = .FALSE.
  end subroutine uw_conv_end

!#####################################################################
!#####################################################################

  SUBROUTINE uw_conv(is, js, Time, tb, qv, ub, vb, pmid,pint,zmid,zint, & !input
       qtr, omega, delt, pblht, ustar, bstar, qstar, land,              & !input
       coldT, asol, lat, lon, cush, tkep, do_strat,                     & !input
       skip_calculation, max_available_cf,                              & !input
       tten, qvten, qlten, qiten, qaten, qnten,                         & !output
       uten, vten, rain, snow, cmf, liq_pflx,                           & !output
       ice_pflx, cldql, cldqi, cldqa, cldqn,                            & !output
       tracers, trtend, uw_wetdep, cbmfo, gusto)
!      sflx, lflx, tdt_rad, tdt_dyn, qvdt_dyn, qldt_dyn,                & !inputforMSE
!      qidt_dyn, dgz_dyn, ddp_dyn, hdt_dgz_adj, dgz_phy,                & !inputforMSE
!      tdt_tot, qvdt_dif, qldt_dif, qidt_dif, hfint0 )                  & !inputforMSE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     SHALLOW CONVECTION SCHEME
!     Described in Bretherton et. al (MWR, April 2004)
!     For info contact Ming Zhao: ming.zhao@noaa.gov
!
!     Inputs: see below
!
!     Outputs: see below
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    implicit none

    type(time_type), intent(in)  :: Time
    integer,         intent(in)  :: is, js
    real,            intent(in)  :: delt

    real, intent(in), dimension(:,:,:)   :: ub,vb !wind profile (m/s)
    real, intent(in), dimension(:,:,:)   :: zint  !height@model interfaces(m)
    real, intent(in), dimension(:,:,:)   :: pint  !pressure@model interfaces(pa)
    real, intent(in), dimension(:,:,:)   :: tb    !temperature profile (K)
    real, intent(in), dimension(:,:,:)   :: qv    !specific humidity profile (kg/kg)
    real, intent(in), dimension(:,:,:,:) :: qtr   !cloud tracers (liq_wat, ice_wat, cld_amt, liq_drp)
    real, intent(in), dimension(:,:,:)   :: pmid  !pressure@model mid-levels (pa)
    real, intent(in), dimension(:,:,:)   :: zmid  !height@model mid-levels (m)
    real, intent(in), dimension(:,:,:)   :: omega !omega (Pa/s)
    real, intent(in), dimension(:,:)     :: land  !land fraction
    real, intent(in), dimension(:,:,:)   :: max_available_cf !  largest
                                     ! realizable value for uw cld frac
                                   ! after accounting for deep cld frac
    logical,intent(in), dimension(:,:)   :: skip_calculation ! do not
                                                 ! calculate where .true.
    logical,intent(in)                   :: do_strat !logical flag
    logical,intent(in), dimension(:,:)   :: coldT    !logical flag

    real, intent(in),    dimension(:,:)  :: pblht, ustar, bstar, qstar,  lat, lon !pbl height...
    real, intent(inout), dimension(:,:)  :: cush  ! convective scale height (m)

!required when do_mse_budget=T
!    real, intent(in),  dimension(:,:,:) :: tdt_rad, tdt_dyn, qvdt_dyn, dgz_dyn, ddp_dyn, tdt_tot
!    real, intent(in),  dimension(:,:,:) :: qvdt_dif, qldt_dyn, qidt_dyn, qldt_dif, qidt_dif, dgz_phy
!    real, intent(inout), dimension(:,:) :: hdt_dgz_adj, sflx, lflx, hfint0 !column integrated total MSE (J/m2)
!required when do_mse_budget=T

    type(aerosol_type),  intent (in)     :: asol

    real, intent(out), dimension(:,:,:)  :: tten,qvten              ! T,qv tendencies
    real, intent(out), dimension(:,:,:)  :: qlten,qiten,qaten,qnten ! q tendencies
    real, intent(out), dimension(:,:,:)  :: uten,vten               ! u,v tendencies

    real, intent(out), dimension(:,:,:)  :: cldql,cldqi,cldqa, cldqn!in-updraft q
    real, intent(out), dimension(:,:,:)  :: cmf    ! mass flux at level above layer (kg/m2/s)
    real, intent(out), dimension(:,:,:)  :: liq_pflx   ! liq precipitation flux removed from a layer
    real, intent(out), dimension(:,:,:)  :: ice_pflx   ! solid precipitation flux removed from a layer
    real, intent(out), dimension(:,:)    :: rain, snow
    real, intent(inout), dimension(:,:)  :: cbmfo, gusto, tkep! cloud-base mass flux
    real, intent(in),  dimension(:,:,:,:)  :: tracers         ! env. tracers
    real, intent(out), dimension(:,:,:,:)  :: trtend          ! calculated tracer tendencies
    real, intent(out), dimension(:,:,:)  :: uw_wetdep       ! calculated wet depostion for tracers

    integer i, j, k, kl, klm, nk, naer, na, n, ksrc, kinv

    real rhos0j
    real hlsrc, thcsrc, qctsrc, tmp, tmp1, lofactor, crh_th, tvs, qvs, gust_new, gust_dis
    real zsrc, psrc, cbmf_shallow, cbmf_old, cbmf_deep, rkm_shallow, rkm_dp
    real del_crh, dcrh, dcrh0, dpsum
    real pblfact, numx
    real, dimension(size(tb,1),size(tb,2)) :: &
         plcl,       &     ! pressure of lifting condensation level (Pa)
         zlcl,       &     ! depth of lifting condensation level (m)
         plfc,       &     ! pressure of level of free convection (Pa)
         plnb,       &     ! pressure of level of neutral buoyancy (Pa)
         cino,       &     ! cin (J/kg)
         capeo,      &     ! cape(J/kg)
         tkeo,       &     ! tke (m2/s2)
         wrelo,      &     ! release level vertical velocity (m/s)
         ufrco,      &     ! cloud-base updraft fraction
         zinvo,      &     ! surface driven mixed-layer height
         einso,      &     ! estimated inversion strength (K)
         denth,      &
         dqtmp,      &
         dting,      &
         dcino,      &
         dcapeo,     &
         cwfno,      &
         ocode,      &
         xpsrc,      &
         xhlsrc,     &
         xqtsrc,     &
         crho,       &
         rkm_s,      &
         rkm_d,      &
         cush_s,     &
         cush_d,     &
         feq_c,      &
         feq_s,      &
         feq_d,      &
         rhos,       &
         lhflx,      &
         shflx,      &
         cltc,       &
         hfint0_old,  pblht_avg, omg_avg, hlsrc_avg, qtsrc_avg, cape_avg, cin_avg, &
         hten, lts, eis, hdt_tot_int, hdt_ver_int,                                 &
         cpool, bflux, nbuo_s, nbuo_d, pcb_s, pcb_d, pcb_c, pct_s, pct_d, pct_c, omgavg

!commented out when do_mse_buget=T
    real, dimension(size(tb,1),size(tb,2),size(tb,3))   :: tdt_rad, tdt_dyn, qvdt_dyn, dgz_dyn, ddp_dyn, tdt_tot
    real, dimension(size(tb,1),size(tb,2),size(tb,3))   :: qvdt_dif, qldt_dyn, qidt_dyn, qldt_dif, qidt_dif, dgz_phy
    real, dimension(size(tb,1),size(tb,2))              :: hdt_dgz_adj, sflx, lflx, hfint0
!commented out when do_mse_buget=T

    real, dimension(size(tb,1),size(tb,2),size(tb,3)+1) :: hlflx, qtflx, pflx, qtflx_up, qtflx_dn
    real, dimension(size(tb,1),size(tb,2),size(tb,3)+1) :: omgmc_up, nqtflx, omega_up, omega_dn
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: fero_s,fdro_s,fdrso_s, tten_pevap, qvten_pevap, temp
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: qldet_s, qidet_s, qadet_s, qndet_s, peo, hmo, hms, abu, buo_s
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: qldet_d, qidet_d, qadet_d, qndet_d, dbuodp_s, dbuodp_d

    real, dimension(size(tb,1),size(tb,2),size(tb,3)+1) :: cfq_s, cqa_s, cql_s, cqi_s, cqn_s, wuo_s, cmf_s
    real, dimension(size(tb,1),size(tb,2),size(tb,3)+1) :: cfq_d, cqa_d, cql_d, cqi_d, cqn_d, wuo_d, cmf_d
    real, dimension(size(tb,3)+1)                       :: cqa, cql, cqi, cqn

    real, dimension(size(tb,1),size(tb,2),size(tb,3))   :: hdt_vadv_int, hdt_hadv_int, hdt_forc_int,           &
                                                           hdt_adv_int, hdp_dyn_int, hdt_dyn_int, hdt_sum_int, &
                                                           tdt_rad_int, tdt_dyn_int, tdt_dif_int, dgz_dyn_int, &
                                                           qvdt_dyn_int,qvdt_dif_int,qidt_dyn_int,qidt_dif_int,&
                                                           dgz_phy_int, dpint, hfint, hfintn, tdt_conv_int,    &
                                                           qdt_conv_int, hdt_conv_int, dgz_conv_int, ddp_dyn_int
    real, dimension(size(tb,1),size(tb,2),size(tb,3))   :: tv, tvn, dz, dz_n, zmid_n, dthvdp
    real, dimension(size(tb,1),size(tb,2),size(tb,3)+1) :: zint_n

    real, dimension(size(tb,1),size(tb,2))            :: scale_uw, scale_tr
    real :: tnew, qtin, dqt, temp_1, temp_max, temp_min

    !f1p
    real, dimension(size(tracers,1), size(tracers,2), size(tracers,3), size(tracers,4)) :: trtend_nc, rn_diag

!========Option for deep convection=======================================
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: uten_d, vten_d, tten_d,    &
         qvten_d, qlten_d, qiten_d, qaten_d, qnten_d, buo_d, qtten_d,               &
         fero_d, fdro_d, fdrso_d, tten_pevap_d, qvten_pevap_d
    real, dimension(size(tb,1),size(tb,2),size(tb,3)+1) :: hlflx_d, qtflx_d, pflx_d, nqtflx_d
    real, dimension(size(tb,1),size(tb,2)) :: rain_d, snow_d, cwfn_d
    real, dimension(size(tb,1),size(tb,2)) :: dcapedm_d, dcwfndm_d, denth_d, dting_d, dqtmp_d, cbmf_d
    real, dimension(size(tracers,1),size(tracers,2),size(tracers,3),size(tracers,4)) :: trevp_d, trevp_s
    real, dimension(size(tracers,3),size(tracers,4)) :: trtend_t, trwet_t
!f1p
    real, dimension(size(tracers,3),size(tracers,4)) :: trtend_t_nc, trwet_t_nc, rn
!
    type(randomNumberStream), dimension(size(tb,1), size(tb,2)) :: streams
    real, dimension(size(tb,1),size(tb,2)) :: rand_s, rand_d
    integer :: iseed
    real    :: seedwts(8) = (/3000.,1000.,300.,100.,30.,10.,3.,1./)
    real    :: tempseed(8)
    integer :: thisseed(8)
    integer :: seedperm = 0
 !========Option for deep convection=======================================

    real, dimension(size(tb,3)) :: am1, am2, am3, am4, am5, qntmp
    real, dimension(size(tb,3),5) :: amx

    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: pmass    ! layer mass (kg/m2)
    real, dimension(size(tb,1),size(tb,2))            :: tempdiag ! temporary diagnostic variable
    real, dimension(size(tracers,1),size(tracers,2),size(tracers,3),size(tracers,4))  :: trwet
    ! calculated tracer wet deposition tendencies

    integer imax, jmax, kmax
    integer kd, ntracers
    integer ktop_tmp, kbot_tmp
    real :: tten_intg, qvten_intg, cp_inv, half_delt
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: tten_forc, qten_forc
    real, dimension(size(tb,3)) :: tten_tmp, qvten_tmp !f1p
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: tten_rad
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: dissipative_heat
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: tdt_dif
    real  :: latx, lonx, lat1b, lat1e, lon1b, lon1e, lat2b, lat2e, lon2b, lon2e

    logical used
    type(sounding)          :: sd, sd1
    type(adicloud)          :: ac, ac1
    type(cclosure)          :: cc, cc1
    type(cplume)            :: cp, cp1
    type(ctend)             :: ct, ct1
    type(cpnlist)           :: cpn,dpn
    type(deepc)             :: dpc
    integer ::  ier
    character(len=256) :: ermesg

    kd = size(tracers,3)
    ntracers = size(tracers,4)
    call sd_init_k(kd,ntracers,sd);
    call sd_init_k(kd,ntracers,sd1);
    call ac_init_k(kd,ac);
    call ac_init_k(kd,ac1);
    call cp_init_k(kd,ntracers,cp)
    call cp_init_k(kd,ntracers,cp1)
    call ct_init_k(kd,ntracers,ct)
    call ct_init_k(kd,ntracers,ct1)
    !pack namelist parameters into plume and closure structure
    cpn % do_umf_pbl = do_umf_pbl
    cpn % do_qctflx_zero = do_qctflx_zero
    cpn % do_hlflx_zero  = do_hlflx_zero
    cpn % do_new_qnact   = do_new_qnact
    cpn % do_2nd_act     = do_2nd_act
    cpn % do_downdraft   = do_downdraft
    cpn % do_conv_micro_N= do_conv_micro_N
    cpn % do_varying_rpen= do_varying_rpen
    cpn % rpen_choice    = rpen_choice
    cpn % do_subcloud_flx= do_subcloud_flx
    cpn % do_new_subflx  = do_new_subflx
    cpn % do_detran_zero = do_detran_zero
    cpn % N0        = N0
    cpn % rle       = rle
    cpn % rpen      = rpen
    cpn % eis_max   = eis_max
    cpn % eis_min   = eis_min
    cpn % rmaxfrac  = rmaxfrac
    cpn % wmin      = wmin
    cpn % wmax      = wmax
    cpn % rbuoy     = rbuoy
    cpn % rdrag     = rdrag
    cpn % frac_drs  = frac_drs
    cpn % frac_dr0  = frac_dr0
    cpn % bigc      = bigc
    cpn % auto_th0  = auto_th0
    cpn % deltaqc0  = deltaqc0
    cpn % do_pdfpcp = do_pdfpcp
    cpn % do_pmadjt = do_pmadjt
    cpn % do_emmax  = do_emmax
    cpn % do_pnqv   = do_pnqv
    cpn % do_tten_max   = do_tten_max
    cpn % tten_max  = tten_max
    cpn % rkm_max   = rkm_max
    cpn % rkm_min   = rkm_min
    cpn % scaleh0   = scaleh0
    cpn % beta      = beta
    cpn % emfrac_max= emfrac_max
    cpn % auto_rate = auto_rate
    cpn % tcrit     = tcrit
    cpn % cldhgt_max= cldhgt_max
    cpn % cldhgt_max_shallow= cldhgt_max_shallow
    cpn % do_ice    = do_ice
    cpn % do_ppen   = do_ppen
    cpn % do_pevap  = do_pevap
    cpn % hcevap    = hcevap
    cpn % hcevappbl = hcevappbl
    cpn % cfrac     = cfrac
    cpn % pblfac    = pblfac
    cpn % mixing_assumption= mixing_assumption
    cpn % mp_choice = mp_choice
    cpn % de_choice = de_choice
    cpn % Nl_land   = Nl_land
    cpn % Nl_ocean  = Nl_ocean
    cpn % qi_thresh = qi_thresh
    cpn % r_thresh  = r_thresh
    cpn % peff_l    = peff_l
    cpn % peff_i    = peff_i
    cpn % do_forcedlifting= do_forcedlifting
    cpn % do_minmse = do_minmse
    cpn % atopevap  = atopevap
    cpn % wtwmin_ratio = wmin_ratio*wmin_ratio
    cpn % do_auto_aero = do_auto_aero
    cpn % rad_crit = rad_crit
    cpn % wrel_min = wrel_min
    cpn % do_weffect = do_weffect
    cpn % do_new_pblfac = do_new_pblfac
    cpn % weffect    = weffect
    cpn % use_online_aerosol = use_online_aerosol
    cpn % use_new_let = use_new_let
    cpn % use_lcl_only= use_lcl_only
    cpn % do_new_pevap= do_new_pevap
    cpn % stop_at_let = stop_at_let
    cpn % do_limit_wmax= do_limit_wmax
    cpn % do_limit_fdr = do_limit_fdr
    cpn % plev_for = plev_for
    cpn % plev_umf = plev_umf
    cpn % nom_ratio= nom_ratio
    if (ntracers > 0) then
      allocate ( cpn%tracername   (ntracers) )
      allocate ( cpn%tracer_units (ntracers) )
      allocate ( cpn%wetdep       (ntracers) )
      cpn%tracername(:) = tracername(:)
      cpn%tracer_units(:) = tracer_units(:)
      cpn%wetdep(:)%scheme = wetdep(:)%scheme
      cpn%wetdep(:)%Henry_constant = wetdep(:)%Henry_constant
      cpn%wetdep(:)%Henry_variable = wetdep(:)%Henry_variable
      cpn%wetdep(:)%frac_in_cloud = wetdep(:)%frac_in_cloud
      cpn%wetdep(:)%alpha_r = wetdep(:)%alpha_r
      cpn%wetdep(:)%alpha_s = wetdep(:)%alpha_s
      cpn%wetdep(:)%Lwetdep = wetdep(:)%Lwetdep
      cpn%wetdep(:)%Lgas = wetdep(:)%Lgas
      cpn%wetdep(:)%Laerosol = wetdep(:)%Laerosol
      cpn%wetdep(:)%Lice = wetdep(:)%Lice
      allocate ( dpn%tracername   (ntracers) )
      allocate ( dpn%tracer_units (ntracers) )
      allocate ( dpn%wetdep       (ntracers) )
      dpn%tracername(:) = tracername(:)
      dpn%tracer_units(:) = tracer_units(:)
      dpn%wetdep(:)%scheme = wetdep(:)%scheme
      dpn%wetdep(:)%Henry_constant = wetdep(:)%Henry_constant
      dpn%wetdep(:)%Henry_variable = wetdep(:)%Henry_variable
      dpn%wetdep(:)%frac_in_cloud = wetdep(:)%frac_in_cloud
      dpn%wetdep(:)%alpha_r = wetdep(:)%alpha_r
      dpn%wetdep(:)%alpha_s = wetdep(:)%alpha_s
      dpn%wetdep(:)%Lwetdep = wetdep(:)%Lwetdep
      dpn%wetdep(:)%Lgas = wetdep(:)%Lgas
      dpn%wetdep(:)%Laerosol = wetdep(:)%Laerosol
      dpn%wetdep(:)%Lice = wetdep(:)%Lice
    endif
    call cpn_copy(cpn, dpn)

    cc  % igauss    = igauss
    cc  % rkfre     = rkfre
    cc  % rmaxfrac  = rmaxfrac
    cc  % wcrit_min = wcrit_min
    cc  % mass_fact = mass_fact
    cc  % rbuoy     = rbuoy
    cc  % tau_sh    = tau_sh
    cc  % do_old_cbmfmax = do_old_cbmfmax
!========Option for deep convection=======================================
    dpc % cbmf0               = cbmf0
    dpc % rkm_dp1             = rkm_dp1
    dpc % rkm_dp2             = rkm_dp2
    dpc % crh_th_ocean        = crh_th_ocean
    dpc % crh_th_land         = crh_th_land
    dpc % cape_th             = cape_th
    dpc % cin_th              = cin_th
    dpc % cwfn_th             = cwfn_th
    dpc % tau_dp              = tau_dp
    dpc % mixing_assumption_d = mixing_assumption_d
    dpc % do_ppen_d           = do_ppen_d
    dpc % rpen_d              = rpen_d
    dpc % do_pevap_d          = do_pevap_d
    dpc % cfrac_d             = cfrac_d
    dpc % hcevap_d            = hcevap_d
    dpc % hcevappbl_d         = hcevappbl_d
    dpc % pblfac_d            = pblfac_d
    dpc % frac_limit_d        = frac_limit_d
    dpc % dcapedm_th          = dcapedm_th
    dpc % dcwfndm_th          = dcwfndm_th
    dpc % do_forcedlifting_d  = do_forcedlifting_d
    dpc % lofactor_d          = lofactor_d
    dpc % auto_th0_d          = auto_th0_d
    dpc % tcrit_d             = tcrit_d
    dpc % peff_l_d            = peff_l_d
    dpc % peff_i_d            = peff_i_d
    dpc % do_cgust_dp         = do_cgust_dp
    dpc % cgust_choice        = cgust_choice
    dpc % tau_dp_fact         = tau_dp_fact
    dpc % cin_fact            = cin_fact
    dpc % wcrit_min_gust      = wcrit_min_gust
    dpc % src_choice_d        = src_choice_d
    numx                      = duration/delt
!========Option for deep convection=======================================
    imax  = size( tb, 1 )
    jmax  = size( tb, 2 )
    kmax  = size( tb, 3 )
    sd % kmax=kmax
    sd % plev_cin=plev_cin
    sd % plev_omg=plev_omg

    kl=kmax-1
    klm=kl-1

   !initialize 3D variables outside the loop

!commented out when do_mse_budget=T
    tdt_rad=0.; tdt_dyn=0.; qvdt_dyn=0.; dgz_dyn=0.; ddp_dyn=0.;
    tdt_tot=0.; qvdt_dif=0.; qldt_dyn=0.; qldt_dif=0.;
    qidt_dyn=0.; qidt_dif=0.; hdt_dgz_adj=0.; sflx=0.; lflx=0.; dgz_phy=0.; hfint0=0.;
!commented out when do_mse_budget=T

    tten=0.; qvten=0.; qlten=0.; qiten=0.; qaten=0.; qnten=0.;
    uten=0.; vten =0.; rain =0.; snow =0.; plcl =0.; plfc=0.; plnb=0.;
    cldqa=0.; cldql=0.; cldqi=0.; cldqn=0.; zlcl=0.; cltc=0.;

    cqa  =0.; cql  =0.; cqi  =0.; cqn  =0.;
    cqa_s=0.; cql_s=0.; cqi_s=0.; cqn_s=0.;
    hlflx=0.; qtflx=0.; nqtflx=0.; pflx=0.; am1=0.; am2=0.; am3=0.; am4=0.; amx=0.;
    tten_pevap=0.; qvten_pevap=0.; temp=0.;
    ice_pflx = 0. ; liq_pflx = 0.; qtflx_up=0.; qtflx_dn=0.;
    omega_up=0.; omega_dn=0.; omgmc_up=0.;

    hdt_ver_int=0.; hdt_tot_int=0.;
    hdt_vadv_int=0.; hdt_hadv_int=0.; hdt_adv_int=0.; hdp_dyn_int=0.; hdt_dyn_int=0.; hdt_sum_int=0.;
    tdt_rad_int =0.; tdt_dyn_int =0.; tdt_dif_int=0.; dgz_dyn_int=0.; dgz_phy_int=0.; hdt_forc_int=0.;
    qvdt_dyn_int=0.; qidt_dyn_int=0.; qvdt_dif_int=0.; qidt_dif_int=0.; ddp_dyn_int=0.;

    tdt_dif=0.;

    cpool=0.; bflux=0.; scale_uw=1.; scale_tr=1.;

    cino=0.; capeo=0.; tkeo=0.; wrelo=0.; ufrco=0.; zinvo=0.; einso=0.; wuo_s=0.; peo=0.;
    fero_s=0.; fdro_s=0.; fdrso_s=0.; cmf=0.; denth=0.;  dqtmp=0.; ocode=0; cmf_s=0.; cfq_s=0.;
    dcapeo=0.; dcino=0.; xpsrc=0.; xhlsrc=0.; xqtsrc=0.; feq_s=0.; feq_d=0.; feq_c=0; rkm_s=0.;
    trtend=0.; trwet=0.; crho=0.; hmo=0.; hms=0.; abu=0.; dbuodp_s=0.; dbuodp_d=0.;
    pblht_avg=0.; omg_avg=0.; hlsrc_avg=0.; qtsrc_avg=0.; cape_avg=0.; cin_avg=0.;
    qldet_s=0.; qidet_s=0.; qadet_s=0.; qndet_s=0.;
    qldet_d=0.; qidet_d=0.; qadet_d=0.; qndet_d=0.;
    dting = 0.; cush_s=-1.;
    pcb_s=0.; pcb_d=0.; pcb_c=0.;
    pct_s=0.; pct_d=0.; pct_c=0.;
    dissipative_heat = 0.; rhos=0; lhflx=0; shflx=0;
    hfint0_old=hfint0; hfint0=0; dpint=0.; hfint=0.; hfintn=0.;
    tdt_conv_int=0.; qdt_conv_int=0.; hdt_conv_int=0.;
    lts=0.; eis=0.; nbuo_s=0.; nbuo_d=0.; lofactor=1.; omgavg=0.;
    naer = size(asol%aerosol,4)

!========Option for deep convection=======================================
    tten_d=0.; qvten_d=0.; qlten_d=0.; qiten_d=0.; qaten_d=0.; qnten_d=0.;
    uten_d=0.; vten_d =0.; rain_d =0.; snow_d =0.; qtten_d=0.; cfq_d=0.;
    trevp_d=0.; trevp_s=0.; cush_d=-1.;
    cqa_d=0.; cql_d=0.; cqi_d=0.; cqn_d=0.;
    hlflx_d=0.; qtflx_d=0.; nqtflx_d=0.; pflx_d=0.;
    wuo_d=0.; fero_d=0.; fdro_d=0.; fdrso_d=0.;
    cmf_d=0.; buo_d=0.; buo_s=0.;
    denth_d=0.; dting_d=0.; dqtmp_d=0.; cbmf_d=0.; cwfn_d=0.;
    dcapedm_d=0.; dcwfndm_d=0.;
    dcino=0.;
    tten_pevap_d=0.; qvten_pevap_d=0.; rkm_d=0.;

!========Option for deep convection=======================================
    if (do_stochastic_rkm) then
      do j = 1, jmax
        do i=1, imax
          tempseed = tb(i,j,kmax)*seedwts
          thisseed = nint(tempseed)
          streams(i,j) = initializeRandomNumberStream(ishftc(thisseed, seedperm))
          call getRandomNumbers( streams(i,j), rand_s(i,j) )
        enddo
      enddo
      do j = 1, jmax
        do i=1, imax
          tempseed = tb(i,j,1)*seedwts
          thisseed = nint(tempseed)
          streams(i,j) = initializeRandomNumberStream(ishftc(thisseed, seedperm))
          call getRandomNumbers( streams(i,j), rand_d(i,j) )
        enddo
      enddo
    end if

    do j = 1, jmax
       do i=1, imax

         trtend_t=0.; trwet_t=0.;
         do k=1,kmax
           pmass(i,j,k) = (pint(i,j,k+1) - pint(i,j,k))/GRAV
         enddo
    !relaxation TKE back to 0 with time-scale of disscale
    !tkeavg = ustar(i,j)*bstar(i,j)*disscale
    !dissipate tke with length-scale of disscale
    !tkeavg=(ustar(i,j)*bstar(i,j)*disscale)**(2./3.)
    !below following Holtslag and Boville 1993

         if (pblht(i,j).lt.0.) then
            temp_1=0.0
         elseif (pblht(i,j).gt.5000.) then
            temp_1=5000.
         else
            temp_1=pblht(i,j)
         endif

!         bflux(i,j) = 0.5*(0.6*ustar(i,j)*bstar(i,j)*temp_1)**(2./3.)
         temp_1=ustar(i,j)**3.+0.6*ustar(i,j)*bstar(i,j)*temp_1
         if (temp_1 .gt. 0.) temp_1 = 0.5*temp_1**(2./3.)
         tkeo(i,j) = MAX (tkemin, temp_1)

         cbmf_shallow=0. ! Set cbmf_shallow to avoid usage before assignment.
         if (skip_calculation(i,j)) then
           ocode(i,j) = 6
           go to 100
         endif
         call clearit(ac, cc, cp, ct, cp1, ct1);

! restrict grid-box area available to shallow convection to that which
! is not involved with deep convection; added for running with Donner deep
          cp%maxcldfrac = minval(max_available_cf(i,j,:))
          cc%maxcldfrac = cp%maxcldfrac

          cc%scaleh = cush(i,j);
          cush(i,j) = -1.;
          if(cc%scaleh.le.0.0) cc%scaleh=1000.

          am1(:) = 0.; am2(:) = 0.; am3(:) = 0.; am4(:) = 0.; am5(:) = 0.;
          amx(:,:)=0. !amx contains mixing ratio (kg/kg)

          do k=1,kmax
            tmp=1. / (zint(i,j,k)-zint(i,j,k+1)) * 1.0e9 * 1.0e-12
            tmp1=1./pmass(i,j,k)
            if(use_online_aerosol) then
              do na = 1,naer
                if(asol%aerosol_names(na) == 'so4' .or. &
                   asol%aerosol_names(na) == 'so4_anthro' .or. &
                   asol%aerosol_names(na) == 'so4_natural') then     !aerosol unit: kg/m2
                           am1(k)=am1(k)+asol%aerosol(i,j,k,na)*tmp  !am1 unit: g/cm3
                           amx(k,1)=amx(k,1)+asol%aerosol(i,j,k,na)*tmp1  !am1 unit: kg/kg
                else if(asol%aerosol_names(na) == 'omphilic' .or. &
                        asol%aerosol_names(na) == 'omphobic') then
                           am4(k)=am4(k)+asol%aerosol(i,j,k,na)*tmp
                           amx(k,4)=amx(k,4)+asol%aerosol(i,j,k,na)*tmp1
                else if(asol%aerosol_names(na) == 'bcphilic' .or. &
                        asol%aerosol_names(na) == 'bcphobic' .or. &
                        asol%aerosol_names(na) == 'dust1' .or. &
                        asol%aerosol_names(na) == 'dust2' .or. &
                        asol%aerosol_names(na) == 'dust3' .or. &
                        asol%aerosol_names(na) == 'dust_mode1_of_2') then   !h1g, 2015-09-19
                           am2(k)=am2(k)+asol%aerosol(i,j,k,na)*tmp
                           amx(k,2)=amx(k,2)+asol%aerosol(i,j,k,na)*tmp1
                else if(asol%aerosol_names(na) == 'seasalt1' .or. &
                        asol%aerosol_names(na) == 'seasalt2' .or. &
                        asol%aerosol_names(na) == 'seasalt_aitken' .or. &   !h1g, 2015-09-19
                        asol%aerosol_names(na) == 'seasalt_fine'  ) then    !h1g, 2015-09-19
                           am3(k)=am3(k)+asol%aerosol(i,j,k,na)*tmp
                           amx(k,3)=amx(k,3)+asol%aerosol(i,j,k,na)*tmp1
                else if(asol%aerosol_names(na) == 'seasalt3' .or. &
                        asol%aerosol_names(na) == 'seasalt4' .or. &
                        asol%aerosol_names(na) == 'seasalt5' .or. &
                        asol%aerosol_names(na) == 'seasalt_coarse') then    !h1g, 2015-09-19
                           am5(k)=am5(k)+asol%aerosol(i,j,k,na)*tmp
                           amx(k,5)=amx(k,5)+asol%aerosol(i,j,k,na)*tmp1
                end if
              end do
              am2(k)=am2(k)+am3(k)+am4(k)
              amx(k,2)=amx(k,2)+amx(k,3)+amx(k,4)
              if(.not. use_sub_seasalt) then
                am3(k)=am3(k)+am5(k)
                amx(k,3)=amx(k,3)+amx(k,5)
              end if
            else
              am1(k)= asol%aerosol(i,j,k,2)*tmp
              am2(k)= asol%aerosol(i,j,k,1)*tmp
              am3(k)= sea_salt_scale*asol%aerosol(i,j,k,5)*tmp
              am4(k)= om_to_oc*asol%aerosol(i,j,k,3)*tmp
              amx(k,1)= asol%aerosol(i,j,k,2)*tmp1
              amx(k,2)= asol%aerosol(i,j,k,1)*tmp1
              amx(k,3)= sea_salt_scale*asol%aerosol(i,j,k,5)*tmp1
              amx(k,4)= om_to_oc*asol%aerosol(i,j,k,3)*tmp1
            endif
          end do

!========Pack column properties into a sounding structure====================

          if (do_qn) then
             qntmp(:)=qtr(i,j,:,nqn)
          else
             qntmp(:)=0.
          end if
          ksrc=1
          tdt_dif(i,j,:)=tdt_tot(i,j,:)-tdt_rad(i,j,:);
          call pack_sd_k(land(i,j), coldT(i,j), delt, pmid(i,j,:), pint(i,j,:),    &
          zmid(i,j,:), zint(i,j,:), ub(i,j,:), vb(i,j,:), omega(i,j,:), tb(i,j,:), &
          qv(i,j,:), qtr(i,j,:,nql), qtr(i,j,:,nqi), qtr(i,j,:,nqa), qntmp,        &
          am1(:), am2(:), am3(:), am4(:), amx(:,:), tracers(i,j,:,:), src_choice,  &
          tdt_rad(i,j,:), tdt_dyn(i,j,:), qvdt_dyn(i,j,:), qidt_dyn(i,j,:),        &
          dgz_dyn(i,j,:), ddp_dyn(i,j,:), tdt_dif(i,j,:), dgz_phy(i,j,:),          &
          qvdt_dif(i,j,:), qidt_dif(i,j,:), sd, Uw_p)

!========Finite volume intepolation==========================================

          gusto(i,j)=max(gusto(i,j),tkemin)
          sd%do_mse_budget = do_mse_budget
          sd%do_gust_qt = do_gust_qt
          sd%gqt_choice = gqt_choice
          sd%cgust     = gusto(i,j)
          sd%cgust0    = cgust0
          sd%cgust_max = cgust_max
          sd%sigma0    = sigma0
          sd%tke       = tkeo(i,j)
          sd%lat       = lat(i,j)*180/3.1415926
          sd%lon       = lon(i,j)*180/3.1415926

          sd%numx      = numx

          if (use_turb_tke ) sd%tke = tkep(i,j)   !h1g, 2015-08-11

          call extend_sd_k(sd, pblht(i,j), do_ice, Uw_p)

          xpsrc (i,j)= sd%psrc
          xhlsrc(i,j)= sd%hlsrc
          xqtsrc(i,j)= sd%qctsrc
          crho  (i,j)= sd%crh
          lts   (i,j)= sd%lts
          eis   (i,j)= sd%eis
          zinvo(i,j) = sd%zinv
          kinv       = sd%kinv
          einso(i,j) = (sd%hl(kinv+1)-sd%hl(kinv))/(sd%z(kinv+1)-sd%z(kinv))/cp_air

        if (do_mse_budget) then
          tmp=sd%thvbot(1)*sd%exners(1)
          rhos(i,j)= sd%ps(0)/(rdgas*tmp)
!         lhflx(i,j)=rhos(i,j)*ustar(i,j)*qstar(i,j)
!         qvs=sd%qv(1)+sd%ssqct(1)*(sd%ps(0)-sd%p(1))
!         tvs=tmp*(1+0.608*qvs)
!         shflx(i,j)=Uw_p%cp_air*(rhos(i,j)*ustar(i,j)*bstar(i,j)*tvs/Uw_p%grav-0.608*tmp*lhflx(i,j))&
!                       / (1+0.608*qvs);
!         lhflx(i,j)=lhflx(i,j)*Uw_p%hlv

          shflx(i,j) = sflx(i,j)
          lhflx(i,j) = lflx(i,j)*Uw_p%hlv
          hfint0(i,j) =sd%hfint (1);

          do k = 1,kmax
             nk = kmax+1-k
             dpint       (i,j,nk) = sd%dpint   (k);
             hfint       (i,j,nk) = sd%hfint   (k);
             hfintn      (i,j,nk) = sd%hfintn  (k);
             tdt_rad_int (i,j,nk) = sd%tdt_rad (k);
             tdt_dyn_int (i,j,nk) = sd%tdt_dyn (k);
             tdt_dif_int (i,j,nk) = sd%tdt_dif (k);
             qvdt_dyn_int(i,j,nk) = sd%qvdt_dyn(k);
             qidt_dyn_int(i,j,nk) = sd%qidt_dyn(k);
             qvdt_dif_int(i,j,nk) = sd%qvdt_dif(k);
             dgz_dyn_int (i,j,nk) = sd%dgz_dyn (k);
             dgz_phy_int (i,j,nk) = sd%dgz_phy (k);
             hdp_dyn_int (i,j,nk) = sd%hdp_dyn (k);
             ddp_dyn_int (i,j,nk) = sd%ddp_dyn (k);
             hdt_vadv_int(i,j,nk) = sd%hdt_vadv(k);
             hdt_forc_int(i,j,nk) = sd%hdt_forc(k);

             qtflx_up(i,j,nk) = sd%qtflx_up(k)
             qtflx_dn(i,j,nk) = sd%qtflx_dn(k)
             omega_up(i,j,nk) = sd%omega_up(k)
             omega_dn(i,j,nk) = sd%omega_dn(k)
          enddo

          hdt_adv_int (i,j,:)=tdt_dyn_int(i,j,:)+qvdt_dyn_int(i,j,:)-qidt_dyn_int(i,j,:)+dgz_dyn_int(i,j,:);
          hdt_hadv_int(i,j,:)=hdt_adv_int(i,j,:)-hdt_vadv_int(i,j,:)
          hdt_dyn_int (i,j,:)=hdt_adv_int(i,j,:)+hdp_dyn_int (i,j,:)

          hdt_sum_int(i,j,:)=hdt_dyn_int(i,j,:)+tdt_rad_int(i,j,:)+tdt_dif_int(i,j,:)+qvdt_dif_int(i,j,:)
!note:    qvdt_dif_int = lhflx; tdt_dif_int = shflx
!         hdt_sum_int(i,j)=hdt_dyn_int(i,j,kmax)+tdt_rad_int(i,j,kmax)+shflx(i,j)+lhflx(i,j)

!for verification purpose, the current value of hdt_ver_int should equal hdt_sum_int at previous time-step
          hdt_tot_int(i,j) = (hfint0(i,j)-hfint0_old(i,j))/delt
          hdt_ver_int(i,j) = hdt_tot_int(i,j) - hdt_dgz_adj(i,j)
        endif !do_mse_budget

!========Find source air, and do adiabatic cloud lifting======================
          ksrc  =sd%ksrc
          zsrc  =sd%zsrc
          psrc  =sd%psrc
          thcsrc=sd%thcsrc
          qctsrc=sd%qctsrc
          hlsrc =sd%hlsrc

          rkm_shallow=rkm_sh

          if (do_stochastic_rkm) then
             rkm_shallow = rkm_shallow * frac_rkm_pert * rand_s(i,j)
          endif

          if (do_lands) then
             cpn % auto_th0 = auto_th0 * (1. + landfact_m * sd%land)
             call qt_parcel_k (sd%qs(1), qstar(i,j), pblht(i,j), sd%tke, sd%land, 0.0, &
                  pblht0, 1.0, lofactor0, lochoice, qctsrc, lofactor)
             rkm_shallow = rkm_shallow  * lofactor
          end if
          if (do_peff_land) then
             lofactor= 1.- sd%land*(1.- lofactor0)
             cpn % peff_l = peff_l  * lofactor
             cpn % peff_i = peff_i  * lofactor
          end if

          call adi_cloud_k(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, Uw_p, do_fast, do_ice, ac)
          !ac % usrc = sd%u(sd%ktoppbl)
          !ac % vsrc = sd%v(sd%ktoppbl)
          !if (ac%plfc.eq.0) ac%plfc=psrc
          !if (ac%plnb.eq.0) ac%plnb=psrc
          cino (i,j) = ac%cin
          capeo(i,j) = ac%cape
          plcl (i,j) = ac%plcl
          zlcl (i,j) = ac%zlcl-sd%zs(0)
          plfc (i,j) = ac%plfc
          plnb (i,j) = ac%plnb

          do k = 1,kmax
             nk = kmax+1-k
             hmo  (i,j,nk) = sd%hm(k);
             hms  (i,j,nk) = sd%hms(k);
             abu  (i,j,nk) = ac%buo(k);
             dthvdp(i,j,nk)= sd%dthvdp(k)
          end do

          rkm_s(i,j) = rkm_shallow

          if (do_fast) then
             if (ac%klcl.eq.0 .or. ac%plcl.eq.sd%ps(1) .or. ac%plcl.lt.20000.) then
                ocode(i,j)=1; cbmf_shallow=0.; goto 100
             end if
             if (ac%plfc.lt.50000.) then
                ocode(i,j)=2; cbmf_shallow=0.; goto 100
             end if
          end if

!========Cumulus closure to determine cloud base mass flux===================

          cbmf_old=cbmfo(i,j); cc%cbmf=cbmf_old;

          if (iclosure.eq.0) then
             call cclosure_bretherton(sd%tke, cpn, sd, Uw_p, ac, cc)
          else if (iclosure.eq.1) then
             call cclosure_implicit(sd%tke, cpn, sd, Uw_p, ac, cc, delt, rkm_shallow, &
                  do_coldT, sd1, ac1, cc1, cp, ct, ier, ermesg)
             if (ier /= 0) then
               call error_mesg ('subroutine uw_conv iclosure=1 ', ermesg, FATAL)
             endif
          else if (iclosure.eq.2) then
             call cclosure_relaxwfn(sd%tke, cpn, sd, Uw_p, ac, cc, cp, ct, delt,  &
                  rkm_shallow, do_coldT, sd1, ac1, cc1, cp1, ct1, ier, ermesg)
             if (ier /= 0) then
               call error_mesg ('subroutine uw_conv iclosure=2 ', ermesg, FATAL)
             endif
          end if

          cbmfo(i,j) = cc%cbmf
          wrelo(i,j) = cc%wrel
          ufrco(i,j) = cc%ufrc

          if (ac%klcl.eq.0 .or. ac%plcl.eq.sd%ps(1) .or. ac%plcl.lt.20000.) then
             ocode(i,j)=1; cbmf_shallow=0.; goto 100
          end if
          if (ac%plfc.lt.50000.) then
             ocode(i,j)=2; cbmf_shallow=0.; goto 100
          end if
          if(cc%cbmf.lt.1.e-6 .or. cc%wrel.eq.0.) then
             ocode(i,j)=3; cbmf_shallow=0.; goto 100
          end if
          if (do_eis_limit) then
             if (sd%eis .gt. eis_max) then
               ocode(i,j)=10; cbmf_shallow=0.; goto 100
             end if
          end if
          if (do_eis_limitn) then
             if (sd%eis .gt. eis_max) then
               ocode(i,j)=10; cbmf_shallow=0.; goto 200
             end if
          end if
          if (do_lts_limit) then
             if (sd%lts .gt. eis_max) then
               ocode(i,j)=10; cbmf_shallow=0.; goto 100
             end if
          end if
          if (do_lts_limitn) then
             if (sd%lts .gt. eis_max) then
               ocode(i,j)=10; cbmf_shallow=0.; goto 200
             end if
          end if

!========Do shallow cumulus plume calculation================================

          cbmf_shallow = cc%cbmf
          cpn%do_ppen=do_ppen
          cpn%rpen   =rpen
          call cumulus_plume_k(cpn, sd, ac, cp, rkm_shallow, cbmf_shallow, cc%wrel, cc%scaleh, Uw_p, ier, ermesg)
          if (ier /= 0) then
            call error_mesg ('subroutine uw_conv', ermesg, FATAL)
          endif
          if(cp%ltop.lt.cp%krel+2 .or. cp%let.le.cp%krel+1) then
             ocode(i,j)=4; cbmf_shallow=0.; goto 100
          end if
          if(cp%cldhgt.ge.cldhgt_max) then
             ocode(i,j)=5; cbmf_shallow=0.; goto 100
          end if

!========Calculate cumulus produced tendencies===============================

          call cumulus_tend_k(cpn, sd, Uw_p, cp, ct, do_coldT)

!========Unpack convective tendencies========================================
          do k = 1,cp%ltop
             nk = kmax+1-k
             uten  (i,j,nk) = ct%uten (k)
             vten  (i,j,nk) = ct%vten (k)
             qlten (i,j,nk) = ct%qlten(k)
             qiten (i,j,nk) = ct%qiten(k)
             qaten (i,j,nk) = ct%qaten(k)
             qnten (i,j,nk) = ct%qnten(k)
             qvten (i,j,nk) = ct%qvten(k)
             ice_pflx(i,j,nk) = cp%ppti(k)
             liq_pflx(i,j,nk) = cp%pptr(k)
             tten  (i,j,nk) = ct%tten (k)
             rhos0j = sd%ps(k)/(rdgas*0.5*(cp%thvbot(k+1)+cp%thvtop(k))*sd%exners(k))
             pflx  (i,j,nk) = ct%pflx (k)
             hlflx (i,j,nk) = ct%hlflx (k)
             qtflx (i,j,nk) = ct%qctflx(k)
             nqtflx(i,j,nk) = ct%nqtflx(k)
             tten_pevap (i,j,nk) = ct%tevap (k)
             qvten_pevap(i,j,nk) = ct%qevap (k)

             qldet_s(i,j,nk)= ct%qldet(k)
             qidet_s(i,j,nk)= ct%qidet(k)
             qadet_s(i,j,nk)= ct%qadet(k)
             qndet_s(i,j,nk)= ct%qndet(k)

             cqa_s (i,j,nk) = cp%ufrc(k)
             cql_s (i,j,nk) = cp%qlu(k)
             cqi_s (i,j,nk) = cp%qiu(k)
             cqn_s (i,j,nk) = cp%qnu(k)
             cldqa (i,j,nk) = cqa_s (i,j,nk)
             cldql (i,j,nk) = cql_s (i,j,nk)
             cldqi (i,j,nk) = cqi_s (i,j,nk)
             cldqn (i,j,nk) = cqn_s (i,j,nk)

             if (include_emf_s) then
                cmf_s (i,j,nk) = cp%umf(k) + cp%emf(k)
                cmf   (i,j,nk) = cp%umf(k) + cp%emf(k)
             else
                cmf_s (i,j,nk) = cp%umf(k)
                cmf   (i,j,nk) = cp%umf(k)
             end if
             dbuodp_s(i,j,nk)=cp%dbuodp(k);
             buo_s (i,j,nk) = cp%buo(k)
             wuo_s (i,j,nk) = cp%wu (k)
             peo   (i,j,nk) = cp%peff(k)
             fero_s(i,j,nk) = cp%fer(k)
             fdro_s(i,j,nk) = cp%fdr(k)
             fdrso_s(i,j,nk)= cp%fdrsat(k)*cp%fdr(k)!*cp%umf(k)

             do n = 1, size(trtend,4)
              trevp_s(i,j,nk,n) = ct%trevp(k,n)
             enddo
          enddo
          cush_s(i,j)  = cp%cush
          snow  (i,j)  = ct%snow
          rain  (i,j)  = ct%rain
          denth (i,j)  = ct%denth
          dqtmp (i,j)  = ct%dqtmp
          dting (i,j)  = ct%dting
          cpool (i,j)  = ct%cpool
          nbuo_s(i,j)  = cp%nbuo
          pcb_s (i,j)  = cp%prel
          pct_s (i,j)  = cp%ptop

! make sure the predicted tracer tendencies do not produce negative
! tracers due to convective tendencies. if necessary, adjust the
! tendencies.
          if (do_deep) then
            trtend_t = ct%trten
            trwet_t  = ct%trwet
          else
!f1p
!            call check_tracer_realizability (kmax, size(trtend,4), delt, &
!                                           cp%tr, ct%trten, ct%trwet)
             call check_tracer_realizability (kmax, size(trtend,4), delt, &
                                           cp%tr, ct%trten, ct%trwet, pmass(i,j,:), tracer_check_type, rn = rn )
            do k = 1,cp%ltop
              nk = kmax+1-k
              do n = 1, size(trtend,4)
                trtend(i,j,nk,n) = ct%trten(k,n) + ct%trwet(k,n)
                trwet(i,j,nk,n)  = ct%trwet(k,n)
                rn_diag(i,j,nk,n) = rn(k,n)
              enddo
            enddo
          end if

!========Option for deep convection=======================================
100       if (do_deep) then
             cbmf_deep = 0.
             rkm_dp = dpc%rkm_dp1
             crh_th = sd%land*dpc%crh_th_land+(1.-sd%land)*dpc%crh_th_ocean
             tmp = max(min (sd%crh, 1.0), 0.0)
             del_crh = tmp - crh_th
             dcrh0  = crh_max-crh_th
             if (del_crh .gt. 0.) then
                cbmf_deep = 0.0001 !first assuming existence of deep convective cloud base mass flux
                dcrh = del_crh/dcrh0; dcrh = dcrh**(norder) !dcrh = dcrh**(1./norder)
                if (dcrh.gt.1) then
                   rkm_dp = dpc%rkm_dp2
                else
                   rkm_dp  = dpc%rkm_dp1 + dcrh * (dpc%rkm_dp2 - dpc%rkm_dp1)
                end if

                lofactor= 1.- sd%land*(1.- dpc%lofactor_d) !option for introducing land difference
                if (do_lod_rkm) then
                   rkm_dp       = rkm_dp  * lofactor
                elseif (do_lod_cfrac) then
                   dpc % cfrac_d= cfrac_d * lofactor
                elseif (do_lod_tcrit) then
                   dpc % tcrit_d= tcrit_d * lofactor
                elseif (do_lod_cape) then
                   dpc % cape_th= cape_th * lofactor
                elseif (do_lod_tau) then
                   dpc % tau_dp = tau_dp  * lofactor
                end if

                if (do_stochastic_rkm) then
                 rkm_dp = rkm_dp * frac_rkm_pert * rand_d(i,j)
                end if

             end if

             if (do_plev_umf) then
                if (cp%umf_plev .lt. shallow_umf_thresh) then
                   cbmf_deep = 0
                endif
             endif

             dpn % do_ppen  = dpc % do_ppen_d
             dpn % rpen     = dpc % rpen_d
             dpn % do_pevap = dpc % do_pevap_d
             dpn % cfrac    = dpc % cfrac_d
             dpn % hcevap   = dpc % hcevap_d
             dpn % hcevappbl= dpc % hcevappbl_d
             dpn % pblfac   = dpc % pblfac_d
             dpn % tcrit    = dpc % tcrit_d
             dpn % auto_th0 = dpc % auto_th0_d
             dpn % peff_l   = dpc % peff_l_d
             dpn % peff_i   = dpc % peff_i_d
             dpn % mixing_assumption = dpc % mixing_assumption_d
             dpn % do_forcedlifting  = dpc % do_forcedlifting_d

             if (idpchoice.eq.0) then
                call  dpconv0(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                      rkm_dp, cbmf_deep, cp1, ct1, ocode(i,j), ier, ermesg)
             else if (idpchoice.eq.1) then
                call  dpconv1(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                      rkm_dp, cbmf_deep, sd1, ac1, cp1, ct1, ocode(i,j), dcapedm_d(i,j), ier, ermesg)
             else if (idpchoice.eq.2) then
                call  dpconv2(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                      rkm_dp, cbmf_deep, sd1, ac1, cp1, ct1, ocode(i,j),            &
                      dcwfndm_d(i,j), dcapedm_d(i,j), cwfn_d(i,j), lat(i,j), lon(i,j), ier, ermesg)
             else if (idpchoice.eq.3) then
                call  dpconv3(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                      rkm_dp, cbmf_deep, sd1, ac1, cp1, ct1, ocode(i,j),            &
                      dcwfndm_d(i,j), dcapedm_d(i,j), cwfn_d(i,j), lat(i,j), lon(i,j), ier, ermesg)
             end if
             if (ier /= 0) then
                call error_mesg ('uw_conv calling dpconv', ermesg, FATAL)
             endif

             if (cbmf_deep.eq.0) then
                call cp_clear_k(cp1); cp1%maxcldfrac=1.; cp1%cush=-1;
                call ct_clear_k(ct1);
             end if

             do k = 1,kmax !cp1%ltop
                nk = kmax+1-k
                uten_d  (i,j,nk) = ct1%uten (k)
                vten_d  (i,j,nk) = ct1%vten (k)
                qlten_d (i,j,nk) = ct1%qlten(k)
                qiten_d (i,j,nk) = ct1%qiten(k)
                qaten_d (i,j,nk) = ct1%qaten(k)
                qnten_d (i,j,nk) = ct1%qnten(k)
                qvten_d (i,j,nk) = ct1%qvten(k)
                qtten_d (i,j,nk) = ct1%qctten(k)
                tten_d  (i,j,nk) = ct1%tten (k)
                qldet_d (i,j,nk) = ct1%qldet(k)
                qidet_d (i,j,nk) = ct1%qidet(k)
                qadet_d (i,j,nk) = ct1%qadet(k)
                qndet_d (i,j,nk) = ct1%qndet(k)
                pflx_d  (i,j,nk) = ct1%pflx (k)
                hlflx_d (i,j,nk) = ct1%hlflx (k)
                qtflx_d (i,j,nk) = ct1%qctflx(k)
                nqtflx_d(i,j,nk) = ct1%nqtflx(k)
                tten_pevap_d (i,j,nk) = ct1%tevap (k)
                qvten_pevap_d(i,j,nk) = ct1%qevap (k)

                cqa_d   (i,j,nk) = cp1%ufrc(k)
                cql_d   (i,j,nk) = cp1%qlu(k)
                cqi_d   (i,j,nk) = cp1%qiu(k)
                cqn_d   (i,j,nk) = cp1%qnu(k)

                if (include_emf_d) then
                   cmf_d (i,j,nk) = cp1%umf(k) + cp1%emf(k)
                else
                   cmf_d (i,j,nk) = cp1%umf(k)
                end if
                dbuodp_d(i,j,nk) = cp1%dbuodp(k);
                buo_d   (i,j,nk) = cp1%buo(k)
                wuo_d   (i,j,nk) = cp1%wu (k)
                fero_d  (i,j,nk) = cp1%fer(k)
                fdro_d  (i,j,nk) = cp1%fdr(k)
                fdrso_d (i,j,nk) = cp1%fdrsat(k)*cp1%fdr(k)!*cp1%umf(k)
                do n = 1, size(trtend,4)
                   trevp_d(i,j,nk,n) = ct1%trevp(k,n)
                enddo
             enddo
             snow_d  (i,j)  = ct1%snow
             rain_d  (i,j)  = ct1%rain
             cbmf_d  (i,j)  = cbmf_deep
             denth_d (i,j)  = ct1%denth
             dting_d (i,j)  = ct1%dting
             dqtmp_d (i,j)  = ct1%dqtmp
             cush_d  (i,j)  = cp1%cush
             rkm_d   (i,j)  = rkm_dp;
             nbuo_d  (i,j)  = cp1%nbuo
             pcb_d   (i,j)  = cp1%prel
             pct_d   (i,j)  = cp1%ptop

             trtend_t = trtend_t+ct1%trten
             trwet_t  = trwet_t +ct1%trwet
!<f1p
             trtend_t_nc = trtend_t
             trwet_t_nc  = trwet_t
!>
!
!             call check_tracer_realizability (kmax, size(trtend,4), delt, &
!                                              cp1%tr, trtend_t, trwet_t)
             call check_tracer_realizability (kmax, size(trtend,4), delt, &
                                              cp1%tr, trtend_t, trwet_t, pmass(i,j,:), tracer_check_type, rn = rn    )
             do k = 1,kmax!cp1%ltop
               nk = kmax+1-k
               do n = 1, size(trtend,4)
                 trtend(i,j,nk,n) = trtend_t(k,n) + trwet_t(k,n)
                 trwet(i,j,nk,n)  = trwet_t(k,n)
!f1p
                 trtend_nc(i,j,nk,n) = trtend_t_nc(k,n) + trwet_t_nc(k,n)
                 rn_diag(i,j,nk,n) = rn(k,n)
!>
               enddo
             enddo

             uten  (i,j,:) = uten  (i,j,:) + uten_d  (i,j,:)
             vten  (i,j,:) = vten  (i,j,:) + vten_d  (i,j,:)
             qlten (i,j,:) = qlten (i,j,:) + qlten_d (i,j,:)
             qiten (i,j,:) = qiten (i,j,:) + qiten_d (i,j,:)
             qaten (i,j,:) = qaten (i,j,:) + qaten_d (i,j,:)
             qnten (i,j,:) = qnten (i,j,:) + qnten_d (i,j,:)
             qvten (i,j,:) = qvten (i,j,:) + qvten_d (i,j,:)
             tten  (i,j,:) = tten  (i,j,:) + tten_d  (i,j,:)
             pflx  (i,j,:) = pflx  (i,j,:) + pflx_d  (i,j,:)
             hlflx (i,j,:) = hlflx (i,j,:) + hlflx_d (i,j,:)
             qtflx (i,j,:) = qtflx (i,j,:) + qtflx_d (i,j,:)
             nqtflx(i,j,:) = nqtflx(i,j,:) + nqtflx_d(i,j,:)
             cmf   (i,j,:) = cmf_s (i,j,:) + cmf_d   (i,j,:)
             tten_pevap (i,j,:)=tten_pevap (i,j,:) + tten_pevap_d (i,j,:)
             qvten_pevap(i,j,:)=qvten_pevap(i,j,:) + qvten_pevap_d(i,j,:)

             if (do_new_convcld) then
                cqa(kmax+1)=0.; cql(kmax+1)=0.; cqi(kmax+1)=0.; cqn(kmax+1)=0.;
                do k = 1,kmax
                   cqa(k) =cqa_s(i,j,k)+cqa_d(i,j,k)
                   if (cqa(k).ne.0.) then
                      cql(k)=(cql_s(i,j,k)*cqa_s(i,j,k)+cql_d(i,j,k)*cqa_d(i,j,k))/cqa(k)
                      cqi(k)=(cqi_s(i,j,k)*cqa_s(i,j,k)+cqi_d(i,j,k)*cqa_d(i,j,k))/cqa(k)
                      cqn(k)=(cqn_s(i,j,k)*cqa_s(i,j,k)+cqn_d(i,j,k)*cqa_d(i,j,k))/cqa(k)
                   else
                      cql(k)=0.
                      cqi(k)=0.
                      cqn(k)=0.
                   end if
                   cqa(k) = min(cqa(k),1.0)
                end do
                do k = 1,kmax
                   cldqa(i,j,k)=(cqa(k)+cqa(k+1))*0.5
                   cldql(i,j,k)=(cql(k)+cql(k+1))*0.5
                   cldqi(i,j,k)=(cqi(k)+cqi(k+1))*0.5
                   cldqn(i,j,k)=(cqn(k)+cqn(k+1))*0.5
                end do

                do k = 1,kmax
                   cltc (i,j) = max(cltc(i,j),cldqa(i,j,k)) !assuming maximum overlap
                end do
             end if

             snow  (i,j)  = snow  (i,j) + snow_d  (i,j)
             rain  (i,j)  = rain  (i,j) + rain_d  (i,j)
             denth (i,j)  = denth (i,j) + denth_d (i,j)
             dting (i,j)  = dting (i,j) + dting_d (i,j)
             dqtmp (i,j)  = dqtmp (i,j) + dqtmp_d (i,j)
             cpool (i,j)  = cpool (i,j) + ct1%cpool

             feq_c (i,j)  = max(feq_s(i,j), feq_d(i,j))
             pcb_c (i,j)  = max(pcb_s(i,j), pcb_d(i,j))
             pct_c (i,j)  = min(max(pct_s(i,j),0.), max(pct_d(i,j),0.))
             !cbmfo (i,j)  = cc%cbmf
             !cwfno (i,j)  = cc%cwfn
          end if
!========End of do_deep, Option for deep convection=======================================

200       if (do_prog_gust) then
             gusto(i,j)=(gusto(i,j)+geff*cpool(i,j)*delt)/(1+delt/tau_gust)
          endif

!subtract parameterized convective mass flux
          do k = 1,kmax
             tmp=cmf(i,j,k)*Uw_p%grav;
             if ((-omega_up(i,j,k).gt.tmp) .and. (tmp.gt.0)) then
                omgmc_up(i,j,k) = omega_up(i,j,k)+tmp;
             else
                omgmc_up(i,j,k) = omega_up(i,j,k);
             endif
          enddo

       enddo
    enddo

    call sd_end_k(sd)
    call sd_end_k(sd1)
    call ac_end_k(ac)
    call ac_end_k(ac1)
    call cp_end_k(cp)
    call cp_end_k(cp1)
    call ct_end_k(ct)
    call ct_end_k(ct1)
    if (_ALLOCATED ( cpn%tracername    ))  deallocate ( cpn%tracername    )
    if (_ALLOCATED ( cpn%tracer_units  ))  deallocate ( cpn%tracer_units  )
    if (_ALLOCATED ( cpn%wetdep        ))  deallocate ( cpn%wetdep        )
    if (_ALLOCATED ( dpn%tracername    ))  deallocate ( dpn%tracername    )
    if (_ALLOCATED ( dpn%tracer_units  ))  deallocate ( dpn%tracer_units  )
    if (_ALLOCATED ( dpn%wetdep        ))  deallocate ( dpn%wetdep        )

    if (do_uwcmt) then
        half_delt = delt*0.5
        cp_inv    = 1. / Uw_p%cp_air
        dissipative_heat(:,:,:) = -((ub(:,:,:) + half_delt*uten(:,:,:))*uten(:,:,:) + &
                                    (vb(:,:,:) + half_delt*vten(:,:,:))*vten(:,:,:))*cp_inv
        !tten(:,:,:) = tten(:,:,:) + dissipative_heat(:,:,:)
    else
        uten=0.;
        vten=0.;
    end if

    if ( prevent_unreasonable ) then
       temp = qtr(:,:,:,nql)/delt + qlten(:,:,:)
       where (temp(:,:,:) .lt. 0. .and. qlten(:,:,:) .ne. 0.)
         tten (:,:,:) = tten (:,:,:) - temp(:,:,:)*HLV/CP_AIR
         qvten(:,:,:) = qvten(:,:,:) + temp(:,:,:)
         qlten(:,:,:) = qlten(:,:,:) - temp(:,:,:)
       end where

       temp = qtr(:,:,:,nqi)/delt + qiten(:,:,:)
       where (temp(:,:,:) .lt. 0. .and. qiten(:,:,:) .ne. 0.)
         tten (:,:,:) = tten (:,:,:) - temp(:,:,:)*HLS/CP_AIR
         qvten(:,:,:) = qvten(:,:,:) + temp(:,:,:)
         qiten(:,:,:) = qiten(:,:,:) - temp(:,:,:)
       end where

!The following 3 lines are in do_limit_uw; commented out here for reproducing earlier am4
       where (abs(qlten(:,:,:)+qiten(:,:,:))*delt .lt. 1.e-10 )
         qaten(:,:,:) = 0.0
       end where

       temp = qtr(:,:,:,nqa)/delt + qaten(:,:,:)
       where (temp(:,:,:) .lt. 0. .and. qaten(:,:,:) .ne. 0.)
         qaten(:,:,:) = qaten(:,:,:) - temp(:,:,:)
       end where
       where (temp(:,:,:)*delt .gt. 1. .and. qaten(:,:,:) .ne. 0.)
         qaten(:,:,:) = (1. - qtr(:,:,:,nqa))/delt
       end where

       if (do_qn) then
          temp = qtr(:,:,:,nqn)/delt + qnten(:,:,:)
          where (temp(:,:,:) .lt. 0. .and. qnten(:,:,:) .ne. 0.)
            qnten(:,:,:) = qnten(:,:,:) - temp(:,:,:)
          end where
       end if

      !rescaling to prevent negative specific humidity for each grid point
      scale_uw=1.0
      if (do_rescale) then
         temp=1.0
         do k=1,kmax
            do j=1,jmax
               do i=1,imax
                  qtin =  qv(i,j,k)
                  dqt  =  qvten(i,j,k) * delt
                  if ( dqt.lt.0 .and. qtin+dqt.lt.1.e-10 ) then
                     temp(i,j,k) = max( 0.0, -(qtin-1.e-10)/dqt )
                  endif
               enddo
            enddo
         enddo
!scaling factor for each column is the minimum value within that column
         scale_uw = minval(temp,dim=3)
      endif
      !rescaling to prevent excessive temperature tendencies
      if (do_rescale_t) then
         !temp_min=300.; temp_max=200.
         do k=1,kmax
            do j=1,jmax
               do i=1,imax
                  !temp_max = max(temp_max,tb(i,j,k))
                  !temp_min = min(temp_min,tb(i,j,k))
                  tnew  =  tb(i,j,k) + tten(i,j,k) * delt
                  if ( tnew > tmax0 ) then
                     temp_1 = 0.0
                     !print *, 'WARNING: scale_uw = 0 to prevent large T tendencies in UW'
                     !print *, i,j,'lev=',k,'pressure=',pmid(i,j,k),'tb=',tb(i,j,k),'tten=',tten(i,j,k)*delt
                     !print *, 'lat=', sd%lat, 'lon=', sd%lon, 'land=',sd%land
                  else
                     temp_1 = 1.0
                  endif
                  !scaling factor for each column is the minimum value within that column
                  scale_uw(i,j) = min( temp_1, scale_uw(i,j))
               enddo
            enddo
         enddo
      endif

      if (do_rescale .or. do_rescale_t) then
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              uten (i,j,k)  = scale_uw(i,j) * uten (i,j,k)
              vten (i,j,k)  = scale_uw(i,j) * vten (i,j,k)
              tten (i,j,k)  = scale_uw(i,j) * tten (i,j,k)
              qvten(i,j,k)  = scale_uw(i,j) * qvten(i,j,k)
              qlten(i,j,k)  = scale_uw(i,j) * qlten(i,j,k)
              qiten(i,j,k)  = scale_uw(i,j) * qiten(i,j,k)
              qaten(i,j,k)  = scale_uw(i,j) * qaten(i,j,k)
              if (do_qn) qnten(i,j,k) = scale_uw(i,j) * qnten(i,j,k)
              if (k.eq.kmax) then
                rain(i,j) = scale_uw(i,j) * rain(i,j)
                snow(i,j) = scale_uw(i,j) * snow(i,j)
                rain_d(i,j)=scale_uw(i,j) * rain_d(i,j)
                snow_d(i,j)=scale_uw(i,j) * snow_d(i,j)
                cpool(i,j)= scale_uw(i,j) * cpool(i,j)
              endif
            end do
          end do
        end do
      end if
    endif  !end of prevent_unreasonable

    if (do_mse_budget) then
      tv=0.; tvn=0.; dz=0.; dz_n=0.; zint_n=0.; zmid_n=0.; dgz_conv_int=0.;
      tv  = tb(:,:,:) *(1 + 0.608*qv (:,:,:) - qtr (:,:,:,nql)-qtr (:,:,:,nqi))
      tvn = (tb(:,:,:)+tten(:,:,:)*delt) * (1 + 0.608*(qv(:,:,:)+qvten(:,:,:)*delt) &
            - (qtr(:,:,:,nql)+qlten(:,:,:)*delt)-(qtr(:,:,:,nqi)+qiten(:,:,:)*delt))
      zint_n(:,:,kmax+1)=zint(:,:,kmax+1);
      do k=kmax,1,-1
         dz     (:,:,k)=zint(:,:,k)-zint(:,:,k+1)
         dz_n   (:,:,k)=dz(:,:,k)*tvn(:,:,k)/tv(:,:,k)
         zint_n(:,:,k)=zint_n(:,:,k+1)+dz_n(:,:,k)
         zmid_n(:,:,k)=(zint_n(:,:,k)+zint_n(:,:,k+1))*0.5
      enddo
      dgz_conv_int(:,:,:)= Uw_p%grav*(zmid_n(:,:,:)-zmid(:,:,:))/delt * pmass(:,:,:)
      tdt_conv_int(:,:,:)= Uw_p%cp_air * tten(:,:,:) * pmass(:,:,:)
      qdt_conv_int(:,:,:)=(Uw_p%HLv*qvten(:,:,:) - Uw_p%HLf*qiten(:,:,:)) * pmass(:,:,:)
      do k=2, kmax
        dgz_conv_int(:,:,k)= dgz_conv_int(:,:,k)+dgz_conv_int(:,:,k-1)
        tdt_conv_int(:,:,k)= tdt_conv_int(:,:,k)+tdt_conv_int(:,:,k-1)
        qdt_conv_int(:,:,k)= qdt_conv_int(:,:,k)+qdt_conv_int(:,:,k-1)
      end do
      hdt_conv_int(:,:,:) = tdt_conv_int(:,:,:)+qdt_conv_int(:,:,:)+dgz_conv_int(:,:,:)
    endif !do_mse_budget

    do i=1,imax
       do j=1,jmax
          cush(i,j) = cush_s(i,j)
          if (cush_s(i,j) .gt. 0) feq_s(i,j)=1.;
          if (cush_d(i,j) .gt. 0) feq_d(i,j)=1.;
          do k=1,kmax
             cfq_s(i,j,k) = 0.
             cfq_d(i,j,k) = 0.
             if (cmf_s(i,j,k) .gt. 0.) cfq_s(i,j,k) = 1.
             if (cmf_d(i,j,k) .gt. 0.) cfq_d(i,j,k) = 1.
          enddo
       enddo
    enddo

    if (do_imposing_rad_cooling) then
      tten_rad(:,:,:) = 0.0
      do j = 1,jmax
         do i=1,imax
           do k = 1,kmax
              if (tb(i,j,k) .gt. t_thresh) then
                tten_rad (i,j,k) = cooling_rate/86400.
              else
                tten_rad (i,j,k) = (t_strato-tb(i,j,k))/(tau_rad*86400.)
              end if
           enddo
         enddo
      enddo
!      used = send_data( id_tten_rad_uwc,tten_rad*aday,Time, is, js, 1)
    end if

    if (do_imposing_forcing) then
      tten_forc(:,:,:) = 0.0
      qten_forc(:,:,:) = 0.0
      do j = 1,jmax
        do i=1,imax
           tten_forc(i,j,:)=0;
           qten_forc(i,j,:)=0;
           if (use_klevel) then
             k = klevel
             tten_forc(i,j,k)=tdt_rate/86400.
             qten_forc(i,j,k)=qdt_rate/86400.
           else
             kbot_tmp=1;
             ktop_tmp=kmax;
             do k=1,kmax
               if (pmid(i,j,k)>=7500) then
                 ktop_tmp=k
               end if
               if (pmid(i,j,k)>=85000) then
                 kbot_tmp=k
               end if
             enddo
             do k = kbot_tmp,ktop_tmp
                if (pmid(i,j,k)>pres_min .and. pmid(i,j,k)<=pres_max) then
                   tten_forc(i,j,k)=tdt_rate/86400.
                   qten_forc(i,j,k)=qdt_rate/86400.
                end if
             enddo
           end if
        enddo
      enddo
      used = send_data( id_tdt_forc_uwc,tten_forc*aday,Time, is, js, 1)
      used = send_data( id_qdt_forc_uwc,qten_forc*aday,Time, is, js, 1)
    end if

    !diagnostic output
    used = send_data( id_xpsrc_uwc,        xpsrc,              Time, is, js)
    used = send_data( id_xhlsrc_uwc,       xhlsrc,             Time, is, js)
    used = send_data( id_xqtsrc_uwc,       xqtsrc,             Time, is, js)
    used = send_data( id_hlsrc_avg,        hlsrc_avg,          Time, is, js)
    used = send_data( id_qtsrc_avg,        qtsrc_avg,          Time, is, js)
    used = send_data( id_cape_avg,         cape_avg,           Time, is, js)
    used = send_data( id_cin_avg,          cin_avg,            Time, is, js)
    used = send_data( id_omg_avg,          omg_avg,            Time, is, js)
    used = send_data( id_pblht_avg,        pblht_avg,          Time, is, js)

    used = send_data( id_tdt_pevap_uwc,    tten_pevap,         Time, is, js, 1)
    used = send_data( id_qdt_pevap_uwc,    qvten_pevap,        Time, is, js, 1)

    used = send_data( id_tdt_uwc,    tten ,        Time, is, js, 1)
    used = send_data( id_qdt_uwc,    qvten,        Time, is, js, 1)
    used = send_data( id_udt_uwc,    uten,         Time, is, js, 1)
    used = send_data( id_vdt_uwc,    vten,         Time, is, js, 1)
    used = send_data( id_cmf_uwc,    cmf,          Time, is, js, 1)
    used = send_data( id_cqa_uwc,    cldqa,        Time, is, js, 1)
    used = send_data( id_cql_uwc,    cldql,        Time, is, js, 1)
    used = send_data( id_cqi_uwc,    cldqi,        Time, is, js, 1)
    used = send_data( id_cqn_uwc,    cldqn,        Time, is, js, 1)
    used = send_data( id_cltc_uwc,   cltc*100.,    Time, is, js)
    used = send_data( id_pcb_uwc,    pcb_c*0.01,   Time, is, js)
    used = send_data( id_pct_uwc,    pct_c*0.01,   Time, is, js)

    used = send_data( id_cqa_uws,    cqa_s,        Time, is, js, 1)
    used = send_data( id_cql_uws,    cql_s,        Time, is, js, 1)
    used = send_data( id_cqi_uws,    cqi_s,        Time, is, js, 1)
    used = send_data( id_cqn_uws,    cqn_s,        Time, is, js, 1)
    used = send_data( id_pcb_uws,    pcb_s*0.01,   Time, is, js)
    used = send_data( id_pct_uws,    pct_s*0.01,   Time, is, js)

    used = send_data( id_tdt_uws, (tten-tten_d),   Time, is, js, 1)
    used = send_data( id_qdt_uws, (qvten-qvten_d), Time, is, js, 1)
    used = send_data( id_cmf_uws,    cmf_s,        Time, is, js, 1)
    used = send_data( id_wuo_uws,    wuo_s,        Time, is, js, 1)

    used = send_data( id_cfq_uws,    cfq_s,        Time, is, js, 1)
    used = send_data( id_fwu_uws,    wuo_s*cfq_s,  Time, is, js, 1)
    used = send_data( id_fqa_uws,    cqa_s*cfq_s,  Time, is, js, 1)
    used = send_data( id_fql_uws,    cql_s*cfq_s,  Time, is, js, 1)
    used = send_data( id_fqi_uws,    cqi_s*cfq_s,  Time, is, js, 1)
    used = send_data( id_fqn_uws,    cqn_s*cfq_s,  Time, is, js, 1)

    used = send_data( id_peo_uwc,    peo,          Time, is, js, 1)
    used = send_data( id_fer_uws,    fero_s,       Time, is, js, 1)
    used = send_data( id_fdr_uws,    fdro_s,       Time, is, js, 1)
    used = send_data( id_fdrs_uws,   fdrso_s,      Time, is, js, 1)

    used = send_data( id_hlflx_uwc,  hlflx,        Time, is, js, 1)
    used = send_data( id_qtflx_uwc,  qtflx,        Time, is, js, 1)

    used = send_data( id_hmo_uwc,    hmo,          Time, is, js, 1)
    used = send_data( id_hms_uwc,    hms,          Time, is, js, 1)
    used = send_data( id_abu_uwc,    abu,          Time, is, js, 1)
    used = send_data( id_buo_uws,    buo_s,        Time, is, js, 1)
    used = send_data( id_dbuodp_uws, dbuodp_s,     Time, is, js, 1)
    used = send_data( id_dbuodp_uwd, dbuodp_d,     Time, is, js, 1)
    used = send_data( id_dthvdp_uwc, dthvdp,       Time, is, js, 1)

    used = send_data( id_prec_uws, (rain+snow-rain_d-snow_d), Time, is, js )
    used = send_data( id_snow_uws, (snow-snow_d),      Time, is, js )
    used = send_data( id_prec_uwc, (rain+snow),        Time, is, js )
    used = send_data( id_snow_uwc, (snow),             Time, is, js )
    used = send_data( id_cin_uwc,  (cino),             Time, is, js )
    used = send_data( id_cape_uwc, (capeo),            Time, is, js )
    used = send_data( id_gust_uwc, (gusto),            Time, is, js )
    used = send_data( id_cpool_uwc,(cpool),            Time, is, js )
    used = send_data( id_crh_uwc,  (crho),             Time, is, js )
    used = send_data( id_fcrh_uws, (crho*feq_s),       Time, is, js )
    used = send_data( id_fcrh_uwd, (crho*feq_d),       Time, is, js )
    used = send_data( id_pblht_uwc,(pblht),            Time, is, js )
    used = send_data( id_tke_uwc,  (tkeo),             Time, is, js )
    used = send_data( id_tkep_uwc, (tkep),             Time, is, js )
    used = send_data( id_cbmf_uwc, (cbmfo),            Time, is, js )
    used = send_data( id_wrel_uwc, (wrelo),            Time, is, js )
    used = send_data( id_ufrc_uwc, (ufrco),            Time, is, js )
    used = send_data( id_plcl_uwc, (plcl*0.01),        Time, is, js )
    used = send_data( id_zlcl_uwc, (zlcl),             Time, is, js )
    used = send_data( id_plfc_uwc, (plfc*0.01),        Time, is, js )
    used = send_data( id_plnb_uwc, (plnb*0.01),        Time, is, js )
    used = send_data( id_zinv_uwc, (zinvo),            Time, is, js )
    used = send_data( id_cush_uws, (cush_s),           Time, is, js )
    used = send_data( id_dcin_uwc, (dcino),            Time, is, js )
    used = send_data( id_dcape_uwc,(dcapeo),           Time, is, js )
!    used = send_data( id_dwfn_uwc, (dwfno),            Time, is, js )
    used = send_data( id_enth_uwc, (denth),            Time, is, js )
    used = send_data( id_qtmp_uwc, (dqtmp),            Time, is, js )
    used = send_data( id_dting_uwc,(dting),            Time, is, js )
    used = send_data( id_ocode_uwc,(ocode),            Time, is, js )
    used = send_data( id_feq_uwc,  (feq_c),            Time, is, js )
    used = send_data( id_feq_uws,  (feq_s),            Time, is, js )
    used = send_data( id_rkm_uws,  (rkm_s),            Time, is, js )
    used = send_data( id_frkm_uws, (rkm_s*feq_s),      Time, is, js )
    used = send_data( id_scale_uwc,(scale_uw),         Time, is, js )
    used = send_data( id_scaletr_uwc,(scale_tr),       Time, is, js )

    used = send_data( id_lts_uwc,      lts,            Time, is, js )
    used = send_data( id_eis_uwc,      eis,            Time, is, js )

  if (do_mse_budget) then
    used = send_data( id_hfint0,        hfint0,                 Time, is, js)
    used = send_data( id_hfintn0,       hfintn      (:,:,kmax), Time, is, js)
    used = send_data( id_dpint0,        dpint       (:,:,kmax), Time, is, js)
    used = send_data( id_tdt_rad_int0,  tdt_rad_int (:,:,kmax), Time, is, js)
    used = send_data( id_tdt_dyn_int0,  tdt_dyn_int (:,:,kmax), Time, is, js)
    used = send_data( id_tdt_dif_int0,  tdt_dif_int (:,:,kmax), Time, is, js)
    used = send_data( id_qvdt_dyn_int0, qvdt_dyn_int(:,:,kmax), Time, is, js)
    used = send_data( id_qvdt_dif_int0, qvdt_dif_int(:,:,kmax), Time, is, js)
    used = send_data( id_dgz_dyn_int0,  dgz_dyn_int (:,:,kmax), Time, is, js)
    used = send_data( id_dgz_phy_int0,  dgz_phy_int (:,:,kmax), Time, is, js)
    used = send_data( id_hdp_dyn_int0,  hdp_dyn_int (:,:,kmax), Time, is, js)
    used = send_data( id_ddp_dyn_int0,  ddp_dyn_int (:,:,kmax), Time, is, js)
    used = send_data( id_hdt_vadv_int0, hdt_vadv_int(:,:,kmax), Time, is, js)
    used = send_data( id_hdt_hadv_int0, hdt_hadv_int(:,:,kmax), Time, is, js)
    used = send_data( id_hdt_sum_int0,  hdt_sum_int (:,:,kmax), Time, is, js)

    used = send_data( id_hdt_tot_int,  hdt_tot_int,  Time, is, js)
    used = send_data( id_hdt_ver_int,  hdt_ver_int,  Time, is, js)
    used = send_data( id_lhflx_uwc,    lhflx,        Time, is, js)
    used = send_data( id_shflx_uwc,    shflx,        Time, is, js)

    used = send_data( id_dpint,        dpint,        Time, is, js, 1)
    used = send_data( id_hfint,        hfint,        Time, is, js, 1)
    used = send_data( id_hfintn,       hfintn,       Time, is, js, 1)
    used = send_data( id_dgz_conv_int, dgz_conv_int, Time, is, js, 1)
    used = send_data( id_tdt_conv_int, tdt_conv_int, Time, is, js, 1)
    used = send_data( id_qdt_conv_int, qdt_conv_int, Time, is, js, 1)
    used = send_data( id_hdt_conv_int, hdt_conv_int, Time, is, js, 1)

    used = send_data( id_hdt_forc_int, hdt_forc_int, Time, is, js, 1)
    used = send_data( id_tdt_rad_int,  tdt_rad_int,  Time, is, js, 1)
    used = send_data( id_tdt_dyn_int,  tdt_dyn_int,  Time, is, js, 1)
    used = send_data( id_tdt_dif_int,  tdt_dif_int,  Time, is, js, 1)
    used = send_data( id_qvdt_dyn_int, qvdt_dyn_int, Time, is, js, 1)
    used = send_data( id_qvdt_dif_int, qvdt_dif_int, Time, is, js, 1)
    used = send_data( id_dgz_dyn_int,  dgz_dyn_int,  Time, is, js, 1)
    used = send_data( id_dgz_phy_int,  dgz_phy_int,  Time, is, js, 1)
    used = send_data( id_hdp_dyn_int,  hdp_dyn_int,  Time, is, js, 1)
    used = send_data( id_ddp_dyn_int,  ddp_dyn_int,  Time, is, js, 1)
    used = send_data( id_hdt_vadv_int, hdt_vadv_int, Time, is, js, 1)
    used = send_data( id_hdt_hadv_int, hdt_hadv_int, Time, is, js, 1)
    used = send_data( id_hdt_sum_int,  hdt_sum_int,  Time, is, js, 1)

    used = send_data( id_tdt_rad_uwc,  tdt_rad,      Time, is, js, 1)
    used = send_data( id_tdt_dyn_uwc,  tdt_dyn,      Time, is, js, 1)
    used = send_data( id_tdt_dif_uwc,  tdt_dif,      Time, is, js, 1)
    used = send_data( id_qvdt_dyn_uwc, qvdt_dyn,     Time, is, js, 1)
    used = send_data( id_qvdt_dif_uwc, qvdt_dif,     Time, is, js, 1)
    used = send_data( id_dgz_dyn_uwc,  dgz_dyn,      Time, is, js, 1)

    used = send_data( id_nqtflx_uwc,   nqtflx,       Time, is, js, 1)
    used = send_data( id_qtflx_up_uwc, qtflx_up,     Time, is, js, 1)
    used = send_data( id_qtflx_dn_uwc, qtflx_dn,     Time, is, js, 1)
    used = send_data( id_omgmc_up_uwc, omgmc_up,     Time, is, js, 1)
    used = send_data( id_omega_up_uwc, omega_up,     Time, is, js, 1)
    used = send_data( id_omega_dn_uwc, omega_dn,     Time, is, js, 1)
    used = send_data( id_pflx_uwc,     pflx,         Time, is, js, 1)

    used = send_data( id_nbuo_uws,     nbuo_s,       Time, is, js )
   endif !do_mse_budget

    if ( do_uwcmt ) then
      used = send_data( id_tdt_diss_uwc,  dissipative_heat, Time, is, js, 1)
    end if

    if ( do_strat ) then
       used = send_data( id_qldt_uwc,   qlten,  Time, is, js, 1)
       used = send_data( id_qidt_uwc,   qiten,  Time, is, js, 1)
       used = send_data( id_qadt_uwc,   qaten,  Time, is, js, 1)
       used = send_data( id_qndt_uwc,   qnten,  Time, is, js, 1)
       used = send_data( id_qtdt_uwc,  (qvten+qlten+qiten),Time, is, js, 1)

       used = send_data( id_qldt_uws,   qlten-qlten_d,  Time, is, js, 1)
       used = send_data( id_qidt_uws,   qiten-qiten_d,  Time, is, js, 1)
       used = send_data( id_qadt_uws,   qaten-qaten_d,  Time, is, js, 1)
       used = send_data( id_qndt_uws,   qnten-qnten_d,  Time, is, js, 1)
       used = send_data( id_qtdt_uws,  (qvten+qlten+qiten-qvten_d-qlten_d-qiten_d),Time, is, js, 1)

       used = send_data( id_qldet_uws,  qldet_s, Time, is, js, 1)
       used = send_data( id_qidet_uws,  qidet_s, Time, is, js, 1)
       used = send_data( id_qadet_uws,  qadet_s, Time, is, js, 1)
       used = send_data( id_qndet_uws,  qndet_s, Time, is, js, 1)
    end if
!f1p
    if ( allocated(id_rn) ) then
       do n = 1,size(id_rn)
          used = send_data( id_rn(n), rn_diag(:,:,:,n), Time, is, js, 1)
       end do
    end if

    if ( allocated(id_tracerdt_uwc) ) then
       do n = 1,size(id_tracerdt_uwc)
          used = send_data( id_tracerdt_uwc(n), trtend(:,:,:,n), Time, is, js, 1)
       end do
    end if
!f1p
    if ( allocated(id_tracerdt_uwc_nc) ) then
       do n = 1,size(id_tracerdt_uwc_nc)
          used = send_data( id_tracerdt_uwc_nc(n), trtend_nc(:,:,:,n), Time, is, js, 1)
       end do
    end if

    if ( allocated(id_tracerdt_uwc_col) ) then
       do n = 1,size(id_tracerdt_uwc_col)
          if ( id_tracerdt_uwc_col(n) > 0 ) then
            tempdiag = 0.
            do k = 1,kmax
               tempdiag(:,:) = tempdiag(:,:) + trtend(:,:,k,n) * pmass(:,:,k)
            end do
            used = send_data( id_tracerdt_uwc_col(n), tempdiag(:,:), Time, is, js)
          end if
       end do
    end if
!f1p
    if ( allocated(id_tracerdt_uwc_col_nc) ) then
       do n = 1,size(id_tracerdt_uwc_col_nc)
          if ( id_tracerdt_uwc_col_nc(n) > 0 ) then
            tempdiag = 0.
            do k = 1,kmax
               tempdiag(:,:) = tempdiag(:,:) + trtend_nc(:,:,k,n) * pmass(:,:,k)
            end do
            used = send_data( id_tracerdt_uwc_col_nc(n), tempdiag(:,:), Time, is, js)
          end if
       end do
    end if

    if ( allocated(id_tracerdtwet_uwc) ) then
       do n = 1,size(id_tracerdtwet_uwc)
          used = send_data( id_tracerdtwet_uwc(n), trwet(:,:,:,n), Time, is, js, 1)
       end do
    end if

!<<<fp
!this means that uw_wetdep = 0 if wet_uwc_col is not requested for tracer n

!     if ( allocated(id_tracerdtwet_uwc_col) ) then
!        uw_wetdep = 0.
!        do n = 1,size(id_tracerdtwet_uwc_col)
!           if ( id_tracerdtwet_uwc_col(n) > 0 ) then
!              tempdiag = 0.
!              do k = 1,kmax
!                tempdiag(:,:) = tempdiag(:,:) + trwet(:,:,k,n) * pmass(:,:,k)
!             end do
!             used = send_data( id_tracerdtwet_uwc_col(n), tempdiag(:,:), Time, is, js)
!             uw_wetdep(:,:,n) = tempdiag(:,:)
!           end if
!        end do
!     end if

!now calculated uw for all tracers

    uw_wetdep = 0.
    do n=1,ntracers
       tempdiag = 0.
       do k = 1,kmax
          tempdiag(:,:) = tempdiag(:,:) + trwet(:,:,k,n) * pmass(:,:,k)
       end do
       uw_wetdep(:,:,n) = tempdiag(:,:)
    end do

     if ( allocated(id_tracerdtwet_uwc_col) ) then
        do n = 1,size(id_tracerdtwet_uwc_col)
           if ( id_tracerdtwet_uwc_col(n) > 0 ) then
             used = send_data( id_tracerdtwet_uwc_col(n), uw_wetdep(:,:,n), Time, is, js)
           end if
        end do
     end if

!>>>

    if ( allocated(id_trevp_uwc) ) then
       do n = 1,size(id_trevp_uwc)
          used = send_data( id_trevp_uwc(n), trevp_s(:,:,:,n)+trevp_d(:,:,:,n), Time, is, js, 1)
       end do
    end if

!========Option for deep convection=======================================
    if (do_deep) then
       used=send_data( id_tdt_pevap_uwd, tten_pevap_d, Time, is, js, 1)
       used=send_data( id_qdt_pevap_uwd, qvten_pevap_d,Time, is, js, 1)
       used=send_data( id_tdt_uwd,   tten_d,         Time, is, js, 1)
       used=send_data( id_qdt_uwd,   qvten_d,        Time, is, js, 1)
       used=send_data( id_qtdt_uwd,  qtten_d,        Time, is, js, 1)
       used=send_data( id_cmf_uwd,   cmf_d,          Time, is, js, 1)
       used=send_data( id_buo_uwd,   buo_d,          Time, is, js, 1)
       used=send_data( id_wuo_uwd,   wuo_d,          Time, is, js, 1)

       used=send_data( id_cfq_uwd,   cfq_d,          Time, is, js, 1)
       used=send_data( id_fwu_uwd,   wuo_d*cfq_d,    Time, is, js, 1)
       used=send_data( id_fqa_uwd,   cqa_d*cfq_d,    Time, is, js, 1)
       used=send_data( id_fql_uwd,   cql_d*cfq_d,    Time, is, js, 1)
       used=send_data( id_fqi_uwd,   cqi_d*cfq_d,    Time, is, js, 1)
       used=send_data( id_fqn_uwd,   cqn_d*cfq_d,    Time, is, js, 1)

       used=send_data( id_fer_uwd,   fero_d,         Time, is, js, 1)
       used=send_data( id_fdr_uwd,   fdro_d,         Time, is, js, 1)
       used=send_data( id_fdrs_uwd,  fdrso_d,        Time, is, js, 1)
       used=send_data( id_cqa_uwd,   cqa_d,          Time, is, js, 1)
       used=send_data( id_cql_uwd,   cql_d,          Time, is, js, 1)
       used=send_data( id_cqi_uwd,   cqi_d,          Time, is, js, 1)
       used=send_data( id_cqn_uwd,   cqn_d,          Time, is, js, 1)
       used=send_data( id_feq_uwd,   feq_d,          Time, is, js )
       used=send_data( id_pcb_uwd,   pcb_d*0.01,     Time, is, js)
       used=send_data( id_pct_uwd,   pct_d*0.01,     Time, is, js)

       used=send_data( id_hlflx_uwd, hlflx_d,        Time, is, js, 1)
       used=send_data( id_qtflx_uwd, qtflx_d,        Time, is, js, 1)
       used=send_data( id_nqtflx_uwd,nqtflx_d,       Time, is, js, 1)
!       used=send_data( id_trtend_uwd, trtend,        Time, is, js, 1)
!       used=send_data( id_trwet_uwd,  trwet,         Time, is, js, 1)

       used=send_data( id_prec_uwd, (rain_d+snow_d),       Time, is, js )
       used=send_data( id_snow_uwd, (snow_d),              Time, is, js )
       used=send_data( id_cbmf_uwd, (cbmf_d),              Time, is, js )
       used=send_data( id_dcapedm_uwd,(dcapedm_d),         Time, is, js )
       used=send_data( id_dcwfndm_uwd,(dcwfndm_d),         Time, is, js )
       used=send_data( id_cwfn_uwd, (cwfn_d),              Time, is, js )
       used=send_data( id_cush_uwd, (cush_d),              Time, is, js )
       used=send_data( id_enth_uwd, (denth_d),             Time, is, js )
       used=send_data( id_rkm_uwd,  (rkm_d),               Time, is, js )
       used=send_data( id_frkm_uwd, (rkm_d*feq_d),         Time, is, js )
       used=send_data( id_rand_uwd, (rand_d),              Time, is, js )
       used=send_data( id_nbuo_uwd, nbuo_d,                Time, is, js )

       if ( do_strat ) then
          used=send_data( id_qldt_uwd,  qlten_d, Time, is, js, 1)
          used=send_data( id_qidt_uwd,  qiten_d, Time, is, js, 1)
          used=send_data( id_qadt_uwd,  qaten_d, Time, is, js, 1)
          used=send_data( id_qtdt_uwd,  (qvten_d+qlten_d+qiten_d),Time, is, js, 1)
          used=send_data( id_qldet_uwd, qldet_d, Time, is, js, 1)
          used=send_data( id_qidet_uwd, qidet_d, Time, is, js, 1)
          used=send_data( id_qadet_uwd, qadet_d, Time, is, js, 1)
          used=send_data( id_qndet_uwd, qndet_d, Time, is, js, 1)
       end if
       if ( allocated(id_trevp_uwd) ) then
         do n = 1,size(id_trevp_uwd)
           used = send_data( id_trevp_uwd(n), trevp_d(:,:,:,n), Time, is, js, 1)
         end do
       end if
    end if
!========Option for deep convection=======================================

    if (.not.apply_tendency) then
       uten=0.; vten=0.; tten=0.; qvten=0.; cmf=0.; rain=0.; snow=0.;
       qlten=0.; qiten=0.; qaten=0.; qnten=0.;
    end if

    if (zero_out_conv_area) then
       cldql=0.; cldqi=0.; cldqa=0.; cldqn=0.;
    end if

    if (do_imposing_rad_cooling) then
      do j = 1,jmax
         do i=1,imax
           tten(i,j,:) = tten (i,j,:) + tten_rad (i,j,:)
         enddo
       enddo
    end if

    if (do_imposing_forcing) then
       do j = 1,jmax
          do i=1,imax
             tten (i,j,:) = tten (i,j,:) + tten_forc(i,j,:)
             qvten(i,j,:) = qvten(i,j,:) + qten_forc(i,j,:)
          enddo
        enddo
    end if

    if (do_mse_budget) then
      do j = 1,jmax
        do i=1,imax
          hten(i,j)=0
          do k=1,kmax
             hten(i,j)=hten(i,j)+(cp_air*tten(i,j,k)+HLv*qvten(i,j,k)-HLv*qiten(i,j,k))*pmass(i,j,k)
          enddo
        enddo
      enddo
    endif

  END SUBROUTINE UW_CONV

!#####################################################################
!#####################################################################

  subroutine clearit(ac, cc, cp, ct, cp1, ct1)

    type(adicloud), intent(inout) :: ac
    type(cclosure), intent(inout) :: cc
    type(cplume),   intent(inout) :: cp,cp1
    type(ctend),    intent(inout) :: ct,ct1

    call ac_clear_k(ac);
    ac%klcl =0;  ac%klfc =0;  ac%klnb =0; ac%zlcl=0;

    cc%wrel=0.; cc%ufrc=0.; cc%scaleh=0.;

    call cp_clear_k(cp);  cp%maxcldfrac =1.;
    call ct_clear_k(ct);
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call ct_clear_k(ct1);

  end subroutine clearit


!#####################################################################

end MODULE UW_CONV_MOD
