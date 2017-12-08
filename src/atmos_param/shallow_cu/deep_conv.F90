MODULE DEEP_CONV_MOD

  use      fms_mod,         only : write_version_number
  use      Constants_Mod,   ONLY : tfreeze,HLv,HLf,HLs,CP_AIR,GRAV,Kappa,rdgas,rvgas
  use  conv_utilities_mod,  only : uw_params_init
  use  conv_utilities_k_mod,only : sd_init_k, sd_copy_k, sd_end_k,  &
                                   ac_init_k, ac_clear_k, ac_end_k, &
                                   pack_sd_k, adi_cloud_k, extend_sd_k,&
                                   exn_init_k, exn_end_k, findt_init_k,&
                                   findt_end_k, qt_parcel_deep_k, &
                                   adicloud, sounding, uw_params, findt_k, &
                                   exn_k, qsat_k, erfccc

  use  conv_plumes_k_mod,   only : cp_init_k, cp_end_k, cp_clear_k, &
                                   ct_init_k, ct_end_k, ct_clear_k, &
                                   cumulus_tend_k, cumulus_plume_k, &
                                   cplume, ctend, cpnlist

  use  conv_closures_mod,   only : cclosure_bretherton,   &
                                   cclosure_relaxcbmf, &
                                   cclosure_relaxwfn,  &
                                   cclosure_implicit, cclosure

!---------------------------------------------------------------------
  implicit none
  private
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

!-------  interfaces --------

  public  :: cpn_copy, dpconv0, dpconv1, dpconv2, dpconv3

  logical         :: module_is_initialized = .false.

  character(len=7) :: mod_name = 'deep_conv'

  public deepc
  type deepc
     real    :: cbmf0
     real    :: rkm_dp1
     real    :: rkm_dp2
     real    :: crh_th_land
     real    :: crh_th_ocean
     real    :: cape_th
     real    :: cin_th
     real    :: tau_dp
     real    :: rpen_d
     integer :: mixing_assumption_d
     logical :: do_ppen_d
     logical :: do_pevap_d
     real    :: cfrac_d
     real    :: hcevap_d
     real    :: hcevappbl_d
     real    :: pblfac_d
     real    :: dcapedm_th
     real    :: dcwfndm_th
     real    :: frac_limit_d
     real    :: lofactor_d
     real    :: tcrit_d
     real    :: auto_th0_d
     real    :: peff_l_d
     real    :: peff_i_d
     integer :: src_choice_d
     real    :: cwfn_th
     real    :: cin_fact
     real    :: wcrit_min_gust
     real    :: tau_dp_fact
     integer :: cgust_choice
     logical :: do_cgust_dp
     logical :: do_forcedlifting_d
  end type deepc

contains

!#####################################################################
!#####################################################################

  subroutine cpn_copy(cpn, dpn)
    type(cpnlist), intent(in)    :: cpn
    type(cpnlist), intent(inout) :: dpn

    dpn % do_umf_pbl         = cpn % do_umf_pbl
    dpn % do_qctflx_zero     = cpn % do_qctflx_zero
    dpn % do_hlflx_zero      = cpn % do_hlflx_zero
    dpn % do_new_qnact       = cpn % do_new_qnact
    dpn % do_2nd_act         = cpn % do_2nd_act
    dpn % do_downdraft       = cpn % do_downdraft
    dpn % do_conv_micro_N    = cpn % do_conv_micro_N
    dpn % do_varying_rpen    = cpn % do_varying_rpen
    dpn % rpen_choice        = cpn % rpen_choice
    dpn % do_subcloud_flx    = cpn % do_subcloud_flx
    dpn % do_new_subflx      = cpn % do_new_subflx
    dpn % use_lcl_only       = cpn % use_lcl_only
    dpn % do_new_pevap       = cpn % do_new_pevap
    dpn % stop_at_let        = cpn % stop_at_let
    dpn % do_limit_wmax      = cpn % do_limit_wmax
    dpn % do_limit_fdr       = cpn % do_limit_fdr
    dpn % plev_for           = cpn % plev_for
    dpn % plev_umf           = cpn % plev_umf
    dpn % do_detran_zero     = cpn % do_detran_zero
    dpn % N0                 = cpn % N0
    dpn % rle                = cpn % rle
    dpn % rpen               = cpn % rpen
    dpn % eis_max            = cpn % eis_max
    dpn % eis_min            = cpn % eis_min
    dpn % rmaxfrac           = cpn % rmaxfrac
    dpn % wmin               = cpn % wmin
    dpn % wmax               = cpn % wmax
    dpn % rbuoy              = cpn % rbuoy
    dpn % rdrag              = cpn % rdrag
    dpn % frac_drs           = cpn % frac_drs
    dpn % frac_dr0           = cpn % frac_dr0
    dpn % bigc               = cpn % bigc
    dpn % auto_th0           = cpn % auto_th0
    dpn % deltaqc0           = cpn % deltaqc0
    dpn % do_pdfpcp          = cpn % do_pdfpcp
    dpn % do_pmadjt          = cpn % do_pmadjt
    dpn % do_emmax           = cpn % do_emmax
    dpn % do_pnqv            = cpn % do_pnqv
    dpn % emfrac_max         = cpn % emfrac_max
    dpn % auto_rate          = cpn % auto_rate
    dpn % tcrit              = cpn % tcrit
    dpn % cldhgt_max         = cpn % cldhgt_max
    dpn % cldhgt_max_shallow = cpn % cldhgt_max_shallow
    dpn % do_ice             = cpn % do_ice
    dpn % do_ppen            = cpn % do_ppen
    dpn % do_pevap           = cpn % do_pevap
    dpn % hcevap             = cpn % hcevap
    dpn % hcevappbl          = cpn % hcevappbl
    dpn % cfrac              = cpn % cfrac
    dpn % pblfac             = cpn % pblfac
    dpn % do_minmse          = cpn % do_minmse
    dpn % mixing_assumption  = cpn % mixing_assumption
    dpn % mp_choice          = cpn % mp_choice
    dpn % de_choice          = cpn % de_choice
    dpn % Nl_land            = cpn % Nl_land
    dpn % Nl_ocean           = cpn % Nl_ocean
    dpn % qi_thresh          = cpn % qi_thresh
    dpn % r_thresh           = cpn % r_thresh
    dpn % peff_l             = cpn % peff_l
    dpn % peff_i             = cpn % peff_i
    dpn % do_forcedlifting   = cpn % do_forcedlifting
    dpn % atopevap           = cpn % atopevap
    dpn % wtwmin_ratio       = cpn % wtwmin_ratio
    dpn % do_auto_aero       = cpn % do_auto_aero
    dpn % rad_crit           = cpn % rad_crit
    dpn % wrel_min           = cpn % wrel_min
    dpn % do_weffect         = cpn % do_weffect
    dpn % do_new_pblfac      = cpn % do_new_pblfac
    dpn % weffect            = cpn % weffect
    dpn % use_online_aerosol = cpn % use_online_aerosol
    dpn % use_new_let        = cpn % use_new_let
    dpn % do_tten_max        = cpn % do_tten_max
    dpn % tten_max           = cpn % tten_max
    dpn % rkm_max            = cpn % rkm_max
    dpn % rkm_min            = cpn % rkm_min
    dpn % scaleh0            = cpn % scaleh0
    dpn % beta               = cpn % beta
    dpn % nom_ratio          = cpn % nom_ratio

  end subroutine cpn_copy

!#####################################################################
!#####################################################################

  subroutine dpconv0(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                     rkm_dp, cbmf_deep, cp1, ct1, ocode, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(inout)  :: dpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(in)     :: sd
    type(adicloud),  intent(in)     :: ac
    type(cclosure),  intent(in)     :: cc
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    type(cplume),    intent(inout)  :: cp
    type(ctend),     intent(inout)  :: ct
    real,            intent(inout)  :: rkm_dp, cbmf_deep
    type(cplume),    intent(inout)  :: cp1
    type(ctend),     intent(inout)  :: ct1
    real,            intent(inout)  :: ocode
    integer,            intent(out) :: ier
    character(len=256), intent(out) :: ermesg

    real          :: zcldtop, dcrh

    ier = 0
    ermesg = ' '
    if ( (ocode.ne.0 .and. ocode.ne.4) .or. (cbmf_deep.eq.0) ) then
       ocode=6; return
    end if

    zcldtop = sd%z(cp%ltop)
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call ct_clear_k(ct1);
    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf_deep, cc%wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
       ocode=6; return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
    end if

  end subroutine dpconv0


!#####################################################################
!#####################################################################

  subroutine dpconv1(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                     rkm_dp, cbmf_deep, sd1, ac1, cp1, ct1, ocode, dcapedm, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(inout)  :: dpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(in)     :: sd
    type(adicloud),  intent(inout)  :: ac
    type(cclosure),  intent(in)     :: cc
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    type(cplume),    intent(inout)  :: cp
    type(ctend),     intent(inout)  :: ct
    real,            intent(inout)  :: rkm_dp, cbmf_deep
    type(sounding),  intent(inout)  :: sd1
    type(adicloud),  intent(inout)  :: ac1
    type(cplume),    intent(inout)  :: cp1
    type(ctend),     intent(inout)  :: ct1
    real,            intent(inout)  :: ocode, dcapedm
    integer,            intent(out) :: ier
    character(len=256), intent(out) :: ermesg

    real          :: zcldtop, wrel, cbmf0, cbmf_max, tmp
    real          :: zs0, ps0, hl0, thc0
    integer       :: ksrc
    real          :: zsrc, psrc, thcsrc, qctsrc, hlsrc

    ier = 0
    ermesg = ' '
    zcldtop = 2000 !sd%z(cp%ltop)
    wrel = max(cc%wrel, 0.1)

   if (dpc%do_cgust_dp) then
       !if ((ac%cape.gt.ac%cin) .and. (sd%land.gt.0.5)) then
       if ((ac%cape.gt.0) .and. (sd%land.gt.0.5)) then
          dpn%do_forcedlifting = .true.
          dpn%do_ppen = .false.
       endif
   endif

!    if ( (ocode.ne.0 .and. ocode.ne.4) .or. (cbmf_deep.eq.0) ) then
    if ( cbmf_deep.eq.0 ) then
       ocode=6; return
    end if
    if (ac%cape .lt. dpc%cape_th) then
       cbmf_deep = 0.; ocode=7; return
    end if

    cbmf0 = 0.0001
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call ct_clear_k(ct1);
    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
       cbmf_deep = 0.; ocode=8; return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
    end if

    call sd_copy_k(sd, sd1)
    tmp      = 1.
    sd1 % t  = sd1 % t  + ct1%tten  * sd%delt * tmp
    sd1 % qv = sd1 % qv + ct1%qvten * sd%delt * tmp
    sd1 % ql = sd1 % ql + ct1%qlten * sd%delt * tmp
    sd1 % qi = sd1 % qi + ct1%qiten * sd%delt * tmp
    sd1 % qa = sd1 % qa + ct1%qaten * sd%delt * tmp
    sd1 % qn = sd1 % qn + ct1%qnten * sd%delt * tmp
    sd1 % u  = sd1 % u  + ct1%uten  * sd%delt * tmp
    sd1 % v  = sd1 % v  + ct1%vten  * sd%delt * tmp

    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)

    call ac_clear_k(ac1);

    ksrc  =sd1%ksrc
    zsrc  =sd1%zs (ksrc)
    psrc  =sd1%ps (ksrc)
    thcsrc=sd1%thc(ksrc)
    qctsrc=sd1%qct(ksrc)
    hlsrc =sd1%hl (ksrc)
    call adi_cloud_k(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd1, Uw_p, .false., do_ice, ac1)

    dcapedm=(ac%cape-ac1%cape)/cbmf0

    if (dcapedm .lt. dpc%dcapedm_th) then
       cbmf_deep=0.; ocode=9;
       call ct_clear_k(ct1);
       return
    else
       cbmf_deep= (ac%cape - dpc%cape_th) / dcapedm / (dpc%tau_dp/sd%delt)
    end if

    cbmf_max  = (sd%ps(0) - sd%ps(cp1%krel))*(dpc%frac_limit_d/sd%delt)/Grav
    cbmf_deep = max(min(cbmf_deep, cbmf_max), 0.)

    if(cbmf_deep.lt.1.e-10) then
       cbmf_deep=0.; ocode=10;
       call ct_clear_k(ct1);
       return
    end if

    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf_deep, wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
       cbmf_deep = 0.; ocode=11;
       call ct_clear_k(ct1);
       return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
    end if

  end subroutine dpconv1

!#####################################################################
!#####################################################################

  subroutine dpconv2(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                     rkm_dp, cbmf_deep, sd1, ac1, cp1, ct1, ocode, &
                     dcwfndm, dcapedm, cwfn, lat, lon, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(inout)  :: dpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(inout)  :: sd
    type(adicloud),  intent(inout)  :: ac
    type(cclosure),  intent(in)     :: cc
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    type(cplume),    intent(inout)  :: cp
    type(ctend),     intent(inout)  :: ct
    real,            intent(inout)  :: rkm_dp, cbmf_deep
    type(sounding),  intent(inout)  :: sd1
    type(adicloud),  intent(inout)  :: ac1
    type(cplume),    intent(inout)  :: cp1
    type(ctend),     intent(inout)  :: ct1
    real,            intent(inout)  :: ocode, dcwfndm, dcapedm
    real,            intent(out)    :: cwfn
    real,            intent(in)     :: lat, lon
    integer,            intent(out) :: ier
    character(len=256), intent(out) :: ermesg

    real          :: zcldtop, wrel, wrel0, cbmf0, cwfn_new, cwfn_th, cbmf_max, tmp, ufrc
    integer       :: ksrc, k, let, krel
    real          :: latx, lonx, lat1b, lat1e, lon1b, lon1e, lat2b, lat2e, lon2b, lon2e
    real          :: zsrc, psrc, thcsrc, qctsrc, hlsrc, lofactor, taudp

    ier = 0
    ermesg = ' '
    cwfn = 0.; dcwfndm=0.; dcapedm=0.; wrel0=0.1; lofactor=1.0;
    taudp = dpc%tau_dp;

    if ( cbmf_deep.eq.0 ) then
       ocode=6; return
    end if

    cwfn_th = dpc%cwfn_th;
    if (dpc%do_cgust_dp) then
       if (dpc%cgust_choice==0) then
          if (ac%cape.gt.dpc%cape_th .and. cc%wrel.le.0. .and.     &
              ac%cape.gt.dpc%cin_fact*ac%cin .and. sd%cgust.gt.sd%cgust0) then
              dpn%do_forcedlifting = .true.
              dpn%do_ppen  = .false.;
              lofactor     = dpc%lofactor_d
              rkm_dp       = rkm_dp        * lofactor
              dpn % peff_l = dpn % peff_l  / lofactor
              dpn % peff_i = dpn % peff_i  / lofactor
          endif
       else if (dpc%cgust_choice==1 .and. sd%land.gt.0.5) then
          if (ac%cape.gt.dpc%cape_th .and. cc%wrel.le.0. .and.     &
              ac%cape.gt.dpc%cin_fact*ac%cin .and. sd%cgust.gt.sd%cgust0) then
              dpn%do_forcedlifting = .true.
              dpn%do_ppen  = .false.;
              lofactor     = 1.- sd%land*(1.- dpc%lofactor_d)
              rkm_dp       = rkm_dp        * lofactor
              !dpn % peff_l = dpn % peff_l  / lofactor
              !dpn % peff_i = dpn % peff_i  / lofactor
          endif
       else if (dpc%cgust_choice==2 .and. sd%land.gt.0.5) then
          if (ac%cape.gt.dpc%cape_th .and. cc%wrel.le.0. .and.     &
              ac%cape.gt.dpc%cin_fact*ac%cin .and. sd%cgust.gt.sd%cgust0) then
              dpn%do_forcedlifting = .true.
              dpn%do_ppen  = .false.;
              lofactor     = 1.- sd%land*(1.- dpc%lofactor_d)
              rkm_dp       = rkm_dp        * lofactor
              dpn % peff_l = dpn % peff_l  * lofactor
              dpn % peff_i = dpn % peff_i  * lofactor
          endif
       else if (dpc%cgust_choice==3 .and. sd%land.gt.0.5) then
          if (ac%cape.gt.dpc%cape_th .and. cc%wrel.le.0. .and.     &
              ac%cape.gt.dpc%cin_fact*ac%cin .and. sd%cgust.gt.sd%cgust0) then
              dpn%do_forcedlifting = .true.
              dpn%do_ppen  = .false.;
              lofactor     = 1.- sd%land*(1.- dpc%lofactor_d)
              rkm_dp       = rkm_dp        * lofactor
              dpn % peff_l = dpn % peff_l  / lofactor
              dpn % peff_i = dpn % peff_i  / lofactor
          endif
       else if (dpc%cgust_choice==5) then !saved only for testing purpose
          lat1b= 30.; lat1e=40.; lon1b=260.; lon1e=270.
          lat2b=-20.; lat2e=10.; lon2b=285.; lon2e=305.
          latx=lat*180/3.1415926; lonx=lon*180/3.1415926
          if (latx.gt.lat1b .and. latx.lt.lat1e .and. lonx.gt.lon1b .and. lonx.lt.lon1e) then
             tmp=1
          elseif (latx.gt.lat2b .and. latx.lt.lat2e .and. lonx.gt.lon2b .and. lonx.lt.lon2e) then
             tmp=2
          endif
          sd%do_gust_qt = .true.
          sd%gqt_choice = 0
          call extend_sd_k(sd,sd%pblht, do_ice, Uw_p)
          ksrc  =sd%ksrc
          zsrc  =sd%zsrc
          psrc  =sd%psrc
          thcsrc=sd%thcsrc
          qctsrc=sd%qctsrc
          hlsrc =sd%hlsrc
          call ac_clear_k(ac);
          call adi_cloud_k(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, Uw_p, .false., do_ice, ac)
       endif
    end if

    if (cp%ltop == 0) then
       zcldtop = 1000.
    else
        zcldtop = sd%z(cp%ltop)
    endif

    wrel  = max(cc%wrel, wrel0)
    cbmf0 = dpc%cbmf0
    call cp_clear_k(cp);  cp %maxcldfrac=1.;
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call ct_clear_k(ct1);
    call cumulus_plume_k(dpn, sd, ac, cp, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp%ltop.lt.cp%krel+2 .or. cp%let.le.cp%krel+1) then
       cbmf_deep = 0.; ocode=8; return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp, ct1, do_coldT)
    end if

    call sd_copy_k(sd, sd1)
    tmp      = 1.
    sd1 % t  = sd1 % t  + ct1%tten  * sd%delt * tmp
    sd1 % qv = sd1 % qv + ct1%qvten * sd%delt * tmp
    sd1 % ql = sd1 % ql + ct1%qlten * sd%delt * tmp
    sd1 % qi = sd1 % qi + ct1%qiten * sd%delt * tmp
    sd1 % qa = sd1 % qa + ct1%qaten * sd%delt * tmp
    sd1 % qn = sd1 % qn + ct1%qnten * sd%delt * tmp
    sd1 % u  = sd1 % u  + ct1%uten  * sd%delt * tmp
    sd1 % v  = sd1 % v  + ct1%vten  * sd%delt * tmp

    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)

    ksrc  =sd1%ksrc
    zsrc  =sd1%zsrc
    psrc  =sd1%psrc
    thcsrc=sd1%thcsrc
    qctsrc=sd1%qctsrc
    hlsrc =sd1%hlsrc

    call ac_clear_k(ac1);
    call adi_cloud_k(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd1, Uw_p, .false., do_ice, ac1)

    call cumulus_plume_k(dpn, sd1, ac1, cp1, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)

    krel=max(cp%krel,cp1%krel)
    let =min(cp%let, cp1%let)

    cwfn=0.; !unit: m2/s2 or J/kg
    do k = krel,let-1
       cwfn = cwfn + cp%buo(k)/cp%thvtop(k)*cp%dp(k)/sd%rho(k)
    end do

    cwfn_new=0.;
    do k = krel,let-1
       cwfn_new = cwfn_new + cp1%buo(k)/cp1%thvtop(k)*cp1%dp(k)/sd1%rho(k)
    end do

    dcapedm=(ac%cape-ac1%cape)/cbmf0
    dcwfndm =(cwfn - cwfn_new)/cbmf0

    if (dcwfndm .lt. dpc%dcwfndm_th .or. cwfn.lt.0) then
       cbmf_deep=0.; ocode=9; return
    else
       cbmf_deep= (cwfn-cwfn_th) / dcwfndm / (taudp/sd%delt)
    end if

    tmp = sd%ps(0)-sd%ps(krel)

!    if (sd%ksrc .gt. sd%kinv) then !elevated convection
!       tmp = sd%dp(krel)
!    end if

    cbmf_max = tmp*dpc%frac_limit_d/sd%delt/Grav

    cbmf_deep = max(min(cbmf_deep, cbmf_max), 0.)

    if(cbmf_deep.lt.1.e-10) then
       cbmf_deep=0.; ocode=10; return
    end if

    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf_deep, wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
       cbmf_deep = 0.; ocode=11; return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
    end if

  end subroutine dpconv2
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
  subroutine dpconv3(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                     rkm_dp, cbmf_deep, sd1, ac1, cp1, ct1, ocode, &
                     dcwfndm, dcapedm, cwfn, lat, lon, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(inout)  :: dpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(inout)  :: sd
    type(adicloud),  intent(inout)  :: ac
    type(cclosure),  intent(in)     :: cc
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    type(cplume),    intent(inout)  :: cp
    type(ctend),     intent(inout)  :: ct
    real,            intent(inout)  :: rkm_dp, cbmf_deep
    type(sounding),  intent(inout)  :: sd1
    type(adicloud),  intent(inout)  :: ac1
    type(cplume),    intent(inout)  :: cp1
    type(ctend),     intent(inout)  :: ct1
    real,            intent(inout)  :: ocode, dcwfndm, dcapedm
    real,            intent(out)    :: cwfn
    real,            intent(in)     :: lat, lon
    integer,            intent(out) :: ier
    character(len=256), intent(out) :: ermesg

    real          :: zcldtop, wrel, wrel0, cbmf0, cwfn_new, cwfn_th, cbmf_max, tmp, ufrc
    integer       :: ksrc, k, let, krel
    real          :: latx, lonx, lat1b, lat1e, lon1b, lon1e, lat2b, lat2e, lon2b, lon2e
    real          :: zsrc, psrc, thcsrc, qctsrc, hlsrc, lofactor, cape, taudp

    ier = 0
    ermesg = ' '
    cwfn = 0.; dcwfndm=0.; dcapedm=0.; wrel0=0.1; lofactor=1.0;
    taudp = dpc%tau_dp;

    if ( cbmf_deep.eq.0 ) then
       ocode=6; return
    end if

    cwfn_th = dpc%cwfn_th;
    if (dpc%do_cgust_dp) then
       if (dpc%cgust_choice==1 .and. sd%land.gt.0.5) then
          if (ac%cape.gt.dpc%cape_th .and. cc%wrel.le.0. .and.     &
              ac%cape.gt.dpc%cin_fact*ac%cin .and. sd%cgust.gt.sd%cgust0) then
              dpn%do_forcedlifting = .true.
              dpn%do_ppen  = .false.;
              lofactor     = 1.- sd%land*(1.- dpc%lofactor_d)
              rkm_dp       = rkm_dp        * lofactor
              !dpn % peff_l = dpn % peff_l  / lofactor
              !dpn % peff_i = dpn % peff_i  / lofactor
          endif
       else if (dpc%cgust_choice==5) then !saved only for testing purpose
          lat1b= 30.; lat1e=40.; lon1b=260.; lon1e=270.
          lat2b=-20.; lat2e=10.; lon2b=285.; lon2e=305.
          latx=lat*180/3.1415926; lonx=lon*180/3.1415926
          if (latx.gt.lat1b .and. latx.lt.lat1e .and. lonx.gt.lon1b .and. lonx.lt.lon1e) then
             tmp=1
          elseif (latx.gt.lat2b .and. latx.lt.lat2e .and. lonx.gt.lon2b .and. lonx.lt.lon2e) then
             tmp=2
          endif
          sd%do_gust_qt = .true.
          sd%gqt_choice = 0
          call extend_sd_k(sd,sd%pblht, do_ice, Uw_p)
          ksrc  =sd%ksrc
          zsrc  =sd%zsrc
          psrc  =sd%psrc
          thcsrc=sd%thcsrc
          qctsrc=sd%qctsrc
          hlsrc =sd%hlsrc
          call ac_clear_k(ac);
          call adi_cloud_k(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, Uw_p, .false., do_ice, ac)
       endif
    end if

    zcldtop = sd%z(cp%ltop)
    wrel  = max(cc%wrel, wrel0)
    cbmf0 = dpc%cbmf0
    call cp_clear_k(cp);  cp %maxcldfrac=1.;
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call ct_clear_k(ct1);
    call cumulus_plume_k(dpn, sd, ac, cp, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp%ltop.lt.cp%krel+2 .or. cp%let.le.cp%krel+1) then
       cbmf_deep = 0.; ocode=8; return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp, ct1, do_coldT)
    end if

    call sd_copy_k(sd, sd1)
    tmp      = 1.
    sd1 % t  = sd1 % t  + ct1%tten  * sd%delt * tmp
    sd1 % qv = sd1 % qv + ct1%qvten * sd%delt * tmp
    sd1 % ql = sd1 % ql + ct1%qlten * sd%delt * tmp
    sd1 % qi = sd1 % qi + ct1%qiten * sd%delt * tmp
    sd1 % qa = sd1 % qa + ct1%qaten * sd%delt * tmp
    sd1 % qn = sd1 % qn + ct1%qnten * sd%delt * tmp
    sd1 % u  = sd1 % u  + ct1%uten  * sd%delt * tmp
    sd1 % v  = sd1 % v  + ct1%vten  * sd%delt * tmp

    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)

    ksrc  =sd1%ksrc
    zsrc  =sd1%zsrc
    psrc  =sd1%psrc
    thcsrc=sd1%thcsrc
    qctsrc=sd1%qctsrc
    hlsrc =sd1%hlsrc

    call ac_clear_k(ac1);
    call adi_cloud_k(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd1, Uw_p, .false., do_ice, ac1)

    dcapedm=(ac%cape-ac1%cape)/cbmf0

    if (sd%use_capecin_avg) then
       cape = sd%cape_avg
    else
       cape = ac%cape
    endif

    if (dcapedm .lt. dpc%dcapedm_th .or. cape.lt.0.) then
       cbmf_deep=0.; ocode=9; return
    else
       cbmf_deep= (cape - dpc%cape_th) / dcapedm / (taudp/sd%delt)
    end if

    krel=cp%krel
    tmp = sd%ps(0)-sd%ps(krel)
    cbmf_max = tmp*dpc%frac_limit_d/sd%delt/Grav
    cbmf_deep = max(min(cbmf_deep, cbmf_max), 0.)

    if(cbmf_deep.lt.1.e-10) then
       cbmf_deep=0.; ocode=10;
       return
    end if

    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf_deep, wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
       cbmf_deep = 0.; ocode=11;
       return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
    end if

  end subroutine dpconv3
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
  subroutine scale_ct (ct, sd, Uw_p, scale)
    implicit none
    type(ctend),     intent(inout) :: ct
    type(sounding),  intent(in)    :: sd
    type(uw_params), intent(inout) :: Uw_p
    real,            intent(in)    :: scale
    integer  :: i, k

    ct%uten   =ct%uten  *scale;
    ct%vten   =ct%vten  *scale;
    ct%tten   =ct%tten  *scale;
    ct%qvten  =ct%qvten *scale;
    ct%qlten  =ct%qlten *scale;
    ct%qiten  =ct%qiten *scale;
    ct%qaten  =ct%qaten *scale;
    ct%qnten  =ct%qnten *scale;
    ct%qldet  =ct%qldet *scale;
    ct%qidet  =ct%qidet *scale;
    ct%qadet  =ct%qadet *scale;
    ct%qndet  =ct%qndet *scale;
    ct%hlten  =ct%hlten *scale;
    ct%thcten =ct%thcten*scale;
    ct%qctten =ct%qctten*scale;
    ct%qvdiv  =ct%qvdiv *scale;
    ct%qldiv  =ct%qldiv *scale;
    ct%qidiv  =ct%qidiv *scale;
    ct%hlflx  =ct%hlflx *scale;
    ct%thcflx =ct%thcflx*scale;
    ct%qctflx =ct%qctflx*scale;
    ct%qvflx  =ct%qvflx *scale;
    ct%qlflx  =ct%qlflx *scale;
    ct%qiflx  =ct%qiflx *scale;
    ct%qaflx  =ct%qaflx *scale;
    ct%qnflx  =ct%qnflx *scale;
    ct%umflx  =ct%umflx *scale;
    ct%vmflx  =ct%vmflx *scale;
    ct%pflx   =ct%pflx  *scale;
    ct%tevap  =ct%tevap *scale;
    ct%qevap  =ct%qevap *scale;
    ct%trflx  =ct%trflx *scale;
    ct%trten  =ct%trten *scale;
    ct%trwet  =ct%trwet *scale;
    ct%trevp  =ct%trevp *scale;
    ct%dtring =ct%dtring*scale;
    ct%snow   =ct%snow  *scale;
    ct%rain   =ct%rain  *scale;
    ct%denth  =ct%denth *scale;
    ct%dting  =ct%dting *scale;
    ct%dqtmp  =ct%dqtmp *scale;

    ct%dting=0.;
    ct%denth=0.; ct%dqtmp=0.; ct%uav=0.;ct%vav=0.;
    do k = 1,sd%kmax
       ct%denth = ct%denth + (Uw_p%cp_air*ct%tten(k)-Uw_p%HLv*ct%qlten(k)-Uw_p%HLs*ct%qiten(k))*sd%dp(k)/Uw_p%grav
       ct%dqtmp = ct%dqtmp + (ct%qvten(k) + ct%qlten(k) + ct%qiten(k))*sd%dp(k)/Uw_p%grav
       ct%uav   = ct%uav   + ct%uten(k)*sd%dp(k)/Uw_p%grav
       ct%vav   = ct%vav   + ct%vten(k)*sd%dp(k)/Uw_p%grav
       ct%dting = ct%dting + Uw_p%cp_air*ct%tten(k)*sd%dp(k)/Uw_p%grav
    end do
    ct%denth = ct%denth - Uw_p%HLv*ct%rain - Uw_p%HLs*ct%snow
    ct%dqtmp = ct%dqtmp + ct%rain + ct%snow

    do i = 1,size(sd%tr,2)
       ct%dtring(i)=0.;
       do k = 1,sd%kmax
         ct%dtring(i) = ct%dtring(i) + Uw_p%cp_air*ct%trten(k,i)*sd%dp(k)/Uw_p%grav
       enddo
     end do

  end subroutine scale_ct

!#####################################################################
!#####################################################################
  subroutine compute_cwfn3d(dpc, dpn, Uw_p, sd, ac1, cp1, do_coldT, do_ice, &
                     rkm_dp, cwfn3d, cape3d, ier, ermesg)
    implicit none

    type(deepc),     intent(in)     :: dpc
    type(cpnlist),   intent(in)     :: dpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(inout)  :: sd
    type(adicloud),  intent(inout)  :: ac1
    type(cplume),    intent(inout)  :: cp1
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    real,            intent(in)     :: rkm_dp
    real,            intent(inout),  dimension(:)  :: cwfn3d, cape3d
    integer,            intent(inout) :: ier
    character(len=256), intent(inout) :: ermesg

    integer       :: ksrc, i, k, km
    real          :: zcldtop, wrel, cbmf0, tmp, zsrc, psrc, hlsrc, thcsrc, qctsrc,cwfnmax

    ier = 0
    ermesg = ' '
    zcldtop = 2000
    wrel = 0.1
    cbmf0=0.001
    km=15
    cape3d(:)=0.; cwfn3d(:)=0.;
    do k=1,km
       call ac_clear_k(ac1);
       call cp_clear_k(cp1); cp1%maxcldfrac=1.;
       ksrc=k
       zsrc  =sd%zs (ksrc)
       psrc  =sd%ps (ksrc)
       hlsrc =sd%hl (ksrc)
       thcsrc=sd%thc(ksrc)
       qctsrc=sd%qct(ksrc)
       call adi_cloud_k(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, Uw_p, .false., do_ice, ac1)
       cape3d(k)=ac1%cape
       call cumulus_plume_k(dpn, sd, ac1, cp1, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)
       cwfn3d(k)=0.
       do i = cp1%krel,cp1%let
          cwfn3d(k) = cwfn3d(k) + cp1%buo(i)/cp1%thvtop(i)*cp1%dp(i)
       end do
    end do

!    ksrc=1
!    cwfnmax=cwfn3d(1)
!    do k=1,km
!     if (cwfn3d(k) .gt. cwfnmax) then
!       cwfnmax=cwfn3d(k)
!       ksrc=k
!     end if
!    end do
!    sd%ksrc  =ksrc;

  end subroutine compute_cwfn3d
!#####################################################################
!#####################################################################

  subroutine conv_forced(dpc, dpn, Uw_p, sd, ac, do_coldT, do_ice, &
                     rkm_dp, cbmf_deep, cp1, ct1, lat, lon, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(inout)  :: dpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(in)     :: sd
    type(adicloud),  intent(in)     :: ac
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    real,            intent(inout)  :: rkm_dp, cbmf_deep
    type(cplume),    intent(inout)  :: cp1
    type(ctend),     intent(inout)  :: ct1
    real,            intent(in)     :: lat, lon
    integer,            intent(out) :: ier
    character(len=256), intent(out) :: ermesg

    real          :: latx, lonx, lat1b, lat1e, lon1b, lon1e, lat2b, lat2e, lon2b, lon2e
    real          :: zcldtop, wrel, ufrc, tmp
    integer       :: k

    ier = 0
    ermesg = ' '
    zcldtop = 2000

    if (sd%land.gt.0.5) then
      if (ac%cape.gt.dpc%cape_th .and. sd%cgust.gt.sd%cgust0) then
         lat1b= 30.; lat1e=40.; lon1b=260.; lon1e=270.
         lat2b=-20.; lat2e=10.; lon2b=285.; lon2e=305.
         latx=lat*180/3.1415926; lonx=lon*180/3.1415926
         if (latx.gt.lat1b .and. latx.lt.lat1e .and. lonx.gt.lon1b .and. lonx.lt.lon1e) then
            tmp=1
         elseif (latx.gt.lat2b .and. latx.lt.lat2e .and. lonx.gt.lon2b .and. lonx.lt.lon2e) then
            tmp=2
         endif
         dpn%do_forcedlifting = .true.
         dpn%do_ppen = .false. !dpn%do_pevap = .false.
         call cclosure_gust(sd%cgust, sd, Uw_p, ac, dpc%wcrit_min_gust, cbmf_deep, wrel, ufrc)
         if(cbmf_deep.lt.1.e-10) then 
            cbmf_deep=0.;
            call cp_clear_k(cp1); cp1%maxcldfrac=1.; cp1%cush=-1;
            call ct_clear_k(ct1);
            return
         end if
         call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf_deep, wrel, zcldtop, Uw_p, ier, ermesg)
         if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
            cbmf_deep = 0.; 
            call cp_clear_k(cp1); cp1%maxcldfrac=1.; cp1%cush=-1;
            call ct_clear_k(ct1);
            return
         else
            call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
            return
         end if
      endif
    endif
  end subroutine conv_forced
!#####################################################################
!#####################################################################


!#####################################################################
!#####################################################################
  subroutine cclosure_gust(cgust, sd, Uw_p, ac, wcrit_min, cbmf, wrel, ufrc)

    implicit none
    real,            intent(in)    :: cgust
    type(sounding),  intent(in)    :: sd
    type(uw_params), intent(in)    :: Uw_p
    type(adicloud),  intent(in)    :: ac
    real,            intent(in)    :: wcrit_min !wcrit_min=0.2
    real,            intent(out)   :: cbmf, wrel, ufrc

    integer :: kkk
    real    :: sigmaw, wcrit, erfarg, wexp, wtw, cbmf0, tmp

    cbmf=0.; wrel=0.; ufrc=0.;
    wcrit  = sqrt(2. * ac % cin)
    sigmaw = sqrt(cgust)
    wcrit = max(wcrit, wcrit_min*sigmaw)
    cbmf   = ac % rho0lcl * sigmaw / 2.5066 * exp(-0.5*((wcrit/sigmaw)**2.))

    !Diagnose updraft fraction sqrt(2.) = 1.4142
    erfarg=wcrit / (1.4142 * sigmaw)
    if(erfarg.lt.20.)then
       ufrc = min(0.15, 0.5*erfccc(erfarg))
    else
       ufrc = 0.
    endif

    if(ufrc.gt.0.0) then !Diagnose expected value of cloud base vertical velocity
       wexp = cbmf / ac % rho0lcl / ufrc
    else
       wexp = 0.
       cbmf = 0.
    endif

    wtw = wexp * wexp - 2 * ac % cin
    if(wtw.le.0.) then
       wrel=0.;
    else
       wrel=sqrt(wtw)
    end if
    wrel=min(wrel, 25.)

    cbmf0 = sd%dp(sd%ksrc) / sd%delt / Uw_p%GRAV
    kkk = max(max(ac%klcl, sd%kinv),1)
    tmp = (sd%ps(0) - sd%ps(kkk)) * 0.25
    cbmf0 = tmp / sd%delt / Uw_p%GRAV

    cbmf = min (cbmf, cbmf0)
    if (wrel .gt. 0.) then
      ufrc=cbmf / wrel /ac % rho0lcl
    else
      ufrc=0.
    end if

    return

  end subroutine cclosure_gust
!#####################################################################
!#####################################################################


end MODULE DEEP_CONV_MOD
