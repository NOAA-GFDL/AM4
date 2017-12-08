                  module moistproc_kernels_mod

use time_manager_mod,             only: time_type
use constants_mod,                only: CP_AIR, HLV, HLS, SECONDS_PER_DAY
use field_manager_mod,            only: MODEL_ATMOS
use tracer_manager_mod,           only: get_tracer_index
use moist_conv_mod,               only: moist_conv
use uw_conv_mod,                  only: uw_conv
use ras_mod,                      only: ras
use cu_mo_trans_mod,              only: cu_mo_trans
use aerosol_types_mod,            only: aerosol_type
use detr_ice_num_mod ,            only: detr_ice_num
use physics_radiation_exch_mod,   only: cloud_scheme_data_type

implicit none
private
public   moistproc_mca, moistproc_ras, moistproc_cmt, &
         moistproc_uw_conv, moistproc_scale_uw, moistproc_scale_donner


!--------------------- version number ----------------------------------
character(len=128) :: &
version = '$Id$'
character(len=128) :: tagname = '$Name$'

!-----------------------------------------------------------------------




                          contains



!#######################################################################

subroutine moistproc_cmt ( Time, is, js, t, u, v, tracer, pfull, phalf, &
                           zfull, zhalf, pmass, tdt, udt, vdt, rdt,     &
                           ttnd_conv, dt, mc_cmt, det_cmt, diff_cu_mo,  &
                           num_tracers)

type(time_type),          intent(in)    :: Time
integer,                  intent(in)    :: is, js, num_tracers
real,                     intent(in)    :: dt
real, dimension(:,:,:),   intent(in)    :: pfull, phalf, zfull, zhalf,  &
                                           pmass, mc_cmt, det_cmt
real, dimension(:,:,:),   intent(inout) :: t, u, v, tdt, udt, vdt,   &
                                           ttnd_conv, diff_cu_mo
real, dimension(:,:,:,:), intent(inout) :: rdt, tracer

      real, dimension(size(t,1), size(t,2), size(t,3)) :: ttnd, utnd, vtnd
      real, dimension(size(rdt,1), size(rdt,2),      &
                              size(rdt,3),num_tracers) :: qtr
      integer :: n

!------------------------------------------------------------------------
!    call cu_mo_trans to calculate cumulus momentum transport.
!------------------------------------------------------------------------
      call cu_mo_trans (is, js, Time, mc_cmt, t, phalf, pfull, &
                        zhalf, zfull, dt, u, v, tracer,        &
                        pmass, det_cmt, utnd, vtnd, ttnd,      &
                        qtr, diff_cu_mo  )

!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    just obtained from cu_mo_trans.
!---------------------------------------------------------------------
      do n=1, num_tracers
        rdt(:,:,:,n) = rdt(:,:,:,n) + qtr(:,:,:,n)
      end do

!----------------------------------------------------------------------
!    add the temperature, specific humidity and momentum tendencies 
!    due to cumulus transfer (ttnd, qtnd, utnd, vtnd) to the arrays 
!    accumulating these tendencies from all physics processes (tdt, qdt, 
!    udt, vdt).
!----------------------------------------------------------------------
      tdt = tdt + ttnd 
      udt = udt + utnd
      vdt = vdt + vtnd
      ttnd_conv = ttnd_conv + ttnd

!-----------------------------------------------------------------------


end subroutine moistproc_cmt



!#######################################################################

subroutine moistproc_mca      &
         (Time, is, js, t, q, tracer, pfull, phalf, coldT, dtinv, &
          tdt, qdt, rdt, q_tnd, ttnd_conv, qtnd_conv, lprec, fprec,  &
          doing_prog_clouds, num_tracers, tracers_in_mca, num_mca_tracers)

type(time_type),          intent(in)    :: Time
integer,                  intent(in)    :: is, js, num_tracers,  &
                                           num_mca_tracers
real,                     intent(in)    :: dtinv
logical,                  intent(in)    :: doing_prog_clouds
logical, dimension(:),    intent(in)    :: tracers_in_mca
logical, dimension(:,:),  intent(in)    :: coldT
real, dimension(:,:,:),   intent(in)    :: pfull, phalf 
real, dimension(:,:),     intent(inout) :: lprec, fprec
real, dimension(:,:,:),   intent(inout) :: t, q, tdt, qdt, ttnd_conv,  &
                                           qtnd_conv
real, dimension(:,:,:,:), intent(inout) :: rdt, q_tnd
real, dimension(:,:,:,:), intent(out)   :: tracer

      real, dimension(size(t,1), size(t,2))      :: rain, snow
      real, dimension(size(t,1), size(t,2), &
                                 size(t,3))      :: ttnd, qtnd
      real, dimension(size(rdt,1), size(rdt,2),&
               size(rdt,3), num_mca_tracers)     :: trcr, qtr

      integer :: nn, n, nql, nqa, nqi, nqn

      nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
      nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
      nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
      nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )

!---------------------------------------------------------------------
!    check each active tracer to find any that are to be transported 
!    by moist convective adjustment and fill the mca_tracers array with
!    these fields.
!---------------------------------------------------------------------
      nn = 1
      do n=1, num_tracers
        if (tracers_in_mca(n)) then
          trcr(:,:,:,nn) = tracer(:,:,:,n)
          nn = nn + 1
        endif
      end do

!---------------------------------------------------------------------
!    call subroutine moist_conv to obtain the temperature, moisture
!    precipitation and tracer tendencies due to the moist convective
!    adjustment parameterization. currently there is no tracer tendency
!    due to this parameterization.
!---------------------------------------------------------------------
!++++yim Should also account for change in qn dut to moist convective adjustment.

      if (doing_prog_clouds) then
        call moist_conv (t, q, pfull, phalf, coldT, ttnd, qtnd,          &
                         rain, snow, dtinv, Time, is, js, trcr, qtr,     &
                         ql=tracer(:,:,:,nql), qi=tracer(:,:,:,nqi),     &
                         cf=tracer(:,:,:,nqa), qldel=q_tnd(:,:,:,nql),   &
                         qidel=q_tnd(:,:,:,nqi), cfdel=q_tnd(:,:,:,nqa))
      else
        call moist_conv (t, q, pfull, phalf, coldT, ttnd, qtnd,          &
                         rain, snow, dtinv, Time, is, js,                &
                         trcr, qtr )
      endif

!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    just obtained from moist convective adjustment. currently there
!    is no tracer transport by this process.
!    NOTE : the stratcloud tracers are updated within moist_conv.
!---------------------------------------------------------------------
      nn = 1
      do n=1, num_tracers
        if (tracers_in_mca(n)) then
          rdt(:,:,:,n) = rdt(:,:,:,n) + qtr(:,:,:,nn)
          nn = nn + 1
        endif
      end do

!----------------------------------------------------------------------
!    add the temperature and specific humidity tendencies from moist
!    convective adjustment (ttnd, qtnd) to the arrays accumulating 
!    these tendencies from all physics processes (tdt, qdt).
!----------------------------------------------------------------------
      tdt = tdt + ttnd 
      qdt = qdt + qtnd
      ttnd_conv = ttnd_conv + ttnd
      qtnd_conv = qtnd_conv + qtnd

!----------------------------------------------------------------------
!    increment the liquid, solid and total precipitation fields with 
!    the contribution from moist convective adjustment.
!----------------------------------------------------------------------
      lprec  = lprec  + rain
      fprec  = fprec  + snow

!----------------------------------------------------------------------
!    if prognostic clouds are active, add the cloud liquid, ice and area
!    tendencies from moist convective adjustment to the 
!    arrays accumulating these tendencies from all physics processes 
!    (rdt).
!----------------------------------------------------------------------
      if (doing_prog_clouds) then
        rdt(:,:,:,nql) = rdt(:,:,:,nql) + q_tnd(:,:,:,nql)
        rdt(:,:,:,nqi) = rdt(:,:,:,nqi) + q_tnd(:,:,:,nqi)
        rdt(:,:,:,nqa) = rdt(:,:,:,nqa) + q_tnd(:,:,:,nqa)
      endif


end subroutine moistproc_mca


!#######################################################################

subroutine moistproc_ras      &
              (Time, is, js, dt, coldT, t, q, u, v, tracer, pfull, phalf, &
               zhalf, tdt, qdt, udt, vdt, rdt, q_tnd, ttnd, qtnd,    &
               ttnd_conv, qtnd_conv, mc, det0, lprec, fprec, rain, snow,  &
               rain3d, snow3d, Aerosol, doing_prog_clouds, do_liq_num,   &
               num_tracers, tracers_in_ras, num_ras_tracers,             & 
               do_ice_num, detrain_ice_num)

type(time_type),          intent(in)           :: Time
integer,                  intent(in)           :: is, js, num_tracers,  &
                                                  num_ras_tracers
logical,                  intent(in)           :: doing_prog_clouds,   &
                                                  do_liq_num, do_ice_num,&
                                                  detrain_ice_num
real,                     intent(in)           :: dt
logical, dimension(:),    intent(in)           :: tracers_in_ras
logical, dimension(:,:),  intent(in)           :: coldT
real, dimension(:,:,:),   intent(in)           :: pfull, phalf, zhalf
real, dimension(:,:),     intent(inout)        :: lprec, fprec
real, dimension(:,:,:),   intent(inout)        :: t, q, u, v, tdt, qdt,   &
                                                  udt, vdt, ttnd, qtnd,   &
                                                  ttnd_conv, qtnd_conv
real, dimension(:,:,:,:), intent(inout)        :: rdt, tracer, q_tnd
real, dimension(:,:),     intent(out)          :: rain, snow
real, dimension(:,:,:),   intent(out)          :: rain3d,  snow3d, mc, det0

type(aerosol_type),       intent(in), optional :: Aerosol


      real, dimension(size(t,1), size(t,2),    &
                                 size(t,3))           :: utnd, vtnd
      real, dimension(size(rdt,1), size(rdt,2),  &
                         size(rdt,3),num_ras_tracers) :: trcr, qtr

      integer :: nn, n, nql, nqa, nqi, nqn, nqni

!----------------------------------------------------------------------
!    define tracer indices of the prognostic cloud tracers.
!----------------------------------------------------------------------
      nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
      nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
      nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
      nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
      nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )

!----------------------------------------------------------------------
!    if any tracers are to be transported by ras convection, check each
!    active tracer to find those to be transported and fill the 
!    ras_tracers array with these fields.
!---------------------------------------------------------------------
      nn = 1
      do n=1, num_tracers
        if (tracers_in_ras(n)) then
          trcr(:,:,:,nn) = tracer(:,:,:,n)
          nn = nn + 1
        endif
      end do

!----------------------------------------------------------------------
!    call subroutine ras to obtain the temperature, specific humidity,
!    velocity, precipitation and tracer tendencies and mass flux 
!    associated with the relaxed arakawa-schubert parameterization.
!----------------------------------------------------------------------
      if (doing_prog_clouds .and. (.not.do_liq_num)) then
        call ras (is, js, Time, t, q, u, v, pfull, phalf, zhalf, coldT, &
                  dt, ttnd, qtnd, utnd, vtnd, rain3d, snow3d, rain, snow, &
                  trcr, qtr, mc0=mc, det0=det0,        &
                  ql0=tracer(:,:,:,nql), qi0=tracer(:,:,:,nqi),    &
                  qa0=tracer(:,:,:,nqa), dl0=q_tnd(:,:,:,nql),     &
                  di0=q_tnd(:,:,:,nqi), da0=q_tnd(:,:,:,nqa))       

      elseif (doing_prog_clouds .and. do_liq_num) then
        call ras (is, js, Time, t, q, u, v, pfull, phalf, zhalf, coldT, &
                  dt, ttnd, qtnd, utnd, vtnd, rain3d, snow3d, rain, snow, &
                  trcr, qtr, mc0=mc, det0=det0,        &
                  ql0=tracer(:,:,:,nql), qi0=tracer(:,:,:,nqi),    &
                  qa0=tracer(:,:,:,nqa), dl0=q_tnd(:,:,:,nql),     &
                  di0=q_tnd(:,:,:,nqi), da0=q_tnd(:,:,:,nqa),  &      
                  qn0=tracer(:,:,:,nqn), dn0=q_tnd(:,:,:,nqn),     &
                  do_strat=doing_prog_clouds, Aerosol=Aerosol)
      else
        call ras (is, js, Time, t, q, u, v, pfull, phalf, zhalf, coldT, &
                  dt, ttnd, qtnd, utnd, vtnd, rain3d, snow3d, rain, snow, &
                  trcr, qtr, mc0=mc, det0=det0)
      endif

!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    just obtained from ras transport.
!    NOTE : the prognostic cloud tracers are updated within ras.        
!---------------------------------------------------------------------
      nn = 1
      do n=1, num_tracers
        if (tracers_in_ras(n)) then
          rdt(:,:,:,n) = rdt(:,:,:,n) + qtr (:,:,:,nn)
          nn = nn + 1
        endif
      end do

!----------------------------------------------------------------------
!    add the temperature, specific humidity and momentum tendencies 
!    from ras (ttnd, qtnd, utnd, vtnd) to the arrays accumulating 
!    these tendencies from all physics processes (tdt, qdt, udt, vdt).
!----------------------------------------------------------------------
      tdt = tdt + ttnd 
      qdt = qdt + qtnd
      udt = udt + utnd
      vdt = vdt + vtnd

!---------------------------------------------------------------------
!    update the total time tendency of temperature and specific humidity due
!    to convective processes. 
!---------------------------------------------------------------------
      ttnd_conv = ttnd_conv + ttnd
      qtnd_conv = qtnd_conv + qtnd

!----------------------------------------------------------------------
!    if prognostic clouds are activated, add the cloud liquid, ice, area
!    and droplet number tendencies from ras to the arrays accumulating 
!    these tendencies from all physics processes (rdt).
!----------------------------------------------------------------------
      if (doing_prog_clouds) then
        rdt(:,:,:,nql) = rdt(:,:,:,nql) + q_tnd(:,:,:,nql)
        rdt(:,:,:,nqi) = rdt(:,:,:,nqi) + q_tnd(:,:,:,nqi)
        rdt(:,:,:,nqa) = rdt(:,:,:,nqa) + q_tnd(:,:,:,nqa)
        if (do_liq_num) rdt(:,:,:,nqn) = rdt(:,:,:,nqn) + q_tnd(:,:,:,nqn)

!------------------------------------------------------------------------
!    if prognostic ice number is activated, call detr_ice_num to update
!    its value after ice detrainment is calculated (proportional to ice 
!    mass).
!------------------------------------------------------------------------
        IF (do_ice_num .AND. detrain_ice_num) THEN
          CALL detr_ice_num (t, q_tnd(:,:,:,nqi), q_tnd(:,:,:,nqni))   
          rdt(:,:,:,nqni) = rdt(:,:,:,nqni) + q_tnd(:,:,:,nqni)
        END IF
      endif

!----------------------------------------------------------------------
!    increment the liquid and frozen precipitation fields with 
!    the contribution from ras.
!----------------------------------------------------------------------
      lprec  = lprec  + rain
      fprec  = fprec  + snow

!------------------------------------------------------------------------



end subroutine moistproc_ras



!#######################################################################

subroutine moistproc_uw_conv  (           &
              Time, is, ie, js, je, dt, t, q, u, v, tracer,            &
              pfull, phalf, zfull, zhalf, omega, pblht,        &
              ustar, bstar, qstar, shflx, lhflx, land, coldT, Aerosol, &
              tdt_rad, tdt_dyn, qdt_dyn, dgz_dyn, ddp_dyn, tdt_dif,   &
              qdt_dif, hmint, lat, lon, cush, cbmf, cgust, tke, pblhto,  &
              rkmo, taudpo, exist_shconv, exist_dpconv,   & 
              pblht_prev, hlsrc_prev, qtsrc_prev, cape_prev, cin_prev,  &
              tke_prev, &!miz
              cmf, conv_calc_completed,            &
              available_cf_for_uw, tdt, qdt, udt, vdt, rdt,    &
              ttnd_conv, qtnd_conv, lprec, fprec, precip,      &
              liq_precflx, ice_precflx, rain_uw, snow_uw, ttnd_uw,   &
              qtnd_uw, utnd_uw, vtnd_uw, qtruw, qltnd_uw, qitnd_uw, &
              qatnd_uw, qntnd_uw, qnitnd_uw, doing_prog_clouds,   &
              do_limit_uw, do_liq_num, num_tracers, tracers_in_uw,  &
              num_uw_tracers, Cld_props, uw_wetdep, do_ice_num,  &
              detrain_ice_num)

type(time_type),              intent(in)    :: Time
type(aerosol_type),           intent(in)    :: Aerosol
integer,                      intent(in)    :: is, ie,js, je, num_tracers,&
                                               num_uw_tracers
real,                         intent(in)    :: dt
logical,                      intent(in)    :: doing_prog_clouds,  &
                                               do_limit_uw, do_liq_num
logical, dimension(:),        intent(in)    :: tracers_in_uw
logical, dimension(:,:),      intent(in)    :: coldT, conv_calc_completed
real, dimension(:,:),         intent(in)    :: land, ustar, bstar, qstar, &
                                               pblht, shflx, lhflx, lat, &
                                               lon
real, dimension(:,:,:),       intent(in)    :: pfull, phalf, zfull, zhalf, &
                                               omega, t, q, u, v,  &
                                               available_cf_for_uw
real, dimension(:,:,:,:),     intent(in)    :: tracer
real, dimension(:,:),         intent(inout) :: lprec, fprec, precip, cush, &
                                               cbmf, hmint, cgust
real, dimension(:,:),         intent(inout) :: tke, pblhto, rkmo, taudpo
integer, dimension(:,:,:),    intent(inout) :: exist_shconv, exist_dpconv
real, dimension(:,:,:),       intent(inout) :: pblht_prev, hlsrc_prev,  &
                                               qtsrc_prev, cape_prev,   &
                                               cin_prev, tke_prev !miz
real, dimension(:,:,:),       intent(in)    :: tdt_rad, tdt_dyn, qdt_dyn, &
                                               dgz_dyn, ddp_dyn, tdt_dif, &
                                               qdt_dif 
real, dimension(:,:,:),       intent(inout) :: tdt, qdt, udt, vdt,   &
                                               ttnd_conv, qtnd_conv, cmf
real, dimension(:,:,:,:),     intent(inout) :: rdt
type(cloud_scheme_data_type), intent(inout) :: Cld_props
logical, intent(in )                        :: do_ice_num, detrain_ice_num
real, dimension(:,:,:),       intent(out)   :: liq_precflx, ice_precflx
real, dimension(:,:,:),       intent(out)   :: uw_wetdep
real, dimension(:,:),         intent(inout) :: rain_uw, snow_uw
real, dimension(:,:,:),       intent(inout) :: ttnd_uw, qtnd_uw, utnd_uw, &
                                               vtnd_uw, qltnd_uw, qitnd_uw,&
                                               qatnd_uw, qntnd_uw, qnitnd_uw
real, dimension(:,:,:,:),     intent(inout) :: qtruw
                                           

      real, dimension(size(rdt,1), size(rdt,2),      &
                    size(rdt,3), num_uw_tracers) :: trcr

      integer :: n, nn, nql, nqi, nqa, nqn, nqni

!------------------------------------------------------------------------
!    define indices of cloud tracers into tracer array.
!------------------------------------------------------------------------
      nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
      nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
      nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
      nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
      nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )

!----------------------------------------------------------------------
!    if any tracers are to be transported by UW convection, check each
!    active tracer to find those to be transported and fill the 
!    trcr array with these fields.
!---------------------------------------------------------------------
      nn = 1
      do n=1, num_tracers
        if (tracers_in_uw(n)) then
          trcr(:,:,:,nn) = tracer(:,:,:,n)
          nn = nn + 1
        endif
      end do

!-------------------------------------------------------------------------
!    call uw_conv to calculate the effects of shallow convection.
!-------------------------------------------------------------------------
      call uw_conv (is, js, Time, t, q, u, v, pfull, phalf, zfull, zhalf, &
                    tracer, omega, dt, pblht, ustar, bstar, qstar, land,  &
                    coldT, Aerosol,lat, lon, cush, tke, doing_prog_clouds,&
                    conv_calc_completed, available_cf_for_uw, ttnd_uw,    &
                    qtnd_uw, qltnd_uw, qitnd_uw, qatnd_uw, qntnd_uw,      &
                    utnd_uw, vtnd_uw, rain_uw, snow_uw, cmf, liq_precflx, &
                    ice_precflx, Cld_props%liquid_amt, Cld_props%ice_amt, &
                    Cld_props%cloud_area, Cld_props%droplet_number,       &
                    trcr, qtruw, uw_wetdep, cbmf, cgust)

!-------------------------------------------------------------------------
!    call detr_ice_num to calculate the ice number tendency due to 
!    detrainment, which is proportional to the ice mass.
!-------------------------------------------------------------------------
      IF ( do_ice_num .AND. detrain_ice_num ) THEN
        CALL detr_ice_num (t, qitnd_uw(:,:,:),  &
                                              qnitnd_uw(:,:,:) )
      END IF

!-----------------------------------------------------------------------
!    if the subroutine to enforce non-negative cloud tracers is not being
!    used, update the physics tendency, convection tendency and tracer
!    tendency arrays with the contributions from uw_shallow convection.
!    if such enforcement is desired, the updates will occur later (this
!    should be further examined).
!-----------------------------------------------------------------------
      if (.not. do_limit_uw) then
        tdt = tdt + ttnd_uw
        qdt = qdt + qtnd_uw
        udt = udt + utnd_uw
        vdt = vdt + vtnd_uw
        ttnd_conv = ttnd_conv + ttnd_uw
        qtnd_conv = qtnd_conv + qtnd_uw
        lprec = lprec + rain_uw
        fprec = fprec + snow_uw
        precip = precip + rain_uw + snow_uw

        if (doing_prog_clouds) then
          rdt(:,:,:,nql) = rdt(:,:,:,nql) + qltnd_uw
          rdt(:,:,:,nqi) = rdt(:,:,:,nqi) + qitnd_uw
          rdt(:,:,:,nqa) = rdt(:,:,:,nqa) + qatnd_uw
          if (do_liq_num) rdt(:,:,:,nqn) = rdt(:,:,:,nqn) + qntnd_uw
          if (do_ice_num) rdt(:,:,:,nqni) = rdt(:,:,:,nqni) + qnitnd_uw
        endif

!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions
!    just obtained from uw transport.
!---------------------------------------------------------------------
        nn = 1
        do n=1, num_tracers
          if (tracers_in_uw(n)) then
            rdt(:,:,:,n) = rdt(:,:,:,n) + qtruw(:,:,:,nn)
            nn = nn + 1
          endif
        end do
      endif  !(.not. do_limit_uw)

!------------------------------------------------------------------------


end subroutine moistproc_uw_conv





!#######################################################################

subroutine moistproc_scale_donner     &
           (is, ie, js, je, dt, q, delta_temp, delta_q, precip_returned,  &
            total_precip, lheat_precip, liquid_precip, frozen_precip,   &
            pmass, num_tracers, tracers_in_donner, delta_ql, delta_qi,  &
            delta_qa, qlin, qiin, qtr, scale)

!--------------------------------------------------------------------------
!     tendencies coming out of Donner deep are adjusted to prevent
!     the formation of negative water vapor, liquid or ice.
!--------------------------------------------------------------------------

!-------------------------------------------------------------------------
integer,                  intent(in)    :: is, ie, js, je, num_tracers
real ,                    intent(in)    :: dt  ! Timestep of model.
logical, dimension(:),    intent(in)    :: tracers_in_donner
real, dimension(:,:),     intent(inout) :: precip_returned, total_precip, &
                                           lheat_precip
real, dimension(:,:,:),   intent(inout) :: q, delta_temp, delta_q,    &
                                           liquid_precip, frozen_precip, &
                                           pmass, delta_ql, delta_qi,   &
                                           delta_qa, qlin, qiin
real, dimension(:,:,:,:), intent(inout) :: qtr
real, dimension(:,:),     intent(out)   :: scale

!-----------------------------------------------------------------------

      real, dimension(size(q,1), size(q,2), size(q,3)) :: temp

      integer :: n, nn, i, j, k, ix, jx, kx
      real    :: qvin, dqv

!------------------------------------------------------------------------
!    define array sizes.
!------------------------------------------------------------------------
      ix = size(q,1)
      jx = size(q,2)
      kx = size(q,3)

!--------------------------------------------------------------------------
!     (1) Prevent negative liquid and ice specific humidities after
!     tendencies are applied
!--------------------------------------------------------------------------

      where ((qlin + delta_ql) .lt. 0.)
        delta_temp  = delta_temp - (qlin + delta_ql)*HLV/CP_AIR
        delta_q     = delta_q + (qlin + delta_ql)
        delta_ql    = delta_ql - (qlin + delta_ql)
      end where

      where ((qiin + delta_qi) .lt. 0.)
        delta_temp  = delta_temp - (qiin + delta_qi)*HLS/CP_AIR
        delta_q     = delta_q + (qiin + delta_qi)
        delta_qi    = delta_qi - (qiin + delta_qi)
      end where

!------------------------------------------------------------------------
!RSH NOTE: 
!    the 1.e-10 which are used below could be changed to be qmin, as 
!    defined in the physics_driver nml, and available as Exch_ctrl%qmin.
!------------------------------------------------------------------------
      where (abs(delta_ql + delta_qi) .lt. 1.e-10 )
        delta_qa = 0.0
      end where

!------------------------------------------------------------------------
!    (2) Compute limit on Donner tendencies to prevent water vapor
!    from going below 1.e-10, again the value of qmin as defined in
!    physics_driver_nml. 
!------------------------------------------------------------------------

!-----------------------------------------------------------------------
!    compute a scaling factor for each grid point.
!-----------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix
            qvin = q(i,j,k)
            dqv  = delta_q(i,j,k)
            if ( dqv.lt.0 .and. qvin + dqv .lt. 1.e-10 ) then
              temp(i,j,k) = max( 0.0, -(qvin - 1.e-10)/dqv )
            else
              temp(i,j,k) = 1.0
            endif
          end do
        end do
      end do

!-------------------------------------------------------------------------
!    define the scaling factor for each column as the minimum value found
!    within that column.
!-------------------------------------------------------------------------
      scale = minval( temp, dim=3 )

!-------------------------------------------------------------------------
!    scale the convective tendencies of temperature, cloud tracers and
!    tracers transported by the donner convection scheme.
!-------------------------------------------------------------------------
      do k=1,kx
        delta_temp(:,:,k)  = scale*delta_temp(:,:,k)
        delta_q(:,:,k)     = scale*delta_q (:,:,k)
        delta_qa(:,:,k)    = scale*delta_qa(:,:,k)
        delta_ql(:,:,k)    = scale*delta_ql(:,:,k)
        delta_qi(:,:,k)    = scale*delta_qi(:,:,k)
      end do

      nn = 1
      do n=1, num_tracers
        if (tracers_in_donner(n)) then
          do k=1,kx
            qtr(:,:,k,nn) = scale(:,:) * qtr(:,:,k,nn)
          end do
          nn = nn + 1
        endif
      end do

!-------------------------------------------------------------------------
!    scale the precipitation fields and associated enthalpy terms.
!    Precip returned from Donner scheme is recalculated below.
!-------------------------------------------------------------------------
!     precip_returned = scale*precip_returned

      total_precip = scale*total_precip
      lheat_precip = scale*lheat_precip
      do k=1, kx
        liquid_precip(:,:,k) = scale(:,:)*liquid_precip(:,:,k)
        frozen_precip(:,:,k) = scale(:,:)*frozen_precip(:,:,k)
      end do

!-------------------------------------------------------------------------
!    Limit liquid and frozen precip to not have negative values.
!-------------------------------------------------------------------------
      where ( liquid_precip(:,:,:) .lt. 0.)
        liquid_precip(:,:,:) = 0.0
      end where

      where ( frozen_precip(:,:,:) .lt. 0.)
        frozen_precip(:,:,:) = 0.0
      end where

!-------------------------------------------------------------------------
!    dimensions of liquid_precip is [kg(H20)/(kg s)]*(SECONDS_PER_DAY)
!    dimensions of frozen_precip is [kg(H20)/(kg s)]*(SECONDS_PER_DAY)

!    Note that (dt/seconds_per_day) * sum of (liquid_precip(k) + 
!    frozen_precip( k) *pmass(k)) gives precip_returned.
!-------------------------------------------------------------------------
      precip_returned(:,:) = 0.0
      do k=1, kx
        precip_returned(:,:) = precip_returned(:,:) +   &
             (liquid_precip(:,:,k) + frozen_precip(:,:,k))*pmass(:,:,k) *&
                                                   dt/SECONDS_PER_DAY
      end do

!-------------------------------------------------------------------------


end subroutine moistproc_scale_donner



!#######################################################################

subroutine moistproc_scale_uw   &
        (is, ie, js, je, dt, q, tracer, tdt, qdt, udt, vdt, rdt,    &
         ttnd_conv, qtnd_conv, lprec, fprec, precip, qtruw, rain_uw,  &
         snow_uw, ttnd_uw, qtnd_uw, utnd_uw, vtnd_uw, qltnd_uw, qitnd_uw, &
         qatnd_uw, qntnd_uw, qnitnd_uw, doing_prog_clouds, do_liq_num,  &
         num_tracers, tracers_in_uw, scale, do_ice_num)

!------------------------------------------------------------------------
!    tendencies coming out of UW shallow are adjusted to prevent
!    the formation of negative water vapor, liquid or ice.
!------------------------------------------------------------------------

integer,                  intent(in)    :: is, ie, js, je, num_tracers
real,                     intent(in)    :: dt
logical,                  intent(in)    :: doing_prog_clouds, do_liq_num,  &
                                           do_ice_num
logical, dimension(:),    intent(in)    :: tracers_in_uw
real, dimension(:,:,:),   intent(in)    :: q
real, dimension(:,:,:,:), intent(in)    :: tracer
real, dimension(:,:),     intent(inout) :: lprec, fprec, precip
real, dimension(:,:,:),   intent(inout) :: tdt, qdt, udt, vdt,  &
                                             ttnd_conv, qtnd_conv
real, dimension(:,:,:,:), intent(inout) :: rdt
real, dimension(:,:),     intent(out)   :: scale

real, dimension(:,:),     intent(inout) :: rain_uw, snow_uw
real, dimension(:,:,:),   intent(inout) :: ttnd_uw, qtnd_uw, utnd_uw,  &
                                           vtnd_uw, qltnd_uw, qitnd_uw,  &
                                           qatnd_uw, qntnd_uw, qnitnd_uw
real, dimension(:,:,:,:), intent(inout) :: qtruw

!------------------------------------------------------------------------
      real, dimension(size(q,1), size(q,2), size(q,3)) :: temp

      integer :: n, nn, i, j, k, ix, jx, kx, nql, nqi, nqa, nqn, nqni
      real    :: qvin, dqv

!-----------------------------------------------------------------------
!    define array dimensions and the tracer indices for the active cloud 
!    tracers.
!-----------------------------------------------------------------------
      ix = size(q,1)
      jx = size(q,2)
      kx = size(q,3)
      nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
      nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
      nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
      nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
      nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )

!----------------------------------------------------------------------
!    (1) Prevent negative liquid and ice specific humidities after 
!        tendencies are applied.
!-----------------------------------------------------------------------
      temp = tracer(:,:,:,nql)/dt + qltnd_uw
      where (temp(:,:,:) .lt. 0.)
        ttnd_uw  = ttnd_uw  - temp*HLV/CP_AIR
        qtnd_uw  = qtnd_uw  + temp
        qltnd_uw = qltnd_uw - temp
      end where

      temp = tracer(:,:,:,nqi)/dt + qitnd_uw
      where (temp .lt. 0.)
        ttnd_uw  = ttnd_uw  - temp*HLS/CP_AIR
        qtnd_uw  = qtnd_uw  + temp
        qitnd_uw = qitnd_uw - temp
      end where

!------------------------------------------------------------------------
!RSH NOTE: 
!    the 1.e-10 which are used below could be changed to be qmin, as 
!    defined in the physics_driver nml, and available as Exch_ctrl%qmin.
!------------------------------------------------------------------------
      where (abs(qltnd_uw + qitnd_uw)*dt .lt. 1.e-10 )
        qatnd_uw = 0.0
      end where

!------------------------------------------------------------------------
!    (2) Compute limit on UW tendencies to prevent water vapor
!    from going below 1.e-10. The value of 1.e-10 is consistent with qmin
!    defined in physics_driver_nml, and available as Exch_ctrl%qmin.
!------------------------------------------------------------------------

!-----------------------------------------------------------------------
!    compute a scaling factor for each grid point.
!-----------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix
            qvin = q(i,j,k) + tracer(i,j,k,nql) + tracer(i,j,k,nqi)
            dqv  = ( qtnd_uw(i,j,k) + qltnd_uw(i,j,k) + qitnd_uw(i,j,k) )*dt
            if ( dqv.lt.0 .and. qvin + dqv .lt. 1.e-10 ) then
              temp(i,j,k) = max( 0.0, -(qvin - 1.e-10)/dqv )
            else
              temp(i,j,k) = 1.0
            endif
          end do
        end do
      end do

!-----------------------------------------------------------------------
!    the scaling factor for each column is the minimum value found 
!    within that column.
!-----------------------------------------------------------------------
      scale = minval( temp, dim=3 )

!-----------------------------------------------------------------------
!    scale the cloud tracer, momentum, temperature, precipitation  and 
!    transported tracer tendencies returned from uw shallow convection.
!-----------------------------------------------------------------------
      do k=1,kx
        utnd_uw(:,:,k)  = scale*utnd_uw(:,:,k)
        vtnd_uw(:,:,k)  = scale*vtnd_uw(:,:,k)
        ttnd_uw(:,:,k)  = scale*ttnd_uw(:,:,k)
        qtnd_uw(:,:,k)  = scale*qtnd_uw(:,:,k)
        qltnd_uw(:,:,k) = scale*qltnd_uw(:,:,k)
        qitnd_uw(:,:,k) = scale*qitnd_uw(:,:,k)
        qatnd_uw(:,:,k) = scale*qatnd_uw(:,:,k)
      end do

      if (do_liq_num) then
         do k=1,kx
           qntnd_uw(:,:,k) = scale*qntnd_uw(:,:,k)
         end do
      end if

      if (do_ice_num) then
        do k=1,kx
          qnitnd_uw(:,:,k) = scale*qnitnd_uw(:,:,k)
        end do
      end if

      rain_uw(:,:) = scale*rain_uw(:,:)
      snow_uw(:,:) = scale*snow_uw(:,:)

!-----------------------------------------------------------------------
!    update  total-physics and total-convection tendencies, precipitation
!    fields and the cloud tracer fields and tendencies.
!-----------------------------------------------------------------------
      tdt = tdt + ttnd_uw
      qdt = qdt + qtnd_uw
      udt = udt + utnd_uw
      vdt = vdt + vtnd_uw
      ttnd_conv = ttnd_conv + ttnd_uw
      qtnd_conv = qtnd_conv + qtnd_uw

      lprec = lprec + rain_uw
      fprec = fprec + snow_uw
      precip = precip + rain_uw + snow_uw

      if (doing_prog_clouds) then
        rdt(:,:,:,nql) = rdt(:,:,:,nql) + qltnd_uw
        rdt(:,:,:,nqi) = rdt(:,:,:,nqi) + qitnd_uw
        rdt(:,:,:,nqa) = rdt(:,:,:,nqa) + qatnd_uw
        if (do_liq_num) rdt(:,:,:,nqn) = rdt(:,:,:,nqn) + qntnd_uw
        if (do_ice_num) rdt(:,:,:,nqni) = rdt(:,:,:,nqni) + qnitnd_uw
      endif

!------------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    obtained from uw transport.
!------------------------------------------------------------------------
      nn = 1
      do n=1, num_tracers
        if (tracers_in_uw(n)) then
          rdt(:,:,:,n) = rdt(:,:,:,n) + qtruw (:,:,:,nn)
          nn = nn + 1
        endif
      end do

!-------------------------------------------------------------------------


end subroutine moistproc_scale_uw



!#########################################################################



                     end module moistproc_kernels_mod
