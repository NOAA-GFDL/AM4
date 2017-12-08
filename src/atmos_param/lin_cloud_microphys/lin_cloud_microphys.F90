!
! Cloud micro-physics package for GFDL global cloud resolving model
! The algorithms are originally based on Lin et al 1983. Many key
! elements have been changed/improved based on several other publications
! Developer: Shian-Jiann Lin
!
module lin_cld_microphys_mod
 use mpp_mod,           only: stdlog, mpp_pe, mpp_root_pe, mpp_clock_id, &
                              mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE, &
                              input_nml_file, mpp_max
 use diag_manager_mod,  only: register_diag_field, send_data
 use time_manager_mod,  only: time_type, get_time
 use constants_mod,     only: grav, rdgas, rvgas, cp_air, cp_vapor, hlv, hlf, kappa, pi
 use fms_mod,           only: write_version_number, open_namelist_file, &
                              check_nml_error, file_exist, close_file,  &
                              error_mesg, FATAL

 implicit none
 private

 public  lin_cld_microphys_driver, lin_cld_microphys_init, lin_cld_microphys_end, wqs1, wqs2
 public  qsmith_init, qsmith, es2_table1d, es3_table1d, esw_table1d, wqsat_moist, wqsat2_moist
 public  setup_con, qs_blend
 real             :: missing_value = -1.e10
 logical          :: module_is_initialized = .false.
 logical          :: qsmith_tables_initialized = .false.
 character(len=17) :: mod_name = 'lin_cld_microphys'

!==== fms constants ====================
!!! real, parameter :: latv  = hlv             ! = 2.500e6
!!!  real, parameter:: cv_air = 717.56   ! Satoh value
!!! real, parameter :: lati  = hlf             ! = 3.34e5
!!! real, parameter :: lats  = latv+lati       ! = 2.834E6
! rdgas = 287.04;   rvgas = 461.50
! cp_air =rdgas * 7./2. = 1006.64           ! heat capacity at constant pressure (j/kg/k)
! The following two are from Emanuel's book "Atmospheric Convection"
!!! real, parameter :: c_liq = 4190.           ! heat capacity of water at 0C
                                            !
 real, parameter :: eps   = rdgas/rvgas     ! = 0.621971831
 real, parameter :: zvir  = rvgas/rdgas-1.  ! = 0.607789855
 real, parameter :: table_ice  = 273.16  ! freezing point for qs table
 real, parameter :: cv_air =  cp_air - rdgas ! = rdgas * (7/2-1) = 2.5*rdgas=717.68

   real, parameter:: c_liq = 4190.       ! heat capacity of water at 0C
   real, parameter:: c_ice = 2106.           ! heat capacity of ice at 0C: c=c_ice+7.3*(T-Tice) 
   real, parameter:: cv_vap = 1410.0     ! Emanuel value
   real, parameter:: cp_vap = cp_vapor   !  1850.
   real, parameter:: dc_vap =  cp_vap - c_liq     ! = -2368.
   real, parameter:: dc_ice =  c_liq - c_ice      ! = 2112.

   real, parameter:: hlv0 = 2.501e6   ! Emanuel Appendix-2
   real, parameter:: hlf0 = 3.337e5   ! Emanuel
   real, parameter:: t_ice = 273.15
   real, parameter:: Lv0 =  hlv0 - dc_vap*t_ice   ! = 3.147782e6
   real, parameter:: Li0 =  hlf0 - dc_ice*t_ice   ! = -2.431928e5 
! Lv(T) = Lv0 + dc_vap*T, where Lv0  = hlv + (c_liq-cp_vap)*273.15 = 3.147619e6
! Lv decreasing with increasing T; Lv(273.15) = hlv0; Lv(303) ~ 2.43e6
! Li(T) = Li0 + dc_ice*T
! Li increasing with increasing T
! Lv + Li = Lv0+Li0 + (cp_vap - c_ice)*T = Lv0+Li0 - 256*T
!
!==== fms constants ====================

 real, parameter :: qrmin  = 1.e-9
 real, parameter :: qvmin  = 1.e-22      ! min value for water vapor (treated as zero)
 real, parameter :: qcmin  = 1.e-12      ! min value for cloud condensates
 real, parameter :: sfcrho = 1.20        ! surface air density
 real, parameter :: vmin   = 1.e-2       ! minimum fall speed for rain/graupel
 real, parameter :: rhor   = 1.0e3  ! LFO83
 real, parameter :: dz_min = 1.e-2
 real :: cracs, csacr, cgacr, cgacs, acco(3,4), csacw,          &
         craci, csaci, cgacw, cgaci, cracw, cssub(5), cgsub(5), &
         crevp(5), cgfr(2), csmlt(5), cgmlt(5)
 real :: es0, ces0
 real :: pie, rgrav, fac_rc
 real :: lcp, icp, tcp
 real :: h_var

 logical :: sedi_transport = .true.     !
 logical :: do_sedi_w = .false.
 logical :: do_sedi_heat = .true.     !
 logical :: prog_ccn = .false.     ! do prognostic CCN (Yi Ming's method)
 logical :: do_qa   = .false.      ! do inline cloud fraction
 logical :: rad_snow =.false.
 logical :: rad_graupel =.false.
 logical :: rad_rain =.false.
 logical :: fix_negative =.false.
 logical :: do_setup=.true.
 logical :: master
 logical :: p_nonhydro = .false.

 real, allocatable:: table(:), table2(:), table3(:), tablew(:), des(:), des2(:), des3(:), desw(:)
 logical :: tables_are_initialized = .false.

 integer:: id_rh, id_vtr, id_vts,  id_vtg, id_vti, id_rain, id_snow, id_graupel, &
           id_ice, id_prec, id_cond, id_var, id_droplets, id_dbz, id_maxdbz
 real:: lati, latv, lats

 real, parameter :: dt_fr = 8.       ! homogeneous freezing of all cloud water at t_wfr - dt_fr
                                     ! minimum temperature water can exist (Moore & Molinero Nov. 2011, Nature)
                                     ! dt_fr can be considered as the error bar
 integer :: lin_cld_mp_clock   ! clock for timing of driver routine

 real :: t_snow_melt = 12.      ! snow melt tempearture scale factor
 real :: t_grau_melt = 15.      ! graupel melt tempearture scale factor

! The defaults are good for 25-50 km simulation
! For cloud-resolving: 1-5 km
!                       qi0_crt = 0.8E-4
!                       qs0_crt = 0.6E-3
!                       c_psaci = 0.1
!                       c_pgacs = 0.01
!----------------------
! namelist  parameters:
!----------------------
 real :: cld_min = 0.05
 real :: tice  = 273.16  ! set tice = 165. to trun off ice-phase phys (Kessler emulator)

 real :: qc_crt  = 5.0e-8  ! minimum condensate mixing ratio to allow partial cloudiness
 real :: t_min   = 161.  ! Min temperature for ice-phase micro phys
 real :: mp_time = 120.  ! maximum micro-physics time step (sec)

 real :: rh_inc = 0.10   ! rh increment for complete evap of ql and qi
 real :: rh_inr = 0.25
 real :: rh_ins = 0.25   ! rh increment for sublimation of snow

! The following 3 time scales are for melting during terminal falls
 real :: tau_s  = 120.    ! snow melt
 real :: tau_g  = 180.    ! graupel melt
 real :: tau_mlt = 10.    ! ice melting time-scale

! cloud water
 real :: tau_v2l = 600.   ! vapor --> cloud water (condensation)  time scale
 real :: tau_l2v = 600.   ! cloud water --> vapor (evaporation)  time scale
! Graupel
 real :: tau_g2v = 900.   ! Grapuel sublimation time scale
 real :: tau_v2g =21600. ! Grapuel deposition -- make it a slow process

 real :: dw_land  = 0.20  ! base value for subgrid deviation/variability over land
 real :: dw_ocean = 0.15  ! base value for ocean
 real :: ccn_o = 100.
 real :: ccn_l = 250.
 real :: rthresh = 10.0e-6     ! critical cloud drop radius (micro m)

!-------------------------------------------------------------
! WRF/WSM6 scheme: qi_gen = 4.92e-11 * (1.e3*exp(0.1*tmp))**1.33
!  optimized:      qi_gen = 4.92e-11 * exp( 1.33*log(1.e3*exp(0.1*tmp)) )
! qi_gen ~ 4.808e-7 at 0 C; 1.818e-6 at -10 C, 9.82679e-5 at -40C
! the following value is constructed such that qc_crt = 0 at zero C and @ -10C matches
! WRF/WSM6 ice initiation scheme; qi_crt = qi_gen*min(qi_lim, 0.1*tmp) / den
 real :: qi_gen  = 1.818E-6
 real :: qi_lim  = 1.
 real :: ql_mlt  = 4.0e-3    ! max value of cloud water allowed from melted cloud ice
 real :: ql_gen  = 6.0e-4    ! max ql generation during remapping step if fast_sat_adj = .T.
 real :: c0_max  = 1.0e-4    ! maximum total condensates from subgrid saturation
 real :: sat_adj0 = 0.99     ! adjustment factor (0: no, 1: full) during fast_sat_adj

! Cloud condensate upper bounds: "safety valves" for ql & qi
 real :: ql0_max = 4.0e-3    ! max ql value (converted to rain)
 real :: qi0_max = 5.0e-3    ! max qi value (converted to snow)

 real :: qi0_crt = 1.0e-4    ! ice  --> snow autocon threshold (was 1.E-4)
                             ! qi0_crt is highly dependent on horizontal resolution
 real :: qr0_crt = 1.0e-4    ! rain --> snow or graupel/hail threshold
                             ! LFO used *mixing ratio* = 1.E-4 (hail in LFO)
 real :: c_psaut = 1.0e-3   ! autoconversion rate: cloud_ice -> snow
 real :: c_psaci = 0.01     ! accretion: cloud ice --> snow (was 0.1 in Zetac)
 real :: c_piacr = 5.0      ! accretion: rain --> ice:
 real :: c_cracw = 0.9      ! rain accretion efficiency

! Decreasing  clin to reduce csacw (so as to reduce cloud water ---> snow)
 real:: alin = 842.0
 real:: clin = 4.8      ! 4.8 --> 6. (to ehance ql--> qs)

!-----------------
! Graupel control:
!-----------------
 real :: qs0_crt = 2.0e-3   ! snow --> graupel density threshold (0.6e-3 in Purdue Lin scheme)
 real :: c_pgacs = 2.0e-3   ! snow --> graupel "accretion" eff. (was 0.1 in Zetac)

! fall velocity tuning constants:
 real :: den_ref = sfcrho   ! Reference (surface) density for fall speed
                            ! Larger value produce larger fall speed
 real :: vr_fac = 1.
 real :: vs_fac = 1.
 real :: vg_fac = 1.
 real :: vi_fac = 1.

 logical :: fast_sat_adj  = .false.
 logical :: z_slope_liq  = .true.          !  use linear mono slope for autocconversions
 logical :: z_slope_ice  = .false.          !  use linear mono slope for autocconversions
 logical :: use_deng_mace = .true.       ! Helmfield-Donner ice speed
 logical :: do_subgrid_z = .false.       ! 2X resolution sub-grid saturation/cloud scheme
 logical :: use_ccn      = .false.
 logical :: use_ppm      = .true.
 logical :: ppm_rain_fall  = .true.
 logical :: mono_prof = .true.          ! perform terminal fall with mono ppm scheme
 logical :: mp_debug = .false.
 logical :: mp_print = .false.

 real:: global_area = -1.

 real:: tice0, t_wfr
 real:: p_crt   = 100.E2   !
 integer:: k_moist = 100

 namelist /lin_cld_microphys_nml/mp_time, t_min, tau_s, tau_g, dw_land, dw_ocean,  &
                      vr_fac, vs_fac, vg_fac, vi_fac, ql_mlt, do_qa, fix_negative, &
                      qs0_crt, qi_gen, ql0_max, qi0_max, qi0_crt, qr0_crt, fast_sat_adj, &
                      rh_inc, rh_ins, rh_inr, use_deng_mace, use_ccn, do_subgrid_z,  &
                      rthresh, ccn_l, ccn_o, qc_crt, tau_g2v, tau_v2g, sat_adj0,    &
                      c_piacr, tau_mlt, tau_v2l, tau_l2v, qi_lim, ql_gen, c0_max,     &
                      c_psaut, c_psaci, c_pgacs, z_slope_liq, z_slope_ice, prog_ccn,  &
                      c_cracw, alin, clin, p_crt, tice, k_moist, rad_snow, rad_graupel, rad_rain,   &
                      cld_min, use_ppm, ppm_rain_fall, mono_prof, do_sedi_heat,sedi_transport,   &
                      do_sedi_w,  mp_debug, mp_print

!---- version number -----
 character(len=128) :: version = '$Id$'
 character(len=128) :: tagname = '$Name$'

 logical  :: hydrostatic, phys_hydrostatic

 contains


  subroutine lin_cld_microphys_driver(qv, ql, qr, qi, qs, qg, qa, qn,                &
                               qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt,      &
                               pt_dt, pt, w, uin, vin, udt, vdt, dz, delp, area, dt_in, &
                               land,  rain, snow, ice, graupel,                      &
                               iis,iie, jjs,jje, kks,kke, ktop, kbot, time)
! kks == 1; kke == kbot == npz
  type(time_type), intent(in):: time
  integer,         intent(in):: iis,iie, jjs,jje  ! physics window
  integer,         intent(in):: kks,kke           ! vertical dimension
  integer,         intent(in):: ktop, kbot        ! vertical compute domain
  real,            intent(in):: dt_in

  real, intent(in   ), dimension(:,:)  :: area
  real, intent(in   ), dimension(:,:)  :: land  !land fraction
  real, intent(out  ), dimension(:,:)  :: rain, snow, ice, graupel
  real, intent(in   ), dimension(:,:,:):: delp, dz, uin, vin
  real, intent(in   ), dimension(:,:,:):: pt, qv, ql, qr, qi, qs, qg, qa, qn 
  real, intent(inout), dimension(:,:,:):: pt_dt,  qa_dt, udt, vdt, w
  real, intent(inout), dimension(:,:,:):: qv_dt, ql_dt, qr_dt, qi_dt,  &
                                          qs_dt, qg_dt

   call error_mesg('lin_cld_microphys_driver','Lin cloud microphysics is not available.',FATAL)

 end subroutine lin_cld_microphys_driver



 subroutine check_column(ktop, kbot, q, no_fall)
 integer, intent(in):: ktop, kbot
 real,    intent(in):: q(ktop:kbot)
 logical, intent(out):: no_fall
! local:
 integer k

 no_fall = .true.
 do k=ktop, kbot
    if ( q(k) > qrmin ) then
         no_fall = .false.
         exit
    endif
 enddo

 end subroutine check_column


 subroutine lin_cld_microphys_init(id, jd, kd, axes, time, hydrostatic_in,&
                  phys_hydrostatic_in)

    integer,         intent(in) :: id, jd, kd
    integer,         intent(in) :: axes(4)
    type(time_type), intent(in) :: time
    logical,         intent(in) :: hydrostatic_in, phys_hydrostatic_in

   call error_mesg('lin_cld_microphys_init','Lin cloud microphysics is not available.',FATAL)

 end subroutine lin_cld_microphys_init



 subroutine lin_cld_microphys_end

   deallocate ( table  )
   deallocate ( table2 )
   deallocate ( table3 )
   deallocate ( tablew )
   deallocate ( des )
   deallocate ( des2 )
   deallocate ( des3 )
   deallocate ( desw )

   tables_are_initialized = .false.

   call error_mesg('lin_cld_microphys_end','Lin cloud microphysics is not available.',FATAL)

 end subroutine lin_cld_microphys_end



 subroutine setup_con

  master = (mpp_pe().eq.mpp_root_pe())
  rgrav = 1./ grav

  if ( .not. qsmith_tables_initialized ) call qsmith_init
  qsmith_tables_initialized = .true.

 end subroutine setup_con



 subroutine qsmith_init
  integer, parameter:: length=2621
  integer i

  if( .not. tables_are_initialized ) then

    master = (mpp_pe().eq.mpp_root_pe())
    if (master) print*, ' lin MP: initializing qs tables'
!!! DEBUG CODE
!    print*, mpp_pe(), allocated(table), allocated(table2), allocated(table3), allocated(tablew), allocated(des), allocated(des2), allocated(des3), allocated(desw)
!!! END DEBUG CODE

!                            generate es table (dt = 0.1 deg. c)
       allocate ( table( length) )
       allocate ( table2(length) )
       allocate ( table3(length) )
       allocate ( tablew(length) )
       allocate (   des (length) )
       allocate (   des2(length) )
       allocate (   des3(length) )
       allocate (   desw(length) )

       call qs_table (length )
       call qs_table2(length )
       call qs_table3(length )
       call qs_tablew(length )

       do i=1,length-1
           des(i) = max(0.,  table(i+1) -  table(i))
          des2(i) = max(0., table2(i+1) - table2(i))
          des3(i) = max(0., table3(i+1) - table3(i))
          desw(i) = max(0., tablew(i+1) - tablew(i))
       enddo
        des(length) =  des(length-1)
       des2(length) = des2(length-1)
       des3(length) = des3(length-1)
       desw(length) = desw(length-1)

       tables_are_initialized = .true.
  endif

 end subroutine qsmith_init

 real function wqs1(ta, den)
! Pure water phase; universal dry/moist formular using air density
! Input "den" can be either dry or moist air density
  real, intent(in):: ta, den
! local:
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = tablew(it) + (ap1-it)*desw(it)
      wqs1 = es / (rvgas*ta*den)

 end function wqs1

 real function wqs2(ta, den, dqdt)
! Pure water phase; universal dry/moist formular using air density
! Input "den" can be either dry or moist air density
  real, intent(in):: ta, den
  real, intent(out):: dqdt
! local:
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  integer it

  if (.not. tables_are_initialized) call qsmith_init

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = tablew(it) + (ap1-it)*desw(it)
      wqs2 = es / (rvgas*ta*den)
        it = ap1 - 0.5
! Finite diff, del_T = 0.1:
      dqdt = 10.*(desw(it) + (ap1-it)*(desw(it+1)-desw(it))) / (rvgas*ta*den)

 end function wqs2

 real function iqs1(ta, den)
! water-ice phase; universal dry/moist formular using air density
! Input "den" can be either dry or moist air density
  real, intent(in):: ta, den
! local:
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table2(it) + (ap1-it)*des2(it)
      iqs1 = es / (rvgas*ta*den)

 end function iqs1

 real function iqs2(ta, den, dqdt)
! water-ice phase; universal dry/moist formular using air density
! Input "den" can be either dry or moist air density
  real, intent(in):: ta, den
  real, intent(out):: dqdt
! local:
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table2(it) + (ap1-it)*des2(it)
      iqs2 = es / (rvgas*ta*den)
        it = ap1 - 0.5
      dqdt = 10.*(des2(it) + (ap1-it)*(des2(it+1)-des2(it))) / (rvgas*ta*den)

 end function iqs2

 real function qs1d_moist(ta, qv, pa, dqdt)
! 2-phase tabel
  real, intent(in):: ta, pa, qv
  real, intent(out):: dqdt
! local:
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  real, parameter:: eps10 = 10.*eps
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table2(it) + (ap1-it)*des2(it)
      qs1d_moist = eps*es*(1.+zvir*qv)/pa
        it = ap1 - 0.5
      dqdt = eps10*(des2(it) + (ap1-it)*(des2(it+1)-des2(it)))*(1.+zvir*qv)/pa

 end function qs1d_moist

 real function wqsat2_moist(ta, qv, pa, dqdt)
! Pure water phase
  real, intent(in):: ta, pa, qv
  real, intent(out):: dqdt
! local:
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  real, parameter:: eps10 = 10.*eps
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = tablew(it) + (ap1-it)*desw(it)
     wqsat2_moist = eps*es*(1.+zvir*qv)/pa
     dqdt = eps10*(desw(it) + (ap1-it)*(desw(it+1)-desw(it)))*(1.+zvir*qv)/pa

 end function wqsat2_moist

 real function wqsat_moist(ta, qv, pa)
! Pure water phase
  real, intent(in):: ta, pa, qv
! local:
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = tablew(it) + (ap1-it)*desw(it)
     wqsat_moist = eps*es*(1.+zvir*qv)/pa

 end function wqsat_moist

 real function qs1d_m(ta, qv, pa)
! 2-phase tabel
  real, intent(in):: ta, pa, qv
! local:
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  real, parameter:: eps10 = 10.*eps
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table2(it) + (ap1-it)*des2(it)
      qs1d_m = eps*es*(1.+zvir*qv)/pa

 end function qs1d_m

 real function d_sat(ta)
! Computes the difference in saturation vapor *density* between water and ice
  real, intent(in):: ta
  real, parameter:: tmin=table_ice - 160.
  real es_w, es_i, ap1
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
! over Water:
       es_w = tablew(it) + (ap1-it)*desw(it)
! over Ice:
       es_i = table2(it) + (ap1-it)*des2(it)
      d_sat = dim(es_w, es_i)/(rvgas*ta)  ! Take positive difference

 end function d_sat


 real function esw_table(ta)
! pure water phase table
  real, intent(in):: ta
  real, parameter:: tmin=table_ice - 160.
  real  ap1
  integer it
       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
      esw_table = tablew(it) + (ap1-it)*desw(it)
 end function esw_table


 real function es2_table(ta)
! two-phase table
  real, intent(in):: ta
  real, parameter:: tmin=table_ice - 160.
  real  ap1
  integer it
       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
      es2_table = table2(it) + (ap1-it)*des2(it)
 end function es2_table


 subroutine esw_table1d(ta, es, n)
  integer, intent(in):: n
! For waterphase only
  real, intent(in)::  ta(n)
  real, intent(out):: es(n)
  real, parameter:: tmin=table_ice - 160.
  real  ap1
  integer i, it

  do i=1, n
       ap1 = 10.*dim(ta(i), tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
     es(i) = tablew(it) + (ap1-it)*desw(it)
  enddo
 end subroutine esw_table1d



 subroutine es2_table1d(ta, es, n)
  integer, intent(in):: n
! two-phase table with -2C as the transition point for ice-water phase
! For sea ice model
  real, intent(in)::  ta(n)
  real, intent(out):: es(n)
  real, parameter:: tmin=table_ice - 160.
  real  ap1
  integer i, it

  do i=1, n
       ap1 = 10.*dim(ta(i), tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
     es(i) = table2(it) + (ap1-it)*des2(it)
  enddo
 end subroutine es2_table1d


 subroutine es3_table1d(ta, es, n)
  integer, intent(in):: n
! two-phase table with -2C as the transition point for ice-water phase
  real, intent(in)::  ta(n)
  real, intent(out):: es(n)
  real, parameter:: tmin=table_ice - 160.
  real  ap1
  integer i, it

  do i=1, n
       ap1 = 10.*dim(ta(i), tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
     es(i) = table3(it) + (ap1-it)*des3(it)
  enddo
 end subroutine es3_table1d



 subroutine qs_tablew(n)
! 2-phase table
      integer, intent(in):: n
      real:: delt=0.1
      real esbasw, tbasw, esbasi, tbasi, tmin, tem, aa, b, c, d, e
      integer i

! constants
      esbasw = 1013246.0
       tbasw =     373.16
      esbasi =    6107.1
       tbasi =     273.16
        tmin = tbasi - 160.

     do i=1,n
        tem = tmin+delt*real(i-1)
!  compute es over water
!  see smithsonian meteorological tables page 350.
        aa  = -7.90298*(tbasw/tem-1.)
        b   =  5.02808*alog10(tbasw/tem)
        c   = -1.3816e-07*(10**((1.-tem/tbasw)*11.344)-1.)
        d   =  8.1328e-03*(10**((tbasw/tem-1.)*(-3.49149))-1.)
        e   = alog10(esbasw)
        tablew(i) = 0.1 * 10**(aa+b+c+d+e)
     enddo

 end subroutine qs_tablew


 subroutine qs_table2(n)
! 2-phase table
  integer, intent(in):: n
  real:: delt=0.1
  real esbasw, tbasw, esbasi, tbasi, tmin, tem, aa, b, c, d, e
  integer :: i0, i1
  real :: tem0, tem1
  integer i

! constants
      esbasw = 1013246.0
       tbasw =     373.16
      esbasi =    6107.1
       tbasi =     273.16
      tmin = tbasi - 160.

     do i=1,n
        tem = tmin+delt*real(i-1)
        if ( i<= 1600 ) then
!  compute es over ice between -160c and 0 c.
!  see smithsonian meteorological tables page 350.
              aa  = -9.09718 *(tbasi/tem-1.)
              b   = -3.56654 *alog10(tbasi/tem)
              c   =  0.876793*(1.-tem/tbasi)
              e   = alog10(esbasi)
             table2(i) = 0.1 * 10**(aa+b+c+e)
        else
!  compute es over water between 0c and 102c.
!  see smithsonian meteorological tables page 350.
             aa  = -7.90298*(tbasw/tem-1.)
             b   =  5.02808*alog10(tbasw/tem)
             c   = -1.3816e-07*(10**((1.-tem/tbasw)*11.344)-1.)
             d   =  8.1328e-03*(10**((tbasw/tem-1.)*(-3.49149))-1.)
             e   = alog10(esbasw)
             table2(i) = 0.1 * 10**(aa+b+c+d+e)
        endif
     enddo

!----------
! smoother
!----------
      i0 = 1600;  i1 = 1601
      tem0 = 0.25*(table2(i0-1) + 2.*table(i0) + table2(i0+1))
      tem1 = 0.25*(table2(i1-1) + 2.*table(i1) + table2(i1+1))
      table2(i0) = tem0
      table2(i1) = tem1

 end subroutine qs_table2



 subroutine qs_table3(n)
! 2-phase table with "-2 C" as the transition point
  integer, intent(in):: n
  real:: delt=0.1
  real esbasw, tbasw, esbasi, tbasi, tmin, tem, aa, b, c, d, e
  integer :: i0, i1
  real :: tem0, tem1
  integer i

! constants
      esbasw = 1013246.0
       tbasw =     373.16
      esbasi =    6107.1
       tbasi =     273.16
      tmin = tbasi - 160.

     do i=1,n
        tem = tmin+delt*real(i-1)
!       if ( i<= 1600 ) then
        if ( i<= 1580 ) then  ! to -2 C
!  compute es over ice between -160c and 0 c.
!  see smithsonian meteorological tables page 350.
              aa  = -9.09718 *(tbasi/tem-1.)
              b   = -3.56654 *alog10(tbasi/tem)
              c   =  0.876793*(1.-tem/tbasi)
              e   = alog10(esbasi)
             table3(i) = 0.1 * 10**(aa+b+c+e)
        else
!  compute es over water between -2c and 102c.
!  see smithsonian meteorological tables page 350.
             aa  = -7.90298*(tbasw/tem-1.)
             b   =  5.02808*alog10(tbasw/tem)
             c   = -1.3816e-07*(10**((1.-tem/tbasw)*11.344)-1.)
             d   =  8.1328e-03*(10**((tbasw/tem-1.)*(-3.49149))-1.)
             e   = alog10(esbasw)
             table3(i) = 0.1 * 10**(aa+b+c+d+e)
        endif
     enddo

!----------
! smoother
!----------
      i0 = 1580
      tem0 = 0.25*(table3(i0-1) + 2.*table(i0) + table3(i0+1))
      i1 = 1581
      tem1 = 0.25*(table3(i1-1) + 2.*table(i1) + table3(i1+1))
      table3(i0) = tem0
      table3(i1) = tem1

 end subroutine qs_table3


 real function qs_blend(t, p, q)
! Note: this routine is based on "moist" mixing ratio
! Blended mixed phase table
  real, intent(in):: t, p, q
  real es, ap1
  real, parameter:: tmin=table_ice - 160.
  integer it

       ap1 = 10.*dim(t, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table(it) + (ap1-it)*des(it)
      qs_blend = eps*es*(1.+zvir*q)/p

 end function qs_blend

 subroutine qs_table(n)
      integer, intent(in):: n
      real esupc(200)
      real:: delt=0.1
      real esbasw, tbasw, esbasi, tbasi, tmin, tem, aa, b, c, d, e, esh20
      real wice, wh2o
      integer i

! constants
      esbasw = 1013246.0
       tbasw =     373.16
      esbasi =    6107.1
       tbasi =     273.16

!  compute es over ice between -160c and 0 c.
      tmin = tbasi - 160.
!  see smithsonian meteorological tables page 350.
      do i=1,1600
         tem = tmin+delt*real(i-1)
         aa  = -9.09718 *(tbasi/tem-1.)
         b   = -3.56654 *alog10(tbasi/tem)
         c   =  0.876793*(1.-tem/tbasi)
         e   = alog10(esbasi)
         table(i)=10**(aa+b+c+e)
      enddo

!  compute es over water between -20c and 102c.
!  see smithsonian meteorological tables page 350.
      do  i=1,1221
          tem = 253.16+delt*real(i-1)
          aa  = -7.90298*(tbasw/tem-1.)
          b   =  5.02808*alog10(tbasw/tem)
          c   = -1.3816e-07*(10**((1.-tem/tbasw)*11.344)-1.)
          d   =  8.1328e-03*(10**((tbasw/tem-1.)*(-3.49149))-1.)
          e   = alog10(esbasw)
          esh20  = 10**(aa+b+c+d+e)
          if (i <= 200) then
              esupc(i) = esh20
          else
              table(i+1400) = esh20
          endif
      enddo

!  derive blended es over ice and supercooled water between -20c and 0c
      do i=1,200
         tem  = 253.16+delt*real(i-1)
         wice = 0.05*(273.16-tem)
         wh2o = 0.05*(tem-253.16)
         table(i+1400) = wice*table(i+1400)+wh2o*esupc(i)
      enddo

      do i=1,n
         table(i) = table(i)*0.1
      enddo

 end subroutine qs_table


 subroutine qsmith(im, km, ks, t, p, q, qs, dqdt)
! input t in deg k; p (pa) : moist pressure
  integer, intent(in):: im, km, ks
  real, intent(in),dimension(im,km):: t, p, q
  real, intent(out),dimension(im,km):: qs
  real, intent(out), optional:: dqdt(im,km)
! local:
  real, parameter:: eps10 = 10.*eps
  real es(im,km)
  real ap1
  real, parameter:: tmin=table_ice - 160.
  integer i, k, it

  if( .not. tables_are_initialized ) then
       call  qsmith_init
  endif

      do k=ks,km
         do i=1,im
            ap1 = 10.*dim(t(i,k), tmin) + 1.
            ap1 = min(2621., ap1)
            it = ap1
            es(i,k) = table(it) + (ap1-it)*des(it)
            qs(i,k) = eps*es(i,k)*(1.+zvir*q(i,k))/p(i,k)
         enddo
      enddo

      if ( present(dqdt) ) then
      do k=ks,km
           do i=1,im
              ap1 = 10.*dim(t(i,k), tmin) + 1.
              ap1 = min(2621., ap1) - 0.5
              it  = ap1
              dqdt(i,k) = eps10*(des(it)+(ap1-it)*(des(it+1)-des(it)))*(1.+zvir*q(i,k))/p(i,k)
           enddo
      enddo
      endif

 end subroutine qsmith

end module lin_cld_microphys_mod
