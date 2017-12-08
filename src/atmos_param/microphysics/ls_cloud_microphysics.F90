                    module ls_cloud_microphysics_mod

!-----------------------------------------------------------------------
!
!         interface module for cloud microphysics
!         ---------------------------------------
!         OPTIONS AVAILABLE:
!             Rotstayn-Klein microphysics (this used in old strat_cloud)
!             Lin microphysics         
!             MG microphysics (as developed by M. Salzmann)
!             MG-NCAR microphysics (an early version of NCAR microphysics)
!             NCAR microphysics version 1.5 (became available in 2012)
!
!-----------------------------------------------------------------------

! fms modules
use time_manager_mod,      only: time_type, get_time, set_date
use mpp_mod,               only: input_nml_file
use fms_mod,               only: error_mesg, FATAL, NOTE,        &
                                 file_exist, check_nml_error,    &
                                 open_namelist_file, close_file, &
                                 write_version_number, stdlog,   &
                                 mpp_pe, mpp_root_pe, stdlog,    &
                                 mpp_clock_id, mpp_clock_begin,  &
                                 mpp_clock_end, CLOCK_MODULE,    &
                                 CLOCK_MODULE_DRIVER, &
                                 MPP_CLOCK_SYNC, read_data, write_data
use field_manager_mod,     only: MODEL_ATMOS
use tracer_manager_mod,    only: get_tracer_index,&
                                 get_number_tracers, &
                                 get_tracer_names, &
                                 query_method, &
                                 NO_TRACER
use constants_mod,         only: CP_AIR, GRAV, HLV, HLS, HLF, &
                                 RDGAS, RVGAS, TFREEZE, WTMAIR, &
                                 SECONDS_PER_DAY, KAPPA

!  shared physics and physics utilities modules

use physics_types_mod,     only: physics_control_type
use lscloud_types_mod,     only: lscloud_types_init, atmos_state_type, &
                                 diag_id_type, diag_pt_type, &
                                 lsc_constants_type, lscloud_nml_type, &
                                 lscloud_debug_type, particles_type,  &
                                 cloud_state_type, cloud_processes_type,  &
                                 precip_state_type
use aerosol_types_mod,     only: aerosol_type
use physics_radiation_exch_mod,      &
                           only: exchange_control_type
use moist_proc_utils_mod,  only: mp_input_type, mp_output_type,  &
                                 mp_lsdiag_type, mp_nml_type,  &
                                 mp_lsdiag_control_type,  &
                                 mp_conv2ls_type, mp_tendency_type,  &
                                 mp_removal_type

! physics modules

use lin_cld_microphys_mod, only: lin_cld_microphys_init, &
                                 setup_con,   &
                                 lin_cld_microphys_driver, &
                                 lin_cld_microphys_end
use lscloud_debug_mod,     only: write_debug_output
use rotstayn_klein_mp_mod, only: rotstayn_klein_microp, &
                                 rotstayn_klein_microp_init,  &
                                 rotstayn_klein_microp_end
use morrison_gettelman_microp_mod,     &
                           only: morrison_gettelman_microp, &
                                 morrison_gettelman_microp_init, &
                                 morrison_gettelman_microp_end
use cldwat2m_micro_mod,    only: ini_micro, mmicro_pcond, mmicro_end
use micro_mg_mod,          only: micro_mg_init, micro_mg_get_cols,&
                                 micro_mg_tend

implicit none
private

!-----------------------------------------------------------------------
!-------------------- public data/interfaces ---------------------------

public   ls_cloud_microphysics_init, ls_cloud_microphysics, &
         ls_cloud_microphysics_end, ls_cloud_microphysics_time_vary
  
!-----------------------------------------------------------------------
!-------------------- private data -------------------------------------

private  adjust_precip_fields, adjust_for_supersaturation_removal,  &
         destroy_tiny_clouds, destroy_tiny_clouds_clubb


!--------------------- version number ----------------------------------
character(len=128) :: version = '$Id: $'
character(len=128) :: tagname = '$Name: $'

!--------------------------------------------------------------------------
!---namelist---------------------------------------------------------------
 
real    :: lin_microphys_top_press = 10.E2
                                      ! top pressure level at which Lin
                                      ! microphysics will be calculated
logical :: mass_cons = .true.         ! should we ensure water mass 
                                      ! conservation by adjusting precip 
                                      ! to balance column water 
                                      ! mass change ?
integer :: override_liq_num = 0       ! override model predicted droplet
                                      ! number ? 0 = no, otherwise = y
integer :: override_ice_num = 0       ! override model predicted ice 
                                      ! number ? 0 = no, otherwise = y
logical :: use_Meyers = .false.       ! use Meyers formula when overriding
                                      ! model ice particle number ?
logical :: use_Cooper = .false.       ! use Cooper formula when overriding
                                      ! model ice particle number ?
integer, dimension(6) :: init_date = (/ 1980, 1, 1, 0, 0, 0 /)  
                                      ! date to use as base for  
                                      ! defining microphysics start time
real    :: micro_begin_sec  = 0.0     ! begin microphysics this many 
                                      ! seconds after init_date
integer :: top_lev = 1                ! topmost level for ncar microphysics
real    :: min_precip_needing_adjustment     = 0.0      
real    :: lowest_allowed_precip = 0.0
logical :: use_ndust = .false.


namelist / ls_cloud_microphysics_nml /   &
                               lin_microphys_top_press, mass_cons, &
                               override_liq_num, override_ice_num, &
                               use_Meyers, use_Cooper, init_date, &
                               micro_begin_sec, top_lev, &
                               min_precip_needing_adjustment, &
                               lowest_allowed_precip, use_ndust


!-------------------- clock definitions --------------------------------

integer  :: rk_micro_clock, lin_micro_clock, ncar_micro_clock

!----------------------------------------------------------------------
!    module variables retrieved from other modules
!----------------------------------------------------------------------
real    :: qmin
integer :: do_clubb
logical :: limit_conv_cloud_frac
logical :: do_lin_cld_microphys
integer :: super_ice_opt
logical :: do_pdf_clouds
logical :: doing_prog_clouds
real    :: dtcloud, inv_dtcloud
logical :: do_rk_microphys, do_mg_microphys, do_mg_ncar_microphys, &
           do_ncar_microphys
logical :: tiedtke_macrophysics
logical :: dqa_activation, total_activation
integer :: nsphum, nql, nqi, nqa, nqn, nqni, nqr, nqs, nqg

!--------------------------------------------------------------------
!    other module variables
!--------------------------------------------------------------------

integer, parameter          :: r8 = selected_real_kind(12)   
                                  ! 8 byte real
integer                     :: current_days0, current_sec0   
                                  ! variables related to delayed initiation
                                  ! of microphysics
integer                     :: ktop   
                                  ! top layer index for Lin Micro-Physics



logical            :: module_is_initialized = .false.




                             contains



!#######################################################################

subroutine ls_cloud_microphysics_init  (    &
                     Nml_mp, Constants_lsc, Physics_control, &
                     id, jd, kd, Time, axes, pref, Nml_lsc, Exch_ctrl)

!------------------------------------------------------------------------

type(mp_nml_type),           intent(in)    :: Nml_mp
type(lsc_constants_type),    intent(in)    :: Constants_lsc
type(physics_control_type),  intent(in)    :: Physics_control
integer,                     intent(in)    :: id, jd, kd
integer,                     intent(in)    :: axes(4)
type(time_type),             intent(in)    :: Time
real, dimension(:),          intent(in)    :: pref
type(lscloud_nml_type),      intent(in)    :: Nml_lsc
type(exchange_control_type), intent(in)    :: Exch_ctrl

!------------------------------------------------------------------------
! local variables:    
      integer              :: logunit, io, ierr
      type(time_type)      :: Time_init
      character(len=128)   :: errstring ! Output status: non-blank for 
                                        ! error return
      integer              :: k
      integer              :: rk_micro_init_clock, lin_micro_init_clock, &
                              ncar_micro_init_clock

!-----------------------------------------------------------------------

      if (module_is_initialized) return

!-----------------------------------------------------------------------
!    save variables needed from other modules as module variables
!-----------------------------------------------------------------------
      qmin = Exch_ctrl%qmin
      do_clubb = Exch_ctrl%do_clubb
      limit_conv_cloud_frac = Nml_mp%limit_conv_cloud_frac
      do_lin_cld_microphys = Constants_lsc%do_lin_cld_microphys
      super_ice_opt = Nml_lsc%super_ice_opt
      do_pdf_clouds = Nml_lsc%do_pdf_clouds
      doing_prog_clouds = Exch_ctrl%doing_prog_clouds
      do_rk_microphys = Constants_lsc%do_rk_microphys
      do_mg_microphys = Constants_lsc%do_mg_microphys
      do_mg_ncar_microphys = Constants_lsc%do_mg_ncar_microphys
      do_ncar_microphys = Constants_lsc%do_ncar_microphys
      tiedtke_macrophysics = Constants_lsc%tiedtke_macrophysics
      dqa_activation = Constants_lsc%dqa_activation
      total_activation = Constants_lsc%total_activation
      nql = Physics_control%nql
      nqi = Physics_control%nqi
      nqa = Physics_control%nqa
      nqn = Physics_control%nqn
      nqni = Physics_control%nqni
      nqg = Physics_control%nqg
      nqr = Physics_control%nqr
      nqs = Physics_control%nqs

!------------------------------------------------------------------------
!    define clocks for large-scale cloud schemes initialization (local
!    variables) and for the prognostic loop clocks (module variables).
!------------------------------------------------------------------------
      if (do_rk_microphys) then
        rk_micro_init_clock = mpp_clock_id(   &
               '   Ls_cld_micro: rk_micro:Initialization' , &
                                                grain=CLOCK_MODULE_DRIVER )
      else if (do_lin_cld_microphys) then
        lin_micro_init_clock = mpp_clock_id(     &
               '   Ls_cld_micro: lin_micro:Initialization' , &
                                                grain=CLOCK_MODULE_DRIVER )
      else if (do_mg_microphys) then
        ncar_micro_init_clock = mpp_clock_id(    &
               '   Ls_cld_micro: mg_micro:Initialization' , &
                                                grain=CLOCK_MODULE_DRIVER)
      else if (do_mg_ncar_microphys) then
        ncar_micro_init_clock = mpp_clock_id(    &
               '   Ls_cld_micro: mg_ncar_micro:Initialization' , &
                                                grain=CLOCK_MODULE_DRIVER)
      else if (do_ncar_microphys) then
        ncar_micro_init_clock = mpp_clock_id(    &
               '   Ls_cld_micro: ncar_micro:Initialization' , &
                                                grain=CLOCK_MODULE_DRIVER)
      endif

      if (do_rk_microphys) then
        rk_micro_clock = mpp_clock_id(   &
               '   Ls_cld_micro: rk_micro' , &
                                                grain=CLOCK_MODULE_DRIVER )
      else if (do_lin_cld_microphys) then
        lin_micro_clock = mpp_clock_id(     &
               '   Ls_cld_micro: lin_micro' , &
                                                grain=CLOCK_MODULE_DRIVER )
      else if (do_mg_microphys) then
        ncar_micro_clock = mpp_clock_id(    &
               '   Ls_cld_micro: mg_micro' , &
                                                grain=CLOCK_MODULE_DRIVER )
      else if (do_mg_ncar_microphys) then
        ncar_micro_clock = mpp_clock_id(    &
               '   Ls_cld_micro: mg_ncar_micro' , &
                                                grain=CLOCK_MODULE_DRIVER)
      else if (do_ncar_microphys) then
        ncar_micro_clock = mpp_clock_id(    &
               '   Ls_cld_micro: ncar_micro' , &
                                                grain=CLOCK_MODULE_DRIVER )
      endif

!-------------------------------------------------------------------------
!    process namelist.
!-------------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=ls_cloud_microphysics_nml, iostat=io)
      ierr = check_nml_error(io,'ls_cloud_microphysics_nml')
#else
      if ( file_exist('input.nml')) then
        unit = open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=ls_cloud_microphysics_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'ls_cloud_microphysics_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!-------------------------------------------------------------------------
!    write version number and namelist to standard log.
!-------------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe()) &
                          write (logunit, nml=ls_cloud_microphysics_nml)

!-------------------------------------------------------------------------
!    make sure needed modules have been initialized.
!-------------------------------------------------------------------------
      call lscloud_types_init

!-------------------------------------------------------------------------
!     perform consistency / realizability checks.
!-------------------------------------------------------------------------
      if (do_clubb > 0 .and.   &
             (.not. (do_ncar_microphys .or. do_mg_ncar_microphys) ) ) then
        call error_mesg ('ls_cloud_microphys_mod/ls_cloud_microphys_init',&
        'when clubb is activated, must use microphys_scheme = "ncar" or &
                      & "mg_ncar"', FATAL)
      endif

      if (override_ice_num == 1) then
        if ( use_Meyers .and. use_Cooper ) &
          call error_mesg (   &
               'ls_cloud_microphysics/ls_cloud_microphysics_init',&
                  'use_Meyers and use_Cooper cannot both be true',FATAL)
      endif

!-----------------------------------------------------------------------
!    initialize the active microphysics scheme module.
!-----------------------------------------------------------------------
      if (doing_prog_clouds) then

!-----------------------------------------------------------------------
!  rotstayn-klein microphysics
!-----------------------------------------------------------------------
        if (do_rk_microphys) then
          call mpp_clock_begin (rk_micro_init_clock)
          call rotstayn_klein_microp_init (Nml_lsc, Exch_ctrl)
          call mpp_clock_end   (rk_micro_init_clock)

!-----------------------------------------------------------------------
!  lin cloud microphysics
!-----------------------------------------------------------------------
        else if (do_lin_cld_microphys) then
          call mpp_clock_begin (lin_micro_init_clock)
          if (Exch_ctrl%do_liq_num) call error_mesg   &
              ('ls_cloud_microphysics/ls_cloud_microphysics_init',  &
               'do_lin_cld_microphys cannot be active with prognostic &
                                     &droplet scheme (do_liq_num)', FATAL)
          call lin_cld_microphys_init    &
                 (id, jd, kd, axes, Time, Physics_control%hydrostatic,  &
                                         Physics_control%phys_hydrostatic)

!------------------------------------------------------------------------
!    define top model level at which Lin microphysics is active (10 hPa).
!------------------------------------------------------------------------
          ktop     = 1
          do k = 1, kd
            if (pref(k) > lin_microphys_top_press) then
              ktop    = k
              exit
            endif
          enddo
          if (mpp_pe() == mpp_root_pe()) &
                write(*,*) 'Top layer for lin_cld_microphys=',   &
                                                  ktop, pref(ktop)
          call mpp_clock_end   (lin_micro_init_clock)

!-----------------------------------------------------------------------
!  morrison-gettelman microphysics (as done by M. Salzmann)
!-----------------------------------------------------------------------
        else if (do_mg_microphys) then
          call mpp_clock_begin (ncar_micro_init_clock)
          call morrison_gettelman_microp_init (Nml_lsc, Exch_ctrl)
          call mpp_clock_end (ncar_micro_init_clock)

!-----------------------------------------------------------------------
!   an early version of the ncar microphysics 
!-----------------------------------------------------------------------
        else if (do_mg_ncar_microphys) then
          call mpp_clock_begin (ncar_micro_init_clock)
          call ini_micro (GRAV, RDGAS, RVGAS, CP_AIR, TFREEZE, HLV, HLF, &
                          Nml_lsc, Exch_ctrl)
          call mpp_clock_end (ncar_micro_init_clock)

!-----------------------------------------------------------------------
!  latest available ncar microphysics  (ncar v1.5)
!-----------------------------------------------------------------------
        else if (do_ncar_microphys) then
          call mpp_clock_begin (ncar_micro_init_clock)
          call micro_mg_init (r8, GRAV, RDGAS, RVGAS, CP_AIR, TFREEZE, &
                              HLV, HLF, Nml_lsc%do_ice_nucl_wpdf,   &
                              errstring, Exch_ctrl)
          if (trim(errstring) /= '') then
            call error_mesg    &
                  ('ls_cloud_microphysics/ls_cloud_microphysics_init', &
                                                         errstring, FATAL)
          endif
          call mpp_clock_end (ncar_micro_init_clock)

!-----------------------------------------------------------------------
!  no valid microphys scheme chosen
!-----------------------------------------------------------------------
        else
           call error_mesg   &
              ('ls_cloud_microphysics/ls_cloud_microphysics_init', &
               'invalid microphys_scheme option in lscloud_driver nml',  &
                                                                     FATAL)
        endif
      endif  ! (doing_prog_clouds)

!-----------------------------------------------------------------------
!    even if Lin cloud microphysics is not active, still need to 
!    initialize its tables.
!-----------------------------------------------------------------------
      if (.not. do_lin_cld_microphys) then
        call setup_con
      endif

!-------------------------------------------------------------------------
!    get namelist initial time from namelist to determine whether 
!    it is time for microphysics to be active.
!-------------------------------------------------------------------------
      Time_init = set_date( init_date(1), init_date(2), init_date(3),  &
                            init_date(4), init_date(5), init_date(6) )
      call get_time( Time_init, current_sec0, current_days0)

!------------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------------


end subroutine ls_cloud_microphysics_init     





!########################################################################

subroutine ls_cloud_microphysics_time_vary (dtcloud_in)

real, intent(in) :: dtcloud_in

!-----------------------------------------------------------------------
!    define current timestep and its inverse.
!-----------------------------------------------------------------------
      dtcloud = dtcloud_in
      inv_dtcloud = 1.0/dtcloud

!----------------------------------------------------------------------

end subroutine ls_cloud_microphysics_time_vary 


!########################################################################

subroutine ls_cloud_microphysics (    &
                 is, ie, js, je, Time, dt, Input_mp, Output_mp, C2ls_mp,&
                 Tend_mp, Lsdiag_mp, Lsdiag_mp_control, Atmos_state,   &
                 Cloud_state, Particles, Precip_state, Cloud_processes, &
                                                      Removal_mp, Aerosol)

!-----------------------------------------------------------------------

type(mp_input_type),        intent(inout)        :: Input_mp
type(mp_output_type),       intent(inout)        :: Output_mp
type(mp_removal_type),      intent(inout)        :: Removal_mp
type(mp_conv2ls_type),      intent(inout)        :: C2ls_mp
type(mp_tendency_type),     intent(inout)        :: Tend_mp
type(mp_lsdiag_type),       intent(inout)        :: Lsdiag_mp
type(mp_lsdiag_control_type), intent(inout)      :: Lsdiag_mp_control
type(atmos_state_type),     intent(inout)        :: Atmos_state
type(cloud_state_type),     intent(inout)        :: Cloud_state
type(particles_type),       intent(inout)        :: Particles
type(precip_state_type),    intent(inout)        :: Precip_state
type(cloud_processes_type), intent(inout)        :: Cloud_processes
type(time_type),            intent(in)           :: Time
integer,                    intent(in)           :: is, ie, js, je
real,                       intent(in)           :: dt
type(aerosol_type),         intent(in), optional :: Aerosol

!------------------------------------------------------------------------
!   local variables:

      real, dimension(size(Input_mp%tin,1),    &
                            size(Input_mp%tin,2)) ::    &
                               ice_lin, graupel_lin,   &
                               enth_micro_col,  wat_micro_col

      real, dimension( size(Input_mp%tin,1), size(Input_mp%tin,2),   &
                            size(Input_mp%tin,3))   ::   &
                               delp, delz, &
                               ST_micro, SQ_micro, SL_micro, SI_micro, &
                               SN_micro, SNI_micro, D_eros_l, D_eros_i,  &
                               nerosc, nerosi, dqcdt, dqidt, qa_new, &
                               ssat_disposal, ql_new,  qi_new,           &
                               nctend, nitend, qn_new, qni_new, &
                               rho, liqcldf, icecldf, tmp2s,  &
                               accre_enhann, tnd_qsnown, &
                               tnd_nsnown, re_icen, relvarn, &           
                               crystal1, rho_air,  &
                               aerosols_concen, droplets_concen

      real, dimension( size(Input_mp%tin,1), size(Input_mp%tin,2),   &
                            size(Input_mp%tin,3),4) ::  &
                               rbar_dust_4bin, ndust_4bin

      integer,dimension(:),allocatable  :: mgcols         

      integer               :: mgncol
      integer               :: i, j, k, n
      integer               :: ix, jx, kx
      integer               :: nlev
      character(len=128)    :: errstring
      real                  :: current_total_sec
      integer               :: current_sec, current_days
      real                  :: depth

!-------------------------------------------------------------------------
!   define array dimensions
!-------------------------------------------------------------------------
      ix = size(Input_mp%tin,1) 
      jx = size(Input_mp%tin,2) 
      kx = size(Input_mp%tin,3)

!------------------------------------------------------------------------
!   call selected microphysics scheme.
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!    Rotstayn-Klein microphysics
!------------------------------------------------------------------------
      if (do_rk_microphys) then
        call mpp_clock_begin (rk_micro_clock)
        call rotstayn_klein_microp ( &
                      ix, jx, kx, Particles%N3D, total_activation,  &  
                      dtcloud, inv_dtcloud, Input_mp%pfull,&
                      Input_mp%pmass, Atmos_state%airdens,     &
                      Atmos_state%esat0, Cloud_state%ql_in,  &
                      Cloud_state%qi_in, Cloud_state%qa_in,   &
                      Cloud_state%ql_mean, Cloud_state%qa_mean, &
                      Cloud_state%qn_mean, Input_mp%omega,  &
                      Input_mp%tin, Atmos_state%U_ca, &
                      Input_mp%qin, Atmos_state%qs,  &
                      Cloud_processes%D_eros, Cloud_processes%dcond_ls, &
                      Cloud_processes% dcond_ls_ice,         &
                      Cloud_processes%qvg, Atmos_state%gamma,   &
                      Cloud_processes%delta_cf, Particles%drop1,    &
                      Particles%concen_dust_sub, Cloud_state%ql_upd,   &
                      Cloud_state%qi_upd, Cloud_state%qn_upd,       & 
                      Cloud_state%qi_mean, Cloud_state%qa_upd,   &
                      C2ls_mp%convective_humidity_area,   &
                      Lsdiag_mp_control%n_diag_4d, Lsdiag_mp%diag_4d,   &
                      Lsdiag_mp_control%diag_id,     &
                      Lsdiag_mp_control%diag_pt,  &
                      Lsdiag_mp_control%n_diag_4d_kp1,   &
                      Lsdiag_mp%diag_4d_kp1,  &
                      limit_conv_cloud_frac, &
                      Cloud_state%SA_out, Cloud_state%SN_out,        & 
                      Tend_mp%ttnd, Tend_mp%qtnd, Cloud_state%SL_out,  &
                      Cloud_state%SI_out, Removal_mp%rain3d,  &
                      Removal_mp%snow3d, Removal_mp%snowclr3d,    &
                      Precip_state%surfrain, Precip_state%surfsnow,  &
                      Cloud_processes%f_snow_berg )  

!-----------------------------------------------------------------------
!   define the output tendency fields.
!-----------------------------------------------------------------------
        Tend_mp%q_tnd(:,:,:,nql) = Cloud_state%SL_out
        Tend_mp%q_tnd(:,:,:,nqi) = Cloud_state%SI_out
        Tend_mp%q_tnd(:,:,:,nqa) = Cloud_state%SA_out
        if (nqn /= NO_TRACER) &
          Tend_mp%q_tnd(:,:,:,nqn) = Cloud_state%SN_out(:,:,:)
        if (nqni /= NO_TRACER) &
          Tend_mp%q_tnd(:,:,:,nqni) = Cloud_state%SNI_out(:,:,:)
        call mpp_clock_end   (rk_micro_clock)

!-----------------------------------------------------------------------
!   lin cld microphysics is activated
!-----------------------------------------------------------------------
      else if (do_lin_cld_microphys ) then
        call mpp_clock_begin (lin_micro_clock)
        do k=1,kx
          delp(:,:,k) =  &
                   Input_mp%phalf(:,:,k+1) - Input_mp%phalf(:,:,k)
          delz(:,:,k) =  &
                  (Input_mp%zhalf(:,:,k+1) - Input_mp%zhalf(:,:,k))*&
                                  Input_mp%tin(:,:,k)/Input_mp%tm(:,:,k)
        end do

!------------------------------------------------------------------------
!   droplet number concentration in #/cc following Boucher and Lohmann 1995
!   the lin microphysics uses #/cc                                        
!------------------------------------------------------------------------
        do k=1,kx
          do j=1,jx
            do i=1,ix
              depth = (Input_mp%phalf(i,j,k+1) - Input_mp%phalf(i,j,k)) / &
                     ( 9.8*0.029*Input_mp%pfull(i,j,k     )/(8.314*  &
                                                   Input_mp%tin(i,j,k)))
!sulfate concentration in ug/m3 
              aerosols_concen(i,j,k)=(Aerosol%aerosol(i,j,k,1))/depth*1.e9
              if ( Input_mp%lat(i,j) < -1.0472 ) then    ! South of ~60S are treated as ocean
                droplets_concen(i,j,k)= 10.**2.06*  &
                                     (0.7273*aerosols_concen(i,j,k))**0.48
              else
                droplets_concen(i,j,k)= Input_mp%land(i,j) *(10.**2.24*   &
                                 (0.7273*aerosols_concen(i,j,k))**0.257)+ &
                                    (1.-Input_mp%land(i,j))* (10.**2.06*  &
                                     (0.7273*aerosols_concen(i,j,k))**0.48)
              endif
            enddo
          enddo
        end do

        call lin_cld_microphys_driver(       &
                Input_mp%qin,  Input_mp%tracer(:,:,:,nql),   &
                Input_mp%tracer(:,:,:,nqr), Input_mp%tracer(:,:,:,nqi), &
                Input_mp%tracer(:,:,:,nqs), Input_mp%tracer(:,:,:,nqg), &
                Input_mp%tracer(:,:,:,nqa), droplets_concen,Tend_mp%qtnd, &
                Tend_mp%q_tnd(:,:,:,nql), Tend_mp%q_tnd(:,:,:,nqr), &
                Tend_mp%q_tnd(:,:,:,nqi), Tend_mp%q_tnd(:,:,:,nqs), &
                Tend_mp%q_tnd(:,:,:,nqg), tend_mp%q_tnd(:,:,:,nqa),    &
                Tend_mp%ttnd, Input_mp%tin,  Input_mp%w, Input_mp%uin, &
                Input_mp%vin, Output_mp%udt, Output_mp%vdt, delz, delp, &
                Input_mp%area, dt, Input_mp%land, Precip_state%surfrain, &
                Precip_state%surfsnow, ice_lin, graupel_lin, &
                is, ie, js, je, 1, kx, ktop, kx, Time)

!-----------------------------------------------------------------------
! Add all "solid" form of precipitation into surf_snow
!-----------------------------------------------------------------------
        Precip_state%surfsnow = (Precip_state%surfsnow + ice_lin +   &
                                                 graupel_lin) * dt/86400.
        Precip_state%surfrain =  Precip_state%surfrain * dt/86400.

!-----------------------------------------------------------------------
! Update tendencies:
!-----------------------------------------------------------------------
        Output_mp%rdt(:,:,:,nqr) = Output_mp%rdt(:,:,:,nqr) +    &
                                               Tend_mp%q_tnd(:,:,:,nqr)
        Output_mp%rdt(:,:,:,nqs) = Output_mp%rdt(:,:,:,nqs) +   &
                                               Tend_mp%q_tnd(:,:,:,nqs)
        Output_mp%rdt(:,:,:,nqg) = Output_mp%rdt(:,:,:,nqg) +   &
                                               Tend_mp%q_tnd(:,:,:,nqg)
        Tend_mp%ttnd             = Tend_mp%ttnd * dt
        Tend_mp%qtnd             = Tend_mp%qtnd * dt
        Tend_mp%q_tnd(:,:,:,nql) = Tend_mp%q_tnd(:,:,:,nql) * dt
        Tend_mp%q_tnd(:,:,:,nqi) = Tend_mp%q_tnd(:,:,:,nqi) * dt
        Tend_mp%q_tnd(:,:,:,nqa) = Tend_mp%q_tnd(:,:,:,nqa) * dt

!-----------------------------------------------------------------------
! Update rain_wat, snow_wat, graupel_wat
!-----------------------------------------------------------------------
        Input_mp%tracer(:,:,:,nqr) = Input_mp%tracer(:,:,:,nqr) +   &
                                              Tend_mp%q_tnd(:,:,:,nqr)*dt
        Input_mp%tracer(:,:,:,nqs) = Input_mp%tracer(:,:,:,nqs) +   &
                                              Tend_mp%q_tnd(:,:,:,nqs)*dt
        Input_mp%tracer(:,:,:,nqg) = Input_mp%tracer(:,:,:,nqg) +   &
                                              Tend_mp%q_tnd(:,:,:,nqg)*dt
        call mpp_clock_end   (lin_micro_clock)

!-----------------------------------------------------------------------
!  NCAR microphysics (currently 3 flavors)
!-----------------------------------------------------------------------
      else if (do_mg_microphys   .or. &    
               do_mg_ncar_microphys   .or. &    
               do_ncar_microphys) then      
        call mpp_clock_begin (ncar_micro_clock)

!-----------------------------------------------------------------------
!   determine whether NCAR microphysics are active at the current time
!-----------------------------------------------------------------------
        call get_time( time, current_sec, current_days)
        current_total_sec = real(current_sec - current_sec0) +   &
                                    86400.0*(current_days - current_days0)
   
        if (current_total_sec >= micro_begin_sec) then
          Atmos_state%tn = Input_mp%tin + Tend_mp%ttnd
          Atmos_state%qvn = Input_mp%qin + Tend_mp%qtnd

!--------------------------------------------------------------------------
!     define some input fields related to ls condensation and the  cloud 
!     erosion process that are needed when tiedtke macrophysics are active,
!     since the magnitude of these processes is still subject to change
!     based on what the microphysics does. for the non-tiedtke case, the 
!     magnitude of these processes have been locked in before the micro-
!     physics tendencies are calculated, and so these input fields are 
!     set to 0.0.
!--------------------------------------------------------------------------
          do k=1,kx
            do j=1,jx
              do i=1,ix
                Cloud_processes%dcond_ls_tot(i,j,k) =   &
                               Cloud_processes%dcond_ls(i,j,k) +   &
                                      Cloud_processes%dcond_ls_ice(i,j,k) 
                if (tiedtke_macrophysics) then
                  D_eros_i(i,j,k) = -Cloud_state%qi_upd(i,j,k)* &
                                        Cloud_processes%D_eros(i,j,k)/ &
                                                                  dtcloud
                  D_eros_l(i,j,k) = -Cloud_state%ql_upd(i,j,k)* &
                                        Cloud_processes%D_eros(i,j,k)/ &
                                                                  dtcloud
                  if (Cloud_state%ql_upd(i,j,k) >= qmin) then
                    nerosc(i,j,k) = D_eros_l(i,j,k)/  &
                                      Cloud_state%ql_upd(i,j,k)* &
                               Cloud_state%qn_upd(i,j,k)/MAX(0.0001, &
                                               Cloud_state%qa_upd(i,j,k))
                  else
                    nerosc(i,j,k) = 0.
                  endif
                  if (Cloud_state%qi_upd(i,j,k) >= qmin) then
                    nerosi(i,j,k) = D_eros_i(i,j,k)/   &
                                          Cloud_state%qi_upd(i,j,k)* &
                               Cloud_state%qni_upd(i,j,k)/MAX(0.0001, &
                                              Cloud_state%qa_upd(i,j,k))
                  else
                    nerosi(i,j,k) = 0.
                  endif
                  if (Cloud_processes%dcond_ls_tot(i,j,k) > 0.) then
                    if (Atmos_state%tn(i,j,k) <= (tfreeze - 40.) ) then
                      dqcdt (i,j,k) = 0.
                      dqidt(i,j,k) = Cloud_processes%dcond_ls_tot(i,j,k)* &
                                                              inv_dtcloud
                    else
                      dqidt (i,j,k) = 0.
                      dqcdt(i,j,k) = Cloud_processes%dcond_ls_tot(i,j,k)* &
                                                              inv_dtcloud
                    endif
                  else
                    if (Atmos_state%tn(i,j,k) <= tfreeze) then
                      dqcdt(i,j,k) = MAX(  &
                                     Cloud_processes%dcond_ls_tot(i,j,k),&
                                          -Cloud_state%ql_upd(i,j,k))
                      dqidt(i,j,k) = MAX(   &
                                     Cloud_processes%dcond_ls_tot(i,j,k) -&
                                                        dqcdt(i,j,k),   &
                                                -Cloud_state%qi_upd(i,j,k))
                      dqcdt(i,j,k) = dqcdt(i,j,k)*inv_dtcloud
                      dqidt(i,j,k) = dqidt(i,j,k)*inv_dtcloud
                    else
                      dqidt(i,j,k) = 0.
                      dqcdt(i,j,k) = MAX(   &
                                      Cloud_processes%dcond_ls_tot(i,j,k),&
                                           -Cloud_state%ql_upd(i,j,k))* &
                                                              inv_dtcloud
                    endif
                  endif
                else ! (tiedtke)
                  dqidt(i,j,k) = 0.
                  dqcdt(i,j,k) = 0.
                  nerosi(i,j,k) = 0.
                  nerosc(i,j,k) = 0.
                  D_eros_l(i,j,k) = 0.                             
                  D_eros_i(i,j,k) = 0.                             
                endif  ! (tiedtke)
              end do
            end do   
          end do   

!------------------------------------------------------------------------
!    if the 'mg' microphysics (the original,produced by M. Salzmann, with
!    tweaks by H. Guo and R. Hemler) is activated, execute the following:
!------------------------------------------------------------------------
          if (do_mg_microphys) then

!-----------------------------------------------------------------------
!   define activated droplets in units of #/kg (drop1 is #/cc).
!-----------------------------------------------------------------------
            Particles%drop2 = Particles%drop1*1.e6/Atmos_state%airdens

!-----------------------------------------------------------------------
!   execute the microphysics, 1 jrow at a time.
!-----------------------------------------------------------------------
            do j=1,jx

!------------------------------------------------------------------------
!    if debugging is activated, output the temp tendency prior to
!    microphysics.
!------------------------------------------------------------------------
              call write_debug_output (" ST samp bef mg ",   &
                                                       Tend_mp%ttnd, j=j)

!-------------------------------------------------------------------------
!    call morrison-gettelman (mg) microphysics package.
!-------------------------------------------------------------------------
              call morrison_gettelman_microp( &
                   tiedtke_macrophysics, total_activation,   &
                   dqa_activation, j ,ix, jx, kx, dtcloud,   &
                   Input_mp%pfull(:,j,:),  Atmos_state%delp(:,j,:),  &
                   Atmos_state%tn(:,j,:),  Input_mp%tin(:,j,:),    &
                   Atmos_state%qvn(:,j,:), Input_mp%qin(:,j,:),  &
                   Cloud_state%ql_upd(:,j,:), Cloud_state%qi_upd(:,j,:),&
                   Cloud_state%qn_upd(:,j,:), Cloud_state%qni_upd(:,j,:), &
                   Cloud_state%qa_upd(:,j,:), dqcdt(:,j,:), dqidt(:,j,:), &
                   Particles%drop2(:,j,:), Particles%crystal1(:,j,:), &
                   Particles%rbar_dust(:,j,:), Particles%ndust(:,j,:),  &
                   Cloud_processes%delta_cf(:,j,:),   &
                   Cloud_state%qa_upd(:,j,:), Cloud_state%qa_upd_0(:,j,:),&
                   Cloud_state%SA_0(:,j,:), D_eros_l(:,j,:),  &
                   nerosc(:,j,:),  D_eros_i(:,j,:), nerosi(:,j,:), &
                   Atmos_state%gamma(:,j,:), inv_dtcloud,    &
                   Cloud_state%qa_in(:,j,:), Tend_mp%ttnd(:,j,:),   &
                   Tend_mp%qtnd(:,j,:), ssat_disposal(:,j,:), &
                   ST_micro(:,j,:), SQ_micro(:,j,:),&
                   SL_micro(:,j,:), SI_micro(:,j,:),  &
                   SN_micro(:,j,:), SNI_micro(:,j,:),&
                   Cloud_state%SA_out(:,j,:), Removal_mp%rain3d,   &
                   Removal_mp%snow3d, Precip_state%surfrain(:,j),   &
                   Precip_state%surfsnow(:,j), &
                   Precip_state%lsc_rain(:,j,:),   &
                   Precip_state%lsc_snow(:,j,:), &
                   Precip_state%lsc_rain_size(:,j,:),   &
                   Precip_state%lsc_snow_size(:,j,:), &
                   Cloud_processes%f_snow_berg(:,j,:), &
                   Lsdiag_mp_control%n_diag_4d, Lsdiag_mp%diag_4d,   &
                   Lsdiag_mp_control%diag_id, Lsdiag_mp_control%diag_pt)    

!------------------------------------------------------------------------
!    if debugging is activated, output the temp tendency after 
!    microphysics is completed.
!------------------------------------------------------------------------
              call write_debug_output  &
                                  (" ST samp aft mg ", Tend_mp%ttnd, j=j)
            end do   ! j loop

!-------------------------------------------------------------------------
!    executed for mg_ncar or ncar microphysics:
!-------------------------------------------------------------------------
          else if (do_mg_ncar_microphys .or. do_ncar_microphys) then   
          
            rho = Input_mp%pfull/(RDGAS*Atmos_state%tn)

!------------------------------------------------------------------------
!   define amount of activated aerosol to be passed to microphysics.
!------------------------------------------------------------------------
            if (do_clubb > 0 ) then

!------------------------------------------------------------------------
!   for CLUBB, activated ice crystals are supplied by 
!   Particles%Icedrop_act_CLUBB, in units of #/kg. since droplet activation
!   has been considered and resultant droplet number has been updated 
!   in CLUBB, we do not need to include activation within MG microphysics,
!   i.e Particles%drop2 = 0.0.
!------------------------------------------------------------------------
              crystal1 = Particles%Icedrop_act_CLUBB  
              Particles%drop2 = 0.0

!------------------------------------------------------------------------
!   Options to override ice crystal number concentrations from CLUBB:
!------------------------------------------------------------------------
              if (override_ice_num == 1) then
                rho_air = Input_mp%pfull/Atmos_state%tn/RDGAS

                if (use_Meyers) then
!------------------------------------------------------------------------
!   Meyers formula as in original MG microphysics
!-------------------------------------------------------------------------
                  crystal1 = 1000.0*exp(  &
                             (12.96*0.0125*(273.15-Atmos_state%tn))-0.639)
                  crystal1 = crystal1/rho_air

                elseif( use_Cooper ) then
!--------------------------------------------------------------------
! cooper curve (factor of 1000 is to convert from L-1 to m-3)
! put limit on number of nucleated crystals, set to number at T=-35 C
! then convert from m-3 to kg-1.
!--------------------------------------------------------------------
                  crystal1 = 0.005*exp(0.304*(273.15-Atmos_state%tn))*1000.
                  crystal1 = min( crystal1, 208.9e3)/rho_air
                else 

!-----------------------------------------------------------------------
! use a constant 0.5 /kg
!-----------------------------------------------------------------------
                  crystal1 = 1.0e6 * 0.5
                end if
              end if
  
            else ! (do_clubb)
!-----------------------------------------------------------------------
!   for the non-CLUBB case, use the values previously calculated and 
!   input to this routine. convert to units of #/kg.
!-----------------------------------------------------------------------
              crystal1 = Particles%crystal1/rho

!-----------------------------------------------------------------------
!   define activated droplets in units of #/kg (drop1 is #/cc).
!-----------------------------------------------------------------------
              Particles%drop2 = Particles%drop1*1.e6/Atmos_state%airdens
            endif  ! (do_clubb)


!------------------------------------------------------------------------
!   set liquid and ice cloud fraction to be the same as total large-scale 
!   cloud fraction.
!------------------------------------------------------------------------
            liqcldf  = Cloud_state%qa_upd
            icecldf  = Cloud_state%qa_upd

!------------------------------------------------------------------------
!    are these the actual bin centers, or arbitrary values ??
!    radius of 4 dust bins for contact freezing (in the unit of m)
!------------------------------------------------------------------------
            do k = 1,kx
              do j = 1,jx
                do i = 1,ix
                  rbar_dust_4bin(i,j,k,1) = 5.e-6
                  rbar_dust_4bin(i,j,k,2) = 10.e-6
                  rbar_dust_4bin(i,j,k,3) = 15.e-6
                  rbar_dust_4bin(i,j,k,4) = 20.e-6

!------------------------------------------------------------------------
!    define the number of particles in each of the 4 dust bins, if 
!    contact freezing is to be done. the active code below assigns
!    the total number to each size bin, as is done with CLUBB.
!    is this OK, or should the total number be distributed across all the
!    bins? (commented code distributes total equally across bins)
!------------------------------------------------------------------------
                  if ( use_ndust ) then
                    ndust_4bin(i,j,k,1)     = Particles%ndust(i, j, k)
                    ndust_4bin(i,j,k,2)     = Particles%ndust(i, j, k)
                    ndust_4bin(i,j,k,3)     = Particles%ndust(i, j, k)
                    ndust_4bin(i,j,k,4)     = Particles%ndust(i, j, k)
!                   ndust_4bin(i,j,k,1) = 0.25*Particles%ndust(i,j,k) 
!                   ndust_4bin(i,j,k,2) = 0.25*Particles%ndust(i,j,k) 
!                   ndust_4bin(i,j,k,3) = 0.25*Particles%ndust(i,j,k) 
!                   ndust_4bin(i,j,k,4) = 0.25*Particles%ndust(i,j,k) 

!------------------------------------------------------------------------
!    if contact freezing not desired, set ndust = 0. in each bin.
!------------------------------------------------------------------------
                  else
                    ndust_4bin(i,j,k,1)     = 0.0
                    ndust_4bin(i,j,k,2)     = 0.0
                    ndust_4bin(i,j,k,3)     = 0.0
                    ndust_4bin(i,j,k,4)     = 0.0
                  endif
                enddo
              enddo
            enddo

!------------------------------------------------------------------------
!    define the relative variance of the cloud water within each gridbox.
!    when CLUBB is active, spatially-dependent values are returned from 
!    CLUBB; with Tiedtke macrophysics up to this time only a constant 
!    value has been used, though spatial dependence could be introduced.
!------------------------------------------------------------------------
            if (do_clubb > 0 ) then
              relvarn(:,:,:) = Cloud_state%qcvar_clubb(:,:,:)
            else
              relvarn(:,:,:) = Cloud_state%relvarn(:,:,:)
            endif

!------------------------------------------------------------------------
!    if the 'mg_ncar' microphysics (a version of the NCAR microphysics 
!    developed by H. Guo and R. Hemler, based upon a newer release than 
!    that used by M. Salzmann, but following his adaptations for use in
!    GFDL models) is activated, execute the following:
!------------------------------------------------------------------------
            if (do_mg_ncar_microphys) then   

!-----------------------------------------------------------------------
!    execute the microphysics, 1 jrow at a time.
!-----------------------------------------------------------------------
              do j=1,jx
                call mmicro_pcond( &
                     dqa_activation, total_activation,    &
                     tiedtke_macrophysics, .false., j ,jx, kx, ix, ix,  &
                     dtcloud, relvarn(:,j,:), Atmos_state%tn(:,j,:),     &
                     Atmos_state%qvn(:,j,:),  Cloud_state%ql_upd(:,j,:), &
                     Cloud_state%qi_upd(:,j,:), Cloud_state%qn_upd(:,j,:),&
                     Cloud_state%qni_upd(:,j,:), Input_mp%pfull(:,j,:),  &
                     Atmos_state%delp(:,j,:), Input_mp%phalf(:,j,:), &
                     Cloud_state%qa_upd(:,j,:), liqcldf(:,j,:),   &
                     icecldf(:,j,:), Cloud_processes%delta_cf(:,j,:), &
                     D_eros_l(:,j,:), nerosc(:,j,:), &
                     D_eros_i(:,j,:), nerosi(:,j,:), &
                     dqcdt(:,j,:), dqidt(:,j,:), crystal1(:,j,:), &
                     Particles%drop2(:,j,:), rbar_dust_4bin(:,j,:,:), &
                     ndust_4bin(:,j,:,:), &
                     ST_micro(:,j,:), SQ_micro(:,j,:), SL_micro(:,j,:), &
                     SI_micro(:,j,:), SN_micro(:,j,:), SNI_micro(:,j,:), &
                     Precip_state%surfrain(:,j),   &
                     Precip_state%surfsnow(:,j),   &
                     Removal_mp%rain3d(:,j,:), Removal_mp%snow3d(:,j,:), &
                     Precip_state%lsc_rain(:,j,:),   &
                     Precip_state%lsc_snow(:,j,:), &
                     Precip_state%lsc_rain_size(:,j,:),   &
                     Precip_state%lsc_snow_size(:,j,:), &
                     Cloud_processes%f_snow_berg(:,j,:), &
                     Cloud_state%qa_in(:,j,:), Atmos_state%gamma(:,j,:),&
                     Cloud_state%SA_0(:,j,:), Cloud_state%SA_out(:,j,:),  &
                     ssat_disposal (:,j,:), Lsdiag_mp_control%n_diag_4d,  &
                     Lsdiag_mp%diag_4d, Lsdiag_mp_control%diag_id,    &
                     Lsdiag_mp_control%diag_pt)
              end do

!------------------------------------------------------------------------
!    if the 'ncar' microphysics (the newest available version of the NCAR 
!    microphysics, adapted for use in FMS by H. Guo and R. Hemler) 
!    is activated, execute the following:
!------------------------------------------------------------------------
            else if (do_ncar_microphys) then 

!------------------------------------------------------------------------
!    define the topmost model level at which microphysics is to be 
!    calculated (top_lev).  define the number of levels (nlev) over which 
!    microphysics will be calculated (from top_lev to the surface).
!------------------------------------------------------------------------
              nlev = kx - top_lev + 1

!-----------------------------------------------------------------------
!    define additional input fields:
!    accre_enhann -- accretion enhancement factor
!    the following are used if an external cirrus microphysics model 
!    is active (eg, NCAR CARMA model)
!    tnd_qsnown --  snow mass tendency (kg/kg/s)
!    tnd_nsnown(:,:) ! snow number tendency (#/kg/s)
!    re_icen(:,:)    ! ice effective radius (m)
!-----------------------------------------------------------------------
              accre_enhann(:,:,:) = 1.0 
              tnd_qsnown(:,:,:) = 0.     
              tnd_nsnown(:,:,:) = 0.
              re_icen(:,:,:) = 0.

!-------------------------------------------------------------------------
!    execute the microphysics, 1 jrow at a time.
!-------------------------------------------------------------------------
              do j=1,jx

!------------------------------------------------------------------------
!    call subroutine micro_mg_get_cols to identify the columns in which
!    microphysics will be calculated. only those columns meeting certain 
!    criteria below the specified top_lev will be flagged, and experience
!    microphysics.  
!------------------------------------------------------------------------
                call micro_mg_get_cols (   &
                    ix, nlev, top_lev, Atmos_state%qvn(:,j,:), &
                    Cloud_state%ql_upd(:,j,:) + dqcdt(:,j,:)*dtcloud, &
                    Cloud_state%qi_upd(:,j,:) + dqidt(:,j,:)*dtcloud, &
                    mgncol, mgcols, do_clubb > 0)

!------------------------------------------------------------------------
!    if debugging is activated, output the temp tendency prior to
!    microphysics.
!------------------------------------------------------------------------
                call write_debug_output (" ST samp bef mg ",   &
                                                      Tend_mp%ttnd, j=j)

!-------------------------------------------------------------------------
!    call the ncar microphysics routine micro_mg_tend.
!-------------------------------------------------------------------------
                call  micro_mg_tend ( &
                       dqa_activation, total_activation, &
                       tiedtke_macrophysics, j, jx, mgncol, mgcols,   &
                       nlev, top_lev, dtcloud, Atmos_state%tn(:,j,:), &
                       Atmos_state%qvn(:,j,:), Cloud_state%ql_upd(:,j,:), &
                       Cloud_state%qi_upd(:,j,:),   &
                       Cloud_state%qn_upd(:,j,:),   &
                       Cloud_state%qni_upd(:,j,:), relvarn(:,j,:),  &
                       accre_enhann(:,j,:), Input_mp%pfull(:,j,:),  &
                       Atmos_state%delp(:,j,:), Input_mp%phalf(:,j,:), &
                       Cloud_state%qa_upd(:,j,:), liqcldf(:,j,:),   &
                       icecldf(:,j,:), Cloud_processes%delta_cf(:,j,:), &
                       D_eros_l(:,j,:), nerosc(:,j,:), D_eros_i(:,j,:),  &
                       nerosi(:,j,:), dqcdt(:,j,:), dqidt(:,j,:),   &
                       crystal1(:,j,:), Particles%drop2(:,j,:),&
                       rbar_dust_4bin(:,j,:,:),  ndust_4bin(:,j,:,:),   &
                       ST_micro(:,j,:), SQ_micro(:,j,:), SL_micro(:,j,:), &
                       SI_micro(:,j,:), SN_micro(:,j,:), SNI_micro(:,j,:),&
                       Precip_state%surfrain(:,j),    &
                       Precip_state%surfsnow(:,j),&
                       Precip_state%lsc_snow(:,j,:), &
                       Removal_mp%rain3d(:,j,:), Removal_mp%snow3d(:,j,:),&
                       Precip_state%lsc_rain(:,j,:),  &
                       Precip_state%lsc_rain_size(:,j,:),  &
                       Precip_state%lsc_snow_size(:,j,:),  &
                       tnd_qsnown(:,j,:), tnd_nsnown(:,j,:),    &
                       re_icen(:,j,:), errstring,    &
                       Cloud_processes%f_snow_berg(:,j,:), &
                       ssat_disposal(:,j,:), &
                       Lsdiag_mp_control%n_diag_4d, Lsdiag_mp%diag_4d,  &
                       Lsdiag_mp_control%diag_id,    &
                                                 Lsdiag_mp_control%diag_pt)

!------------------------------------------------------------------------
!   convert from effective radius to diameter for use in radiation.
!   in mg and mg-ncar, diameter is returned from microphysics routine,
!   so this step is unneeded.
!------------------------------------------------------------------------
                Precip_state%lsc_rain_size(:,j,:) =   &
                                   2.0*Precip_state%lsc_rain_size(:,j,:)
                Precip_state%lsc_snow_size(:,j,:) =   &
                                   2.0*Precip_state%lsc_snow_size(:,j,:)

!-------------------------------------------------------------------------
!   if an error message was returned from micro_mg_tend output it and
!   stop execution.
!-------------------------------------------------------------------------
                if (trim(errstring) /= '') then
                  call error_mesg (  &
                          'moist_processes/ls_cloud_microphysics', &
                                                         errstring, FATAL)
                endif

!------------------------------------------------------------------------
!   if debugging is activated, output the temp tendency after microphysics.
!------------------------------------------------------------------------
                call write_debug_output (" ST samp aft mg ",  &
                                                       Tend_mp%ttnd, j=j)

              end do  ! end of j loop
     
!------------------------------------------------------------------------
!    calculate column enthalpy and total water changes
!    Note: in MG-microphys, temperature tendency is multiplied by Cp_air.
!------------------------------------------------------------------------
              enth_micro_col(:,:) = 0.0
              wat_micro_col(:,:)  = 0.0
              do j=1,jx
                do i=1,ix
                  do k=1,kx
                    enth_micro_col(i,j) = enth_micro_col(i,j)   +         &
                        ( ST_micro(i,j,k) - HLV*SL_micro(i,j,k) -   &
                                         HLS*SI_micro(i,j,k) )*    &
                                             Atmos_state%delp(i,j,k)/grav

                    wat_micro_col(i,j) = wat_micro_col(i,j)  +            &
                         ( SQ_micro(i,j,k) + SL_micro(i,j,k) +  &
                                        SI_micro(i,j,k) )*   &
                                             Atmos_state%delp(i,j,k)/grav
                  enddo
  
                  enth_micro_col(i,j) = enth_micro_col(i,j) +   &
                                                 (-HLV*1000.0* &
                                        (Precip_state%surfrain(i,j) -  &
                                         Precip_state%surfsnow(i,j)) -  &
                                 HLS*1000.0 * Precip_state%surfsnow(i,j) )

                  wat_micro_col(i,j) = wat_micro_col(i,j) +   &
                                       Precip_state%surfrain(i,j) *1000.0
                enddo
              enddo
            endif ! (do_mg)
          endif  ! (do_mg)

!------------------------------------------------------------------------
!    adjust precip fields to assure mass conservation and realizable
!    values.
!------------------------------------------------------------------------
          call adjust_precip_fields (   &
                              ix, jx, kx, SQ_micro, SL_micro, SI_micro,  &
                                  Atmos_state, Precip_state, Lsdiag_mp, &
                                                       Lsdiag_mp_control ) 

!-----------------------------------------------------------------------
!    update prognostic tendencies due to microphysics terms.
!-----------------------------------------------------------------------
            Tend_mp%qtnd = Tend_mp%qtnd + SQ_micro*dtcloud
            Tend_mp%ttnd = Tend_mp%ttnd + ST_micro/cp_air*dtcloud
            Tend_mp%q_tnd(:,:,:,nql) = Tend_mp%q_tnd(:,:,:,nql) +   &
                                                SL_micro(:,:,:)*dtcloud
            Tend_mp%q_tnd(:,:,:,nqi) = Tend_mp%q_tnd(:,:,:,nqi) +   &
                                                SI_micro(:,:,:)*dtcloud
            Tend_mp%q_tnd(:,:,:,nqn) = Tend_mp%q_tnd(:,:,:,nqn) +   &
                                                SN_micro(:,:,:)*dtcloud
            Tend_mp%q_tnd(:,:,:,nqni)= Tend_mp%q_tnd(:,:,:,nqni) +  &
                                                SNI_micro(:,:,:)*dtcloud

            Cloud_state%SL_out = Cloud_state%SL_out + SL_micro*dtcloud
            Cloud_state%SI_out = Cloud_state%SI_out + SI_micro*dtcloud
            Cloud_state%SN_out = Cloud_state%SN_out + SN_micro*dtcloud
            Cloud_state%SNI_out = Cloud_state%SNI_out + SNI_micro*dtcloud

!------------------------------------------------------------------------
!    adjustment to fields needed after removing supersaturation (only
!    for non-CLUBB case, as originally coded).
!------------------------------------------------------------------------
          if (do_clubb <= 0 ) then
            call adjust_for_supersaturation_removal (  &
                    ix, jx, kx, C2ls_mp, Input_mp, Atmos_state,  &
                         ssat_disposal, Particles, Cloud_state, Lsdiag_mp,&
                                                    Lsdiag_mp_control ) 
          endif

!------------------------------------------------------------------------
!    process output fields after microphysics is completed 
!    for the CLUBB case.
!------------------------------------------------------------------------
          if (do_clubb > 0 ) then

!-------------------------------------------------------------------------
!    remove any clouds with less condensate present than the specified 
!    allowable minimum.
!-------------------------------------------------------------------------
            call destroy_tiny_clouds_clubb (   &
                  ix, jx, kx, Cloud_state, Tend_mp, Lsdiag_mp, &
                        Lsdiag_mp_control, C2ls_mp, Input_mp, Atmos_state)


!-------------------------------------------------------------------------
!     ---> h1g, 06-14-2013, in order to reproduce bit-wise identical 
!     results as AM3-CLUBB
!-------------------------------------------------------------------------
            Tend_mp%ttnd = Tend_mp%ttnd/dtcloud
            Tend_mp%ttnd = Tend_mp%ttnd*dtcloud
     
            Tend_mp%qtnd = Tend_mp%qtnd/dtcloud
            Tend_mp%qtnd = Tend_mp%qtnd*dtcloud

            Tend_mp%q_tnd(:,:,:,nql) = Tend_mp%q_tnd(:,:,:,nql)/dtcloud
            Tend_mp%q_tnd(:,:,:,nql) = Tend_mp%q_tnd(:,:,:,nql)*dtcloud

            Tend_mp%q_tnd(:,:,:,nqn) = Tend_mp%q_tnd(:,:,:,nqn)/dtcloud
            Tend_mp%q_tnd(:,:,:,nqn) = Tend_mp%q_tnd(:,:,:,nqn)*dtcloud

            Tend_mp%q_tnd(:,:,:,nqni) = Tend_mp%q_tnd(:,:,:,nqni)/dtcloud
            Tend_mp%q_tnd(:,:,:,nqni) = Tend_mp%q_tnd(:,:,:,nqni)*dtcloud
    
            Tend_mp%q_tnd(:,:,:,nqi) = Tend_mp%q_tnd(:,:,:,nqi)/dtcloud
            Tend_mp%q_tnd(:,:,:,nqi) = Tend_mp%q_tnd(:,:,:,nqi)*dtcloud

!------------------------------------------------------------------------
!    process fields after microphysics is completed for the non-CLUBB case.
!------------------------------------------------------------------------
          else  !(do_clubb)

!-------------------------------------------------------------------------
!    remove any clouds with less condensate present than the specified 
!    allowable minimum.
!-------------------------------------------------------------------------
            call destroy_tiny_clouds (    &
                   ix, jx, kx, Cloud_state, Tend_mp, Lsdiag_mp, &
                      Lsdiag_mp_control, C2ls_mp, Input_mp, Atmos_state)

!-----------------------------------------------------------------------
!    define output fields.
!-----------------------------------------------------------------------
            Tend_mp%q_tnd(:,:,:,nql) = Cloud_state%SL_out
            Tend_mp%q_tnd(:,:,:,nqi) = Cloud_state%SI_out
            Tend_mp%q_tnd(:,:,:,nqa) = Cloud_state%SA_out
            if (nqn /= NO_TRACER) &
                     Tend_mp%q_tnd(:,:,:,nqn) = Cloud_state%SN_out(:,:,:)
            if (nqni /= NO_TRACER) &
                     Tend_mp%q_tnd(:,:,:,nqni) = Cloud_state%SNI_out(:,:,:)

!-----------------------------------------------------------------------
!    define the total precipitating ice field for use in COSP (stored in
!    Removal_mp%snowclr3d).it is moved to this field so that total cloud
!    ice plus snow is contained in (Removal_mp%snowclr3d + the cloud ice
!    field), as in the R-K microphysics case.
!-----------------------------------------------------------------------
            Removal_mp%snowclr3d = Removal_mp%snow3d

          endif !  do_clubb
        endif !   current_total_sec >= micro_begin_sec

        call mpp_clock_end   (ncar_micro_clock)
!-------------------------------------------------------------------------
!    exit with error if no valid microphysics scheme was specified.
!-------------------------------------------------------------------------
      else    ! do rk
        call error_mesg ('ls_cloud_microphysics/ls_cloud_microphysics', &
              'invalid lscloud_driver_nml microphys_scheme option', FATAL)
      endif    ! (do_rk)

!---------------------------------------------------------------------


end subroutine ls_cloud_microphysics    



!#######################################################################

subroutine ls_cloud_microphysics_end                                

!------------------------------------------------------------------------

      integer   :: rk_micro_term_clock, lin_micro_term_clock, &
                   ncar_micro_term_clock

!------------------------------------------------------------------------

      if (.not. module_is_initialized) return

!------------------------------------------------------------------------
!    define clocks for each microphysics scheme.
!------------------------------------------------------------------------
      if (do_rk_microphys) then
        rk_micro_term_clock = mpp_clock_id(   &
               '   Ls_cld_micro: rk_micro:Termination' , &
                                                grain=CLOCK_MODULE_DRIVER )
      else if (do_lin_cld_microphys) then
        lin_micro_term_clock = mpp_clock_id(     &
               '   Ls_cld_micro: lin_micro:Termination' , &
                                                grain=CLOCK_MODULE_DRIVER )
      else if (do_mg_microphys) then
        ncar_micro_term_clock = mpp_clock_id(    &
               '   Ls_cld_micro: mg_micro:Termination' , &
                                                grain=CLOCK_MODULE_DRIVER)
      else if (do_mg_ncar_microphys) then
        ncar_micro_term_clock = mpp_clock_id(    &
               '   Ls_cld_micro: mg_ncar_micro:Termination' , &
                                                grain=CLOCK_MODULE_DRIVER)
      else if (do_ncar_microphys) then
        ncar_micro_term_clock = mpp_clock_id(    &
               '   Ls_cld_micro: ncar_micro:Termination' , &
                                                grain=CLOCK_MODULE_DRIVER)
      endif

!-------------------------------------------------------------------------
!    call the termination routines for each of the available
!    microphysics schemes which has one.
!-------------------------------------------------------------------------
      if (do_rk_microphys ) then
        call mpp_clock_begin (rk_micro_term_clock)
        call rotstayn_klein_microp_end
        call mpp_clock_end   (rk_micro_term_clock)
      else if (do_mg_microphys ) then
        call mpp_clock_begin (ncar_micro_term_clock)
        call morrison_gettelman_microp_end 
        call mpp_clock_end   (ncar_micro_term_clock)
      else if (do_lin_cld_microphys) then
        call mpp_clock_begin (lin_micro_term_clock)
        call lin_cld_microphys_end
        call mpp_clock_end   (lin_micro_term_clock)
      else if (do_mg_ncar_microphys ) then
        call mpp_clock_begin (ncar_micro_term_clock)
        call mmicro_end
        call mpp_clock_end   (ncar_micro_term_clock)
      endif
      
      module_is_initialized = .false.

!----------------------------------------------------------------------

end subroutine ls_cloud_microphysics_end


!########################################################################

subroutine adjust_precip_fields (    &
              ix, jx, kx, SQ_micro, SL_micro, SI_micro,     &
                                   Atmos_state, Precip_state, Lsdiag_mp, &
                                                 Lsdiag_mp_control )

!------------------------------------------------------------------------
!    subroutine adjust_precip_fields modifies the surface precipitation to
!    balance the atmospheric tendencies of water, and thus conserve water
!    mass. Any needed adjustments are available for examination as
!    model netcdf diagnostics. 
!------------------------------------------------------------------------

integer,                    intent(in)    :: ix, jx, kx
type(mp_lsdiag_type),       intent(inout) :: Lsdiag_mp
type(mp_lsdiag_control_type), intent(inout) :: Lsdiag_mp_control
type(atmos_state_type),     intent(inout) :: Atmos_state
type(precip_state_type),    intent(inout) :: Precip_state
real, dimension (:,:,:),    intent(in)    :: SL_micro, SI_micro, SQ_micro


!----------------------------------------------------------------------
!   local variables:

      real, dimension (ix, jx)     :: m1, m2, scalef
      integer                      :: i, j, k

!----------------------------------------------------------------------
!    if enforcement of water mass conservation is desired, compute
!    the precip reaching the ground and the net change in vapor,
!    cloud water and cloud ice.
!----------------------------------------------------------------------
      if (mass_cons) then
        do j=1,jx
          do i=1,ix
            m1(i,j) = 0.
            do k=1,kx
              m1(i,j) = m1(i,j) +   &
                   (SQ_micro(i,j,k) + SL_micro(i,j,k) + SI_micro(i,j,k))* &
                                    dtcloud*Atmos_state%delp(i,j,k)/grav
            end do
            m2(i,j) = 1.e3*Precip_state%surfrain(i,j)*dtcloud

!------------------------------------------------------------------------
!    for small precip, adjustment for conservation may be ignored. other-
!    wise, compute the ratio of condensate loss to precip at the surface.
!------------------------------------------------------------------------
            if ( (do_clubb > 0 .and.   &
                       m2(i,j) .GT. min_precip_needing_adjustment) .or. &
                   (do_clubb == 0 .and. m2(i,j) .ne. 0.0)) THEN
              scalef(i,j) = -m1(i,j)/m2(i,j)
 
!-----------------------------------------------------------------------
!   define diagnostics capturing the rate (kg/m2/s) that the precip
!   field is adjusted to balance the loss of atmospheric water mass.
!-----------------------------------------------------------------------
              if (Lsdiag_mp_control%diag_id%rain_mass_conv > 0   ) &
              Lsdiag_mp%diag_4d(i,j,1,   &
                            Lsdiag_mp_control%diag_pt%rain_mass_conv) = &
                           (scalef(i,j)*Precip_state%surfrain(i,j) -    &
                                        Precip_state%surfrain(i,j))*1.0e3
              if (Lsdiag_mp_control%diag_id%snow_mass_conv > 0   ) &
              Lsdiag_mp%diag_4d(i,j,1,   &
                           Lsdiag_mp_control%diag_pt%snow_mass_conv) = &
                              (scalef(i,j)*Precip_state%surfsnow(i,j) -  &
                                         Precip_state%surfsnow(i,j))*1.0e3
 
!------------------------------------------------------------------------
!    modify the output rain and snow precip fields.
!------------------------------------------------------------------------
              Precip_state%surfrain(i,j) =    &
                                   scalef(i,j)*Precip_state%surfrain(i,j)
              Precip_state%surfsnow(i,j) =    &
                                   scalef(i,j)*Precip_state%surfsnow(i,j)
            end if
          end do
        end do
      end if 

!------------------------------------------------------------------------
!    save the rain and snow precipitation fields before any lower limit
!    is imposed (usually 0.0).
!------------------------------------------------------------------------
      if (Lsdiag_mp_control%diag_id%neg_rain > 0) &
        Lsdiag_mp%diag_4d(:,:,1,    &
                      Lsdiag_mp_control%diag_pt%neg_rain) = 1.0e3*    &
          (Precip_state%surfrain(:,:) - Precip_state%surfsnow(:,:))* &
                                                                   dtcloud
      if (Lsdiag_mp_control%diag_id%neg_snow > 0) &
        Lsdiag_mp%diag_4d(:,:,1,    &
                        Lsdiag_mp_control%diag_pt%neg_snow) = 1.0e3*    &
          (Precip_state%surfsnow(:,:))*dtcloud 

!-----------------------------------------------------------------------
!    impose lower limit.
!-----------------------------------------------------------------------
      Precip_state%surfrain = max(     &
             1.e3*(Precip_state%surfrain - Precip_state%surfsnow)*   &
                                        dtcloud , lowest_allowed_precip) 
      Precip_state%surfsnow = max(    &
             1.e3*Precip_state%surfsnow*dtcloud,   &
                                                  lowest_allowed_precip)

!-----------------------------------------------------------------------
!    compute amount of precip which has been eliminated by this
!    adjustment.
!-----------------------------------------------------------------------
      if (Lsdiag_mp_control%diag_id%neg_rain > 0) &
          Lsdiag_mp%diag_4d(:,:,1,Lsdiag_mp_control%diag_pt%neg_rain) =   &
          -1.0*( (Precip_state%surfrain(:,:))  -   &
                  Lsdiag_mp%diag_4d(:,:,1,   &
                              Lsdiag_mp_control%diag_pt%neg_rain))/dtcloud 
      if (Lsdiag_mp_control%diag_id%neg_snow > 0) &
          Lsdiag_mp%diag_4d(:,:,1,Lsdiag_mp_control%diag_pt%neg_snow) =   &
          -1.0*( (Precip_state%surfsnow(:,:))  -   &
                  Lsdiag_mp%diag_4d(:,:,1,    &
                             Lsdiag_mp_control%diag_pt%neg_snow))/dtcloud

!-------------------------------------------------------------------------

end subroutine adjust_precip_fields



!#########################################################################

subroutine adjust_for_supersaturation_removal (  &
                      ix, jx, kx, C2ls_mp, Input_mp, Atmos_state, &
                       ssat_disposal, Particles, Cloud_state, Lsdiag_mp, &
                                                       lsdiag_mp_control ) 

!-----------------------------------------------------------------------
!    with tiedtke macrophysics, supersaturation removal results in an 
!    increase in cloudiness to the max allowable cloudiness in the grid 
!    box and a consequent increase in activated aerosols due to this 
!    increase in coverage when the Ming dqa activation is being used.
!-----------------------------------------------------------------------

integer,                    intent(in)    :: ix, jx, kx
type(mp_lsdiag_type),       intent(inout) :: Lsdiag_mp
type(mp_lsdiag_control_type), intent(inout) :: Lsdiag_mp_control
type(atmos_state_type),     intent(inout) :: Atmos_state
type(mp_input_type),        intent(inout) :: Input_mp    
type(mp_conv2ls_type),      intent(inout) :: C2ls_mp     
type(cloud_state_type),     intent(inout) :: Cloud_state
type(particles_type),       intent(inout) :: Particles
real, dimension(:,:,:),     intent(in)    :: ssat_disposal


!-----------------------------------------------------------------------
!   local variables:

      real, dimension(ix, jx, kx)  :: rho, tmp2s
      integer                      :: i, j, k

!----------------------------------------------------------------------
!    process only occurs if using tiedtke macrophysics and not doing
!    pdf clouds.
!----------------------------------------------------------------------
      if (tiedtke_macrophysics .and.  .not. do_pdf_clouds) then

!-----------------------------------------------------------------------
!    where supersaturation is present, define the effects of removing it 
!    on the cloud area and cloud particle / ice crystal number.
!-----------------------------------------------------------------------
        do k=1,kx
          do j=1,jx
            do i=1,ix
              if (ssat_disposal(i,j,k) > 0.0) then
 
!-----------------------------------------------------------------------
!    define the density (rho).
!-----------------------------------------------------------------------
                rho(i,j,k) = Input_mp%pfull(i,j,k)/   &
                                        (RDGAS*Atmos_state%tn(i,j,k))

!-----------------------------------------------------------------------
!    define the area unavailable for large-scale clouds due to it 
!    containing convective cloud (tmp2s).
!-----------------------------------------------------------------------
                if (limit_conv_cloud_frac) then
                  tmp2s(i,j,k) = C2ls_mp%convective_humidity_area(i,j,k)
                else
                  tmp2s(i,j,k) = 0.
                endif

!-----------------------------------------------------------------------
!    when dqa activation is being used, the increase in cloud area results
!    in an increase in activated ice particles and cloud nuclei,
!    proportional to the cloud area increase. save the incremental 
!    increase due to removing superstauration as diagnostics.
!-----------------------------------------------------------------------
                if (dqa_activation) then
                  if (ssat_disposal(i,j,k) == 2.) then
                    Cloud_state%SNi_out(i,j,k) =     &
                        Cloud_state%SNi_out(i,j,k) + &
                          (Particles%crystal1(i,j,k)/rho(i,j,k)*  &
                               (1. - Cloud_state%qa_upd(i,j,k) -  &
                                   tmp2s(i,j,k))/dtcloud)*dtcloud
                    if (Lsdiag_mp_control%diag_id%qnidt_super +   &
                            Lsdiag_mp_control%diag_id%qni_super_col > 0 ) &
                        Lsdiag_mp%diag_4d(i,j,k,  &
                              Lsdiag_mp_control%diag_pt%qnidt_super) =    &
                        Particles%crystal1(i,j,k)/rho(i,j,k)*  &
                     (1. - Cloud_state%qa_upd(i,j,k) - tmp2s(i,j,k))/  &
                                                                  dtcloud  

                  else if (ssat_disposal(i,j,k) == 1.) then
                    Cloud_state%SN_out(i,j,k) =    &
                        Cloud_state%SN_out(i,j,k) +  &
                          (Particles%drop2(i,j,k)*    &
                               (1. - Cloud_state%qa_upd(i,j,k) -   &
                                  tmp2s(i,j,k))/dtcloud)*dtcloud
                    if (Lsdiag_mp_control%diag_id%qndt_super +   &
                            Lsdiag_mp_control%diag_id%qn_super_col > 0 ) &
                         Lsdiag_mp%diag_4d(i,j,k,   &
                             Lsdiag_mp_control%diag_pt%qndt_super ) =    &
                           Particles%drop2(i,j,k)*               &
                     (1. - Cloud_state%qa_upd(i,j,k) - tmp2s(i,j,k))/  &
                                                                  dtcloud  
                  endif
                end if ! dqa_activation
                if (max(Lsdiag_mp_control%diag_id%qadt_super,  &
                            Lsdiag_mp_control%diag_id%qa_super_col) > 0) then
                  Lsdiag_mp%diag_4d(i,j,k,  &
                        Lsdiag_mp_control%diag_pt%qadt_super ) = &              
                     (1. - Cloud_state%qa_upd(i,j,k) - tmp2s(i,j,k))/  &
                                                                  dtcloud  
                endif

!-------------------------------------------------------------------------
!    add the change to the cloud area increment resulting from this 
!    process (SA_out), and update the model cloud area after this process 
!    is completed (Cloud_state%qa_upd).
!-------------------------------------------------------------------------
                Cloud_state%SA_out(i,j,k) =   &
                      Cloud_state%SA_out(i,j,k) + &
                          (1. - Cloud_state%qa_upd(i,j,k) - tmp2s(i,j,k))  
                Cloud_state%qa_upd(i,j,k) = 1. - tmp2s(i,j,k)     
              endif ! ssat_disposal > 0.0
            end do
          end do
        end do
      end if 

!-------------------------------------------------------------------------


end subroutine adjust_for_supersaturation_removal


!######################################################################

subroutine destroy_tiny_clouds (   &
                  ix, jx, kx, Cloud_state, Tend_mp, Lsdiag_mp,   &
                        Lsdiag_mp_control, C2ls_mp, Input_mp, Atmos_state)

!-----------------------------------------------------------------------
!    routine to conservatively remove unacceptably small clouds.
!-----------------------------------------------------------------------

integer,                    intent(in)    :: ix, jx, kx
type(mp_lsdiag_type),       intent(inout) :: Lsdiag_mp
type(mp_lsdiag_control_type), intent(inout) :: Lsdiag_mp_control
type(mp_tendency_type),     intent(inout) :: Tend_mp
type(atmos_state_type),     intent(inout) :: Atmos_state
type(mp_input_type),        intent(inout) :: Input_mp    
type(mp_conv2ls_type),      intent(inout) :: C2ls_mp    
type(cloud_state_type),     intent(inout) :: Cloud_state

!-----------------------------------------------------------------------
!   local variables:

      real, dimension (ix,jx,kx)   :: ql_new, qi_new, qn_new, qni_new, &
                                      qa_new
      integer :: i,j,k


!-----------------------------------------------------------------------
!    define current cloud and particle values.
!----------------------------------------------------------------------
      ql_new  = Cloud_state%ql_in  + Cloud_state%SL_out 
      qi_new  = Cloud_state%qi_in  + Cloud_state%SI_out
      qn_new  = Cloud_state%qn_in  + Cloud_state%SN_out         
      qni_new = Cloud_state%qni_in + Cloud_state%SNi_out         

!-----------------------------------------------------------------------
!    if these values are lower than acceptable, or if the new cloud area
!    is lower than acceptable, set the tendency to balance the input value,
!    so that the field is 0. upon exiting this routine. include 
!    adjustments to temp and vapor to conserve energy and water mass.
!-----------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix
            if ((ql_new(i,j,k) <= qmin  .and.    &
                 qi_new(i,j,k) <= qmin)   .or.   &
                (Cloud_state%qa_upd(i,j,k) <=        qmin)) then
              Cloud_state%SL_out(i,j,k) = Cloud_state%SL_out(i,j,k) -  &
                                                              ql_new(i,j,k)
              Cloud_state%SI_out(i,j,k) = Cloud_state%SI_out(i,j,k) -  &
                                                              qi_new(i,j,k)
              Cloud_state%SA_out(i,j,k) = Cloud_state%SA_out(i,j,k) -  &
                                                  Cloud_state%qa_upd(i,j,k)
              Tend_mp%ttnd(i,j,k) = Tend_mp%ttnd(i,j,k) -    &
                            (hlv*ql_new(i,j,k) + hls*qi_new(i,j,k))/cp_air
              Tend_mp%qtnd(i,j,k) = Tend_mp%qtnd(i,j,k) +    &
                                           (ql_new(i,j,k) + qi_new(i,j,k))
              Cloud_state%SN_out(i,j,k) = Cloud_state%SN_out(i,j,k) -  &
                                                             qn_new(i,j,k)
              Cloud_state%SNi_out(i,j,k) = Cloud_state%SNi_out(i,j,k) -  &
                                                            qni_new(i,j,k)

!------------------------------------------------------------------------
!    save diagnostics defining the adjustments made here to destroy the
!    clouds.
!------------------------------------------------------------------------
              if (Lsdiag_mp_control%diag_id%qldt_destr > 0 .or.   &
                             Lsdiag_mp_control%diag_id%ql_destr_col > 0) &
                  Lsdiag_mp%diag_4d(i,j,k,  &
                             Lsdiag_mp_control%diag_pt%qldt_destr) =  &
                                      - ql_new(i,j,k)/dtcloud
              if (Lsdiag_mp_control%diag_id%qidt_destr > 0 .or.   &
                             Lsdiag_mp_control%diag_id%qi_destr_col > 0) &
                Lsdiag_mp%diag_4d(i,j,k,     &
                                Lsdiag_mp_control%diag_pt%qidt_destr) =  &
                                      - qi_new(i,j,k)/dtcloud
              if (Lsdiag_mp_control%diag_id%qadt_destr > 0 .or.   &
                             Lsdiag_mp_control%diag_id%qa_destr_col > 0) &
                   Lsdiag_mp%diag_4d(i,j,k,    &
                                Lsdiag_mp_control%diag_pt%qadt_destr) =  &
                          - Cloud_state%qa_upd(i,j,k)/dtcloud
              if (Lsdiag_mp_control%diag_id%qndt_destr > 0 .or.   &
                             Lsdiag_mp_control%diag_id%qn_destr_col > 0) &
                Lsdiag_mp%diag_4d(i,j,k,   &
                              Lsdiag_mp_control%diag_pt%qndt_destr) =  &
                                      - qn_new(i,j,k)/dtcloud
              if (Lsdiag_mp_control%diag_id%qnidt_destr +    &
                            Lsdiag_mp_control%diag_id%qni_destr_col > 0) &
                Lsdiag_mp%diag_4d(i,j,k,   &
                            Lsdiag_mp_control%diag_pt%qnidt_destr) = &
                                    - qni_new(i,j,k)/dtcloud
              if (Lsdiag_mp_control%diag_id%qdt_destr +    &
                              Lsdiag_mp_control%diag_id%q_destr_col > 0) &
                Lsdiag_mp%diag_4d(i,j,k,    &
                           Lsdiag_mp_control%diag_pt%qdt_destr) =   &
                      (ql_new(i,j,k) + qi_new(i,j,k))/dtcloud
            endif
          end do
        end do
      end do

!-----------------------------------------------------------------------
!    redefine the new cloud tracer values.
!-----------------------------------------------------------------------
      ql_new  =  Cloud_state%ql_in  + Cloud_state%SL_out 
      qi_new  =  Cloud_state%qi_in  + Cloud_state%SI_out
      qn_new  =  Cloud_state%qn_in  + Cloud_state%SN_out
      qni_new =  Cloud_state%qni_in + Cloud_state%SNI_out

!-----------------------------------------------------------------------
!    if the new value of cloud water is too small (including negative 
!    roundoff values), and the vapor will remain positive when 
!    conservatively adjusted, eliminate the cloudwater by adjusting the
!    vapor.  
!-----------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix
            if (abs(ql_new(i,j,k)) .le. qmin  .and.  &
                 Input_mp%qin(i,j,k) + Tend_mp%qtnd(i,j,k) +   &
                                                 ql_new(i,j,k) > 0.0) then
              Cloud_state%SL_out(i,j,k) =  - Cloud_state%ql_in(i,j,k)
              Tend_mp%qtnd(i,j,k) = Tend_mp%qtnd(i,j,k) + ql_new(i,j,k)
              Tend_mp%ttnd(i,j,k) = Tend_mp%ttnd(i,j,k) -    &
                                               (hlv*ql_new(i,j,k))/cp_air

!------------------------------------------------------------------------
!    compute diagnostic for this liquid loss due to this cleanup.
!------------------------------------------------------------------------
              if (Lsdiag_mp_control%diag_id%qdt_cleanup_liquid +    &
                    Lsdiag_mp_control%diag_id%q_cleanup_liquid_col > 0) &
                 Lsdiag_mp%diag_4d(i,j,k,   &
                       Lsdiag_mp_control%diag_pt%qdt_cleanup_liquid) =   &
                                                   ql_new(i,j,k)/dtcloud

!------------------------------------------------------------------------
!    with the removal of all liquid, the cloud droplet number must also
!    be set to 0.0. define diagnostic for droplet loss due to this cleanup.
!------------------------------------------------------------------------
              Cloud_state%SN_out(i,j,k) = Cloud_state%SN_out(i,j,k) -   &
                                                           qn_new(i,j,k) 
              if (Lsdiag_mp_control%diag_id%qndt_cleanup +   &
                          Lsdiag_mp_control%diag_id%qn_cleanup_col > 0) &
                Lsdiag_mp%diag_4d(i,j,k,   &
                         Lsdiag_mp_control%diag_pt%qndt_cleanup) = &
                                              - qn_new(i,j,k)/dtcloud
            endif
          end do
        end do
      end do

!-----------------------------------------------------------------------
!    if the new value of cloud ice is too small (including negative 
!    roundoff values), and the vapor will remain positive when 
!    conservatively adjusted, eliminate the cloudice by adjusting the
!    vapor.  
!-----------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix
            if (abs(qi_new(i,j,k)) .le. qmin  .and.  &
                  Input_mp%qin(i,j,k) + Tend_mp%qtnd(i,j,k) +   &
                                                qi_new(i,j,k) > 0.0) then
              Cloud_state%SI_out(i,j,k) =  - Cloud_state%qi_in(i,j,k)
              Tend_mp%qtnd(i,j,k) = Tend_mp%qtnd(i,j,k) + qi_new(i,j,k)
              Tend_mp%ttnd(i,j,k) = Tend_mp%ttnd(i,j,k) -    &
                                                (hls*qi_new(i,j,k))/cp_air

!------------------------------------------------------------------------
!    compute diagnostic for this liquid loss due to this cleanup.
!------------------------------------------------------------------------
              if (Lsdiag_mp_control%diag_id%qdt_cleanup_ice +    &
                        Lsdiag_mp_control%diag_id%q_cleanup_ice_col > 0) &
                Lsdiag_mp%diag_4d(i,j,k,   &
                         Lsdiag_mp_control%diag_pt%qdt_cleanup_ice) =   &
                                       qi_new(i,j,k)/dtcloud

!------------------------------------------------------------------------
!    with the removal of all ice, the ice particle number must also
!    be set to 0.0. define diagnostic for crystal loss due to this cleanup.
!------------------------------------------------------------------------
              Cloud_state%SNI_out(i,j,k) = Cloud_state%SNI_out(i,j,k) -  &
                                                             qni_new(i,j,k) 
              if (Lsdiag_mp_control%diag_id%qnidt_cleanup +    &
                    Lsdiag_mp_control%diag_id%qni_cleanup_col > 0) &
               Lsdiag_mp%diag_4d(i,j,k,     &
                         Lsdiag_mp_control%diag_pt%qnidt_cleanup) = &
                                  - qni_new(i,j,k)/dtcloud
            endif
          end do
        end do
      end do
    
!-----------------------------------------------------------------------
!    force the change in ice crystal number to not be so large as to 
!    eliminate more crystals than were present initially. save a diagnostic
!    if desired.
!-----------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix
            if (Lsdiag_mp_control%diag_id%qnidt_cleanup2 +    &
                   Lsdiag_mp_control%diag_id%qni_cleanup2_col > 0) &
              Lsdiag_mp%diag_4d(i,j,k,    &
                              Lsdiag_mp_control%diag_pt%qnidt_cleanup2) = &
                                           Cloud_state%SNi_out(i,j,k)

            Cloud_state%SNi_out(i,j,k) = MAX(Cloud_state%SNi_out(i,j,k), &
                                            - Cloud_state%qni_in(i,j,k))
            if (Lsdiag_mp_control%diag_id%qnidt_cleanup2 +    &
                        Lsdiag_mp_control%diag_id%qni_cleanup2_col > 0) &
                  Lsdiag_mp%diag_4d(i,j,k,   &
                          Lsdiag_mp_control%diag_pt%qnidt_cleanup2) =    &
            (Lsdiag_mp%diag_4d(i,j,k,    &
                    Lsdiag_mp_control%diag_pt%qnidt_cleanup2) - &
                     Cloud_state%SNi_out(i,j,k))*inv_dtcloud
          end do
        end do
      end do

    
!-----------------------------------------------------------------------
!    force the change in cloud droplet number to not be so large as to 
!    eliminate more droplets than were present initially. save a diagnostic
!    if desired.
!-----------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix
            if (Lsdiag_mp_control%diag_id%qndt_cleanup2 +     &
                     Lsdiag_mp_control%diag_id%qn_cleanup2_col > 0) &
                   Lsdiag_mp%diag_4d(i,j,k,    &
                          Lsdiag_mp_control%diag_pt%qndt_cleanup2) = &
                                            Cloud_state%SN_out(i,j,k)
            Cloud_state%SN_out(i,j,k) = MAX(Cloud_state%SN_out(i,j,k),  &
                                              - Cloud_state%qn_in(i,j,k))
            if (Lsdiag_mp_control%diag_id%qndt_cleanup2 +    &
                          Lsdiag_mp_control%diag_id%qn_cleanup2_col > 0) &
               Lsdiag_mp%diag_4d(i,j,k,    &
                           Lsdiag_mp_control%diag_pt%qndt_cleanup2) = &
                 (Lsdiag_mp%diag_4d(i,j,k,   &
                            Lsdiag_mp_control%diag_pt%qndt_cleanup2) -   &
                     Cloud_state%SN_out(i,j,k))*inv_dtcloud
          end do
        end do
      end do

!----------------------------------------------------------------------
!    make sure the new cloud area is not smaller than the minimum 
!    allowable. if not set the tendency so that cloud area is reduced to 
!    zero after the step. save a diagnostic if desired.
!----------------------------------------------------------------------
      if (Lsdiag_mp_control%diag_id%qadt_destr +    &
                          Lsdiag_mp_control%diag_id%qa_destr_col > 0)    &
      Lsdiag_mp%diag_4d(:,:,:,Lsdiag_mp_control%diag_pt%qadt_destr) =    &
         Lsdiag_mp%diag_4d(:,:,:,Lsdiag_mp_control%diag_pt%qadt_destr) +  &
                              Cloud_state%SA_out*inv_dtcloud

      qa_new = Cloud_state%qa_in + Cloud_state%SA_out

      where ( abs(qa_new) .le. qmin )
        Cloud_state%SA_out  = -Cloud_state%qa_in
      endwhere

      if (Lsdiag_mp_control%diag_id%qadt_destr +   &
                          Lsdiag_mp_control%diag_id%qa_destr_col > 0)    &
      Lsdiag_mp%diag_4d(:,:,:,Lsdiag_mp_control%diag_pt%qadt_destr) =    &
      Lsdiag_mp%diag_4d(:,:,:,Lsdiag_mp_control%diag_pt%qadt_destr) -     &
                              Cloud_state%SA_out*inv_dtcloud

!------------------------------------------------------------------------
!    enforce constraints on the cloud area. the change must not be so
!    large as to eliminate more cloud than is present. also it must not
!    be so large as to more than fill the available area in the grid box
!    (some area may have been taken up by the convective system, so the
!    max available area is (1 - conv area). Include a diagnostic if
!    desired. this constraint has already been imposed with r-k 
!    microphysics, as part of the destruction diagnostic.
!------------------------------------------------------------------------
      if (do_mg_ncar_microphys .or. do_ncar_microphys .or. &
                                                  do_mg_microphys) then
        if (Lsdiag_mp_control%diag_id%qadt_limits +    &
                         Lsdiag_mp_control%diag_id%qa_limits_col > 0)    &
       Lsdiag_mp%diag_4d(:,:,:,Lsdiag_mp_control%diag_pt%qadt_limits) =   &
                                  Cloud_state%SA_out(:,:,:)

        Cloud_state%SA_out = MAX(Cloud_state%SA_out,-Cloud_state%qa_in)
        Cloud_state%SA_out = MIN(Cloud_state%SA_out,   &
                             1. - C2ls_mp%convective_humidity_area - &
                                                       Cloud_state%qa_in)
        if (Lsdiag_mp_control%diag_id%qadt_limits +   &
                    Lsdiag_mp_control%diag_id%qa_limits_col      > 0)    &
         Lsdiag_mp%diag_4d(:,:,:,Lsdiag_mp_control%diag_pt%qadt_limits) = &
                      (Cloud_state%SA_out(:,:,:) - &
        Lsdiag_mp%diag_4d(:,:,:,Lsdiag_mp_control%diag_pt%qadt_limits))* &
                                                               inv_dtcloud
      endif

!-----------------------------------------------------------------------


end subroutine destroy_tiny_clouds



!#######################################################################

subroutine destroy_tiny_clouds_clubb (    &
                   ix, jx, kx, Cloud_state, Tend_mp, Lsdiag_mp, &
                        Lsdiag_mp_control, C2ls_mp, Input_mp, Atmos_state)

integer,                    intent(in)    :: ix,jx,kx
type(atmos_state_type),     intent(inout) :: Atmos_state
type(cloud_state_type),     intent(inout) :: Cloud_state
type(mp_tendency_type),     intent(inout) :: Tend_mp
type(mp_conv2ls_type),      intent(inout) :: C2ls_mp
type(mp_input_type),        intent(in   ) :: Input_mp
type(mp_lsdiag_type),       intent(inout) :: Lsdiag_mp
type(mp_lsdiag_control_type), intent(inout) :: Lsdiag_mp_control

!----------------------------------------------------------------------
!   local variables:

      integer :: i, j, k
      real, dimension (ix,jx,kx) :: ql_new, qi_new, qn_new, &
                                    qni_new, qa_new

!-----------------------------------------------------------------------
!    define current cloud and particle values.
!----------------------------------------------------------------------
     ql_new = Input_mp%tracer(:,:,:,nql)+Tend_mp%q_tnd(:,:,:,nql)
     qi_new = Input_mp%tracer(:,:,:,nqi)+Tend_mp%q_tnd(:,:,:,nqi)
     qn_new = Input_mp%tracer(:,:,:,nqn)+Tend_mp%q_tnd(:,:,:,nqn)
     qni_new = Input_mp%tracer(:,:,:,nqni)+Tend_mp%q_tnd(:,:,:,nqni)
     qa_new = Input_mp%tracer(:,:,:,nqa) + Tend_mp%q_tnd(:,:,:,nqa)

!-----------------------------------------------------------------------
!    if these values are lower than acceptable, set the tendency to 
!    balance the input value, so that the field is 0. upon exiting this 
!    loop. include adjustments to temp and vapor to conserve energy and 
!    water mass.
!-----------------------------------------------------------------------
     do k=1,kx
       do j=1,jx
         do i=1,ix
           if ( (ql_new(i,j,k) .le. qmin) .and.    &
                                       (qi_new(i,j,k) .le. qmin) ) then
             Tend_mp%qtnd(i,j,k) = Tend_mp%qtnd(i,j,k) +   &
                                           ql_new(i,j,k) + qi_new(i,j,k) 
             Tend_mp%ttnd(i,j,k) = Tend_mp%ttnd(i,j,k) -   &
                     (hlv*(ql_new(i,j,k) ) + hls*(qi_new(i,j,k) ) )/cp_air
             Tend_mp%q_tnd(i,j,k,nql) = Tend_mp%q_tnd(i,j,k,nql) -  &
                                                          (ql_new(i,j,k))
             Tend_mp%q_tnd(i,j,k,nqi) = Tend_mp%q_tnd(i,j,k,nqi) -   &
                                                          (qi_new(i,j,k))
             Tend_mp%q_tnd(i,j,k,nqa) = Tend_mp%q_tnd(i,j,k,nqa) -   &
                                                             qa_new(i,j,k) 
             Tend_mp%q_tnd(i,j,k,nqn) = Tend_mp%q_tnd(i,j,k,nqn) -   &
                                                           (qn_new(i,j,k) )
             Tend_mp%q_tnd(i,j,k,nqni)= Tend_mp%q_tnd(i,j,k,nqni) -   &
                                                          (qni_new(i,j,k) )

!------------------------------------------------------------------------
!    save diagnostics defining the adjustments made here to destroy the
!    clouds.
!------------------------------------------------------------------------
             if ( Lsdiag_mp_control%diag_id%qldt_destr > 0  .or.   &
                           Lsdiag_mp_control%diag_id%ql_destr_col > 0 )  &
                Lsdiag_mp%diag_4d(i,j,k,   &
                        Lsdiag_mp_control%diag_pt%qldt_destr) =   &
                                  - (ql_new(i,j,k) )/dtcloud
             if ( Lsdiag_mp_control%diag_id%qidt_destr > 0  .or.   &
                          Lsdiag_mp_control%diag_id%qi_destr_col > 0   ) &
               Lsdiag_mp%diag_4d(i,j,k,  &
                             Lsdiag_mp_control%diag_pt%qidt_destr) =   &
                                   - (qi_new(i,j,k) )/dtcloud
             if ( Lsdiag_mp_control%diag_id%qadt_destr > 0  .or.   &
                          Lsdiag_mp_control%diag_id%qa_destr_col > 0   ) &
               Lsdiag_mp%diag_4d(i,j,k,   &
                              Lsdiag_mp_control%diag_pt%qadt_destr) =   &
                                     - qa_new(i,j,k)/dtcloud
             if ( Lsdiag_mp_control%diag_id%qndt_destr > 0  .or.   &
                           Lsdiag_mp_control%diag_id%qn_destr_col > 0  ) &
               Lsdiag_mp%diag_4d(i,j,k,    &
                               Lsdiag_mp_control%diag_pt%qndt_destr) =  &
                                   - (qn_new(i,j,k) )/dtcloud
             if ( Lsdiag_mp_control%diag_id%qnidt_destr > 0   )          &
               Lsdiag_mp%diag_4d(i,j,k,    &
                              Lsdiag_mp_control%diag_pt%qnidt_destr) =  &
                                  - (qni_new(i,j,k) )/dtcloud
             if ( Lsdiag_mp_control%diag_id%qdt_destr > 0 )             &
              Lsdiag_mp%diag_4d(i,j,k,   &
                            Lsdiag_mp_control%diag_pt%qdt_destr) =  &
                      (ql_new(i,j,k) + qi_new(i,j,k))/dtcloud

           endif
         end do
       end do
     end do

!-----------------------------------------------------------------------
!    redefine the new cloud tracer values.
!-----------------------------------------------------------------------
      ql_new = Input_mp%tracer(:,:,:,nql)+Tend_mp%q_tnd(:,:,:,nql)
      qi_new = Input_mp%tracer(:,:,:,nqi)+Tend_mp%q_tnd(:,:,:,nqi)
      qn_new = Input_mp%tracer(:,:,:,nqn)+Tend_mp%q_tnd(:,:,:,nqn)
      qni_new = Input_mp%tracer(:,:,:,nqni)+Tend_mp%q_tnd(:,:,:,nqni)

!-----------------------------------------------------------------------
!    if only the new value of cloud water is too small (including negative 
!    values), eliminate the cloudwater by conservatively adjusting 
!    the vapor.  
!-----------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix
            if ( (ql_new(i,j,k) .le. qmin) ) then
              Tend_mp%qtnd(i,j,k) = Tend_mp%qtnd(i,j,k) + (ql_new(i,j,k))
              Tend_mp%ttnd(i,j,k) = Tend_mp%ttnd(i,j,k) -   &
                                              (hlv*(ql_new(i,j,k)))/cp_air
              Tend_mp%q_tnd(i,j,k,nql) = Tend_mp%q_tnd(i,j,k,nql) -  &
                                                           (ql_new(i,j,k))

!------------------------------------------------------------------------
!    compute diagnostic for this liquid loss due to this cleanup.
!------------------------------------------------------------------------
              if ( Lsdiag_mp_control%diag_id%qdt_cleanup_liquid > 0 ) &
                Lsdiag_mp%diag_4d(i,j,k,    &
                      Lsdiag_mp_control%diag_pt%qdt_cleanup_liquid) =   &
                                      (ql_new(i,j,k))/dtcloud

!------------------------------------------------------------------------
!    with the removal of all liquid, the cloud droplet number must also
!    be set to 0.0. define diagnostic for droplet loss due to this cleanup.
!------------------------------------------------------------------------
              Tend_mp%q_tnd(i,j,k,nqn) = Tend_mp%q_tnd(i,j,k,nqn) -   &
                                                          (qn_new(i,j,k))
              IF ( Lsdiag_mp_control%diag_id%qndt_cleanup > 0 ) &
                Lsdiag_mp%diag_4d(i,j,k,    &
                           Lsdiag_mp_control%diag_pt%qndt_cleanup) = &
                                     -(qn_new(i,j,k))/dtcloud
            endif
          end do
        end do
      end do

!-----------------------------------------------------------------------
!    if only the new value of cloud ice is too small (including negative 
!    values), eliminate the cloud ice by conservatively adjusting the
!    vapor.  
!-----------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix
            if ( (qi_new(i,j,k) .le. qmin) ) then
              Tend_mp%qtnd(i,j,k) = Tend_mp%qtnd(i,j,k) + (qi_new(i,j,k))
              Tend_mp%ttnd(i,j,k) = Tend_mp%ttnd(i,j,k) -   &
                                              (hls*(qi_new(i,j,k)))/cp_air
              Tend_mp%q_tnd(i,j,k,nqi) = Tend_mp%q_tnd(i,j,k,nqi) -  &
                                                         (qi_new(i,j,k))

!------------------------------------------------------------------------
!    compute diagnostic for this ice loss due to this cleanup.
!------------------------------------------------------------------------
              if ( Lsdiag_mp_control%diag_id%qdt_cleanup_ice > 0 ) &
                Lsdiag_mp%diag_4d(i,j,k,   &
                           Lsdiag_mp_control%diag_pt%qdt_cleanup_ice) =  &
                                      (qi_new(i,j,k))/dtcloud

!------------------------------------------------------------------------
!    with the removal of all ice, the ice crystal number must also
!    be set to 0.0. define diagnostic for crystal loss due to this cleanup.
!------------------------------------------------------------------------
              Tend_mp%q_tnd(i,j,k,nqni) = Tend_mp%q_tnd(i,j,k,nqni) -  &
                                                         (qni_new(i,j,k))
              if ( Lsdiag_mp_control%diag_id%qnidt_cleanup > 0 ) &
                 Lsdiag_mp%diag_4d(i,j,k,   &
                             Lsdiag_mp_control%diag_pt%qnidt_cleanup) =   &
                                   - (qni_new(i,j,k))/dtcloud
            endif
          end do
        end do
      end do

!------------------------------------------------------------------------

end subroutine destroy_tiny_clouds_clubb



!########################################################################


          end module ls_cloud_microphysics_mod

