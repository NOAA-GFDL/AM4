                module radiation_driver_diag_mod

! <OVERVIEW>
!    radiation_driver_diag_mod provides diagnostics and
!    solar interpolator to the radiation driver.
! </OVERVIEW>
!  <DIAGFIELDS>
!  Diagnostic fields may be output to a netcdf file by specifying the
!  module name radiation and the desired field names (given below)
!  in file diag_table. See the documentation for diag_manager.
!  
!  Diagnostic fields for module name: radiation
!  
!     field name      field description
!     ----------      -----------------
!  
!     alb_sfc         surface albedo (percent)
!     coszen          cosine of the solar zenith angle
!  
!     tdt_sw          temperature tendency for SW radiation (deg_K/sec)
!     tdt_lw          Temperature tendency for LW radiation (deg_K/sec)
!     swdn_toa        SW flux down at TOA (watts/m2)
!     rsdt            SW flux down at TOA (watts/m2)
!     swup_toa        SW flux up at TOA (watts/m2)
!     rsut            SW flux up at TOA (watts/m2)
!     olr             outgoing longwave radiation (watts/m2)
!     rlut            outgoing longwave radiation (watts/m2)
!     swup_sfc        SW flux up at surface (watts/m2)
!     rsus            SW flux up at surface (watts/m2)
!     swdn_sfc        SW flux down at surface (watts/m2)
!     rsds            SW flux down at surface (watts/m2)
!     lwup_sfc        LW flux up at surface  (watts/m2)
!     rlus            LW flux up at surface  (watts/m2)
!     lwdn_sfc        LW flux down at surface (watts/m2)
!     rlds            LW flux down at surface (watts/m2)
!  
!  NOTE: When namelist variable do_clear_sky_pass = .true. an additional clear sky
!        diagnostic fields may be saved.
!  
!     tdt_sw_clr      clear sky temperature tendency for SW radiation (deg_K/sec)
!     tdt_lw_clr      clear sky Temperature tendency for LW radiation (deg_K/sec)
!     swdn_toa_clr    clear sky SW flux down at TOA (watts/m2)
!     swup_toa_clr    clear sky SW flux up at TOA (watts/m2)
!     rsutcs          clear sky SW flux up at TOA (watts/m2)
!     olr_clr         clear sky outgoing longwave radiation (watts/m2)
!     rlutcs          clear sky outgoing longwave radiation (watts/m2)
!     swup_sfc_clr    clear sky SW flux up at surface (watts/m2)
!     rsuscs          clear sky SW flux up at surface (watts/m2)
!     swdn_sfc_clr    clear sky SW flux down at surface (watts/m2)
!     rsdscs          clear sky SW flux down at surface (watts/m2)
!     lwup_sfc_clr    clear sky LW flux up at surface  (watts/m2)
!     lwdn_sfc_clr    clear sky LW flux down at surface (watts/m2)
!     rldscs          clear sky LW flux down at surface (watts/m2)
!  </DIAGFIELDS>


use mpp_mod,               only: input_nml_file
use fms_mod,               only: fms_init, &
                                 mpp_pe, mpp_root_pe, &
                                 open_namelist_file, stdlog, stdout, &
                                 file_exist, FATAL, WARNING, NOTE, &
                                 close_file, &
                                 write_version_number, check_nml_error,&
                                 error_mesg
use fms_io_mod,            only: restore_state, &
                                 register_restart_field, restart_file_type, &
                                 save_restart
use diag_manager_mod,      only: register_diag_field, send_data, &
                                 diag_manager_init, get_base_time
use diag_data_mod,         only: CMOR_MISSING_VALUE
use time_manager_mod,      only: time_manager_init, time_type, operator(>)
use constants_mod,         only: constants_init, STEFAN, SECONDS_PER_DAY, &
                                 CP_AIR, RADIAN, WTMCO2, WTMAIR, GRAV

use atmos_cmip_diag_mod,   only: register_cmip_diag_field_3d, &
                                 register_cmip_diag_field_2d, &
                                 send_cmip_data_3d, &
                                 cmip_diag_id_type, &
                                 query_cmip_diag_id
use atmos_global_diag_mod, only: register_global_diag_field, &
                                 buffer_global_diag, &
                                 send_global_diag

! atmos physics modules

use diag_integral_mod,     only: diag_integral_init, &
                                 diag_integral_field_init, &
                                 sum_diag_integral_field

! radiation modules

use radiation_driver_types_mod, only: radiation_control_type, &
                                      astronomy_type, &
                                      rad_output_type

use aerosolrad_types_mod,  only: aerosolrad_control_type

use shortwave_types_mod,   only: sw_output_type
use longwave_types_mod,    only: lw_output_type
use radiative_gases_types_mod, only: radiative_gases_type

use shortwave_driver_mod, only: get_solar_constant

!--------------------------------------------------------------------

implicit none 
private 

!------------ version number for this module --------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

!----------------------------------------------------------------------
!    radiation_driver_diag_mod is an interface between
!    the radiation_driver_mod and diag_manager_mod
!----------------------------------------------------------------------

public  radiation_driver_diag_init, &
        update_rad_fields, &
        produce_radiation_diagnostics, &
        radiation_driver_diag_end, &
        radiation_driver_diag_endts, &
        write_solar_interp_restart_nc

private initialize_diagnostic_integrals,   &
        diag_field_init, &
        solar_flux_save_init, &
        solar_interp_register_restart

!-----------------------------------------------------------------------
!------- namelist ---------

logical ::     &    
        all_step_diagnostics = .false.!  are lw and sw radiative bdy
                                      !  fluxes and atmospheric heating 
                                      !  rates to be output on physics 
                                      !  steps ?

real    :: trop_ht_at_poles = 30000.  !  assumed height of tropoause at
                                      !  poles for case of tropause
                                      !  linearly varying with latitude
                                      !  [ Pa ]
real    :: trop_ht_at_eq    = 10000.  !  assumed height of tropoause at
                                      !  equator for case of tropause
                                      !  linearly varying with latitude
                                      !  [ Pa ]
real    :: trop_ht_constant = 20000.  !  assumed height of tropoause   
                                      !  when assumed constant       
                                      !  [ Pa ]
logical :: constant_tropo = .true.    !  generate tropopause fluxes when
                                      !  tropopause ht assumed constant?
logical :: linear_tropo   = .true.    !  generate tropopause fluxes when
                                      !  tropopause assumed to vary
                                      !  linearly with latitude?
logical :: thermo_tropo   = .false.   !  generate tropopause fluxes when
                                      !  tropopause determined thermo-
                                      !  dynamically ?

namelist /radiation_driver_diag_nml/ all_step_diagnostics, &
                                     trop_ht_at_poles, trop_ht_at_eq, &
                                     trop_ht_constant, constant_tropo, &
                                     linear_tropo, thermo_tropo


!-----------------------------------------------------------------------
type sw_flux_save_type
   real, dimension(:,:),   pointer :: flux_sw_down_total_dir_clr=>NULL(), &
                                      flux_sw_down_total_dif_clr=>NULL(), &
                                      flux_sw_down_vis_clr=>NULL()
   real, dimension(:,:),   pointer :: flux_sw_surf=>NULL(), &
                                      flux_sw_surf_dir=>NULL(), &
                                      flux_sw_surf_refl_dir=>NULL(), &
                                      flux_sw_surf_dif=>NULL(), &
                                      flux_sw_down_vis_dir=>NULL(), &
                                      flux_sw_down_vis_dif=>NULL(), &
                                      flux_sw_down_total_dir=>NULL(), &
                                      flux_sw_down_total_dif=>NULL(), &
                                      flux_sw_vis=>NULL(), &
                                      flux_sw_vis_dir=>NULL(), &
                                      flux_sw_refl_vis_dir=>NULL(), &
                                      flux_sw_vis_dif=>NULL()
   real, dimension(:,:,:), pointer :: sw_heating=>NULL(),    &
                                      tot_heating=>NULL(), &
                                      dfsw=>NULL(), ufsw=>NULL(), &
                                      fsw=>NULL(), hsw=>NULL()
   real, dimension(:,:,:), pointer :: sw_heating_clr=>NULL(), &
                                      tot_heating_clr=>NULL(), &
                                      dfswcf=>NULL(), ufswcf=>NULL(), &
                                      fswcf=>NULL(), hswcf=>NULL()
   real, dimension(:,:,:), pointer :: dfsw_ad=>NULL(), ufsw_ad=>NULL()
   real, dimension(:,:,:), pointer :: dfswcf_ad=>NULL(), ufswcf_ad=>NULL()
end type sw_flux_save_type

type lw_flux_save_type
   real, dimension(:,:,:),   pointer :: tdtlw=>NULL(), tdtlw_clr=>NULL()
   real, dimension(:,:,:),   pointer :: flxnet=>NULL(), flxnetcf=>NULL()
   real, dimension(:,:),     pointer :: olr=>NULL(), lwups=>NULL(), &
                                        lwdns=>NULL(), olr_clr=>NULL(), &
                                        lwups_clr=>NULL(), lwdns_clr=>NULL()
   real, dimension(:,:),     pointer :: olr_ad=>NULL(), lwups_ad=>NULL(), &
                                        lwdns_ad=>NULL(), olr_ad_clr=>NULL(), &
                                        lwups_ad_clr=>NULL(), lwdns_ad_clr=>NULL()
end type lw_flux_save_type

type diag_special_type
  real, dimension(:,:,:), pointer :: swdn_trop, swup_trop
  real, dimension(:,:,:), pointer :: netlw_trop
  real, dimension(:,:,:), pointer :: swdn_trop_clr, swup_trop_clr
  real, dimension(:,:,:), pointer :: netlw_trop_clr
end type diag_special_type

!-----------------------------------------------------------------------
!---- private data ----
!-- for netcdf restart
type(restart_file_type), pointer, save :: Solar_restart => NULL()
type(restart_file_type), pointer, save :: Tile_restart => NULL()
logical :: doing_netcdf_restart = .false.

!    solar_save is used when renormalize_sw_fluxes is active, to save
!    the solar factor (fracday*cosz/r**2) from the previous radiation
!    step so that the radiative forcing terms may be adjusted on each
!    timestep to reflect the current solar forcing.
!
!    Sw_flux_save contains radiative forcing terms from the previous
!    radiation time step which must be saved when renormalization is
!    activated.
!
!    Lw_flux_save contains longwave radiative forcing terms from the
!    previous radiation time step so that they may be output in the
!    diagnostics file on every physics step, if desired, so that when 
!    renormalize_sw_fluxes is active, total radiative terms may be easily
!    generated.

real, allocatable, dimension(:,:)   ::  solar_save
type(sw_flux_save_type), save :: Sw_flux_save
type(lw_flux_save_type), save :: Lw_flux_save
type(diag_special_type), save :: Diag_special
!-----------------------------------------------------------------------
!    diagnostics variables
integer, parameter :: MX_SPEC_LEVS = 4 
                             ! number of special levels at
                             ! which radiative fluxes are to be 
                             ! calculated for diagnostic purposes

character(len=16)            :: mod_name = 'radiation'
integer                      :: id_alb_sfc, id_cosz, id_fracday, &
                                id_alb_sfc_avg, &
                                id_alb_sfc_vis_dir, id_alb_sfc_nir_dir,&
                                id_alb_sfc_vis_dif, id_alb_sfc_nir_dif
integer                      :: id_flux_sw_dir, id_flux_sw_dif, &
                                id_flux_sw_refl_dir,  &
                                id_flux_sw_refl_vis_dir, id_flux_sw, &
                                id_flux_sw_down_vis_dir, &
                                id_flux_sw_down_vis_dif, &
                                id_flux_sw_down_total_dir, &
                                id_flux_sw_down_total_dif, &
                                id_flux_sw_down_total_dir_clr, &
                                id_flux_sw_down_total_dif_clr, &
                                id_flux_sw_down_vis_clr, &
                                id_flux_sw_vis, &
                                id_flux_sw_vis_dir, &
                                id_flux_sw_vis_dif, &
                                id_co2mass, id_ch4global, id_n2oglobal, &
                                id_cfc11global, id_cfc12global, &
                                id_cfc113global, id_hcfc22global, &
                                id_rrvco2, id_rrvf11, id_rrvf12, &
                                id_rrvf113, id_rrvf22, id_rrvch4, &
                                id_rrvn2o, id_co2_tf, id_ch4_tf, &
                                id_n2o_tf, id_solar_constant
integer                      :: id_conc_drop, id_conc_ice

integer                      :: id_allradp, id_heat2d, id_heat2d_sw
integer, dimension(2)        :: id_tdt_sw,   id_tdt_lw,  &
                                id_ufsw, id_dfsw,  &
                                id_flxnet, &
                                id_swdn_toa, id_swup_toa, id_olr, &
                                id_netrad_toa,  id_netrad_1_Pa,  &
                                id_swup_sfc, id_swdn_sfc,         &
                                id_lwup_sfc, id_lwdn_sfc
integer                      :: id_rlds, id_rldscs, id_rlus, id_rsds,   &
                                id_rsdscs, id_rsus, id_rsuscs, id_rsdt, &
                                id_rsut, id_rsutcs, id_rlut, id_rlutcs, &
                                id_rtmt, id_rsdsdiff, id_rsdscsdiff
integer                      :: id_rsdsaf, id_rsusaf, id_rsutaf, &
                                id_rsdscsaf, id_rsuscsaf, id_rsutcsaf, &
                                id_rldsaf, id_rlutaf, id_rldscsaf, id_rlutcsaf
type(cmip_diag_id_type)      :: ID_tntr, ID_tntrs, ID_tntrscs, ID_tntrl, ID_tntrlcs, &
                                ID_rsu, ID_rsucs, ID_rsd, ID_rsdcs, &
                                ID_rsuaf, ID_rsucsaf, ID_rsdaf, ID_rsdcsaf
integer, dimension(MX_SPEC_LEVS,2)   :: id_swdn_special,   &
                                        id_swup_special,  &
                                        id_netlw_special
integer, dimension(2)        :: id_swtoa, id_swsfc,               &
                                id_lwsfc,                         &
                                id_swtoa_ad, id_swsfc_ad,         &
                                id_swdn_sfc_ad,                   &
                                id_swup_sfc_ad,                   &
                                id_swup_toa_ad,                   &
                                id_olr_ad, id_lwsfc_ad

! globally averaged diagnostics
integer :: id_rlut_g, id_rlutcs_g, id_rsut_g, id_rsutcs_g, id_rsdt_g
integer :: id_rss_g

real, dimension(:,:), allocatable :: swups_acc, swdns_acc
real, dimension(:,:), allocatable :: olr_intgl, swabs_intgl

real                         :: missing_value = -999.
character(len=8)             :: std_digits   = 'f8.3'
character(len=8)             :: extra_digits = 'f16.11'

logical  ::  do_swaerosol_forcing
logical  ::  do_lwaerosol_forcing
integer  ::  indx_swaf
integer  ::  indx_lwaf

!-----------------------------------------------------------------------

logical :: module_is_initialized = .false.

!-----------------------------------------------------------------------

                         contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine radiation_driver_diag_init (Time, id, jd, kmax, axes, &
                                       Rad_control, Aerosolrad_control)

!--------------------------------------------------------------------
type(time_type),               intent(in) :: Time
integer,                       intent(in) :: id, jd, kmax
integer, dimension(4),         intent(in) :: axes
type(radiation_control_type),  intent(in) :: Rad_control
type(aerosolrad_control_type), intent(in) :: Aerosolrad_control
!---------------------------------------------------------------------
!   local variables
      integer           ::   unit, io, ierr, logunit

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized. note that data_override_init cannot
!    be called successfully from here (a data_override_mod feature);
!    instead it relies upon a check for previous initialization when
!    subroutine data_override is called.
!---------------------------------------------------------------------
      call fms_init
      call constants_init
      call diag_manager_init
      call time_manager_init
      call diag_integral_init

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=radiation_driver_diag_nml, iostat=io)
      ierr = check_nml_error(io,'radiation_driver_diag_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=radiation_driver_diag_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'radiation_driver_diag_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
       if (mpp_pe() == mpp_root_pe() ) &
            write (logunit, nml=radiation_driver_diag_nml)

!---------------------------------------------------------------------
!    save aerosol forcing flags as module variables
!---------------------------------------------------------------------
      do_swaerosol_forcing = Aerosolrad_control%do_swaerosol_forcing
      do_lwaerosol_forcing = Aerosolrad_control%do_lwaerosol_forcing

!    indexing for aerosol forcing output
      indx_swaf = Aerosolrad_control%indx_swaf
      indx_lwaf = Aerosolrad_control%indx_lwaf

!---------------------------------------------------------------------
!    allocate space for variables which must be saved when sw fluxes
!    are renormalized or diagnostics are desired to be output on every
!    physics step.
!---------------------------------------------------------------------
      if (Rad_control%renormalize_sw_fluxes .or. all_step_diagnostics) then
          allocate (solar_save                       (id,jd))
          allocate (Sw_flux_save%flux_sw_surf          (id,jd))
          allocate (Sw_flux_save%flux_sw_surf_dir      (id,jd))
          allocate (Sw_flux_save%flux_sw_surf_refl_dir (id,jd))
          allocate (Sw_flux_save%flux_sw_surf_dif      (id,jd))
          allocate (Sw_flux_save%flux_sw_down_vis_dir  (id,jd))
          allocate (Sw_flux_save%flux_sw_down_vis_dif  (id,jd))
          allocate (Sw_flux_save%flux_sw_down_total_dir(id,jd))
          allocate (Sw_flux_save%flux_sw_down_total_dif(id,jd))
          allocate (Sw_flux_save%flux_sw_vis           (id,jd))
          allocate (Sw_flux_save%flux_sw_vis_dir       (id,jd))
          allocate (Sw_flux_save%flux_sw_refl_vis_dir  (id,jd))
          allocate (Sw_flux_save%flux_sw_vis_dif       (id,jd))
          allocate (Sw_flux_save%sw_heating            (id,jd,kmax))
          allocate (Sw_flux_save%tot_heating           (id,jd,kmax))
          allocate (Sw_flux_save%dfsw                  (id,jd,kmax+1))
          allocate (Sw_flux_save%ufsw                  (id,jd,kmax+1))
          allocate (Sw_flux_save% fsw                  (id,jd,kmax+1))
          allocate (Sw_flux_save% hsw                  (id,jd,kmax))
          if (do_swaerosol_forcing) then
              allocate (Sw_flux_save%dfsw_ad             (id,jd,kmax+1))
              allocate (Sw_flux_save%ufsw_ad             (id,jd,kmax+1))
          endif
          if (Rad_control%do_totcld_forcing) then 
              allocate (Sw_flux_save%sw_heating_clr            (id,jd,kmax))
              allocate (Sw_flux_save%tot_heating_clr           (id,jd,kmax))
              allocate (Sw_flux_save%dfswcf                    (id,jd,kmax+1))
              allocate (Sw_flux_save%ufswcf                    (id,jd,kmax+1))
              allocate (Sw_flux_save% fswcf                    (id,jd,kmax+1))
              allocate (Sw_flux_save% hswcf                    (id,jd,kmax))
              allocate (Sw_flux_save%flux_sw_down_total_dir_clr(id,jd))
              allocate (Sw_flux_save%flux_sw_down_total_dif_clr(id,jd))
              allocate (Sw_flux_save%flux_sw_down_vis_clr      (id,jd)) 
              if (do_swaerosol_forcing) then
                  allocate (Sw_flux_save%dfswcf_ad               (id,jd,kmax+1))
                  allocate (Sw_flux_save%ufswcf_ad               (id,jd,kmax+1))
              endif
          endif

          ! allocate space for special shortwave diagnostics
          allocate(Diag_special%swdn_trop (id,jd,MX_SPEC_LEVS))
          allocate(Diag_special%swup_trop (id,jd,MX_SPEC_LEVS))
          if (Rad_control%do_totcld_forcing) then
              allocate(Diag_special%swdn_trop_clr (id,jd,MX_SPEC_LEVS))
              allocate(Diag_special%swup_trop_clr (id,jd,MX_SPEC_LEVS))
          endif
      endif

!---------------------------------------------------------------------
!    allocate space for variables which must be saved when lw fluxes
!    are to be output on every physics step.
!---------------------------------------------------------------------
      if (all_step_diagnostics) then
          allocate (Lw_flux_save%olr             (id,jd))
          allocate (Lw_flux_save%lwups           (id,jd))
          allocate (Lw_flux_save%lwdns           (id,jd))
          allocate (Lw_flux_save%tdtlw           (id,jd,kmax))
          allocate (Lw_flux_save%flxnet          (id,jd,kmax+1))
          if (do_lwaerosol_forcing) then
              allocate (Lw_flux_save%olr_ad        (id,jd))
              allocate (Lw_flux_save%lwups_ad      (id,jd))
              allocate (Lw_flux_save%lwdns_ad      (id,jd))
          endif
          if (Rad_control%do_totcld_forcing) then
              allocate (Lw_flux_save%olr_clr           (id,jd))
              allocate (Lw_flux_save%lwups_clr         (id,jd))
              allocate (Lw_flux_save%lwdns_clr         (id,jd))
              allocate (Lw_flux_save%tdtlw_clr         (id,jd,kmax))
              allocate (Lw_flux_save%flxnetcf          (id,jd,kmax+1))
              if (do_lwaerosol_forcing) then
                  allocate (Lw_flux_save%olr_ad_clr      (id,jd))
                  allocate (Lw_flux_save%lwups_ad_clr    (id,jd))
                  allocate (Lw_flux_save%lwdns_ad_clr    (id,jd))
              endif
          endif

          ! allocate space for special longwave diagnostics
          allocate(Diag_special%netlw_trop(id,jd,MX_SPEC_LEVS))
          if (Rad_control%do_totcld_forcing) then
              allocate(Diag_special%netlw_trop_clr(id,jd,MX_SPEC_LEVS))
          endif
      endif

!----------------------------------------------------------------------
!    define characteristics of desired diagnostic integrals. 
!----------------------------------------------------------------------
      call initialize_diagnostic_integrals (id, jd)

!----------------------------------------------------------------------
!    register the desired netcdf output variables with the 
!    diagnostics_manager.
!----------------------------------------------------------------------
      call diag_field_init (id, jd, Time, axes, &
                            Rad_control%do_totcld_forcing, &
                            Aerosolrad_control%do_swaerosol, &
                            Aerosolrad_control%do_lwaerosol)

!----------------------------------------------------------------------
      if  (Rad_control%using_restart_file) then
!----------------------------------------------------------------------
!    Register fields to be written out to restart file.
!    Add the solar fields needed to restart the solar interplator
!----------------------------------------------------------------------
        if (Rad_control%renormalize_sw_fluxes) then
          call solar_interp_register_restart('solar_interp.res.nc', &
                                             Rad_control%do_totcld_forcing)

!-----------------------------------------------------------------------
!    if a valid restart file exists, then read by calling restore_state
!-----------------------------------------------------------------------
          if ( file_exist('INPUT/solar_interp.res.nc')) then
            call restore_state(Tile_restart)
          endif
        endif
      endif ! (using_restart_file)

!-----------------------------------------------------------------------

      module_is_initialized = .true.

!-----------------------------------------------------------------------


end subroutine radiation_driver_diag_init


!######################################################################
! <SUBROUTINE NAME="update_rad_fields">
!  <OVERVIEW>
!    update_rad_fields defines the current radiative heating rate, 
!    surface long and short wave fluxes and cosine of zenith angle
!    to be returned to physics_driver, including renormalization 
!    effects when that option is activated.
!  </OVERVIEW>
!  <DESCRIPTION>
!    update_rad_fields defines the current radiative heating rate, 
!    surface long and short wave fluxes and cosine of zenith angle
!    to be returned to physics_driver, including renormalization 
!    effects when that option is activated.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call update_rad_fields (is, ie, js, je, Time_diag, Astro_phys,   &
!                              Sw_output, Astro, Rad_output, flux_ratio)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending i,j indices in global storage arrays
!  </IN>
!  <IN NAME="Time_diag" TYPE="time_type">
!      Time on next timestep, used as stamp for diag-
!                        nostic output  [ time_type  (days, seconds) ] 
!  </IN>
!  <INOUT NAME="Astro" TYPE="astronomy_type">
!   astronomical properties on model grid, usually
!                   valid over radiation timestep on entry, on exit are 
!                   valid over model timestep when renormalizing
!  </INOUT>
!  <IN NAME="Astro_phys" TYPE="astronomy_type">
!   astronomical properties on model grid, valid over 
!                   physics timestep, used when renormalizing sw fluxes
!  </IN>
!  <INOUT NAME="Rad_output" TYPE="rad_output_type">
!   Radiation output from radiation package, contains variables
!                     which are output from radiation_driver to the 
!                     calling routine, and then used elsewhere within
!                     the component models.
!  </INOUT>
!  <IN NAME="Sw_output" TYPE="sw_output_type">
!   shortwave radiation output data from the 
!                        sea_esf_rad radiation package  when that 
!                        package is active
!  </IN>
!  <OUT NAME="flux_ratio" TYPE="real">
!   factor to multiply the radiation step values of 
!                   sw fluxes and heating rates by in order to get
!                   current physics timestep values
!  </OUT>
! </SUBROUTINE>
!

subroutine update_rad_fields (is, ie, js, je, Time_diag, Astro, Astro_phys,   &
                              Rad_control, Sw_output, Rad_output, flux_ratio)

!---------------------------------------------------------------------
!    update_rad_fields defines the current radiative heating rate, 
!    surface long and short wave fluxes and cosine of zenith angle
!    to be returned to physics_driver, including renormalization 
!    effects when that option is activated.
!--------------------------------------------------------------------

integer,                 intent(in)    ::  is, ie, js, je
type(time_type),         intent(in)    ::  Time_diag
type(astronomy_type),    intent(in)    ::  Astro
type(astronomy_type),    intent(in)    ::  Astro_phys
type(radiation_control_type),       intent(in)    ::  Rad_control
type(sw_output_type), dimension(:), intent(in)    ::  Sw_output
type(rad_output_type),              intent(inout) ::  Rad_output
real,  dimension(:,:),              intent(out)   ::  flux_ratio

!-------------------------------------------------------------------
!  intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data 
!                   in the physics_window being integrated
!      Time_diag    time on next timestep, used as stamp for diag-
!                   nostic output  [ time_type  (days, seconds) ]  
!      Astro        astronomical properties on model grid,
!                   valid over radiation timestep
!                   [astronomy_type]
!      Astro_phys   astronomical properties on model grid, valid over 
!                   physics timestep, used when renormalizing sw fluxes
!                   [astronomy_type]
!      Sw_output    shortwave output variables on model grid,
!                   [sw_output_type]     
!
!  intent(inout) variables:
!
!      Rad_output   radiation output variables on model grid, valid
!                   on entry over either physics or radiation timestep, 
!                   on exit are valid over physics step when renormal-
!                   izing sw fluxes
!                   [rad_output_type]     
!
!  intent(out) variables:
!
!      flux_ratio   factor to multiply the radiation step values of 
!                   sw fluxes and heating rates by in order to get
!                   current physics timestep values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:
!
      real, dimension (is:ie, js:je, &
                       size(Rad_output%tdt_rad,3))  ::  tdtlw, tdtlw_clr
      integer   :: i, j, k

!---------------------------------------------------------------------
!  local variables:
!
!     tdtlw              longwave heating rate
!                        [ deg K sec(-1) ]
!     tdtlw_clr          longwave heating rate under clear sky 
!                        conditions
!                        [ deg K sec(-1) ]
!     i,j,k              do-loop indices
!
!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_diag_mod',  &
               'module has not been initialized', FATAL)

!---------------------------------------------------------------------

      if (Rad_control%renormalize_sw_fluxes) then

!----------------------------------------------------------------------
!    if sw fluxes are to be renormalized, save the heating rates, fluxes
!    and solar factor calculated on radiation steps.
!---------------------------------------------------------------------
        if (Rad_control%do_sw_rad) then
          solar_save(is:ie,js:je)  = Astro%solar(:,:)
          call solar_flux_save_init (is, ie,js, je, Sw_output, Rad_output, &
                                     Rad_control%do_totcld_forcing)

!---------------------------------------------------------------------
!    define the ratio of the solar factor valid over this physics step
!    to that valid over the current radiation timestep.
!---------------------------------------------------------------------
          do j=1,je-js+1
            do i=1,ie-is+1
              if (solar_save(i+is-1,j+js-1) /= 0.0) then
                flux_ratio(i, j) = Astro_phys%solar(i,j)/   &
                                            solar_save(i+is-1,j+js-1)
              else
                flux_ratio(i,j) = 0.0
              endif
            end do
          end do

!----------------------------------------------------------------------
!    on non-radiation steps define the ratio of the current solar factor
!    valid for this physics step to that valid for the last radiation 
!    step. 
!----------------------------------------------------------------------
        else 
          do j=1,je-js+1
            do i=1,ie-is+1
              if (solar_save(i+is-1,j+js-1) /= 0.0) then
                flux_ratio(i, j) = Astro_phys%solar(i,j)/   &
                                            solar_save(i+is-1,j+js-1)
              else
                flux_ratio(i,j) = 0.0
              endif
            end do
          end do
        endif  ! (do_sw_rad)

!---------------------------------------------------------------------
!    redefine the total and shortwave heating rates, along with surface
!    sw fluxes, as a result of the difference in solar factor (the 
!    relative earth-sun motion) between the current physics and current
!    radiation timesteps.
!---------------------------------------------------------------------
        tdtlw(:,:,:) = Sw_flux_save%tot_heating(is:ie,js:je,:) - Sw_flux_save%sw_heating(is:ie,js:je,:)
        do k=1, size(Rad_output%tdt_rad,3)
          Rad_output%tdtsw(is:ie,js:je,k) = Sw_flux_save%sw_heating(is:ie,js:je,k)*flux_ratio(:,:)
        end do
        do k=1, size(Rad_output%tdt_rad,3)+1
          Rad_output%ufsw(is:ie,js:je,k) = Sw_flux_save%ufsw(is:ie,js:je,k)*flux_ratio(:,:)
          Rad_output%dfsw(is:ie,js:je,k) = Sw_flux_save%dfsw(is:ie,js:je,k)*flux_ratio(:,:)
        end do
        Rad_output%tdt_rad(is:ie,js:je,:) = tdtlw(:,:,:) + Rad_output%tdtsw(is:ie,js:je,:)
        Rad_output%flux_sw_surf          (is:ie,js:je) = flux_ratio(:,:)*Sw_flux_save%flux_sw_surf          (is:ie,js:je)
        Rad_output%flux_sw_surf_dir      (is:ie,js:je) = flux_ratio(:,:)*Sw_flux_save%flux_sw_surf_dir      (is:ie,js:je)
        Rad_output%flux_sw_surf_refl_dir (is:ie,js:je) = flux_ratio(:,:)*Sw_flux_save%flux_sw_surf_refl_dir (is:ie,js:je)
        Rad_output%flux_sw_surf_dif      (is:ie,js:je) = flux_ratio(:,:)*Sw_flux_save%flux_sw_surf_dif      (is:ie,js:je)
        Rad_output%flux_sw_down_vis_dir  (is:ie,js:je) = flux_ratio(:,:)*Sw_flux_save%flux_sw_down_vis_dir  (is:ie,js:je)
        Rad_output%flux_sw_down_vis_dif  (is:ie,js:je) = flux_ratio(:,:)*Sw_flux_save%flux_sw_down_vis_dif  (is:ie,js:je)
        Rad_output%flux_sw_down_total_dir(is:ie,js:je) = flux_ratio(:,:)*Sw_flux_save%flux_sw_down_total_dir(is:ie,js:je)
        Rad_output%flux_sw_down_total_dif(is:ie,js:je) = flux_ratio(:,:)*Sw_flux_save%flux_sw_down_total_dif(is:ie,js:je)
        Rad_output%flux_sw_vis           (is:ie,js:je) = flux_ratio(:,:)*Sw_flux_save%flux_sw_vis           (is:ie,js:je)
        Rad_output%flux_sw_vis_dir       (is:ie,js:je) = flux_ratio(:,:)*Sw_flux_save%flux_sw_vis_dir       (is:ie,js:je)
        Rad_output%flux_sw_refl_vis_dir  (is:ie,js:je) = flux_ratio(:,:)*Sw_flux_save%flux_sw_refl_vis_dir  (is:ie,js:je)
        Rad_output%flux_sw_vis_dif       (is:ie,js:je) = flux_ratio(:,:)*Sw_flux_save%flux_sw_vis_dif       (is:ie,js:je)
        if (Rad_control%do_totcld_forcing) then
          tdtlw_clr(:,:,:) = Sw_flux_save%tot_heating_clr(is:ie,js:je,:) - Sw_flux_save%sw_heating_clr (is:ie,js:je,:)
          Rad_output%flux_sw_down_total_dir_clr(is:ie,js:je) = flux_ratio(:,:)*Sw_flux_save%flux_sw_down_total_dir_clr(is:ie,js:je)
          Rad_output%flux_sw_down_total_dif_clr(is:ie,js:je) = flux_ratio(:,:)*Sw_flux_save%flux_sw_down_total_dif_clr(is:ie,js:je)
          Rad_output%flux_sw_down_vis_clr      (is:ie,js:je) = flux_ratio(:,:)*Sw_flux_save%flux_sw_down_vis_clr      (is:ie,js:je)
          do k=1, size(Rad_output%tdt_rad,3)
            Rad_output%tdtsw_clr(is:ie,js:je,k) = Sw_flux_save%sw_heating_clr (is:ie,js:je,k)*flux_ratio(:,:)
          end do
          do k=1, size(Rad_output%tdt_rad,3)+1
            Rad_output%ufsw_clr(is:ie,js:je,k) = Sw_flux_save%ufswcf(is:ie,js:je,k)*flux_ratio(:,:)
            Rad_output%dfsw_clr(is:ie,js:je,k) = Sw_flux_save%dfswcf(is:ie,js:je,k)*flux_ratio(:,:)
          end do
          Rad_output%tdt_rad_clr(is:ie,js:je,:) = tdtlw_clr(:,:,:) + Rad_output%tdtsw_clr(is:ie,js:je,:)
        endif
      else if (all_step_diagnostics) then

!----------------------------------------------------------------------
!    if sw fluxes are to be output on every physics step, save the 
!    heating rates and fluxes calculated on radiation steps.
!---------------------------------------------------------------------
        if (Rad_control%do_sw_rad) then
          call solar_flux_save_init (is,ie,js,je, Sw_output, Rad_output, &
                                     Rad_control%do_totcld_forcing)
        endif
      else
        flux_ratio(:,:) = 1.0
      endif  ! (renormalize_sw_fluxes)

!--------------------------------------------------------------------

end subroutine update_rad_fields 

!#####################################################################
! <SUBROUTINE NAME="initialize_diagnostic_integrals">
!  <OVERVIEW>
!    initialize_diagnostic_integrals registers the desired integrals 
!    with diag_integral_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    initialize_diagnostic_integrals registers the desired integrals 
!    with diag_integral_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call initialize_diagnostic_integrals
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine initialize_diagnostic_integrals (id, jd)
integer, intent(in) :: id, jd

!---------------------------------------------------------------------
!    initialize_diagnostic_integrals registers the desired integrals 
!    with diag_integral_mod.
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    initialize standard global quantities for integral package. 
!----------------------------------------------------------------------
      call diag_integral_field_init ('olr',    std_digits)
      call diag_integral_field_init ('abs_sw', std_digits)
!     call diag_integral_field_init ('olr_clr',    std_digits)
!     call diag_integral_field_init ('abs_sw_clr', std_digits)

!---------------------------------------------------------------------
!    allocate space for the global integrals being accumulated in 
!    this module.
!---------------------------------------------------------------------
      allocate (olr_intgl(id,jd))
      allocate (swabs_intgl(id,jd))

!--------------------------------------------------------------------

end subroutine initialize_diagnostic_integrals

!######################################################################
! <SUBROUTINE NAME="diag_field_init">
!  <OVERVIEW>
!    diag_field_init registers the desired diagnostic fields with the
!    diagnostics manager.
!  </OVERVIEW>
!  <DESCRIPTION>
!    diag_field_init registers the desired diagnostic fields with the
!    diagnostics manager.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diag_field_init ( Time, axes )
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   Current time
!  </IN>
!  <IN NAME="axes" TYPE="integer">
!   diagnostic variable axes for netcdf files
!  </IN>
! </SUBROUTINE>
!
subroutine diag_field_init ( id, jd, Time, axes, &
                             do_totcld_forcing, do_swaerosol, do_lwaerosol )

!---------------------------------------------------------------------
!    diag_field_init registers the desired diagnostic fields with the
!    diagnostics manager.
!---------------------------------------------------------------------

integer,         intent(in) :: id, jd
type(time_type), intent(in) :: Time
integer        , intent(in) :: axes(4)
logical,         intent(in) :: do_totcld_forcing
logical,         intent(in) :: do_swaerosol
logical,         intent(in) :: do_lwaerosol

!--------------------------------------------------------------------
!  intent(in) variables
!
!      Time        current time
!      axes        data axes for use with diagnostic fields
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables

      character(len=8)  ::   clr
      character(len=16) ::   clr2, lwaer_prep, swaer_prep
      integer           ::   bxes(4)
      integer           ::   i, k, n
      character(len=16) ::   spec_names(MX_SPEC_LEVS)
      character(len=24) ::   spec_long_names(MX_SPEC_LEVS)
!BW   integer           ::   id_lev, cmip_axes(4)

!--------------------------------------------------------------------
!  local variables:
!
!       clr          character string used in netcdf variable short name
!       clr2         character string used in netcdf variable long name
!       n            number of passes through name generation loop
!       i            do-loop index
!
!--------------------------------------------------------------------
!    names for special diagnostics fields
!-------------------------------------------------------------------
     spec_names = (/ character(len=16) :: '200hPa', 'lin_trop', 'therm_trop', '1_Pa' /)
     spec_long_names = (/ character(len=24) :: '200 hPa', 'linear tropopause', &
                          'thermo tropopause', '1 hPa' /)

!-------------------------------------------------------------------
!    define variable axis array with elements (1:3) valid for variables
!    defined at flux levels.
!-------------------------------------------------------------------
      bxes(1:2) = axes(1:2)
      bxes(3) = axes(4)
      bxes(4) = axes(4)

!---------------------------------------------------------------------
!    determine how many passes are needed through the name generation 
!    loop. 
!---------------------------------------------------------------------
      if (do_totcld_forcing) then
        n= 2
      else
        n= 1
      endif

      if (do_swaerosol ) then
        swaer_prep = 'without'
      else
        swaer_prep = 'with'
      endif
      if (do_lwaerosol ) then
        lwaer_prep = 'without'
      else
        lwaer_prep = 'with'
      endif

!---------------------------------------------------------------------
!    generate names for standard and clear sky diagnostic fields. if 
!    clear sky values being generated, generate the clear sky names
!    on pass 1, followed by the standard names.
!---------------------------------------------------------------------
      do i = 1, n
        if ( i == n) then
          clr  = "    "
          clr2 = "          "
        else
          clr  = "_clr"
          clr2 = "clear sky "
        endif

        ! register special fluxes
        do k = 1, MX_SPEC_LEVS
          id_swdn_special(k,i) = register_diag_field (mod_name,   &
              'swdn_'//trim(spec_names(k))//trim(clr), axes(1:2), Time, &
              trim(clr2)//'SW flux down at '//trim(spec_long_names(k)), &
              'watts/m2', missing_value=missing_value)

          id_swup_special(k,i) = register_diag_field (mod_name,   &
              'swup_'//trim(spec_names(k))//trim(clr), axes(1:2), Time, &
              trim(clr2)//'SW flux up at '//trim(spec_long_names(k)), &
              'watts/m2', missing_value=missing_value)

          id_netlw_special(1,i) = register_diag_field (mod_name,   &
              'netlw_'//trim(spec_names(k))//trim(clr), axes(1:2), Time, &
              trim(clr2)//'net LW flux at '//trim(spec_long_names(k)), &
              'watts/m2', missing_value=missing_value)
        end do


        id_tdt_sw(i) = register_diag_field (mod_name,   &
                'tdt_sw'//trim(clr), axes(1:3), Time, & 
                trim(clr2)//'temperature tendency for SW radiation', &
                'deg_K/sec', missing_value=missing_value) 

        id_ufsw(i) = register_diag_field (mod_name,   &
               'allufsw'//trim(clr), bxes(1:3), Time, &
               trim(clr2)//'upward sw flux', &
               'watts/m2', missing_value=missing_value)

        id_dfsw(i) = register_diag_field (mod_name,   &
               'alldfsw'//trim(clr), bxes(1:3), Time, &
               trim(clr2)//'downward sw flux', &
               'watts/m2', missing_value=missing_value)

        id_flxnet(i) = register_diag_field (mod_name,   &
               'allnetlw'//trim(clr), bxes(1:3), Time, &
               trim(clr2)//'net lw flux', &
               'watts/m2', missing_value=missing_value)

        id_tdt_lw(i) = register_diag_field (mod_name,    &
                'tdt_lw'//trim(clr), axes(1:3), Time, &
                trim(clr2)//'temperature tendency for LW radiation', &
                'deg_K/sec', missing_value=missing_value)

        id_swdn_toa(i) = register_diag_field (mod_name,   &
                'swdn_toa'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux down at TOA', &
                'watts/m2', missing_value=missing_value)

        id_swup_toa(i) = register_diag_field (mod_name,    &
                'swup_toa'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux up at TOA', &
                'watts/m2', missing_value=missing_value)

        id_olr(i) = register_diag_field (mod_name,   &
                'olr'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'outgoing longwave radiation', &
                'watts/m2', missing_value=missing_value)

        id_netrad_toa(i) = register_diag_field (mod_name,   &
                'netrad_toa'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'net radiation (lw + sw) at toa', &
                'watts/m2', missing_value=missing_value)

        id_netrad_1_Pa(i) = register_diag_field (mod_name,   &
                'netrad_1_Pa'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'net radiation (lw + sw) at 1 Pa', &
                'watts/m2', missing_value=missing_value)

        id_swup_sfc(i) = register_diag_field (mod_name,    &
                'swup_sfc'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux up at surface', &
                'watts/m2', missing_value=missing_value)

        id_swdn_sfc(i) = register_diag_field (mod_name,     &
                'swdn_sfc'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux down at surface', &
                'watts/m2', missing_value=missing_value)

        id_lwup_sfc(i) = register_diag_field (mod_name,   &
                'lwup_sfc'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'LW flux up at surface', &
                'watts/m2', missing_value=missing_value)

        id_lwdn_sfc(i) = register_diag_field (mod_name,    &
                'lwdn_sfc'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'LW flux down at surface', &
                'watts/m2', missing_value=missing_value)

        id_swtoa(i) = register_diag_field (mod_name,    &
                 'swtoa'//trim(clr), axes(1:2), Time, &
                  trim(clr2)//' Net SW flux at TOA ', &
                  'watts/m2', missing_value=missing_value)

        id_swsfc(i) = register_diag_field (mod_name,    &
                  'swsfc'//trim(clr), axes(1:2), Time, &
                  trim(clr2)//' Net SW flux at surface', &
                  'watts/m2', missing_value=missing_value)
 
        id_lwsfc(i) = register_diag_field (mod_name,    &
                  'lwsfc'//trim(clr), axes(1:2), Time, &
                 trim(clr2)//' Net LW flux at surface', &
                  'watts/m2', missing_value=missing_value)
 
        id_swtoa_ad(i) = register_diag_field (mod_name,    &
                  'swtoa_ad'//trim(clr), axes(1:2), Time, &
                  trim(clr2)//' Net SW flux at TOA '// trim(swaer_prep) &
                                 // ' aerosol', &
                  'watts/m2', missing_value=missing_value)
 
        id_swsfc_ad(i) = register_diag_field (mod_name,    &
                 'swsfc_ad'//trim(clr), axes(1:2), Time, &
            trim(clr2)//' Net SW flux at surface '// trim(swaer_prep) &
           // ' aerosol', &
                'watts/m2', missing_value=missing_value)

       id_swdn_sfc_ad(i) = register_diag_field (mod_name,    &
                 'swdn_sfc_ad'//trim(clr), axes(1:2), Time, &
                 trim(clr2)//' SW flux down at surface '// &
                   trim(swaer_prep) // ' aerosol', &
                 'watts/m2', missing_value=missing_value)
 
        id_swup_sfc_ad(i) = register_diag_field (mod_name,    &
                 'swup_sfc_ad'//trim(clr), axes(1:2), Time, &
                 trim(clr2)//' SW flux up at surface ' //   &
                   trim(swaer_prep) // ' aerosol', &
                  'watts/m2', missing_value=missing_value)

        id_swup_toa_ad(i) = register_diag_field (mod_name,    &
                 'swup_toa_ad'//trim(clr), axes(1:2), Time, &
                 trim(clr2)//' SW flux up at TOA '  //  &
                   trim(swaer_prep) // ' aerosol', &
                 'watts/m2', missing_value=missing_value)
 
         id_olr_ad(i) = register_diag_field (mod_name,    &
                  'lwtoa_ad'//trim(clr), axes(1:2), Time, &
                  trim(clr2)//' Net LW flux at TOA (olr) ' //  &
                   trim(lwaer_prep) // ' aerosol', &
                  'watts/m2', missing_value=missing_value)

         id_lwsfc_ad(i) = register_diag_field (mod_name,    &
                  'lwsfc_ad'//trim(clr), axes(1:2), Time, &
                  trim(clr2)//' Net LW flux at surface  ' //   &
                   trim(lwaer_prep) // ' aerosol', &
                 'watts/m2', missing_value=missing_value)
 
       end do

        ID_tntr = register_cmip_diag_field_3d ( mod_name, 'tntr', Time, &
           'Tendency of Air Temperature due to Radiative Heating', 'K s-1', &
            standard_name = 'tendency_of_air_temperature_due_to_radiative_heating' )

        ID_tntrs = register_cmip_diag_field_3d ( mod_name, 'tntrs', Time, &
           'Tendency of Air Temperature due to Shortwave Radiative Heating', 'K s-1', &
            standard_name = 'tendency_of_air_temperature_due_to_shortwave_heating' )

        ID_tntrscs = register_cmip_diag_field_3d (mod_name, 'tntrscs', Time, & 
           'Tendency of Air Temperature due to Clear Sky Shortwave Radiative Heating', 'K s-1', &
            standard_name = 'tendency_of_air_temperature_due_to_clear_sky_shortwave_heating' )

        ID_rsu = register_cmip_diag_field_3d (mod_name, 'rsu', Time, & 
                             'Upwelling Shortwave Radiation', 'W m-2', &
            standard_name = 'upwelling_shortwave_flux_in_air', axis='half' )

        ID_rsucs = register_cmip_diag_field_3d (mod_name, 'rsucs', Time, & 
                             'Upwelling Clear-Sky Shortwave Radiation', 'W m-2', &
            standard_name = 'upwelling_shortwave_flux_in_air_assuming_clear_sky', axis='half' )

        ID_rsuaf = register_cmip_diag_field_3d (mod_name, 'rsuaf', Time, & 
                             'Upwelling Clean-Sky Shortwave Radiation at each level', 'W m-2', &
            standard_name = 'upwelling_shortwave_flux_in_air_assuming_clean_sky', axis='half' )

        ID_rsucsaf = register_cmip_diag_field_3d (mod_name, 'rsucsaf', Time, & 
                             'Upwelling Clean-Clear-Sky Shortwave Radiation at each level', 'W m-2', &
            standard_name = 'upwelling_shortwave_flux_in_air_assuming_clean_clear_sky', axis='half' )

        ID_rsd = register_cmip_diag_field_3d (mod_name, 'rsd', Time, & 
                             'Downwelling Shortwave Radiation', 'W m-2', &
            standard_name = 'downwelling_shortwave_flux_in_air', axis='half' )

        ID_rsdcs = register_cmip_diag_field_3d (mod_name, 'rsdcs', Time, & 
                             'Downwelling Clear-Sky Shortwave Radiation', 'W m-2', &
            standard_name = 'downwelling_shortwave_flux_in_air_assuming_clear_sky', axis='half' )

        ID_rsdaf = register_cmip_diag_field_3d (mod_name, 'rsdaf', Time, & 
                             'Downwelling Clean-Sky Shortwave Radiation at each level', 'W m-2', &
            standard_name = 'downwelling_shortwave_flux_in_air_assuming_clean_sky', axis='half' )

        ID_rsdcsaf = register_cmip_diag_field_3d (mod_name, 'rsdcsaf', Time, & 
                             'Downwelling Clean-Clear-Sky Shortwave Radiation at each level', 'W m-2', &
            standard_name = 'downwelling_shortwave_flux_in_air_assuming_clean_clear_sky', axis='half' )

        ID_tntrl = register_cmip_diag_field_3d (mod_name, 'tntrl', Time, & 
           'Tendency of Air Temperature due to Longwave Radiative Heating', 'K s-1', &
            standard_name = 'tendency_of_air_temperature_due_to_longwave_heating' )

        ID_tntrlcs = register_cmip_diag_field_3d (mod_name, 'tntrlcs', Time, & 
           'Tendency of Air Temperature due to Clear Sky Longwave Radiative Heating', 'K s-1', &
            standard_name = 'tendency_of_air_temperature_due_to_clear_sky_longwave_heating' )

        id_rlds = register_cmip_diag_field_2d (mod_name, 'rlds', Time, &
                'Surface Downwelling Longwave Radiation', 'W m-2', &
            standard_name = 'surface_downwelling_longwave_flux_in_air')

        id_rldscs = register_cmip_diag_field_2d (mod_name, 'rldscs', Time, &
                'Surface Downwelling Clear-Sky Longwave Radiation',  'W m-2', &
           standard_name = 'surface_downwelling_longwave_flux_in_air_assuming_clear_sky')

        id_rldsaf = register_cmip_diag_field_2d (mod_name, 'rldsaf', Time, &
                'Surface Downwelling Aerosol-Free Longwave Radiation',  'W m-2', &
           standard_name = 'surface_downwelling_longwave_flux_in_air_assuming_clean_sky')

        id_rldscsaf = register_cmip_diag_field_2d (mod_name, 'rldscsaf', Time, &
                'Surface Downwelling Clear-Sky, Aerosol-Free Longwave Radiation',  'W m-2', &
           standard_name = 'surface_downwelling_longwave_flux_in_air_assuming_clean_clear_sky')

        id_rlus = register_cmip_diag_field_2d (mod_name, 'rlus', Time, &
                'Surface Upwelling Longwave Radiation', 'W m-2', &
            standard_name = 'surface_upwelling_longwave_flux_in_air')

        id_rsds = register_cmip_diag_field_2d (mod_name, 'rsds', Time, &
                'Surface Downwelling Shortwave Radiation', 'W m-2', &
            standard_name = 'surface_downwelling_shortwave_flux_in_air')

        id_rsdscs = register_cmip_diag_field_2d (mod_name, 'rsdscs', Time, &
                'Surface Downwelling Clear-Sky Shortwave Radiation', 'W m-2', &
           standard_name = 'surface_downwelling_shortwave_flux_in_air_assuming_clear_sky')

        id_rsdsaf = register_cmip_diag_field_2d (mod_name, 'rsdsaf', Time, &
                'Surface Downwelling Aerosol-Free Shortwave Radiation', 'W m-2', &
           standard_name = 'surface_downwelling_shortwave_flux_in_air_assuming_clean_sky')

        id_rsdscsaf = register_cmip_diag_field_2d (mod_name, 'rsdscsaf', Time, &
                'Surface Downwelling Clear-Sky, Aerosol-Free Shortwave Radiation', 'W m-2', &
           standard_name = 'surface_downwelling_shortwave_flux_in_air_assuming_clean_clear_sky')

        id_rsus = register_cmip_diag_field_2d (mod_name, 'rsus', Time, &
                'Surface Upwelling Shortwave Radiation', 'W m-2', &
            standard_name = 'surface_upwelling_shortwave_flux_in_air')

        id_rsuscs = register_cmip_diag_field_2d (mod_name, 'rsuscs', Time, &
                'Surface Upwelling Clear-Sky Shortwave Radiation', 'W m-2', &
           standard_name = 'surface_upwelling_shortwave_flux_in_air_assuming_clear_sky')

        id_rsusaf = register_cmip_diag_field_2d (mod_name, 'rsusaf', Time, &
                'Surface Upwelling Aerosol-Free Shortwave Radiation', 'W m-2', &
           standard_name = 'surface_upwelling_shortwave_flux_in_air_assuming_clean_sky')

        id_rsuscsaf = register_cmip_diag_field_2d (mod_name, 'rsuscsaf', Time, &
                'Surface Upwelling Clear-Sky, Aerosol-Free Shortwave Radiation', 'W m-2', &
           standard_name = 'surface_upwelling_shortwave_flux_in_air_assuming_clean_clear_sky')

        id_rsdt = register_cmip_diag_field_2d (mod_name, 'rsdt', Time, &
                        'TOA Incident Shortwave Radiation', 'W m-2',   &
                         standard_name = 'toa_incoming_shortwave_flux')

        id_rsut = register_cmip_diag_field_2d (mod_name, 'rsut', Time, &
                       'TOA Outgoing Shortwave Radiation', 'W m-2',    &
                        standard_name = 'toa_outgoing_shortwave_flux')

        id_rsutcs = register_cmip_diag_field_2d (mod_name, 'rsutcs', Time, &
                 'TOA Outgoing Clear-Sky Shortwave Radiation', 'W m-2',    &
                 standard_name = 'toa_outgoing_shortwave_flux_assuming_clear_sky')

        id_rsutaf = register_cmip_diag_field_2d (mod_name, 'rsutaf', Time, &
                 'TOA Outgoing Aerosol-Free Shortwave Radiation', 'W m-2',    &
                 standard_name = 'toa_outgoing_shortwave_flux_in_air_assuming_clean_sky')

        id_rsutcsaf = register_cmip_diag_field_2d (mod_name, 'rsutcsaf', Time, &
                 'TOA Outgoing Clear-Sky, Aerosol-Free Shortwave Radiation', 'W m-2',    &
                 standard_name = 'toa_outgoing_shortwave_flux_in_air_assuming_clean_clear_sky')

        id_rlut = register_cmip_diag_field_2d (mod_name, 'rlut', Time, &
                        'TOA Outgoing Longwave Radiation', 'W m-2',    &
                         standard_name = 'toa_outgoing_longwave_flux')

        id_rlutcs = register_cmip_diag_field_2d (mod_name, 'rlutcs', Time, &
                   'TOA Outgoing Clear-Sky Longwave Radiation', 'W m-2',   &
                 standard_name = 'toa_outgoing_longwave_flux_assuming_clear_sky')

        id_rlutaf = register_cmip_diag_field_2d (mod_name, 'rlutaf', Time, &
                   'TOA Outgoing Aerosol-Free Longwave Radiation', 'W m-2',   &
                 standard_name = 'toa_outgoing_longwave_flux_in_air_assuming_clean_sky')

        id_rlutcsaf = register_cmip_diag_field_2d (mod_name, 'rlutcsaf', Time, &
                   'TOA Outgoing Clear-Sky, Aerosol-Free Longwave Radiation', 'W m-2',   &
                 standard_name = 'toa_outgoing_longwave_flux_in_air_assuming_clean_clear_sky')

        id_rtmt = register_cmip_diag_field_2d (mod_name, 'rtmt', Time, &
                        'Net Downward Flux at Top of Model', 'W m-2',  &
                standard_name = 'net_downward_radiative_flux_at_top_of_atmosphere_model')
  
        id_rsdsdiff = register_cmip_diag_field_2d (mod_name, 'rsdsdiff', Time, &
               'Surface Diffuse Downwelling Shortwave Radiation', 'W m-2',  &
               standard_name = 'surface_diffuse_downwelling_shortwave_flux_in_air')
  
        id_rsdscsdiff = register_cmip_diag_field_2d (mod_name, 'rsdscsdiff', Time, &
               'Surface Diffuse Downwelling Clear Sky Shortwave Radiation', 'W m-2',  &
               standard_name = 'surface_diffuse_downwelling_shortwave_flux_in_air_assuming_clear_sky')


         id_allradp   = register_diag_field (mod_name,   &
                 'allradp', axes(1:3), Time, &
                 'temperature tendency for SW + LW radiation', &
                 'deg_K/sec', missing_value=missing_value)

         id_heat2d   = register_diag_field (mod_name,   &
                 'heat2d_rad', axes(1:2), Time, &
                 'integrated net radiative heating', 'watts/m2')

         id_heat2d_sw = register_diag_field (mod_name,   &
                 'heat2d_sw', axes(1:2), Time, &
                 'integrated shortwave radiative heating', 'watts/m2')

!----- initialize global integrals for netCDF output -----

        id_rsdt_g = register_global_diag_field ('rsdt', Time, &
                        'TOA Incident Shortwave Radiation', 'W m-2',   &
                       standard_name='toa_incoming_shortwave_flux', buffer=.true.)

        id_rsut_g = register_global_diag_field ('rsut', Time, &
                       'TOA Outgoing Shortwave Radiation', 'W m-2',    &
                      standard_name='toa_outgoing_shortwave_flux', buffer=.true.)

        id_rsutcs_g = register_global_diag_field ('rsutcs', Time, &
                 'TOA Outgoing Clear-Sky Shortwave Radiation', 'W m-2',    &
                 standard_name='toa_outgoing_shortwave_flux_assuming_clear_sky', buffer=.true.)

        id_rlut_g = register_global_diag_field ('rlut', Time, &
                        'TOA Outgoing Longwave Radiation', 'W m-2',    &
                      standard_name='toa_outgoing_longwave_flux', buffer=.true.)

        id_rlutcs_g = register_global_diag_field ('rlutcs', Time, &
                   'TOA Outgoing Clear-Sky Longwave Radiation', 'W m-2',   &
                 standard_name='toa_outgoing_longwave_flux_assuming_clear_sky', buffer=.true.)

        id_rss_g = register_global_diag_field ('rss', Time, &
                        'Net shortwave radiation at the surface', 'W m-2',   &
                       standard_name='surface_net_shortwave_flux', buffer=.true.)


!----------------------------------------------------------------------
!    register fields that are not clear-sky depedent.
!----------------------------------------------------------------------
        id_conc_drop = register_diag_field (mod_name,   &
                   'conc_drop', axes(1:3), Time, & 
                   'drop concentration ', &
                   'g/m^3', missing_value=missing_value) 
    
        id_conc_ice = register_diag_field (mod_name,   &
                   'conc_ice', axes(1:3), Time, & 
                   'ice concentration ', &
                   'g/m^3', missing_value=missing_value) 
 
      
      id_flux_sw_dir = register_diag_field (mod_name,    &
                'flux_sw_dir', axes(1:2), Time, &
                'net direct sfc sw flux', 'watts/m2', &
                missing_value=missing_value)

      id_flux_sw_refl_dir = register_diag_field (mod_name,    &
                'flux_sw_refl_dir', axes(1:2), Time, &
                'refl sw dir from sfc', 'watts/m2', &
                missing_value=missing_value)

      id_flux_sw_refl_vis_dir = register_diag_field (mod_name,    &
                'flux_sw_refl_vis_dir', axes(1:2), Time, &
                'refl sw vis dir from sfc', 'watts/m2', &
                missing_value=missing_value)

      id_flux_sw_dif = register_diag_field (mod_name,    &
                'flux_sw_dif', axes(1:2), Time, &
                'net diffuse sfc sw flux', 'watts/m2', &
                missing_value=missing_value)

      id_flux_sw     = register_diag_field (mod_name,    &
                'flux_sw', axes(1:2), Time, &
                'net sfc sw flux', 'watts/m2', &
                missing_value=missing_value)

      id_flux_sw_down_vis_dir = register_diag_field (mod_name,    &
                'flux_sw_down_vis_dir', axes(1:2), Time, &
                'downward direct visible sfc sw flux', 'watts/m2', &
                 missing_value=missing_value)

      id_flux_sw_down_vis_dif = register_diag_field (mod_name,    &
                'flux_sw_down_vis_dif', axes(1:2), Time, &
                'downward diffuse visible sfc sw flux', 'watts/m2', &
                 missing_value=missing_value)

      id_flux_sw_down_total_dir = register_diag_field (mod_name,    &
                'flux_sw_down_total_dir', axes(1:2), Time, &
                'downward direct total sfc sw flux', 'watts/m2', &
                 missing_value=missing_value)

      id_flux_sw_down_total_dif = register_diag_field (mod_name,    &
               'flux_sw_down_total_dif', axes(1:2), Time, &
               'downward diffuse total sfc sw flux', 'watts/m2', &
                missing_value=missing_value)

    if (do_totcld_forcing) then

      id_flux_sw_down_total_dir_clr = register_diag_field (mod_name,  &
               'flux_sw_down_total_dir_clr', axes(1:2), Time, &
               'downward clearsky direct total sfc sw flux',  &
               'watts/m2',  missing_value=missing_value)
  
      id_flux_sw_down_total_dif_clr = register_diag_field (mod_name,  &
               'flux_sw_down_total_dif_clr', axes(1:2), Time, &
               'downward clearsky diffuse total sfc sw flux',  &
               'watts/m2', missing_value=missing_value)

      id_flux_sw_down_vis_clr = register_diag_field (mod_name,    &
                'flux_sw_down_vis_clr', axes(1:2), Time, &
                'downward visible sfc sw flux clear sky', 'watts/m2', &
                 missing_value=missing_value)

    endif 

      id_flux_sw_vis = register_diag_field (mod_name,    &
               'flux_sw_vis', axes(1:2), Time, &
               'net visible sfc sw flux', 'watts/m2', &
                 missing_value=missing_value)

      id_flux_sw_vis_dir = register_diag_field (mod_name,    &
               'flux_sw_vis_dir', axes(1:2), Time, &
               'net direct visible sfc sw flux', 'watts/m2', &
                 missing_value=missing_value)

      id_flux_sw_vis_dif = register_diag_field (mod_name,    &
                'flux_sw_vis_dif', axes(1:2), Time, &
                'net diffuse visible sfc sw flux', 'watts/m2', &
                  missing_value=missing_value)


      id_solar_constant = register_diag_field (mod_name,    &
                  'solar_constant', Time, &
                  'solar constant', 'watts/m2', &
                  missing_value=missing_value)      
                           
      id_co2_tf = register_diag_field (mod_name,    &
                  'co2_tf', Time, &
                  'co2 mixing ratio used for tf calculation', 'ppmv', &
                  missing_value=missing_value)      
                           
      id_ch4_tf = register_diag_field (mod_name,    &
                  'ch4_tf', Time, &
                  'ch4 mixing ratio used for tf calculation', 'ppbv', &
                  missing_value=missing_value)      
                           
      id_n2o_tf = register_diag_field (mod_name,    &
                  'n2o_tf', Time, &
                  'n2o mixing ratio used for tf calculation', 'ppbv', &
                  missing_value=missing_value)      
                           
      id_rrvco2 = register_diag_field (mod_name,    &
                  'rrvco2', Time, &
                  'co2 mixing ratio', 'ppmv', &
                  missing_value=missing_value)      
                           
      id_rrvf11 = register_diag_field (mod_name,    &
                  'rrvf11', Time, &
                  'f11 mixing ratio', 'pptv', &
                  missing_value=missing_value)
        
      id_rrvf12 = register_diag_field (mod_name,    &
                  'rrvf12', Time, &
                  'f12 mixing ratio', 'pptv', &
                  missing_value=missing_value)

      id_rrvf113 = register_diag_field (mod_name,    &
                   'rrvf113', Time, &
                   'f113 mixing ratio', 'pptv', &
                   missing_value=missing_value)
 
       id_rrvf22 = register_diag_field (mod_name,    &
                   'rrvf22', Time, &
                   'f22 mixing ratio', 'pptv', &
                   missing_value=missing_value)

       id_rrvch4 = register_diag_field (mod_name,    &
                   'rrvch4', Time, &
                   'ch4 mixing ratio', 'ppbv', &
                   missing_value=missing_value)

       id_rrvn2o = register_diag_field (mod_name,    &
                   'rrvn2o', Time, &
                   'n2o mixing ratio', 'ppbv', &
                   missing_value=missing_value)

      id_co2mass = register_diag_field (mod_name, 'co2mass', Time, &
                  'Total Atmospheric Mass of CO2', 'kg', &
                  standard_name = 'atmosphere_mass_of_carbon_dioxide', &
                  missing_value=CMOR_MISSING_VALUE) 
                           
      id_cfc11global = register_diag_field (mod_name,    &
                  'cfc11global', Time, &
                  'Global Mean Mole Fraction of CFC11', '1e-12', &
                  standard_name = 'mole_fraction_of_cfc11_in_air', &
                  missing_value=CMOR_MISSING_VALUE) 
        
      id_cfc12global = register_diag_field (mod_name,    &
                  'cfc12global', Time, &
                  'Global Mean Mole Fraction of CFC12', '1e-12', &
                  standard_name = 'mole_fraction_of_cfc12_in_air', &
                  missing_value=CMOR_MISSING_VALUE) 

      id_cfc113global = register_diag_field (mod_name,    &
                   'cfc113global', Time, &
                  'Global Mean Mole Fraction of CFC113', '1e-12', &
                  standard_name = 'mole_fraction_of_cfc113_in_air', &
                  missing_value=CMOR_MISSING_VALUE) 
 
       id_hcfc22global = register_diag_field (mod_name,    &
                   'hcfc22global', Time, &
                  'Global Mean Mole Fraction of HCFC22', '1e-12', &
                  standard_name = 'mole_fraction_of_hcfc22_in_air', &
                  missing_value=CMOR_MISSING_VALUE) 

       id_ch4global = register_diag_field (mod_name,    &
                   'ch4global', Time, &
                  'Global Mean Mole Fraction of CH4', '1e-09', &
                  standard_name = 'mole_fraction_of_methane_in_air', &
                  missing_value=CMOR_MISSING_VALUE) 

       id_n2oglobal = register_diag_field (mod_name,    &
                   'n2oglobal', Time, &
                  'Global Mean Mole Fraction of N2O', '1e-09', &
                standard_name = 'mole_fraction_of_nitrous_oxide_in_air', &
                  missing_value=CMOR_MISSING_VALUE)

         id_alb_sfc_avg = register_diag_field (mod_name,    &
                 'averaged_alb_sfc', axes(1:2), Time, &
                 'surface albedo', 'percent', &
                 missing_value=missing_value)
         if (id_alb_sfc_avg > 0) then
           allocate (swdns_acc(id,jd))
           allocate (swups_acc(id,jd))
           swups_acc = 0.0
           swdns_acc = 1.0e-35
         endif
      id_alb_sfc = register_diag_field (mod_name,    &
                'alb_sfc', axes(1:2), Time, &
                'surface albedo', 'percent', &
                  missing_value=missing_value) 

      id_alb_sfc_vis_dir = register_diag_field (mod_name,    &
                'alb_sfc_vis_dir', axes(1:2), Time, &
!               'surface albedo_vis_dir', 'percent')
! BUGFIX
                'surface albedo_vis_dir', 'percent', &
                 missing_value=missing_value)
      id_alb_sfc_nir_dir = register_diag_field (mod_name,    &
                'alb_sfc_nir_dir', axes(1:2), Time, &
!               'surface albedo_nir', 'percent')
! BUGFIX
                'surface albedo_nir_dir', 'percent', &
                  missing_value=missing_value)
 
      id_alb_sfc_vis_dif = register_diag_field (mod_name,    &
                 'alb_sfc_vis_dif', axes(1:2), Time, &
!               'surface albedo_vis', 'percent')
! BUGFIX
                 'surface albedo_vis_dif', 'percent', &
                  missing_value=missing_value)
      id_alb_sfc_nir_dif = register_diag_field (mod_name,    &
                 'alb_sfc_nir_dif', axes(1:2), Time, &
!               'surface albedo_nir', 'percent')
! BUGFIX
                 'surface albedo_nir_dif', 'percent', &
                   missing_value=missing_value)
      id_cosz = register_diag_field (mod_name,    &
                'cosz',axes(1:2),  Time,    &
                'cosine of zenith angle',    &
                'none', missing_value=missing_value)

      id_fracday = register_diag_field (mod_name,   &
                'fracday',axes(1:2), Time,   &
                'daylight fraction of radiation timestep',   &
                'percent', missing_value=missing_value)
        

!-----------------------------------------------------------------------

end subroutine diag_field_init

!######################################################################
! <SUBROUTINE NAME="produce_radiation_diagnostics">
!  <OVERVIEW>
!    produce_radiation_diagnostics produces netcdf output and global 
!    and hemispheric integrals of radiation package variables.
!  </OVERVIEW>
!  <DESCRIPTION>
!    produce_radiation_diagnostics produces netcdf output and global 
!    and hemispheric integrals of radiation package variables.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call produce_radiation_diagnostics          &
!                            (is, ie, js, je, Time_diag, lat, ts, asfc, &
!                             flux_ratio, Astro, Rad_output, Lw_output, &
!                             Sw_output, Lsc_microphys)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending i,j indices in global storage arrays
!  </IN> 

!  <IN NAME="Time_diag" TYPE="time_type">
!      Time_diag         time on next timestep, used as stamp for diag-
!                        nostic output  [ time_type  (days, seconds) ] 
!  </IN>

!  <IN NAME="lat" TYPE="real">
!    lat        mean latitude (in radians) of all grid boxes processed by this
!               call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <IN NAME="ts" TYPE="real">
!   Surface skin temperature
!  </IN>
!  <IN NAME="asfc" TYPE="real">
!   surface albedo
!  </IN>
!  <IN NAME="flux_ratio" TYPE="real">
!   renormalization factor for sw fluxes and heating 
!                   rates
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
!   astronomical input data for the radiation package
!  </IN>
!  <IN NAME="Rad_output" TYPE="rad_output_type">
!   Radiation output from radiation package, contains variables
!                     which are output from radiation_driver to the 
!                     calling routine, and then used elsewhere within
!                     the component models.
!  </IN>
!  <IN NAME="Lw_output" TYPE="lw_output_type">
!      longwave radiation output data from the 
!                        sea_esf_rad radiation package, when that 
!                        package is active
!  </IN>
!  <IN NAME="Sw_output" TYPE="sw_output_type">
!   shortwave radiation output data from the 
!                        sea_esf_rad radiation package  when that 
!                        package is active
!  </IN>
!  <IN NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale clouds,
!                        when the microphysical package is active
!  </IN>
! </SUBROUTINE>
!
subroutine produce_radiation_diagnostics          &
                 (is, ie, js, je, Time_diag, Time, lat, atm_mass, ts, pflux, phalf, &
                  asfc_vis_dir, asfc_nir_dir, asfc_vis_dif, asfc_nir_dif, &
                  flux_ratio, Astro, Astro_phys, Rad_output, Rad_gases,&
                  Rad_control, Lw_output, Sw_output)
!BW               Lsc_microphys)

!--------------------------------------------------------------------
!    produce_radiation_diagnostics produces netcdf output and global 
!    and hemispheric integrals of radiation package variables.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
integer,                 intent(in)             :: is, ie, js, je
type(time_type),         intent(in)             :: Time_diag
type(time_type),         intent(in)             :: Time
real,dimension(:,:),     intent(in)             :: lat, ts
real,                    intent(in)             :: atm_mass
real,dimension(:,:,:),   intent(in)             :: pflux, phalf
real, dimension(:,:),    intent(in)             :: asfc_vis_dir, &
                                                   asfc_nir_dir, &
                                                   asfc_vis_dif, &
                                                   asfc_nir_dif
real,dimension(:,:),     intent(in)             :: flux_ratio
type(astronomy_type),    intent(in)             :: Astro, Astro_phys
type(rad_output_type),   intent(in)             :: Rad_output
type(radiative_gases_type), intent(in)          :: Rad_gases
type(radiation_control_type),  intent(in)       :: Rad_control
type(lw_output_type), dimension(:), intent(in), optional :: Lw_output
type(sw_output_type), dimension(:), intent(in), optional :: Sw_output
!BW type(microphysics_type), intent(in), optional   :: Lsc_microphys
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data 
!                   in the physics_window being integrated
!      Time_diag    time on next timestep, used as stamp for diagnostic 
!                   output  [ time_type  (days, seconds) ]  
!      lat          latitude of model points  [ radians ]
!      ts           surface temperature  [ deg K ]
!      asfc         surface albedo  [ dimensionless ]
!      flux_ratio   renormalization factor for sw fluxes and heating 
!                   rates [ dimensionless ]
!      Astro        astronomical  variables input to the radiation
!                   package [ dimensionless ]
!      Rad_output   rad_output_type variable containing radiation 
!                   output fields
!      Rad_gases    radiative_gases_type variable containing co2 mixing
!                   ratio
!
!
!    intent(in) optional variables:
!
!      Lw_output    lw_output_type variable containing output from 
!                   the longwave radiation code of the
!                   sea_esf_rad package, on the model grid
!      Sw_output    sw_output_type variable containing output from 
!                   the shortwave radiation code of the
!                   sea_esf_rad package, on the model grid
!      Lsc_microphys  microphysics_type structure, contains variables
!                     describing the microphysical properties of the
!                     large-scale clouds, passed through to lower
!                     level routines
!        
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  local variables

      real, dimension (ie-is+1,je-js+1) ::           & 
                                                swin, swout, olr, &
                                                swups, swdns, lwups, &
                                                lwdns, swin_clr,   &
                                                swout_clr, olr_clr, &
                                                swups_clr, swdns_clr,&
                                                lwups_clr, lwdns_clr   

      real, dimension (ie-is+1,je-js+1, MX_SPEC_LEVS) ::    & 
                                                swdn_trop,  &
                                                swdn_trop_clr, &
                                                swup_trop, &
                                                swup_trop_clr, &
                                                netlw_trop, &
                                                netlw_trop_clr

      real, dimension (ie-is+1,je-js+1, size(Rad_output%tdtsw,3)) ::  &
                                                tdtlw, tdtlw_clr,&
                                                hsw, hswcf, pmass

      real, dimension (ie-is+1,je-js+1, size(Rad_output%tdtsw,3)+1) :: &
                                                dfsw, ufsw,  &
                                                dfswcf, ufswcf,&
                                                flxnet, flxnetcf, &
                                                fsw, fswcf
      real, dimension (ie-is+1,je-js+1) ::      &
                                         swin_ad,     swout_ad, olr_ad,&
                              swups_ad,    swdns_ad, lwups_ad,lwdns_ad,&
                                 swin_ad_clr, swout_ad_clr, olr_ad_clr,&
                  swups_ad_clr, swdns_ad_clr, lwups_ad_clr, lwdns_ad_clr, &
                                                               heat2d
      real, dimension (ie-is+1,je-js+1, size(Rad_output%tdtsw,3)+1) :: &
                                               dfsw_ad, ufsw_ad,  &
                                               dfswcf_ad, ufswcf_ad

      integer           :: j, k
      integer           :: ipass
      logical           :: used
      integer           :: iind, jind
      integer           :: kmax
      real              :: solar_constant

!      asfc         surface albedo  [ dimensionless ]
!      asfc_vis_dir surface visible albedo  [ dimensionless ]
!      asfc_nir_dir surface nir albedo  [ dimensionless ]
!      asfc_vis_dif surface visible albedo  [ dimensionless ]
!      asfc_nir_dif surface nir albedo  [ dimensionless ]

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_diag_mod',  &
               'module has not been initialized', FATAL)

!-------------------------------------------------------------------
!    define tropopause fluxes for diagnostic use later
!-------------------------------------------------------------------
      if (Rad_control%do_sw_rad .or. Rad_control%do_lw_rad) then
         call flux_trop_calc (lat, pflux, phalf, Rad_control,  &
                              Lw_output(1), Sw_output(1), &
                              swdn_trop, swup_trop, &
                              swdn_trop_clr, swup_trop_clr, &
                              netlw_trop, netlw_trop_clr)
      endif

!-------------------------------------------------------------------
!    save special shortwave diagnostics
!-------------------------------------------------------------------
      if (Rad_control%do_sw_rad) then
        if (Rad_control%renormalize_sw_fluxes .or. all_step_diagnostics) then
          Diag_special%swdn_trop(is:ie,js:je,:) = swdn_trop(:,:,:)
          Diag_special%swup_trop(is:ie,js:je,:) = swup_trop(:,:,:)
          if (Rad_control%do_totcld_forcing) then
            Diag_special%swdn_trop_clr(is:ie,js:je,:) = swdn_trop_clr(:,:,:)
            Diag_special%swup_trop_clr(is:ie,js:je,:) = swup_trop_clr(:,:,:)
          endif
        endif
      endif
      
!---------------------------------------------------------------------
!    if sw flux renormalization is active, modify the fluxes calculated
!    on the last radiation step by the normalization factor based on
!    the difference in solar factor between the current model step and
!    the current radiation step.
!----------------------------------------------------------------------
      kmax = size (Rad_output%tdtsw,3)
      if (Rad_control%renormalize_sw_fluxes) then
        do k=1, kmax         
          hsw(:,:,k) = Sw_flux_save%hsw(is:ie,js:je,k)*flux_ratio(:,:)
        end do
        do k=1, kmax+1             
          if (do_swaerosol_forcing) then
            dfsw_ad(:,:,k) = Sw_flux_save%dfsw_ad(is:ie,js:je,k)*flux_ratio(:,:)
            ufsw_ad(:,:,k) = Sw_flux_save%ufsw_ad(is:ie,js:je,k)*flux_ratio(:,:)
          endif
          dfsw(:,:,k) = Sw_flux_save%dfsw(is:ie,js:je,k)*flux_ratio(:,:)
          ufsw(:,:,k) = Sw_flux_save%ufsw(is:ie,js:je,k)*flux_ratio(:,:)
          fsw (:,:,k) = Sw_flux_save%fsw (is:ie,js:je,k)*flux_ratio(:,:)
        end do
        do k=1,MX_SPEC_LEVS
          swdn_trop(:,:,k) = Diag_special%swdn_trop(is:ie,js:je,k)*flux_ratio(:,:)
          swup_trop(:,:,k) = Diag_special%swup_trop(is:ie,js:je,k)*flux_ratio(:,:)
        end do
        if (Rad_control%do_totcld_forcing) then
          do k=1, kmax            
            hswcf(:,:,k) = Sw_flux_save%hswcf(is:ie,js:je,k)*flux_ratio(:,:)
          end do
          do k=1, kmax+1            
            if (do_swaerosol_forcing) then
              dfswcf_ad(:,:,k) = Sw_flux_save%dfswcf_ad(is:ie,js:je,k)*flux_ratio(:,:)
              ufswcf_ad(:,:,k) = Sw_flux_save%ufswcf_ad(is:ie,js:je,k)*flux_ratio(:,:)
            endif
            dfswcf(:,:,k) = Sw_flux_save%dfswcf(is:ie,js:je,k)*flux_ratio(:,:)
            ufswcf(:,:,k) = Sw_flux_save%ufswcf(is:ie,js:je,k)*flux_ratio(:,:)
            fswcf (:,:,k) = Sw_flux_save%fswcf (is:ie,js:je,k)*flux_ratio(:,:)
          end do
          do k=1,MX_SPEC_LEVS
            swdn_trop_clr(:,:,k) = Diag_special%swdn_trop_clr(is:ie,js:je,k)*flux_ratio(:,:)
            swup_trop_clr(:,:,k) = Diag_special%swup_trop_clr(is:ie,js:je,k)*flux_ratio(:,:)
          end do
        endif

!----------------------------------------------------------------------
!    if renormalization is not active and this is a radiation step
!    (i.e., diagnostics desired), define the variables to be output as
!    the values present in Sw_output.
!---------------------------------------------------------------------
      else if (Rad_control%do_sw_rad) then
        do k=1, kmax            
          hsw(:,:,k) = Sw_output(1)%hsw(:,:,k)
        end do
        do k=1, kmax+1             
          if (do_swaerosol_forcing) then
            dfsw_ad(:,:,k) = Sw_output(indx_swaf)%dfsw(:,:,k)
            ufsw_ad(:,:,k) = Sw_output(indx_swaf)%ufsw(:,:,k)
          endif
          dfsw(:,:,k) = Sw_output(1)%dfsw(:,:,k)
          ufsw(:,:,k) = Sw_output(1)%ufsw(:,:,k)
          fsw(:,:,k) = Sw_output(1)%fsw(:,:,k)
        end do
        if (Rad_control%do_totcld_forcing) then
          do k=1, kmax             
            hswcf(:,:,k) = Sw_output(1)%hswcf(:,:,k)
          end do
          do k=1, kmax+1            
            if (do_swaerosol_forcing) then
              dfswcf_ad(:,:,k) = Sw_output(indx_swaf)%dfswcf(:,:,k)
              ufswcf_ad(:,:,k) = Sw_output(indx_swaf)%ufswcf(:,:,k)  
            endif
            dfswcf(:,:,k) = Sw_output(1)%dfswcf(:,:,k)
            ufswcf(:,:,k) = Sw_output(1)%ufswcf(:,:,k)
            fswcf(:,:,k) = Sw_output(1)%fswcf(:,:,k)
          end do
        endif

!----------------------------------------------------------------------
!    if renormalization is not active and this is not a radiation step
!    but all_step_diagnostics is activated (i.e., diagnostics desired),
!    define the variables to be output as the values previously saved
!    in the xxx_save variables.
!---------------------------------------------------------------------
      else if (all_step_diagnostics) then
        do k=1, kmax
          hsw(:,:,k) = Sw_flux_save%hsw(is:ie,js:je,k)
        end do
        do k=1, kmax+1
          if (do_swaerosol_forcing) then
            dfsw_ad(:,:,k) = Sw_flux_save%dfsw_ad(is:ie,js:je,k)
            ufsw_ad(:,:,k) = Sw_flux_save%ufsw_ad(is:ie,js:je,k)
          endif
          dfsw(:,:,k) = Sw_flux_save%dfsw(is:ie,js:je,k)
          ufsw(:,:,k) = Sw_flux_save%ufsw(is:ie,js:je,k)
          fsw (:,:,k) = Sw_flux_save%fsw (is:ie,js:je,k)
        end do
        swdn_trop(:,:,:) = Diag_special%swdn_trop(is:ie,js:je,:)
        swup_trop(:,:,:) = Diag_special%swup_trop(is:ie,js:je,:)
        if (Rad_control%do_totcld_forcing) then
          do k=1, kmax
            hswcf(:,:,k) = Sw_flux_save%hswcf(is:ie,js:je,k)
          end do
          do k=1, kmax+1
            if (do_swaerosol_forcing) then
              dfswcf_ad(:,:,k) = Sw_flux_save%dfswcf_ad(is:ie,js:je,k)
              ufswcf_ad(:,:,k) = Sw_flux_save%ufswcf_ad(is:ie,js:je,k)
            endif
            dfswcf(:,:,k) = Sw_flux_save%dfswcf(is:ie,js:je,k)
            ufswcf(:,:,k) = Sw_flux_save%ufswcf(is:ie,js:je,k)
            fswcf (:,:,k) = Sw_flux_save%fswcf (is:ie,js:je,k)
          end do
          swdn_trop_clr(:,:,:) = Diag_special%swdn_trop_clr(is:ie,js:je,:)
          swup_trop_clr(:,:,:) = Diag_special%swup_trop_clr(is:ie,js:je,:)
        endif
      endif

!---------------------------------------------------------------------
!    define the sw diagnostic arrays.
!---------------------------------------------------------------------
      if (Rad_control%renormalize_sw_fluxes .or. Rad_control%do_sw_rad .or.    &
          all_step_diagnostics) then
        if (do_swaerosol_forcing) then
          swin_ad (:,:) = dfsw_ad(:,:,1)
          swout_ad(:,:) = ufsw_ad(:,:,1)
          swups_ad(:,:) = ufsw_ad(:,:,kmax+1)
          swdns_ad(:,:) = dfsw_ad(:,:,kmax+1)
        endif
        swin (:,:) = dfsw(:,:,1)
        swout(:,:) = ufsw(:,:,1)
        swups(:,:) = ufsw(:,:,kmax+1)
        swdns(:,:) = dfsw(:,:,kmax+1)
        if (Rad_control%do_totcld_forcing) then
          if (do_swaerosol_forcing) then
            swin_ad_clr (:,:) = dfswcf_ad(:,:,1)
            swout_ad_clr(:,:) = ufswcf_ad(:,:,1)
            swups_ad_clr(:,:) = ufswcf_ad(:,:,kmax+1)
            swdns_ad_clr(:,:) = dfswcf_ad(:,:,kmax+1)
          endif
          swin_clr (:,:) = dfswcf(:,:,1)
          swout_clr(:,:) = ufswcf(:,:,1)
          swups_clr(:,:) = ufswcf(:,:,kmax+1)
          swdns_clr(:,:) = dfswcf(:,:,kmax+1)
        endif

        if (id_alb_sfc_avg > 0) then
          swups_acc(is:ie,js:je) = swups_acc(is:ie, js:je) + swups(:,:)
          swdns_acc(is:ie,js:je) = swdns_acc(is:ie, js:je) + swdns(:,:)
        endif
 
!---------------------------------------------------------------------
!   send standard sw diagnostics to diag_manager.
!---------------------------------------------------------------------
      if (Time_diag > Time) then
        if (Rad_control%do_totcld_forcing) then
          ipass = 2
        else
          ipass = 1
        endif

!------- sw tendency -----------
        if (id_tdt_sw(ipass) > 0 ) then
          used = send_data (id_tdt_sw(ipass),    &
                            Rad_output%tdtsw(is:ie,js:je,:),   &
                            Time_diag, is, js, 1)
        endif

!---- 3d upward sw flux ---------
        if (id_ufsw(ipass) > 0 ) then
          used = send_data (id_ufsw(ipass),    &
                            Rad_output%ufsw(is:ie,js:je,:),   &
                            Time_diag, is, js, 1)
        endif

!---- 3d downward sw flux ---------
        if (id_dfsw(ipass) > 0 ) then
          used = send_data (id_dfsw(ipass),    &
                            Rad_output%dfsw(is:ie,js:je,:),   &
                            Time_diag, is, js, 1)
        endif


!------- incoming sw flux toa -------
        if (id_swdn_toa(ipass) > 0 ) then
          used = send_data (id_swdn_toa(ipass), swin,   &
                            Time_diag, is, js )
        endif

!------- outgoing sw flux toa -------
        if (id_swup_toa(ipass) > 0 ) then
          used = send_data (id_swup_toa(ipass), swout,    &
                            Time_diag, is, js )
        endif

        do k = 1, MX_SPEC_LEVS
!------- incoming sw flux trop -------
          if (id_swdn_special(k,ipass) > 0 ) then
            used = send_data (id_swdn_special(k,ipass), swdn_trop(:,:,k), Time_diag, is, js )
          endif

!------- outgoing sw flux trop -------
          if (id_swdn_special(k,ipass) > 0 ) then
            used = send_data (id_swup_special(k,ipass), swup_trop(:,:,k), Time_diag, is, js )
          endif
        end do

!------- upward sw flux surface -------
        if (id_swup_sfc(ipass) > 0 ) then
          used = send_data (id_swup_sfc(ipass), swups,    &
                            Time_diag, is, js )
        endif

!------- downward sw flux surface -------
        if (id_swdn_sfc(ipass) > 0 ) then
          used = send_data (id_swdn_sfc(ipass), swdns,   &
                            Time_diag, is, js )
        endif
        
!------- net sw flux at toa -------
       if (id_swtoa(ipass) > 0 ) then
          used = send_data (id_swtoa(ipass), swin-swout,   &
                             Time_diag, is, js )
         endif

!------- net sw flux at surface -------
         if (id_swsfc(ipass) > 0 ) then
           used = send_data (id_swsfc(ipass), swdns-swups,   &
                             Time_diag, is, js )
         endif

      if (do_swaerosol_forcing) then

!------- net sw flux at toa -------
         if (id_swtoa_ad(ipass) > 0 ) then
           used = send_data (id_swtoa_ad(ipass), swin_ad-swout_ad,   &
                             Time_diag, is, js )
        endif
 
!------- net sw flux at surface -------
         if (id_swsfc_ad(ipass) > 0 ) then
           used = send_data (id_swsfc_ad(ipass), swdns_ad-swups_ad,   &
                           Time_diag, is, js )
         endif
 
!------- sw flux down at surface -------
         if (id_swdn_sfc_ad(ipass) > 0 ) then
           used = send_data (id_swdn_sfc_ad(ipass), swdns_ad,   &
                             Time_diag, is, js )
         endif

!------- sw flux up at surface -------
         if (id_swup_sfc_ad(ipass) > 0 ) then
           used = send_data (id_swup_sfc_ad(ipass), swups_ad,   &
                            Time_diag, is, js )
         endif

!------- outgoing sw flux toa -------
         if (id_swup_toa_ad(ipass) > 0 ) then
           used = send_data (id_swup_toa_ad(ipass), swout_ad,    &
                            Time_diag, is, js )
        endif
     endif

!-----------------------------------------------
!-------  cmip diagnostics (full-sky)  ---------
!-----------------------------------------------

!------- sw tendency -----------
        if (query_cmip_diag_id(ID_tntrs)) then
          used = send_cmip_data_3d (ID_tntrs, Rad_output%tdtsw(is:ie,js:je,:),  &
                            Time_diag, is, js, 1)
        endif

!------- 3d upward sw flux -------
        if (query_cmip_diag_id(ID_rsu)) then
          used = send_cmip_data_3d (ID_rsu, ufsw(is:ie,js:je,:), Time_diag, is, js, 1)
        endif

!------- 3d downward sw flux -------
        if (query_cmip_diag_id(ID_rsd)) then
          used = send_cmip_data_3d (ID_rsd, dfsw(is:ie,js:je,:), Time_diag, is, js, 1)
        endif

!------- downward sw flux surface -------
        if (id_rsds > 0 ) then
          used = send_data (id_rsds, swdns, Time_diag, is, js )
        endif

!------- upward sw flux surface -------
        if (id_rsus > 0 ) then
          used = send_data (id_rsus, swups, Time_diag, is, js )
        endif

!------- incoming sw flux toa -------
        if (id_rsdt > 0 ) then
          used = send_data (id_rsdt, swin, Time_diag, is, js )
        endif

!------- outgoing sw flux toa -------
        if (id_rsut > 0 ) then
          used = send_data (id_rsut, swout, Time_diag, is, js )
        endif

        if (do_swaerosol_forcing) then

!------- downward sw flux surface (without aerosols) -------
          if (id_rsdsaf > 0 ) then
            used = send_data (id_rsdsaf, swdns_ad, Time_diag, is, js )
          endif

!------- upward sw flux surface (without aerosols) -------
          if (id_rsusaf > 0 ) then
            used = send_data (id_rsusaf, swups_ad, Time_diag, is, js )
          endif

!------- outgoing sw flux toa -------
          if (id_rsutaf > 0 ) then
            used = send_data (id_rsutaf, swout_ad, Time_diag, is, js )
          endif

!------- 3d upward sw flux -------
          if (query_cmip_diag_id(ID_rsuaf)) then
            used = send_cmip_data_3d (ID_rsuaf, ufsw_ad(is:ie,js:je,:), Time_diag, is, js, 1)
          endif

!------- 3d downward sw flux -------
          if (query_cmip_diag_id(ID_rsdaf)) then
            used = send_cmip_data_3d (ID_rsdaf, ufsw_ad(is:ie,js:je,:), Time_diag, is, js, 1)
          endif

        endif

!----------------------------------------------------------------------
!    now pass clear-sky diagnostics, if they have been calculated.
!----------------------------------------------------------------------
        if (Rad_control%do_totcld_forcing) then
          ipass = 1

!------- sw tendency -----------
          if (id_tdt_sw(ipass) > 0 ) then
            used = send_data (id_tdt_sw(ipass),   &
                              Rad_output%tdtsw_clr(is:ie,js:je,:),  &
                              Time_diag, is, js, 1)
          endif

!---- 3d upward sw flux ---------
         if (id_ufsw(ipass) > 0 ) then
           used = send_data (id_ufsw(ipass), Rad_output%ufsw_clr(is:ie,js:je,:), &
                             Time_diag, is, js, 1)
         endif
 
!---- 3d downward sw flux ---------
         if (id_dfsw(ipass) > 0 ) then
           used = send_data (id_dfsw(ipass), Rad_output%dfsw_clr(is:ie,js:je,:), &
                             Time_diag, is, js, 1)
         endif

!------- incoming sw flux toa -------
          if (id_swdn_toa(ipass) > 0 ) then
            used = send_data (id_swdn_toa(ipass), swin_clr,    &
                              Time_diag, is, js )
          endif

!------- outgoing sw flux toa -------
          if (id_swup_toa(ipass) > 0 ) then
            used = send_data (id_swup_toa(ipass), swout_clr,  &
                              Time_diag, is, js )
          endif

          do k = 1, MX_SPEC_LEVS
!------- incoming sw flux trop -------
            if (id_swdn_special(1,ipass) > 0 ) then
              used = send_data (id_swdn_special(k,ipass), swdn_trop_clr(:,:,k), &
                                Time_diag, is, js )
            endif

!------- outgoing sw flux trop -------
            if (id_swup_special(1,ipass) > 0 ) then
              used = send_data (id_swup_special(k,ipass), swup_trop_clr(:,:,k), &
                                Time_diag, is, js )
            endif
          end do

!------- upward sw flux surface -------
          if (id_swup_sfc(ipass) > 0 ) then
            used = send_data (id_swup_sfc(ipass), swups_clr,   &
                              Time_diag, is, js )
          endif

!------- downward sw flux surface -------
          if (id_swdn_sfc(ipass) > 0 ) then
            used = send_data (id_swdn_sfc(ipass), swdns_clr,    &
                              Time_diag, is, js )
          endif

!------- net sw flux at toa -------
        if (id_swtoa(ipass) > 0 ) then
           used = send_data (id_swtoa(ipass), swin_clr-swout_clr,   &
                             Time_diag, is, js )
         endif
 
!------- net sw flux at surface -------
         if (id_swsfc(ipass) > 0 ) then
           used = send_data (id_swsfc(ipass), swdns_clr-swups_clr,   &
                            Time_diag, is, js )
         endif
     if (do_swaerosol_forcing) then

!------- net sw flux at toa -------
        if (id_swtoa_ad(ipass) > 0 ) then
          used = send_data (id_swtoa_ad(ipass), swin_ad_clr-swout_ad_clr,   &
                            Time_diag, is, js )
        endif
 
!------- net sw flux at surface -------
         if (id_swsfc_ad(ipass) > 0 ) then
           used = send_data (id_swsfc_ad(ipass), swdns_ad_clr-swups_ad_clr,   &
                             Time_diag, is, js )
         endif
 
!------- sw flux down at surface -------
         if (id_swdn_sfc_ad(ipass) > 0 ) then
           used = send_data (id_swdn_sfc_ad(ipass), swdns_ad_clr,   &
                             Time_diag, is, js )
        endif
 
!------- sw flux up at surface -------
        if (id_swup_sfc_ad(ipass) > 0 ) then
          used = send_data (id_swup_sfc_ad(ipass), swups_ad_clr,   &
                            Time_diag, is, js )
        endif
 
!------- outgoing sw flux toa -------
        if (id_swup_toa_ad(ipass) > 0 ) then
           used = send_data (id_swup_toa_ad(ipass), swout_ad_clr,    &
                            Time_diag, is, js )
        endif

     endif

!-----------------------------------------------
!------  cmip diagnostics (clear-sky)  ---------
!-----------------------------------------------

!------- sw tendency -----------
          if (query_cmip_diag_id(ID_tntrscs)) then
            used = send_cmip_data_3d (ID_tntrscs, Rad_output%tdtsw_clr(is:ie,js:je,:),  &
                              Time_diag, is, js, 1)
          endif

!------- 3d upward sw flux -------
        if (query_cmip_diag_id(ID_rsucs)) then
          used = send_cmip_data_3d (ID_rsucs, ufswcf(is:ie,js:je,:), Time_diag, is, js, 1)
        endif

!------- 3d downward sw flux -------
        if (query_cmip_diag_id(ID_rsdcs)) then
          used = send_cmip_data_3d (ID_rsdcs, dfswcf(is:ie,js:je,:), Time_diag, is, js, 1)
        endif

!------- downward sw flux surface -------
          if (id_rsdscs > 0 ) then
            used = send_data (id_rsdscs, swdns_clr, Time_diag, is, js )
          endif

!------- upward sw flux surface -------
          if (id_rsuscs > 0 ) then
            used = send_data (id_rsuscs, swups_clr, Time_diag, is, js )
          endif

!------- outgoing sw flux toa -------
          if (id_rsutcs > 0 ) then
            used = send_data (id_rsutcs, swout_clr, Time_diag, is, js )
          endif

          if (do_swaerosol_forcing) then

!------- downward sw flux surface (clean-clear-sky) -------
            if (id_rsdscsaf > 0 ) then
              used = send_data (id_rsdscsaf, swdns_ad_clr, Time_diag, is, js )
            endif

!------- upward sw flux surface (clean-clear-sky) -------
            if (id_rsuscsaf > 0 ) then
              used = send_data (id_rsuscsaf, swups_ad_clr, Time_diag, is, js )
            endif

!------- outgoing sw flux toa (clean-clear-sky) -------
            if (id_rsutcsaf > 0 ) then
              used = send_data (id_rsutcsaf, swout_ad_clr, Time_diag, is, js )
            endif

!------- 3d upward sw flux ---------
            if (query_cmip_diag_id(ID_rsucsaf)) then
              used = send_cmip_data_3d (ID_rsucsaf, ufswcf_ad(is:ie,js:je,:), Time_diag, is, js, 1)
            endif

!------- 3d downward sw flux ---------
            if (query_cmip_diag_id(ID_rsdcsaf)) then
              used = send_cmip_data_3d (ID_rsdcsaf, dfswcf_ad(is:ie,js:je,:), Time_diag, is, js, 1)
            endif

          endif

         endif  ! (Rad_control%do_totcld_forcing)

!-----------------------------------------------------------------------
!    send cloud-forcing-independent diagnostics to diagnostics manager.
!-----------------------------------------------------------------------

!---- 3d total radiative heating ---------
        if (id_allradp > 0 ) then
          used = send_data (id_allradp    ,    &
                            Rad_output%tdt_rad(is:ie,js:je,:),   &
                            Time_diag, is, js, 1)
        endif

        if (query_cmip_diag_id(ID_tntr)) then
          used = send_cmip_data_3d (ID_tntr, Rad_output%tdt_rad(is:ie,js:je,:),  &
                            Time_diag, is, js, 1)
        endif

! integrated radiative heating (can be used to check energy conservation)
! NOTE: phalf is now pass directly from radiation_driver argument list
!       so this diagnostic can be computed every timestep
  !BW if (Rad_control%do_sw_rad .or. Rad_control%do_lw_rad) then
        if (id_heat2d > 0 .or. id_heat2d_sw > 0) then
          do k=1,kmax
            pmass(:,:,k) = phalf(:,:,k+1)-phalf(:,:,k)
          enddo
        endif
        if (id_heat2d > 0) then
          heat2d = CP_AIR/GRAV * sum(Rad_output%tdt_rad(is:ie,js:je,:)*pmass,3)  
          used = send_data (id_heat2d, heat2d, Time_diag, is, js )
        endif
        if (id_heat2d_sw > 0) then
          heat2d = CP_AIR/GRAV * sum(Rad_output%tdtsw(is:ie,js:je,:)*pmass,3)  
          used = send_data (id_heat2d_sw, heat2d, Time_diag, is, js )
        endif
  !BW endif

!------- conc_drop  -------------------------
!BW       if (Rad_control%do_sw_rad .or. Rad_control%do_lw_rad) then
!BW         if ( id_conc_drop > 0 ) then
!BW           used = send_data (id_conc_drop, Lsc_microphys%conc_drop, &
!BW                             Time_diag, is, js, 1)
!BW         endif
!BW       endif
  
!------- conc_ice  -------------------------
!BW       if (Rad_control%do_sw_rad .or. Rad_control%do_lw_rad) then
!BW         if (id_conc_ice > 0 ) then
!BW           used = send_data (id_conc_ice, Lsc_microphys%conc_ice, &
!BW                             Time_diag, is, js, 1)
!BW         endif
!BW       endif

!------- solar constant  -------------------------
        if (Rad_control%do_sw_rad .or. Rad_control%do_lw_rad) then
          if ( id_solar_constant > 0 ) then
            call get_solar_constant (solar_constant)
            used = send_data ( id_solar_constant, solar_constant, Time_diag )
          endif
        endif

!------- co2 mixing ratio used for tf calculation  -------------------
        if (Rad_control%do_sw_rad .or. Rad_control%do_lw_rad) then
          if ( id_co2_tf > 0 ) then
            used = send_data ( id_co2_tf,   &
                               1.0E6*Rad_gases%co2_for_last_tf_calc,  &
                               Time_diag )
          endif
        endif
 
!------- ch4 mixing ratio used for tf calculation   ---------------
        if (Rad_control%do_sw_rad .or. Rad_control%do_lw_rad) then
          if ( id_ch4_tf > 0 ) then
            used = send_data ( id_ch4_tf,  &
                               1.0E9*Rad_gases%ch4_for_last_tf_calc,  &
                               Time_diag )
          endif
        endif
 
!------- n2o mixing ratio used for tf calculation  ---------------
        if (Rad_control%do_sw_rad .or. Rad_control%do_lw_rad) then
          if ( id_n2o_tf > 0 ) then
            used = send_data ( id_n2o_tf,   &
                               1.0E9*Rad_gases%n2o_for_last_tf_calc,  &
                               Time_diag )
          endif
        endif
 
!------- co2 mixing ratio  -------------------------
        if (Rad_control%do_sw_rad .or. Rad_control%do_lw_rad) then
          if ( id_rrvco2 > 0 ) then
            used = send_data ( id_rrvco2, 1.0E6*Rad_gases%rrvco2,  &
                               Time_diag )
          endif
          if ( id_co2mass > 0 ) then
            used = send_data ( id_co2mass, atm_mass*WTMCO2/(WTMAIR*GRAV)* &
                                           Rad_gases%rrvco2,  &
                               Time_diag )
          endif
        endif
 
!------- f11 mixing ratio  -------------------------
        if (Rad_control%do_sw_rad .or. Rad_control%do_lw_rad) then
          if ( id_rrvf11 > 0 ) then
            used = send_data ( id_rrvf11, 1.0E12*Rad_gases%rrvf11,  &
                               Time_diag )
          endif
          if ( id_cfc11global > 0 ) then
            used = send_data ( id_cfc11global, 1.0E12*Rad_gases%rrvf11,  &
                               Time_diag )
          endif
        endif
 
!------- f12 mixing ratio  -------------------------
        if (Rad_control%do_sw_rad .or. Rad_control%do_lw_rad) then
          if ( id_rrvf12 > 0 ) then
            used = send_data ( id_rrvf12, 1.0E12*Rad_gases%rrvf12,  &
                               Time_diag )
          endif
          if ( id_cfc12global > 0 ) then
            used = send_data ( id_cfc12global, 1.0E12*Rad_gases%rrvf12,  &
                               Time_diag )
          endif
        endif
 
!------- f113 mixing ratio  -------------------------
        if (Rad_control%do_sw_rad .or. Rad_control%do_lw_rad) then
          if ( id_rrvf113 > 0 ) then
            used = send_data ( id_rrvf113, 1.0E12*Rad_gases%rrvf113,  &
                               Time_diag )
          endif
          if ( id_cfc113global > 0 ) then
            used = send_data ( id_cfc113global, 1.0E12*Rad_gases%rrvf113,  &
                               Time_diag )
          endif
        endif

!------- f22 mixing ratio  -------------------------
        if (Rad_control%do_sw_rad .or. Rad_control%do_lw_rad) then
          if ( id_rrvf22 > 0 ) then
            used = send_data ( id_rrvf22, 1.0E12*Rad_gases%rrvf22,  &
                               Time_diag )
          endif
          if ( id_hcfc22global > 0 ) then
            used = send_data ( id_hcfc22global, 1.0E12*Rad_gases%rrvf22,  &
                               Time_diag )
          endif
        endif

!------- ch4 mixing ratio  -------------------------
        if (Rad_control%do_sw_rad .or. Rad_control%do_lw_rad) then
          if ( id_rrvch4 > 0 ) then
            used = send_data ( id_rrvch4, 1.0E9*Rad_gases%rrvch4,  &
                               Time_diag )
          endif
          if ( id_ch4global > 0 ) then
            used = send_data ( id_ch4global, 1.0E9*Rad_gases%rrvch4,  &
                               Time_diag )
          endif
        endif

!------- n2o mixing ratio  -------------------------
        if (Rad_control%do_sw_rad .or. Rad_control%do_lw_rad) then
          if ( id_rrvn2o > 0 ) then
            used = send_data ( id_rrvn2o, 1.0E9*Rad_gases%rrvn2o,  &
                               Time_diag )
          endif
          if ( id_n2oglobal > 0 ) then
            used = send_data ( id_n2oglobal, 1.0E9*Rad_gases%rrvn2o,  &
                               Time_diag )
          endif
        endif

!------- surface albedo  -------------------------
        if ( id_alb_sfc_avg > 0 ) then
          used = send_data ( id_alb_sfc_avg,  &
                  100.*swups_acc(is:ie,js:je)/swdns_acc(is:ie,js:je), &
                     Time_diag, is, js )
        endif
        if ( id_alb_sfc > 0 ) then
          used = send_data ( id_alb_sfc, 100.*swups/(1.0e-35 + swdns), &
                     Time_diag, is, js )
        endif

!------- surface visible albedo  -------------------------
        if ( id_alb_sfc_vis_dir > 0 ) then
          used = send_data ( id_alb_sfc_vis_dir, &
                         100.*asfc_vis_dir, Time_diag, is, js )
        endif
        if ( id_alb_sfc_vis_dif > 0 ) then
          used = send_data ( id_alb_sfc_vis_dif, &
                         100.*asfc_vis_dif, Time_diag, is, js )
        endif
 
!------- surface nir albedo  -------------------------
        if ( id_alb_sfc_nir_dir > 0 ) then
          used = send_data ( id_alb_sfc_nir_dir, &
                         100.*asfc_nir_dir, Time_diag, is, js )
        endif
        if ( id_alb_sfc_nir_dif > 0 ) then
           used = send_data ( id_alb_sfc_nir_dif, &
                         100.*asfc_nir_dif, Time_diag, is, js )
        endif
 
!------- surface net sw flux, direct and diffuse  --------------------
        if ( id_flux_sw > 0 ) then
         used = send_data ( id_flux_sw, &
                (Rad_output%flux_sw_surf_dir(is:ie,js:je) + &
                 Rad_output%flux_sw_surf_dif(is:ie,js:je)), &
                                          Time_diag, is, js )
        endif

        if ( id_flux_sw_dir > 0 ) then
         used = send_data ( id_flux_sw_dir, &
          Rad_output%flux_sw_surf_dir( is:ie,js:je), Time_diag,  &
                                                              is, js )
        endif
        
        if ( id_flux_sw_refl_dir > 0 ) then
         used = send_data ( id_flux_sw_refl_dir, &
          Rad_output%flux_sw_surf_refl_dir( is:ie,js:je), Time_diag,  &
                                                              is, js )
        endif

        if ( id_flux_sw_refl_vis_dir > 0 ) then
         used = send_data ( id_flux_sw_refl_vis_dir, &
          Rad_output%flux_sw_refl_vis_dir( is:ie,js:je), Time_diag,  &
                                                              is, js )
        endif
        if ( id_flux_sw_dif > 0 ) then
          used = send_data ( id_flux_sw_dif, &
           Rad_output%flux_sw_surf_dif(is:ie,js:je), Time_diag, &
                                                              is, js )
        endif

!------- surface downward visible sw flux, direct and diffuse ----------
        if ( id_flux_sw_down_vis_dir > 0 ) then
          used = send_data ( id_flux_sw_down_vis_dir, &
                     Rad_output%flux_sw_down_vis_dir(is:ie,js:je), &
                     Time_diag, is, js )
        endif
        if ( id_flux_sw_down_vis_dif > 0 ) then
          used = send_data ( id_flux_sw_down_vis_dif, &
                     Rad_output%flux_sw_down_vis_dif(is:ie, js:je), &
                     Time_diag, is, js )
        endif
 
!------- surface downward total sw flux, direct and diffuse  ----------
        if ( id_flux_sw_down_total_dir > 0 ) then
          used = send_data ( id_flux_sw_down_total_dir,  &
                     Rad_output%flux_sw_down_total_dir(is:ie,js:je),  &
                     Time_diag, is, js )
        endif
        if ( id_flux_sw_down_total_dif > 0 ) then
         used = send_data ( id_flux_sw_down_total_dif,  &
                    Rad_output%flux_sw_down_total_dif(is:ie,js:je),  &
                    Time_diag, is, js )
        endif
        if ( id_rsdsdiff > 0 ) then
         used = send_data ( id_rsdsdiff,  &
                    Rad_output%flux_sw_down_total_dif(is:ie,js:je),  &
                    Time_diag, is, js )
        endif

      if (Rad_control%do_totcld_forcing) then
 
!------- surface downward total sw flux, direct and diffuse  ----------
        if ( id_flux_sw_down_total_dir_clr > 0 ) then
          used = send_data ( id_flux_sw_down_total_dir_clr,  &
                 Rad_output%flux_sw_down_total_dir_clr(is:ie,js:je),&
                 Time_diag, is, js )
        endif
        if ( id_flux_sw_down_total_dif_clr > 0 ) then
          used = send_data ( id_flux_sw_down_total_dif_clr,  &
                  Rad_output%flux_sw_down_total_dif_clr(is:ie,js:je),  &
                  Time_diag, is, js )
        endif
        if ( id_rsdscsdiff > 0 ) then
         used = send_data ( id_rsdscsdiff,  &
                    Rad_output%flux_sw_down_total_dif_clr(is:ie,js:je),  &
                    Time_diag, is, js )
        endif
        if ( id_flux_sw_down_vis_clr > 0 ) then
          used = send_data ( id_flux_sw_down_vis_clr, &
                     Rad_output%flux_sw_down_vis_clr(is:ie, js:je), &
                     Time_diag, is, js )
        endif
      endif
 
!------- surface net visible sw flux, total, direct and diffuse -------
        if ( id_flux_sw_vis > 0 ) then
          used = send_data ( id_flux_sw_vis,   &
                Rad_output%flux_sw_vis(is:ie,js:je), Time_diag, is, js )
        endif
        if ( id_flux_sw_vis_dir > 0 ) then
          used = send_data ( id_flux_sw_vis_dir,   &
            Rad_output%flux_sw_vis_dir(is:ie,js:je), Time_diag, is, js )
        endif
        if ( id_flux_sw_vis_dif > 0 ) then
          used = send_data ( id_flux_sw_vis_dif,  &
            Rad_output%flux_sw_vis_dif(is:ie,js:je), Time_diag, is, js )
        endif

! use the correct astrononmy values (NEED TO CHECK THIS)
! could do this with pointer

        if (Rad_control%renormalize_sw_fluxes) then
!------- cosine of zenith angle ----------------
            if ( id_cosz > 0 ) then
              used = send_data ( id_cosz, Astro_phys%cosz, Time_diag, is, js )
            endif

!------- daylight fraction  --------------
            if ( id_fracday > 0 ) then
              used = send_data (id_fracday, Astro_phys%fracday, Time_diag, is, js )
            end if
         else
            ! only on radiation steps (not available for all_step_diagnostics)
            if (Rad_control%do_sw_rad) then
!------- cosine of zenith angle ----------------
               if ( id_cosz > 0 ) then
                 used = send_data ( id_cosz, Astro%cosz, Time_diag, is, js )
               endif

!------- daylight fraction  --------------
               if ( id_fracday > 0 ) then
                 used = send_data (id_fracday, Astro%fracday, Time_diag, is, js )
               end if
            endif
         endif


      endif
      endif   ! (renormalize_sw_fluxes .or. do_rad .or.   
              !  all_step_diagnostics)

!---------------------------------------------------------------------
!    define the longwave diagnostic arrays for the sea-esf radiation 
!    package.  convert to mks units.
!---------------------------------------------------------------------
      if (Rad_control%do_lw_rad) then
        olr  (:,:)   = Lw_output(1)%flxnet(:,:,1)
        lwups(:,:)   =   STEFAN*ts(:,:  )**4
        lwdns(:,:)   = lwups(:,:) - Lw_output(1)%flxnet(:,:,kmax+1)
        tdtlw(:,:,:) = Lw_output(1)%heatra(:,:,:)/ SECONDS_PER_DAY
        flxnet(:,:,:) = Lw_output(1)%flxnet(:,:,:)
        if (do_lwaerosol_forcing) then
          olr_ad  (:,:)   = Lw_output(indx_lwaf)%flxnet(:,:,1)
          lwups_ad(:,:)   = STEFAN*ts(:,:  )**4
          lwdns_ad(:,:)   = lwups_ad(:,:) -    &
                                Lw_output(indx_lwaf)%flxnet(:,:,kmax+1)
        endif

        if (Rad_control%do_totcld_forcing) then
          olr_clr  (:,:)   = Lw_output(1)%flxnetcf(:,:,1)
          lwups_clr(:,:)   =              STEFAN*ts(:,:  )**4
          lwdns_clr(:,:)   = lwups_clr(:,:) -    & 
                             Lw_output(1)%flxnetcf(:,:,kmax+1)
          tdtlw_clr(:,:,:) = Lw_output(1)%heatracf(:,:,:)/SECONDS_PER_DAY
          flxnetcf(:,:,:) = Lw_output(1)%flxnetcf(:,:,:)
          if (do_lwaerosol_forcing) then
            olr_ad_clr  (:,:)   = Lw_output(indx_lwaf)%flxnetcf(:,:,1)
            lwups_ad_clr(:,:)   = STEFAN*ts(:,:  )**4
            lwdns_ad_clr(:,:)   = lwups_ad_clr(:,:) -    &
                             Lw_output(indx_lwaf)%flxnetcf(:,:,kmax+1)
          endif
        endif

!---------------------------------------------------------------------
!    if diagnostics are desired on all physics steps, save the arrays 
!    for later use.
!---------------------------------------------------------------------
        if (all_step_diagnostics) then
          if (do_lwaerosol_forcing) then
            Lw_flux_save%olr_ad  (is:ie,js:je)   = olr_ad(:,:)
            Lw_flux_save%lwups_ad(is:ie,js:je)   = lwups_ad(:,:)
            Lw_flux_save%lwdns_ad(is:ie,js:je)   = lwdns_ad(:,:)
          endif
          Lw_flux_save%olr  (is:ie,js:je)   = olr(:,:)
          Lw_flux_save%lwups(is:ie,js:je)   = lwups(:,:)
          Lw_flux_save%lwdns(is:ie,js:je)   = lwdns(:,:)
          Lw_flux_save%tdtlw(is:ie,js:je,:) = tdtlw(:,:,:)
          Lw_flux_save%flxnet(is:ie,js:je,:) = Lw_output(1)%flxnet(:,:,:)
          Diag_special%netlw_trop(is:ie,js:je,:) = netlw_trop(:,:,:)

          if (Rad_control%do_totcld_forcing) then
            if (do_lwaerosol_forcing) then 
              Lw_flux_save%olr_ad_clr  (is:ie,js:je)   = olr_ad_clr(:,:)
              Lw_flux_save%lwups_ad_clr(is:ie,js:je)   = lwups_ad_clr(:,:)
              Lw_flux_save%lwdns_ad_clr(is:ie,js:je)   = lwdns_ad_clr(:,:)
            endif
            Lw_flux_save%olr_clr  (is:ie,js:je)   = olr_clr(:,:)
            Lw_flux_save%flxnetcf(is:ie,js:je,:)  = Lw_output(1)%flxnetcf(:,:,:)
            Lw_flux_save%lwups_clr(is:ie,js:je)   = lwups_clr(:,:)
            Lw_flux_save%lwdns_clr(is:ie,js:je)   = lwdns_clr(:,:)
            Lw_flux_save%tdtlw_clr(is:ie,js:je,:) = tdtlw_clr(:,:,:)
            Diag_special%netlw_trop_clr(is:ie,js:je,:) = netlw_trop_clr(:,:,:)
          endif
        endif

!---------------------------------------------------------------------
!    if this is not a radiation step, but diagnostics are desired,
!    define the fields from the xxx_save variables.
!---------------------------------------------------------------------
      else if (all_step_diagnostics) then  ! (do_lw_rad)
        if (do_lwaerosol_forcing) then 
           olr_ad(:,:)     = Lw_flux_save%olr_ad  (is:ie,js:je)
           lwups_ad(:,:)   = Lw_flux_save%lwups_ad(is:ie,js:je)
           lwdns_ad(:,:)   = Lw_flux_save%lwdns_ad(is:ie,js:je) 
        endif
        olr(:,:)     = Lw_flux_save%olr  (is:ie,js:je)
        lwups(:,:)   = Lw_flux_save%lwups(is:ie,js:je)
        lwdns(:,:)   = Lw_flux_save%lwdns(is:ie,js:je)
        tdtlw(:,:,:) = Lw_flux_save%tdtlw(is:ie,js:je,:)
        flxnet(:,:,:) = Lw_flux_save%flxnet(is:ie,js:je,:)
        netlw_trop(:,:,:) = Diag_special%netlw_trop(is:ie,js:je,:)

        if (Rad_control%do_totcld_forcing) then
          if (do_lwaerosol_forcing) then
            olr_ad_clr(:,:)     = Lw_flux_save%olr_ad_clr  (is:ie,js:je)
            lwups_ad_clr(:,:)   = Lw_flux_save%lwups_ad_clr(is:ie,js:je)
            lwdns_ad_clr(:,:)   = Lw_flux_save%lwdns_ad_clr(is:ie,js:je)
          endif
          olr_clr(:,:)     = Lw_flux_save%olr_clr  (is:ie,js:je)
          lwups_clr(:,:)   = Lw_flux_save%lwups_clr(is:ie,js:je)
          lwdns_clr(:,:)   = Lw_flux_save%lwdns_clr(is:ie,js:je)
          tdtlw_clr(:,:,:) = Lw_flux_save%tdtlw_clr(is:ie,js:je,:)
          flxnetcf (:,:,:) = Lw_flux_save%flxnetcf(is:ie,js:je,:)
          netlw_trop_clr(:,:,:) =   &
                             Diag_special%netlw_trop_clr(is:ie,js:je,:)
        endif
      endif  ! (all_step_diagnostics)


      if (Rad_control%do_lw_rad .or. all_step_diagnostics) then
      if (Time_diag > Time) then
!---------------------------------------------------------------------
!   send standard lw diagnostics to diag_manager.
!---------------------------------------------------------------------
        if (Rad_control%do_totcld_forcing) then
          ipass = 2
        else
          ipass = 1
        endif

!---- net lw flux ---------
        if (id_flxnet(ipass) > 0 ) then
          used = send_data (id_flxnet(ipass),    &
                           flxnet,   &
                           Time_diag, is, js, 1)
        endif

!------- lw tendency -----------
        if (id_tdt_lw(ipass) > 0 ) then
          used = send_data (id_tdt_lw(ipass), tdtlw,    &
                            Time_diag, is, js, 1)
        endif

!------- outgoing lw flux toa (olr) -------
        if (id_olr(ipass) > 0 ) then
          used = send_data (id_olr(ipass), olr,    &
                            Time_diag, is, js )
        endif

!------- net radiation (lw + sw) at toa -------
        if (id_netrad_toa(ipass) > 0 ) then
          used = send_data (id_netrad_toa(ipass),   &
                            swin - swout - olr, &
                            Time_diag, is, js )
        endif

!------- net radiation (lw + sw) at 1 Pa-------
        if (id_netrad_1_Pa(ipass) > 0 ) then
          used = send_data (id_netrad_1_Pa(ipass),   &
               swdn_trop(:,:,4) -swup_trop(:,:,4) -netlw_trop(:,:,4), &
                             Time_diag, is, js )
        endif

        do k = 1, MX_SPEC_LEVS
!------- net lw flux trop (netlw_trop) -------
          if (id_netlw_special(k,ipass) > 0 ) then
            used = send_data (id_netlw_special(k,ipass), netlw_trop(:,:,k),  &
                              Time_diag, is, js )
          endif
        end do

!------- upward lw flux surface -------
        if ( id_lwup_sfc(ipass) > 0 ) then
          used = send_data (id_lwup_sfc(ipass), lwups,    &
                            Time_diag, is, js )
        endif

!------- downward lw flux surface -------
        if (id_lwdn_sfc(ipass) > 0 ) then
          used = send_data (id_lwdn_sfc(ipass), lwdns,    &
                            Time_diag, is, js )
        endif

!------- net lw flux surface -------
         if ( id_lwsfc(ipass) > 0 ) then
           used = send_data (id_lwsfc(ipass), lwups-lwdns,    &
                             Time_diag, is, js )
         endif
 
     if (do_lwaerosol_forcing) then

!------- outgoing lw flux toa (olr) with aerosols-------
        if (id_olr_ad(ipass) > 0 ) then
          used = send_data (id_olr_ad(ipass), olr_ad,    &
                             Time_diag, is, js )
        endif

!------- net lw flux surface -------
        if ( id_lwsfc_ad(ipass) > 0 ) then
           used = send_data (id_lwsfc_ad(ipass), lwups_ad-lwdns_ad,    &
                             Time_diag, is, js )
        endif
     endif

!----------------------------------------
!  longwave data for cmip names
!----------------------------------------

!------- lw tendency -----------
        if (query_cmip_diag_id(ID_tntrl)) then
          used = send_cmip_data_3d (ID_tntrl, tdtlw,    &
                            Time_diag, is, js, 1)
        endif

!------- downward lw flux surface -------
        if (id_rlds > 0 ) then
          used = send_data (id_rlds, lwdns,    &
                            Time_diag, is, js )
        endif

!------- upward lw flux surface -------
        if ( id_rlus > 0 ) then
          used = send_data (id_rlus, lwups, Time_diag, is, js )
        endif

!------- outgoing lw flux toa (olr) -------
        if (id_rlut > 0 ) then
          used = send_data (id_rlut, olr, Time_diag, is, js )
        endif

!------- net radiation (lw + sw) at top of atmos model-------
        if (id_rtmt > 0 ) then
          used = send_data (id_rtmt,   &
               swdn_trop(:,:,4) -swup_trop(:,:,4) -netlw_trop(:,:,4), &
                             Time_diag, is, js )
        endif

        if (do_lwaerosol_forcing) then

!------- downward lw flux surface (without aerosols) -------
          if ( id_rldsaf > 0 ) then
            used = send_data (id_rldsaf, lwdns_ad, Time_diag, is, js )
          endif

!------- outgoing lw flux toa (without aerosols) -------
          if (id_rlutaf > 0 ) then
            used = send_data (id_rlutaf, olr_ad, Time_diag, is, js )
          endif

        endif

!----------------------------------------------------------------------
!    now pass clear-sky diagnostics, if they have been calculated.
!----------------------------------------------------------------------
        if (Rad_control%do_totcld_forcing) then
          ipass = 1

!---- net lw flux ---------
        if (id_flxnet(ipass) > 0 ) then
          used = send_data (id_flxnet(ipass),    &
                            flxnetcf,   &
                            Time_diag, is, js, 1)
        endif

!------- lw tendency -----------
          if (id_tdt_lw(ipass) > 0 ) then
            used = send_data (id_tdt_lw(ipass), tdtlw_clr,    &
                              Time_diag, is, js, 1)
          endif

!------- outgoing lw flux toa (olr) -------
          if (id_olr(ipass) > 0 ) then
            used = send_data (id_olr(ipass), olr_clr,   &
                              Time_diag, is, js )
          endif

!------- net radiation (lw + sw) toa -------
          if (id_netrad_toa(ipass) > 0 ) then
            used = send_data (id_netrad_toa(ipass),   &
                              swin_clr - swout_clr - olr_clr,   &
                              Time_diag, is, js )
          endif

          do k = 1, MX_SPEC_LEVS
!------- net lw flux trop (netlw_trop) -------
            if (id_netlw_special(k,ipass) > 0 ) then
              used = send_data (id_netlw_special(k,ipass), netlw_trop_clr(:,:,k), &
                                Time_diag, is, js )
            endif
          end do

!------- upward lw flux surface -------
          if (id_lwup_sfc(ipass) > 0 ) then
            used = send_data (id_lwup_sfc(ipass), lwups_clr,   &
                              Time_diag, is, js )
          endif

!------- downward lw flux surface -------
          if (id_lwdn_sfc(ipass) > 0 ) then
            used = send_data (id_lwdn_sfc(ipass), lwdns_clr,   &
                              Time_diag, is, js )
          endif

!------- net lw flux surface -------
         if ( id_lwsfc(ipass) > 0 ) then
           used = send_data (id_lwsfc(ipass), lwups_clr-lwdns_clr,    &
                             Time_diag, is, js )
        endif   
   
     if (do_lwaerosol_forcing) then

!------- outgoing lw flux toa (olr) with aerosols-------
         if (id_olr_ad(ipass) > 0 ) then
          used = send_data (id_olr_ad(ipass), olr_ad_clr,    &
                            Time_diag, is, js )
         endif   
   
!------- net lw flux surface -------
         if ( id_lwsfc_ad(ipass) > 0 ) then
           used = send_data (id_lwsfc_ad(ipass), lwups_ad_clr-lwdns_ad_clr,    &
                           Time_diag, is, js )
        endif
      endif

!------------------------------------------
!  longwave clear-sky data for cmip names
!------------------------------------------

!------- lw tendency -----------
        if (query_cmip_diag_id(ID_tntrlcs)) then
          used = send_cmip_data_3d (ID_tntrlcs, tdtlw_clr,    &
                            Time_diag, is, js, 1)
        endif

!------- downward lw flux surface -------
          if (id_rldscs > 0 ) then
            used = send_data (id_rldscs, lwdns_clr, Time_diag, is, js )
          endif

!------- outgoing lw flux toa (olr) -------
          if (id_rlutcs > 0 ) then
            used = send_data (id_rlutcs, olr_clr, Time_diag, is, js )
          endif

          if (do_lwaerosol_forcing) then

!------- downward lw flux surface (clean-clear-sky) -------
            if (id_rldscsaf > 0 ) then
              used = send_data (id_rldscsaf, lwdns_ad_clr, Time_diag, is, js )
            endif

!------- outgoing lw flux toa (clean-clear-sky) -------
            if (id_rlutcsaf > 0 ) then
              used = send_data (id_rlutcsaf, olr_ad_clr, Time_diag, is, js )
            endif

          endif

        endif  ! (Rad_control%do_totcld_forcing)
        endif
      endif  ! (do_lw_rad .or. all_step_diagnostics)

!--------------------------------------------------------------------
!    now define various diagnostic integrals.
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!    accumulate global integral quantities 
!--------------------------------------------------------------------
      if (Rad_control%do_lw_rad .or. all_step_diagnostics) then
        olr_intgl(is:ie,js:je) = olr(:,:)

        if (id_rlut_g   > 0) call buffer_global_diag (id_rlut_g,   olr,     Time_diag, is, js)
        if (id_rlutcs_g > 0) call buffer_global_diag (id_rlutcs_g, olr_clr, Time_diag, is, js)
      endif

      if (Rad_control%renormalize_sw_fluxes .or. Rad_control%do_sw_rad .or.    &
          all_step_diagnostics) then
        swabs_intgl(is:ie,js:je) = swin(:,:) - swout(:,:)

        if (id_rsdt_g   > 0) call buffer_global_diag (id_rsdt_g,   swin,      Time_diag, is, js)
        if (id_rsut_g   > 0) call buffer_global_diag (id_rsut_g,   swout,     Time_diag, is, js)
        if (id_rsutcs_g > 0) call buffer_global_diag (id_rsutcs_g, swout_clr, Time_diag, is, js)
        if (id_rss_g    > 0) call buffer_global_diag (id_rss_g,  swdns-swups, Time_diag, is, js)
      endif

!--------------------------------------------------------------------

end subroutine produce_radiation_diagnostics

!###################################################################

subroutine radiation_driver_diag_end (Rad_control)

type(radiation_control_type),  intent(in) :: Rad_control

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_diag_mod',  &
               'module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    write restart file if desired; the file is not necessary if job 
!    ends on step prior to radiation ts, or if restart seamlessness 
!    is not required.
!---------------------------------------------------------------------
      if (Rad_control%using_restart_file) then
        if (Rad_control%renormalize_sw_fluxes) then
          call write_solar_interp_restart_nc
        endif
      endif

!---------------------------------------------------------------------
!    release space for renormalization arrays, if that option is active.
!---------------------------------------------------------------------
      if (Rad_control%renormalize_sw_fluxes .or. all_step_diagnostics)  then
         deallocate (solar_save)
         deallocate (Sw_flux_save%flux_sw_surf, &
                     Sw_flux_save%sw_heating, &
                     Sw_flux_save%flux_sw_surf_dir,   &
                     Sw_flux_save%flux_sw_surf_refl_dir,   &
                     Sw_flux_save%flux_sw_surf_dif,   &
                     Sw_flux_save%flux_sw_down_vis_dir,   &
                     Sw_flux_save%flux_sw_down_vis_dif,   &
                     Sw_flux_save%flux_sw_down_total_dir, &
                     Sw_flux_save%flux_sw_down_total_dif, &
                     Sw_flux_save%flux_sw_vis, &
                     Sw_flux_save%flux_sw_vis_dir, &
                     Sw_flux_save%flux_sw_refl_vis_dir, &
                     Sw_flux_save%flux_sw_vis_dif, &
                     Sw_flux_save%tot_heating, &
                     Sw_flux_save%dfsw, Sw_flux_save%ufsw,   &
                     Sw_flux_save%fsw, Sw_flux_save%hsw)
        if (do_swaerosol_forcing) then
          deallocate (Sw_flux_save%dfsw_ad, Sw_flux_save%ufsw_ad)
        endif
        if (Rad_control%do_totcld_forcing) then
          deallocate (Sw_flux_save%sw_heating_clr, &
                      Sw_flux_save%tot_heating_clr,  &
                      Sw_flux_save%dfswcf, Sw_flux_save%ufswcf, &
                      Sw_flux_save%fswcf,  Sw_flux_save%hswcf, &
                      Sw_flux_save%flux_sw_down_total_dir_clr, &
                      Sw_flux_save%flux_sw_down_total_dif_clr, &
                      Sw_flux_save%flux_sw_down_vis_clr )
          if (do_swaerosol_forcing) then
            deallocate (Sw_flux_save%dfswcf_ad, Sw_flux_save%ufswcf_ad)
          endif
        endif
      endif

!---------------------------------------------------------------------
!    release space needed when all_step_diagnostics is active.
!---------------------------------------------------------------------
      if (all_step_diagnostics)  then
        deallocate (Lw_flux_save%olr, &
                    Lw_flux_save%lwups, &
                    Lw_flux_save%lwdns, &
                    Lw_flux_save%flxnet, &
                    Lw_flux_save%tdtlw)
        if (do_lwaerosol_forcing) then
          deallocate (Lw_flux_save%olr_ad, &
                      Lw_flux_save%lwups_ad, &
                      Lw_flux_save%lwdns_ad)
        endif
        if (Rad_control%do_totcld_forcing) then
          deallocate (Lw_flux_save%olr_clr, &
                      Lw_flux_save%lwups_clr, &
                      Lw_flux_save%lwdns_clr, &
                      Lw_flux_save%flxnetcf, &
                      Lw_flux_save%tdtlw_clr)
          if (do_lwaerosol_forcing) then
            deallocate (Lw_flux_save%olr_ad_clr, &
                        Lw_flux_save%lwups_ad_clr,  &
                        Lw_flux_save%lwdns_ad_clr)
          endif
        endif
      endif

      module_is_initialized = .false.

end subroutine radiation_driver_diag_end

!###################################################################

subroutine solar_flux_save_init (is, ie, js, je, Sw_output, &
                                 Rad_output, do_totcld_forcing)

integer,                            intent(in) :: is, js, ie, je
type(sw_output_type), dimension(:), intent(in) :: Sw_output
type(rad_output_type),              intent(in) :: Rad_output
logical,                            intent(in) :: do_totcld_forcing

!----------------------------------------
!  module variables (initialized in init)
!     do_swaerosol_forcing
!     indx_swaf

    Sw_flux_save%flux_sw_surf         (is:ie,js:je) = Rad_output%flux_sw_surf         (is:ie,js:je)
    Sw_flux_save%flux_sw_surf_dir     (is:ie,js:je) = Rad_output%flux_sw_surf_dir     (is:ie,js:je)
    Sw_flux_save%flux_sw_surf_refl_dir(is:ie,js:je) = Rad_output%flux_sw_surf_refl_dir(is:ie,js:je)
    Sw_flux_save%flux_sw_surf_dif     (is:ie,js:je) = Rad_output%flux_sw_surf_dif     (is:ie,js:je)
    Sw_flux_save%flux_sw_down_vis_dir (is:ie,js:je) = Rad_output%flux_sw_down_vis_dir (is:ie,js:je)

    Sw_flux_save%flux_sw_down_vis_dif  (is:ie,js:je) = Rad_output%flux_sw_down_vis_dif  (is:ie,js:je)
    Sw_flux_save%flux_sw_down_total_dir(is:ie,js:je) = Rad_output%flux_sw_down_total_dir(is:ie,js:je)
    Sw_flux_save%flux_sw_down_total_dif(is:ie,js:je) = Rad_output%flux_sw_down_total_dif(is:ie,js:je)
    Sw_flux_save%flux_sw_vis           (is:ie,js:je) = Rad_output%flux_sw_vis           (is:ie,js:je)
    Sw_flux_save%flux_sw_vis_dir       (is:ie,js:je) = Rad_output%flux_sw_vis_dir       (is:ie,js:je)
    Sw_flux_save%flux_sw_refl_vis_dir  (is:ie,js:je) = Rad_output%flux_sw_refl_vis_dir  (is:ie,js:je)
    Sw_flux_save%flux_sw_vis_dif       (is:ie,js:je) = Rad_output%flux_sw_vis_dif       (is:ie,js:je)

    Sw_flux_save%sw_heating (is:ie,js:je,:) = Rad_output%tdtsw  (is:ie,js:je,:)
    Sw_flux_save%tot_heating(is:ie,js:je,:) = Rad_output%tdt_rad(is:ie,js:je,:)

    Sw_flux_save%dfsw       (is:ie,js:je,:) = Sw_output(1)%dfsw(:,:,:)
    Sw_flux_save%ufsw       (is:ie,js:je,:) = Sw_output(1)%ufsw(:,:,:)
    Sw_flux_save%fsw        (is:ie,js:je,:) = Sw_output(1)%fsw (:,:,:)
    Sw_flux_save%hsw        (is:ie,js:je,:) = Sw_output(1)%hsw (:,:,:)

    if (do_swaerosol_forcing) then
      Sw_flux_save%dfsw_ad  (is:ie,js:je,:) = Sw_output(indx_swaf)%dfsw(:,:,:)
      Sw_flux_save%ufsw_ad  (is:ie,js:je,:) = Sw_output(indx_swaf)%ufsw(:,:,:)
    endif

    if (do_totcld_forcing) then 
      Sw_flux_save%sw_heating_clr (is:ie,js:je,:) = Rad_output%tdtsw_clr(is:ie,js:je,:)
      Sw_flux_save%tot_heating_clr(is:ie,js:je,:) = Rad_output%tdt_rad_clr(is:ie,js:je,:)
      Sw_flux_save%dfswcf         (is:ie,js:je,:) = Sw_output(1)%dfswcf(:,:,:)
      Sw_flux_save%ufswcf         (is:ie,js:je,:) = Sw_output(1)%ufswcf(:,:,:)
      Sw_flux_save% fswcf         (is:ie,js:je,:) = Sw_output(1)%fswcf (:,:,:)
      Sw_flux_save% hswcf         (is:ie,js:je,:) = Sw_output(1)%hswcf (:,:,:)
      Sw_flux_save%flux_sw_down_total_dir_clr(is:ie,js:je) = Rad_output%flux_sw_down_total_dir_clr(is:ie,js:je)
      Sw_flux_save%flux_sw_down_total_dif_clr(is:ie,js:je) = Rad_output%flux_sw_down_total_dif_clr(is:ie,js:je)
      Sw_flux_save%flux_sw_down_vis_clr      (is:ie,js:je) = Rad_output%flux_sw_down_vis_clr      (is:ie,js:je)
      if (do_swaerosol_forcing) then
        Sw_flux_save%dfswcf_ad  (is:ie,js:je,:) = Sw_output(indx_swaf)%dfswcf(:,:,:)
        Sw_flux_save%ufswcf_ad  (is:ie,js:je,:) = Sw_output(indx_swaf)%ufswcf(:,:,:)
      endif
    endif

end subroutine solar_flux_save_init

!###################################################################

subroutine solar_interp_register_restart(fname, do_totcld_forcing)
  character(len=*), intent(in) :: fname
  logical,          intent(in) :: do_totcld_forcing

  character(len=64)            :: fname2
  integer                      :: id_restart

   ! all data is distributed on tile files
  !call get_mosaic_tile_file(fname, fname2, .false. ) 
   allocate(Tile_restart)
   doing_netcdf_restart = .true.

!  NOTE: there could a problem when restarting a model with renormalize_sw_fluxes = true
!  if the restart file has int_renormalize_sw_fluxes = 0 (i.e., no sw fluxes saved)

   id_restart = register_restart_field(Tile_restart, fname, 'solar_save', solar_save)
   id_restart = register_restart_field(Tile_restart, fname, 'flux_sw_surf_save', Sw_flux_save%flux_sw_surf)
   id_restart = register_restart_field(Tile_restart, fname, 'flux_sw_surf_dir_save', Sw_flux_save%flux_sw_surf_dir)
   id_restart = register_restart_field(Tile_restart, fname, 'flux_sw_surf_refl_dir_save', Sw_flux_save%flux_sw_surf_refl_dir)
   id_restart = register_restart_field(Tile_restart, fname, 'flux_sw_surf_dif_save', Sw_flux_save%flux_sw_surf_dif)
   id_restart = register_restart_field(Tile_restart, fname, 'flux_sw_down_vis_dir_save', Sw_flux_save%flux_sw_down_vis_dir)
   id_restart = register_restart_field(Tile_restart, fname, 'flux_sw_down_vis_dif_save', Sw_flux_save%flux_sw_down_vis_dif)
   id_restart = register_restart_field(Tile_restart, fname, 'flux_sw_down_total_dir_save', Sw_flux_save%flux_sw_down_total_dir)
   id_restart = register_restart_field(Tile_restart, fname, 'flux_sw_down_total_dif_save', Sw_flux_save%flux_sw_down_total_dif)
   id_restart = register_restart_field(Tile_restart, fname, 'flux_sw_vis_save', Sw_flux_save%flux_sw_vis)
   id_restart = register_restart_field(Tile_restart, fname, 'flux_sw_vis_dir_save', Sw_flux_save%flux_sw_vis_dir)
   id_restart = register_restart_field(Tile_restart, fname, 'flux_sw_refl_vis_dir_save', Sw_flux_save%flux_sw_refl_vis_dir)
   id_restart = register_restart_field(Tile_restart, fname, 'flux_sw_vis_dif_save', Sw_flux_save%flux_sw_vis_dif)
   id_restart = register_restart_field(Tile_restart, fname, 'sw_heating_save', Sw_flux_save%sw_heating)
   id_restart = register_restart_field(Tile_restart, fname, 'tot_heating_save', Sw_flux_save%tot_heating)
   id_restart = register_restart_field(Tile_restart, fname, 'dfsw_save', Sw_flux_save%dfsw)
   id_restart = register_restart_field(Tile_restart, fname, 'ufsw_save', Sw_flux_save%ufsw)
   id_restart = register_restart_field(Tile_restart, fname, 'fsw_save', Sw_flux_save%fsw)
   id_restart = register_restart_field(Tile_restart, fname, 'hsw_save', Sw_flux_save%hsw)
   if (do_totcld_forcing) then
      id_restart = register_restart_field(Tile_restart, fname, 'sw_heating_clr_save', Sw_flux_save%sw_heating_clr)
      id_restart = register_restart_field(Tile_restart, fname, 'tot_heating_clr_save', Sw_flux_save%tot_heating_clr)
      id_restart = register_restart_field(Tile_restart, fname, 'dfswcf_save', Sw_flux_save%dfswcf)
      id_restart = register_restart_field(Tile_restart, fname, 'ufswcf_save', Sw_flux_save%ufswcf)
      id_restart = register_restart_field(Tile_restart, fname, 'fswcf_save', Sw_flux_save%fswcf)
      id_restart = register_restart_field(Tile_restart, fname, 'hswcf_save', Sw_flux_save%hswcf)
      id_restart = register_restart_field(Tile_restart, fname, 'flux_sw_down_total_dir_clr_save', Sw_flux_save%flux_sw_down_total_dir_clr)
      id_restart = register_restart_field(Tile_restart, fname, 'flux_sw_down_total_dif_clr_save', Sw_flux_save%flux_sw_down_total_dif_clr)
      id_restart = register_restart_field(Tile_restart, fname, 'flux_sw_down_vis_clr_save', Sw_flux_save%flux_sw_down_vis_clr)
   endif

end subroutine solar_interp_register_restart

!###################################################################

subroutine write_solar_interp_restart_nc (timestamp)
character(len=*), intent(in), optional :: timestamp

      if (.not.doing_netcdf_restart) return
!---------------------------------------------------------------------
!    only the root pe will write control information -- the last value 
!    in the list of restart versions and the alarm information.
!---------------------------------------------------------------------
      if (mpp_pe() == mpp_root_pe() ) then
        call error_mesg('radiation_driver_diag_mod', 'Writing netCDF formatted restart file:'//&
                        '  RESTART/radiation_driver.res.nc (solar interp data)', NOTE)
      endif

!---------------------------------------------------------------------
!    write out the optional time average restart data. note that 
!    do_average and renormalize_sw_fluxes may not both be true.
!---------------------------------------------------------------------

! Make sure that the restart_versions variable is up to date.
      call save_restart(Tile_restart, timestamp)

end subroutine write_solar_interp_restart_nc

!###################################################################

subroutine radiation_driver_diag_endts (Rad_control)

type(radiation_control_type),  intent(in) :: Rad_control

logical :: used
!---------------------------------------------------------------------
!    compute and write out global integrals
!---------------------------------------------------------------------

      if (Rad_control%do_lw_rad .or. all_step_diagnostics) then
        call sum_diag_integral_field ('olr',    olr_intgl)

        if (id_rlut_g   > 0) used = send_global_diag (id_rlut_g)
        if (id_rlutcs_g > 0) used = send_global_diag (id_rlutcs_g)
      endif

      if (Rad_control%renormalize_sw_fluxes .or. Rad_control%do_sw_rad .or.    &
          all_step_diagnostics) then
        call sum_diag_integral_field ('abs_sw', swabs_intgl )

        if (id_rsdt_g   > 0) used = send_global_diag (id_rsdt_g)
        if (id_rsut_g   > 0) used = send_global_diag (id_rsut_g)
        if (id_rsutcs_g > 0) used = send_global_diag (id_rsutcs_g)
        if (id_rss_g    > 0) used = send_global_diag (id_rss_g)
      endif

!---------------------------------------------------------------------

end subroutine radiation_driver_diag_endts

!###################################################################
! <SUBROUTINE NAME="flux_trop_calc">
!  <OVERVIEW>
!    flux_trop_calc defines the shortwave and longwave fluxes at the
!    tropopause immediately after the computation of fluxes at model
!    levels by the radiation algorithms (invoked by radiation_calc).
!  </OVERVIEW>
!  <DESCRIPTION>
!    flux_trop_calc defines the shortwave and longwave fluxes at the
!    tropopause immediately after the computation of fluxes at model
!    levels by the radiation algorithms (invoked by radiation_calc).
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_trop_calc          (lat, pflux, &
!                                 Lw_output, Sw_output, &
!                                 swdn_trop, swup_trop, &
!                                 swdn_trop_clr, swup_trop_clr, &
!                                 netlw_trop, netlw_trop_clr )
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!    mean latitude (in radians) of all grid boxes processed by this
!    call to flux_trop_calc   [real, dimension(:,:)]
!  </IN>
!  <IN NAME="pflux" TYPE="real">
!   pressure (in pascals) at layer boundaries [real, dimension(:,:,:)]
!  </IN>
!  <INOUT NAME="Lw_output" TYPE="lw_output_type">
!      longwave radiation output data from the 
!                        sea_esf_rad radiation package, when that 
!                        package is active
!  </INOUT>
!  <INOUT NAME="Sw_output" TYPE="sw_output_type">
!   shortwave radiation output data from the 
!                        sea_esf_rad radiation package  when that 
!                        package is active
!  </INOUT>
! </SUBROUTINE>


subroutine flux_trop_calc (lat, pflux, phalf, &
                           Rad_control, Lw_output, Sw_output, &
                           swdn_trop, swup_trop, &
                           swdn_trop_clr, swup_trop_clr, &
                           netlw_trop, netlw_trop_clr)

real,dimension(:,:),          intent(in)  :: lat
real,dimension(:,:,:),        intent(in)  :: pflux, phalf
type(radiation_control_type), intent(in)  :: Rad_control
type(lw_output_type),         intent(in)  :: Lw_output
type(sw_output_type),         intent(in)  :: Sw_output
real, dimension(:,:,:),       intent(out) :: swdn_trop, swup_trop, &
                                             swdn_trop_clr, swup_trop_clr
real, dimension(:,:,:),       intent(out) :: netlw_trop, netlw_trop_clr

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      lat          latitude of model points  [ radians ]
!      pflux        pressure at layer boundaries [ Pa ]
!      phalf        pressure at dynamical layer boundaries [ Pa ]
!      Lw_output    lw_output_type variable containing output from 
!                   the longwave radiation code of the
!                   sea_esf_rad package, on the model grid
!      Sw_output    sw_output_type variable containing output from 
!                   the shortwave radiation code of the
!                   sea_esf_rad package, on the model grid

      real, dimension (size(pflux,1),size(pflux,2)) :: &
                                                lat_deg, tropo_ht


      integer           :: j, k
      integer           :: ki, i
      integer           :: kmax
      real              :: wtlo, wthi

      kmax = size(pflux,3) - 1

!---------------------------------------------------------------------
!    compute net downward flux at 1 Pa (top of dynmamical model)
!    here dynamical pressure top is hard-wired to 1 Pa.
!---------------------------------------------------------------------
 
      do j = 1, size(pflux,2)
        do i = 1, size(pflux,1)
         !wtlo = (1.0          - pflux(i,j,1))/ &
          wtlo = (phalf(i,j,1) - pflux(i,j,1))/ &
                 (pflux(i,j,2) - pflux(i,j,1))
          wthi = 1.0 - wtlo
          if (Rad_control%do_lw_rad) then
            netlw_trop(i,j,4) = wthi*Lw_output%flxnet(i,j,1) + &
                                wtlo*Lw_output%flxnet(i,j,2)
            if (Rad_control%do_totcld_forcing) then
              netlw_trop_clr(i,j,4) = wthi*Lw_output%flxnetcf(i,j,1) + &
                                      wtlo*Lw_output%flxnetcf(i,j,2)
            endif
          endif
          if (Rad_control%do_sw_rad) then
            swdn_trop(i,j,4) = wthi*Sw_output%dfsw(i,j,1) + &
                               wtlo*Sw_output%dfsw(i,j,2)
            swup_trop(i,j,4) = wthi*Sw_output%ufsw(i,j,1) + &
                               wtlo*Sw_output%ufsw(i,j,2)
            if (Rad_control%do_totcld_forcing) then
              swdn_trop_clr(i,j,4) = wthi*Sw_output%dfswcf(i,j,1) +&
                                     wtlo*Sw_output%dfswcf(i,j,2)
              swup_trop_clr(i,j,4) = wthi*Sw_output%ufswcf(i,j,1) +&
                                     wtlo*Sw_output%ufswcf(i,j,2)
            endif
          endif
        enddo
      enddo


      if (constant_tropo) then
        tropo_ht(:,:) = trop_ht_constant
! interpolate the fluxes between the appropriate pressures bracketing
! (trop) 

        do j = 1, size(pflux,2)
          do i = 1, size(pflux,1)
            do k = kmax+1,2,-1
              if (pflux(i,j,k) >= tropo_ht(i,j) .and.     &
                  pflux(i,j,k-1) < tropo_ht(i,j))      then
                ki = k
!   the indices for high,low pressure bracketing "tropo_ht" are ki, ki-1
                wtlo = (tropo_ht(i,j) - pflux(i,j,ki-1))/ &
                       (pflux(i,j,ki) - pflux(i,j,ki-1))
                wthi = 1.0 - wtlo
                if (Rad_control%do_lw_rad) then
                  netlw_trop(i,j,1) = wtlo*Lw_output%flxnet(i,j,ki) + &
                                      wthi*Lw_output%flxnet(i,j,ki-1)
                  if (Rad_control%do_totcld_forcing) then
                    netlw_trop_clr(i,j,1) = wtlo*Lw_output%flxnetcf(i,j,ki) + &
                                            wthi*Lw_output%flxnetcf(i,j,ki-1)
                  endif
                endif
                if (Rad_control%do_sw_rad) then
                  swdn_trop(i,j,1) = wtlo*Sw_output%dfsw(i,j,ki) + &
                                     wthi*Sw_output%dfsw(i,j,ki-1)
                  swup_trop(i,j,1) = wtlo*Sw_output%ufsw(i,j,ki) + &
                                     wthi*Sw_output%ufsw(i,j,ki-1)
                  if (Rad_control%do_totcld_forcing) then
                    swdn_trop_clr(i,j,1) = wtlo*Sw_output%dfswcf(i,j,ki) +&
                                           wthi*Sw_output%dfswcf(i,j,ki-1)
                    swup_trop_clr(i,j,1) = wtlo*Sw_output%ufswcf(i,j,ki) +&
                                           wthi*Sw_output%ufswcf(i,j,ki-1)
                  endif
                endif
                exit
              endif
            enddo !k
          enddo   !i
        enddo     !j
      endif


      if (linear_tropo) then
        lat_deg(:,:) = lat(:,:)*RADIAN
        tropo_ht(:,:) = trop_ht_at_eq + ABS(lat_deg(:,:))*  &
                        (trop_ht_at_poles - trop_ht_at_eq)/90.
! interpolate the fluxes between the appropriate pressures bracketing

        do i = 1,size(pflux,1)
          do j = 1,size(pflux,2)
            do k = kmax+1,2,-1
              if (pflux(i,j,k) >= tropo_ht(i,j) .and.     &
                  pflux(i,j,k-1) < tropo_ht(i,j))      then
                ki = k
!   the indices for high,low pressure bracketing "tropo_ht" are ki, ki-1
                wtlo = (tropo_ht(i,j) - pflux(i,j,ki-1))/ &
                       (pflux(i,j,ki) - pflux(i,j,ki-1))
                wthi = 1.0 - wtlo
                if (Rad_control%do_lw_rad) then
                  netlw_trop(i,j,2) = wtlo*Lw_output%flxnet(i,j,ki) + &
                                      wthi*Lw_output%flxnet(i,j,ki-1)
                  if (Rad_control%do_totcld_forcing) then
                    netlw_trop_clr(i,j,2) = wtlo*Lw_output%flxnetcf(i,j,ki) + &
                                            wthi*Lw_output%flxnetcf(i,j,ki-1)
                  endif
                endif
                if (Rad_control%do_sw_rad) then
                  swdn_trop(i,j,2) = wtlo*Sw_output%dfsw(i,j,ki) + &
                                     wthi*Sw_output%dfsw(i,j,ki-1)
                  swup_trop(i,j,2) = wtlo*Sw_output%ufsw(i,j,ki) + &
                                     wthi*Sw_output%ufsw(i,j,ki-1)
                  if (Rad_control%do_totcld_forcing) then
                    swdn_trop_clr(i,j,2) = wtlo*Sw_output%dfswcf(i,j,ki) + &
                                           wthi*Sw_output%dfswcf(i,j,ki-1)
                    swup_trop_clr(i,j,2) = wtlo*Sw_output%ufswcf(i,j,ki) + &
                                           wthi*Sw_output%ufswcf(i,j,ki-1)
                  endif
                endif
                exit
              endif
            enddo !k
          enddo !j
        enddo !i
      endif


      if (thermo_tropo) then
        call error_mesg ( 'radiation_driver_mod', &
              'thermo_tropo option not yet available', FATAL)
! interpolate the fluxes between the appropriate pressures bracketing

        do i = 1, size(pflux,1)
          do j = 1, size(pflux,2)
            do k = kmax+1,2,-1
              if (pflux(i,j,k) >= tropo_ht(i,j) .and.     &
                  pflux(i,j,k-1) < tropo_ht(i,j))      then
                ki = k
!   the indices for high,low pressure bracketing "tropo_ht" are ki, ki-1
                wtlo = (tropo_ht(i,j) - pflux(i,j,ki-1))/ &
                       (pflux(i,j,ki) - pflux(i,j,ki-1))
                wthi = 1.0 - wtlo
                if (Rad_control%do_lw_rad) then
                  netlw_trop(i,j,3) = wtlo*Lw_output%flxnet(i,j,ki) + &
                                      wthi*Lw_output%flxnet(i,j,ki-1)
                  if (Rad_control%do_totcld_forcing) then
                    netlw_trop_clr(i,j,3) = wtlo*Lw_output%flxnetcf(i,j,ki) + &
                                          wthi*Lw_output%flxnetcf(i,j,ki-1)
                  endif
                endif
                if (Rad_control%do_sw_rad) then
                  swdn_trop(i,j,3) = wtlo*Sw_output%dfsw(i,j,ki) + &
                                     wthi*Sw_output%dfsw(i,j,ki-1)
                  swup_trop(i,j,3) = wtlo*Sw_output%ufsw(i,j,ki) + &
                                     wthi*Sw_output%ufsw(i,j,ki-1)
                  if (Rad_control%do_totcld_forcing) then
                    swdn_trop_clr(i,j,3) = wtlo*Sw_output%dfswcf(i,j,ki) + &
                                           wthi*Sw_output%dfswcf(i,j,ki-1)
                    swup_trop_clr(i,j,3) = wtlo*Sw_output%ufswcf(i,j,ki) + &
                                           wthi*Sw_output%ufswcf(i,j,ki-1)
                  endif
                endif
                exit
              endif
            enddo !k
          enddo !j
        enddo !i
      endif

!----------------------------------------------------------------------

end subroutine flux_trop_calc
                                       
!####################################################################

                end module radiation_driver_diag_mod

