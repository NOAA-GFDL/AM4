                module radiation_driver_mod
!> <CONTACT EMAIL="Fei.Liu@noaa.gov">
!!  fil
!! </CONTACT>
!! <REVIEWER EMAIL="">
!!  
!! </REVIEWER>
!! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
!! <OVERVIEW>
!!    radiation_driver_mod is the interface between physics_driver_mod
!!    and the radiation parameterizations, currently only the
!!    sea_esf_rad radiation package. it provides 
!!    radiative heating rates, boundary radiative fluxes, and any other 
!!    radiation package output fields to other component models of the
!!    modeling system.
!! </OVERVIEW>
!! <DESCRIPTION>
!! The following modules are called from this driver module:
!!
!!   1) astronomy
!!
!!   2) cloud properties
!!
!!   3) prescribed zonal ozone
!!
!!   4) longwave and shortwave radiation driver
!! </DESCRIPTION>
!! <INFO>
!
!!   <REFERENCE>  For a specific list of radiation references see the
!!     longwave and shortwave documentation.          </REFERENCE>
!!   <BUG>
!!For some of the diagnostics fields that represent fractional amounts,
!!    such as reflectivity and absorptivity, the units are incorrectly
!!    given as percent.
!!</BUG>
!! </INFO>
!    shared modules:

use block_control_mod,     only: block_control_type
use mpp_mod,               only: input_nml_file
use fms_mod,               only: fms_init, mpp_clock_id, &
                                 mpp_clock_begin, mpp_clock_end, &
                                 CLOCK_MODULE_DRIVER, CLOCK_MODULE, &
                                 CLOCK_ROUTINE, &
                                 field_exist, field_size, &
                                 mpp_pe, mpp_root_pe, &
                                 open_namelist_file, stdlog, stdout, &
                                 file_exist, FATAL, WARNING, NOTE, &
                                 close_file, read_data, write_data, &
                                 write_version_number, check_nml_error,&
                                 error_mesg, mpp_chksum, &
                                 read_data, mpp_error
use fms_io_mod,            only: restore_state, &
                                 register_restart_field, restart_file_type, &
                                 save_restart, get_mosaic_tile_file
use time_manager_mod,      only: time_type, set_date, set_time,  &
                                 get_time,    operator(+),       &
                                 print_date, time_manager_init, &
                                 assignment(=), &
                                 operator(-), operator(/=), get_date,&
                                 operator(<), operator(>=), operator(>)
use sat_vapor_pres_mod,    only: sat_vapor_pres_init, compute_qs
use constants_mod,         only: constants_init, RDGAS, RVGAS,   &
                                 STEFAN, GRAV, SECONDS_PER_DAY,  &
                                 RADIAN, diffac
use data_override_mod,     only: data_override

use field_manager_mod,     only: MODEL_ATMOS
use tracer_manager_mod,    only: tracer_manager_init, &
                                 get_number_tracers, &
                                 get_tracer_index, NO_TRACER

! shared radiation package modules:

use radiation_types_mod,   only: radiation_type, &
                                 radiation_input_control_type, &
                                 radiation_input_block_type, &
                                 radiation_input_glblqty_type
use physics_radiation_exch_mod,only: exchange_control_type, &
                                     clouds_from_moist_block_type, &
                                     radiation_flux_control_type, &
                                     radiation_flux_block_type, &
                                     radiation_flux_type, &
                                     cosp_from_rad_type, &
                                     cosp_from_rad_control_type, &
                                     cosp_from_rad_block_type, &
                                     alloc_radiation_flux_type, &
                                     dealloc_radiation_flux_type

use radiation_driver_types_mod, only: radiation_control_type, &
                                      astronomy_type, surface_type, &
                                      atmos_input_type, &
                                      rad_output_type, astronomy_inp_type, &
                                      rad_output_init

use     aerosol_types_mod, only: aerosol_type, &
                                 aerosol_time_vary_type

use  aerosolrad_types_mod, only: aerosolrad_control_type, &
                                 aerosolrad_diag_type

use    cloudrad_types_mod, only: cld_specification_type, &
                                 microphysics_type, &
                                 cloudrad_control_type

!  physics support modules:

use astronomy_mod,         only: astronomy_init, annual_mean_solar, &
                                 daily_mean_solar, diurnal_solar, &
                                 astronomy_end

!  radiation component modules:

use longwave_driver_mod,  only: longwave_driver_init,   &   
                                longwave_driver_time_vary, &
                                longwave_driver, &
                                longwave_driver_endts, &
                                longwave_driver_end, &
                                longwave_number_of_bands, &
                                lw_table_type, &
                                longwave_get_tables

use longwave_types_mod,    only: lw_output_type, lw_diagnostics_type, assignment(=)

use shortwave_driver_mod, only: shortwave_driver_init,  &
                                shortwave_driver,  &
                                shortwave_driver_end, &
                                shortwave_driver_time_vary, &
                                shortwave_number_of_bands, &
                                get_solar_constant

use shortwave_types_mod,   only: sw_output_type, assignment(=)

use rad_output_file_mod,   only: rad_output_file_init, &
                                 write_rad_output_file,    &
                                 rad_output_file_end

use aerosolrad_driver_mod, only: aerosolrad_driver_init, &
                                 aerosolrad_driver_time_vary, &
                                 aerosolrad_driver, &
                                 aerosolrad_driver_endts, &
                                 aerosolrad_driver_end, &
                                 aerosolrad_driver_dealloc

use cloudrad_driver_mod,   only: cloudrad_driver_init, &
                                 cloudrad_driver, &
                                 cloudrad_driver_end

use cloudrad_diagnostics_mod,      &
                           only: obtain_cloud_tau_and_em, &
                                 modis_yim, modis_cmip

use radiative_gases_mod,   only: radiative_gases_init,   &
                                 radiative_gases_time_vary, &
                                 radiative_gases_endts, &
                                 define_radiative_gases, &
                                 radiative_gases_end,     &    
                                 radiative_gases_restart, &
                                 get_longwave_gas_flag

use radiative_gases_types_mod, only: radiative_gases_type, &
                                     assignment(=)

use radiation_driver_diag_mod, only: radiation_driver_diag_init, &
                                     radiation_driver_diag_end, &
                                     radiation_driver_diag_endts, &
                                     update_rad_fields, &
                                     produce_radiation_diagnostics, &
                                     write_solar_interp_restart_nc

!--------------------------------------------------------------------

implicit none 
private 

!----------------------------------------------------------------------
!    radiation_driver_mod is the interface between physics_driver_mod
!    and the radiation parameterizations, currently only the
!    sea_esf_rad radiation package. it provides 
!    radiative heating rates, boundary radiative fluxes, and any other 
!    radiation package output fields to other component models of the
!    modeling system.
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!------------ version number for this module --------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'


!---------------------------------------------------------------------
!------ interfaces -----
! <PUBLIC>
!use radiation_driver_mod [,only: radiation_driver_init,
!                                 radiation_driver,
!                                 radiation_driver_end]
!   radiation_driver_init
!       Must be called once before subroutine radiation_driver to
!       initialize the module (read namelist input and restart file).
!       Also calls the initialization routines for other modules used.
!   radiation_driver
!       Called every time step (not on the radiation time step)
!       to compute the longwave and shortwave radiative tendencies.
!   radiation_driver_end
!       Called once at the end of a model run to terminate the module (write
!       a restart file). Also calls the termination routines for other
!       modules used.
!Notes:
! 1) A namelist interface controls runtime options.
! 3) A restart file radiation_driver.res is generated by this module.
!</PUBLIC>

public    radiation_driver_init, radiation_driver, &
          radiation_driver_time_vary, radiation_driver_endts, &
          radiation_driver_end, radiation_driver_restart, &
          radiation_diag_type_dealloc

! inherited from other modules
public    lw_output_type, sw_output_type


private  & 


! called from radiation_driver_end:
          write_restart_nc, &

! called from radiation_driver:
          obtain_astronomy_variables, radiation_calc,    &
          deallocate_arrays, &
          return_cosp_inputs, &
          define_atmos_input_fields, atmos_input_dealloc, &
          define_surface, surface_dealloc, &

! called from radiation_driver_time_vary:
          define_rad_times, &

! called from define_atmos_input_fields:
          calculate_auxiliary_variables


!-----------------------------------------------------------------------
!------- namelist ---------
logical :: using_restart_file = .true. ! if set to .false, restart file
                                       ! will NOT be written by this 
                                       ! module; this will not affect
                                       ! answers as long as job is 
                                       ! restarted on a radiation
                                       ! timestep
logical ::  do_clear_sky_pass= .false.!  are the clear-sky radiation
                                      !  diagnostics to be calculated ?
character(len=24) ::    &
            zenith_spec = '      '    !  string defining how zenith 
                                      !  angle is computed. acceptable
                                      !  values: 'daily_mean', 'annual_
                                      !  mean', 'diurnally_varying'
character(len=16) ::   &
                rad_package='sea_esf' !  string defining the radiation
                                      !  package being used. acceptable
                                      !  values : 'sea_esf'
logical ::     &
         renormalize_sw_fluxes=.false.!  should sw fluxes and the zenith
                                      !  angle be renormalized on each 
                                      !  timestep because of the 
                                      !  movement of earth wrt the sun ?
integer, dimension(6) ::    &
    rad_date = (/ 0, 0, 0, 0, 0, 0 /) !  fixed date for which radiation
                                      !  is to be valid (applies to
                                      !  solar info, ozone, clouds)
                                      !  [yr, mo, day, hr, min, sec]
logical :: rsd=.false.                !  (repeat same day) - call 
                                      !  radiation for the specified 
                                      !  rad_date (yr,mo,day), but run 
                                      !  through the diurnal cycle (hr,
                                      !  min,sec)

logical :: use_mixing_ratio = .false. !  assumes q is mixing ratio
                                      !  rather than specific humidity
logical :: doing_data_override = .false.  
                                      !  input fields to the radiation
                                      !  package are being overriden
                                      !  using data_override_mod ?
logical :: overriding_temps = .false. !  temperature and ts fields are
                                      !  overriden ?
logical :: overriding_sphum = .false. !  specific humidity field is
                                      !  overriden ?
logical :: overriding_clouds = .false.!  cloud specification fields are
                                      !  overriden ?
logical :: overriding_albedo = .false.!  surface albedo field is
                                      !  overriden ?
logical :: overriding_aerosol = .false.
                                      !  aerosol fields are overriden ?
logical :: use_co2_tracer_field = .false.
                                      !  obtain co2 field for use by 
                                      !  radiation package from co2
                                      !  tracer field ?
logical :: use_uniform_solar_input = .false.
                                      !  the (lat,lon) values used to
                                      !  calculate zenith angle are
                                      !  uniform across the grid ?
real    :: lat_for_solar_input = 100. !  latitude to be used when uni-
                                      !  form solar input is activated
                                      !  [ degrees ]
real    :: lon_for_solar_input = 500. !  longitude to be used when uni-
                                      !  form solar input is activated
                                      !  [ degrees ]

logical :: always_calculate = .false. !  radiation calculation is done
                                      !  on every call to 
                                      !  radiation_driver ?

logical :: treat_sfc_refl_dir_as_dif = .true.
                                      ! when true, solar direct  beam
                                      ! radiation reflected from the
                                      ! surface is seen as diffuse by
                                      ! the exchange grid. when false, it 
                                      ! is seen as direct, changing solar
                                      ! input to the sfc, and eliminating
                                      ! negative diffuse sw fluxes at the
                                      ! sfc, which cause problems in ESM.

integer ::  rad_time_step = 0         !  radiative time step in seconds


integer ::  sw_rad_time_step = 0      !  radiative time step in seconds
logical :: use_single_lw_sw_ts = .true. ! lw and sw are integrated
                                        ! using rad_time_step ? if 
                                       ! false, then lw uses 
                                       ! rad_time_step, sw uses 
                                       ! sw_rad_time_step

logical :: apply_temp_limits = .true.  ! upper and lower limits are
                                       ! are applied to temperature
logical :: apply_vapor_limits = .true. ! lower limit is applied to
                                       ! water vapor

logical :: nonzero_rad_flux_init = .false.

logical :: do_radiation = .true.

logical :: do_conserve_energy = .false.
                                      ! when true, the actually model layer
                                      ! interface pressures (p_half) are used
                                      ! for the radiation calculation thus
                                      ! conserving energy

! <NAMELIST NAME="radiation_driver_nml">
!  <DATA NAME="rad_time_step" UNITS="" TYPE="integer" DIM="" DEFAULT="14400">
!The radiative time step in seconds.
!  </DATA>
!  <DATA NAME="do_clear_sky_pass" UNITS="" TYPE="logical" DIM="" DEFAULT="">
! are the clear-sky radiation
!  diagnostics to be calculated ?
!  </DATA>
!  <DATA NAME="zenith_spec" UNITS="" TYPE="character" DIM="" DEFAULT="">
!string defining how zenith 
!  angle is computed. acceptable
!  values: 'daily_mean', 'annual_
!  mean', 'diurnally_varying'
!  </DATA>
!  <DATA NAME="rad_package" UNITS="" TYPE="character" DIM="" DEFAULT="">
!string defining the radiation
!  package being used. acceptable
!  values : 'sea_esf'
!  </DATA>
!  </DATA>
!  <DATA NAME="renormalize_sw_fluxes" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!should sw fluxes and the zenith
!  angle be renormalized on each 
!  timestep because of the 
!  movement of earth wrt the sun ?
!  </DATA>
!  <DATA NAME="rad_date" UNITS="" TYPE="integer" DIM="" DEFAULT="">
!fixed date for which radiation
!  is to be valid (applies to
!  solar info, ozone, clouds)
!  [yr, mo, day, hr, min, sec]
!  </DATA>
!  <DATA NAME="rsd" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!(repeat same day) - call 
!  radiation for the specified 
!  rad_date (yr,mo,day), but run 
!  through the diurnal cycle (hr,
!  min,sec)
!  </DATA>
!  <DATA NAME="use_mixing_ratio" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!assumes q is mixing ratio
!  rather than specific humidity
!  </DATA>
!  <DATA NAME="doing_data_override" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!input fields to the radiation
!  package are being overriden
!  using data_override_mod ?
!  </DATA>
!  <DATA NAME="overriding_temps" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!temperature and ts fields are
!  overriden ?
!  </DATA>
!  <DATA NAME="overriding_sphum" UNITS="" TYPE="logical" DIM="" DEFAULT="">
! specific humidity field is
!  overriden ?
!  </DATA>
!  <DATA NAME="overriding_clouds" UNITS="" TYPE="logical" DIM="" DEFAULT="">
! cloud specification fields are
!  overriden ?
!  </DATA>
!  <DATA NAME="overriding_albedo" UNITS="" TYPE="logical" DIM="" DEFAULT="">
! surface albedo field is
!  overriden ?
!  </DATA>
!  <DATA NAME="overriding_aerosol" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!aerosol fields are overriden ?
!  </DATA>
!  <DATA NAME="use_co2_tracer_field" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!use co2 value from co2 tracer field?
!  </DATA>
!  <DATA NAME="always_calculate" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!  calculate radiative fluxes and heating rates on every call to 
!  radiation_driver ?
!  </DATA>
!  <DATA NAME="always_calculate" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!fluxes and heating rates should
! be calculatd on each call to
! radiation_driver ? (true for
! standalone applications)
!  </DATA>
!  <DATA NAME="use_uniform_solar_input" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!  the (lat,lon) values used to
!  calculate zenith angle are
!  uniform across the grid ?
!  </DATA>
!  <DATA NAME="lat_for_solar_input" UNITS="" TYPE="real" DIM="" DEFAULT="">
!  latitude to be used when uni-
!  form solar input is activated
!  [ degrees ]
!  </DATA>
!  <DATA NAME="lon_for_solar_input" UNITS="" TYPE="real" DIM="" DEFAULT="">
!  longitude to be used when uni-
!  form solar input is activated
!  [ degrees ]
!  </DATA>
! </NAMELIST>
!
namelist /radiation_driver_nml/ do_radiation, &
                                do_clear_sky_pass, &
                                using_restart_file, &
                                zenith_spec, rad_package,    &
                                renormalize_sw_fluxes, &
                                rad_date, rsd,    &
                                use_mixing_ratio, &
                                doing_data_override, &
                                overriding_temps, overriding_sphum, &
                                overriding_clouds, overriding_albedo, &
                                overriding_aerosol, &
                                use_co2_tracer_field, &
                                treat_sfc_refl_dir_as_dif, &
                                always_calculate,  &
                                use_uniform_solar_input, &
                                lat_for_solar_input, lon_for_solar_input, &
                                rad_time_step, sw_rad_time_step, use_single_lw_sw_ts, &
                                nonzero_rad_flux_init, &
                                do_conserve_energy
!---------------------------------------------------------------------
!---- public data ----

!------- data structures -------

! radiation_diag_type stores data used by the radiation_diag_mod
! this is typically only used for the standalone radiation code

public radiation_diag_type
type radiation_diag_type
  type(sw_output_type), dimension(:), pointer :: Sw_output=>NULL()
  type(lw_output_type), dimension(:), pointer :: Lw_output=>NULL()
  type(lw_diagnostics_type)                   :: Lw_diagnostics
  type(lw_table_type),                pointer :: Lw_tables
  real, dimension(:,:,:,:,:), pointer :: emrndlw=>NULL(), &
                                         emmxolw=>NULL(), &
                                         cldsct=>NULL(),  &
                                         cldext=>NULL(),  &
                                         cldasymm=>NULL()
  real, dimension(:,:,:,:), pointer   :: camtsw=>NULL(), &
                                         crndlw=>NULL()
  real, dimension(:,:,:), pointer     :: cmxolw=>NULL(), &
                                         press=>NULL(),  &
                                         pflux=>NULL(),  &
                                         temp=>NULL(),   &
                                         rh2o=>NULL(),   &
                                         qo3=>NULL()
  real, dimension(:,:), pointer       :: asfc_vis_dir=>NULL(), &
                                         asfc_nir_dir=>NULL(), &
                                         asfc_vis_dif=>NULL(), &
                                         asfc_nir_dif=>NULL(), &
                                         cosz=>NULL(), &
                                         fracday=>NULL()
  real :: rrsun, rrvco2, rrvf11, rrvf12, rrvf113, rrvf22, rrvch4, rrvn2o
  real :: solar_constant
  integer :: indx_swaf, indx_lwaf
  logical :: do_totcld_forcing
  logical :: do_swaerosol, do_lwaerosol
  logical :: do_swaerosol_forcing, do_lwaerosol_forcing
  logical :: do_diurnal, do_annual, do_daily_mean
  logical :: do_radiation_diag = .false.
end type radiation_diag_type

!     Sw_output  sw_output_type variable containing shortwave radiation output data
!     Lw_output  lw_output_type variable containing longwave radiation output data
!     Lw_diagnostics  lw_diagnostics_type variable containing diagnostic
!                     longwave output used by the radiation diagnostics module
!     emmxolw    lw cloud emissivity for maximally overlapped clouds.
!                [ dimensionless ]
!     emrndlw    lw cloud emissivity for randomly overlapped clouds.
!                [ dimensionless ]
!     cmxolw     amount of maximally overlapped longwave cloud.
!                [ dimensionless ]
!     crndlw     amount of randomly overlapped longwave cloud.
!                [ dimensionless ]
!     camtsw     shortwave cloud amount [ dimensionless ]
!     cldext      cloud extinction coefficient [ km -1 ]  
!     cldsct      cloud scattering coefficient [ dimensionless ]
!     cldsasymm   cloud asymmetry factor  [ dimensionless ]

!---------------------------------------------------------------------
!---- private data ----
!-- for netcdf restart
type(restart_file_type), pointer, save :: Rad_restart => NULL()
type(restart_file_type), pointer, save :: Til_restart => NULL()
logical                                :: in_different_file = .false.
integer                                :: int_renormalize_sw_fluxes
integer                                :: int_do_clear_sky_pass

!--- data for concurrent radiation restarts ---
integer :: ido_conc_rad = 0
integer :: idonner_meso = 0
integer :: idoing_donner = 0
integer :: idoing_uw_conv = 0
type(restart_file_type), pointer, save :: Rad_restart_conc => NULL()
type(restart_file_type), pointer, save :: Til_restart_conc => NULL()
logical                                :: in_different_file_conc = .false.

type(radiation_flux_block_type) :: Restart

logical :: do_concurrent_radiation = .false.

!--- miscellaneous data ---

type(radiative_gases_type), save :: Rad_gases_tv
type(aerosolrad_control_type) ::  Aerosolrad_control
type(cloudrad_control_type) :: Cldrad_control
type(radiation_control_type) :: Rad_control
integer :: radiation_clock, longwave_clock, shortwave_clock

!    miscellaneous control variables:

integer   :: nt                       ! total no. of tracers
integer   :: ntp                      ! total no. of prognostic tracers

!---------------------------------------------------------------------
!    logical  flags.

logical ::  module_is_initialized = .false. ! module initialized?
logical ::  do_rad                          ! is this a radiation step ?
logical ::  do_lw_rad, do_sw_rad            ! is this a radiation step ?
logical ::  use_rad_date                    ! specify time of radiation
                                            ! independent of model time?

!---------------------------------------------------------------------
!    list of restart files readable by this module.
!
!                 sea_esf_rad.res:
!
!     version 1:  sea_esf_rad.res file version used initially in 
!                 AM2 model series (through galway code, AM2p8). this
!                 is the only version of sea_esf_rad.res ever produced.
!
!                 radiation_driver.res:
!
!     version 1:  not readable by this module.
!     version 2:  added cosine of zenith angle as an output to
!                 radiation_driver.res  (6/27/00)
!     version 3:  added restart variables needed when sw renormalization
!                 is active. (3/21/02)
!     version 4:  added longwave heating rate as separate output 
!                 variable, since it is needed as input to edt_mod
!                 and entrain_mod. (7/17/02)
!     version 5:  removed variables associated with the former 
!                 do_average namelist option (7/23/03)
!     version 6:  added writing of sw tropospheric fluxes (up and
!                 down) so that they are available for the renormal-
!                 ization case (developed by ds, 10/03; added to 
!                 trunk code 01/14/04).
!     version 7:  added swdn to saved variables (developed by slm 
!                 11/23/03, added to trunk code 01/14/04).
!     version 8:  includes additional sw fluxes at sfc, used with
!                 land model (11/13/03).
!     version 9:  consolidation of version 6 and version 8. (version 7
!                 replaced by version 8.)
!     version 10: adds 2 clr sky sw down diffuse and direct sfc flux
!                 diagnostic variables (10/18/04)
!     version 11: adds flux_sw_down_vis_clr diagnostic variable for use
!                 in assessing polar ice maintainability (6/19/07)
!---------------------------------------------------------------------
integer, dimension(10) :: restart_versions     = (/ 2, 3, 4, 5, 6,  &
                                                   7, 8, 9, 10, 11 /)
integer                :: vers ! version number of the restart file being read

!-----------------------------------------------------------------------
!    these arrays must be preserved across timesteps:
!
!    Rad_output is a rad_output_type variable with the following 
!    components:
!          tdt_rad        radiative (sw + lw) heating rate
!          flux_sw_surf   net (down-up) sw flux at surface
!          flux_sw_surf_dir   net (down-up) sw flux at surface
!          flux_sw_surf_refl_dir   dir sw flux reflected at surface
!          flux_sw_surf_dif   net (down-up) sw flux at surface
!          flux_sw_down_vis_dir  downward visible sw flux at surface
!          flux_sw_down_vis_dif  downward visible sw flux at surface
!          flux_sw_down_total_dir  downward total sw flux at surface
!          flux_sw_down_total_dif  downward total sw flux at surface
!          flux_sw_down_total_dir_clr  downward total direct sw flux at 
!                                      surface  (clear sky)
!          flux_sw_down_total_dif_clr  downward total diffuse sw flux 
!                                      at surface   (clear sky)
!          flux_sw_down_vis_clr  downward visible sw flux at surface
!                                       (clear sky)
!          flux_sw_vis    net visible sw flux at surface
!          flux_sw_vis_dir    net visible sw flux at surface
!          flux_sw_refl_vis_dir reflected direct visible sw flux at surface
!          flux_sw_vis_dif net visible sw flux at surface
!          flux_lw_surf   downward lw flux at surface
!          coszen_angle   cosine of the zenith angle (used for the 
!                         last radiation calculation)
!          tdt_rad_clr    net radiative heating rate in the absence of
!                         cloud
!          tdtsw          shortwave heating rate
!          tdtsw_clr      shortwave heating rate in he absence of cloud
!          tdtlw_clr       longwave heating rate in he absence of cloud
!          tdtlw          longwave heating rate
!          ufsw          upward sw flux
!          dfsw          downward sw flux
!          ufsw_clr      upward sw flux
!          dfsw_clr      downward sw flux
!          flxnet        net lw flux
!          flxnetcf      net lw flux, cloud free
!          extinction    SW extinction (band 4 centered on 1 micron) for volcanoes

type(rad_output_type),save          ::  Rad_output

!-----------------------------------------------------------------------
!    time-step-related constants
 
integer    :: lwrad_alarm    !  time interval until the next radiation 
                             !  calculation (seconds)
integer    :: swrad_alarm    !  time interval until the next radiation 
                             !  calculation (seconds)
type(time_type) :: Rad_time  !  time at which the climatologically-
                             !  determined, time-varying input fields to
                             !  radiation should apply 
                             !  [ time_type (days, seconds)]
integer    :: dt             !  physics time step (frequency of calling 
                             !  radiation_driver)  [ seconds ]
integer :: lw_rad_time_step  !  longwave radiative time step in seconds

!-----------------------------------------------------------------------
!    timing clocks       

integer                      :: misc_clock, clouds_clock, calc_clock

!--------------------------------------------------------------------
! miscellaneous variables and indices

integer        ::  ks         !  model grid coordinate of top level
                              !  at which radiation is calculated 
integer        ::  ke         !  model grid coordinate of bottommost
                              !  level at which radiation is calculated

integer        ::  ksrad=1    !  always set to 1
integer        ::  kerad      !  number of layers in radiation grid

real           ::  rh2o_lower_limit_seaesf=2.0E-07
                              !  smallest value of h2o mixing ratio 
                              !  allowed with sea_esf_rad package
real           ::  rh2o_lower_limit
                              !  smallest value of h2o mixing ratio 
                              !  allowed in the current experiment
real           ::  temp_lower_limit=100.0  ! [ K ]
                              !  smallest value of temperature      
                              !  allowed in the current experiment
real           ::  temp_upper_limit=370.00  ! [ K ]
                              !  largest value of temperature 
                              !  allowed in the current experiment

real           ::  surf_flx_init=50.0  ! [w / m^2 ]
                              !  value to which surface lw and sw fluxes
                              !  are set in the absence of a .res file
                              !  containing them

real           ::  coszen_angle_init=0.50
                              !  value to which cosine of zenith angle  
                              !  is set in the absence of a .res file
                              !  containing it

real           ::  log_p_at_top=2.0
                              !  assumed value of ln of ratio of pres-
                              !  sure at flux level 2 to that at model
                              !  top (needed for deltaz calculation,
                              !  is infinite for model top at p = 0.0,
                              !  this value is used to give a reasonable
                              !  deltaz)
real,parameter ::  D608 = (RVGAS-RDGAS)/RDGAS
                              !  virtual temperature factor  
real,parameter ::  D622 = RDGAS/RVGAS
                              ! ratio of gas constants - dry air to 
                              ! water vapor
real,parameter ::  D378 = 1.0 - D622  
                              ! 1 - gas constant ratio
integer :: id, jd

! number of stochastic columns if do_stochastic_clouds = true
integer :: num_lw_stoch_columns=1, num_sw_stoch_columns=1
! number of ica profiles if do_ica_calcs = true
integer :: num_lw_ica_profiles=1, num_sw_ica_profiles=1
! flag for stochastic/ica calcs
! do_stochastic => 1; do_ica_calcs => 2
integer :: flag_stoch

! <DATASET NAME="Restart file">
! A restart data set called radiation_driver.res(.nc) saves the
!     global fields for the current radiative tendency, net shortwave
!     surface flux, downward longwave surface flux, and cosine of the
!     zenith angle. If the namelist variable do_average=true,
!     then additional time averaged global data is written.
!     If the restart file is not present when initializing then the
!     radiative tendency is set to zero, the SW and LW surface fluxes
!     to 50 watts/m2, and the cosine of the zenith angle to 0.50.
!     Since radiation is usually computed on the first time step when
!     restarting, these values may have little or no effect.  If the
!     restart file is not present time average data is also set to zero.
! </DATASET>
!<REFERENCE> For a specific list of radiation references see the
!     longwave and shortwave documentation.</REFERENCE>
!---------------------------------------------------------------------


                         contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!######################################################################
! <SUBROUTINE NAME="radiation_driver_init">
!  <OVERVIEW>
!   radiation_driver_init is the constructor for radiation_driver_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   radiation_driver_init is the constructor for radiation_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiation_driver_init (lonb, latb, pref, axes, Time, &
!                                  aerosol_names)
!
!  </TEMPLATE>
!  <IN NAME="lonb" TYPE="real">
!    lonb      Longitude in radians for all (i.e., the global size)
!              grid box corners, the size of lonb should be one more
!              than the number of points along the x-axis and y-axis.
!                 [real, dimension(:,:)]
!  </IN>
!  <IN NAME="latb" TYPE="real">
!    latb      Latitude in radians for all (i.e., the global size)
!              grid box corners, the size of latb should be one more
!              than the number of latitude points along the x-axis and y-axis.
!                 [real, dimension(:,:)]
!  </IN>
!  <IN NAME="pref" TYPE="real">
!    pref      Two reference profiles of pressure at full model levels
!              plus the surface (nlev+1). The first profile assumes a surface
!              pressure of 101325 pa, and the second profile assumes 
!              81060 pa.  [real, dimension(nlev+1,2)]
!  </IN>
!  <IN NAME="axes" TYPE="integer">
!    axes      The axis indices that are returned by previous calls to
!              diag_axis_init. The values of this array correspond to the
!              x, y, full (p)level, and half (p)level axes. These are the
!              axes that diagnostic fields are output on.
!                 [integer, dimension(4)]
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!    Time      The current time.  [time_type]
!  </IN>
!   <ERROR MSG="must have two reference pressure profile" STATUS="FATAL">
!     The input argument pref must have a second dimension size of 2.
!   </ERROR>
!   <ERROR MSG="restart version ## cannot be read by this module version" STATUS="FATAL">
!     You have attempted to read a radiation_driver.res file with either
!       no restart version number or an incorrect restart version number.
!   </ERROR>
!
!   <NOTE>
!    radiation time step has changed, next radiation time also changed
!       The radiation time step from the namelist input did not match
!       the radiation time step from the radiation restart file.
!       The next time for radiation will be adjusted for the new  namelist
!       input) value.
!   </NOTE>
! </SUBROUTINE>
!
subroutine radiation_driver_init (Time, lonb, latb, axes, &
                                  Exch_ctrl, Atm_block, Radiation, Rad_flux)

!---------------------------------------------------------------------
!   radiation_driver_init is the constructor for radiation_driver_mod.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
type(time_type),             intent(in)    :: Time
real, dimension(:,:),        intent(in)    :: lonb, latb
integer, dimension(4),       intent(in)    :: axes
type(exchange_control_type), intent(inout) :: Exch_ctrl
type(block_control_type),    intent(in)    :: Atm_block
type(radiation_type),        intent(inout) :: Radiation
type(radiation_flux_type),   intent(inout) :: Rad_flux(:)

!----------------------------------------------------------------------
!   intent(in) variables:
!
!       lonb           2d array of model longitudes on cell corners 
!                      [ radians ]
!       latb           2d array of model latitudes at cell corners 
!                      [ radians ]
!       axes           diagnostic variable axes
!       Time           current time [time_type(days, seconds)]
!
!---------------------------------------------------------------------
!   local variables

      character(len=64), dimension(:), pointer :: aerosol_names => NULL()
      character(len=64), dimension(:), pointer :: aerosol_family_names => NULL()

      integer          ::  radiation_init_clock, &
                           radiative_gases_init_clock, &
                           cloud_spec_init_clock, &
                           aerosol_init_clock

      integer           ::   unit, io, ierr, logunit, outunit, n
      integer           ::   ibs, ibe, jbs, jbe
      integer           ::   kmax 
      integer           ::   nyr, nv, nband
      integer           ::   yr, month, year, dum
      integer           ::   ico2
      integer           ::   num_sw_bands, num_lw_bands


!---------------------------------------------------------------------
!   local variables
! 
!        unit    io unit number for namelist file
!        io      error status returned from io operation
!        ierr    error code
!        id      number of grid points in x direction (on processor)
!        jd      number of grid points in y direction (on processor)
!        kmax    number of model layers
!
!   aerosol_names        = names associated with activated aerosol species
!   aerosol_family_names = names associated with activated aerosol families
!                
!---------------------------------------------------------------------

      
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
      call time_manager_init
      call sat_vapor_pres_init
      call constants_init
      call tracer_manager_init

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=radiation_driver_nml, iostat=io)
      ierr = check_nml_error(io,'radiation_driver_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=radiation_driver_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'radiation_driver_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!--------------------------------------------------------------------
!    make sure other namelist variables are consistent with 
!    doing_data_override. Validate here to prevent potentially mis-
!    leading values from going into the stdlog file.
!--------------------------------------------------------------------
      if (.not. doing_data_override) then
        overriding_temps   = .false.
        overriding_sphum   = .false.
        overriding_albedo  = .false.
        overriding_clouds  = .false.
        overriding_aerosol = .false.
      endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
           write (logunit, nml=radiation_driver_nml)

!--------------------------------------------------------------------
!    define the model dimensions on the local processor.
!---------------------------------------------------------------------
      id = size(lonb,1)-1
      jd = size(latb,2)-1
      kmax  = Atm_block%npz
      call get_number_tracers (MODEL_ATMOS, num_tracers=nt, &
                                            num_prog=ntp)

!---------------------------------------------------------------------
!    verify that a valid radiation time step has been specified.
!---------------------------------------------------------------------
      if (rad_time_step <= 0) then
        call error_mesg ('radiation_driver_mod', &
            ' radiation timestep must be set to a positive integer', &
              FATAL)
      endif
      if (.not. use_single_lw_sw_ts) then
        if (sw_rad_time_step <= 0) then
          call error_mesg ('radiation_driver_mod', &
           ' sw radiation timestep must be set to a positive integer', &
              FATAL)
        endif
      endif

      if (use_single_lw_sw_ts .and. (sw_rad_time_step /= 0.0 .and. &
          sw_rad_time_step /= rad_time_step) ) then
        call error_mesg ('radiation_driver_mod', &
         'to avoid confusion, sw_rad_time_step must either remain at &
                  &default value of 0.0, or be same as rad_time_step &
                        &when use_single_lw_sw_ts is .true.', FATAL)
      endif
      if (use_single_lw_sw_ts) then
        sw_rad_time_step = rad_time_step
      endif
      lw_rad_time_step = rad_time_step

      if (MOD(INT(SECONDS_PER_DAY), lw_rad_time_step) /= 0) then
        call error_mesg ('radiation_drive_mod', &
             'lw radiation timestep currently restricted to be an &
                       &integral factor of seconds in a day', FATAL)
      endif
      if (MOD(INT(SECONDS_PER_DAY), sw_rad_time_step) /= 0) then
        call error_mesg ('radiation_driver_mod', &
             'sw radiation timestep currently restricted to be an &
                       &integral factor of seconds in a day', FATAL)
      endif

!---------------------------------------------------------------------
!    set logical variable defining the radiation scheme desired from the
!    namelist-input character string. set lower limit to water vapor 
!    mixing ratio that the radiation code will see, to assure keeping 
!    within radiation lookup tables. exit if value is invalid.
!---------------------------------------------------------------------
      rh2o_lower_limit = rh2o_lower_limit_seaesf

!---------------------------------------------------------------------
!    set control variable indicating whether ozone effects are to be
!    included in the radiative calculation
!---------------------------------------------------------------------
!BW   Rad_control%do_o3 = do_o3 

!---------------------------------------------------------------------
!    set control variables indicating whether the effects of other
!    radiative gases are to be included in the longwave radiation
!    calculation
!---------------------------------------------------------------------
!BW   Rad_control%do_ch4_lw = do_ch4_lw ! moved to radiative gases
!BW   Rad_control%do_n2o_lw = do_n2o_lw 
!BW   Rad_control%do_co2_lw = do_co2_lw 
!BW   Rad_control%do_cfc_lw = do_cfc_lw 
      
!---------------------------------------------------------------------
!    stop execution if overriding of aerosol data has been requested.
!    code to do so has not yet been written.
!---------------------------------------------------------------------
      if (overriding_aerosol) then
        call error_mesg ('radiation_driver_mod', &
                'overriding of aerosol data not yet implemented', FATAL)
      endif

!RSH:
!RSH    if use_co2_tracer_field is .true., verify here that there is
!RSH   in fact a co2 field included within the tracer array. if not,
!RSH   call error_mesg and abort execution.
!RSH
      if(use_co2_tracer_field) then
         ico2 = get_tracer_index(MODEL_ATMOS, 'co2')
         if(ico2 == NO_TRACER) then
            call error_mesg('radiation_driver_mod', &
                 'co2 must be present as a tracer when use_co2_tracer_field is .true.', FATAL)
         endif
      endif

!--------------------------------------------------------------------
!    set logical variables defining how the solar zenith angle is to
!    be  defined from the namelist-input character string.  exit if the
!    character string is invalid.
!--------------------------------------------------------------------
      if (zenith_spec == 'diurnally_varying') then
        Rad_control%do_diurnal = .true.
        Rad_control%do_annual = .false.
        Rad_control%do_daily_mean = .false.
      else if (zenith_spec == 'daily_mean') then
        Rad_control%do_diurnal = .false.
        Rad_control%do_annual = .false.
        Rad_control%do_daily_mean = .true.
      else if (zenith_spec == 'annual_mean') then
        Rad_control%do_diurnal = .false.
        Rad_control%do_annual = .true.
        Rad_control%do_daily_mean = .false.
      else
        call error_mesg ('radiation_driver_mod', &    
            'string provided for zenith_spec is invalid', FATAL)
      endif

!--------------------------------------------------------------------
!    check if spacially-uniform solar input has been requested. if it
!    has, verify that the requested lat and lon are valid, and convert
!    them to radians.
!--------------------------------------------------------------------
      if (use_uniform_solar_input) then
        if (lat_for_solar_input < -90. .or. &
            lat_for_solar_input >  90. ) then
          call error_mesg ('radiation_driver_mod', &
            'specified latitude for uniform solar input is invalid', &
                                                            FATAL)
        else
          lat_for_solar_input = lat_for_solar_input/RADIAN
        endif
        if (lon_for_solar_input < 0. .or. &
            lon_for_solar_input > 360. ) then
          call error_mesg ('radiation_driver_mod', &
            'specified longitude for uniform solar input is invalid', &
                                                             FATAL)
        else
          lon_for_solar_input = lon_for_solar_input/RADIAN
        endif
      endif

!---------------------------------------------------------------------
!    can only renormalize shortwave fluxes when diurnally_varying
!    radiation is used.
!---------------------------------------------------------------------
     if (renormalize_sw_fluxes .and. .not. Rad_control%do_diurnal) then
       call error_mesg ('radiation_driver_mod',  &
       ' can only renormalize sw fluxes when using diurnally-varying'//&
                       ' solar radiation', FATAL)
     endif


!----------------------------------------------------------------------
!    store the controls for solar interpolator and all step diagnostics
!----------------------------------------------------------------------
      Rad_control%renormalize_sw_fluxes = renormalize_sw_fluxes

!---------------------------------------------------------------------
!    define the starting and ending vertical indices of the radiation
!    grid. 
!---------------------------------------------------------------------
      ks = 1
      ke = kmax
      kerad = kmax
       
!---------------------------------------------------------------------
!    be sure both reference pressure profiles have been provided.
!----------------------------------------------------------------------
      if (size(Radiation%glbl_qty%pref,2) /= 2)    &
        call error_mesg ('radiation_driver_mod', &
         'must provide two reference pressure profiles (pref).', FATAL)

!--------------------------------------------------------------------
!    do the initialization specific to the sea_esf_rad radiation
!    package.
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    define control variables indicating whether the clear-sky forcing
!    should be calculated. set a flag to indicate that the variable
!    has been defined.
!---------------------------------------------------------------------
      Rad_control%do_totcld_forcing = do_clear_sky_pass

!---------------------------------------------------------------------
!    allocate space for module variables to contain values which must
!    be saved between timesteps (these are used on every timestep,
!    but only calculated on radiation steps).
!---------------------------------------------------------------------
      call Rad_output%alloc (id, jd, kmax, Rad_control%do_totcld_forcing) 

!-----------------------------------------------------------------------
!    should stardard radiation restart files be written?
!    they are only needed when when the model does not terminate
!    at a radiation time step 
!---------------------------------------------------------------------

   Rad_control%using_restart_file = using_restart_file

   if  (Rad_control%using_restart_file) then

!----------------------------------------------------------------------
!    Register fields to be written out to restart file.
     call rad_driver_register_restart('radiation_driver.res.nc')

!-----------------------------------------------------------------------
!    if a valid restart file exists, call read_restart_file to read it.
!-----------------------------------------------------------------------
        if ( file_exist('INPUT/radiation_driver.res.nc')) then
          call restore_state(Rad_restart)
          if(in_different_file) call restore_state(Til_restart)
        else if ( file_exist('INPUT/radiation_driver.res')) then
          call error_mesg ('radiation_driver_mod', &
              'Native restarts no longer supported', FATAL)
!----------------------------------------------------------------------
!    if no restart file is present, initialize the needed fields until
!    the radiation package may be called. initial surface flux is set 
!    to 100 wm-2, and is only used for initial guess of sea ice temp.
!    set rad_alarm to be 1 second from now, ie., on the first step of 
!    the job.
!-----------------------------------------------------------------------
        else
          lwrad_alarm                = 1
          swrad_alarm                = 1
          if (mpp_pe() == mpp_root_pe() ) then
          call error_mesg ('radiation_driver_mod', &
           'radiation to be calculated on first step: no restart file&
                                                 & present', NOTE)
          endif
          call Rad_output%initvalues
         !call rad_output_init (Rad_output)
          if (mpp_pe() == mpp_root_pe() ) then
            call error_mesg ('radiation_driver_mod', &
           'no acceptable radiation restart file present; therefore'//&
           ' will initialize input fields', NOTE)
          endif
        endif

!---------------------------------------------------------------------
!    if not using restart file, then initialize fields it would contain.
!    it is the responsibility of the user to assure restart is on a
!    radiation timestep so that restart seamlessness is maintained. if
!    restart is done on a non-radiation step, restart seamlessness will 
!    be lost if a restart file is not available.
!---------------------------------------------------------------------
   else  ! (using_restart_file)
     lwrad_alarm                = 1
     swrad_alarm                = 1
     if (mpp_pe() == mpp_root_pe() ) then
       call error_mesg ('radiation_driver_mod', &
          'radiation to be calculated on first step: user asserts that&
           & this is a scheduled radiation step;  if it is not, &
                           &restart seamlessness will be lost ', NOTE)
     endif
     call Rad_output%initvalues
    !call rad_output_init (Rad_output)
   endif ! (using_restart_file)

!--------------------------------------------------------------------
!    do the initialization for the radiation package.
!--------------------------------------------------------------------
      radiation_clock       =       &
                mpp_clock_id( '   Physics_down: Radiation',    &
                   grain=CLOCK_MODULE_DRIVER )
      cloud_spec_init_clock       =       &
        mpp_clock_id( '   Physics_driver_init: Cloud spec: Initialization', &
                       grain=CLOCK_MODULE_DRIVER )
      aerosol_init_clock       =       &
        mpp_clock_id( '   Physics_driver_init: Aerosol: Initialization', &
                       grain=CLOCK_MODULE_DRIVER )
      radiative_gases_init_clock       =       &
        mpp_clock_id( '   Physics_driver_init: Radiative gases: Initialization', &
                       grain=CLOCK_MODULE_DRIVER )
      radiation_init_clock       =       &
        mpp_clock_id( '   Physics_driver_init: Radiation: Initialization', &
                       grain=CLOCK_MODULE_DRIVER )

!---------------------------------------------------------------------
!    initialize clocks to time LW/SW core routines
!---------------------------------------------------------------------
      longwave_clock =      &
                  mpp_clock_id ('   Physics_down: Radiation: lw', &
                        grain=CLOCK_ROUTINE)
      shortwave_clock =     &
                  mpp_clock_id ('   Physics_down: Radiation: sw', &
                        grain=CLOCK_ROUTINE)

!-----------------------------------------------------------------------
!    initialize radiative_gases_mod.
!-----------------------------------------------------------------------
      call mpp_clock_begin ( radiative_gases_init_clock )
      call radiative_gases_init (lw_rad_time_step, &
                                 Radiation%glbl_qty%pref, latb, lonb)
      call mpp_clock_end ( radiative_gases_init_clock )

!---------------------------------------------------------------------
!    if h2o effects are not to be included in the longwave radiative
!    calculations, then set the lower limit for h2o to zero
!---------------------------------------------------------------------
      if (.not. get_longwave_gas_flag('h2o')) then
        rh2o_lower_limit = 0.0
      endif

!---------------------------------------------------------------------
!    initialize the modules that are accessed from radiation_driver_mod.
!---------------------------------------------------------------------
      call longwave_driver_init  (Radiation%glbl_qty%pref(ks:ke+1,:))
      call shortwave_driver_init (Rad_control)

!---------------------------------------------------------------------
!    get the number of SW and LW bands
!    this will be needed by clouds and aerosols
!---------------------------------------------------------------------
      call shortwave_number_of_bands (num_sw_bands)
      call longwave_number_of_bands  (num_lw_bands)

!-----------------------------------------------------------------------
!    initialize clouds.     
!-----------------------------------------------------------------------
      call mpp_clock_begin ( cloud_spec_init_clock )
      call cloudrad_driver_init (Time, rad_time_step, &
                                 lonb, latb, axes, &
                                 Radiation%glbl_qty%pref, &
                                 num_sw_bands, num_lw_bands, &
                                 Cldrad_control, &
                                 Exch_ctrl)
      call mpp_clock_end ( cloud_spec_init_clock )

!-----------------------------------------------------------------------
!    initialize aerosols
!-----------------------------------------------------------------------
      call mpp_clock_begin ( aerosol_init_clock )
      call aerosolrad_driver_init (lonb, latb, kmax,   &
                                   num_sw_bands, num_lw_bands, &
                                   Aerosolrad_control, &
                                   aerosol_names, aerosol_family_names)
      call mpp_clock_end ( aerosol_init_clock )

!---------------------------------------------------------------------
!    initialize the diagnostics and solar interpolator module
!---------------------------------------------------------------------
      call radiation_driver_diag_init (Time, id, jd, kmax, axes, &
                                       Rad_control, Aerosolrad_control)

!-----------------------------------------------------------------------

      call rad_output_file_init    (axes, Time, aerosol_names, &
                                    aerosol_family_names, &
                                    Rad_control, Aerosolrad_control)

!--------------------------------------------------------------------
!    initialize the astronomy_package.
!--------------------------------------------------------------------
      if (Rad_control%do_annual) then
        call astronomy_init (latb, lonb)
      else
        call astronomy_init
      endif

!-----------------------------------------------------------------------
!    check if optional radiative date should be used.
!-----------------------------------------------------------------------
      if (rad_date(1) > 1900 .and.                        &
          rad_date(2) >   0  .and. rad_date(2) < 13 .and. &
          rad_date(3) >   0  .and. rad_date(3) < 32 ) then
        use_rad_date = .true.
      else
        use_rad_date = .false.
      endif

!---------------------------------------------------------------------
!    Rad_flux type
!---------------------------------------------------------------------
      if (size(Rad_flux,1) .gt. 1) then 
        do_concurrent_radiation = .true.
      endif
      call alloc_radiation_flux_type (Rad_flux, nonzero_rad_flux_init, Atm_block)

!--------------------------------------------------------------------
!    initialize clocks to time portions of the code called from 
!    radiation_driver.
!--------------------------------------------------------------------
      misc_clock =    &
            mpp_clock_id ('   Physics_down: Radiation: misc', &
                grain = CLOCK_MODULE)
      clouds_clock =   &
            mpp_clock_id ('   Physics_down: Radiation: clds', &
               grain = CLOCK_MODULE)
      calc_clock =    &
            mpp_clock_id ('   Physics_down: Radiation: calc', &
                grain = CLOCK_MODULE)

!---------------------------------------------------------------------
!     if Rad_time is unchanging between timesteps, or the same day is being
!     repeated, switch to the alternative seed generation procedure to
!     assure unique temporal and spatial seeds for the stochastic cloud
!     parameterization.
!---------------------------------------------------------------------
! NEED TO MOVE ???
      if ( (rsd .or. use_rad_date)) then  
        Cldrad_control%use_temp_for_seed = .true.
        call error_mesg ('radiation_driver_init', &
             'Will use temp as basis for stochastic cloud seed; &
                                &Rad_time is not monotonic', NOTE)
      endif

!---------------------------------------------------------------------
!    verify that stochastic clouds have been activated if the COSP 
!    simulator output has been requested.
!---------------------------------------------------------------------
      if (Exch_ctrl%do_cosp .and.   &
          (.not. Cldrad_control%do_stochastic_clouds) ) then
        call error_mesg ('radiation_driver_init', &
         'cannot call COSP simulator unless stochastic clouds are &
           &activated (do_stochastic_clouds in strat_clouds_W_nml)', &
                                                                  FATAL)
      endif

!--------------------------------------------------------------------
!    return the potential number of stochastic columns.
!--------------------------------------------------------------------
      Exch_ctrl%ncol = Cldrad_control%num_sw_cloud_bands + Cldrad_control%num_lw_cloud_bands
      flag_stoch = 0

     if (Cldrad_control%do_stochastic_clouds) then
        num_lw_stoch_columns = Cldrad_control%num_lw_cloud_bands
        num_sw_stoch_columns = Cldrad_control%num_sw_cloud_bands
        flag_stoch = 1
     endif

!--------------------------------------------------------------------
!    set the number of ica profiles
!--------------------------------------------------------------------
     if (Cldrad_control%do_ica_calcs) then
        num_lw_ica_profiles = Cldrad_control%num_lw_cloud_bands
        num_sw_ica_profiles = Cldrad_control%num_sw_cloud_bands
        flag_stoch = 2
     endif

!---------------------------------------------------------------------
!    deallocate space for local pointers.
!---------------------------------------------------------------------
        deallocate (aerosol_names, aerosol_family_names)

!----------------------------------------------------------------------
! restart logic needs to be handled for many different cases
! 1) concurrent radiation where radiation restart data must
!    be read in when the model restarts
! 2) serial radiation where restart behavior is governed
!    by 'using_restart_file' namelist parameter (handled in radiation_driver_init)
!----------------------------------------------------------------------
   if (do_concurrent_radiation) then
     if (mpp_pe() == mpp_root_pe()) call error_mesg ('radiation_driver_mod:', &
        'concurrent radiation restart active', NOTE)
     call conc_rad_register_restart('conc_radiation_driver.res.nc', &
                                    Rad_flux(size(Rad_flux,1)), Exch_ctrl, Atm_block)
     if (file_exist('INPUT/conc_radiation_driver.res.nc')) then
       call restore_state(Rad_restart_conc)
       if (in_different_file_conc) call restore_state(Til_restart_conc)
     else
       call error_mesg ('radiation_driver_mod', 'restart file conc_radiation_driver.res.nc not found',NOTE)
     endif
100 FORMAT("CHECKSUM::",A32," = ",Z20)
     outunit = stdout()
     write(outunit,*) 'BEGIN CHECKSUM(radiation_driver_init - concurrent):: '
     write(outunit,100) 'tdt_rad                ', mpp_chksum(Restart%tdt_rad                )
     write(outunit,100) 'tdt_lw                 ', mpp_chksum(Restart%tdt_lw                 )
     write(outunit,100) 'flux_sw                ', mpp_chksum(Restart%flux_sw                )
     write(outunit,100) 'flux_sw_dir            ', mpp_chksum(Restart%flux_sw_dir            )
     write(outunit,100) 'flux_sw_dif            ', mpp_chksum(Restart%flux_sw_dif            )
     write(outunit,100) 'flux_sw_down_vis_dir   ', mpp_chksum(Restart%flux_sw_down_vis_dir   )
     write(outunit,100) 'flux_sw_down_vis_dif   ', mpp_chksum(Restart%flux_sw_down_vis_dif   )
     write(outunit,100) 'flux_sw_down_total_dir ', mpp_chksum(Restart%flux_sw_down_total_dir )
     write(outunit,100) 'flux_sw_down_total_dif ', mpp_chksum(Restart%flux_sw_down_total_dif )
     write(outunit,100) 'flux_sw_vis            ', mpp_chksum(Restart%flux_sw_vis            )
     write(outunit,100) 'flux_sw_vis_dir        ', mpp_chksum(Restart%flux_sw_vis_dir        )
     write(outunit,100) 'flux_sw_vis_dif        ', mpp_chksum(Restart%flux_sw_vis_dif        )
     write(outunit,100) 'flux_lw                ', mpp_chksum(Restart%flux_lw                )
     write(outunit,100) 'coszen                 ', mpp_chksum(Restart%coszen                 )
     write(outunit,100) 'extinction             ', mpp_chksum(Restart%extinction             )
     do n = 1,size(Rad_flux(size(Rad_flux,1))%block,1)
       ibs = Atm_block%ibs(n)-Atm_block%isc+1
       ibe = Atm_block%ibe(n)-Atm_block%isc+1
       jbs = Atm_block%jbs(n)-Atm_block%jsc+1
       jbe = Atm_block%jbe(n)-Atm_block%jsc+1

       Rad_flux(size(Rad_flux,1))%block(n)%tdt_rad                = Restart%tdt_rad(ibs:ibe,jbs:jbe,:)
       Rad_flux(size(Rad_flux,1))%block(n)%tdt_lw                 = Restart%tdt_lw(ibs:ibe,jbs:jbe,:)
       Rad_flux(size(Rad_flux,1))%block(n)%flux_sw                = Restart%flux_sw(ibs:ibe,jbs:jbe)
       Rad_flux(size(Rad_flux,1))%block(n)%flux_sw_dir            = Restart%flux_sw_dir(ibs:ibe,jbs:jbe)
       Rad_flux(size(Rad_flux,1))%block(n)%flux_sw_dif            = Restart%flux_sw_dif(ibs:ibe,jbs:jbe)
       Rad_flux(size(Rad_flux,1))%block(n)%flux_sw_down_vis_dir   = Restart%flux_sw_down_vis_dir(ibs:ibe,jbs:jbe)
       Rad_flux(size(Rad_flux,1))%block(n)%flux_sw_down_vis_dif   = Restart%flux_sw_down_vis_dif(ibs:ibe,jbs:jbe)
       Rad_flux(size(Rad_flux,1))%block(n)%flux_sw_down_total_dir = Restart%flux_sw_down_total_dir(ibs:ibe,jbs:jbe)
       Rad_flux(size(Rad_flux,1))%block(n)%flux_sw_down_total_dif = Restart%flux_sw_down_total_dif(ibs:ibe,jbs:jbe)
       Rad_flux(size(Rad_flux,1))%block(n)%flux_sw_vis            = Restart%flux_sw_vis(ibs:ibe,jbs:jbe)
       Rad_flux(size(Rad_flux,1))%block(n)%flux_sw_vis_dir        = Restart%flux_sw_vis_dir(ibs:ibe,jbs:jbe)
       Rad_flux(size(Rad_flux,1))%block(n)%flux_sw_vis_dif        = Restart%flux_sw_vis_dif(ibs:ibe,jbs:jbe)
       Rad_flux(size(Rad_flux,1))%block(n)%flux_lw                = Restart%flux_lw(ibs:ibe,jbs:jbe)
       Rad_flux(size(Rad_flux,1))%block(n)%coszen                 = Restart%coszen(ibs:ibe,jbs:jbe)
       Rad_flux(size(Rad_flux,1))%block(n)%extinction             = Restart%extinction(ibs:ibe,jbs:jbe,:)
     end do
   endif

!---run the time vary routines in order to pre-load interpolation data files
!rab   call aerosol_time_vary (Time, Aerosol_rad)
!rab   call radiative_gases_time_vary (Time, Radiation%glbl_qty%gavg_q, Rad_gases_tv)
!rab   call radiation_driver_time_vary (Time, Rad_gases_tv)

!---------------------------------------------------------------------
!    set flag to indicate that module has been successfully initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!--------------------------------------------------------------------

end subroutine radiation_driver_init


!######################################################################
 
subroutine radiation_driver_time_vary (Time, Time_next, gavg_rrv, &
                                       Rad_flux_control)
 
!---------------------------------------------------------------------
!    radiation_driver_time_vary calculates time-dependent, 
!    space-independent quantities needed within the modules of the 
!    radiation package.
!---------------------------------------------------------------------
 
type(time_type),         intent(in)     :: Time, Time_next
real, dimension(:),      intent(in)     :: gavg_rrv
type(radiation_flux_control_type), intent(inout) :: Rad_flux_control
!----------------------------------------------------------------------
!  Call define_rad_times to obtain the time to be used in the
!  radiation calculation (Rad_time) and to determine which, if any,
!  externally-supplied inputs to radiation_driver must be obtained on
!  this timestep.  Logical flags are returned indicating the need or 
!  lack of need for the aerosol fields, the cloud fields, the
!  radiative gas fields, and the basic atmospheric variable fields.
!----------------------------------------------------------------------
 
    call define_rad_times (Time, Time_next, Rad_time)

    if (do_rad) then

!----------------------------------------------------------------------
!    call radiative_gases_time_vary to retrieve appropriate gas fields from
!    the climatology, if that source for radiative gases is being used.
!----------------------------------------------------------------------
       call radiative_gases_time_vary (Rad_time, gavg_rrv, &
                                       Rad_control%do_lw_rad, Rad_gases_tv)

       if (Aerosolrad_control%do_aerosol) &
           call aerosolrad_driver_time_vary (Rad_time, Aerosolrad_control)

       call shortwave_driver_time_vary (Rad_time)
       call longwave_driver_time_vary  (Rad_gases_tv)
    endif

    Rad_flux_control%do_rad = do_rad

!----------------------------------------------------------------------
 
end subroutine radiation_driver_time_vary
 
!####################################################################

subroutine radiation_driver_endts (is, js)

integer, intent(in)  :: is,js

!---------------------------------------------------------------------

      if (do_rad) then
         call radiative_gases_endts (Rad_gases_tv)
         if (Aerosolrad_control%do_aerosol) call aerosolrad_driver_endts
         call longwave_driver_endts (Rad_gases_tv)
      endif

      call radiation_driver_diag_endts (Rad_control)

!---------------------------------------------------------------------
!    complete radiation step. if this was a radiation step, set the 
!    radiation alarm to go off rad_time_step seconds from now, and
!    set do_rad to false, so that radiation will not be calculated 
!    again until the alarm goes off.
!--------------------------------------------------------------------
      if (.not. always_calculate) then
        if (do_lw_rad) then
          lwrad_alarm = lwrad_alarm + lw_rad_time_step
          do_lw_rad = .false.
        endif
        if (do_sw_rad) then
          swrad_alarm = swrad_alarm + sw_rad_time_step
          do_sw_rad = .false.
        endif

        if (.not. do_lw_rad .and. .not. do_sw_rad)  then
          do_rad = .false.
        else
          do_rad = .true.
        endif

      endif  ! (always_calculate)

      Rad_control%do_lw_rad = do_lw_rad
      Rad_control%do_sw_rad = do_sw_rad

!---------------------------------------------------------------------
 
end subroutine radiation_driver_endts


!#####################################################################
! <SUBROUTINE NAME="radiation_driver">
!  <OVERVIEW>
!    radiation_driver adds the radiative heating rate to the temperature
!    tendency and obtains the radiative boundary fluxes and cosine of 
!    the solar zenith angle to be used in the other component models.
!  </OVERVIEW>
!  <DESCRIPTION>
!    radiation_driver adds the radiative heating rate to the temperature
!    tendency and obtains the radiative boundary fluxes and cosine of 
!    the solar zenith angle to be used in the other component models.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiation_driver (is, ie, js, je, Time, Time_next,  &
!                             lat, lon, Surface, Atmos_input, &
!                             Aerosol, Cld_spec, Rad_gases, &
!                             Lsc_microphys, Meso_microphys,    &
!                             Cell_microphys, Radiation)
!
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending i,j indices in global storage arrays
!  </IN> 
!  <IN NAME="Time" TYPE="time_type">
!   current model time 
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   The time used for diagnostic output
!  </IN>
!  <IN NAME="lon" TYPE="real">
!    lon        mean longitude (in radians) of all grid boxes processed by
!               this call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <IN NAME="lat" TYPE="real">
!    lat        mean latitude (in radians) of all grid boxes processed by this
!               call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <INOUT NAME="Surface" TYPE="surface_type">
!   Surface input data to radiation package
!  </INOUT>
!  <INOUT NAME="Atmos_input" TYPE="atmos_input_type">
!   Atmospheric input data to radiation package
!  </INOUT>
!  <INOUT NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol climatological input data to radiation package
!  </INOUT>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!   Cloud microphysical and physical parameters to radiation package, 
!                     contains var-
!                     iables defining the cloud distribution, passed 
!                     through to lower level routines
!  </INOUT>
!  <INOUT NAME="Rad_gases" TYPE="radiative_gases_type">
!   Radiative gases properties to radiation package, , contains var-
!                     iables defining the radiatively active gases, 
!                     passed through to lower level routines
!  </INOUT>
!  <INOUT NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale
!                      clouds
!  </INOUT>
! <ERROR MSG="radiation_driver_init must first be called" STATUS="FALTA">
! You have not called radiation_driver_init before calling
!       radiation_driver.
! </ERROR>
! <ERROR MSG="Time_next <= Time" STATUS="FALTA">
! Time arguments to radiation_driver are producing a time step <= 0.
!       Check that the time argumnets passed to the physics_driver are
!       correct.
! </ERROR>
! </SUBROUTINE>
!

subroutine radiation_driver (is, ie, js, je, npz,            &
                             Time, Time_next,                &
                             lat, lon,                       &
                             Radiation_control,              &
                             Radiation_input_block,          &
                             Radiation_glbl_qty,             &
                             frac_land, albedo,              &
                             albedo_vis_dir, albedo_nir_dir, &
                             albedo_vis_dif, albedo_nir_dif, &
                             t_surf_rad,                     &
                             Exch_ctrl,                      &
                             Rad_flux_block,                 &
                             Cosp_rad_control,               &
                             Cosp_rad_block,                 &
                             Moist_clouds_block,             &
                             cloudtemp,   cloudvapor,        &
                             aerosoltemp, aerosolvapor,      &
                             aerosolpress,                   &
                             Astronomy_inp, Rad_diag )

!---------------------------------------------------------------------
!    radiation_driver adds the radiative heating rate to the temperature
!    tendency and obtains the radiative boundary fluxes and cosine of 
!    the solar zenith angle to be used in the other component models.
!---------------------------------------------------------------------

integer,                 intent(in)         :: is, ie, js, je, npz
type(time_type),         intent(in)         :: Time, Time_next
real,dimension(:,:),     intent(in)         :: lat, lon

type(radiation_input_control_type), intent(in) :: Radiation_control
type(radiation_input_block_type),   intent(in) :: Radiation_input_block
type(radiation_input_glblqty_type), intent(in) :: Radiation_glbl_qty

real,dimension(:,:),     intent(in), target :: frac_land, &
                                               albedo, &
                                               albedo_vis_dir, albedo_nir_dir, &
                                               albedo_vis_dif, albedo_nir_dif
real,dimension(:,:),     intent(in)         :: t_surf_rad

type(exchange_control_type),        intent(in)    :: Exch_ctrl
type(radiation_flux_block_type),    intent(inout) :: Rad_flux_block
type(cosp_from_rad_control_type),   intent(inout) :: Cosp_rad_control
type(cosp_from_rad_block_type),     intent(inout) :: Cosp_rad_block
type(clouds_from_moist_block_type), intent(in)    :: Moist_clouds_block

! optional inputs typically used by the standalone radiation code

real, dimension(:,:,:),  intent(in), optional    :: cloudtemp,    &
                                                    cloudvapor, &
                                                    aerosoltemp, &
                                                    aerosolvapor, &
                                                    aerosolpress
type(astronomy_inp_type),  intent(inout), optional :: Astronomy_inp
type(radiation_diag_type), intent(out),   optional :: Rad_diag

!----------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      Time           current model time [ time_type (days, seconds) ] 
!      Time_next      time on next timestep, used as stamp for diagnos-
!                     tic output  [ time_type  (days, seconds) ]  
!      lat            latitude of model points  [ radians ]
!      lon            longitude of model points [ radians ]
!
!   intent(inout) variables:
!
!      Surface        surface_type structure, contains variables 
!                     defining the surface characteristics, including
!                     the following component referenced in this 
!                     routine:
!
!         asfc          surface albedo  [ dimensionless ]
!
!      Atmos_input    atmos_input_type structure, contains variables
!                     defining atmospheric state, including the follow-
!                     ing component referenced in this routine
!
!         tsfc          surface temperature [ deg K ]
!
!      Aerosol        aerosol_type structure, contains variables
!                     defining aerosol fields, passed through to
!                     lower level routines
!      Cld_spec       cld_specification_type structure, contains var-
!                     iables defining the cloud distribution, passed 
!                     through to lower level routines
!      Rad_gases      radiative_gases_type structure, contains var-
!                     iables defining the radiatively active gases, 
!                     passed through to lower level routines
!      Lsc_microphys  microphysics_type structure, contains variables
!                     describing the microphysical properties of the
!                     large-scale clouds, passed through to lower
!                     level routines
!
!   intent(inout), optional variables:
!
!        tdt_rad         radiative (sw + lw) heating rate
!                        [ deg K / sec ]
!        flux_sw_surf    net (down-up) sw surface flux 
!                        [ watts / m^^2 ]
!        flux_lw_surf    downward lw surface flux 
!                        [ watts / m^^2 ]
!        coszen_angle    cosine of the zenith angle which will be used 
!                        for the next ocean_albedo calculation 
!                        [ dimensionless ]
!        tdtlw           longwave heating rate
!                        [ deg K / sec ]
!      Astronomy_inp  astronomy_input_type structure, optionally used
!                     to input astronomical forcings, when it is desired
!                     to specify them rather than use astronomy_mod.
!                     Used in various standalone applications.
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      type(cld_specification_type)                   :: Cld_spec
      type(radiative_gases_type)                     :: Rad_gases
      type(atmos_input_type)                         :: Atmos_input
      type(surface_type)                             :: Surface
      type(aerosol_type)                             :: Aerosol
      type(microphysics_type)                        :: Model_microphys


      type(astronomy_type)               :: Astro, Astro_phys
      type(lw_output_type), dimension(Aerosolrad_control%indx_lwaf) :: Lw_output
      type(sw_output_type), dimension(Aerosolrad_control%indx_swaf) :: Sw_output
      type(lw_diagnostics_type)                         :: Lw_diagnostics
      type(aerosolrad_diag_type)     :: Aerosolrad_diags

      real, dimension (ie-is+1, je-js+1) :: flux_ratio, &
                                            lat_uniform, lon_uniform
      integer :: nextinct

! local arrays for aerosol properties

real, dimension(size(Radiation_input_block%t,1), &
                size(Radiation_input_block%t,2), &
                size(Radiation_input_block%t,3), &
                Aerosolrad_control%num_lw_aerosol_bands) :: aerooptdep, aerooptdep_volc
real, dimension(size(Radiation_input_block%t,1), &
                size(Radiation_input_block%t,2), &
                size(Radiation_input_block%t,3), &
                Aerosolrad_control%num_sw_aerosol_bands) :: aeroasymfac, aerosctopdep, aeroextopdep

! local arrays for longwave cloud properties

real, dimension(size(Radiation_input_block%t,1), &
                size(Radiation_input_block%t,2), &
                size(Radiation_input_block%t,3)) :: cmxolw

real, dimension(size(Radiation_input_block%t,1), &
                size(Radiation_input_block%t,2), &
                size(Radiation_input_block%t,3), &
                num_lw_stoch_columns) :: crndlw

real, dimension(size(Radiation_input_block%t,1), &
                size(Radiation_input_block%t,2), &
                size(Radiation_input_block%t,3), &
                Cldrad_control%num_lw_cloud_bands, &
                num_lw_ica_profiles) :: emrndlw, emmxolw

! local arrays for shortwave cloud properties

real, dimension(size(Radiation_input_block%t,1), &
                size(Radiation_input_block%t,2), &
                size(Radiation_input_block%t,3), &
                num_sw_stoch_columns) :: camtsw

real, dimension(size(Radiation_input_block%t,1), &
                size(Radiation_input_block%t,2), &
                size(Radiation_input_block%t,3), &
                Cldrad_control%num_sw_cloud_bands, &
                num_sw_ica_profiles) :: cldsctsw, cldextsw, cldasymmsw

! local pointers for input data
real, dimension(:),       pointer :: gavg_rrv
real, dimension(:,:,:),   pointer :: t, p_full, p_half, z_full, z_half, q
real, dimension(:,:,:,:), pointer :: r, rm

                
!-------------------------------------------------------------------
!   local variables:
!
!      Astro             astronomical properties on model grid, usually
!                        valid over radiation timestep
!                        [astronomy_type]
!      Astro_phys        astronomical properties on model grid, valid 
!                        over current physics timestep
!                        [astronomy_type]
!      Lw_output         sea longwave output fields on model grid,
!                        [lw_output_type]
!      Sw_output         esf shortwave output fields on model grid,
!                        [sw_output_type]
!      flux_ratio        value  used to renormalize sw fluxes and 
!                        heating rates to account for earth-sun motion
!                        during the radiation timestep
!
!----------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!----------------------------------------------------------------------

      if (.not. do_radiation) then
        Rad_flux_block%tdt_rad = 0.0
        Rad_flux_block%tdt_lw = 0.0
        Rad_flux_block%flux_sw = 0.0
        Rad_flux_block%flux_sw_dir = 0.0
        Rad_flux_block%flux_sw_dif = 0.0
        Rad_flux_block%flux_sw_down_vis_dir = 0.0
        Rad_flux_block%flux_sw_down_vis_dif = 0.0
        Rad_flux_block%flux_sw_down_total_dir = 0.0
        Rad_flux_block%flux_sw_down_total_dif = 0.0
        Rad_flux_block%flux_sw_vis = 0.0
        Rad_flux_block%flux_sw_vis_dir = 0.0
        Rad_flux_block%flux_sw_vis_dif = 0.0
        Rad_flux_block%flux_lw = 0.0
        Rad_flux_block%coszen  = 0.0
        Rad_flux_block%extinction = 0.0

        return
      endif

! input data
      t => Radiation_input_block%t
      r => Radiation_input_block%q
      q => Radiation_input_block%q(:,:,:,Radiation_control%sphum)
      p_full => Radiation_input_block%p_full
      p_half => Radiation_input_block%p_half
      z_full => Radiation_input_block%z_full
      z_half => Radiation_input_block%z_half
      gavg_rrv => Radiation_glbl_qty%gavg_q

!----------------------------------------------------------------------
!    prepare to calculate radiative fluxes and heating rates. obtain the 
!    needed atmospheric fields, and needed inputs from other physics 
!    modules.
!----------------------------------------------------------------------
      call mpp_clock_begin ( radiation_clock )

!---------------------------------------------------------------------
!    call define_surface to define a surface_type variable containing
!    the surface albedoes and land fractions for each grid box. this
!    variable must be provided on all timesteps for use in generating
!    netcdf output.
!---------------------------------------------------------------------
!BW   call define_surface (is, ie, js, je, albedo, albedo_vis_dir,  &
!BW                        albedo_nir_dir, albedo_vis_dif, &
!BW                        albedo_nir_dif, frac_land, Surface)
      Surface%land => frac_land
      Surface%asfc_vis_dir => albedo_vis_dir
      Surface%asfc_nir_dir => albedo_nir_dir
      Surface%asfc_vis_dif => albedo_vis_dif
      Surface%asfc_nir_dif => albedo_nir_dif


!---------------------------------------------------------------------
!    if the basic atmospheric input variables to the radiation package
!    are needed, pass the model pressure (p_full, p_half), temperature 
!    (t, t_surf_rad) and specific humidity (q) to subroutine
!    define_atmos_input_fields, which will put these fields and some
!    additional auxiliary fields into the form desired by the radiation
!    package and store them as components of the derived-type variable 
!    Atmos_input.
!---------------------------------------------------------------------
      if (do_rad) then
        call define_atmos_input_fields     &
                              (is, ie, js, je, p_full, p_half, t, q, &
                               t_surf_rad, r, gavg_rrv, Atmos_input, &
                               cloudtemp, cloudvapor, aerosoltemp,   &
                               aerosolvapor, aerosolpress)
      endif

!---------------------------------------------------------------------
!    if this is a radiation step, or if the astronomical inputs to
!    radiation (solar, cosz, fracday, rrsun) need to be obtained 
!    because of time averaging or renormalization, call 
!    obtain_astronomy_variables to do so.
!---------------------------------------------------------------------
      if (do_rad .or. renormalize_sw_fluxes .or.  present(Astronomy_inp)) then 
        if (use_uniform_solar_input) then
          if (present (Astronomy_inp)) then
            call error_mesg ('radiation_driver_mod', &
              'cannot specify both use_uniform_solar_input AND use&
              & Astronomy_inp to specify astronomical variables', &
                                                               FATAL)
          endif
          lat_uniform(:,:) = lat_for_solar_input
          lon_uniform(:,:) = lon_for_solar_input
          call obtain_astronomy_variables (is, ie, js, je,  &
                                           lat_uniform, lon_uniform,  &
                                           Astro, Astro_phys)
        else
          if (present (Astronomy_inp)) then
            Rad_control%do_diurnal = .false.
            Rad_control%do_annual = .false.
            Rad_control%do_daily_mean = .false.
          endif
          call obtain_astronomy_variables (is, ie, js, je, lat, lon,  &
                                           Astro, Astro_phys, &
                                           Astronomy_inp =  &
                                                          Astronomy_inp)
        endif
      endif

!---------------------------------------------------------------------
!    if the radiative gases are needed, call define_radiative_gases to 
!    obtain the values to be used for the radiatively-active gases and 
!    place them in radiative_gases_type derived-type variable Rad_gases.
!---------------------------------------------------------------------
      if (do_rad) then

!--------------------------------------------------------------------
!    fill the contents of the radiative_gases_type variable which
!    will be passed to the radiation package. 
!---------------------------------------------------------------------
        Rad_gases = Rad_gases_tv
        call define_radiative_gases (is, ie, js, je, Rad_time, lat, &
                                     Atmos_input%pflux, r, Rad_gases)
      endif

!---------------------------------------------------------------------
!    if the aerosol fields are needed as input to the radiation_package,
!    call aerosol_driver to access the aerosol data and place it into 
!    an aerosol_type derived-type variable Aerosol.
!---------------------------------------------------------------------
      call mpp_clock_begin (misc_clock)
      if (do_rad .and. Aerosolrad_control%do_aerosol) then
        call aerosolrad_driver (is, ie, js, je, Rad_time, &
                                Astro%fracday, &
                                Atmos_input%phalf, &
                                Atmos_input%pflux, &
                                Atmos_input%aerosolrelhum, &
                                Atmos_input%deltaz, r, &
                                Rad_control%do_cmip_sw_diagnostics, &
                                Aerosolrad_control, &
                                Aerosolrad_diags, Aerosol, &
                                Rad_output%extinction(is:ie,js:je,:), &
                                aerooptdep, aerooptdep_volc, &
                                aeroasymfac, aerosctopdep, aeroextopdep)
      endif
      call mpp_clock_end (misc_clock)

!---------------------------------------------------------------------
!    if the cloud fields are needed, call cloudrad_driver     
!    to retrieve microphysical data which is returned in 
!    microphysics_type variable Model_microphys.
!    also returned are cloud-radiative properties.
!---------------------------------------------------------------------
      call mpp_clock_begin (clouds_clock)
      if (do_rad) then
         call cloudrad_driver (is, ie, js, je,               &
                               Time, Time_next, Rad_time,    &
                               lat, Surface%land, Astro%cosz, &
                               z_half, z_full, r, Atmos_input%tsfc, &
                               Atmos_input%press, Atmos_input%pflux, &
                               Atmos_input%temp, Atmos_input%deltaz, &
                               Atmos_input%cloudtemp, Atmos_input%cloudvapor, &
                               Atmos_input%clouddeltaz, &
                               Cldrad_control, &
                               Aerosol, Cld_spec, Model_microphys,  &
                               Moist_clouds_block, &
                               crndlw, cmxolw, emrndlw, emmxolw, &
                               camtsw, cldsctsw, cldextsw, cldasymmsw )
      endif ! (do_rad)
      call mpp_clock_end (clouds_clock)

!--------------------------------------------------------------------
!    on radiation timesteps, call radiation_calc to determine new radia-
!    tive fluxes and heating rates.
!---------------------------------------------------------------------
      call mpp_clock_begin (calc_clock)
      if (do_rad) then
        call radiation_calc (is, ie, js, je, lat, lon, &
                             Atmos_input%press, Atmos_input%pflux,  &
                             Atmos_input%temp,  Atmos_input%tflux,  &
                             Atmos_input%rh2o,  Atmos_input%deltaz, &
                             Atmos_input%tsfc,  &
                             Surface%asfc_vis_dir, Surface%asfc_nir_dir, &
                             Surface%asfc_vis_dif, Surface%asfc_nir_dif, &
                             Astro, Rad_gases, &
                             aerooptdep, aerooptdep_volc, &
                             aeroasymfac, aerosctopdep, aeroextopdep, &
                             crndlw, cmxolw, emrndlw, emmxolw, &
                             camtsw, cldsctsw, cldextsw, cldasymmsw, &
                             Rad_output, Lw_output, Sw_output, Lw_diagnostics)
      endif
      call mpp_clock_end (calc_clock)

!-------------------------------------------------------------------
!    on all timesteps, call update_rad_fields to update the temperature 
!    tendency and define the fluxes needed by other component models.
!    if the shortwave fluxes are to be renormalized because of the 
!    change in zenith angle since the last radiation timestep, that also
!    is done in this subroutine. 
!-------------------------------------------------------------------
      call mpp_clock_begin (misc_clock)
      call update_rad_fields (is, ie, js, je, Time_next, Astro, &
                              Astro_phys, Rad_control, &
                              Sw_output, Rad_output, flux_ratio)
                                

!-------------------------------------------------------------------
!    call produce_radiation_diagnostics to produce radiation 
!    diagnostics, both fields and integrals.
!-------------------------------------------------------------------
      call produce_radiation_diagnostics        &
                          (is, ie, js, je, Time_next, Time, lat, &
                           Radiation_glbl_qty%atm_mass, Atmos_input%tsfc, &
                           Atmos_input%pflux, p_half,  &
                      !BW  Atmos_input%pflux, Atmos_input%phalf,  &
                           Surface%asfc_vis_dir, Surface%asfc_nir_dir, &
                           Surface%asfc_vis_dif, Surface%asfc_nir_dif, &
                           flux_ratio,  Astro, Astro_phys, &
                           Rad_output, Rad_gases, Rad_control, &
                           Lw_output=Lw_output,&
                           Sw_output=Sw_output)

!---------------------------------------------------------------------
!    call write_rad_output_file to produce a netcdf output file of 
!    radiation-package-relevant variables. note that this is called
!    only on radiation steps, so that the effects of sw renormalization
!    will not be seen in the variables of the data file written by
!    write_rad_output_file.
!---------------------------------------------------------------------
      if (do_lw_rad .and. do_sw_rad) then
        if (Aerosolrad_control%do_aerosol) then
          call write_rad_output_file (is, ie, js, je,   &
                                      Atmos_input%tsfc, &
                                      Surface%asfc_vis_dir, &
                                      Surface%asfc_nir_dir, &
                                      Surface%asfc_vis_dif, &
                                      Surface%asfc_nir_dif, &
                                      Atmos_input%press, Atmos_input%pflux, &
                                      Atmos_input%phalf, Atmos_input%temp,  &
                                      Atmos_input%rh2o, Rad_gases%qo3, &
                                      Atmos_input%deltaz, &
                                      Rad_control, Aerosolrad_control, &
                                      Rad_output, Sw_output(1),  &
                                      Lw_output(1), Cld_spec, & 
                                      Time_next, Time, &
                                      Aerosol=Aerosol, &
                                     !Aerosol_props=Aerosol_props, &
                                      Aerosolrad_diags=Aerosolrad_diags)
        else
          call write_rad_output_file (is, ie, js, je,  &
                                      Atmos_input%tsfc, &
                                      Surface%asfc_vis_dir, &
                                      Surface%asfc_nir_dir, &
                                      Surface%asfc_vis_dif, &
                                      Surface%asfc_nir_dif, &
                                      Atmos_input%press, Atmos_input%pflux, &
                                      Atmos_input%phalf, Atmos_input%temp,  &
                                      Atmos_input%rh2o, Rad_gases%qo3, &
                                      Atmos_input%deltaz, &
                                      Rad_control, Aerosolrad_control, &
                                      Rad_output, Sw_output(1),   &
                                      Lw_output(1), Cld_spec, &
                                      Time_next, Time)
        endif
      endif ! (do_rad)


!---------------------------------------------------------------------
!    allocate and set data for radiation_diag_mod, if needed
!    typically used for diagnostics with the standalone radiation
!---------------------------------------------------------------------

      if (present(Rad_diag)) then
        call set_radiation_diag_type ( Atmos_input%press, &
                       Atmos_input%pflux, Atmos_input%temp, Atmos_input%rh2o, &
                       Surface%asfc_vis_dir, Surface%asfc_nir_dir, &
                       Surface%asfc_vis_dif, Surface%asfc_nir_dif, Astro, &
                       Rad_gases, crndlw, cmxolw, emrndlw, emmxolw, &
                       camtsw, cldsctsw, cldextsw, cldasymmsw, &
                       Sw_output, Lw_output, Lw_diagnostics, &
                       Rad_diag )
      end if

!---------------------------------------------------------------------
!    call deallocate_arrays to deallocate the array space associated 
!    with stack-resident derived-type variables.
!---------------------------------------------------------------------
        call deallocate_arrays (Astro, Astro_phys,    &
                                Sw_output, Lw_output, &
                                Lw_diagnostics, &
                                Aerosol, Aerosolrad_diags)

!--------------------------------------------------------------------

      call mpp_clock_end (misc_clock)

!---------------------------------------------------------------------
!    if COSP is activated and this is a step upon which the cosp
!    simulator is to be called, call return_cosp_inputs to obtain the 
!    radiative inputs that COSP will need.
!---------------------------------------------------------------------
      if (Exch_ctrl%do_cosp .or. Exch_ctrl%do_modis_yim) then
!       if (Cosp_rad_control%step_to_call_cosp) then
        if (Cosp_rad_control%step_to_call_cosp .or. do_rad) then
          call return_cosp_inputs        &
                (is, ie, js, je,  &
                 Time_next, Atmos_input%psfc, Atmos_input%press, Atmos_input%pflux, &
                 Atmos_input%temp, Atmos_input%deltaz, Cosp_rad_block%stoch_cloud_type, &
                 Cosp_rad_block%stoch_conc_drop, Cosp_rad_block%stoch_conc_ice, &
                 Cosp_rad_block%stoch_size_drop, Cosp_rad_block%stoch_size_ice, &
                 Cosp_rad_block%tau_stoch, Cosp_rad_block%lwem_stoch, &
                 Model_microphys, Cldrad_control, Exch_ctrl)

          Cosp_rad_block%mr_ozone(:,:,:) = Rad_gases%qo3(:,:,:)
          where (Rad_output%flux_sw_surf(is:ie,js:je) > 0.0)
             Cosp_rad_block%daytime(:,:) = 1.0
          elsewhere
             Cosp_rad_block%daytime(:,:) = 0.0
          endwhere

        endif ! (step_to_call_cosp)

! 11/15/15 define here rather than in physics_driver_down
         Cosp_rad_block%tsurf_save(:,:) = t_surf_rad(:,:)

        if (Exch_ctrl%do_modis_yim .and. do_rad) then
          call  modis_yim (is, js, Time_next, Atmos_input%psfc,  &
                           Atmos_input%press,  &
                           Cosp_rad_block%tau_stoch(:,:,:,:), &
                           Model_microphys)
        endif
      endif ! (Exch_ctrl%do_cosp)

!-------------------------------------------------------------------
!    process the variables returned from radiation_driver_mod. the 
!    radiative heating rate is added to the accumulated physics heating
!    rate (tdt). net surface lw and sw fluxes and the cosine of the 
!    zenith angle are placed in locations where they can be exported
!    for use in other component models. the lw heating rate is stored
!    in a module variable for potential use in other physics modules.
!    the radiative heating rate is also added to a variable which is
!    accumulating the radiative and turbulent heating rates.
!-------------------------------------------------------------------
      Rad_flux_block%tdt_rad = Rad_output%tdt_rad(is:ie,js:je,:)
      Rad_flux_block%tdt_lw  = Rad_output%tdtlw(is:ie,js:je,:)
      Rad_flux_block%flux_sw = Rad_output%flux_sw_surf(is:ie,js:je)
      if (treat_sfc_refl_dir_as_dif) then
         Rad_flux_block%flux_sw_dir     = Rad_output%flux_sw_surf_dir(is:ie,js:je)
         Rad_flux_block%flux_sw_dif     = Rad_output%flux_sw_surf_dif(is:ie,js:je)
         Rad_flux_block%flux_sw_vis_dir = Rad_output%flux_sw_vis_dir (is:ie,js:je)
         Rad_flux_block%flux_sw_vis_dif = Rad_output%flux_sw_vis_dif (is:ie,js:je)
      else
         Rad_flux_block%flux_sw_dir     = Rad_output%flux_sw_surf_dir(is:ie,js:je) - Rad_output%flux_sw_surf_refl_dir(is:ie,js:je)
         Rad_flux_block%flux_sw_dif     = Rad_output%flux_sw_surf_dif(is:ie,js:je) + Rad_output%flux_sw_surf_refl_dir(is:ie,js:je)
         Rad_flux_block%flux_sw_vis_dir = Rad_output%flux_sw_vis_dir (is:ie,js:je) - Rad_output%flux_sw_refl_vis_dir (is:ie,js:je)
         Rad_flux_block%flux_sw_vis_dif = Rad_output%flux_sw_vis_dif (is:ie,js:je) + Rad_output%flux_sw_refl_vis_dir (is:ie,js:je)
      endif

      Rad_flux_block%flux_sw_down_vis_dir   = Rad_output%flux_sw_down_vis_dir(is:ie,js:je)
      Rad_flux_block%flux_sw_down_vis_dif   = Rad_output%flux_sw_down_vis_dif(is:ie,js:je)
      Rad_flux_block%flux_sw_down_total_dir = Rad_output%flux_sw_down_total_dir(is:ie,js:je)
      Rad_flux_block%flux_sw_down_total_dif = Rad_output%flux_sw_down_total_dif(is:ie,js:je)
      Rad_flux_block%flux_sw_vis            = Rad_output%flux_sw_vis (is:ie,js:je)

      Rad_flux_block%flux_lw = Rad_output%flux_lw_surf(is:ie,js:je)
      Rad_flux_block%coszen  = Rad_output%coszen_angle(is:ie,js:je)
      if (do_rad) Rad_flux_block%extinction = Rad_output%extinction(is:ie,js:je,:)

!---------------------------------------------------------------------
!    call routines to deallocate the components of the derived type 
!    arrays input to radiation_driver.
!---------------------------------------------------------------------
      if (do_rad) then
        call Rad_gases%dealloc
        call Cld_spec%dealloc (Cldrad_control)   !BW, Lsc_microphys)
        call atmos_input_dealloc (Atmos_input)
        call Model_microphys%dealloc (Cldrad_control)
      endif

!BW   call surface_dealloc (Surface)
      Surface%land => null()
      Surface%asfc_vis_dir => null()
      Surface%asfc_nir_dir => null()
      Surface%asfc_vis_dif => null()
      Surface%asfc_nir_dif => null()

!   nullify pointers to input data
      t => null()
      r => null()
      q => null()
      p_full => null()
      p_half => null()
      z_full => null()
      z_half => null()
      gavg_rrv => null()

      call mpp_clock_end ( radiation_clock )

!---------------------------------------------------------------------


end subroutine radiation_driver


!#####################################################################
! <SUBROUTINE NAME="define_rad_times">
!  <OVERVIEW>
!    subroutine define_rad_times determines whether radiation is to be 
!    calculated on the current timestep, and defines logical variables 
!    which determine whether various input fields to radiation_driver 
!    need to be retrieved on the current step.
!  </OVERVIEW>
!  <DESCRIPTION>
!    subroutine define_rad_times determines whether radiation is to be 
!    calculated on the current timestep, and defines logical variables 
!    which determine whether various input fields to radiation_driver 
!    need to be retrieved on the current step.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  define_rad_times (Time, Time_next, Rad_time_out)
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current model time
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   The time used for diagnostic output
!  </IN>
!  <INOUT NAME="Rad_time_out" TYPE="time_type">
!   time at which the climatologically-determined,
!                     time-varying input fields to radiation should 
!                     apply    
!  </INOUT>
! </SUBROUTINE>
!
subroutine define_rad_times (Time, Time_next, Rad_time_out)

!--------------------------------------------------------------------
!    subroutine define_rad_times determines whether radiation is to be 
!    calculated on the current timestep, and defines logical variables 
!    which determine whether various input fields to radiation_driver 
!    need to be retrieved on the current step.
!-------------------------------------------------------------------- 

!---------------------------------------------------------------------
type(time_type), intent(in)     ::  Time, Time_next
type(time_type), intent(inout)  ::  Rad_time_out
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     Time            current model time  
!                     [ time_type, days and seconds]
!     Time_next       model time on the next atmospheric timestep
!                     [ time_type, days and seconds]
!     
!   intent(inout) variables:
!
!     Rad_time_out    time at which the climatologically-determined, 
!                     time-varying input fields to radiation should 
!                     apply    
!                     [ time_type, days and seconds]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer        :: year, month, day, sec
      integer        :: dum, tod(3)
      integer        :: nband

!---------------------------------------------------------------------
!   local variables:
!
!      day            day component of atmospheric timestep
!                     [ days ]
!      sec            seconds component of atmospheric timestep
!                     [ seconds ]
!      dum            dummy variable
!      tod            hours, minutes and seconds components of current
!                     time
!                     [ hours, minutes, seconds ]
!
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    store the atmospheric timestep into a module variable for later
!    use.
!--------------------------------------------------------------------
      call get_time (Time_next-Time, sec, day)    
      dt = day*SECONDS_PER_DAY + sec

!--------------------------------------------------------------------
!    verify that the radiation timestep is an even multiple of the 
!    physics timestep.
!---------------------------------------------------------------------
      if (MOD(lw_rad_time_step, dt) /= 0) then
        call error_mesg ('radiation_driver_mod',  &
    ' lw radiation timestep is not integral multiple of physics step', &
                                                           FATAL)
      endif
      if (MOD(sw_rad_time_step, dt) /= 0) then
        call error_mesg ('radiation_driver_mod',  &
       ' sw radiation timestep is not integral multiple of physics step', &
                                                           FATAL)
      endif

!-------------------------------------------------------------------
!    for the standalone case, new radiation outputs are calculated on 
!    every step, using climatological variable values at the time spec-
!    ified by the input argument Time. 
!-------------------------------------------------------------------
      if (always_calculate) then
        do_rad = .true.
        do_sw_rad = .true.
        do_lw_rad = .true.
        Rad_time = Time
        Rad_control%do_lw_rad = do_lw_rad
        Rad_control%do_sw_rad = do_sw_rad

!--------------------------------------------------------------------
!    determine if this is a radiation time step by decrementing the time
!    to alarm by the current model timestep.  if the alarm "goes off", 
!    i.e., is .le. 0, set do_rad to true, indicating this is a radiation
!    step. otherwise set it to .false. . 
!--------------------------------------------------------------------
      else
        lwrad_alarm = lwrad_alarm -  dt
        swrad_alarm = swrad_alarm -  dt
        if (lwrad_alarm <= 0) then
          do_lw_rad = .true.
        else
          do_lw_rad = .false.
        endif
        if (swrad_alarm <= 0) then
          do_sw_rad = .true.
        else
          do_sw_rad = .false.
        endif
        if (do_sw_rad .or. do_lw_rad) then
           do_rad = .true.
        else
          do_rad = .false.
        endif
      Rad_control%do_lw_rad = do_lw_rad
      Rad_control%do_sw_rad = do_sw_rad

!-------------------------------------------------------------------
!    define the time to be used in defining the time-varying input 
!    fields for the radiation calculation (Rad_time). 
!-------------------------------------------------------------------
        if (rsd) then

!--------------------------------------------------------------------
!    if this is a repeat-same-day (rsd) experiment, define Rad_time
!    as the specified year-month-day (rad_date(1:3)), and the 
!    hr-min-sec of the current time (Time).
!---------------------------------------------------------------------
          if (.not. use_rad_date)   &
            call error_mesg ('radiation_driver_mod', &  
              'if (rsd), must set rad_date(1:3) to valid date', FATAL)
            call get_date (Time, dum, dum, dum, tod(1), tod(2), tod(3))
            Rad_time = set_date (rad_date(1), rad_date(2),& 
                                 rad_date(3), tod(1), tod(2), &
                                 tod(3))

!---------------------------------------------------------------------
!    if the specified date option is active, define Rad_time to be that
!    date and time.
!----------------------------------------------------------------------
        else if (use_rad_date) then
          Rad_time = set_date (rad_date(1), rad_date(2), rad_date(3),  &
                               rad_date(4), rad_date(5), rad_date(6))

!---------------------------------------------------------------------
!    if neither of these special cases is active, define Rad_time as
!    the current time (Time).
!---------------------------------------------------------------------
        else
          Rad_time = Time
        endif  ! (rsd)
      endif  ! (always_calculate)

!---------------------------------------------------------------------
!    place the time at which radiation is to be applied into an output
!    variable.
!---------------------------------------------------------------------
      Rad_time_out = Rad_time

!---------------------------------------------------------------------


end subroutine define_rad_times


!######################################################################
! <SUBROUTINE NAME="define_atmos_input_fields">
!  <OVERVIEW>
!    define_atmos_input_fields converts the atmospheric input fields 
!    (pfull, phalf, t, q, ts) to the form needed by the radiation 
!    modules, and when needed returns radiation-ready fields of pressure
!    (press, psfc), temperature (temp, tsfc), water vapor mixing ratio 
!    (rh2o) and several auxiliary variables in the derived type 
!    structure Atmos_input. the optional input variables are present
!    when running radiative feedback studies (sa_model), and are needed
!    to allow variation of temperature and vapor fields while holding 
!    the aerosol and cloud amounts fixed.
!  </OVERVIEW>
!  <DESCRIPTION>
!    define_atmos_input_fields converts the atmospheric input fields 
!    (pfull, phalf, t, q, ts) to the form needed by the radiation 
!    modules, and when needed returns radiation-ready fields of pressure
!    (press, psfc), temperature (temp, tsfc), water vapor mixing ratio 
!    (rh2o) and several auxiliary variables in the derived type 
!    structure Atmos_input. the optional input variables are present
!    when running radiative feedback studies (sa_model), and are needed
!    to allow variation of temperature and vapor fields while holding 
!    the aerosol and cloud amounts fixed.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  define_atmos_input_fields (is, ie, js, je, pfull, phalf, &
!                                      t, q, ts, r, gavg_rrv, Atmos_input, &
!                                      cloudtemp, cloudvapor, &
!                                      aerosoltemp, aerosolvapor, &
!                                      aerosolpress)  
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!    starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="pfull" TYPE="real">
!   pressure at full levels
!  </IN>
!  <IN NAME="phalf" TYPE="real">
!   pressure at half levels
!  </IN>
!  <IN NAME="t" TYPE="real">
!   temperature at full levels
!  </IN>
!  <IN NAME="q" TYPE="real">
!   specific humidity of water vapor at full levels
!  </IN>
!  <IN NAME="ts" TYPE="real">
!   surface temperature
!  </IN>
!  <IN NAME="r" TYPE="real">
!   tracer array
!  </IN>
!  <IN NAME="gavg_rrv" TYPE="real">
!   global average array of tracer volume mixxing ratio
!  </IN>
!  <INOUT NAME="Atmos_input" TYPE="atmos_input_type">
!   atmos_input type structure, contains the 
!                    following components defined in this subroutine
!  </INOUT>
!  <IN NAME="cloudtemp" TYPE="real">
!    temperature to be seen by clouds (used in 
!                         sa_gcm feedback studies) 
!  </IN>
!  <IN NAME="cloudvapor" TYPE="real">
!   water vapor to be seen by clouds (used in 
!                         sa_gcm feedback studies) 
!  </IN>
!  <IN NAME="aerosoltemp" TYPE="real">
!   required in sa_gcm mode, absent otherwise:
!                         temperature field to be used by aerosol param-
!                         eterization 
!  </IN>
!  <IN NAME="aerosolvapor" TYPE="real">
!   required in sa_gcm mode, absent otherwise:
!                         water vapor field to be used by aerosol param-
!                         eterization 
!  </IN>
!  <IN NAME="aerosolpress" TYPE="real">
!   required in sa_gcm mode, absent otherwise:
!                         pressure field to be used by aerosol param-
!                         eterization
!  </IN>
! </SUBROUTINE>
!
subroutine define_atmos_input_fields (is, ie, js, je, pfull, phalf, &
                                      t, q, ts, r, gavg_rrv, Atmos_input, &
                                      cloudtemp, cloudvapor, &
                                      aerosoltemp, aerosolvapor, &
                                      aerosolpress)     

!---------------------------------------------------------------------
!    define_atmos_input_fields converts the atmospheric input fields 
!    (pfull, phalf, t, q, ts) to the form needed by the radiation 
!    modules, and when needed returns radiation-ready fields of pressure
!    (press, psfc), temperature (temp, tsfc), water vapor mixing ratio 
!    (rh2o) and several auxiliary variables in the derived type 
!    structure Atmos_input. the optional input variables are present
!    when running radiative feedback studies (sa_model), and are needed
!    to allow variation of temperature and vapor fields while holding 
!    the aerosol and cloud amounts fixed.
!---------------------------------------------------------------------
     
integer,                 intent(in)              :: is, ie, js, je
real, dimension(:,:,:),  intent(in), target      :: pfull, phalf, t, q
real, dimension(:,:),    intent(in), target      :: ts
real, dimension(:),      intent(in)              :: gavg_rrv
real, dimension(:,:,:,:),intent(in)              :: r
type(atmos_input_type),  intent(inout)           :: Atmos_input
real, dimension(:,:,:),  intent(in), optional, target :: cloudtemp,    &
                                                         cloudvapor, &
                                                         aerosoltemp, &
                                                         aerosolvapor, &
                                                         aerosolpress

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      pfull        pressure at full levels [ kg / (m s^2) ]
!      phalf        pressure at half levels [ kg / (m s^2) ]
!      t            temperature at full levels [ deg K]
!      q            specific humidity of water vapor at full levels
!                   [ dimensionless ]
!      ts           surface temperature  [ deg K ]
!
!   intent(out) variables:
!
!      Atmos_input   atmos_input type structure, contains the 
!                    following components defined in this subroutine
!         psfc          surface pressure 
!                       [ (kg /( m s^2) ] 
!         tsfc          surface temperature
!                       [ deg K ]
!         temp          temperature at model levels (1:nlev), surface
!                       temperature is stored at value nlev+1; if eta
!                       coordinates, surface value stored in below 
!                       ground points
!                       [ deg K ]
!         press         pressure at model levels (1:nlev), surface 
!                       pressure is stored at index value nlev+1
!                       [ (kg /( m s^2) ] 
!         rh2o          mixing ratio of water vapor at model full levels
!                       [ non-dimensional ]
!         deltaz        model vertical grid separation
!                       [meters]
!         pflux         average of pressure at adjacent model levels
!                       [ (kg /( m s^2) ] 
!         tflux         average of temperature at adjacent model levels
!                       [ deg K ]
!         rel_hum       relative humidity
!                       [ dimensionless ]
!         cloudtemp     temperature to be seen by clouds (used in 
!                       sa_gcm feedback studies) 
!                       [ degrees K ]
!         cloudvapor    water vapor to be seen by clouds (used in 
!                       sa_gcm feedback studies) 
!                       [ nondimensional ]
!         clouddeltaz   deltaz to be used in defining cloud paths (used
!                       in sa_gcm feedback studies)
!                       [ meters ]
!         aerosoltemp   temperature to be seen by aerosols (used in 
!                       sa_gcm feedback studies) 
!                       [ degrees K ]
!         aerosolvapor  water vapor to be seen by aerosols (used in 
!                       sa_gcm feedback studies) 
!                       [ nondimensional ]
!         aerosolpress  pressure field to be seen by aerosols (used in 
!                       sa_gcm feedback studies) 
!                       [ Pa ]
!         aerosolrelhum relative humidity seen by aerosol package,
!                       used in sa_gcm feedback studies
!                       [ dimensionless ]
!
!   intent(in), optional variables:
!
!      cloudtemp          temperature to be seen by clouds (used in 
!                         sa_gcm feedback studies) 
!                         [ degrees K ]
!      cloudvapor         water vapor to be seen by clouds (used in 
!                         sa_gcm feedback studies) 
!                         [ nondimensional ]
!      aerosoltemp        required in sa_gcm mode, absent otherwise:
!                         temperature field to be used by aerosol param-
!                         eterization 
!      aerosolvapor       required in sa_gcm mode, absent otherwise:
!                         water vapor field to be used by aerosol param-
!                         eterization 
!      aerosolpress       required in sa_gcm mode, absent otherwise:
!                         pressure field to be used by aerosol param-
!                         eterization
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables
 
      integer :: i, j, k, kb
      integer :: kmax
      logical :: override
      type(time_type)  :: Data_time
      integer                  :: ico2

!---------------------------------------------------------------------
!  local variables
!
!     i, j, k      do loop indices
!     kb           vertical index of lowest atmospheric level (when
!                  using eta coordinates)
!     kmax         number of model layers
!
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!----------------------------------------------------------------------
!    define the number of model layers.
!----------------------------------------------------------------------
      kmax = size(t,3)

!---------------------------------------------------------------------
!    if the temperature, cloud, or aerosol input data is to be
!    overriden, define the time slice of data which is to be used. allocate
!    storage for the temperature data which will be needed for these
!    cases.
!    THIS OPTION HAS BEEN DIABLED
!---------------------------------------------------------------------
      if (doing_data_override) then
        Data_time = Rad_time                                  
        if (overriding_temps .or. overriding_aerosol .or. &
            overriding_clouds .or. overriding_sphum) then
          call error_mesg ('radiation_driver_mod', &
            'data override for radiation has been disabled in this version', FATAL)
        endif
      endif

!---------------------------------------------------------------------
!    allocate atmos input data arrays
!---------------------------------------------------------------------
      allocate(Atmos_input%temp(size(t,1),size(t,2),kmax))
      allocate(Atmos_input%rh2o(size(t,1),size(t,2),kmax))

      Atmos_input%temp = t
      Atmos_input%rh2o = q

!---------------------------------------------------------------------
!    humidity, convert to mixing ratio
!---------------------------------------------------------------------
      if (.not.use_mixing_ratio) then
         Atmos_input%rh2o = Atmos_input%rh2o/(1.0-Atmos_input%rh2o)
      endif

!---------------------------------------------------------------------
!    surface temp & pressure at model levels and layer interfaces
!---------------------------------------------------------------------
      Atmos_input%tsfc => ts
      Atmos_input%press => pfull
      Atmos_input%phalf => phalf

!---------------------------------------------------------------------
!    define values of surface pressure and temperature.
!--------------------------------------------------------------------
      Atmos_input%psfc => phalf(:,:,kmax+1)
        
!---------------------------------------------------------------------
!    define the cloudtemp component of Atmos_input. 
!---------------------------------------------------------------------
      allocate(Atmos_input%cloudtemp(size(t,1),size(t,2),kmax))
      if (present(cloudtemp)) then
         Atmos_input%cloudtemp = cloudtemp
      else
         Atmos_input%cloudtemp = Atmos_input%temp
      endif

!---------------------------------------------------------------------
!    define the cloudvapor component of Atmos_input.
!---------------------------------------------------------------------
      allocate(Atmos_input%cloudvapor(size(t,1),size(t,2),kmax))
      if (present(cloudvapor)) then
         Atmos_input%cloudvapor = cloudvapor
         if (.not.use_mixing_ratio) then
            Atmos_input%cloudvapor = Atmos_input%cloudvapor/(1.0-Atmos_input%cloudvapor)
         endif
      else
         Atmos_input%cloudvapor = Atmos_input%rh2o
      endif

!---------------------------------------------------------------------
!    define the aerosoltemp component of Atmos_input.
!---------------------------------------------------------------------
      allocate(Atmos_input%aerosoltemp(size(t,1),size(t,2),kmax))
      if (present(aerosoltemp)) then
         Atmos_input%aerosoltemp = aerosoltemp
      else
         Atmos_input%aerosoltemp = Atmos_input%temp
      endif

!---------------------------------------------------------------------
!    define the aerosolvapor component of Atmos_input.
!---------------------------------------------------------------------
      allocate(Atmos_input%aerosolvapor(size(t,1),size(t,2),kmax))
      if (present(aerosolvapor)) then
         Atmos_input%aerosolvapor = aerosolvapor
         if (.not.use_mixing_ratio) then
            Atmos_input%aerosolvapor = Atmos_input%aerosolvapor/(1.0-Atmos_input%aerosolvapor)
         endif
      else
         Atmos_input%aerosolvapor = Atmos_input%rh2o
      endif

!------------------------------------------------------------------
!    define the aerosolpress component of Atmos_input.
!---------------------------------------------------------------------
      if (present(aerosolpress)) then
         Atmos_input%aerosolpress => aerosolpress
      else
         Atmos_input%aerosolpress => Atmos_input%press
      endif
 
!------------------------------------------------------------------
!    be sure that the magnitude of the water vapor mixing ratio field 
!    to be input to the radiation code is no smaller than the value of 
!    rh2o_lower_limit, which is 2.0E-07 when running the sea_esf
!    radiation code and 3.0e-06 when running the original radiation
!    code. Likewise, the temperature that the radiation code sees is
!    constrained to lie between 100K and 370K. these are the limits of
!    the tables referenced within the radiation package.
!      exception:
!    if do_h2o is false, the lower limit of h2o is zero, and radiation
!    tables will not be called.
!-----------------------------------------------------------------------
!BW   if (do_rad) then
      if (apply_vapor_limits) then
        Atmos_input%rh2o(:,:,ks:ke) = MAX(Atmos_input%rh2o(:,:,ks:ke), rh2o_lower_limit)
        Atmos_input%cloudvapor(:,:,ks:ke) = MAX(Atmos_input%cloudvapor(:,:,ks:ke), rh2o_lower_limit)
        Atmos_input%aerosolvapor(:,:,ks:ke) = MAX(Atmos_input%aerosolvapor(:,:,ks:ke), rh2o_lower_limit)
      endif

      if (apply_temp_limits) then
        Atmos_input%temp(:,:,ks:ke) = MIN(MAX(Atmos_input%temp(:,:,ks:ke), temp_lower_limit), temp_upper_limit)
        Atmos_input%cloudtemp(:,:,ks:ke) = MIN(MAX(Atmos_input%cloudtemp(:,:,ks:ke), temp_lower_limit), temp_upper_limit)
        Atmos_input%aerosoltemp(:,:,ks:ke) = MIN(MAX(Atmos_input%aerosoltemp(:,:,ks:ke), temp_lower_limit), temp_upper_limit)
      endif
!BW   endif

!--------------------------------------------------------------------
!    call calculate_aulixiary_variables to compute pressure and 
!    temperature arrays at flux levels and an array of model deltaz.
!--------------------------------------------------------------------
!BW   if (do_rad) then
        call calculate_auxiliary_variables (Atmos_input)
!BW   endif

!RSH
!RSH   define here the values for Atmos_input%tracer_co2.
!RSH
!fil   the error message should never be printed as that code should never
!      be executed, it's an extra guard against user error.
      if (use_co2_tracer_field ) then
         allocate ( Atmos_input%tracer_co2(size(t,1), size(t,2), kmax) )
         ico2 = get_tracer_index(MODEL_ATMOS, 'co2')
         if(ico2 /= NO_TRACER) then
            Atmos_input%tracer_co2(:,:,:) = r(:,:,:,ico2)
            Atmos_input%g_rrvco2 = gavg_rrv(ico2)
         else
            call error_mesg('radiation_driver', &
              'ico2 cannot be NO_TRACER when use_co2_tracer_field is .true.', FATAL)
         endif
      endif


!----------------------------------------------------------------------


end subroutine define_atmos_input_fields 


!#####################################################################
! <SUBROUTINE NAME="define_surface">
!  <OVERVIEW>
!    define_surface stores the input values of land fraction and 
!    surface albedo in a surface_type structure Surface. 
!  </OVERVIEW>
!  <DESCRIPTION>
!    define_surface stores the input values of land fraction and 
!    surface albedo in a surface_type structure Surface. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call define_surface (is, ie, js, je, albedo, land, Surface)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!    starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="albedo" TYPE="real">
!   surface albedo
!  </IN>
!  <IN NAME="land" TYPE="real">
!   fraction of grid box which is land 
!  </IN>
!  <INOUT NAME="Surface" TYPE="surface_type">
!   surface_type structure to be valued
!  </INOUT>
! </SUBROUTINE>
!
subroutine define_surface (is, ie, js, je, albedo, albedo_vis_dir,   &
                           albedo_nir_dir, albedo_vis_dif, &
                           albedo_nir_dif, land, Surface)

!---------------------------------------------------------------------
!    define_surface stores the input values of land fraction and 
!    surface albedo in a surface_type structure Surface.  
!---------------------------------------------------------------------
     
integer,                 intent(in)              :: is, ie, js, je
real, dimension(:,:),    intent(in)              :: albedo, land, &
                                                    albedo_vis_dir,    &
                                                    albedo_nir_dir, &
                                                    albedo_vis_dif,    &
                                                    albedo_nir_dif
type(surface_type),      intent(inout)           :: Surface     

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      albedo       surface albedo  [ dimensionless ]
!      albedo_vis_dir surface visible direct albedo  [ dimensionless ]
!      albedo_nir_dir surface nir direct albedo  [ dimensionless ]
!      albedo_vis_dif surface visible diffuse albedo  [ dimensionless ]
!      albedo_nir_dif surface nir diffuse albedo  [ dimensionless ]
!      land         fraction of grid box which is land [ dimensionless ]
!
!   intent(out) variables:
!
!      Surface       surface_type structure, contains the 
!                    following components defined in this subroutine
!         asfc          surface albedo
!                       [ non-dimensional ]
!         asfc_vis_dir  surface direct visible albedo
!                       [ non-dimensional ]
!         asfc_nir_dir  surface direct nir albedo
!                       [ non-dimensional ]
!         asfc_vis_dif  surface diffuse visible albedo
!                       [ non-dimensional ]
!         asfc_nir_dif  surface diffuse nir albedo
!                       [ non-dimensional ]
!         land          fraction of grid box covered by land
!                       [ non-dimensional ]
!
!---------------------------------------------------------------------

     logical :: override
     type(time_type)  :: Data_time
     real, dimension (size(albedo,1), size(albedo,2)) :: albedo_vis_dir2,   &
                                                         albedo_nir_dir2, &
                                                         albedo_vis_dif2,  &
                                                         albedo_nir_dif2

!-------------------------------------------------------------------
!    verify that the module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

      if (do_rad) then
        if (doing_data_override) then        
!---------------------------------------------------------------------
!    if the albedo data is to be overriden, define the time from which
!    the data is to be retrieved.
!---------------------------------------------------------------------
          if (overriding_albedo) then
            Data_time = Rad_time 


!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's surface
!    albedo data from the override file. if the process fails,
!    write an error message; if it succeeds move the data from the 
!    current window into array albedo2.
!---------------------------------------------------------------------
!           call data_override ('ATM', 'albedonew', albedo_proc,   &
!                             Data_time, override=override)
!           if ( .not. override) then
!             call error_mesg ('radiation_driver_mod', &
!             'cvisrfgd => albedo not overridden successfully', FATAL)
!           else
!             albedo2(:,:) =      albedo_proc(is:ie,js:je)
!           endif
            
            call data_override ('ATM', 'albedo_nir_dir_new',   &
                                albedo_nir_dir2,   &
                                Data_time, override=override, &
                                is_in=is, ie_in=ie, js_in=js, je_in=je)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
                'nirdir => albedo not overridden successfully', FATAL)
            endif

            call data_override ('ATM', 'albedo_nir_dif_new',   &
                                albedo_nir_dif2,   &
                                Data_time, override=override, &
                                is_in=is, ie_in=ie, js_in=js, je_in=je)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
                'nirdif => albedo not overridden successfully', FATAL)
            endif

            call data_override ('ATM', 'albedo_vis_dir_new',   &
                                albedo_vis_dir2,   &
                                Data_time, override=override, &
                                is_in=is, ie_in=ie, js_in=js, je_in=je)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
               'visdir => albedo not overridden successfully', FATAL)
            endif

            call data_override ('ATM', 'albedo_vis_dif_new',   &
                                albedo_vis_dif2,   &
                                Data_time, override=override, &
                                is_in=is, ie_in=ie, js_in=js, je_in=je)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
              'visdif => albedo not overridden successfully', FATAL)
           endif

!--------------------------------------------------------------------
!    if albedo data is not being overriden, define albedo2 to be the 
!    model value of albedo.
!--------------------------------------------------------------------
          else
!           albedo2 = albedo
            albedo_vis_dir2 = albedo_vis_dir
            albedo_nir_dir2 = albedo_nir_dir
            albedo_vis_dif2 = albedo_vis_dif
            albedo_nir_dif2 = albedo_nir_dif
          endif
        else ! (doing data_override)       
!         albedo2 = albedo
          albedo_vis_dir2 = albedo_vis_dir
          albedo_nir_dir2 = albedo_nir_dir
          albedo_vis_dif2 = albedo_vis_dif
          albedo_nir_dif2 = albedo_nir_dif
        endif
      else ! (do_rad)
!       albedo2 = albedo
        albedo_vis_dir2 = albedo_vis_dir
        albedo_nir_dir2 = albedo_nir_dir
        albedo_vis_dif2 = albedo_vis_dif
        albedo_nir_dif2 = albedo_nir_dif
      endif ! (do_rad)

!---------------------------------------------------------------------
!    allocate space for the components of the derived type variable
!    Surface.     
!---------------------------------------------------------------------
      allocate (Surface%asfc (size(albedo,1), size(albedo,2)) )
      allocate (Surface%asfc_vis_dir (size(albedo,1), size(albedo,2) ) )
      allocate (Surface%asfc_nir_dir (size(albedo,1), size(albedo,2) ) )
      allocate (Surface%asfc_vis_dif (size(albedo,1), size(albedo,2) ) )
      allocate (Surface%asfc_nir_dif (size(albedo,1), size(albedo,2) ) )
      allocate (Surface%land (size(albedo,1), size(albedo,2)) )

 
!------------------------------------------------------------------
!    define the fractional land area of each grid box and the surface
!    albedo from the input argument values.
!------------------------------------------------------------------
      Surface%land(:,:) = land(:,:)
      Surface%asfc(:,:) = albedo (:,:)

!pjp  Should the albedos below all be set to albedo2,
!pjp  or should they be included in the override data,
!pjp  or should it not be changed?

      Surface%asfc_vis_dir(:,:) = albedo_vis_dir2(:,:)
      Surface%asfc_nir_dir(:,:) = albedo_nir_dir2(:,:)
      Surface%asfc_vis_dif(:,:) = albedo_vis_dif2(:,:)
      Surface%asfc_nir_dif(:,:) = albedo_nir_dif2(:,:)
     

!----------------------------------------------------------------------


end subroutine define_surface    


!#####################################################################
! <SUBROUTINE NAME="surface_dealloc">
!  <OVERVIEW>
!    surface_dealloc deallocates the array components of the
!    surface_type structure Surface.
!  </OVERVIEW>
!  <DESCRIPTION>
!    surface_dealloc deallocates the array components of the
!    surface_type structure Surface.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call surface_dealloc (Surface)
!  </TEMPLATE>
!  <INOUT NAME="Surface" TYPE="surface_type">
!   surface_type structure to be deallocated
!  </INOUT>
! </SUBROUTINE>
!
subroutine surface_dealloc (Surface)

!----------------------------------------------------------------------
!    surface_dealloc deallocates the array components of the
!    surface_type structure Surface.
!----------------------------------------------------------------------

type(surface_type), intent(inout) :: Surface

!--------------------------------------------------------------------
!   intent(inout) variable:
!
!      Surface        surface_type structure, contains variables 
!                     defining the surface albedo and land fraction
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!-------------------------------------------------------------------
!    deallocate components of surface_type structure.
!-------------------------------------------------------------------
      deallocate (Surface%asfc)
      deallocate (Surface%asfc_vis_dir )
      deallocate (Surface%asfc_nir_dir )
      deallocate (Surface%asfc_vis_dif )
      deallocate (Surface%asfc_nir_dif )
      deallocate (Surface%land)

!--------------------------------------------------------------------


end subroutine surface_dealloc 


!#####################################################################
! <SUBROUTINE NAME="atmos_input_dealloc">
!  <OVERVIEW>
!    atmos_input_dealloc deallocates the array components of the
!    atmos_input_type structure Atmos_input.
!  </OVERVIEW>
!  <DESCRIPTION>
!    atmos_input_dealloc deallocates the array components of the
!    atmos_input_type structure Atmos_input.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call atmos_input_dealloc (Atmos_input)
!  </TEMPLATE>
!  <INOUT NAME="Atmos_input" TYPE="atmos_input_type">
!      atmos_input_type structure, contains variables 
!                     defining the atmospheric pressure, temperature
!                     and moisture distribution.
!  </INOUT>
! </SUBROUTINE>
!
subroutine atmos_input_dealloc (Atmos_input)

!----------------------------------------------------------------------
!    atmos_input_dealloc deallocates the array components of the
!    atmos_input_type structure Atmos_input.
!----------------------------------------------------------------------

type(atmos_input_type), intent(inout) :: Atmos_input

!--------------------------------------------------------------------
!   intent(inout) variable:
!
!      Atmos_input    atmos_input_type structure, contains variables 
!                     defining the atmospheric pressure, temperature
!                     and moisture distribution.
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    deallocate components of atmos_input_type structure.
!---------------------------------------------------------------------

      ! variables always allocated
      deallocate(Atmos_input%temp)
      deallocate(Atmos_input%rh2o)
      deallocate(Atmos_input%cloudtemp)
      deallocate(Atmos_input%cloudvapor)
      deallocate(Atmos_input%aerosoltemp)
      deallocate(Atmos_input%aerosolvapor)

      ! variables always using pointers
      Atmos_input%press => null()
      Atmos_input%phalf => null()
      Atmos_input%psfc  => null()
      Atmos_input%tsfc  => null()
      Atmos_input%aerosolpress => null()

      ! flux pressure
      if (do_conserve_energy) then
        Atmos_input%pflux => null()
      else
        deallocate (Atmos_input%pflux)
      endif

      ! auxiliary variables always allocated
      deallocate (Atmos_input%relhum     )
      deallocate (Atmos_input%tflux      )
      deallocate (Atmos_input%deltaz     )
      deallocate (Atmos_input%clouddeltaz)
      deallocate (Atmos_input%aerosolrelhum)

      if (ASSOCIATED(Atmos_input%tracer_co2)) deallocate(Atmos_input%tracer_co2)
!--------------------------------------------------------------------


end subroutine atmos_input_dealloc 


!#####################################################################
! <SUBROUTINE NAME="radiation_driver_end">
!  <OVERVIEW>
!   radiation_driver_end is the destructor for radiation_driver_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   radiation_driver_end is the destructor for radiation_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiation_driver_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine radiation_driver_end (Rad_flux, Atm_block)

type(radiation_flux_type), intent(in) :: Rad_flux
type(block_control_type),  intent(in) :: Atm_block

!----------------------------------------------------------------------
!    radiation_driver_end is the destructor for radiation_driver_mod.
!----------------------------------------------------------------------

integer :: cloud_spec_term_clock, aerosol_term_clock, &
           radiative_gases_term_clock, radiation_term_clock
integer :: outunit

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!---------------------------------------------------------------------

      cloud_spec_term_clock       =       &
        mpp_clock_id( '   Phys_driver_term: Cloud spec: Termination', &
                       grain=CLOCK_MODULE_DRIVER )
      aerosol_term_clock       =       &
        mpp_clock_id( '   Phys_driver_term: Aerosol: Termination', &
                       grain=CLOCK_MODULE_DRIVER )
      radiative_gases_term_clock       =       &
        mpp_clock_id( '   Phys_driver_term: Radiative gases: Termination', &
                       grain=CLOCK_MODULE_DRIVER )
      radiation_term_clock       =       &
        mpp_clock_id( '   Phys_driver_term: Radiation: Termination', &
                       grain=CLOCK_MODULE_DRIVER )

      call mpp_clock_begin ( radiation_term_clock )

!---------------------------------------------------------------------
!    save restarts if necessary
!---------------------------------------------------------------------
      call radiation_driver_restart_nc (Rad_flux, Atm_block)

!---------------------------------------------------------------------
!    wrap up modules initialized by this module.
!---------------------------------------------------------------------
      call astronomy_end

!---------------------------------------------------------------------
!    wrap up modules specific to the radiation package in use.
!---------------------------------------------------------------------

      call radiation_driver_diag_end (Rad_control)
      call rad_output_file_end
      call longwave_driver_end
      call shortwave_driver_end

      call mpp_clock_begin ( cloud_spec_term_clock )
      call cloudrad_driver_end (Cldrad_control)
      call mpp_clock_end ( cloud_spec_term_clock )

      call mpp_clock_begin ( radiative_gases_term_clock )
      call radiative_gases_end
      call mpp_clock_end ( radiative_gases_term_clock )

      call mpp_clock_begin ( aerosol_term_clock )
      call aerosolrad_driver_end (Aerosolrad_control)
      call mpp_clock_end ( aerosol_term_clock )

!---------------------------------------------------------------------
!    release space used for module variables that hold data between
!    timesteps.
!---------------------------------------------------------------------
      call Rad_output%dealloc

!----------------------------------------------------------------------------
!    print out checksum info for Rad_flux when concurrent radiation is active
!----------------------------------------------------------------------------
      if (do_concurrent_radiation) then
        outunit = stdout()
        write(outunit,'(a,i1,a)') 'ENDING CHECKSUM(radiation_driver_restart - concurrent):: '
        write(outunit,100) 'tdt_rad                ', mpp_chksum(Restart%tdt_rad                )
        write(outunit,100) 'tdt_lw                 ', mpp_chksum(Restart%tdt_lw                 )
        write(outunit,100) 'flux_sw                ', mpp_chksum(Restart%flux_sw                )
        write(outunit,100) 'flux_sw_dir            ', mpp_chksum(Restart%flux_sw_dir            )
        write(outunit,100) 'flux_sw_dif            ', mpp_chksum(Restart%flux_sw_dif            )
        write(outunit,100) 'flux_sw_down_vis_dir   ', mpp_chksum(Restart%flux_sw_down_vis_dir   )
        write(outunit,100) 'flux_sw_down_vis_dif   ', mpp_chksum(Restart%flux_sw_down_vis_dif   )
        write(outunit,100) 'flux_sw_down_total_dir ', mpp_chksum(Restart%flux_sw_down_total_dir )
        write(outunit,100) 'flux_sw_down_total_dif ', mpp_chksum(Restart%flux_sw_down_total_dif )
        write(outunit,100) 'flux_sw_vis            ', mpp_chksum(Restart%flux_sw_vis            )
        write(outunit,100) 'flux_sw_vis_dir        ', mpp_chksum(Restart%flux_sw_vis_dir        )
        write(outunit,100) 'flux_sw_vis_dif        ', mpp_chksum(Restart%flux_sw_vis_dif        )
        write(outunit,100) 'flux_lw                ', mpp_chksum(Restart%flux_lw                )
        write(outunit,100) 'coszen                 ', mpp_chksum(Restart%coszen                 )
        write(outunit,100) 'extinction             ', mpp_chksum(Restart%extinction             )
100 FORMAT("CHECKSUM::",A32," = ",Z20)
      endif

!----------------------------------------------------------------------
!    set initialization status flag.
!----------------------------------------------------------------------
      module_is_initialized = .false.

!----------------------------------------------------------------------

end subroutine radiation_driver_end

!######################################################################

subroutine return_cosp_inputs (  &
                      is, ie, js, je, Time_diag, &
                      psfc, press, pflux, temp, deltaz, stoch_cloud_type, &
                      stoch_conc_drop, stoch_conc_ice, stoch_size_drop,&
                      stoch_size_ice, tau_stoch, lwem_stoch, &
                      Model_microphys, Cldrad_control, Exch_ctrl)


!---------------------------------------------------------------------
!    subroutine return_cosp_inputs calculates and returns the fields 
!    needed as input by the COSP simulator.
!---------------------------------------------------------------------

integer, intent(in)                     :: is,ie, js, je
type(time_type),        intent(in)      :: Time_diag
real, dimension(:,:),   intent(in)      :: psfc
real, dimension(:,:,:), intent(in)      :: press, pflux, temp, deltaz
type(microphysics_type), intent(inout)  :: Model_microphys
type(cloudrad_control_type), intent(in) :: Cldrad_control
type(exchange_control_type), intent(in) :: Exch_ctrl
real, dimension(:,:,:,:), intent(out)   ::    &
                                stoch_cloud_type, stoch_conc_drop, &
                                stoch_conc_ice, stoch_size_drop,  &
                                stoch_size_ice, tau_stoch, lwem_stoch

!-------------------------------------------------------------------
!   local variables

    integer :: ncld, ic1, ic2, iclast
    integer :: strat_index, shallow_index, donner_meso_index, donner_cell_index

!-------------------------------------------------------------------

!-------------------------------------------------------------------
      call obtain_cloud_tau_and_em (is, js, deltaz, Cldrad_control, &
                                    Model_microphys, &
                                    tau_stoch(:,:,:,:),  &
                                    lwem_stoch(:,:,:,:) )

!-------------------------------------------------------------------
      if (Exch_ctrl%do_cosp) then

!-------------------------------------------------------------------
!    define indexing of cloud schemes
!    Model_microphys%scheme_name contains a comma-separated list of all schemes
!-------------------------------------------------------------------
        strat_index = 0
        shallow_index = 0
        donner_meso_index = 0
        donner_cell_index = 0

        ncld = 0
        ic1 = 1
        iclast = len_trim(Model_microphys%scheme_name)
        do while (ic1 .le. iclast)
          ic2 = index( Model_microphys%scheme_name(ic1:iclast), ',' ) ! parse string
          if ( ic2 > 0 ) then
            ic2 = ic1 + ic2 - 2
          else
            ic2 = iclast
          endif
          ncld = ncld + 1
          if (Model_microphys%scheme_name(ic1:ic2) == 'strat_cloud') strat_index = ncld
          if (Model_microphys%scheme_name(ic1:ic2) == 'uw_conv')     shallow_index = ncld
          if (Model_microphys%scheme_name(ic1:ic2) == 'donner_meso') donner_meso_index = ncld
          if (Model_microphys%scheme_name(ic1:ic2) == 'donner_cell') donner_cell_index = ncld
          ic1 = ic2 + 2
        enddo

!---------------------------------------------------------------------
!    save the stochastic cloud type in each subcolumn.
!    output values of 0 --> no cloud
!           values of 1 --> stratiform cloud
!           values of 2 --> convective cloud
!    input values are 0(none), 1(strat), 2(donnermeso), 3(donnercell), 
!    4(uw)
!---------------------------------------------------------------------
        stoch_cloud_type(:,:,:,:) =   &
                         Model_microphys%stoch_cloud_type(:,:,:,:)
         
!---------------------------------------------------------------------
!    donner meso clouds may be treated either as large-scale or
!    convective clouds, dependent on donner_meso_is_largescale.
!---------------------------------------------------------------------
        if (Exch_ctrl%donner_meso_is_largescale) then
          where (stoch_cloud_type(:,:,:,:) == donner_meso_index)  ! == 2)
            stoch_cloud_type(:,:,:,:) = strat_index               ! = 1
          end where
          where (stoch_cloud_type(:,:,:,:) == donner_cell_index .or. &
                 stoch_cloud_type(:,:,:,:) == shallow_index)      ! >= 3)
            stoch_cloud_type(:,:,:,:) = donner_meso_index         ! = 2)
          end where
        else
          where (stoch_cloud_type(:,:,:,:) == donner_meso_index .or. &
                 stoch_cloud_type(:,:,:,:) == donner_cell_index .or. &
                 stoch_cloud_type(:,:,:,:) == shallow_index)      ! >= 2)
            stoch_cloud_type(:,:,:,:) = donner_meso_index         ! = 2)
          end where
        endif    

!---------------------------------------------------------------------
!    save the particle concentrations and sizes seen by the radiation
!    package in each stochastic column.
!---------------------------------------------------------------------
        stoch_conc_drop(:,:,:,:) =  &
                             Model_microphys%stoch_conc_drop(:,:,:,:)
        stoch_conc_ice (:,:,:,:) =  &
                             Model_microphys%stoch_conc_ice (:,:,:,:)
        stoch_size_drop(:,:,:,:) =  &
                             Model_microphys%stoch_size_drop(:,:,:,:)
        stoch_size_ice (:,:,:,:) =  &
                             Model_microphys%stoch_size_ice (:,:,:,:)

      endif

!-------------------------------------------------------------------
!     if (Exch_ctrl%do_modis_yim) then
!       call  modis_yim (is, js, Time_diag, psfc, press, &
!                        Tau_stoch(:,:,:,:), Model_microphys)
!     endif
      call modis_cmip (is, js, Time_diag, press, pflux, temp, &
                       Model_microphys)

!-------------------------------------------------------------------

end subroutine return_cosp_inputs 

!#######################################################################
! <SUBROUTINE NAME="radiation_driver_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine radiation_driver_restart(Rad_flux, Atm_block, timestamp)
  type(radiation_flux_type), intent(in)  :: Rad_flux
  type(block_control_type),  intent(in)  :: Atm_block
  character(len=*), intent(in), optional :: timestamp

      call radiation_driver_restart_nc (Rad_flux, Atm_block, timestamp)
      call write_solar_interp_restart_nc (timestamp)
      call radiative_gases_restart (timestamp)

 end subroutine radiation_driver_restart

!#######################################################################
! <SUBROUTINE NAME="radiation_driver_restart_nc">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine radiation_driver_restart_nc(Rad_flux, Atm_block, timestamp)
  type(radiation_flux_type), intent(in)  :: Rad_flux
  type(block_control_type),  intent(in)  :: Atm_block
  character(len=*), intent(in), optional :: timestamp

!---local variables
  integer :: n, ibs, ibe, jbs, jbe

!---------------------------------------------------------------------
!    write restart file if desired; the file is not necessary if job 
!    ends on step prior to radiation ts, or if restart seamlessness 
!    is not required.
!---------------------------------------------------------------------
  if (do_concurrent_radiation) then
    do n = 1,size(Rad_flux%block,1)
      ibs = Atm_block%ibs(n)-Atm_block%isc+1
      ibe = Atm_block%ibe(n)-Atm_block%isc+1
      jbs = Atm_block%jbs(n)-Atm_block%jsc+1
      jbe = Atm_block%jbe(n)-Atm_block%jsc+1

      Restart%tdt_rad(ibs:ibe,jbs:jbe,:)              = Rad_flux%block(n)%tdt_rad
      Restart%tdt_lw(ibs:ibe,jbs:jbe,:)               = Rad_flux%block(n)%tdt_lw
      Restart%flux_sw(ibs:ibe,jbs:jbe)                = Rad_flux%block(n)%flux_sw
      Restart%flux_sw_dir(ibs:ibe,jbs:jbe)            = Rad_flux%block(n)%flux_sw_dir
      Restart%flux_sw_dif(ibs:ibe,jbs:jbe)            = Rad_flux%block(n)%flux_sw_dif
      Restart%flux_sw_down_vis_dir(ibs:ibe,jbs:jbe)   = Rad_flux%block(n)%flux_sw_down_vis_dir
      Restart%flux_sw_down_vis_dif(ibs:ibe,jbs:jbe)   = Rad_flux%block(n)%flux_sw_down_vis_dif
      Restart%flux_sw_down_total_dir(ibs:ibe,jbs:jbe) = Rad_flux%block(n)%flux_sw_down_total_dir
      Restart%flux_sw_down_total_dif(ibs:ibe,jbs:jbe) = Rad_flux%block(n)%flux_sw_down_total_dif
      Restart%flux_sw_vis(ibs:ibe,jbs:jbe)            = Rad_flux%block(n)%flux_sw_vis
      Restart%flux_sw_vis_dir(ibs:ibe,jbs:jbe)        = Rad_flux%block(n)%flux_sw_vis_dir
      Restart%flux_sw_vis_dif(ibs:ibe,jbs:jbe)        = Rad_flux%block(n)%flux_sw_vis_dif
      Restart%flux_lw(ibs:ibe,jbs:jbe)                = Rad_flux%block(n)%flux_lw
      Restart%coszen(ibs:ibe,jbs:jbe)                 = Rad_flux%block(n)%coszen
      Restart%extinction(ibs:ibe,jbs:jbe,:)           = Rad_flux%block(n)%extinction
    end do

    if (mpp_pe() == mpp_root_pe() ) then 
       call error_mesg('radiation_driver_mod', 'Writing netCDF formatted restart file: '//&
                       'RESTART/conc_radiation_driver.res.nc', NOTE)
    endif

    call save_restart(Rad_restart_conc,timestamp)
    call save_restart(Til_restart_conc,timestamp)
  endif

  if (using_restart_file) then
!---------------------------------------------------------------------
! Make sure that the restart_versions variable is up to date.
!---------------------------------------------------------------------
    vers = restart_versions(size(restart_versions(:)))
    call write_restart_nc(timestamp)
  endif

end subroutine radiation_driver_restart_nc
! </SUBROUTINE> NAME="radiation_driver_restart_nc"

!#######################################################################

subroutine radiation_diag_type_dealloc ( Rad_diag )
type(radiation_diag_type), intent(inout) :: Rad_diag

integer :: n

    Rad_diag%do_radiation_diag = .false.

    deallocate(Rad_diag%emrndlw)
    deallocate(Rad_diag%emmxolw)
    deallocate(Rad_diag%cldsct)
    deallocate(Rad_diag%cldext)
    deallocate(Rad_diag%cldasymm)
    deallocate(Rad_diag%camtsw)
    deallocate(Rad_diag%crndlw)
    deallocate(Rad_diag%cmxolw)

    deallocate(Rad_diag%press)
    deallocate(Rad_diag%pflux)
    deallocate(Rad_diag%temp)
    deallocate(Rad_diag%rh2o)
    deallocate(Rad_diag%qo3)

    deallocate(Rad_diag%asfc_vis_dir)
    deallocate(Rad_diag%asfc_nir_dir)
    deallocate(Rad_diag%asfc_vis_dif)
    deallocate(Rad_diag%asfc_nir_dif)

    deallocate(Rad_diag%cosz)
    deallocate(Rad_diag%fracday)

    nullify(Rad_diag%Lw_tables)

    do n = 1, size(Rad_diag%Sw_output)
      call Rad_diag%Sw_output(n)%dealloc
    end do
    do n = 1, size(Rad_diag%Lw_output)
      call Rad_diag%Lw_output(n)%dealloc
    end do
    call Rad_diag%Lw_diagnostics%dealloc

!-------------------------------------------------------------------

end subroutine radiation_diag_type_dealloc

!-------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!---------------------------------------------------------------------

subroutine write_restart_nc(timestamp)
  character(len=*), intent(in), optional :: timestamp


  if( .not. Rad_control%using_restart_file ) return
!---------------------------------------------------------------------
!    only the root pe will write control information -- the last value 
!    in the list of restart versions and the alarm information.
!---------------------------------------------------------------------
        if (mpp_pe() == mpp_root_pe() ) then
         call error_mesg('radiation_driver_mod', 'Writing netCDF formatted restart file: RESTART/radiation_driver.res.nc', NOTE)
        endif

!---------------------------------------------------------------------
!    write out the optional time average restart data. note that 
!    do_average and renormalize_sw_fluxes may not both be true.
!---------------------------------------------------------------------
        int_renormalize_sw_fluxes = 0
        int_do_clear_sky_pass = 0
        if(renormalize_sw_fluxes) int_renormalize_sw_fluxes = 1
        if(do_clear_sky_pass) int_do_clear_sky_pass = 1

! Make sure that the restart_versions variable is up to date.
        vers = restart_versions(size(restart_versions(:)))
        call save_restart(Rad_restart, timestamp)
        if(in_different_file) call save_restart(Til_restart, timestamp)
!---------------------------------------------------------------------

end subroutine write_restart_nc

!#####################################################################

subroutine conc_rad_register_restart(fname, Rad_flux, Exch_ctrl, Atm_block)
  character(len=*),            intent(in) :: fname
  type(radiation_flux_type),   intent(in) :: Rad_flux
  type(exchange_control_type), intent(in) :: Exch_ctrl
  type(block_control_type),    intent(in) :: Atm_block

!---local variables
  character(len=64)                     :: fname2
  integer :: ix, jx, npz, id_restart

   call get_mosaic_tile_file(fname, fname2, .false. )
   allocate(Rad_restart_conc)
   if(trim(fname2) == trim(fname)) then
      Til_restart_conc => Rad_restart_conc
      in_different_file_conc = .false.
   else
      in_different_file_conc = .true.
      allocate(Til_restart_conc)
   endif

   if (do_concurrent_radiation) ido_conc_rad = 1
   if (Exch_ctrl%donner_meso_is_largescale) idonner_meso = 1
   if (Exch_ctrl%doing_donner) idoing_donner = 1
   if (Exch_ctrl%doing_uw_conv) idoing_uw_conv = 1

   id_restart = register_restart_field(Rad_restart_conc, fname, 'do_concurrent_radiation', ido_conc_rad, no_domain=.true.)
   id_restart = register_restart_field(Rad_restart_conc, fname, 'donner_meso_is_largescale', idonner_meso, no_domain=.true.)
   id_restart = register_restart_field(Rad_restart_conc, fname, 'doing_donner', idoing_donner, no_domain=.true.)
   id_restart = register_restart_field(Rad_restart_conc, fname, 'doing_uw_conv', idoing_uw_conv, no_domain=.true.)

   ix = Atm_block%iec-Atm_block%isc+1
   jx = Atm_block%jec-Atm_block%jsc+1
   npz = Atm_block%npz
   call Restart%alloc (ix,jx,npz)

   id_restart = register_restart_field(Til_restart_conc, fname, 'tdt_rad',                Restart%tdt_rad)
   id_restart = register_restart_field(Til_restart_conc, fname, 'tdt_lw',                 Restart%tdt_lw)
   id_restart = register_restart_field(Til_restart_conc, fname, 'flux_sw',                Restart%flux_sw)
   id_restart = register_restart_field(Til_restart_conc, fname, 'flux_sw_dir',            Restart%flux_sw_dir)
   id_restart = register_restart_field(Til_restart_conc, fname, 'flux_sw_dif',            Restart%flux_sw_dif)
   id_restart = register_restart_field(Til_restart_conc, fname, 'flux_sw_down_vis_dir',   Restart%flux_sw_down_vis_dir)
   id_restart = register_restart_field(Til_restart_conc, fname, 'flux_sw_down_vis_dif',   Restart%flux_sw_down_vis_dif)
   id_restart = register_restart_field(Til_restart_conc, fname, 'flux_sw_down_total_dir', Restart%flux_sw_down_total_dir)
   id_restart = register_restart_field(Til_restart_conc, fname, 'flux_sw_down_total_dif', Restart%flux_sw_down_total_dif)
   id_restart = register_restart_field(Til_restart_conc, fname, 'flux_sw_vis',            Restart%flux_sw_vis)
   id_restart = register_restart_field(Til_restart_conc, fname, 'flux_sw_vis_dir',        Restart%flux_sw_vis_dir)
   id_restart = register_restart_field(Til_restart_conc, fname, 'flux_sw_vis_dif',        Restart%flux_sw_vis_dif)
   id_restart = register_restart_field(Til_restart_conc, fname, 'flux_lw',                Restart%flux_lw)
   id_restart = register_restart_field(Til_restart_conc, fname, 'coszen',                 Restart%coszen)
   id_restart = register_restart_field(Til_restart_conc, fname, 'extinction',             Restart%extinction)

end subroutine conc_rad_register_restart

!#####################################################################

subroutine rad_driver_register_restart(fname)
  character(len=*), intent(in) :: fname
!---------------------------------------------------------------------
  character(len=64)            :: fname2
  integer                      :: id_restart
!---------------------------------------------------------------------

   call get_mosaic_tile_file(fname, fname2, .false. ) 
   allocate(Rad_restart)
   if(trim(fname2) == trim(fname)) then
      Til_restart => Rad_restart
      in_different_file = .false.
   else
      in_different_file = .true.
      allocate(Til_restart)
   endif

  id_restart = register_restart_field(Rad_restart, fname, 'vers', vers, no_domain=.true.)
  id_restart = register_restart_field(Rad_restart, fname, 'lwrad_alarm', lwrad_alarm, mandatory=.false.,no_domain=.true.)
  id_restart = register_restart_field(Rad_restart, fname, 'swrad_alarm', swrad_alarm, mandatory=.false.,no_domain=.true.)
  id_restart = register_restart_field(Rad_restart, fname, 'lw_rad_time_step', lw_rad_time_step, mandatory=.false.,no_domain=.true.)
  id_restart = register_restart_field(Rad_restart, fname, 'sw_rad_time_step', sw_rad_time_step, mandatory=.false.,no_domain=.true.)
  id_restart = register_restart_field(Rad_restart, fname, 'renormalize_sw_fluxes', int_renormalize_sw_fluxes,no_domain=.true.)
  id_restart = register_restart_field(Rad_restart, fname, 'do_clear_sky_pass', int_do_clear_sky_pass,no_domain=.true.)
  id_restart = register_restart_field(Til_restart, fname, 'tdt_rad', Rad_output%tdt_rad)
  id_restart = register_restart_field(Til_restart, fname, 'tdtlw', Rad_output%tdtlw)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_surf', Rad_output%flux_sw_surf)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_surf_dir', Rad_output%flux_sw_surf_dir)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_surf_refl_dir', Rad_output%flux_sw_surf_refl_dir)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_surf_dif', Rad_output%flux_sw_surf_dif)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_down_vis_dir', Rad_output%flux_sw_down_vis_dir)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_down_vis_dif', Rad_output%flux_sw_down_vis_dif)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_down_total_dir', Rad_output%flux_sw_down_total_dir)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_down_total_dif', Rad_output%flux_sw_down_total_dif)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_vis', Rad_output%flux_sw_vis)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_vis_dir', Rad_output%flux_sw_vis_dir)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_refl_vis_dir', Rad_output%flux_sw_refl_vis_dir)
  id_restart = register_restart_field(Til_restart, fname, 'flux_sw_vis_dif', Rad_output%flux_sw_vis_dif)
  id_restart = register_restart_field(Til_restart, fname, 'flux_lw_surf', Rad_output%flux_lw_surf)
  id_restart = register_restart_field(Til_restart, fname, 'coszen_angle', Rad_output%coszen_angle)

!---------------------------------------------------------------------

end subroutine rad_driver_register_restart

!######################################################################
! <SUBROUTINE NAME="obtain_astronomy_variables">
!  <OVERVIEW>
!    obtain_astronomy_variables retrieves astronomical variables, valid 
!    at the requested time and over the requested time intervals.
!  </OVERVIEW>
!  <DESCRIPTION>
!    obtain_astronomy_variables retrieves astronomical variables, valid 
!    at the requested time and over the requested time intervals.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_astronomy_variables (is, ie, js, je, lat, lon,     &
!                                       Astro, Astro_phys)  
!
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending i,j indices in global storage arrays
!  </IN> 
!  <INOUT NAME="Astro" TYPE="astronomy_type">
!     astronomy_type structure; It will
!                    be used to determine the insolation at toa seen
!                    by the shortwave radiation code
!  </INOUT>
!  <INOUT NAME="Astro_phys" TYPE="astronomy_type">
!     astronomy_type structure, defined when renormal-
!                    ization is active. the same components are defined
!                    as for Astro, but they are valid over the current
!                    physics timestep.
!  </INOUT>
!  <INOUT NAME="Astronomy_inp" TYPE="astronomy_inp_type">
!     astronomy_inp_type structure, optionally used to input astronom-
!     ical forcings, when it is desired to specify them rather than use
!     astronomy_mod. Used in various standalone applications.
!  </INOUT>
!  <IN NAME="lon" TYPE="real">
!    lon        mean longitude (in radians) of all grid boxes processed by
!               this call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <IN NAME="lat" TYPE="real">
!    lat        mean latitude (in radians) of all grid boxes processed by this
!               call to radiation_driver   [real, dimension(:,:)]
!  </IN>
! </SUBROUTINE>
!
subroutine obtain_astronomy_variables (is, ie, js, je, lat, lon,     &
                                       Astro, Astro_phys, Astronomy_inp)  

!---------------------------------------------------------------------
!    obtain_astronomy_variables retrieves astronomical variables, valid 
!    at the requested time and over the requested time intervals.
!---------------------------------------------------------------------
integer,                     intent(in)    ::  is, ie, js, je
real, dimension(:,:),        intent(in)    ::  lat, lon
type(astronomy_type),        intent(inout) ::  Astro, Astro_phys
type(astronomy_inp_type),   intent(inout), optional ::  &
                                               Astronomy_inp

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      lat          latitude of model points  
!                   [ radians ]
!      lon          longitude of model points 
!                   [ radians ]
!
!   intent(inout) variables:
!
!      Astro         astronomy_type structure; contains the following
!                    components defined in this subroutine that will
!                    be used to determine the insolation at toa seen
!                    by the shortwave radiation code 
!         solar         shortwave flux factor: cosine of zenith angle *
!                       daylight fraction / (earth-sun distance squared)
!                       [ non-dimensional ]
!         cosz          cosine of zenith angle --  mean value over
!                       appropriate averaging interval
!                       [ non-dimensional ]
!         fracday       fraction of timestep during which the sun is 
!                       shining
!                       [ non-dimensional ]
!         rrsun         inverse of square of earth-sun distance, 
!                       relative to the mean square of earth-sun 
!                       distance
!                       [ non-dimensional ]
!
!      Astro_phys    astronomy_type structure, defined when renormal-
!                    ization is active. the same components are defined
!                    as for Astro, but they are valid over the current
!                    physics timestep.
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      type(time_type)                   :: Dt_zen, Dt_zen2
      type(time_type)                   :: Rad1    
      real, dimension(ie-is+1, je-js+1) ::                            &
                                           cosz_r, solar_r, fracday_r, &
                                           cosz_p, solar_p, fracday_p, &
                                           cosz_a, fracday_a
      real                              :: rrsun_r, rrsun_p, rrsun_a
     !integer                           :: nz
      

!--------------------------------------------------------------------
!  local variables:
!
!     Dt_zen        time-type variable containing the components of the
!                   radiation time step, needed unless do_average is
!                   true or this is not a radiation step and renormal-
!                   ize_sw_fluxes is true
!     Dt_zen2       time-type variable containing the components of the
!                   physics time step, needed when renormalize_sw_fluxes
!                   or do_average is true
!     cosz_r        cosine of zenith angle --  mean value over
!                   radiation time step            
!                   [ non-dimensional ]
!     solar_r       shortwave flux factor relevant over radiation time
!                   step: cosine of zenith angle * daylight fraction / 
!                   (earth-sun distance squared)
!                   [ non-dimensional ]
!     fracday_r     fraction of timestep during which the sun is 
!                   shining over radiation time step
!                   [ non-dimensional ]
!     cosz_p        cosine of zenith angle --  mean value over
!                   physics time step            
!                   [ non-dimensional ]
!     solar_p       shortwave flux factor relevant over physics time
!                   step: cosine of zenith angle * daylight fraction / 
!                   (earth-sun distance squared)
!                   [ non-dimensional ]
!     fracday_p     fraction of timestep during which the sun is 
!                   shining over physics time step
!                   [ non-dimensional ]
!     cosz_a        cosine of zenith angle --  mean value over
!                   next radiation time step            
!                   [ non-dimensional ]
!     solar_a       shortwave flux factor relevant over next radiation 
!                   time step: cosine of zenith angle * daylight 
!                   fraction / (earth-sun distance squared)
!                   [ non-dimensional ]
!     fracday_a     fraction of timestep during which the sun is 
!                   shining over next radiation time step
!                   [ non-dimensional ]
!     rrsun_r       inverse of square of earth-sun distance, 
!                   relative to the mean square of earth-sun 
!                   distance, valid over radiation time step
!                   [ non-dimensional ]
!     rrsun_p       inverse of square of earth-sun distance, 
!                   relative to the mean square of earth-sun 
!                   distance, valid over physics time step
!                   [ non-dimensional ]
!     rrsun_a       inverse of square of earth-sun distance, 
!                   relative to the mean square of earth-sun 
!                   distance, valid over next radiation time step
!                   [ non-dimensional ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    allocate the components of the astronomy_type structure which will
!    return the astronomical inputs to radiation (cosine of zenith 
!    angle, daylight fraction, solar flux factor and earth-sun distance)
!    that are to be used on the current step.
!---------------------------------------------------------------------
      call Astro%alloc (size(lat,1), size(lat,2))

!---------------------------------------------------------------------
!    case 0: input parameters.
!---------------------------------------------------------------------
      if (present (Astronomy_inp)) then
        Astro%rrsun        = Astronomy_inp%rrsun
        Astro%fracday(:,:) = Astronomy_inp%fracday(is:ie,js:je)
        Astro%cosz (:,:) = cos(Astronomy_inp%zenith_angle(is:ie,js:je)/RADIAN)
        Astro%solar(:,:) = Astro%cosz(:,:)*Astro%fracday(:,:)*Astro%rrsun
        Rad_output%coszen_angle(is:ie,js:je) = Astro%cosz(:,:)

!---------------------------------------------------------------------
!    case 1: diurnally-varying shortwave radiation.
!---------------------------------------------------------------------
      else if (Rad_control%do_diurnal) then

!-------------------------------------------------------------------
!    convert the radiation timestep and the model physics timestep
!    to time_type variables.
!-------------------------------------------------------------------
        Dt_zen  = set_time (sw_rad_time_step, 0)
        Dt_zen2 = set_time (dt, 0)
        
!---------------------------------------------------------------------
!    calculate the astronomical factors averaged over the radiation time
!    step between Rad_time and Rad_time + Dt_zen. these values are 
!    needed on radiation steps. output is stored in Astro_rad.
!---------------------------------------------------------------------
        if (do_sw_rad) then
!  calculation for full radiation step:
          call diurnal_solar (lat, lon, Rad_time, cosz_r, fracday_r, &
                              rrsun_r, dt_time=Dt_zen)
          fracday_r = MIN (fracday_r, 1.00)
          solar_r = cosz_r*fracday_r*rrsun_r
        endif

!---------------------------------------------------------------------
!    calculate the astronomical factors averaged over the physics time
!    step between Rad_time and Rad_time + Dt_zen2. these values are
!    needed if either renormalization or time-averaging is active. store
!    the astronomical outputs in Astro_phys.
!---------------------------------------------------------------------
        if (renormalize_sw_fluxes) then
          call diurnal_solar (lat, lon, Rad_time, cosz_p, fracday_p, &
                              rrsun_p, dt_time=Dt_zen2)
          fracday_p = MIN (fracday_p, 1.00)
          solar_p = cosz_p*fracday_p*rrsun_p
        endif

!--------------------------------------------------------------------
!    define the astronomy_type variable(s) to be returned and used in 
!    the radiation calculation. Astro contains the values to be used
!    in the radiation calculation, Astro_phys contains values relevant 
!    over the current physics timestep and is used for renormalization.
!    when renormalization is active, the physics step set is always 
!    needed, and in addition on radiation steps, the radiation step
!    values are needed. 
!---------------------------------------------------------------------
        if (renormalize_sw_fluxes) then
          if (do_sw_rad) then
            Astro%cosz    = cosz_r
            Astro%fracday = fracday_r
            Astro%solar   = solar_r
            Astro%rrsun   = rrsun_r
          endif
          call Astro_phys%alloc (size(lat,1), size(lat,2))
          Astro_phys%cosz    = cosz_p
          Astro_phys%fracday = fracday_p
          Astro_phys%solar   = solar_p
          Astro_phys%rrsun   = rrsun_p

!---------------------------------------------------------------------
!    if renormalization is active, then only the values applicable over
!    radiation steps are needed. 
!---------------------------------------------------------------------
        else                 
          if (do_sw_rad) then
            Astro%cosz    = cosz_r
            Astro%fracday = fracday_r
            Astro%solar   = solar_r
            Astro%rrsun   = rrsun_r
          endif
        endif

!---------------------------------------------------------------------
!    when in the gcm and on a radiation calculation step, define cosine
!    of zenith angle valid over the next radiation step. this is needed 
!    so that the ocean albedo (function of zenith angle) may be properly
!    defined and provided as input to the radiation package on the next
!    timestep.
!----------------------------------------------------------------------
        if (do_sw_rad) then
          call diurnal_solar (lat, lon, Rad_time+Dt_zen, cosz_a,   &
                              fracday_a, rrsun_a, dt_time=Dt_zen)
          Rad_output%coszen_angle(is:ie,js:je) = cosz_a(:,:)
        endif  ! (do_sw_rad)

!---------------------------------------------------------------------
!    case 2: annual-mean shortwave radiation.
!---------------------------------------------------------------------
      else if (Rad_control%do_annual) then
        call annual_mean_solar (js, je, lat, Astro%cosz, Astro%solar,&
                                Astro%fracday, Astro%rrsun)

!---------------------------------------------------------------------
!    save the cosine of zenith angle on the current step to be used to 
!    calculate ocean albedo for use on the next radiation timestep.
!---------------------------------------------------------------------
        Rad_output%coszen_angle(is:ie,js:je) = Astro%cosz(:,:)

!---------------------------------------------------------------------
!    case 3: daily-mean shortwave radiation.
!---------------------------------------------------------------------
      else if (Rad_control%do_daily_mean) then
        call daily_mean_solar (lat, Rad_time, Astro%cosz,  &
                               Astro%fracday, Astro%rrsun)
        Astro%solar = Astro%cosz*Astro%rrsun*Astro%fracday

!---------------------------------------------------------------------
!    save the cosine of zenith angle on the current step to be used to 
!    calculate ocean albedo for use on the next radiation timestep.
!---------------------------------------------------------------------
        Rad_output%coszen_angle(is:ie,js:je) = Astro%cosz(:,:)

!----------------------------------------------------------------------
!    if none of the above options are active, write an error message and
!    stop execution.
!----------------------------------------------------------------------
      else
        call error_mesg('radiation_driver_mod', &
             ' no valid zenith angle specification', FATAL)
      endif

!-------------------------------------------------------------------


end subroutine obtain_astronomy_variables 


!####################################################################
! <SUBROUTINE NAME="radiation_calc">
!  <OVERVIEW>
!    radiation_calc is called on radiation timesteps and calculates
!    the long- and short-wave radiative fluxes and heating rates, and
!    obtains the radiation output fields needed in other portions of
!    the model.
!  </OVERVIEW>
!  <DESCRIPTION>
!    radiation_calc is called on radiation timesteps and calculates
!    the long- and short-wave radiative fluxes and heating rates, and
!    obtains the radiation output fields needed in other portions of
!    the model.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiation_calc (is, ie, js, je,  &
!                           lat, lon, Atmos_input, Surface, &
!                           Astro, Rad_gases, aerooptdep, &
!                           aeroasymfac, aerosctopdep, aeroextopdep, &
!                           crndlw, cmxolw, emrndlw, emmxolw, &
!                           camtsw, cldsct, cldext, cldasymm, &
!                           Rad_output, Lw_output, Sw_output)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending i,j indices in global storage arrays
!  </IN> 
!  <IN NAME="Time_diag" TYPE="time_type">
!      Time_diag         time on next timestep, used as stamp for diag-
!                        nostic output  [ time_type  (days, seconds) ] 
!  </IN>
!  <IN NAME="lon" TYPE="real">
!    lon        mean longitude (in radians) of all grid boxes processed by
!               this call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <IN NAME="lat" TYPE="real">
!    lat        mean latitude (in radians) of all grid boxes processed by this
!               call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   Atmospheric input data to radiation package
!  </IN>
!  <IN NAME="Surface" TYPE="surface_type">
!   Surface input data to radiation package
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
!   astronomical input data for the radiation package
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!   Radiative gases properties to radiation package, , contains var-
!                     iables defining the radiatively active gases, 
!                     passed through to lower level routines
!  </IN>
!  <IN NAME="aerooptdep" TYPE="real">
!   Longwave aerosol optical depth (mean in model layers)
!  </IN>
!  <IN NAME="crndlw" TYPE="real">
!   Longwave cloud amount for random overlap clouds
!  </IN>
!  <IN NAME="cmxolw" TYPE="real">
!   Longwave cloud amount for maximum overlap clouds
!  </IN>
!  <IN NAME="emrndlw" TYPE="real">
!   cloud emissivity for random overlap clouds
!   by longwave band and profile
!  </IN>
!  <IN NAME="emmxolw" TYPE="real">
!   cloud emissivity for maximum overlap clouds
!   by longwave band and profile
!  </IN>
!  <IN NAME="camtsw" TYPE="real">
!   Cloud amount for shortwave clouds. If stochastic clouds is implemented
!   then cloud amount by band.
!  </IN>
!  <IN NAME="cldext" TYPE="real">
!   Cloud extinction parameter
!  </IN>
!  <IN NAME="cldsct" TYPE="real">
!   Cloud single scattering albedo
!  </IN>
!  <IN NAME="cldasymm" TYPE="real">
!   Cloud asymmetric parameter
!  </IN>
!  <INOUT NAME="Rad_output" TYPE="rad_output_type">
!   Radiation output from radiation package, contains variables
!                     which are output from radiation_driver to the 
!                     calling routine, and then used elsewhere within
!                     the component models.
!  </INOUT>
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
!
subroutine radiation_calc (is, ie, js, je, lat, lon, &
                           press, pflux, temp, tflux, rh2o, deltaz, tsfc, &
                           asfc_vis_dir, asfc_nir_dir, &
                           asfc_vis_dif, asfc_nir_dif, Astro, Rad_gases, &
                           aerooptdep, aerooptdep_volc, &
                           aeroasymfac, aerosctopdep, aeroextopdep, &
                           crndlw, cmxolw, emrndlw, emmxolw, &
                           camtsw, cldsct, cldext, cldasymm, &
                           Rad_output, Lw_output, Sw_output, Lw_diagnostics)

!--------------------------------------------------------------------
!    radiation_calc is called on radiation timesteps and calculates
!    the long- and short-wave radiative fluxes and heating rates, and
!    obtains the radiation output fields needed in other portions of
!    the model.
!-----------------------------------------------------------------------

!--------------------------------------------------------------------
integer,                      intent(in)             :: is, ie, js, je
real, dimension(:,:),         intent(in)             :: lat, lon
real, dimension(:,:,:),       intent(in)             :: press, pflux, temp, &
                                                        tflux, rh2o, deltaz
real, dimension(:,:),         intent(in)             :: tsfc
real, dimension(:,:),         intent(in)             :: asfc_vis_dir, &
                                                        asfc_nir_dir, &
                                                        asfc_vis_dif, &
                                                        asfc_nir_dif
type(astronomy_type),         intent(in)             :: Astro
type(radiative_gases_type),   intent(inout)          :: Rad_gases
real, dimension(:,:,:,:),     intent(in)             :: aerooptdep, &
                                                        aerooptdep_volc, &
                                                        aeroasymfac, &
                                                        aerosctopdep, &
                                                        aeroextopdep
real, dimension(:,:,:,:),     intent(in)             :: crndlw
real, dimension(:,:,:),       intent(in)             :: cmxolw
real, dimension(:,:,:,:,:),   intent(in)             :: emrndlw, emmxolw
real, dimension(:,:,:,:),     intent(in)             :: camtsw
real, dimension(:,:,:,:,:),   intent(in)             :: cldsct, cldext, cldasymm
type(rad_output_type),         intent(inout)         :: Rad_output
type(lw_output_type), dimension(:), intent(inout)    :: Lw_output
type(sw_output_type), dimension(:), intent(inout)    :: Sw_output
type(lw_diagnostics_type),          intent(inout)    :: Lw_diagnostics

!-----------------------------------------------------------------------
!    intent(in) variables:
!
!      is,ie,js,je       starting/ending subdomain i,j indices of data 
!                        in the physics_window being integrated
!      lat               latitude of model points on model grid 
!                        [ radians ]
!      lon               longitude of model points on model grid 
!                        [ radians ]
!      Atmos_input       atmospheric input data for the radiation 
!                        package 
!                        [ atmos_input_type ]
!      Surface           surface input data to the radiation package
!                        [ surface_type ]
!      Rad_gases         radiative gas input data for the radiation 
!                        package
!                        [ radiative_gases_type ]
!      aerooptdep        longwave aerosol optical depth in model layers
!      aeroasymfac, aerosctopdep, aeroextopdep
!                        shortwave aerosol radiative properties
!      crndlw            longwave cloud amount for random overlap clouds
!      cmxolw            longwave cloud amount for maximum overlap clouds
!      emrndlw           longwave cloud emissivity for random overlap clouds by band
!      emmxolw           longwave cloud emissivity for maximum overlap clouds by band
!      camtsw            shortwave cloud amount,
!                        if stochastic clouds is implemented then cloud amount by band
!      cldsct, cldext, cldasymm
!                        shortwave cloud radiative properties
!      Astro             astronomical input data for the radiation 
!                        package 
!                        [ astronomy_type ]
!
!
!    intent(out) variables:
!
!      Rad_output        radiation output data needed by other modules
!                        [ rad_output_type ]
!      Lw_output         longwave radiation output data from the 
!                        sea_esf_rad radiation package, when that 
!                        package is active
!                        [ lw_output_type ]
!          The following are the components of Lw_output:
!                 flxnet    net longwave flux at model flux levels 
!                           (including the ground and the top of the 
!                           atmosphere).
!                 heatra    longwave heating rates in model layers.
!                 flxnetcf  net longwave flux at model flux levels 
!                           (including the ground and the top of the 
!                           atmosphere) computed for cloud-free case.
!                 heatra    longwave heating rates in model layers 
!                           computed for cloud-free case.
!      Sw_output         shortwave radiation output data from the 
!                        sea_esf_rad radiation package  when that 
!                        package is active
!                        [ sw_output_type ]
!
!----------------------------------------------------------------------

      integer :: kmax

!----------------------------------------------------------------------
!    compute longwave radiation
!----------------------------------------------------------------------
      if (do_lw_rad) then
        call mpp_clock_begin (longwave_clock)
        call longwave_driver (press, pflux, temp, tflux, rh2o, deltaz,  &
                              Rad_gases, emrndlw, emmxolw, crndlw, cmxolw, &
                              aerooptdep, aerooptdep_volc, &
                              flag_stoch, Rad_control, &
                              Aerosolrad_control%do_lwaerosol, &
                              Aerosolrad_control%volcanic_lw_aerosols, &
                              Lw_output, Lw_diagnostics)
        call mpp_clock_end (longwave_clock)
      endif

!----------------------------------------------------------------------
!    compute shortwave radiation
!----------------------------------------------------------------------
      if (do_sw_rad) then
        call mpp_clock_begin (shortwave_clock)
        call shortwave_driver (press, pflux, temp, rh2o, deltaz, &
                               asfc_vis_dir, asfc_nir_dir, &
                               asfc_vis_dif, asfc_nir_dif, Astro, &
                               aeroasymfac, aerosctopdep, aeroextopdep, &
                               Rad_gases, camtsw, cldsct, cldext, cldasymm, &
                               flag_stoch, Rad_control, &
                               Aerosolrad_control%do_swaerosol, Sw_output)
        call mpp_clock_end (shortwave_clock)
      endif

!---------------------------------------------------------------------
!    define the components of Rad_output to be passed back to 
!    radiation_driver --  total and shortwave radiative heating rates 
!    for standard and clear-sky case (if desired), and surface long- 
!    and short-wave fluxes.
!---------------------------------------------------------------------
      if (do_sw_rad) then
        Rad_output%tdtsw(is:ie,js:je,:) =    &
                          Sw_output(1)%hsw(:,:,:)/SECONDS_PER_DAY
        Rad_output%ufsw(is:ie,js:je,:) =    &
                          Sw_output(1)%ufsw(:,:,:)
        Rad_output%dfsw(is:ie,js:je,:) =    &
                          Sw_output(1)%dfsw(:,:,:)
      endif
      if (do_lw_rad) then
         Rad_output%tdtlw(is:ie,js:je,:) =   &
                       Lw_output(1)%heatra(:,:,:)/SECONDS_PER_DAY
         Rad_output%flxnet(is:ie,js:je,:) =  &
                    Lw_output(1)%flxnet(:,:,:)
      endif

      Rad_output%tdt_rad (is:ie,js:je,:) =  &
                     (Rad_output%tdtsw(is:ie,js:je,:) +   &
                               Rad_output%tdtlw(is:ie,js:je,:))

      if (do_clear_sky_pass) then
        if (do_sw_rad) then
          Rad_output%tdtsw_clr(is:ie,js:je,:) =   &
                     Sw_output(1)%hswcf(:,:,:)/SECONDS_PER_DAY
          Rad_output%ufsw_clr(is:ie,js:je,:) =   &
                     Sw_output(1)%ufswcf(:,:,:)
          Rad_output%dfsw_clr(is:ie,js:je,:) =   &
                     Sw_output(1)%dfswcf(:,:,:)
          Rad_output%flux_sw_down_total_dir_clr(is:ie,js:je) =&
                         Sw_output(1)%dfsw_dir_sfc_clr(:,:)
          Rad_output%flux_sw_down_total_dif_clr(is:ie,js:je) =&
                          Sw_output(1)%dfsw_dif_sfc_clr(:,:)
          Rad_output%flux_sw_down_vis_clr(is:ie,js:je) =   &
                           Sw_output(1)%dfsw_vis_sfc_clr(:,:)
        endif
        if (do_lw_rad) then
          Rad_output%tdtlw_clr(is:ie,js:je,:) =   &
                  Lw_output(1)%heatracf(:,:,:)/SECONDS_PER_DAY
          Rad_output%flxnetcf(is:ie,js:je,:) =  &
                        Lw_output(1)%flxnet(:,:,:)
        endif
        Rad_output%tdt_rad_clr(is:ie,js:je,:) =    &
                   (Rad_output%tdtsw_clr(is:ie,js:je,:) +  &
                          Rad_output%tdtlw_clr(is:ie,js:je,:))
      endif

      kmax = size (Rad_output%tdtsw,3)
      if (do_sw_rad) then
        Rad_output%flux_sw_surf(is:ie,js:je) =   &
                        Sw_output(1)%dfsw(:,:,kmax+1) - &
                                Sw_output(1)%ufsw(:,:,kmax+1)
        Rad_output%flux_sw_surf_dir(is:ie,js:je) =   &
                                Sw_output(1)%dfsw_dir_sfc(:,:)
        Rad_output%flux_sw_surf_refl_dir(is:ie,js:je) =   &
                                Sw_output(1)%ufsw_dir_sfc(:,:)
        Rad_output%flux_sw_surf_dif(is:ie,js:je) =   &
                             Sw_output(1)%dfsw_dif_sfc(:,:) - &
                                  Sw_output(1)%ufsw_dif_sfc(:,:)
        Rad_output%flux_sw_down_vis_dir(is:ie,js:je) =   &
                               Sw_output(1)%dfsw_vis_sfc_dir(:,:)
        Rad_output%flux_sw_down_vis_dif(is:ie,js:je) =   &
                               Sw_output(1)%dfsw_vis_sfc_dif(:,:)
        Rad_output%flux_sw_down_total_dir(is:ie,js:je) =   &
                                   Sw_output(1)%dfsw_dir_sfc(:,:)
        Rad_output%flux_sw_down_total_dif(is:ie,js:je) =   &
                                  Sw_output(1)%dfsw_dif_sfc(:,:)
        Rad_output%flux_sw_vis (is:ie,js:je) =   &
                           Sw_output(1)%dfsw_vis_sfc(:,:) - &
                                    Sw_output(1)%ufsw_vis_sfc(:,:)
        Rad_output%flux_sw_vis_dir (is:ie,js:je) =   &
                               Sw_output(1)%dfsw_vis_sfc_dir(:,:)
        Rad_output%flux_sw_refl_vis_dir (is:ie,js:je) =   &
                               Sw_output(1)%ufsw_vis_sfc_dir(:,:)
        Rad_output%flux_sw_vis_dif (is:ie,js:je) =   &
                        Sw_output(1)%dfsw_vis_sfc_dif(:,:) - &
                               Sw_output(1)%ufsw_vis_sfc_dif(:,:)
      endif
      if (do_lw_rad) then
        Rad_output%flux_lw_surf(is:ie,js:je) =    &
                     STEFAN*tsfc(:,:)**4 -   &
                                   Lw_output(1)%flxnet(:,:,kmax+1)
      endif

!---------------------------------------------------------------------


end subroutine radiation_calc


!######################################################################
! <SUBROUTINE NAME="deallocate_arrays">
!  <OVERVIEW>
!    deallocate_arrays deallocates the array space of local 
!    derived-type variables.
!  </OVERVIEW>
!  <DESCRIPTION>
!    deallocate_arrays deallocates the array space of local 
!    derived-type variables.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call deallocate_arrays (Astro, Astro_phys, Lw_output, Sw_output)
!  </TEMPLATE>
!  <INOUT NAME="Astro" TYPE="astronomy_type">
!   astronomical data for the radiation package
!  </INOUT>
!  <INOUT NAME="Astro_phys" TYPE="astronomy_type">
!   astronomical data for the radiation package
!  </INOUT>
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
!
subroutine deallocate_arrays (Astro, Astro_phys,  &
                              Sw_output, Lw_output, &
                              Lw_diagnostics, &
                              Aerosol, Aerosolrad_diags)

!---------------------------------------------------------------------
!    deallocate_arrays deallocates the array space of local 
!    derived-type variables.
!---------------------------------------------------------------------

type(astronomy_type),              intent(inout)  :: Astro, Astro_phys
type(sw_output_type),dimension(:), intent(inout)  :: Sw_output
type(lw_output_type),dimension(:), intent(inout)  :: Lw_output
type(lw_diagnostics_type),         intent(inout)  :: Lw_diagnostics
type(aerosol_type),                intent(inout)  :: Aerosol
type(aerosolrad_diag_type),        intent(inout)  :: Aerosolrad_diags

      integer  ::  n

!--------------------------------------------------------------------
!    deallocate the variables in Astro and Astro_phys.
!--------------------------------------------------------------------
      if ( do_rad .or. renormalize_sw_fluxes ) then 
        call Astro%dealloc
        if ( do_sw_rad .and. renormalize_sw_fluxes .and. Rad_control%do_diurnal ) then 
            call Astro_phys%dealloc
        endif
      endif

!--------------------------------------------------------------------
!    deallocate the variables in Lw_output.
!--------------------------------------------------------------------
      if (do_lw_rad) then
        do n = 1, Aerosolrad_control%indx_lwaf
          call Lw_output(n)%dealloc
        end do
        call Lw_diagnostics%dealloc
      endif

!--------------------------------------------------------------------
!    deallocate the variables in Sw_output.
!--------------------------------------------------------------------
      if (do_sw_rad) then
        do n = 1, Aerosolrad_control%indx_swaf
          call Sw_output(n)%dealloc
        end do
      endif

!--------------------------------------------------------------------
!    deallocate the window-resident variables in Aerosolrad_diags
!--------------------------------------------------------------------
      if (do_rad .and. Aerosolrad_control%do_aerosol) then
        call aerosolrad_driver_dealloc (Aerosol, Aerosolrad_diags)
      endif

!---------------------------------------------------------------------


end subroutine deallocate_arrays 


!#####################################################################
! <SUBROUTINE NAME="calculate_auxiliary_variables">
!  <OVERVIEW>
!    calculate_auxiliary_variables defines values of model delta z and
!    relative humidity, and the values of pressure and temperature at
!    the grid box vertical interfaces.
!  </OVERVIEW>
!  <DESCRIPTION>
!    calculate_auxiliary_variables defines values of model delta z and
!    relative humidity, and the values of pressure and temperature at
!    the grid box vertical interfaces.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call calculate_auxiliary_variables (Atmos_input)
!  </TEMPLATE>
!  <INOUT NAME="Atmos_input" TYPE="atmos_input_type">
!   atmos_input_type variable, its press and temp
!                   components are input, and its deltaz, relhum, 
!                   pflux, tflux and aerosolrelhum components are 
!                   calculated here and output.
!  </INOUT>
! </SUBROUTINE>
!
subroutine calculate_auxiliary_variables (Atmos_input)

!----------------------------------------------------------------------
!    calculate_auxiliary_variables defines values of model delta z and
!    relative humidity, and the values of pressure and temperature at
!    the grid box vertical interfaces.
!---------------------------------------------------------------------

type(atmos_input_type), intent(inout)  :: Atmos_input

!--------------------------------------------------------------------
!   intent(inout) variables
!
!      Atmos_input  atmos_input_type variable, its press and temp
!                   components are input, and its deltaz, relhum, 
!                   pflux, tflux and aerosolrelhum components are 
!                   calculated here and output.
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables

      real, dimension (size(Atmos_input%temp, 1), &
                       size(Atmos_input%temp, 2), &
                       size(Atmos_input%temp, 3)) :: &
                                                   qsat, qv, tv
      integer   ::  k
      integer   ::  id, jd, kmax



!---------------------------------------------------------------------
!    allocate space for the auxillary components of the 
!    atmos_input_type variable
!---------------------------------------------------------------------
      id   = size(Atmos_input%temp,1)
      jd   = size(Atmos_input%temp,2)
      kmax = size(Atmos_input%temp,3)
      allocate ( Atmos_input%tflux (id,jd,kmax+1) )
      allocate ( Atmos_input%deltaz(id, jd,kmax) )

      allocate ( Atmos_input%relhum(id,jd,kmax) )
      allocate ( Atmos_input%clouddeltaz(id,jd,kmax) )
      allocate ( Atmos_input%aerosolrelhum(id,jd,kmax) )

      if (do_conserve_energy) then
        Atmos_input%pflux => Atmos_input%phalf
      else
!--------------------------------------------------------------------
!    define flux level pressures (pflux) as midway between data level
!    (layer-mean) pressures.
!--------------------------------------------------------------------
        allocate ( Atmos_input%pflux (id,jd,kmax+1) )
        do k=ks+1,ke
          Atmos_input%pflux(:,:,k) = 0.5E+00*  &
                  (Atmos_input%press(:,:,k-1) + Atmos_input%press(:,:,k))
        enddo
        Atmos_input%pflux(:,:,ks  ) = 0.0E+00
        Atmos_input%pflux(:,:,ke+1) = Atmos_input%psfc(:,:)
      endif
!--------------------------------------------------------------------
!    specify temperatures at flux levels (tflux)
!--------------------------------------------------------------------
      do k=ks+1,ke
        Atmos_input%tflux(:,:,k) = 0.5E+00*  &
                (Atmos_input%temp (:,:,k-1) + Atmos_input%temp (:,:,k))
      enddo
      Atmos_input%tflux(:,:,ks  ) = Atmos_input%temp(:,:,ks)
      Atmos_input%tflux(:,:,ke+1) = Atmos_input%tsfc(:,:)

!-------------------------------------------------------------------
!    define deltaz in meters.
!-------------------------------------------------------------------
      tv(:,:,:) = Atmos_input%temp(:,:,ks:ke)*    &
                  (1.0 + D608*Atmos_input%rh2o(:,:,:))
      Atmos_input%deltaz(:,:,ks) = log_p_at_top*RDGAS*tv(:,:,ks)/GRAV
      do k =ks+1,ke   
        Atmos_input%deltaz(:,:,k) = alog(Atmos_input%pflux(:,:,k+1)/  &
                                         Atmos_input%pflux(:,:,k))*   &
                                         RDGAS*tv(:,:,k)/GRAV
      enddo

!-------------------------------------------------------------------
!    define deltaz in meters to be used in cloud feedback analysis.
!-------------------------------------------------------------------
      tv(:,:,:) = Atmos_input%cloudtemp(:,:,ks:ke)*    &
                  (1.0 + D608*Atmos_input%cloudvapor(:,:,:))
      Atmos_input%clouddeltaz(:,:,ks) = log_p_at_top*RDGAS*  &
                                        tv(:,:,ks)/GRAV
      do k =ks+1,ke   
        Atmos_input%clouddeltaz(:,:,k) =    &
                            alog(Atmos_input%pflux(:,:,k+1)/  &
                                         Atmos_input%pflux(:,:,k))*   &
                                         RDGAS*tv(:,:,k)/GRAV
      end do

!------------------------------------------------------------------
!    define the relative humidity.
!------------------------------------------------------------------
      qv(:,:,1:kmax) = Atmos_input%rh2o(:,:,1:kmax) /    &
                                   (1.0 + Atmos_input%rh2o(:,:,1:kmax))
      call compute_qs (Atmos_input%temp(:,:,1:kmax),  &
                       Atmos_input%press(:,:,1:kmax),  &
                       qsat(:,:,1:kmax), q = qv(:,:,1:kmax))
      do k=1,kmax
        Atmos_input%relhum(:,:,k) = qv(:,:,k) / qsat(:,:,k)
        Atmos_input%relhum(:,:,k) = MIN (Atmos_input%relhum(:,:,k), 1.0)
      end do

!------------------------------------------------------------------
!    define the relative humidity seen by the aerosol code.
!------------------------------------------------------------------
        qv(:,:,1:kmax) = Atmos_input%aerosolvapor(:,:,1:kmax) /    &
                         (1.0 + Atmos_input%aerosolvapor(:,:,1:kmax))
        call compute_qs (Atmos_input%aerosoltemp(:,:,1:kmax), &
                         Atmos_input%aerosolpress(:,:,1:kmax),  &
                         qsat(:,:,1:kmax), q = qv(:,:,1:kmax))
      do k=1,kmax
         Atmos_input%aerosolrelhum(:,:,k) = qv(:,:,k) / qsat(:,:,k)
         Atmos_input%aerosolrelhum(:,:,k) =    &
                            MIN (Atmos_input%aerosolrelhum(:,:,k), 1.0)
      end do
 
!----------------------------------------------------------------------


end subroutine calculate_auxiliary_variables


!#######################################################################

subroutine set_radiation_diag_type ( press, pflux, temp, rh2o, &
             asfc_vis_dir, asfc_nir_dir, &
             asfc_vis_dif, asfc_nir_dif, Astro, &
             Rad_gases, crndlw, cmxolw, emrndlw, emmxolw, &
             camtsw, cldsct, cldext, cldasymm, &
             Sw_output, Lw_output, Lw_diagnostics, &
             Rad_diag)

!--------------------------------------------------------------------
!  allocates and assigns data to the radiation_diag_type
!  this data type is used by the radiation_diag module
!--------------------------------------------------------------------

real, dimension(:,:,:),          intent(in) :: press, pflux, temp, rh2o
real, dimension(:,:),            intent(in) :: asfc_vis_dir, asfc_nir_dir, &
                                               asfc_vis_dif, asfc_nir_dif
type(astronomy_type),            intent(in) :: Astro
type(radiative_gases_type),      intent(in) :: Rad_gases
real, dimension(:,:,:,:),        intent(in) :: crndlw
real, dimension(:,:,:),          intent(in) :: cmxolw
real, dimension(:,:,:,:,:),      intent(in) :: emrndlw, emmxolw
real, dimension(:,:,:,:),        intent(in) :: camtsw
real, dimension(:,:,:,:,:),      intent(in) :: cldsct, cldext, cldasymm
type(sw_output_type), dimension(:), intent(in)  :: Sw_output
type(lw_output_type), dimension(:), intent(in)  :: Lw_output
type(lw_diagnostics_type),          intent(in)  :: Lw_diagnostics
type(radiation_diag_type),          intent(out) :: Rad_diag

!--------------------------------------------------------------------
! local variables

integer :: id, jd, kd, n

!--------------------------------------------------------------------

    Rad_diag%do_radiation_diag = .true.

    id = size(rh2o,1)
    jd = size(rh2o,2)
    kd = size(rh2o,3)
      
  ! aerosol and cloud properties
    allocate(Rad_diag%emrndlw (id,jd,kd,size(emrndlw,4),size(emrndlw,5)))
    allocate(Rad_diag%emmxolw (id,jd,kd,size(emmxolw,4),size(emmxolw,5)))
    allocate(Rad_diag%cldsct  (id,jd,kd,size(cldsct ,4),size(cldsct ,5)))
    allocate(Rad_diag%cldext  (id,jd,kd,size(cldext ,4),size(cldext ,5)))
    allocate(Rad_diag%cldasymm(id,jd,kd,size(cldasymm,4),size(cldasymm,5)))
    allocate(Rad_diag%camtsw  (id,jd,kd,size(camtsw,4)))
    allocate(Rad_diag%crndlw  (id,jd,kd,size(crndlw,4)))
    allocate(Rad_diag%cmxolw  (id,jd,kd))
    Rad_diag%emrndlw  = emrndlw
    Rad_diag%emmxolw  = emmxolw
    Rad_diag%cldsct   = cldsct
    Rad_diag%cldext   = cldext
    Rad_diag%cldasymm = cldasymm
    Rad_diag%camtsw   = camtsw
    Rad_diag%crndlw   = crndlw
    Rad_diag%cmxolw   = cmxolw

  ! atmospheric state variables
     allocate(Rad_diag%press(id,jd,size(press,3)))
     allocate(Rad_diag%pflux(id,jd,size(pflux,3)))
     allocate(Rad_diag%temp (id,jd,size(temp,3)))
     allocate(Rad_diag%rh2o (id,jd,kd))
     allocate(Rad_diag%qo3  (id,jd,kd))
     Rad_diag%press = press
     Rad_diag%pflux = pflux
     Rad_diag%temp  = temp
     Rad_diag%rh2o  = rh2o
     Rad_diag%qo3   = Rad_gases%qo3

  ! surface albedo
    allocate(Rad_diag%asfc_vis_dir(id,jd))
    allocate(Rad_diag%asfc_nir_dir(id,jd))
    allocate(Rad_diag%asfc_vis_dif(id,jd))
    allocate(Rad_diag%asfc_nir_dif(id,jd))
    Rad_diag%asfc_vis_dir = asfc_vis_dir
    Rad_diag%asfc_nir_dir = asfc_nir_dir
    Rad_diag%asfc_vis_dif = asfc_vis_dif
    Rad_diag%asfc_nir_dif = asfc_nir_dif

  ! astronomy
    allocate(Rad_diag%cosz   (id,jd))
    allocate(Rad_diag%fracday(id,jd))
    Rad_diag%cosz    = Astro%cosz
    Rad_diag%fracday = Astro%fracday
    Rad_diag%rrsun   = Astro%rrsun
    call get_solar_constant (Rad_diag%solar_constant)
    Rad_diag%do_diurnal = Rad_control%do_diurnal
    Rad_diag%do_annual = Rad_control%do_annual
    Rad_diag%do_daily_mean = Rad_control%do_daily_mean

  ! radiative gases
    Rad_diag%rrvco2  = Rad_gases%rrvco2
    Rad_diag%rrvf11  = Rad_gases%rrvf11
    Rad_diag%rrvf12  = Rad_gases%rrvf12
    Rad_diag%rrvf113 = Rad_gases%rrvf113
    Rad_diag%rrvf22  = Rad_gases%rrvf22
    Rad_diag%rrvch4  = Rad_gases%rrvch4
    Rad_diag%rrvn2o  = Rad_gases%rrvn2o

  ! copy sw annd lw output/diagnostics
    allocate(Rad_diag%Sw_output(size(Sw_output)))
    allocate(Rad_diag%Lw_output(size(Lw_output)))
    do n = 1, size(Sw_output)
      call Rad_diag%Sw_output(n)%alloc (id, jd, kd, Rad_control%do_totcld_forcing)
      Rad_diag%Sw_output(n) = Sw_output(n)
    end do
    do n = 1, size(Lw_output)
      call Rad_diag%Lw_output(n)%alloc (id, jd, kd, Rad_control%do_totcld_forcing)
      Rad_diag%Lw_output(n) = Lw_output(n)
    end do
    call Rad_diag%Lw_diagnostics%alloc (id, jd, kd, &
                                        size(Lw_diagnostics%flx1e1f,3), &
                                        size(Lw_diagnostics%exctsn,4),  &
                                        Rad_control%do_totcld_forcing)
    Rad_diag%Lw_diagnostics = Lw_diagnostics

  ! pointer to longwave tables
    call longwave_get_tables(Rad_diag%Lw_tables)


  ! miscellaneous control variables
    Rad_diag%do_totcld_forcing = Rad_control%do_totcld_forcing
    Rad_diag%do_swaerosol = Aerosolrad_control%do_swaerosol
    Rad_diag%do_lwaerosol = Aerosolrad_control%do_lwaerosol
    Rad_diag%do_swaerosol_forcing = Aerosolrad_control%do_swaerosol_forcing
    Rad_diag%do_lwaerosol_forcing = Aerosolrad_control%do_lwaerosol_forcing
    Rad_diag%indx_swaf = Aerosolrad_control%indx_swaf
    Rad_diag%indx_lwaf = Aerosolrad_control%indx_lwaf

!--------------------------------------------------------------------

end subroutine set_radiation_diag_type

!#######################################################################

                 end module radiation_driver_mod

