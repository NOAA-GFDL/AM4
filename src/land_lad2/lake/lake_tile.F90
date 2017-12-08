#include <fms_platform.h>

module lake_tile_mod

use mpp_domains_mod, only : &
     domain2d, mpp_get_compute_domain, mpp_pass_sg_to_ug

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : file_exist, check_nml_error, read_data, close_file, stdlog
use constants_mod, only : PI, tfreeze, hlf
use land_constants_mod, only : NBANDS
use land_data_mod, only : log_version, lnd_sg, lnd
use land_io_mod, only : init_cover_field
use land_tile_selectors_mod, only : tile_selector_type, SEL_LAKE, register_tile_selector
use tiling_input_types_mod, only : lake_predefined_type
use land_debug_mod, only : is_watch_point

implicit none
private

! ==== public interfaces =====================================================
public :: lake_pars_type
public :: lake_tile_type

public :: new_lake_tile, delete_lake_tile
public :: new_lake_tile_predefined
public :: lake_tiles_can_be_merged, merge_lake_tiles
public :: lake_is_selected
public :: get_lake_tile_tag
public :: lake_tile_stock_pe
public :: lake_tile_heat

public :: read_lake_data_namelist
public :: lake_cover_cold_start

public :: lake_radiation
public :: lake_roughness
public :: lake_data_thermodynamics

public :: max_lev
public :: lake_width_inside_lake
public :: large_lake_sill_width
! =====end of public interfaces ==============================================
interface new_lake_tile
   module procedure lake_tile_ctor
   module procedure lake_tile_copy_ctor
end interface
interface new_lake_tile_predefined
   module procedure lake_tile_ctor_predefined
   module procedure lake_tile_copy_ctor
end interface

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'lake_tile_mod'
#include "../shared/version_variable.inc"

integer, parameter :: max_lev          = 80
integer, parameter :: n_dim_lake_types = 1  ! size of lookup table
real,    parameter :: psi_wilt         = -150.  ! matric head at wilting
real,    parameter :: comp             = 0.001  ! m^-1

! from the modis brdf/albedo product user's guide:
real, parameter :: g_iso  = 1.
real, parameter :: g_vol  = 0.189184
real, parameter :: g_geo  = -1.377622
real, parameter :: g0_iso = 1.0
real, parameter :: g1_iso = 0.0
real, parameter :: g2_iso = 0.0
real, parameter :: g0_vol = -0.007574
real, parameter :: g1_vol = -0.070987
real, parameter :: g2_vol =  0.307588
real, parameter :: g0_geo = -1.284909
real, parameter :: g1_geo = -0.166314
real, parameter :: g2_geo =  0.041840

! ==== types =================================================================
type :: lake_pars_type
  real w_sat
  real awc_lm2
  real k_sat_ref
  real psi_sat_ref
  real chb
  real alpha              ! *** REPLACE LATER BY alpha(layer)
  real heat_capacity_ref
  real thermal_cond_ref
  real refl_dry_dir(NBANDS)
  real refl_dry_dif(NBANDS)
  real refl_sat_dir(NBANDS)
  real refl_sat_dif(NBANDS)
  real emis_dry
  real emis_sat
  real z0_momentum
  real z0_momentum_ice
  real depth_sill
  real width_sill
  real whole_area
  real connected_to_next
  real backwater
  real backwater_1
  real rsa_exp         ! riparian source-area exponent
end type lake_pars_type

type :: lake_tile_type
   integer :: tag ! kind of the lake
   type(lake_pars_type)          :: pars

   real, allocatable :: dz(:)
   real, allocatable :: wl(:)
   real, allocatable :: ws(:)
   real, allocatable :: T(:)
   real, allocatable :: K_z(:)

   real, allocatable :: w_fc(:)
   real, allocatable :: w_wilt(:)
   real :: Eg_part_ref
   real :: z0_scalar
   real :: z0_scalar_ice
   real :: geothermal_heat_flux
   real, allocatable :: e(:),f(:)
   real, allocatable :: heat_capacity_dry(:)
end type lake_tile_type

! ==== module data ===========================================================
real, public :: &
     cpw = 1952.0, & ! specific heat of water vapor at constant pressure
     clw = 4218.0, & ! specific heat of water (liquid)
     csw = 2106.0    ! specific heat of water (ice)
logical, public :: use_brdf

!---- namelist ---------------------------------------------------------------
real    :: lake_width_inside_lake = 1.e5
real    :: large_lake_sill_width = 200.
real    :: min_lake_frac         = 0.
real    :: max_lake_rh           = 1.
real    :: lake_rh_exp           = 1.
real    :: dry_lake_depth_frac   = 0.
real    :: k_over_B              = 0.25      ! reset to 0 for MCM
real    :: k_over_B_ice          = 0.25
real    :: rate_fc               = 0.1/86400 ! 0.1 mm/d drainage rate at FC
real    :: sfc_heat_factor       = 1
real    :: geothermal_heat_flux_constant = 0.0  ! true continental average is ~0.065 W/m2
integer, public :: num_l         = 18           ! number of lake levels
integer, public :: n_outlet      = 0
logical :: use_lm2_awc           = .false.
logical, public :: lake_specific_width    = .false.
integer :: n_map_1st_lake_type = 10

integer, public :: outlet_face(100)
integer, public :: outlet_i(100)
integer, public :: outlet_j(100)
real, public :: outlet_width(100)

! from analysis of modis data (ignoring temperature dependence):
  real :: f_iso_ice(NBANDS) = (/ 0.056, 0.131 /)
  real :: f_vol_ice(NBANDS) = (/ 0.017, 0.053 /)
  real :: f_geo_ice(NBANDS) = (/ 0.004, 0.010 /)
  real :: f_iso_liq(NBANDS) = (/ 0.056, 0.131 /)
  real :: f_vol_liq(NBANDS) = (/ 0.017, 0.053 /)
  real :: f_geo_liq(NBANDS) = (/ 0.004, 0.010 /)

! ---- remainder are used only for cold start ---------
logical :: round_frac_down       = .false.  ! when false, any lake_frac < min_lake_frac
                                            ! is set to min_lake_frac.
                                            ! when true, any lake_frac < min_lake_frac
                                            ! is set to 0.

character(len=16):: lake_to_use     = 'single-tile'
       ! 'multi-tile' for tiled soil [default]
       ! 'single-tile' for geographically varying soil with single type per
       !     model grid cell
       ! 'uniform' for global constant soil, e.g., to reproduce MCM
       ! 'from-rivers' to get the fraction of lakes from river module
logical :: use_single_lake       = .false.   ! true for single global lake,
                                             ! e.g., to recover MCM
logical :: use_mcm_albedo        = .false.   ! .true. for CLIMAP albedo inputs
logical :: use_single_geo        = .false.   ! .true. for global gw res time,
                                             ! e.g., to recover MCM
integer :: lake_index_constant   = 1         ! index of global constant lake,
                                             ! used when use_single_lake
real    :: rsa_exp_global        = 1.5
real, dimension(n_dim_lake_types) :: &
  dat_w_sat             =(/ 1.000   /),&
  dat_awc_lm2           =(/ 1.000   /),&
  dat_k_sat_ref         =(/ 0.021   /),&
  dat_psi_sat_ref       =(/ -.059   /),&
  dat_chb               =(/   3.5   /),&
  dat_heat_capacity_ref =(/ 8.4e7   /),&
  dat_thermal_cond_ref  =(/ 8.4e7   /),&
  dat_emis_dry          =(/ 0.950   /),&
  dat_emis_sat          =(/ 0.980   /),&
  dat_z0_momentum       =(/ 1.4e-4  /),&
  dat_z0_momentum_ice   =(/ 1.4e-4  /),&
  dat_tf_depr           =(/  0.00   /)
real, dimension(n_dim_lake_types, NBANDS) :: &
     dat_refl_dry_dif, dat_refl_dry_dir, &
     dat_refl_sat_dif, dat_refl_sat_dir
data &
                   !  VIS    NIR
  dat_refl_dry_dif / 0.060, 0.060 /, &
  dat_refl_dry_dir / 0.060, 0.060 /, &
  dat_refl_sat_dir / 0.060, 0.060 /, &
  dat_refl_sat_dif / 0.060, 0.060 /
integer, dimension(n_dim_lake_types) :: &
  input_cover_types     =(/ 10 /)
character(len=4), dimension(n_dim_lake_types) :: &
  tile_names            =(/ 'lake' /)

namelist /lake_data_nml/ lake_width_inside_lake, &
     large_lake_sill_width, &
     lake_specific_width, n_outlet, outlet_face, outlet_i, outlet_j, outlet_width, &
     min_lake_frac, round_frac_down, max_lake_rh, lake_rh_exp, &
     dry_lake_depth_frac, &
     lake_to_use,input_cover_types, tile_names, &
     k_over_B, k_over_B_ice,         &
     rate_fc, sfc_heat_factor, geothermal_heat_flux_constant,  &
     num_l, &
     use_lm2_awc,    n_map_1st_lake_type, &
     use_single_lake,           use_mcm_albedo,            &
     use_single_geo,            lake_index_constant,         &
     rsa_exp_global,      &
     dat_w_sat,               dat_awc_lm2,     &
     dat_k_sat_ref,            &
     dat_psi_sat_ref,               dat_chb,          &
     dat_heat_capacity_ref,         dat_thermal_cond_ref,   &
     dat_refl_dry_dir,            dat_refl_sat_dir,              &
     dat_refl_dry_dif,            dat_refl_sat_dif,              &
     dat_emis_dry,              dat_emis_sat,                &
     dat_z0_momentum,   dat_z0_momentum_ice,        dat_tf_depr, &
     f_iso_ice, f_vol_ice, f_geo_ice, f_iso_liq, f_vol_liq, f_geo_liq


!---- end of namelist --------------------------------------------------------

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine read_lake_data_namelist(lake_n_lev)
  integer, intent(out) :: lake_n_lev
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: i
  real    :: z

  call log_version(version, module_name, &
  __FILE__)
#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=lake_data_nml, iostat=io)
     ierr = check_nml_error(io, 'lake_data_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=lake_data_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'lake_data_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  unit=stdlog()
  write(unit, nml=lake_data_nml)

  ! initialize global module data here

  ! register selectors for tile-specific diagnostics
  do i=1, n_dim_lake_types
     call register_tile_selector(tile_names(i), long_name='',&
          tag = SEL_LAKE, idata1 = i, area_depends_on_time=.FALSE. )
  enddo

  ! set up output arguments
  lake_n_lev = num_l
end subroutine read_lake_data_namelist

! ============================================================================
function lake_tile_ctor_predefined(tag,lake_predefined,itile) result(ptr)
  type(lake_tile_type), pointer :: ptr ! return value
  integer, intent(in)           :: tag ! kind of lake
  type(lake_predefined_type), intent(in) :: lake_predefined
  integer, intent(in) :: itile

  allocate(ptr)
  ptr%tag = tag
  ! allocate storage for tile data
  allocate(ptr%dz     (num_l), &
           ptr%wl     (num_l), &
           ptr%ws     (num_l), &
           ptr%T      (num_l), &
           ptr%K_z    (num_l), &
           ptr%w_fc   (num_l),  &
           ptr%w_wilt (num_l),  &
           ptr%heat_capacity_dry(num_l),  &
           ptr%e      (num_l),  &
           ptr%f      (num_l)   )
  call init_lake_data_0d_predefined(ptr,lake_predefined,itile)

end function lake_tile_ctor_predefined

! ============================================================================
function lake_tile_ctor(tag) result(ptr)
  type(lake_tile_type), pointer :: ptr ! return value
  integer, intent(in)           :: tag ! kind of lake

  allocate(ptr)
  ptr%tag = tag
  ! allocate storage for tile data
  allocate(ptr%dz     (num_l), &
           ptr%wl     (num_l), &
           ptr%ws     (num_l), &
           ptr%T      (num_l), &
           ptr%K_z    (num_l), &
           ptr%w_fc   (num_l),  &
           ptr%w_wilt (num_l),  &
           ptr%heat_capacity_dry(num_l),  &
           ptr%e      (num_l),  &
           ptr%f      (num_l)   )
  call init_lake_data_0d(ptr)

end function lake_tile_ctor


! ============================================================================
function lake_tile_copy_ctor(lake) result(ptr)
  type(lake_tile_type), pointer    :: ptr  ! return value
  type(lake_tile_type), intent(in) :: lake ! tile to copy

  allocate(ptr)
  ptr=lake ! copy all non-pointer data
  ! no need to allocate storage for allocatable components of the tile, because
  ! F2003 takes care of that, and also takes care of copying data
end function lake_tile_copy_ctor


! ============================================================================
subroutine delete_lake_tile(ptr)
  type(lake_tile_type), pointer :: ptr

  ! no need to deallocate components of the tile, because F2003 takes care of
  ! allocatable components deallocation when tile is deallocated
  deallocate(ptr)
end subroutine delete_lake_tile

subroutine init_lake_data_0d_predefined(lake,lake_predefined,itile)
  type(lake_tile_type), intent(inout) :: lake
  type(lake_predefined_type), intent(in) :: lake_predefined
  integer, intent(in) :: itile

!  real rsa_exp         ! riparian source-area exponent

  integer :: k
  k = lake%tag

  !Variables previously defined in gridded data
  lake%pars%connected_to_next = lake_predefined%connected_to_next(itile)
  lake%pars%whole_area = lake_predefined%whole_area(itile)
  lake%pars%depth_sill = lake_predefined%depth_sill(itile)
  lake%pars%backwater = lake_predefined%backwater(itile)
  lake%pars%backwater_1 = lake_predefined%backwater_1(itile)
  lake%pars%width_sill= lake_predefined%width_sill(itile)

  !Variables previosuly defined in look-up tables
  lake%pars%w_sat             = lake_predefined%w_sat(itile)
  lake%pars%awc_lm2           = lake_predefined%awc_lm2(itile)
  lake%pars%k_sat_ref         = lake_predefined%k_sat_ref(itile)
  lake%pars%psi_sat_ref       = lake_predefined%psi_sat_ref(itile)
  lake%pars%chb               = lake_predefined%chb(itile)
  lake%pars%alpha             = lake_predefined%alpha(itile)
  lake%pars%heat_capacity_ref = lake_predefined%heat_capacity_ref(itile)
  lake%pars%thermal_cond_ref  = lake_predefined%thermal_cond_ref(itile)
  lake%pars%refl_dry_dir      = lake_predefined%refl_dry_dir(itile,:)
  lake%pars%refl_dry_dif      = lake_predefined%refl_dry_dif(itile,:)
  lake%pars%refl_sat_dir      = lake_predefined%refl_sat_dir(itile,:)
  lake%pars%refl_sat_dif      = lake_predefined%refl_sat_dif(itile,:)
  lake%pars%emis_dry          = lake_predefined%emis_dry(itile)
  lake%pars%emis_sat          = lake_predefined%emis_sat(itile)
  lake%pars%z0_momentum       = lake_predefined%z0_momentum(itile)
  lake%pars%z0_momentum_ice   = lake_predefined%z0_momentum_ice(itile)

  if (is_watch_point()) then
   print*,'lake parameters'
   print*,'connected_to_next: ',lake%pars%connected_to_next
   print*,'whole_area: ',lake%pars%whole_area
   print*,'depth_sill: ',lake%pars%depth_sill
   print*,'backwater: ',lake%pars%backwater
   print*,'backwater_1: ',lake%pars%backwater_1
   print*,'width_sill: ',lake%pars%width_sill
   print*,'w_sat o:',dat_w_sat(k),' n:',lake%pars%w_sat
   print*,'awc_lm2 o:',dat_awc_lm2(k),' n:',lake%pars%awc_lm2
   print*,'k_sat_ref o:',dat_k_sat_ref(k),' n:',lake%pars%k_sat_ref
   print*,'psi_sat_ref o:',dat_psi_sat_ref(k),' n:',lake%pars%psi_sat_ref
   print*,'chb o:',dat_chb(k),' n:',lake%pars%chb
   print*,'alpha o:',1.0,' n:',lake%pars%alpha
   print*,'heat_capacity_ref o:',dat_heat_capacity_ref(k),' n:',lake%pars%heat_capacity_ref
   print*,'thermal_cond_ref o:',dat_thermal_cond_ref(k),' n:',lake%pars%thermal_cond_ref
   print*,'refl_dry_dir o:',dat_refl_dry_dir(k,:),' n:',lake%pars%refl_dry_dir
   print*,'refl_dry_dif o:',dat_refl_dry_dif(k,:),' n:',lake%pars%refl_dry_dif
   print*,'refl_sat_dir o:',dat_refl_sat_dir(k,:),' n:',lake%pars%refl_sat_dir
   print*,'refl_sat_dif o:',dat_refl_sat_dif(k,:),' n:',lake%pars%refl_sat_dif
   print*,'emis_dry o:',dat_emis_dry(k),' n:',lake%pars%emis_dry 
   print*,'emis_sat o:',dat_emis_sat(k),' n:',lake%pars%emis_sat
   print*,'z0_momentum o:',dat_z0_momentum(k),' n:',lake%pars%z0_momentum
   print*,'z0_momentum_ice o:',dat_z0_momentum_ice(k),'n:',lake%pars%z0_momentum_ice
  endif

  lake%pars%rsa_exp           = rsa_exp_global

  ! -------- derived constant lake parameters --------
  ! w_fc (field capacity) set to w at which hydraulic conductivity equals
  ! a nominal drainage rate "rate_fc"
  ! w_wilt set to w at which psi is psi_wilt
  if (use_lm2_awc) then
     lake%w_wilt(:) = 0.15
     lake%w_fc  (:) = 0.15 + lake%pars%awc_lm2
  else
     lake%w_wilt(:) = lake%pars%w_sat &
          *(lake%pars%psi_sat_ref/(psi_wilt*lake%pars%alpha))**(1/lake%pars%chb)
     lake%w_fc  (:) = lake%pars%w_sat &
          *(rate_fc/(lake%pars%k_sat_ref*lake%pars%alpha**2))**(1/(3+2*lake%pars%chb))
  endif

  ! below made use of phi_e from parlange via entekhabi
  lake%Eg_part_ref = (-4*lake%w_fc(1)**2*lake%pars%k_sat_ref*lake%pars%psi_sat_ref*lake%pars%chb &
       /(PI*lake%pars%w_sat)) * (lake%w_fc(1)/lake%pars%w_sat)**(2+lake%pars%chb)   &
       *(2*PI/(3*lake%pars%chb**2*(1+3/lake%pars%chb)*(1+4/lake%pars%chb)))/2

  lake%z0_scalar = lake%pars%z0_momentum * exp(-k_over_B)
  lake%z0_scalar_ice = lake%pars%z0_momentum_ice * exp(-k_over_B_ice)
  lake%geothermal_heat_flux = geothermal_heat_flux_constant

end subroutine 

subroutine init_lake_data_0d(lake)
  type(lake_tile_type), intent(inout) :: lake

!  real rsa_exp         ! riparian source-area exponent

  integer :: k
  k = lake%tag

  lake%pars%w_sat             = dat_w_sat            (k)
  lake%pars%awc_lm2           = dat_awc_lm2          (k)
  lake%pars%k_sat_ref         = dat_k_sat_ref        (k)
  lake%pars%psi_sat_ref       = dat_psi_sat_ref      (k)
  lake%pars%chb               = dat_chb              (k)
  lake%pars%alpha             = 1
  lake%pars%heat_capacity_ref = dat_heat_capacity_ref(k)
  lake%pars%thermal_cond_ref  = dat_thermal_cond_ref (k)
  lake%pars%refl_dry_dir      = dat_refl_dry_dir     (k,:)
  lake%pars%refl_dry_dif      = dat_refl_dry_dif     (k,:)
  lake%pars%refl_sat_dir      = dat_refl_sat_dir     (k,:)
  lake%pars%refl_sat_dif      = dat_refl_sat_dif     (k,:)
  lake%pars%emis_dry          = dat_emis_dry         (k)
  lake%pars%emis_sat          = dat_emis_sat         (k)
  lake%pars%z0_momentum       = dat_z0_momentum      (k)
  lake%pars%z0_momentum_ice   = dat_z0_momentum_ice  (k)

  !!!
!  print*,'lake parameters'
!  print*,'connected_to_next: ',lake%pars%connected_to_next
!  print*,'whole_area: ',lake%pars%whole_area
!  print*,'depth_sill: ',lake%pars%depth_sill
!  print*,'backwater: ',lake%pars%backwater
!  print*,'backwater_1: ',lake%pars%backwater_1
!  print*,'width_sill: ',lake%pars%width_sill
!  print*,'w_sat o:',dat_w_sat(k),' n:',lake%pars%w_sat
!  print*,'awc_lm2 o:',dat_awc_lm2(k),' n:',lake%pars%awc_lm2
!  print*,'k_sat_ref o:',dat_k_sat_ref(k),' n:',lake%pars%k_sat_ref
!  print*,'psi_sat_ref o:',dat_psi_sat_ref(k),' n:',lake%pars%psi_sat_ref
!  print*,'chb o:',dat_chb(k),' n:',lake%pars%chb
!  print*,'alpha o:',1.0,' n:',lake%pars%alpha
!  print*,'heat_capacity_ref o:',dat_heat_capacity_ref(k),' n:',lake%pars%heat_capacity_ref
!  print*,'thermal_cond_ref o:',dat_thermal_cond_ref(k),' n:',lake%pars%thermal_cond_ref
!  print*,'refl_dry_dir o:',dat_refl_dry_dir(k,:),' n:',lake%pars%refl_dry_dir
!  print*,'refl_dry_dif o:',dat_refl_dry_dif(k,:),' n:',lake%pars%refl_dry_dif
!  print*,'refl_sat_dir o:',dat_refl_sat_dir(k,:),' n:',lake%pars%refl_sat_dir
!  print*,'refl_sat_dif o:',dat_refl_sat_dif(k,:),' n:',lake%pars%refl_sat_dif
!  print*,'emis_dry o:',dat_emis_dry(k),' n:',lake%pars%emis_dry 
!  print*,'emis_sat o:',dat_emis_sat(k),' n:',lake%pars%emis_sat
!  print*,'z0_momentum o:',dat_z0_momentum(k),' n:',lake%pars%z0_momentum
!  print*,'z0_momentum_ice o:',dat_z0_momentum_ice(k),'n:',lake%pars%z0_momentum_ice
  !!!

  lake%pars%rsa_exp           = rsa_exp_global

  ! -------- derived constant lake parameters --------
  ! w_fc (field capacity) set to w at which hydraulic conductivity equals
  ! a nominal drainage rate "rate_fc"
  ! w_wilt set to w at which psi is psi_wilt
  if (use_lm2_awc) then
     lake%w_wilt(:) = 0.15
     lake%w_fc  (:) = 0.15 + lake%pars%awc_lm2
  else
     lake%w_wilt(:) = lake%pars%w_sat &
          *(lake%pars%psi_sat_ref/(psi_wilt*lake%pars%alpha))**(1/lake%pars%chb)
     lake%w_fc  (:) = lake%pars%w_sat &
          *(rate_fc/(lake%pars%k_sat_ref*lake%pars%alpha**2))**(1/(3+2*lake%pars%chb))
  endif

  ! below made use of phi_e from parlange via entekhabi
  lake%Eg_part_ref = (-4*lake%w_fc(1)**2*lake%pars%k_sat_ref*lake%pars%psi_sat_ref*lake%pars%chb &
       /(PI*lake%pars%w_sat)) * (lake%w_fc(1)/lake%pars%w_sat)**(2+lake%pars%chb)   &
       *(2*PI/(3*lake%pars%chb**2*(1+3/lake%pars%chb)*(1+4/lake%pars%chb)))/2

  lake%z0_scalar = lake%pars%z0_momentum * exp(-k_over_B)
  lake%z0_scalar_ice = lake%pars%z0_momentum_ice * exp(-k_over_B_ice)
  lake%geothermal_heat_flux = geothermal_heat_flux_constant

end subroutine init_lake_data_0d


! ============================================================================
function lake_cover_cold_start(land_mask, lonb, latb, domain) result (lake_frac)
! creates and initializes a field of fractional lake coverage
  logical, intent(in) :: land_mask(:)    ! land mask
  real,    intent(in) :: lonb(:,:), latb(:,:)! boundaries of the grid cells
  real,    pointer    :: lake_frac (:,:) ! output: map of lake fractional coverage
  type(domain2d), intent(in) :: domain
  real :: lake_frac_sg(lnd_sg%is:lnd_sg%ie,lnd_sg%js:lnd_sg%je)
  allocate( lake_frac(size(land_mask(:)),n_dim_lake_types))

  if (trim(lake_to_use)=='from-rivers') then
     lake_frac = 0.0
     if (file_exist('INPUT/river_data.nc', domain)) then
         call read_data('INPUT/river_data.nc', 'lake_frac', lake_frac_sg, &
                        domain=domain)
         call mpp_pass_sg_to_ug(lnd%domain, lake_frac_sg, lake_frac(:,1))
     endif
     ! make sure 'missing values' don't get into the result
     where (lake_frac < 0) lake_frac = 0
     where (lake_frac > 1) lake_frac = 1
  else
     call init_cover_field(lake_to_use, 'INPUT/ground_type.nc', 'cover','frac', &
          lonb, latb, lake_index_constant, input_cover_types, lake_frac)
  endif

  if (round_frac_down) then
      where (lake_frac.gt.0. .and. lake_frac.lt.min_lake_frac) lake_frac = 0.
    else
      where (lake_frac.gt.0. .and. lake_frac.lt.min_lake_frac) lake_frac = min_lake_frac
    endif

end function lake_cover_cold_start

! =============================================================================
function lake_tiles_can_be_merged(lake1,lake2) result(response)
  logical :: response
  type(lake_tile_type), intent(in) :: lake1,lake2

  response = (lake1%tag==lake2%tag)
end function lake_tiles_can_be_merged

! =============================================================================
! combine two lake tiles with specified weights; the results goes into the
! second one
! THIS NEEDS TO BE REVISED FOR TILE-DEPENDENT DZ
subroutine merge_lake_tiles(t1,w1,t2,w2)
  type(lake_tile_type), intent(in)    :: t1
  type(lake_tile_type), intent(inout) :: t2
  real, intent(in) :: w1, w2 ! relative weights

  ! ---- local vars
  real    :: x1, x2 ! normalized relative weights
  real    :: HEAT1, HEAT2 ! temporaries for heat
  real    :: C1, C2 ! heat capacities
  real    :: gw
  integer :: i

WRITE (*,*) 'SORRY, BUT merge_lake_tiles NEEDS TO BE REVISED TO ALLOW FOR ', &
            'HORIZONTALLY VARYING VERTICAL DISCRETIZATION'

  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1.0 - x1

  ! combine state variables
  do i = 1,num_l
    ! calculate "dry" heat capacities:
    C1 = sfc_heat_factor*t1%pars%heat_capacity_ref
    C2 = sfc_heat_factor*t2%pars%heat_capacity_ref
    ! calculate heat content at this level for both source tiles
    HEAT1 = &
    0.
!    (C1*dz(i)+clw*t1%wl(i)+csw*t1%ws(i)) * (t1%T(i)-tfreeze)
    HEAT2 = &
    0.
!    (C2*dz(i)+clw*t2%wl(i)+csw*t2%ws(i)) * (t2%T(i)-tfreeze)
    ! calculate (and assign) combined water mass
    t2%wl(i) = t1%wl(i)*x1 + t2%wl(i)*x2
    t2%ws(i) = t1%ws(i)*x1 + t2%ws(i)*x2
    ! if dry heat capacity of combined lake is to be changed, update it here
    ! ...
    ! calculate combined temperature, based on total heat content and combined
    ! heat capacity
!    t2%T(i) = (HEAT1*x1+HEAT2*x2) / &
!      (C2*dz(i)+clw*t2%wl(i)+csw*t2%ws(i)) + tfreeze
    t2%T(i) = 0.
  enddo

end subroutine merge_lake_tiles

! =============================================================================
! returns true if tile fits the specified selector
function lake_is_selected(lake, sel)
  logical lake_is_selected
  type(tile_selector_type),  intent(in) :: sel
  type(lake_tile_type),      intent(in) :: lake

  lake_is_selected = (sel%idata1 == lake%tag)
end function lake_is_selected


! ============================================================================
! returns tag of the tile
function get_lake_tile_tag(lake) result(tag)
  integer :: tag
  type(lake_tile_type), intent(in) :: lake

  tag = lake%tag
end function get_lake_tile_tag


! ============================================================================
! compute bare-lake albedos and bare-lake emissivity
subroutine lake_radiation ( lake, cosz, &
     lake_refl_dir, lake_refl_dif, lake_refl_lw, lake_emis )
  type(lake_tile_type), intent(in) :: lake
  real, intent(in) :: cosz
  real, intent(out) :: lake_refl_dir(NBANDS), lake_refl_dif(NBANDS), lake_refl_lw, lake_emis

  ! ---- local vars
  real :: lake_sfc_vlc, blend
  real :: liq_value_dir(NBANDS), ice_value_dir(NBANDS)
  real :: liq_value_dif(NBANDS), ice_value_dif(NBANDS)
  real :: zenith_angle, zsq, zcu

  ! ---- radiation properties
  lake_sfc_vlc = lake%wl(1)/lake%dz(1)
  blend        = lake_sfc_vlc/lake%pars%w_sat
  if (use_brdf) then
     zenith_angle = acos(cosz)
     zsq = zenith_angle*zenith_angle
     zcu = zenith_angle*zsq
     liq_value_dir =  f_iso_liq*(g0_iso+g1_iso*zsq+g2_iso*zcu) &
                    + f_vol_liq*(g0_vol+g1_vol*zsq+g2_vol*zcu) &
                    + f_geo_liq*(g0_geo+g1_geo*zsq+g2_geo*zcu)
     ice_value_dir =  f_iso_ice*(g0_iso+g1_iso*zsq+g2_iso*zcu) &
                    + f_vol_ice*(g0_vol+g1_vol*zsq+g2_vol*zcu) &
                    + f_geo_ice*(g0_geo+g1_geo*zsq+g2_geo*zcu)
     liq_value_dif = g_iso*f_iso_liq + g_vol*f_vol_liq + g_geo*f_geo_liq
     ice_value_dif = g_iso*f_iso_ice + g_vol*f_vol_ice + g_geo*f_geo_ice
  else
     liq_value_dir = lake%pars%refl_sat_dir
     ice_value_dir = lake%pars%refl_dry_dir
     liq_value_dif = lake%pars%refl_sat_dif
     ice_value_dif = lake%pars%refl_dry_dif
  endif
  lake_refl_dir = ice_value_dir + blend*(liq_value_dir-ice_value_dir)
  lake_refl_dif = ice_value_dif + blend*(liq_value_dif-ice_value_dif)
  lake_emis     = lake%pars%emis_dry   + blend*(lake%pars%emis_sat-lake%pars%emis_dry  )
  lake_refl_lw  = 1 - lake_emis
end subroutine lake_radiation

! ============================================================================
! compute bare-lake roughness
subroutine lake_roughness ( lake,lake_z0s, lake_z0m )
  type(lake_tile_type), intent(in)  :: lake
  real,                 intent(out) :: lake_z0s, lake_z0m

  if (lake%ws(1).le.0.) then
      lake_z0s = lake%z0_scalar
      lake_z0m = lake%pars%z0_momentum
    else
      lake_z0s = lake%z0_scalar_ice
      lake_z0m = lake%pars%z0_momentum_ice
    endif
end subroutine lake_roughness

! ============================================================================
! compute lake thermodynamic properties.
subroutine lake_data_thermodynamics ( lake_pars, lake_depth, &
     lake_rh, heat_capacity_dry, thermal_cond)
  type(lake_pars_type), intent(in)  :: lake_pars
  real,                 intent(in)  :: lake_depth
  real,                 intent(out) :: lake_rh
  real,                 intent(out) :: heat_capacity_dry(:)
  real,                 intent(out) :: thermal_cond(:)

  ! ---- local vars
  integer l
  real lake_depth_frac, lake_rh_base

! ----------------------------------------------------------------------------

  lake_depth_frac = (lake_depth/lake_pars%depth_sill)
  lake_rh_base = (lake_depth_frac-dry_lake_depth_frac)/(1.-dry_lake_depth_frac)
  if (lake_rh_base.gt.0.) then
      lake_rh = min(max_lake_rh, lake_rh_base**lake_rh_exp )
    else
      lake_rh = 0.
    endif

  do l = 1, num_l
     heat_capacity_dry(l) = lake_pars%heat_capacity_ref
     thermal_cond(l)  = lake_pars%thermal_cond_ref
  enddo

end subroutine lake_data_thermodynamics

! ============================================================================
subroutine lake_tile_stock_pe (lake, twd_liq, twd_sol  )
  type(lake_tile_type),  intent(in)    :: lake
  real,                  intent(out)   :: twd_liq, twd_sol
  integer n

  twd_liq = 0.
  twd_sol = 0.
  do n=1, size(lake%wl)
    twd_liq = twd_liq + lake%wl(n)
    twd_sol = twd_sol + lake%ws(n)
    enddo

end subroutine lake_tile_stock_pe


! ============================================================================
! returns lake tile heat content, J/m2
function lake_tile_heat (lake) result(heat) ; real heat
  type(lake_tile_type),  intent(in)  :: lake

  integer :: i

  heat = 0
  do i = 1, num_l
     heat = heat + &
          (lake%heat_capacity_dry(i)*lake%dz(i) + clw*lake%wl(i) &
	     + csw*lake%ws(i))*(lake%T(i)-tfreeze) + &
          hlf*lake%ws(i)
  enddo
end function lake_tile_heat

end module lake_tile_mod
