! ============================================================================
! hillslope model module
! ============================================================================
module hillslope_mod

#include "../shared/debug.inc"

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use mpp_mod, only: mpp_pe, mpp_root_pe
use fms_mod, only: error_mesg, file_exist, close_file, check_nml_error, &
     stdlog, FATAL, NOTE, WARNING

use fms_io_mod, only: restart_file_type, free_restart_type, &
      set_domain, nullify_domain
use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_enum_type, &
     first_elmt, loop_over_tiles, get_elmt_indices
use land_utils_mod, only : put_to_tiles_r0d_fptr
use land_tile_diag_mod, only : diag_buff_type, &
     register_tiled_static_field, set_default_diag_filter, &
     send_tile_data_r0d_fptr, &
     send_tile_data_i0d_fptr, OP_SUM
use land_data_mod, only : lnd_sg, log_version, lnd
use land_io_mod, only : read_field
use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_restart_axis, add_int_tile_data, get_int_tile_data, &
     print_netcdf_error
use nf_utils_mod,  only : nfu_inq_dim
use land_debug_mod, only : is_watch_point, is_watch_cell, set_current_point
use land_transitions_mod, only : do_landuse_change
use vegn_harvesting_mod , only : do_harvesting
use hillslope_tile_mod , only : register_hlsp_selectors
use constants_mod, only : tfreeze
use soil_tile_mod, only : gw_option, GW_TILED, initval, soil_tile_type, &
     gw_scale_length, gw_scale_relief

implicit none
private

! ==== public method interfaces =====================================================
public :: read_hlsp_namelist
public :: hlsp_init      ! read surface parameters, read restart file, set
                         ! hillslope-position-dependent parameters
public :: hlsp_init_predefined ! Initialize hillslope using predefined tile
                               ! parameters
!public :: get_max_hidx    ! evaluate maximum hillslope indices for the gridcell
public :: hlsp_coldfracs ! determine # and fractions of tiles within hillslopes for cold start
                          ! and determine hillslope indices
public :: retrieve_hlsp_indices ! returns hillslope position and parent to land_cover_cold_start.
                                ! All indices are 0 if .not. do_hillslope_model.
public :: horiz_wt_depth_to_init ! returns water table depth to init, when soil:horiz_init_wt == .true.
public :: hlsp_end
public :: save_hlsp_restart
public :: hlsp_config_check     ! Check configuration for errors, at the end of land_model_init sequence.
                                ! Also deallocate any module variables used during cold start.
public :: calculate_wt_init  ! Calculates water table depth for initialization to be returned by
                             ! horiz_wt_depth_to_init.

! =====end of public interfaces ==============================================

! =====private methods
private :: read_hillslope_surfdat ! shared function used by hlsp_init and hlsp_cold_fracs
private :: hlsp_diag_init
private :: meanelev ! used by calculate_wt_init


! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'hillslope'
#include "../shared/version_variable.inc"

integer, parameter :: max_vc = 30 ! Max num_vertclusters that can be input from namelist for
                                  ! tile horizontal grid.
! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
logical, protected, public :: do_hillslope_model = .false. ! Activates model
logical             :: fixed_num_vertclusters    = .true.  ! False ==> adaptive hillslope "grid"
                                                           ! Hardwired for now
integer, protected, public :: num_vertclusters   = 5       ! Number of vertical clusters (initial, if
                                                           ! fixed_num_vertclusters == .false.)
integer, protected, public :: max_num_topo_hlsps = 3       ! Global maximum # allowed actively modeled
                                                           ! hillslopes in each gridcell.
logical             :: hillslope_horz_subdiv     = .true.  ! Allow tiles in the same hillslope position
                                                           ! to be subdivided for landuse change, etc.
logical             :: hillslope_topo_subdiv     = .false. ! Allow multiple instances of each
                                                           ! topographic hillslope type based on land
                                                           ! use history, etc.
                                                           ! Hardwired
logical             :: soil_types_by_hlsp        = .false. ! Assign soil types to each hillslope
                                                           ! based on hillslope input datafile.
! ZMS Could add option to have equal area rather than equal length tiles.
real, protected, public :: strm_depth_penetration = 1.5    ! (m) depth to which streams communicate
                                                           ! with lowland tiles
logical, protected, public :: use_hlsp_aspect_in_gwflow = .false. ! True ==> head gradients effectively
                       ! reduced by sin(hillslope angle).
logical, protected, public :: use_geohydrodata   = .true.  ! True ==> input some gridcell-level data
                       ! from the "geohydrology.nc" dataset to initialize hillslope soil properties,
                       ! and keep some of the diagnostic hydrology calculations based on the gridcell
                       ! values.
!real, protected, public   :: pond               = 0.      ! [mm] water / ice allowed to pool on
!                       ! surface before run-off
logical, protected, public :: stiff_do_explicit = .false.  ! update water profile explicitly
                       ! due to inter-tile flows in case becomes stiff during timestep.
                       ! Hard-wired for now.
logical             :: diagnostics_by_cluster    = .false.  ! True ==> create diagnostics for
                                                            ! every tile cluster from 1:num_vertclusters
logical             :: init_wt_strmelev      = .true.    ! initialize water table at stream depth
                                                         ! else use requested soil:init_wtdep
logical             :: equal_length_tiles    = .true.    ! .false. ==> tile lengths input from namelist
real, dimension(max_vc) :: dl = 0.0  ! [-] vector of length fractions input from namelist for tile horizontal grid
                               ! Ordered from stream-most tile to hilltop.
logical, protected, public :: dammed_strm_bc = .false.  ! True ==> Hydraulic head increases with stream
                       !depth, limiting outflow with large strm_depth_penetration
logical, protected, public :: simple_inundation = .false. ! True ==> simple inundation scheme based on
                       ! microtopography
logical, protected, public :: exp_inundation = .false. ! Use inun_frac based on exp{-zwt_b/microtopo} even
                                                ! if .not. simple_inundation
real, protected, public    :: surf_flow_velocity = 1.  ! [m/s] Assumed nominal surface runoff velocity at a slope
                    ! (tangent) of 1. (for simple_inundation)
logical, protected, public :: limit_intertile_flow = .false. ! True ==> Limit explicit inter-tile flows
                       ! to improve numerical stability
real, protected, public    :: flow_ratio_limit = 1.    ! max delta psi to length ratio allowed, if limit_intertile_flow
logical, protected, public :: tiled_DOC_flux = .false. ! True ==> Calculate DOC fluxes for soil carbon model

character(len=256)  :: hillslope_surfdata = 'INPUT/hillslope.nc'
character(len=24)   :: hlsp_interpmethod = 'nearest'
character(*), parameter  :: hlsp_rst_ifname = 'INPUT/hlsp.res.nc'
character(*), parameter  :: hlsp_rst_ofname = 'hlsp.res.nc'

! Removed do_hillslope_model from namelist.  Tie to gw_option in soil_tile.
namelist /hlsp_nml/ num_vertclusters, max_num_topo_hlsps, hillslope_horz_subdiv, &
                    hillslope_surfdata, hlsp_interpmethod, soil_types_by_hlsp, &
                    strm_depth_penetration, use_hlsp_aspect_in_gwflow, use_geohydrodata, &
                    diagnostics_by_cluster, init_wt_strmelev, dammed_strm_bc, &
                    simple_inundation, surf_flow_velocity, dl, equal_length_tiles, &
                    limit_intertile_flow, flow_ratio_limit, exp_inundation, tiled_DOC_flux
! hardwired: fixed_num_vertclusters, hillslope_topo_subdiv, stiff_do_explicit
!---- end of namelist --------------------------------------------------------

logical             :: module_is_initialized = .false.
character(len=24)   :: hlsp_surf_dimname = 'nhlsps'

! ---- diagnostic field IDs
integer, protected, public :: &!id_soil_e_depth,
                    id_microtopo, id_tile_hlsp_length, id_tile_hlsp_slope, &
                    id_tile_hlsp_elev, id_tile_hlsp_hpos, id_tile_hlsp_width, & !id_transm_bedrock, &
                    id_hidx_j, id_hidx_k

! vars to save between subroutine calls during cold start
logical  :: cold_start = .false.  ! flag for cold start (==true) or initial conditions (==false)
! dims i,j, tile
integer, private, allocatable :: hidx_j_loc (:,:) ! hillslope position index
integer, private, allocatable :: hidx_k_loc (:,:) ! hillslope parent index

! dims i,j, hk, hj
! ZMS For now make available only during cold start. Perhaps later make available until model end.
real, private, allocatable :: fracs_loc(:,:,:) ! fractions by topo bin
real, private, allocatable :: elev_loc(:,:,:) ! elev by topo bin [m]
real, private, allocatable :: init_wt_loc (:,:,:) ! water table depth to init by hillslope cluster


! ==== end of module variables ===============================================

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

contains

! ============================================================================
subroutine read_hlsp_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call log_version(version, module_name, &
  __FILE__)
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=hlsp_nml, iostat=io)
  ierr = check_nml_error(io, 'hlsp_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=hlsp_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'hlsp_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=hlsp_nml)
  endif

  ! Set overall model switch based on soil_tile_mod:gw_option.  This requires that call to
  ! read_soil_namelist occurs before call to read_hlsp_namelist.
  if (gw_option == GW_TILED) then
     do_hillslope_model = .true.
  else
     do_hillslope_model = .false.
  end if

  if (mpp_pe() == mpp_root_pe()) then
     if (do_hillslope_model) then
        call error_mesg(module_name, 'Tiled hillslope model activated based on namelist option '// &
                        'gw_option == "tiled".', NOTE)
     else
        call error_mesg(module_name, 'Tiled hillslope model not activated based on namelist option'// &
                        'gw_option /= "tiled".', NOTE)
     end if
  end if

  if (do_hillslope_model .and. (.not. equal_length_tiles) .and. num_vertclusters > max_vc) &
     call error_mesg(module_name, 'Error: num_vertclusters requested is larger than max_vc, the '// &
                     'maximum size of the input tile length fraction vector.', FATAL)

  if (do_hillslope_model)  call register_hlsp_selectors(max_num_topo_hlsps, num_vertclusters, &
                                                          diagnostics_by_cluster)

end subroutine read_hlsp_namelist


! Read surface data for model initialization. Called by hlsp_init and hlsp_coldfracs.
! ============================================================================
subroutine read_hillslope_surfdat ( ls, le, num_topo_hlsps, frac_topo_hlsps, soil_e_depth, &
                                    microtopo, hlsp_length, hlsp_slope, hlsp_slope_exp, &
                                    hlsp_top_width, k_sat_gw, soiltype)
                                    ! Removed hlsp_stream_width

  integer, intent(in)   :: ls, le  ! unstructure domain index bounds
  integer, intent(out)  :: num_topo_hlsps(ls:le)  ! number of modeled hillslope topo types
  real, intent(out)     :: frac_topo_hlsps(ls:le, max_num_topo_hlsps)
                           ! prevalence of modeled hillslope topo types
  real, intent(out)     :: soil_e_depth(ls:le, max_num_topo_hlsps), &
                           microtopo(ls:le, max_num_topo_hlsps), &
                           hlsp_length(ls:le, max_num_topo_hlsps), &
                           hlsp_slope(ls:le, max_num_topo_hlsps), &
                           hlsp_slope_exp(ls:le, max_num_topo_hlsps), & !hlsp_stream_width(ls:le, max_num_topo_hlsps), &
                           hlsp_top_width(ls:le, max_num_topo_hlsps), &
                           k_sat_gw(ls:le, max_num_topo_hlsps)
  integer, optional, intent(out) :: soiltype(ls:le, max_num_topo_hlsps)
  integer, allocatable :: ibuffer(:,:)
  real   , allocatable :: rbuffer(:,:)
  integer :: ierr, ncid, totnumhlsps ! err, file id and tot # of hillslopes per gcell on surfdata

  ! Check length of nhlsps dimension on input file.
  __NF_ASRT__(nf_open(hillslope_surfdata, NF_NOWRITE,ncid))
  ierr = nfu_inq_dim(ncid, hlsp_surf_dimname, totnumhlsps)
  __NF_ASRT__(nf_close(ncid))
  !write(*,*)'totnumhlsps = ',totnumhlsps,', max_num_topo_hlsps = ', &
  !      max_num_topo_hlsps
  if (ierr > 0) call error_mesg(module_name, 'Error reading file '// hillslope_surfdata // ','// &
                                'dimension '// hlsp_surf_dimname // '.', FATAL)

  if (totnumhlsps < max_num_topo_hlsps) &
     call error_mesg(module_name, 'Number of active modeled hillslopes per gridcell requested '// &
                     'in namelist exceeds that available in hillslope inputdata.', &
                     FATAL)

  allocate(ibuffer(ls:le, totnumhlsps), rbuffer(ls:le, totnumhlsps))

  if (mpp_pe() == mpp_root_pe()) &
     call error_mesg(module_name, 'read_hillslope_surfdat: WARNING: Using nearest interpolation '// &
                     'which currently may not be robust for cubic sphere grid if hillslope dataset '// &
                     'is not on the native grid.', NOTE)
  ! Note: this function is not currently a robust "nearest" interpolation for cubic-sphere
  ! grids and will need to be updated.

  call read_field( hillslope_surfdata, 'NUM_TOPO_HLSPS', lnd%lon, lnd%lat, num_topo_hlsps, &
                   interp='nearest' )

  call read_field( hillslope_surfdata, 'FRAC_TOPO_HLSPS', lnd%lon, lnd%lat, rbuffer, &
                   interp='nearest' )
  frac_topo_hlsps(:,:) = rbuffer(:,1:max_num_topo_hlsps)

  if (.not. use_geohydrodata) then
     call read_field( hillslope_surfdata, 'SOIL_E_DEPTH', lnd%lon, lnd%lat, rbuffer, &
                      interp=hlsp_interpmethod )
     soil_e_depth(:,:) = rbuffer(:,1:max_num_topo_hlsps)
  else
     soil_e_depth(:,:) = initval ! will not be used
  end if
  call read_field( hillslope_surfdata, 'MICROTOPO', lnd%lon, lnd%lat, rbuffer, &
                   interp=hlsp_interpmethod )
  microtopo(:,:) = rbuffer(:,1:max_num_topo_hlsps)

  call read_field( hillslope_surfdata, 'HLSP_LENGTH', lnd%lon, lnd%lat, rbuffer, &
                   interp=hlsp_interpmethod )
  hlsp_length(:,:) = rbuffer(:,1:max_num_topo_hlsps)
  if (use_geohydrodata) hlsp_length(:,:) = hlsp_length(:,:)*gw_scale_length

  call read_field( hillslope_surfdata, 'HLSP_SLOPE', lnd%lon, lnd%lat, rbuffer, &
                   interp=hlsp_interpmethod )
  ! Hillslope elevation at top divided by hillslope length
  hlsp_slope(:,:) = rbuffer(:,1:max_num_topo_hlsps)
  if (use_geohydrodata) hlsp_slope(:,:) = hlsp_slope(:,:)*gw_scale_relief

  call read_field( hillslope_surfdata, 'HLSP_SLOPE_EXP', lnd%lon, lnd%lat, rbuffer, &
                   interp=hlsp_interpmethod )
  ! Hillslope profile will follow equation z = H(x/L)^a, where a=hlsp_slope_exp,
  ! H = max elevation, and L = max length from stream.
  hlsp_slope_exp(:,:) = rbuffer(:,1:max_num_topo_hlsps)

!  call read_field( hillslope_surfdata, 'HLSP_STREAM_WIDTH', lnd_sg%lon, lnd_sg%lat, rbuffer, &
!                   interp=hlsp_interpmethod )
!  hlsp_stream_width(:,:) = rbuffer(:,1:max_num_topo_hlsps)

  call read_field( hillslope_surfdata, 'HLSP_TOP_WIDTH', lnd%lon, lnd%lat, rbuffer, &
                   interp=hlsp_interpmethod )
  hlsp_top_width(:,:) = rbuffer(:,1:max_num_topo_hlsps)

  if (.not. use_geohydrodata) then
     call read_field( hillslope_surfdata, 'BEDROCK_KSAT', lnd%lon, lnd%lat, rbuffer, &
                      interp=hlsp_interpmethod )
     k_sat_gw(:,:) = rbuffer(:,1:max_num_topo_hlsps)
  else
     k_sat_gw(:,:) = initval ! will not be used
  end if

  if (present(soiltype)) then
     call read_field( hillslope_surfdata, 'SOILTYPE', lnd%lon, lnd%lat, ibuffer, &
                     interp='nearest')
     soiltype(:,:) = ibuffer(:,1:max_num_topo_hlsps)
  end if

  deallocate(rbuffer, ibuffer)

end subroutine read_hillslope_surfdat


! ============================================================================
! Calculate fractional areas of tiles in each gridcell for cold start, and determine
! hillslope indices. Called from soil_cover_cold_start, only if do_hillslope_model.
subroutine hlsp_coldfracs ( tile_frac, n_dim_soil_types )
  real, pointer  :: tile_frac (:,:)
  integer, intent(in) :: n_dim_soil_types
  integer, allocatable :: num_topo_hlsps(:)
  real, allocatable :: frac_topo_hlsps(:,:)
  real, allocatable, dimension(:,:) :: soil_e_depth, microtopo, hlsp_length, &
                                       hlsp_slope, hlsp_slope_exp, &
                                       hlsp_top_width, k_sat_gw ! Removed hlsp_stream_width
  integer, allocatable, dimension(:,:) :: soiltype
  integer :: num_lnd, num_tiles, lls, lle, l, t, h, n, &
             ll, st ! index ranges and bounds
  real, allocatable :: soilfracs0(:) ! temporary
  real :: sumfracs, hillfrac, convergence, cfact
  real :: smallnumber = 1.e-12  ! Should be several orders of magnitude greater than machine
                                ! precision for error checking.
  real :: NN ! num_vertclusters, real
  real :: sum_l ! sum of dl
  character(len=256) :: mesg
  NN = real(num_vertclusters)
  num_lnd = size(tile_frac,1)
  num_tiles = size(tile_frac,2)
  lls = lnd%ls
  lle = lnd%le
  if ( (num_lnd /= lle - lls + 1) ) &
     call error_mesg(module_name, 'Land mask used in soil_cover_cold_start and passed to hlsp_coldfracs'// &
                     'is NOT equal in size to local unstructured grid bounds.', FATAL)

  if (.not. equal_length_tiles) then ! Error checking on input dl vector of tile length fracs
     sum_l = sum(dl(1:num_vertclusters))
     if (abs(sum_l - 1.) > smallnumber) then
        write(mesg,*)'Error: sum(dl) /= 1. sum(dl) = ', sum_l, &
                     '. Bad namelist dl in hlsp_nml. dl must match num_vertclusters.'
        call error_mesg(module_name, mesg, FATAL)
     end if
     if (.not. all(dl(1:num_vertclusters) > 0.) .or. .not. all(dl(1:num_vertclusters) < 1.)) &
        call error_mesg(module_name, 'Error: bad input dl in hlsp_nml. dl must match num_vertclusters.', &
                        FATAL)
  end if

  allocate(num_topo_hlsps(lls:lle), frac_topo_hlsps(lls:lle, max_num_topo_hlsps), &
            soil_e_depth(lls:lle, max_num_topo_hlsps), &
            microtopo(lls:lle, max_num_topo_hlsps), &
            hlsp_length(lls:lle, max_num_topo_hlsps), &
            hlsp_slope(lls:lle, max_num_topo_hlsps), &
            hlsp_slope_exp(lls:lle, max_num_topo_hlsps), &
            hlsp_top_width(lls:lle, max_num_topo_hlsps), &
            k_sat_gw(lls:lle, max_num_topo_hlsps) )
  if (soil_types_by_hlsp) &
     allocate ( soiltype(lls:lle, max_num_topo_hlsps) )

  allocate(hidx_j_loc(num_lnd, num_tiles), hidx_k_loc(num_lnd, num_tiles))

  ! For use in other subroutines
  cold_start = .true. ! if this subroutine is called then model is cold starting
  allocate(fracs_loc(lls:lle, max_num_topo_hlsps, num_vertclusters), &
           elev_loc(lls:lle, max_num_topo_hlsps, num_vertclusters))
  fracs_loc(:,:,:) = initval
  elev_loc(:,:,:) = initval

  ! Retrieve surface data
  if (soil_types_by_hlsp) then
     call read_hillslope_surfdat ( lls, lle, num_topo_hlsps, frac_topo_hlsps, &
                                   soil_e_depth, microtopo, hlsp_length, hlsp_slope, &
                                   hlsp_slope_exp, & !hlsp_stream_width, &
                                   hlsp_top_width, k_sat_gw, soiltype=soiltype)
  else
     call read_hillslope_surfdat ( lls, lle, num_topo_hlsps, frac_topo_hlsps, &
                                   soil_e_depth, microtopo, hlsp_length, hlsp_slope, &
                                   hlsp_slope_exp, & !hlsp_stream_width, &
                                   hlsp_top_width, k_sat_gw)
  end if

  ! soil_cover_cold_start now configured so that tile_frac must be sent back with
  ! clumps of repeating sets of topographically adjacent tiles. (So if there are multiple
  ! soil types, these will be repeated in each topographic bin. See how "soiltags" is set
  ! in soil_tile_mod.)

  allocate(soilfracs0(n_dim_soil_types))

  ! Walk through land model grid and reassign fractions. Distribute the n_dim_soil_types
  ! across the up to max_num_topo_hlsps * num_vertclusters topo bins.

  do l = 1, num_lnd
     ll = lls + l - 1

     soilfracs0(:) = tile_frac(l, 1:n_dim_soil_types) ! Copy existing fractions

     num_topo_hlsps(ll) = min(num_topo_hlsps(ll), max_num_topo_hlsps)
                              ! surface data could potentially allow for more than namelist
                              ! (but not less)

     sumfracs = sum(frac_topo_hlsps(ll, 1:num_topo_hlsps(ll)))
     ! Only the active num_topo_hlsps in each gcell will be used; relative hillslope
     ! fraction will be normalized to 1 (to conserve total tile fraction returned by
     ! soil_cover_cold_start).

     if (sumfracs > (1. + smallnumber) .or. sumfracs <= 0.) then
        write(mesg,*) 'hlsp_curfracs: Invalid hillslope fractions at land '// &
                      'point ll:', ll, '.'
        call error_mesg(module_name, mesg, FATAL)
     end if

     ! Allocate to tiles
     t = 1
     tile_frac(l, 1:num_tiles) = 0.
     do h = 1, max_num_topo_hlsps
        if (h <= num_topo_hlsps(ll)) then

           convergence = hlsp_top_width(ll, h) !/ 1. hlsp_stream_width(li, lj, h)
           if (convergence <= 0.) then
              write(mesg, *)'hlsp_coldfracs: land point l, hillslope h:', ll, h, '.', &
                            'hlsp_top_width = ', convergence, ' <= 0!'
              call error_mesg(module_name, mesg, FATAL)
           end if
           ! Converging hillslopes have value > 1; diverging < 1.

           sum_l = 0 ! running sum of dl
           do n = 1, num_vertclusters ! Order from lowest elevation to highest in each hillslope.
              ! Calculate fractions for topo bin
              ! Initial tiles will have assumed trapezoid shapes as projected onto horizontal
              ! surface.
              if (equal_length_tiles) then
                 hillfrac = frac_topo_hlsps(ll, h) / sumfracs &  ! hillslope overall area
                            * ( (2*NN-1 - (n-1)*2)/(2*NN) + (1+(n-1)*2)/(2*NN)*convergence ) & ! hillslope relative area
                            / real(num_vertclusters) / (1. + convergence) * 2. ! normalization
                 fracs_loc(ll, h, n) = hillfrac
              else ! length fraction assigned by input
                 sum_l = sum_l + dl(n)
                 !Convergence factor will be width at dl divided by average width
                 cfact = (convergence*(sum_l - dl(n)/2.) + 1. - (sum_l - dl(n)/2.)) &
                         *2./(1. + convergence)
                 hillfrac = frac_topo_hlsps(ll, h) / sumfracs & ! hillslope overall area
                            * dl(n) * cfact ! hillslope relative area
!!!!
                 fracs_loc(ll, h, n) = hillfrac
              end if
              ! distribute among soil types
              do st = 1, n_dim_soil_types
                 if (.not. soil_types_by_hlsp) then
                    tile_frac(l, t) = hillfrac * soilfracs0(st)
                 else if (soiltype(ll,h) == st) then
                    tile_frac(l, t) = hillfrac * sum(soilfracs0(:))
                 ! else weight remains zero
                 end if
                 ! Assign tile indices (which will be returned in retrieve_hillslope_indices)
                 ! Indices for zero-weight tiles will be ignored in land_cover_cold_start_0d
                 hidx_j_loc(l, t) = n
                 hidx_k_loc(l, t) = h
                 ! Advance tile index
                 t = t + 1

                 ! Debug
#ifdef ZMSDEBUG_INIT
                 !if (is_watch_point()) then
                    write(*,*)'Assigned tile fraction in hillslope: ll, hidx_j,', &
                              ' hidx_k, tile_frac:', ll, n, h, tile_frac(l, t-1)
                 !end if
#endif
              end do
           end do
        end if
     end do


     ! Error checking
     if (abs(sum(tile_frac(l, :)) - sum(soilfracs0(:))) > smallnumber) then
        write(mesg,*) 'hlsp_curfracs: Error in calculating tile fractions at '// &
                      'point ll:', ll, ', error:', &
                      sum(tile_frac(l, :)) - sum(soilfracs0(:))
        call error_mesg(module_name, mesg, FATAL)
     end if

  end do ! num_lnd

  deallocate(num_topo_hlsps, frac_topo_hlsps, soil_e_depth, &
             microtopo, hlsp_length, &
             hlsp_slope, hlsp_slope_exp, &!hlsp_stream_width, &
             hlsp_top_width, k_sat_gw, soilfracs0)

  if (soil_types_by_hlsp) deallocate(soiltype)
end subroutine hlsp_coldfracs



! ============================================================================
! return hillslope indices
subroutine retrieve_hlsp_indices (hlsp_pos, hlsp_par)
   integer, pointer :: hlsp_pos(:,:), hlsp_par(:,:)

   if (do_hillslope_model) then
      if (size(hlsp_pos,1) /= size(hidx_j_loc,1) .or. size(hlsp_pos,2) /= size(hidx_j_loc,2) ) &
         call error_mesg(module_name, &
               'Wrong dimension size in "hillslope_mod:retrieve_hillslope_indices.',FATAL)
      hlsp_pos(:,:) = hidx_j_loc(:,:)
      hlsp_par(:,:) = hidx_k_loc(:,:)
      deallocate(hidx_j_loc, hidx_k_loc)
   else
      hlsp_pos(:,:) = 0
      hlsp_par(:,:) = 0
   end if

end subroutine retrieve_hlsp_indices


! ============================================================================
! initialize hillslope model
subroutine hlsp_init(id_ug)
  integer,intent(in) :: id_ug !<Unstructured axis id.

  ! ---- local vars
  integer :: unit         ! unit for various i/o
  type(land_tile_enum_type)     :: ce  ! tail and current tile list elements
  type(land_tile_type), pointer :: tile   ! pointer to current tile
  real, allocatable, dimension(:,:) :: frac_topo_hlsps
  integer, allocatable, dimension(:) :: num_topo_hlsps
  real, allocatable, dimension(:,:) :: soil_e_depth, microtopo, hlsp_length, &
                                     hlsp_slope, hlsp_slope_exp, &
                                     hlsp_top_width, k_sat_gw ! Removed hlsp_stream_width
  integer :: lis, lie ! unstruct grid bounds
  integer :: li, lj, lk, ll ! lat, lon, tile indices
  integer :: hj, hk ! hillslope pos, par indices
  real    :: rhj, NN ! real hj, num_vertclusters
  real    :: c, a ! convergence, hlsp_slope
  character(len=256) :: mesg
  logical :: restart_exists
  real :: tfreeze_diff
  real :: hpos ! local horizontal position (-)
  real :: hbdu ! local upstream bdy (-)
  real :: hbdd ! local downstream bdy (-)
  type(land_restart_type) :: restart

  module_is_initialized = .TRUE.

  if (.not. do_hillslope_model) then
     ! Set hillslope indices all to 0 and return
     ce = first_elmt(land_tile_map)
     do while(loop_over_tiles(ce,tile))
        if (associated(tile%soil)) then
           tile%soil%hidx_j = 0
           tile%soil%hidx_k = 0
        end if
     end do
     return
  end if

  ! initialize hillslope-dependent diagnostic fields
  call hlsp_diag_init(id_ug)

  ! -------- initialize state --------
  lis = lnd%ls
  lie = lnd%le
  NN = real(num_vertclusters)

  allocate(num_topo_hlsps(lis:lie), frac_topo_hlsps(lis:lie, max_num_topo_hlsps), &
            soil_e_depth(lis:lie, max_num_topo_hlsps), &
            microtopo(lis:lie, max_num_topo_hlsps), &
            hlsp_length(lis:lie, max_num_topo_hlsps), &
            hlsp_slope(lis:lie, max_num_topo_hlsps), &
            hlsp_slope_exp(lis:lie, max_num_topo_hlsps), &
            hlsp_top_width(lis:lie, max_num_topo_hlsps), &
            k_sat_gw(lis:lie, max_num_topo_hlsps) )

  ! Retrieve surface data
  call read_hillslope_surfdat ( lis, lie, num_topo_hlsps, frac_topo_hlsps, &
                                soil_e_depth, microtopo, hlsp_length, hlsp_slope, hlsp_slope_exp, &
                                hlsp_top_width, k_sat_gw)

  call open_land_restart(restart,hlsp_rst_ifname,restart_exists)
  if (restart_exists) then
     call error_mesg(module_name, 'hlsp_init, '// &
          'reading NetCDF restart "'//trim(hlsp_rst_ifname)//'"', &
          NOTE)
     if (cold_start) &
        call error_mesg(module_name, 'hlsp_init: coldfracs subroutine called even though restart file '// &
                             'exists! Inconsistency of "cold_start" in hillslope_mod.', FATAL)
     call get_int_tile_data(restart, 'HIDX_J', soil_hidx_j_ptr)
     call get_int_tile_data(restart, 'HIDX_K', soil_hidx_k_ptr)
  else
     call error_mesg(module_name, 'hlsp_init: '// &
          'cold-starting hillslope model',&
          NOTE)
     ! These indices will have been set in land_cover_cold_start_0d
     if (.not. cold_start) &
        call error_mesg(module_name, 'hlsp_init: coldfracs subroutine not called even though restart file '// &
                             'does not exist! Inconsistency of "cold_start" in hillslope_mod.', FATAL)

  endif
  call free_land_restart(restart)

  tfreeze_diff = 0. ! Initialize before tile loop.

  ! Assign hillslope-topographic-position dependent parameters to tiles.
  ce = first_elmt(land_tile_map, ls=lis)
  do while(loop_over_tiles(ce,tile,ll,lk))
     call set_current_point(ll,lk)

     if(.not.associated(tile)) then
        call error_mesg(module_name, "tile%soil is not associated", FATAL)
     endif

     if (.not.associated(tile%soil)) cycle

     hj = tile%soil%hidx_j
     hk = tile%soil%hidx_k
     ! ZMS Note: To allow multiple instances of each topo hillslope, add check here to see if
     ! hk > max_num_topo_hillslopes.  If so, use mod(hk, max_num_topo_hillslopes), or
     ! max_num_topo_hillslopes if mod == 0.

     ! ZMS Note: If code to be used to allow different properties based on vertical position
     ! WITHIN hillslopes, it would be done here...
     if (.not. use_geohydrodata) then
        tile%soil%pars%soil_e_depth = soil_e_depth(ll,hk)
        tile%soil%pars%k_sat_gw = k_sat_gw(ll,hk)

        ! Set LM3.1 variables based on overall hillslope values:
        tile%soil%pars%hillslope_length = hlsp_length(ll,hk)
        tile%soil%pars%hillslope_relief = hlsp_length(ll,hk)*hlsp_slope(ll,hk)
        c = hlsp_top_width(ll,hk)
        a = hlsp_slope_exp(ll,hk)
        if (a > -1.) then ! else abort below
           tile%soil%pars%hillslope_zeta_bar = 2.*(1. + c + a*c)/(2.+3.*a+a**2.)/(1.+c) ! ZMS Check this
        end if
     end if
     tile%soil%pars%microtopo = microtopo(ll,hk)

     ! Check valid hlsp_slope_exp
     if (hlsp_slope_exp(ll,hk) <= 0.) then
        write(mesg,*)'hlsp_init: invalid hillslope slope exponent i,j, hillslope k = ',li,lj,hk, &
                     '. hlsp_slope_exp = ', hlsp_slope_exp(ll,hk), ' <= 0!'
        call error_mesg(module_name, mesg, FATAL)
     end if

     if (fixed_num_vertclusters .and. equal_length_tiles) then
        tile%soil%pars%tile_hlsp_length = hlsp_length(ll,hk) / NN
        rhj = real(hj)
        tile%soil%pars%tile_hlsp_slope  = hlsp_slope(ll,hk) * &
                ( (rhj/NN)**hlsp_slope_exp(ll,hk) - ( (rhj-1)/NN)**hlsp_slope_exp(ll,hk) ) / &
                ( rhj/NN - (rhj-1)/NN)
        tile%soil%pars%tile_hlsp_elev   = hlsp_slope(ll,hk) * hlsp_length(ll,hk) * &
                                          0.5*( (rhj/NN)**hlsp_slope_exp(ll,hk) + &
                                                ((rhj-1)/NN)**hlsp_slope_exp(ll,hk) )
                                      ! use mean elevation of endpoints to match tile_hlsp_slope
        tile%soil%pars%tile_hlsp_hpos   = hlsp_length(ll,hk) * (rhj-0.5)/NN
!        tile%soil%pars%tile_hlsp_width  = hlsp_stream_width(li,lj,hk) + (rhj-0.5)/NN * &
!                                          (hlsp_top_width(li,lj,hk) - hlsp_stream_width(li,lj,hk))
        tile%soil%pars%tile_hlsp_width  = 1. + (rhj-0.5)/NN * &
                                          (hlsp_top_width(ll,hk) - 1.)
     else if (fixed_num_vertclusters) then ! tile_hlsp_length set from input namelist dl
        tile%soil%pars%tile_hlsp_length = hlsp_length(ll,hk) * dl(hj)
        hpos = dl(hj)/2.
        if (hj > 1) hpos = hpos + sum(dl(1:hj-1))
        hbdd = hpos - dl(hj)/2.
        hbdu = hpos + dl(hj)/2.
        tile%soil%pars%tile_hlsp_slope  = hlsp_slope(ll,hk) * &
                 (hbdu**hlsp_slope_exp(ll,hk) - hbdd**hlsp_slope_exp(ll,hk)) / dl(hj)
        tile%soil%pars%tile_hlsp_elev   = hlsp_slope(ll,hk) * hlsp_length(ll,hk) * &
                                          0.5*( hbdu**hlsp_slope_exp(ll,hk) + &
                                                hbdd**hlsp_slope_exp(ll,hk) )
                                      ! use mean elevation of endpoints to match tile_hlsp_slope
        tile%soil%pars%tile_hlsp_hpos   = hlsp_length(ll,hk) * hpos
        tile%soil%pars%tile_hlsp_width  = 1. + hpos * &
                                          (hlsp_top_width(ll,hk) - 1.)

     ! else will need to save these variables in restart file
     end if
     ! Set elev_loc for later use
     if (.not. restart_exists) then
        elev_loc(ll,hk,hj) = tile%soil%pars%tile_hlsp_elev
     end if

     ! Debug
     if (is_watch_cell()) then
        write(*,*)'use_geohydrodata = ', use_geohydrodata
        write(*,*)'hlsp_init: li,lj,hj,hk: ',li,lj,hj,hk
        write(*,*)'soil_e_depth: ', tile%soil%pars%soil_e_depth
        write(*,*)'microtopo: ', tile%soil%pars%microtopo
        write(*,*)'k_sat_gw: ', tile%soil%pars%k_sat_gw
        write(*,*)'tile_hlsp_length: ', tile%soil%pars%tile_hlsp_length
        write(*,*)'tile_hlsp_slope: ', tile%soil%pars%tile_hlsp_slope
        write(*,*)'tile_hlsp_elev: ', tile%soil%pars%tile_hlsp_elev
        write(*,*)'tile_hlsp_hpos: ', tile%soil%pars%tile_hlsp_hpos
        write(*,*)'tile_hlsp_width: ', tile%soil%pars%tile_hlsp_width
     end if

     ! Check for variable freezing point depression, not currently implemented to be consistent
     ! with energy conservation
     if (tile%soil%pars%tfreeze /= tfreeze) then
        if (tfreeze_diff == 0.) then
           tfreeze_diff = (tile%soil%pars%tfreeze - tfreeze)
        else if (tfreeze_diff /= tile%soil%pars%tfreeze - tfreeze) then
           call error_mesg(module_name, 'Freezing point depression appears to be initialized with' // &
                           ' spatially variable values.  This is not currently implemented to be ' // &
                           'consistent with energy conservation with the Hillslope Model.', FATAL)
        end if
     end if

  end do

  ! ---- static [for now] diagnostic section
  ! List of fields:
   !id_soil_e_depth, id_microtopo, id_tile_hlsp_length, id_tile_hlsp_slope, &
   !id_tile_hlsp_elev, id_tile_hlsp_hpos, id_tile_hlsp_width, id_transm_bedrock, &
   !id_hidx_j, id_hidx_k
   ! soil_e_depth and k_sat_gw will be done in soil_init.
   call send_tile_data_i0d_fptr(id_hidx_j,              soil_hidx_j_ptr)
   call send_tile_data_i0d_fptr(id_hidx_k,              soil_hidx_k_ptr)
!   call send_tile_data_r0d_fptr(id_soil_e_depth,        soil_soil_e_depth_ptr)
   call send_tile_data_r0d_fptr(id_microtopo,           soil_microtopo_ptr)
   call send_tile_data_r0d_fptr(id_tile_hlsp_length,    soil_tile_hlsp_length_ptr)
   call send_tile_data_r0d_fptr(id_tile_hlsp_slope,     soil_tile_hlsp_slope_ptr)
   call send_tile_data_r0d_fptr(id_tile_hlsp_elev,      soil_tile_hlsp_elev_ptr)
   call send_tile_data_r0d_fptr(id_tile_hlsp_hpos,      soil_tile_hlsp_hpos_ptr)
   call send_tile_data_r0d_fptr(id_tile_hlsp_width,     soil_tile_hlsp_width_ptr)
!   call send_tile_data_r0d_fptr(id_transm_bedrock,      soil_transm_bedrock_ptr)

   deallocate(num_topo_hlsps, frac_topo_hlsps, &
            soil_e_depth, &
            microtopo, &
            hlsp_length, &
            hlsp_slope, &
            hlsp_slope_exp, &
            hlsp_top_width, &
            k_sat_gw)

end subroutine hlsp_init

! ============================================================================
! initialize hillslope model
subroutine hlsp_init_predefined(id_ug)
  integer,intent(in) :: id_ug !<Unstructured axis id.

  ! ---- local vars
  integer :: unit         ! unit for various i/o
  type(land_tile_enum_type)     :: ce  ! tail and current tile list elements
  type(land_tile_type), pointer :: tile   ! pointer to current tile
  character(len=256) :: mesg
  logical :: restart_exists
  logical :: found
  integer :: siz(4), tsize
  type(land_restart_type) :: restart

  ! change the initialization flag to true
  module_is_initialized = .TRUE.

  ! initialize hillslope-dependent diagnostic fields
  call hlsp_diag_init(id_ug)

  call open_land_restart(restart,hlsp_rst_ifname,restart_exists)
  if (restart_exists) then
     call error_mesg(module_name, 'hlsp_init, '// &
          'reading NetCDF restart "'//trim(hlsp_rst_ifname)//'"', &
          NOTE)
     if (cold_start) &
        call error_mesg(module_name, 'hlsp_init: coldfracs subroutine called even though restart file '// &
                             'exists! Inconsistency of "cold_start" in hillslope_mod.', FATAL)
     call get_int_tile_data(restart, 'HIDX_J', soil_hidx_j_ptr)
     call get_int_tile_data(restart, 'HIDX_K', soil_hidx_k_ptr)
  else
     call error_mesg(module_name, 'hlsp_init: '// &
          'cold-starting hillslope model',&
          NOTE)

  endif
  call free_land_restart(restart)

  ! ---- static [for now] diagnostic section
  ! List of fields:
   !id_soil_e_depth, id_microtopo, id_tile_hlsp_length, id_tile_hlsp_slope, &
   !id_tile_hlsp_elev, id_tile_hlsp_hpos, id_tile_hlsp_width, id_transm_bedrock, &
   !id_hidx_j, id_hidx_k
   ! soil_e_depth and k_sat_gw will be done in soil_init.
   call send_tile_data_i0d_fptr(id_hidx_j,         soil_hidx_j_ptr)
   call send_tile_data_i0d_fptr(id_hidx_k,         soil_hidx_k_ptr)
!   call send_tile_data_r0d_fptr(id_soil_e_depth,  soil_soil_e_depth_ptr)
   call send_tile_data_r0d_fptr(id_microtopo,      soil_microtopo_ptr)
   call send_tile_data_r0d_fptr(id_tile_hlsp_length, soil_tile_hlsp_length_ptr)
   call send_tile_data_r0d_fptr(id_tile_hlsp_slope,  soil_tile_hlsp_slope_ptr)
   call send_tile_data_r0d_fptr(id_tile_hlsp_elev,   soil_tile_hlsp_elev_ptr)
   call send_tile_data_r0d_fptr(id_tile_hlsp_hpos,   soil_tile_hlsp_hpos_ptr)
   call send_tile_data_r0d_fptr(id_tile_hlsp_width,  soil_tile_hlsp_width_ptr)
!   call send_tile_data_r0d_fptr(id_transm_bedrock,  soil_transm_bedrock_ptr)

 !Print out parameters
 if (is_watch_point()) then
   ce = first_elmt(land_tile_map)
   do while (loop_over_tiles(ce,tile))
      if (.not.associated(tile%soil))cycle
      print*,'soil_e_depth',tile%soil%pars%soil_e_depth
      print*,'microtopo',tile%soil%pars%microtopo
      print*,'k_sat_gw',tile%soil%pars%k_sat_gw
      print*,'hlsp_elev',tile%soil%pars%tile_hlsp_elev
      print*,'hlsp_hpos',tile%soil%pars%tile_hlsp_hpos
      print*,'hlsp_width',tile%soil%pars%tile_hlsp_width
      print*,'hlsp_slope',tile%soil%pars%tile_hlsp_slope
      print*,'hlsp_length',tile%soil%pars%tile_hlsp_length
      print*,'microtopo',tile%soil%pars%microtopo
      print*,'hidx_j',tile%soil%hidx_j
      print*,'hidx_k',tile%soil%hidx_k
   end do
 endif

end subroutine hlsp_init_predefined

! ============================================================================
subroutine hlsp_diag_init(id_ug)
   integer,intent(in) :: id_ug !<Unstructured axis id.

   ! ---- local vars
   integer :: axes(1)

   ! define array of axis indices
   axes = (/id_ug/)

   ! set the default sub-sampling filter for the fields below
   call set_default_diag_filter('soil')

   ! define static [ZMS: FOR TESTING] fields

   ! List of fields:
   !id_soil_e_depth, id_microtopo, id_tile_hlsp_length, id_tile_hlsp_slope, &
   !id_tile_hlsp_elev, id_tile_hlsp_hpos, id_tile_hlsp_width, id_transm_bedrock, &
   !id_hidx_j, id_hidx_k
!   id_soil_e_depth = register_tiled_static_field ( module_name, 'soil_e_depth',  &
!      axes, 'soil hydraulic depth-scale', 'm', missing_value=-100.0 )
   ! soil_e_depth and k_sat_gw done in soil_init
   id_microtopo = register_tiled_static_field ( module_name, 'microtopo', &
      axes, 'microtopographic roughness length-scale', 'm', missing_value=-100.0 )
   id_tile_hlsp_length = register_tiled_static_field ( module_name, 'tile_hlsp_length', &
      axes, 'horizontal length of tiles along the direction of hillslope', 'm', &
      missing_value=-100.0 )
   id_tile_hlsp_slope = register_tiled_static_field ( module_name, 'tile_hlsp_slope', &
      axes, 'topographic slope of tiles (tangent of slope angle)', '-', &
      missing_value=-100.0 )
   id_tile_hlsp_elev = register_tiled_static_field ( module_name, 'tile_hlsp_elev', &
      axes, 'vertical elevation of tiles in hillslope with respect to stream', 'm', &
      missing_value=-100.0 )
   id_tile_hlsp_hpos = register_tiled_static_field ( module_name, 'tile_hlsp_hposition', &
      axes, 'horizontal position of tile along the direction of hillslope', 'm', missing_value=-100.0 )
   id_tile_hlsp_width = register_tiled_static_field ( module_name, 'tile_hlsp_width', &
      axes, 'representative horizontal width of tile perpendicular to hillslope', '-', &
      missing_value=-100.0 )
!   id_transm_bedrock = register_tiled_static_field ( module_name, 'bedrock_transmissivity', &
!      axes, 'bedrock hydraulic transmissivity', 'm^2/s', missing_value=-100.0 )
   id_hidx_j = register_tiled_static_field ( module_name, 'hillslope_position', &
      axes, 'horizontal position index along hillslope', missing_value=0., op=OP_SUM )
   id_hidx_k = register_tiled_static_field ( module_name, 'hillslope_parent', &
      axes, 'index of hillslope parent', missing_value=0., op=OP_SUM )

end subroutine hlsp_diag_init


! ============================================================================
subroutine hlsp_end ()

  module_is_initialized =.FALSE.

end subroutine hlsp_end


! ============================================================================
subroutine save_hlsp_restart (tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars ----------------------------------------------------------
  character(267) :: filename
  type(land_restart_type) :: restart ! restart file i/o object

  if (.not. do_hillslope_model) return

  call error_mesg(module_name,'writing hillslope restart',NOTE)
! must set domain so that io_domain is available
! Note that filename is updated for tile & rank numbers during file creation
  filename = trim(timestamp)//hlsp_rst_ofname
  call init_land_restart(restart, filename, soil_tile_exists, tile_dim_length)

  call add_int_tile_data(restart,'HIDX_J',soil_hidx_j_ptr,'hillslope position index','-')
  call add_int_tile_data(restart,'HIDX_K',soil_hidx_k_ptr,'hillslope parent index','-')
  ! save performs io domain aggregation through mpp_io as with regular domain data
  call save_land_restart(restart)
  call free_land_restart(restart)

end subroutine save_hlsp_restart

! ============================================================================
subroutine hlsp_config_check()

  if (.not. do_hillslope_model) return

  ! Error checking
  if ( (do_landuse_change .or. do_harvesting) .and. (.not. hillslope_horz_subdiv) .and. &
        (.not. hillslope_topo_subdiv) ) &
     call error_mesg(module_name, 'Land use change and harvesting require '// &
                     'the hillslope model to allow horizontal or topographic hillslope subdivision.',&
                     FATAL)

  if ((.not. fixed_num_vertclusters) .or. hillslope_topo_subdiv) &
     call error_mesg(module_name, 'hlsp_init: fixed_num_vertclusters == .true.,'// &
                     'hillslope_horz_subdiv == .true., and hillslope_topo_subdiv == .false'// &
                     'is currently required.', &
                     FATAL)

  ! ZMS fill in this function
  ! NWC - Why is was this ever conditional??? Eventually set externally
  if ((do_landuse_change .or. do_harvesting) .and. hillslope_horz_subdiv) then
      call transitions_disturbance_length_init()
  end if

  ! Deallocate variables used during init, as this function is called at end of land_model init
  ! sequence.
  if (cold_start) then
     deallocate(fracs_loc, elev_loc)
     if (allocated(init_wt_loc)) deallocate(init_wt_loc)
  end if

end subroutine hlsp_config_check

! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function soil_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   soil_tile_exists = associated(tile%soil)
end function soil_tile_exists


! ============================================================================
! For now just set disturb_scale to 1/2 hillslope_length. Later this will be set externally
! in transitions code, or in vegetation / peat code.
subroutine transitions_disturbance_length_init()
   type(land_tile_enum_type)     :: ce   ! tile list enumerator
   type(land_tile_type), pointer :: tile ! pointer to current tile

   ce = first_elmt(land_tile_map)
   do while(loop_over_tiles(ce,tile))
      if (associated(tile%soil)) then
         tile%soil%pars%disturb_scale = tile%soil%pars%tile_hlsp_length * num_vertclusters / 2.
      end if
   end do

end subroutine transitions_disturbance_length_init


! ============================================================================
! Return water table depth init for soil tile to soil_init.
subroutine horiz_wt_depth_to_init (soil, ce, init_wt, dry)
   ! Arguments
   type(soil_tile_type), intent(in)       :: soil
   type(land_tile_enum_type), intent(in)  :: ce  ! enumerator for soil, for retrieving indices
   real, intent(out)                      :: init_wt  ! [m] local wt depth to prescribe
   logical, intent(in), optional          :: dry ! flag indicating dry point from cold start data
   ! Local variables
   integer :: hk,hj,li,lj,t,ll ! coordinates of current tile
   logical :: drylocal = .false.

   if (present(dry)) drylocal = dry

   if (cold_start) then
      call get_elmt_indices(ce,i=li,j=lj,k=t,l=ll)
      call set_current_point(ll,t)
      if (init_wt_strmelev .or. drylocal) then
         init_wt = soil%pars%tile_hlsp_elev
      else
         hk = soil%hidx_k
         hj = soil%hidx_j
         init_wt = init_wt_loc(ll, hk, hj) ! Calculated previously in calculate_wt_init
      end if
      ! Debug
      if (is_watch_cell()) then
         write(*,*)'hlsp_mod: horiz_wt_depth_to_init: Assigning water table depth at watch cell i,j:', &
                    li,lj, '.'
         write(*,*)'hillslope indices hk, hj:',hk, hj, '.'
         write(*,*)'tile elev = ', soil%pars%tile_hlsp_elev
         write(*,*)'water table depth to init = ', init_wt
      end if
   else
      init_wt = 0. ! will be written over by restart soil state variables
   end if

end subroutine horiz_wt_depth_to_init


! ============================================================================
! If cold start and .not. init_wt_strmelev, calculate water table depth to init. Modifies
! module variable init_wt_local. Called from soil_init.
subroutine calculate_wt_init(wtdep_requested)
   real, intent(in)   :: wtdep_requested ! [m] average or horizontal water table depth
   integer   :: lls, lle  ! lat, lon index bounds
   integer   :: ll, li, lj, hk, hj ! indices
   real, allocatable :: depth(:) ! [m] prognostic depth for location of water table, relative to tile elev
   real, allocatable :: elev(:) ! [m] tile elevation in this hillslope
   real, allocatable :: fracs(:) ! gridcell area fractions in this hillslope
   real :: melev            ! [m] hillslope or partial-hillslope mean elevation
   real :: wtelev           ! [m] elevation of (non-seepage) water table in hillslope
   real :: fracs_sat, fracs_unsat ! sum of seepage and non-seepage area fractions
   real :: wteff            ! [m] calculated water table depth
   real, parameter :: tolerance = 1.e-9 ! [m] allowed roundoff error
   logical :: lerror = .false. ! error in calculation
!   character(len=512) :: mesg
   lls = lnd%ls
   lle = lnd%le

   if (cold_start) then ! init_wt_local will be used
      allocate(init_wt_loc(lls:lle, max_num_topo_hlsps, num_vertclusters))
      init_wt_loc(:,:,:) = 0.
      if (.not. init_wt_strmelev) then ! extra-tile information will be required.
         allocate(depth(1:num_vertclusters), elev(1:num_vertclusters), fracs(1:num_vertclusters))
         do ll = lls, lle
            do hk=1,max_num_topo_hlsps
               if (fracs_loc(ll, hk, 1) > 0 .and. elev_loc(ll, hk, 1) /= initval) then
               ! this hillslope is active
               ! Retrieve elev and fracs
                  fracs(:) = fracs_loc(ll, hk, :)
                  elev(:) = elev_loc(ll, hk, :)
                  ! Initialize wtelev to meanelev - wtdep_requested.
                  melev = meanelev(elev, fracs, 1, num_vertclusters)
                  wtelev = melev - wtdep_requested
                  ! Initialize depth
                  depth(:) = elev(:) - wtelev

                  ! Proceed in loop to catch seepage faces where depth < 0
                  hj = 1
                  do while(hj <= num_vertclusters .and. depth(hj) < 0.)
                     depth(hj) = 0.
                     if (hj < num_vertclusters) then
                        melev = meanelev(elev, fracs, hj + 1, num_vertclusters)
                        fracs_sat = sum(fracs(1:hj))
                        fracs_unsat = sum(fracs(hj+1:num_vertclusters))
                        wtelev = melev - wtdep_requested * (fracs_sat+fracs_unsat)/fracs_unsat
                        depth(hj+1:num_vertclusters) = elev(hj+1:num_vertclusters) - wtelev
                     end if
                     hj = hj + 1
                  end do

                  ! Check solution
                  if (wtdep_requested >= 0.) then
                     wteff = sum(fracs(:)*depth(:))/sum(fracs(:))
                     if (abs(wteff - wtdep_requested) > tolerance) lerror = .true.
                     if (.not. all(depth >= 0.)) lerror = .true.
                  else if (.not. all(depth == 0.)) then
                     lerror = .true.
                  end if
                  if (lerror) then
                     li = lnd%i_index(ll)
                     lj = lnd%j_index(ll)
                     write(*,*)'Error in calculation of water table depth at li, lj, hk = ', &
                                  li, lj, hk, '; wtdep_requested = ', wtdep_requested, '.'
                     write(*,*)'fracs = ', fracs(:)
                     write(*,*)'elev = ', elev(:)
                     write(*,*)'depth = ', depth(:)
                     call error_mesg(module_name, 'Error in calculating water table to init.', &
                                     FATAL)
                  end if

                  ! Set init_wt_loc
                  init_wt_loc(ll, hk, :) = depth(:)

               end if ! fracs_loc
            end do !hk
         end do !ll

         deallocate(depth, elev, fracs)
      end if ! init_wt_strmelev
   end if ! cold_start

end subroutine calculate_wt_init

! Calculate the mean of elev weighted by area from lowbound to upbound
! ============================================================================
function meanelev(elev, area, lowbound, upbound) result(melev)
   real :: melev ! mean elevation
   real, intent(in) :: elev(:), area(:)
   integer, intent(in) :: lowbound, upbound

   if (sum(area(lowbound:upbound)) > 0) then
      melev = sum(elev(lowbound:upbound)*area(lowbound:upbound))/ &
              sum(area(lowbound:upbound))
   else
      call error_mesg(module_name, ': meanelev called with zero area. Aborting.', FATAL)
   end if

end function meanelev

! ============================================================================
! cohort accessor functions: given a pointer to cohort, return a pointer to a
! specific member of the cohort structure
#define DEFINE_SOIL_ACCESSOR_0D(xtype,x) subroutine soil_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%x;endif;\
end subroutine
#define DEFINE_SOIL_COMPONENT_ACCESSOR_0D(xtype,component,x) subroutine soil_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%soil))p=>t%soil%component%x;endif;\
end subroutine

DEFINE_SOIL_ACCESSOR_0D(integer,hidx_j)
DEFINE_SOIL_ACCESSOR_0D(integer,hidx_k)

!DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars, soil_e_depth)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars, microtopo)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars, tile_hlsp_length)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars, tile_hlsp_slope)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars, tile_hlsp_elev)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars, tile_hlsp_hpos)
DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars, tile_hlsp_width)
!DEFINE_SOIL_COMPONENT_ACCESSOR_0D(real,pars, transm_bedrock)


end module hillslope_mod



