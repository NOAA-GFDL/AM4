! ============================================================================
! snow model module
! ============================================================================
module snow_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : error_mesg, file_exist, check_nml_error, &
     stdlog, close_file, mpp_pe, mpp_root_pe, FATAL, NOTE
use time_manager_mod,   only: time_type_to_real
use constants_mod,      only: tfreeze, hlv, hlf, PI

use land_constants_mod, only : NBANDS
use snow_tile_mod, only : &
     snow_tile_type, read_snow_data_namelist, &
     snow_data_thermodynamics, snow_data_area, &
     snow_active, &
     snow_data_hydraulics, max_lev, cpw, clw, csw, use_brdf

use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_enum_type, &
     first_elmt, loop_over_tiles
use land_data_mod, only : lnd, log_version
use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_restart_axis, add_tile_data, get_tile_data
use land_debug_mod, only : is_watch_point

implicit none
private

! ==== public interfaces =====================================================
public :: read_snow_namelist
public :: snow_init
public :: snow_end
public :: save_snow_restart
public :: snow_get_depth_area
public :: snow_sfc_water
public :: sweep_tiny_snow
public :: snow_step_1
public :: snow_step_2
! =====end of public interfaces ==============================================


! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'snow_mod'
#include "../shared/version_variable.inc"

! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
logical :: retro_heat_capacity  = .false.
logical :: lm2  = .false.
logical :: steal = .false.
character(len=16):: albedo_to_use = ''  ! or 'brdf-params'
real :: max_snow       = 1000.
real :: wet_max        = 0.0  ! TEMP, move to snow_data
real :: snow_density   = 300. ! TEMP, move to snow_data and generalize
real :: init_temp = 260.   ! cold-start snow T
real :: init_pack_ws   =   0.
real :: init_pack_wl   =   0.
real :: min_snow_mass = 0.
logical :: prevent_tiny_snow = .FALSE. ! if true, tiny snow is removed at the 
   ! beginning of fast time step to avoid numerical issues. There is no harm
   ! in doing that, but it changes answers, so for compatibility with older code 
   ! turn it off.

namelist /snow_nml/ retro_heat_capacity, lm2, steal, albedo_to_use, &
                    max_snow, wet_max, snow_density, &
                    init_temp, init_pack_ws, init_pack_wl, &
                    min_snow_mass, prevent_tiny_snow
!---- end of namelist --------------------------------------------------------

logical         :: module_is_initialized =.FALSE.
real            :: delta_time
integer         :: num_l    ! # of snow layers
! next three 'z' variables are all normalized by total snow pack depth
real            :: dz (max_lev) ! relative thicknesses of layers
real            :: z  (max_lev) ! relative depths of layer bounds
real            :: zz (max_lev) ! relative depths of layer centers
real            :: heat_capacity_retro = 1.6e6
real            :: mc_fict

! ==== end of module variables ===============================================

contains

! ============================================================================
subroutine read_snow_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: l            ! layer iterator

  call read_snow_data_namelist(num_l,dz,mc_fict)

  call log_version(version, module_name, &
  __FILE__)
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=snow_nml, iostat=io)
  ierr = check_nml_error(io, 'snow_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=snow_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'snow_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=snow_nml)
  endif

  ! -------- set up vertical discretization --------
  zz(1) = 0
  do l = 1, num_l
     zz(l+1) = zz(l) + dz(l)
     z(l)    = 0.5*(zz(l+1) + zz(l))
  enddo

end subroutine read_snow_namelist


! ============================================================================
! initialize snow model
subroutine snow_init()

  ! ---- local vars ----------------------------------------------------------
  integer :: k
  type(land_tile_enum_type)     :: ce    ! tile list enumerator
  type(land_tile_type), pointer :: tile  ! pointer to current tile
  character(*), parameter :: restart_file_name='INPUT/snow.res.nc'
  type(land_restart_type) :: restart
  logical :: restart_exists

  module_is_initialized = .TRUE.
  delta_time = time_type_to_real(lnd%dt_fast)

  ! -------- initialize snow state --------
  call open_land_restart(restart,restart_file_name,restart_exists)
  if (restart_exists) then
     call error_mesg('snow_init', 'reading NetCDF restart "'//trim(restart_file_name)//'"', NOTE)
     call get_tile_data(restart, 'temp', 'zfull', snow_temp_ptr)
     call get_tile_data(restart, 'wl'  , 'zfull', snow_wl_ptr)
     call get_tile_data(restart, 'ws'  , 'zfull', snow_ws_ptr)
  else
     call error_mesg('snow_init', 'cold-starting snow', NOTE)
     ce = first_elmt(land_tile_map)
     do while(loop_over_tiles(ce, tile))
        if (.not.associated(tile%snow)) cycle
        do k = 1,num_l
           tile%snow%wl(k) = init_pack_wl * dz(k)
           tile%snow%ws(k) = init_pack_ws * dz(k)
           tile%snow%T(k)  = init_temp
        enddo
     enddo
  endif
  call free_land_restart(restart)

  if (trim(albedo_to_use)=='') then
     use_brdf = .false.
  elseif (trim(albedo_to_use)=='brdf-params') then
     use_brdf = .true.
  else
     call error_mesg('snow_init',&
          'option albedo_to_use="'//&
          trim(albedo_to_use)//'" is invalid, use "" or "brdf-params"',&
          FATAL)
  endif

end subroutine snow_init


! ============================================================================
subroutine snow_end ()

  module_is_initialized =.FALSE.

end subroutine snow_end


! ============================================================================
subroutine save_snow_restart (tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars
  character(267) :: filename
  type(land_restart_type) :: restart ! restart file i/o object

  call error_mesg('snow_end','writing NetCDF restart',NOTE)
! Note that filename is updated for tile & rank numbers during file creation
  filename = trim(timestamp)//'snow.res.nc'
  call init_land_restart(restart, filename, snow_tile_exists, tile_dim_length)
  call add_restart_axis(restart,'zfull',zz(1:num_l),'Z',longname='depth of level centers',sense=-1)

  call add_tile_data(restart,'temp','zfull', snow_temp_ptr, 'snow temperature','degrees_K')
  call add_tile_data(restart,'wl'  ,'zfull', snow_wl_ptr,   'snow liquid water content','kg/m2')
  call add_tile_data(restart,'ws'  ,'zfull', snow_ws_ptr,   'snow solid water content','kg/m2')

  call save_land_restart(restart)
  call free_land_restart(restart)

end subroutine save_snow_restart

! ============================================================================
subroutine snow_get_depth_area(snow, snow_depth, snow_area)
  type(snow_tile_type), intent(in) :: snow
  real, intent(out) :: snow_depth, snow_area

  integer :: l

  snow_depth= 0.0
  do l = 1, num_l
     snow_depth = snow_depth + snow%ws(l)
  enddo
  snow_depth = snow_depth / snow_density
  call snow_data_area (snow_depth, snow_area )
end subroutine


! ============================================================================
! if snow amount is below specified limit, sweeps it into runoff
subroutine sweep_tiny_snow(snow, lrunf, frunf, hlrunf, hfrunf)
  type(snow_tile_type), intent(inout) :: snow
  real, intent(out) :: lrunf, frunf, hlrunf, hfrunf

  real :: snow_mass
  integer :: l

  lrunf=0 ; frunf=0 ; hlrunf=0 ; hfrunf=0
  if (.not.prevent_tiny_snow) return ! do nothing, return zeros

  snow_mass  = sum(snow%ws)
  ! check if the snow is small enough to warrant sweeping
  if ( snow_mass<0 .or. snow_mass >= min_snow_mass ) return

  lrunf  = sum(snow%wl) ; frunf  = snow_mass
  hlrunf = 0.0 ; hfrunf = 0.0
  do l = 1, num_l
     hlrunf = hlrunf + clw*snow%wl(l)*(snow%T(l)-tfreeze)
     hfrunf = hfrunf + csw*snow%ws(l)*(snow%T(l)-tfreeze)
  enddo
  snow%ws = 0 ; snow%wl = 0
end subroutine sweep_tiny_snow

! ============================================================================
! update snow properties explicitly for time step.
! integrate snow-heat conduction equation upward from bottom of snow
! to surface, delivering linearization of surface ground heat flux.
subroutine snow_step_1 ( snow, snow_G_Z, snow_G_TZ, &
                         snow_rh, &
                         snow_area, snow_G0, snow_DGDT )
  type(snow_tile_type), intent(inout) :: snow
  real,                 intent(in) :: snow_G_Z
  real,                 intent(in) :: snow_G_TZ
  real,                 intent(out):: &
       snow_rh, snow_area, snow_G0, snow_DGDT

  ! ---- local vars
  real :: snow_depth, bbb, denom, dt_e
  real, dimension(num_l):: aaa, ccc, thermal_cond, dz_phys, heat_capacity
  integer :: l

! ----------------------------------------------------------------------------
! in preparation for implicit energy balance, determine various measures
! of water availability, so that vapor fluxes will not exceed mass limits
! ----------------------------------------------------------------------------

  call snow_data_thermodynamics ( snow_rh, thermal_cond )
  snow_depth= 0.0
  do l = 1, num_l
     snow_depth = snow_depth + snow%ws(l)
  enddo
  snow_depth = snow_depth / snow_density
  call snow_data_area (snow_depth, snow_area )

  do l = 1, num_l
     dz_phys(l) = dz(l)*snow_depth
  enddo

  if (retro_heat_capacity) then
     do l = 1, num_l
        heat_capacity(l) = heat_capacity_retro*dz_phys(l)
     enddo
  else
     do l = 1, num_l
        heat_capacity(l) = mc_fict*dz(l) + &
             clw*snow%wl(l) + csw*snow%ws(l)
     enddo
  endif

!  if(num_l > 1) then
  if (snow_depth > 0) then
     do l = 1, num_l-1
        dt_e = 2 / ( dz_phys(l+1)/thermal_cond(l+1) &
                     + dz_phys(l)/thermal_cond(l)   )
        aaa(l+1) = - dt_e * delta_time / heat_capacity(l+1)
        ccc(l)   = - dt_e * delta_time / heat_capacity(l)
     enddo

     bbb = 1.0 - aaa(num_l) + delta_time*snow_G_TZ/heat_capacity(num_l)
     denom = bbb
     dt_e = aaa(num_l)*(snow%T(num_l) - snow%T(num_l-1)) &
          - delta_time*snow_G_Z/heat_capacity(num_l)
     snow%e(num_l-1) = -aaa(num_l)/denom
     snow%f(num_l-1) = dt_e/denom

     do l = num_l-1, 2, -1
        bbb = 1.0 - aaa(l) - ccc(l)
        denom = bbb + ccc(l)*snow%e(l)
        dt_e = - ( ccc(l)*(snow%T(l+1) - snow%T(l)  ) &
                  -aaa(l)*(snow%T(l)   - snow%T(l-1)) )
        snow%e(l-1) = -aaa(l)/denom
        snow%f(l-1) = (dt_e - ccc(l)*snow%f(l))/denom
     enddo

     denom = delta_time/heat_capacity(1)
     snow_G0    = ccc(1)*(snow%T(2)- snow%T(1) &
          + snow%f(1)) / denom
     snow_DGDT  = (1 - ccc(1)*(1-snow%e(1))) / denom
  endif

!    else  ! one-level case
!      denom = delta_time/heat_capacity(1)
!      snow_G0    = 0.
!      snow_DGDT  = 1. / denom
!    end if

  if (snow_depth <= 0) then
     snow_G0   = snow_G_Z
     snow_DGDT = snow_G_TZ
  endif

  if(is_watch_point()) then
     write(*,*) 'snow_depth', snow_depth
     write(*,*) '############ snow_step_1 output'
     write(*,*) 'mask      ', .true.
     write(*,*) 'snow_rh   ', snow_rh
     write(*,*) 'snow_area ', snow_area
     write(*,*) 'snow_G_Z  ', snow_G_Z
     write(*,*) 'snow_G_TZ ', snow_G_TZ
     write(*,*) 'snow_G0   ', snow_G0
     write(*,*) 'snow_DGDT ', snow_DGDT
     write(*,*) '############ end of snow_step_1 output'
  endif

end subroutine snow_step_1


! ============================================================================
! apply boundary flows to snow water and move snow water vertically.
  subroutine snow_step_2 ( snow,                     &
                           vegn_lprec, vegn_fprec, vegn_hlprec, vegn_hfprec, &
                           DTg,  Mg_imp,  evapg,  fswg,  flwg,  sensg,  &
                           use_tfreeze_in_grnd_latent, subs_DT, &
                           subs_M_imp, subs_evap, subs_fsw, subs_flw, subs_sens,  &
                           snow_fsw, snow_flw, snow_sens, &
                           snow_levap, snow_fevap, snow_melt, &
                           snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, &
                           snow_hlrunf, snow_hfrunf, snow_Tbot, snow_Cbot, snow_C, &
                           snow_avrg_T )
  type(snow_tile_type), intent(inout) :: snow
  real, intent(in) :: &
     vegn_lprec, vegn_fprec, vegn_hlprec, vegn_hfprec
  real, intent(in) :: &
     DTg, Mg_imp, evapg, fswg, flwg, sensg
  logical, intent(in) :: use_tfreeze_in_grnd_latent
  real, intent(out) :: &
         subs_DT, subs_M_imp, subs_evap, subs_fsw, subs_flw, subs_sens, &
         snow_fsw, snow_flw, snow_sens, &
         snow_levap, snow_fevap, snow_melt, &
         snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, &
         snow_hlrunf, snow_hfrunf, snow_Tbot, snow_Cbot, snow_C, snow_avrg_T

  ! ---- local vars
  real, dimension(num_l) :: del_t, M_layer
  real :: depth, snow_subl, &
         cap0, dW_l, dW_s, dcap, dheat,&
         melt, melt_per_deg, drain,       &
         snow_mass, sum_liq, &
         sum_heat, sum_sno, &
         snow_transfer, frac,&
         liq_rate, hliq_rate,&
         sno_rate, hsno_rate, fict_heat, &
         evapg_lm2, vegn_fprec_lm2, &
         snow_LMASS, snow_FMASS, snow_HEAT
  integer :: l, l_old
  real :: new_ws(num_l)
  real :: new_wl(num_l)
  real :: new_T(num_l)
  ! --------------------------------------------------------------------------

  snow_subl = snow_subl_frac(snow)
  depth= 0.
  do l = 1, num_l
    depth = depth + snow%ws(l)
  enddo
  depth = depth / snow_density

  if(is_watch_point()) then
     write(*,*) '############ snow_step_2 input'
     write(*,*) 'mask       ', .TRUE.
     write(*,*) 'snow_subl  ', snow_subl
     write(*,*) 'vegn_lprec ', vegn_lprec
     write(*,*) 'vegn_fprec ', vegn_fprec
     write(*,*) 'vegn_hlprec', vegn_hlprec
     write(*,*) 'vegn_hfprec', vegn_hfprec
     write(*,*) 'DTg        ', DTg
     write(*,*) 'Mg_imp     ', Mg_imp
     write(*,*) 'evapg      ', evapg
     write(*,*) 'fswg       ', fswg
     write(*,*) 'flwg       ', flwg
     write(*,*) 'sensg      ', sensg
     write(*,*) '############ end of snow_step_2 input'

     write(*,*) 'depth   ', depth
     do l = 1, num_l
        write(*,'(i2,3(x,a,g23.16))') l,&
             ' wl=', snow%wl(l),&
             ' ws=', snow%ws(l),&
             ' T =', snow%T(l)
     enddo
  endif

  snow_LMASS = 0; snow_FMASS = 0; snow_HEAT = 0
  do l = 1, num_l;
        snow_LMASS = snow_LMASS + snow%wl(l)
        snow_FMASS = snow_FMASS + snow%ws(l)
        snow_HEAT = snow_HEAT + &
          (mc_fict*dz(l) + clw*snow%wl(l) + csw*snow%ws(l))  &
                                                * (snow%T(l)-tfreeze)
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 1.01 ***** '
     write(*,*) 'LMASS         ', snow_LMASS
     write(*,*) 'FMASS         ', snow_FMASS
     write(*,*) 'HEAT          ', snow_HEAT
  endif

  ! ---- record fluxes -------------------------------------------------------
  if (lm2.and.steal) then
    if (snow_FMASS-Mg_imp > 0.) then
      if (evapg <= (snow_FMASS-Mg_imp)/delta_time) then
          evapg_lm2 = evapg
          vegn_fprec_lm2 = vegn_fprec
      else if (evapg <= (snow_FMASS-Mg_imp)/delta_time+vegn_fprec) then
          evapg_lm2 = evapg
          vegn_fprec_lm2 = vegn_fprec - evapg + (snow_FMASS-Mg_imp)/delta_time
      else
          evapg_lm2 = (snow_FMASS-Mg_imp)/delta_time+vegn_fprec
          vegn_fprec_lm2 = 0.
      endif
    else
      evapg_lm2 = 0.
      vegn_fprec_lm2 = vegn_fprec
    endif
  else
     evapg_lm2 = evapg
     vegn_fprec_lm2 = vegn_fprec
  endif
  vegn_fprec_lm2 = vegn_fprec
  if (depth>0) then
        snow_fsw   = fswg
        snow_flw   = flwg
        snow_sens  = sensg
        snow_levap = evapg_lm2*(1-snow_subl)
        snow_fevap = evapg_lm2*   snow_subl
  else
        snow_fsw    = 0
        snow_flw    = 0
        snow_sens   = 0
        snow_levap  = 0
        snow_fevap  = 0
  endif
  subs_fsw = fswg - snow_fsw
  subs_flw = flwg - snow_flw
  subs_evap = evapg - snow_levap - snow_fevap
  subs_sens = sensg - snow_sens

  ! ---- load surface temp change and perform back substitution --------------
  if (depth>0) then
      del_t(1) = DTg
      snow%T(1)  = snow%T(1) + del_t(1)
  endif
  if ( num_l > 1) then
     do l = 1, num_l-1
        if (depth>0) then
            del_t(l+1) = snow%e(l) * del_t(l) + snow%f(l)
            snow%T(l+1) = snow%T(l+1) + del_t(l+1)
        endif
     enddo
  endif
  if (depth>0) then
    subs_DT = del_t(num_l)
  else
    subs_DT = DTg
  endif

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 2 ***** '
     do l = 1, num_l
        write(*,'(i2,a,g23.16)') l,' T =', snow%T(l)
     enddo
  endif

  snow_LMASS = 0; snow_FMASS = 0; snow_HEAT = 0
  do l = 1, num_l
        snow_LMASS = snow_LMASS + snow%wl(l)
        snow_FMASS = snow_FMASS + snow%ws(l)
        snow_HEAT = snow_HEAT + &
          (mc_fict*dz(l) + clw*snow%wl(l) + csw*snow%ws(l))  &
                                                * (snow%T(l)-tfreeze)
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 2.01 ***** '
     write(*,*) 'LMASS         ', snow_LMASS
     write(*,*) 'FMASS         ', snow_FMASS
     write(*,*) 'HEAT          ', snow_HEAT
  endif

  ! ---- evaporation and sublimation -----------------------------------------
  if (depth>0) then
      snow%wl(1) = snow%wl(1) - snow_levap*delta_time
      snow%ws(1) = snow%ws(1) - snow_fevap*delta_time
      cap0 = mc_fict*dz(1) + clw*snow%wl(1) + csw*snow%ws(1)
      ! T adjustment for nonlinear terms (del_T)*(del_W)
      dheat = delta_time*(clw*snow_levap+csw*snow_fevap)*del_T(1)
      ! take out extra heat not claimed in advance for evaporation
      if (use_tfreeze_in_grnd_latent) dheat = dheat &
            - delta_time*((cpw-clw)*snow_levap+(cpw-csw)*snow_fevap) &
                               *(snow%T(1)-del_T(1)-tfreeze)
      snow%T(1)  = snow%T(1)  + dheat/cap0
    endif

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 2.5 ***** '
     do l = 1, num_l
        write(*,'(i2,3(a,g23.16))')l,&
             ' wl=', snow%wl(l),&
             ' ws=', snow%ws(l),&
             ' T =', snow%T(l)
     enddo
  endif

  snow_LMASS = 0; snow_FMASS = 0; snow_HEAT = 0
  do l = 1, num_l
        snow_LMASS = snow_LMASS + snow%wl(l)
        snow_FMASS = snow_FMASS + snow%ws(l)
        snow_HEAT = snow_HEAT + &
          (mc_fict*dz(l) + clw*snow%wl(l) + csw*snow%ws(l))  &
                                                * (snow%T(l)-tfreeze)
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 2.51 ***** '
     write(*,*) 'LMASS         ', snow_LMASS
     write(*,*) 'FMASS         ', snow_FMASS
     write(*,*) 'HEAT          ', snow_HEAT
  endif

  ! ---- distribute implicit phase change downward through snow layers -------
  if (depth>0) then
      snow_melt = Mg_imp/delta_time
  else
      snow_melt = 0
  endif
  M_layer = 0.
  subs_M_imp = Mg_imp
  do l = 1, num_l
    if (depth>0 .and. subs_M_imp.gt.0) then
        M_layer(l) =  min( subs_M_imp, max(0.,snow%ws(l)) )
        subs_M_imp = subs_M_imp - M_layer(l)
    endif
  enddo
  if (depth>0) then
      M_layer(1) = M_layer(1) + subs_M_imp
      subs_M_imp = 0.
  endif
  do l = 1, num_l
    if (depth>0) then
          cap0 = mc_fict*dz(l) + clw*snow%wl(l) + csw*snow%ws(l)
          snow%wl(l) = snow%wl(l) + M_layer(l)
          snow%ws(l) = snow%ws(l) - M_layer(l)
          snow%T(l)  = tfreeze + (cap0*(snow%T(l)-tfreeze) ) &
                                                          / ( cap0 + (clw-csw)*M_layer(l) )
    endif
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 3 ***** '
     do l = 1, num_l
        write(*,'(i2,3(a,g23.16))') l,&
             ' wl=', snow%wl(l),&
             ' ws=', snow%ws(l),&
             '  T=', snow%T(l)
     enddo
  endif

  snow_LMASS = 0; snow_FMASS = 0; snow_HEAT = 0
  do l = 1, num_l
    snow_LMASS = snow_LMASS + snow%wl(l)
    snow_FMASS = snow_FMASS + snow%ws(l)
        snow_HEAT = snow_HEAT + &
          (mc_fict*dz(l) + clw*snow%wl(l) + csw*snow%ws(l))  &
                                                * (snow%T(l)-tfreeze)
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 3.01 ***** '
     write(*,*) 'LMASS         ', snow_LMASS
     write(*,*) 'FMASS         ', snow_FMASS
     write(*,*) 'HEAT          ', snow_HEAT
  endif

! ----------------------------------------------------------------------------
!  call snow_data_hydraulics (pars, snow%wl, psi, hyd_cond )

! ---- remainder of mass fluxes and associated sensible heat fluxes ----------
  liq_rate = vegn_lprec
  sno_rate = vegn_fprec_lm2
  hliq_rate = vegn_hlprec
  if (vegn_fprec.ne.0.) then
          hsno_rate = vegn_hfprec*(vegn_fprec_lm2/vegn_fprec)
  else
          hsno_rate = 0.
  endif

  do l = 1, num_l
    if(depth>0 .or. vegn_fprec_lm2>0) then
    ! ---- mix inflow with existing snow and water ---------------------------
          cap0 = mc_fict*dz(l) + clw*snow%wl(l) + csw*snow%ws(l)
          dW_l = liq_rate*delta_time
          dW_s = sno_rate*delta_time
          dcap = clw*dW_l + csw*dW_s
          snow%ws(l) = snow%ws(l) + dW_s
          snow%wl(l) = snow%wl(l) + dW_l
          snow%T(l)  = tfreeze + (cap0*(snow%T(l)-tfreeze) &
                               + (hsno_rate+hliq_rate)*delta_time) /(cap0 + dcap)
    endif

    if(is_watch_point()) then
       write(*,*) ' ***** snow_step_2 checkpoint 4a ***** '
       write(*,'(i2,3(a,g23.16))') l,&
            ' wl=', snow%wl(l),&
            ' ws=', snow%ws(l),&
            '  T=', snow%T(l)
    endif

    if (depth>0 .or. vegn_fprec_lm2>0) then
    ! ---- compute explicit melt/freeze --------------------------------------
          melt_per_deg = (cap0+dcap)/hlf
          if (snow%ws(l)>0 .and. snow%T(l)>tfreeze) then
                  melt =  min(snow%ws(l), (snow%T(l)-tfreeze)*melt_per_deg)
      elseif (snow%wl(l)>0 .and. snow%T(l)<tfreeze) then
                  melt = -min(snow%wl(l), (tfreeze-snow%T(l))*melt_per_deg)
      else
                  melt = 0
          endif
          snow_melt = snow_melt + melt/delta_time
          snow%wl(l) = snow%wl(l) + melt
          snow%ws(l) = snow%ws(l) - melt
!        where (cap0+dcap.ne.0.) &
!        snow%T(l)  = snow%T(l)  - melt/melt_per_deg
          snow%T(l) = tfreeze &
                 + ((cap0+dcap)*(snow%T(l)-tfreeze) - hlf*melt) &
                                                          / ( cap0+dcap + (clw-csw)*melt )
    endif

   if(is_watch_point()) then
      write(*,*) ' ***** snow_step_2 checkpoint 4b ***** '
      write(*,'(i2,3(a,g23.16))')l,&
            ' wl=', snow%wl(l),&
            ' ws=', snow%ws(l),&
            '  T=', snow%T(l)
   endif

   if (depth>0 .or. vegn_fprec_lm2>0) then
    ! ---- compute drainage from this layer to next --------------------------
        drain = max (0., snow%wl(l) - wet_max*snow%ws(l))
        snow%wl(l) = snow%wl(l) - drain
        liq_rate = drain / delta_time
        hliq_rate = clw*liq_rate*(snow%T(l)-tfreeze)
        sno_rate = 0
        hsno_rate = 0
    endif
  enddo

  snow_lprec  = liq_rate
  snow_hlprec = hliq_rate

  snow_LMASS = 0; snow_FMASS = 0; snow_HEAT = 0
  do l = 1, num_l
        snow_LMASS = snow_LMASS + snow%wl(l)
        snow_FMASS = snow_FMASS + snow%ws(l)
        snow_HEAT = snow_HEAT + &
          (mc_fict*dz(l) + clw*snow%wl(l) + csw*snow%ws(l))  &
                                                * (snow%T(l)-tfreeze)
  enddo


  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 4c ***** '
     write(*,*) 'LMASS         ', snow_LMASS
     write(*,*) 'FMASS         ', snow_FMASS
     write(*,*) 'HEAT          ', snow_HEAT
  endif

! ---- conceptually remove fictitious mass/heat for the moment ---------------
  fict_heat = 0.
  do l = 1, num_l
    fict_heat = fict_heat + dz(l)*snow%T(l)     ! (*mc_fict)
  enddo

  snow_mass  = sum(snow%ws)
  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 4d ***** '
     write(*,*) 'max_snow    ', max_snow
     write(*,*) 'snow_mass   ', snow_mass
     do l=1, num_l
       write(*,'(i2,3(a,g23.16))') l,&
            ' wl=', snow%wl(l),&
            ' ws=', snow%ws(l),&
            '  T=', snow%T(l)
     enddo
  endif


! ---- remove any isolated snow molecules (!) or sweep any excess snow from top of pack ----

  snow_lrunf  = 0.
  snow_frunf  = 0.
  snow_hlrunf = 0.
  snow_hfrunf = 0.
  if (0. < snow_mass .and. snow_mass < min_snow_mass ) then
        do l = 1, num_l
          snow_hlrunf = snow_hlrunf  &
            + clw*snow%wl(l)*(snow%T(l)-tfreeze)
          snow_hfrunf = snow_hfrunf  &
            + csw*snow%ws(l)*(snow%T(l)-tfreeze)
        enddo
        snow_lrunf  = sum(snow%wl)
        snow_frunf  = snow_mass
        snow_mass   = 0.
        snow%ws = 0.
        snow%wl = 0.
        if(is_watch_point()) then
           write(*,*) ' ***** snow_step_2 checkpoint 4e ***** '
           write(*,*) 'snow_lrunf    ', snow_lrunf
           write(*,*) 'snow_frunf    ', snow_frunf
        endif
  else if (max_snow < snow_mass) then
        snow_frunf  = snow_mass - max_snow
        snow_mass  = max_snow
        sum_sno  = 0
        snow_transfer = 0
        do l = 1, num_l
          if (sum_sno + snow%ws(l) > snow_frunf) then
              snow_transfer = snow_frunf - sum_sno
          else
              snow_transfer = snow%ws(l)
          endif
          if (snow%ws(l) > 0) then
              frac = snow_transfer / snow%ws(l)
          else
              frac = 1.
          endif
          sum_sno  = sum_sno  + snow_transfer
          snow_lrunf  = snow_lrunf  +     frac*snow%wl(l)
          snow_hlrunf = snow_hlrunf + clw*frac*snow%wl(l)*(snow%T(l)-tfreeze)
          snow_hfrunf = snow_hfrunf + csw*frac*snow%ws(l)*(snow%T(l)-tfreeze)
          snow%ws(l) = (1-frac)*snow%ws(l)
          snow%wl(l) = (1-frac)*snow%wl(l)
        enddo
        if(is_watch_point()) then
           write(*,*) ' ***** snow_step_2 checkpoint 4f ***** '
           write(*,*) 'snow_lrunf    ', snow_lrunf
           write(*,*) 'snow_frunf    ', snow_frunf
        endif
  endif
  snow_lrunf  = snow_lrunf  / delta_time
  snow_frunf  = snow_frunf  / delta_time
  snow_hlrunf = snow_hlrunf / delta_time
  snow_hfrunf = snow_hfrunf / delta_time

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 5 ***** '
     write(*,*) 'fict_heat         ', fict_heat
     do l = 1, num_l
        write(*,'(i2,3(a,g23.16))')l,&
             ' wl=', snow%wl(l),&
             ' ws=', snow%ws(l),&
             ' T =', snow%T(l)
     enddo
  endif

  snow_LMASS = 0; snow_FMASS = 0; snow_HEAT = 0
  do l = 1, num_l
        snow_LMASS = snow_LMASS + snow%wl(l)
        snow_FMASS = snow_FMASS + snow%ws(l)
        snow_HEAT = snow_HEAT + &
          (mc_fict*dz(l) + clw*snow%wl(l) + csw*snow%ws(l))  &
                                                * (snow%T(l)-tfreeze)
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 5.01 ***** '
     write(*,*) 'LMASS         ', snow_LMASS
     write(*,*) 'FMASS         ', snow_FMASS
     write(*,*) 'HEAT          ', snow_HEAT
  endif

  depth= 0.
  new_ws=0
  new_wl=0
  new_T=0
  do l = 1, num_l
    depth = depth + snow%ws(l)
  enddo
  depth = depth / snow_density

!************************** fudge to avoid T=NaN from too-small mass **
!   if(depth*snow_density < min_snow_mass .and. depth>0.) then
!       depth = 0
!       snow%ws = 0
!       snow%wl = 0
!     endif

! ---- re-layer the snowpack ------------------------------------------------
  do l = 1, num_l
     if (depth > 0) then
        new_ws(l) = snow_mass*dz(l)
        sum_sno = 0
        sum_liq = 0
        sum_heat = 0
     endif
     do l_old = 1, num_l
        if (depth > 0) then
           if (sum_sno + snow%ws(l_old) > new_ws(l)) then
              snow_transfer = new_ws(l) - sum_sno
           else
              snow_transfer = snow%ws(l_old)
           endif
           if (snow%ws(l_old) .ne. 0.) then
              frac = snow_transfer / snow%ws(l_old)
           else
              frac = 1
           endif
           sum_sno  = sum_sno  + snow_transfer
           sum_liq  = sum_liq  + frac*     snow%wl(l_old)
           sum_heat = sum_heat + frac*&
                (clw*snow%wl(l_old) + csw*snow%ws(l_old))&
                *snow%T(l_old)
           snow%ws(l_old) = (1.-frac)*snow%ws(l_old)
           snow%wl(l_old) = (1.-frac)*snow%wl(l_old)
           if(is_watch_point()) then
              write(*,*) 'l=',l, ' l_old=',l_old,snow_transfer,frac,&
                   sum_sno,sum_liq, sum_heat
           endif
        endif

     enddo
     if (depth > 0) then
        new_wl(l) = sum_liq
        new_T(l)  = sum_heat / (clw*new_wl(l) + csw*new_ws(l))
     endif
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 5.1 ***** '
     write(*,*) 'depth             ', depth
     write(*,*) 'fict_heat         ', fict_heat
     do l = 1, num_l
        write(*,'(i2,3(a,g23.16))')l,&
             ' new_wl=', new_wl(l),&
             ' new_ws=', new_ws(l),&
             ' new_T =', new_T(l)
     enddo
  endif

! add back fictional mass/heat
  do l = 1, num_l
    if (depth > 0) &
    new_T(l) = ( &
    (clw*new_wl(l) + csw*new_ws(l))*new_T(l)  &
      + mc_fict*dz(l)*fict_heat ) &
      / (clw*new_wl(l) + csw*new_ws(l) + dz(l)*mc_fict)
  enddo

  do l = 1, num_l
    if (depth > 0) then
      snow%ws(l) = new_ws(l)
      snow%wl(l) = new_wl(l)
      snow%T(l)  = new_T(l)
    endif
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 6 ***** '
     write(*,*) 'evap         ', subs_evap
     write(*,*) 'snow_lprec', snow_lprec
     write(*,*) 'depth        ', depth
     do l = 1, num_l
        write(*,'(i2,3(a,g23.16))')l,&
             ' wl=', snow%wl(l),&
             ' ws=', snow%ws(l),&
             ' T =', snow%T(l)
     enddo
  endif

  snow_LMASS = 0
  snow_FMASS = 0
  snow_HEAT = 0
  do l = 1, num_l
        snow_LMASS = snow_LMASS + snow%wl(l)
        snow_FMASS = snow_FMASS + snow%ws(l)
        snow_HEAT = snow_HEAT + &
          (mc_fict*dz(l) + clw*snow%wl(l) + csw*snow%ws(l))  &
                            * (snow%T(l)-tfreeze)
  enddo
  snow_Tbot = snow%T(num_l)
  snow_Cbot = mc_fict*dz(num_l) &
        + clw*snow%wl(num_l) + csw*snow%ws(num_l)
  snow_C = sum(mc_fict*dz(1:num_l) &
        + clw*snow%wl(1:num_l) + csw*snow%ws(1:num_l))
  snow_avrg_T = snow_HEAT/snow_C+tfreeze

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 7 ***** '
     write(*,*) 'LMASS         ', snow_LMASS
     write(*,*) 'FMASS         ', snow_FMASS
     write(*,*) 'HEAT          ', snow_HEAT
  endif

end subroutine snow_step_2

! ============================================================================
subroutine snow_sfc_water(snow, grnd_liq, grnd_ice, grnd_subl, grnd_tf)
  type(snow_tile_type), intent(in) :: snow
  real, intent(out) :: &
     grnd_liq, grnd_ice, & ! surface liquid and ice, respectively, kg/m2
     grnd_subl, &          ! fraction of vapor flux that sublimates
     grnd_tf               ! freezing temperature at the surface, degK

  grnd_liq = snow%wl(1) ! only liquid in the top snow layer is available to freeze implicitly
  grnd_ice = sum(snow%ws(:)) ! snow in any layer can be melted implicitly
  if (snow_active(snow)) then
     grnd_subl = 1.0
  else
     grnd_subl = 0
  endif
  grnd_tf = tfreeze
end subroutine snow_sfc_water

! ============================================================================
real function snow_subl_frac(snow)
  type(snow_tile_type), intent(in) :: snow

  snow_subl_frac = 0
  if (snow_active(snow)) snow_subl_frac = 1.0
end function snow_subl_frac

! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function snow_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   snow_tile_exists = associated(tile%snow)
end function snow_tile_exists

! ============================================================================
! accessor functions: given a pointer to a land tile, they return pointer
! to the desired member of the land tile, of NULL if this member does not
! exist.
subroutine snow_temp_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%T(i)
   endif
end subroutine snow_temp_ptr

subroutine snow_wl_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%wl(i)
   endif
end subroutine snow_wl_ptr

subroutine snow_ws_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr => tile%snow%ws(i)
   endif
end subroutine snow_ws_ptr

end module snow_mod



