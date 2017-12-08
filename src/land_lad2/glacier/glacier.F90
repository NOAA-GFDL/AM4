! ============================================================================
! glac model module
! ============================================================================
module glacier_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : error_mesg, file_exist, check_nml_error, stdlog, close_file, &
     mpp_pe, mpp_root_pe, FATAL, NOTE

use time_manager_mod,   only: time_type_to_real
use diag_manager_mod,   only: diag_axis_init
use constants_mod,      only: tfreeze, hlv, hlf, dens_h2o

use glac_tile_mod,      only: glac_tile_type, &
     read_glac_data_namelist, glac_data_thermodynamics, glac_data_hydraulics, &
     max_lev, cpw, clw, csw, use_brdf

use land_tile_mod, only : land_tile_map, land_tile_type, land_tile_enum_type, &
     first_elmt, loop_over_tiles
use land_tile_diag_mod, only : register_tiled_diag_field, send_tile_data, diag_buff_type, &
     set_default_diag_filter
use land_data_mod, only : lnd, log_version
use land_tile_io_mod, only: land_restart_type, &
     init_land_restart, open_land_restart, save_land_restart, free_land_restart, &
     add_restart_axis, add_tile_data, get_tile_data
use land_debug_mod, only : is_watch_point

implicit none
private

! ==== public interfaces =====================================================
public :: read_glac_namelist
public :: glac_init
public :: glac_end
public :: save_glac_restart
public :: glac_sfc_water
public :: glac_step_1
public :: glac_step_2
! =====end of public interfaces ==============================================


! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'glacier'
#include "../shared/version_variable.inc"

! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
logical :: lm2                   = .true.  ! *** CODE WORKS ONLY FOR .TRUE. !!! ****
logical :: conserve_glacier_mass = .true.
character(len=16):: albedo_to_use = ''  ! or 'brdf-params'
real    :: init_temp            = 260.       ! cold-start glac T
real    :: init_w               = 150.       ! cold-start w(l)/dz(l)
namelist /glac_nml/ lm2, conserve_glacier_mass,  albedo_to_use, &
                    init_temp, init_w, cpw, clw, csw
!---- end of namelist --------------------------------------------------------

logical         :: module_is_initialized =.FALSE.
real            :: delta_time       ! fast time step

integer         :: num_l            ! # of water layers
real            :: dz    (max_lev)  ! thicknesses of layers
real            :: zfull (max_lev)
real            :: zhalf (max_lev+1)


! ---- diagnostic field IDs
integer :: id_zhalf, id_zfull
integer :: id_lwc, id_swc, id_temp, id_ie, id_sn, id_bf, id_hie, id_hsn, id_hbf

! ==== end of module variables ===============================================

contains

! ============================================================================
subroutine read_glac_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: l            ! level iterator

  call read_glac_data_namelist(num_l, dz)

  call log_version(version, module_name, &
  __FILE__)
#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=glac_nml, iostat=io)
     ierr = check_nml_error(io, 'glac_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=glac_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'glac_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=glac_nml)
  endif

  ! ---- set up vertical discretization
  zhalf(1) = 0
  do l = 1, num_l;
     zhalf(l+1) = zhalf(l) + dz(l)
     zfull(l) = 0.5*(zhalf(l+1) + zhalf(l))
  enddo

end subroutine read_glac_namelist


! ============================================================================
! initialize glacier model
subroutine glac_init (id_ug)
  integer, intent(in)  :: id_ug !<Unstructured axis id.

  ! ---- local vars
  type(land_tile_enum_type)     :: ce   ! tile list enumerator
  type(land_tile_type), pointer :: tile ! pointer to current tile
  type(land_restart_type) :: restart
  logical :: restart_exists
  character(*), parameter :: restart_file_name='INPUT/glac.res.nc'

  module_is_initialized = .TRUE.
  delta_time = time_type_to_real(lnd%dt_fast)

  ! -------- initialize glac state --------
  call open_land_restart(restart,restart_file_name,restart_exists)
  if (restart_exists) then
     call error_mesg('glac_init',&
          'reading NetCDF restart "'//trim(restart_file_name)//'"',&
          NOTE)
     call get_tile_data(restart, 'temp', 'zfull', glac_temp_ptr)
     call get_tile_data(restart, 'wl',   'zfull', glac_wl_ptr)
     call get_tile_data(restart, 'ws',   'zfull', glac_ws_ptr)
  else
     call error_mesg('glac_init',&
          'cold-starting glacier',&
          NOTE)
     ce = first_elmt(land_tile_map)
     do while(loop_over_tiles(ce,tile))
        if (.not.associated(tile%glac)) cycle

        if (init_temp.ge.tfreeze.or.lm2) then      ! USE glac TFREEZE HERE
           tile%glac%wl(1:num_l) = init_w*dz(1:num_l)
           tile%glac%ws(1:num_l) = 0
        else
           tile%glac%wl(1:num_l) = 0
           tile%glac%ws(1:num_l) = init_w*dz(1:num_l)
        endif
        tile%glac%T             = init_temp
     enddo
  endif
  call free_land_restart(restart)

  if (trim(albedo_to_use)=='brdf-params') then
     use_brdf = .true.
  else if (trim(albedo_to_use)=='') then
     use_brdf = .false.
  else
     call error_mesg('glac_init',&
          'option albedo_to_use="'// trim(albedo_to_use)//&
          '" is invalid, use "brdf-params", or nothing ("")',&
          FATAL)
  endif

  if (.not.lm2) then
     call error_mesg('glac_init',&
          'currently only lm2=.TRUE. is supported',&
          FATAL)
  endif

  call glac_diag_init (id_ug, zfull(1:num_l), zhalf(1:num_l+1) )

end subroutine glac_init


! ============================================================================
subroutine glac_end ()

  module_is_initialized =.FALSE.

end subroutine glac_end

! ============================================================================
subroutine save_glac_restart (tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars
  character(267) :: filename
  type(land_restart_type) :: restart ! restart file i/o object

  call error_mesg('glac_end','writing NetCDF restart',NOTE)
! must set domain so that io_domain is available
! Note that filename is updated for tile & rank numbers during file creation
  filename = trim(timestamp)//'glac.res.nc'
  call init_land_restart(restart, filename, glac_tile_exists, tile_dim_length)
  call add_restart_axis(restart,'zfull',zfull(1:num_l),'Z','m','full level',sense=-1)

  ! Output data provides signature
  call add_tile_data(restart,'temp', 'zfull', glac_temp_ptr, longname='glacier temperature',  units='degrees_K')
  call add_tile_data(restart,'wl',   'zfull', glac_wl_ptr,   longname='liquid water content', units='kg/m2')
  call add_tile_data(restart,'ws',   'zfull', glac_ws_ptr,   longname='solid water content',  units='kg/m2')

  ! save performs io domain aggregation through mpp_io as with regular domain data
  call save_land_restart(restart)
  call free_land_restart(restart)
end subroutine save_glac_restart

! ============================================================================
subroutine glac_sfc_water(glac, grnd_liq, grnd_ice, grnd_subl, grnd_tf)
  type(glac_tile_type), intent(in)  :: glac
  real, intent(out) :: &
     grnd_liq, grnd_ice, & ! surface liquid and ice, respectively, kg/m2
     grnd_subl, &          ! fraction of vapor flux that sublimates
     grnd_tf               ! freezing temperature at the surface, degK

  if (lm2) then
     grnd_liq = 0
     grnd_ice = 1.e6
  else
     grnd_liq  = max(glac%wl(1), 0.0)
     grnd_ice  = max(glac%ws(1), 0.0)
  endif
  if (grnd_liq + grnd_ice > 0 ) then
     grnd_subl = grnd_ice / (grnd_liq + grnd_ice)
  else
     grnd_subl = 0
  endif

  grnd_tf = glac%pars%tfreeze
end subroutine glac_sfc_water

! ============================================================================
function glac_subl_frac(glac); real glac_subl_frac
  type(glac_tile_type), intent(in)  :: glac

  real :: grnd_liq, grnd_ice, grnd_subl, grnd_tf
  call glac_sfc_water(glac, grnd_liq, grnd_ice, grnd_subl, grnd_tf)
  glac_subl_frac = grnd_subl
end function glac_subl_frac

! ============================================================================
! update glac properties explicitly for time step.
! integrate glac-heat conduction equation upward from bottom of glac
! to surface, delivering linearization of surface ground heat flux.
subroutine glac_step_1 ( glac, &
                         glac_rh, glac_G0, &
                         glac_DGDT, conserve_glacier_mass_out )
  type(glac_tile_type),intent(inout) :: glac
  real, intent(out) :: &
       glac_rh, &
       glac_G0, &
       glac_DGDT
  logical, intent(out) :: conserve_glacier_mass_out

  ! ---- local vars
  real                   :: bbb, denom, dt_e
  real, dimension(num_l) :: aaa, ccc, thermal_cond, heat_capacity, vlc, vsc
  integer :: l

! ----------------------------------------------------------------------------
! in preparation for implicit energy balance, determine various measures
! of water availability, so that vapor fluxes will not exceed mass limits
! ----------------------------------------------------------------------------

  conserve_glacier_mass_out = conserve_glacier_mass

  if(is_watch_point()) then
    write(*,*) 'checkpoint gs1 a'
    write(*,*) 'mask    ',  .TRUE.
    write(*,*) 'T       ', glac%T(1)
  endif

  do l = 1, num_l
     vlc(l) = max(0.0, glac%wl(l) / (dens_h2o * dz(l)))
     vsc(l) = max(0.0, glac%ws(l) / (dens_h2o * dz(l)))
  enddo

  call glac_data_thermodynamics ( glac%pars, vlc(1), vsc(1),  &
       glac_rh, glac%heat_capacity_dry, thermal_cond )

  do l = 1, num_l
     heat_capacity(l) = glac%heat_capacity_dry(l)*dz(l) &
          + clw*glac%wl(l) + csw*glac%ws(l)
  enddo

  if(num_l > 1) then
     do l = 1, num_l-1
        dt_e = 2 / ( dz(l+1)/thermal_cond(l+1) &
                + dz(l)/thermal_cond(l)   )
        aaa(l+1) = - dt_e * delta_time / heat_capacity(l+1)
        ccc(l)   = - dt_e * delta_time / heat_capacity(l)
     enddo

     bbb = 1.0 - aaa(num_l)
     denom = bbb
     dt_e = aaa(num_l)*(glac%T(num_l) - glac%T(num_l-1)) &
               + glac%geothermal_heat_flux * delta_time / heat_capacity(num_l)
     glac%e(num_l-1) = -aaa(num_l)/denom
     glac%f(num_l-1) = dt_e/denom
     do l = num_l-1, 2, -1
        bbb = 1.0 - aaa(l) - ccc(l)
        denom = bbb + ccc(l)*glac%e(l)
        dt_e = - ( ccc(l)*(glac%T(l+1) - glac%T(l)  ) &
             -aaa(l)*(glac%T(l)   - glac%T(l-1)) )
        glac%e(l-1) = -aaa(l)/denom
        glac%f(l-1) = (dt_e - ccc(l)*glac%f(l))/denom
     enddo
     denom = delta_time/(heat_capacity(1) )
     glac_G0    = ccc(1)*(glac%T(2)- glac%T(1) &
          + glac%f(1)) / denom
     glac_DGDT  = (1 - ccc(1)*(1-glac%e(1))) / denom
  else  ! one-level case
     denom = delta_time/heat_capacity(1)
     glac_G0    = 0.
     glac_DGDT  = 1. / denom
  endif

  if(is_watch_point())then
     write(*,*) 'checkpoint gs1 c'
     write(*,*) 'mask    ', .TRUE.
     write(*,*) 'rh      ', glac_rh
     write(*,*) 'G0      ', glac_G0
     write(*,*) 'DGDT    ', glac_DGDT
     do l = 1, num_l
        write(*,*) 'T(dbg,l)', glac%T(l)
     enddo

  endif
end subroutine glac_step_1


! ============================================================================
! apply boundary flows to glac water and move glac water vertically.
  subroutine glac_step_2 ( glac, diag, snow_lprec, snow_hlprec,  &
                           subs_DT, subs_M_imp, subs_evap, &
                           glac_levap, glac_fevap, glac_melt, &
                           glac_lrunf, glac_hlrunf, glac_Ttop, glac_Ctop )
! *** WARNING!!! MOST OF THIS CODE IS SIMPLY COPIED FROM SOIL FOR POSSIBLE
! FUTURE DEVELOPMENT. (AND SIMILAR CODE IN SOIL MOD HAS BEEN FURTHER DEVELOPED,
! SO THIS IS MAINLY JUNK.) ONLY THE LM2 BRANCHES WORK. FOR LM2, THE SURFACE OF THE
! GLACIER IS EFFECTIVELY SEALED W.R.T. MASS TRANSER:
! NO LIQUID CAN INFILTRATE, AND NO GLACIER MASS
! CAN ENTER THE ATMOSPHERE. ONLY SUPERFICIAL SNOW PARTICIPATES IN THE WATER
! CYCLE. HOWEVER, SENSIBLE HEAT TRANSFER AND GLACIER MELT CAN OCCUR,
! TO AN UNLIMITED EXTENT, AS NEEDED
! TO KEEP GLACIER AT OR BELOW FREEZING. MELT WATER STAYS IN GLACIER.

  type(glac_tile_type), intent(inout) :: glac
  type(diag_buff_type), intent(inout) :: diag
  real, intent(in) :: &
     snow_lprec, &
     snow_hlprec, &
     subs_DT,       &!
     subs_M_imp,       &! rate of phase change of non-evaporated glac water
     subs_evap
  real, intent(out) :: &
     glac_levap, glac_fevap, glac_melt, &
     glac_lrunf, glac_hlrunf, glac_Ttop, glac_Ctop

  ! ---- local vars ----------------------------------------------------------
  real, dimension(num_l) :: del_t, eee, fff, &
             psi, DThDP, hyd_cond, DKDP, K, DKDPm, DKDPp, grad, &
             vlc, vsc, dW_l, u_minus, u_plus, DPsi, glac_w_fc
  real, dimension(num_l+1) :: flow
  real, dimension(num_l  ) :: div
  real :: &
     lprec_eff, hlprec_eff, tflow, hcap,cap_flow, &
     melt_per_deg, melt, glac_subl, &
     lrunf_sn,lrunf_ie,lrunf_bf, hlrunf_sn,hlrunf_ie,hlrunf_bf, &
     Qout, DQoutDP,&
     tau_gw, c0, c1, c2, x, aaa, bbb, ccc, ddd, xxx, Dpsi_min, Dpsi_max
  logical :: stiff
  real, dimension(num_l-1) :: del_z
  integer :: l
  real :: jj
  ! --------------------------------------------------------------------------

  jj = 1.
  DPsi = 0.0
  c1   = 0.0

  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 1 ***** '
     write(*,*) 'mask    ', .TRUE.
     write(*,*) 'subs_evap    ', subs_evap
     write(*,*) 'snow_lprec   ', snow_lprec
     write(*,*) 'subs_M_imp   ', subs_M_imp
     write(*,*) 'theta_s ', glac%pars%w_sat
     do l = 1, num_l
        write(*,'(i2.2,99(a,g23.16))')l,&
             ' T =', glac%T(l),&
             ' Th=', (glac%ws(l)+glac%wl(l))/(dens_h2o*dz(l)),&
             ' wl=', glac%wl(l),&
             ' ws=', glac%ws(l)
     enddo

  endif

  ! ---- record fluxes ---------
  glac_subl = glac_subl_frac(glac)
  IF (LM2) THEN ! EVAP SHOULD BE ZERO ANYWAY, BUT THIS IS JUST TO BE SURE...
     glac_levap  = 0.
     glac_fevap  = 0.
  ELSE
     glac_levap  = subs_evap*(1-glac_subl)
     glac_fevap  = subs_evap*   glac_subl
  ENDIF
  glac_melt   = subs_M_imp / delta_time

  ! ---- load surface temp change and perform back substitution --------------
  del_t(1) = subs_DT
  glac%T(1) = glac%T(1) + del_t(1)
  if ( num_l > 1) then
    do l = 1, num_l-1
      del_t(l+1) = glac%e(l) * del_t(l) + glac%f(l)
      glac%T(l+1) = glac%T(l+1) + del_t(l+1)
    end do
  end if

  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 2 ***** '
     write(*,*) 'levap=',glac_levap
     write(*,*) 'fevap=',glac_fevap
     write(*,*) 'subs_M_imp=',subs_M_imp
     do l = 1, num_l
        write(*,'(i2.2,x,a,g23.16)') l, 'T', glac%T(l)
     enddo
  endif

IF (LM2) THEN ! *********************************************************
    glac_lrunf  = snow_lprec
    glac_hlrunf = snow_hlprec
ELSE   ! ****************************************************************
  ! ---- extract evap from glac and do implicit melt --------------------
    glac%wl(1) = glac%wl(1) - glac_levap*delta_time
    glac%ws(1) = glac%ws(1) - glac_fevap*delta_time
    hcap = glac%heat_capacity_dry(1)*dz(1) &
                       + clw*glac%wl(1) + csw*glac%ws(1)
    glac%T(1) = glac%T(1) + (   &
                  +((clw-cpw)*glac_levap                              &
                  + (csw-cpw)*glac_fevap)*(glac%T(1)  -tfreeze) &
                                               )*delta_time/ hcap
    glac%wl(1) = glac%wl(1) + subs_M_imp
    glac%ws(1) = glac%ws(1) - subs_M_imp
    glac%T(1)  = tfreeze + (hcap*(glac%T(1)-tfreeze) ) &
                              / ( hcap + (clw-csw)*subs_M_imp )
  ! ---- remainder of mass fluxes and associated sensible heat fluxes --------
  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 3 ***** '
     do l = 1, num_l
        write(*,'(i2.2,99(a,g23.16))') l,&
             ' T =', glac%T(l),&
             ' wl=', glac%wl(l),&
             ' ws=', glac%ws(l)
     enddo
  endif

  ! ---- fetch glac hydraulic properties -------------------------------------
  vlc=1;vsc=0
  do l = 1, num_l
     vlc(l) = max(0., glac%wl(l) / (dens_h2o*dz(l)))
     vsc(l) = max(0., glac%ws(l) / (dens_h2o*dz(l)))
  enddo
  call glac_data_hydraulics (glac, vlc, vsc, &
                   psi, DThDP, hyd_cond, DKDP, Dpsi_min, Dpsi_max, tau_gw, &
                   glac_w_fc )
     if(is_watch_point()) then
        write(*,*) ' ***** glac_step_2 checkpoint 3.1 ***** '
        do l = 1, num_l
           write(*,'(i2.2,99(x,a,g23.16))') l, 'vlc', vlc(l),&
                'K  ', hyd_cond(l)
        enddo

     endif
    div = 0.
    do l = 1, num_l
      div(l) = 0.15*dens_h2o*dz(l)/tau_gw
    enddo
    lrunf_bf = sum(div)

  ! ---- glac-water flow ----------------------------------------------------
    stiff = all(DThDP.eq.0)
    if (snow_lprec/=0 .and. psi(num_l)>0) then
      lrunf_sn = snow_lprec*min((psi(num_l)/zhalf(num_l))**glac%pars%rsa_exp,1.)
      hlrunf_sn = lrunf_sn*snow_hlprec/snow_lprec
    else
      lrunf_sn = 0.
      hlrunf_sn = 0.
    endif
    lprec_eff = snow_lprec - lrunf_sn
    hlprec_eff = snow_hlprec - hlrunf_sn
    flow(1) = delta_time*lprec_eff
    do l = 1, num_l-1
      del_z(l) = zfull(l+1)-zfull(l)
      K(l) = 0.5*(hyd_cond(l)+hyd_cond(l+1))
      DKDPm(l) = 0. !0.5*DKDP(l)
      DKDPp(l) = 0. ! 0.5*DKDP(l+1)
!        K(l) = hyd_cond(l)
!        DKDPm(l) = DKDP(l)
!        DKDPp(l) = 0
      grad(l)  = jj*(psi(l+1)-psi(l))/del_z(l) - 1
    enddo


    if(is_watch_point()) then
       write(*,*) ' ***** glac_step_2 checkpoint 3.1 ***** '
       do l = 1, num_l
          write(*,'(i2.2,x,a,99g23.16)') l, 'DThDP,hyd_cond,psi,DKDP', &
               DThDP(l), hyd_cond(l), psi(l), DKDP(l)
       enddo
       do l = 1, num_l-1
          write(*,'(i2.2,x,a,99g23.16)') l, 'K,DKDPm,DKDPp,grad,del_z', &
               K(l), DKDPm(l), DKDPp(l), grad(l)
      enddo
    endif

    l = num_l
    xxx = dens_h2o*dz(l)*DThDP(l)/delta_time
    aaa =     - ( jj* K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
        bbb = xxx - (- jj*K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1) )
        ddd = - K(l-1) *grad(l-1) - div(l)
    eee(l-1) = -aaa/bbb
    fff(l-1) =  ddd/bbb

    if(is_watch_point()) then
       write(*,'(a,i4,99g23.16)') 'l,a,b, ,d', l,aaa, bbb,ddd
    endif


    do l = num_l-1, 2, -1
      xxx = dens_h2o*dz(l)*DThDP(l)/delta_time
      aaa = - ( jj*K(l-1)/del_z(l-1) - DKDPm(l-1)*grad(l-1))
      bbb = xxx-( -jj*K(l-1)/del_z(l-1) - DKDPp(l-1)*grad(l-1)&
                  -jj*K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
      ccc =   - (  jj*K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
      ddd =       K(l)*grad(l) - K(l-1)*grad(l-1) &
                            - div(l)
      eee(l-1) =                    -aaa/(bbb+ccc*eee(l))
      fff(l-1) =  (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
      if(is_watch_point()) then
         write(*,'(a,i4,99g23.16)') 'l,a,b,c,d', l,aaa, bbb,ccc,ddd
      endif
    enddo

    l = 1
    xxx = dens_h2o*dz(l)*DThDP(l)/delta_time
    bbb = xxx - ( -jj*K(l  )/del_z(l  ) + DKDPm(l  )*grad(l  ))
    ccc =     - (  jj*K(l  )/del_z(l  ) + DKDPp(l  )*grad(l  ))
    ddd =          flow(1)/delta_time +    K(l)     *grad(l) &
                            - div(l)
    if (stiff) then
      dPsi(l) =  - psi(l)
    else
      dPsi(l) = (ddd-ccc*fff(l))/(bbb+ccc*eee(l))
      dPsi(l) = min (dPsi(l), Dpsi_max)
      dPsi(l) = max (dPsi(l), Dpsi_min)
    endif
    flow(l) = (dPsi(l)*(bbb+ccc*eee(l))+ccc*fff(l) &
                      - K(l)*grad(l))*delta_time
    lrunf_ie         = lprec_eff - flow(l)/delta_time

    if(is_watch_point()) then
       write(*,'(a,i4,99g23.16)') 'l,  b,c,d', l, bbb,ccc,ddd

       write(*,*) ' ***** glac_step_2 checkpoint 3.2 ***** '
       write(*,*) 'ie,sn,bf:', lrunf_ie,lrunf_sn,lrunf_bf
       do l = 1, num_l-1
          write(*,'(a,i4,99g23.16)') 'l,eee(l),fff(l)', l,eee(l), fff(l)
       enddo
       write(*,*) 'DThDP(1)', DThDP(1)
       write(*,*) 'ddd(1)', ddd
       write(*,*) 'ccc(1)', ccc
       write(*,*) 'bbb(1)', bbb
       write(*,*) 'dPsi(1)', dPsi(1)
       write(*,*) 'Psi(1)', Psi(1)
    endif

    do l = 2, num_l
      dPsi(l) = eee(l-1)*dPsi(l-1) + fff(l-1)
    enddo

    do l = 1, num_l-1
      flow(l+1) = delta_time*( &
           -K(l)*(grad(l)&
           +jj*(DPsi(l+1)-DPsi(l))/ del_z(l)) &
           -grad(l)*(DKDPp(l)*Dpsi(l+1)+ &
                           DKDPm(l)*Dpsi(l) )  )
      dW_l(l) = flow(l) - flow(l+1) - div(l)*delta_time
      glac%wl(l) = glac%wl(l) + dW_l(l)
    enddo
    flow(num_l+1) = 0.
    dW_l(num_l) = flow(num_l) - flow(num_l+1) &
                          - div(num_l)*delta_time
    glac%wl(num_l) = glac%wl(num_l) + dW_l(num_l)

  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 3.3 ***** '
     write(*,*) 'psi_sat',glac%pars%psi_sat_ref
     write(*,*) 'Dpsi_max',Dpsi_max
     do l = 1, num_l
        write(*,'(i2.2,99(a,g23.16))')l,&
             ' Th=', (glac%ws(l)+glac%wl(l))/(dens_h2o*dz(l)),&
             ' wl=', glac%wl(l),&
             ' ws=', glac%ws(l),&
             'Dpsi=', dPsi(l), &
             'flow=', flow(l)
     enddo
  endif

  if  (snow_hlprec.ne.0.) then
    tflow = tfreeze + snow_hlprec/(clw*snow_lprec)
  else
    tflow = tfreeze
  endif

  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 3.4 ***** '
     write(*,*) ' tfreeze', tfreeze
     write(*,*) ' snow_hlprec', snow_hlprec
  endif

! For initial testing, use top-down-flow weights to advect heat.
  u_minus = 1.
  u_plus  = 0.
  if (flow(1).lt.0.) u_minus(1) = 0.
  hcap = (glac%heat_capacity_dry(num_l)*dz(num_l) &
                              + csw*glac%ws(num_l))/clw
  aaa = -flow(num_l) * u_minus(num_l)
  bbb =  hcap + glac%wl(num_l) - dW_l(num_l) - aaa
  eee(num_l-1) = -aaa/bbb
  fff(num_l-1) = aaa*(glac%T(num_l)-glac%T(num_l-1)) / bbb

  do l = num_l-1, 2, -1
    hcap = (glac%heat_capacity_dry(l)*dz(l) &
                              + csw*glac%ws(l))/clw
    aaa = -flow(l)   * u_minus(l)
    ccc =  flow(l+1) * u_plus (l)
    bbb =  hcap + glac%wl(l) - dW_l(l) - aaa - ccc
    eee(l-1) = -aaa / ( bbb +ccc*eee(l) )
    fff(l-1) = (   aaa*(glac%T(l)-glac%T(l-1))    &
                       + ccc*(glac%T(l)-glac%T(l+1))    &
                       - ccc*fff(l) ) / ( bbb +ccc*eee(l) )
  enddo

  hcap = (glac%heat_capacity_dry(1)*dz(1) + csw*glac%ws(1))/clw
  aaa = -flow(1) * u_minus(1)
  ccc =  flow(2) * u_plus (1)
  bbb =  hcap + glac%wl(1) - dW_l(1) - aaa - ccc

  del_t(1) =  (  aaa*(glac%T(1)-tflow          ) &
                     + ccc*(glac%T(1)-glac%T(2)) &
                     - ccc*fff(1) ) / (bbb+ccc*eee(1))
  glac%T(1) = glac%T(1) + del_t(1)

  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 3.4.1 ***** '
     write(*,*) 'hcap', hcap
     write(*,*) 'aaa', aaa
     write(*,*) 'bbb', bbb
     write(*,*) 'ccc', ccc
     write(*,*) 'del_t(1)', del_t(1)
     write(*,*) ' T(1)', glac%T(1)
  endif

  do l = 1, num_l-1
    del_t(l+1) = eee(l)*del_t(l) + fff(l)
    glac%T(l+1) = glac%T(l+1) + del_t(l+1)
  enddo

  tflow = glac%T(num_l)

  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 3.5 ***** '
     write(*,*) 'hcap', hcap
     write(*,*) 'cap_flow', cap_flow
     do l = 1, num_l
        write(*,'(i2.2,99(a,g23.16))')l, ' T', glac%T(l), ' flow ',flow(l)
     enddo
     write(*,*) 'delta_time,tau_gw,c0,c1,c2,x', delta_time,tau_gw,c0,&
          c1,c2,x
     write(*,*) 'level=', num_l+1, ' flow ',flow(num_l+1)
  endif

  ! THIS T AVERAGING IS WRONG, BECAUSE IT NEGLECTS THE MEDIUM  ***
  ! ALSO, FREEZE-THAW IS NEEDED!
  ! PROBABLY THIS SECTION WILL BE DELETED ANYWAY, WITH GW TREATED ABOVE.
    if (lprec_eff.ne.0. .and. flow(1).ge.0. ) then
      hlrunf_ie = lrunf_ie*hlprec_eff/lprec_eff
    else if (flow(1).lt.0. ) then
      hlrunf_ie = hlprec_eff - (flow(1)/delta_time)*clw &
                         *(glac%T(1)-tfreeze)
    else
      hlrunf_ie = 0.
    endif
    hlrunf_bf = clw*sum(div*(glac%T-tfreeze))
    glac_lrunf  = lrunf_sn + lrunf_ie + lrunf_bf
    glac_hlrunf = hlrunf_sn + hlrunf_bf + hlrunf_ie
  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 3.7 ***** '
     write(*,*) 'hcap', hcap
     write(*,*) 'cap_flow', cap_flow
     do l = 1, num_l
        write(*,'(i2.2,99(a,g23.16))')l, ' T', glac%T(l)
     enddo
  endif
    do l = 1, num_l
    ! ---- compute explicit melt/freeze --------------------------------------
      hcap = glac%heat_capacity_dry(l)*dz(l) &
               + clw*glac%wl(l) + csw*glac%ws(l)
      melt_per_deg = hcap/hlf
      if (glac%ws(l)>0 .and. glac%T(l)>glac%pars%tfreeze) then
        melt =  min(glac%ws(l), (glac%T(l)-glac%pars%tfreeze)*melt_per_deg)
      else if (glac%wl(l)>0 .and. glac%T(l)<glac%pars%tfreeze) then
        melt = -min(glac%wl(l), (glac%pars%tfreeze-glac%T(l))*melt_per_deg)
      else
        melt = 0
      endif

      if(is_watch_point()) then
         write(*,'(a,i4,99g23.16)') 'l,T,wl(1),ws(1),melt:', l,glac%T(l), glac%wl(l), &
              glac%ws(l), melt
      endif

      glac%wl(l) = glac%wl(l) + melt
      glac%ws(l) = glac%ws(l) - melt
      glac%T(l) = tfreeze &
         + (hcap*(glac%T(l)-tfreeze) - hlf*melt) &
                              / ( hcap + (clw-csw)*melt )
      if(is_watch_point()) then
         write(*,'(a,i4,99g23.16)') 'l,T,wl(1),ws(1):', l,glac%T(l), glac%wl(l), &
              glac%ws(l)
      endif

      glac_melt = glac_melt + melt / delta_time
    enddo
  if(is_watch_point()) then
     write(*,*) ' ***** glac_step_2 checkpoint 5 ***** '
     write(*,*) 'i,j,k,melt:',&
          glac_melt*delta_time
     do l = 1, num_l
        write(*,'(i2.2,99(a,g23.16))')l, &
             ' T =', glac%T(l), &
             ' Th=', (glac%ws(l)+glac%wl(l))/(dens_h2o*dz(l)),&
             ' wl=', glac%wl(l),&
             ' ws=', glac%ws(l)
     enddo
  endif

ENDIF  !*****************************************************************************

  glac_Ttop = glac%T(1)
  glac_Ctop = glac%heat_capacity_dry(1)*dz(1) &
    + clw*glac%wl(1) + csw*glac%ws(1)

! ----------------------------------------------------------------------------
! given solution for surface energy balance, write diagnostic output.
!

  ! ---- diagnostic section
  call send_tile_data (id_temp, glac%T,     diag )
  call send_tile_data (id_lwc,  glac%wl(1:num_l)/dz(1:num_l), diag )
  call send_tile_data (id_swc,  glac%ws(1:num_l)/dz(1:num_l), diag )
  if (.not.lm2) then
     call send_tile_data (id_ie,   lrunf_ie,        diag )
     call send_tile_data (id_sn,   lrunf_sn,        diag )
     call send_tile_data (id_bf,   lrunf_bf,        diag )
     call send_tile_data (id_hie,  hlrunf_ie,       diag )
     call send_tile_data (id_hsn,  hlrunf_sn,       diag )
     call send_tile_data (id_hbf,  hlrunf_bf,       diag )
  endif

end subroutine glac_step_2

! ============================================================================
subroutine glac_diag_init (id_ug, zfull, zhalf )
  integer,         intent(in) :: id_ug   !<Unstructured axis id.
  real,            intent(in) :: zfull(:)! Full levels, m
  real,            intent(in) :: zhalf(:)! Half levels, m

  ! ---- local vars ----------------------------------------------------------
  integer :: axes(2)

  ! define vertical axis
  id_zhalf = diag_axis_init ( &
       'glac_zhalf', zhalf, 'meters', 'z', 'half level',  -1, set_name='glac' )
  id_zfull = diag_axis_init ( &
       'glac_zfull', zfull, 'meters', 'z', 'full level',  -1, set_name='glac', &
       edges=id_zhalf )

  ! define array of axis indices
  axes = (/id_ug,id_zfull/)

  ! set the default sub-sampling filter for the fields below
  call set_default_diag_filter('glac')

  ! define diagnostic fields
  id_lwc = register_tiled_diag_field ( module_name, 'glac_liq', axes,        &
       lnd%time, 'bulk density of liquid water', 'kg/m3', missing_value=-100.0 )
  id_swc  = register_tiled_diag_field ( module_name, 'glac_ice',  axes,      &
       lnd%time, 'bulk density of solid water', 'kg/m3',  missing_value=-100.0 )
  id_temp  = register_tiled_diag_field ( module_name, 'glac_T',  axes,       &
       lnd%time, 'temperature',            'degK',  missing_value=-100.0 )
  if (.not.lm2) then
     id_ie  = register_tiled_diag_field ( module_name, 'glac_rie',  axes(1:1),  &
          lnd%time, 'inf exc runf',            'kg/(m2 s)',  missing_value=-100.0 )
     id_sn  = register_tiled_diag_field ( module_name, 'glac_rsn',  axes(1:1),  &
          lnd%time, 'satn runf',            'kg/(m2 s)',  missing_value=-100.0 )
     id_bf  = register_tiled_diag_field ( module_name, 'glac_rbf',  axes(1:1),  &
          lnd%time, 'baseflow',            'kg/(m2 s)',  missing_value=-100.0 )
     id_hie  = register_tiled_diag_field ( module_name, 'glac_hie',  axes(1:1), &
          lnd%time, 'heat ie runf',            'W/m2',  missing_value=-100.0 )
     id_hsn  = register_tiled_diag_field ( module_name, 'glac_hsn',  axes(1:1), &
          lnd%time, 'heat sn runf',            'W/m2',  missing_value=-100.0 )
     id_hbf  = register_tiled_diag_field ( module_name, 'glac_hbf',  axes(1:1), &
          lnd%time, 'heat bf runf',            'W/m2',  missing_value=-100.0 )
  endif

end subroutine glac_diag_init

! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function glac_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   glac_tile_exists = associated(tile%glac)
end function glac_tile_exists


! ============================================================================
! accessor functions: given a pointer to a land tile, they return pointer
! to the desired member of the land tile, of NULL if this member does not
! exist.
subroutine glac_temp_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%glac)) ptr => tile%glac%T(i)
   endif
end subroutine glac_temp_ptr

subroutine glac_wl_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%glac)) ptr => tile%glac%wl(i)
   endif
end subroutine glac_wl_ptr

subroutine glac_ws_ptr(tile, i, ptr)
   type(land_tile_type), pointer :: tile ! input
   integer             , intent(in) :: i ! index in the array
   real                , pointer :: ptr  ! returned pointer to the data
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%glac)) ptr => tile%glac%ws(i)
   endif
end subroutine glac_ws_ptr

end module glacier_mod
