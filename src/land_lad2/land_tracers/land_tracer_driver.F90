module land_tracer_driver_mod

#include "../shared/debug.inc"

use constants_mod, only : rdgas
use time_manager_mod, only : time_type, time_type_to_real
use fms_mod, only : lowercase, stdout, stdlog
use field_manager_mod , only : MODEL_ATMOS, MODEL_LAND, parse
use tracer_manager_mod, only : NO_TRACER, get_tracer_index, get_tracer_names, query_method
use table_printer_mod

use land_constants_mod, only : d608, kBoltz
use land_debug_mod, only : is_watch_point
use land_data_mod, only : lnd, log_version
use land_tracers_mod, only : ntcana, isphum, ico2
use land_tile_mod, only : land_tile_type, land_tile_grnd_T
use land_tile_diag_mod, only : set_default_diag_filter, &
        register_tiled_diag_field, send_tile_data

use cana_tile_mod, only : canopy_air_mass_for_tracers
use vegn_data_mod, only : spdata
use vegn_tile_mod, only : vegn_tile_LAI
use vegn_cohort_mod, only : get_vegn_wet_frac

! import interfaces from non-generic tracer modules, e.g.:
! use land_dust_mod, only : land_dust_init, land_dust_end, update_land_dust

implicit none
private

! ==== public interfaces =====================================================
public :: land_tracer_driver_init
public :: land_tracer_driver_end
public :: update_cana_tracers
! ==== end of public interfaces ==============================================

! ---- module constants ------------------------------------------------------
character(len=*), parameter :: module_name = 'land_tracer_driver_mod'
#include "../shared/version_variable.inc"
character(len=*), parameter :: diag_name   = 'land_tracers'

real :: diffusivity_h2o = 0.282e-4 ! diffusivity of water vapor m2/s,
    ! Cussler, E. L. (1997). Diffusion: Mass Transfer in Fluid Systems (2nd ed.).
    ! New York: Cambridge University Press. ISBN 0-521-45078-0.

! ---- data types -----------------------------------------------------------
type :: tracer_data_type
   character(32) :: name = ''  ! tracer name
   integer :: tr_atm  = NO_TRACER ! index of this tracer in atmos tracer array

   logical :: is_generic    = .TRUE. ! flag of generic tracer; initialization of non-generic tracers should turn it to FALSE
   logical :: do_deposition = .TRUE. ! if true, generic dry deposition is used
   ! dry deposition parameters. The default values are set as O3 parameters from (Wesely, 1989)
   real    :: diff_ratio    = 1.6    ! ratio of water vapor molecular diffusivity in the air to that of the tracer, unitless
   real    :: Henry_const   = 0.01   ! Henry's law constant for the gas, M atm-1
   real    :: reactivity    = 1.0    ! normalized reactivity factor
   real    :: r_lw          = 1000.0 ! resistance of wet portion of the leaf, s/m
   real    :: r_ls          = 1000.0 ! resistance of snow covered portion of the leaf, s/m
   real    :: r_gsS         = 100.0  ! parameter of ground resistance
   real    :: r_gsO         = 300.0  ! parameter of ground resistance

   integer :: & ! diag field IDs
     id_emis,      id_ddep,  &
     id_flux_atm,  id_dfdtr, &
     id_con_v_lam, id_con_g_lam, &
     id_con_v,     id_con_g, &
     id_conc
end type tracer_data_type

! ---- private module variables ----------------------------------------------
logical :: module_is_initialized = .FALSE.
real, save :: dt ! fast time step, s
type(tracer_data_type), allocatable :: trdata(:)

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
subroutine land_tracer_driver_init(id_ug)
  integer,intent(in) :: id_ug !<Unstructured axis id.
!----------

  integer :: tr ! tracer index
  real    :: value ! temporary storage for parsing input
  character(32)  :: name, units, funits ! name and units of the tracer and flux
  character(128) :: longname ! long name of the tracer
  character(32)  :: method
  character(512) :: parameters
  type(table_printer_type) :: table

  ! write the version and tag name to the logfile
  call log_version(version, module_name, &
  __FILE__)

  ! calculate time step
  dt  = time_type_to_real(lnd%dt_fast) ! store in a module variable for convenience

  ! allocate storage for tracer informations
  allocate(trdata(ntcana))
  trdata(isphum)%is_generic = .FALSE.
  trdata(ico2)%is_generic   = .FALSE.

  ! initialize non-generic tracers, e.g.:
  ! call land_dust_init(id_ug,trdata(:)%is_generic)

  ! NOTE that (1) the non-generic tracer init must skip all non-generic tracers
  ! that have already been initialized (in case there is a conflict), and
  ! (2) it must set trdata(:)%is_generic to FALSE for the tracers it claims

  ! initialize generic tracer parameters
  do tr = 1, ntcana
     if (.not.trdata(tr)%is_generic) cycle ! skip all non-generic tracers

     call get_tracer_names(MODEL_LAND, tr, trdata(tr)%name)
     trdata(tr)%tr_atm = get_tracer_index (MODEL_ATMOS, trdata(tr)%name)

     ! set up deposition flag
     trdata(tr)%do_deposition = .FALSE.
     method = ''; parameters = ''
     if (trdata(tr)%tr_atm >0) then
        if (query_method('dry_deposition', MODEL_ATMOS, trdata(tr)%tr_atm, method)) then
           trdata(tr)%do_deposition = (index(lowercase(method),'land:lm3')>0)
        endif
     endif
     ! set up deposition parameters
     if(query_method('dry_deposition', MODEL_LAND, tr, method, parameters)) then
        if ( parse(parameters, 'diff_ratio',  value) > 0 ) trdata(tr)%diff_ratio  = value
        if ( parse(parameters, 'Henry_const', value) > 0 ) trdata(tr)%Henry_const = value
        if ( parse(parameters, 'reactivity',  value) > 0 ) trdata(tr)%reactivity  = value
        if ( parse(parameters, 'r_lw',        value) > 0 ) trdata(tr)%r_lw        = value
        if ( parse(parameters, 'r_ls',        value) > 0 ) trdata(tr)%r_ls        = value
        if ( parse(parameters, 'r_gsS',       value) > 0 ) trdata(tr)%r_gsS       = value
        if ( parse(parameters, 'r_gsO',       value) > 0 ) trdata(tr)%r_gsO       = value
     endif
  enddo

  call init_with_headers(table, trdata(:)%name)
  call add_row(table, 'do_deposition',    trdata(:)%do_deposition)
  call add_row(table, 'diff_ratio',       trdata(:)%diff_ratio)
  call add_row(table, 'Henry_const',      trdata(:)%Henry_const)
  call add_row(table, 'reactivity',       trdata(:)%reactivity)
  call add_row(table, 'r_lw',             trdata(:)%r_lw)
  call add_row(table, 'r_ls',             trdata(:)%r_ls)
  call add_row(table, 'r_gsS',            trdata(:)%r_gsS)
  call add_row(table, 'r_gsO',            trdata(:)%r_gsO)

  call print(table,stdlog())
  call print(table,stdout())

  ! register diag fields for generic tracers
  call set_default_diag_filter('land')
  do tr = 1, ntcana
     call get_tracer_names(MODEL_LAND, tr, name, longname, units)

     funits = flux_units(units)
     trdata(tr)%id_flux_atm = &
       register_tiled_diag_field(diag_name, trim(name)//'_flux_atm', &
       (/id_ug/),  lnd%time, trim(name)//' flux to the atmosphere', &
       trim(funits), missing_value=-1.0)
     ! TODO: verify units of dfdtr
     trdata(tr)%id_dfdtr = &
       register_tiled_diag_field(diag_name, trim(name)//'_dfdtr', &
       (/id_ug/),  lnd%time,'derivative of '//trim(name)//' flux to the atmosphere', &
       trim(funits), missing_value=-1.0)
     if (trdata(tr)%do_deposition) then
        ! TODO: initialize parameters of generic dry deposition here

        trdata(tr)%id_ddep = &
          register_tiled_diag_field(diag_name, trim(name)//'_ddep', &
          (/id_ug/),  lnd%time, trim(name)//' dry deposition', 'kg/(m2 s)', &
          missing_value=-1.0)
        trdata(tr)%id_con_v_lam = &
          register_tiled_diag_field(diag_name, trim(name)//'_con_v_lam', &
          (/id_ug/),  lnd%time, 'quasi-laminar conductance between canopy and canopy air for '//trim(name), &
          'm/s', missing_value=-1.0)
        trdata(tr)%id_con_g_lam = &
          register_tiled_diag_field(diag_name, trim(name)//'_con_g_lam', &
          (/id_ug/),  lnd%time, 'quasi-laminar conductance between ground and canopy air for '//trim(name), &
          'm/s', missing_value=-1.0)
        trdata(tr)%id_con_v = &
          register_tiled_diag_field(diag_name, trim(name)//'_con_v', &
          (/id_ug/),  lnd%time, 'total conductance between canopy and canopy air for'//trim(name), &
          'm/s', missing_value=-1.0)
        trdata(tr)%id_con_g = &
          register_tiled_diag_field(diag_name, trim(name)//'_con_g', &
          (/id_ug/),  lnd%time, 'total conductance between ground and canopy air for '//trim(name), &
          'm/s', missing_value=-1.0)
        trdata(tr)%id_conc = &
          register_tiled_diag_field(diag_name, trim(name), &
          (/id_ug/),  lnd%time, 'concentration or '//trim(name)//' in canopy air', &
          units, missing_value=-1.0)
     endif
  enddo
  module_is_initialized = .TRUE.
end subroutine land_tracer_driver_init


! ============================================================================
subroutine land_tracer_driver_end()
  ! call finalizers for all non-generic tracers, e.g.:
  ! call land_dust_end()

  deallocate(trdata)
  module_is_initialized = .FALSE.
end subroutine land_tracer_driver_end


! ============================================================================
! calculates conductance of laminar sublayer for given tracer.
! Seinfeld and Pandis (1989) Atmospheric Chemistry and Physics, Wiley, 1998, p963
elemental real function laminar_conductance(T,p,sphum,ustar,diff) result(cond)
   real, intent(in) :: T       ! air temperature, deg K
   real, intent(in) :: p       ! pressure, Pa
   real, intent(in) :: sphum   ! specific humidity, kg/kg
   real, intent(in) :: ustar   ! friction velocity, m/s
   real, intent(in) :: diff    ! molecular diffusivity of the species (m2/s)

   real :: rho_air ! density of air, kg/m3
   ! ---- local vars
   real :: viscosity ! dynamic viscosity, kg/(m s)
   real :: Sc        ! Schmidt number, unitless

   viscosity = 1.458E-6 * T**1.5/(T+110.4)  ! Dynamic viscosity [kg/(m s)]
   rho_air = p/(rdgas*T*(1+d608*sphum))
   if (diff > 0) then
      Sc = viscosity/(rho_air*diff)
      cond = 0.2*ustar*Sc**(-2.0/3.0)
   else
      cond = 0.0
   endif
end function laminar_conductance


! ============================================================================
! updates concentration of tracers in the canopy air, taking into account dry
! deposition and exchange with the atmosphere
subroutine update_cana_tracers(tile, tr_flux, dfdtr, &
     precip_l, precip_s, pressure, ustar, con_g, con_v, stomatal_cond )
  type(land_tile_type), intent(inout) :: tile
  real, intent(in) :: tr_flux(:) ! fluxes of tracers
  real, intent(in) :: dfdtr  (:) ! derivatives of tracer fluxes w.r.t. concentrations
  real, intent(in) :: pressure   ! atmospheric pressure, N/m2
  real, intent(in) :: precip_l   ! liquid precipitation, kg/(m2 s)
  real, intent(in) :: precip_s   ! solid precipitation, kg/(m2 s)
  real, intent(in) :: ustar ! friction velocity, m/s
  ! real, intent(in) :: wind10 ! wind at 10 m above displacement height, m/s
  real, intent(in) :: con_g ! aerodynamic conductance between canopy air and ground for tracers
  real, intent(in) :: con_v ! aerodynamic conductance between canopy air and canopy for tracers
  real, intent(in) :: stomatal_cond ! integral stomatal conductance of canopy (that is, multiplied by LAI), for water vapor, m/s

  integer :: tr      ! tracer index
  integer :: species ! vegetation species index
  real    :: rho     ! density of canopy air
  real    :: dq      ! canopy air tracer tendency per time step
  real    :: con_v_lam ! quasi-laminar layer conductance for vegetation, m/s
  real    :: con_g_lam ! quasi-laminar layer conductance for ground, m/s
  real    :: con_st  ! total stomatal conductance, scaled by dry leaf area, m/s
  real    :: con_mx  ! "mesophyll conductance", m/s
  real    :: con_cu  ! total cuticular conductance, including dry and wet areas, m/s
  real    :: con_gr  ! "ground conductance", m/s
  real    :: cv0, cv1, cv2, cg0 ! intermediate values for total conductance calculations, m/s
  real    :: cv, cg  ! total conductances for vegetation and ground surface, m/s
  real    :: f_atm   ! flux of the tracer to the atmosphere, kg/(m2 s)
  real    :: ddep    ! dry deposition of the tracer, kg/(m2 s)
  real    :: LAI     ! leaf area index, m2/m2
  real    :: ustar_s ! friction velocity at the top of quasi-laminar layer, m/s
  real    :: emis(ntcana) ! tracer sources
  real    ::  ft, & ! fraction of canopy not covered by intercepted water/snow
              fw, & ! fraction of canopy covered by intercepted water
              fs    ! fraction of canopy covered by intercepted snow

  if (is_watch_point()) then
     write(*,*) 'update_cana_tracers input'
     __DEBUG1__(tr_flux)
     __DEBUG1__(dfdtr)
     __DEBUG1__(pressure)
     __DEBUG2__(precip_l, precip_s)
     __DEBUG1__(ustar)
     __DEBUG3__(con_g, con_v, stomatal_cond)
     write(*,*) 'end of update_cana_tracers input'
  endif

  ! update non-generic tracers, e.g.:
  ! call update_land_dust(tile, i, j, tr_flux, dfdtr, &
  !   precip_l, precip_s, pressure, ustar, wind10, con_g, con_v )
  ! wind10 is not passed to this subroutine yet

  ! update generic tracers
  ! calculate tracers sources
  emis(:) = 0.0
  ! TODO: add non-zero sources for generic tracers

  if (associated(tile%vegn)) then
     LAI = vegn_tile_LAI(tile%vegn)
  else
     LAI = 0.0
  endif
  ! ustar passed into this subroutine corresponds to the momentum dissipation
  ! (tau = rho*ustar**2) on the entire land surface. If vegetation is present,
  ! momentum dissipates on both leaves and ground; for the momentum conservation
  ! we need to reduce dissipation per unit surface proportionally to the total
  ! surface area. Factor of 2 takes into account that the leaves are 2-sided.
   ustar_s = ustar/sqrt(2*LAI+1.0)

  rho = pressure/(rdgas*tile%cana%T*(1+d608*tile%cana%tr(isphum))) ! air density
  if (associated(tile%vegn)) then
     ! our species
     species = tile%vegn%cohorts(1)%species
     ! adjust the stomatal conductance for the fraction covered by intercepted water/snow
     call get_vegn_wet_frac ( tile%vegn%cohorts(1), fw=fw, fs=fs )
     ft = 1 - fw - fs ! fraction of canopy not covered by water or snow
     con_st = ft*stomatal_cond
  else
     con_st = 0.0
  endif

  if (is_watch_point()) then
     __DEBUG4__(LAI,ustar,ustar_s,con_st)
  endif

  do tr = 1, ntcana
     if (.not.trdata(tr)%is_generic) cycle

     if (trdata(tr)%do_deposition) then
        if (associated(tile%vegn)) then
           ! conductance of quasi-laminar layers
           con_v_lam = LAI*laminar_conductance(tile%vegn%cohorts(1)%Tv, pressure, &
                           tile%cana%tr(isphum), ustar_s, diffusivity_h2o/trdata(tr)%diff_ratio)
           ! cuticular conductance: includes deposition od dry and wet portions of the
           ! leaves
           con_cu = LAI * ( ft*spdata(species)%tracer_cuticular_cond &
                              *(trdata(tr)%Henry_const*1e-5 + trdata(tr)%reactivity) &
                            + fw / trdata(tr)%r_lw + fs / trdata(tr)%r_ls )
        else
           con_v_lam = 0.0
           con_cu    = 0.0
        endif

        ! estimate conductances
        cv0 = 0.0; cv1 = 0.0; cv = 0.0
        ! aerodynamic (turbulent) + quasi-laminar conductance
        if (con_v+con_v_lam > 0) cv0 = con_v*con_v_lam/(con_v+con_v_lam)

        if (con_st > 0) then
           ! mesophyll conductance
           con_mx = trdata(tr)%Henry_const/3000.0 + 100.0*trdata(tr)%reactivity
           ! combined mesophyll + stomatal conductance
           cv1 = 1/(trdata(tr)%diff_ratio/con_st+1/con_mx)
        endif

        cv2 = cv1 + con_cu ! sum of stomatal/mesophyll and cuticular conductances

        ! total conductance between canopy air and canopy
        if (cv0+cv2>0) cv = cv0*cv2/(cv0+cv2)

        ! TODO: reactivity of the species on the ground must be taken into account
        !       this reactivity may be different for snow/water/soil
        con_g_lam = laminar_conductance(land_tile_grnd_T(tile), pressure, &
                               tile%cana%tr(isphum), ustar_s, diffusivity_h2o/trdata(tr)%diff_ratio)
        cg0 = con_g*con_g_lam/(con_g+con_g_lam)
        ! ground conductance
        con_gr = trdata(tr)%Henry_const*1e-5/trdata(tr)%r_gsS &
               + trdata(tr)%reactivity/trdata(tr)%r_gsO
        cg = cg0*con_gr/(cg0+con_gr)
     else
        cv = 0 ; cg = 0; con_v_lam = 0 ; con_g_lam = 0
     endif
     ! update tracer concentration -- this needs to be done even when do_deposition
     ! is FALSE, because source might be non-zero
     dq = (-tr_flux(tr)-rho*(cv+cg)*tile%cana%tr(tr)+emis(tr)) &
        / (canopy_air_mass_for_tracers/dt + dfdtr(tr) + rho*(cv+cg))
     if (is_watch_point()) then
        write(*,*) 'update_cana_tracers : ', trim(trdata(tr)%name)
        __DEBUG4__(tile%cana%tr(tr), tr_flux(tr), dfdtr(tr), emis(tr))
        __DEBUG3__(rho,cv,cg)
        __DEBUG2__(dq,tile%cana%tr(tr) + dq)
     endif
     tile%cana%tr(tr) = tile%cana%tr(tr) + dq
     ! ---- final values of the fluxes, for diagnostics
     ddep  = rho*(cv+cg)*tile%cana%tr(tr)
     f_atm = tr_flux(tr)+dfdtr(tr)*dq
     ! ---- diagnostic section
     call send_tile_data(trdata(tr)%id_con_v_lam,  con_v_lam,  tile%diag)
     call send_tile_data(trdata(tr)%id_con_g_lam,  con_g_lam,  tile%diag)
     call send_tile_data(trdata(tr)%id_con_v,      cv,         tile%diag)
     call send_tile_data(trdata(tr)%id_con_g,      cg,         tile%diag)
     call send_tile_data(trdata(tr)%id_emis,       emis(tr),   tile%diag)
     call send_tile_data(trdata(tr)%id_ddep,       ddep,       tile%diag)
     call send_tile_data(trdata(tr)%id_flux_atm,   f_atm,      tile%diag)
     call send_tile_data(trdata(tr)%id_dfdtr,      dfdtr(tr),  tile%diag)
     call send_tile_data(trdata(tr)%id_conc,       tile%cana%tr(tr), tile%diag)
  enddo

end subroutine update_cana_tracers


! ============================================================================
function flux_units(tracer_units) result (units)
   character(32) :: units
   character(*), intent(in) :: tracer_units

   select case (trim(lowercase(tracer_units)))
   case ('mmr')
     units = 'kg/m2/s'
   case ('kg/kg')
     units = 'kg/m2/s'
   case ('vmr')
     units = 'mole/m2/s'
   case ('mol/mol')
     units = 'mole/m2/s'
   case ('mole/mole')
     units = 'mole/m2/s'
   case default
     units = trim(tracer_units)//' kg/(m2 s)'
   end select
end function flux_units

end module land_tracer_driver_mod
