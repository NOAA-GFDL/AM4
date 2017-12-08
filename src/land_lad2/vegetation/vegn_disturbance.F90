! ============================================================================
! vegetation disturbances
! ============================================================================
module vegn_disturbance_mod

use land_constants_mod, only : seconds_per_year
use land_data_mod,   only : log_version
use vegn_data_mod,   only : spdata, fsc_wood, fsc_liv, fsc_froot, agf_bs, LEAF_OFF
use vegn_tile_mod,   only : vegn_tile_type
use soil_tile_mod,   only : soil_tile_type
use vegn_cohort_mod, only : vegn_cohort_type, height_from_biomass, lai_from_biomass, &
     update_biomass_pools
use soil_mod,        only : add_root_litter
use soil_carbon_mod, only : add_litter, soil_carbon_option, &
     SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE
use fms_mod, only : error_mesg, FATAL

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_disturbance_init
public :: vegn_nat_mortality
public :: vegn_disturbance
public :: update_fuel
! =====end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'vegn_disturbance_mod'
#include "../shared/version_variable.inc"

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! =============================================================================
subroutine vegn_disturbance_init()
  call log_version(version, module_name, &
  __FILE__)
end subroutine vegn_disturbance_init

subroutine vegn_disturbance(vegn, soil, dt)
  type(vegn_tile_type), intent(inout) :: vegn ! vegetation data
  type(soil_tile_type), intent(inout) :: soil ! soil data
  real, intent(in) :: dt ! time since last disturbance calculations, s


  real, parameter :: BMIN = 1e-10; ! should be the same as in growth function
  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  real :: precip;
  real :: delta;
  real :: fraction_lost;
  real :: drought_month;
  real :: deltat
  integer :: i
  integer :: sp ! shorhand for cohort species
  real :: new_fast_C_ag, new_fast_C_bg, new_slow_C_ag, new_slow_C_bg

  deltat = dt/seconds_per_year ! convert time interval to years

  !  Disturbance Rates
  precip=vegn%p_ann*86400*365;
  drought_month = vegn%lambda;
  vegn%disturbance_rate(0) = 0.0;
  vegn%disturbance_rate(1) = 0.0;

  call calculate_patch_disturbance_rates(vegn)

  ! Fire disturbance implicitly, i.e.  not patch creating
  vegn%area_disturbed_by_fire = (1.0-exp(-vegn%disturbance_rate(1)*deltat));

  do i = 1,vegn%n_cohorts
     cc => vegn%cohorts(i)
     sp = cc%species

     fraction_lost = 1.0-exp(-vegn%disturbance_rate(1)*deltat);

     ! "dead" biomass : wood + sapwood
     delta = (cc%bwood+cc%bsw)*fraction_lost;
     select case (soil_carbon_option)
     case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
        soil%slow_soil_C(1) = soil%slow_soil_C(1) + (1.0-spdata(sp)%smoke_fraction)*delta*(1-fsc_wood);
        soil%fast_soil_C(1) = soil%fast_soil_C(1) + (1.0-spdata(sp)%smoke_fraction)*delta*   fsc_wood;
     case (SOILC_CORPSE)
        call add_litter(soil%coarseWoodLitter,(/(1.0-spdata(sp)%smoke_fraction)*delta*fsc_wood*agf_bs,&
                       (1.0-spdata(sp)%smoke_fraction)*delta*(1.0-fsc_wood)*agf_bs,0.0/))
        soil%coarsewoodlitter_fsc_in=soil%coarsewoodlitter_fsc_in+(1.0-spdata(sp)%smoke_fraction)*delta*fsc_wood*agf_bs
        soil%coarsewoodlitter_ssc_in=soil%coarsewoodlitter_ssc_in+(1.0-spdata(sp)%smoke_fraction)*delta*(1.0-fsc_wood)*agf_bs

        !fsc_in and ssc_in are updated in add_root_litter
        call add_root_litter(soil,vegn,(/(1.0-spdata(sp)%smoke_fraction)*delta*fsc_wood*(1.0-agf_bs),&
                       (1.0-spdata(sp)%smoke_fraction)*delta*(1.0-fsc_wood)*(1.0-agf_bs),0.0/))
     case default
        call error_mesg('vegn_disturbance','soil_carbon_option is invalid. This should never happen. Contact developer.', FATAL)
     end select

     cc%bwood = cc%bwood * (1-fraction_lost);
     cc%bsw   = cc%bsw   * (1-fraction_lost);

     vegn%csmoke_pool = vegn%csmoke_pool + spdata(sp)%smoke_fraction*delta;

     ! for budget tracking - temporarily not keeping wood and the rest separately,ens
     !      soil%ssc_in(1)+=delta*(1.0-spdata(sp)%smoke_fraction)*(1-fsc_wood); */
     !      soil%fsc_in(1)+=delta*(1.0-spdata(sp)%smoke_fraction)*fsc_wood; */
     if (soil_carbon_option==SOILC_CENTURY.or.soil_carbon_option==SOILC_CENTURY_BY_LAYER) then
        soil%ssc_in(1) = soil%ssc_in(1)+(cc%bwood+cc%bsw)*fraction_lost *(1.0-spdata(sp)%smoke_fraction)
        !     soil%fsc_in(1)+=cc%bsw*fraction_lost *(1.0-spdata(sp)%smoke_fraction);
     endif
     vegn%veg_out = vegn%veg_out+delta;

     !"alive" biomass: leaves, roots, and virtual pool
     delta = (cc%bl+cc%blv+cc%br)*fraction_lost;
     select case (soil_carbon_option)
     case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
        soil%fast_soil_C(1) = soil%fast_soil_C(1) + (1.0-spdata(sp)%smoke_fraction)*delta*    fsc_liv ;
        soil%slow_soil_C(1) = soil%slow_soil_C(1) + (1.0-spdata(sp)%smoke_fraction)*delta*(1- fsc_liv);
        soil%fsc_in(1) = soil%fsc_in(1)+delta*(1.0-spdata(sp)%smoke_fraction);
     case (SOILC_CORPSE)
        new_fast_C_ag=(cc%bl+cc%blv)*fraction_lost*fsc_liv*(1.0-spdata(sp)%smoke_fraction)
        new_slow_C_ag=(cc%bl+cc%blv)*fraction_lost*(1.0-fsc_liv)*(1.0-spdata(sp)%smoke_fraction)
        new_fast_C_bg=cc%br*fraction_lost*fsc_froot*(1.0-spdata(sp)%smoke_fraction)
        new_slow_C_bg=cc%br*fraction_lost*(1.0-fsc_froot)*(1.0-spdata(sp)%smoke_fraction)
        call add_litter(soil%leafLitter,(/new_fast_C_ag,new_slow_C_ag,0.0/))
        soil%leaflitter_fsc_in=soil%leaflitter_fsc_in+new_fast_C_ag
        soil%leaflitter_ssc_in=soil%leaflitter_ssc_in+new_slow_C_ag
        !fsc_in and ssc_in updated in add_root_litter
        call add_root_litter(soil,vegn,(/new_fast_C_bg,new_slow_C_bg,0.0/))
     case default
        call error_mesg('vegn_disturbance','soil_carbon_option is invalid. This should never happen. Contact developer.', FATAL)
     end select

     cc%bl  = cc%bl  * (1-fraction_lost);
     cc%blv = cc%blv * (1-fraction_lost);
     cc%br  = cc%br  * (1-fraction_lost);

     vegn%csmoke_pool = vegn%csmoke_pool + spdata(sp)%smoke_fraction*delta;

     vegn%veg_out = vegn%veg_out+delta;

     !"living" biomass:leaves, roots and sapwood
     delta = cc%bliving*fraction_lost;
     cc%bliving = cc%bliving - delta;

     if(cc%bliving < BMIN) then
        ! remove vegetaion competely
        select case (soil_carbon_option)
        case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
           soil%fast_soil_C(1) = soil%fast_soil_C(1) + fsc_liv*cc%bliving+ fsc_wood*cc%bwood;
           soil%slow_soil_C(1) = soil%slow_soil_C(1) + (1.- fsc_liv)*cc%bliving+ (1-fsc_wood)*cc%bwood;

           soil%fsc_in(1) = soil%fsc_in(1) + cc%bwood+cc%bliving;
        case (SOILC_CORPSE)
           ! remove vegetation competely
           new_fast_C_ag=fsc_liv*(cc%bl+cc%blv+cc%bsw*agf_bs) + fsc_wood*(cc%bwood)*agf_bs
           new_slow_C_ag=(1-fsc_liv)*(cc%bl+cc%blv+cc%bsw*agf_bs) + (1-fsc_wood)*(cc%bwood)*agf_bs
           new_fast_C_bg=fsc_froot*cc%br+fsc_liv*cc%bsw*(1-agf_bs) + fsc_wood*cc%bwood*(1-agf_bs)
           new_slow_C_bg=(1-fsc_froot)*cc%br + (1-fsc_liv)*cc%bsw*(1-agf_bs) + (1-fsc_wood)*cc%bwood*(1-agf_bs)

           call add_litter(soil%leafLitter,(/fsc_liv*(cc%bl),(1-fsc_liv)*(cc%bl),0.0/))
           call add_litter(soil%coarseWoodLitter,(/fsc_liv*(cc%blv+cc%bsw*agf_bs) + fsc_wood*(cc%bwood)*agf_bs,&
                                  (1-fsc_liv)*(cc%blv+cc%bsw*agf_bs) + (1-fsc_wood)*(cc%bwood)*agf_bs,0.0/))
           soil%leaflitter_fsc_in=soil%leaflitter_fsc_in+fsc_liv*(cc%bl)
           soil%leaflitter_ssc_in=soil%leaflitter_ssc_in+(1-fsc_liv)*(cc%bl)
           soil%coarsewoodlitter_fsc_in=soil%coarsewoodlitter_fsc_in+&
                  fsc_liv*(cc%blv+cc%bsw*agf_bs) + fsc_wood*(cc%bwood)*agf_bs
           soil%coarsewoodlitter_ssc_in=soil%coarsewoodlitter_ssc_in+&
                  (1-fsc_liv)*(cc%blv+cc%bsw*agf_bs) + (1-fsc_wood)*(cc%bwood)*agf_bs
           call add_root_litter(soil, vegn, (/new_fast_C_bg,new_slow_C_bg,0.0/))
        case default
           call error_mesg('vegn_disturbance','soil_carbon_option is invalid. This should never happen. Contact developer.', FATAL)
        end select
        vegn%veg_out = vegn%veg_out + cc%bwood+cc%bliving;

        cc%bliving = 0.;
        cc%bwood   = 0.;
     endif
     call update_biomass_pools(cc)
  enddo

  vegn%csmoke_rate = vegn%csmoke_pool; ! kg C/(m2 yr)
end subroutine vegn_disturbance

! ============================================================================
subroutine calculate_patch_disturbance_rates(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  real :: fire_prob;
  real :: fuel;

  fuel = vegn%fuel

#if SIMPLE_FIRE
  ! CALCULATE FIRE DISTURBANCE RATES
  vegn%disturbance_rate(1)=fire(vegn);
#else

  ! lambda is the number of drought months;
  fire_prob = vegn%lambda/(1.+vegn%lambda);
  ! compute average fuel during fire months
  if (vegn%lambda > 0.00001 ) fuel = fuel/vegn%lambda;
  vegn%disturbance_rate(1) = fuel * fire_prob;

  ! put a threshold for very dry years for warm places
  if (vegn%t_ann > 273.16 .and. vegn%lambda > 3.)  vegn%disturbance_rate(1)=0.33;
#endif

  if(vegn%disturbance_rate(1) > 0.33) vegn%disturbance_rate(1)=0.33;

  ! this is only true for the one cohort per patch case
  vegn%disturbance_rate(0) = spdata(vegn%cohorts(1)%species)%treefall_disturbance_rate;
  vegn%total_disturbance_rate = vegn%disturbance_rate(1)+vegn%disturbance_rate(0);

  vegn%fuel = fuel;
end subroutine calculate_patch_disturbance_rates


! ============================================================================
function fire(vegn) result(fireterm)
  real :: fireterm; ! return value
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  real :: precip_av ! average precipitation, mm/year

  fireterm = 0
  precip_av = vegn%p_ann * seconds_per_year;
!!$  vegn%ignition_rate = 0.00;
  vegn%fuel = vegn%total_biomass;

  if(vegn%fuel>0.0) then
     if(precip_av < 400.+40.*(vegn%t_ann-273.16)) then
        fireterm = vegn%fuel*(400. + 40.*(vegn%t_ann-273.16) - precip_av);
     endif
  endif
end function


! ============================================================================
subroutine update_fuel(vegn, wilt)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: wilt ! ratio of wilting to saturated water content

  ! ---- local constants
  !  these three used to be in data
  real, parameter :: fire_height_threashold = 100;
  real, parameter :: fp1 = 1.; ! disturbance rate per kgC/m2 of fuel
  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  real :: theta_crit; ! critical ratio of average soil water to sat. water
  real :: ignition_rate;
  real ::  babove;
  integer :: i

  do i = 1,vegn%n_cohorts
     cc => vegn%cohorts(i)
     ! calculate theta_crit: actually either fact_crit_fire or cnst_crit_fire
     ! is zero, enforced by logic in the vegn_data.F90
     theta_crit = spdata(cc%species)%cnst_crit_fire &
           + wilt*spdata(cc%species)%fact_crit_fire
     theta_crit = max(0.0,min(1.0, theta_crit))
     if((cc%height < fire_height_threashold) &
          .and.(vegn%theta_av_fire < theta_crit)  &
          .and.(vegn%tsoil_av > 278.16)) then
        babove = cc%bl + agf_bs * (cc%bsw + cc%bwood + cc%blv);
        ! this is fuel available durng the drought months only
        vegn%fuel = vegn%fuel + spdata(cc%species)%fuel_intensity*babove;
     endif
  enddo

  ! note that the ignition rate calculation based on the value of theta_crit for
  ! the last cohort -- currently it doesn't matter since we have just one cohort,
  ! but something needs to be done about that in the future
  ignition_rate = 0.;
  if ( (vegn%theta_av_fire < theta_crit) &
       .and. (vegn%tsoil_av>278.16)) ignition_rate = 1.;
  vegn%lambda = vegn%lambda + ignition_rate;

end subroutine update_fuel




! ============================================================================
subroutine vegn_nat_mortality(vegn, soil, deltat)
  type(vegn_tile_type), intent(inout) :: vegn  ! vegetation data
  type(soil_tile_type), intent(inout) :: soil  ! soil data
  real, intent(in) :: deltat ! time since last mortality calculations, s

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  real :: delta;
  real :: fraction_lost;
  real :: bdead, balive; ! combined biomass pools
  integer :: i

  vegn%disturbance_rate(0)        = 0.0;
  vegn%area_disturbed_by_treefall = 0.0;

  do i = 1,vegn%n_cohorts
     cc => vegn%cohorts(i)
     ! Treat treefall disturbance implicitly, i.e. not creating a new tile.
     ! note that this disturbance rate calculation only works for the one cohort per
     ! tile case -- in case of multiple cohort disturbance rate perhaps needs to be
     ! accumulated (or averaged? or something else?) over the cohorts.
     vegn%disturbance_rate(0) = spdata(cc%species)%treefall_disturbance_rate;
     vegn%area_disturbed_by_treefall = &
          1.0-exp(-vegn%disturbance_rate(0)*deltat/seconds_per_year);

     ! calculate combined biomass pools
     balive = cc%bl + cc%blv + cc%br;
     bdead  = cc%bsw + cc%bwood;
     ! ens need a daily PATCH_FREQ here, for now it is set to 48
     fraction_lost = 1.0-exp(-vegn%disturbance_rate(0)*deltat/seconds_per_year);

     ! "dead" biomass : wood + sapwood
     delta = bdead*fraction_lost;

     select case (soil_carbon_option)
     case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
        soil%slow_soil_C(1) = soil%slow_soil_C(1) + (1-fsc_wood)*delta;
        soil%fast_soil_C(1) = soil%fast_soil_C(1) +    fsc_wood *delta;
     case (SOILC_CORPSE)
        !Add above ground fraction to top soil layer, and the rest to the soil profile
        call add_litter(soil%coarseWoodLitter,(/fsc_wood *delta*agf_bs,(1-fsc_wood)*delta*agf_bs,0.0/))
        soil%coarseWoodlitter_fsc_in=soil%coarseWoodlitter_fsc_in+fsc_wood *delta*agf_bs
        soil%coarseWoodlitter_ssc_in=soil%coarseWoodlitter_ssc_in+(1-fsc_wood)*delta*agf_bs
        !fsc_in and ssc_in updated in add_root_litter
        call add_root_litter(soil,vegn,  (/fsc_wood *delta*(1.0-agf_bs),(1-fsc_wood)*delta*(1.0-agf_bs),0.0/))
     case default
        call error_mesg('vegn_nat_mortality','soil_carbon_option is invalid. This should never happen. Contact developer.', FATAL)
     end select

     cc%bwood = cc%bwood * (1-fraction_lost);
     cc%bsw   = cc%bsw   * (1-fraction_lost);

     ! for budget tracking -temporarily
     ! It doesn't look correct to me: ssc_in should probably include factor
     ! (1-fsc_wood) and the whole calculations should be moved up in front
     ! of bwood and bsw modification
     ! soil%fsc_in(1)+= cc%bsw*fraction_lost;
     if (soil_carbon_option==SOILC_CENTURY.or.soil_carbon_option==SOILC_CENTURY_BY_LAYER) then
         soil%ssc_in(1)  = soil%ssc_in(1)  + (cc%bwood+cc%bsw)*fraction_lost;
     endif
     vegn%veg_out = vegn%veg_out + delta;

     ! kill the live biomass if mortality is set to affect it
     if (spdata(cc%species)%mortality_kills_balive) then
        delta = (cc%bl + cc%blv + cc%br)*fraction_lost

        select case (soil_carbon_option)
        case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
           soil%slow_soil_C(1) = soil%slow_soil_C(1) + (1-fsc_liv)*delta;
           soil%fast_soil_C(1) = soil%fast_soil_C(1) +    fsc_liv *delta;
        case (SOILC_CORPSE)
           call add_litter(soil%leafLitter,(/(cc%bl+cc%blv)*fraction_lost*fsc_liv,(cc%bl+cc%blv)*fraction_lost*(1.0-fsc_liv),0.0/))
           call add_root_litter(soil,vegn,(/cc%br*fraction_lost*fsc_froot,cc%br*fraction_lost*(1.0-fsc_froot),0.0/))
        case default
           call error_mesg('vegn_nat_mortality','soil_carbon_option is invalid. This should never happen. Contact developer.', FATAL)
        end select

        cc%br  = cc%br  * (1-fraction_lost);
        cc%bl  = cc%bl  * (1-fraction_lost);
        cc%blv = cc%blv * (1-fraction_lost);

        ! for budget tracking
        soil%ssc_in(1)  = soil%ssc_in(1)  + (cc%bwood+cc%bsw)*fraction_lost;
        vegn%veg_out = vegn%veg_out + delta;
     endif

     cc%bliving = cc%bsw + cc%bl + cc%br + cc%blv;
     call update_biomass_pools(cc);
  enddo

end subroutine vegn_nat_mortality


end module vegn_disturbance_mod
