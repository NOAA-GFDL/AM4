 
module aerosol_types_mod

!--------------------------------------------------------------------

use time_manager_mod,  only: time_type
use interpolator_mod,  only: interpolate_type

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!------ public data structures ------
!--------------------------------------------------------------------

public   aerosol_type

type aerosol_type
     real,       dimension(:,:,:,:), pointer :: aerosol=>NULL()
     logical,    dimension(:,:),     pointer :: family_members=>NULL()
     character(len=64), dimension(:), pointer :: aerosol_names=>NULL()
end type aerosol_type

!--------------------------------------------------------------------

public aerosol_time_vary_type

!  Interp                interpolate_type variable containing the
!                        information about the aerosol species
!  Time                  time for which data is obtained from
!                        aerosol timeseries
!  being_overridden      is a given aerosol field to be overridden
!                        based on the model data_table?
!  output_override_info  should override info about each
!                        aerosol field be output (will be
!                        set to .false. after first time step)
!  override_counter      used to count calls to aerosol_endts
!                        so that output_override_info may
!                        set to .false. after physics_up
!  nfields               number of active aerosol species
!  nfamilies             number of active aerosol families

type aerosol_time_vary_type
     type(interpolate_type), dimension(:), pointer :: Interp=>NULL()
     type(time_type),        dimension(:), pointer :: Time=>NULL()
     logical, dimension(:), pointer :: being_overridden=>NULL()
     integer  :: nfields=0
     integer  :: nfamilies=0
     logical :: output_override_info = .true.
     integer :: override_counter = 0
     logical :: variable_is_initialized = .false.
end type aerosol_time_vary_type

!------------------------------------------------------------------

!------------------------------------------------------------------

public   aerosol_properties_type

!    aerextband
!    aerssalbband
!    aerasymmband
!    aerextbandlw
!    aerssalbbandlw
!    sulfate_index
!    optical_index

type aerosol_properties_type
     integer, dimension(:,:,:), pointer  :: ivol=>NULL()
     real, dimension(:,:), pointer  :: aerextband=>NULL(),   &
                                       aerssalbband=>NULL(), &
                                       aerasymmband=>NULL(), &
                                       aerextbandlw=>NULL(), &
                                       aerssalbbandlw=>NULL(), &
                                       aerextbandlw_cn=>NULL(), &
                                       aerssalbbandlw_cn=>NULL()
     real, dimension(:,:,:,:), pointer :: sw_ext=>NULL(), &
                                          sw_ssa=>NULL(), &
                                          sw_asy=>NULL(), &
                                          lw_ext=>NULL(), &
                                          lw_ssa=>NULL(), &
                                          lw_asy=>NULL()
!yim
     integer, dimension(:,:), pointer :: sulfate_index=>NULL()
     integer, dimension(:), pointer :: optical_index=>NULL()
     integer, dimension(:), pointer :: omphilic_index=>NULL()
     integer, dimension(:), pointer :: bcphilic_index=>NULL()
     integer, dimension(:), pointer :: seasalt1_index=>NULL()
     integer, dimension(:), pointer :: seasalt2_index=>NULL()
     integer, dimension(:), pointer :: seasalt3_index=>NULL()
     integer, dimension(:), pointer :: seasalt4_index=>NULL()
     integer, dimension(:), pointer :: seasalt5_index=>NULL()
     integer, dimension(:), pointer :: seasalt_aitken_index=>NULL()
     integer, dimension(:), pointer :: seasalt_fine_index=>NULL()
     integer, dimension(:), pointer :: seasalt_coarse_index=>NULL()
     integer                        :: sulfate_flag
     integer                        :: omphilic_flag
     integer                        :: bcphilic_flag
     integer                        :: seasalt1_flag
     integer                        :: seasalt2_flag
     integer                        :: seasalt3_flag
     integer                        :: seasalt4_flag
     integer                        :: seasalt5_flag
     integer                        :: seasalta_flag
     integer                        :: seasaltf_flag
     integer                        :: seasaltc_flag
!yim
     integer                        :: bc_flag
end type aerosol_properties_type



!------------------------------------------------------------------

public   aerosol_diagnostics_type

!    extopdep
!    absopdep
!    extopdep_vlcno
!    absopdep_vlcno
!    lw_extopdep_vlcno
!    lw_absopdep_vlcno
!    sw_heating_vlcno

type aerosol_diagnostics_type
     real, dimension(:,:,:,:),   pointer :: sw_heating_vlcno=>NULL()
     real, dimension(:,:,:,:,:), pointer :: extopdep=>NULL(), &
                                            absopdep=>NULL()
     real, dimension(:,:,:,:,:), pointer :: asymdep=>NULL()

     real, dimension(:,:,:,:),   pointer :: extopdep_vlcno=>NULL(), &
                                            absopdep_vlcno=>NULL(), &
                                            lw_extopdep_vlcno=>NULL(), &
                                            lw_absopdep_vlcno=>NULL()
end type aerosol_diagnostics_type

!------------------------------------------------------------------


end module aerosol_types_mod

