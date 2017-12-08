module tropchem_types_mod

!f1p
!store nml
!facilitate diagnostics

  public

  real :: small_value

  type tropchem_diag

!id for some diagnostics
     integer :: nb_diag
     integer :: ind_pso4_h2o2
     integer :: ind_pso4_o3
     integer :: ind_pH     
     integer :: ind_cH     
     integer :: ind_lwc    
     integer :: ind_enh4,ind_ehcoo,ind_ech3coo,ind_ehco3,ind_eco3,ind_eoh,ind_eno3,ind_eso4,ind_ehso3,ind_eso3,ind_ealk
     integer :: ind_phno3_d(5), ind_phno3_g_d, ind_ghno3_d
     integer :: ind_pso4_d(5), ind_pso4_g_d, ind_surfc_nh4no3, ind_surfc_isorropia

  end type tropchem_diag

  type tropchem_opt

!to store nml. Add options here
     logical               :: retain_cm3_bugs
     logical               :: do_fastjx_photo 
     integer               :: aerosol_thermo
     real                  :: gN2O5
     real                  :: gNO3
     real                  :: gNO2
     real                  :: gSO2
     real                  :: gH2SO4_dust
     real                  :: gSO2_dust
     real                  :: gNH3
     real                  :: gHNO3_dust
     real                  :: gNO3_dust
     real                  :: gN2O5_dust
     integer               :: gHNO3_dust_dynamic 
     real                  :: gSO2_dynamic
     real                  :: gHO2
     real                  :: min_lwc_for_cloud_chem
     logical               :: check_convergence
     logical               :: use_lsc_in_fastjx
     integer               :: cloud_chem, cloud_chem_ph_solver
     logical               :: het_chem_fine_aerosol_only
     real                  :: cloud_H
     logical               :: cloud_ho2_h2o2
!     logical               :: do_h2so4_nucleation 
!     real                  :: frac_dust_incloud
     real                  :: frac_aerosol_incloud
     real                  :: max_rh_aerosol
     logical               :: limit_no3
     integer               :: het_chem
     character(len=128)    :: sim_data_flsp
     logical               :: time_varying_solarflux
  end type tropchem_opt

  CONTAINS

  subroutine tropchem_types_init(trop_diag,ismall_value)
!initialize all indices for diagnostic to 0
!setup small value

    real, intent(in)                   :: ismall_value
    type(tropchem_diag), intent(inout) :: trop_diag

    trop_diag%nb_diag        = 0
    trop_diag%ind_pso4_h2o2  = 0
    trop_diag%ind_pso4_o3    = 0
    trop_diag%ind_cH         = 0
    trop_diag%ind_lwc        = 0
    trop_diag%ind_phno3_d    = 0
    trop_diag%ind_phno3_g_d  = 0
    trop_diag%ind_ghno3_d    = 0

    small_value = ismall_value
    
  end subroutine tropchem_types_init

end module tropchem_types_mod
