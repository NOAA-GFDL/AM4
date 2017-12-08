module tiling_input_types_mod

type :: lake_predefined_type
 integer :: nc_grpid,nlake,nband
 real,pointer,dimension(:) :: frac
 real,pointer,dimension(:) :: w_sat
 real,pointer,dimension(:) :: awc_lm2
 real,pointer,dimension(:) :: k_sat_ref
 real,pointer,dimension(:) :: psi_sat_ref
 real,pointer,dimension(:) :: chb
 real,pointer,dimension(:) :: alpha
 real,pointer,dimension(:) :: heat_capacity_ref
 real,pointer,dimension(:) :: thermal_cond_ref
 real,pointer,dimension(:,:) :: refl_dry_dir
 real,pointer,dimension(:,:) :: refl_dry_dif
 real,pointer,dimension(:,:) :: refl_sat_dir
 real,pointer,dimension(:,:) :: refl_sat_dif
 real,pointer,dimension(:) :: emis_dry
 real,pointer,dimension(:) :: emis_sat
 real,pointer,dimension(:) :: z0_momentum
 real,pointer,dimension(:) :: z0_momentum_ice
 real,pointer,dimension(:) :: depth_sill
 real,pointer,dimension(:) :: width_sill
 real,pointer,dimension(:) :: whole_area
 real,pointer,dimension(:) :: connected_to_next
 real,pointer,dimension(:) :: backwater
 real,pointer,dimension(:) :: backwater_1
 real,pointer,dimension(:) :: rsa_exp         ! riparian source-area exponent
end type lake_predefined_type

type :: glacier_predefined_type

 integer :: nc_grpid,nglacier,nband
 real,pointer,dimension(:) :: frac
 real,pointer,dimension(:) :: w_sat
 real,pointer,dimension(:) :: awc_lm2
 real,pointer,dimension(:) :: k_sat_ref
 real,pointer,dimension(:) :: psi_sat_ref
 real,pointer,dimension(:) :: chb
 real,pointer,dimension(:) :: alpha
 real,pointer,dimension(:) :: heat_capacity_ref
 real,pointer,dimension(:) :: thermal_cond_ref
 real,pointer,dimension(:,:) :: refl_max_dir
 real,pointer,dimension(:,:) :: refl_max_dif
 real,pointer,dimension(:,:) :: refl_min_dir
 real,pointer,dimension(:,:) :: refl_min_dif
 real,pointer,dimension(:) :: emis_dry
 real,pointer,dimension(:) :: emis_sat
 real,pointer,dimension(:) :: z0_momentum
 real,pointer,dimension(:) :: tfreeze

end type glacier_predefined_type

type :: soil_predefined_type

 !miscellanous
 integer :: nsoil,nc_grpid,nband
 real,pointer,dimension(:) :: frac
 !soil
 real,pointer,dimension(:) :: dat_w_sat
 real,pointer,dimension(:) :: dat_awc_lm2
 real,pointer,dimension(:) :: dat_k_sat_ref
 real,pointer,dimension(:) :: dat_psi_sat_ref
 real,pointer,dimension(:) :: dat_chb
 real,pointer,dimension(:) :: dat_heat_capacity_dry
 real,pointer,dimension(:) :: dat_thermal_cond_dry
 real,pointer,dimension(:) :: dat_thermal_cond_sat
 real,pointer,dimension(:) :: dat_thermal_cond_exp
 real,pointer,dimension(:) :: dat_thermal_cond_scale
 real,pointer,dimension(:) :: dat_thermal_cond_weight
 real,pointer,dimension(:,:) :: dat_refl_dry_dir 
 real,pointer,dimension(:,:) :: dat_refl_dry_dif
 real,pointer,dimension(:,:) :: dat_refl_sat_dir 
 real,pointer,dimension(:,:) :: dat_refl_sat_dif 
 real,pointer,dimension(:) :: dat_emis_dry
 real,pointer,dimension(:) :: dat_emis_sat
 real,pointer,dimension(:) :: dat_z0_momentum
 real,pointer,dimension(:) :: dat_tf_depr
 real,pointer,dimension(:) :: rsa_exp_global
 real,pointer,dimension(:) :: gw_res_time
 real,pointer,dimension(:) :: gw_hillslope_length
 real,pointer,dimension(:) :: gw_scale_length
 real,pointer,dimension(:) :: gw_hillslope_zeta_bar
 real,pointer,dimension(:) :: gw_hillslope_relief
 real,pointer,dimension(:) :: gw_scale_relief
 real,pointer,dimension(:) :: gw_soil_e_depth
 real,pointer,dimension(:) :: gw_scale_soil_depth
 real,pointer,dimension(:) :: gw_hillslope_a
 real,pointer,dimension(:) :: gw_hillslope_n
 real,pointer,dimension(:) :: gw_perm
 real,pointer,dimension(:) :: gw_scale_perm
 !hillslope tiling
 real,pointer,dimension(:) :: microtopo
 real,pointer,dimension(:) :: tile_hlsp_length
 real,pointer,dimension(:) :: tile_hlsp_slope
 real,pointer,dimension(:) :: tile_hlsp_elev
 real,pointer,dimension(:) :: tile_hlsp_hpos
 real,pointer,dimension(:) :: tile_hlsp_width
 integer,pointer,dimension(:) :: hidx_k
 integer,pointer,dimension(:) :: hidx_j
 !vegetation
 integer,pointer,dimension(:) :: vegn

end type soil_predefined_type

type :: metadata_predefined_type

 integer :: ntile,nband
 integer,pointer,dimension(:) :: tid
 integer,pointer,dimension(:) :: tile
 integer,pointer,dimension(:) :: ttype
 real,pointer,dimension(:) :: frac

end type metadata_predefined_type

! public type -- used to temporarily store all the input data for a cell
type :: tile_parameters_type

 !metadata
 type(metadata_predefined_type),pointer :: metadata
 !soil
 type(soil_predefined_type),pointer :: soil
 !lake
 type(lake_predefined_type),pointer :: lake
 !glacier
 type(glacier_predefined_type),pointer :: glacier

end type tile_parameters_type

end module
