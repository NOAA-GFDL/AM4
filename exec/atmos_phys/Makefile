# Makefile created by mkmf 19.3.0

CPPDEFS = -DINTERNAL_FILE_NML -g -DCLUBB

OTHERFLAGS = -I$(BUILDROOT)fms -I$(SRCROOT)shared/include

TMPFILES = *.f90 *.mod

include $(MK_TEMPLATE)

.DEFAULT:
	-echo $@ does not exist.
all: libatmos_phys.a
AM3.mat.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.mat.F90 AM3.mods.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.mat.F90
AM3.mods.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.mods.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.mods.F90
AM3.subs.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.subs.F90 AM3.mods.o mo_jpl.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.subs.F90
AM3_fastjx.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3_fastjx.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3_fastjx.F90
AM3_fphoto.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3_fphoto.F90 AM3_fastjx.o mo_chem_utls.o moz.mods.o AM3.mods.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3_fphoto.F90
CLUBB_3D_var.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_3D_var.F90 parameters_model.o constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_3D_var.F90
CLUBB_driver_SCM.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_driver_SCM.F90 aerosol_types.o aer_ccn_act.o aer_ccn_act_k.o ice_nucl.o alt_cloud.o clubb_core.o clubb_precision.o constants_clubb.o error_code.o grid_class.o parameter_indices.o parameters_model.o parameters_tunable.o stats_subs.o stats_variables.o T_in_K_module.o variables_prognostic_module.o polysvp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_driver_SCM.F90
MISR_simulator.o: $(SRCROOT)atmos_phys/atmos_param/cosp/MISR_simulator/MISR_simulator.F $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/MISR_simulator/MISR_simulator.F
Skw_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/Skw_module.F90 constants_clubb.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/Skw_module.F90
T_in_K_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/T_in_K_module.F90 constants_clubb.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/T_in_K_module.F90
advance_helper_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_helper_module.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_helper_module.F90
advance_sclrm_Nd_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/advance_sclrm_Nd_module.F90 clubb_precision.o grid_class.o parameters_model.o lapack_wrap.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/advance_sclrm_Nd_module.F90
advance_windm_edsclrm_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_windm_edsclrm_module.F90 grid_class.o parameters_model.o parameters_tunable.o model_flags.o clubb_precision.o stats_type.o stats_variables.o clip_explicit.o error_code.o constants_clubb.o sponge_layer_damping.o lapack_wrap.o diffusion.o mean_adv.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_windm_edsclrm_module.F90
advance_wp2_wp3_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_wp2_wp3_module.F90 grid_class.o parameters_tunable.o stats_type.o stats_variables.o constants_clubb.o model_flags.o clubb_precision.o error_code.o lapack_wrap.o fill_holes.o clip_explicit.o csr_matrix_class_3array.o gmres_wrap.o gmres_cache.o diffusion.o mean_adv.o hyper_diffusion_4th_ord.o advance_helper_module.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_wp2_wp3_module.F90
advance_xm_wpxp_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_xm_wpxp_module.F90 parameters_tunable.o constants_clubb.o parameters_model.o grid_class.o model_flags.o mono_flux_limiter.o pdf_parameter_module.o clubb_precision.o error_code.o stats_type.o stats_variables.o sponge_layer_damping.o diffusion.o mean_adv.o clip_semi_implicit.o advance_helper_module.o lapack_wrap.o pos_definite_module.o clip_explicit.o fill_holes.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_xm_wpxp_module.F90
advance_xp2_xpyp_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_xp2_xpyp_module.F90 constants_clubb.o model_flags.o parameters_tunable.o parameters_model.o grid_class.o clubb_precision.o clip_explicit.o stats_type.o error_code.o stats_variables.o array_index.o diffusion.o mean_adv.o advance_helper_module.o lapack_wrap.o interpolation.o fill_holes.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_xp2_xpyp_module.F90
aer_ccn_act.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90 aer_ccn_act_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90
aer_ccn_act_k.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90
aer_in_act.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90
aerosol.o: $(SRCROOT)atmos_phys/atmos_shared/aerosol/aerosol.F90 aerosol_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/aerosol/aerosol.F90
aerosol_cloud.o: $(SRCROOT)atmos_phys/atmos_param/aerosol_cloud/aerosol_cloud.F90 aerosol_types.o polysvp.o aerosol_params.o aer_ccn_act.o ice_nucl.o lscloud_types.o lscloud_debug.o moist_processes_utils.o physics_radiation_exch.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/aerosol_cloud/aerosol_cloud.F90
aerosol_params.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90
aerosol_thermodynamics.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.F90 tropchem_types_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.F90
aerosol_types.o: $(SRCROOT)atmos_phys/atmos_shared/aerosol/aerosol_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/aerosol/aerosol_types.F90
aerosolrad_driver.o: $(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_driver.F90 aerosol_types.o aerosol.o aerosolrad_types.o aerosolrad_package.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_driver.F90
aerosolrad_package.o: $(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_package.F90 aerosol_types.o aerosolrad_types.o sealw99.o esfsw_driver.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_package.F90
aerosolrad_types.o: $(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_types.F90
alt_cloud.o: $(SRCROOT)atmos_phys/atmos_param/clubb/alt_cloud.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/alt_cloud.F90
anl_erf.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/anl_erf.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/anl_erf.F90
array_index.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/array_index.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/array_index.F90
array_lib.o: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/array_lib.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H mrgrnk.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/array_lib.F90
atmos_age_tracer.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_age_tracer.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_age_tracer.F90
atmos_carbon_aerosol.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_carbon_aerosol.F90 atmos_cmip_diag.o atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_carbon_aerosol.F90
atmos_ch3i.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_ch3i.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_ch3i.F90
atmos_ch4.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_ch4.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_ch4.F90
atmos_cmip_diag.o: $(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.F90 $(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_interp.inc
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag	$(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.F90
atmos_co2.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_co2.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_co2.F90
atmos_convection_tracer.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_convection_tracer.F90 atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_convection_tracer.F90
atmos_dust.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_dust.F90 atmos_cmip_diag.o atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_dust.F90
atmos_global_diag.o: $(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag/atmos_global_diag.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag/atmos_global_diag.F90
atmos_lib.o: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/atmos_lib.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/atmos_lib.F90
atmos_nudge.o: $(SRCROOT)atmos_phys/atmos_shared/atmos_nudge/atmos_nudge.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/atmos_nudge/atmos_nudge.F90
atmos_radon.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_radon.F90 atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_radon.F90
atmos_regional_tracer_driver.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_regional_tracer_driver.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_regional_tracer_driver.F90
atmos_sea_salt.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sea_salt.F90 atmos_cmip_diag.o atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sea_salt.F90
atmos_soa.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_soa.F90 atmos_cmip_diag.o xactive_bvoc.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_soa.F90
atmos_sulfate.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sulfate.F90 atmos_cmip_diag.o cloud_chem.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sulfate.F90
atmos_sulfur_hex.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sulfur_hex.F90 atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sulfur_hex.F90
atmos_tracer_driver.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tracer_driver.F90 atmos_cmip_diag.o atmos_tracer_utilities.o atmos_radon.o atmos_carbon_aerosol.o atmos_convection_tracer.o atmos_sulfur_hex.o atmos_ch3i.o atmos_sea_salt.o atmos_dust.o atmos_sulfate.o atmos_soa.o tropchem_driver.o atmos_regional_tracer_driver.o strat_chem_driver.o atmos_age_tracer.o atmos_co2.o atmos_ch4.o atmos_tropopause.o xactive_bvoc.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tracer_driver.F90
atmos_tracer_utilities.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tracer_utilities.F90 atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tracer_utilities.F90
atmos_tropopause.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tropopause.F90 atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tropopause.F90
betaDistribution.o: $(SRCROOT)atmos_phys/atmos_param/cloud_generator/betaDistribution.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/cloud_generator/betaDistribution.F90
betts_miller.o: $(SRCROOT)atmos_phys/atmos_param/betts_miller/betts_miller.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/betts_miller/betts_miller.F90
bm_massflux.o: $(SRCROOT)atmos_phys/atmos_param/betts_miller/bm_massflux.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/betts_miller/bm_massflux.F90
bm_omp.o: $(SRCROOT)atmos_phys/atmos_param/betts_miller/bm_omp.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/betts_miller/bm_omp.F90
bulkphys_rad.o: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/bulkphys_rad.F90 cloudrad_types.o strat_clouds_W.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/clouds/bulkphys_rad.F90
calc_Re.o: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/calc_Re.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H math_lib.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/calc_Re.F90
calendar.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/calendar.F90 clubb_precision.o constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/calendar.F90
cg_drag.o: $(SRCROOT)atmos_phys/atmos_param/cg_drag/cg_drag.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/cg_drag/cg_drag.F90
check_nan.o: $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/check_nan.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/lscloud_driver/check_nan.F90
cldwat2m_micro.o: $(SRCROOT)atmos_phys/atmos_param/microphysics/cldwat2m_micro.F90 gamma_mg.o lscloud_types.o simple_pdf.o physics_radiation_exch.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/microphysics/cldwat2m_micro.F90
clip_explicit.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clip_explicit.F90 grid_class.o parameters_model.o model_flags.o clubb_precision.o stats_type.o stats_variables.o constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clip_explicit.F90
clip_semi_implicit.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clip_semi_implicit.F90 clubb_precision.o constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clip_semi_implicit.F90
cloud_chem.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/cloud_chem.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/cloud_chem.F90
cloud_generator.o: $(SRCROOT)atmos_phys/atmos_param/cloud_generator/cloud_generator.F90 betaDistribution.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/cloud_generator/cloud_generator.F90
cloud_obs.o: $(SRCROOT)atmos_phys/atmos_param/cloud_obs/cloud_obs.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/cloud_obs/cloud_obs.F90
cloud_rad.o: $(SRCROOT)atmos_phys/atmos_param/cloud_rad/cloud_rad.F90 gamma_mg.o lscloud_constants.o aerosol_types.o physics_radiation_exch.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/cloud_rad/cloud_rad.F90
cloud_spec.o: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloud_spec.F90 get_random_number_stream.o aerosol_types.o physics_radiation_exch.o cloudrad_types.o strat_clouds_W.o donner_deep_clouds_W.o uw_clouds_W.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloud_spec.F90
cloud_zonal.o: $(SRCROOT)atmos_phys/atmos_param/cloud_zonal/cloud_zonal.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/cloud_zonal/cloud_zonal.F90
cloudrad_diagnostics.o: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_diagnostics.F90 atmos_cmip_diag.o cloudrad_types.o microphys_rad.o isccp_clouds.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_diagnostics.F90
cloudrad_driver.o: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_driver.F90 aerosol_types.o physics_radiation_exch.o cloudrad_types.o cloud_spec.o cloudrad_package.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_driver.F90
cloudrad_package.o: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_package.F90 cloudrad_types.o cloudrad_diagnostics.o bulkphys_rad.o microphys_rad.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_package.F90
cloudrad_types.o: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_types.F90
clouds.o: $(SRCROOT)atmos_phys/atmos_param/clouds/clouds.F90 cloud_rad.o cloud_zonal.o cloud_obs.o rh_clouds.o diag_cloud.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clouds/clouds.F90
clubb_core.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clubb_core.F90 constants_clubb.o parameters_tunable.o parameters_model.o model_flags.o grid_class.o numerical_check.o variables_diagnostic_module.o pdf_parameter_module.o advance_sclrm_Nd_module.o advance_xm_wpxp_module.o advance_xp2_xpyp_module.o surface_varnce_module.o pdf_closure_module.o mixing_length.o advance_windm_edsclrm_module.o saturation.o advance_wp2_wp3_module.o clubb_precision.o error_code.o Skw_module.o clip_explicit.o T_in_K_module.o stats_subs.o stats_type.o stats_variables.o fill_holes.o sigma_sqd_w_module.o parameter_indices.o variables_prognostic_module.o csr_matrix_class_3array.o gmres_wrap.o gmres_cache.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clubb_core.F90
clubb_precision.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clubb_precision.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clubb_precision.F90
constants_clubb.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/constants_clubb.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/constants_clubb.F90
conv_closures.o: $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_closures.F90 conv_utilities_k.o conv_plumes_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_closures.F90
conv_plumes.o: $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_plumes.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_plumes.F90
conv_plumes_k.o: $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_plumes_k.F90 $(SRCROOT)FMS/include/fms_platform.h aer_ccn_act_k.o conv_utilities_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_plumes_k.F90
conv_utilities.o: $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_utilities.F90 conv_utilities_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_utilities.F90
conv_utilities_k.o: $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_utilities_k.F90 $(SRCROOT)FMS/include/fms_platform.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_utilities_k.F90
convection_driver.o: $(SRCROOT)atmos_phys/atmos_param/convection_driver/convection_driver.F90 physics_types.o vert_diff_driver.o physics_radiation_exch.o betts_miller.o bm_massflux.o bm_omp.o donner_deep.o moist_conv.o uw_conv.o ras.o dry_adj.o detr_ice_num.o rh_clouds.o cu_mo_trans.o mo_hook.o aerosol_types.o moist_processes_utils.o moistproc_kernels.o atmos_tracer_utilities.o atmos_global_diag.o atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/convection_driver/convection_driver.F90
corr_matrix_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/corr_matrix_module.F90 input_reader.o matrix_operations.o constants_clubb.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/corr_matrix_module.F90
cosp.o: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H cosp_types.o cosp_simulator.o cosp_modis_simulator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/cosp.F90
cosp_constants.o: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_constants.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/cosp_constants.F90
cosp_diagnostics.o: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_diagnostics.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H cosp_types.o cosp_io.o cosp_constants.o lmd_ipsl_stats.o cosp_modis_simulator.o modis_simulator.o cosp_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/cosp_diagnostics.F90
cosp_driver.o: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_driver.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H physics_radiation_exch.o physics_types.o cosp_types.o cosp.o cosp_io.o cosp_constants.o cosp_modis_simulator.o cosp_diagnostics.o cosp_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/cosp_driver.F90
cosp_io.o: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_io.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o cosp_modis_simulator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/cosp_io.F90
cosp_isccp_simulator.o: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_isccp_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/cosp_isccp_simulator.F90
cosp_lidar.o: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_lidar.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/cosp_lidar.F90
cosp_misr_simulator.o: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_misr_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/cosp_misr_simulator.F90
cosp_modis_simulator.o: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_modis_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o modis_simulator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/cosp_modis_simulator.F90
cosp_radar.o: $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/cosp_radar.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o cosp_utils.o radar_simulator_types.o array_lib.o atmos_lib.o format_input.o mrgrnk.o math_lib.o optics_lib.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/llnl/cosp_radar.F90
cosp_rttov_simulator.o: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_rttov_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/cosp_rttov_simulator.F90
cosp_simulator.o: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o cosp_radar.o cosp_lidar.o cosp_isccp_simulator.o cosp_modis_simulator.o cosp_misr_simulator.o cosp_rttov_simulator.o cosp_stats.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/cosp_simulator.F90
cosp_stats.o: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_stats.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o llnl_stats.o lmd_ipsl_stats.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/cosp_stats.F90
cosp_types.o: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_types.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_utils.o radar_simulator_types.o scale_LUTs_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/cosp_types.F90
cosp_utils.o: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_utils.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/cosp_utils.F90
csr_matrix_class_3array.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/csr_matrix_class_3array.F90 constants_clubb.o grid_class.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/csr_matrix_class_3array.F90
cu_mo_trans.o: $(SRCROOT)atmos_phys/atmos_param/cu_mo_trans/cu_mo_trans.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/cu_mo_trans/cu_mo_trans.F90
cumulus_closure_k.o: $(SRCROOT)atmos_phys/atmos_param/donner_deep/cumulus_closure_k.F90 donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/donner_deep/cumulus_closure_k.F90
damping_driver.o: $(SRCROOT)atmos_phys/atmos_param/damping_driver/damping_driver.F90 mg_drag.o cg_drag.o topo_drag.o atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/damping_driver/damping_driver.F90
deep_conv.o: $(SRCROOT)atmos_phys/atmos_param/shallow_cu/deep_conv.F90 conv_utilities.o conv_utilities_k.o conv_plumes_k.o conv_closures.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/shallow_cu/deep_conv.F90
detr_ice_num.o: $(SRCROOT)atmos_phys/atmos_param/convection_driver/detr_ice_num.F90 lscloud_constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/convection_driver/detr_ice_num.F90
diag_cloud.o: $(SRCROOT)atmos_phys/atmos_param/diag_cloud/diag_cloud.F90 cloud_zonal.o diag_cloud_rad.o shallow_conv.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/diag_cloud/diag_cloud.F90
diag_cloud_rad.o: $(SRCROOT)atmos_phys/atmos_param/diag_cloud_rad/diag_cloud_rad.F90 cloud_rad.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/diag_cloud_rad/diag_cloud_rad.F90
diffusion.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/diffusion.F90 grid_class.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/diffusion.F90
diffusivity.o: $(SRCROOT)atmos_phys/atmos_param/diffusivity/diffusivity.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/diffusivity/diffusivity.F90
donner_cape_k.o: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_cape_k.F90 donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_cape_k.F90
donner_cloud_model_k.o: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_cloud_model_k.F90 donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_cloud_model_k.F90
donner_deep.o: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep.F90 donner_types.o conv_utilities_k.o conv_plumes_k.o fms_donner.o nonfms_donner.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep.F90
donner_deep_clouds_W.o: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/donner_deep_clouds_W.F90 cloudrad_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/clouds/donner_deep_clouds_W.F90
donner_deep_k.o: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep_k.F90 donner_types.o conv_utilities_k.o conv_plumes_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep_k.F90
donner_deep_miz.o: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep_miz.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep_miz.F90
donner_lite_k.o: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_lite_k.F90 donner_types.o conv_utilities_k.o conv_plumes_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_lite_k.F90
donner_lscloud_k.o: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_lscloud_k.F90 donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_lscloud_k.F90
donner_meso_k.o: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_meso_k.F90 donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_meso_k.F90
donner_rad_k.o: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_rad_k.F90 donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_rad_k.F90
donner_types.o: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_types.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_types.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/donner_deep	$(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_types.F90
donner_utilities_k.o: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_utilities_k.F90 donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_utilities_k.F90
dry_adj.o: $(SRCROOT)atmos_phys/atmos_param/dry_adj/dry_adj.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/dry_adj/dry_adj.F90
dsd.o: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/dsd.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H array_lib.o math_lib.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/dsd.F90
edt.o: $(SRCROOT)atmos_phys/atmos_param/edt/edt.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/edt/edt.F90
endian.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/endian.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/endian.F90
entrain.o: $(SRCROOT)atmos_phys/atmos_param/entrain/entrain.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/entrain/entrain.F90
error_code.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/error_code.F90 constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/error_code.F90
esfsw_bands.o: $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_bands.F90 esfsw_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_bands.F90
esfsw_driver.o: $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_driver.F90 esfsw_parameters.o esfsw_bands.o esfsw_utilities.o shortwave_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_driver.F90
esfsw_parameters.o: $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_parameters.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_parameters.F90
esfsw_utilities.o: $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_utilities.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_utilities.F90
extrapolation.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/extrapolation.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/extrapolation.F90
file_functions.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/file_functions.F90 clubb_precision.o constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/file_functions.F90
fill_holes.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/fill_holes.F90 grid_class.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/fill_holes.F90
fms_donner.o: $(SRCROOT)atmos_phys/atmos_param/donner_deep/fms_donner.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_nml.h atmos_tracer_utilities.o donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/donner_deep	$(SRCROOT)atmos_phys/atmos_param/donner_deep/fms_donner.F90
format_input.o: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/format_input.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H array_lib.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/format_input.F90
gamma_mg.o: $(SRCROOT)atmos_phys/atmos_param/microphysics/gamma_mg.F90 lscloud_constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/microphysics/gamma_mg.F90
gas_tf.o: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/gas_tf.F90 longwave_utilities.o longwave_params.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/longwave/gas_tf.F90
gases.o: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/gases.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/gases.F90
get_random_number_stream.o: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/get_random_number_stream.F90 cloudrad_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/clouds/get_random_number_stream.F90
gmres_cache.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/gmres_cache.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/gmres_cache.F90
gmres_wrap.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/gmres_wrap.F90 gmres_cache.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/gmres_wrap.F90
grey_radiation.o: $(SRCROOT)atmos_phys/atmos_param/grey_radiation/grey_radiation.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/grey_radiation/grey_radiation.F90
grid_class.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/grid_class.F90 clubb_precision.o constants_clubb.o error_code.o file_functions.o interpolation.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/grid_class.F90
hydrostatic_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/hydrostatic_module.F90 constants_clubb.o grid_class.o clubb_precision.o interpolation.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/hydrostatic_module.F90
hyper_diffusion_4th_ord.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/hyper_diffusion_4th_ord.F90 clubb_precision.o grid_class.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/hyper_diffusion_4th_ord.F90
icarus.o: $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/icarus.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/icarus.F90
ice_nucl.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90 aer_ccn_act_k.o aerosol_params.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90
input_names.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/input_names.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/input_names.F90
input_reader.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/input_reader.F90 clubb_precision.o constants_clubb.o input_names.o interpolation.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/input_reader.F90
interpolation.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/interpolation.F90 clubb_precision.o constants_clubb.o model_flags.o error_code.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/interpolation.F90
isccp_cloud_types.o: $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F
isccp_clouds.o: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/isccp_clouds.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/clouds/isccp_clouds.F90
isoropiaIIcode.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/isoropiaIIcode.F $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/isrpia.inc
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/isoropiaIIcode.F
lapack_wrap.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/lapack_wrap.F90 constants_clubb.o error_code.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/lapack_wrap.F90
lidar_simulator.o: $(SRCROOT)atmos_phys/atmos_param/cosp/actsim/lidar_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/actsim/lidar_simulator.F90
lin_cloud_microphys.o: $(SRCROOT)atmos_phys/atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90
llnl_stats.o: $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/llnl_stats.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/llnl/llnl_stats.F90
lmd_ipsl_stats.o: $(SRCROOT)atmos_phys/atmos_param/cosp/actsim/lmd_ipsl_stats.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H llnl_stats.o cosp_constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/actsim/lmd_ipsl_stats.F90
longwave_clouds.o: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_clouds.F90 longwave_utilities.o longwave_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_clouds.F90
longwave_driver.o: $(SRCROOT)atmos_phys/atmos_param/radiation/driver/longwave_driver.F90 radiation_driver_types.o aerosolrad_types.o sealw99.o longwave_types.o radiative_gases.o radiative_gases_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/driver/longwave_driver.F90
longwave_fluxes.o: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_fluxes.F90 longwave_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_fluxes.F90
longwave_params.o: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_params.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_params.F90
longwave_tables.o: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_tables.F90 longwave_utilities.o longwave_params.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_tables.F90
longwave_types.o: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_types.F90
longwave_utilities.o: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_utilities.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_utilities.F90
ls_cloud_macrophysics.o: $(SRCROOT)atmos_phys/atmos_param/macrophysics/ls_cloud_macrophysics.F90 physics_types.o tiedtke_macro.o aerosol_cloud.o lscloud_types.o CLUBB_driver_SCM.o aerosol_types.o moist_processes_utils.o physics_radiation_exch.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/macrophysics/ls_cloud_macrophysics.F90
ls_cloud_microphysics.o: $(SRCROOT)atmos_phys/atmos_param/microphysics/ls_cloud_microphysics.F90 physics_types.o lscloud_types.o aerosol_types.o physics_radiation_exch.o moist_processes_utils.o lin_cloud_microphys.o lscloud_debug.o rotstayn_klein_mp.o morrison_gettelman_microp.o cldwat2m_micro.o micro_mg.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/microphysics/ls_cloud_microphysics.F90
lscale_cond.o: $(SRCROOT)atmos_phys/atmos_param/lscale_cond/lscale_cond.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/lscale_cond/lscale_cond.F90
lscloud_constants.o: $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_constants.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_constants.F90
lscloud_debug.o: $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_debug.F90 lscloud_types.o moist_processes_utils.o check_nan.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_debug.F90
lscloud_driver.o: $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_driver.F90 physics_types.o lscale_cond.o lscloud_debug.o lscloud_types.o polysvp.o lscloud_netcdf.o rh_clouds.o physics_radiation_exch.o ls_cloud_macrophysics.o ls_cloud_microphysics.o aerosol_cloud.o strat_cloud.o moist_processes_utils.o aerosol_types.o atmos_tracer_utilities.o atmos_global_diag.o atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_driver.F90
lscloud_netcdf.o: $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_netcdf.F90 lscloud_types.o atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_netcdf.F90
lscloud_types.o: $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_types.F90 atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_types.F90
lw_gases_stdtf.o: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/lw_gases_stdtf.F90 longwave_utilities.o gas_tf.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/longwave/lw_gases_stdtf.F90
m_tracname.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/m_tracname.F90 moz.mods.o AM3.mods.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/m_tracname.F90
math_lib.o: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/math_lib.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H mrgrnk.o array_lib.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/math_lib.F90
matrix_operations.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/matrix_operations.F90 clubb_precision.o error_code.o constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/matrix_operations.F90
mean_adv.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mean_adv.F90 grid_class.o model_flags.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mean_adv.F90
mg_drag.o: $(SRCROOT)atmos_phys/atmos_param/mg_drag/mg_drag.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/mg_drag/mg_drag.F90
micro_mg.o: $(SRCROOT)atmos_phys/atmos_param/microphysics/micro_mg.F90 gamma_mg.o lscloud_types.o simple_pdf.o physics_radiation_exch.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/microphysics/micro_mg.F90
microphys_rad.o: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/microphys_rad.F90 cloudrad_types.o sealw99.o esfsw_driver.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/clouds/microphys_rad.F90
mixing_length.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mixing_length.F90 constants_clubb.o parameters_tunable.o parameters_model.o grid_class.o numerical_check.o saturation.o error_code.o model_flags.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mixing_length.F90
mo_chem_utls.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90 moz.mods.o AM3.mods.o m_tracname.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90
mo_chemdr.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chemdr.F90 tropchem_types_mod.o mo_chem_utls.o moz.mods.o AM3.mods.o mo_photo.o mo_exp_slv.o mo_imp_slv.o mo_rodas_slv.o mo_usrrxt.o mo_setinv.o moz.subs.o AM3.subs.o mo_setsox.o mo_fphoto.o AM3_fphoto.o strat_chem_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chemdr.F90
mo_chemini.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chemini.F90 mo_photo.o mo_chem_utls.o mo_usrrxt.o moz.mods.o AM3.mods.o mo_exp_slv.o mo_imp_slv.o mo_rodas_slv.o mo_read_sim_chm.o mo_fphoto.o AM3_fphoto.o tropchem_types_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chemini.F90
mo_exp_slv.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90 moz.mods.o AM3.mods.o mo_chem_utls.o moz.mat.o AM3.mat.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90
mo_fastjx.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fastjx.F90 tropchem_types_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fastjx.F90
mo_fphoto.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.F90 mo_fastjx.o mo_chem_utls.o moz.mods.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.F90
mo_hook.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_hook.F90 atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_hook.F90
mo_imp_slv.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90 moz.mods.o AM3.mods.o mo_chem_utls.o m_tracname.o moz.mat.o AM3.mat.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90
mo_jpl.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_jpl.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_jpl.F90
mo_photo.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_photo.F90 mo_chem_utls.o moz.mods.o AM3.mods.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_photo.F90
mo_read_sim_chm.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90 moz.mods.o AM3.mods.o m_tracname.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90
mo_rodas_slv.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90 moz.mods.o AM3.mods.o mo_chem_utls.o moz.mat.o AM3.mat.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90
mo_setinv.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_setinv.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_setinv.F90
mo_setsox.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/mo_setsox.F90 mo_chem_utls.o tropchem_types_mod.o cloud_chem.o aerosol_thermodynamics.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/mo_setsox.F90
mo_usrrxt.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90 aerosol_thermodynamics.o strat_chem_utilities.o tropchem_types_mod.o mo_chem_utls.o moz.mods.o AM3.mods.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90
model_flags.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/model_flags.F90 constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/model_flags.F90
modis_simulator.o: $(SRCROOT)atmos_phys/atmos_param/cosp/MODIS_simulator/modis_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H cosp_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/MODIS_simulator/modis_simulator.F90
moist_conv.o: $(SRCROOT)atmos_phys/atmos_param/moist_conv/moist_conv.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/moist_conv/moist_conv.F90
moist_processes.o: $(SRCROOT)atmos_phys/atmos_param/moist_processes/moist_processes.F90 physics_types.o physics_radiation_exch.o lscloud_driver.o convection_driver.o atmos_global_diag.o vert_diff_driver.o aerosol_types.o moist_processes_utils.o atmos_tracer_utilities.o atmos_dust.o atmos_tracer_driver.o atmos_sea_salt.o atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/moist_processes/moist_processes.F90
moist_processes_utils.o: $(SRCROOT)atmos_phys/atmos_param/moist_processes/moist_processes_utils.F90 lscloud_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/moist_processes/moist_processes_utils.F90
moistproc_kernels.o: $(SRCROOT)atmos_phys/atmos_param/convection_driver/moistproc_kernels.F90 moist_conv.o uw_conv.o ras.o cu_mo_trans.o aerosol_types.o detr_ice_num.o physics_radiation_exch.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/convection_driver/moistproc_kernels.F90
mono_flux_limiter.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mono_flux_limiter.F90 grid_class.o constants_clubb.o clubb_precision.o error_code.o fill_holes.o stats_type.o stats_variables.o mean_adv.o lapack_wrap.o anl_erf.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mono_flux_limiter.F90
morrison_gettelman_microp.o: $(SRCROOT)atmos_phys/atmos_param/microphysics/morrison_gettelman_microp.F90 polysvp.o lscloud_debug.o lscloud_types.o gamma_mg.o lscloud_constants.o simple_pdf.o physics_radiation_exch.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/microphysics/morrison_gettelman_microp.F90
moz.mat.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.F90 moz.mods.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.F90
moz.mods.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.F90
moz.subs.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.F90 moz.mods.o mo_jpl.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.F90
mrgrnk.o: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/mrgrnk.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/mrgrnk.F90
mt95.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mt95.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mt95.f90
my25_turb.o: $(SRCROOT)atmos_phys/atmos_param/my25_turb/my25_turb.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/my25_turb/my25_turb.F90
nonfms_donner.o: $(SRCROOT)atmos_phys/atmos_param/donner_deep/nonfms_donner.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_nml.h donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/donner_deep	$(SRCROOT)atmos_phys/atmos_param/donner_deep/nonfms_donner.F90
numerical_check.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/numerical_check.F90 grid_class.o clubb_precision.o parameters_model.o pdf_parameter_module.o stats_variables.o constants_clubb.o variables_diagnostic_module.o variables_prognostic_module.o parameters_microphys.o error_code.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/numerical_check.F90
optical_path.o: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/optical_path.F90 longwave_utilities.o longwave_params.o lw_gases_stdtf.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/longwave/optical_path.F90
optics_lib.o: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/optics_lib.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/optics_lib.F90
output_grads.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/output_grads.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/recl.inc clubb_precision.o constants_clubb.o stat_file_module.o model_flags.o endian.o calendar.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/output_grads.F90
output_netcdf.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/output_netcdf.F90 stat_file_module.o clubb_precision.o constants_clubb.o parameters_model.o parameters_tunable.o parameter_indices.o model_flags.o parameters_microphys.o parameters_radiation.o calendar.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/output_netcdf.F90
ozone.o: $(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/ozone.F90 radiative_gases_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/ozone.F90
parameter_indices.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameter_indices.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameter_indices.F90
parameters_microphys.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_microphys.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_microphys.F90
parameters_model.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_model.F90 clubb_precision.o constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_model.F90
parameters_radiation.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_radiation.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_radiation.F90
parameters_tunable.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_tunable.F90 parameter_indices.o grid_class.o clubb_precision.o constants_clubb.o error_code.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_tunable.F90
pdf_closure_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pdf_closure_module.F90 constants_clubb.o parameters_model.o parameters_tunable.o pdf_parameter_module.o anl_erf.o numerical_check.o saturation.o error_code.o stats_variables.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pdf_closure_module.F90
pdf_parameter_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pdf_parameter_module.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pdf_parameter_module.F90
pf_to_mr.o: $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/pf_to_mr.F $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/llnl/pf_to_mr.F
physics_driver.o: $(SRCROOT)atmos_phys/atmos_param/physics_driver/physics_driver.F90 atmos_tracer_driver.o atmos_cmip_diag.o aerosol_types.o physics_radiation_exch.o physics_types.o moist_processes_utils.o aerosol.o cosp_driver.o moist_processes.o vert_turb_driver.o vert_diff_driver.o damping_driver.o grey_radiation.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/physics_driver/physics_driver.F90
physics_radiation_exch.o: $(SRCROOT)atmos_phys/atmos_param/physics_radiation_exch/physics_radiation_exch.F90 $(SRCROOT)FMS/include/fms_platform.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_phys/atmos_param/physics_radiation_exch/physics_radiation_exch.F90
physics_types.o: $(SRCROOT)atmos_phys/atmos_param/physics_driver/physics_types.F90 $(SRCROOT)FMS/include/fms_platform.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_phys/atmos_param/physics_driver/physics_types.F90
polysvp.o: $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/polysvp.F90 lscloud_constants.o lscloud_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/lscloud_driver/polysvp.F90
pos_definite_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pos_definite_module.F90 grid_class.o constants_clubb.o clubb_precision.o error_code.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pos_definite_module.F90
prec_scops.o: $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/prec_scops.F $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/llnl/prec_scops.F
qe_moist_convection.o: $(SRCROOT)atmos_phys/atmos_param/qe_moist_convection/qe_moist_convection.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/qe_moist_convection/qe_moist_convection.F90
rad_output_file.o: $(SRCROOT)atmos_phys/atmos_param/radiation/driver/rad_output_file.F90 aerosol_types.o atmos_cmip_diag.o radiation_driver_types.o aerosolrad_types.o cloudrad_types.o esfsw_driver.o shortwave_types.o longwave_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/driver/rad_output_file.F90
radar_simulator.o: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H mrgrnk.o array_lib.o math_lib.o optics_lib.o radar_simulator_types.o scale_LUTs_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator.F90
radar_simulator_init.o: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator_init.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H radar_simulator_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator_init.F90
radar_simulator_types.o: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator_types.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator_types.F90
radiation_diag.o: $(SRCROOT)atmos_phys/atmos_param/radiation/util/radiation_diag.F90 radiation_driver.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/util/radiation_diag.F90
radiation_driver.o: $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver.F90 radiation_types.o physics_radiation_exch.o radiation_driver_types.o aerosol_types.o aerosolrad_types.o cloudrad_types.o longwave_driver.o longwave_types.o shortwave_driver.o shortwave_types.o rad_output_file.o aerosolrad_driver.o cloudrad_driver.o cloudrad_diagnostics.o radiative_gases.o radiative_gases_types.o radiation_driver_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver.F90
radiation_driver_diag.o: $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver_diag.F90 atmos_cmip_diag.o atmos_global_diag.o radiation_driver_types.o aerosolrad_types.o shortwave_types.o longwave_types.o radiative_gases_types.o shortwave_driver.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver_diag.F90
radiation_driver_types.o: $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver_types.F90
radiation_types.o: $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_types.F90 $(SRCROOT)FMS/include/fms_platform.h atmos_co2.o atmos_ch4.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_types.F90
radiative_gases.o: $(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/radiative_gases.F90 sealw99.o radiative_gases_types.o ozone.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/radiative_gases.F90
radiative_gases_types.o: $(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/radiative_gases_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/radiative_gases_types.F90
ras.o: $(SRCROOT)atmos_phys/atmos_param/ras/ras.F90 aerosol_types.o aer_ccn_act.o moist_processes_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/ras/ras.F90
rh_clouds.o: $(SRCROOT)atmos_phys/atmos_param/rh_clouds/rh_clouds.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/rh_clouds/rh_clouds.F90
rotstayn_klein_mp.o: $(SRCROOT)atmos_phys/atmos_param/microphysics/rotstayn_klein_mp.F90 cloud_generator.o aer_in_act.o polysvp.o lscloud_types.o lscloud_debug.o physics_radiation_exch.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/microphysics/rotstayn_klein_mp.F90
saturation.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/saturation.F90 model_flags.o clubb_precision.o constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/saturation.F90
scale_LUTs_io.o: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/scale_LUTs_io.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H radar_simulator_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/scale_LUTs_io.F90
scops.o: $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/scops.F $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/congvec.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp -I$(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7	$(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/scops.F
sealw99.o: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/sealw99.F90 longwave_utilities.o longwave_types.o longwave_params.o longwave_clouds.o longwave_fluxes.o longwave_tables.o optical_path.o gas_tf.o lw_gases_stdtf.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/longwave/sealw99.F90
shallow_conv.o: $(SRCROOT)atmos_phys/atmos_param/shallow_conv/shallow_conv.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/shallow_conv/shallow_conv.F90
shallow_physics.o: $(SRCROOT)atmos_phys/atmos_param/shallow_physics/shallow_physics.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/shallow_physics/shallow_physics.F90
shortwave_driver.o: $(SRCROOT)atmos_phys/atmos_param/radiation/driver/shortwave_driver.F90 radiation_driver_types.o esfsw_driver.o shortwave_types.o radiative_gases_types.o solar_data_driver.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/driver/shortwave_driver.F90
shortwave_types.o: $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/shortwave_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/shortwave_types.F90
sigma_sqd_w_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/sigma_sqd_w_module.F90 constants_clubb.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/sigma_sqd_w_module.F90
simple_pdf.o: $(SRCROOT)atmos_phys/atmos_param/microphysics/simple_pdf.F90 betaDistribution.o lscloud_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/microphysics/simple_pdf.F90
solar_data_driver.o: $(SRCROOT)atmos_phys/atmos_param/radiation/driver/solar_data_driver.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/driver/solar_data_driver.F90
sponge_layer_damping.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/sponge_layer_damping.F90 clubb_precision.o grid_class.o constants_clubb.o interpolation.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/sponge_layer_damping.F90
stable_bl_turb.o: $(SRCROOT)atmos_phys/atmos_param/stable_bl_turb/stable_bl_turb.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/stable_bl_turb/stable_bl_turb.F90
stat_file_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stat_file_module.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stat_file_module.F90
stats_LH_zt.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_LH_zt.F90 constants_clubb.o stats_variables.o stats_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_LH_zt.F90
stats_rad_zm.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_rad_zm.F90 constants_clubb.o stats_variables.o stats_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_rad_zm.F90
stats_rad_zt.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_rad_zt.F90 constants_clubb.o stats_variables.o stats_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_rad_zt.F90
stats_sfc.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_sfc.F90 constants_clubb.o stats_variables.o stats_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_sfc.F90
stats_subs.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_subs.F90 stats_variables.o clubb_precision.o output_grads.o output_netcdf.o stats_zm.o stats_zt.o stats_LH_zt.o stats_rad_zt.o stats_rad_zm.o stats_sfc.o error_code.o constants_clubb.o parameters_microphys.o grid_class.o variables_diagnostic_module.o pdf_parameter_module.o T_in_K_module.o parameters_model.o stats_type.o fill_holes.o interpolation.o saturation.o array_index.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_subs.F90
stats_type.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_type.F90 stat_file_module.o clubb_precision.o grid_class.o error_code.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_type.F90
stats_variables.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_variables.F90 stats_type.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_variables.F90
stats_zm.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_zm.F90 constants_clubb.o stats_variables.o stats_type.o parameters_model.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_zm.F90
stats_zt.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_zt.F90 constants_clubb.o stats_variables.o stats_type.o parameters_model.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_zt.F90
strat_chem_driver.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90 strat_chem_model.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90
strat_chem_model.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/stratchem/strat_chem_model.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/stratchem/strat_chem_model.F90
strat_chem_utilities.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90
strat_cloud.o: $(SRCROOT)atmos_phys/atmos_param/strat_cloud/strat_cloud.F90 aerosol_types.o lscloud_types.o moist_processes_utils.o physics_types.o physics_radiation_exch.o lscloud_netcdf.o aerosol_cloud.o betaDistribution.o aer_in_act.o aer_ccn_act.o cloud_generator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/strat_cloud/strat_cloud.F90
strat_clouds_W.o: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/strat_clouds_W.F90 aerosol_types.o physics_radiation_exch.o cloud_rad.o cloudrad_types.o get_random_number_stream.o cloud_generator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/clouds/strat_clouds_W.F90
surface_varnce_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/surface_varnce_module.F90 parameters_model.o constants_clubb.o numerical_check.o error_code.o array_index.o stats_type.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/surface_varnce_module.F90
tiedtke_macro.o: $(SRCROOT)atmos_phys/atmos_param/macrophysics/tiedtke_macro.F90 betaDistribution.o lscloud_types.o lscloud_debug.o polysvp.o moist_processes_utils.o physics_radiation_exch.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/macrophysics/tiedtke_macro.F90
tke_turb.o: $(SRCROOT)atmos_phys/atmos_param/tke_turb/tke_turb.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/tke_turb/tke_turb.F90
topo_drag.o: $(SRCROOT)atmos_phys/atmos_param/topo_drag/topo_drag.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/topo_drag/topo_drag.F90
tropchem_driver.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem_driver.F90 tropchem_types_mod.o atmos_cmip_diag.o atmos_tracer_utilities.o mo_chemdr.o mo_setsox.o mo_chemini.o m_tracname.o moz.mods.o AM3.mods.o mo_hook.o strat_chem_utilities.o mo_chem_utls.o atmos_sulfate.o shortwave_driver.o cloud_chem.o aerosol_thermodynamics.o mo_usrrxt.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem_driver.F90
tropchem_types_mod.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.F90
two_stream_gray_rad.o: $(SRCROOT)atmos_phys/atmos_param/two_stream_gray_rad/two_stream_gray_rad.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/two_stream_gray_rad/two_stream_gray_rad.F90
uw_clouds_W.o: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/uw_clouds_W.F90 physics_radiation_exch.o cloudrad_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/radiation/clouds/uw_clouds_W.F90
uw_conv.o: $(SRCROOT)atmos_phys/atmos_param/shallow_cu/uw_conv.F90 $(SRCROOT)FMS/include/fms_platform.h atmos_cmip_diag.o atmos_tracer_utilities.o moist_processes_utils.o aerosol_types.o aer_ccn_act.o conv_utilities.o conv_utilities_k.o conv_plumes_k.o conv_closures.o deep_conv.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_phys/atmos_param/shallow_cu/uw_conv.F90
variables_diagnostic_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_diagnostic_module.F90 pdf_parameter_module.o clubb_precision.o constants_clubb.o parameters_model.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_diagnostic_module.F90
variables_prognostic_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_prognostic_module.F90 pdf_parameter_module.o clubb_precision.o constants_clubb.o parameters_model.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_prognostic_module.F90
variables_radiation_module.o: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_radiation_module.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_radiation_module.F90
vert_advection.o: $(SRCROOT)atmos_phys/atmos_shared/vert_advection/vert_advection.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/vert_advection/vert_advection.F90
vert_diff.o: $(SRCROOT)atmos_phys/atmos_param/vert_diff/vert_diff.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/vert_diff/vert_diff.F90
vert_diff_driver.o: $(SRCROOT)atmos_phys/atmos_param/vert_diff_driver/vert_diff_driver.F90 vert_diff.o atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/vert_diff_driver/vert_diff_driver.F90
vert_turb_driver.o: $(SRCROOT)atmos_phys/atmos_param/vert_turb_driver/vert_turb_driver.F90 my25_turb.o tke_turb.o diffusivity.o edt.o physics_radiation_exch.o physics_types.o shallow_conv.o stable_bl_turb.o entrain.o moist_processes_utils.o atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/vert_turb_driver/vert_turb_driver.F90
wet_deposition_0D.o: $(SRCROOT)atmos_phys/atmos_param/donner_deep/wet_deposition_0D.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_param/donner_deep/wet_deposition_0D.F90
xactive_bvoc.o: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/xactive_bvoc.F90 moz.mods.o m_tracname.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_phys/atmos_shared/tracer_driver/xactive_bvoc.F90
zeff.o: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/zeff.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H math_lib.o optics_lib.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_phys/atmos_param/cosp	$(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/zeff.F90
./atmos_cmip_interp.inc: $(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_interp.inc
	cp $(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_interp.inc .
./atmos_convection_tracer.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_convection_tracer.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_convection_tracer.F90 .
./T_in_K_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/T_in_K_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/T_in_K_module.F90 .
./check_nan.F90: $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/check_nan.F90
	cp $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/check_nan.F90 .
./AM3_fphoto.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3_fphoto.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3_fphoto.F90 .
./mo_hook.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_hook.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_hook.F90 .
./rh_clouds.F90: $(SRCROOT)atmos_phys/atmos_param/rh_clouds/rh_clouds.F90
	cp $(SRCROOT)atmos_phys/atmos_param/rh_clouds/rh_clouds.F90 .
./sealw99.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/sealw99.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/sealw99.F90 .
./uw_clouds_W.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/uw_clouds_W.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/uw_clouds_W.F90 .
./recl.inc: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/recl.inc
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/recl.inc .
./mo_read_sim_chm.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90 .
./strat_chem_model.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/stratchem/strat_chem_model.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/stratchem/strat_chem_model.F90 .
./array_index.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/array_index.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/array_index.F90 .
./csr_matrix_class_3array.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/csr_matrix_class_3array.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/csr_matrix_class_3array.F90 .
./constants_clubb.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/constants_clubb.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/constants_clubb.F90 .
./betts_miller.F90: $(SRCROOT)atmos_phys/atmos_param/betts_miller/betts_miller.F90
	cp $(SRCROOT)atmos_phys/atmos_param/betts_miller/betts_miller.F90 .
./parameters_tunable.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_tunable.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_tunable.F90 .
./variables_radiation_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_radiation_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_radiation_module.F90 .
./clubb_precision.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clubb_precision.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clubb_precision.F90 .
./cosp_stats.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_stats.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_stats.F90 .
./aerosol_types.F90: $(SRCROOT)atmos_phys/atmos_shared/aerosol/aerosol_types.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/aerosol/aerosol_types.F90 .
./parameters_radiation.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_radiation.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_radiation.F90 .
./aer_ccn_act.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90 .
./mean_adv.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mean_adv.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mean_adv.F90 .
./atmos_radon.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_radon.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_radon.F90 .
./lidar_simulator.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/actsim/lidar_simulator.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/actsim/lidar_simulator.F90 .
./cloud_zonal.F90: $(SRCROOT)atmos_phys/atmos_param/cloud_zonal/cloud_zonal.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cloud_zonal/cloud_zonal.F90 .
./entrain.F90: $(SRCROOT)atmos_phys/atmos_param/entrain/entrain.F90
	cp $(SRCROOT)atmos_phys/atmos_param/entrain/entrain.F90 .
./variables_diagnostic_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_diagnostic_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_diagnostic_module.F90 .
./mt95.f90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mt95.f90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mt95.f90 .
./radar_simulator.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator.F90 .
./strat_chem_utilities.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90 .
./pf_to_mr.F: $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/pf_to_mr.F
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/pf_to_mr.F .
./sigma_sqd_w_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/sigma_sqd_w_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/sigma_sqd_w_module.F90 .
./morrison_gettelman_microp.F90: $(SRCROOT)atmos_phys/atmos_param/microphysics/morrison_gettelman_microp.F90
	cp $(SRCROOT)atmos_phys/atmos_param/microphysics/morrison_gettelman_microp.F90 .
./format_input.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/format_input.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/format_input.F90 .
./cloud_chem.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/cloud_chem.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/cloud_chem.F90 .
./shortwave_types.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/shortwave_types.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/shortwave_types.F90 .
./rotstayn_klein_mp.F90: $(SRCROOT)atmos_phys/atmos_param/microphysics/rotstayn_klein_mp.F90
	cp $(SRCROOT)atmos_phys/atmos_param/microphysics/rotstayn_klein_mp.F90 .
./parameters_microphys.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_microphys.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_microphys.F90 .
./lscloud_driver.F90: $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_driver.F90 .
./mo_exp_slv.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90 .
./bm_massflux.F90: $(SRCROOT)atmos_phys/atmos_param/betts_miller/bm_massflux.F90
	cp $(SRCROOT)atmos_phys/atmos_param/betts_miller/bm_massflux.F90 .
./conv_utilities_k.F90: $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_utilities_k.F90
	cp $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_utilities_k.F90 .
./moist_processes.F90: $(SRCROOT)atmos_phys/atmos_param/moist_processes/moist_processes.F90
	cp $(SRCROOT)atmos_phys/atmos_param/moist_processes/moist_processes.F90 .
./calendar.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/calendar.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/calendar.F90 .
./gamma_mg.F90: $(SRCROOT)atmos_phys/atmos_param/microphysics/gamma_mg.F90
	cp $(SRCROOT)atmos_phys/atmos_param/microphysics/gamma_mg.F90 .
./atmos_regional_tracer_driver.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_regional_tracer_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_regional_tracer_driver.F90 .
./stats_subs.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_subs.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_subs.F90 .
./calc_Re.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/calc_Re.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/calc_Re.F90 .
./longwave_types.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_types.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_types.F90 .
./ice_nucl.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90 .
./moz.mods.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.F90 .
./pdf_closure_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pdf_closure_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pdf_closure_module.F90 .
./mo_imp_slv.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90 .
./donner_deep_k.F90: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep_k.F90
	cp $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep_k.F90 .
./advance_sclrm_Nd_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/advance_sclrm_Nd_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/advance_sclrm_Nd_module.F90 .
./conv_plumes_k.F90: $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_plumes_k.F90
	cp $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_plumes_k.F90 .
./interpolation.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/interpolation.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/interpolation.F90 .
./congvec.H: $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/congvec.H
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/congvec.H .
./mo_chemini.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chemini.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chemini.F90 .
./mg_drag.F90: $(SRCROOT)atmos_phys/atmos_param/mg_drag/mg_drag.F90
	cp $(SRCROOT)atmos_phys/atmos_param/mg_drag/mg_drag.F90 .
./pos_definite_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pos_definite_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pos_definite_module.F90 .
./conv_plumes.F90: $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_plumes.F90
	cp $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_plumes.F90 .
./xactive_bvoc.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/xactive_bvoc.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/xactive_bvoc.F90 .
./advance_xm_wpxp_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_xm_wpxp_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_xm_wpxp_module.F90 .
./atmos_tracer_utilities.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tracer_utilities.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tracer_utilities.F90 .
./radiative_gases_types.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/radiative_gases_types.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/radiative_gases_types.F90 .
./physics_types.F90: $(SRCROOT)atmos_phys/atmos_param/physics_driver/physics_types.F90
	cp $(SRCROOT)atmos_phys/atmos_param/physics_driver/physics_types.F90 .
./stats_type.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_type.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_type.F90 .
./diag_cloud.F90: $(SRCROOT)atmos_phys/atmos_param/diag_cloud/diag_cloud.F90
	cp $(SRCROOT)atmos_phys/atmos_param/diag_cloud/diag_cloud.F90 .
./atmos_sulfur_hex.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sulfur_hex.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sulfur_hex.F90 .
./donner_lscloud_k.F90: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_lscloud_k.F90
	cp $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_lscloud_k.F90 .
./donner_deep_miz.F90: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep_miz.F90
	cp $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep_miz.F90 .
./stats_LH_zt.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_LH_zt.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_LH_zt.F90 .
./cloudrad_package.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_package.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_package.F90 .
./atmos_global_diag.F90: $(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag/atmos_global_diag.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag/atmos_global_diag.F90 .
./hyper_diffusion_4th_ord.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/hyper_diffusion_4th_ord.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/hyper_diffusion_4th_ord.F90 .
./aerosol_thermodynamics.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.F90 .
./donner_nml.h: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_nml.h
	cp $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_nml.h .
./cg_drag.F90: $(SRCROOT)atmos_phys/atmos_param/cg_drag/cg_drag.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cg_drag/cg_drag.F90 .
./alt_cloud.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/alt_cloud.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/alt_cloud.F90 .
./ras.F90: $(SRCROOT)atmos_phys/atmos_param/ras/ras.F90
	cp $(SRCROOT)atmos_phys/atmos_param/ras/ras.F90 .
./aerosol_cloud.F90: $(SRCROOT)atmos_phys/atmos_param/aerosol_cloud/aerosol_cloud.F90
	cp $(SRCROOT)atmos_phys/atmos_param/aerosol_cloud/aerosol_cloud.F90 .
./lscloud_netcdf.F90: $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_netcdf.F90
	cp $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_netcdf.F90 .
./two_stream_gray_rad.F90: $(SRCROOT)atmos_phys/atmos_param/two_stream_gray_rad/two_stream_gray_rad.F90
	cp $(SRCROOT)atmos_phys/atmos_param/two_stream_gray_rad/two_stream_gray_rad.F90 .
./atmos_co2.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_co2.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_co2.F90 .
./moz.mat.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.F90 .
./fill_holes.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/fill_holes.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/fill_holes.F90 .
./ls_cloud_macrophysics.F90: $(SRCROOT)atmos_phys/atmos_param/macrophysics/ls_cloud_macrophysics.F90
	cp $(SRCROOT)atmos_phys/atmos_param/macrophysics/ls_cloud_macrophysics.F90 .
./esfsw_driver.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_driver.F90 .
./longwave_clouds.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_clouds.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_clouds.F90 .
./cloud_spec.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloud_spec.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloud_spec.F90 .
./clip_semi_implicit.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clip_semi_implicit.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clip_semi_implicit.F90 .
./vert_advection.F90: $(SRCROOT)atmos_phys/atmos_shared/vert_advection/vert_advection.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/vert_advection/vert_advection.F90 .
./stats_zt.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_zt.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_zt.F90 .
./aer_ccn_act_k.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90 .
./donner_types.F90: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_types.F90
	cp $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_types.F90 .
./radiation_types.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_types.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_types.F90 .
./bulkphys_rad.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/bulkphys_rad.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/bulkphys_rad.F90 .
./anl_erf.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/anl_erf.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/anl_erf.F90 .
./rad_output_file.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/driver/rad_output_file.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/driver/rad_output_file.F90 .
./mo_photo.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_photo.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_photo.F90 .
./scale_LUTs_io.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/scale_LUTs_io.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/scale_LUTs_io.F90 .
./lscale_cond.F90: $(SRCROOT)atmos_phys/atmos_param/lscale_cond/lscale_cond.F90
	cp $(SRCROOT)atmos_phys/atmos_param/lscale_cond/lscale_cond.F90 .
./cloudrad_diagnostics.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_diagnostics.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_diagnostics.F90 .
./cosp_defs.H: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H .
./extrapolation.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/extrapolation.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/extrapolation.F90 .
./surface_varnce_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/surface_varnce_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/surface_varnce_module.F90 .
./array_lib.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/array_lib.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/array_lib.F90 .
./radiation_diag.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/util/radiation_diag.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/util/radiation_diag.F90 .
./sponge_layer_damping.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/sponge_layer_damping.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/sponge_layer_damping.F90 .
./cumulus_closure_k.F90: $(SRCROOT)atmos_phys/atmos_param/donner_deep/cumulus_closure_k.F90
	cp $(SRCROOT)atmos_phys/atmos_param/donner_deep/cumulus_closure_k.F90 .
./radiation_driver.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver.F90 .
./lapack_wrap.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/lapack_wrap.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/lapack_wrap.F90 .
./cloud_obs.F90: $(SRCROOT)atmos_phys/atmos_param/cloud_obs/cloud_obs.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cloud_obs/cloud_obs.F90 .
./cosp_constants.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_constants.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_constants.F90 .
./donner_deep_clouds_W.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/donner_deep_clouds_W.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/donner_deep_clouds_W.F90 .
./deep_conv.F90: $(SRCROOT)atmos_phys/atmos_param/shallow_cu/deep_conv.F90
	cp $(SRCROOT)atmos_phys/atmos_param/shallow_cu/deep_conv.F90 .
./polysvp.F90: $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/polysvp.F90
	cp $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/polysvp.F90 .
./cloud_rad.F90: $(SRCROOT)atmos_phys/atmos_param/cloud_rad/cloud_rad.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cloud_rad/cloud_rad.F90 .
./prec_scops.F: $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/prec_scops.F
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/prec_scops.F .
./longwave_utilities.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_utilities.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_utilities.F90 .
./model_flags.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/model_flags.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/model_flags.F90 .
./clouds.F90: $(SRCROOT)atmos_phys/atmos_param/clouds/clouds.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clouds/clouds.F90 .
./conv_utilities.F90: $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_utilities.F90
	cp $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_utilities.F90 .
./lscloud_debug.F90: $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_debug.F90
	cp $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_debug.F90 .
./vert_turb_driver.F90: $(SRCROOT)atmos_phys/atmos_param/vert_turb_driver/vert_turb_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_param/vert_turb_driver/vert_turb_driver.F90 .
./MISR_simulator.F: $(SRCROOT)atmos_phys/atmos_param/cosp/MISR_simulator/MISR_simulator.F
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/MISR_simulator/MISR_simulator.F .
./parameter_indices.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameter_indices.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameter_indices.F90 .
./CLUBB_driver_SCM.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_driver_SCM.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_driver_SCM.F90 .
./lw_gases_stdtf.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/lw_gases_stdtf.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/lw_gases_stdtf.F90 .
./shallow_conv.F90: $(SRCROOT)atmos_phys/atmos_param/shallow_conv/shallow_conv.F90
	cp $(SRCROOT)atmos_phys/atmos_param/shallow_conv/shallow_conv.F90 .
./mono_flux_limiter.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mono_flux_limiter.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mono_flux_limiter.F90 .
./mo_fastjx.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fastjx.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fastjx.F90 .
./mo_usrrxt.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90 .
./wet_deposition_0D.F90: $(SRCROOT)atmos_phys/atmos_param/donner_deep/wet_deposition_0D.F90
	cp $(SRCROOT)atmos_phys/atmos_param/donner_deep/wet_deposition_0D.F90 .
./donner_lite_k.F90: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_lite_k.F90
	cp $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_lite_k.F90 .
./moist_processes_utils.F90: $(SRCROOT)atmos_phys/atmos_param/moist_processes/moist_processes_utils.F90
	cp $(SRCROOT)atmos_phys/atmos_param/moist_processes/moist_processes_utils.F90 .
./math_lib.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/math_lib.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/math_lib.F90 .
./grid_class.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/grid_class.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/grid_class.F90 .
./atmos_soa.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_soa.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_soa.F90 .
./cosp_radar.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/cosp_radar.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/cosp_radar.F90 .
./cosp_simulator.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_simulator.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_simulator.F90 .
./endian.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/endian.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/endian.F90 .
./conv_closures.F90: $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_closures.F90
	cp $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_closures.F90 .
./tiedtke_macro.F90: $(SRCROOT)atmos_phys/atmos_param/macrophysics/tiedtke_macro.F90
	cp $(SRCROOT)atmos_phys/atmos_param/macrophysics/tiedtke_macro.F90 .
./esfsw_parameters.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_parameters.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_parameters.F90 .
./mo_setinv.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_setinv.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_setinv.F90 .
./cloudrad_driver.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_driver.F90 .
./cu_mo_trans.F90: $(SRCROOT)atmos_phys/atmos_param/cu_mo_trans/cu_mo_trans.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cu_mo_trans/cu_mo_trans.F90 .
./stats_sfc.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_sfc.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_sfc.F90 .
./isoropiaIIcode.F: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/isoropiaIIcode.F
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/isoropiaIIcode.F .
./m_tracname.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/m_tracname.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/m_tracname.F90 .
./mrgrnk.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/mrgrnk.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/mrgrnk.F90 .
./mo_jpl.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_jpl.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_jpl.F90 .
./aerosolrad_package.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_package.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_package.F90 .
./advance_helper_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_helper_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_helper_module.F90 .
./atmos_tropopause.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tropopause.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tropopause.F90 .
./AM3.subs.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.subs.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.subs.F90 .
./diag_cloud_rad.F90: $(SRCROOT)atmos_phys/atmos_param/diag_cloud_rad/diag_cloud_rad.F90
	cp $(SRCROOT)atmos_phys/atmos_param/diag_cloud_rad/diag_cloud_rad.F90 .
./atmos_tracer_driver.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tracer_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tracer_driver.F90 .
./lmd_ipsl_stats.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/actsim/lmd_ipsl_stats.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/actsim/lmd_ipsl_stats.F90 .
./micro_mg.F90: $(SRCROOT)atmos_phys/atmos_param/microphysics/micro_mg.F90
	cp $(SRCROOT)atmos_phys/atmos_param/microphysics/micro_mg.F90 .
./saturation.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/saturation.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/saturation.F90 .
./AM3.mods.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.mods.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.mods.F90 .
./cosp_isccp_simulator.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_isccp_simulator.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_isccp_simulator.F90 .
./hydrostatic_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/hydrostatic_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/hydrostatic_module.F90 .
./stable_bl_turb.F90: $(SRCROOT)atmos_phys/atmos_param/stable_bl_turb/stable_bl_turb.F90
	cp $(SRCROOT)atmos_phys/atmos_param/stable_bl_turb/stable_bl_turb.F90 .
./radiation_driver_diag.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver_diag.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver_diag.F90 .
./donner_utilities_k.F90: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_utilities_k.F90
	cp $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_utilities_k.F90 .
./mo_chemdr.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chemdr.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chemdr.F90 .
./edt.F90: $(SRCROOT)atmos_phys/atmos_param/edt/edt.F90
	cp $(SRCROOT)atmos_phys/atmos_param/edt/edt.F90 .
./microphys_rad.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/microphys_rad.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/microphys_rad.F90 .
./radar_simulator_types.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator_types.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator_types.F90 .
./CLUBB_3D_var.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_3D_var.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_3D_var.F90 .
./mo_rodas_slv.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90 .
./tropchem_driver.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem_driver.F90 .
./simple_pdf.F90: $(SRCROOT)atmos_phys/atmos_param/microphysics/simple_pdf.F90
	cp $(SRCROOT)atmos_phys/atmos_param/microphysics/simple_pdf.F90 .
./radiative_gases.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/radiative_gases.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/radiative_gases.F90 .
./cloudrad_types.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_types.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_types.F90 .
./aer_in_act.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90 .
./fms_donner.F90: $(SRCROOT)atmos_phys/atmos_param/donner_deep/fms_donner.F90
	cp $(SRCROOT)atmos_phys/atmos_param/donner_deep/fms_donner.F90 .
./esfsw_bands.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_bands.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_bands.F90 .
./dsd.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/dsd.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/dsd.F90 .
./atmos_sea_salt.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sea_salt.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sea_salt.F90 .
./cldwat2m_micro.F90: $(SRCROOT)atmos_phys/atmos_param/microphysics/cldwat2m_micro.F90
	cp $(SRCROOT)atmos_phys/atmos_param/microphysics/cldwat2m_micro.F90 .
./shallow_physics.F90: $(SRCROOT)atmos_phys/atmos_param/shallow_physics/shallow_physics.F90
	cp $(SRCROOT)atmos_phys/atmos_param/shallow_physics/shallow_physics.F90 .
./detr_ice_num.F90: $(SRCROOT)atmos_phys/atmos_param/convection_driver/detr_ice_num.F90
	cp $(SRCROOT)atmos_phys/atmos_param/convection_driver/detr_ice_num.F90 .
./qe_moist_convection.F90: $(SRCROOT)atmos_phys/atmos_param/qe_moist_convection/qe_moist_convection.F90
	cp $(SRCROOT)atmos_phys/atmos_param/qe_moist_convection/qe_moist_convection.F90 .
./atmos_ch3i.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_ch3i.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_ch3i.F90 .
./betaDistribution.F90: $(SRCROOT)atmos_phys/atmos_param/cloud_generator/betaDistribution.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cloud_generator/betaDistribution.F90 .
./tke_turb.F90: $(SRCROOT)atmos_phys/atmos_param/tke_turb/tke_turb.F90
	cp $(SRCROOT)atmos_phys/atmos_param/tke_turb/tke_turb.F90 .
./vert_diff_driver.F90: $(SRCROOT)atmos_phys/atmos_param/vert_diff_driver/vert_diff_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_param/vert_diff_driver/vert_diff_driver.F90 .
./icarus.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/icarus.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/icarus.F90 .
./optics_lib.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/optics_lib.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/optics_lib.F90 .
./bm_omp.F90: $(SRCROOT)atmos_phys/atmos_param/betts_miller/bm_omp.F90
	cp $(SRCROOT)atmos_phys/atmos_param/betts_miller/bm_omp.F90 .
./output_grads.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/output_grads.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/output_grads.F90 .
./strat_cloud.F90: $(SRCROOT)atmos_phys/atmos_param/strat_cloud/strat_cloud.F90
	cp $(SRCROOT)atmos_phys/atmos_param/strat_cloud/strat_cloud.F90 .
./lscloud_constants.F90: $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_constants.F90
	cp $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_constants.F90 .
./aerosolrad_types.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_types.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_types.F90 .
./mo_chem_utls.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90 .
./optical_path.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/optical_path.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/optical_path.F90 .
./get_random_number_stream.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/get_random_number_stream.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/get_random_number_stream.F90 .
./aerosol_params.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90 .
./gmres_cache.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/gmres_cache.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/gmres_cache.F90 .
./atmos_lib.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/atmos_lib.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/atmos_lib.F90 .
./isrpia.inc: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/isrpia.inc
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/isrpia.inc .
./dry_adj.F90: $(SRCROOT)atmos_phys/atmos_param/dry_adj/dry_adj.F90
	cp $(SRCROOT)atmos_phys/atmos_param/dry_adj/dry_adj.F90 .
./radar_simulator_init.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator_init.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator_init.F90 .
./input_reader.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/input_reader.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/input_reader.F90 .
./nonfms_donner.F90: $(SRCROOT)atmos_phys/atmos_param/donner_deep/nonfms_donner.F90
	cp $(SRCROOT)atmos_phys/atmos_param/donner_deep/nonfms_donner.F90 .
./pdf_parameter_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pdf_parameter_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pdf_parameter_module.F90 .
./isccp_clouds.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/isccp_clouds.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/isccp_clouds.F90 .
./advance_windm_edsclrm_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_windm_edsclrm_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_windm_edsclrm_module.F90 .
./longwave_fluxes.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_fluxes.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_fluxes.F90 .
./file_functions.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/file_functions.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/file_functions.F90 .
./advance_wp2_wp3_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_wp2_wp3_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_wp2_wp3_module.F90 .
./atmos_age_tracer.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_age_tracer.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_age_tracer.F90 .
./gas_tf.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/gas_tf.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/gas_tf.F90 .
./error_code.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/error_code.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/error_code.F90 .
./lin_cloud_microphys.F90: $(SRCROOT)atmos_phys/atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90
	cp $(SRCROOT)atmos_phys/atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90 .
./cosp_utils.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_utils.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_utils.F90 .
./donner_meso_k.F90: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_meso_k.F90
	cp $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_meso_k.F90 .
./lscloud_types.F90: $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_types.F90
	cp $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_types.F90 .
./clip_explicit.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clip_explicit.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clip_explicit.F90 .
./AM3_fastjx.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3_fastjx.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3_fastjx.F90 .
./stats_variables.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_variables.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_variables.F90 .
./moistproc_kernels.F90: $(SRCROOT)atmos_phys/atmos_param/convection_driver/moistproc_kernels.F90
	cp $(SRCROOT)atmos_phys/atmos_param/convection_driver/moistproc_kernels.F90 .
./donner_cape_k.F90: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_cape_k.F90
	cp $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_cape_k.F90 .
./cosp_driver.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_driver.F90 .
./donner_cloud_model_k.F90: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_cloud_model_k.F90
	cp $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_cloud_model_k.F90 .
./numerical_check.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/numerical_check.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/numerical_check.F90 .
./stat_file_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stat_file_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stat_file_module.F90 .
./input_names.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/input_names.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/input_names.F90 .
./atmos_sulfate.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sulfate.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sulfate.F90 .
./strat_clouds_W.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/strat_clouds_W.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/strat_clouds_W.F90 .
./cosp_diagnostics.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_diagnostics.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_diagnostics.F90 .
./AM3.mat.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.mat.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.mat.F90 .
./gases.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/gases.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/gases.F90 .
./grey_radiation.F90: $(SRCROOT)atmos_phys/atmos_param/grey_radiation/grey_radiation.F90
	cp $(SRCROOT)atmos_phys/atmos_param/grey_radiation/grey_radiation.F90 .
./longwave_tables.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_tables.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_tables.F90 .
./strat_chem_driver.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90 .
./moist_conv.F90: $(SRCROOT)atmos_phys/atmos_param/moist_conv/moist_conv.F90
	cp $(SRCROOT)atmos_phys/atmos_param/moist_conv/moist_conv.F90 .
./atmos_dust.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_dust.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_dust.F90 .
./ozone.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/ozone.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/ozone.F90 .
./ls_cloud_microphysics.F90: $(SRCROOT)atmos_phys/atmos_param/microphysics/ls_cloud_microphysics.F90
	cp $(SRCROOT)atmos_phys/atmos_param/microphysics/ls_cloud_microphysics.F90 .
./cosp_lidar.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_lidar.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_lidar.F90 .
./corr_matrix_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/corr_matrix_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/corr_matrix_module.F90 .
./advance_xp2_xpyp_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_xp2_xpyp_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_xp2_xpyp_module.F90 .
./longwave_driver.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/driver/longwave_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/driver/longwave_driver.F90 .
./cosp_misr_simulator.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_misr_simulator.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_misr_simulator.F90 .
./mo_setsox.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/mo_setsox.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/mo_setsox.F90 .
./donner_deep.F90: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep.F90
	cp $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep.F90 .
./isccp_cloud_types.F: $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F .
./llnl_stats.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/llnl_stats.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/llnl_stats.F90 .
./cosp.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/cosp.F90 .
./mixing_length.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mixing_length.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mixing_length.F90 .
./cosp_rttov_simulator.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_rttov_simulator.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_rttov_simulator.F90 .
./modis_simulator.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/MODIS_simulator/modis_simulator.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/MODIS_simulator/modis_simulator.F90 .
./gmres_wrap.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/gmres_wrap.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/gmres_wrap.F90 .
./solar_data_driver.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/driver/solar_data_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/driver/solar_data_driver.F90 .
./diffusion.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/diffusion.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/diffusion.F90 .
./zeff.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/zeff.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/zeff.F90 .
./atmos_cmip_diag.F90: $(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.F90 .
./fms_platform.h: $(SRCROOT)FMS/include/fms_platform.h
	cp $(SRCROOT)FMS/include/fms_platform.h .
./topo_drag.F90: $(SRCROOT)atmos_phys/atmos_param/topo_drag/topo_drag.F90
	cp $(SRCROOT)atmos_phys/atmos_param/topo_drag/topo_drag.F90 .
./atmos_carbon_aerosol.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_carbon_aerosol.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_carbon_aerosol.F90 .
./physics_radiation_exch.F90: $(SRCROOT)atmos_phys/atmos_param/physics_radiation_exch/physics_radiation_exch.F90
	cp $(SRCROOT)atmos_phys/atmos_param/physics_radiation_exch/physics_radiation_exch.F90 .
./aerosolrad_driver.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_driver.F90 .
./uw_conv.F90: $(SRCROOT)atmos_phys/atmos_param/shallow_cu/uw_conv.F90
	cp $(SRCROOT)atmos_phys/atmos_param/shallow_cu/uw_conv.F90 .
./stats_rad_zt.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_rad_zt.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_rad_zt.F90 .
./physics_driver.F90: $(SRCROOT)atmos_phys/atmos_param/physics_driver/physics_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_param/physics_driver/physics_driver.F90 .
./damping_driver.F90: $(SRCROOT)atmos_phys/atmos_param/damping_driver/damping_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_param/damping_driver/damping_driver.F90 .
./donner_rad_k.F90: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_rad_k.F90
	cp $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_rad_k.F90 .
./scops.F: $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/scops.F
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/scops.F .
./radiation_driver_types.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver_types.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver_types.F90 .
./moz.subs.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.F90 .
./vert_diff.F90: $(SRCROOT)atmos_phys/atmos_param/vert_diff/vert_diff.F90
	cp $(SRCROOT)atmos_phys/atmos_param/vert_diff/vert_diff.F90 .
./mo_fphoto.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.F90 .
./stats_zm.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_zm.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_zm.F90 .
./esfsw_utilities.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_utilities.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_utilities.F90 .
./shortwave_driver.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/driver/shortwave_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/driver/shortwave_driver.F90 .
./my25_turb.F90: $(SRCROOT)atmos_phys/atmos_param/my25_turb/my25_turb.F90
	cp $(SRCROOT)atmos_phys/atmos_param/my25_turb/my25_turb.F90 .
./clubb_core.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clubb_core.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clubb_core.F90 .
./variables_prognostic_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_prognostic_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_prognostic_module.F90 .
./atmos_ch4.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_ch4.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_ch4.F90 .
./aerosol.F90: $(SRCROOT)atmos_phys/atmos_shared/aerosol/aerosol.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/aerosol/aerosol.F90 .
./tropchem_types_mod.F90: $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.F90 .
./Skw_module.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/Skw_module.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/Skw_module.F90 .
./longwave_params.F90: $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_params.F90
	cp $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_params.F90 .
./diffusivity.F90: $(SRCROOT)atmos_phys/atmos_param/diffusivity/diffusivity.F90
	cp $(SRCROOT)atmos_phys/atmos_param/diffusivity/diffusivity.F90 .
./cloud_generator.F90: $(SRCROOT)atmos_phys/atmos_param/cloud_generator/cloud_generator.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cloud_generator/cloud_generator.F90 .
./cosp_io.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_io.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_io.F90 .
./cosp_modis_simulator.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_modis_simulator.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_modis_simulator.F90 .
./matrix_operations.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/matrix_operations.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/matrix_operations.F90 .
./stats_rad_zm.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_rad_zm.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_rad_zm.F90 .
./atmos_nudge.F90: $(SRCROOT)atmos_phys/atmos_shared/atmos_nudge/atmos_nudge.F90
	cp $(SRCROOT)atmos_phys/atmos_shared/atmos_nudge/atmos_nudge.F90 .
./convection_driver.F90: $(SRCROOT)atmos_phys/atmos_param/convection_driver/convection_driver.F90
	cp $(SRCROOT)atmos_phys/atmos_param/convection_driver/convection_driver.F90 .
./cosp_types.F90: $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_types.F90
	cp $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_types.F90 .
./output_netcdf.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/output_netcdf.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/output_netcdf.F90 .
./donner_types.h: $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_types.h
	cp $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_types.h .
./parameters_model.F90: $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_model.F90
	cp $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_model.F90 .
SRC = $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_diagnostic_module.F90 $(SRCROOT)atmos_phys/atmos_param/microphysics/gamma_mg.F90 $(SRCROOT)atmos_phys/atmos_shared/aerosol/aerosol_types.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/Skw_module.F90 $(SRCROOT)atmos_phys/atmos_param/clouds/clouds.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_types.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_age_tracer.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_tunable.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_rad_zt.F90 $(SRCROOT)atmos_phys/atmos_param/vert_diff/vert_diff.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/hydrostatic_module.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/mo_setsox.F90 $(SRCROOT)atmos_phys/atmos_param/cloud_generator/betaDistribution.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_zt.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tropopause.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_xm_wpxp_module.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/corr_matrix_module.F90 $(SRCROOT)atmos_phys/atmos_param/betts_miller/bm_massflux.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/util/radiation_diag.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_utilities_k.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/calc_Re.F90 $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_types.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_driver.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/matrix_operations.F90 $(SRCROOT)atmos_phys/atmos_shared/aerosol/aerosol.F90 $(SRCROOT)atmos_phys/atmos_param/moist_processes/moist_processes_utils.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_types.F90 $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_constants.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sea_salt.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_utilities.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/grid_class.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pdf_parameter_module.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_cape_k.F90 $(SRCROOT)atmos_phys/atmos_param/diag_cloud/diag_cloud.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_jpl.F90 $(SRCROOT)atmos_phys/atmos_shared/vert_advection/vert_advection.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_meso_k.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/gmres_cache.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/scops.F $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/optics_lib.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/array_index.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_cu/deep_conv.F90 $(SRCROOT)atmos_phys/atmos_param/convection_driver/moistproc_kernels.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/dsd.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sulfate.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/csr_matrix_class_3array.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_plumes_k.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/constants_clubb.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_conv/shallow_conv.F90 $(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag/atmos_global_diag.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_regional_tracer_driver.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_prognostic_module.F90 $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/polysvp.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/get_random_number_stream.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_bands.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mt95.f90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_diagnostics.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3_fastjx.F90 $(SRCROOT)atmos_phys/atmos_param/stable_bl_turb/stable_bl_turb.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tracer_driver.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_LH_zt.F90 $(SRCROOT)atmos_phys/atmos_param/convection_driver/convection_driver.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/xactive_bvoc.F90 $(SRCROOT)atmos_phys/atmos_param/diag_cloud_rad/diag_cloud_rad.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/strat_clouds_W.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_io.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_cu/uw_conv.F90 $(SRCROOT)atmos_phys/atmos_param/topo_drag/topo_drag.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_package.F90 $(SRCROOT)atmos_phys/atmos_param/physics_driver/physics_types.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_variables.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/driver/rad_output_file.F90 $(SRCROOT)atmos_phys/atmos_param/moist_processes/moist_processes.F90 $(SRCROOT)atmos_phys/atmos_param/my25_turb/my25_turb.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/math_lib.F90 $(SRCROOT)atmos_phys/atmos_param/microphysics/simple_pdf.F90 $(SRCROOT)atmos_phys/atmos_param/dry_adj/dry_adj.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/output_netcdf.F90 $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/check_nan.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90 $(SRCROOT)atmos_phys/atmos_param/microphysics/ls_cloud_microphysics.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/output_grads.F90 $(SRCROOT)atmos_phys/atmos_param/tke_turb/tke_turb.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clip_explicit.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_radiation.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_dust.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fastjx.F90 $(SRCROOT)atmos_phys/atmos_param/physics_radiation_exch/physics_radiation_exch.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clip_semi_implicit.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/sealw99.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_wp2_wp3_module.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/stratchem/strat_chem_model.F90 $(SRCROOT)atmos_phys/atmos_param/entrain/entrain.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/driver/longwave_driver.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clubb_precision.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_types.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_radiation_module.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stat_file_module.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_model.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_constants.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mono_flux_limiter.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chemini.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem_driver.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/sponge_layer_damping.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/alt_cloud.F90 $(SRCROOT)atmos_phys/atmos_param/betts_miller/bm_omp.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_package.F90 $(SRCROOT)atmos_phys/atmos_param/moist_conv/moist_conv.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_xp2_xpyp_module.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/interpolation.F90 $(SRCROOT)atmos_phys/atmos_param/ras/ras.F90 $(SRCROOT)atmos_phys/atmos_param/microphysics/morrison_gettelman_microp.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/bulkphys_rad.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_physics/shallow_physics.F90 $(SRCROOT)atmos_phys/atmos_param/lscale_cond/lscale_cond.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/atmos_lib.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_microphys.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_lscloud_k.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/sigma_sqd_w_module.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/hyper_diffusion_4th_ord.F90 $(SRCROOT)atmos_phys/atmos_param/damping_driver/damping_driver.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameter_indices.F90 $(SRCROOT)atmos_phys/atmos_param/mg_drag/mg_drag.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator_types.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/pf_to_mr.F $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_ch3i.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/gases.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver_diag.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/gmres_wrap.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_tables.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_cloud_model_k.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/radiative_gases_types.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/microphys_rad.F90 $(SRCROOT)atmos_phys/atmos_param/microphysics/micro_mg.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90 $(SRCROOT)atmos_phys/atmos_param/aerosol_cloud/aerosol_cloud.F90 $(SRCROOT)atmos_phys/atmos_param/edt/edt.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_driver.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_hook.F90 $(SRCROOT)atmos_phys/atmos_param/microphysics/rotstayn_klein_mp.F90 $(SRCROOT)atmos_phys/atmos_param/cloud_rad/cloud_rad.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/radiative_gases.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/isoropiaIIcode.F $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/zeff.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_soa.F90 $(SRCROOT)atmos_phys/atmos_param/vert_turb_driver/vert_turb_driver.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/cloud_chem.F90 $(SRCROOT)atmos_phys/atmos_param/rh_clouds/rh_clouds.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_driver_SCM.F90 $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_netcdf.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/saturation.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pdf_closure_module.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_types.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/mrgrnk.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/wet_deposition_0D.F90 $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_debug.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_fluxes.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_driver.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_type.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_parameters.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/fms_donner.F90 $(SRCROOT)atmos_phys/atmos_param/cloud_zonal/cloud_zonal.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_convection_tracer.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_closures.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/cosp_radar.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_ch4.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/optical_path.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_utilities.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/fill_holes.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.F90 $(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/anl_erf.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_isccp_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/cu_mo_trans/cu_mo_trans.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mixing_length.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/isccp_clouds.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/array_lib.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/m_tracname.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator_init.F90 $(SRCROOT)atmos_phys/atmos_param/betts_miller/betts_miller.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_carbon_aerosol.F90 $(SRCROOT)atmos_phys/atmos_param/macrophysics/ls_cloud_macrophysics.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_co2.F90 $(SRCROOT)atmos_phys/atmos_param/two_stream_gray_rad/two_stream_gray_rad.F90 $(SRCROOT)atmos_phys/atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/numerical_check.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/gas_tf.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/lw_gases_stdtf.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_lidar.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/endian.F90 $(SRCROOT)atmos_phys/atmos_param/qe_moist_convection/qe_moist_convection.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_utils.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/MODIS_simulator/modis_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/MISR_simulator/MISR_simulator.F $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_params.F90 $(SRCROOT)atmos_phys/atmos_param/cg_drag/cg_drag.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mean_adv.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clubb_core.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep_miz.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_modis_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/donner_deep_clouds_W.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_rad_k.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/ozone.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_simulator.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90 $(SRCROOT)atmos_phys/atmos_param/cloud_generator/cloud_generator.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/advance_sclrm_Nd_module.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/diffusion.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/llnl_stats.F90 $(SRCROOT)atmos_phys/atmos_param/grey_radiation/grey_radiation.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_misr_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/uw_clouds_W.F90 $(SRCROOT)atmos_phys/atmos_param/macrophysics/tiedtke_macro.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pos_definite_module.F90 $(SRCROOT)atmos_phys/atmos_param/cloud_obs/cloud_obs.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/calendar.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/T_in_K_module.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_diagnostics.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/error_code.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tracer_utilities.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/surface_varnce_module.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_utilities.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/lapack_wrap.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/input_reader.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/file_functions.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_lite_k.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.mat.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep_k.F90 $(SRCROOT)atmos_phys/atmos_param/microphysics/cldwat2m_micro.F90 $(SRCROOT)atmos_phys/atmos_param/physics_driver/physics_driver.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_3D_var.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/cumulus_closure_k.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/prec_scops.F $(SRCROOT)atmos_phys/atmos_param/radiation/driver/solar_data_driver.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_sfc.F90 $(SRCROOT)atmos_phys/atmos_param/convection_driver/detr_ice_num.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_subs.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloud_spec.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/driver/shortwave_driver.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver_types.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/model_flags.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chemdr.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_driver.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_helper_module.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_rttov_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/nonfms_donner.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_setinv.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_types.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/scale_LUTs_io.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/shortwave_types.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_rad_zm.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_clouds.F90 $(SRCROOT)atmos_phys/atmos_shared/atmos_nudge/atmos_nudge.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_photo.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90 $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_driver.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.subs.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/extrapolation.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/input_names.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3_fphoto.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/icarus.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_zm.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_radon.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_types.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_stats.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.mods.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sulfur_hex.F90 $(SRCROOT)atmos_phys/atmos_param/diffusivity/diffusivity.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_windm_edsclrm_module.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_utilities_k.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_plumes.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/format_input.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/actsim/lmd_ipsl_stats.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/actsim/lidar_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/strat_cloud/strat_cloud.F90 $(SRCROOT)atmos_phys/atmos_param/vert_diff_driver/vert_diff_driver.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_nml.h $(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_interp.inc $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/recl.inc $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_types.h $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/congvec.H $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/isrpia.inc
OBJ = variables_diagnostic_module.o gamma_mg.o aerosol_types.o Skw_module.o clouds.o cosp_types.o atmos_age_tracer.o parameters_tunable.o stats_rad_zt.o vert_diff.o hydrostatic_module.o mo_setsox.o betaDistribution.o stats_zt.o atmos_tropopause.o advance_xm_wpxp_module.o corr_matrix_module.o bm_massflux.o radiation_diag.o conv_utilities_k.o calc_Re.o lscloud_types.o cloudrad_driver.o matrix_operations.o aerosol.o moist_processes_utils.o cloudrad_types.o lscloud_constants.o atmos_sea_salt.o longwave_utilities.o grid_class.o pdf_parameter_module.o donner_deep.o donner_cape_k.o diag_cloud.o mo_jpl.o vert_advection.o tropchem_types_mod.o donner_meso_k.o gmres_cache.o scops.o optics_lib.o array_index.o deep_conv.o moistproc_kernels.o dsd.o atmos_sulfate.o csr_matrix_class_3array.o conv_plumes_k.o constants_clubb.o shallow_conv.o atmos_global_diag.o atmos_regional_tracer_driver.o variables_prognostic_module.o polysvp.o get_random_number_stream.o esfsw_bands.o aerosol_thermodynamics.o mt95.o cloudrad_diagnostics.o AM3_fastjx.o stable_bl_turb.o atmos_tracer_driver.o stats_LH_zt.o convection_driver.o isccp_cloud_types.o xactive_bvoc.o diag_cloud_rad.o strat_clouds_W.o cosp_io.o uw_conv.o topo_drag.o aerosolrad_package.o physics_types.o stats_variables.o rad_output_file.o moist_processes.o my25_turb.o math_lib.o simple_pdf.o dry_adj.o output_netcdf.o check_nan.o ice_nucl.o ls_cloud_microphysics.o output_grads.o tke_turb.o clip_explicit.o parameters_radiation.o atmos_dust.o mo_fastjx.o physics_radiation_exch.o clip_semi_implicit.o sealw99.o advance_wp2_wp3_module.o strat_chem_model.o entrain.o aerosol_params.o longwave_driver.o clubb_precision.o longwave_types.o variables_radiation_module.o stat_file_module.o parameters_model.o cosp_constants.o mono_flux_limiter.o mo_chemini.o tropchem_driver.o sponge_layer_damping.o alt_cloud.o bm_omp.o mo_chem_utls.o cloudrad_package.o moist_conv.o aer_ccn_act_k.o advance_xp2_xpyp_module.o interpolation.o ras.o morrison_gettelman_microp.o bulkphys_rad.o shallow_physics.o lscale_cond.o atmos_lib.o parameters_microphys.o donner_lscloud_k.o sigma_sqd_w_module.o hyper_diffusion_4th_ord.o damping_driver.o parameter_indices.o mg_drag.o radar_simulator_types.o pf_to_mr.o atmos_ch3i.o gases.o radiation_driver_diag.o gmres_wrap.o longwave_tables.o donner_cloud_model_k.o radiative_gases_types.o microphys_rad.o micro_mg.o mo_read_sim_chm.o aerosol_cloud.o edt.o esfsw_driver.o mo_hook.o rotstayn_klein_mp.o cloud_rad.o cosp.o moz.subs.o radiative_gases.o mo_rodas_slv.o isoropiaIIcode.o zeff.o atmos_soa.o vert_turb_driver.o cloud_chem.o rh_clouds.o CLUBB_driver_SCM.o lscloud_netcdf.o saturation.o pdf_closure_module.o mo_usrrxt.o donner_types.o mrgrnk.o wet_deposition_0D.o lscloud_debug.o longwave_fluxes.o aerosolrad_driver.o stats_type.o esfsw_parameters.o fms_donner.o cloud_zonal.o atmos_convection_tracer.o conv_closures.o radar_simulator.o cosp_radar.o atmos_ch4.o optical_path.o conv_utilities.o mo_fphoto.o fill_holes.o moz.mods.o atmos_cmip_diag.o anl_erf.o cosp_isccp_simulator.o cu_mo_trans.o mixing_length.o isccp_clouds.o array_lib.o m_tracname.o radar_simulator_init.o betts_miller.o atmos_carbon_aerosol.o ls_cloud_macrophysics.o atmos_co2.o two_stream_gray_rad.o lin_cloud_microphys.o numerical_check.o aer_in_act.o gas_tf.o lw_gases_stdtf.o cosp_lidar.o moz.mat.o endian.o qe_moist_convection.o cosp_utils.o modis_simulator.o MISR_simulator.o longwave_params.o cg_drag.o strat_chem_driver.o mean_adv.o clubb_core.o donner_deep_miz.o cosp_modis_simulator.o donner_deep_clouds_W.o donner_rad_k.o ozone.o cosp_simulator.o aer_ccn_act.o cloud_generator.o advance_sclrm_Nd_module.o diffusion.o llnl_stats.o grey_radiation.o cosp_misr_simulator.o uw_clouds_W.o tiedtke_macro.o pos_definite_module.o cloud_obs.o calendar.o T_in_K_module.o cosp_diagnostics.o error_code.o atmos_tracer_utilities.o surface_varnce_module.o esfsw_utilities.o lapack_wrap.o input_reader.o file_functions.o mo_imp_slv.o donner_lite_k.o AM3.mat.o donner_deep_k.o cldwat2m_micro.o physics_driver.o CLUBB_3D_var.o cumulus_closure_k.o prec_scops.o solar_data_driver.o stats_sfc.o detr_ice_num.o stats_subs.o cloud_spec.o shortwave_driver.o radiation_driver_types.o model_flags.o mo_chemdr.o cosp_driver.o mo_exp_slv.o advance_helper_module.o cosp_rttov_simulator.o nonfms_donner.o mo_setinv.o aerosolrad_types.o scale_LUTs_io.o shortwave_types.o stats_rad_zm.o longwave_clouds.o atmos_nudge.o radiation_driver.o mo_photo.o strat_chem_utilities.o lscloud_driver.o AM3.subs.o extrapolation.o input_names.o AM3_fphoto.o icarus.o stats_zm.o atmos_radon.o radiation_types.o cosp_stats.o AM3.mods.o atmos_sulfur_hex.o diffusivity.o advance_windm_edsclrm_module.o donner_utilities_k.o conv_plumes.o format_input.o lmd_ipsl_stats.o lidar_simulator.o strat_cloud.o vert_diff_driver.o
OFF = $(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_interp.inc $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_convection_tracer.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/T_in_K_module.F90 $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/check_nan.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3_fphoto.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_hook.F90 $(SRCROOT)atmos_phys/atmos_param/rh_clouds/rh_clouds.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/sealw99.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/uw_clouds_W.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/recl.inc $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/stratchem/strat_chem_model.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/array_index.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/csr_matrix_class_3array.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/constants_clubb.F90 $(SRCROOT)atmos_phys/atmos_param/betts_miller/betts_miller.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_tunable.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_radiation_module.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clubb_precision.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_stats.F90 $(SRCROOT)atmos_phys/atmos_shared/aerosol/aerosol_types.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_radiation.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mean_adv.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_radon.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/actsim/lidar_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/cloud_zonal/cloud_zonal.F90 $(SRCROOT)atmos_phys/atmos_param/entrain/entrain.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_diagnostic_module.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mt95.f90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/pf_to_mr.F $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/sigma_sqd_w_module.F90 $(SRCROOT)atmos_phys/atmos_param/microphysics/morrison_gettelman_microp.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/format_input.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/cloud_chem.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/shortwave_types.F90 $(SRCROOT)atmos_phys/atmos_param/microphysics/rotstayn_klein_mp.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_microphys.F90 $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_driver.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90 $(SRCROOT)atmos_phys/atmos_param/betts_miller/bm_massflux.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_utilities_k.F90 $(SRCROOT)atmos_phys/atmos_param/moist_processes/moist_processes.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/calendar.F90 $(SRCROOT)atmos_phys/atmos_param/microphysics/gamma_mg.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_regional_tracer_driver.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_subs.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/calc_Re.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_types.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pdf_closure_module.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep_k.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/advance_sclrm_Nd_module.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_plumes_k.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/interpolation.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/congvec.H $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chemini.F90 $(SRCROOT)atmos_phys/atmos_param/mg_drag/mg_drag.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pos_definite_module.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_plumes.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/xactive_bvoc.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_xm_wpxp_module.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tracer_utilities.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/radiative_gases_types.F90 $(SRCROOT)atmos_phys/atmos_param/physics_driver/physics_types.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_type.F90 $(SRCROOT)atmos_phys/atmos_param/diag_cloud/diag_cloud.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sulfur_hex.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_lscloud_k.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep_miz.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_LH_zt.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_package.F90 $(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag/atmos_global_diag.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/hyper_diffusion_4th_ord.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_nml.h $(SRCROOT)atmos_phys/atmos_param/cg_drag/cg_drag.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/alt_cloud.F90 $(SRCROOT)atmos_phys/atmos_param/ras/ras.F90 $(SRCROOT)atmos_phys/atmos_param/aerosol_cloud/aerosol_cloud.F90 $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_netcdf.F90 $(SRCROOT)atmos_phys/atmos_param/two_stream_gray_rad/two_stream_gray_rad.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_co2.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/fill_holes.F90 $(SRCROOT)atmos_phys/atmos_param/macrophysics/ls_cloud_macrophysics.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_driver.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_clouds.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloud_spec.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clip_semi_implicit.F90 $(SRCROOT)atmos_phys/atmos_shared/vert_advection/vert_advection.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_zt.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_types.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_types.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/bulkphys_rad.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/anl_erf.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/driver/rad_output_file.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_photo.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/scale_LUTs_io.F90 $(SRCROOT)atmos_phys/atmos_param/lscale_cond/lscale_cond.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_diagnostics.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_defs.H $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/extrapolation.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/surface_varnce_module.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/array_lib.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/util/radiation_diag.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/sponge_layer_damping.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/cumulus_closure_k.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/lapack_wrap.F90 $(SRCROOT)atmos_phys/atmos_param/cloud_obs/cloud_obs.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_constants.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/donner_deep_clouds_W.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_cu/deep_conv.F90 $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/polysvp.F90 $(SRCROOT)atmos_phys/atmos_param/cloud_rad/cloud_rad.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/prec_scops.F $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_utilities.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/model_flags.F90 $(SRCROOT)atmos_phys/atmos_param/clouds/clouds.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_utilities.F90 $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_debug.F90 $(SRCROOT)atmos_phys/atmos_param/vert_turb_driver/vert_turb_driver.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/MISR_simulator/MISR_simulator.F $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameter_indices.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_driver_SCM.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/lw_gases_stdtf.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_conv/shallow_conv.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mono_flux_limiter.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fastjx.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/wet_deposition_0D.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_lite_k.F90 $(SRCROOT)atmos_phys/atmos_param/moist_processes/moist_processes_utils.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/math_lib.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/grid_class.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_soa.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/cosp_radar.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/endian.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_cu/conv_closures.F90 $(SRCROOT)atmos_phys/atmos_param/macrophysics/tiedtke_macro.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_parameters.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_setinv.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_driver.F90 $(SRCROOT)atmos_phys/atmos_param/cu_mo_trans/cu_mo_trans.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_sfc.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/isoropiaIIcode.F $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/m_tracname.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/mrgrnk.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_jpl.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_package.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_helper_module.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tropopause.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.subs.F90 $(SRCROOT)atmos_phys/atmos_param/diag_cloud_rad/diag_cloud_rad.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_tracer_driver.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/actsim/lmd_ipsl_stats.F90 $(SRCROOT)atmos_phys/atmos_param/microphysics/micro_mg.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/saturation.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.mods.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_isccp_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/hydrostatic_module.F90 $(SRCROOT)atmos_phys/atmos_param/stable_bl_turb/stable_bl_turb.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver_diag.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_utilities_k.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chemdr.F90 $(SRCROOT)atmos_phys/atmos_param/edt/edt.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/microphys_rad.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator_types.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_3D_var.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem_driver.F90 $(SRCROOT)atmos_phys/atmos_param/microphysics/simple_pdf.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/radiative_gases.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/cloudrad_types.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/fms_donner.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_bands.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/dsd.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sea_salt.F90 $(SRCROOT)atmos_phys/atmos_param/microphysics/cldwat2m_micro.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_physics/shallow_physics.F90 $(SRCROOT)atmos_phys/atmos_param/convection_driver/detr_ice_num.F90 $(SRCROOT)atmos_phys/atmos_param/qe_moist_convection/qe_moist_convection.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_ch3i.F90 $(SRCROOT)atmos_phys/atmos_param/cloud_generator/betaDistribution.F90 $(SRCROOT)atmos_phys/atmos_param/tke_turb/tke_turb.F90 $(SRCROOT)atmos_phys/atmos_param/vert_diff_driver/vert_diff_driver.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/icarus.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/optics_lib.F90 $(SRCROOT)atmos_phys/atmos_param/betts_miller/bm_omp.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/output_grads.F90 $(SRCROOT)atmos_phys/atmos_param/strat_cloud/strat_cloud.F90 $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_constants.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_types.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/optical_path.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/get_random_number_stream.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/gmres_cache.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/atmos_lib.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/isrpia.inc $(SRCROOT)atmos_phys/atmos_param/dry_adj/dry_adj.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/radar_simulator_init.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/input_reader.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/nonfms_donner.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/pdf_parameter_module.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/isccp_clouds.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_windm_edsclrm_module.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_fluxes.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/file_functions.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_wp2_wp3_module.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_age_tracer.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/gas_tf.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/error_code.F90 $(SRCROOT)atmos_phys/atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_utils.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_meso_k.F90 $(SRCROOT)atmos_phys/atmos_param/lscloud_driver/lscloud_types.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clip_explicit.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3_fastjx.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_variables.F90 $(SRCROOT)atmos_phys/atmos_param/convection_driver/moistproc_kernels.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_cape_k.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_driver.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_cloud_model_k.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/numerical_check.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stat_file_module.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/input_names.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_sulfate.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/clouds/strat_clouds_W.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_diagnostics.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM3_chem/AM3.mat.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/gases.F90 $(SRCROOT)atmos_phys/atmos_param/grey_radiation/grey_radiation.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_tables.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90 $(SRCROOT)atmos_phys/atmos_param/moist_conv/moist_conv.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_dust.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/radiative_gases/ozone.F90 $(SRCROOT)atmos_phys/atmos_param/microphysics/ls_cloud_microphysics.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_lidar.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/corr_matrix_module.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/advance_xp2_xpyp_module.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/driver/longwave_driver.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_misr_simulator.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/mo_setsox.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_deep.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F $(SRCROOT)atmos_phys/atmos_param/cosp/llnl/llnl_stats.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/mixing_length.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_rttov_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/MODIS_simulator/modis_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/gmres_wrap.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/driver/solar_data_driver.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/diffusion.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/quickbeam/zeff.F90 $(SRCROOT)atmos_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)atmos_phys/atmos_param/topo_drag/topo_drag.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_carbon_aerosol.F90 $(SRCROOT)atmos_phys/atmos_param/physics_radiation_exch/physics_radiation_exch.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/aerosols/aerosolrad_driver.F90 $(SRCROOT)atmos_phys/atmos_param/shallow_cu/uw_conv.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_rad_zt.F90 $(SRCROOT)atmos_phys/atmos_param/physics_driver/physics_driver.F90 $(SRCROOT)atmos_phys/atmos_param/damping_driver/damping_driver.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_rad_k.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/icarus-scops-3.7/scops.F $(SRCROOT)atmos_phys/atmos_param/radiation/driver/radiation_driver_types.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.F90 $(SRCROOT)atmos_phys/atmos_param/vert_diff/vert_diff.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_zm.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/shortwave/esfsw_utilities.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/driver/shortwave_driver.F90 $(SRCROOT)atmos_phys/atmos_param/my25_turb/my25_turb.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/clubb_core.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/variables_prognostic_module.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/atmos_ch4.F90 $(SRCROOT)atmos_phys/atmos_shared/aerosol/aerosol.F90 $(SRCROOT)atmos_phys/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/Skw_module.F90 $(SRCROOT)atmos_phys/atmos_param/radiation/longwave/longwave_params.F90 $(SRCROOT)atmos_phys/atmos_param/diffusivity/diffusivity.F90 $(SRCROOT)atmos_phys/atmos_param/cloud_generator/cloud_generator.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_io.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_modis_simulator.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/matrix_operations.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/stats_rad_zm.F90 $(SRCROOT)atmos_phys/atmos_shared/atmos_nudge/atmos_nudge.F90 $(SRCROOT)atmos_phys/atmos_param/convection_driver/convection_driver.F90 $(SRCROOT)atmos_phys/atmos_param/cosp/cosp_types.F90 $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/output_netcdf.F90 $(SRCROOT)atmos_phys/atmos_param/donner_deep/donner_types.h $(SRCROOT)atmos_phys/atmos_param/clubb/CLUBB_core/parameters_model.F90
clean: neat
	-rm -f .libatmos_phys.a.cppdefs $(OBJ) libatmos_phys.a

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libatmos_phys.a: $(OBJ)
	$(AR) $(ARFLAGS) libatmos_phys.a $(OBJ)

