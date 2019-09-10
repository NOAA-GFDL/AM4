# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/clubb/CLUBB_core.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libparameters_radiation.la libparameter_indices.la libextrapolation.la \
libmono_flux_limiter.la libgrid_class.la libclubb_precision.la \
libadvance_helper_module.la libadvance_xm_wpxp_module.la \
libdiffusion.la libendian.la libadvance_wp2_wp3_module.la \
libmixing_length.la libconstants_clubb.la libstats_type.la \
libvariables_prognostic_module.la libhyper_diffusion_4th_ord.la \
libinput_names.la libmatrix_operations.la libanl_erf.la libcalendar.la \
libpos_definite_module.la libstats_LH_zt.la libparameters_model.la \
libsaturation.la libSkw_module.la libsigma_sqd_w_module.la \
libhydrostatic_module.la libstats_zt.la libstats_zm.la \
libstats_subs.la libclip_semi_implicit.la libparameters_microphys.la \
liboutput_grads.la libmodel_flags.la libinterpolation.la \
libadvance_windm_edsclrm_module.la libfill_holes.la libgmres_wrap.la \
libT_in_K_module.la libmean_adv.la libstats_variables.la \
libcorr_matrix_module.la libparameters_tunable.la liboutput_netcdf.la \
libpdf_closure_module.la libarray_index.la libclip_explicit.la \
libsurface_varnce_module.la libgmres_cache.la libstats_rad_zm.la \
libstats_sfc.la libstat_file_module.la liblapack_wrap.la \
libinput_reader.la liberror_code.la libadvance_xp2_xpyp_module.la \
libcsr_matrix_class_3array.la libvariables_diagnostic_module.la \
libnumerical_check.la libstats_rad_zt.la libpdf_parameter_module.la \
libsponge_layer_damping.la libfile_functions.la \
libvariables_radiation_module.la libclubb_core.la

# Each convenience library depends on its source.
libparameters_radiation_la_SOURCES = parameters_radiation.F90
libparameter_indices_la_SOURCES = parameter_indices.F90
libextrapolation_la_SOURCES = extrapolation.F90
libmono_flux_limiter_la_SOURCES = mono_flux_limiter.F90
libgrid_class_la_SOURCES = grid_class.F90
libclubb_precision_la_SOURCES = clubb_precision.F90
libadvance_helper_module_la_SOURCES = advance_helper_module.F90
libadvance_xm_wpxp_module_la_SOURCES = advance_xm_wpxp_module.F90
libdiffusion_la_SOURCES = diffusion.F90
libendian_la_SOURCES = endian.F90
libadvance_wp2_wp3_module_la_SOURCES = advance_wp2_wp3_module.F90
libmixing_length_la_SOURCES = mixing_length.F90
libconstants_clubb_la_SOURCES = constants_clubb.F90
libstats_type_la_SOURCES = stats_type.F90
libvariables_prognostic_module_la_SOURCES = variables_prognostic_module.F90
libhyper_diffusion_4th_ord_la_SOURCES = hyper_diffusion_4th_ord.F90
libinput_names_la_SOURCES = input_names.F90
libmatrix_operations_la_SOURCES = matrix_operations.F90
libanl_erf_la_SOURCES = anl_erf.F90
libcalendar_la_SOURCES = calendar.F90
libpos_definite_module_la_SOURCES = pos_definite_module.F90
libstats_LH_zt_la_SOURCES = stats_LH_zt.F90
libparameters_model_la_SOURCES = parameters_model.F90
libsaturation_la_SOURCES = saturation.F90
libSkw_module_la_SOURCES = Skw_module.F90
libsigma_sqd_w_module_la_SOURCES = sigma_sqd_w_module.F90
libhydrostatic_module_la_SOURCES = hydrostatic_module.F90
libstats_zt_la_SOURCES = stats_zt.F90
libstats_zm_la_SOURCES = stats_zm.F90
libstats_subs_la_SOURCES = stats_subs.F90
libclip_semi_implicit_la_SOURCES = clip_semi_implicit.F90
libparameters_microphys_la_SOURCES = parameters_microphys.F90
liboutput_grads_la_SOURCES = output_grads.F90
libmodel_flags_la_SOURCES = model_flags.F90
libinterpolation_la_SOURCES = interpolation.F90
libadvance_windm_edsclrm_module_la_SOURCES = advance_windm_edsclrm_module.F90
libfill_holes_la_SOURCES = fill_holes.F90
libgmres_wrap_la_SOURCES = gmres_wrap.F90
libT_in_K_module_la_SOURCES = T_in_K_module.F90
libmean_adv_la_SOURCES = mean_adv.F90
libstats_variables_la_SOURCES = stats_variables.F90
libcorr_matrix_module_la_SOURCES = corr_matrix_module.F90
libparameters_tunable_la_SOURCES = parameters_tunable.F90
liboutput_netcdf_la_SOURCES = output_netcdf.F90
libpdf_closure_module_la_SOURCES = pdf_closure_module.F90
libarray_index_la_SOURCES = array_index.F90
libclip_explicit_la_SOURCES = clip_explicit.F90
libsurface_varnce_module_la_SOURCES = surface_varnce_module.F90
libgmres_cache_la_SOURCES = gmres_cache.F90
libstats_rad_zm_la_SOURCES = stats_rad_zm.F90
libstats_sfc_la_SOURCES = stats_sfc.F90
libstat_file_module_la_SOURCES = stat_file_module.F90
liblapack_wrap_la_SOURCES = lapack_wrap.F90
libinput_reader_la_SOURCES = input_reader.F90
liberror_code_la_SOURCES = error_code.F90
libadvance_xp2_xpyp_module_la_SOURCES = advance_xp2_xpyp_module.F90
libcsr_matrix_class_3array_la_SOURCES = csr_matrix_class_3array.F90
libvariables_diagnostic_module_la_SOURCES = variables_diagnostic_module.F90
libnumerical_check_la_SOURCES = numerical_check.F90
libstats_rad_zt_la_SOURCES = stats_rad_zt.F90
libpdf_parameter_module_la_SOURCES = pdf_parameter_module.F90
libsponge_layer_damping_la_SOURCES = sponge_layer_damping.F90
libfile_functions_la_SOURCES = file_functions.F90
libvariables_radiation_module_la_SOURCES = variables_radiation_module.F90
libclubb_core_la_SOURCES = clubb_core.F90

# Each mod file depends on the .o file.
parameters_radiation_mod.mod: parameters_radiation.$(OBJEXT)
parameter_indices_mod.mod: parameter_indices.$(OBJEXT)
extrapolation_mod.mod: extrapolation.$(OBJEXT)
mono_flux_limiter_mod.mod: mono_flux_limiter.$(OBJEXT)
grid_class_mod.mod: grid_class.$(OBJEXT)
clubb_precision_mod.mod: clubb_precision.$(OBJEXT)
advance_helper_module_mod.mod: advance_helper_module.$(OBJEXT)
advance_xm_wpxp_module_mod.mod: advance_xm_wpxp_module.$(OBJEXT)
diffusion_mod.mod: diffusion.$(OBJEXT)
endian_mod.mod: endian.$(OBJEXT)
advance_wp2_wp3_module_mod.mod: advance_wp2_wp3_module.$(OBJEXT)
mixing_length_mod.mod: mixing_length.$(OBJEXT)
constants_clubb_mod.mod: constants_clubb.$(OBJEXT)
stats_type_mod.mod: stats_type.$(OBJEXT)
variables_prognostic_module_mod.mod: variables_prognostic_module.$(OBJEXT)
hyper_diffusion_4th_ord_mod.mod: hyper_diffusion_4th_ord.$(OBJEXT)
input_names_mod.mod: input_names.$(OBJEXT)
matrix_operations_mod.mod: matrix_operations.$(OBJEXT)
anl_erf_mod.mod: anl_erf.$(OBJEXT)
calendar_mod.mod: calendar.$(OBJEXT)
pos_definite_module_mod.mod: pos_definite_module.$(OBJEXT)
stats_LH_zt_mod.mod: stats_LH_zt.$(OBJEXT)
parameters_model_mod.mod: parameters_model.$(OBJEXT)
saturation_mod.mod: saturation.$(OBJEXT)
Skw_module_mod.mod: Skw_module.$(OBJEXT)
sigma_sqd_w_module_mod.mod: sigma_sqd_w_module.$(OBJEXT)
hydrostatic_module_mod.mod: hydrostatic_module.$(OBJEXT)
stats_zt_mod.mod: stats_zt.$(OBJEXT)
stats_zm_mod.mod: stats_zm.$(OBJEXT)
stats_subs_mod.mod: stats_subs.$(OBJEXT)
clip_semi_implicit_mod.mod: clip_semi_implicit.$(OBJEXT)
parameters_microphys_mod.mod: parameters_microphys.$(OBJEXT)
output_grads_mod.mod: output_grads.$(OBJEXT)
model_flags_mod.mod: model_flags.$(OBJEXT)
interpolation_mod.mod: interpolation.$(OBJEXT)
advance_windm_edsclrm_module_mod.mod: advance_windm_edsclrm_module.$(OBJEXT)
fill_holes_mod.mod: fill_holes.$(OBJEXT)
gmres_wrap_mod.mod: gmres_wrap.$(OBJEXT)
T_in_K_module_mod.mod: T_in_K_module.$(OBJEXT)
mean_adv_mod.mod: mean_adv.$(OBJEXT)
stats_variables_mod.mod: stats_variables.$(OBJEXT)
corr_matrix_module_mod.mod: corr_matrix_module.$(OBJEXT)
parameters_tunable_mod.mod: parameters_tunable.$(OBJEXT)
output_netcdf_mod.mod: output_netcdf.$(OBJEXT)
pdf_closure_module_mod.mod: pdf_closure_module.$(OBJEXT)
array_index_mod.mod: array_index.$(OBJEXT)
clip_explicit_mod.mod: clip_explicit.$(OBJEXT)
surface_varnce_module_mod.mod: surface_varnce_module.$(OBJEXT)
gmres_cache_mod.mod: gmres_cache.$(OBJEXT)
stats_rad_zm_mod.mod: stats_rad_zm.$(OBJEXT)
stats_sfc_mod.mod: stats_sfc.$(OBJEXT)
stat_file_module_mod.mod: stat_file_module.$(OBJEXT)
lapack_wrap_mod.mod: lapack_wrap.$(OBJEXT)
input_reader_mod.mod: input_reader.$(OBJEXT)
error_code_mod.mod: error_code.$(OBJEXT)
advance_xp2_xpyp_module_mod.mod: advance_xp2_xpyp_module.$(OBJEXT)
csr_matrix_class_3array_mod.mod: csr_matrix_class_3array.$(OBJEXT)
variables_diagnostic_module_mod.mod: variables_diagnostic_module.$(OBJEXT)
numerical_check_mod.mod: numerical_check.$(OBJEXT)
stats_rad_zt_mod.mod: stats_rad_zt.$(OBJEXT)
pdf_parameter_module_mod.mod: pdf_parameter_module.$(OBJEXT)
sponge_layer_damping_mod.mod: sponge_layer_damping.$(OBJEXT)
file_functions_mod.mod: file_functions.$(OBJEXT)
variables_radiation_module_mod.mod: variables_radiation_module.$(OBJEXT)
clubb_core_mod.mod: clubb_core.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
parameters_radiation.$(OBJEXT): clubb_precision_mod.mod
extrapolation.$(OBJEXT): clubb_precision_mod.mod
mono_flux_limiter.$(OBJEXT): grid_class_mod.mod clubb_precision_mod.mod constants_clubb_mod.mod stats_type_mod.mod anl_erf_mod.mod fill_holes_mod.mod mean_adv_mod.mod stats_variables_mod.mod lapack_wrap_mod.mod error_code_mod.mod
grid_class.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod interpolation_mod.mod error_code_mod.mod file_functions_mod.mod
advance_helper_module.$(OBJEXT): clubb_precision_mod.mod
advance_xm_wpxp_module.$(OBJEXT): mono_flux_limiter_mod.mod grid_class_mod.mod clubb_precision_mod.mod advance_helper_module_mod.mod diffusion_mod.mod constants_clubb_mod.mod stats_type_mod.mod pos_definite_module_mod.mod parameters_model_mod.mod clip_semi_implicit_mod.mod model_flags_mod.mod fill_holes_mod.mod mean_adv_mod.mod stats_variables_mod.mod parameters_tunable_mod.mod clip_explicit_mod.mod lapack_wrap_mod.mod error_code_mod.mod pdf_parameter_module_mod.mod sponge_layer_damping_mod.mod
diffusion.$(OBJEXT): grid_class_mod.mod clubb_precision_mod.mod
advance_wp2_wp3_module.$(OBJEXT): grid_class_mod.mod clubb_precision_mod.mod advance_helper_module_mod.mod diffusion_mod.mod constants_clubb_mod.mod stats_type_mod.mod hyper_diffusion_4th_ord_mod.mod model_flags_mod.mod fill_holes_mod.mod gmres_wrap_mod.mod mean_adv_mod.mod stats_variables_mod.mod parameters_tunable_mod.mod clip_explicit_mod.mod gmres_cache_mod.mod lapack_wrap_mod.mod error_code_mod.mod
mixing_length.$(OBJEXT): grid_class_mod.mod clubb_precision_mod.mod constants_clubb_mod.mod parameters_model_mod.mod saturation_mod.mod model_flags_mod.mod parameters_tunable_mod.mod error_code_mod.mod numerical_check_mod.mod
constants_clubb.$(OBJEXT): clubb_precision_mod.mod
stats_type.$(OBJEXT): grid_class_mod.mod clubb_precision_mod.mod stat_file_module_mod.mod error_code_mod.mod
variables_prognostic_module.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod parameters_model_mod.mod pdf_parameter_module_mod.mod
hyper_diffusion_4th_ord.$(OBJEXT): grid_class_mod.mod clubb_precision_mod.mod
matrix_operations.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod error_code_mod.mod
calendar.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod
pos_definite_module.$(OBJEXT): grid_class_mod.mod clubb_precision_mod.mod constants_clubb_mod.mod error_code_mod.mod
stats_LH_zt.$(OBJEXT): constants_clubb_mod.mod stats_type_mod.mod stats_variables_mod.mod
parameters_model.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod
saturation.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod model_flags_mod.mod
Skw_module.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod
sigma_sqd_w_module.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod
hydrostatic_module.$(OBJEXT): grid_class_mod.mod clubb_precision_mod.mod constants_clubb_mod.mod interpolation_mod.mod
stats_zt.$(OBJEXT): constants_clubb_mod.mod stats_type_mod.mod parameters_model_mod.mod stats_variables_mod.mod error_code_mod.mod
stats_zm.$(OBJEXT): constants_clubb_mod.mod stats_type_mod.mod parameters_model_mod.mod stats_variables_mod.mod error_code_mod.mod
stats_subs.$(OBJEXT): grid_class_mod.mod clubb_precision_mod.mod constants_clubb_mod.mod stats_type_mod.mod stats_LH_zt_mod.mod parameters_model_mod.mod saturation_mod.mod stats_zt_mod.mod stats_zm_mod.mod parameters_microphys_mod.mod output_grads_mod.mod interpolation_mod.mod fill_holes_mod.mod T_in_K_module_mod.mod stats_variables_mod.mod output_netcdf_mod.mod array_index_mod.mod stats_rad_zm_mod.mod stats_sfc_mod.mod error_code_mod.mod variables_diagnostic_module_mod.mod stats_rad_zt_mod.mod pdf_parameter_module_mod.mod
clip_semi_implicit.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod
parameters_microphys.$(OBJEXT): clubb_precision_mod.mod
output_grads.$(OBJEXT): clubb_precision_mod.mod endian_mod.mod constants_clubb_mod.mod calendar_mod.mod model_flags_mod.mod stat_file_module_mod.mod
model_flags.$(OBJEXT): constants_clubb_mod.mod
interpolation.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod model_flags_mod.mod error_code_mod.mod
advance_windm_edsclrm_module.$(OBJEXT): grid_class_mod.mod clubb_precision_mod.mod diffusion_mod.mod constants_clubb_mod.mod stats_type_mod.mod parameters_model_mod.mod model_flags_mod.mod mean_adv_mod.mod stats_variables_mod.mod parameters_tunable_mod.mod clip_explicit_mod.mod lapack_wrap_mod.mod error_code_mod.mod sponge_layer_damping_mod.mod
fill_holes.$(OBJEXT): grid_class_mod.mod clubb_precision_mod.mod
gmres_wrap.$(OBJEXT): clubb_precision_mod.mod gmres_cache_mod.mod
T_in_K_module.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod
mean_adv.$(OBJEXT): grid_class_mod.mod clubb_precision_mod.mod model_flags_mod.mod
stats_variables.$(OBJEXT): clubb_precision_mod.mod stats_type_mod.mod
corr_matrix_module.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod matrix_operations_mod.mod input_reader_mod.mod
parameters_tunable.$(OBJEXT): parameter_indices_mod.mod grid_class_mod.mod clubb_precision_mod.mod constants_clubb_mod.mod error_code_mod.mod
output_netcdf.$(OBJEXT): parameters_radiation_mod.mod parameter_indices_mod.mod clubb_precision_mod.mod constants_clubb_mod.mod calendar_mod.mod parameters_model_mod.mod parameters_microphys_mod.mod model_flags_mod.mod parameters_tunable_mod.mod stat_file_module_mod.mod
pdf_closure_module.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod anl_erf_mod.mod parameters_model_mod.mod saturation_mod.mod stats_variables_mod.mod parameters_tunable_mod.mod error_code_mod.mod numerical_check_mod.mod pdf_parameter_module_mod.mod
clip_explicit.$(OBJEXT): grid_class_mod.mod clubb_precision_mod.mod constants_clubb_mod.mod stats_type_mod.mod parameters_model_mod.mod model_flags_mod.mod stats_variables_mod.mod
surface_varnce_module.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod stats_type_mod.mod parameters_model_mod.mod array_index_mod.mod error_code_mod.mod numerical_check_mod.mod
gmres_cache.$(OBJEXT): clubb_precision_mod.mod
stats_rad_zm.$(OBJEXT): constants_clubb_mod.mod stats_type_mod.mod stats_variables_mod.mod
stats_sfc.$(OBJEXT): constants_clubb_mod.mod stats_type_mod.mod stats_variables_mod.mod
stat_file_module.$(OBJEXT): clubb_precision_mod.mod
lapack_wrap.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod error_code_mod.mod
input_reader.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod input_names_mod.mod interpolation_mod.mod
error_code.$(OBJEXT): constants_clubb_mod.mod
advance_xp2_xpyp_module.$(OBJEXT): grid_class_mod.mod clubb_precision_mod.mod advance_helper_module_mod.mod diffusion_mod.mod constants_clubb_mod.mod stats_type_mod.mod parameters_model_mod.mod model_flags_mod.mod interpolation_mod.mod fill_holes_mod.mod mean_adv_mod.mod stats_variables_mod.mod parameters_tunable_mod.mod array_index_mod.mod clip_explicit_mod.mod lapack_wrap_mod.mod error_code_mod.mod
csr_matrix_class_3array.$(OBJEXT): grid_class_mod.mod constants_clubb_mod.mod
variables_diagnostic_module.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod parameters_model_mod.mod pdf_parameter_module_mod.mod
numerical_check.$(OBJEXT): grid_class_mod.mod clubb_precision_mod.mod constants_clubb_mod.mod variables_prognostic_module_mod.mod parameters_model_mod.mod parameters_microphys_mod.mod stats_variables_mod.mod error_code_mod.mod variables_diagnostic_module_mod.mod pdf_parameter_module_mod.mod
stats_rad_zt.$(OBJEXT): constants_clubb_mod.mod stats_type_mod.mod stats_variables_mod.mod
pdf_parameter_module.$(OBJEXT): clubb_precision_mod.mod
sponge_layer_damping.$(OBJEXT): grid_class_mod.mod clubb_precision_mod.mod constants_clubb_mod.mod interpolation_mod.mod
file_functions.$(OBJEXT): clubb_precision_mod.mod constants_clubb_mod.mod
variables_radiation_module.$(OBJEXT): clubb_precision_mod.mod
clubb_core.$(OBJEXT): parameter_indices_mod.mod grid_class_mod.mod clubb_precision_mod.mod advance_xm_wpxp_module_mod.mod advance_wp2_wp3_module_mod.mod mixing_length_mod.mod constants_clubb_mod.mod stats_type_mod.mod variables_prognostic_module_mod.mod parameters_model_mod.mod saturation_mod.mod Skw_module_mod.mod sigma_sqd_w_module_mod.mod stats_subs_mod.mod model_flags_mod.mod advance_windm_edsclrm_module_mod.mod fill_holes_mod.mod gmres_wrap_mod.mod T_in_K_module_mod.mod stats_variables_mod.mod parameters_tunable_mod.mod pdf_closure_module_mod.mod clip_explicit_mod.mod surface_varnce_module_mod.mod gmres_cache_mod.mod error_code_mod.mod advance_xp2_xpyp_module_mod.mod variables_diagnostic_module_mod.mod numerical_check_mod.mod pdf_parameter_module_mod.mod

# Mod files are built and then installed as headers.
MODFILES = parameters_radiation_mod.mod parameter_indices_mod.mod \
extrapolation_mod.mod mono_flux_limiter_mod.mod grid_class_mod.mod \
clubb_precision_mod.mod advance_helper_module_mod.mod \
advance_xm_wpxp_module_mod.mod diffusion_mod.mod endian_mod.mod \
advance_wp2_wp3_module_mod.mod mixing_length_mod.mod \
constants_clubb_mod.mod stats_type_mod.mod \
variables_prognostic_module_mod.mod hyper_diffusion_4th_ord_mod.mod \
input_names_mod.mod matrix_operations_mod.mod anl_erf_mod.mod \
calendar_mod.mod pos_definite_module_mod.mod stats_LH_zt_mod.mod \
parameters_model_mod.mod saturation_mod.mod Skw_module_mod.mod \
sigma_sqd_w_module_mod.mod hydrostatic_module_mod.mod stats_zt_mod.mod \
stats_zm_mod.mod stats_subs_mod.mod clip_semi_implicit_mod.mod \
parameters_microphys_mod.mod output_grads_mod.mod model_flags_mod.mod \
interpolation_mod.mod advance_windm_edsclrm_module_mod.mod \
fill_holes_mod.mod gmres_wrap_mod.mod T_in_K_module_mod.mod \
mean_adv_mod.mod stats_variables_mod.mod corr_matrix_module_mod.mod \
parameters_tunable_mod.mod output_netcdf_mod.mod \
pdf_closure_module_mod.mod array_index_mod.mod clip_explicit_mod.mod \
surface_varnce_module_mod.mod gmres_cache_mod.mod stats_rad_zm_mod.mod \
stats_sfc_mod.mod stat_file_module_mod.mod lapack_wrap_mod.mod \
input_reader_mod.mod error_code_mod.mod \
advance_xp2_xpyp_module_mod.mod csr_matrix_class_3array_mod.mod \
variables_diagnostic_module_mod.mod numerical_check_mod.mod \
stats_rad_zt_mod.mod pdf_parameter_module_mod.mod \
sponge_layer_damping_mod.mod file_functions_mod.mod \
variables_radiation_module_mod.mod clubb_core_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = recl.inc

CLEANFILES = *.mod
