# Makefile created by mkmf 19.3.1

CPPDEFS = -DINTERNAL_FILE_NML -Duse_libMPI -Duse_netCDF

TMPFILES = *genmod.f90

include $(MK_TEMPLATE)

.DEFAULT:
	-echo $@ does not exist.
all: libfms.a
affinity.o: $(SRCROOT)FMS/affinity/affinity.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)FMS/affinity/affinity.c
amip_interp.o: $(SRCROOT)FMS/amip_interp/amip_interp.F90 $(SRCROOT)FMS/include/file_version.h time_interp.o time_manager.o get_cal_time.o horiz_interp.o fms.o constants.o platform.o mpp.o fms2_io.o fms_io.o mpp_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/amip_interp/amip_interp.F90
astronomy.o: $(SRCROOT)FMS/astronomy/astronomy.F90 $(SRCROOT)FMS/include/file_version.h fms.o time_manager.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/astronomy/astronomy.F90
atmos_ocean_fluxes.o: $(SRCROOT)FMS/coupler/atmos_ocean_fluxes.F90 $(SRCROOT)FMS/include/file_version.h mpp.o fms.o coupler_types.o field_manager.o fm_util.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/coupler/atmos_ocean_fluxes.F90
axis_utils.o: $(SRCROOT)FMS/axis_utils/axis_utils.F90 $(SRCROOT)FMS/include/file_version.h mpp_io.o mpp.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/axis_utils/axis_utils.F90
axis_utils2.o: $(SRCROOT)FMS/axis_utils/axis_utils2.F90 $(SRCROOT)FMS/include/file_version.h mpp.o fms.o fms2_io.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/axis_utils/axis_utils2.F90
blackboxio.o: $(SRCROOT)FMS/fms2_io/blackboxio.F90 mpp_domains.o fms_io_utils.o netcdf_io.o fms_netcdf_domain_io.o fms_netcdf_unstructured_domain_io.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/fms2_io/blackboxio.F90
block_control.o: $(SRCROOT)FMS/block_control/block_control.F90 mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/block_control/block_control.F90
cloud_interpolator.o: $(SRCROOT)FMS/drifters/cloud_interpolator.F90 $(SRCROOT)FMS/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/drifters/cloud_interpolator.F90
column_diagnostics.o: $(SRCROOT)FMS/column_diagnostics/column_diagnostics.F90 $(SRCROOT)FMS/include/file_version.h mpp_io.o fms.o time_manager.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/column_diagnostics/column_diagnostics.F90
constants.o: $(SRCROOT)FMS/constants/constants.F90 $(SRCROOT)FMS/include/file_version.h platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/constants/constants.F90
coupler_types.o: $(SRCROOT)FMS/coupler/coupler_types.F90 $(SRCROOT)FMS/include/file_version.h fms.o fms_io.o time_manager.o diag_manager.o data_override.o mpp_domains.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/coupler/coupler_types.F90
create_xgrid.o: $(SRCROOT)FMS/mosaic/create_xgrid.c $(SRCROOT)FMS/mosaic/mosaic_util.h $(SRCROOT)FMS/mosaic/create_xgrid.h $(SRCROOT)FMS/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)FMS/mosaic	$(SRCROOT)FMS/mosaic/create_xgrid.c
data_override.o: $(SRCROOT)FMS/data_override/data_override.F90 $(SRCROOT)FMS/include/file_version.h constants.o mpp_io.o mpp.o horiz_interp.o time_interp_external.o time_interp_external2.o fms.o fms_io.o axis_utils.o axis_utils2.o mpp_domains.o time_manager.o fms2_io.o get_grid_version_mpp.o get_grid_version_fms2io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/data_override/data_override.F90
diag_axis.o: $(SRCROOT)FMS/diag_manager/diag_axis.F90 $(SRCROOT)FMS/include/file_version.h platform.o mpp_domains.o fms.o diag_data.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/diag_manager/diag_axis.F90
diag_data.o: $(SRCROOT)FMS/diag_manager/diag_data.F90 $(SRCROOT)FMS/include/file_version.h platform.o time_manager.o mpp_domains.o mpp_io.o fms.o fms2_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/diag_manager/diag_data.F90
diag_grid.o: $(SRCROOT)FMS/diag_manager/diag_grid.F90 $(SRCROOT)FMS/include/file_version.h platform.o constants.o fms.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/diag_manager/diag_grid.F90
diag_integral.o: $(SRCROOT)FMS/diag_integral/diag_integral.F90 $(SRCROOT)FMS/include/file_version.h platform.o time_manager.o mpp.o fms.o fms2_io.o constants.o ensemble_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/diag_integral/diag_integral.F90
diag_manager.o: $(SRCROOT)FMS/diag_manager/diag_manager.F90 $(SRCROOT)FMS/include/file_version.h platform.o time_manager.o mpp_io.o mpp.o fms.o fms_io.o diag_axis.o diag_util.o diag_data.o diag_table.o diag_output.o diag_grid.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/diag_manager/diag_manager.F90
diag_output.o: $(SRCROOT)FMS/diag_manager/diag_output.F90 platform.o mpp_io.o mpp_domains.o mpp.o diag_axis.o diag_data.o time_manager.o fms.o fms2_io.o axis_utils2.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/diag_manager/diag_output.F90
diag_table.o: $(SRCROOT)FMS/diag_manager/diag_table.F90 mpp_io.o mpp.o fms.o time_manager.o constants.o diag_data.o diag_util.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/diag_manager/diag_table.F90
diag_util.o: $(SRCROOT)FMS/diag_manager/diag_util.F90 $(SRCROOT)FMS/include/file_version.h platform.o diag_data.o diag_axis.o diag_output.o diag_grid.o fms.o fms_io.o mpp_domains.o time_manager.o mpp_io.o mpp.o constants.o fms2_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/diag_manager/diag_util.F90
drifters.o: $(SRCROOT)FMS/drifters/drifters.F90 $(SRCROOT)FMS/drifters/fms_switches.h $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/drifters/drifters_push.h $(SRCROOT)FMS/drifters/drifters_set_field.h $(SRCROOT)FMS/drifters/drifters_compute_k.h mpp.o mpp_domains.o drifters_core.o drifters_input.o drifters_io.o drifters_comm.o cloud_interpolator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/drifters -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/drifters/drifters.F90
drifters_comm.o: $(SRCROOT)FMS/drifters/drifters_comm.F90 $(SRCROOT)FMS/drifters/fms_switches.h mpp.o mpp_domains.o drifters_core.o drifters_input.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/drifters	$(SRCROOT)FMS/drifters/drifters_comm.F90
drifters_core.o: $(SRCROOT)FMS/drifters/drifters_core.F90 $(SRCROOT)FMS/include/file_version.h platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/drifters/drifters_core.F90
drifters_input.o: $(SRCROOT)FMS/drifters/drifters_input.F90 $(SRCROOT)FMS/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/drifters/drifters_input.F90
drifters_io.o: $(SRCROOT)FMS/drifters/drifters_io.F90 $(SRCROOT)FMS/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/drifters/drifters_io.F90
ensemble_manager.o: $(SRCROOT)FMS/coupler/ensemble_manager.F90 fms.o mpp.o fms_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/coupler/ensemble_manager.F90
field_manager.o: $(SRCROOT)FMS/field_manager/field_manager.F90 $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/field_manager/parse.inc mpp.o mpp_io.o fms.o fms2_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include -I$(SRCROOT)FMS/field_manager	$(SRCROOT)FMS/field_manager/field_manager.F90
fm_util.o: $(SRCROOT)FMS/field_manager/fm_util.F90 $(SRCROOT)FMS/include/file_version.h field_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/field_manager/fm_util.F90
fms.o: $(SRCROOT)FMS/fms/fms.F90 $(SRCROOT)FMS/include/file_version.h mpp.o mpp_domains.o mpp_io.o fms_io.o fms2_io.o memutils.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/fms/fms.F90
fms2_io.o: $(SRCROOT)FMS/fms2_io/fms2_io.F90 fms_io_utils.o netcdf_io.o fms_netcdf_domain_io.o fms_netcdf_unstructured_domain_io.o blackboxio.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/fms2_io/fms2_io.F90
fms_affinity.o: $(SRCROOT)FMS/affinity/fms_affinity.F90 $(SRCROOT)FMS/include/file_version.h mpp.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/affinity/fms_affinity.F90
fms_io.o: $(SRCROOT)FMS/fms/fms_io.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/fms/read_data_2d.inc $(SRCROOT)FMS/fms/read_data_3d.inc $(SRCROOT)FMS/fms/read_data_4d.inc $(SRCROOT)FMS/fms/write_data.inc $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_axis.inc $(SRCROOT)FMS/fms/fms_io_unstructured_setup_one_field.inc $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_field.inc $(SRCROOT)FMS/fms/fms_io_unstructured_save_restart.inc $(SRCROOT)FMS/fms/fms_io_unstructured_read.inc $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_name.inc $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_unit.inc $(SRCROOT)FMS/fms/fms_io_unstructured_file_unit.inc $(SRCROOT)FMS/fms/fms_io_unstructured_get_field_size.inc $(SRCROOT)FMS/fms/fms_io_unstructured_field_exist.inc mpp_io.o mpp_domains.o mpp.o platform.o mpp_parameter.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include -I$(SRCROOT)FMS/fms	$(SRCROOT)FMS/fms/fms_io.F90
fms_io_utils.o: $(SRCROOT)FMS/fms2_io/fms_io_utils.F90 $(SRCROOT)FMS/fms2_io/include/array_utils.inc $(SRCROOT)FMS/fms2_io/include/array_utils_char.inc $(SRCROOT)FMS/fms2_io/include/get_data_type_string.inc $(SRCROOT)FMS/fms2_io/include/get_checksum.inc mpp.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/fms2_io/include	$(SRCROOT)FMS/fms2_io/fms_io_utils.F90
fms_netcdf_domain_io.o: $(SRCROOT)FMS/fms2_io/fms_netcdf_domain_io.F90 $(SRCROOT)FMS/fms2_io/include/register_domain_restart_variable.inc $(SRCROOT)FMS/fms2_io/include/domain_read.inc $(SRCROOT)FMS/fms2_io/include/domain_write.inc $(SRCROOT)FMS/fms2_io/include/compute_global_checksum.inc mpp.o mpp_domains.o fms_io_utils.o netcdf_io.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/fms2_io/include	$(SRCROOT)FMS/fms2_io/fms_netcdf_domain_io.F90
fms_netcdf_unstructured_domain_io.o: $(SRCROOT)FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90 $(SRCROOT)FMS/fms2_io/include/register_unstructured_domain_restart_variable.inc $(SRCROOT)FMS/fms2_io/include/unstructured_domain_read.inc $(SRCROOT)FMS/fms2_io/include/unstructured_domain_write.inc mpp_domains.o fms_io_utils.o netcdf_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/fms2_io/include	$(SRCROOT)FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90
gaussian_topog.o: $(SRCROOT)FMS/topography/gaussian_topog.F90 $(SRCROOT)FMS/include/file_version.h fms.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/topography/gaussian_topog.F90
get_cal_time.o: $(SRCROOT)FMS/time_manager/get_cal_time.F90 $(SRCROOT)FMS/include/file_version.h fms.o time_manager.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/time_manager/get_cal_time.F90
get_grid_version_fms2io.o: $(SRCROOT)FMS/data_override/get_grid_version_fms2io.F90 constants.o mpp.o mpp_domains.o fms2_io.o mosaic2.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/data_override/get_grid_version_fms2io.F90
get_grid_version_mpp.o: $(SRCROOT)FMS/data_override/get_grid_version_mpp.F90 constants.o mpp.o fms_io.o fms.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/data_override/get_grid_version_mpp.F90
gradient.o: $(SRCROOT)FMS/mosaic/gradient.F90 $(SRCROOT)FMS/include/file_version.h mpp.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/mosaic/gradient.F90
gradient_c2l.o: $(SRCROOT)FMS/mosaic/gradient_c2l.c $(SRCROOT)FMS/mosaic/constant.h $(SRCROOT)FMS/mosaic/mosaic_util.h $(SRCROOT)FMS/mosaic/gradient_c2l.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)FMS/mosaic	$(SRCROOT)FMS/mosaic/gradient_c2l.c
grid.o: $(SRCROOT)FMS/mosaic/grid.F90 $(SRCROOT)FMS/include/file_version.h mpp.o constants.o fms.o fms_io.o mosaic.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/mosaic/grid.F90
horiz_interp.o: $(SRCROOT)FMS/horiz_interp/horiz_interp.F90 $(SRCROOT)FMS/include/file_version.h fms.o mpp.o constants.o horiz_interp_type.o horiz_interp_conserve.o horiz_interp_bilinear.o horiz_interp_bicubic.o horiz_interp_spherical.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/horiz_interp/horiz_interp.F90
horiz_interp_bicubic.o: $(SRCROOT)FMS/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)FMS/include/file_version.h mpp.o fms.o horiz_interp_type.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/horiz_interp/horiz_interp_bicubic.F90
horiz_interp_bilinear.o: $(SRCROOT)FMS/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)FMS/include/file_version.h mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/horiz_interp/horiz_interp_bilinear.F90
horiz_interp_conserve.o: $(SRCROOT)FMS/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h mpp.o fms.o fms_io.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/horiz_interp/horiz_interp_conserve.F90
horiz_interp_spherical.o: $(SRCROOT)FMS/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)FMS/include/file_version.h mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/horiz_interp/horiz_interp_spherical.F90
horiz_interp_type.o: $(SRCROOT)FMS/horiz_interp/horiz_interp_type.F90 mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/horiz_interp/horiz_interp_type.F90
interp.o: $(SRCROOT)FMS/mosaic/interp.c $(SRCROOT)FMS/mosaic/mosaic_util.h $(SRCROOT)FMS/mosaic/interp.h $(SRCROOT)FMS/mosaic/create_xgrid.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)FMS/mosaic	$(SRCROOT)FMS/mosaic/interp.c
interpolator.o: $(SRCROOT)FMS/interpolator/interpolator.F90 $(SRCROOT)FMS/include/file_version.h mpp.o mpp_io.o mpp_domains.o diag_manager.o fms.o fms2_io.o horiz_interp.o time_manager.o time_interp.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/interpolator/interpolator.F90
memutils.o: $(SRCROOT)FMS/memutils/memutils.F90 mpp.o mpp_memutils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/memutils/memutils.F90
mersennetwister.o: $(SRCROOT)FMS/random_numbers/mersennetwister.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/random_numbers/mersennetwister.F90
monin_obukhov.o: $(SRCROOT)FMS/monin_obukhov/monin_obukhov.F90 $(SRCROOT)FMS/include/file_version.h constants.o mpp.o fms.o monin_obukhov_inter.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/monin_obukhov/monin_obukhov.F90
monin_obukhov_inter.o: $(SRCROOT)FMS/monin_obukhov/monin_obukhov_inter.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/monin_obukhov/monin_obukhov_inter.F90
mosaic.o: $(SRCROOT)FMS/mosaic/mosaic.F90 $(SRCROOT)FMS/include/file_version.h fms.o mpp.o mpp_io.o fms_io.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/mosaic/mosaic.F90
mosaic2.o: $(SRCROOT)FMS/mosaic/mosaic2.F90 $(SRCROOT)FMS/include/file_version.h fms.o mpp.o mpp_domains.o constants.o fms2_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/mosaic/mosaic2.F90
mosaic_util.o: $(SRCROOT)FMS/mosaic/mosaic_util.c $(SRCROOT)FMS/mosaic/mosaic_util.h $(SRCROOT)FMS/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)FMS/mosaic	$(SRCROOT)FMS/mosaic/mosaic_util.c
mpp.o: $(SRCROOT)FMS/mpp/mpp.F90 $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/mpp/include/system_clock.h $(SRCROOT)FMS/mpp/include/mpp_util.inc $(SRCROOT)FMS/mpp/include/mpp_util_mpi.inc $(SRCROOT)FMS/mpp/include/mpp_util_nocomm.inc $(SRCROOT)FMS/mpp/include/mpp_error_a_a.h $(SRCROOT)FMS/mpp/include/mpp_error_a_s.h $(SRCROOT)FMS/mpp/include/mpp_error_s_a.h $(SRCROOT)FMS/mpp/include/mpp_error_s_s.h $(SRCROOT)FMS/mpp/include/mpp_comm.inc $(SRCROOT)FMS/mpp/include/mpp_comm_mpi.inc $(SRCROOT)FMS/mpp/include/mpp_transmit_mpi.h $(SRCROOT)FMS/mpp/include/mpp_transmit.inc $(SRCROOT)FMS/mpp/include/mpp_reduce_mpi.h $(SRCROOT)FMS/mpp/include/mpp_sum_mpi.h $(SRCROOT)FMS/mpp/include/mpp_sum.inc $(SRCROOT)FMS/mpp/include/mpp_sum_mpi_ad.h $(SRCROOT)FMS/mpp/include/mpp_sum_ad.inc $(SRCROOT)FMS/mpp/include/mpp_alltoall_mpi.h $(SRCROOT)FMS/mpp/include/mpp_type_mpi.h $(SRCROOT)FMS/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)FMS/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_reduce_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm_ad.h $(SRCROOT)FMS/mpp/include/mpp_alltoall_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_type_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_chksum_int.h $(SRCROOT)FMS/mpp/include/mpp_chksum_scalar.h $(SRCROOT)FMS/mpp/include/mpp_chksum.h $(SRCROOT)FMS/mpp/include/mpp_gather.h $(SRCROOT)FMS/mpp/include/mpp_scatter.h mpp_parameter.o mpp_data.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include -I$(SRCROOT)FMS/mpp/include	$(SRCROOT)FMS/mpp/mpp.F90
mpp_data.o: $(SRCROOT)FMS/mpp/mpp_data.F90 $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/mpp/include/mpp_data_mpi.inc $(SRCROOT)FMS/mpp/include/mpp_data_nocomm.inc mpp_parameter.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include -I$(SRCROOT)FMS/mpp/include	$(SRCROOT)FMS/mpp/mpp_data.F90
mpp_domains.o: $(SRCROOT)FMS/mpp/mpp_domains.F90 $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/mpp/include/mpp_define_nest_domains.inc $(SRCROOT)FMS/mpp/include/mpp_domains_util.inc $(SRCROOT)FMS/mpp/include/mpp_domains_comm.inc $(SRCROOT)FMS/mpp/include/mpp_domains_define.inc $(SRCROOT)FMS/mpp/include/mpp_domains_misc.inc $(SRCROOT)FMS/mpp/include/mpp_update_domains2D.h $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_nonblock.h $(SRCROOT)FMS/mpp/include/mpp_do_update_nonblock.h $(SRCROOT)FMS/mpp/include/mpp_do_updateV_nonblock.h $(SRCROOT)FMS/mpp/include/mpp_do_update.h $(SRCROOT)FMS/mpp/include/mpp_do_updateV.h $(SRCROOT)FMS/mpp/include/mpp_do_check.h $(SRCROOT)FMS/mpp/include/mpp_do_checkV.h $(SRCROOT)FMS/mpp/include/mpp_update_nest_domains.h $(SRCROOT)FMS/mpp/include/mpp_do_update_nest.h $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_ad.h $(SRCROOT)FMS/mpp/include/mpp_do_update_ad.h $(SRCROOT)FMS/mpp/include/mpp_do_updateV_ad.h $(SRCROOT)FMS/mpp/include/mpp_do_redistribute.h $(SRCROOT)FMS/mpp/include/mpp_get_boundary.h $(SRCROOT)FMS/mpp/include/mpp_get_boundary_ad.h $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary.h $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary_ad.h $(SRCROOT)FMS/mpp/include/mpp_group_update.h $(SRCROOT)FMS/mpp/include/group_update_pack.inc $(SRCROOT)FMS/mpp/include/group_update_unpack.inc $(SRCROOT)FMS/mpp/include/mpp_domains_reduce.inc $(SRCROOT)FMS/mpp/include/mpp_global_reduce.h $(SRCROOT)FMS/mpp/include/mpp_global_sum.h $(SRCROOT)FMS/mpp/include/mpp_global_sum_tl.h $(SRCROOT)FMS/mpp/include/mpp_global_sum_ad.h $(SRCROOT)FMS/mpp/include/mpp_global_field.h $(SRCROOT)FMS/mpp/include/mpp_global_field_ad.h $(SRCROOT)FMS/mpp/include/mpp_do_global_field.h $(SRCROOT)FMS/mpp/include/mpp_do_global_field_ad.h $(SRCROOT)FMS/mpp/include/mpp_unstruct_domain.inc $(SRCROOT)FMS/mpp/include/mpp_unstruct_pass_data.h $(SRCROOT)FMS/mpp/include/mpp_global_field_ug.h mpp_parameter.o mpp_data.o mpp.o mpp_memutils.o mpp_efp.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include -I$(SRCROOT)FMS/mpp/include	$(SRCROOT)FMS/mpp/mpp_domains.F90
mpp_efp.o: $(SRCROOT)FMS/mpp/mpp_efp.F90 mpp.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/mpp/mpp_efp.F90
mpp_io.o: $(SRCROOT)FMS/mpp/mpp_io.F90 $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/mpp/include/mpp_io_util.inc $(SRCROOT)FMS/mpp/include/mpp_io_misc.inc $(SRCROOT)FMS/mpp/include/mpp_io_connect.inc $(SRCROOT)FMS/mpp/include/mpp_io_read.inc $(SRCROOT)FMS/mpp/include/mpp_read_2Ddecomp.h $(SRCROOT)FMS/mpp/include/mpp_read_compressed.h $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.inc $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.h $(SRCROOT)FMS/mpp/include/mpp_io_write.inc $(SRCROOT)FMS/mpp/include/mpp_write_2Ddecomp.h $(SRCROOT)FMS/mpp/include/mpp_write_compressed.h $(SRCROOT)FMS/mpp/include/mpp_write_unlimited_axis.h $(SRCROOT)FMS/mpp/include/mpp_write.h $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_write.inc $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_read.inc mpp_parameter.o mpp.o mpp_domains.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include -I$(SRCROOT)FMS/mpp/include	$(SRCROOT)FMS/mpp/mpp_io.F90
mpp_memuse.o: $(SRCROOT)FMS/mpp/mpp_memuse.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)FMS/mpp/mpp_memuse.c
mpp_memutils.o: $(SRCROOT)FMS/mpp/mpp_memutils.F90 mpp.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/mpp/mpp_memutils.F90
mpp_parameter.o: $(SRCROOT)FMS/mpp/mpp_parameter.F90 $(SRCROOT)FMS/include/file_version.h platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/mpp/mpp_parameter.F90
mpp_utilities.o: $(SRCROOT)FMS/mpp/mpp_utilities.F90 $(SRCROOT)FMS/include/file_version.h mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/mpp/mpp_utilities.F90
netcdf_io.o: $(SRCROOT)FMS/fms2_io/netcdf_io.F90 $(SRCROOT)FMS/fms2_io/include/netcdf_add_restart_variable.inc $(SRCROOT)FMS/fms2_io/include/netcdf_read_data.inc $(SRCROOT)FMS/fms2_io/include/netcdf_write_data.inc $(SRCROOT)FMS/fms2_io/include/register_global_attribute.inc $(SRCROOT)FMS/fms2_io/include/register_variable_attribute.inc $(SRCROOT)FMS/fms2_io/include/get_global_attribute.inc $(SRCROOT)FMS/fms2_io/include/get_variable_attribute.inc $(SRCROOT)FMS/fms2_io/include/compressed_write.inc $(SRCROOT)FMS/fms2_io/include/compressed_read.inc mpp.o fms_io_utils.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/fms2_io/include	$(SRCROOT)FMS/fms2_io/netcdf_io.F90
platform.o: $(SRCROOT)FMS/platform/platform.F90 $(SRCROOT)FMS/include/fms_platform.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/platform/platform.F90
quicksort.o: $(SRCROOT)FMS/drifters/quicksort.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/drifters/quicksort.F90
random_numbers.o: $(SRCROOT)FMS/random_numbers/random_numbers.F90 mersennetwister.o time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/random_numbers/random_numbers.F90
read_mosaic.o: $(SRCROOT)FMS/mosaic/read_mosaic.c $(SRCROOT)FMS/mosaic/read_mosaic.h $(SRCROOT)FMS/mosaic/constant.h $(SRCROOT)FMS/mosaic/mosaic_util.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)FMS/mosaic	$(SRCROOT)FMS/mosaic/read_mosaic.c
sat_vapor_pres.o: $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)FMS/include/file_version.h constants.o fms.o mpp.o sat_vapor_pres_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres.F90
sat_vapor_pres_k.o: $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)FMS/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres_k.F90
stock_constants.o: $(SRCROOT)FMS/exchange/stock_constants.F90 $(SRCROOT)FMS/include/file_version.h mpp.o fms.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/exchange/stock_constants.F90
time_interp.o: $(SRCROOT)FMS/time_interp/time_interp.F90 $(SRCROOT)FMS/include/file_version.h time_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/time_interp/time_interp.F90
time_interp_external.o: $(SRCROOT)FMS/time_interp/time_interp_external.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h fms.o mpp.o mpp_io.o time_manager.o get_cal_time.o mpp_domains.o time_interp.o axis_utils.o platform.o horiz_interp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/time_interp/time_interp_external.F90
time_interp_external2.o: $(SRCROOT)FMS/time_interp/time_interp_external2.F90 $(SRCROOT)FMS/include/file_version.h platform.o fms.o mpp.o time_manager.o get_cal_time.o mpp_domains.o time_interp.o axis_utils2.o horiz_interp.o fms2_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/time_interp/time_interp_external2.F90
time_manager.o: $(SRCROOT)FMS/time_manager/time_manager.F90 $(SRCROOT)FMS/include/file_version.h platform.o constants.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/time_manager/time_manager.F90
topography.o: $(SRCROOT)FMS/topography/topography.F90 $(SRCROOT)FMS/include/file_version.h gaussian_topog.o horiz_interp.o fms.o fms2_io.o fms_io.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/topography/topography.F90
tracer_manager.o: $(SRCROOT)FMS/tracer_manager/tracer_manager.F90 $(SRCROOT)FMS/include/file_version.h mpp.o fms.o field_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/tracer_manager/tracer_manager.F90
tridiagonal.o: $(SRCROOT)FMS/tridiagonal/tridiagonal.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/tridiagonal/tridiagonal.F90
xgrid.o: $(SRCROOT)FMS/exchange/xgrid.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h fms.o mpp.o mpp_domains.o constants.o mosaic2.o stock_constants.o gradient.o fms2_io.o mpp_io.o fms_io.o mosaic.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/exchange/xgrid.F90
./gradient_c2l.h: $(SRCROOT)FMS/mosaic/gradient_c2l.h
	cp $(SRCROOT)FMS/mosaic/gradient_c2l.h .
./fms_switches.h: $(SRCROOT)FMS/drifters/fms_switches.h
	cp $(SRCROOT)FMS/drifters/fms_switches.h .
./compressed_read.inc: $(SRCROOT)FMS/fms2_io/include/compressed_read.inc
	cp $(SRCROOT)FMS/fms2_io/include/compressed_read.inc .
./diag_output.F90: $(SRCROOT)FMS/diag_manager/diag_output.F90
	cp $(SRCROOT)FMS/diag_manager/diag_output.F90 .
./mpp_util_mpi.inc: $(SRCROOT)FMS/mpp/include/mpp_util_mpi.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_util_mpi.inc .
./diag_data.F90: $(SRCROOT)FMS/diag_manager/diag_data.F90
	cp $(SRCROOT)FMS/diag_manager/diag_data.F90 .
./mpp_transmit.inc: $(SRCROOT)FMS/mpp/include/mpp_transmit.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_transmit.inc .
./get_data_type_string.inc: $(SRCROOT)FMS/fms2_io/include/get_data_type_string.inc
	cp $(SRCROOT)FMS/fms2_io/include/get_data_type_string.inc .
./fms_io_unstructured_get_file_name.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_name.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_name.inc .
./mpp_do_updateV_ad.h: $(SRCROOT)FMS/mpp/include/mpp_do_updateV_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_updateV_ad.h .
./mpp_get_boundary_ad.h: $(SRCROOT)FMS/mpp/include/mpp_get_boundary_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_get_boundary_ad.h .
./mpp_global_field_ad.h: $(SRCROOT)FMS/mpp/include/mpp_global_field_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_global_field_ad.h .
./diag_util.F90: $(SRCROOT)FMS/diag_manager/diag_util.F90
	cp $(SRCROOT)FMS/diag_manager/diag_util.F90 .
./mpp_unstruct_pass_data.h: $(SRCROOT)FMS/mpp/include/mpp_unstruct_pass_data.h
	cp $(SRCROOT)FMS/mpp/include/mpp_unstruct_pass_data.h .
./data_override.F90: $(SRCROOT)FMS/data_override/data_override.F90
	cp $(SRCROOT)FMS/data_override/data_override.F90 .
./mpp_error_a_s.h: $(SRCROOT)FMS/mpp/include/mpp_error_a_s.h
	cp $(SRCROOT)FMS/mpp/include/mpp_error_a_s.h .
./mpp_alltoall_nocomm.h: $(SRCROOT)FMS/mpp/include/mpp_alltoall_nocomm.h
	cp $(SRCROOT)FMS/mpp/include/mpp_alltoall_nocomm.h .
./fms_affinity.F90: $(SRCROOT)FMS/affinity/fms_affinity.F90
	cp $(SRCROOT)FMS/affinity/fms_affinity.F90 .
./mpp_chksum_int.h: $(SRCROOT)FMS/mpp/include/mpp_chksum_int.h
	cp $(SRCROOT)FMS/mpp/include/mpp_chksum_int.h .
./mpp_write.h: $(SRCROOT)FMS/mpp/include/mpp_write.h
	cp $(SRCROOT)FMS/mpp/include/mpp_write.h .
./read_data_4d.inc: $(SRCROOT)FMS/fms/read_data_4d.inc
	cp $(SRCROOT)FMS/fms/read_data_4d.inc .
./read_mosaic.c: $(SRCROOT)FMS/mosaic/read_mosaic.c
	cp $(SRCROOT)FMS/mosaic/read_mosaic.c .
./tracer_manager.F90: $(SRCROOT)FMS/tracer_manager/tracer_manager.F90
	cp $(SRCROOT)FMS/tracer_manager/tracer_manager.F90 .
./mpp_io_misc.inc: $(SRCROOT)FMS/mpp/include/mpp_io_misc.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_io_misc.inc .
./mpp_sum_mpi.h: $(SRCROOT)FMS/mpp/include/mpp_sum_mpi.h
	cp $(SRCROOT)FMS/mpp/include/mpp_sum_mpi.h .
./sat_vapor_pres.F90: $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres.F90
	cp $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres.F90 .
./mpp_do_updateV.h: $(SRCROOT)FMS/mpp/include/mpp_do_updateV.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_updateV.h .
./mpp_group_update.h: $(SRCROOT)FMS/mpp/include/mpp_group_update.h
	cp $(SRCROOT)FMS/mpp/include/mpp_group_update.h .
./mpp_write_2Ddecomp.h: $(SRCROOT)FMS/mpp/include/mpp_write_2Ddecomp.h
	cp $(SRCROOT)FMS/mpp/include/mpp_write_2Ddecomp.h .
./mpp_do_check.h: $(SRCROOT)FMS/mpp/include/mpp_do_check.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_check.h .
./read_data_3d.inc: $(SRCROOT)FMS/fms/read_data_3d.inc
	cp $(SRCROOT)FMS/fms/read_data_3d.inc .
./ensemble_manager.F90: $(SRCROOT)FMS/coupler/ensemble_manager.F90
	cp $(SRCROOT)FMS/coupler/ensemble_manager.F90 .
./drifters_comm.F90: $(SRCROOT)FMS/drifters/drifters_comm.F90
	cp $(SRCROOT)FMS/drifters/drifters_comm.F90 .
./write_data.inc: $(SRCROOT)FMS/fms/write_data.inc
	cp $(SRCROOT)FMS/fms/write_data.inc .
./mosaic2.F90: $(SRCROOT)FMS/mosaic/mosaic2.F90
	cp $(SRCROOT)FMS/mosaic/mosaic2.F90 .
./mpp_do_updateV_nonblock.h: $(SRCROOT)FMS/mpp/include/mpp_do_updateV_nonblock.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_updateV_nonblock.h .
./mersennetwister.F90: $(SRCROOT)FMS/random_numbers/mersennetwister.F90
	cp $(SRCROOT)FMS/random_numbers/mersennetwister.F90 .
./mpp_io_connect.inc: $(SRCROOT)FMS/mpp/include/mpp_io_connect.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_io_connect.inc .
./astronomy.F90: $(SRCROOT)FMS/astronomy/astronomy.F90
	cp $(SRCROOT)FMS/astronomy/astronomy.F90 .
./horiz_interp_spherical.F90: $(SRCROOT)FMS/horiz_interp/horiz_interp_spherical.F90
	cp $(SRCROOT)FMS/horiz_interp/horiz_interp_spherical.F90 .
./mosaic.F90: $(SRCROOT)FMS/mosaic/mosaic.F90
	cp $(SRCROOT)FMS/mosaic/mosaic.F90 .
./cloud_interpolator.F90: $(SRCROOT)FMS/drifters/cloud_interpolator.F90
	cp $(SRCROOT)FMS/drifters/cloud_interpolator.F90 .
./mpp_sum_mpi_ad.h: $(SRCROOT)FMS/mpp/include/mpp_sum_mpi_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_sum_mpi_ad.h .
./mpp_domains_comm.inc: $(SRCROOT)FMS/mpp/include/mpp_domains_comm.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_domains_comm.inc .
./mpp_error_s_a.h: $(SRCROOT)FMS/mpp/include/mpp_error_s_a.h
	cp $(SRCROOT)FMS/mpp/include/mpp_error_s_a.h .
./axis_utils2.F90: $(SRCROOT)FMS/axis_utils/axis_utils2.F90
	cp $(SRCROOT)FMS/axis_utils/axis_utils2.F90 .
./drifters_compute_k.h: $(SRCROOT)FMS/drifters/drifters_compute_k.h
	cp $(SRCROOT)FMS/drifters/drifters_compute_k.h .
./mpp_update_domains2D_ad.h: $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_ad.h .
./monin_obukhov.F90: $(SRCROOT)FMS/monin_obukhov/monin_obukhov.F90
	cp $(SRCROOT)FMS/monin_obukhov/monin_obukhov.F90 .
./fms_io_unstructured_save_restart.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_save_restart.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_save_restart.inc .
./horiz_interp_bicubic.F90: $(SRCROOT)FMS/horiz_interp/horiz_interp_bicubic.F90
	cp $(SRCROOT)FMS/horiz_interp/horiz_interp_bicubic.F90 .
./topography.F90: $(SRCROOT)FMS/topography/topography.F90
	cp $(SRCROOT)FMS/topography/topography.F90 .
./netcdf_io.F90: $(SRCROOT)FMS/fms2_io/netcdf_io.F90
	cp $(SRCROOT)FMS/fms2_io/netcdf_io.F90 .
./unstructured_domain_read.inc: $(SRCROOT)FMS/fms2_io/include/unstructured_domain_read.inc
	cp $(SRCROOT)FMS/fms2_io/include/unstructured_domain_read.inc .
./register_domain_restart_variable.inc: $(SRCROOT)FMS/fms2_io/include/register_domain_restart_variable.inc
	cp $(SRCROOT)FMS/fms2_io/include/register_domain_restart_variable.inc .
./fms_io_unstructured_read.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_read.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_read.inc .
./horiz_interp_conserve.F90: $(SRCROOT)FMS/horiz_interp/horiz_interp_conserve.F90
	cp $(SRCROOT)FMS/horiz_interp/horiz_interp_conserve.F90 .
./mpp_chksum.h: $(SRCROOT)FMS/mpp/include/mpp_chksum.h
	cp $(SRCROOT)FMS/mpp/include/mpp_chksum.h .
./mpp_io_unstructured_write.inc: $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_write.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_write.inc .
./amip_interp.F90: $(SRCROOT)FMS/amip_interp/amip_interp.F90
	cp $(SRCROOT)FMS/amip_interp/amip_interp.F90 .
./mpp_read_compressed.h: $(SRCROOT)FMS/mpp/include/mpp_read_compressed.h
	cp $(SRCROOT)FMS/mpp/include/mpp_read_compressed.h .
./drifters_input.F90: $(SRCROOT)FMS/drifters/drifters_input.F90
	cp $(SRCROOT)FMS/drifters/drifters_input.F90 .
./netcdf_read_data.inc: $(SRCROOT)FMS/fms2_io/include/netcdf_read_data.inc
	cp $(SRCROOT)FMS/fms2_io/include/netcdf_read_data.inc .
./mpp_scatter.h: $(SRCROOT)FMS/mpp/include/mpp_scatter.h
	cp $(SRCROOT)FMS/mpp/include/mpp_scatter.h .
./mpp_data_mpi.inc: $(SRCROOT)FMS/mpp/include/mpp_data_mpi.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_data_mpi.inc .
./gaussian_topog.F90: $(SRCROOT)FMS/topography/gaussian_topog.F90
	cp $(SRCROOT)FMS/topography/gaussian_topog.F90 .
./get_variable_attribute.inc: $(SRCROOT)FMS/fms2_io/include/get_variable_attribute.inc
	cp $(SRCROOT)FMS/fms2_io/include/get_variable_attribute.inc .
./fms_io_unstructured_file_unit.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_file_unit.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_file_unit.inc .
./mpp_do_update_nest.h: $(SRCROOT)FMS/mpp/include/mpp_do_update_nest.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_update_nest.h .
./get_global_attribute.inc: $(SRCROOT)FMS/fms2_io/include/get_global_attribute.inc
	cp $(SRCROOT)FMS/fms2_io/include/get_global_attribute.inc .
./mpp_util_nocomm.inc: $(SRCROOT)FMS/mpp/include/mpp_util_nocomm.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_util_nocomm.inc .
./mpp_reduce_nocomm.h: $(SRCROOT)FMS/mpp/include/mpp_reduce_nocomm.h
	cp $(SRCROOT)FMS/mpp/include/mpp_reduce_nocomm.h .
./create_xgrid.c: $(SRCROOT)FMS/mosaic/create_xgrid.c
	cp $(SRCROOT)FMS/mosaic/create_xgrid.c .
./mpp_utilities.F90: $(SRCROOT)FMS/mpp/mpp_utilities.F90
	cp $(SRCROOT)FMS/mpp/mpp_utilities.F90 .
./fms_io_utils.F90: $(SRCROOT)FMS/fms2_io/fms_io_utils.F90
	cp $(SRCROOT)FMS/fms2_io/fms_io_utils.F90 .
./group_update_unpack.inc: $(SRCROOT)FMS/mpp/include/group_update_unpack.inc
	cp $(SRCROOT)FMS/mpp/include/group_update_unpack.inc .
./mpp_do_checkV.h: $(SRCROOT)FMS/mpp/include/mpp_do_checkV.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_checkV.h .
./drifters_push.h: $(SRCROOT)FMS/drifters/drifters_push.h
	cp $(SRCROOT)FMS/drifters/drifters_push.h .
./fms_io_unstructured_register_restart_axis.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_axis.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_axis.inc .
./mpp_sum_ad.inc: $(SRCROOT)FMS/mpp/include/mpp_sum_ad.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_sum_ad.inc .
./mpp_data.F90: $(SRCROOT)FMS/mpp/mpp_data.F90
	cp $(SRCROOT)FMS/mpp/mpp_data.F90 .
./fms_io_unstructured_get_field_size.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_get_field_size.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_get_field_size.inc .
./get_cal_time.F90: $(SRCROOT)FMS/time_manager/get_cal_time.F90
	cp $(SRCROOT)FMS/time_manager/get_cal_time.F90 .
./mpp_domains.F90: $(SRCROOT)FMS/mpp/mpp_domains.F90
	cp $(SRCROOT)FMS/mpp/mpp_domains.F90 .
./mpp_do_update_nonblock.h: $(SRCROOT)FMS/mpp/include/mpp_do_update_nonblock.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_update_nonblock.h .
./column_diagnostics.F90: $(SRCROOT)FMS/column_diagnostics/column_diagnostics.F90
	cp $(SRCROOT)FMS/column_diagnostics/column_diagnostics.F90 .
./mpp_update_domains2D_nonblock.h: $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_nonblock.h
	cp $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_nonblock.h .
./mpp_comm_mpi.inc: $(SRCROOT)FMS/mpp/include/mpp_comm_mpi.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_comm_mpi.inc .
./constants.F90: $(SRCROOT)FMS/constants/constants.F90
	cp $(SRCROOT)FMS/constants/constants.F90 .
./mpp_sum.inc: $(SRCROOT)FMS/mpp/include/mpp_sum.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_sum.inc .
./register_variable_attribute.inc: $(SRCROOT)FMS/fms2_io/include/register_variable_attribute.inc
	cp $(SRCROOT)FMS/fms2_io/include/register_variable_attribute.inc .
./mpp_io_unstructured_read.inc: $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_read.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_read.inc .
./mpp_do_global_field.h: $(SRCROOT)FMS/mpp/include/mpp_do_global_field.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_global_field.h .
./mpp_gather.h: $(SRCROOT)FMS/mpp/include/mpp_gather.h
	cp $(SRCROOT)FMS/mpp/include/mpp_gather.h .
./system_clock.h: $(SRCROOT)FMS/mpp/include/system_clock.h
	cp $(SRCROOT)FMS/mpp/include/system_clock.h .
./mpp_get_boundary.h: $(SRCROOT)FMS/mpp/include/mpp_get_boundary.h
	cp $(SRCROOT)FMS/mpp/include/mpp_get_boundary.h .
./time_interp_external.F90: $(SRCROOT)FMS/time_interp/time_interp_external.F90
	cp $(SRCROOT)FMS/time_interp/time_interp_external.F90 .
./mpp_read_2Ddecomp.h: $(SRCROOT)FMS/mpp/include/mpp_read_2Ddecomp.h
	cp $(SRCROOT)FMS/mpp/include/mpp_read_2Ddecomp.h .
./mpp_read_distributed_ascii.inc: $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.inc .
./get_grid_version_mpp.F90: $(SRCROOT)FMS/data_override/get_grid_version_mpp.F90
	cp $(SRCROOT)FMS/data_override/get_grid_version_mpp.F90 .
./array_utils_char.inc: $(SRCROOT)FMS/fms2_io/include/array_utils_char.inc
	cp $(SRCROOT)FMS/fms2_io/include/array_utils_char.inc .
./fm_util.F90: $(SRCROOT)FMS/field_manager/fm_util.F90
	cp $(SRCROOT)FMS/field_manager/fm_util.F90 .
./mpp_io_write.inc: $(SRCROOT)FMS/mpp/include/mpp_io_write.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_io_write.inc .
./netcdf_write_data.inc: $(SRCROOT)FMS/fms2_io/include/netcdf_write_data.inc
	cp $(SRCROOT)FMS/fms2_io/include/netcdf_write_data.inc .
./compressed_write.inc: $(SRCROOT)FMS/fms2_io/include/compressed_write.inc
	cp $(SRCROOT)FMS/fms2_io/include/compressed_write.inc .
./mpp_io.F90: $(SRCROOT)FMS/mpp/mpp_io.F90
	cp $(SRCROOT)FMS/mpp/mpp_io.F90 .
./create_xgrid.h: $(SRCROOT)FMS/mosaic/create_xgrid.h
	cp $(SRCROOT)FMS/mosaic/create_xgrid.h .
./mpp_read_distributed_ascii.h: $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.h
	cp $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.h .
./axis_utils.F90: $(SRCROOT)FMS/axis_utils/axis_utils.F90
	cp $(SRCROOT)FMS/axis_utils/axis_utils.F90 .
./mpp_io_util.inc: $(SRCROOT)FMS/mpp/include/mpp_io_util.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_io_util.inc .
./gradient_c2l.c: $(SRCROOT)FMS/mosaic/gradient_c2l.c
	cp $(SRCROOT)FMS/mosaic/gradient_c2l.c .
./interp.c: $(SRCROOT)FMS/mosaic/interp.c
	cp $(SRCROOT)FMS/mosaic/interp.c .
./horiz_interp_type.F90: $(SRCROOT)FMS/horiz_interp/horiz_interp_type.F90
	cp $(SRCROOT)FMS/horiz_interp/horiz_interp_type.F90 .
./platform.F90: $(SRCROOT)FMS/platform/platform.F90
	cp $(SRCROOT)FMS/platform/platform.F90 .
./mpp_global_sum.h: $(SRCROOT)FMS/mpp/include/mpp_global_sum.h
	cp $(SRCROOT)FMS/mpp/include/mpp_global_sum.h .
./unstructured_domain_write.inc: $(SRCROOT)FMS/fms2_io/include/unstructured_domain_write.inc
	cp $(SRCROOT)FMS/fms2_io/include/unstructured_domain_write.inc .
./tridiagonal.F90: $(SRCROOT)FMS/tridiagonal/tridiagonal.F90
	cp $(SRCROOT)FMS/tridiagonal/tridiagonal.F90 .
./mpp_io_read.inc: $(SRCROOT)FMS/mpp/include/mpp_io_read.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_io_read.inc .
./diag_axis.F90: $(SRCROOT)FMS/diag_manager/diag_axis.F90
	cp $(SRCROOT)FMS/diag_manager/diag_axis.F90 .
./mpp_do_update_ad.h: $(SRCROOT)FMS/mpp/include/mpp_do_update_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_update_ad.h .
./domain_read.inc: $(SRCROOT)FMS/fms2_io/include/domain_read.inc
	cp $(SRCROOT)FMS/fms2_io/include/domain_read.inc .
./quicksort.F90: $(SRCROOT)FMS/drifters/quicksort.F90
	cp $(SRCROOT)FMS/drifters/quicksort.F90 .
./mpp_comm.inc: $(SRCROOT)FMS/mpp/include/mpp_comm.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_comm.inc .
./time_interp_external2.F90: $(SRCROOT)FMS/time_interp/time_interp_external2.F90
	cp $(SRCROOT)FMS/time_interp/time_interp_external2.F90 .
./fms_platform.h: $(SRCROOT)FMS/include/fms_platform.h
	cp $(SRCROOT)FMS/include/fms_platform.h .
./grid.F90: $(SRCROOT)FMS/mosaic/grid.F90
	cp $(SRCROOT)FMS/mosaic/grid.F90 .
./mpp_update_nest_domains.h: $(SRCROOT)FMS/mpp/include/mpp_update_nest_domains.h
	cp $(SRCROOT)FMS/mpp/include/mpp_update_nest_domains.h .
./drifters_core.F90: $(SRCROOT)FMS/drifters/drifters_core.F90
	cp $(SRCROOT)FMS/drifters/drifters_core.F90 .
./mpp_chksum_scalar.h: $(SRCROOT)FMS/mpp/include/mpp_chksum_scalar.h
	cp $(SRCROOT)FMS/mpp/include/mpp_chksum_scalar.h .
./memutils.F90: $(SRCROOT)FMS/memutils/memutils.F90
	cp $(SRCROOT)FMS/memutils/memutils.F90 .
./mpp_transmit_mpi.h: $(SRCROOT)FMS/mpp/include/mpp_transmit_mpi.h
	cp $(SRCROOT)FMS/mpp/include/mpp_transmit_mpi.h .
./mpp_do_get_boundary.h: $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary.h .
./drifters.F90: $(SRCROOT)FMS/drifters/drifters.F90
	cp $(SRCROOT)FMS/drifters/drifters.F90 .
./diag_integral.F90: $(SRCROOT)FMS/diag_integral/diag_integral.F90
	cp $(SRCROOT)FMS/diag_integral/diag_integral.F90 .
./register_unstructured_domain_restart_variable.inc: $(SRCROOT)FMS/fms2_io/include/register_unstructured_domain_restart_variable.inc
	cp $(SRCROOT)FMS/fms2_io/include/register_unstructured_domain_restart_variable.inc .
./mpp_error_a_a.h: $(SRCROOT)FMS/mpp/include/mpp_error_a_a.h
	cp $(SRCROOT)FMS/mpp/include/mpp_error_a_a.h .
./atmos_ocean_fluxes.F90: $(SRCROOT)FMS/coupler/atmos_ocean_fluxes.F90
	cp $(SRCROOT)FMS/coupler/atmos_ocean_fluxes.F90 .
./fms_io_unstructured_get_file_unit.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_unit.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_unit.inc .
./horiz_interp.F90: $(SRCROOT)FMS/horiz_interp/horiz_interp.F90
	cp $(SRCROOT)FMS/horiz_interp/horiz_interp.F90 .
./interp.h: $(SRCROOT)FMS/mosaic/interp.h
	cp $(SRCROOT)FMS/mosaic/interp.h .
./diag_manager.F90: $(SRCROOT)FMS/diag_manager/diag_manager.F90
	cp $(SRCROOT)FMS/diag_manager/diag_manager.F90 .
./time_interp.F90: $(SRCROOT)FMS/time_interp/time_interp.F90
	cp $(SRCROOT)FMS/time_interp/time_interp.F90 .
./read_mosaic.h: $(SRCROOT)FMS/mosaic/read_mosaic.h
	cp $(SRCROOT)FMS/mosaic/read_mosaic.h .
./mpp_unstruct_domain.inc: $(SRCROOT)FMS/mpp/include/mpp_unstruct_domain.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_unstruct_domain.inc .
./read_data_2d.inc: $(SRCROOT)FMS/fms/read_data_2d.inc
	cp $(SRCROOT)FMS/fms/read_data_2d.inc .
./mpp_sum_nocomm_ad.h: $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm_ad.h .
./file_version.h: $(SRCROOT)FMS/include/file_version.h
	cp $(SRCROOT)FMS/include/file_version.h .
./sat_vapor_pres_k.F90: $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres_k.F90
	cp $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres_k.F90 .
./diag_table.F90: $(SRCROOT)FMS/diag_manager/diag_table.F90
	cp $(SRCROOT)FMS/diag_manager/diag_table.F90 .
./mpp.F90: $(SRCROOT)FMS/mpp/mpp.F90
	cp $(SRCROOT)FMS/mpp/mpp.F90 .
./fms_io_unstructured_field_exist.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_field_exist.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_field_exist.inc .
./array_utils.inc: $(SRCROOT)FMS/fms2_io/include/array_utils.inc
	cp $(SRCROOT)FMS/fms2_io/include/array_utils.inc .
./fms_netcdf_domain_io.F90: $(SRCROOT)FMS/fms2_io/fms_netcdf_domain_io.F90
	cp $(SRCROOT)FMS/fms2_io/fms_netcdf_domain_io.F90 .
./get_grid_version_fms2io.F90: $(SRCROOT)FMS/data_override/get_grid_version_fms2io.F90
	cp $(SRCROOT)FMS/data_override/get_grid_version_fms2io.F90 .
./interpolator.F90: $(SRCROOT)FMS/interpolator/interpolator.F90
	cp $(SRCROOT)FMS/interpolator/interpolator.F90 .
./field_manager.F90: $(SRCROOT)FMS/field_manager/field_manager.F90
	cp $(SRCROOT)FMS/field_manager/field_manager.F90 .
./fms_io_unstructured_register_restart_field.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_field.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_field.inc .
./fms.F90: $(SRCROOT)FMS/fms/fms.F90
	cp $(SRCROOT)FMS/fms/fms.F90 .
./fms_io.F90: $(SRCROOT)FMS/fms/fms_io.F90
	cp $(SRCROOT)FMS/fms/fms_io.F90 .
./fms_io_unstructured_setup_one_field.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_setup_one_field.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_setup_one_field.inc .
./mpp_error_s_s.h: $(SRCROOT)FMS/mpp/include/mpp_error_s_s.h
	cp $(SRCROOT)FMS/mpp/include/mpp_error_s_s.h .
./mpp_domains_define.inc: $(SRCROOT)FMS/mpp/include/mpp_domains_define.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_domains_define.inc .
./gradient.F90: $(SRCROOT)FMS/mosaic/gradient.F90
	cp $(SRCROOT)FMS/mosaic/gradient.F90 .
./mpp_domains_reduce.inc: $(SRCROOT)FMS/mpp/include/mpp_domains_reduce.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_domains_reduce.inc .
./fms2_io.F90: $(SRCROOT)FMS/fms2_io/fms2_io.F90
	cp $(SRCROOT)FMS/fms2_io/fms2_io.F90 .
./fms_netcdf_unstructured_domain_io.F90: $(SRCROOT)FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90
	cp $(SRCROOT)FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90 .
./drifters_io.F90: $(SRCROOT)FMS/drifters/drifters_io.F90
	cp $(SRCROOT)FMS/drifters/drifters_io.F90 .
./mpp_domains_misc.inc: $(SRCROOT)FMS/mpp/include/mpp_domains_misc.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_domains_misc.inc .
./mpp_sum_nocomm.h: $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm.h
	cp $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm.h .
./mosaic_util.c: $(SRCROOT)FMS/mosaic/mosaic_util.c
	cp $(SRCROOT)FMS/mosaic/mosaic_util.c .
./group_update_pack.inc: $(SRCROOT)FMS/mpp/include/group_update_pack.inc
	cp $(SRCROOT)FMS/mpp/include/group_update_pack.inc .
./mpp_reduce_mpi.h: $(SRCROOT)FMS/mpp/include/mpp_reduce_mpi.h
	cp $(SRCROOT)FMS/mpp/include/mpp_reduce_mpi.h .
./random_numbers.F90: $(SRCROOT)FMS/random_numbers/random_numbers.F90
	cp $(SRCROOT)FMS/random_numbers/random_numbers.F90 .
./mpp_do_redistribute.h: $(SRCROOT)FMS/mpp/include/mpp_do_redistribute.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_redistribute.h .
./monin_obukhov_inter.F90: $(SRCROOT)FMS/monin_obukhov/monin_obukhov_inter.F90
	cp $(SRCROOT)FMS/monin_obukhov/monin_obukhov_inter.F90 .
./mpp_global_field.h: $(SRCROOT)FMS/mpp/include/mpp_global_field.h
	cp $(SRCROOT)FMS/mpp/include/mpp_global_field.h .
./mpp_transmit_nocomm.h: $(SRCROOT)FMS/mpp/include/mpp_transmit_nocomm.h
	cp $(SRCROOT)FMS/mpp/include/mpp_transmit_nocomm.h .
./mpp_global_field_ug.h: $(SRCROOT)FMS/mpp/include/mpp_global_field_ug.h
	cp $(SRCROOT)FMS/mpp/include/mpp_global_field_ug.h .
./coupler_types.F90: $(SRCROOT)FMS/coupler/coupler_types.F90
	cp $(SRCROOT)FMS/coupler/coupler_types.F90 .
./mpp_type_mpi.h: $(SRCROOT)FMS/mpp/include/mpp_type_mpi.h
	cp $(SRCROOT)FMS/mpp/include/mpp_type_mpi.h .
./register_global_attribute.inc: $(SRCROOT)FMS/fms2_io/include/register_global_attribute.inc
	cp $(SRCROOT)FMS/fms2_io/include/register_global_attribute.inc .
./mpp_memutils.F90: $(SRCROOT)FMS/mpp/mpp_memutils.F90
	cp $(SRCROOT)FMS/mpp/mpp_memutils.F90 .
./mpp_do_get_boundary_ad.h: $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary_ad.h .
./constant.h: $(SRCROOT)FMS/mosaic/constant.h
	cp $(SRCROOT)FMS/mosaic/constant.h .
./mpp_do_global_field_ad.h: $(SRCROOT)FMS/mpp/include/mpp_do_global_field_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_global_field_ad.h .
./domain_write.inc: $(SRCROOT)FMS/fms2_io/include/domain_write.inc
	cp $(SRCROOT)FMS/fms2_io/include/domain_write.inc .
./mpp_data_nocomm.inc: $(SRCROOT)FMS/mpp/include/mpp_data_nocomm.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_data_nocomm.inc .
./stock_constants.F90: $(SRCROOT)FMS/exchange/stock_constants.F90
	cp $(SRCROOT)FMS/exchange/stock_constants.F90 .
./parse.inc: $(SRCROOT)FMS/field_manager/parse.inc
	cp $(SRCROOT)FMS/field_manager/parse.inc .
./mpp_parameter.F90: $(SRCROOT)FMS/mpp/mpp_parameter.F90
	cp $(SRCROOT)FMS/mpp/mpp_parameter.F90 .
./netcdf_add_restart_variable.inc: $(SRCROOT)FMS/fms2_io/include/netcdf_add_restart_variable.inc
	cp $(SRCROOT)FMS/fms2_io/include/netcdf_add_restart_variable.inc .
./block_control.F90: $(SRCROOT)FMS/block_control/block_control.F90
	cp $(SRCROOT)FMS/block_control/block_control.F90 .
./mpp_type_nocomm.h: $(SRCROOT)FMS/mpp/include/mpp_type_nocomm.h
	cp $(SRCROOT)FMS/mpp/include/mpp_type_nocomm.h .
./mpp_comm_nocomm.inc: $(SRCROOT)FMS/mpp/include/mpp_comm_nocomm.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_comm_nocomm.inc .
./drifters_set_field.h: $(SRCROOT)FMS/drifters/drifters_set_field.h
	cp $(SRCROOT)FMS/drifters/drifters_set_field.h .
./mpp_define_nest_domains.inc: $(SRCROOT)FMS/mpp/include/mpp_define_nest_domains.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_define_nest_domains.inc .
./mpp_write_unlimited_axis.h: $(SRCROOT)FMS/mpp/include/mpp_write_unlimited_axis.h
	cp $(SRCROOT)FMS/mpp/include/mpp_write_unlimited_axis.h .
./get_checksum.inc: $(SRCROOT)FMS/fms2_io/include/get_checksum.inc
	cp $(SRCROOT)FMS/fms2_io/include/get_checksum.inc .
./diag_grid.F90: $(SRCROOT)FMS/diag_manager/diag_grid.F90
	cp $(SRCROOT)FMS/diag_manager/diag_grid.F90 .
./xgrid.F90: $(SRCROOT)FMS/exchange/xgrid.F90
	cp $(SRCROOT)FMS/exchange/xgrid.F90 .
./time_manager.F90: $(SRCROOT)FMS/time_manager/time_manager.F90
	cp $(SRCROOT)FMS/time_manager/time_manager.F90 .
./compute_global_checksum.inc: $(SRCROOT)FMS/fms2_io/include/compute_global_checksum.inc
	cp $(SRCROOT)FMS/fms2_io/include/compute_global_checksum.inc .
./mpp_domains_util.inc: $(SRCROOT)FMS/mpp/include/mpp_domains_util.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_domains_util.inc .
./affinity.c: $(SRCROOT)FMS/affinity/affinity.c
	cp $(SRCROOT)FMS/affinity/affinity.c .
./blackboxio.F90: $(SRCROOT)FMS/fms2_io/blackboxio.F90
	cp $(SRCROOT)FMS/fms2_io/blackboxio.F90 .
./mpp_alltoall_mpi.h: $(SRCROOT)FMS/mpp/include/mpp_alltoall_mpi.h
	cp $(SRCROOT)FMS/mpp/include/mpp_alltoall_mpi.h .
./mpp_do_update.h: $(SRCROOT)FMS/mpp/include/mpp_do_update.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_update.h .
./mpp_update_domains2D.h: $(SRCROOT)FMS/mpp/include/mpp_update_domains2D.h
	cp $(SRCROOT)FMS/mpp/include/mpp_update_domains2D.h .
./mpp_memuse.c: $(SRCROOT)FMS/mpp/mpp_memuse.c
	cp $(SRCROOT)FMS/mpp/mpp_memuse.c .
./mpp_util.inc: $(SRCROOT)FMS/mpp/include/mpp_util.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_util.inc .
./mpp_efp.F90: $(SRCROOT)FMS/mpp/mpp_efp.F90
	cp $(SRCROOT)FMS/mpp/mpp_efp.F90 .
./mpp_global_sum_tl.h: $(SRCROOT)FMS/mpp/include/mpp_global_sum_tl.h
	cp $(SRCROOT)FMS/mpp/include/mpp_global_sum_tl.h .
./mosaic_util.h: $(SRCROOT)FMS/mosaic/mosaic_util.h
	cp $(SRCROOT)FMS/mosaic/mosaic_util.h .
./mpp_global_reduce.h: $(SRCROOT)FMS/mpp/include/mpp_global_reduce.h
	cp $(SRCROOT)FMS/mpp/include/mpp_global_reduce.h .
./mpp_write_compressed.h: $(SRCROOT)FMS/mpp/include/mpp_write_compressed.h
	cp $(SRCROOT)FMS/mpp/include/mpp_write_compressed.h .
./horiz_interp_bilinear.F90: $(SRCROOT)FMS/horiz_interp/horiz_interp_bilinear.F90
	cp $(SRCROOT)FMS/horiz_interp/horiz_interp_bilinear.F90 .
./mpp_global_sum_ad.h: $(SRCROOT)FMS/mpp/include/mpp_global_sum_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_global_sum_ad.h .
SRC = $(SRCROOT)FMS/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)FMS/mpp/mpp_data.F90 $(SRCROOT)FMS/mosaic/mosaic2.F90 $(SRCROOT)FMS/diag_manager/diag_axis.F90 $(SRCROOT)FMS/fms2_io/fms2_io.F90 $(SRCROOT)FMS/amip_interp/amip_interp.F90 $(SRCROOT)FMS/axis_utils/axis_utils.F90 $(SRCROOT)FMS/column_diagnostics/column_diagnostics.F90 $(SRCROOT)FMS/platform/platform.F90 $(SRCROOT)FMS/mosaic/mosaic_util.c $(SRCROOT)FMS/data_override/get_grid_version_fms2io.F90 $(SRCROOT)FMS/fms2_io/fms_netcdf_domain_io.F90 $(SRCROOT)FMS/exchange/xgrid.F90 $(SRCROOT)FMS/diag_manager/diag_manager.F90 $(SRCROOT)FMS/constants/constants.F90 $(SRCROOT)FMS/diag_manager/diag_util.F90 $(SRCROOT)FMS/mpp/mpp_memuse.c $(SRCROOT)FMS/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)FMS/field_manager/field_manager.F90 $(SRCROOT)FMS/data_override/get_grid_version_mpp.F90 $(SRCROOT)FMS/diag_manager/diag_grid.F90 $(SRCROOT)FMS/time_manager/get_cal_time.F90 $(SRCROOT)FMS/coupler/ensemble_manager.F90 $(SRCROOT)FMS/drifters/drifters_input.F90 $(SRCROOT)FMS/fms2_io/blackboxio.F90 $(SRCROOT)FMS/random_numbers/random_numbers.F90 $(SRCROOT)FMS/fms2_io/fms_io_utils.F90 $(SRCROOT)FMS/mpp/mpp_utilities.F90 $(SRCROOT)FMS/block_control/block_control.F90 $(SRCROOT)FMS/diag_manager/diag_output.F90 $(SRCROOT)FMS/time_manager/time_manager.F90 $(SRCROOT)FMS/tridiagonal/tridiagonal.F90 $(SRCROOT)FMS/random_numbers/mersennetwister.F90 $(SRCROOT)FMS/mosaic/read_mosaic.c $(SRCROOT)FMS/drifters/drifters_core.F90 $(SRCROOT)FMS/drifters/quicksort.F90 $(SRCROOT)FMS/fms2_io/netcdf_io.F90 $(SRCROOT)FMS/diag_integral/diag_integral.F90 $(SRCROOT)FMS/fms/fms.F90 $(SRCROOT)FMS/astronomy/astronomy.F90 $(SRCROOT)FMS/mosaic/gradient.F90 $(SRCROOT)FMS/mpp/mpp_memutils.F90 $(SRCROOT)FMS/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)FMS/affinity/fms_affinity.F90 $(SRCROOT)FMS/mosaic/create_xgrid.c $(SRCROOT)FMS/drifters/drifters_io.F90 $(SRCROOT)FMS/monin_obukhov/monin_obukhov.F90 $(SRCROOT)FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90 $(SRCROOT)FMS/diag_manager/diag_data.F90 $(SRCROOT)FMS/field_manager/fm_util.F90 $(SRCROOT)FMS/exchange/stock_constants.F90 $(SRCROOT)FMS/fms/fms_io.F90 $(SRCROOT)FMS/drifters/drifters.F90 $(SRCROOT)FMS/mosaic/interp.c $(SRCROOT)FMS/topography/gaussian_topog.F90 $(SRCROOT)FMS/mosaic/gradient_c2l.c $(SRCROOT)FMS/mosaic/mosaic.F90 $(SRCROOT)FMS/horiz_interp/horiz_interp.F90 $(SRCROOT)FMS/time_interp/time_interp_external.F90 $(SRCROOT)FMS/mpp/mpp_domains.F90 $(SRCROOT)FMS/data_override/data_override.F90 $(SRCROOT)FMS/coupler/coupler_types.F90 $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)FMS/axis_utils/axis_utils2.F90 $(SRCROOT)FMS/time_interp/time_interp_external2.F90 $(SRCROOT)FMS/mpp/mpp_parameter.F90 $(SRCROOT)FMS/monin_obukhov/monin_obukhov_inter.F90 $(SRCROOT)FMS/mpp/mpp_efp.F90 $(SRCROOT)FMS/mpp/mpp.F90 $(SRCROOT)FMS/tracer_manager/tracer_manager.F90 $(SRCROOT)FMS/drifters/drifters_comm.F90 $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)FMS/memutils/memutils.F90 $(SRCROOT)FMS/diag_manager/diag_table.F90 $(SRCROOT)FMS/topography/topography.F90 $(SRCROOT)FMS/interpolator/interpolator.F90 $(SRCROOT)FMS/horiz_interp/horiz_interp_type.F90 $(SRCROOT)FMS/mosaic/grid.F90 $(SRCROOT)FMS/mpp/mpp_io.F90 $(SRCROOT)FMS/time_interp/time_interp.F90 $(SRCROOT)FMS/drifters/cloud_interpolator.F90 $(SRCROOT)FMS/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)FMS/affinity/affinity.c $(SRCROOT)FMS/coupler/atmos_ocean_fluxes.F90 $(SRCROOT)FMS/mpp/include/mpp_gather.h $(SRCROOT)FMS/mpp/include/mpp_global_sum_ad.h $(SRCROOT)FMS/mpp/include/mpp_get_boundary.h $(SRCROOT)FMS/mpp/include/system_clock.h $(SRCROOT)FMS/mpp/include/mpp_global_sum_tl.h $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_read.inc $(SRCROOT)FMS/mpp/include/mpp_util.inc $(SRCROOT)FMS/mpp/include/mpp_do_global_field.h $(SRCROOT)FMS/mpp/include/mpp_write_compressed.h $(SRCROOT)FMS/mpp/include/mpp_global_reduce.h $(SRCROOT)FMS/mosaic/mosaic_util.h $(SRCROOT)FMS/fms2_io/include/register_variable_attribute.inc $(SRCROOT)FMS/mpp/include/mpp_do_update.h $(SRCROOT)FMS/mpp/include/mpp_alltoall_mpi.h $(SRCROOT)FMS/mpp/include/mpp_sum.inc $(SRCROOT)FMS/mpp/include/mpp_update_domains2D.h $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_nonblock.h $(SRCROOT)FMS/fms2_io/include/compute_global_checksum.inc $(SRCROOT)FMS/mpp/include/mpp_domains_util.inc $(SRCROOT)FMS/mpp/include/mpp_comm_mpi.inc $(SRCROOT)FMS/mpp/include/mpp_do_update_nonblock.h $(SRCROOT)FMS/mpp/include/mpp_sum_ad.inc $(SRCROOT)FMS/fms/fms_io_unstructured_get_field_size.inc $(SRCROOT)FMS/fms2_io/include/get_checksum.inc $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_axis.inc $(SRCROOT)FMS/drifters/drifters_push.h $(SRCROOT)FMS/mpp/include/mpp_do_checkV.h $(SRCROOT)FMS/mpp/include/group_update_unpack.inc $(SRCROOT)FMS/mpp/include/mpp_write_unlimited_axis.h $(SRCROOT)FMS/mpp/include/mpp_reduce_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_util_nocomm.inc $(SRCROOT)FMS/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)FMS/mpp/include/mpp_type_nocomm.h $(SRCROOT)FMS/drifters/drifters_set_field.h $(SRCROOT)FMS/mpp/include/mpp_define_nest_domains.inc $(SRCROOT)FMS/field_manager/parse.inc $(SRCROOT)FMS/fms2_io/include/get_global_attribute.inc $(SRCROOT)FMS/mpp/include/mpp_do_update_nest.h $(SRCROOT)FMS/fms2_io/include/netcdf_add_restart_variable.inc $(SRCROOT)FMS/mpp/include/mpp_data_nocomm.inc $(SRCROOT)FMS/fms2_io/include/get_variable_attribute.inc $(SRCROOT)FMS/fms/fms_io_unstructured_file_unit.inc $(SRCROOT)FMS/mpp/include/mpp_data_mpi.inc $(SRCROOT)FMS/mosaic/constant.h $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary_ad.h $(SRCROOT)FMS/mpp/include/mpp_scatter.h $(SRCROOT)FMS/mpp/include/mpp_do_global_field_ad.h $(SRCROOT)FMS/fms2_io/include/domain_write.inc $(SRCROOT)FMS/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_read_compressed.h $(SRCROOT)FMS/mpp/include/mpp_do_redistribute.h $(SRCROOT)FMS/mpp/include/mpp_global_field.h $(SRCROOT)FMS/fms2_io/include/netcdf_read_data.inc $(SRCROOT)FMS/mpp/include/mpp_type_mpi.h $(SRCROOT)FMS/fms2_io/include/register_global_attribute.inc $(SRCROOT)FMS/mpp/include/mpp_global_field_ug.h $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_write.inc $(SRCROOT)FMS/mpp/include/mpp_reduce_mpi.h $(SRCROOT)FMS/mpp/include/group_update_pack.inc $(SRCROOT)FMS/fms2_io/include/register_domain_restart_variable.inc $(SRCROOT)FMS/fms2_io/include/unstructured_domain_read.inc $(SRCROOT)FMS/mpp/include/mpp_chksum.h $(SRCROOT)FMS/mpp/include/mpp_domains_misc.inc $(SRCROOT)FMS/fms/fms_io_unstructured_read.inc $(SRCROOT)FMS/mpp/include/mpp_error_s_s.h $(SRCROOT)FMS/mpp/include/mpp_domains_reduce.inc $(SRCROOT)FMS/mpp/include/mpp_domains_define.inc $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_field.inc $(SRCROOT)FMS/fms/fms_io_unstructured_save_restart.inc $(SRCROOT)FMS/fms/fms_io_unstructured_setup_one_field.inc $(SRCROOT)FMS/drifters/drifters_compute_k.h $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_ad.h $(SRCROOT)FMS/mpp/include/mpp_domains_comm.inc $(SRCROOT)FMS/mpp/include/mpp_sum_mpi_ad.h $(SRCROOT)FMS/fms2_io/include/array_utils.inc $(SRCROOT)FMS/fms/fms_io_unstructured_field_exist.inc $(SRCROOT)FMS/mpp/include/mpp_error_s_a.h $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm_ad.h $(SRCROOT)FMS/fms/read_data_2d.inc $(SRCROOT)FMS/mpp/include/mpp_unstruct_domain.inc $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/mosaic/read_mosaic.h $(SRCROOT)FMS/mpp/include/mpp_io_connect.inc $(SRCROOT)FMS/fms/write_data.inc $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_unit.inc $(SRCROOT)FMS/fms/read_data_3d.inc $(SRCROOT)FMS/mpp/include/mpp_do_updateV_nonblock.h $(SRCROOT)FMS/mosaic/interp.h $(SRCROOT)FMS/mpp/include/mpp_error_a_a.h $(SRCROOT)FMS/fms2_io/include/register_unstructured_domain_restart_variable.inc $(SRCROOT)FMS/mpp/include/mpp_write_2Ddecomp.h $(SRCROOT)FMS/mpp/include/mpp_group_update.h $(SRCROOT)FMS/mpp/include/mpp_do_check.h $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary.h $(SRCROOT)FMS/mpp/include/mpp_sum_mpi.h $(SRCROOT)FMS/mpp/include/mpp_chksum_scalar.h $(SRCROOT)FMS/mpp/include/mpp_io_misc.inc $(SRCROOT)FMS/mpp/include/mpp_do_updateV.h $(SRCROOT)FMS/mpp/include/mpp_transmit_mpi.h $(SRCROOT)FMS/mpp/include/mpp_comm.inc $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/mpp/include/mpp_update_nest_domains.h $(SRCROOT)FMS/fms/read_data_4d.inc $(SRCROOT)FMS/mpp/include/mpp_io_read.inc $(SRCROOT)FMS/mpp/include/mpp_write.h $(SRCROOT)FMS/fms2_io/include/domain_read.inc $(SRCROOT)FMS/mpp/include/mpp_chksum_int.h $(SRCROOT)FMS/mpp/include/mpp_do_update_ad.h $(SRCROOT)FMS/mpp/include/mpp_error_a_s.h $(SRCROOT)FMS/mpp/include/mpp_unstruct_pass_data.h $(SRCROOT)FMS/fms2_io/include/unstructured_domain_write.inc $(SRCROOT)FMS/mpp/include/mpp_alltoall_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_global_sum.h $(SRCROOT)FMS/mpp/include/mpp_io_util.inc $(SRCROOT)FMS/mpp/include/mpp_get_boundary_ad.h $(SRCROOT)FMS/mpp/include/mpp_global_field_ad.h $(SRCROOT)FMS/fms2_io/include/compressed_write.inc $(SRCROOT)FMS/fms2_io/include/netcdf_write_data.inc $(SRCROOT)FMS/mpp/include/mpp_io_write.inc $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.h $(SRCROOT)FMS/mpp/include/mpp_do_updateV_ad.h $(SRCROOT)FMS/mosaic/create_xgrid.h $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_name.inc $(SRCROOT)FMS/fms2_io/include/get_data_type_string.inc $(SRCROOT)FMS/mpp/include/mpp_util_mpi.inc $(SRCROOT)FMS/mpp/include/mpp_transmit.inc $(SRCROOT)FMS/mosaic/gradient_c2l.h $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.inc $(SRCROOT)FMS/mpp/include/mpp_read_2Ddecomp.h $(SRCROOT)FMS/fms2_io/include/compressed_read.inc $(SRCROOT)FMS/fms2_io/include/array_utils_char.inc $(SRCROOT)FMS/drifters/fms_switches.h
OBJ = horiz_interp_spherical.o mpp_data.o mosaic2.o diag_axis.o fms2_io.o amip_interp.o axis_utils.o column_diagnostics.o platform.o mosaic_util.o get_grid_version_fms2io.o fms_netcdf_domain_io.o xgrid.o diag_manager.o constants.o diag_util.o mpp_memuse.o horiz_interp_bilinear.o field_manager.o get_grid_version_mpp.o diag_grid.o get_cal_time.o ensemble_manager.o drifters_input.o blackboxio.o random_numbers.o fms_io_utils.o mpp_utilities.o block_control.o diag_output.o time_manager.o tridiagonal.o mersennetwister.o read_mosaic.o drifters_core.o quicksort.o netcdf_io.o diag_integral.o fms.o astronomy.o gradient.o mpp_memutils.o horiz_interp_bicubic.o fms_affinity.o create_xgrid.o drifters_io.o monin_obukhov.o fms_netcdf_unstructured_domain_io.o diag_data.o fm_util.o stock_constants.o fms_io.o drifters.o interp.o gaussian_topog.o gradient_c2l.o mosaic.o horiz_interp.o time_interp_external.o mpp_domains.o data_override.o coupler_types.o sat_vapor_pres.o axis_utils2.o time_interp_external2.o mpp_parameter.o monin_obukhov_inter.o mpp_efp.o mpp.o tracer_manager.o drifters_comm.o sat_vapor_pres_k.o memutils.o diag_table.o topography.o interpolator.o horiz_interp_type.o grid.o mpp_io.o time_interp.o cloud_interpolator.o horiz_interp_conserve.o affinity.o atmos_ocean_fluxes.o
OFF = $(SRCROOT)FMS/mosaic/gradient_c2l.h $(SRCROOT)FMS/drifters/fms_switches.h $(SRCROOT)FMS/fms2_io/include/compressed_read.inc $(SRCROOT)FMS/diag_manager/diag_output.F90 $(SRCROOT)FMS/mpp/include/mpp_util_mpi.inc $(SRCROOT)FMS/diag_manager/diag_data.F90 $(SRCROOT)FMS/mpp/include/mpp_transmit.inc $(SRCROOT)FMS/fms2_io/include/get_data_type_string.inc $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_name.inc $(SRCROOT)FMS/mpp/include/mpp_do_updateV_ad.h $(SRCROOT)FMS/mpp/include/mpp_get_boundary_ad.h $(SRCROOT)FMS/mpp/include/mpp_global_field_ad.h $(SRCROOT)FMS/diag_manager/diag_util.F90 $(SRCROOT)FMS/mpp/include/mpp_unstruct_pass_data.h $(SRCROOT)FMS/data_override/data_override.F90 $(SRCROOT)FMS/mpp/include/mpp_error_a_s.h $(SRCROOT)FMS/mpp/include/mpp_alltoall_nocomm.h $(SRCROOT)FMS/affinity/fms_affinity.F90 $(SRCROOT)FMS/mpp/include/mpp_chksum_int.h $(SRCROOT)FMS/mpp/include/mpp_write.h $(SRCROOT)FMS/fms/read_data_4d.inc $(SRCROOT)FMS/mosaic/read_mosaic.c $(SRCROOT)FMS/tracer_manager/tracer_manager.F90 $(SRCROOT)FMS/mpp/include/mpp_io_misc.inc $(SRCROOT)FMS/mpp/include/mpp_sum_mpi.h $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)FMS/mpp/include/mpp_do_updateV.h $(SRCROOT)FMS/mpp/include/mpp_group_update.h $(SRCROOT)FMS/mpp/include/mpp_write_2Ddecomp.h $(SRCROOT)FMS/mpp/include/mpp_do_check.h $(SRCROOT)FMS/fms/read_data_3d.inc $(SRCROOT)FMS/coupler/ensemble_manager.F90 $(SRCROOT)FMS/drifters/drifters_comm.F90 $(SRCROOT)FMS/fms/write_data.inc $(SRCROOT)FMS/mosaic/mosaic2.F90 $(SRCROOT)FMS/mpp/include/mpp_do_updateV_nonblock.h $(SRCROOT)FMS/random_numbers/mersennetwister.F90 $(SRCROOT)FMS/mpp/include/mpp_io_connect.inc $(SRCROOT)FMS/astronomy/astronomy.F90 $(SRCROOT)FMS/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)FMS/mosaic/mosaic.F90 $(SRCROOT)FMS/drifters/cloud_interpolator.F90 $(SRCROOT)FMS/mpp/include/mpp_sum_mpi_ad.h $(SRCROOT)FMS/mpp/include/mpp_domains_comm.inc $(SRCROOT)FMS/mpp/include/mpp_error_s_a.h $(SRCROOT)FMS/axis_utils/axis_utils2.F90 $(SRCROOT)FMS/drifters/drifters_compute_k.h $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_ad.h $(SRCROOT)FMS/monin_obukhov/monin_obukhov.F90 $(SRCROOT)FMS/fms/fms_io_unstructured_save_restart.inc $(SRCROOT)FMS/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)FMS/topography/topography.F90 $(SRCROOT)FMS/fms2_io/netcdf_io.F90 $(SRCROOT)FMS/fms2_io/include/unstructured_domain_read.inc $(SRCROOT)FMS/fms2_io/include/register_domain_restart_variable.inc $(SRCROOT)FMS/fms/fms_io_unstructured_read.inc $(SRCROOT)FMS/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)FMS/mpp/include/mpp_chksum.h $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_write.inc $(SRCROOT)FMS/amip_interp/amip_interp.F90 $(SRCROOT)FMS/mpp/include/mpp_read_compressed.h $(SRCROOT)FMS/drifters/drifters_input.F90 $(SRCROOT)FMS/fms2_io/include/netcdf_read_data.inc $(SRCROOT)FMS/mpp/include/mpp_scatter.h $(SRCROOT)FMS/mpp/include/mpp_data_mpi.inc $(SRCROOT)FMS/topography/gaussian_topog.F90 $(SRCROOT)FMS/fms2_io/include/get_variable_attribute.inc $(SRCROOT)FMS/fms/fms_io_unstructured_file_unit.inc $(SRCROOT)FMS/mpp/include/mpp_do_update_nest.h $(SRCROOT)FMS/fms2_io/include/get_global_attribute.inc $(SRCROOT)FMS/mpp/include/mpp_util_nocomm.inc $(SRCROOT)FMS/mpp/include/mpp_reduce_nocomm.h $(SRCROOT)FMS/mosaic/create_xgrid.c $(SRCROOT)FMS/mpp/mpp_utilities.F90 $(SRCROOT)FMS/fms2_io/fms_io_utils.F90 $(SRCROOT)FMS/mpp/include/group_update_unpack.inc $(SRCROOT)FMS/mpp/include/mpp_do_checkV.h $(SRCROOT)FMS/drifters/drifters_push.h $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_axis.inc $(SRCROOT)FMS/mpp/include/mpp_sum_ad.inc $(SRCROOT)FMS/mpp/mpp_data.F90 $(SRCROOT)FMS/fms/fms_io_unstructured_get_field_size.inc $(SRCROOT)FMS/time_manager/get_cal_time.F90 $(SRCROOT)FMS/mpp/mpp_domains.F90 $(SRCROOT)FMS/mpp/include/mpp_do_update_nonblock.h $(SRCROOT)FMS/column_diagnostics/column_diagnostics.F90 $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_nonblock.h $(SRCROOT)FMS/mpp/include/mpp_comm_mpi.inc $(SRCROOT)FMS/constants/constants.F90 $(SRCROOT)FMS/mpp/include/mpp_sum.inc $(SRCROOT)FMS/fms2_io/include/register_variable_attribute.inc $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_read.inc $(SRCROOT)FMS/mpp/include/mpp_do_global_field.h $(SRCROOT)FMS/mpp/include/mpp_gather.h $(SRCROOT)FMS/mpp/include/system_clock.h $(SRCROOT)FMS/mpp/include/mpp_get_boundary.h $(SRCROOT)FMS/time_interp/time_interp_external.F90 $(SRCROOT)FMS/mpp/include/mpp_read_2Ddecomp.h $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.inc $(SRCROOT)FMS/data_override/get_grid_version_mpp.F90 $(SRCROOT)FMS/fms2_io/include/array_utils_char.inc $(SRCROOT)FMS/field_manager/fm_util.F90 $(SRCROOT)FMS/mpp/include/mpp_io_write.inc $(SRCROOT)FMS/fms2_io/include/netcdf_write_data.inc $(SRCROOT)FMS/fms2_io/include/compressed_write.inc $(SRCROOT)FMS/mpp/mpp_io.F90 $(SRCROOT)FMS/mosaic/create_xgrid.h $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.h $(SRCROOT)FMS/axis_utils/axis_utils.F90 $(SRCROOT)FMS/mpp/include/mpp_io_util.inc $(SRCROOT)FMS/mosaic/gradient_c2l.c $(SRCROOT)FMS/mosaic/interp.c $(SRCROOT)FMS/horiz_interp/horiz_interp_type.F90 $(SRCROOT)FMS/platform/platform.F90 $(SRCROOT)FMS/mpp/include/mpp_global_sum.h $(SRCROOT)FMS/fms2_io/include/unstructured_domain_write.inc $(SRCROOT)FMS/tridiagonal/tridiagonal.F90 $(SRCROOT)FMS/mpp/include/mpp_io_read.inc $(SRCROOT)FMS/diag_manager/diag_axis.F90 $(SRCROOT)FMS/mpp/include/mpp_do_update_ad.h $(SRCROOT)FMS/fms2_io/include/domain_read.inc $(SRCROOT)FMS/drifters/quicksort.F90 $(SRCROOT)FMS/mpp/include/mpp_comm.inc $(SRCROOT)FMS/time_interp/time_interp_external2.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/mosaic/grid.F90 $(SRCROOT)FMS/mpp/include/mpp_update_nest_domains.h $(SRCROOT)FMS/drifters/drifters_core.F90 $(SRCROOT)FMS/mpp/include/mpp_chksum_scalar.h $(SRCROOT)FMS/memutils/memutils.F90 $(SRCROOT)FMS/mpp/include/mpp_transmit_mpi.h $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary.h $(SRCROOT)FMS/drifters/drifters.F90 $(SRCROOT)FMS/diag_integral/diag_integral.F90 $(SRCROOT)FMS/fms2_io/include/register_unstructured_domain_restart_variable.inc $(SRCROOT)FMS/mpp/include/mpp_error_a_a.h $(SRCROOT)FMS/coupler/atmos_ocean_fluxes.F90 $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_unit.inc $(SRCROOT)FMS/horiz_interp/horiz_interp.F90 $(SRCROOT)FMS/mosaic/interp.h $(SRCROOT)FMS/diag_manager/diag_manager.F90 $(SRCROOT)FMS/time_interp/time_interp.F90 $(SRCROOT)FMS/mosaic/read_mosaic.h $(SRCROOT)FMS/mpp/include/mpp_unstruct_domain.inc $(SRCROOT)FMS/fms/read_data_2d.inc $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm_ad.h $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)FMS/diag_manager/diag_table.F90 $(SRCROOT)FMS/mpp/mpp.F90 $(SRCROOT)FMS/fms/fms_io_unstructured_field_exist.inc $(SRCROOT)FMS/fms2_io/include/array_utils.inc $(SRCROOT)FMS/fms2_io/fms_netcdf_domain_io.F90 $(SRCROOT)FMS/data_override/get_grid_version_fms2io.F90 $(SRCROOT)FMS/interpolator/interpolator.F90 $(SRCROOT)FMS/field_manager/field_manager.F90 $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_field.inc $(SRCROOT)FMS/fms/fms.F90 $(SRCROOT)FMS/fms/fms_io.F90 $(SRCROOT)FMS/fms/fms_io_unstructured_setup_one_field.inc $(SRCROOT)FMS/mpp/include/mpp_error_s_s.h $(SRCROOT)FMS/mpp/include/mpp_domains_define.inc $(SRCROOT)FMS/mosaic/gradient.F90 $(SRCROOT)FMS/mpp/include/mpp_domains_reduce.inc $(SRCROOT)FMS/fms2_io/fms2_io.F90 $(SRCROOT)FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90 $(SRCROOT)FMS/drifters/drifters_io.F90 $(SRCROOT)FMS/mpp/include/mpp_domains_misc.inc $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm.h $(SRCROOT)FMS/mosaic/mosaic_util.c $(SRCROOT)FMS/mpp/include/group_update_pack.inc $(SRCROOT)FMS/mpp/include/mpp_reduce_mpi.h $(SRCROOT)FMS/random_numbers/random_numbers.F90 $(SRCROOT)FMS/mpp/include/mpp_do_redistribute.h $(SRCROOT)FMS/monin_obukhov/monin_obukhov_inter.F90 $(SRCROOT)FMS/mpp/include/mpp_global_field.h $(SRCROOT)FMS/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_global_field_ug.h $(SRCROOT)FMS/coupler/coupler_types.F90 $(SRCROOT)FMS/mpp/include/mpp_type_mpi.h $(SRCROOT)FMS/fms2_io/include/register_global_attribute.inc $(SRCROOT)FMS/mpp/mpp_memutils.F90 $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary_ad.h $(SRCROOT)FMS/mosaic/constant.h $(SRCROOT)FMS/mpp/include/mpp_do_global_field_ad.h $(SRCROOT)FMS/fms2_io/include/domain_write.inc $(SRCROOT)FMS/mpp/include/mpp_data_nocomm.inc $(SRCROOT)FMS/exchange/stock_constants.F90 $(SRCROOT)FMS/field_manager/parse.inc $(SRCROOT)FMS/mpp/mpp_parameter.F90 $(SRCROOT)FMS/fms2_io/include/netcdf_add_restart_variable.inc $(SRCROOT)FMS/block_control/block_control.F90 $(SRCROOT)FMS/mpp/include/mpp_type_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)FMS/drifters/drifters_set_field.h $(SRCROOT)FMS/mpp/include/mpp_define_nest_domains.inc $(SRCROOT)FMS/mpp/include/mpp_write_unlimited_axis.h $(SRCROOT)FMS/fms2_io/include/get_checksum.inc $(SRCROOT)FMS/diag_manager/diag_grid.F90 $(SRCROOT)FMS/exchange/xgrid.F90 $(SRCROOT)FMS/time_manager/time_manager.F90 $(SRCROOT)FMS/fms2_io/include/compute_global_checksum.inc $(SRCROOT)FMS/mpp/include/mpp_domains_util.inc $(SRCROOT)FMS/affinity/affinity.c $(SRCROOT)FMS/fms2_io/blackboxio.F90 $(SRCROOT)FMS/mpp/include/mpp_alltoall_mpi.h $(SRCROOT)FMS/mpp/include/mpp_do_update.h $(SRCROOT)FMS/mpp/include/mpp_update_domains2D.h $(SRCROOT)FMS/mpp/mpp_memuse.c $(SRCROOT)FMS/mpp/include/mpp_util.inc $(SRCROOT)FMS/mpp/mpp_efp.F90 $(SRCROOT)FMS/mpp/include/mpp_global_sum_tl.h $(SRCROOT)FMS/mosaic/mosaic_util.h $(SRCROOT)FMS/mpp/include/mpp_global_reduce.h $(SRCROOT)FMS/mpp/include/mpp_write_compressed.h $(SRCROOT)FMS/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)FMS/mpp/include/mpp_global_sum_ad.h
clean: neat
	-rm -f .libfms.a.cppdefs $(OBJ) libfms.a

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libfms.a: $(OBJ)
	$(AR) $(ARFLAGS) libfms.a $(OBJ)


