# Makefile created by mkmf 19.3.1

CPPDEFS = -DINTERNAL_FILE_NML -Duse_libMPI -Duse_netCDF

TMPFILES = *genmod.f90

include $(MK_TEMPLATE)

.DEFAULT:
	-echo $@ does not exist.
all: libfms.a
MersenneTwister.o: $(SRCROOT)FMS/random_numbers/MersenneTwister.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/random_numbers/MersenneTwister.F90
affinity.o: $(SRCROOT)FMS/affinity/affinity.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)FMS/affinity/affinity.c
amip_interp.o: $(SRCROOT)FMS/amip_interp/amip_interp.F90 $(SRCROOT)FMS/include/file_version.h time_interp.o time_manager.o get_cal_time.o horiz_interp.o fms.o constants.o platform.o mpp.o fms2_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/amip_interp/amip_interp.F90
astronomy.o: $(SRCROOT)FMS/astronomy/astronomy.F90 $(SRCROOT)FMS/include/file_version.h fms.o time_manager.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/astronomy/astronomy.F90
atmos_ocean_fluxes.o: $(SRCROOT)FMS/coupler/atmos_ocean_fluxes.F90 $(SRCROOT)FMS/include/file_version.h mpp.o fms.o coupler_types.o field_manager.o fm_util.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/coupler/atmos_ocean_fluxes.F90
axis_utils.o: $(SRCROOT)FMS/axis_utils/axis_utils.F90 $(SRCROOT)FMS/include/file_version.h mpp_io.o mpp.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/axis_utils/axis_utils.F90
axis_utils2.o: $(SRCROOT)FMS/axis_utils/axis_utils2.F90 $(SRCROOT)FMS/include/file_version.h mpp.o fms.o fms2_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/axis_utils/axis_utils2.F90
blackboxio.o: $(SRCROOT)FMS/fms2_io/blackboxio.F90 mpp_domains.o fms_io_utils.o netcdf_io.o fms_netcdf_domain_io.o fms_netcdf_unstructured_domain_io.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/fms2_io/blackboxio.F90
block_control.o: $(SRCROOT)FMS/block_control/block_control.F90 $(SRCROOT)FMS/include/fms_platform.h mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/block_control/block_control.F90
cloud_interpolator.o: $(SRCROOT)FMS/drifters/cloud_interpolator.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/drifters/cloud_interpolator.F90
column_diagnostics.o: $(SRCROOT)FMS/column_diagnostics/column_diagnostics.F90 $(SRCROOT)FMS/include/file_version.h mpp_io.o fms.o time_manager.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/column_diagnostics/column_diagnostics.F90
constants.o: $(SRCROOT)FMS/constants/constants.F90 $(SRCROOT)FMS/include/file_version.h platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/constants/constants.F90
coupler_types.o: $(SRCROOT)FMS/coupler/coupler_types.F90 $(SRCROOT)FMS/include/file_version.h fms.o fms_io.o time_manager.o diag_manager.o data_override.o mpp_domains.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/coupler/coupler_types.F90
create_xgrid.o: $(SRCROOT)FMS/mosaic/create_xgrid.c $(SRCROOT)FMS/mosaic/mosaic_util.h $(SRCROOT)FMS/mosaic/create_xgrid.h $(SRCROOT)FMS/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)FMS/mosaic	$(SRCROOT)FMS/mosaic/create_xgrid.c
data_override.o: $(SRCROOT)FMS/data_override/data_override.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h platform.o constants.o mpp.o horiz_interp.o time_interp_external2.o fms_io.o fms.o axis_utils2.o mpp_domains.o time_manager.o fms2_io.o mosaic2.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/data_override/data_override.F90
diag_axis.o: $(SRCROOT)FMS/diag_manager/diag_axis.F90 $(SRCROOT)FMS/include/file_version.h platform.o mpp_domains.o fms.o diag_data.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/diag_manager/diag_axis.F90
diag_data.o: $(SRCROOT)FMS/diag_manager/diag_data.F90 $(SRCROOT)FMS/include/file_version.h platform.o time_manager.o mpp_domains.o mpp_io.o fms.o fms2_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/diag_manager/diag_data.F90
diag_grid.o: $(SRCROOT)FMS/diag_manager/diag_grid.F90 $(SRCROOT)FMS/include/file_version.h platform.o constants.o fms.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/diag_manager/diag_grid.F90
diag_integral.o: $(SRCROOT)FMS/diag_integral/diag_integral.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h time_manager.o mpp.o fms.o fms2_io.o constants.o ensemble_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/diag_integral/diag_integral.F90
diag_manager.o: $(SRCROOT)FMS/diag_manager/diag_manager.F90 $(SRCROOT)FMS/include/file_version.h platform.o time_manager.o mpp_io.o mpp.o fms.o fms_io.o diag_axis.o diag_util.o diag_data.o diag_table.o diag_output.o diag_grid.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/diag_manager/diag_manager.F90
diag_output.o: $(SRCROOT)FMS/diag_manager/diag_output.F90 platform.o mpp_io.o mpp_domains.o mpp.o diag_axis.o diag_data.o time_manager.o fms.o fms2_io.o axis_utils2.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/diag_manager/diag_output.F90
diag_table.o: $(SRCROOT)FMS/diag_manager/diag_table.F90 mpp_io.o mpp.o fms.o time_manager.o constants.o diag_data.o diag_util.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/diag_manager/diag_table.F90
diag_util.o: $(SRCROOT)FMS/diag_manager/diag_util.F90 $(SRCROOT)FMS/include/file_version.h platform.o diag_data.o diag_axis.o diag_output.o diag_grid.o fms.o fms_io.o mpp_domains.o time_manager.o mpp_io.o mpp.o constants.o fms2_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/diag_manager/diag_util.F90
drifters.o: $(SRCROOT)FMS/drifters/drifters.F90 $(SRCROOT)FMS/drifters/fms_switches.h $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/drifters/drifters_push.h $(SRCROOT)FMS/drifters/drifters_set_field.h $(SRCROOT)FMS/drifters/drifters_compute_k.h mpp.o mpp_domains.o drifters_core.o drifters_input.o drifters_io.o drifters_comm.o cloud_interpolator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/drifters -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/drifters/drifters.F90
drifters_comm.o: $(SRCROOT)FMS/drifters/drifters_comm.F90 $(SRCROOT)FMS/drifters/fms_switches.h $(SRCROOT)FMS/include/fms_platform.h mpp.o mpp_domains.o drifters_core.o drifters_input.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/drifters -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/drifters/drifters_comm.F90
drifters_core.o: $(SRCROOT)FMS/drifters/drifters_core.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/drifters/drifters_core.F90
drifters_input.o: $(SRCROOT)FMS/drifters/drifters_input.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h
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
fms_io_utils.o: $(SRCROOT)FMS/fms2_io/fms_io_utils.F90 $(SRCROOT)FMS/fms2_io/include/array_utils.inc $(SRCROOT)FMS/fms2_io/include/array_utils_char.inc $(SRCROOT)FMS/fms2_io/include/get_data_type_string.inc $(SRCROOT)FMS/fms2_io/include/get_checksum.inc mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/fms2_io/include	$(SRCROOT)FMS/fms2_io/fms_io_utils.F90
fms_netcdf_domain_io.o: $(SRCROOT)FMS/fms2_io/fms_netcdf_domain_io.F90 $(SRCROOT)FMS/fms2_io/include/register_domain_restart_variable.inc $(SRCROOT)FMS/fms2_io/include/domain_read.inc $(SRCROOT)FMS/fms2_io/include/domain_write.inc $(SRCROOT)FMS/fms2_io/include/compute_global_checksum.inc mpp.o mpp_domains.o fms_io_utils.o netcdf_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/fms2_io/include	$(SRCROOT)FMS/fms2_io/fms_netcdf_domain_io.F90
fms_netcdf_unstructured_domain_io.o: $(SRCROOT)FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90 $(SRCROOT)FMS/fms2_io/include/register_unstructured_domain_restart_variable.inc $(SRCROOT)FMS/fms2_io/include/unstructured_domain_read.inc $(SRCROOT)FMS/fms2_io/include/unstructured_domain_write.inc mpp_domains.o fms_io_utils.o netcdf_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/fms2_io/include	$(SRCROOT)FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90
gaussian_topog.o: $(SRCROOT)FMS/topography/gaussian_topog.F90 $(SRCROOT)FMS/include/file_version.h fms.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/topography/gaussian_topog.F90
get_cal_time.o: $(SRCROOT)FMS/time_manager/get_cal_time.F90 $(SRCROOT)FMS/include/file_version.h fms.o time_manager.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/time_manager/get_cal_time.F90
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
interpolator.o: $(SRCROOT)FMS/interpolator/interpolator.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h mpp.o mpp_domains.o diag_manager.o fms.o fms2_io.o horiz_interp.o time_manager.o time_interp.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/interpolator/interpolator.F90
memutils.o: $(SRCROOT)FMS/memutils/memutils.F90 mpp.o mpp_memutils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/memutils/memutils.F90
monin_obukhov.o: $(SRCROOT)FMS/monin_obukhov/monin_obukhov.F90 $(SRCROOT)FMS/include/file_version.h constants.o mpp.o fms.o monin_obukhov_kernel.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/monin_obukhov/monin_obukhov.F90
monin_obukhov_kernel.o: $(SRCROOT)FMS/monin_obukhov/monin_obukhov_kernel.F90 $(SRCROOT)FMS/include/fms_platform.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/monin_obukhov/monin_obukhov_kernel.F90
mosaic.o: $(SRCROOT)FMS/mosaic/mosaic.F90 $(SRCROOT)FMS/include/file_version.h fms.o mpp.o mpp_io.o fms_io.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/mosaic/mosaic.F90
mosaic2.o: $(SRCROOT)FMS/mosaic/mosaic2.F90 $(SRCROOT)FMS/include/file_version.h fms.o mpp.o mpp_domains.o constants.o fms2_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/mosaic/mosaic2.F90
mosaic_util.o: $(SRCROOT)FMS/mosaic/mosaic_util.c $(SRCROOT)FMS/mosaic/mosaic_util.h $(SRCROOT)FMS/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)FMS/mosaic	$(SRCROOT)FMS/mosaic/mosaic_util.c
mpp.o: $(SRCROOT)FMS/mpp/mpp.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/mpp/include/system_clock.h $(SRCROOT)FMS/mpp/include/mpp_util.inc $(SRCROOT)FMS/mpp/include/mpp_util_sma.inc $(SRCROOT)FMS/mpp/include/mpp_util_mpi.inc $(SRCROOT)FMS/mpp/include/mpp_util_nocomm.inc $(SRCROOT)FMS/mpp/include/mpp_error_a_a.h $(SRCROOT)FMS/mpp/include/mpp_error_a_s.h $(SRCROOT)FMS/mpp/include/mpp_error_s_a.h $(SRCROOT)FMS/mpp/include/mpp_error_s_s.h $(SRCROOT)FMS/mpp/include/mpp_comm.inc $(SRCROOT)FMS/mpp/include/mpp_comm_sma.inc $(SRCROOT)FMS/mpp/include/mpp_transmit_sma.h $(SRCROOT)FMS/mpp/include/mpp_transmit.inc $(SRCROOT)FMS/mpp/include/mpp_reduce_sma.h $(SRCROOT)FMS/mpp/include/mpp_sum_sma.h $(SRCROOT)FMS/mpp/include/mpp_sum.inc $(SRCROOT)FMS/mpp/include/mpp_alltoall_sma.h $(SRCROOT)FMS/mpp/include/mpp_type_sma.h $(SRCROOT)FMS/mpp/include/mpp_comm_mpi.inc $(SRCROOT)FMS/mpp/include/mpp_transmit_mpi.h $(SRCROOT)FMS/mpp/include/mpp_reduce_mpi.h $(SRCROOT)FMS/mpp/include/mpp_sum_mpi.h $(SRCROOT)FMS/mpp/include/mpp_sum_mpi_ad.h $(SRCROOT)FMS/mpp/include/mpp_sum_ad.inc $(SRCROOT)FMS/mpp/include/mpp_alltoall_mpi.h $(SRCROOT)FMS/mpp/include/mpp_type_mpi.h $(SRCROOT)FMS/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)FMS/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_reduce_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm_ad.h $(SRCROOT)FMS/mpp/include/mpp_alltoall_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_type_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_chksum_int.h $(SRCROOT)FMS/mpp/include/mpp_chksum_scalar.h $(SRCROOT)FMS/mpp/include/mpp_chksum.h $(SRCROOT)FMS/mpp/include/mpp_gather.h $(SRCROOT)FMS/mpp/include/mpp_scatter.h mpp_parameter.o mpp_data.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include -I$(SRCROOT)FMS/mpp/include	$(SRCROOT)FMS/mpp/mpp.F90
mpp_data.o: $(SRCROOT)FMS/mpp/mpp_data.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/mpp/include/mpp_data_sma.inc $(SRCROOT)FMS/mpp/include/mpp_data_mpi.inc $(SRCROOT)FMS/mpp/include/mpp_data_nocomm.inc mpp_parameter.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include -I$(SRCROOT)FMS/mpp/include	$(SRCROOT)FMS/mpp/mpp_data.F90
mpp_domains.o: $(SRCROOT)FMS/mpp/mpp_domains.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/mpp/include/mpp_define_nest_domains.inc $(SRCROOT)FMS/mpp/include/mpp_domains_util.inc $(SRCROOT)FMS/mpp/include/mpp_domains_comm.inc $(SRCROOT)FMS/mpp/include/mpp_domains_define.inc $(SRCROOT)FMS/mpp/include/mpp_domains_misc.inc $(SRCROOT)FMS/mpp/include/mpp_update_domains2D.h $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_nonblock.h $(SRCROOT)FMS/mpp/include/mpp_do_update_nonblock.h $(SRCROOT)FMS/mpp/include/mpp_do_updateV_nonblock.h $(SRCROOT)FMS/mpp/include/mpp_do_update.h $(SRCROOT)FMS/mpp/include/mpp_do_updateV.h $(SRCROOT)FMS/mpp/include/mpp_do_check.h $(SRCROOT)FMS/mpp/include/mpp_do_checkV.h $(SRCROOT)FMS/mpp/include/mpp_update_nest_domains.h $(SRCROOT)FMS/mpp/include/mpp_do_update_nest.h $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_ad.h $(SRCROOT)FMS/mpp/include/mpp_do_update_ad.h $(SRCROOT)FMS/mpp/include/mpp_do_updateV_ad.h $(SRCROOT)FMS/mpp/include/mpp_do_redistribute.h $(SRCROOT)FMS/mpp/include/mpp_get_boundary.h $(SRCROOT)FMS/mpp/include/mpp_get_boundary_ad.h $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary.h $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary_ad.h $(SRCROOT)FMS/mpp/include/mpp_group_update.h $(SRCROOT)FMS/mpp/include/group_update_pack.inc $(SRCROOT)FMS/mpp/include/group_update_unpack.inc $(SRCROOT)FMS/mpp/include/mpp_domains_reduce.inc $(SRCROOT)FMS/mpp/include/mpp_global_reduce.h $(SRCROOT)FMS/mpp/include/mpp_global_sum.h $(SRCROOT)FMS/mpp/include/mpp_global_sum_tl.h $(SRCROOT)FMS/mpp/include/mpp_global_sum_ad.h $(SRCROOT)FMS/mpp/include/mpp_global_field.h $(SRCROOT)FMS/mpp/include/mpp_global_field_ad.h $(SRCROOT)FMS/mpp/include/mpp_do_global_field.h $(SRCROOT)FMS/mpp/include/mpp_do_global_field_ad.h $(SRCROOT)FMS/mpp/include/mpp_unstruct_domain.inc $(SRCROOT)FMS/mpp/include/mpp_unstruct_pass_data.h $(SRCROOT)FMS/mpp/include/mpp_global_field_ug.h mpp_parameter.o mpp_data.o mpp.o mpp_memutils.o mpp_pset.o mpp_efp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include -I$(SRCROOT)FMS/mpp/include	$(SRCROOT)FMS/mpp/mpp_domains.F90
mpp_efp.o: $(SRCROOT)FMS/mpp/mpp_efp.F90 $(SRCROOT)FMS/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/mpp/mpp_efp.F90
mpp_io.o: $(SRCROOT)FMS/mpp/mpp_io.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/mpp/include/mpp_io_util.inc $(SRCROOT)FMS/mpp/include/mpp_io_misc.inc $(SRCROOT)FMS/mpp/include/mpp_io_connect.inc $(SRCROOT)FMS/mpp/include/mpp_io_read.inc $(SRCROOT)FMS/mpp/include/mpp_read_2Ddecomp.h $(SRCROOT)FMS/mpp/include/mpp_read_compressed.h $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.inc $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.h $(SRCROOT)FMS/mpp/include/mpp_io_write.inc $(SRCROOT)FMS/mpp/include/mpp_write_2Ddecomp.h $(SRCROOT)FMS/mpp/include/mpp_write_compressed.h $(SRCROOT)FMS/mpp/include/mpp_write_unlimited_axis.h $(SRCROOT)FMS/mpp/include/mpp_write.h $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_write.inc $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_read.inc mpp_parameter.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include -I$(SRCROOT)FMS/mpp/include	$(SRCROOT)FMS/mpp/mpp_io.F90
mpp_memuse.o: $(SRCROOT)FMS/mpp/mpp_memuse.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)FMS/mpp/mpp_memuse.c
mpp_memutils.o: $(SRCROOT)FMS/mpp/mpp_memutils.F90 $(SRCROOT)FMS/fms/../include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/fms/../include	$(SRCROOT)FMS/mpp/mpp_memutils.F90
mpp_parameter.o: $(SRCROOT)FMS/mpp/mpp_parameter.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/mpp/mpp_parameter.F90
mpp_pset.o: $(SRCROOT)FMS/mpp/mpp_pset.F90 $(SRCROOT)FMS/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/mpp/mpp_pset.F90
mpp_utilities.o: $(SRCROOT)FMS/mpp/mpp_utilities.F90 $(SRCROOT)FMS/include/file_version.h mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/mpp/mpp_utilities.F90
netcdf_io.o: $(SRCROOT)FMS/fms2_io/netcdf_io.F90 $(SRCROOT)FMS/fms2_io/include/netcdf_add_restart_variable.inc $(SRCROOT)FMS/fms2_io/include/netcdf_read_data.inc $(SRCROOT)FMS/fms2_io/include/netcdf_write_data.inc $(SRCROOT)FMS/fms2_io/include/register_global_attribute.inc $(SRCROOT)FMS/fms2_io/include/register_variable_attribute.inc $(SRCROOT)FMS/fms2_io/include/get_global_attribute.inc $(SRCROOT)FMS/fms2_io/include/get_variable_attribute.inc $(SRCROOT)FMS/fms2_io/include/compressed_write.inc $(SRCROOT)FMS/fms2_io/include/compressed_read.inc mpp.o fms_io_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/fms2_io/include	$(SRCROOT)FMS/fms2_io/netcdf_io.F90
nsclock.o: $(SRCROOT)FMS/mpp/nsclock.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)FMS/mpp/nsclock.c
platform.o: $(SRCROOT)FMS/platform/platform.F90 $(SRCROOT)FMS/include/fms_platform.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/platform/platform.F90
quicksort.o: $(SRCROOT)FMS/drifters/quicksort.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/drifters/quicksort.F90
random_numbers.o: $(SRCROOT)FMS/random_numbers/random_numbers.F90 MersenneTwister.o time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/random_numbers/random_numbers.F90
read_mosaic.o: $(SRCROOT)FMS/mosaic/read_mosaic.c $(SRCROOT)FMS/mosaic/read_mosaic.h $(SRCROOT)FMS/mosaic/constant.h $(SRCROOT)FMS/mosaic/mosaic_util.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)FMS/mosaic	$(SRCROOT)FMS/mosaic/read_mosaic.c
sat_vapor_pres.o: $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)FMS/include/file_version.h constants.o fms.o mpp.o sat_vapor_pres_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres.F90
sat_vapor_pres_k.o: $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)FMS/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres_k.F90
stock_constants.o: $(SRCROOT)FMS/exchange/stock_constants.F90 $(SRCROOT)FMS/include/file_version.h mpp.o fms.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/exchange/stock_constants.F90
threadloc.o: $(SRCROOT)FMS/mpp/threadloc.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)FMS/mpp/threadloc.c
time_interp.o: $(SRCROOT)FMS/time_interp/time_interp.F90 $(SRCROOT)FMS/include/file_version.h time_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/time_interp/time_interp.F90
time_interp_external.o: $(SRCROOT)FMS/time_interp/time_interp_external.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h fms.o mpp.o mpp_io.o time_manager.o get_cal_time.o mpp_domains.o time_interp.o axis_utils.o platform.o horiz_interp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/time_interp/time_interp_external.F90
time_interp_external2.o: $(SRCROOT)FMS/time_interp/time_interp_external2.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h fms.o mpp.o time_manager.o get_cal_time.o mpp_domains.o time_interp.o axis_utils2.o platform.o horiz_interp.o fms2_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/time_interp/time_interp_external2.F90
time_manager.o: $(SRCROOT)FMS/time_manager/time_manager.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h constants.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/time_manager/time_manager.F90
topography.o: $(SRCROOT)FMS/topography/topography.F90 $(SRCROOT)FMS/include/file_version.h gaussian_topog.o horiz_interp.o fms.o fms2_io.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/topography/topography.F90
tracer_manager.o: $(SRCROOT)FMS/tracer_manager/tracer_manager.F90 $(SRCROOT)FMS/include/file_version.h mpp.o fms.o field_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/tracer_manager/tracer_manager.F90
tridiagonal.o: $(SRCROOT)FMS/tridiagonal/tridiagonal.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMS/tridiagonal/tridiagonal.F90
xgrid.o: $(SRCROOT)FMS/exchange/xgrid.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h fms.o mpp.o mpp_domains.o constants.o mosaic2.o stock_constants.o gradient.o fms2_io.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)FMS/exchange/xgrid.F90
./create_xgrid.c: $(SRCROOT)FMS/mosaic/create_xgrid.c
	cp $(SRCROOT)FMS/mosaic/create_xgrid.c .
./mpp_util_nocomm.inc: $(SRCROOT)FMS/mpp/include/mpp_util_nocomm.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_util_nocomm.inc .
./mpp_global_field_ug.h: $(SRCROOT)FMS/mpp/include/mpp_global_field_ug.h
	cp $(SRCROOT)FMS/mpp/include/mpp_global_field_ug.h .
./tracer_manager.F90: $(SRCROOT)FMS/tracer_manager/tracer_manager.F90
	cp $(SRCROOT)FMS/tracer_manager/tracer_manager.F90 .
./drifters.F90: $(SRCROOT)FMS/drifters/drifters.F90
	cp $(SRCROOT)FMS/drifters/drifters.F90 .
./mpp_define_nest_domains.inc: $(SRCROOT)FMS/mpp/include/mpp_define_nest_domains.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_define_nest_domains.inc .
./drifters_core.F90: $(SRCROOT)FMS/drifters/drifters_core.F90
	cp $(SRCROOT)FMS/drifters/drifters_core.F90 .
./parse.inc: $(SRCROOT)FMS/field_manager/parse.inc
	cp $(SRCROOT)FMS/field_manager/parse.inc .
./mpp_update_domains2D.h: $(SRCROOT)FMS/mpp/include/mpp_update_domains2D.h
	cp $(SRCROOT)FMS/mpp/include/mpp_update_domains2D.h .
./fms_io_unstructured_get_field_size.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_get_field_size.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_get_field_size.inc .
./read_mosaic.h: $(SRCROOT)FMS/mosaic/read_mosaic.h
	cp $(SRCROOT)FMS/mosaic/read_mosaic.h .
./mpp_io_connect.inc: $(SRCROOT)FMS/mpp/include/mpp_io_connect.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_io_connect.inc .
./diag_integral.F90: $(SRCROOT)FMS/diag_integral/diag_integral.F90
	cp $(SRCROOT)FMS/diag_integral/diag_integral.F90 .
./get_global_attribute.inc: $(SRCROOT)FMS/fms2_io/include/get_global_attribute.inc
	cp $(SRCROOT)FMS/fms2_io/include/get_global_attribute.inc .
./cloud_interpolator.F90: $(SRCROOT)FMS/drifters/cloud_interpolator.F90
	cp $(SRCROOT)FMS/drifters/cloud_interpolator.F90 .
./mpp_do_global_field.h: $(SRCROOT)FMS/mpp/include/mpp_do_global_field.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_global_field.h .
./fms_platform.h: $(SRCROOT)FMS/fms/../include/fms_platform.h
	cp $(SRCROOT)FMS/fms/../include/fms_platform.h .
./mpp_chksum.h: $(SRCROOT)FMS/mpp/include/mpp_chksum.h
	cp $(SRCROOT)FMS/mpp/include/mpp_chksum.h .
./netcdf_read_data.inc: $(SRCROOT)FMS/fms2_io/include/netcdf_read_data.inc
	cp $(SRCROOT)FMS/fms2_io/include/netcdf_read_data.inc .
./platform.F90: $(SRCROOT)FMS/platform/platform.F90
	cp $(SRCROOT)FMS/platform/platform.F90 .
./register_unstructured_domain_restart_variable.inc: $(SRCROOT)FMS/fms2_io/include/register_unstructured_domain_restart_variable.inc
	cp $(SRCROOT)FMS/fms2_io/include/register_unstructured_domain_restart_variable.inc .
./fms_platform.h: $(SRCROOT)FMS/include/fms_platform.h
	cp $(SRCROOT)FMS/include/fms_platform.h .
./fms_io.F90: $(SRCROOT)FMS/fms/fms_io.F90
	cp $(SRCROOT)FMS/fms/fms_io.F90 .
./mpp_memuse.c: $(SRCROOT)FMS/mpp/mpp_memuse.c
	cp $(SRCROOT)FMS/mpp/mpp_memuse.c .
./interp.h: $(SRCROOT)FMS/mosaic/interp.h
	cp $(SRCROOT)FMS/mosaic/interp.h .
./mpp_global_sum_ad.h: $(SRCROOT)FMS/mpp/include/mpp_global_sum_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_global_sum_ad.h .
./mosaic.F90: $(SRCROOT)FMS/mosaic/mosaic.F90
	cp $(SRCROOT)FMS/mosaic/mosaic.F90 .
./mpp_do_check.h: $(SRCROOT)FMS/mpp/include/mpp_do_check.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_check.h .
./file_version.h: $(SRCROOT)FMS/include/file_version.h
	cp $(SRCROOT)FMS/include/file_version.h .
./mpp_do_update.h: $(SRCROOT)FMS/mpp/include/mpp_do_update.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_update.h .
./fms_io_unstructured_save_restart.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_save_restart.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_save_restart.inc .
./mpp_reduce_nocomm.h: $(SRCROOT)FMS/mpp/include/mpp_reduce_nocomm.h
	cp $(SRCROOT)FMS/mpp/include/mpp_reduce_nocomm.h .
./mpp_global_sum.h: $(SRCROOT)FMS/mpp/include/mpp_global_sum.h
	cp $(SRCROOT)FMS/mpp/include/mpp_global_sum.h .
./fms_netcdf_domain_io.F90: $(SRCROOT)FMS/fms2_io/fms_netcdf_domain_io.F90
	cp $(SRCROOT)FMS/fms2_io/fms_netcdf_domain_io.F90 .
./mpp_do_get_boundary_ad.h: $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary_ad.h .
./amip_interp.F90: $(SRCROOT)FMS/amip_interp/amip_interp.F90
	cp $(SRCROOT)FMS/amip_interp/amip_interp.F90 .
./diag_table.F90: $(SRCROOT)FMS/diag_manager/diag_table.F90
	cp $(SRCROOT)FMS/diag_manager/diag_table.F90 .
./nsclock.c: $(SRCROOT)FMS/mpp/nsclock.c
	cp $(SRCROOT)FMS/mpp/nsclock.c .
./unstructured_domain_write.inc: $(SRCROOT)FMS/fms2_io/include/unstructured_domain_write.inc
	cp $(SRCROOT)FMS/fms2_io/include/unstructured_domain_write.inc .
./mpp_data_mpi.inc: $(SRCROOT)FMS/mpp/include/mpp_data_mpi.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_data_mpi.inc .
./tridiagonal.F90: $(SRCROOT)FMS/tridiagonal/tridiagonal.F90
	cp $(SRCROOT)FMS/tridiagonal/tridiagonal.F90 .
./mpp_global_sum_tl.h: $(SRCROOT)FMS/mpp/include/mpp_global_sum_tl.h
	cp $(SRCROOT)FMS/mpp/include/mpp_global_sum_tl.h .
./mpp_read_compressed.h: $(SRCROOT)FMS/mpp/include/mpp_read_compressed.h
	cp $(SRCROOT)FMS/mpp/include/mpp_read_compressed.h .
./system_clock.h: $(SRCROOT)FMS/mpp/include/system_clock.h
	cp $(SRCROOT)FMS/mpp/include/system_clock.h .
./horiz_interp_conserve.F90: $(SRCROOT)FMS/horiz_interp/horiz_interp_conserve.F90
	cp $(SRCROOT)FMS/horiz_interp/horiz_interp_conserve.F90 .
./constants.F90: $(SRCROOT)FMS/constants/constants.F90
	cp $(SRCROOT)FMS/constants/constants.F90 .
./mpp_transmit.inc: $(SRCROOT)FMS/mpp/include/mpp_transmit.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_transmit.inc .
./mpp_global_field_ad.h: $(SRCROOT)FMS/mpp/include/mpp_global_field_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_global_field_ad.h .
./mpp_util_sma.inc: $(SRCROOT)FMS/mpp/include/mpp_util_sma.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_util_sma.inc .
./coupler_types.F90: $(SRCROOT)FMS/coupler/coupler_types.F90
	cp $(SRCROOT)FMS/coupler/coupler_types.F90 .
./mpp_io_unstructured_read.inc: $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_read.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_read.inc .
./blackboxio.F90: $(SRCROOT)FMS/fms2_io/blackboxio.F90
	cp $(SRCROOT)FMS/fms2_io/blackboxio.F90 .
./diag_data.F90: $(SRCROOT)FMS/diag_manager/diag_data.F90
	cp $(SRCROOT)FMS/diag_manager/diag_data.F90 .
./fms_io_unstructured_get_file_unit.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_unit.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_unit.inc .
./get_cal_time.F90: $(SRCROOT)FMS/time_manager/get_cal_time.F90
	cp $(SRCROOT)FMS/time_manager/get_cal_time.F90 .
./time_interp_external2.F90: $(SRCROOT)FMS/time_interp/time_interp_external2.F90
	cp $(SRCROOT)FMS/time_interp/time_interp_external2.F90 .
./mpp_do_updateV.h: $(SRCROOT)FMS/mpp/include/mpp_do_updateV.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_updateV.h .
./mpp_do_global_field_ad.h: $(SRCROOT)FMS/mpp/include/mpp_do_global_field_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_global_field_ad.h .
./mpp_data_sma.inc: $(SRCROOT)FMS/mpp/include/mpp_data_sma.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_data_sma.inc .
./fms_switches.h: $(SRCROOT)FMS/drifters/fms_switches.h
	cp $(SRCROOT)FMS/drifters/fms_switches.h .
./create_xgrid.h: $(SRCROOT)FMS/mosaic/create_xgrid.h
	cp $(SRCROOT)FMS/mosaic/create_xgrid.h .
./mpp_memutils.F90: $(SRCROOT)FMS/mpp/mpp_memutils.F90
	cp $(SRCROOT)FMS/mpp/mpp_memutils.F90 .
./ensemble_manager.F90: $(SRCROOT)FMS/coupler/ensemble_manager.F90
	cp $(SRCROOT)FMS/coupler/ensemble_manager.F90 .
./mpp_data.F90: $(SRCROOT)FMS/mpp/mpp_data.F90
	cp $(SRCROOT)FMS/mpp/mpp_data.F90 .
./mpp_io.F90: $(SRCROOT)FMS/mpp/mpp_io.F90
	cp $(SRCROOT)FMS/mpp/mpp_io.F90 .
./stock_constants.F90: $(SRCROOT)FMS/exchange/stock_constants.F90
	cp $(SRCROOT)FMS/exchange/stock_constants.F90 .
./register_global_attribute.inc: $(SRCROOT)FMS/fms2_io/include/register_global_attribute.inc
	cp $(SRCROOT)FMS/fms2_io/include/register_global_attribute.inc .
./fms.F90: $(SRCROOT)FMS/fms/fms.F90
	cp $(SRCROOT)FMS/fms/fms.F90 .
./mpp_domains_define.inc: $(SRCROOT)FMS/mpp/include/mpp_domains_define.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_domains_define.inc .
./mpp_reduce_mpi.h: $(SRCROOT)FMS/mpp/include/mpp_reduce_mpi.h
	cp $(SRCROOT)FMS/mpp/include/mpp_reduce_mpi.h .
./xgrid.F90: $(SRCROOT)FMS/exchange/xgrid.F90
	cp $(SRCROOT)FMS/exchange/xgrid.F90 .
./block_control.F90: $(SRCROOT)FMS/block_control/block_control.F90
	cp $(SRCROOT)FMS/block_control/block_control.F90 .
./group_update_pack.inc: $(SRCROOT)FMS/mpp/include/group_update_pack.inc
	cp $(SRCROOT)FMS/mpp/include/group_update_pack.inc .
./mpp_domains_util.inc: $(SRCROOT)FMS/mpp/include/mpp_domains_util.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_domains_util.inc .
./mpp_do_update_ad.h: $(SRCROOT)FMS/mpp/include/mpp_do_update_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_update_ad.h .
./mpp_pset.F90: $(SRCROOT)FMS/mpp/mpp_pset.F90
	cp $(SRCROOT)FMS/mpp/mpp_pset.F90 .
./mpp_gather.h: $(SRCROOT)FMS/mpp/include/mpp_gather.h
	cp $(SRCROOT)FMS/mpp/include/mpp_gather.h .
./mpp_transmit_nocomm.h: $(SRCROOT)FMS/mpp/include/mpp_transmit_nocomm.h
	cp $(SRCROOT)FMS/mpp/include/mpp_transmit_nocomm.h .
./mpp_sum_mpi_ad.h: $(SRCROOT)FMS/mpp/include/mpp_sum_mpi_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_sum_mpi_ad.h .
./gradient_c2l.h: $(SRCROOT)FMS/mosaic/gradient_c2l.h
	cp $(SRCROOT)FMS/mosaic/gradient_c2l.h .
./mpp_io_read.inc: $(SRCROOT)FMS/mpp/include/mpp_io_read.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_io_read.inc .
./atmos_ocean_fluxes.F90: $(SRCROOT)FMS/coupler/atmos_ocean_fluxes.F90
	cp $(SRCROOT)FMS/coupler/atmos_ocean_fluxes.F90 .
./mpp_comm_sma.inc: $(SRCROOT)FMS/mpp/include/mpp_comm_sma.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_comm_sma.inc .
./netcdf_io.F90: $(SRCROOT)FMS/fms2_io/netcdf_io.F90
	cp $(SRCROOT)FMS/fms2_io/netcdf_io.F90 .
./horiz_interp_type.F90: $(SRCROOT)FMS/horiz_interp/horiz_interp_type.F90
	cp $(SRCROOT)FMS/horiz_interp/horiz_interp_type.F90 .
./field_manager.F90: $(SRCROOT)FMS/field_manager/field_manager.F90
	cp $(SRCROOT)FMS/field_manager/field_manager.F90 .
./threadloc.c: $(SRCROOT)FMS/mpp/threadloc.c
	cp $(SRCROOT)FMS/mpp/threadloc.c .
./compute_global_checksum.inc: $(SRCROOT)FMS/fms2_io/include/compute_global_checksum.inc
	cp $(SRCROOT)FMS/fms2_io/include/compute_global_checksum.inc .
./astronomy.F90: $(SRCROOT)FMS/astronomy/astronomy.F90
	cp $(SRCROOT)FMS/astronomy/astronomy.F90 .
./mpp_update_domains2D_ad.h: $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_ad.h .
./monin_obukhov.F90: $(SRCROOT)FMS/monin_obukhov/monin_obukhov.F90
	cp $(SRCROOT)FMS/monin_obukhov/monin_obukhov.F90 .
./MersenneTwister.F90: $(SRCROOT)FMS/random_numbers/MersenneTwister.F90
	cp $(SRCROOT)FMS/random_numbers/MersenneTwister.F90 .
./mpp_write_2Ddecomp.h: $(SRCROOT)FMS/mpp/include/mpp_write_2Ddecomp.h
	cp $(SRCROOT)FMS/mpp/include/mpp_write_2Ddecomp.h .
./mpp_update_domains2D_nonblock.h: $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_nonblock.h
	cp $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_nonblock.h .
./mpp_chksum_int.h: $(SRCROOT)FMS/mpp/include/mpp_chksum_int.h
	cp $(SRCROOT)FMS/mpp/include/mpp_chksum_int.h .
./compressed_read.inc: $(SRCROOT)FMS/fms2_io/include/compressed_read.inc
	cp $(SRCROOT)FMS/fms2_io/include/compressed_read.inc .
./mpp_read_distributed_ascii.inc: $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.inc .
./mpp_type_nocomm.h: $(SRCROOT)FMS/mpp/include/mpp_type_nocomm.h
	cp $(SRCROOT)FMS/mpp/include/mpp_type_nocomm.h .
./diag_output.F90: $(SRCROOT)FMS/diag_manager/diag_output.F90
	cp $(SRCROOT)FMS/diag_manager/diag_output.F90 .
./mpp_unstruct_pass_data.h: $(SRCROOT)FMS/mpp/include/mpp_unstruct_pass_data.h
	cp $(SRCROOT)FMS/mpp/include/mpp_unstruct_pass_data.h .
./array_utils_char.inc: $(SRCROOT)FMS/fms2_io/include/array_utils_char.inc
	cp $(SRCROOT)FMS/fms2_io/include/array_utils_char.inc .
./fms_io_unstructured_register_restart_axis.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_axis.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_axis.inc .
./netcdf_add_restart_variable.inc: $(SRCROOT)FMS/fms2_io/include/netcdf_add_restart_variable.inc
	cp $(SRCROOT)FMS/fms2_io/include/netcdf_add_restart_variable.inc .
./topography.F90: $(SRCROOT)FMS/topography/topography.F90
	cp $(SRCROOT)FMS/topography/topography.F90 .
./horiz_interp.F90: $(SRCROOT)FMS/horiz_interp/horiz_interp.F90
	cp $(SRCROOT)FMS/horiz_interp/horiz_interp.F90 .
./drifters_compute_k.h: $(SRCROOT)FMS/drifters/drifters_compute_k.h
	cp $(SRCROOT)FMS/drifters/drifters_compute_k.h .
./fms_netcdf_unstructured_domain_io.F90: $(SRCROOT)FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90
	cp $(SRCROOT)FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90 .
./time_interp.F90: $(SRCROOT)FMS/time_interp/time_interp.F90
	cp $(SRCROOT)FMS/time_interp/time_interp.F90 .
./netcdf_write_data.inc: $(SRCROOT)FMS/fms2_io/include/netcdf_write_data.inc
	cp $(SRCROOT)FMS/fms2_io/include/netcdf_write_data.inc .
./group_update_unpack.inc: $(SRCROOT)FMS/mpp/include/group_update_unpack.inc
	cp $(SRCROOT)FMS/mpp/include/group_update_unpack.inc .
./mosaic_util.h: $(SRCROOT)FMS/mosaic/mosaic_util.h
	cp $(SRCROOT)FMS/mosaic/mosaic_util.h .
./mpp_write_unlimited_axis.h: $(SRCROOT)FMS/mpp/include/mpp_write_unlimited_axis.h
	cp $(SRCROOT)FMS/mpp/include/mpp_write_unlimited_axis.h .
./mpp_update_nest_domains.h: $(SRCROOT)FMS/mpp/include/mpp_update_nest_domains.h
	cp $(SRCROOT)FMS/mpp/include/mpp_update_nest_domains.h .
./mpp_error_s_s.h: $(SRCROOT)FMS/mpp/include/mpp_error_s_s.h
	cp $(SRCROOT)FMS/mpp/include/mpp_error_s_s.h .
./mpp_io_unstructured_write.inc: $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_write.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_write.inc .
./grid.F90: $(SRCROOT)FMS/mosaic/grid.F90
	cp $(SRCROOT)FMS/mosaic/grid.F90 .
./mpp_utilities.F90: $(SRCROOT)FMS/mpp/mpp_utilities.F90
	cp $(SRCROOT)FMS/mpp/mpp_utilities.F90 .
./mpp_alltoall_sma.h: $(SRCROOT)FMS/mpp/include/mpp_alltoall_sma.h
	cp $(SRCROOT)FMS/mpp/include/mpp_alltoall_sma.h .
./mpp_comm_nocomm.inc: $(SRCROOT)FMS/mpp/include/mpp_comm_nocomm.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_comm_nocomm.inc .
./monin_obukhov_kernel.F90: $(SRCROOT)FMS/monin_obukhov/monin_obukhov_kernel.F90
	cp $(SRCROOT)FMS/monin_obukhov/monin_obukhov_kernel.F90 .
./mpp_transmit_sma.h: $(SRCROOT)FMS/mpp/include/mpp_transmit_sma.h
	cp $(SRCROOT)FMS/mpp/include/mpp_transmit_sma.h .
./axis_utils2.F90: $(SRCROOT)FMS/axis_utils/axis_utils2.F90
	cp $(SRCROOT)FMS/axis_utils/axis_utils2.F90 .
./read_data_3d.inc: $(SRCROOT)FMS/fms/read_data_3d.inc
	cp $(SRCROOT)FMS/fms/read_data_3d.inc .
./mpp_get_boundary_ad.h: $(SRCROOT)FMS/mpp/include/mpp_get_boundary_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_get_boundary_ad.h .
./drifters_io.F90: $(SRCROOT)FMS/drifters/drifters_io.F90
	cp $(SRCROOT)FMS/drifters/drifters_io.F90 .
./mpp_unstruct_domain.inc: $(SRCROOT)FMS/mpp/include/mpp_unstruct_domain.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_unstruct_domain.inc .
./mpp_alltoall_nocomm.h: $(SRCROOT)FMS/mpp/include/mpp_alltoall_nocomm.h
	cp $(SRCROOT)FMS/mpp/include/mpp_alltoall_nocomm.h .
./fms_io_unstructured_read.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_read.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_read.inc .
./horiz_interp_bilinear.F90: $(SRCROOT)FMS/horiz_interp/horiz_interp_bilinear.F90
	cp $(SRCROOT)FMS/horiz_interp/horiz_interp_bilinear.F90 .
./random_numbers.F90: $(SRCROOT)FMS/random_numbers/random_numbers.F90
	cp $(SRCROOT)FMS/random_numbers/random_numbers.F90 .
./drifters_input.F90: $(SRCROOT)FMS/drifters/drifters_input.F90
	cp $(SRCROOT)FMS/drifters/drifters_input.F90 .
./axis_utils.F90: $(SRCROOT)FMS/axis_utils/axis_utils.F90
	cp $(SRCROOT)FMS/axis_utils/axis_utils.F90 .
./mpp_domains_comm.inc: $(SRCROOT)FMS/mpp/include/mpp_domains_comm.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_domains_comm.inc .
./mpp_write.h: $(SRCROOT)FMS/mpp/include/mpp_write.h
	cp $(SRCROOT)FMS/mpp/include/mpp_write.h .
./mpp_data_nocomm.inc: $(SRCROOT)FMS/mpp/include/mpp_data_nocomm.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_data_nocomm.inc .
./register_domain_restart_variable.inc: $(SRCROOT)FMS/fms2_io/include/register_domain_restart_variable.inc
	cp $(SRCROOT)FMS/fms2_io/include/register_domain_restart_variable.inc .
./mpp_chksum_scalar.h: $(SRCROOT)FMS/mpp/include/mpp_chksum_scalar.h
	cp $(SRCROOT)FMS/mpp/include/mpp_chksum_scalar.h .
./mpp_comm.inc: $(SRCROOT)FMS/mpp/include/mpp_comm.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_comm.inc .
./mpp_alltoall_mpi.h: $(SRCROOT)FMS/mpp/include/mpp_alltoall_mpi.h
	cp $(SRCROOT)FMS/mpp/include/mpp_alltoall_mpi.h .
./get_data_type_string.inc: $(SRCROOT)FMS/fms2_io/include/get_data_type_string.inc
	cp $(SRCROOT)FMS/fms2_io/include/get_data_type_string.inc .
./quicksort.F90: $(SRCROOT)FMS/drifters/quicksort.F90
	cp $(SRCROOT)FMS/drifters/quicksort.F90 .
./mpp_comm_mpi.inc: $(SRCROOT)FMS/mpp/include/mpp_comm_mpi.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_comm_mpi.inc .
./read_data_2d.inc: $(SRCROOT)FMS/fms/read_data_2d.inc
	cp $(SRCROOT)FMS/fms/read_data_2d.inc .
./fms_io_unstructured_file_unit.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_file_unit.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_file_unit.inc .
./mpp.F90: $(SRCROOT)FMS/mpp/mpp.F90
	cp $(SRCROOT)FMS/mpp/mpp.F90 .
./mpp_io_write.inc: $(SRCROOT)FMS/mpp/include/mpp_io_write.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_io_write.inc .
./sat_vapor_pres_k.F90: $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres_k.F90
	cp $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres_k.F90 .
./get_checksum.inc: $(SRCROOT)FMS/fms2_io/include/get_checksum.inc
	cp $(SRCROOT)FMS/fms2_io/include/get_checksum.inc .
./affinity.c: $(SRCROOT)FMS/affinity/affinity.c
	cp $(SRCROOT)FMS/affinity/affinity.c .
./fms_io_unstructured_field_exist.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_field_exist.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_field_exist.inc .
./constant.h: $(SRCROOT)FMS/mosaic/constant.h
	cp $(SRCROOT)FMS/mosaic/constant.h .
./mpp_reduce_sma.h: $(SRCROOT)FMS/mpp/include/mpp_reduce_sma.h
	cp $(SRCROOT)FMS/mpp/include/mpp_reduce_sma.h .
./register_variable_attribute.inc: $(SRCROOT)FMS/fms2_io/include/register_variable_attribute.inc
	cp $(SRCROOT)FMS/fms2_io/include/register_variable_attribute.inc .
./mpp_sum_nocomm.h: $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm.h
	cp $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm.h .
./mpp_sum_ad.inc: $(SRCROOT)FMS/mpp/include/mpp_sum_ad.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_sum_ad.inc .
./mpp_efp.F90: $(SRCROOT)FMS/mpp/mpp_efp.F90
	cp $(SRCROOT)FMS/mpp/mpp_efp.F90 .
./fms_io_unstructured_register_restart_field.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_field.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_field.inc .
./mpp_do_get_boundary.h: $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary.h .
./mpp_read_distributed_ascii.h: $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.h
	cp $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.h .
./unstructured_domain_read.inc: $(SRCROOT)FMS/fms2_io/include/unstructured_domain_read.inc
	cp $(SRCROOT)FMS/fms2_io/include/unstructured_domain_read.inc .
./drifters_push.h: $(SRCROOT)FMS/drifters/drifters_push.h
	cp $(SRCROOT)FMS/drifters/drifters_push.h .
./mpp_type_mpi.h: $(SRCROOT)FMS/mpp/include/mpp_type_mpi.h
	cp $(SRCROOT)FMS/mpp/include/mpp_type_mpi.h .
./mpp_group_update.h: $(SRCROOT)FMS/mpp/include/mpp_group_update.h
	cp $(SRCROOT)FMS/mpp/include/mpp_group_update.h .
./fms_io_unstructured_get_file_name.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_name.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_name.inc .
./mpp_scatter.h: $(SRCROOT)FMS/mpp/include/mpp_scatter.h
	cp $(SRCROOT)FMS/mpp/include/mpp_scatter.h .
./gradient.F90: $(SRCROOT)FMS/mosaic/gradient.F90
	cp $(SRCROOT)FMS/mosaic/gradient.F90 .
./mpp_get_boundary.h: $(SRCROOT)FMS/mpp/include/mpp_get_boundary.h
	cp $(SRCROOT)FMS/mpp/include/mpp_get_boundary.h .
./mpp_domains.F90: $(SRCROOT)FMS/mpp/mpp_domains.F90
	cp $(SRCROOT)FMS/mpp/mpp_domains.F90 .
./mpp_global_reduce.h: $(SRCROOT)FMS/mpp/include/mpp_global_reduce.h
	cp $(SRCROOT)FMS/mpp/include/mpp_global_reduce.h .
./mpp_error_s_a.h: $(SRCROOT)FMS/mpp/include/mpp_error_s_a.h
	cp $(SRCROOT)FMS/mpp/include/mpp_error_s_a.h .
./mosaic_util.c: $(SRCROOT)FMS/mosaic/mosaic_util.c
	cp $(SRCROOT)FMS/mosaic/mosaic_util.c .
./array_utils.inc: $(SRCROOT)FMS/fms2_io/include/array_utils.inc
	cp $(SRCROOT)FMS/fms2_io/include/array_utils.inc .
./mpp_error_a_a.h: $(SRCROOT)FMS/mpp/include/mpp_error_a_a.h
	cp $(SRCROOT)FMS/mpp/include/mpp_error_a_a.h .
./write_data.inc: $(SRCROOT)FMS/fms/write_data.inc
	cp $(SRCROOT)FMS/fms/write_data.inc .
./mpp_io_misc.inc: $(SRCROOT)FMS/mpp/include/mpp_io_misc.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_io_misc.inc .
./mpp_do_updateV_ad.h: $(SRCROOT)FMS/mpp/include/mpp_do_updateV_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_updateV_ad.h .
./horiz_interp_spherical.F90: $(SRCROOT)FMS/horiz_interp/horiz_interp_spherical.F90
	cp $(SRCROOT)FMS/horiz_interp/horiz_interp_spherical.F90 .
./fms_io_unstructured_setup_one_field.inc: $(SRCROOT)FMS/fms/fms_io_unstructured_setup_one_field.inc
	cp $(SRCROOT)FMS/fms/fms_io_unstructured_setup_one_field.inc .
./interpolator.F90: $(SRCROOT)FMS/interpolator/interpolator.F90
	cp $(SRCROOT)FMS/interpolator/interpolator.F90 .
./compressed_write.inc: $(SRCROOT)FMS/fms2_io/include/compressed_write.inc
	cp $(SRCROOT)FMS/fms2_io/include/compressed_write.inc .
./sat_vapor_pres.F90: $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres.F90
	cp $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres.F90 .
./mosaic2.F90: $(SRCROOT)FMS/mosaic/mosaic2.F90
	cp $(SRCROOT)FMS/mosaic/mosaic2.F90 .
./mpp_do_redistribute.h: $(SRCROOT)FMS/mpp/include/mpp_do_redistribute.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_redistribute.h .
./mpp_read_2Ddecomp.h: $(SRCROOT)FMS/mpp/include/mpp_read_2Ddecomp.h
	cp $(SRCROOT)FMS/mpp/include/mpp_read_2Ddecomp.h .
./domain_write.inc: $(SRCROOT)FMS/fms2_io/include/domain_write.inc
	cp $(SRCROOT)FMS/fms2_io/include/domain_write.inc .
./mpp_do_update_nest.h: $(SRCROOT)FMS/mpp/include/mpp_do_update_nest.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_update_nest.h .
./mpp_do_update_nonblock.h: $(SRCROOT)FMS/mpp/include/mpp_do_update_nonblock.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_update_nonblock.h .
./horiz_interp_bicubic.F90: $(SRCROOT)FMS/horiz_interp/horiz_interp_bicubic.F90
	cp $(SRCROOT)FMS/horiz_interp/horiz_interp_bicubic.F90 .
./drifters_comm.F90: $(SRCROOT)FMS/drifters/drifters_comm.F90
	cp $(SRCROOT)FMS/drifters/drifters_comm.F90 .
./mpp_domains_reduce.inc: $(SRCROOT)FMS/mpp/include/mpp_domains_reduce.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_domains_reduce.inc .
./read_mosaic.c: $(SRCROOT)FMS/mosaic/read_mosaic.c
	cp $(SRCROOT)FMS/mosaic/read_mosaic.c .
./mpp_global_field.h: $(SRCROOT)FMS/mpp/include/mpp_global_field.h
	cp $(SRCROOT)FMS/mpp/include/mpp_global_field.h .
./mpp_transmit_mpi.h: $(SRCROOT)FMS/mpp/include/mpp_transmit_mpi.h
	cp $(SRCROOT)FMS/mpp/include/mpp_transmit_mpi.h .
./mpp_write_compressed.h: $(SRCROOT)FMS/mpp/include/mpp_write_compressed.h
	cp $(SRCROOT)FMS/mpp/include/mpp_write_compressed.h .
./fms_io_utils.F90: $(SRCROOT)FMS/fms2_io/fms_io_utils.F90
	cp $(SRCROOT)FMS/fms2_io/fms_io_utils.F90 .
./mpp_type_sma.h: $(SRCROOT)FMS/mpp/include/mpp_type_sma.h
	cp $(SRCROOT)FMS/mpp/include/mpp_type_sma.h .
./fm_util.F90: $(SRCROOT)FMS/field_manager/fm_util.F90
	cp $(SRCROOT)FMS/field_manager/fm_util.F90 .
./mpp_sum_nocomm_ad.h: $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm_ad.h
	cp $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm_ad.h .
./mpp_sum_sma.h: $(SRCROOT)FMS/mpp/include/mpp_sum_sma.h
	cp $(SRCROOT)FMS/mpp/include/mpp_sum_sma.h .
./drifters_set_field.h: $(SRCROOT)FMS/drifters/drifters_set_field.h
	cp $(SRCROOT)FMS/drifters/drifters_set_field.h .
./diag_util.F90: $(SRCROOT)FMS/diag_manager/diag_util.F90
	cp $(SRCROOT)FMS/diag_manager/diag_util.F90 .
./mpp_domains_misc.inc: $(SRCROOT)FMS/mpp/include/mpp_domains_misc.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_domains_misc.inc .
./data_override.F90: $(SRCROOT)FMS/data_override/data_override.F90
	cp $(SRCROOT)FMS/data_override/data_override.F90 .
./diag_manager.F90: $(SRCROOT)FMS/diag_manager/diag_manager.F90
	cp $(SRCROOT)FMS/diag_manager/diag_manager.F90 .
./gradient_c2l.c: $(SRCROOT)FMS/mosaic/gradient_c2l.c
	cp $(SRCROOT)FMS/mosaic/gradient_c2l.c .
./mpp_error_a_s.h: $(SRCROOT)FMS/mpp/include/mpp_error_a_s.h
	cp $(SRCROOT)FMS/mpp/include/mpp_error_a_s.h .
./get_variable_attribute.inc: $(SRCROOT)FMS/fms2_io/include/get_variable_attribute.inc
	cp $(SRCROOT)FMS/fms2_io/include/get_variable_attribute.inc .
./diag_grid.F90: $(SRCROOT)FMS/diag_manager/diag_grid.F90
	cp $(SRCROOT)FMS/diag_manager/diag_grid.F90 .
./mpp_io_util.inc: $(SRCROOT)FMS/mpp/include/mpp_io_util.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_io_util.inc .
./gaussian_topog.F90: $(SRCROOT)FMS/topography/gaussian_topog.F90
	cp $(SRCROOT)FMS/topography/gaussian_topog.F90 .
./mpp_do_updateV_nonblock.h: $(SRCROOT)FMS/mpp/include/mpp_do_updateV_nonblock.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_updateV_nonblock.h .
./mpp_util_mpi.inc: $(SRCROOT)FMS/mpp/include/mpp_util_mpi.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_util_mpi.inc .
./mpp_sum.inc: $(SRCROOT)FMS/mpp/include/mpp_sum.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_sum.inc .
./column_diagnostics.F90: $(SRCROOT)FMS/column_diagnostics/column_diagnostics.F90
	cp $(SRCROOT)FMS/column_diagnostics/column_diagnostics.F90 .
./memutils.F90: $(SRCROOT)FMS/memutils/memutils.F90
	cp $(SRCROOT)FMS/memutils/memutils.F90 .
./diag_axis.F90: $(SRCROOT)FMS/diag_manager/diag_axis.F90
	cp $(SRCROOT)FMS/diag_manager/diag_axis.F90 .
./mpp_do_checkV.h: $(SRCROOT)FMS/mpp/include/mpp_do_checkV.h
	cp $(SRCROOT)FMS/mpp/include/mpp_do_checkV.h .
./domain_read.inc: $(SRCROOT)FMS/fms2_io/include/domain_read.inc
	cp $(SRCROOT)FMS/fms2_io/include/domain_read.inc .
./interp.c: $(SRCROOT)FMS/mosaic/interp.c
	cp $(SRCROOT)FMS/mosaic/interp.c .
./mpp_util.inc: $(SRCROOT)FMS/mpp/include/mpp_util.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_util.inc .
./fms2_io.F90: $(SRCROOT)FMS/fms2_io/fms2_io.F90
	cp $(SRCROOT)FMS/fms2_io/fms2_io.F90 .
./time_manager.F90: $(SRCROOT)FMS/time_manager/time_manager.F90
	cp $(SRCROOT)FMS/time_manager/time_manager.F90 .
./mpp_parameter.F90: $(SRCROOT)FMS/mpp/mpp_parameter.F90
	cp $(SRCROOT)FMS/mpp/mpp_parameter.F90 .
./read_data_4d.inc: $(SRCROOT)FMS/fms/read_data_4d.inc
	cp $(SRCROOT)FMS/fms/read_data_4d.inc .
./mpp_sum_mpi.h: $(SRCROOT)FMS/mpp/include/mpp_sum_mpi.h
	cp $(SRCROOT)FMS/mpp/include/mpp_sum_mpi.h .
./time_interp_external.F90: $(SRCROOT)FMS/time_interp/time_interp_external.F90
	cp $(SRCROOT)FMS/time_interp/time_interp_external.F90 .
./fms_affinity.F90: $(SRCROOT)FMS/affinity/fms_affinity.F90
	cp $(SRCROOT)FMS/affinity/fms_affinity.F90 .
SRC = $(SRCROOT)FMS/fms/fms_io.F90 $(SRCROOT)FMS/mosaic/gradient_c2l.c $(SRCROOT)FMS/drifters/quicksort.F90 $(SRCROOT)FMS/mosaic/interp.c $(SRCROOT)FMS/time_interp/time_interp_external2.F90 $(SRCROOT)FMS/exchange/xgrid.F90 $(SRCROOT)FMS/topography/gaussian_topog.F90 $(SRCROOT)FMS/mpp/mpp_io.F90 $(SRCROOT)FMS/mpp/threadloc.c $(SRCROOT)FMS/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)FMS/mpp/mpp_memutils.F90 $(SRCROOT)FMS/fms2_io/fms2_io.F90 $(SRCROOT)FMS/mosaic/mosaic.F90 $(SRCROOT)FMS/drifters/drifters_core.F90 $(SRCROOT)FMS/time_interp/time_interp.F90 $(SRCROOT)FMS/diag_manager/diag_manager.F90 $(SRCROOT)FMS/platform/platform.F90 $(SRCROOT)FMS/mpp/nsclock.c $(SRCROOT)FMS/time_interp/time_interp_external.F90 $(SRCROOT)FMS/fms2_io/fms_io_utils.F90 $(SRCROOT)FMS/mpp/mpp_pset.F90 $(SRCROOT)FMS/horiz_interp/horiz_interp_type.F90 $(SRCROOT)FMS/mosaic/create_xgrid.c $(SRCROOT)FMS/monin_obukhov/monin_obukhov.F90 $(SRCROOT)FMS/diag_integral/diag_integral.F90 $(SRCROOT)FMS/exchange/stock_constants.F90 $(SRCROOT)FMS/time_manager/get_cal_time.F90 $(SRCROOT)FMS/mosaic/gradient.F90 $(SRCROOT)FMS/affinity/fms_affinity.F90 $(SRCROOT)FMS/diag_manager/diag_table.F90 $(SRCROOT)FMS/mpp/mpp.F90 $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)FMS/drifters/drifters_comm.F90 $(SRCROOT)FMS/tridiagonal/tridiagonal.F90 $(SRCROOT)FMS/axis_utils/axis_utils.F90 $(SRCROOT)FMS/memutils/memutils.F90 $(SRCROOT)FMS/fms/fms.F90 $(SRCROOT)FMS/drifters/drifters_input.F90 $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)FMS/mpp/mpp_parameter.F90 $(SRCROOT)FMS/mosaic/read_mosaic.c $(SRCROOT)FMS/data_override/data_override.F90 $(SRCROOT)FMS/drifters/cloud_interpolator.F90 $(SRCROOT)FMS/mpp/mpp_memuse.c $(SRCROOT)FMS/amip_interp/amip_interp.F90 $(SRCROOT)FMS/mpp/mpp_data.F90 $(SRCROOT)FMS/random_numbers/random_numbers.F90 $(SRCROOT)FMS/topography/topography.F90 $(SRCROOT)FMS/constants/constants.F90 $(SRCROOT)FMS/diag_manager/diag_data.F90 $(SRCROOT)FMS/mosaic/mosaic_util.c $(SRCROOT)FMS/mpp/mpp_efp.F90 $(SRCROOT)FMS/field_manager/fm_util.F90 $(SRCROOT)FMS/block_control/block_control.F90 $(SRCROOT)FMS/fms2_io/fms_netcdf_domain_io.F90 $(SRCROOT)FMS/diag_manager/diag_grid.F90 $(SRCROOT)FMS/diag_manager/diag_output.F90 $(SRCROOT)FMS/mosaic/grid.F90 $(SRCROOT)FMS/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)FMS/mosaic/mosaic2.F90 $(SRCROOT)FMS/coupler/atmos_ocean_fluxes.F90 $(SRCROOT)FMS/diag_manager/diag_axis.F90 $(SRCROOT)FMS/astronomy/astronomy.F90 $(SRCROOT)FMS/diag_manager/diag_util.F90 $(SRCROOT)FMS/coupler/coupler_types.F90 $(SRCROOT)FMS/fms2_io/blackboxio.F90 $(SRCROOT)FMS/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)FMS/random_numbers/MersenneTwister.F90 $(SRCROOT)FMS/interpolator/interpolator.F90 $(SRCROOT)FMS/coupler/ensemble_manager.F90 $(SRCROOT)FMS/tracer_manager/tracer_manager.F90 $(SRCROOT)FMS/fms2_io/netcdf_io.F90 $(SRCROOT)FMS/column_diagnostics/column_diagnostics.F90 $(SRCROOT)FMS/drifters/drifters_io.F90 $(SRCROOT)FMS/mpp/mpp_utilities.F90 $(SRCROOT)FMS/time_manager/time_manager.F90 $(SRCROOT)FMS/mpp/mpp_domains.F90 $(SRCROOT)FMS/axis_utils/axis_utils2.F90 $(SRCROOT)FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90 $(SRCROOT)FMS/horiz_interp/horiz_interp.F90 $(SRCROOT)FMS/affinity/affinity.c $(SRCROOT)FMS/field_manager/field_manager.F90 $(SRCROOT)FMS/drifters/drifters.F90 $(SRCROOT)FMS/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)FMS/monin_obukhov/monin_obukhov_kernel.F90 $(SRCROOT)FMS/mpp/include/mpp_util_sma.inc $(SRCROOT)FMS/mpp/include/mpp_type_mpi.h $(SRCROOT)FMS/mpp/include/mpp_global_field_ad.h $(SRCROOT)FMS/mpp/include/mpp_scatter.h $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_name.inc $(SRCROOT)FMS/mpp/include/mpp_group_update.h $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_read.inc $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.h $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary.h $(SRCROOT)FMS/mpp/include/system_clock.h $(SRCROOT)FMS/drifters/drifters_push.h $(SRCROOT)FMS/fms2_io/include/unstructured_domain_read.inc $(SRCROOT)FMS/mpp/include/mpp_transmit.inc $(SRCROOT)FMS/mpp/include/mpp_error_s_a.h $(SRCROOT)FMS/mpp/include/mpp_error_a_a.h $(SRCROOT)FMS/fms2_io/include/array_utils.inc $(SRCROOT)FMS/mpp/include/mpp_get_boundary.h $(SRCROOT)FMS/mpp/include/mpp_global_reduce.h $(SRCROOT)FMS/fms2_io/include/register_variable_attribute.inc $(SRCROOT)FMS/mpp/include/mpp_reduce_sma.h $(SRCROOT)FMS/mosaic/constant.h $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_field.inc $(SRCROOT)FMS/mpp/include/mpp_read_compressed.h $(SRCROOT)FMS/mpp/include/mpp_global_sum_tl.h $(SRCROOT)FMS/mpp/include/mpp_data_mpi.inc $(SRCROOT)FMS/fms2_io/include/unstructured_domain_write.inc $(SRCROOT)FMS/mpp/include/mpp_sum_ad.inc $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm.h $(SRCROOT)FMS/fms/fms_io_unstructured_file_unit.inc $(SRCROOT)FMS/mosaic/interp.h $(SRCROOT)FMS/mpp/include/mpp_global_sum_ad.h $(SRCROOT)FMS/fms/read_data_2d.inc $(SRCROOT)FMS/fms/fms_io_unstructured_save_restart.inc $(SRCROOT)FMS/mpp/include/mpp_do_update.h $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/mpp/include/mpp_do_check.h $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/mpp/include/mpp_comm_mpi.inc $(SRCROOT)FMS/fms2_io/include/register_unstructured_domain_restart_variable.inc $(SRCROOT)FMS/fms/fms_io_unstructured_field_exist.inc $(SRCROOT)FMS/fms2_io/include/get_checksum.inc $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary_ad.h $(SRCROOT)FMS/mpp/include/mpp_io_write.inc $(SRCROOT)FMS/mpp/include/mpp_global_sum.h $(SRCROOT)FMS/mpp/include/mpp_reduce_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_write.h $(SRCROOT)FMS/field_manager/parse.inc $(SRCROOT)FMS/mpp/include/mpp_define_nest_domains.inc $(SRCROOT)FMS/fms2_io/include/register_domain_restart_variable.inc $(SRCROOT)FMS/mosaic/read_mosaic.h $(SRCROOT)FMS/mpp/include/mpp_data_nocomm.inc $(SRCROOT)FMS/fms/fms_io_unstructured_get_field_size.inc $(SRCROOT)FMS/mpp/include/mpp_update_domains2D.h $(SRCROOT)FMS/mpp/include/mpp_util_nocomm.inc $(SRCROOT)FMS/mpp/include/mpp_domains_comm.inc $(SRCROOT)FMS/mpp/include/mpp_global_field_ug.h $(SRCROOT)FMS/fms/../include/fms_platform.h $(SRCROOT)FMS/fms2_io/include/netcdf_read_data.inc $(SRCROOT)FMS/fms2_io/include/get_data_type_string.inc $(SRCROOT)FMS/mpp/include/mpp_chksum.h $(SRCROOT)FMS/mpp/include/mpp_do_global_field.h $(SRCROOT)FMS/fms2_io/include/get_global_attribute.inc $(SRCROOT)FMS/mpp/include/mpp_chksum_scalar.h $(SRCROOT)FMS/mpp/include/mpp_io_connect.inc $(SRCROOT)FMS/mpp/include/mpp_alltoall_mpi.h $(SRCROOT)FMS/mpp/include/mpp_comm.inc $(SRCROOT)FMS/mpp/include/mpp_util.inc $(SRCROOT)FMS/mosaic/mosaic_util.h $(SRCROOT)FMS/mpp/include/group_update_unpack.inc $(SRCROOT)FMS/mpp/include/mpp_error_s_s.h $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_write.inc $(SRCROOT)FMS/mpp/include/mpp_write_unlimited_axis.h $(SRCROOT)FMS/mpp/include/mpp_update_nest_domains.h $(SRCROOT)FMS/drifters/drifters_compute_k.h $(SRCROOT)FMS/fms2_io/include/netcdf_write_data.inc $(SRCROOT)FMS/mpp/include/mpp_get_boundary_ad.h $(SRCROOT)FMS/mpp/include/mpp_sum_mpi.h $(SRCROOT)FMS/mpp/include/mpp_alltoall_nocomm.h $(SRCROOT)FMS/fms/fms_io_unstructured_read.inc $(SRCROOT)FMS/mpp/include/mpp_unstruct_domain.inc $(SRCROOT)FMS/fms/read_data_4d.inc $(SRCROOT)FMS/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)FMS/mpp/include/mpp_alltoall_sma.h $(SRCROOT)FMS/fms/read_data_3d.inc $(SRCROOT)FMS/mpp/include/mpp_transmit_sma.h $(SRCROOT)FMS/fms2_io/include/get_variable_attribute.inc $(SRCROOT)FMS/mpp/include/mpp_unstruct_pass_data.h $(SRCROOT)FMS/mpp/include/mpp_error_a_s.h $(SRCROOT)FMS/fms2_io/include/compressed_read.inc $(SRCROOT)FMS/mpp/include/mpp_chksum_int.h $(SRCROOT)FMS/mpp/include/mpp_type_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.inc $(SRCROOT)FMS/mpp/include/mpp_domains_misc.inc $(SRCROOT)FMS/fms2_io/include/netcdf_add_restart_variable.inc $(SRCROOT)FMS/mpp/include/mpp_sum.inc $(SRCROOT)FMS/fms2_io/include/domain_read.inc $(SRCROOT)FMS/mpp/include/mpp_do_checkV.h $(SRCROOT)FMS/mpp/include/mpp_do_updateV_nonblock.h $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_axis.inc $(SRCROOT)FMS/mpp/include/mpp_util_mpi.inc $(SRCROOT)FMS/fms2_io/include/array_utils_char.inc $(SRCROOT)FMS/mpp/include/mpp_io_util.inc $(SRCROOT)FMS/mpp/include/mpp_transmit_mpi.h $(SRCROOT)FMS/mpp/include/mpp_comm_sma.inc $(SRCROOT)FMS/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_sum_mpi_ad.h $(SRCROOT)FMS/mosaic/gradient_c2l.h $(SRCROOT)FMS/mpp/include/mpp_domains_reduce.inc $(SRCROOT)FMS/mpp/include/mpp_gather.h $(SRCROOT)FMS/mpp/include/mpp_io_read.inc $(SRCROOT)FMS/mpp/include/mpp_global_field.h $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm_ad.h $(SRCROOT)FMS/mpp/include/mpp_type_sma.h $(SRCROOT)FMS/drifters/drifters_set_field.h $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_nonblock.h $(SRCROOT)FMS/mpp/include/mpp_write_2Ddecomp.h $(SRCROOT)FMS/mpp/include/mpp_sum_sma.h $(SRCROOT)FMS/mpp/include/mpp_write_compressed.h $(SRCROOT)FMS/fms2_io/include/compute_global_checksum.inc $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_ad.h $(SRCROOT)FMS/mpp/include/mpp_do_redistribute.h $(SRCROOT)FMS/fms2_io/include/compressed_write.inc $(SRCROOT)FMS/mpp/include/mpp_domains_define.inc $(SRCROOT)FMS/fms2_io/include/domain_write.inc $(SRCROOT)FMS/fms2_io/include/register_global_attribute.inc $(SRCROOT)FMS/mpp/include/mpp_read_2Ddecomp.h $(SRCROOT)FMS/mpp/include/mpp_do_global_field_ad.h $(SRCROOT)FMS/mpp/include/mpp_do_updateV.h $(SRCROOT)FMS/mpp/include/mpp_do_updateV_ad.h $(SRCROOT)FMS/mpp/include/mpp_io_misc.inc $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_unit.inc $(SRCROOT)FMS/fms/write_data.inc $(SRCROOT)FMS/mosaic/create_xgrid.h $(SRCROOT)FMS/fms/fms_io_unstructured_setup_one_field.inc $(SRCROOT)FMS/drifters/fms_switches.h $(SRCROOT)FMS/mpp/include/mpp_data_sma.inc $(SRCROOT)FMS/mpp/include/mpp_do_update_nonblock.h $(SRCROOT)FMS/mpp/include/mpp_do_update_nest.h $(SRCROOT)FMS/mpp/include/mpp_do_update_ad.h $(SRCROOT)FMS/mpp/include/mpp_domains_util.inc $(SRCROOT)FMS/mpp/include/group_update_pack.inc $(SRCROOT)FMS/mpp/include/mpp_reduce_mpi.h
OBJ = fms_io.o gradient_c2l.o quicksort.o interp.o time_interp_external2.o xgrid.o gaussian_topog.o mpp_io.o threadloc.o horiz_interp_bicubic.o mpp_memutils.o fms2_io.o mosaic.o drifters_core.o time_interp.o diag_manager.o platform.o nsclock.o time_interp_external.o fms_io_utils.o mpp_pset.o horiz_interp_type.o create_xgrid.o monin_obukhov.o diag_integral.o stock_constants.o get_cal_time.o gradient.o fms_affinity.o diag_table.o mpp.o sat_vapor_pres_k.o drifters_comm.o tridiagonal.o axis_utils.o memutils.o fms.o drifters_input.o sat_vapor_pres.o mpp_parameter.o read_mosaic.o data_override.o cloud_interpolator.o mpp_memuse.o amip_interp.o mpp_data.o random_numbers.o topography.o constants.o diag_data.o mosaic_util.o mpp_efp.o fm_util.o block_control.o fms_netcdf_domain_io.o diag_grid.o diag_output.o grid.o horiz_interp_spherical.o mosaic2.o atmos_ocean_fluxes.o diag_axis.o astronomy.o diag_util.o coupler_types.o blackboxio.o horiz_interp_bilinear.o MersenneTwister.o interpolator.o ensemble_manager.o tracer_manager.o netcdf_io.o column_diagnostics.o drifters_io.o mpp_utilities.o time_manager.o mpp_domains.o axis_utils2.o fms_netcdf_unstructured_domain_io.o horiz_interp.o affinity.o field_manager.o drifters.o horiz_interp_conserve.o monin_obukhov_kernel.o
OFF = $(SRCROOT)FMS/mosaic/create_xgrid.c $(SRCROOT)FMS/mpp/include/mpp_util_nocomm.inc $(SRCROOT)FMS/mpp/include/mpp_global_field_ug.h $(SRCROOT)FMS/tracer_manager/tracer_manager.F90 $(SRCROOT)FMS/drifters/drifters.F90 $(SRCROOT)FMS/mpp/include/mpp_define_nest_domains.inc $(SRCROOT)FMS/drifters/drifters_core.F90 $(SRCROOT)FMS/field_manager/parse.inc $(SRCROOT)FMS/mpp/include/mpp_update_domains2D.h $(SRCROOT)FMS/fms/fms_io_unstructured_get_field_size.inc $(SRCROOT)FMS/mosaic/read_mosaic.h $(SRCROOT)FMS/mpp/include/mpp_io_connect.inc $(SRCROOT)FMS/diag_integral/diag_integral.F90 $(SRCROOT)FMS/fms2_io/include/get_global_attribute.inc $(SRCROOT)FMS/drifters/cloud_interpolator.F90 $(SRCROOT)FMS/mpp/include/mpp_do_global_field.h $(SRCROOT)FMS/fms/../include/fms_platform.h $(SRCROOT)FMS/mpp/include/mpp_chksum.h $(SRCROOT)FMS/fms2_io/include/netcdf_read_data.inc $(SRCROOT)FMS/platform/platform.F90 $(SRCROOT)FMS/fms2_io/include/register_unstructured_domain_restart_variable.inc $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/fms/fms_io.F90 $(SRCROOT)FMS/mpp/mpp_memuse.c $(SRCROOT)FMS/mosaic/interp.h $(SRCROOT)FMS/mpp/include/mpp_global_sum_ad.h $(SRCROOT)FMS/mosaic/mosaic.F90 $(SRCROOT)FMS/mpp/include/mpp_do_check.h $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/mpp/include/mpp_do_update.h $(SRCROOT)FMS/fms/fms_io_unstructured_save_restart.inc $(SRCROOT)FMS/mpp/include/mpp_reduce_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_global_sum.h $(SRCROOT)FMS/fms2_io/fms_netcdf_domain_io.F90 $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary_ad.h $(SRCROOT)FMS/amip_interp/amip_interp.F90 $(SRCROOT)FMS/diag_manager/diag_table.F90 $(SRCROOT)FMS/mpp/nsclock.c $(SRCROOT)FMS/fms2_io/include/unstructured_domain_write.inc $(SRCROOT)FMS/mpp/include/mpp_data_mpi.inc $(SRCROOT)FMS/tridiagonal/tridiagonal.F90 $(SRCROOT)FMS/mpp/include/mpp_global_sum_tl.h $(SRCROOT)FMS/mpp/include/mpp_read_compressed.h $(SRCROOT)FMS/mpp/include/system_clock.h $(SRCROOT)FMS/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)FMS/constants/constants.F90 $(SRCROOT)FMS/mpp/include/mpp_transmit.inc $(SRCROOT)FMS/mpp/include/mpp_global_field_ad.h $(SRCROOT)FMS/mpp/include/mpp_util_sma.inc $(SRCROOT)FMS/coupler/coupler_types.F90 $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_read.inc $(SRCROOT)FMS/fms2_io/blackboxio.F90 $(SRCROOT)FMS/diag_manager/diag_data.F90 $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_unit.inc $(SRCROOT)FMS/time_manager/get_cal_time.F90 $(SRCROOT)FMS/time_interp/time_interp_external2.F90 $(SRCROOT)FMS/mpp/include/mpp_do_updateV.h $(SRCROOT)FMS/mpp/include/mpp_do_global_field_ad.h $(SRCROOT)FMS/mpp/include/mpp_data_sma.inc $(SRCROOT)FMS/drifters/fms_switches.h $(SRCROOT)FMS/mosaic/create_xgrid.h $(SRCROOT)FMS/mpp/mpp_memutils.F90 $(SRCROOT)FMS/coupler/ensemble_manager.F90 $(SRCROOT)FMS/mpp/mpp_data.F90 $(SRCROOT)FMS/mpp/mpp_io.F90 $(SRCROOT)FMS/exchange/stock_constants.F90 $(SRCROOT)FMS/fms2_io/include/register_global_attribute.inc $(SRCROOT)FMS/fms/fms.F90 $(SRCROOT)FMS/mpp/include/mpp_domains_define.inc $(SRCROOT)FMS/mpp/include/mpp_reduce_mpi.h $(SRCROOT)FMS/exchange/xgrid.F90 $(SRCROOT)FMS/block_control/block_control.F90 $(SRCROOT)FMS/mpp/include/group_update_pack.inc $(SRCROOT)FMS/mpp/include/mpp_domains_util.inc $(SRCROOT)FMS/mpp/include/mpp_do_update_ad.h $(SRCROOT)FMS/mpp/mpp_pset.F90 $(SRCROOT)FMS/mpp/include/mpp_gather.h $(SRCROOT)FMS/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_sum_mpi_ad.h $(SRCROOT)FMS/mosaic/gradient_c2l.h $(SRCROOT)FMS/mpp/include/mpp_io_read.inc $(SRCROOT)FMS/coupler/atmos_ocean_fluxes.F90 $(SRCROOT)FMS/mpp/include/mpp_comm_sma.inc $(SRCROOT)FMS/fms2_io/netcdf_io.F90 $(SRCROOT)FMS/horiz_interp/horiz_interp_type.F90 $(SRCROOT)FMS/field_manager/field_manager.F90 $(SRCROOT)FMS/mpp/threadloc.c $(SRCROOT)FMS/fms2_io/include/compute_global_checksum.inc $(SRCROOT)FMS/astronomy/astronomy.F90 $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_ad.h $(SRCROOT)FMS/monin_obukhov/monin_obukhov.F90 $(SRCROOT)FMS/random_numbers/MersenneTwister.F90 $(SRCROOT)FMS/mpp/include/mpp_write_2Ddecomp.h $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_nonblock.h $(SRCROOT)FMS/mpp/include/mpp_chksum_int.h $(SRCROOT)FMS/fms2_io/include/compressed_read.inc $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.inc $(SRCROOT)FMS/mpp/include/mpp_type_nocomm.h $(SRCROOT)FMS/diag_manager/diag_output.F90 $(SRCROOT)FMS/mpp/include/mpp_unstruct_pass_data.h $(SRCROOT)FMS/fms2_io/include/array_utils_char.inc $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_axis.inc $(SRCROOT)FMS/fms2_io/include/netcdf_add_restart_variable.inc $(SRCROOT)FMS/topography/topography.F90 $(SRCROOT)FMS/horiz_interp/horiz_interp.F90 $(SRCROOT)FMS/drifters/drifters_compute_k.h $(SRCROOT)FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90 $(SRCROOT)FMS/time_interp/time_interp.F90 $(SRCROOT)FMS/fms2_io/include/netcdf_write_data.inc $(SRCROOT)FMS/mpp/include/group_update_unpack.inc $(SRCROOT)FMS/mosaic/mosaic_util.h $(SRCROOT)FMS/mpp/include/mpp_write_unlimited_axis.h $(SRCROOT)FMS/mpp/include/mpp_update_nest_domains.h $(SRCROOT)FMS/mpp/include/mpp_error_s_s.h $(SRCROOT)FMS/mpp/include/mpp_io_unstructured_write.inc $(SRCROOT)FMS/mosaic/grid.F90 $(SRCROOT)FMS/mpp/mpp_utilities.F90 $(SRCROOT)FMS/mpp/include/mpp_alltoall_sma.h $(SRCROOT)FMS/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)FMS/monin_obukhov/monin_obukhov_kernel.F90 $(SRCROOT)FMS/mpp/include/mpp_transmit_sma.h $(SRCROOT)FMS/axis_utils/axis_utils2.F90 $(SRCROOT)FMS/fms/read_data_3d.inc $(SRCROOT)FMS/mpp/include/mpp_get_boundary_ad.h $(SRCROOT)FMS/drifters/drifters_io.F90 $(SRCROOT)FMS/mpp/include/mpp_unstruct_domain.inc $(SRCROOT)FMS/mpp/include/mpp_alltoall_nocomm.h $(SRCROOT)FMS/fms/fms_io_unstructured_read.inc $(SRCROOT)FMS/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)FMS/random_numbers/random_numbers.F90 $(SRCROOT)FMS/drifters/drifters_input.F90 $(SRCROOT)FMS/axis_utils/axis_utils.F90 $(SRCROOT)FMS/mpp/include/mpp_domains_comm.inc $(SRCROOT)FMS/mpp/include/mpp_write.h $(SRCROOT)FMS/mpp/include/mpp_data_nocomm.inc $(SRCROOT)FMS/fms2_io/include/register_domain_restart_variable.inc $(SRCROOT)FMS/mpp/include/mpp_chksum_scalar.h $(SRCROOT)FMS/mpp/include/mpp_comm.inc $(SRCROOT)FMS/mpp/include/mpp_alltoall_mpi.h $(SRCROOT)FMS/fms2_io/include/get_data_type_string.inc $(SRCROOT)FMS/drifters/quicksort.F90 $(SRCROOT)FMS/mpp/include/mpp_comm_mpi.inc $(SRCROOT)FMS/fms/read_data_2d.inc $(SRCROOT)FMS/fms/fms_io_unstructured_file_unit.inc $(SRCROOT)FMS/mpp/mpp.F90 $(SRCROOT)FMS/mpp/include/mpp_io_write.inc $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)FMS/fms2_io/include/get_checksum.inc $(SRCROOT)FMS/affinity/affinity.c $(SRCROOT)FMS/fms/fms_io_unstructured_field_exist.inc $(SRCROOT)FMS/mosaic/constant.h $(SRCROOT)FMS/mpp/include/mpp_reduce_sma.h $(SRCROOT)FMS/fms2_io/include/register_variable_attribute.inc $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm.h $(SRCROOT)FMS/mpp/include/mpp_sum_ad.inc $(SRCROOT)FMS/mpp/mpp_efp.F90 $(SRCROOT)FMS/fms/fms_io_unstructured_register_restart_field.inc $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary.h $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.h $(SRCROOT)FMS/fms2_io/include/unstructured_domain_read.inc $(SRCROOT)FMS/drifters/drifters_push.h $(SRCROOT)FMS/mpp/include/mpp_type_mpi.h $(SRCROOT)FMS/mpp/include/mpp_group_update.h $(SRCROOT)FMS/fms/fms_io_unstructured_get_file_name.inc $(SRCROOT)FMS/mpp/include/mpp_scatter.h $(SRCROOT)FMS/mosaic/gradient.F90 $(SRCROOT)FMS/mpp/include/mpp_get_boundary.h $(SRCROOT)FMS/mpp/mpp_domains.F90 $(SRCROOT)FMS/mpp/include/mpp_global_reduce.h $(SRCROOT)FMS/mpp/include/mpp_error_s_a.h $(SRCROOT)FMS/mosaic/mosaic_util.c $(SRCROOT)FMS/fms2_io/include/array_utils.inc $(SRCROOT)FMS/mpp/include/mpp_error_a_a.h $(SRCROOT)FMS/fms/write_data.inc $(SRCROOT)FMS/mpp/include/mpp_io_misc.inc $(SRCROOT)FMS/mpp/include/mpp_do_updateV_ad.h $(SRCROOT)FMS/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)FMS/fms/fms_io_unstructured_setup_one_field.inc $(SRCROOT)FMS/interpolator/interpolator.F90 $(SRCROOT)FMS/fms2_io/include/compressed_write.inc $(SRCROOT)FMS/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)FMS/mosaic/mosaic2.F90 $(SRCROOT)FMS/mpp/include/mpp_do_redistribute.h $(SRCROOT)FMS/mpp/include/mpp_read_2Ddecomp.h $(SRCROOT)FMS/fms2_io/include/domain_write.inc $(SRCROOT)FMS/mpp/include/mpp_do_update_nest.h $(SRCROOT)FMS/mpp/include/mpp_do_update_nonblock.h $(SRCROOT)FMS/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)FMS/drifters/drifters_comm.F90 $(SRCROOT)FMS/mpp/include/mpp_domains_reduce.inc $(SRCROOT)FMS/mosaic/read_mosaic.c $(SRCROOT)FMS/mpp/include/mpp_global_field.h $(SRCROOT)FMS/mpp/include/mpp_transmit_mpi.h $(SRCROOT)FMS/mpp/include/mpp_write_compressed.h $(SRCROOT)FMS/fms2_io/fms_io_utils.F90 $(SRCROOT)FMS/mpp/include/mpp_type_sma.h $(SRCROOT)FMS/field_manager/fm_util.F90 $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm_ad.h $(SRCROOT)FMS/mpp/include/mpp_sum_sma.h $(SRCROOT)FMS/drifters/drifters_set_field.h $(SRCROOT)FMS/diag_manager/diag_util.F90 $(SRCROOT)FMS/mpp/include/mpp_domains_misc.inc $(SRCROOT)FMS/data_override/data_override.F90 $(SRCROOT)FMS/diag_manager/diag_manager.F90 $(SRCROOT)FMS/mosaic/gradient_c2l.c $(SRCROOT)FMS/mpp/include/mpp_error_a_s.h $(SRCROOT)FMS/fms2_io/include/get_variable_attribute.inc $(SRCROOT)FMS/diag_manager/diag_grid.F90 $(SRCROOT)FMS/mpp/include/mpp_io_util.inc $(SRCROOT)FMS/topography/gaussian_topog.F90 $(SRCROOT)FMS/mpp/include/mpp_do_updateV_nonblock.h $(SRCROOT)FMS/mpp/include/mpp_util_mpi.inc $(SRCROOT)FMS/mpp/include/mpp_sum.inc $(SRCROOT)FMS/column_diagnostics/column_diagnostics.F90 $(SRCROOT)FMS/memutils/memutils.F90 $(SRCROOT)FMS/diag_manager/diag_axis.F90 $(SRCROOT)FMS/mpp/include/mpp_do_checkV.h $(SRCROOT)FMS/fms2_io/include/domain_read.inc $(SRCROOT)FMS/mosaic/interp.c $(SRCROOT)FMS/mpp/include/mpp_util.inc $(SRCROOT)FMS/fms2_io/fms2_io.F90 $(SRCROOT)FMS/time_manager/time_manager.F90 $(SRCROOT)FMS/mpp/mpp_parameter.F90 $(SRCROOT)FMS/fms/read_data_4d.inc $(SRCROOT)FMS/mpp/include/mpp_sum_mpi.h $(SRCROOT)FMS/time_interp/time_interp_external.F90 $(SRCROOT)FMS/affinity/fms_affinity.F90
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

