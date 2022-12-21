# Makefile created by mkmf 19.3.1
CPPDEFS = -DINTERNAL_FILE_NML -DSPMD -DCLIMATE_NUDGE

OTHERFLAGS = -I$(BUILDROOT)atmos_phys -I$(BUILDROOT)fms/build/.mods

TMPFILES = *.o *.mod

include $(MK_TEMPLATE)

.DEFAULT:
	-echo $@ does not exist.
all: libatmos_cubed_sphere.a
a2b_edge.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/a2b_edge.F90 fv_grid_utils.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/model/a2b_edge.F90
atmosphere.o: $(SRCROOT)GFDL_atmos_cubed_sphere/driver/GFDL/atmosphere.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h fv_arrays.o fv_control.o fv_eta.o fv_dynamics.o fv_nesting.o fv_diagnostics.o fv_cmip_diag.o fv_restart.o fv_timing.o fv_mp_mod.o fv_sg.o fv_update_phys.o fv_io.o fv_regional_bc.o fv_climate_nudge.o fv_ada_nudge.o fv_nudge.o gfdl_mp.o cloud_diagnosis.o coarse_graining.o coarse_grained_diagnostics.o coarse_grained_restart_files.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)GFDL_atmos_cubed_sphere/driver/GFDL/atmosphere.F90
boundary.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/boundary.F90 fv_mp_mod.o fv_arrays.o fv_timing.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/model/boundary.F90
cloud_diagnosis.o: $(SRCROOT)GFDL_atmos_cubed_sphere/driver/SHiELD/cloud_diagnosis.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/driver/SHiELD/cloud_diagnosis.F90
coarse_grained_diagnostics.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_grained_diagnostics.F90 fv_arrays.o fv_diagnostics.o fv_mapz.o coarse_graining.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_grained_diagnostics.F90
coarse_grained_restart_files.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_grained_restart_files.F90 coarse_graining.o fv_arrays.o fv_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_grained_restart_files.F90
coarse_graining.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_graining.F90 fv_mapz.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_graining.F90
dyn_core.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/dyn_core.F90 fv_mp_mod.o sw_core.o a2b_edge.o nh_core.o tp_core.o fv_timing.o fv_diagnostics.o fv_diag_column.o fv_update_phys.o fv_ada_nudge.o fv_nudge.o fv_arrays.o boundary.o fv_regional_bc.o test_cases.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/model/dyn_core.F90
external_ic.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/external_ic.F90 $(SRCROOT)FMS/include/file_version.h external_sst.o fv_arrays.o fv_diagnostics.o fv_grid_utils.o fv_io.o fv_mapz.o fv_regional_bc.o fv_mp_mod.o fv_surf_map.o fv_timing.o init_hydro.o fv_fill.o fv_eta.o sim_nc_mod.o fv_nudge.o test_cases.o boundary.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/external_ic.F90
external_sst.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/external_sst.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/external_sst.F90
fv_ada_nudge.o: $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_ada_nudge.F90 $(SRCROOT)FMS/include/file_version.h external_sst.o fv_grid_utils.o fv_diagnostics.o tp_core.o fv_mapz.o fv_mp_mod.o fv_timing.o sim_nc_mod.o fv_arrays.o fv_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_ada_nudge.F90
fv_arrays.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_arrays.F90 $(SRCROOT)FMS/include/fms_platform.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_arrays.F90
fv_climate_nudge.o: $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_climate_nudge.F90 $(SRCROOT)FMS/include/file_version.h read_climate_nudge_data.o fv_mapz.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_climate_nudge.F90
fv_cmip_diag.o: $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_cmip_diag.F90 $(SRCROOT)FMS/include/file_version.h fv_mapz.o fv_arrays.o fv_diagnostics.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_cmip_diag.F90
fv_cmp.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_cmp.F90 fv_arrays.o gfdl_cloud_microphys.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_cmp.F90
fv_control.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_control.F90 $(SRCROOT)FMS/include/file_version.h fv_io.o fv_restart.o fv_arrays.o fv_grid_utils.o fv_eta.o fv_grid_tools.o fv_mp_mod.o test_cases.o fv_timing.o fv_diagnostics.o coarse_grained_restart_files.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_control.F90
fv_diag_column.o: $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_diag_column.F90 $(SRCROOT)FMS/include/file_version.h fv_arrays.o fv_grid_utils.o fv_sg.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_diag_column.F90
fv_diagnostics.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_diagnostics.F90 $(SRCROOT)FMS/include/file_version.h $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_diagnostics.h fv_arrays.o fv_mapz.o fv_mp_mod.o fv_eta.o fv_grid_utils.o a2b_edge.o fv_surf_map.o fv_sg.o gfdl_mp.o fv_diag_column.o gfdl_cloud_microphys.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include -I$(SRCROOT)GFDL_atmos_cubed_sphere/tools	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_diagnostics.F90
fv_dynamics.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_dynamics.F90 dyn_core.o fv_mapz.o fv_tracer2d.o fv_grid_utils.o fv_fill.o fv_mp_mod.o fv_timing.o fv_diagnostics.o fv_sg.o fv_nesting.o fv_regional_bc.o boundary.o fv_arrays.o fv_nudge.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_dynamics.F90
fv_eta.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_eta.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_eta.h fv_mp_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)GFDL_atmos_cubed_sphere/tools	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_eta.F90
fv_fill.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_fill.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_fill.F90
fv_grid_tools.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_grid_tools.F90 fv_arrays.o fv_grid_utils.o fv_timing.o fv_mp_mod.o sorted_index.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_grid_tools.F90
fv_grid_utils.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_grid_utils.F90 $(SRCROOT)FMS/include/fms_platform.h external_sst.o fv_arrays.o fv_eta.o fv_mp_mod.o fv_timing.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_grid_utils.F90
fv_io.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_io.F90 external_sst.o fv_arrays.o fv_eta.o fv_mp_mod.o fv_treat_da_inc.o fv_mapz.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_io.F90
fv_mapz.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_mapz.F90 fv_grid_utils.o fv_fill.o fv_arrays.o fv_timing.o fv_mp_mod.o fv_cmp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_mapz.F90
fv_mp_mod.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_mp_mod.F90 fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_mp_mod.F90
fv_nesting.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_nesting.F90 fv_sg.o boundary.o fv_arrays.o fv_grid_utils.o init_hydro.o fv_mapz.o fv_timing.o fv_mp_mod.o fv_diagnostics.o sw_core.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_nesting.F90
fv_nggps_diag.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_nggps_diag.F90 fv_diagnostics.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_nggps_diag.F90
fv_nudge.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_nudge.F90 $(SRCROOT)FMS/include/file_version.h external_sst.o fv_grid_utils.o fv_diagnostics.o tp_core.o fv_mapz.o fv_mp_mod.o fv_timing.o sim_nc_mod.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_nudge.F90
fv_regional_bc.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_regional_bc.F90 fv_arrays.o fv_diagnostics.o fv_grid_utils.o fv_mapz.o fv_mp_mod.o fv_fill.o fv_eta.o boundary.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_regional_bc.F90
fv_restart.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_restart.F90 fv_arrays.o fv_io.o fv_grid_utils.o fv_diagnostics.o init_hydro.o test_cases.o fv_mp_mod.o fv_surf_map.o external_ic.o fv_eta.o boundary.o fv_timing.o fv_treat_da_inc.o coarse_grained_restart_files.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_restart.F90
fv_sg.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_sg.F90 gfdl_cloud_microphys.o fv_mp_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_sg.F90
fv_surf_map.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_surf_map.F90 fv_grid_utils.o fv_mp_mod.o fv_timing.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_surf_map.F90
fv_timing.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_timing.F90 fv_mp_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_timing.F90
fv_tracer2d.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_tracer2d.F90 tp_core.o fv_mp_mod.o fv_timing.o boundary.o fv_regional_bc.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_tracer2d.F90
fv_treat_da_inc.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_treat_da_inc.F90 fv_arrays.o fv_grid_utils.o fv_mp_mod.o sim_nc_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_treat_da_inc.F90
fv_update_phys.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_update_phys.F90 fv_mp_mod.o fv_arrays.o boundary.o fv_eta.o fv_timing.o fv_diagnostics.o fv_mapz.o fv_climate_nudge.o fv_ada_nudge.o fv_nudge.o fv_grid_utils.o fv_nesting.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_update_phys.F90
gfdl_cloud_microphys.o: $(SRCROOT)GFDL_atmos_cubed_sphere/driver/SHiELD/gfdl_cloud_microphys.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/driver/SHiELD/gfdl_cloud_microphys.F90
gfdl_mp.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/gfdl_mp.F90 fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/model/gfdl_mp.F90
init_hydro.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/init_hydro.F90 fv_grid_utils.o fv_mp_mod.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/init_hydro.F90
nh_core.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/nh_core.F90 tp_core.o nh_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/model/nh_core.F90
nh_utils.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/nh_utils.F90 tp_core.o sw_core.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/model/nh_utils.F90
read_climate_nudge_data.o: $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/read_climate_nudge_data.F90 $(SRCROOT)FMS/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/read_climate_nudge_data.F90
sim_nc_mod.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/sim_nc_mod.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/sim_nc_mod.F90
sorted_index.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/sorted_index.F90 fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/sorted_index.F90
sw_core.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/sw_core.F90 tp_core.o fv_mp_mod.o fv_arrays.o a2b_edge.o test_cases.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/model/sw_core.F90
test_cases.o: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/test_cases.F90 $(SRCROOT)FMS/include/file_version.h init_hydro.o fv_mp_mod.o fv_grid_utils.o fv_surf_map.o fv_grid_tools.o fv_eta.o fv_sg.o fv_diagnostics.o fv_arrays.o gfdl_cloud_microphys.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)GFDL_atmos_cubed_sphere/tools/test_cases.F90
tp_core.o: $(SRCROOT)GFDL_atmos_cubed_sphere/model/tp_core.F90 fv_grid_utils.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)GFDL_atmos_cubed_sphere/model/tp_core.F90
./test_cases.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/test_cases.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/test_cases.F90 .
./gfdl_cloud_microphys.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/driver/SHiELD/gfdl_cloud_microphys.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/driver/SHiELD/gfdl_cloud_microphys.F90 .
./fv_climate_nudge.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_climate_nudge.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_climate_nudge.F90 .
./tp_core.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/tp_core.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/tp_core.F90 .
./coarse_grained_diagnostics.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_grained_diagnostics.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_grained_diagnostics.F90 .
./fv_fill.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_fill.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_fill.F90 .
./fv_update_phys.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_update_phys.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_update_phys.F90 .
./sim_nc_mod.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/sim_nc_mod.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/sim_nc_mod.F90 .
./fv_tracer2d.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_tracer2d.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_tracer2d.F90 .
./fv_control.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_control.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_control.F90 .
./fv_diag_column.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_diag_column.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_diag_column.F90 .
./fv_mapz.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_mapz.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_mapz.F90 .
./fv_diagnostics.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_diagnostics.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_diagnostics.F90 .
./file_version.h: $(SRCROOT)FMS/include/file_version.h
	cp $(SRCROOT)FMS/include/file_version.h .
./gfdl_mp.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/gfdl_mp.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/gfdl_mp.F90 .
./read_climate_nudge_data.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/read_climate_nudge_data.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/read_climate_nudge_data.F90 .
./fv_restart.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_restart.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_restart.F90 .
./fv_treat_da_inc.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_treat_da_inc.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_treat_da_inc.F90 .
./fv_io.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_io.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_io.F90 .
./fv_timing.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_timing.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_timing.F90 .
./fv_sg.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_sg.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_sg.F90 .
./cloud_diagnosis.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/driver/SHiELD/cloud_diagnosis.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/driver/SHiELD/cloud_diagnosis.F90 .
./fv_nggps_diag.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_nggps_diag.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_nggps_diag.F90 .
./boundary.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/boundary.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/boundary.F90 .
./fv_cmip_diag.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_cmip_diag.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_cmip_diag.F90 .
./atmosphere.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/driver/GFDL/atmosphere.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/driver/GFDL/atmosphere.F90 .
./init_hydro.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/init_hydro.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/init_hydro.F90 .
./fv_eta.h: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_eta.h
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_eta.h .
./dyn_core.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/dyn_core.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/dyn_core.F90 .
./sorted_index.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/sorted_index.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/sorted_index.F90 .
./fv_grid_tools.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_grid_tools.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_grid_tools.F90 .
./fv_cmp.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_cmp.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_cmp.F90 .
./fv_arrays.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_arrays.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_arrays.F90 .
./fv_eta.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_eta.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_eta.F90 .
./fv_diagnostics.h: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_diagnostics.h
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_diagnostics.h .
./fv_nudge.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_nudge.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_nudge.F90 .
./nh_core.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/nh_core.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/nh_core.F90 .
./external_sst.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/external_sst.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/external_sst.F90 .
./nh_utils.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/nh_utils.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/nh_utils.F90 .
./fv_surf_map.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_surf_map.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_surf_map.F90 .
./fms_platform.h: $(SRCROOT)FMS/include/fms_platform.h
	cp $(SRCROOT)FMS/include/fms_platform.h .
./coarse_grained_restart_files.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_grained_restart_files.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_grained_restart_files.F90 .
./fv_mp_mod.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_mp_mod.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_mp_mod.F90 .
./sw_core.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/sw_core.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/sw_core.F90 .
./fv_ada_nudge.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_ada_nudge.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_ada_nudge.F90 .
./coarse_graining.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_graining.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_graining.F90 .
./fv_dynamics.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_dynamics.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_dynamics.F90 .
./a2b_edge.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/a2b_edge.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/a2b_edge.F90 .
./fv_regional_bc.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_regional_bc.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_regional_bc.F90 .
./fv_grid_utils.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_grid_utils.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_grid_utils.F90 .
./fv_nesting.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_nesting.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_nesting.F90 .
./external_ic.F90: $(SRCROOT)GFDL_atmos_cubed_sphere/tools/external_ic.F90
	cp $(SRCROOT)GFDL_atmos_cubed_sphere/tools/external_ic.F90 .
SRC = $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_grid_tools.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/driver/GFDL/atmosphere.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_mapz.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_nggps_diag.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_tracer2d.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/tp_core.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/init_hydro.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_nudge.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/boundary.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/nh_core.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_regional_bc.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_grained_restart_files.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_diagnostics.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_grained_diagnostics.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_control.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_timing.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_sg.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_treat_da_inc.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/a2b_edge.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_dynamics.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/gfdl_mp.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_restart.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/external_sst.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_arrays.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_fill.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_update_phys.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_grid_utils.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_cmp.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/sim_nc_mod.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_diag_column.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/sorted_index.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/external_ic.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_io.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_eta.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_ada_nudge.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/driver/SHiELD/gfdl_cloud_microphys.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_graining.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/nh_utils.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_cmip_diag.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/test_cases.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/read_climate_nudge_data.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_climate_nudge.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/dyn_core.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/sw_core.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_mp_mod.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_surf_map.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_nesting.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/driver/SHiELD/cloud_diagnosis.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_diagnostics.h $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_eta.h
OBJ = fv_grid_tools.o atmosphere.o fv_mapz.o fv_nggps_diag.o fv_tracer2d.o tp_core.o init_hydro.o fv_nudge.o boundary.o nh_core.o fv_regional_bc.o coarse_grained_restart_files.o fv_diagnostics.o coarse_grained_diagnostics.o fv_control.o fv_timing.o fv_sg.o fv_treat_da_inc.o a2b_edge.o fv_dynamics.o gfdl_mp.o fv_restart.o external_sst.o fv_arrays.o fv_fill.o fv_update_phys.o fv_grid_utils.o fv_cmp.o sim_nc_mod.o fv_diag_column.o sorted_index.o external_ic.o fv_io.o fv_eta.o fv_ada_nudge.o gfdl_cloud_microphys.o coarse_graining.o nh_utils.o fv_cmip_diag.o test_cases.o read_climate_nudge_data.o fv_climate_nudge.o dyn_core.o sw_core.o fv_mp_mod.o fv_surf_map.o fv_nesting.o cloud_diagnosis.o
OFF = $(SRCROOT)GFDL_atmos_cubed_sphere/tools/test_cases.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/driver/SHiELD/gfdl_cloud_microphys.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_climate_nudge.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/tp_core.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_grained_diagnostics.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_fill.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_update_phys.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/sim_nc_mod.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_tracer2d.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_control.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_diag_column.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_mapz.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_diagnostics.F90 $(SRCROOT)FMS/include/file_version.h $(SRCROOT)GFDL_atmos_cubed_sphere/model/gfdl_mp.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/read_climate_nudge_data.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_restart.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_treat_da_inc.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_io.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_timing.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_sg.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/driver/SHiELD/cloud_diagnosis.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_nggps_diag.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/boundary.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_cmip_diag.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/driver/GFDL/atmosphere.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/init_hydro.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_eta.h $(SRCROOT)GFDL_atmos_cubed_sphere/model/dyn_core.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/sorted_index.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_grid_tools.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_cmp.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_arrays.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_eta.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_diagnostics.h $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_nudge.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/nh_core.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/external_sst.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/nh_utils.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_surf_map.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_grained_restart_files.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/fv_mp_mod.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/sw_core.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/GFDL_tools/fv_ada_nudge.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/coarse_graining.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_dynamics.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/a2b_edge.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_regional_bc.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_grid_utils.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/model/fv_nesting.F90 $(SRCROOT)GFDL_atmos_cubed_sphere/tools/external_ic.F90
clean: neat
	-rm -f .libatmos_cubed_sphere.a.cppdefs $(OBJ) libatmos_cubed_sphere.a

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libatmos_cubed_sphere.a: $(OBJ)
	$(AR) $(ARFLAGS) libatmos_cubed_sphere.a $(OBJ)

