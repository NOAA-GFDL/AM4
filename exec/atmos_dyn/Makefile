# Makefile created by mkmf 19.3.0

CPPDEFS = -DINTERNAL_FILE_NML -g -DCLIMATE_NUDGE -DSPMD 

OTHERFLAGS = -I$(BUILDROOT)atmos_phys -I$(BUILDROOT)fms -I$(SRCROOT)shared/include

include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libatmos_dyn.a
a2b_edge.o: $(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90 fv_grid_utils.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90
atmos_model.o: $(SRCROOT)atmos_drivers/coupled/atmos_model.F90 atmosphere.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_drivers/coupled/atmos_model.F90
atmosphere.o: $(SRCROOT)atmos_cubed_sphere/driver/GFDL/atmosphere.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)FMS/include/file_version.h fv_arrays.o fv_control.o fv_eta.o fv_io.o fv_dynamics.o fv_nesting.o fv_diagnostics.o fv_cmip_diag.o fv_restart.o fv_timing.o fv_mp_mod.o fv_sg.o fv_update_phys.o fv_climate_nudge.o fv_ada_nudge.o fv_nudge.o boundary.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_cubed_sphere/driver/GFDL/atmosphere.F90
boundary.o: $(SRCROOT)atmos_cubed_sphere/model/boundary.F90 fv_mp_mod.o fv_arrays.o fv_timing.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/boundary.F90
dyn_core.o: $(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90 fv_mp_mod.o sw_core.o a2b_edge.o nh_core.o tp_core.o fv_timing.o fv_diagnostics.o fv_update_phys.o fv_ada_nudge.o fv_nudge.o fv_arrays.o boundary.o test_cases.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90
external_ic.o: $(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90 $(SRCROOT)FMS/include/file_version.h external_sst.o fv_arrays.o fv_diagnostics.o fv_grid_utils.o fv_io.o fv_mapz.o fv_mp_mod.o fv_surf_map.o fv_timing.o init_hydro.o fv_fill.o fv_eta.o sim_nc_mod.o fv_nudge.o test_cases.o boundary.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90
external_sst.o: $(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90
fv_ada_nudge.o: $(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_ada_nudge.F90 $(SRCROOT)FMS/include/file_version.h external_sst.o fv_grid_utils.o fv_diagnostics.o tp_core.o fv_mapz.o fv_mp_mod.o fv_timing.o sim_nc_mod.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_ada_nudge.F90
fv_arrays.o: $(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90 $(SRCROOT)FMS/include/fms_platform.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90
fv_climate_nudge.o: $(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_climate_nudge.F90 $(SRCROOT)FMS/include/file_version.h read_climate_nudge_data.o fv_mapz.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_climate_nudge.F90
fv_cmip_diag.o: $(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_cmip_diag.F90 $(SRCROOT)FMS/include/file_version.h fv_mapz.o fv_arrays.o fv_diagnostics.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_cmip_diag.F90
fv_cmp.o: $(SRCROOT)atmos_cubed_sphere/model/fv_cmp.F90 fv_mp_mod.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_cmp.F90
fv_control.o: $(SRCROOT)atmos_cubed_sphere/model/fv_control.F90 $(SRCROOT)FMS/include/file_version.h fv_io.o fv_restart.o fv_arrays.o fv_grid_utils.o fv_eta.o fv_grid_tools.o fv_mp_mod.o test_cases.o fv_timing.o fv_diagnostics.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_cubed_sphere/model/fv_control.F90
fv_current_grid.o: $(SRCROOT)atmos_cubed_sphere/model/fv_current_grid.F90 $(SRCROOT)FMS/include/fms_platform.h fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_cubed_sphere/model/fv_current_grid.F90
fv_diagnostics.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90 $(SRCROOT)FMS/include/file_version.h fv_arrays.o fv_mapz.o fv_mp_mod.o fv_eta.o fv_grid_utils.o a2b_edge.o fv_surf_map.o fv_sg.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90
fv_dynamics.o: $(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90 dyn_core.o fv_mapz.o fv_tracer2d.o fv_grid_utils.o fv_fill.o fv_mp_mod.o fv_timing.o fv_diagnostics.o fv_sg.o fv_nesting.o boundary.o fv_arrays.o fv_nudge.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90
fv_eta.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90 fv_mp_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90
fv_fill.o: $(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90
fv_grid_tools.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90 fv_arrays.o fv_grid_utils.o fv_timing.o fv_mp_mod.o sorted_index.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90
fv_grid_utils.o: $(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90 $(SRCROOT)FMS/include/fms_platform.h external_sst.o fv_arrays.o fv_eta.o fv_mp_mod.o fv_timing.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90
fv_io.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90 external_sst.o fv_arrays.o fv_eta.o fv_mp_mod.o fv_mapz.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90
fv_mapz.o: $(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90 fv_grid_utils.o fv_fill.o fv_arrays.o fv_timing.o fv_mp_mod.o fv_cmp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90
fv_mp_mod.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90 fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90
fv_nesting.o: $(SRCROOT)atmos_cubed_sphere/model/fv_nesting.F90 fv_sg.o fv_restart.o boundary.o fv_mp_mod.o fv_arrays.o fv_grid_utils.o init_hydro.o fv_mapz.o fv_timing.o fv_diagnostics.o sw_core.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_nesting.F90
fv_nggps_diag.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_nggps_diag.F90 fv_diagnostics.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_nggps_diag.F90
fv_nudge.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90 $(SRCROOT)FMS/include/file_version.h external_sst.o fv_grid_utils.o fv_diagnostics.o tp_core.o fv_mapz.o fv_mp_mod.o fv_timing.o sim_nc_mod.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90
fv_restart.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90 fv_arrays.o fv_io.o fv_grid_utils.o fv_diagnostics.o init_hydro.o test_cases.o fv_mp_mod.o fv_surf_map.o external_ic.o fv_eta.o boundary.o fv_timing.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90
fv_sg.o: $(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90 fv_mp_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90
fv_surf_map.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90 fv_grid_utils.o fv_mp_mod.o fv_timing.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90
fv_timing.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90 fv_mp_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90
fv_tracer2d.o: $(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90 tp_core.o fv_mp_mod.o fv_timing.o boundary.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90
fv_update_phys.o: $(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90 fv_mp_mod.o fv_arrays.o boundary.o fv_eta.o fv_timing.o fv_diagnostics.o fv_mapz.o fv_climate_nudge.o fv_ada_nudge.o fv_nudge.o fv_grid_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90
init_hydro.o: $(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90 fv_grid_utils.o fv_mp_mod.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90
nh_core.o: $(SRCROOT)atmos_cubed_sphere/model_nh_null/nh_core.F90 fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model_nh_null/nh_core.F90
read_climate_nudge_data.o: $(SRCROOT)atmos_cubed_sphere/GFDL_tools/read_climate_nudge_data.F90 $(SRCROOT)FMS/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_cubed_sphere/GFDL_tools/read_climate_nudge_data.F90
sim_nc_mod.o: $(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90
sorted_index.o: $(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90 fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90
sw_core.o: $(SRCROOT)atmos_cubed_sphere/model/sw_core.F90 fv_mp_mod.o tp_core.o fv_arrays.o a2b_edge.o test_cases.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/sw_core.F90
test_cases.o: $(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90 init_hydro.o fv_mp_mod.o fv_grid_utils.o fv_surf_map.o fv_grid_tools.o fv_eta.o fv_sg.o fv_diagnostics.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90
tp_core.o: $(SRCROOT)atmos_cubed_sphere/model/tp_core.F90 fv_mp_mod.o fv_grid_utils.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/tp_core.F90
./fv_restart.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90 .
./fv_tracer2d.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90 .
./test_cases.F90: $(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90 .
./fv_nudge.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90 .
./file_version.h: $(SRCROOT)FMS/include/file_version.h
	cp $(SRCROOT)FMS/include/file_version.h .
./fv_current_grid.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_current_grid.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_current_grid.F90 .
./sim_nc_mod.F90: $(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90 .
./tp_core.F90: $(SRCROOT)atmos_cubed_sphere/model/tp_core.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/tp_core.F90 .
./atmos_model.F90: $(SRCROOT)atmos_drivers/coupled/atmos_model.F90
	cp $(SRCROOT)atmos_drivers/coupled/atmos_model.F90 .
./external_sst.F90: $(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90 .
./fv_sg.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90 .
./sw_core.F90: $(SRCROOT)atmos_cubed_sphere/model/sw_core.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/sw_core.F90 .
./dyn_core.F90: $(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90 .
./fv_cmp.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_cmp.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_cmp.F90 .
./fv_eta.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90 .
./nh_core.F90: $(SRCROOT)atmos_cubed_sphere/model_nh_null/nh_core.F90
	cp $(SRCROOT)atmos_cubed_sphere/model_nh_null/nh_core.F90 .
./fv_grid_utils.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90 .
./fv_nggps_diag.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_nggps_diag.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_nggps_diag.F90 .
./boundary.F90: $(SRCROOT)atmos_cubed_sphere/model/boundary.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/boundary.F90 .
./read_climate_nudge_data.F90: $(SRCROOT)atmos_cubed_sphere/GFDL_tools/read_climate_nudge_data.F90
	cp $(SRCROOT)atmos_cubed_sphere/GFDL_tools/read_climate_nudge_data.F90 .
./fv_diagnostics.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90 .
./fv_ada_nudge.F90: $(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_ada_nudge.F90
	cp $(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_ada_nudge.F90 .
./init_hydro.F90: $(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90 .
./fv_control.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_control.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_control.F90 .
./fv_timing.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90 .
./fv_surf_map.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90 .
./fv_arrays.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90 .
./a2b_edge.F90: $(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90 .
./fv_climate_nudge.F90: $(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_climate_nudge.F90
	cp $(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_climate_nudge.F90 .
./fv_grid_tools.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90 .
./fv_update_phys.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90 .
./sorted_index.F90: $(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90 .
./fv_mp_mod.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90 .
./fv_nesting.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_nesting.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_nesting.F90 .
./fv_io.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90 .
./atmosphere.F90: $(SRCROOT)atmos_cubed_sphere/driver/GFDL/atmosphere.F90
	cp $(SRCROOT)atmos_cubed_sphere/driver/GFDL/atmosphere.F90 .
./fv_cmip_diag.F90: $(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_cmip_diag.F90
	cp $(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_cmip_diag.F90 .
./external_ic.F90: $(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90 .
./fv_mapz.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90 .
./fv_dynamics.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90 .
./fms_platform.h: $(SRCROOT)FMS/include/fms_platform.h
	cp $(SRCROOT)FMS/include/fms_platform.h .
./fv_fill.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90 .
SRC = $(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90 $(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90 $(SRCROOT)atmos_cubed_sphere/model/tp_core.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_cmp.F90 $(SRCROOT)atmos_cubed_sphere/model/sw_core.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90 $(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90 $(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_cmip_diag.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90 $(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_nesting.F90 $(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90 $(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_current_grid.F90 $(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90 $(SRCROOT)atmos_drivers/coupled/atmos_model.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_control.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90 $(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_climate_nudge.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_nggps_diag.F90 $(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_ada_nudge.F90 $(SRCROOT)atmos_cubed_sphere/GFDL_tools/read_climate_nudge_data.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90 $(SRCROOT)atmos_cubed_sphere/model_nh_null/nh_core.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90 $(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90 $(SRCROOT)atmos_cubed_sphere/driver/GFDL/atmosphere.F90 $(SRCROOT)atmos_cubed_sphere/model/boundary.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90 $(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90 $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/include/fms_platform.h
OBJ = fv_mapz.o fv_restart.o a2b_edge.o tp_core.o fv_cmp.o sw_core.o fv_tracer2d.o dyn_core.o fv_cmip_diag.o fv_mp_mod.o fv_sg.o test_cases.o fv_fill.o fv_update_phys.o fv_nesting.o sim_nc_mod.o external_ic.o fv_arrays.o fv_current_grid.o init_hydro.o fv_io.o atmos_model.o fv_control.o fv_timing.o fv_climate_nudge.o fv_nggps_diag.o fv_ada_nudge.o read_climate_nudge_data.o fv_nudge.o nh_core.o fv_dynamics.o external_sst.o fv_diagnostics.o fv_eta.o atmosphere.o boundary.o fv_grid_tools.o fv_surf_map.o sorted_index.o fv_grid_utils.o
OFF = $(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90 $(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90 $(SRCROOT)FMS/include/file_version.h $(SRCROOT)atmos_cubed_sphere/model/fv_current_grid.F90 $(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90 $(SRCROOT)atmos_cubed_sphere/model/tp_core.F90 $(SRCROOT)atmos_drivers/coupled/atmos_model.F90 $(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90 $(SRCROOT)atmos_cubed_sphere/model/sw_core.F90 $(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_cmp.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90 $(SRCROOT)atmos_cubed_sphere/model_nh_null/nh_core.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_nggps_diag.F90 $(SRCROOT)atmos_cubed_sphere/model/boundary.F90 $(SRCROOT)atmos_cubed_sphere/GFDL_tools/read_climate_nudge_data.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90 $(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_ada_nudge.F90 $(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_control.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90 $(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90 $(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_climate_nudge.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90 $(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_nesting.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90 $(SRCROOT)atmos_cubed_sphere/driver/GFDL/atmosphere.F90 $(SRCROOT)atmos_cubed_sphere/GFDL_tools/fv_cmip_diag.F90 $(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90
clean: neat
	-rm -f .libatmos_dyn.a.cppdefs $(OBJ) libatmos_dyn.a

neat:
	-rm -f $(TMPFILES) *.mod

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libatmos_dyn.a: $(OBJ)
	$(AR) $(ARFLAGS) libatmos_dyn.a $(OBJ)

