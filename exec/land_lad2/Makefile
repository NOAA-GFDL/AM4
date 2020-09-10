# Makefile created by mkmf 19.3.0
CPPDEFS = -DINTERNAL_FILE_NML

OTHERFLAGS = -I$(BUILDROOT)fms

include $(MK_TEMPLATE)

.DEFAULT:
	-echo $@ does not exist.
all: libland_lad2.a
cana_tile.o: $(SRCROOT)land_lad2/canopy_air/cana_tile.F90 land_tracers.o land_tile_diag_sel.o land_constants.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
canopy_air.o: $(SRCROOT)land_lad2/canopy_air/canopy_air.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc sphum.o land_constants.o land_tracers.o cana_tile.o land_tile.o land_data.o land_tile_io.o land_debug.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
glac_tile.o: $(SRCROOT)land_lad2/glacier/glac_tile.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_constants.o land_io.o land_tile_diag_sel.o land_data.o tiling_input_types.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)shared/include -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
glacier.o: $(SRCROOT)land_lad2/glacier/glacier.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc glac_tile.o land_tile.o land_tile_diag.o land_data.o land_tile_io.o land_debug.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
hillslope.o: $(SRCROOT)land_lad2/soil/hillslope.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_tile.o land_utils.o land_tile_diag.o land_data.o land_io.o land_tile_io.o nf_utils.o land_debug.o transitions.o vegn_harvesting.o hillslope_tile.o soil_tile.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
hillslope_hydrology.o: $(SRCROOT)land_lad2/soil/hillslope_hydrology.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc soil_tile.o land_tile.o land_data.o land_debug.o hillslope.o land_tile_diag.o soil_carbon.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
hillslope_tile.o: $(SRCROOT)land_lad2/soil/hillslope_tile.F90 land_tile_diag_sel.o soil_tile.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
lake.o: $(SRCROOT)land_lad2/lake/lake.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc lake_tile.o land_tile.o land_tile_diag.o land_data.o land_tile_io.o land_debug.o land_utils.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
lake_tile.o: $(SRCROOT)land_lad2/lake/lake_tile.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_constants.o land_data.o land_io.o land_tile_diag_sel.o tiling_input_types.o land_debug.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)shared/include -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
land_constants.o: $(SRCROOT)land_lad2/land_constants.F90
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
land_data.o: $(SRCROOT)land_lad2/land_data.F90 $(SRCROOT)land_lad2/shared/version_variable.inc
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
land_debug.o: $(SRCROOT)land_lad2/shared/land_debug.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_data.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
land_dust.o: $(SRCROOT)land_lad2/land_tracers/land_dust.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_constants.o cana_tile.o soil_tile.o snow_tile.o vegn_tile.o vegn_data.o land_tile.o land_tile_diag.o land_data.o land_io.o land_tracers.o land_debug.o table_printer.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
land_io.o: $(SRCROOT)land_lad2/shared/land_io.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_numerics.o nf_utils.o land_data.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
land_model.o: $(SRCROOT)land_lad2/land_model.F90 $(SRCROOT)land_lad2/shared/debug.inc $(SRCROOT)land_lad2/shared/version_variable.inc sphum.o land_constants.o land_tracers.o land_tracer_driver.o glacier.o lake.o soil.o soil_carbon.o snow.o vegetation.o cana_tile.o canopy_air.o river.o topo_rough.o soil_tile.o vegn_tile.o lake_tile.o glac_tile.o snow_tile.o land_numerics.o land_io.o land_tile.o land_data.o nf_utils.o land_tile_io.o land_tile_diag.o land_debug.o vegn_static_override.o transitions.o hillslope.o hillslope_hydrology.o vegn_data.o tiling_input.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
land_numerics.o: $(SRCROOT)land_lad2/shared/land_numerics.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_data.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
land_tile.o: $(SRCROOT)land_lad2/land_tile.F90 land_constants.o glac_tile.o lake_tile.o soil_tile.o hillslope_tile.o cana_tile.o vegn_tile.o snow_tile.o land_tile_diag_sel.o land_tile_diag_buff.o land_data.o tiling_input_types.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
land_tile_diag.o: $(SRCROOT)land_lad2/shared/land_tile_diag.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_tile_diag_sel.o land_tile.o vegn_cohort.o land_data.o land_debug.o land_tile_diag_buff.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
land_tile_diag_buff.o: $(SRCROOT)land_lad2/shared/land_tile_diag_buff.F90
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
land_tile_diag_sel.o: $(SRCROOT)land_lad2/shared/land_tile_diag_sel.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_data.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
land_tile_io.o: $(SRCROOT)land_lad2/shared/land_tile_io.F90 nf_utils.o land_io.o land_tile.o land_data.o land_utils.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
land_tracer_driver.o: $(SRCROOT)land_lad2/land_tracers/land_tracer_driver.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc table_printer.o land_constants.o land_debug.o land_data.o land_tracers.o land_tile.o land_tile_diag.o cana_tile.o vegn_data.o vegn_tile.o vegn_cohort.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
land_tracers.o: $(SRCROOT)land_lad2/land_tracers/land_tracers.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_data.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
land_utils.o: $(SRCROOT)land_lad2/shared/land_utils.F90 land_tile.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
nf_utils.o: $(SRCROOT)land_lad2/shared/nf_utils/nf_utils.F90 nfu.o nfc.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
nfc.o: $(SRCROOT)land_lad2/shared/nf_utils/nfc.F90 $(SRCROOT)land_lad2/shared/nf_utils/getput_compressed.inc nfu.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/shared/nf_utils
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
nfu.o: $(SRCROOT)land_lad2/shared/nf_utils/nfu.F90 $(SRCROOT)land_lad2/shared/nf_utils/getput.inc
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/shared/nf_utils
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
river.o: $(SRCROOT)land_lad2/river/river.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc river_type.o river_physics.o land_tile.o land_data.o lake_tile.o land_io.o table_printer.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
river_physics.o: $(SRCROOT)land_lad2/river/river_physics.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc river_type.o lake.o lake_tile.o land_debug.o land_data.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
river_type.o: $(SRCROOT)land_lad2/river/river_type.F90
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
snow.o: $(SRCROOT)land_lad2/snow/snow.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_constants.o snow_tile.o land_tile.o land_data.o land_tile_io.o land_debug.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
snow_tile.o: $(SRCROOT)land_lad2/snow/snow_tile.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_constants.o land_tile_diag_sel.o land_data.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)shared/include -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
soil.o: $(SRCROOT)land_lad2/soil/soil.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_constants.o soil_tile.o soil_carbon.o land_tile.o land_utils.o land_tile_diag.o land_data.o land_io.o land_tile_io.o vegn_data.o vegn_tile.o land_debug.o uptake.o hillslope.o hillslope_hydrology.o river.o land_numerics.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
soil_carbon.o: $(SRCROOT)land_lad2/soil/soil_carbon.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_constants.o vegn_data.o land_numerics.o land_debug.o land_data.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
soil_tile.o: $(SRCROOT)land_lad2/soil/soil_tile.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_constants.o land_tile_diag_sel.o soil_carbon.o land_data.o land_debug.o tiling_input_types.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)shared/include -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
sphum.o: $(SRCROOT)land_lad2/shared/sphum.F90 land_debug.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
table_printer.o: $(SRCROOT)land_lad2/shared/table_printer.F90
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
tiling_input.o: $(SRCROOT)land_lad2/predefined_tiles/tiling_input.F90 land_data.o vegn_cohort.o land_tile.o tiling_input_types.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
tiling_input_types.o: $(SRCROOT)land_lad2/predefined_tiles/tiling_input_types.F90
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
topo_rough.o: $(SRCROOT)land_lad2/topo_rough/topo_rough.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_data.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
transitions.o: $(SRCROOT)land_lad2/transitions/transitions.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc nfu.o vegn_data.o cana_tile.o snow_tile.o vegn_tile.o soil_tile.o land_tile.o land_tile_diag.o land_tile_io.o land_data.o vegn_harvesting.o land_debug.o land_numerics.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)shared/include -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
uptake.o: $(SRCROOT)land_lad2/soil/uptake.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc soil_tile.o land_debug.o land_data.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
vegetation.o: $(SRCROOT)land_lad2/vegetation/vegetation.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc sphum.o vegn_tile.o soil_tile.o land_constants.o land_tile.o land_tile_diag.o land_data.o land_io.o land_tile_io.o vegn_data.o vegn_cohort.o canopy_air.o vegn_cohort_io.o land_debug.o vegn_radiation.o vegn_photosynthesis.o vegn_static_override.o vegn_dynamics.o vegn_disturbance.o vegn_harvesting.o soil_carbon.o soil.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
vegn_cohort.o: $(SRCROOT)land_lad2/vegetation/vegn_cohort.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc land_constants.o vegn_data.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
vegn_cohort_io.o: $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.F90 $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.inc nf_utils.o land_io.o land_tile.o land_tile_io.o vegn_cohort.o land_data.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/vegetation
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
vegn_data.o: $(SRCROOT)land_lad2/vegetation/vegn_data.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_constants.o land_tile_diag_sel.o land_data.o table_printer.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
vegn_disturbance.o: $(SRCROOT)land_lad2/vegetation/vegn_disturbance.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_constants.o land_data.o vegn_data.o vegn_tile.o soil_tile.o vegn_cohort.o soil.o soil_carbon.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
vegn_dynamics.o: $(SRCROOT)land_lad2/vegetation/vegn_dynamics.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_constants.o land_tile_diag.o vegn_data.o vegn_tile.o soil_tile.o vegn_cohort.o soil_carbon.o soil.o land_debug.o land_data.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
vegn_harvesting.o: $(SRCROOT)land_lad2/vegetation/vegn_harvesting.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc vegn_data.o soil_tile.o vegn_tile.o vegn_cohort.o soil_carbon.o land_data.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
vegn_photosynthesis.o: $(SRCROOT)land_lad2/vegetation/vegn_photosynthesis.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc sphum.o land_constants.o land_debug.o vegn_data.o vegn_tile.o vegn_cohort.o land_data.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
vegn_radiation.o: $(SRCROOT)land_lad2/vegetation/vegn_radiation.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc land_constants.o vegn_data.o vegn_tile.o vegn_cohort.o snow_tile.o land_debug.o land_data.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
vegn_static_override.o: $(SRCROOT)land_lad2/vegetation/vegn_static_override.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc $(SRCROOT)land_lad2/vegetation/read_remap_cohort_data.inc $(SRCROOT)land_lad2/vegetation/read_remap_cohort_data_new.inc nf_utils.o land_data.o land_io.o land_numerics.o land_tile_io.o land_tile.o vegn_cohort.o vegn_cohort_io.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared -I$(SRCROOT)land_lad2/vegetation
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
vegn_tile.o: $(SRCROOT)land_lad2/vegetation/vegn_tile.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc land_constants.o land_io.o land_tile_diag_sel.o vegn_data.o vegn_cohort.o
	$(eval preproc := $(addsuffix .DO_NOT_MODIFY.f90,$(basename $(notdir $<))))
	cpp -nostdinc -C -v $(CPPDEFS) $(CPPFLAGS) $< > $(preproc) -I$(SRCROOT)land_lad2/predefined_tiles/../shared
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -o $@ 	$(preproc)
./uptake.F90: $(SRCROOT)land_lad2/soil/uptake.F90
	cp $(SRCROOT)land_lad2/soil/uptake.F90 .
./land_tile_diag.F90: $(SRCROOT)land_lad2/shared/land_tile_diag.F90
	cp $(SRCROOT)land_lad2/shared/land_tile_diag.F90 .
./soil.F90: $(SRCROOT)land_lad2/soil/soil.F90
	cp $(SRCROOT)land_lad2/soil/soil.F90 .
./vegetation.F90: $(SRCROOT)land_lad2/vegetation/vegetation.F90
	cp $(SRCROOT)land_lad2/vegetation/vegetation.F90 .
./fms_platform.h: $(SRCROOT)FMS/include/fms_platform.h
	cp $(SRCROOT)FMS/include/fms_platform.h .
./vegn_photosynthesis.F90: $(SRCROOT)land_lad2/vegetation/vegn_photosynthesis.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_photosynthesis.F90 .
./debug.inc: $(SRCROOT)land_lad2/shared/debug.inc
	cp $(SRCROOT)land_lad2/shared/debug.inc .
./hillslope.F90: $(SRCROOT)land_lad2/soil/hillslope.F90
	cp $(SRCROOT)land_lad2/soil/hillslope.F90 .
./debug.inc: $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc
	cp $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc .
./land_data.F90: $(SRCROOT)land_lad2/land_data.F90
	cp $(SRCROOT)land_lad2/land_data.F90 .
./tiling_input_types.F90: $(SRCROOT)land_lad2/predefined_tiles/tiling_input_types.F90
	cp $(SRCROOT)land_lad2/predefined_tiles/tiling_input_types.F90 .
./getput_compressed.inc: $(SRCROOT)land_lad2/shared/nf_utils/getput_compressed.inc
	cp $(SRCROOT)land_lad2/shared/nf_utils/getput_compressed.inc .
./snow.F90: $(SRCROOT)land_lad2/snow/snow.F90
	cp $(SRCROOT)land_lad2/snow/snow.F90 .
./land_tracer_driver.F90: $(SRCROOT)land_lad2/land_tracers/land_tracer_driver.F90
	cp $(SRCROOT)land_lad2/land_tracers/land_tracer_driver.F90 .
./land_tile_io.F90: $(SRCROOT)land_lad2/shared/land_tile_io.F90
	cp $(SRCROOT)land_lad2/shared/land_tile_io.F90 .
./vegn_cohort_io.inc: $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.inc
	cp $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.inc .
./vegn_data.F90: $(SRCROOT)land_lad2/vegetation/vegn_data.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_data.F90 .
./lake.F90: $(SRCROOT)land_lad2/lake/lake.F90
	cp $(SRCROOT)land_lad2/lake/lake.F90 .
./vegn_dynamics.F90: $(SRCROOT)land_lad2/vegetation/vegn_dynamics.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_dynamics.F90 .
./land_io.F90: $(SRCROOT)land_lad2/shared/land_io.F90
	cp $(SRCROOT)land_lad2/shared/land_io.F90 .
./land_constants.F90: $(SRCROOT)land_lad2/land_constants.F90
	cp $(SRCROOT)land_lad2/land_constants.F90 .
./land_tile.F90: $(SRCROOT)land_lad2/land_tile.F90
	cp $(SRCROOT)land_lad2/land_tile.F90 .
./land_numerics.F90: $(SRCROOT)land_lad2/shared/land_numerics.F90
	cp $(SRCROOT)land_lad2/shared/land_numerics.F90 .
./tiling_input.F90: $(SRCROOT)land_lad2/predefined_tiles/tiling_input.F90
	cp $(SRCROOT)land_lad2/predefined_tiles/tiling_input.F90 .
./canopy_air.F90: $(SRCROOT)land_lad2/canopy_air/canopy_air.F90
	cp $(SRCROOT)land_lad2/canopy_air/canopy_air.F90 .
./hillslope_tile.F90: $(SRCROOT)land_lad2/soil/hillslope_tile.F90
	cp $(SRCROOT)land_lad2/soil/hillslope_tile.F90 .
./soil_carbon.F90: $(SRCROOT)land_lad2/soil/soil_carbon.F90
	cp $(SRCROOT)land_lad2/soil/soil_carbon.F90 .
./topo_rough.F90: $(SRCROOT)land_lad2/topo_rough/topo_rough.F90
	cp $(SRCROOT)land_lad2/topo_rough/topo_rough.F90 .
./hillslope_hydrology.F90: $(SRCROOT)land_lad2/soil/hillslope_hydrology.F90
	cp $(SRCROOT)land_lad2/soil/hillslope_hydrology.F90 .
./nfu.F90: $(SRCROOT)land_lad2/shared/nf_utils/nfu.F90
	cp $(SRCROOT)land_lad2/shared/nf_utils/nfu.F90 .
./read_remap_cohort_data.inc: $(SRCROOT)land_lad2/vegetation/read_remap_cohort_data.inc
	cp $(SRCROOT)land_lad2/vegetation/read_remap_cohort_data.inc .
./cana_tile.F90: $(SRCROOT)land_lad2/canopy_air/cana_tile.F90
	cp $(SRCROOT)land_lad2/canopy_air/cana_tile.F90 .
./vegn_disturbance.F90: $(SRCROOT)land_lad2/vegetation/vegn_disturbance.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_disturbance.F90 .
./vegn_tile.F90: $(SRCROOT)land_lad2/vegetation/vegn_tile.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_tile.F90 .
./river_type.F90: $(SRCROOT)land_lad2/river/river_type.F90
	cp $(SRCROOT)land_lad2/river/river_type.F90 .
./soil_tile.F90: $(SRCROOT)land_lad2/soil/soil_tile.F90
	cp $(SRCROOT)land_lad2/soil/soil_tile.F90 .
./land_tile_diag_buff.F90: $(SRCROOT)land_lad2/shared/land_tile_diag_buff.F90
	cp $(SRCROOT)land_lad2/shared/land_tile_diag_buff.F90 .
./river_physics.F90: $(SRCROOT)land_lad2/river/river_physics.F90
	cp $(SRCROOT)land_lad2/river/river_physics.F90 .
./snow_tile.F90: $(SRCROOT)land_lad2/snow/snow_tile.F90
	cp $(SRCROOT)land_lad2/snow/snow_tile.F90 .
./glac_tile.F90: $(SRCROOT)land_lad2/glacier/glac_tile.F90
	cp $(SRCROOT)land_lad2/glacier/glac_tile.F90 .
./vegn_harvesting.F90: $(SRCROOT)land_lad2/vegetation/vegn_harvesting.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_harvesting.F90 .
./version_variable.inc: $(SRCROOT)land_lad2/shared/version_variable.inc
	cp $(SRCROOT)land_lad2/shared/version_variable.inc .
./land_dust.F90: $(SRCROOT)land_lad2/land_tracers/land_dust.F90
	cp $(SRCROOT)land_lad2/land_tracers/land_dust.F90 .
./table_printer.F90: $(SRCROOT)land_lad2/shared/table_printer.F90
	cp $(SRCROOT)land_lad2/shared/table_printer.F90 .
./sphum.F90: $(SRCROOT)land_lad2/shared/sphum.F90
	cp $(SRCROOT)land_lad2/shared/sphum.F90 .
./glacier.F90: $(SRCROOT)land_lad2/glacier/glacier.F90
	cp $(SRCROOT)land_lad2/glacier/glacier.F90 .
./vegn_radiation.F90: $(SRCROOT)land_lad2/vegetation/vegn_radiation.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_radiation.F90 .
./land_tracers.F90: $(SRCROOT)land_lad2/land_tracers/land_tracers.F90
	cp $(SRCROOT)land_lad2/land_tracers/land_tracers.F90 .
./vegn_cohort_io.F90: $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.F90 .
./lake_tile.F90: $(SRCROOT)land_lad2/lake/lake_tile.F90
	cp $(SRCROOT)land_lad2/lake/lake_tile.F90 .
./river.F90: $(SRCROOT)land_lad2/river/river.F90
	cp $(SRCROOT)land_lad2/river/river.F90 .
./getput.inc: $(SRCROOT)land_lad2/shared/nf_utils/getput.inc
	cp $(SRCROOT)land_lad2/shared/nf_utils/getput.inc .
./vegn_static_override.F90: $(SRCROOT)land_lad2/vegetation/vegn_static_override.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_static_override.F90 .
./transitions.F90: $(SRCROOT)land_lad2/transitions/transitions.F90
	cp $(SRCROOT)land_lad2/transitions/transitions.F90 .
./read_remap_cohort_data_new.inc: $(SRCROOT)land_lad2/vegetation/read_remap_cohort_data_new.inc
	cp $(SRCROOT)land_lad2/vegetation/read_remap_cohort_data_new.inc .
./land_tile_diag_sel.F90: $(SRCROOT)land_lad2/shared/land_tile_diag_sel.F90
	cp $(SRCROOT)land_lad2/shared/land_tile_diag_sel.F90 .
./vegn_cohort.F90: $(SRCROOT)land_lad2/vegetation/vegn_cohort.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_cohort.F90 .
./land_debug.F90: $(SRCROOT)land_lad2/shared/land_debug.F90
	cp $(SRCROOT)land_lad2/shared/land_debug.F90 .
./land_utils.F90: $(SRCROOT)land_lad2/shared/land_utils.F90
	cp $(SRCROOT)land_lad2/shared/land_utils.F90 .
./nf_utils.F90: $(SRCROOT)land_lad2/shared/nf_utils/nf_utils.F90
	cp $(SRCROOT)land_lad2/shared/nf_utils/nf_utils.F90 .
./nfc.F90: $(SRCROOT)land_lad2/shared/nf_utils/nfc.F90
	cp $(SRCROOT)land_lad2/shared/nf_utils/nfc.F90 .
./version_variable.inc: $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc
	cp $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc .
./land_model.F90: $(SRCROOT)land_lad2/land_model.F90
	cp $(SRCROOT)land_lad2/land_model.F90 .
SRC = $(SRCROOT)land_lad2/vegetation/vegn_tile.F90 $(SRCROOT)land_lad2/soil/soil.F90 $(SRCROOT)land_lad2/land_tracers/land_tracers.F90 $(SRCROOT)land_lad2/canopy_air/cana_tile.F90 $(SRCROOT)land_lad2/soil/hillslope_hydrology.F90 $(SRCROOT)land_lad2/soil/soil_carbon.F90 $(SRCROOT)land_lad2/land_data.F90 $(SRCROOT)land_lad2/vegetation/vegn_cohort.F90 $(SRCROOT)land_lad2/shared/land_numerics.F90 $(SRCROOT)land_lad2/glacier/glacier.F90 $(SRCROOT)land_lad2/shared/nf_utils/nfc.F90 $(SRCROOT)land_lad2/lake/lake.F90 $(SRCROOT)land_lad2/soil/hillslope_tile.F90 $(SRCROOT)land_lad2/shared/table_printer.F90 $(SRCROOT)land_lad2/soil/hillslope.F90 $(SRCROOT)land_lad2/shared/land_tile_diag_buff.F90 $(SRCROOT)land_lad2/snow/snow.F90 $(SRCROOT)land_lad2/shared/nf_utils/nf_utils.F90 $(SRCROOT)land_lad2/canopy_air/canopy_air.F90 $(SRCROOT)land_lad2/shared/land_tile_diag.F90 $(SRCROOT)land_lad2/river/river.F90 $(SRCROOT)land_lad2/vegetation/vegn_photosynthesis.F90 $(SRCROOT)land_lad2/vegetation/vegn_static_override.F90 $(SRCROOT)land_lad2/shared/land_io.F90 $(SRCROOT)land_lad2/lake/lake_tile.F90 $(SRCROOT)land_lad2/land_tile.F90 $(SRCROOT)land_lad2/vegetation/vegn_disturbance.F90 $(SRCROOT)land_lad2/vegetation/vegetation.F90 $(SRCROOT)land_lad2/topo_rough/topo_rough.F90 $(SRCROOT)land_lad2/shared/sphum.F90 $(SRCROOT)land_lad2/land_tracers/land_dust.F90 $(SRCROOT)land_lad2/glacier/glac_tile.F90 $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.F90 $(SRCROOT)land_lad2/predefined_tiles/tiling_input.F90 $(SRCROOT)land_lad2/soil/soil_tile.F90 $(SRCROOT)land_lad2/shared/land_tile_diag_sel.F90 $(SRCROOT)land_lad2/shared/land_tile_io.F90 $(SRCROOT)land_lad2/river/river_physics.F90 $(SRCROOT)land_lad2/land_tracers/land_tracer_driver.F90 $(SRCROOT)land_lad2/vegetation/vegn_harvesting.F90 $(SRCROOT)land_lad2/vegetation/vegn_dynamics.F90 $(SRCROOT)land_lad2/land_constants.F90 $(SRCROOT)land_lad2/vegetation/vegn_data.F90 $(SRCROOT)land_lad2/river/river_type.F90 $(SRCROOT)land_lad2/shared/land_debug.F90 $(SRCROOT)land_lad2/transitions/transitions.F90 $(SRCROOT)land_lad2/shared/nf_utils/nfu.F90 $(SRCROOT)land_lad2/predefined_tiles/tiling_input_types.F90 $(SRCROOT)land_lad2/land_model.F90 $(SRCROOT)land_lad2/snow/snow_tile.F90 $(SRCROOT)land_lad2/shared/land_utils.F90 $(SRCROOT)land_lad2/vegetation/vegn_radiation.F90 $(SRCROOT)land_lad2/soil/uptake.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc $(SRCROOT)land_lad2/shared/debug.inc $(SRCROOT)land_lad2/shared/nf_utils/getput_compressed.inc $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)land_lad2/vegetation/read_remap_cohort_data_new.inc $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.inc $(SRCROOT)land_lad2/vegetation/read_remap_cohort_data.inc $(SRCROOT)land_lad2/shared/version_variable.inc $(SRCROOT)land_lad2/shared/nf_utils/getput.inc $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc
OBJ = vegn_tile.o soil.o land_tracers.o cana_tile.o hillslope_hydrology.o soil_carbon.o land_data.o vegn_cohort.o land_numerics.o glacier.o nfc.o lake.o hillslope_tile.o table_printer.o hillslope.o land_tile_diag_buff.o snow.o nf_utils.o canopy_air.o land_tile_diag.o river.o vegn_photosynthesis.o vegn_static_override.o land_io.o lake_tile.o land_tile.o vegn_disturbance.o vegetation.o topo_rough.o sphum.o land_dust.o glac_tile.o vegn_cohort_io.o tiling_input.o soil_tile.o land_tile_diag_sel.o land_tile_io.o river_physics.o land_tracer_driver.o vegn_harvesting.o vegn_dynamics.o land_constants.o vegn_data.o river_type.o land_debug.o transitions.o nfu.o tiling_input_types.o land_model.o snow_tile.o land_utils.o vegn_radiation.o uptake.o
OFF = $(SRCROOT)land_lad2/soil/uptake.F90 $(SRCROOT)land_lad2/shared/land_tile_diag.F90 $(SRCROOT)land_lad2/soil/soil.F90 $(SRCROOT)land_lad2/vegetation/vegetation.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)land_lad2/vegetation/vegn_photosynthesis.F90 $(SRCROOT)land_lad2/shared/debug.inc $(SRCROOT)land_lad2/soil/hillslope.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/debug.inc $(SRCROOT)land_lad2/land_data.F90 $(SRCROOT)land_lad2/predefined_tiles/tiling_input_types.F90 $(SRCROOT)land_lad2/shared/nf_utils/getput_compressed.inc $(SRCROOT)land_lad2/snow/snow.F90 $(SRCROOT)land_lad2/land_tracers/land_tracer_driver.F90 $(SRCROOT)land_lad2/shared/land_tile_io.F90 $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.inc $(SRCROOT)land_lad2/vegetation/vegn_data.F90 $(SRCROOT)land_lad2/lake/lake.F90 $(SRCROOT)land_lad2/vegetation/vegn_dynamics.F90 $(SRCROOT)land_lad2/shared/land_io.F90 $(SRCROOT)land_lad2/land_constants.F90 $(SRCROOT)land_lad2/land_tile.F90 $(SRCROOT)land_lad2/shared/land_numerics.F90 $(SRCROOT)land_lad2/predefined_tiles/tiling_input.F90 $(SRCROOT)land_lad2/canopy_air/canopy_air.F90 $(SRCROOT)land_lad2/soil/hillslope_tile.F90 $(SRCROOT)land_lad2/soil/soil_carbon.F90 $(SRCROOT)land_lad2/topo_rough/topo_rough.F90 $(SRCROOT)land_lad2/soil/hillslope_hydrology.F90 $(SRCROOT)land_lad2/shared/nf_utils/nfu.F90 $(SRCROOT)land_lad2/vegetation/read_remap_cohort_data.inc $(SRCROOT)land_lad2/canopy_air/cana_tile.F90 $(SRCROOT)land_lad2/vegetation/vegn_disturbance.F90 $(SRCROOT)land_lad2/vegetation/vegn_tile.F90 $(SRCROOT)land_lad2/river/river_type.F90 $(SRCROOT)land_lad2/soil/soil_tile.F90 $(SRCROOT)land_lad2/shared/land_tile_diag_buff.F90 $(SRCROOT)land_lad2/river/river_physics.F90 $(SRCROOT)land_lad2/snow/snow_tile.F90 $(SRCROOT)land_lad2/glacier/glac_tile.F90 $(SRCROOT)land_lad2/vegetation/vegn_harvesting.F90 $(SRCROOT)land_lad2/shared/version_variable.inc $(SRCROOT)land_lad2/land_tracers/land_dust.F90 $(SRCROOT)land_lad2/shared/table_printer.F90 $(SRCROOT)land_lad2/shared/sphum.F90 $(SRCROOT)land_lad2/glacier/glacier.F90 $(SRCROOT)land_lad2/vegetation/vegn_radiation.F90 $(SRCROOT)land_lad2/land_tracers/land_tracers.F90 $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.F90 $(SRCROOT)land_lad2/lake/lake_tile.F90 $(SRCROOT)land_lad2/river/river.F90 $(SRCROOT)land_lad2/shared/nf_utils/getput.inc $(SRCROOT)land_lad2/vegetation/vegn_static_override.F90 $(SRCROOT)land_lad2/transitions/transitions.F90 $(SRCROOT)land_lad2/vegetation/read_remap_cohort_data_new.inc $(SRCROOT)land_lad2/shared/land_tile_diag_sel.F90 $(SRCROOT)land_lad2/vegetation/vegn_cohort.F90 $(SRCROOT)land_lad2/shared/land_debug.F90 $(SRCROOT)land_lad2/shared/land_utils.F90 $(SRCROOT)land_lad2/shared/nf_utils/nf_utils.F90 $(SRCROOT)land_lad2/shared/nf_utils/nfc.F90 $(SRCROOT)land_lad2/predefined_tiles/../shared/version_variable.inc $(SRCROOT)land_lad2/land_model.F90
clean: neat
	$(RM) $(OBJ) libland_lad2.a
	$(RM) *.DO_NOT_MODIFY.f90
	$(RM) *.mod

neat:
	$(RM) $(TMPFILES)

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libland_lad2.a: $(OBJ)
	$(AR) $(ARFLAGS) libland_lad2.a $(OBJ)

