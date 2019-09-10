# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/land_lad2/vegetation.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libvegn_photosynthesis.la libvegn_data.la libvegn_dynamics.la \
libvegn_tile.la libvegn_cohort.la libvegn_static_override.la \
libvegn_harvesting.la libvegetation.la libvegn_radiation.la \
libvegn_disturbance.la libvegn_cohort_io.la

# Each convenience library depends on its source.
libvegn_photosynthesis_la_SOURCES = vegn_photosynthesis.F90
libvegn_data_la_SOURCES = vegn_data.F90
libvegn_dynamics_la_SOURCES = vegn_dynamics.F90
libvegn_tile_la_SOURCES = vegn_tile.F90
libvegn_cohort_la_SOURCES = vegn_cohort.F90
libvegn_static_override_la_SOURCES = vegn_static_override.F90
libvegn_harvesting_la_SOURCES = vegn_harvesting.F90
libvegetation_la_SOURCES = vegetation.F90
libvegn_radiation_la_SOURCES = vegn_radiation.F90
libvegn_disturbance_la_SOURCES = vegn_disturbance.F90
libvegn_cohort_io_la_SOURCES = vegn_cohort_io.F90

# Each mod file depends on the .o file.
vegn_photosynthesis_mod.mod: vegn_photosynthesis.$(OBJEXT)
vegn_data_mod.mod: vegn_data.$(OBJEXT)
vegn_dynamics_mod.mod: vegn_dynamics.$(OBJEXT)
vegn_tile_mod.mod: vegn_tile.$(OBJEXT)
vegn_cohort_mod.mod: vegn_cohort.$(OBJEXT)
vegn_static_override_mod.mod: vegn_static_override.$(OBJEXT)
vegn_harvesting_mod.mod: vegn_harvesting.$(OBJEXT)
vegetation_mod.mod: vegetation.$(OBJEXT)
vegn_radiation_mod.mod: vegn_radiation.$(OBJEXT)
vegn_disturbance_mod.mod: vegn_disturbance.$(OBJEXT)
vegn_cohort_io_mod.mod: vegn_cohort_io.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
vegn_photosynthesis.$(OBJEXT): vegn_data_mod.mod vegn_tile_mod.mod vegn_cohort_mod.mod
vegn_dynamics.$(OBJEXT): vegn_data_mod.mod vegn_tile_mod.mod vegn_cohort_mod.mod
vegn_tile.$(OBJEXT): vegn_data_mod.mod vegn_cohort_mod.mod vegetation_mod.mod
vegn_cohort.$(OBJEXT): vegn_data_mod.mod vegn_tile_mod.mod
vegn_static_override.$(OBJEXT): vegn_cohort_mod.mod
vegn_harvesting.$(OBJEXT): vegn_data_mod.mod vegn_tile_mod.mod vegn_cohort_mod.mod
vegetation.$(OBJEXT): vegn_photosynthesis_mod.mod vegn_data_mod.mod vegn_dynamics_mod.mod vegn_tile_mod.mod vegn_cohort_mod.mod vegn_harvesting_mod.mod vegn_radiation_mod.mod vegn_disturbance_mod.mod
vegn_radiation.$(OBJEXT): vegn_data_mod.mod vegn_tile_mod.mod vegn_cohort_mod.mod
vegn_disturbance.$(OBJEXT): vegn_data_mod.mod vegn_tile_mod.mod vegn_cohort_mod.mod
vegn_cohort_io.$(OBJEXT): vegn_cohort_mod.mod

# Mod files are built and then installed as headers.
MODFILES = vegn_photosynthesis_mod.mod vegn_data_mod.mod vegn_dynamics_mod.mod \
vegn_tile_mod.mod vegn_cohort_mod.mod vegn_static_override_mod.mod \
vegn_harvesting_mod.mod vegetation_mod.mod vegn_radiation_mod.mod \
vegn_disturbance_mod.mod vegn_cohort_io_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = vegn_cohort_io.inc read_remap_cohort_data_new.inc \
read_remap_cohort_data.inc

CLEANFILES = *.mod
