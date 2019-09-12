# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/land_lad2.
#
# Ed Hartnett 2019-07-30

# Run make targets in each of these subdirs
SUBDIRS = transitions canopy_air glacier . vegetation soil snow	\
topo_rough river lake shared predefined_tiles land_tracers

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libland_tile.la libland_constants.la	\
libland_model.la libland_data.la

# Each convenience library depends on its source.
libland_tile_la_SOURCES = land_tile.F90
libland_constants_la_SOURCES = land_constants.F90
libland_model_la_SOURCES = land_model.F90
libland_data_la_SOURCES = land_data.F90

# Each mod file depends on the .o file.
land_tile_mod.mod: land_tile.$(OBJEXT)
land_constants_mod.mod: land_constants.$(OBJEXT)
land_model_mod.mod: land_model.$(OBJEXT)
land_data_mod.mod: land_data.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
land_tile.$(OBJEXT): land_constants_mod.mod land_data_mod.mod
land_model.$(OBJEXT): land_tile_mod.mod land_constants_mod.mod land_data_mod.mod

# Mod files are built and then installed as headers.
MODFILES = land_tile_mod.mod land_constants_mod.mod land_model_mod.mod \
land_data_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = land_version.inc

CLEANFILES = *.mod
