# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/land_lad2/canopy_air.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libcanopy_air.la libcana_tile.la

# Each convenience library depends on its source.
libcanopy_air_la_SOURCES = canopy_air.F90
libcana_tile_la_SOURCES = cana_tile.F90

# Each mod file depends on the .o file.
canopy_air_mod.mod: canopy_air.$(OBJEXT)
cana_tile_mod.mod: cana_tile.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
canopy_air.$(OBJEXT): cana_tile_mod.mod

# Mod files are built and then installed as headers.
MODFILES = canopy_air_mod.mod cana_tile_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
