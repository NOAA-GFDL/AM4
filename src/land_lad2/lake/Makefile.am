# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/land_lad2/lake.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = liblake.la liblake_tile.la

# Each convenience library depends on its source.
liblake_la_SOURCES = lake.F90
liblake_tile_la_SOURCES = lake_tile.F90

# Each mod file depends on the .o file.
lake_mod.mod: lake.$(OBJEXT)
lake_tile_mod.mod: lake_tile.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
lake.$(OBJEXT): lake_tile_mod.mod

# Mod files are built and then installed as headers.
MODFILES = lake_mod.mod lake_tile_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
