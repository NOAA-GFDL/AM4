# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/land_lad2/snow.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libsnow.la libsnow_tile.la

# Each convenience library depends on its source.
libsnow_la_SOURCES = snow.F90
libsnow_tile_la_SOURCES = snow_tile.F90

# Each mod file depends on the .o file.
snow_mod.mod: snow.$(OBJEXT)
snow_tile_mod.mod: snow_tile.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
snow.$(OBJEXT): snow_tile_mod.mod

# Mod files are built and then installed as headers.
MODFILES = snow_mod.mod snow_tile_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
