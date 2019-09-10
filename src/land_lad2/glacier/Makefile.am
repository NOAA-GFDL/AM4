# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/land_lad2/glacier.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libglac_tile.la libglacier.la

# Each convenience library depends on its source.
libglac_tile_la_SOURCES = glac_tile.F90
libglacier_la_SOURCES = glacier.F90

# Each mod file depends on the .o file.
glac_tile_mod.mod: glac_tile.$(OBJEXT)
glacier_mod.mod: glacier.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
glacier.$(OBJEXT): glac_tile_mod.mod

# Mod files are built and then installed as headers.
MODFILES = glac_tile_mod.mod glacier_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
