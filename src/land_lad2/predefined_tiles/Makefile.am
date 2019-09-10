# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/land_lad2/predefined_tiles.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libtiling_input.la libtiling_input_types.la

# Each convenience library depends on its source.
libtiling_input_la_SOURCES = tiling_input.F90
libtiling_input_types_la_SOURCES = tiling_input_types.F90

# Each mod file depends on the .o file.
tiling_input_mod.mod: tiling_input.$(OBJEXT)
tiling_input_types_mod.mod: tiling_input_types.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
tiling_input.$(OBJEXT): tiling_input_types_mod.mod

# Mod files are built and then installed as headers.
MODFILES = tiling_input_mod.mod tiling_input_types_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
