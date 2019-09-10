# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/land_lad2/topo_rough.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libtopo_rough.la

# Each convenience library depends on its source.
libtopo_rough_la_SOURCES = topo_rough.F90

# Each mod file depends on the .o file.
topo_rough_mod.mod: topo_rough.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = topo_rough_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
