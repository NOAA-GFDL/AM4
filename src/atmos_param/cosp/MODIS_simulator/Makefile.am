# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/cosp/MODIS_simulator.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libmodis_simulator.la

# Each convenience library depends on its source.
libmodis_simulator_la_SOURCES = modis_simulator.F90

# Each mod file depends on the .o file.
modis_simulator_mod.mod: modis_simulator.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = modis_simulator_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
