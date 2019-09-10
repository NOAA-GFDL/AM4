# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/grey_radiation.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libgrey_radiation.la

# Each convenience library depends on its source.
libgrey_radiation_la_SOURCES = grey_radiation.F90

# Each mod file depends on the .o file.
grey_radiation_mod.mod: grey_radiation.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = grey_radiation_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
