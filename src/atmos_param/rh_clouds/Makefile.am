# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/rh_clouds.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = librh_clouds.la

# Each convenience library depends on its source.
librh_clouds_la_SOURCES = rh_clouds.F90

# Each mod file depends on the .o file.
rh_clouds_mod.mod: rh_clouds.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = rh_clouds_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = rh_clouds.html

CLEANFILES = *.mod
