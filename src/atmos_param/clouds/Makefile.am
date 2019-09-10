# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/clouds.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libclouds.la

# Each convenience library depends on its source.
libclouds_la_SOURCES = clouds.F90

# Each mod file depends on the .o file.
clouds_mod.mod: clouds.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = clouds_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = clouds.html

CLEANFILES = *.mod
