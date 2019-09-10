# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_shared/vert_advection.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libvert_advection.la

# Each convenience library depends on its source.
libvert_advection_la_SOURCES = vert_advection.F90

# Each mod file depends on the .o file.
vert_advection_mod.mod: vert_advection.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = vert_advection_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = vert_advection.html

CLEANFILES = *.mod
