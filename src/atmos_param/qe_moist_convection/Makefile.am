# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/qe_moist_convection.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libqe_moist_convection.la

# Each convenience library depends on its source.
libqe_moist_convection_la_SOURCES = qe_moist_convection.F90

# Each mod file depends on the .o file.
qe_moist_convection_mod.mod: qe_moist_convection.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = qe_moist_convection_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
