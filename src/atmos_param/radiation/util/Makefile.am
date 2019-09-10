# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/radiation/util.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libradiation_diag.la

# Each convenience library depends on its source.
libradiation_diag_la_SOURCES = radiation_diag.F90

# Each mod file depends on the .o file.
radiation_diag_mod.mod: radiation_diag.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = radiation_diag_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
