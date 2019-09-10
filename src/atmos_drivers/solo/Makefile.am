# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_drivers/solo.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libatmos_model.la

# Each convenience library depends on its source.
libatmos_model_la_SOURCES = atmos_model.F90

# Each mod file depends on the .o file.
atmos_model_mod.mod: atmos_model.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = atmos_model_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = atmos_model.html

CLEANFILES = *.mod
