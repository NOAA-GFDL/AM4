# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/strat_cloud.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libstrat_cloud.la

# Each convenience library depends on its source.
libstrat_cloud_la_SOURCES = strat_cloud.F90

# Each mod file depends on the .o file.
strat_cloud_mod.mod: strat_cloud.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = strat_cloud_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = strat_cloud.tech.ps strat_cloud.html

CLEANFILES = *.mod
