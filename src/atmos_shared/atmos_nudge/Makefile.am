# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_shared/atmos_nudge.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libatmos_nudge.la

# Each convenience library depends on its source.
libatmos_nudge_la_SOURCES = atmos_nudge.F90

# Each mod file depends on the .o file.
atmos_nudge_mod.mod: atmos_nudge.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = atmos_nudge_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
