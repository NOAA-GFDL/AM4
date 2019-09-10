# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/tke_turb.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libtke_turb.la

# Each convenience library depends on its source.
libtke_turb_la_SOURCES = tke_turb.F90

# Each mod file depends on the .o file.
tke_turb_mod.mod: tke_turb.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = tke_turb_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
