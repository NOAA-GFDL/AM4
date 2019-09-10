# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/lscale_cond.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = liblscale_cond.la

# Each convenience library depends on its source.
liblscale_cond_la_SOURCES = lscale_cond.F90

# Each mod file depends on the .o file.
lscale_cond_mod.mod: lscale_cond.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = lscale_cond_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = lscale_cond.html

CLEANFILES = *.mod
