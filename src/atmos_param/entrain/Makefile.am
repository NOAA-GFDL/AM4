# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/entrain.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libentrain.la

# Each convenience library depends on its source.
libentrain_la_SOURCES = entrain.F90

# Each mod file depends on the .o file.
entrain_mod.mod: entrain.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = entrain_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = entrain.html

CLEANFILES = *.mod
