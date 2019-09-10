# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/ras.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libras.la

# Each convenience library depends on its source.
libras_la_SOURCES = ras.F90

# Each mod file depends on the .o file.
ras_mod.mod: ras.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = ras_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = ras.html

CLEANFILES = *.mod
