# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/edt.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libedt.la

# Each convenience library depends on its source.
libedt_la_SOURCES = edt.F90

# Each mod file depends on the .o file.
edt_mod.mod: edt.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = edt_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
