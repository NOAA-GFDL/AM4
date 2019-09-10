# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/cosp/icarus-scops-3.7.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libicarus.la

# Each convenience library depends on its source.
libicarus_la_SOURCES = icarus.F90

# Each mod file depends on the .o file.
icarus_mod.mod: icarus.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = icarus_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = congvec.inc

CLEANFILES = *.mod
