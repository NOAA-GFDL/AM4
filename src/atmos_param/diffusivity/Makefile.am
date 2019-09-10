# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/diffusivity.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libdiffusivity.la

# Each convenience library depends on its source.
libdiffusivity_la_SOURCES = diffusivity.F90

# Each mod file depends on the .o file.
diffusivity_mod.mod: diffusivity.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = diffusivity_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = diffusivity.html diffusivity.tech.ps

CLEANFILES = *.mod
