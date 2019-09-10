# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/diag_integral.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libdiag_integral.la

# Each convenience library depends on its source.
libdiag_integral_la_SOURCES = diag_integral.F90

# Each mod file depends on the .o file.
diag_integral_mod.mod: diag_integral.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = diag_integral_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = diag_integral.html

CLEANFILES = *.mod
