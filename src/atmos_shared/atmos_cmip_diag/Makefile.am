# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_shared/atmos_cmip_diag.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libatmos_cmip_diag.la libatmos_global_diag.la

# Each convenience library depends on its source.
libatmos_cmip_diag_la_SOURCES = atmos_cmip_diag.F90
libatmos_global_diag_la_SOURCES = atmos_global_diag.F90

# Each mod file depends on the .o file.
atmos_cmip_diag_mod.mod: atmos_cmip_diag.$(OBJEXT)
atmos_global_diag_mod.mod: atmos_global_diag.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = atmos_cmip_diag_mod.mod atmos_global_diag_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = atmos_cmip_interp.inc

CLEANFILES = *.mod
