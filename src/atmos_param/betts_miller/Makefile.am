# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/betts_miller.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libbetts_miller.la libbm_omp.la libbm_massflux.la

# Each convenience library depends on its source.
libbetts_miller_la_SOURCES = betts_miller.F90
libbm_omp_la_SOURCES = bm_omp.F90
libbm_massflux_la_SOURCES = bm_massflux.F90

# Each mod file depends on the .o file.
betts_miller_mod.mod: betts_miller.$(OBJEXT)
bm_omp_mod.mod: bm_omp.$(OBJEXT)
bm_massflux_mod.mod: bm_massflux.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = betts_miller_mod.mod bm_omp_mod.mod bm_massflux_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
