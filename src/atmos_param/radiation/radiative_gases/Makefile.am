# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/radiation/radiative_gases.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libradiative_gases_types.la libozone.la libradiative_gases.la

# Each convenience library depends on its source.
libradiative_gases_types_la_SOURCES = radiative_gases_types.F90
libozone_la_SOURCES = ozone.F90
libradiative_gases_la_SOURCES = radiative_gases.F90

# Each mod file depends on the .o file.
radiative_gases_types_mod.mod: radiative_gases_types.$(OBJEXT)
ozone_mod.mod: ozone.$(OBJEXT)
radiative_gases_mod.mod: radiative_gases.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
ozone.$(OBJEXT): radiative_gases_types_mod.mod radiative_gases_mod.mod
radiative_gases.$(OBJEXT): radiative_gases_types_mod.mod ozone_mod.mod

# Mod files are built and then installed as headers.
MODFILES = radiative_gases_types_mod.mod ozone_mod.mod radiative_gases_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
