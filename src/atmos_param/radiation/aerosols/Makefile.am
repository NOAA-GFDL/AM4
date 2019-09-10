# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/radiation/aerosols.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libaerosolrad_driver.la libaerosolrad_package.la \
libaerosolrad_types.la

# Each convenience library depends on its source.
libaerosolrad_driver_la_SOURCES = aerosolrad_driver.F90
libaerosolrad_package_la_SOURCES = aerosolrad_package.F90
libaerosolrad_types_la_SOURCES = aerosolrad_types.F90

# Each mod file depends on the .o file.
aerosolrad_driver_mod.mod: aerosolrad_driver.$(OBJEXT)
aerosolrad_package_mod.mod: aerosolrad_package.$(OBJEXT)
aerosolrad_types_mod.mod: aerosolrad_types.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
aerosolrad_driver.$(OBJEXT): aerosolrad_package_mod.mod aerosolrad_types_mod.mod
aerosolrad_package.$(OBJEXT): aerosolrad_types_mod.mod

# Mod files are built and then installed as headers.
MODFILES = aerosolrad_driver_mod.mod aerosolrad_package_mod.mod \
aerosolrad_types_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
