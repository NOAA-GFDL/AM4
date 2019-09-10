# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/radiation/shortwave.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libshortwave_types.la libesfsw_bands.la libesfsw_parameters.la \
libesfsw_utilities.la libesfsw_driver.la

# Each convenience library depends on its source.
libshortwave_types_la_SOURCES = shortwave_types.F90
libesfsw_bands_la_SOURCES = esfsw_bands.F90
libesfsw_parameters_la_SOURCES = esfsw_parameters.F90
libesfsw_utilities_la_SOURCES = esfsw_utilities.F90
libesfsw_driver_la_SOURCES = esfsw_driver.F90

# Each mod file depends on the .o file.
shortwave_types_mod.mod: shortwave_types.$(OBJEXT)
esfsw_bands_mod.mod: esfsw_bands.$(OBJEXT)
esfsw_parameters_mod.mod: esfsw_parameters.$(OBJEXT)
esfsw_utilities_mod.mod: esfsw_utilities.$(OBJEXT)
esfsw_driver_mod.mod: esfsw_driver.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
esfsw_bands.$(OBJEXT): esfsw_utilities_mod.mod
esfsw_driver.$(OBJEXT): shortwave_types_mod.mod esfsw_bands_mod.mod esfsw_parameters_mod.mod esfsw_utilities_mod.mod

# Mod files are built and then installed as headers.
MODFILES = shortwave_types_mod.mod esfsw_bands_mod.mod esfsw_parameters_mod.mod \
esfsw_utilities_mod.mod esfsw_driver_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
