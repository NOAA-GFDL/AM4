# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/lscloud_driver.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libpolysvp.la liblscloud_debug.la liblscloud_driver.la \
liblscloud_types.la liblscloud_constants.la libcheck_nan.la \
liblscloud_netcdf.la

# Each convenience library depends on its source.
libpolysvp_la_SOURCES = polysvp.F90
liblscloud_debug_la_SOURCES = lscloud_debug.F90
liblscloud_driver_la_SOURCES = lscloud_driver.F90
liblscloud_types_la_SOURCES = lscloud_types.F90
liblscloud_constants_la_SOURCES = lscloud_constants.F90
libcheck_nan_la_SOURCES = check_nan.F90
liblscloud_netcdf_la_SOURCES = lscloud_netcdf.F90

# Each mod file depends on the .o file.
polysvp_mod.mod: polysvp.$(OBJEXT)
lscloud_debug_mod.mod: lscloud_debug.$(OBJEXT)
lscloud_driver_mod.mod: lscloud_driver.$(OBJEXT)
lscloud_types_mod.mod: lscloud_types.$(OBJEXT)
lscloud_constants_mod.mod: lscloud_constants.$(OBJEXT)
check_nan_mod.mod: check_nan.$(OBJEXT)
lscloud_netcdf_mod.mod: lscloud_netcdf.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
polysvp.$(OBJEXT): lscloud_types_mod.mod lscloud_constants_mod.mod
lscloud_debug.$(OBJEXT): lscloud_types_mod.mod check_nan_mod.mod
lscloud_driver.$(OBJEXT): polysvp_mod.mod lscloud_debug_mod.mod lscloud_types_mod.mod lscloud_netcdf_mod.mod
lscloud_netcdf.$(OBJEXT): lscloud_types_mod.mod

# Mod files are built and then installed as headers.
MODFILES = polysvp_mod.mod lscloud_debug_mod.mod lscloud_driver_mod.mod \
lscloud_types_mod.mod lscloud_constants_mod.mod check_nan_mod.mod \
lscloud_netcdf_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
