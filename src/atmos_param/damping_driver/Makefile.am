# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/damping_driver.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libdamping_driver.la

# Each convenience library depends on its source.
libdamping_driver_la_SOURCES = damping_driver.F90

# Each mod file depends on the .o file.
damping_driver_mod.mod: damping_driver.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = damping_driver_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = damping_driver.html

CLEANFILES = *.mod
