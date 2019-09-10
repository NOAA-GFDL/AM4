# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/vert_turb_driver.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libvert_turb_driver.la

# Each convenience library depends on its source.
libvert_turb_driver_la_SOURCES = vert_turb_driver.F90

# Each mod file depends on the .o file.
vert_turb_driver_mod.mod: vert_turb_driver.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = vert_turb_driver_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = vert_turb_driver.html

CLEANFILES = *.mod
