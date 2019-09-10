# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/physics_driver.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libphysics_types.la libphysics_driver.la

# Each convenience library depends on its source.
libphysics_types_la_SOURCES = physics_types.F90
libphysics_driver_la_SOURCES = physics_driver.F90

# Each mod file depends on the .o file.
physics_types_mod.mod: physics_types.$(OBJEXT)
physics_driver_mod.mod: physics_driver.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
physics_driver.$(OBJEXT): physics_types_mod.mod

# Mod files are built and then installed as headers.
MODFILES = physics_types_mod.mod physics_driver_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = physics_driver.html Radiation_overview.tech.ps

CLEANFILES = *.mod
