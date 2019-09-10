# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/shallow_physics.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libshallow_physics.la

# Each convenience library depends on its source.
libshallow_physics_la_SOURCES = shallow_physics.F90

# Each mod file depends on the .o file.
shallow_physics_mod.mod: shallow_physics.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = shallow_physics_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = shallow_physics.html

CLEANFILES = *.mod
