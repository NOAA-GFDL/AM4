# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/topo_drag.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libtopo_drag.la

# Each convenience library depends on its source.
libtopo_drag_la_SOURCES = topo_drag.F90

# Each mod file depends on the .o file.
topo_drag_mod.mod: topo_drag.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = topo_drag_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
