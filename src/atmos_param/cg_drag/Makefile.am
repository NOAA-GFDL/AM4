# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/cg_drag.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libcg_drag.la

# Each convenience library depends on its source.
libcg_drag_la_SOURCES = cg_drag.F90

# Each mod file depends on the .o file.
cg_drag_mod.mod: cg_drag.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = cg_drag_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
