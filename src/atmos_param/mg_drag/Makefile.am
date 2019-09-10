# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/mg_drag.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libmg_drag.la

# Each convenience library depends on its source.
libmg_drag_la_SOURCES = mg_drag.F90

# Each mod file depends on the .o file.
mg_drag_mod.mod: mg_drag.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = mg_drag_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = mg_drag.html

CLEANFILES = *.mod
