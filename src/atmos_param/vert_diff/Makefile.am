# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/vert_diff.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libvert_diff.la

# Each convenience library depends on its source.
libvert_diff_la_SOURCES = vert_diff.F90

# Each mod file depends on the .o file.
vert_diff_mod.mod: vert_diff.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = vert_diff_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = vert_diff.html vert_diff.tech.ps

CLEANFILES = *.mod
