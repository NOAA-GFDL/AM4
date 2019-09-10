# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/dry_adj.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libdry_adj.la

# Each convenience library depends on its source.
libdry_adj_la_SOURCES = dry_adj.F90

# Each mod file depends on the .o file.
dry_adj_mod.mod: dry_adj.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = dry_adj_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = dry_adj.html

CLEANFILES = *.mod
