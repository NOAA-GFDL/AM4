# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/land_lad2/transitions.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libtransitions.la

# Each convenience library depends on its source.
libtransitions_la_SOURCES = transitions.F90

# Each mod file depends on the .o file.
transitions_mod.mod: transitions.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = transitions_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
