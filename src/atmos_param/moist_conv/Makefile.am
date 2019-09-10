# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/moist_conv.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libmoist_conv.la

# Each convenience library depends on its source.
libmoist_conv_la_SOURCES = moist_conv.F90

# Each mod file depends on the .o file.
moist_conv_mod.mod: moist_conv.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = moist_conv_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = moist_conv.html

CLEANFILES = *.mod
