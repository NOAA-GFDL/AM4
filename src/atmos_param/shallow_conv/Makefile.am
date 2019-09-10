# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/shallow_conv.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libshallow_conv.la

# Each convenience library depends on its source.
libshallow_conv_la_SOURCES = shallow_conv.F90

# Each mod file depends on the .o file.
shallow_conv_mod.mod: shallow_conv.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = shallow_conv_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = shallow_conv.html

CLEANFILES = *.mod
