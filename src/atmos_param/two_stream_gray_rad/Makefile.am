# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/two_stream_gray_rad.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libtwo_stream_gray_rad.la

# Each convenience library depends on its source.
libtwo_stream_gray_rad_la_SOURCES = two_stream_gray_rad.F90

# Each mod file depends on the .o file.
two_stream_gray_rad_mod.mod: two_stream_gray_rad.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = two_stream_gray_rad_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
