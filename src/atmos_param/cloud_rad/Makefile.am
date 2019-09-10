# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/cloud_rad.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libcloud_rad.la

# Each convenience library depends on its source.
libcloud_rad_la_SOURCES = cloud_rad.F90

# Each mod file depends on the .o file.
cloud_rad_mod.mod: cloud_rad.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = cloud_rad_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = cloud_rad.html

CLEANFILES = *.mod
