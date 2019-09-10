# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/lin_cloud_microphys.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = liblin_cloud_microphys.la

# Each convenience library depends on its source.
liblin_cloud_microphys_la_SOURCES = lin_cloud_microphys.F90

# Each mod file depends on the .o file.
lin_cloud_microphys_mod.mod: lin_cloud_microphys.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = lin_cloud_microphys_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
