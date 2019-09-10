# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/cloud_obs.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libcloud_obs.la

# Each convenience library depends on its source.
libcloud_obs_la_SOURCES = cloud_obs.F90

# Each mod file depends on the .o file.
cloud_obs_mod.mod: cloud_obs.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = cloud_obs_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = cloud_obs.html

CLEANFILES = *.mod
