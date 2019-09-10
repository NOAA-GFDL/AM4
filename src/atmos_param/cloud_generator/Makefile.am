# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/cloud_generator.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libcloud_generator.la libbetaDistribution.la

# Each convenience library depends on its source.
libcloud_generator_la_SOURCES = cloud_generator.F90
libbetaDistribution_la_SOURCES = betaDistribution.F90

# Each mod file depends on the .o file.
cloud_generator_mod.mod: cloud_generator.$(OBJEXT)
betaDistribution_mod.mod: betaDistribution.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = cloud_generator_mod.mod betaDistribution_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
