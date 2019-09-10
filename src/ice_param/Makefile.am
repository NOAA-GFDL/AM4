# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/ice_param.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libocean_rough.la libocean_albedo.la libice_albedo.la

# Each convenience library depends on its source.
libocean_rough_la_SOURCES = ocean_rough.F90
libocean_albedo_la_SOURCES = ocean_albedo.F90
libice_albedo_la_SOURCES = ice_albedo.F90

# Each mod file depends on the .o file.
ocean_rough_mod.mod: ocean_rough.$(OBJEXT)
ocean_albedo_mod.mod: ocean_albedo.$(OBJEXT)
ice_albedo_mod.mod: ice_albedo.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = ocean_rough_mod.mod ocean_albedo_mod.mod ice_albedo_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = ocean_albedo.html ice_albedo.html ocean_rough.html

CLEANFILES = *.mod
