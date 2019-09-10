# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_shared/aerosol.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libaerosol.la libaerosol_types.la

# Each convenience library depends on its source.
libaerosol_la_SOURCES = aerosol.F90
libaerosol_types_la_SOURCES = aerosol_types.F90

# Each mod file depends on the .o file.
aerosol_mod.mod: aerosol.$(OBJEXT)
aerosol_types_mod.mod: aerosol_types.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
aerosol.$(OBJEXT): aerosol_types_mod.mod

# Mod files are built and then installed as headers.
MODFILES = aerosol_mod.mod aerosol_types_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
