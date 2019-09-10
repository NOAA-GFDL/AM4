# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/macrophysics.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libtiedtke_macro.la libls_cloud_macrophysics.la

# Each convenience library depends on its source.
libtiedtke_macro_la_SOURCES = tiedtke_macro.F90
libls_cloud_macrophysics_la_SOURCES = ls_cloud_macrophysics.F90

# Each mod file depends on the .o file.
tiedtke_macro_mod.mod: tiedtke_macro.$(OBJEXT)
ls_cloud_macrophysics_mod.mod: ls_cloud_macrophysics.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
ls_cloud_macrophysics.$(OBJEXT): tiedtke_macro_mod.mod

# Mod files are built and then installed as headers.
MODFILES = tiedtke_macro_mod.mod ls_cloud_macrophysics_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
