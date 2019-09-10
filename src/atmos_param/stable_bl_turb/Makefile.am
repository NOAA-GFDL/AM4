# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/stable_bl_turb.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libstable_bl_turb.la

# Each convenience library depends on its source.
libstable_bl_turb_la_SOURCES = stable_bl_turb.F90

# Each mod file depends on the .o file.
stable_bl_turb_mod.mod: stable_bl_turb.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = stable_bl_turb_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = stable_bl_turb.html

CLEANFILES = *.mod
