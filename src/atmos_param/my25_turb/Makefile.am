# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/my25_turb.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libmy25_turb.la

# Each convenience library depends on its source.
libmy25_turb_la_SOURCES = my25_turb.F90

# Each mod file depends on the .o file.
my25_turb_mod.mod: my25_turb.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = my25_turb_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = my25_turb.html

CLEANFILES = *.mod
