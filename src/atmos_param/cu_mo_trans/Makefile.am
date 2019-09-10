# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/cu_mo_trans.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libcu_mo_trans.la

# Each convenience library depends on its source.
libcu_mo_trans_la_SOURCES = cu_mo_trans.F90

# Each mod file depends on the .o file.
cu_mo_trans_mod.mod: cu_mo_trans.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = cu_mo_trans_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = cu_mo_trans.html

CLEANFILES = *.mod
