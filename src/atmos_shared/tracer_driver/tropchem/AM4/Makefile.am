# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_shared/tracer_driver/tropchem/AM4.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libmoz_mods.la libmo_fphoto.la libmoz_subs.la libmoz_mat.la \
libmo_fastjx.la

# Each convenience library depends on its source.
libmoz_mods_la_SOURCES = moz.mods.F90
libmo_fphoto_la_SOURCES = mo_fphoto.F90
libmoz_subs_la_SOURCES = moz.subs.F90
libmoz_mat_la_SOURCES = moz.mat.F90
libmo_fastjx_la_SOURCES = mo_fastjx.F90

# Each mod file depends on the .o file.
moz.mods_mod.mod: moz.mods.$(OBJEXT)
mo_fphoto_mod.mod: mo_fphoto.$(OBJEXT)
moz.subs_mod.mod: moz.subs.$(OBJEXT)
moz.mat_mod.mod: moz.mat.$(OBJEXT)
mo_fastjx_mod.mod: mo_fastjx.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = moz_mods_mod.mod mo_fphoto_mod.mod moz_subs_mod.mod moz_mat_mod.mod \
mo_fastjx_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
