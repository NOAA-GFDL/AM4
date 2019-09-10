# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_shared/tracer_driver/tropchem/AM3.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libAM3_mods.la libAM3_fastjx.la libAM3_fphoto.la libAM3_subs.la \
libAM3_mat.la

# Each convenience library depends on its source.
libAM3_mods_la_SOURCES = AM3.mods.F90
libAM3_fastjx_la_SOURCES = AM3_fastjx.F90
libAM3_fphoto_la_SOURCES = AM3_fphoto.F90
libAM3_subs_la_SOURCES = AM3.subs.F90
libAM3_mat_la_SOURCES = AM3.mat.F90

# Each mod file depends on the .o file.
AM3.mods_mod.mod: AM3.mods.$(OBJEXT)
AM3_fastjx_mod.mod: AM3_fastjx.$(OBJEXT)
AM3_fphoto_mod.mod: AM3_fphoto.$(OBJEXT)
AM3.subs_mod.mod: AM3.subs.$(OBJEXT)
AM3.mat_mod.mod: AM3.mat.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = AM3_mods_mod.mod AM3_fastjx_mod.mod AM3_fphoto_mod.mod \
AM3_subs_mod.mod AM3_mat_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
