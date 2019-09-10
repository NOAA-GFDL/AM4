# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/ice_sis.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libice_model.la libice_type.la libice_bergs.la libice_grid.la \
librot.la libice_spec.la libice_dyn.la libmask.la libice_thm.la

# Each convenience library depends on its source.
libice_model_la_SOURCES = ice_model.F90
libice_type_la_SOURCES = ice_type.F90
libice_bergs_la_SOURCES = ice_bergs.F90
libice_grid_la_SOURCES = ice_grid.F90
librot_la_SOURCES = rot.F90
libice_spec_la_SOURCES = ice_spec.F90
libice_dyn_la_SOURCES = ice_dyn.F90
libmask_la_SOURCES = mask.F90
libice_thm_la_SOURCES = ice_thm.F90

# Each mod file depends on the .o file.
ice_model_mod.mod: ice_model.$(OBJEXT)
ice_type_mod.mod: ice_type.$(OBJEXT)
ice_bergs_mod.mod: ice_bergs.$(OBJEXT)
ice_grid_mod.mod: ice_grid.$(OBJEXT)
rot_mod.mod: rot.$(OBJEXT)
ice_spec_mod.mod: ice_spec.$(OBJEXT)
ice_dyn_mod.mod: ice_dyn.$(OBJEXT)
mask_mod.mod: mask.$(OBJEXT)
ice_thm_mod.mod: ice_thm.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
ice_model.$(OBJEXT): ice_type_mod.mod ice_bergs_mod.mod ice_grid_mod.mod ice_spec_mod.mod ice_dyn_mod.mod ice_thm_mod.mod
ice_type.$(OBJEXT): ice_bergs_mod.mod ice_grid_mod.mod ice_dyn_mod.mod ice_thm_mod.mod
ice_dyn.$(OBJEXT): ice_grid_mod.mod ice_thm_mod.mod
mask.$(OBJEXT): rot_mod.mod

# Mod files are built and then installed as headers.
MODFILES = ice_model_mod.mod ice_type_mod.mod ice_bergs_mod.mod ice_grid_mod.mod \
rot_mod.mod ice_spec_mod.mod ice_dyn_mod.mod mask_mod.mod \
ice_thm_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = sis.ps sis.html

CLEANFILES = *.mod
