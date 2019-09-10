# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_shared/tracer_driver/aer_ccn_act.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libaer_ccn_act_k.la libaer_in_act.la libaerosol_params.la \
libice_nucl.la libaer_ccn_act.la

# Each convenience library depends on its source.
libaer_ccn_act_k_la_SOURCES = aer_ccn_act_k.F90
libaer_in_act_la_SOURCES = aer_in_act.F90
libaerosol_params_la_SOURCES = aerosol_params.F90
libice_nucl_la_SOURCES = ice_nucl.F90
libaer_ccn_act_la_SOURCES = aer_ccn_act.F90

# Each mod file depends on the .o file.
aer_ccn_act_k_mod.mod: aer_ccn_act_k.$(OBJEXT)
aer_in_act_mod.mod: aer_in_act.$(OBJEXT)
aerosol_params_mod.mod: aerosol_params.$(OBJEXT)
ice_nucl_mod.mod: ice_nucl.$(OBJEXT)
aer_ccn_act_mod.mod: aer_ccn_act.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
ice_nucl.$(OBJEXT): aer_ccn_act_k_mod.mod aerosol_params_mod.mod aer_ccn_act_mod.mod
aer_ccn_act.$(OBJEXT): aer_ccn_act_k_mod.mod

# Mod files are built and then installed as headers.
MODFILES = aer_ccn_act_k_mod.mod aer_in_act_mod.mod aerosol_params_mod.mod \
ice_nucl_mod.mod aer_ccn_act_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
