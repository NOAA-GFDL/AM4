# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_shared/tracer_driver/stratchem.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libstrat_chem_driver.la libstrat_chem_model.la

# Each convenience library depends on its source.
libstrat_chem_driver_la_SOURCES = strat_chem_driver.F90
libstrat_chem_model_la_SOURCES = strat_chem_model.F90

# Each mod file depends on the .o file.
strat_chem_driver_mod.mod: strat_chem_driver.$(OBJEXT)
strat_chem_model_mod.mod: strat_chem_model.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = strat_chem_driver_mod.mod strat_chem_model_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
