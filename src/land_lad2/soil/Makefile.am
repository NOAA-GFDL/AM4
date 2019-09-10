# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/land_lad2/soil.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libsoil_tile.la libsoil.la libhillslope.la libhillslope_tile.la \
libhillslope_hydrology.la libuptake.la libsoil_carbon.la

# Each convenience library depends on its source.
libsoil_tile_la_SOURCES = soil_tile.F90
libsoil_la_SOURCES = soil.F90
libhillslope_la_SOURCES = hillslope.F90
libhillslope_tile_la_SOURCES = hillslope_tile.F90
libhillslope_hydrology_la_SOURCES = hillslope_hydrology.F90
libuptake_la_SOURCES = uptake.F90
libsoil_carbon_la_SOURCES = soil_carbon.F90

# Each mod file depends on the .o file.
soil_tile_mod.mod: soil_tile.$(OBJEXT)
soil_mod.mod: soil.$(OBJEXT)
hillslope_mod.mod: hillslope.$(OBJEXT)
hillslope_tile_mod.mod: hillslope_tile.$(OBJEXT)
hillslope_hydrology_mod.mod: hillslope_hydrology.$(OBJEXT)
uptake_mod.mod: uptake.$(OBJEXT)
soil_carbon_mod.mod: soil_carbon.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
soil_tile.$(OBJEXT): soil_mod.mod soil_carbon_mod.mod
soil.$(OBJEXT): soil_tile_mod.mod hillslope_mod.mod hillslope_hydrology_mod.mod uptake_mod.mod soil_carbon_mod.mod
hillslope.$(OBJEXT): soil_tile_mod.mod soil_mod.mod hillslope_tile_mod.mod
hillslope_tile.$(OBJEXT): soil_tile_mod.mod soil_mod.mod
hillslope_hydrology.$(OBJEXT): soil_tile_mod.mod soil_mod.mod hillslope_mod.mod soil_carbon_mod.mod
uptake.$(OBJEXT): soil_tile_mod.mod soil_mod.mod

# Mod files are built and then installed as headers.
MODFILES = soil_tile_mod.mod soil_mod.mod hillslope_mod.mod \
hillslope_tile_mod.mod hillslope_hydrology_mod.mod uptake_mod.mod \
soil_carbon_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
