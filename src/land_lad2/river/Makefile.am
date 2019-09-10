# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/land_lad2/river.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libriver.la libriver_type.la libriver_physics.la

# Each convenience library depends on its source.
libriver_la_SOURCES = river.F90
libriver_type_la_SOURCES = river_type.F90
libriver_physics_la_SOURCES = river_physics.F90

# Each mod file depends on the .o file.
river_mod.mod: river.$(OBJEXT)
river_type_mod.mod: river_type.$(OBJEXT)
river_physics_mod.mod: river_physics.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
river.$(OBJEXT): river_type_mod.mod river_physics_mod.mod
river_physics.$(OBJEXT): river_mod.mod river_type_mod.mod

# Mod files are built and then installed as headers.
MODFILES = river_mod.mod river_type_mod.mod river_physics_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
