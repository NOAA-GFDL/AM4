# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/land_lad2/land_tracers.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libland_tracer_driver.la libland_tracers.la libland_dust.la

# Each convenience library depends on its source.
libland_tracer_driver_la_SOURCES = land_tracer_driver.F90
libland_tracers_la_SOURCES = land_tracers.F90
libland_dust_la_SOURCES = land_dust.F90

# Each mod file depends on the .o file.
land_tracer_driver_mod.mod: land_tracer_driver.$(OBJEXT)
land_tracers_mod.mod: land_tracers.$(OBJEXT)
land_dust_mod.mod: land_dust.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
land_tracer_driver.$(OBJEXT): land_tracers_mod.mod land_dust_mod.mod
land_dust.$(OBJEXT): land_tracers_mod.mod

# Mod files are built and then installed as headers.
MODFILES = land_tracer_driver_mod.mod land_tracers_mod.mod land_dust_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
