# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/cosp/actsim.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = liblidar_simulator.la liblmd_ipsl_stats.la

# Each convenience library depends on its source.
liblidar_simulator_la_SOURCES = lidar_simulator.F90
liblmd_ipsl_stats_la_SOURCES = lmd_ipsl_stats.F90

# Each mod file depends on the .o file.
lidar_simulator_mod.mod: lidar_simulator.$(OBJEXT)
lmd_ipsl_stats_mod.mod: lmd_ipsl_stats.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = lidar_simulator_mod.mod lmd_ipsl_stats_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
