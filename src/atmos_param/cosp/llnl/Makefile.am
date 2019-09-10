# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/cosp/llnl.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libcosp_radar.la libllnl_stats.la

# Each convenience library depends on its source.
libcosp_radar_la_SOURCES = cosp_radar.F90
libllnl_stats_la_SOURCES = llnl_stats.F90

# Each mod file depends on the .o file.
cosp_radar_mod.mod: cosp_radar.$(OBJEXT)
llnl_stats_mod.mod: llnl_stats.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = cosp_radar_mod.mod llnl_stats_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
