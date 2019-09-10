# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/cosp.
#
# Ed Hartnett 2019-07-30

# Run make targets in each of these subdirs
SUBDIRS = icarus-scops-3.7 MODIS_simulator llnl actsim MISR_simulator quickbeam

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libcosp_simulator.la libcosp.la libcosp_diagnostics.la \
libcosp_modis_simulator.la libcosp_stats.la libcosp_constants.la \
libcosp_misr_simulator.la libcosp_utils.la libcosp_io.la \
libcosp_types.la libcosp_rttov_simulator.la libcosp_driver.la \
libcosp_lidar.la libcosp_isccp_simulator.la

# Each convenience library depends on its source.
libcosp_simulator_la_SOURCES = cosp_simulator.F90
libcosp_la_SOURCES = cosp.F90
libcosp_diagnostics_la_SOURCES = cosp_diagnostics.F90
libcosp_modis_simulator_la_SOURCES = cosp_modis_simulator.F90
libcosp_stats_la_SOURCES = cosp_stats.F90
libcosp_constants_la_SOURCES = cosp_constants.F90
libcosp_misr_simulator_la_SOURCES = cosp_misr_simulator.F90
libcosp_utils_la_SOURCES = cosp_utils.F90
libcosp_io_la_SOURCES = cosp_io.F90
libcosp_types_la_SOURCES = cosp_types.F90
libcosp_rttov_simulator_la_SOURCES = cosp_rttov_simulator.F90
libcosp_driver_la_SOURCES = cosp_driver.F90
libcosp_lidar_la_SOURCES = cosp_lidar.F90
libcosp_isccp_simulator_la_SOURCES = cosp_isccp_simulator.F90

# Each mod file depends on the .o file.
cosp_simulator_mod.mod: cosp_simulator.$(OBJEXT)
cosp_mod.mod: cosp.$(OBJEXT)
cosp_diagnostics_mod.mod: cosp_diagnostics.$(OBJEXT)
cosp_modis_simulator_mod.mod: cosp_modis_simulator.$(OBJEXT)
cosp_stats_mod.mod: cosp_stats.$(OBJEXT)
cosp_constants_mod.mod: cosp_constants.$(OBJEXT)
cosp_misr_simulator_mod.mod: cosp_misr_simulator.$(OBJEXT)
cosp_utils_mod.mod: cosp_utils.$(OBJEXT)
cosp_io_mod.mod: cosp_io.$(OBJEXT)
cosp_types_mod.mod: cosp_types.$(OBJEXT)
cosp_rttov_simulator_mod.mod: cosp_rttov_simulator.$(OBJEXT)
cosp_driver_mod.mod: cosp_driver.$(OBJEXT)
cosp_lidar_mod.mod: cosp_lidar.$(OBJEXT)
cosp_isccp_simulator_mod.mod: cosp_isccp_simulator.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
cosp_driver.$(OBJEXT): cosp_mod.mod cosp_diagnostics_mod.mod

# Mod files are built and then installed as headers.
MODFILES = cosp_simulator_mod.mod cosp_mod.mod cosp_diagnostics_mod.mod \
cosp_modis_simulator_mod.mod cosp_stats_mod.mod cosp_constants_mod.mod \
cosp_misr_simulator_mod.mod cosp_utils_mod.mod cosp_io_mod.mod \
cosp_types_mod.mod cosp_rttov_simulator_mod.mod cosp_driver_mod.mod \
cosp_lidar_mod.mod cosp_isccp_simulator_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
