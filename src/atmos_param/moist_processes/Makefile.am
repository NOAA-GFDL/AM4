# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/moist_processes.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libmoist_processes.la libmoist_processes_utils.la

# Each convenience library depends on its source.
libmoist_processes_la_SOURCES = moist_processes.F90
libmoist_processes_utils_la_SOURCES = moist_processes_utils.F90

# Each mod file depends on the .o file.
moist_processes_mod.mod: moist_processes.$(OBJEXT)
moist_processes_utils_mod.mod: moist_processes_utils.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = moist_processes_mod.mod moist_processes_utils_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = moist_processes.html

CLEANFILES = *.mod
