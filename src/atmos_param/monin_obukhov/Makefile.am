# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/monin_obukhov.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libmonin_obukhov.la libmonin_obukhov_kernel.la

# Each convenience library depends on its source.
libmonin_obukhov_la_SOURCES = monin_obukhov.F90
libmonin_obukhov_kernel_la_SOURCES = monin_obukhov_kernel.F90

# Each mod file depends on the .o file.
monin_obukhov_mod.mod: monin_obukhov.$(OBJEXT)
monin_obukhov_kernel_mod.mod: monin_obukhov_kernel.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
monin_obukhov_kernel.$(OBJEXT): monin_obukhov_mod.mod

# Mod files are built and then installed as headers.
MODFILES = monin_obukhov_mod.mod monin_obukhov_kernel_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = monin_obukhov_interfaces.h monin_obukhov.tech.ps monin_obukhov.html

CLEANFILES = *.mod
