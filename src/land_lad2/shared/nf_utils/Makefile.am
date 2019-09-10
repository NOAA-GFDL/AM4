# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/land_lad2/shared/nf_utils.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libnf_utils.la libnfc.la libnfu.la

# Each convenience library depends on its source.
libnf_utils_la_SOURCES = nf_utils.F90
libnfc_la_SOURCES = nfc.F90
libnfu_la_SOURCES = nfu.F90

# Each mod file depends on the .o file.
nf_utils_mod.mod: nf_utils.$(OBJEXT)
nfc_mod.mod: nfc.$(OBJEXT)
nfu_mod.mod: nfu.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
nf_utils.$(OBJEXT): nfc_mod.mod nfu_mod.mod
nfc.$(OBJEXT): nfu_mod.mod

# Mod files are built and then installed as headers.
MODFILES = nf_utils_mod.mod nfc_mod.mod nfu_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = getput.inc getput_compressed.inc

CLEANFILES = *.mod
