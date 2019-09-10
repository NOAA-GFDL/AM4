# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/shallow_cu.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libuw_conv.la libconv_plumes_k.la libconv_plumes.la \
libconv_utilities_k.la libdeep_conv.la libconv_closures.la \
libconv_utilities.la

# Each convenience library depends on its source.
libuw_conv_la_SOURCES = uw_conv.F90
libconv_plumes_k_la_SOURCES = conv_plumes_k.F90
libconv_plumes_la_SOURCES = conv_plumes.F90
libconv_utilities_k_la_SOURCES = conv_utilities_k.F90
libdeep_conv_la_SOURCES = deep_conv.F90
libconv_closures_la_SOURCES = conv_closures.F90
libconv_utilities_la_SOURCES = conv_utilities.F90

# Each mod file depends on the .o file.
uw_conv_mod.mod: uw_conv.$(OBJEXT)
conv_plumes_k_mod.mod: conv_plumes_k.$(OBJEXT)
conv_plumes_mod.mod: conv_plumes.$(OBJEXT)
conv_utilities_k_mod.mod: conv_utilities_k.$(OBJEXT)
deep_conv_mod.mod: deep_conv.$(OBJEXT)
conv_closures_mod.mod: conv_closures.$(OBJEXT)
conv_utilities_mod.mod: conv_utilities.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
conv_utilities.$(OBJEXT): conv_utilities_k_mod.mod

# Mod files are built and then installed as headers.
MODFILES = uw_conv_mod.mod conv_plumes_k_mod.mod conv_plumes_mod.mod \
conv_utilities_k_mod.mod deep_conv_mod.mod conv_closures_mod.mod \
conv_utilities_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
