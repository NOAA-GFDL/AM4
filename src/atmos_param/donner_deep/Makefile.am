# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/donner_deep.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libdonner_deep_k.la libdonner_meso_k.la libcumulus_closure_k.la \
libdonner_deep.la libdonner_cape_k.la libdonner_rad_k.la \
libdonner_deep_miz.la libfms_donner.la libwet_deposition_0D.la \
libdonner_cloud_model_k.la libdonner_lite_k.la libdonner_types.la \
libnonfms_donner.la libdonner_lscloud_k.la libdonner_utilities_k.la

# Each convenience library depends on its source.
libdonner_deep_k_la_SOURCES = donner_deep_k.F90
libdonner_meso_k_la_SOURCES = donner_meso_k.F90
libcumulus_closure_k_la_SOURCES = cumulus_closure_k.F90
libdonner_deep_la_SOURCES = donner_deep.F90
libdonner_cape_k_la_SOURCES = donner_cape_k.F90
libdonner_rad_k_la_SOURCES = donner_rad_k.F90
libdonner_deep_miz_la_SOURCES = donner_deep_miz.F90
libfms_donner_la_SOURCES = fms_donner.F90
libwet_deposition_0D_la_SOURCES = wet_deposition_0D.F90
libdonner_cloud_model_k_la_SOURCES = donner_cloud_model_k.F90
libdonner_lite_k_la_SOURCES = donner_lite_k.F90
libdonner_types_la_SOURCES = donner_types.F90
libnonfms_donner_la_SOURCES = nonfms_donner.F90
libdonner_lscloud_k_la_SOURCES = donner_lscloud_k.F90
libdonner_utilities_k_la_SOURCES = donner_utilities_k.F90

# Each mod file depends on the .o file.
donner_deep_k_mod.mod: donner_deep_k.$(OBJEXT)
donner_meso_k_mod.mod: donner_meso_k.$(OBJEXT)
cumulus_closure_k_mod.mod: cumulus_closure_k.$(OBJEXT)
donner_deep_mod.mod: donner_deep.$(OBJEXT)
donner_cape_k_mod.mod: donner_cape_k.$(OBJEXT)
donner_rad_k_mod.mod: donner_rad_k.$(OBJEXT)
donner_deep_miz_mod.mod: donner_deep_miz.$(OBJEXT)
fms_donner_mod.mod: fms_donner.$(OBJEXT)
wet_deposition_0D_mod.mod: wet_deposition_0D.$(OBJEXT)
donner_cloud_model_k_mod.mod: donner_cloud_model_k.$(OBJEXT)
donner_lite_k_mod.mod: donner_lite_k.$(OBJEXT)
donner_types_mod.mod: donner_types.$(OBJEXT)
nonfms_donner_mod.mod: nonfms_donner.$(OBJEXT)
donner_lscloud_k_mod.mod: donner_lscloud_k.$(OBJEXT)
donner_utilities_k_mod.mod: donner_utilities_k.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
donner_deep_k.$(OBJEXT): donner_types_mod.mod
donner_meso_k.$(OBJEXT): donner_types_mod.mod
cumulus_closure_k.$(OBJEXT): donner_types_mod.mod
donner_deep.$(OBJEXT): fms_donner_mod.mod donner_types_mod.mod nonfms_donner_mod.mod
donner_cape_k.$(OBJEXT): donner_types_mod.mod
donner_rad_k.$(OBJEXT): donner_types_mod.mod
fms_donner.$(OBJEXT): donner_types_mod.mod
donner_cloud_model_k.$(OBJEXT): donner_types_mod.mod
donner_lite_k.$(OBJEXT): donner_types_mod.mod
nonfms_donner.$(OBJEXT): donner_types_mod.mod
donner_lscloud_k.$(OBJEXT): donner_types_mod.mod
donner_utilities_k.$(OBJEXT): donner_types_mod.mod

# Mod files are built and then installed as headers.
MODFILES = donner_deep_k_mod.mod donner_meso_k_mod.mod cumulus_closure_k_mod.mod \
donner_deep_mod.mod donner_cape_k_mod.mod donner_rad_k_mod.mod \
donner_deep_miz_mod.mod fms_donner_mod.mod wet_deposition_0D_mod.mod \
donner_cloud_model_k_mod.mod donner_lite_k_mod.mod \
donner_types_mod.mod nonfms_donner_mod.mod donner_lscloud_k_mod.mod \
donner_utilities_k_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = Quickstart.ps donner_types.h donner_nml.h

CLEANFILES = *.mod
