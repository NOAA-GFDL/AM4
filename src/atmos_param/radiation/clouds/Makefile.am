# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/radiation/clouds.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libisccp_clouds.la libcloudrad_types.la libdonner_deep_clouds_W.la \
libget_random_number_stream.la libstrat_clouds_W.la \
libcloudrad_diagnostics.la libuw_clouds_W.la libcloudrad_package.la \
libcloudrad_driver.la libbulkphys_rad.la libcloud_spec.la \
libmicrophys_rad.la

# Each convenience library depends on its source.
libisccp_clouds_la_SOURCES = isccp_clouds.F90
libcloudrad_types_la_SOURCES = cloudrad_types.F90
libdonner_deep_clouds_W_la_SOURCES = donner_deep_clouds_W.F90
libget_random_number_stream_la_SOURCES = get_random_number_stream.F90
libstrat_clouds_W_la_SOURCES = strat_clouds_W.F90
libcloudrad_diagnostics_la_SOURCES = cloudrad_diagnostics.F90
libuw_clouds_W_la_SOURCES = uw_clouds_W.F90
libcloudrad_package_la_SOURCES = cloudrad_package.F90
libcloudrad_driver_la_SOURCES = cloudrad_driver.F90
libbulkphys_rad_la_SOURCES = bulkphys_rad.F90
libcloud_spec_la_SOURCES = cloud_spec.F90
libmicrophys_rad_la_SOURCES = microphys_rad.F90

# Each mod file depends on the .o file.
isccp_clouds_mod.mod: isccp_clouds.$(OBJEXT)
cloudrad_types_mod.mod: cloudrad_types.$(OBJEXT)
donner_deep_clouds_W_mod.mod: donner_deep_clouds_W.$(OBJEXT)
get_random_number_stream_mod.mod: get_random_number_stream.$(OBJEXT)
strat_clouds_W_mod.mod: strat_clouds_W.$(OBJEXT)
cloudrad_diagnostics_mod.mod: cloudrad_diagnostics.$(OBJEXT)
uw_clouds_W_mod.mod: uw_clouds_W.$(OBJEXT)
cloudrad_package_mod.mod: cloudrad_package.$(OBJEXT)
cloudrad_driver_mod.mod: cloudrad_driver.$(OBJEXT)
bulkphys_rad_mod.mod: bulkphys_rad.$(OBJEXT)
cloud_spec_mod.mod: cloud_spec.$(OBJEXT)
microphys_rad_mod.mod: microphys_rad.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
donner_deep_clouds_W.$(OBJEXT): cloudrad_types_mod.mod
get_random_number_stream.$(OBJEXT): cloudrad_types_mod.mod
strat_clouds_W.$(OBJEXT): cloudrad_types_mod.mod
cloudrad_diagnostics.$(OBJEXT): isccp_clouds_mod.mod cloudrad_types_mod.mod microphys_rad_mod.mod
uw_clouds_W.$(OBJEXT): cloudrad_types_mod.mod
cloudrad_package.$(OBJEXT): cloudrad_types_mod.mod cloudrad_diagnostics_mod.mod bulkphys_rad_mod.mod microphys_rad_mod.mod
cloudrad_driver.$(OBJEXT): cloudrad_types_mod.mod cloudrad_package_mod.mod cloud_spec_mod.mod
bulkphys_rad.$(OBJEXT): cloudrad_types_mod.mod strat_clouds_W_mod.mod
cloud_spec.$(OBJEXT): cloudrad_types_mod.mod donner_deep_clouds_W_mod.mod strat_clouds_W_mod.mod uw_clouds_W_mod.mod
microphys_rad.$(OBJEXT): cloudrad_types_mod.mod

# Mod files are built and then installed as headers.
MODFILES = isccp_clouds_mod.mod cloudrad_types_mod.mod \
donner_deep_clouds_W_mod.mod get_random_number_stream_mod.mod \
strat_clouds_W_mod.mod cloudrad_diagnostics_mod.mod \
uw_clouds_W_mod.mod cloudrad_package_mod.mod cloudrad_driver_mod.mod \
bulkphys_rad_mod.mod cloud_spec_mod.mod microphys_rad_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
