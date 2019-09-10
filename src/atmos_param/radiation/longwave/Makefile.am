# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/radiation/longwave.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = liblongwave_fluxes.la liblongwave_clouds.la libsealw99.la \
liblongwave_params.la liblongwave_tables.la libgas_tf.la \
liblongwave_utilities.la liblw_gases_stdtf.la liboptical_path.la \
liblongwave_types.la

# Each convenience library depends on its source.
liblongwave_fluxes_la_SOURCES = longwave_fluxes.F90
liblongwave_clouds_la_SOURCES = longwave_clouds.F90
libsealw99_la_SOURCES = sealw99.F90
liblongwave_params_la_SOURCES = longwave_params.F90
liblongwave_tables_la_SOURCES = longwave_tables.F90
libgas_tf_la_SOURCES = gas_tf.F90
liblongwave_utilities_la_SOURCES = longwave_utilities.F90
liblw_gases_stdtf_la_SOURCES = lw_gases_stdtf.F90
liboptical_path_la_SOURCES = optical_path.F90
liblongwave_types_la_SOURCES = longwave_types.F90

# Each mod file depends on the .o file.
longwave_fluxes_mod.mod: longwave_fluxes.$(OBJEXT)
longwave_clouds_mod.mod: longwave_clouds.$(OBJEXT)
sealw99_mod.mod: sealw99.$(OBJEXT)
longwave_params_mod.mod: longwave_params.$(OBJEXT)
longwave_tables_mod.mod: longwave_tables.$(OBJEXT)
gas_tf_mod.mod: gas_tf.$(OBJEXT)
longwave_utilities_mod.mod: longwave_utilities.$(OBJEXT)
lw_gases_stdtf_mod.mod: lw_gases_stdtf.$(OBJEXT)
optical_path_mod.mod: optical_path.$(OBJEXT)
longwave_types_mod.mod: longwave_types.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
longwave_fluxes.$(OBJEXT): longwave_types_mod.mod
longwave_clouds.$(OBJEXT): longwave_utilities_mod.mod longwave_types_mod.mod
sealw99.$(OBJEXT): longwave_fluxes_mod.mod longwave_clouds_mod.mod longwave_params_mod.mod longwave_tables_mod.mod gas_tf_mod.mod longwave_utilities_mod.mod lw_gases_stdtf_mod.mod optical_path_mod.mod longwave_types_mod.mod
longwave_tables.$(OBJEXT): longwave_params_mod.mod longwave_utilities_mod.mod
gas_tf.$(OBJEXT): longwave_params_mod.mod longwave_utilities_mod.mod
lw_gases_stdtf.$(OBJEXT): gas_tf_mod.mod longwave_utilities_mod.mod
optical_path.$(OBJEXT): longwave_params_mod.mod longwave_utilities_mod.mod lw_gases_stdtf_mod.mod

# Mod files are built and then installed as headers.
MODFILES = longwave_fluxes_mod.mod longwave_clouds_mod.mod sealw99_mod.mod \
longwave_params_mod.mod longwave_tables_mod.mod gas_tf_mod.mod \
longwave_utilities_mod.mod lw_gases_stdtf_mod.mod optical_path_mod.mod \
longwave_types_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
