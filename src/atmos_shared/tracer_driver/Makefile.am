# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_shared/tracer_driver.
#
# Ed Hartnett 2019-07-30

# Run make targets in each of these subdirs
SUBDIRS = stratchem aer_ccn_act tropchem

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libatmos_tracer_utilities.la libatmos_convection_tracer.la \
libatmos_co2.la libcloud_chem.la libatmos_ch3i.la \
libatmos_tracer_driver.la libatmos_sea_salt.la \
libatmos_carbon_aerosol.la libatmos_dust.la libatmos_radon.la \
libatmos_regional_tracer_driver.la libatmos_soa.la \
libatmos_age_tracer.la libatmos_sulfur_hex.la libatmos_sulfate.la

# Each convenience library depends on its source.
libatmos_tracer_utilities_la_SOURCES = atmos_tracer_utilities.F90
libatmos_convection_tracer_la_SOURCES = atmos_convection_tracer.F90
libatmos_co2_la_SOURCES = atmos_co2.F90
libcloud_chem_la_SOURCES = cloud_chem.F90
libatmos_ch3i_la_SOURCES = atmos_ch3i.F90
libatmos_tracer_driver_la_SOURCES = atmos_tracer_driver.F90
libatmos_sea_salt_la_SOURCES = atmos_sea_salt.F90
libatmos_carbon_aerosol_la_SOURCES = atmos_carbon_aerosol.F90
libatmos_dust_la_SOURCES = atmos_dust.F90
libatmos_radon_la_SOURCES = atmos_radon.F90
libatmos_regional_tracer_driver_la_SOURCES = atmos_regional_tracer_driver.F90
libatmos_soa_la_SOURCES = atmos_soa.F90
libatmos_age_tracer_la_SOURCES = atmos_age_tracer.F90
libatmos_sulfur_hex_la_SOURCES = atmos_sulfur_hex.F90
libatmos_sulfate_la_SOURCES = atmos_sulfate.F90

# Each mod file depends on the .o file.
atmos_tracer_utilities_mod.mod: atmos_tracer_utilities.$(OBJEXT)
atmos_convection_tracer_mod.mod: atmos_convection_tracer.$(OBJEXT)
atmos_co2_mod.mod: atmos_co2.$(OBJEXT)
cloud_chem_mod.mod: cloud_chem.$(OBJEXT)
atmos_ch3i_mod.mod: atmos_ch3i.$(OBJEXT)
atmos_tracer_driver_mod.mod: atmos_tracer_driver.$(OBJEXT)
atmos_sea_salt_mod.mod: atmos_sea_salt.$(OBJEXT)
atmos_carbon_aerosol_mod.mod: atmos_carbon_aerosol.$(OBJEXT)
atmos_dust_mod.mod: atmos_dust.$(OBJEXT)
atmos_radon_mod.mod: atmos_radon.$(OBJEXT)
atmos_regional_tracer_driver_mod.mod: atmos_regional_tracer_driver.$(OBJEXT)
atmos_soa_mod.mod: atmos_soa.$(OBJEXT)
atmos_age_tracer_mod.mod: atmos_age_tracer.$(OBJEXT)
atmos_sulfur_hex_mod.mod: atmos_sulfur_hex.$(OBJEXT)
atmos_sulfate_mod.mod: atmos_sulfate.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
atmos_convection_tracer.$(OBJEXT): atmos_tracer_utilities_mod.mod
atmos_tracer_driver.$(OBJEXT): atmos_tracer_utilities_mod.mod atmos_convection_tracer_mod.mod atmos_co2_mod.mod atmos_ch3i_mod.mod atmos_sea_salt_mod.mod atmos_carbon_aerosol_mod.mod atmos_dust_mod.mod atmos_radon_mod.mod atmos_regional_tracer_driver_mod.mod atmos_soa_mod.mod atmos_age_tracer_mod.mod atmos_sulfur_hex_mod.mod atmos_sulfate_mod.mod
atmos_sea_salt.$(OBJEXT): atmos_tracer_utilities_mod.mod
atmos_carbon_aerosol.$(OBJEXT): atmos_tracer_utilities_mod.mod
atmos_dust.$(OBJEXT): atmos_tracer_utilities_mod.mod
atmos_radon.$(OBJEXT): atmos_tracer_utilities_mod.mod
atmos_sulfur_hex.$(OBJEXT): atmos_tracer_utilities_mod.mod
atmos_sulfate.$(OBJEXT): cloud_chem_mod.mod

# Mod files are built and then installed as headers.
MODFILES = atmos_tracer_utilities_mod.mod atmos_convection_tracer_mod.mod \
atmos_co2_mod.mod cloud_chem_mod.mod atmos_ch3i_mod.mod \
atmos_tracer_driver_mod.mod atmos_sea_salt_mod.mod \
atmos_carbon_aerosol_mod.mod atmos_dust_mod.mod atmos_radon_mod.mod \
atmos_regional_tracer_driver_mod.mod atmos_soa_mod.mod \
atmos_age_tracer_mod.mod atmos_sulfur_hex_mod.mod \
atmos_sulfate_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = atmos_convection_tracer.html atmos_soa.html atmos_tracer_driver.html \
atmos_tracer_utilities.html atmos_sulfur_hex.html atmos_co2.html \
atmos_sea_salt.html atmos_sulfate.html atmos_dust.html \
atmos_carbon_aerosol.html atmos_radon.html

CLEANFILES = *.mod
