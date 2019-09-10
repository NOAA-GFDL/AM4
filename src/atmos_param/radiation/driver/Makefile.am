# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/radiation/driver.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libsolar_data_driver.la libradiation_driver_types.la \
libradiation_driver_diag.la librad_output_file.la \
libradiation_driver.la liblongwave_driver.la libradiation_types.la \
libshortwave_driver.la

# Each convenience library depends on its source.
libsolar_data_driver_la_SOURCES = solar_data_driver.F90
libradiation_driver_types_la_SOURCES = radiation_driver_types.F90
libradiation_driver_diag_la_SOURCES = radiation_driver_diag.F90
librad_output_file_la_SOURCES = rad_output_file.F90
libradiation_driver_la_SOURCES = radiation_driver.F90
liblongwave_driver_la_SOURCES = longwave_driver.F90
libradiation_types_la_SOURCES = radiation_types.F90
libshortwave_driver_la_SOURCES = shortwave_driver.F90

# Each mod file depends on the .o file.
solar_data_driver_mod.mod: solar_data_driver.$(OBJEXT)
radiation_driver_types_mod.mod: radiation_driver_types.$(OBJEXT)
radiation_driver_diag_mod.mod: radiation_driver_diag.$(OBJEXT)
rad_output_file_mod.mod: rad_output_file.$(OBJEXT)
radiation_driver_mod.mod: radiation_driver.$(OBJEXT)
longwave_driver_mod.mod: longwave_driver.$(OBJEXT)
radiation_types_mod.mod: radiation_types.$(OBJEXT)
shortwave_driver_mod.mod: shortwave_driver.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
radiation_driver_diag.$(OBJEXT): radiation_driver_types_mod.mod radiation_driver_mod.mod shortwave_driver_mod.mod
rad_output_file.$(OBJEXT): radiation_driver_types_mod.mod radiation_driver_mod.mod
radiation_driver.$(OBJEXT): radiation_driver_types_mod.mod radiation_driver_diag_mod.mod rad_output_file_mod.mod longwave_driver_mod.mod radiation_types_mod.mod shortwave_driver_mod.mod
longwave_driver.$(OBJEXT): radiation_driver_types_mod.mod radiation_driver_mod.mod
shortwave_driver.$(OBJEXT): solar_data_driver_mod.mod radiation_driver_types_mod.mod radiation_driver_mod.mod

# Mod files are built and then installed as headers.
MODFILES = solar_data_driver_mod.mod radiation_driver_types_mod.mod \
radiation_driver_diag_mod.mod rad_output_file_mod.mod \
radiation_driver_mod.mod longwave_driver_mod.mod \
radiation_types_mod.mod shortwave_driver_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
