# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/convection_driver.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libmoistproc_kernels.la libdetr_ice_num.la libconvection_driver.la

# Each convenience library depends on its source.
libmoistproc_kernels_la_SOURCES = moistproc_kernels.F90
libdetr_ice_num_la_SOURCES = detr_ice_num.F90
libconvection_driver_la_SOURCES = convection_driver.F90

# Each mod file depends on the .o file.
moistproc_kernels_mod.mod: moistproc_kernels.$(OBJEXT)
detr_ice_num_mod.mod: detr_ice_num.$(OBJEXT)
convection_driver_mod.mod: convection_driver.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
moistproc_kernels.$(OBJEXT): detr_ice_num_mod.mod
convection_driver.$(OBJEXT): moistproc_kernels_mod.mod detr_ice_num_mod.mod

# Mod files are built and then installed as headers.
MODFILES = moistproc_kernels_mod.mod detr_ice_num_mod.mod \
convection_driver_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
