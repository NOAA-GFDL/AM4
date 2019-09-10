# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/clubb.
#
# Ed Hartnett 2019-07-30

# Run make targets in each of these subdirs
SUBDIRS = CLUBB_core

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libadvance_sclrm_Nd_module.la libCLUBB_driver_SCM.la libalt_cloud.la \
libCLUBB_3D_var.la

# Each convenience library depends on its source.
libadvance_sclrm_Nd_module_la_SOURCES = advance_sclrm_Nd_module.F90
libCLUBB_driver_SCM_la_SOURCES = CLUBB_driver_SCM.F90
libalt_cloud_la_SOURCES = alt_cloud.F90
libCLUBB_3D_var_la_SOURCES = CLUBB_3D_var.F90

# Each mod file depends on the .o file.
advance_sclrm_Nd_module_mod.mod: advance_sclrm_Nd_module.$(OBJEXT)
CLUBB_driver_SCM_mod.mod: CLUBB_driver_SCM.$(OBJEXT)
alt_cloud_mod.mod: alt_cloud.$(OBJEXT)
CLUBB_3D_var_mod.mod: CLUBB_3D_var.$(OBJEXT)

# Some mods are dependant on other mods in this dir.

# Mod files are built and then installed as headers.
MODFILES = advance_sclrm_Nd_module_mod.mod CLUBB_driver_SCM_mod.mod \
alt_cloud_mod.mod CLUBB_3D_var_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
