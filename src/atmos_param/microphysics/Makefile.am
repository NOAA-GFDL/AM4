# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/microphysics.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libsimple_pdf.la librotstayn_klein_mp.la libgamma_mg.la \
libls_cloud_microphysics.la libcldwat2m_micro.la \
libmorrison_gettelman_microp.la libmicro_mg.la

# Each convenience library depends on its source.
libsimple_pdf_la_SOURCES = simple_pdf.F90
librotstayn_klein_mp_la_SOURCES = rotstayn_klein_mp.F90
libgamma_mg_la_SOURCES = gamma_mg.F90
libls_cloud_microphysics_la_SOURCES = ls_cloud_microphysics.F90
libcldwat2m_micro_la_SOURCES = cldwat2m_micro.F90
libmorrison_gettelman_microp_la_SOURCES = morrison_gettelman_microp.F90
libmicro_mg_la_SOURCES = micro_mg.F90

# Each mod file depends on the .o file.
simple_pdf_mod.mod: simple_pdf.$(OBJEXT)
rotstayn_klein_mp_mod.mod: rotstayn_klein_mp.$(OBJEXT)
gamma_mg_mod.mod: gamma_mg.$(OBJEXT)
ls_cloud_microphysics_mod.mod: ls_cloud_microphysics.$(OBJEXT)
cldwat2m_micro_mod.mod: cldwat2m_micro.$(OBJEXT)
morrison_gettelman_microp_mod.mod: morrison_gettelman_microp.$(OBJEXT)
micro_mg_mod.mod: micro_mg.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
ls_cloud_microphysics.$(OBJEXT): rotstayn_klein_mp_mod.mod cldwat2m_micro_mod.mod morrison_gettelman_microp_mod.mod micro_mg_mod.mod
cldwat2m_micro.$(OBJEXT): simple_pdf_mod.mod gamma_mg_mod.mod
morrison_gettelman_microp.$(OBJEXT): simple_pdf_mod.mod gamma_mg_mod.mod
micro_mg.$(OBJEXT): simple_pdf_mod.mod gamma_mg_mod.mod cldwat2m_micro_mod.mod

# Mod files are built and then installed as headers.
MODFILES = simple_pdf_mod.mod rotstayn_klein_mp_mod.mod gamma_mg_mod.mod \
ls_cloud_microphysics_mod.mod cldwat2m_micro_mod.mod \
morrison_gettelman_microp_mod.mod micro_mg_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
