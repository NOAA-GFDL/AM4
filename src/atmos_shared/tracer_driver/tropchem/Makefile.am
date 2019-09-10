# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_shared/tracer_driver/tropchem.
#
# Ed Hartnett 2019-07-30

# Run make targets in each of these subdirs
SUBDIRS = AM3 AM4

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libtropchem_types_mod.la libstrat_chem_utilities.la libmo_photo.la \
libtropchem_driver.la libmo_setinv.la libaerosol_thermodynamics.la \
libmo_usrrxt.la libm_tracname.la libmo_imp_slv.la libmo_exp_slv.la \
libmo_hook.la libmo_setsox.la libmo_read_sim_chm.la libmo_rodas_slv.la \
libmo_chem_utls.la libmo_jpl.la libmo_chemdr.la libmo_chemini.la

# Each convenience library depends on its source.
libtropchem_types_mod_la_SOURCES = tropchem_types_mod.F90
libstrat_chem_utilities_la_SOURCES = strat_chem_utilities.F90
libmo_photo_la_SOURCES = mo_photo.F90
libtropchem_driver_la_SOURCES = tropchem_driver.F90
libmo_setinv_la_SOURCES = mo_setinv.F90
libaerosol_thermodynamics_la_SOURCES = aerosol_thermodynamics.F90
libmo_usrrxt_la_SOURCES = mo_usrrxt.F90
libm_tracname_la_SOURCES = m_tracname.F90
libmo_imp_slv_la_SOURCES = mo_imp_slv.F90
libmo_exp_slv_la_SOURCES = mo_exp_slv.F90
libmo_hook_la_SOURCES = mo_hook.F90
libmo_setsox_la_SOURCES = mo_setsox.F90
libmo_read_sim_chm_la_SOURCES = mo_read_sim_chm.F90
libmo_rodas_slv_la_SOURCES = mo_rodas_slv.F90
libmo_chem_utls_la_SOURCES = mo_chem_utls.F90
libmo_jpl_la_SOURCES = mo_jpl.F90
libmo_chemdr_la_SOURCES = mo_chemdr.F90
libmo_chemini_la_SOURCES = mo_chemini.F90

# Each mod file depends on the .o file.
tropchem_types_mod_mod.mod: tropchem_types_mod.$(OBJEXT)
strat_chem_utilities_mod.mod: strat_chem_utilities.$(OBJEXT)
mo_photo_mod.mod: mo_photo.$(OBJEXT)
tropchem_driver_mod.mod: tropchem_driver.$(OBJEXT)
mo_setinv_mod.mod: mo_setinv.$(OBJEXT)
aerosol_thermodynamics_mod.mod: aerosol_thermodynamics.$(OBJEXT)
mo_usrrxt_mod.mod: mo_usrrxt.$(OBJEXT)
m_tracname_mod.mod: m_tracname.$(OBJEXT)
mo_imp_slv_mod.mod: mo_imp_slv.$(OBJEXT)
mo_exp_slv_mod.mod: mo_exp_slv.$(OBJEXT)
mo_hook_mod.mod: mo_hook.$(OBJEXT)
mo_setsox_mod.mod: mo_setsox.$(OBJEXT)
mo_read_sim_chm_mod.mod: mo_read_sim_chm.$(OBJEXT)
mo_rodas_slv_mod.mod: mo_rodas_slv.$(OBJEXT)
mo_chem_utls_mod.mod: mo_chem_utls.$(OBJEXT)
mo_jpl_mod.mod: mo_jpl.$(OBJEXT)
mo_chemdr_mod.mod: mo_chemdr.$(OBJEXT)
mo_chemini_mod.mod: mo_chemini.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
mo_photo.$(OBJEXT): mo_chem_utls_mod.mod
tropchem_driver.$(OBJEXT): aerosol_thermodynamics_mod.mod mo_usrrxt_mod.mod
aerosol_thermodynamics.$(OBJEXT): tropchem_types_mod_mod.mod
mo_usrrxt.$(OBJEXT): tropchem_types_mod_mod.mod strat_chem_utilities_mod.mod mo_chem_utls_mod.mod
mo_imp_slv.$(OBJEXT): m_tracname_mod.mod mo_chem_utls_mod.mod
mo_exp_slv.$(OBJEXT): mo_chem_utls_mod.mod
mo_setsox.$(OBJEXT): tropchem_types_mod_mod.mod aerosol_thermodynamics_mod.mod mo_chem_utls_mod.mod
mo_rodas_slv.$(OBJEXT): mo_chem_utls_mod.mod
mo_chem_utls.$(OBJEXT): m_tracname_mod.mod
mo_chemdr.$(OBJEXT): tropchem_types_mod_mod.mod strat_chem_utilities_mod.mod mo_photo_mod.mod mo_setinv_mod.mod mo_usrrxt_mod.mod mo_setsox_mod.mod mo_chem_utls_mod.mod
mo_chemini.$(OBJEXT): tropchem_types_mod_mod.mod mo_usrrxt_mod.mod mo_chem_utls_mod.mod

# Mod files are built and then installed as headers.
MODFILES = tropchem_types_mod_mod.mod strat_chem_utilities_mod.mod \
mo_photo_mod.mod tropchem_driver_mod.mod mo_setinv_mod.mod \
aerosol_thermodynamics_mod.mod mo_usrrxt_mod.mod m_tracname_mod.mod \
mo_imp_slv_mod.mod mo_exp_slv_mod.mod mo_hook_mod.mod \
mo_setsox_mod.mod mo_read_sim_chm_mod.mod mo_rodas_slv_mod.mod \
mo_chem_utls_mod.mod mo_jpl_mod.mod mo_chemdr_mod.mod \
mo_chemini_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = tropchem_driver.html isrpia.inc

CLEANFILES = *.mod
