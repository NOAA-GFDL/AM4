# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/atmos_param/cosp/quickbeam.
#
# Ed Hartnett 2019-07-30

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libradar_simulator_types.la libradar_simulator.la libscale_LUTs_io.la \
libradar_simulator_init.la libatmos_lib.la libformat_input.la \
libmrgrnk.la libgases.la libdsd.la libmath_lib.la libarray_lib.la \
liboptics_lib.la libcalc_Re.la libzeff.la

# Each convenience library depends on its source.
libradar_simulator_types_la_SOURCES = radar_simulator_types.F90
libradar_simulator_la_SOURCES = radar_simulator.F90
libscale_LUTs_io_la_SOURCES = scale_LUTs_io.F90
libradar_simulator_init_la_SOURCES = radar_simulator_init.F90
libatmos_lib_la_SOURCES = atmos_lib.F90
libformat_input_la_SOURCES = format_input.F90
libmrgrnk_la_SOURCES = mrgrnk.F90
libgases_la_SOURCES = gases.F90
libdsd_la_SOURCES = dsd.F90
libmath_lib_la_SOURCES = math_lib.F90
libarray_lib_la_SOURCES = array_lib.F90
liboptics_lib_la_SOURCES = optics_lib.F90
libcalc_Re_la_SOURCES = calc_Re.F90
libzeff_la_SOURCES = zeff.F90

# Each mod file depends on the .o file.
radar_simulator_types_mod.mod: radar_simulator_types.$(OBJEXT)
radar_simulator_mod.mod: radar_simulator.$(OBJEXT)
scale_LUTs_io_mod.mod: scale_LUTs_io.$(OBJEXT)
radar_simulator_init_mod.mod: radar_simulator_init.$(OBJEXT)
atmos_lib_mod.mod: atmos_lib.$(OBJEXT)
format_input_mod.mod: format_input.$(OBJEXT)
mrgrnk_mod.mod: mrgrnk.$(OBJEXT)
gases_mod.mod: gases.$(OBJEXT)
dsd_mod.mod: dsd.$(OBJEXT)
math_lib_mod.mod: math_lib.$(OBJEXT)
array_lib_mod.mod: array_lib.$(OBJEXT)
optics_lib_mod.mod: optics_lib.$(OBJEXT)
calc_Re_mod.mod: calc_Re.$(OBJEXT)
zeff_mod.mod: zeff.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
radar_simulator.$(OBJEXT): radar_simulator_types_mod.mod scale_LUTs_io_mod.mod math_lib_mod.mod array_lib_mod.mod optics_lib_mod.mod
scale_LUTs_io.$(OBJEXT): radar_simulator_types_mod.mod radar_simulator_mod.mod
radar_simulator_init.$(OBJEXT): radar_simulator_types_mod.mod radar_simulator_mod.mod
format_input.$(OBJEXT): array_lib_mod.mod
dsd.$(OBJEXT): math_lib_mod.mod array_lib_mod.mod
math_lib.$(OBJEXT): array_lib_mod.mod
calc_Re.$(OBJEXT): math_lib_mod.mod
zeff.$(OBJEXT): math_lib_mod.mod optics_lib_mod.mod

# Mod files are built and then installed as headers.
MODFILES = radar_simulator_types_mod.mod radar_simulator_mod.mod \
scale_LUTs_io_mod.mod radar_simulator_init_mod.mod atmos_lib_mod.mod \
format_input_mod.mod mrgrnk_mod.mod gases_mod.mod dsd_mod.mod \
math_lib_mod.mod array_lib_mod.mod optics_lib_mod.mod calc_Re_mod.mod \
zeff_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
