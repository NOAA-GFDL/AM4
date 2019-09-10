# This is an automake file for the GFDL AM4 model.
# This builds directory ./src/land_lad2/shared.
#
# Ed Hartnett 2019-07-30

# Run make targets in each of these subdirs
SUBDIRS = nf_utils

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libland_utils.la libsphum.la libland_debug.la libland_numerics.la \
libland_tile_diag.la libtable_printer.la libland_tile_diag_buff.la \
libland_tile_io.la libland_tile_diag_sel.la libland_io.la

# Each convenience library depends on its source.
libland_utils_la_SOURCES = land_utils.F90
libsphum_la_SOURCES = sphum.F90
libland_debug_la_SOURCES = land_debug.F90
libland_numerics_la_SOURCES = land_numerics.F90
libland_tile_diag_la_SOURCES = land_tile_diag.F90
libtable_printer_la_SOURCES = table_printer.F90
libland_tile_diag_buff_la_SOURCES = land_tile_diag_buff.F90
libland_tile_io_la_SOURCES = land_tile_io.F90
libland_tile_diag_sel_la_SOURCES = land_tile_diag_sel.F90
libland_io_la_SOURCES = land_io.F90

# Each mod file depends on the .o file.
land_utils_mod.mod: land_utils.$(OBJEXT)
sphum_mod.mod: sphum.$(OBJEXT)
land_debug_mod.mod: land_debug.$(OBJEXT)
land_numerics_mod.mod: land_numerics.$(OBJEXT)
land_tile_diag_mod.mod: land_tile_diag.$(OBJEXT)
table_printer_mod.mod: table_printer.$(OBJEXT)
land_tile_diag_buff_mod.mod: land_tile_diag_buff.$(OBJEXT)
land_tile_io_mod.mod: land_tile_io.$(OBJEXT)
land_tile_diag_sel_mod.mod: land_tile_diag_sel.$(OBJEXT)
land_io_mod.mod: land_io.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
land_tile_diag.$(OBJEXT): land_debug_mod.mod
land_tile_io.$(OBJEXT): land_utils_mod.mod land_io_mod.mod
land_tile_diag_sel.$(OBJEXT): land_tile_diag_mod.mod
land_io.$(OBJEXT): land_numerics_mod.mod

# Mod files are built and then installed as headers.
MODFILES = land_utils_mod.mod sphum_mod.mod land_debug_mod.mod \
land_numerics_mod.mod land_tile_diag_mod.mod table_printer_mod.mod \
land_tile_diag_buff_mod.mod land_tile_io_mod.mod \
land_tile_diag_sel_mod.mod land_io_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = debug.inc version_variable.inc

CLEANFILES = *.mod
