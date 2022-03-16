# Makefile created by mkmf 19.3.2

CPPDEFS = -Duse_libMPI -Duse_netCDF -DINTERNAL_FILE_NML -DSPMD -DLAND_GRID_FROM_ATMOS -DLAND_BND_TRACERS -DGFDL_COMPATIBLE_MICROP

OTHERFLAGS = -I$(BUILDROOT)fms/build/.mods -I$(SRCROOT)FMS/include 

TMPFILES = *.f90 *.mod *.o

include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libice.a
ice_albedo.o: $(SRCROOT)ice_param/ice_albedo.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_param/ice_albedo.F90
ocean_albedo.o: $(SRCROOT)ice_param/ocean_albedo.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_param/ocean_albedo.F90
ocean_rough.o: $(SRCROOT)ice_param/ocean_rough.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_param/ocean_rough.F90
./ice_albedo.F90: $(SRCROOT)ice_param/ice_albedo.F90
	cp $(SRCROOT)ice_param/ice_albedo.F90 .
./ocean_albedo.F90: $(SRCROOT)ice_param/ocean_albedo.F90
	cp $(SRCROOT)ice_param/ocean_albedo.F90 .
./ocean_rough.F90: $(SRCROOT)ice_param/ocean_rough.F90
	cp $(SRCROOT)ice_param/ocean_rough.F90 .
SRC = $(SRCROOT)ice_param/ocean_rough.F90 $(SRCROOT)ice_param/ice_albedo.F90 $(SRCROOT)ice_param/ocean_albedo.F90
OBJ = ocean_rough.o ice_albedo.o ocean_albedo.o
OFF = $(SRCROOT)ice_param/ice_albedo.F90 $(SRCROOT)ice_param/ocean_albedo.F90 $(SRCROOT)ice_param/ocean_rough.F90
clean: neat
	-rm -f .libice.a.cppdefs $(OBJ) libice.a

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libice.a: $(OBJ)
	$(AR) $(ARFLAGS) libice.a $(OBJ)

