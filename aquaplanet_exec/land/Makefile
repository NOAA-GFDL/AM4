# Makefile created by mkmf 19.3.2

CPPDEFS = -DINTERNAL_FILE_NML -Duse_libMPI -Duse_netCDF -DSPMD -DLAND_GRID_FROM_ATMOS -DLAND_BND_TRACERS -DGFDL_COMPATIBLE_MICROP

OTHERFLAGS = -I$(BUILDROOT)fms/build/.mods -I$(SRCROOT)FMS/include 

TMPFILES = *.f90 *.mod *.o

include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libland.a
land_model.o: $(SRCROOT)land_null/land_model.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_null/land_model.F90
./land_model.F90: $(SRCROOT)land_null/land_model.F90
	cp $(SRCROOT)land_null/land_model.F90 .
SRC = $(SRCROOT)land_null/land_model.F90
OBJ = land_model.o
OFF = $(SRCROOT)land_null/land_model.F90
clean: neat
	-rm -f .libland.a.cppdefs $(OBJ) libland.a

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libland.a: $(OBJ)
	$(AR) $(ARFLAGS) libland.a $(OBJ)

