# Makefile created by mkmf 19.3.0

CPPDEFS = -DINTERNAL_FILE_NML -g -DCLIMATE_NUDGE -DSPMD 

OTHERFLAGS = -I$(BUILDROOT)atmos_cubed_sphere -I$(BUILDROOT)atmos_phys -I$(BUILDROOT)fms -I$(SRCROOT)FMS/include

include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libatmos_dyn.a
atmos_model.o: $(SRCROOT)atmos_drivers/coupled/atmos_model.F90 $(SRCROOT)FMS/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)atmos_drivers/coupled/atmos_model.F90
./atmos_model.F90: $(SRCROOT)atmos_drivers/coupled/atmos_model.F90
	cp $(SRCROOT)atmos_drivers/coupled/atmos_model.F90 .
./file_version.h: $(SRCROOT)FMS/include/file_version.h
	cp $(SRCROOT)FMS/include/file_version.h .
SRC = $(SRCROOT)atmos_drivers/coupled/atmos_model.F90 $(SRCROOT)FMS/include/file_version.h
OBJ = atmos_model.o
OFF = $(SRCROOT)atmos_drivers/coupled/atmos_model.F90 $(SRCROOT)FMS/include/file_version.h
clean: neat
	-rm -f .libatmos_dyn.a.cppdefs  libatmos_dyn.a *.mod *.o

neat:
	-rm -f $(TMPFILES) *.mod

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libatmos_dyn.a: $(OBJ)
	$(AR) $(ARFLAGS) libatmos_dyn.a $(OBJ)

