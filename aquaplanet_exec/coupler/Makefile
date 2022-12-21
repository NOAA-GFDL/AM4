# Makefile created by mkmf 19.3.2

OTHERFLAGS = -I$(BUILDROOT)atmos_dyn -I$(BUILDROOT)atmos_cubed_sphere -I$(BUILDROOT)ice -I$(BUILDROOT)land -I$(BUILDROOT)atmos_phys -I$(BUILDROOT)fms/build/.mods -I$(SRCROOT)FMS/include 

TMPFILES = *.f90 *.mod *.o

include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libcoupler.a
coupler_main.o: $(SRCROOT)FMScoupler/simple/coupler_main.F90 ice_model.o flux_exchange.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMScoupler/simple/coupler_main.F90
flux_exchange.o: $(SRCROOT)FMScoupler/simple/flux_exchange.F90 ice_model.o surface_flux.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMScoupler/simple/flux_exchange.F90
ice_model.o: $(SRCROOT)FMScoupler/simple/ice_model.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMScoupler/simple/ice_model.F90
surface_flux.o: $(SRCROOT)FMScoupler/shared/surface_flux.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMScoupler/shared/surface_flux.F90
./surface_flux.F90: $(SRCROOT)FMScoupler/shared/surface_flux.F90
	cp $(SRCROOT)FMScoupler/shared/surface_flux.F90 .
./flux_exchange.F90: $(SRCROOT)FMScoupler/simple/flux_exchange.F90
	cp $(SRCROOT)FMScoupler/simple/flux_exchange.F90 .
./coupler_main.F90: $(SRCROOT)FMScoupler/simple/coupler_main.F90
	cp $(SRCROOT)FMScoupler/simple/coupler_main.F90 .
./ice_model.F90: $(SRCROOT)FMScoupler/simple/ice_model.F90
	cp $(SRCROOT)FMScoupler/simple/ice_model.F90 .
SRC = $(SRCROOT)FMScoupler/simple/ice_model.F90 $(SRCROOT)FMScoupler/simple/flux_exchange.F90 $(SRCROOT)FMScoupler/shared/surface_flux.F90 $(SRCROOT)FMScoupler/simple/coupler_main.F90
OBJ = ice_model.o flux_exchange.o surface_flux.o coupler_main.o
OFF = $(SRCROOT)FMScoupler/shared/surface_flux.F90 $(SRCROOT)FMScoupler/simple/flux_exchange.F90 $(SRCROOT)FMScoupler/simple/coupler_main.F90 $(SRCROOT)FMScoupler/simple/ice_model.F90
clean: neat
	-rm -f .libcoupler.a.cppdefs $(OBJ) libcoupler.a 

neat:
	-rm -f $(TMPFILES) *.mod

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libcoupler.a: $(OBJ)
	$(AR) $(ARFLAGS) libcoupler.a $(OBJ)

