# Makefile created by mkmf 19.3.0

CPPDEFS = -DINTERNAL_FILE_NML

OTHERFLAGS = -I$(BUILDROOT)atmos_dyn -I$(BUILDROOT)ice_sis -I$(BUILDROOT)atmos_phys -I$(BUILDROOT)mom6 -I$(BUILDROOT)land_lad2 -I$(BUILDROOT)fms

include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libcoupler.a
atmos_ocean_dep_fluxes_calc.o: $(SRCROOT)coupler/atmos_ocean_dep_fluxes_calc.f90 surface_flux.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c $(SRCROOT)coupler/atmos_ocean_dep_fluxes_calc.f90
atmos_ocean_fluxes_calc.o: $(SRCROOT)coupler/atmos_ocean_fluxes_calc.f90 surface_flux.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c $(SRCROOT)coupler/atmos_ocean_fluxes_calc.f90
atm_land_ice_flux_exchange.o: $(SRCROOT)coupler/atm_land_ice_flux_exchange.F90 surface_flux.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)coupler/atm_land_ice_flux_exchange.F90
coupler_main.o: $(SRCROOT)coupler/coupler_main.F90 flux_exchange.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)coupler/coupler_main.F90
flux_exchange.o: $(SRCROOT)coupler/flux_exchange.F90 atm_land_ice_flux_exchange.o land_ice_flux_exchange.o ice_ocean_flux_exchange.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)coupler/flux_exchange.F90
ice_ocean_flux_exchange.o: $(SRCROOT)coupler/ice_ocean_flux_exchange.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)coupler/ice_ocean_flux_exchange.F90
land_ice_flux_exchange.o: $(SRCROOT)coupler/land_ice_flux_exchange.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)coupler/land_ice_flux_exchange.F90
surface_flux.o: $(SRCROOT)coupler/surface_flux.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)coupler/surface_flux.F90

SRC = $(SRCROOT)coupler/land_ice_flux_exchange.F90 $(SRCROOT)coupler/surface_flux.F90 $(SRCROOT)coupler/ice_ocean_flux_exchange.F90 $(SRCROOT)coupler/coupler_main.F90 $(SRCROOT)coupler/atm_land_ice_flux_exchange.F90 $(SRCROOT)coupler/flux_exchange.F90 $(SRCROOT)coupler/atmos_ocean_fluxes_calc.f90 $(SRCROOT)coupler/atmos_ocean_dep_fluxes_calc.f90
OBJ = land_ice_flux_exchange.o surface_flux.o ice_ocean_flux_exchange.o coupler_main.o atm_land_ice_flux_exchange.o flux_exchange.o atmos_ocean_fluxes_calc.o atmos_ocean_dep_fluxes_calc.o

clean: neat
	$(RM) $(OBJ) libcoupler.a
	$(RM) *.mod

neat:
	$(RM) $(TMPFILES)

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libcoupler.a: $(OBJ)
	$(AR) $(ARFLAGS) libcoupler.a $(OBJ)

