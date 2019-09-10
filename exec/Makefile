# Makefile for Experiment 'cm4p12_warsaw'

BUILDROOT = $(dir $(abspath $(firstword $(MAKEFILE_LIST))))
SRCROOT = $(abspath $(BUILDROOT)../src)/ # Final '/' added to denote directory

MK_TEMPLATE = $(BUILDROOT)templates/intel.mk

# Override default compile options. See templates/intel.mk for more compile options.
# BLD_TYPE configures the type of build, PROD (production), REPRO (regression) or
# DEBUG (debug options)
BLD_TYPE = PROD
OPENMP = t # Build with OpenMP enabled

include $(MK_TEMPLATE)

SUBMAKEFLAGS = BUILDROOT=$(BUILDROOT) SRCROOT=$(SRCROOT) MK_TEMPLATE=$(MK_TEMPLATE) BLD_TYPE=$(BLD_TYPE)

fms_cm4p12_warsaw.x: coupler/libcoupler.a ice_sis/libice_sis.a atmos_dyn/libatmos_dyn.a land_lad2/libland_lad2.a atmos_phys/libatmos_phys.a mom6/libmom6.a fms/libfms.a
	$(LD) $^ $(LDFLAGS) -o $@ $(STATIC_LIBS)

fms/libfms.a:  FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)fms $(@F)

atmos_phys/libatmos_phys.a: fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)atmos_phys $(@F)

atmos_dyn/libatmos_dyn.a: atmos_phys/libatmos_phys.a fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)atmos_dyn $(@F)

ice_sis/libice_sis.a: mom6/libmom6.a fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)ice_sis $(@F)

land_lad2/libland_lad2.a: fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)land_lad2 $(@F)

# MOM6 is currently unable to run with OpenMP enabled
mom6/libmom6.a: fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) --directory=$(BUILDROOT)mom6 $(@F)

coupler/libcoupler.a: atmos_dyn/libatmos_dyn.a ice_sis/libice_sis.a atmos_phys/libatmos_phys.a mom6/libmom6.a land_lad2/libland_lad2.a fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)coupler $(@F)

FORCE:

clean:
	$(MAKE) --directory=$(BUILDROOT)fms clean
	$(MAKE) --directory=$(BUILDROOT)atmos_phys clean
	$(MAKE) --directory=$(BUILDROOT)atmos_dyn clean
	$(MAKE) --directory=$(BUILDROOT)ice_sis clean
	$(MAKE) --directory=$(BUILDROOT)land_lad2 clean
	$(MAKE) --directory=$(BUILDROOT)mom6 clean
	$(MAKE) --directory=$(BUILDROOT)coupler clean

clean_all: clean
	$(RM) *.x

