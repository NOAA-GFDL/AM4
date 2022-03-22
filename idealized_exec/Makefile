# Makefile for Experiment 'aquaplanet_2022.01'

BUILDROOT = $(dir $(abspath $(firstword $(MAKEFILE_LIST))))
MKTEMPLATEROOT = $(abspath $(BUILDROOT)../exec/templates)/
SRCROOT = $(abspath $(BUILDROOT)../src)/ # Final '/' added to denote directory

# Select the make template (default is intel)
ifeq ($(gcc),on)
  MK_TEMPLATE = $(MKTEMPLATEROOT)gnu.mk
else
  MK_TEMPLATE = $(MKTEMPLATEROOT)intel.mk
endif
# Set shell to sh
ifeq ($(SH),sh)
  SHELL=/bin/sh
endif

# Override default compile options. See templates/intel.mk for more compile options.
# BLD_TYPE configures the type of build, PROD (production), REPRO (regression) or
# DEBUG (debug options)
BLD_TYPE = PROD
OPENMP = t # Build with OpenMP enabled

include $(MK_TEMPLATE)

SUBMAKEFLAGS = BUILDROOT=$(BUILDROOT) SRCROOT=$(SRCROOT) MK_TEMPLATE=$(MK_TEMPLATE) BLD_TYPE=$(BLD_TYPE)

aquaplanet_2022.01.x: coupler/libcoupler.a atmos_dyn/libatmos_dyn.a ice/libice.a atmos_cubed_sphere/libatmos_cubed_sphere.a atmos_phys/libatmos_phys.a land/libland.a fms/build/libFMS/.libs/libFMS.a
	$(LD) $^ $(LDFLAGS) -o $@ $(STATIC_LIBS)

fms/build/libFMS/.libs/libFMS.a:  FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)fms $(@F) 

atmos_phys/libatmos_phys.a: fms/build/libFMS/.libs/libFMS.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)atmos_phys $(@F)

atmos_cubed_sphere/libatmos_cubed_sphere.a: atmos_phys/libatmos_phys.a fms/build/libFMS/.libs/libFMS.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)atmos_cubed_sphere $(@F)

atmos_dyn/libatmos_dyn.a: atmos_cubed_sphere/libatmos_cubed_sphere.a atmos_phys/libatmos_phys.a fms/build/libFMS/.libs/libFMS.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)atmos_dyn $(@F)

land/libland.a: fms/build/libFMS/.libs/libFMS.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)land $(@F)

ice/libice.a: fms/build/libFMS/.libs/libFMS.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)ice $(@F)

coupler/libcoupler.a: atmos_dyn/libatmos_dyn.a ice/libice.a land/libland.a atmos_phys/libatmos_phys.a fms/build/libFMS/.libs/libFMS.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)coupler $(@F)

FORCE:

clean:
	$(MAKE) --directory=$(BUILDROOT)fms clean
	$(MAKE) --directory=$(BUILDROOT)atmos_phys clean
	$(MAKE) --directory=$(BUILDROOT)atmos_dyn clean
	$(MAKE) --directory=$(BUILDROOT)atmos_cubed_sphere clean
	$(MAKE) --directory=$(BUILDROOT)ice clean
	$(MAKE) --directory=$(BUILDROOT)land clean
	$(MAKE) --directory=$(BUILDROOT)coupler clean

clean_all: clean
	$(RM) *.x

