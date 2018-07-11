#FC = ifort
#CC = icc
#CXX = CC
#LD = ifort

 FC = mpiifort
 CC = mpiicc
 CXX = CC
 LD = mpiifort


DEBUG =
REPRO =
VERBOSE =
OPENMP =

need := 3.81
ok := $(filter $(need),$(firstword $(sort $(MAKE_VERSION) $(need))))
ifneq ($(need),$(ok))
$(error Need at least make version $(need).  Load module gmake/3.81)
endif

MAKEFLAGS += --jobs=2

NETCDF_DIR = /opt/netcdf/4.6.1/INTEL
NETCDF_ROOT = $(NETCDF_DIR)
HDF5_DIR = /opt/hdf5/1.10.1/INTEL
HDF5_ROOT = $(HDF5_DIR)
#MPI_ROOT    = /opt/intel/2017_up1/impi/2017.1.132/include64
MPI_ROOT    = /opt/intel/2017_up2/impi/2017.2.174/mic/
INCLUDE = -I$(NETCDF_ROOT)/include -I$(HDF5_ROOT)/include
#-L/usr/lib64 -lhdf5_cpp -lhdf5_fortran -lhdf5_hl -lhdf5_hl_cpp -lhdf5hl_fortran -lhdf5_hl -lhdf5
#-L/home/Thomas.Robinson/hdf5/lib
FPPFLAGS := -fpp -Wp,-w $(INCLUDE)

FFLAGS := -msse2 -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -i4 -r8 -nowarn -sox -traceback $(INCLUDE)
FFLAGS_OPT = -O3 -fp-model source -qoverride-limits
#-debug minimal -fp-model source -qoverride-limits
FFLAGS_DEBUG = -g -O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fpe0 -ftrapuv
FFLAGS_REPRO = -O2 -debug minimal -fp-model source -I$(MPI_ROOT)/include
FFLAGS_OPENMP = -qopenmp
#-L/home/Thomas.Robinson/hdf5/lib -lhdf5_cpp -lhdf5_fortran -lhdf5_hl -lhdf5_hl_cpp -lhdf5hl_fortran -lhdf5_hl -lhdf5 -I/home/Thomas.Robinson/hdf5/include -I/opt/intel/2017_up2/advisor_2017.1.2.501009/include/intel64 -L/opt/intel/2017_up2/advisor_2017.1.2.501009/lib64 -Bdynamic -shared-intel -g
#-I/home/Thomas.Robinson/hdf5/include -L/home/Thomas.Robinson/hdf5/lib
#-L/opt/hdf5/1.10.0-patch1/lib
#-lhfd5 -lhdf5_fortran
FFLAGS_VERBOSE = -v -V -what -warn all

CFLAGS := -D__IFC -msse2 -sox -traceback -g $(INCLUDE) -I$(MPI_ROOT)/include
CFLAGS_OPT = -O2 -debug minimal
CFLAGS_REPRO = -O2 -debug minimal $(INCLUDE) -I$(MPI_ROOT)/include
CFLAGS_OPENMP = -qopenmp -I/opt/intel/2017_up2/advisor_2017.1.2.501009/include -g
CFLAGS_DEBUG = -O0 -g -ftrapuv $(INCLUDE) -I$(MPI_ROOT)/include
CFLAGS_VERBOSE = -w3

FFLAGS_TEST = -O3 -debug minimal -fp-model source -qoverride-limits
CFLAGS_TEST = -O2

LDFLAGS := -L$(HDF5_ROOT)/lib -lhdf5_cpp -lhdf5_fortran -lhdf5_hl -lhdf5_hl_cpp -lhdf5hl_fortran -lhdf5_hl -lhdf5
LDFLAGS_OPENMP := -qopenmp
#-L/home/Thomas.Robinson/hdf5/lib -lhdf5_cpp -lhdf5_fortran -lhdf5_hl -lhdf5_hl_cpp -lhdf5hl_fortran -lhdf5_hl -lhdf5
#-lhfd5 -lhdf5_fortran
LDFLAGS_VERBOSE := -Wl,-V,--verbose,-cref,-M

LIBS :=

ifneq ($(REPRO),)
CFLAGS += $(CFLAGS_REPRO)
FFLAGS += $(FFLAGS_REPRO)
else ifneq ($(DEBUG),)
CFLAGS += $(CFLAGS_DEBUG)
FFLAGS += $(FFLAGS_DEBUG)
else ifneq ($(TEST),)
CFLAGS += $(CFLAGS_TEST)
FFLAGS += $(FFLAGS_TEST)
else
CFLAGS += $(CFLAGS_OPT)
FFLAGS += $(FFLAGS_OPT)
endif

ifneq ($(OPENMP),)
CFLAGS += $(CFLAGS_OPENMP)
FFLAGS += $(FFLAGS_OPENMP)
LDFLAGS += $(LDFLAGS_OPENMP)

LIBS += -L$(INTEL_PATH)/$(INTEL_MAJOR_VERSION)/$(INTEL_MINOR_VERSION)/lib/intel64 -lifcoremt
endif

ifneq ($(VERBOSE),)
CFLAGS += $(CFLAGS_VERBOSE)
FFLAGS += $(FFLAGS_VERBOSE)
LDFLAGS += $(LDFLAGS_VERBOSE)
endif

ifeq ($(NETCDF),3)
  ifneq ($(findstring -Duse_netCDF,$(CPPDEFS)),)
    CPPDEFS += -Duse_LARGEFILE
  endif
endif

ifneq ($(findstring netcdf-4.0.1,$(LOADEDMODULES)),)
  LIBS += -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz
else
  LIBS += -lnetcdf -lnetcdff -lhdf5_hl -lhdf5 -lz
endif

LIBS +=
LDFLAGS += $(LIBS)

RM = rm -f
SHELL = /bin/csh -f
TMPFILES = .*.m *.B *.L *.i *.i90 *.l *.s *.mod *.opt

.SUFFIXES: .F .F90 .H .L .T .f .f90 .h .i .i90 .l .o .s .opt .x

.f.L:
	$(FC) $(FFLAGS) -c -listing $*.f
.f.opt:
	$(FC) $(FFLAGS) -c -opt_report_level max -opt_report_phase all -opt_report_file $*.opt $*.f
.f.l:
	$(FC) $(FFLAGS) -c $(LIST) $*.f
.f.T:
	$(FC) $(FFLAGS) -c -cif $*.f
.f.o:
	$(FC) $(FFLAGS) -c $*.f
.f.s:
	$(FC) $(FFLAGS) -S $*.f
.f.x:
	$(FC) $(FFLAGS) -o $*.x $*.f *.o $(LDFLAGS)
.f90.L:
	$(FC) $(FFLAGS) -c -listing $*.f90
.f90.opt:
	$(FC) $(FFLAGS) -c -opt_report_level max -opt_report_phase all -opt_report_file $*.opt $*.f90
.f90.l:
	$(FC) $(FFLAGS) -c $(LIST) $*.f90
.f90.T:
	$(FC) $(FFLAGS) -c -cif $*.f90
.f90.o:
	$(FC) $(FFLAGS) -c $*.f90
.f90.s:
	$(FC) $(FFLAGS) -c -S $*.f90
.f90.x:
	$(FC) $(FFLAGS) -o $*.x $*.f90 *.o $(LDFLAGS)
.F.L:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -listing $*.F
.F.opt:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -opt_report_level max -opt_report_phase all -opt_report_file $*.opt $*.F
.F.l:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c $(LIST) $*.F
.F.T:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -cif $*.F
.F.f:
	$(FC) $(CPPDEFS) $(FPPFLAGS) -EP $*.F > $*.f
.F.i:
	$(FC) $(CPPDEFS) $(FPPFLAGS) -P $*.F
.F.o:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c $*.F
.F.s:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -S $*.F
.F.x:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -o $*.x $*.F *.o $(LDFLAGS)
.F90.L:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -listing $*.F90
.F90.opt:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -opt_report_level max -opt_report_phase all -opt_report_file $*.opt $*.F90
.F90.l:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c $(LIST) $*.F90
.F90.T:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -cif $*.F90
.F90.f90:
	$(FC) $(CPPDEFS) $(FPPFLAGS) -EP $*.F90 > $*.f90
.F90.i90:
	$(FC) $(CPPDEFS) $(FPPFLAGS) -P $*.F90
.F90.o:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c $*.F90
.F90.s:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -S $*.F90
.F90.x:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -o $*.x $*.F90 *.o $(LDFLAGS)
