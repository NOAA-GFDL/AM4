# Template for the Intel Compilers on Linux systems
#
# Typical use with mkmf
# mkmf -t linux-intel.mk -c"-Duse_libMPI -Duse_netCDF" path_names /usr/local/include

############
# Command Macros
FC = mpifort
CC = mpicc
CXX = mpicpc
LD = mpifort

#######################
# Build target macros
#
# Macros that modify compiler flags used in the build.  Target
# macrose are usually set on the call to make:
#
#    make BLD_TYPE=PROD NETCDF=3
#
# Most target macros are activated when their value is non-blank.
# Some have a single value that is checked.  Others will use the
# value of the macro in the compile command.

# BLD_TYPE
# Determines the type of build.  Values are:
# PROD - Use the production settings (default)
# REPRO - Extra options to guarentee run to run reproducibility.
# DEBUG - Compile with debug options (-O0 -g)
# TEST - Use additional compiler options defined in FFLAGS_TEST
#        and CFLAGS_TEST
ifndef BLD_TYPE
BLD_TYPE = PROD
endif

# NETCDF_FLAGS
# NETCDF_LIBS
# If defined, use the NETCDF compile and link options defined in these
# variables.  If these options are not defined, the makefile will
# attempt to get the correct options from the `nf-config` command.

# MPI_FLAGS
# MPI_LIBS
# If defined, use the MPI compile and link options defined in these
# variables.  If these options are not defined, the makefile will
# attempt to get the correct options from the `pkg-config` for mpich2
# MPI library.

# VERBOSE
# If non-blank, add additional verbosity compiler options

# OPENMP
# If non-blank, compile with openmp enabled

# NO_OVERRIDE_LIMITS
# If non-blank, do not use the -qoverride-limits compiler option.
# Default behavior is to compile with -qoverride-limits.

# NETCDF
# If value is '3' (default) and CPPDEFS contains '-Duse_netCDF', then
# the additional cpp macro '-Duse_LARGEFILE' is added to the CPPDEFS
# macro.
ifndef NETCDF
NETCDF = 3
endif

# INCLUDES
#A list of -I Include directories to be added to the the compile
#command.

# ISA
# The Intel Instruction Set Archetecture (ISA) compile options to use.
# If blank, than use the default ISA settings for the host.
ifndef ISA
ISA = -msse2
endif

# COVERAGE
# If non-blank Add the code coverage compile options.

# Need to use at least GNU Make version 3.81
need := 3.81
ok := $(filter $(need),$(firstword $(sort $(MAKE_VERSION) $(need))))
ifneq ($(need),$(ok))
$(error Need at least make version $(need).  Load module gmake/3.81)
endif

MAKEFLAGS += --jobs=$(shell grep '^processor' /proc/cpuinfo | wc -l)

# Macro for Fortran preprocessor
FPPFLAGS = $(INCLUDES)
# Fortran Compiler flags for the NetCDF library
ifndef NETCDF_FLAGS
FPPFLAGS += $(shell nf-config --fflags)
else
FPPFLAGS += $(NETCDF_FLAGS)
endif
# Fortran Compiler flags for the MPICH MPI library
ifndef MPI_FLAGS
FPPFLAGS += $(shell pkg-config --cflags-only-I mpich2-c)
else
FPPFLAGS += $(MPI_FLAGS)
endif

# Base set of Fortran compiler flags
FFLAGS := -fcray-pointer -fdefault-real-8 -fdefault-double-8 -Waliasing -ffree-line-length-none -fno-range-check

# Flags based on perforance target (production (OPT), reproduction (REPRO), or debug (DEBUG)
FFLAGS_OPT = -O2 -fno-expensive-optimizations
FFLAGS_REPRO =
FFLAGS_DEBUG = -O0 -g -W -fbounds-check -ffpe-trap=invalid,zero,overflow

# Flags to add additional build options
FFLAGS_OPENMP = -fopenmp
FFLAGS_OVERRIDE_LIMITS = 
FFLAGS_VERBOSE = -Wall -Wextra
FFLAGS_COVERAGE =

# Macro for C preprocessor
CPPFLAGS = -D__IFC $(INCLUDES)
# C Compiler flags for the NetCDF library
ifndef NETCDF_FLAGS
CPPFLAGS += $(shell nc-config --cflags)
else
CPPFLAGS += $(NETCDF_FLAGS)
endif
# C Compiler flags for the MPICH MPI library
ifndef MPI_FLAGS
CPPFLAGS += $(shell pkg-config --cflags-only-I mpich2-c)
else
CPPFLAGS += $(MPI_FLAGS)
endif

# Base set of C compiler flags
CFLAGS := 

# Flags based on perforance target (production (OPT), reproduction (REPRO), or debug (DEBUG)
CFLAGS_PROD = -O2
CFLAGS_REPRO = -O2
CFLAGS_DEBUG = -O0 -g 

# Flags to add additional build options
CFLAGS_OPENMP = -fopenmp
CFLAGS_VERBOSE = -Wall -Wextra
CFLAGS_COVERAGE =

# Optional Testing compile flags.  If FFLAGS_TEST or CFLAGS_TEST are not defined, then the PROD
# compile settings will be used
ifndef FFLAGS_TEST
FFLAGS_TEST = $(FFLAGS_PROD)
endif
ifndef CFLAGS_TEST
CFLAGS_TEST = $(CFLAGS_OPT)
endif

# Linking flags
LDFLAGS :=
LDFLAGS_OPENMP := -fopenmp
LDFLAGS_VERBOSE :=
LDFLAGS_COVERAGE :=

# Start with a blank LIBS
LIBS =
# NetCDF library flags
ifndef NETCDF_LIBS
LIBS += $(shell nf-config --flibs)
else
LIBS += $(NETCDF_LIBS)
endif
# MPICH MPI library flags
ifndef MPI_LIBS
LIBS += $(shell pkg-config --libs mpich2-f90)
else
LIBS += $(MPI_LIBS)
endif
# HDF library flags
ifndef HDF_LIBS
LIBS += -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran
else
LIBS += $(HDF_LIBS)
endif
# MKL library flags
ifndef MKL_LIBS
#LIBS += -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
else
LIBS += $(MKL_LIBS)
endif

# Get compile flags based on target macros.
ifeq ($(BLD_TYPE),REPRO)
CFLAGS += $(CFLAGS_REPRO)
FFLAGS += $(FFLAGS_REPRO)
else ifeq ($(BLD_TYPE),DEBUG)
CFLAGS += $(CFLAGS_DEBUG)
FFLAGS += $(FFLAGS_DEBUG)
else ifeq ($(BLD_TYPE),TEST)
CFLAGS += $(CFLAGS_TEST)
FFLAGS += $(FFLAGS_TEST)
else
CFLAGS += $(CFLAGS_PROD)
FFLAGS += $(FFLAGS_PROD)
endif

ifdef OPENMP
CFLAGS += $(CFLAGS_OPENMP)
FFLAGS += $(FFLAGS_OPENMP)
LDFLAGS += $(LDFLAGS_OPENMP)
endif

ifdef ISA
CFLAGS += $(ISA)
FFLAGS += $(ISA)
endif

ifndef NO_OVERRIDE_LIMITS
FFLAGS += $(FFLAGS_OVERRIDE_LIMITS)
endif

ifdef VERBOSE
CFLAGS += $(CFLAGS_VERBOSE)
FFLAGS += $(FFLAGS_VERBOSE)
LDFLAGS += $(LDFLAGS_VERBOSE)
endif

ifeq ($(NETCDF),3)
  # add the use_LARGEFILE cppdef
  ifneq ($(findstring -Duse_netCDF,$(CPPDEFS)),)
    CPPDEFS += -Duse_LARGEFILE
  endif
endif

ifdef COVERAGE
ifdef BUILDROOT
PROF_DIR=-prof-dir=$(BUILDROOT)
endif
CFLAGS += $(CFLAGS_COVERAGE) $(PROF_DIR)
FFLAGS += $(FFLAGS_COVERAGE) $(PROF_DIR)
LDFLAGS += $(LDFLAGS_COVERAGE) $(PROF_DIR)
endif

LDFLAGS += $(LIBS)

#---------------------------------------------------------------------------
# you should never need to change any lines below.

# see the MIPSPro F90 manual for more details on some of the file extensions
# discussed here.
# this makefile template recognizes fortran sourcefiles with extensions
# .f, .f90, .F, .F90. Given a sourcefile <file>.<ext>, where <ext> is one of
# the above, this provides a number of default actions:

# make <file>.opt	create an optimization report
# make <file>.o		create an object file
# make <file>.s		create an assembly listing
# make <file>.x		create an executable file, assuming standalone
#			source
# make <file>.i		create a preprocessed file (for .F)
# make <file>.i90	create a preprocessed file (for .F90)

# The macro TMPFILES is provided to slate files like the above for removal.

RM = rm -f
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
