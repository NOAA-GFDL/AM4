#!/bin/sh

singularity build -f intel_netcdf_ubuntu.sif Singularity.intel_netcdf
singularity build -f am4_2021.02_ubuntu_intel.sif Singularity.intel_am4
