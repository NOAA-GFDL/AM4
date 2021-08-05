#!/bin/sh

singularity build -f intel_netcdf_ubuntu.sif Singularity.intel2021.2_netcdfc4.7.4
singularity build -f am4_2021.03_ubuntu_intel.sif Singularity.intel_am4
