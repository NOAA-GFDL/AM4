#!/bin/sh

singularity build -f intel2021.2_netcdfc4.7.4_ubuntu.sif Singularity.intel_netcdf
singularity build -f am4_2021.03_ubuntu_intel.sif Singularity.intel_am4
