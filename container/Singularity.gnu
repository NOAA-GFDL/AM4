Bootstrap: docker
From: thomasrobinson/centos7-netcdff:4.5.3-c4.7.4-gcc-mpich-slurm
Stage: build
## Singularity def file used to create AM4

%post
## Set up spack
 . /opt/spack/share/spack/setup-env.sh
## Make the AM4 directory
 mkdir -p /opt/AM4
 cd /opt
## Build the AM4 from github
 git clone --recursive https://github.com/NOAA-GFDL/AM4.git -b 2021.02 
     cd AM4/exec 
     make -j 20 gcc=on HDF_INCLUDE=-I/opt/hdf5/include SH=sh CLUBB=off 
     cp am4_xanadu_2021.02.x /opt/AM4 
     make clean_all
 chmod 777 /opt/AM4/am4_xanadu_2021.02.x

## Add the AM4 executable to the path
%environment
ENV PATH=/opt/AM4/:${PATH}

## Run AM4
%runscript
 ulimit -s unlimited
 /opt/AM4/am4_xanadu_2021.02.x

