FROM thomasrobinson/centos7-netcdff:4.5.3-c4.7.4-gcc-mpich-slurm
## Dockerfile used to create AM4

## Set up spack
RUN . /opt/spack/share/spack/setup-env.sh
## Make the AM4 directory
RUN mkdir -p /opt/AM4
## Build the AM4 from github
RUN git clone --recursive https://github.com/NOAA-GFDL/AM4.git -b 2021.02 \
    && cd AM4/exec \ 
    && make gcc=on HDF_INCLUDE=-I/opt/hdf5/include SH=sh CLUBB=off \
    && cp am4_xanadu_2021.02.x /opt/AM4 \
    && make clean_all
## Add the AM4 executable to the path
ENV PATH=/opt/AM4/:${PATH}
## Add permissions to the AM4
RUN chmod 777 /opt/AM4/am4_xanadu_2021.02.x

