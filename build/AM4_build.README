## Requirements
        netcdf          https://www.unidata.ucar.edu/software/netcdf/
        hdf5            https://www.unidata.ucar.edu/software/netcdf/
        mkmf            https://github.com/NOAA-GFDL/mkmf
        list_paths      https://github.com/NOAA-GFDL/mkmf
        c-shell
        compiler

## Quick Compiling Instructions 
1. tcsh
2. cd build
3. ./compile.csh |& tee log.compile

## Important Notes 
The compile script is a csh script.  If you are already in csh or tcsh, you can 
skip 1.

mkmf and list_paths must be in your path.  The line in env.cshrc attempts 
to add these to your path assuming you cloned the mkmf repository one directory 
up 
`set path = ($path ../mkmf )`

If you are not using modules, you can delete the module load lines in end.cshrc

The top of compile.csh should be edited to match where you intend to build, 
where your source (src_dir) is located, and the compile template (mkmf_template).
If you plan on using gcc/gfortran, switch intel.mk to gnu.mk:{
        # ---------------- Set build, src and stage directories

        set src_dir = ../src
        set bld_dir = ${PWD}
        set ptmp_dir = /tmp

        # ---------------- Make template

        set mkmf_template = intel.mk

        # ---------------- set environment

        if ( $echoOn ) unset echo
        source $bld_dir/env.cshrc
        if ( $echoOn ) set echo
}

There are 4 compile options:{
1. Default - uses -O3
2. REPRO=on - Uses -O2
3. DEBUG=on - uses -O0
4. OPENMP=on - uses -qopenmp
You can invoke these options on the make line in the compile script. If you 
want to compile with -O2 and openMP, the make line should look like this:
        make  REPRO=on OPENMP=on NETCDF=3 fms_cm4p12_warsaw.x
NOTE: if compiling with gcc/gfortran, do not compile with openMP.  Remove it 
from the make line in the compile script
}


