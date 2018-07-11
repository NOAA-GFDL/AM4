
module load intel
module load hdf5
module load netcdf

setenv KMP_STACKSIZE 512m
setenv NC_BLKSZ 1M
setenv F_UFMTENDIAN big

set path = ($path ../mkmf ) 
