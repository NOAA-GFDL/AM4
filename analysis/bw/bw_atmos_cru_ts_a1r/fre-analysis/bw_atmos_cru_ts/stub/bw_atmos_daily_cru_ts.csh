#!/bin/csh -f
#------------------------------------
#PBS -N bw_atmos_cru_daily_timeseries
#PBS -l size=1
#PBS -l walltime=12:00:00
#PBS -r y
#PBS -j oe
#PBS -o
#PBS -q batch
#----------------------------------
# Source data: pp/atmos/ts/daily/Xyr

# variables set by frepp
 set in_data_dir
 set out_dir
 set descriptor
 set yr1 
 set yr2 
 set databegyr
 set dataendyr
 set datachunk
 set staticfile
 set fremodule
 set freanalysismodule

# make sure valid platform and required modules are loaded
if (`gfdl_platform` == "hpcs-csc") then
   source $MODULESHOME/init/csh
   module purge
   module use -a /home/fms/local/modulefiles
   module load $fremodule
   module load $freanalysismodule
   module load ncarg/6.2.1
   module load ifort
   module load git
else
   echo "ERROR: invalid platform"
   exit 1
endif

# check again?
if (! $?FRE_ANALYSIS_GIT_URL) then
   echo "ERROR: environment variable FRE_ANALYSIS_GIT_URL not set."
   exit 1
endif

# clone the source code from the repository if it does not exist

set GIT_REPOSITORY = $FRE_ANALYSIS_GIT_URL/bw
set FRE_CODE_TAG = master
set PACKAGE_NAME = bw_atmos_cru_ts
set FRE_CODE_BASE = $TMPDIR/fre-analysis

if (! -e $FRE_CODE_BASE/$PACKAGE_NAME) then
   if (! -e $FRE_CODE_BASE) mkdir $FRE_CODE_BASE
   cd $FRE_CODE_BASE
   git clone -b $FRE_CODE_TAG --recursive $GIT_REPOSITORY/$PACKAGE_NAME.git
endif

##################
# run the script
##################

set options = "-i $in_data_dir -d $descriptor -o $out_dir -y $yr1,$yr2 -c $databegyr,$dataendyr,$datachunk -s $staticfile --statsfile"

# daily temp (tasmean, tasmin, tasmax, tasdiff); and wet days
$FRE_CODE_BASE/$PACKAGE_NAME/runscript.pl $options -t mean,stdev -v tas
#$FRE_CODE_BASE/$PACKAGE_NAME/runscript.pl $options -v wet

