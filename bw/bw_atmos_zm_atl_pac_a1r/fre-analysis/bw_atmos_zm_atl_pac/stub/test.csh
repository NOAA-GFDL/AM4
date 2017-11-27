#!/bin/csh -f
#------------------------------------
#PBS -N bw_atmos_atl_pac
#PBS -l size=1
#PBS -l walltime=04:00:00
#PBS -r y
#PBS -j oe
#PBS -o
#PBS -q batch
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Script: bw_atmos_zm_atl_pac.csh
# Author: Bruce Wyman
# Source data: pp/atmos/av/monthly_Xyr
# Output: Creates figures in $out_dir/atmos_${yr1}_${yr2}
#
# Sample frepp usage (http://www.gfdl.noaa.gov/fms/fre/#analysis):
# <component type="atmos">
#    <timeAverage ... >
#       <analysis script="script_name [options]"/>
#    </timeAverage>
# </component>

# variables set by frepp
 set in_data_dir = "/archive/h1g/awg/ulm_201505/c96L32_am4g5r1/gfdl.ncrc2-intel-prod-openmp/pp/atmos/av/monthly_30yr"
 set in_data_file = "atmos.1981-2010.{01,02,03,04,05,06,07,08,09,10,11,12}.nc"
 set descriptor = "c96L32_am4g5r1"
 set out_dir = "/nbhome/$USER/awg/ulm_201505/c96L32_am4g5r1_test"
 set yr1 = 1981
 set yr2 = 2010
 set fremodule = "fre"
 set freanalysismodule = "fre-analysis"

# make sure valid platform and required modules are loaded
if (`gfdl_platform` == "hpcs-csc") then
   source $MODULESHOME/init/csh
   module purge
   module use -a /home/fms/local/modulefiles
   module load $fremodule
   module load $freanalysismodule
   module load ncarg/6.2.1
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
set PACKAGE_NAME = bw_atmos_zm_atl_pac
set FRE_CODE_BASE = $TMPDIR/fre-analysis

if (! -e $FRE_CODE_BASE/$PACKAGE_NAME) then
   if (! -e $FRE_CODE_BASE) mkdir $FRE_CODE_BASE
   cd $FRE_CODE_BASE
   git clone -b $FRE_CODE_TAG --recursive $GIT_REPOSITORY/$PACKAGE_NAME.git
endif

##################
# run the script
##################

if ($#argv == 0) then
  set fields = (tauu pr curl) # plot all fields by default
else
  set fields = ($argv)
endif

set options = "-i $in_data_dir -d $descriptor -y $yr1,$yr2 -o $out_dir"

foreach field ($fields)
  if (`echo $field | perl -e '$f=<>;chomp$f;if(grep{$f eq $_}qw/tauu pr curl/){print"1"}else{print"0"}'`) then
    $FRE_CODE_BASE/$PACKAGE_NAME/run_zm_atl_pac.pl -v $field $options $in_data_file
  else
    echo "WARNING: invalid field: $field"
  endif
end

