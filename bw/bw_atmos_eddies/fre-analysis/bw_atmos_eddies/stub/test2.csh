#!/bin/csh -f
#------------------------------------
#PBS -N atmos_eddies
#PBS -l size=1
#PBS -l walltime=12:00:00
#PBS -r y
#PBS -j oe
#PBS -o
#PBS -q batch
#------------------------------------
# Source data: pp/atmos/ts/daily/Xxyr

# variables set by frepp
 set in_data_dir = /archive/Bonnie.Samuels/ulm_plus_LowSens_Drag_mom6_2015.03.13/CM4_c96L48_am4c3_2000_KPP_min6m/gfdl.ncrc2-intel-prod-openmp/pp/atmos/ts/daily/5yr
 set descriptor = CM4_test
 set out_dir = /nbhome/$USER/awg/fre-analysis-test/CM4_test
 set yr1 = 0001
 set yr2 = 0010
 set databegyr = 0001
 set dataendyr = 0010
 set datachunk = 5
 set staticfile = /archive/Bonnie.Samuels/ulm_plus_LowSens_Drag_mom6_2015.03.13/CM4_c96L48_am4c3_2000_KPP_min6m/gfdl.ncrc2-intel-prod-openmp/pp/atmos/atmos.static.nc
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
   module load ifort
   module load git
else
   echo "ERROR: invalid platform"
   exit 1
endif

# check again?
if (! $?FRE_ANALYSIS_HOME) then
   echo "ERROR: environment variable FRE_ANALYSIS_HOME not set."
   exit 1
endif

# clone the source code from the repository if it does not exist

set GIT_REPOSITORY = "file:///home/bw/git-repository"
set FRE_CODE_TAG = testing
set PACKAGE_NAME = bw_atmos_eddies
set FRE_CODE_BASE = $TMPDIR/fre-analysis

if (! -e $FRE_CODE_BASE/$PACKAGE_NAME) then
   if (! -e $FRE_CODE_BASE) mkdir $FRE_CODE_BASE
   cd $FRE_CODE_BASE
   git clone -b $FRE_CODE_TAG --recursive $GIT_REPOSITORY/fre-analysis/$PACKAGE_NAME.git
endif

##################
# run the script
##################

set options = "-i $in_data_dir -d $descriptor -o $out_dir -y $yr1,$yr2 -c $databegyr,$dataendyr,$datachunk -s $staticfile"
set filter = 2   # 0 = no filter, 1 = filtered, 2 = both

$FRE_CODE_BASE/$PACKAGE_NAME/run_eddies.pl $options -f $filter -p tt10,hh10,uu10,uu200,vv200,uv200,hh500,qt850,uv850

