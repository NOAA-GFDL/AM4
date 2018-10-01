#!/bin/sh

# Sample run script to run the am4p0 experiment

# ***********************************************************************
# Modify the settings in this section to match your system's environment
# and the directory locations to the executable, input data, initial
# conditions data and work directory.

# Name of the mpiexec program to use
mpiexec_prog=aprun
# Option used to specify number of MPI process to run (usually -n or -np)
mpiexec_nopt=-n
# Option used to specify number of OpenMP threads to run
mpiexec_topt=-d

# Where to perform the run
# If using AM4.tar, this should be AM4_run
workDir=/path/to/run/dir

## Location of data files
# The input files can be skipped if the input is already in $workdir/INPUT
#inputDataRoot=/path/to/input/data
#inputDataTar=${inputDataRoot}/inputData.tar.gz
#initCondTar=${inputDataRoot}/.tar.gz

# Location of executable (run with $mpiexec_prog)
executable=/path/to/executable/fms_cm4p12_warsaw.x


## Run parameters
#total_npes is the number of cores to run on, omp_threads is the number of
# openMP threads
total_npes=432
omp_threads=2

# End of configuration section
# ***********************************************************************

# Enviornment settings for run
export KMP_STACKSIZE=512m
export NC_BLKSZ=1M
export F_UFMTENDIAN=big

# Remember CWD
initialDir=$(pwd)

# check of required programs
if ! hash tar 2> /dev/null
then
  echo "ERROR: Unable to find \`tar\` in PATH." 1>&2
  echo "ERROR: Halting script." 1>&2
fi
if ! hash ${mpiexec_prog} 2> /dev/null
then
  echo "ERROR: Unable to find \`${mpiexec_prog}\` in PATH." 1>&2
  echo "ERROR: Halting script." 1>&2
fi


# Verify work directory exists, if not create it
if [ ! -e ${workDir} ]
then
  mkdir -p ${workDir}
  if [ $? -ne 0 ]
  then
    echo "ERROR: Unable to create work directory \"${workDir}\"." 1>&2
    echo "ERROR: Halting script." 1>&2
    exit 1
  fi
elif [ ! -d ${workDir} ]
then
  echo "ERROR: Work directory \"${workDir}\" is not a directory." 1>&2
  echo "ERROR: Halting script." 1>&2
  exit 1
fi

# Check if work directory is empty, warn if not
if [ $(ls -1qA ${workDir} | wc -l) -gt 0 ]
then
  echo "NOTE: Work directory \"${workDir}\" is not empty." 1>&2
  echo "NOTE: Data in \"${workDir}\" will be overwritten." 1>&2
fi

# Enter working directory, and setup the directory
cd ${workDir}
if [ $? -ne 0 ]
then
  echo "ERROR: Unable \`cd\` into work directory \"${workDir}\"." 1>&2
  echo "ERROR: Halting script." 1>&2
  exit 1
fi

# Create RESTART directory, if it doesn't eixt.
if [ ! -e RESTART ]
then
  mkdir RESTART
  if [ $? -ne 0 ]
  then
    echo "ERROR: Unable to create directory \"${workDir}/RESTART\"." 1>&2
    echo "ERROR: Halting script." 1>&2
    exit 1
  fi
elif [ ! -d RESTART ]
then
  echo "ERROR: Directory \"${workDir}/RESTART\" is not a directory." 1>&2
  echo "ERROR: Halting script." 1>&2
  exit 1
elif [ $(ls -1qA ${workDir}/RESTART | wc -l) -gt 0 ]
then
  echo "WARNING: Directory \"${workDir}/RESTART\" is not empty." 1>&2
  echo "WARNING: Contents will be overwritten." 1>&2
fi

## Use this section if you are untar'ing the input data ##
## Not required if using AM4.tar out of the box          ##
## Extract the input data
#tar xf ${inputDataTar}
#if [ $? -ne 0 ]
#then
#  echo "ERROR: Unable to extract data from \"${inputDataTar}\"." 1>&2
#  echo "ERROR: Halting script." 1>&2
#fi
#
#tar xf ${initCondTar}
#if [ $? -ne 0 ]
#then
#  echo "ERROR: Unable to extract data from \"${initCondTar}\"." 1>&2
#  echo "ERROR: Halting script." 1>&2
#fi

# Run the model
${mpiexec_prog} ${mpiexec_nopt} ${total_npes} ${mpiexec_topt} ${omp_threads} ${executable} 2>&1 | tee ${workDir}/fms.out
if [ $? -ne 0 ]
then
  echo "ERROR: Run failed." 1>&2
  echo "ERROR: Output from run in \"${workDir}/fms.out\"." 1>&2
  exit 1
fi

# Return to the initial directory
cd ${initialDir}
