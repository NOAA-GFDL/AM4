# GFDL AM4 Model

[![DOI](https://zenodo.org/badge/102487636.svg)](https://zenodo.org/badge/latestdoi/102487636)

This repository includes the public release of the GFDL AM4 model
code.  The AM4 model is described in the
[two](https://doi.org/10.1002/2017MS001208)
[articles](https://doi.org/10.1002/2017MS001209) published in the
[Journal of Advances in Modeling Earth Systems
(JAMES)](https://agupubs.onlinelibrary.wiley.com/journal/19422466).
More information on the model and access to the output is available on
the [AM4 data and code
site](http://data1.gfdl.noaa.gov/nomads/forms/am4.0/) at the
[Geophysical Fluid Dynamics Laboratory
(GFDL)](https://www.gfdl.noaa.gov).

The layout of this package includes the following directories:

* src - The source code for the AM4 model
* exec - The build directory with Makefiles for building the model executable
* run - Sample run script
* analysis - Sample analysis scripts

## Cloning Instructions

This repository uses [git
submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to
point to other repositories.  Thus, care should be taken when cloning,
and updating the source to ensure all source.  To obtain all source,
use the following git command

```
git clone --recursive https://github.com/NOAA-GFDL/AM4.git
```

The `--recursive` option to `git clone` instructs git to recursively
clone all submodules.  In the event the repository was not cloned
using the `--recursive` option, the following step must be taken to
obtain all sources:

```
# From within the AM4 parent directory
git submodule update --init --recursive
```

## Source Code

All model source is contained in the [src](src) directory.  GFDL
tracks code using the git version control system.  This package
includes a single version of the following GFDL model components.  The
git hash listed corresponds to the commit hash in the internal GFDL
git repository.

Component | Commit Hash
--------- | -----------
atmos_drivers | 5ee95d6abf0879594551dd7e6635dff4004c4010
atmos_param | 2e94acfd8621e85216bf822c395a8c3f15a511a5
atmos_shared | a557d4d7bab033ef1ad1d400a62fe07a97ccb477
ice_param | 1553c8bc4f9a66791c89367b6f327147523155ed
ice_sis | ccc7328dcd79706dd5c17c8bab660222886fc80b
land_lad2 | a220288ecb289bf9d793d051fc5076072874ce07

The following components are available in the
[NOAA-GFDL](https://github.com/NOAA-GFDL) github organization:

* [MOM6](https://github.com/NOAA-GFDL/MOM6)
* [coupler](https://github.com/NOAA-GFDL/coupler)
* [FMS](https://github.com/NOAA-GFDL/FMS) (as [shared](src/shared))
* [GFDL_atmos_cubed_sphere (tag AM4.0)](https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere) (as [atmos_cubed_sphere](src/atmos_cubed_sphere))

## Building AM4

The [exec](exec) directory contains Makefiles that can be used to
build the AM4 executable.  These Makefiles were generated using the
[Make Makefile (mkmf)](https://github.com/NOAA-GFDL/mkmf) program.
Included in the exec direcgtory is a sample make template file for the
Intel compilers ([intel.mk](exec/templates/intel.mk)).  This make
template can be used on any system with a relatively recent version of
the Intel compilers, the netCDF 4 library and the MPICH2 MPI library.
Included in the [intel.mk](exec/templates/intel.mk) file are
additional settings that can be modified during the build.  


To run the default build (-O3 -msse2), go to the exec directory and
enter the command
```
make
```
If you would like to change some of the compiler options, there are several different
options to add to the make command.  For example
```
make ISA=-xhost BLD_TYPE=REPRO
```
will replace -msse with -xhost and -O3 with -O2.  The three options for
`BLD_TYPE` are  
`PROD` (-O3)  
`REPRO` (-O2)    
`DEBUG` (-O0 and other traps)  
All of the make line options can be
found in the [intel.mk](exec/templates/intel.mk) file.

## Obtaining the input data

The input data required for running the AM4 model can be found on
[GFDL's data
portal](http://data1.gfdl.noaa.gov/nomads/forms/am4.0/) .

The file `AM4.tar.gz` contains a configured run directory to run a
sample experiment of the AM4 model.  Included in the tar file is a
README.AM4_run with more instructions on how to configure the AM4 run
directory.

On Linux systems, the `wget` command is usually sufficient to download the data
file:

```
wget ftp://nomads.gfdl.noaa.gov/users/Ming.Zhao/AM4Documentation/GFDL-AM4.0/inputData/AM4_run.tar.gz
```

To ensure the file downloaded is complete and not corrupted, download one of the two files:

```
wget ftp://nomads.gfdl.noaa.gov/users/Ming.Zhao/AM4Documentation/GFDL-AM4.0/inputData/AM4_run.tar.gz.sha256
wget ftp://nomads.gfdl.noaa.gov/users/Ming.Zhao/AM4Documentation/GFDL-AM4.0/inputData/AM4_run.tar.gz.sig
```

and run the following command that corresponds to the signature file downloaded:

```
sha256sum -c AM4_run.tar.gz.sha256
```

```
gpg --verify AM4_run.tar.gz.sig
```

## Running AM4

Included in the run directory is a sample run script for reference.
To run the AM4 sample experiment, first download the data file
mentioned in [Obtaining the Input data](#obtaining-the-input-data)
section.  Modify the variables in the configuration section in the
sample run script, and then run the script.

The sample data and run script are configured to run on 216
processors.  To run on a different number of processors, or modify the
experiment, refer to the `README.AM4_run` file included in the AM4
data tarball.

Note: The `input.nml` file (found in the AM4 data tarball) contains
Fortran namelists and namelist variables that modify, at run time, the
model.  To learn more about the settings in the `input.nml` file,
please refer to source code where the namelist/variable are defined.

## Analysis Scripts

Some of the climate analysis scripts run at NOAA GFDL and used in the
AM4 documentation papers are located in the analysis directory.
Within each analysis suite, is a [jupyter
notebook](https://jupyter-notebook.readthedocs.io/en/stable/), both
readable and runnable from your local jupyter environment, provided
all dependencies are installed.

E.g.

* [Radiation processor](analysis/cjs1/radiation_atmos_av_mon/radiation_atmos_av_mon.ipynb)
* [Long-term DJF seasonal mean](analysis/bw/bw_atmos_cru_ts_a1r/bw_atmos_monthly_cru_ts.1980-2014.ipynb)
* [Zonal_mean_zonal_wind_stress](analysis/bw/bw_atmos_zm_atl_pac_a1r/bw_atmos_atl_pac.1980-2014.ipynb)
* [PCMDI Metrics Portrait Plot](analysis/pcmdimetrics/portraitPlot-AM4.AMIP.ipynb)

## Model output and Other References

Please refer to the [AM4 data and code
site](http://data1.gfdl.noaa.gov/nomads/forms/am4.0/) for details
about where to find model and OBS data used in the papers.

For all analysis figures and pertaining data, please use the AM4
documentation papers as the original reference.

Please direct your questions and feedback to
gfdl.climate.model.info@noaa.gov

## Disclaimer

The United States Department of Commerce (DOC) GitHub project code is
provided on an 'as is' basis and the user assumes responsibility for
its use.  DOC has relinquished control of the information and no
longer has responsibility to protect the integrity, confidentiality,
or availability of the information.  Any claims against the Department
of Commerce stemming from the use of its GitHub project will be
governed by all applicable Federal law.  Any reference to specific
commercial products, processes, or services by service mark,
trademark, manufacturer, or otherwise, does not constitute or imply
their endorsement, recommendation or favoring by the Department of
Commerce.  The Department of Commerce seal and logo, or the seal and
logo of a DOC bureau, shall not be used in any manner to imply
endorsement of any commercial product or activity by DOC or the United
States Government.

This project code is made available through GitHub but is managed by
NOAA-GFDL at https://gitlab.gfdl.noaa.gov.
