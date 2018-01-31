# AM4 Instructions

To checkout the source code for the AM4 model, use the following git
command <br>

```
git clone -r https://github.com/NOAA-GFDL/AM4.git
```

# Source Code

The source code is located in the src directory. This repository
contains the code for the following folders:

* atmos_cubed_sphere
* atmos_drivers
* atmos_param
* atmos_shared
* ice_param
* ice_sis
* land_lad2

The following folders are available on github and are linked as git
submodules

* MOM6
* coupler
* shared

# Analysis Scripts 
Some of the climate analysis scripts run at NOAA GFDL and used in the
AM4 documentation papers are located in the analysis directory.
Within each analysis suite, there is a jupyter notebook, both readable
and runnable from your local jupyter environment, provided all
dependencies are installed.

E.g.

* [Radiation processor](analysis/cjs1/radiation_atmos_av_mon/radiation_atmos_av_mon.ipynb)
* [Long-term DJF seasonal mean](analysis/bw/bw_atmos_cru_ts_a1r/bw_atmos_monthly_cru_ts.1980-2014.ipynb) 
* [Zonal_mean_zonal_wind_stress](analysis/bw/bw_atmos_zm_atl_pac_a1r/bw_atmos_atl_pac.1980-2014.ipynb)

Please refer https://www.gfdl.noaa.gov/am4.0-model for details about
where to find model and OBS data used in the papers.

For all analysis figures and pertaining data, please use the AM4
documentation papers as the original reference.  Please direct your
questions and feedback to gfdl.climate.model.info@noaa.gov


# Disclaimer

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
