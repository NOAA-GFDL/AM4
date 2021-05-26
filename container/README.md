# AM4 containers
This is a very basic *0.1* version of a README for these containers.  Please feel free to
add to it or open a GitHub issue if there is something missing.
 
## Building with Docker
The Dockerfiles are set up to build an AM4 run using Docker.  There are two Dockerfiles, 
one to build with intel oneAPI 2021.2 and one to build with GCC 10.2.0

## Building with Singularity
The Singularity definition files are included to build using intel oneAPI 2021.2 compilers.
You can build using the singularity_build.sh script
```bash
./singularity_build.sh
```

## Running using singularity
The containers are all using mpich-compatible MPI, so if you run using singularity bind 
or hybrid methods, make sure you are using some flavor of mpich and not openmpi.

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

