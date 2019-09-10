# This is the main configure file for the GFDL AM4
# package.
# Ed Hartnett 7/29/2019

AC_PREREQ([2.59])

# Initialize with name, version, and support email address.
AC_INIT([AM4], [2.0-development], [])

# Find out about the host we're building on.
AC_CANONICAL_HOST

# Find out about the target we're building for.
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE([foreign dist-zip subdir-objects])

# Keep libtool macros in an m4 directory.
AC_CONFIG_MACRO_DIR([m4])

# Set up libtool.
LT_PREREQ([2.4])
LT_INIT()

# Find the C compiler.
AC_PROG_CC
AM_PROG_CC_C_O
AC_C_CONST
AC_PROG_CPP

# Find the Fortran compiler.
AC_PROG_FC
AC_PROG_F77

# Find the install program.
AC_PROG_INSTALL

# Check to see if any macros must be set to enable large (>2GB) files.
AC_SYS_LARGEFILE

# Require MPI.
AC_CHECK_FUNC([MPI_Init], [], [AC_MSG_ERROR([MPI C library required])])

# Require FMS. But this does not work because fms_platform.h is actually a fortran file.
# AC_CHECK_HEADERS([fms_platform.h], [], [AC_MSG_ERROR([Can't find fms_platform.h.])])

dnl # Check for netCDF C library.
dnl AC_SEARCH_LIBS([nc_create], [netcdf], [],
dnl                             [AC_MSG_ERROR([Can't find or link to the netcdf C library, set CPPFLAGS/LDFLAGS.])])

# Check for netCDF Fortran library.
AC_LANG_PUSH(Fortran)
AC_SEARCH_LIBS([nf_create], [netcdff], [],
                            [AC_MSG_ERROR([Can't find or link to the netcdf Fortran library, set CPPFLAGS/LDFLAGS.])])
AC_LANG_POP(Fortran)

# Require netCDF.
#AC_CHECK_FUNC([nf_open], [], [AC_MSG_ERROR([NetCDF Fortran library required])])

# These defines are required for the build.
#AC_DEFINE([use_netCDF], [1])
#AC_DEFINE([use_libMPI], [1])
#AC_DEFINE([SPMD], [1])

# These files will be created when the configure script is run.
AC_CONFIG_FILES([Makefile
        analysis/Makefile
        analysis/cjs1/Makefile
        analysis/cjs1/radiation_atmos_av_mon/Makefile
        analysis/cjs1/radiation_atmos_av_mon/plot_file_scripts/Makefile
        analysis/pcmdimetrics/Makefile
        analysis/bw/Makefile
        analysis/bw/bw_atmos_cru_ts_a1r/Makefile
        analysis/bw/bw_atmos_cru_ts_a1r/fre-analysis/Makefile
        analysis/bw/bw_atmos_cru_ts_a1r/fre-analysis/bw_atmos_cru_ts/Makefile
        analysis/bw/bw_atmos_cru_ts_a1r/fre-analysis/bw_atmos_cru_ts/stub/Makefile
        analysis/bw/bw_atmos_cru_ts_a1r/fre-analysis/bw_atmos_cru_ts/ncl/Makefile
        analysis/bw/bw_atmos_cru_ts_a1r/fre-analysis/bw_atmos_cru_ts/ncl/colormaps/Makefile
        analysis/bw/bw_atmos_cru_ts_a1r/fre-analysis/bw_atmos_cru_ts/ncl/wrapit/Makefile
        analysis/bw/bw_atmos_cru_ts_a1r/fre-analysis/bw_atmos_cru_ts/shared/Makefile
        analysis/bw/bw_atmos_cru_ts_a1r/fre-analysis/bw_atmos_cru_ts/shared/bw/Makefile
        analysis/bw/bw_atmos_cru_ts_a1r/fre-analysis/bw_atmos_cru_ts/shared/bw/ncl/Makefile
        analysis/bw/bw_atmos_cru_ts_a1r/fre-analysis/bw_atmos_cru_ts/shared/bw/ncl/wrapit/Makefile
        analysis/bw/bw_atmos_zm_atl_pac_a1r/Makefile
        analysis/bw/bw_atmos_zm_atl_pac_a1r/fre-analysis/Makefile
        analysis/bw/bw_atmos_zm_atl_pac_a1r/fre-analysis/bw_atmos_zm_atl_pac/Makefile
        analysis/bw/bw_atmos_zm_atl_pac_a1r/fre-analysis/bw_atmos_zm_atl_pac/stub/Makefile
        analysis/bw/bw_atmos_zm_atl_pac_a1r/fre-analysis/bw_atmos_zm_atl_pac/ncl/Makefile
        src/Makefile
        src/ice_param/Makefile
        src/coupler/Makefile
        src/land_lad2/Makefile
        src/land_lad2/transitions/Makefile
        src/land_lad2/canopy_air/Makefile
        src/land_lad2/glacier/Makefile
        src/land_lad2/vegetation/Makefile
        src/land_lad2/soil/Makefile
        src/land_lad2/snow/Makefile
        src/land_lad2/topo_rough/Makefile
        src/land_lad2/river/Makefile
        src/land_lad2/lake/Makefile
        src/land_lad2/shared/Makefile
        src/land_lad2/shared/nf_utils/Makefile
        src/land_lad2/predefined_tiles/Makefile
        src/land_lad2/land_tracers/Makefile
        src/atmos_shared/Makefile
        src/atmos_shared/aerosol/Makefile
        src/atmos_shared/atmos_nudge/Makefile
        src/atmos_shared/vert_advection/Makefile
        src/atmos_shared/atmos_cmip_diag/Makefile
        src/atmos_shared/tracer_driver/Makefile
        src/atmos_shared/tracer_driver/stratchem/Makefile
        src/atmos_shared/tracer_driver/aer_ccn_act/Makefile
        src/atmos_shared/tracer_driver/tropchem/Makefile
        src/atmos_shared/tracer_driver/tropchem/AM3/Makefile
        src/atmos_shared/tracer_driver/tropchem/AM4/Makefile
        src/MOM6/Makefile
        src/atmos_param/Makefile
        src/atmos_param/dry_adj/Makefile
        src/atmos_param/grey_radiation/Makefile
        src/atmos_param/microphysics/Makefile
        src/atmos_param/mg_drag/Makefile
        src/atmos_param/cloud_generator/Makefile
        src/atmos_param/betts_miller/Makefile
        src/atmos_param/physics_driver/Makefile
        src/atmos_param/rh_clouds/Makefile
        src/atmos_param/cloud_rad/Makefile
        src/atmos_param/diffusivity/Makefile
        src/atmos_param/vert_diff/Makefile
        src/atmos_param/strat_cloud/Makefile
        src/atmos_param/lin_cloud_microphys/Makefile
        src/atmos_param/physics_radiation_exch/Makefile
        src/atmos_param/shallow_conv/Makefile
        src/atmos_param/donner_deep/Makefile
        src/atmos_param/aerosol_cloud/Makefile
        src/atmos_param/vert_diff_driver/Makefile
        src/atmos_param/cloud_obs/Makefile
        src/atmos_param/diag_cloud_rad/Makefile
        src/atmos_param/cloud_zonal/Makefile
        src/atmos_param/moist_conv/Makefile
        src/atmos_param/clouds/Makefile
        src/atmos_param/diag_cloud/Makefile
        src/atmos_param/lscale_cond/Makefile
        src/atmos_param/two_stream_gray_rad/Makefile
        src/atmos_param/ras/Makefile
        src/atmos_param/cg_drag/Makefile
        src/atmos_param/cu_mo_trans/Makefile
        src/atmos_param/convection_driver/Makefile
        src/atmos_param/stable_bl_turb/Makefile
        src/atmos_param/cosp/Makefile
        src/atmos_param/cosp/icarus-scops-3.7/Makefile
        src/atmos_param/cosp/MODIS_simulator/Makefile
        src/atmos_param/cosp/llnl/Makefile
        src/atmos_param/cosp/actsim/Makefile
        src/atmos_param/cosp/MISR_simulator/Makefile
        src/atmos_param/cosp/quickbeam/Makefile
        src/atmos_param/shallow_physics/Makefile
        src/atmos_param/monin_obukhov/Makefile
        src/atmos_param/qe_moist_convection/Makefile
        src/atmos_param/diag_integral/Makefile
        src/atmos_param/topo_drag/Makefile
        src/atmos_param/lscloud_driver/Makefile
        src/atmos_param/tke_turb/Makefile
        src/atmos_param/edt/Makefile
        src/atmos_param/entrain/Makefile
        src/atmos_param/macrophysics/Makefile
        src/atmos_param/vert_turb_driver/Makefile
        src/atmos_param/shallow_cu/Makefile
        src/atmos_param/moist_processes/Makefile
        src/atmos_param/damping_driver/Makefile
        src/atmos_param/radiation/Makefile
        src/atmos_param/radiation/longwave/Makefile
        src/atmos_param/radiation/aerosols/Makefile
        src/atmos_param/radiation/shortwave/Makefile
        src/atmos_param/radiation/util/Makefile
        src/atmos_param/radiation/clouds/Makefile
        src/atmos_param/radiation/driver/Makefile
        src/atmos_param/radiation/radiative_gases/Makefile
        src/atmos_param/clubb/Makefile
        src/atmos_param/clubb/CLUBB_core/Makefile
        src/atmos_param/my25_turb/Makefile
        src/atmos_cubed_sphere/Makefile
        src/shared/Makefile
        src/ice_sis/Makefile
        src/atmos_drivers/Makefile
        src/atmos_drivers/solo/Makefile
        src/atmos_drivers/coupled/Makefile
        exec/Makefile
        exec/templates/Makefile
        exec/coupler/Makefile
        exec/fms/Makefile
        exec/land_lad2/Makefile
        exec/mom6/Makefile
        exec/atmos_dyn/Makefile
        exec/atmos_phys/Makefile
        exec/ice_sis/Makefile
        run/Makefile
        ])
AC_OUTPUT()
