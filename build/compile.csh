#!/bin/tcsh -fx

# Compile Script 

set -r echoOn = $?echo

if ( $echoOn ) unset echo
echo "<NOTE> : Starting at $HOST on `date`"
if ( $echoOn ) set echo

#unalias *

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

# ---------------- write main Makefile

sed -e 's/<TAB>/\t/' >$bld_dir/Makefile <<END
# Makefile for Experiment 'cm4p12_warsaw'

SRCROOT = $src_dir/
BUILDROOT = $bld_dir/
STAGEDIR = $ptmp_dir/$bld_dir/

MK_TEMPLATE = $mkmf_template
include \$(MK_TEMPLATE)

fms_cm4p12_warsaw.x: coupler/libcoupler.a ice_sis/libice_sis.a atmos_dyn/libatmos_dyn.a land_lad2/libland_lad2.a atmos_phys/libatmos_phys.a mom6/libmom6.a fms/libfms.a
<TAB>\$(LD) \$^ \$(LDFLAGS) -o \$@ \$(STATIC_LIBS)

fms/libfms.a:  FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=fms \$(@F)

atmos_phys/libatmos_phys.a: fms/libfms.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=atmos_phys \$(@F)

atmos_dyn/libatmos_dyn.a: atmos_phys/libatmos_phys.a fms/libfms.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=atmos_dyn \$(@F)

ice_sis/libice_sis.a: mom6/libmom6.a fms/libfms.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=ice_sis \$(@F)

land_lad2/libland_lad2.a: fms/libfms.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=land_lad2 \$(@F)

mom6/libmom6.a: fms/libfms.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE) OPENMP="" --directory=mom6 \$(@F)

coupler/libcoupler.a: atmos_dyn/libatmos_dyn.a ice_sis/libice_sis.a atmos_phys/libatmos_phys.a mom6/libmom6.a land_lad2/libland_lad2.a fms/libfms.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=coupler \$(@F)

FORCE:

stage:
<TAB>install -d \$(STAGEDIR)
<TAB>install -m 555 fms_cm4p12_warsaw.x \$(STAGEDIR)

clean:
<TAB>\$(MAKE) --directory=fms clean
<TAB>\$(MAKE) --directory=atmos_phys clean
<TAB>\$(MAKE) --directory=atmos_dyn clean
<TAB>\$(MAKE) --directory=ice_sis clean
<TAB>\$(MAKE) --directory=land_lad2 clean
<TAB>\$(MAKE) --directory=mom6 clean
<TAB>\$(MAKE) --directory=coupler clean

localize:
<TAB>\$(MAKE) -f \$(BUILDROOT)fms/Makefile localize
<TAB>\$(MAKE) -f \$(BUILDROOT)atmos_phys/Makefile localize
<TAB>\$(MAKE) -f \$(BUILDROOT)atmos_dyn/Makefile localize
<TAB>\$(MAKE) -f \$(BUILDROOT)ice_sis/Makefile localize
<TAB>\$(MAKE) -f \$(BUILDROOT)land_lad2/Makefile localize
<TAB>\$(MAKE) -f \$(BUILDROOT)mom6/Makefile localize
<TAB>\$(MAKE) -f \$(BUILDROOT)coupler/Makefile localize

distclean:
<TAB>\$(RM) -r fms
<TAB>\$(RM) -r atmos_phys
<TAB>\$(RM) -r atmos_dyn
<TAB>\$(RM) -r ice_sis
<TAB>\$(RM) -r land_lad2
<TAB>\$(RM) -r mom6
<TAB>\$(RM) -r coupler
<TAB>\$(RM) fms_cm4p12_warsaw.x
<TAB>\$(RM) Makefile

END

# ---------------- create component Makefiles

mkdir -p $bld_dir/fms
list_paths -o $bld_dir/fms/pathnames_fms $src_dir/shared
cd $bld_dir
pushd fms
mkmf -m Makefile -a $src_dir -b $bld_dir -p libfms.a -t $mkmf_template -g -c "-DINTERNAL_FILE_NML -g -Duse_libMPI -Duse_netCDF" -Imom6/src/MOM6/pkg/CVMix-src/include -Ishared/include -Ishared/mpp/include $bld_dir/fms/pathnames_fms
popd

mkdir -p $bld_dir/atmos_phys
list_paths -o $bld_dir/atmos_phys/pathnames_atmos_phys $src_dir/atmos_param $src_dir/atmos_shared
cd $bld_dir
pushd atmos_phys
mkmf -m Makefile -a $src_dir -b $bld_dir -p libatmos_phys.a -t $mkmf_template -g -c "-DINTERNAL_FILE_NML -g" -o "-I$bld_dir/fms" -Imom6/src/MOM6/pkg/CVMix-src/include -Ishared/include -Ishared/mpp/include $bld_dir/atmos_phys/pathnames_atmos_phys
popd

mkdir -p $bld_dir/atmos_dyn
list_paths -o $bld_dir/atmos_dyn/pathnames_atmos_dyn $src_dir/atmos_drivers/coupled $src_dir/atmos_cubed_sphere/driver/coupled $src_dir/atmos_cubed_sphere/model $src_dir/atmos_cubed_sphere/model_nh_null $src_dir/atmos_cubed_sphere/tools $src_dir/atmos_cubed_sphere/GFDL_tools
cd $bld_dir
pushd atmos_dyn
mkmf -m Makefile -a $src_dir -b $bld_dir -p libatmos_dyn.a -t $mkmf_template -g -c "-DINTERNAL_FILE_NML -g -DSPMD -DCLIMATE_NUDGE" -o "-I$bld_dir/atmos_phys -I$bld_dir/fms" -Imom6/src/MOM6/pkg/CVMix-src/include -Ishared/include -Ishared/mpp/include $bld_dir/atmos_dyn/pathnames_atmos_dyn
popd

mkdir -p $bld_dir/ice_sis
list_paths -o $bld_dir/ice_sis/pathnames_ice_sis $src_dir/ice_sis $src_dir/ice_param
cd $bld_dir
pushd ice_sis
mkmf -m Makefile -a $src_dir -b $bld_dir -p libice_sis.a -t $mkmf_template -g -c "-DINTERNAL_FILE_NML -g -Duse_netCDF" -o "-I$bld_dir/mom6 -I$bld_dir/fms" -Imom6/src/MOM6/pkg/CVMix-src/include -Ishared/include -Ishared/mpp/include $bld_dir/ice_sis/pathnames_ice_sis
popd

mkdir -p $bld_dir/land_lad2
list_paths -o $bld_dir/land_lad2/pathnames_land_lad2 $src_dir/land_lad2
cd $bld_dir
pushd land_lad2
mkmf -m Makefile -a $src_dir -b $bld_dir -p libland_lad2.a -t $mkmf_template --use-cpp -g -c "-DINTERNAL_FILE_NML -g -nostdinc " -o "-I$bld_dir/fms" -Imom6/src/MOM6/pkg/CVMix-src/include -Ishared/include -Ishared/mpp/include $bld_dir/land_lad2/pathnames_land_lad2
# -I/usr/include -I/usr/lib64/gfortran/modules
popd

mkdir -p $bld_dir/mom6
list_paths -o $bld_dir/mom6/pathnames_mom6 $src_dir/mom6/src/MOM6/config_src/dynamic $src_dir/mom6/src/MOM6/config_src/coupled_driver $src_dir/mom6/src/MOM6/src/*/ $src_dir/mom6/src/MOM6/src/*/*/ $src_dir/ocean_shared/generic_tracers $src_dir/ocean_shared/mocsy/src
cd $bld_dir
pushd mom6
mkmf -m Makefile -a $src_dir -b $bld_dir -p libmom6.a -t $mkmf_template -g -c "-DINTERNAL_FILE_NML -g -DMAX_FIELDS_=100 -DNOT_SET_AFFINITY -D_USE_MOM6_DIAG -D_USE_GENERIC_TRACER -DUSE_PRECISION=2 -D_FILE_VERSION="'"`git-version-string $<`"'"" -o "-I$bld_dir/fms" -Imom6/src/MOM6/pkg/CVMix-src/include -Ishared/include -Ishared/mpp/include $bld_dir/mom6/pathnames_mom6
popd

mkdir -p $bld_dir/coupler
list_paths -o $bld_dir/coupler/pathnames_coupler $src_dir/coupler
cd $bld_dir
pushd coupler
mkmf -m Makefile -a $src_dir -b $bld_dir -p libcoupler.a -t $mkmf_template -g -c "-DINTERNAL_FILE_NML -g" -o "-I$bld_dir/atmos_dyn -I$bld_dir/ice_sis -I$bld_dir/atmos_phys -I$bld_dir/mom6 -I$bld_dir/land_lad2 -I$bld_dir/fms" -Imom6/src/MOM6/pkg/CVMix-src/include -Ishared/include -Ishared/mpp/include $bld_dir/coupler/pathnames_coupler
popd

# ---------------- call make on the main Makefile

make  OPENMP=on NETCDF=3 fms_cm4p12_warsaw.x

if ( $status == 0 ) then
  if ( $?NiNaC_LVL ) then
    if ( $NiNaC_LVL > 0 ) then
      # Run NiNaC
      $NiNaC_BldRx $src_dir $bld_dir
      if ( $status != 0 ) then
        if ( $echoOn ) unset echo
        echo "NiNaC Note: While NiNaC loaded attempt at NiNaC_BldRx failed with exit status $status : FRE continuing as normal."
        if ( $echoOn ) set echo
      endif
    endif
  endif

  if ( $echoOn ) unset echo
  echo "<NOTE> : make succeeded for cm4p12_warsaw."
  if ( $echoOn ) set echo
else
  if ( $echoOn ) unset echo
  echo "*ERROR*: make failed for cm4p12_warsaw."
  if ( $echoOn ) set echo
  exit 1
endif

exit 0
