#!/usr/bin/perl
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
Getopt::Long::Configure("bundling");

# make sure the required platform and module are loaded

my $platform = "desktop"; chomp $platform;
my $public = "yes";
if (($platform eq "desktop") && ($public eq "no")) {
   my $error = required_modules( qw/fre fre-analysis ncarg ifort/ );
   die "Error: Program terminated" if $error;
} else {
if($public eq "no") {  die "Error: Invalid platform";}
}


#-------------------
# argument parsing
#-------------------

my %Opt = (HELP=>0, SETUP_ONLY=>0, VERBOSE=>0, FORCE_COPY=>0, STATSFILE=>0);

my $status = GetOptions ('h|help'           => \$Opt{HELP},
                         'S|setup_only'     => \$Opt{SETUP_ONLY},
                         'F|force_copy'     => \$Opt{FORCE_COPY},
                         'V|verbose'        => \$Opt{VERBOSE},
                         'f|statsfile!'     => \$Opt{STATSFILE},
                         'i|in_data_dir=s'  => \$Opt{in_data_dir},
                         'o|out_dir=s'      => \$Opt{out_dir},
                         'd|descriptor=s'   => \$Opt{descriptor},
                         's|staticfile=s'   => \$Opt{staticfile},
                         'y|year_range=s'   => \$Opt{year_range},
                         'c|data_chunks=s'  => \$Opt{data_chunks},
                         'g|regrid=s'       => \$Opt{regrid},
                         't|stats=s'        => \$Opt{stats},
                         'r|regions=s'      => \$Opt{regions},
                         'v|variable=s'     => \$Opt{variable});

usage() if $Opt{HELP};

#-------------------
# argument checking
#-------------------

# required arguments

if (!$Opt{in_data_dir}) {
   print STDERR "ERROR: no argument given for input directory\n";
   $Opt{HELP} = 1;
}

if (!$Opt{year_range}) {
   print STDERR "ERROR: no argument given for yr1,yr2\n";
   $Opt{HELP} = 1;
}

if (!$Opt{variable}) {
   print STDERR "ERROR: no argument given for variable\n";
   $Opt{HELP} = 1;
}

usage() if $Opt{HELP};


#-------------------------
# process input arguments
#-------------------------

# user supplied year range

if ( $Opt{year_range} =~ /(\d+),(\d+)/ ) {
   $Opt{yr1} = $1; $Opt{yr2} = $2;
} else {
   die "ERROR: invalid entry for yr1,yr2";
}

# determine data chunks (the years in each time series file)

if ($Opt{data_chunks}) {
   my @yrs = split /,/, $Opt{data_chunks};
   if (@yrs == 1) {
      $Opt{databegyr} = $Opt{yr1};
      $Opt{dataendyr} = $Opt{yr2};
      $Opt{datachunk} = $yrs[0];
   } elsif (@yrs == 3) {
      ($Opt{databegyr},$Opt{dataendyr},$Opt{datachunk}) = @yrs;
   } else {
      die "ERROR: invalid entry for data_chunk";
   }
} else {
   $Opt{databegyr} = $Opt{yr1};
   $Opt{dataendyr} = $Opt{yr2};
   $Opt{datachunk} = ($Opt{dataendyr} - $Opt{databegyr}) + 1
}
print "databegyr,dataendyr,datachunk: ".$Opt{databegyr}.", ".$Opt{dataendyr}.", ".$Opt{datachunk}."\n" if $Opt{VERBOSE};


# error check: data years must include analysis years

if ( $Opt{yr1} < $Opt{databegyr} || $Opt{yr2} > $Opt{dataendyr} ) {
   die "ERROR: requested model analysis years (yr1,yr2) not within data years (databegyr,dataendyr)";
}

my $year_range = sprintf("%4.4d-%4.4d",($Opt{yr1},$Opt{yr2}));

# regridding method must be 'bilinear' or 'conserve'
if ( $Opt{regrid} ) {
   if ( $Opt{regrid} ne "bilinear" && $Opt{regrid} ne "conserve" ) {
      die "ERROR: regrid type must be either \'bilinear\' or \'conserve\'";
   }
}

# region file must exist
if ( $Opt{regions} ) {
  die "ERROR: region file \'".$Opt{regions}."\' does not exist" if (! -e $Opt{regions});
}

#-------------------------------
# determine the component name
#-------------------------------
my $component_name;
if ($Opt{in_data_dir} =~ /pp.*\/(.*?)\/ts\//) {
  $component_name = $1;
  print "model component: $component_name\n" if $Opt{VERBOSE};
}

#-------------------------------
# determine the data frequency
#-------------------------------
my $frequency;
if ($Opt{in_data_dir} =~ /\/ts\/(.*?)\//) {
  $frequency = $1;
  print "data frequency: $frequency\n" if $Opt{VERBOSE};
  die "ERROR: invalid data frequency" if ($frequency ne "monthly" && $frequency ne "daily");
} else {
  die "ERROR: data frequency could not be determined";
}

#-------------------------------
# variable name(s)
#-------------------------------
# currently only near surface temperature and wet days for daily data
# all variables for monthly data

my @varlist;
if ( grep /^$Opt{variable}$/, qw/tas pr wet clt/ ) { 
  die "ERROR: Wet days requires daily data" if ($frequency eq "monthly" && $Opt{variable} eq "wet");
  push @varlist, $Opt{variable};
  push @varlist, ("tasmax","tasmin") if ($Opt{variable} eq "tas");
 #push @varlist, $Opt{variable}      if  $frequency eq "monthly";
 #push @varlist, ("tasmax","tasmin") if ($frequency eq "daily" && $Opt{variable} eq "tas");
 #push @varlist, ("wet")             if ($frequency eq "daily" && $Opt{variable} eq "wet");
} else {
  die "ERROR: invalid variable requested: ".$Opt{variable};
}
die "ERROR: unable to determine variables needed, var=".$Opt{variable}.", freq=$frequency" if !@varlist;

#----------------------------------------
# determine static file (for land mask)
#----------------------------------------
if (!$Opt{staticfile}) {
 #if ($Opt{in_data_dir} =~ /^(.*)\/.*\/.*\/.*$/) {
  if ($Opt{in_data_dir} =~ /^(.*\/$component_name)\/.*$/) {
    $Opt{staticfile} = "$1/$component_name.static.nc";
    die "ERROR: static file could not be located" if (!-e $Opt{staticfile});
  } else {
    die "ERROR: could not determine static file name";
  }
  print "static file: ".$Opt{staticfile}."\n" if $Opt{VERBOSE};
} else {
  die "ERROR: static file does not exist" if (!-e $Opt{staticfile});
}

#------------------
# script locations
#------------------
my $package_location = substr(abs_path($0),0,rindex(abs_path($0),"/"));
my $wetday_script = "$package_location/ncl/compute_wet_days.ncl";
my $timeseries_script = "$package_location/ncl/compute_monthly_timeseries.ncl";
my $ncl_script = "$package_location/ncl/plot_ts_data.ncl";
my $date_script = "$package_location/shared/bw/file_dates.pl";
my $label_script = "$package_location/shared/bw/addpagelabels.pl";

# set environment variable for access in NCL code
$ENV{BW_PACKAGE_ROOT} = $package_location;

# put color map directory in search path
$ENV{NCARG_COLORMAPS} = "$package_location/ncl/colormaps:".$ENV{NCARG_ROOT}."/lib/ncarg/colormaps";

# ncl version number
my $nclversion = `ncl -V`; chomp $nclversion;

# directories
die "ERROR: TMPDIR nor defined: invalid platform" if !$ENV{TMPDIR};
my $tmpdir = $ENV{TMPDIR};
my $workdir = "$tmpdir/cruts31";
mkdir $workdir if (!-e $workdir);
chdir $workdir;
my $subdirectory = "Wyman.cru_ts";
# location of observed datasets
my $dataPath = $ENV{FRE_ANALYSIS_ARCHIVE}."/data/monthly/cru_ts_322";

######################
#  model data setup
######################

my %cvarmap = ( tas=>"tas", tasmin=>"tasmin", tasmax=>"tasmax", pr=>"pr", wet=>"pr", clt=>"clt" );
my %gvarmap = ( tas=>"t_ref", tasmin=>"t_ref_min", tasmax=>"t_ref_max", pr=>"precip", wet=>"precip", clt=>"tot_cld_amt" );

# dates for file names (monthly data)
my @date_list = split /\n/, `$date_script -f $frequency -Y $Opt{yr1} -Z $Opt{yr2} $Opt{databegyr} $Opt{dataendyr} $Opt{datachunk}`;
# full_date_range = YYYYMM-YYYYMM
my $full_date_range = $date_list[0] . $date_list[$#date_list];
$full_date_range =~ s/-\d{12}-/-/ if $frequency eq "monthly";
if ($frequency eq "daily") {
  $full_date_range =~ s/-\d{16}-/-/;
  $full_date_range =~ s/01-/-/; $full_date_range =~ s/31$//;  # remove days
}

mkdir "model" if (!-e "model");

# dmget time series files
my @dmget_list = ();
my @missing_list = ();
foreach my $var (@varlist) {
  my $tsfile = "model/$var.$full_date_range.nc";
  if (!-e $tsfile || $Opt{FORCE_COPY}) {
    my $cvar = $cvarmap{$var};
    my $gvar = $gvarmap{$var};
    foreach my $date (@date_list) {
      my $basefile = $Opt{in_data_dir} . "/$component_name.$date";
      if (-e "$basefile.$cvar.nc") {
        push @dmget_list, "$basefile.$cvar.nc";  # cmip-named file
      } elsif (-e "$basefile.$gvar.nc") {
        push @dmget_list, "$basefile.$gvar.nc";  # gfdl-named file
      } else {
        if ( grep /^$var$/, qw/tasmin tasmax tas/ ) { 
          print "WARNING: Missing model file: $basefile.{$var,$gvar}.nc\n";
          push @missing_list, $var;
        } else {
          die "ERROR: Missing model file: $basefile.{$var,$gvar}.nc";
        }
      }
    }
  }
}
# land mask file (add if something to copy)
push @dmget_list, $Opt{staticfile} if (@dmget_list && !-e "model/sftlf.nc");


# exclude variable from the list if missing
print "varlist (before): @varlist\n" if @missing_list && $Opt{VERBOSE};
foreach my $miss (@missing_list) {
  foreach my $index (grep { $varlist[$_] eq $miss } 0 .. $#varlist) { 
    splice(@varlist, $index, 1);
  }
}
print "varlist (after): @varlist\n" if @missing_list && $Opt{VERBOSE};

#----------------------
# dmget and copy files
#----------------------
if (@dmget_list) {
  my @commands;
  push @commands, "# dmgetting ".scalar(@dmget_list)." files";
  push @commands, "dmget @dmget_list";
  push @commands, "# copying files";
  push @commands, "gcp @dmget_list $workdir";
  print  "#------------ model data ------------\n";
  execute_commands(@commands);

  # extract land mask field from static file
  if (!-e "model/sftlf.nc") {
    my @commands;
    push @commands, "# extract land mask";
    if (isfilevar("$component_name.static.nc","sftlf")) {
      push @commands, "ncks -h -v sftlf $component_name.static.nc model/sftlf.nc";
    } elsif (isfilevar("$component_name.static.nc","land_mask")) {
      # rename land_mask to sftlf
      push @commands, "ncks -h -v land_mask $component_name.static.nc model/sftlf.nc";
      push @commands, "ncrename -h -v land_mask,sftlf model/sftlf.nc";
      my @ncatted_opts = ( "-a long_name,sftlf,m,c,\"Land Area Fraction\"",
                           "-a units,sftlf,m,c,\"1.0\"", "-a valid_range,sftlf,d,," );
      push @commands, "ncatted -h @ncatted_opts model/sftlf.nc";
    } else {
      die "ERROR: Missing land-sea mask in file $component_name.static.nc";
    }
    execute_commands(@commands);
  }

  #---------------------
  # loop thru variables
  #---------------------
  foreach my $var (@varlist) {
    my $cvar = $cvarmap{$var};
    my $gvar = $gvarmap{$var};
    my @commands;
    my @files;
    # loop thru consecutive time series files
    foreach my $date (@date_list) {
      my $ofile = "$component_name.$date.$var.nc";
      my ($ifile,$ivar);
      if (-e "$component_name.$date.$cvar.nc") {
        $ifile = "$component_name.$date.$cvar.nc";
        $ivar = $cvar;
      } elsif (-e "$component_name.$date.$gvar.nc") {
        $ifile = "$component_name.$date.$gvar.nc";
        $ivar = $gvar;
      } else {
        die "ERROR: Could not determine input file: $component_name.$date.{$cvar,$gvar}.nc";
      }

      # compute wet day frequency
      if ($var eq "wet") {
         my $file = 
         push @commands, "echo NCL Version $nclversion";
         push @commands, "ncl -Q \'ifile=\"$ifile\"\' \'ofile=\"$ofile\"\' \'var=\"$ivar\"\' \'crit=\"1.0\"\' $wetday_script";
         push @commands, "rm -f $ifile";
      # otherwise change to cmip variable names
      } elsif ($ivar ne $var) { 
        push @commands, "mv $ifile $ofile";
        push @commands, "ncrename -h -v $ivar,$var $ofile";
      }

      # compute monthly time series from daily
      if ($frequency eq "daily") {
        my $mfile = $ofile; $mfile =~ s/01-/-/; $mfile =~ s/31\./\./;
        push @commands, "echo NCL Version $nclversion";
        push @commands, "ncl -Q \'ifile=\"$ofile\"\' \'ofile=\"$mfile\"\' $timeseries_script";
        push @commands, "rm -f $ofile";
        push @files, $mfile;
      } else {
        push @files, $ofile;
      }
    }
    # concatenate monthly timeseries file into long timeseries
    my $tsfile = "model/$var.$full_date_range.nc";
    if (@files > 1 || $files[0] ne $tsfile) {
      push @commands, "ncrcat -h @files $tsfile";
      push @commands, "rm -f @files";
    }
    # append land mask field from static file
    push @commands, "ncks -h -A model/sftlf.nc $tsfile";
    # execute all the commands for this variable
    execute_commands(@commands);
  }
}

# auxilary files (mean & diff)
#if ($frequency eq "daily" && grep(/^tasmax$/, @varlist) && grep(/^tasmin$/, @varlist)) {
if (grep(/^tasmax$/, @varlist) && grep(/^tasmin$/, @varlist)) {
   my @commands = compute_auxilary_files("model",$full_date_range,"model/sftlf.nc");
   if (@commands) {
     print "#------------ auxilary variables -----------\n";
     execute_commands(@commands);
  }
}

# cleanup/remove the land mask files
if (-e "$component_name.static.nc") {
  my @commands = ("# clean up static file");
  push @commands, "rm -f $component_name.static.nc";
  execute_commands(@commands);
}

#-----------------------------------------------
# determine model name (used as label in plot)
#-----------------------------------------------
if (!$Opt{descriptor}) {
  # first, try to read it from the model climatology file
  my $tsfile = "model/".$varlist[0].".$full_date_range.nc";
  if (`ncdump -h $tsfile` =~ /:title = "(.*)" ;/) {
     $Opt{descriptor} = $1;
  } else {
    # second, try to read it from the path name
    if ($Opt{in_data_dir} =~ /.*\/(.*?)\/(.*?)\/pp.*\//) {
       my $desc1 = $1; my $desc2 = $2;
       if ($2 =~ /^gfdl.ncrc/) {
          $Opt{descriptor} = $desc1;
       } else {
          $Opt{descriptor} = $desc2;
       }
    }
  }
   $Opt{descriptor} = "Model" if !$Opt{descriptor};
   print "model descriptor: ".$Opt{descriptor}."\n" if $Opt{VERBOSE};
}


#########################
#  observed data setup
#########################
# need to adjust dates for climo if they are not in years of obs data
my ($oyr1,$oyr2) = obsyears($Opt{yr1},$Opt{yr2});
my $year_range_obs = sprintf("%4.4d-%4.4d",($oyr1,$oyr2));
my $full_date_range_obs = sprintf("%4.4d01-%4.4d31",($oyr1,$oyr2));

my %varmap = ( tasmin=>"tmn", tasmax=>"tmx", tas=>"tmp", pr=>"pre", wet=>"wet", clt=>"cld" );
my $obsdir = "cru";
mkdir $obsdir;

my @dmget_list;
my @commands;
foreach my $var (@varlist) {
   my $obsfile = "$obsdir/$var.$full_date_range_obs.nc";
   if (!-e $obsfile) {
      if ($varmap{$var}) {
         my $ks = ($oyr1-1901)*12;
         my $ke = ($oyr2-1901)*12 + 11;
         my $rawfile = "cru_ts_322.1901-2013.".$varmap{$var}.".nc";
         push @dmget_list, "$dataPath/$rawfile";
         push @commands, "ncks -d time,$ks,$ke $rawfile $obsfile";
         push @commands, "ncrename -h -v ".$varmap{$var}.",$var $obsfile" if ($varmap{$var} ne $var);
         push @commands, "rm -f $rawfile";
      } else {
         die "ERROR: invalid variable name with no CRU Mapping: $var";
      }
   }
}
if (@dmget_list) {
   print  "#------------ observed data ------------\n";
   print  "gcp @dmget_list .\n";
   system("gcp @dmget_list .");

   execute_commands(@commands);
}

# auxilary files (mean & diff)
#if ($frequency eq "daily" && grep(/^tasmax$/, @varlist) && grep(/^tasmin$/, @varlist)) {
if (grep(/^tasmax$/, @varlist) && grep(/^tasmin$/, @varlist)) {
  my @commands = compute_auxilary_files($obsdir,$full_date_range_obs);
  if (@commands) {
    print "#------------ auxilary variables -----------\n";
    execute_commands(@commands);
  }
}

exit if $Opt{SETUP_ONLY};

# add temperature mean and diff 
#if ($frequency eq "daily" && grep(/^tasmax$/, @varlist) && grep(/^tasmin$/, @varlist)) {
if (grep(/^tasmax$/, @varlist) && grep(/^tasmin$/, @varlist)) {
    push @varlist, qw/tasmean tasdiff/;
}

# append atmos_yr1_yr2 to out_dir
my $outdir = $Opt{out_dir} . "/atmos_" . $Opt{yr1} . "_" . $Opt{yr2} if $Opt{out_dir};

#######################
#  plotting section
#######################

# compile code
if (!-e "horiz_interp.so") {
   my @commands = ("mkdir -p wrapit");
   push @commands, "gcp $package_location/shared/bw/ncl/wrapit/horiz_interp{.stub,_float.f90,_double.f90,.inc} wrapit";
   push @commands, "WRAPIT -in -n horiz_interp wrapit/horiz_interp{.stub,_float.f90,_double.f90}";
   execute_commands(@commands);
}

foreach my $var (@varlist) {
   my $mfile = "model/$var.$full_date_range.nc";
   my $ofile = "cru/$var.$full_date_range_obs.nc";
   die "ERROR: file not found: $mfile" if (!-e $mfile);
   die "ERROR: file not found: $ofile" if (!-e $ofile);
   my @options;
   push @options, "\'var=\"$var\"\'";
   push @options, "\'mfile=\"$mfile\"\'";
   push @options, "\'ofile=\"$ofile\"\'";
   push @options, "\'mlab=\"".$Opt{descriptor}."\"\'" if $Opt{descriptor};
   push @options, "\'olab=\"".uc($obsdir)." TS 3.22\"\'";    # <<<==== HARD-CODED CRU DATA VERSION
   push @options, "\'mdate=\"$year_range\"\'";
   push @options, "\'odate=\"$year_range_obs\"\'";
  #push @options, "\'rfile=\"$package_location/ncl/regions_us.txt\"\'" if ($dataset eq "prism");
   push @options, "\'rfile=\"".$Opt{regions}."\"\'" if $Opt{regions};
   push @options, "\'regrid=\"".$Opt{regrid}."\"\'" if $Opt{regrid};
   push @options, "\'stats=\"".$Opt{stats}."\"\'" if $Opt{stats};
   push @options, "\'statsfile=\"statistics.".uc(substr($frequency,0,1)).".txt\"\'" if $Opt{STATSFILE};
   print "NCL Version $nclversion\n";
   print  "ncl -Q @options $ncl_script\n";
   system("ncl -Q @options $ncl_script");

   # add labels
   foreach my $plot ( glob "$var.*.ps" ) {
      print  "$label_script -O -l ANN,DJF,MAM,JJA,SON $plot\n";
      system("$label_script -O -l ANN,DJF,MAM,JJA,SON $plot");
   }

   # move to output directory?
   if ($Opt{out_dir}) {
      my @gcp_files;
      foreach my $plot ( glob "$var.*.ps" ) {
         print  "gzip -f $plot\n";
         system("gzip -f $plot");
         push @gcp_files, $plot.".gz";
      }
      print  "gcp -cd @gcp_files $outdir/$subdirectory/$var/\n";
      system("gcp -cd @gcp_files $outdir/$subdirectory/$var/");
   }
         
}

# save statistics file (if it exists)
if ($Opt{STATSFILE} && $Opt{out_dir}) {
  my $F = uc(substr($frequency,0,1));
  if (-e "statistics.$F.txt") {
    # gcp will overwrite older stats file
    print  "gcp -cd statistics.$F.txt $outdir/$subdirectory/statistics/\n";
    system("gcp -cd statistics.$F.txt $outdir/$subdirectory/statistics/");
  }
}


########################################################################
##################    E N D   O F   S C R I P T    #####################
########################################################################

# input: /aaa/bbb/ccc.xxx
# returns: ccc.xxx
sub tailname {
   my $tail = shift;
   while ($tail =~ s/^.*\///) {}
   return $tail;
}

#----------------------------

sub required_modules {
   my $message = ""; 
   my $error = 0;
   my @current_modules = split/:/, $ENV{"LOADEDMODULES"};
   foreach my $module (@_) {
      my @M = grep /^$module\/.*/, @current_modules;
      if (@M == 0) {
         print STDERR "ERROR: module $module not loaded\n";
         $message .= "\nmodule load $module";
         $message .= "/6.1.2" if ($module eq "ncarg");  # special case
         $error++;
      }   
   }   
   print STDERR "Try running:$message\n" if $error;
   return $error;  # status code
}
    
#----------------------------

sub compute_auxilary_files ($$;$) {
   my ($dir,$yrange,$maskfile)= @_;
   my @commands;
   # auxilary files (mean & diff)
   if (!-e "$dir/tasmean.$yrange.nc" || !-e "$dir/tasdiff.$yrange.nc") {
      if (-e "$dir/tasmax.$yrange.nc" && -e "$dir/tasmin.$yrange.nc") {
         push @commands, "cp $dir/tasmin.$yrange.nc $dir/tasminmax.$yrange.nc";
         push @commands, "ncks -h -A $dir/tasmax.$yrange.nc $dir/tasminmax.$yrange.nc";
         if (!-e "$dir/tasmean.$yrange.nc") {
            push @commands, "ncap2 -hv -s \'tasmean=(tasmax+tasmin)/2\' $dir/tasminmax.$yrange.nc $dir/tasmean.$yrange.nc";
            push @commands, "ncatted -h -a long_name,tasmean,m,c,\"near-surface temperature mean\" $dir/tasmean.$yrange.nc";
            push @commands, "ncks -h -A $maskfile $dir/tasmean.$yrange.nc" if $maskfile;
         }
         if (!-e "$dir/tasdiff.$yrange.nc") {
            push @commands, "ncap2 -hv -s \'tasdiff=tasmax-tasmin\' $dir/tasminmax.$yrange.nc $dir/tasdiff.$yrange.nc";
            push @commands, "ncatted -h -a long_name,tasdiff,m,c,\"near-surface temperature diurnal range\" $dir/tasdiff.$yrange.nc";
            push @commands, "ncks -h -A $maskfile $dir/tasdiff.$yrange.nc" if $maskfile;
         }
         push @commands, "rm -f $dir/tasminmax.$yrange.nc";
      } else {
         die "ERROR: required files do not exist: $dir/tasmax.$yrange.nc && $dir/tasmin.$yrange.nc";
      }
   }
   return @commands;
}

#----------------------------

sub execute_commands {
  foreach my $cmd (@_) {
    print "$cmd\n";
   #print "EXEC> $cmd\n" if ($cmd !~ /^#/);
    system ($cmd) if ($cmd !~ /^#/);
  }
}

#----------------------------

sub obsyears {
   my ($iyrbeg,$iyrend) = @_; 
   my ($yrbeg,$yrend) = ($iyrbeg,$iyrend);
   my ($databeg,$dataend) = (1901,2013);

   if ($yrbeg < $databeg) {
      $yrend = $yrend + ($databeg-$yrbeg);
      $yrbeg = $databeg;
   }   
   if ($yrend > $dataend) {
      my $myrs = $yrbeg-$databeg;
      my $nyrs = $yrend-$dataend;
      $myrs = $nyrs if ($nyrs < $myrs);
      $yrbeg = $yrbeg - $myrs if ($myrs > 0); 
      $yrend = $dataend;
   }   
   print STDERR "NOTE: Adjusted starting date from $iyrbeg to $yrbeg\n" if ($yrbeg != $iyrbeg);
   print STDERR "NOTE: Adjusted ending date from $iyrend to $yrend\n" if ($yrend != $iyrend);
   return ($yrbeg,$yrend);
}

#----------------------------
# return logical whether a variable is in a file

sub isfilevar {
   my ($file,$var) = @_; 
   my $dump = `ncdump -h $file`;
   if ($dump =~ /\t\w+ $var\(.+\)/) {
     return 1;
   } else {
     return 0;
   }
}

#----------------------------

sub usage {
   my $cmdname = substr($0,rindex($0,"/")+1);
   print "

[1mOVERVIEW[0m
   Generates 3-panel comparison plots of Model versus CRU (TS 3.1) data. 
   Global and regional plots are generated.  Currently only precipitation is available.

[1mUSAGE[0m
   $cmdname -i [4mINDIR[0m -y [4mYR1[0m,[4mYR2[0m -v [4mLIST[0m [OPTIONS...]

[1mREQUIRED OPTIONS[0m
   [1m-i[0m, [1m--in_data_dir[0m [4mPATH[0m
         Path for the input directory containing the time series of the daily pressure level fields.

   [1m-y[0m, [1m--year_range[0m [4mYR1[0m,[4mYR2[0m
         The year range for the analyzed data.
         This can be a subset of the years of avaliable data.

   [1m-v[0m, [1m--variable[0m [4mLIST[0m
         Variable to be plotted.
         Possible daily timeseries data variables are: 'tas' or 'wet'.
         Possible monthly timeseries data variables are: 'tas', 'pr', or 'clt'.

[1mOPTIONAL[0m

   [1m-s[0m, [1m--staticfile[0m [4mPATH[0m
         Full path name of the file containing static fields. The land mask will be
         extracted from this file.

   [1m-d[0m, [1m--descriptor[0m [4mNAME[0m
         A descriptive label placed on the figure, usually the experiemnt name.
         If not given then if possible [4mNAME[0m will be determined from the
         input directory path.
 
   [1m-o[0m, [1m--out_dir[0m [4mPATH[0m
         Path for the output directory where the figures are saved.
         The directory \"atmos_YR1_YR2\" will be appended to this path.
         If [1m-o[0m [4mPATH[0m is not present then the figures remain in /ftmp.

   [1m-c[0m, [1m--year_range[0m [4mDATABEGYR[0m,[4mDATAENDYR[0m,[4mDATACHUNK[0m
         Defines the file structure in the the input directory. [4mDATABEGYR[0m is the first year,
         [4mDATAENDYR[0m is the last year, and [4mDATACHUNK[0m is the file size in years.
         If [4mDATABEGYR[0m,[4mDATAENDYR[0m,[4mDATACHUNK[0m are not given, then they are derived from
         [4mYR1[0m,[4mYR2[0m with the assumption that all the data will be analyzed and that there
         is one file. If the optional form of this is used [1m-c[0m [4mDATACHUNK[0m, then
         [4mDATABEGYR[0m and [4mDATAENDYR[0m will be set to [4mYR1[0m and [4mYR2[0m, respectively.

   [1m-V[0m, [1m--verbose[0m
         Output messages about data processing and progress.

[1mEXAMPLE[0m
   $cmdname -i /archive/user/myExperName/pp/atmos_daily_plev/ts/daily/1yr \
            -o /net/user/figures/myExperName/analysis \
            -d myExperName -y 1981,2000 -c 1981,2000,10

";
   exit 1;
}


