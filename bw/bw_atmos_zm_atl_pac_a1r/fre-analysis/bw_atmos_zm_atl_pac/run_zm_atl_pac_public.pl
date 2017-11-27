#!/usr/bin/perl
#------------------------------------------------------------------
# Creates a 3-panel plot of zonally averaged eastward wind stress.
# The 3 panels are for the Pacific Ocean, Atlantic Ocean, and
# global land. In addition to the model curve, CMIP5 enesemble models
# and observations are also plotted for comparison.
#------------------------------------------------------------------
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
Getopt::Long::Configure("bundling");

my @required_modules = qw/fre fre-analysis ncarg/;

# make sure the required modules are loaded from parent shell
my $platform = `gfdl_platform`; chomp $platform;
my $message = ""; 
my $public = "yes";
if (($platform eq "desktop") && ($public eq "no")) {
   my $error = 0;
   foreach my $module (@required_modules) {
      my @M = grep /^$module\/.*/, split/:/, $ENV{"LOADEDMODULES"};
      if (@M == 0) {
         print "ERROR: module $module not loaded\n";
         $message .= "\nmodule load $module";
         $message .= "/6.2.1" if ($module eq "ncarg");
         $error++;
      }   
   }   
   if ($error) {
      print "Try running:$message\n";
      exit 1;
   }   
} else {
   if($public eq "no") {die "Invalid platform";}
}
#-----------------------------------------------
# argument parsing for monthly-mean input data
#-----------------------------------------------

my %Opt = (HELP=>0, SETUP_ONLY=>0, VERBOSE=>0, QUIET=>0, USE_ALL_OBS=>0, variable=>"tauu");
my $status = GetOptions ('h|help'           => \$Opt{HELP},
                         'S|setup_only'     => \$Opt{SETUP_ONLY},
                         'V|verbose'        => \$Opt{VERBOSE},
                         'Q|quiet'          => \$Opt{QUIET},
                         'A|all_obs'        => \$Opt{USE_ALL_OBS},
                         'i|in_data_dir=s'  => \$Opt{in_data_dir},
                         'o|out_dir=s'      => \$Opt{out_dir},
                         'd|descriptor=s'   => \$Opt{descriptor},
                         'y|year_range=s'   => \$Opt{year_range},
                         'v|variable=s'     => \$Opt{variable});

#-------------------------
# argument error checking
#-------------------------

if (!$Opt{in_data_dir} && !$Opt{HELP}) {
   print "ERROR: no argument given for input directory\n";
   $Opt{HELP} = 1;
}

if (!$Opt{year_range} && !$Opt{HELP}) {
   print "ERROR: no argument given for yr1,yr2\n";
   $Opt{HELP} = 1;
} else {
   # user supplied year range 
   my @yrs = split /,/, $Opt{year_range};
   if (@yrs == 2) {
      ($Opt{yr1},$Opt{yr2}) = @yrs;
   } else {
      print "ERROR: invalid entry for yr1,yr2\n";
      $Opt{HELP} = 1;
   }
}
my $yr_range_label = $Opt{yr1}."-".$Opt{yr2};

# input files (must be 12 files)
my @in_files;
if (!$Opt{HELP}) {
   if (@ARGV == 12) {
      @in_files = @ARGV;
   } elsif (@ARGV == 0) {
      # try to determine file names
      my $component_name;
      if ($Opt{in_data_dir} =~ /pp.*\/(.*?)\/av\//) {
         $component_name = $1;
      }
      if ($component_name) {
         foreach my $mon (qw/01 02 03 04 05 06 07 08 09 10 11 12/) {
            my $monfile = join ".", ($component_name, $yr_range_label, $mon, "nc");
            push @in_files, $monfile if (-e $Opt{in_data_dir}."/".$monfile);
         }
      } else {
         print "ERROR: could not determine component name\n";
      }
      if (@in_files != 12) {
         print "ERROR: could not determine number of input files\n";
         $Opt{HELP} = 1;
      }
   } else {
      print "ERROR: incorrect number of files\n";
      $Opt{HELP} = 1;
   }
}


#--------------------------
# usage message then exit
#--------------------------

if ($Opt{HELP}) {
   usage();
   exit 1;
}

#-------------------------------
# check that input files exist
#-------------------------------

my @dmfiles;
foreach my $file (@in_files) {
   my $dmfile = $Opt{in_data_dir}."/".$file;
   push @dmfiles,$dmfile if (-e $dmfile);
}
die "ERROR: some input files do not exist" if (@dmfiles != 12);


#-------------------------------
# variable(s) to be processed
#-------------------------------
die "ERROR: invalid variable: ".$Opt{variable} if ($Opt{variable} ne "tauu" && $Opt{variable} ne "curl" && $Opt{variable} ne "pr");
my @varlist = $Opt{variable};
@varlist = qw/tauu tauv/ if ($Opt{variable} eq "curl");

#-------------------------------
# directory locations
#-------------------------------

# abs_path to determine package location -- more portable
my $package_location = substr(abs_path($0),0,rindex(abs_path($0),"/"));

# set environment variable for access in NCL code
$ENV{BW_PACKAGE_ROOT} = $package_location;

 my $BW_ANALYSIS_ARCHIVE = $ENV{FRE_ANALYSIS_ARCHIVE};
my $archDir = $ENV{FRE_ANALYSIS_ARCHIVE} . "/CMIP5";

my $vortScript = "$package_location/ncl/compute_vorticity.ncl";
my $climScript = "$package_location/ncl/compute_monthly_climatology.ncl";
my $plotScript = "$package_location/ncl/plot_zonal_clim.ncl";
my $label_script = "$package_location/shared/bw/addpagelabels.pl";

my $workDir = $ENV{TMPDIR} . "/bw_".$Opt{variable}."_atl_pac";
mkdir $workDir if (!-e $workDir);
chdir $workDir;

#----------------------------------
# copy data
#----------------------------------

# choose year range for CMIP5 model (several choices)
my $cmip5dates = "197901-200812";
#$cmip5dates = "198101-200012" if ($yr_range_label eq "1981-2000");
#$cmip5dates = "198101-200012" if ($yr_range_label eq "1981-1990");
# starting ending year for plotting script
my $cmip5_yrbeg = substr($cmip5dates,0,4);
my $cmip5_yrend = substr($cmip5dates,7,4);

# check if ensemble mean files have been computed
# if they have then skip
my @climFiles;
if (-d "cmip5") {
  @climFiles = glob "cmip5/".$Opt{variable}."_Amon_*_r0i0p0_$cmip5dates-clim.nc";
}
print "For variable \"".$Opt{variable}."\" there were \"".scalar(@climFiles)."\" local CMIP5 climo files found.\n";
  
# list of climo files for all cmip5 models
my @cmip5files;
foreach my $var (@varlist) {
  if (@climFiles == 0) {
    push @cmip5files, "$archDir/$var\_climo/$var\_Amon_CMIP5_$cmip5dates-clim.nc.tar";
     }
     }

# climatology for cmip5 models (copied to directory "models")
if (@cmip5files) {
  my @commands = ("#--- Copying CMIP5 Climatology ---");
  push @commands, "gcp @cmip5files .";
  push @commands, "mkdir -p models";
  foreach my $cmip5file (@cmip5files) {
    push @commands, "tar -xvf ".tailname($cmip5file)." -C models";
    push @commands, "rm -f ".tailname($cmip5file);
  }
  execute_commands(@commands);

  # compute curl of tau (vorticity)
  if ($Opt{variable} eq "curl") {
    my @commands;
    foreach my $ufile ( glob "models/tauu_Amon_*_$cmip5dates-clim.nc" ) {
      my $vfile = $ufile; $vfile =~ s/\/tauu_/\/tauv_/;
      if (-e $vfile) {
        my $ofile = $ufile; $ofile =~ s/\/tauu_/\/curl_/;
        push @commands, compute_curl_tau($ufile,$vfile,$ofile);
      }
    }
    execute_commands(@commands);
  }

  #-----------------------------------------
  # compute ensemble mean for cmip5 models
  #-----------------------------------------

  my $var = $Opt{variable};

  # gather runs for each model
  my %models;
  opendir(DIR,"models");
  while (my $file = readdir(DIR)) {
    if ($file =~ /$var\_Amon_(.*)_amip_.*\.nc/) {
      my @models = keys %models;
      if (!grep{$_ eq $1} @models) {
        $models{$1} = [$file];
      } else {
        my $aref = $models{$1};
        push @$aref,$file;
        $models{$1} = $aref;
      }
    }
  }

  # compute ensemble means
  my @commands = ("mkdir -p cmip5");
  foreach my $model (keys %models) {
    my $runs = $models{$model};
    my $ename = @$runs[0]; $ename =~ s/_r\d{1,2}i\d{1,2}p\d{1,2}_/_r0i0p0_/;
    push @commands, "ncea -h -p models @$runs cmip5/$ename";
  }
  execute_commands(@commands);
}


#-------------------------------------------------------
# original observed data (copied to directory "obs")
#-------------------------------------------------------

# use the cmip5 starting and ending year
my $yrbeg = $cmip5_yrbeg;
my $yrend = $cmip5_yrend;

# mapping
my %map_var = ( merra => { sftlf=>"land_mask", tauu=>"taux", tauv=>"tauy", pr=>"prectot" },
                ecmwf => { sftlf=>"lsm",       tauu=>"iews", tauv=>"inss" } );
my %map_dates = ( merra=>"197901-201212",
                  ecmwf=>"197901-201212" );

my @dmget_files;
my @commands;
foreach my $obs (qw/merra ecmwf/) {
  # skip precip for ecmwf
  next if ($Opt{variable} eq "pr" && $obs eq "ecmwf");

  my $climVar = $Opt{variable};
 #$climVar = "wind_stress" if ($Opt{variable} eq "tauu" || $Opt{variable} eq "curl");
  my $climFile = "obs/$climVar.$obs.climatology.nc";

  # only create climatology file if it does not exist
  if (!-e $climFile) {
    push @commands, "mkdir -p obs/$obs";

    # create land/sea mask as fraction
    my $sftlf = $map_var{$obs}{sftlf};
    push @dmget_files, "$BW_ANALYSIS_ARCHIVE/data/monthly/$obs/$sftlf.nc.gz";
    push @commands, "gcp $BW_ANALYSIS_ARCHIVE/data/monthly/$obs/$sftlf.nc.gz obs/$obs/";
    push @commands, "gunzip obs/$obs/$sftlf.nc.gz";
    push @commands, "ncrename -h -v $sftlf,sftlf obs/$obs/$sftlf.nc obs/$obs/sftlf.nc";
   #push @commands, "ncap2 -v -s \"sftlf=$sftlf*100.f\" obs/$obs/$sftlf.nc obs/$obs/sftlf.nc";

    # rename variable/files to CMIP5 names
    my $dates = $map_dates{$obs};
    foreach my $var (@varlist) {
      my $dvar = $map_var{$obs}{$var};
      push @dmget_files, "$BW_ANALYSIS_ARCHIVE/data/monthly/$obs/$dvar.$dates.nc.gz";
      push @commands, "gcp $BW_ANALYSIS_ARCHIVE/data/monthly/$obs/$dvar.$dates.nc.gz obs/$obs/";
      push @commands, "gunzip obs/$obs/$dvar.$dates.nc.gz";
      push @commands, "ncrename -h -v $dvar,$var obs/$obs/$dvar.$dates.nc";
      push @commands, "mv obs/$obs/$dvar.$dates.nc obs/$obs/$var.$dates.nc";
    }

    my $tsFile = "obs/$obs/".$Opt{variable}.".$dates.nc";

    # special cases for variables
    if ($Opt{variable} eq "curl") {
      push @commands, "cp obs/$obs/tauu.$dates.nc $tsFile";
      push @commands, "ncks -h -A obs/$obs/tauv.$dates.nc $tsFile";
    }

    # add land/sea mask to timeseries file
    push @commands, "ncks -h -A obs/$obs/sftlf.nc $tsFile";

    # create climatology file
    my @clim_options;
    push @clim_options, "\'ifile=\"$tsFile\"\'";
    push @clim_options, "\'ofile=\"$climFile\"\'";
    push @clim_options, "yr1=$yrbeg";
    push @clim_options, "yr2=$yrend";
    push @commands, "ncl @clim_options $climScript";

    # compute curl of wind stress
    if ($Opt{variable} eq "curl") {
      push @commands, compute_curl_tau($climFile,$climFile,"obs/$obs/vorticity.nc");
      # move vort file to climo
      push @commands, "mv $climFile obs/$obs/wind_stress.$obs.climatology.nc";
      push @commands, "mv obs/$obs/vorticity.nc $climFile";
    }
  }
}

# get the older wind stress observation files (COADS & ERS)
if (($Opt{variable} eq "tauu" || $Opt{variable} eq "curl") && $Opt{USE_ALL_OBS}) {
  my $obsfile = "wind_stress.obs.climatology.nc.tar.gz";
  if (!-e "obs/wind_stress.coads.climatology.nc" || !-e "obs/wind_stress.ers.climatology.nc") {
    push @dmget_files, "$archDir/$obsfile";
    push @commands, "gcp $archDir/$obsfile $obsfile";
    push @commands, "tar -xvf $obsfile obs/wind_stress.{coads,ers}.climatology.nc";
    push @commands, "rm -f $obsfile";
  }
  if ($Opt{variable} eq "curl") {
    foreach my $obs (qw/coads ers/) {
      push @commands, compute_curl_tau("obs/wind_stress.$obs.climatology.nc","obs/wind_stress.$obs.climatology.nc","obs/curl.$obs.climatology.nc");
    }
  }
}
# get the older precip observation files (CMAP & GPCP)
if ($Opt{variable} eq "pr") {
  my $obsfile = "precip.obs.climatology.nc.tar.gz";
  if (!-e "obs/pr.cmap.climatology.nc" || !-e "obs/pr.gpcp.climatology.nc") {
    push @dmget_files, "$archDir/$obsfile";
    push @commands, "gcp $archDir/$obsfile $obsfile";
    push @commands, "tar -xvf $obsfile obs/pr.{cmap,gpcp}.climatology.nc";
    push @commands, "rm -f $obsfile";
  }
}

# execute the commands for observed data
if (@commands) {
  print "#--- Copying Observed Climatology: ---\n" if !$Opt{QUIET};
  execute_commands("dmget @dmget_files") if @dmget_files;
  execute_commands(@commands);
}

###############
# model data
###############

if (!-e "model_climatology.nc") {

  my %mapping = ( tauu=>"tau_x", pr=>"precip", tauv=>"tau_y" );

  # execute the commands to retrieve the model data
  my @commands = ("#--- Copying Model Climatology ---",
                  "dmget @dmfiles",
                  "gcp @dmfiles $workDir");
  execute_commands(@commands);

  #---- extract the variables needed ----
  my @gvars = model_variable_name($in_files[0],\@varlist);
  my $gvarlist = join(",",@gvars) . "," . additional_variables($in_files[0],\@gvars);
  my $lfrac = model_mask_name($in_files[0]); #"land_mask";

  my @commands;
  # extract land-sea mask
  if ($lfrac eq "land_mask") {
   #push @commands, "ncap2 -v -s \"sftlf=$lfrac*100.f\" ".$in_files[0]." sftlf.nc";
    push @commands, "ncks -h -v $lfrac ".$in_files[0]." sftlf.nc";
    push @commands, "ncrename -h -v $lfrac,sftlf sftlf.nc";
    push @commands, "ncatted -h -a units,sftlf,m,c,\"%\" -a valid_range,sftlf,d,, sftlf.nc";
  } else {
    push @commands, "ncks -h -v $lfrac ".$in_files[0]." sftlf.nc";
  }
  # extract all variables/bounds needed (concatenating into an annual climo)
  push @commands, "ncrcat -h -v $gvarlist @in_files model_climatology.nc";
  # flip sign on CMIP-named wind stresses to match GFDL and observed wind stress
  my @ncap_cmds;
  push @ncap_cmds, "tauu=-tauu" if (grep{$_ eq "tauu"} @gvars);
  push @ncap_cmds, "tauv=-tauv" if (grep{$_ eq "tauv"} @gvars);
  if (@ncap_cmds) {
    my $ncap_cmds = join(";",@ncap_cmds);
    push @commands, "mv model_climatology.nc tmp.model_climatology.nc";
    push @commands, "ncap2 -h -s \"$ncap_cmds\" tmp.model_climatology.nc model_climatology.nc";
    push @commands, "rm -f tmp.model_climatology.nc";
  }
  # rename GFDL variables to CMIP names (if needed)
  foreach my $var (@varlist) {
    if (!grep{$_ eq $var} @gvars) {
      my $gvar = $mapping{$var};
      push @commands, "ncrename -h -v $gvar,$var model_climatology.nc";
    }
  }
  # append land-sea mask into climatology
  push @commands, "ncks -A -h sftlf.nc model_climatology.nc";
  push @commands, "rm -f @in_files";
  push @commands, "rm -f sftlf.nc";

  # compute curl of wind stress
  if ($Opt{variable} eq "curl") {
    push @commands, "mv model_climatology.nc model_tau_climatology.nc";
    push @commands, compute_curl_tau("model_tau_climatology.nc","model_tau_climatology.nc","model_climatology.nc");
  }

  # execute the commands for model data
  execute_commands(@commands);
}

#--------------------------------------------
# determine model descriptor if not supplied
# used as label on plot
#--------------------------------------------

if (!$Opt{descriptor}) {
  $Opt{descriptor} = determine_model_descriptor("model_climatology.nc",$Opt{in_data_dir});
  if (!$Opt{descriptor}) {
    print "WARNING: no model descriptor, using generic term \"Model\"\n" if !$Opt{QUIET};
    $Opt{descriptor} = "Model";
  }
}
#-------------------------------
# model descriptor + year range
#-------------------------------
my $model_descriptor = $Opt{descriptor}." (".$yr_range_label.")";


# exit if setting up data only
exit if $Opt{SETUP_ONLY};

#####################
# run the graphics
#####################

print "#--- Running Plotting Routines ---\n" if !$Opt{QUIET};
my $output_file = $Opt{variable}."_atl_pac.ps";
my @options;
push @options, "yr1=$cmip5_yrbeg";
push @options, "yr2=$cmip5_yrend";
push @options, "\'model_descriptor=\"$model_descriptor\"\'";
push @options, "\'psfile=\"$output_file\"\'";
push @options, "\'variable=\"".$Opt{variable}."\"\'";
push @options, "ALL_OBS=True" if $Opt{USE_ALL_OBS};
print  "ncl @options $plotScript\n" if !$Opt{QUIET};
system("ncl @options $plotScript");

# modify the page labels
print  "$label_script -O -l ANN,DJF,MAM,JJA,SON $output_file\n" if !$Opt{QUIET};
system("$label_script -O -l ANN,DJF,MAM,JJA,SON $output_file");

# save the file
if ($Opt{out_dir}) {
   print "--- Saving Postscript File ---\n" if !$Opt{QUIET};
   my $out_dir = $Opt{out_dir} . "/atmos_" . $Opt{yr1} . "_" . $Opt{yr2};
   $out_dir = "gfdl:".$out_dir if ($out_dir =~ /^\/net\// || $out_dir =~ /^\/net2\// || $out_dir =~ /^\/home\// || $out_dir =~ /^\/nbhome\//);
   # compress the file
   print  "gzip  $output_file\n" if !$Opt{QUIET};
   system("gzip  $output_file");
   my $zfile = $output_file; $zfile =~ s/\.ps$/.ps.gz/;
   print  "gcp -cd $zfile $out_dir/Wyman.".$Opt{variable}."_atl_pac/\n" if !$Opt{QUIET};
   system("gcp -cd $zfile $out_dir/Wyman.".$Opt{variable}."_atl_pac/");
}
 

###########################################
###########################################
###########################################

# help message
sub usage {
   my $cmdname = substr($0,rindex($0,"/")+1);
   print "

[1mUSAGE[0m
   $cmdname -i [4mINDIR[0m -y [4mYR1[0m,[4mYR2[0m [OPTIONS...]

";
   return 0;
}

#---------------------------------------------
# return string after last period (i.e., the suffix)
sub tailname {
   my $tail = shift;
   while ($tail =~ s/^.*\///) {}
   return $tail;
}

#---------------------------------------------
# return the string before last period including path (i.e., the trunk)
sub trunkname {
   my $trunk = shift;
   return substr($trunk,0,rindex($trunk,"."));
}

#---------------------------------------------
# execute the commands for model data
sub execute_commands {
  foreach my $cmd (@_) {
    print "$cmd\n" if !$Opt{QUIET};
    system($cmd) if ($cmd !~ /^#/);
    exit 1 if $?;
  }
  return 1;
}

#---------------------------------------------
# determine additional variables in file
# axis bounds and "old" time average info
sub additional_variables {
  my $file = shift;
  my $vars = shift; # reference
  my $dump = `ncdump -h $file`;
  my $DEBUG_LOCAL = 0;
  my @addvars;
  foreach my $var (@$vars) {
    if ($dump =~ /\t\w+ $var\((.+)\) ;/) {
      print STDERR "Found var $var\n" if $DEBUG_LOCAL;

      # search for dimension bounds
      foreach my $dim (split /,\s*/, $1) {
        print STDERR "Found dim $dim\n" if $DEBUG_LOCAL;
        if ($dump =~ /\t\t$dim:bounds = "(.*)"/) {
          my $bnd = $1;
          if ($dump =~ /\t\w+ $bnd\(.+\)/) {
            push @addvars, $bnd if !grep /^$bnd$/, @addvars;
          }
        } elsif ($dump =~ /\t\t$dim:climatology = "(.*)"/) {
          my $bnd = $1;
          if ($dump =~ /\t\w+ $bnd\(.+\)/) {
            push @addvars, $bnd if !grep /^$bnd$/, @addvars;
          }
        }
      }   

      # search for time average info variables
      if ($dump =~ /\t\t$var:time_avg_info = "(.*)"/) {
        foreach my $tavg (split /,/, $1) {
          if ($dump =~ /\t\w+ $tavg\(/) {
            push @addvars, $tavg if !grep /^$tavg$/, @addvars;
          }
        }
      }
    }
  }

  return join ",", @addvars;
}

#---------------------------------------------
# determine model variable name
# input is cmip name
# also check old gfdl name
sub model_variable_name {
  my $file = shift;
  my $vref = shift;
  my @gvars;
  my %mapping = ( tauu=>"tau_x", pr=>"precip", tauv=>"tau_y" );
  my $dump = `ncdump -h $file`;
  foreach my $var (@$vref) {
    if ($dump =~ /\t\w+ $var\((.+)\) ;/) {
      push @gvars, $var;
    } elsif ($mapping{$var}) {
      my $gvar = $mapping{$var};
      push @gvars, $gvar if ($dump =~ /\t\w+ $gvar\((.+)\) ;/);
    }
  }
  die "ERROR in model_variable_name: variable name(s) not recognized" if (scalar(@$vref) ne scalar(@gvars));
  return @gvars;
}

#---------------------------------------------
# determine model land-sea name
sub model_mask_name {
  my $file = shift;
  my $mask;
  my $dump = `ncdump -h $file`;
  foreach my $var (qw/land_mask sftlf/) {
    $mask = $var if ($dump =~ /\t\w+ $var\(.+\) ;/);
  }
  die "ERROR in model_mask_name: neither \"land_mask\" or \"sftlf\" exists in file $file" if !$mask;
  return $mask;
}

#---------------------------------------------
# compute curl of winds

sub compute_curl_tau($$$;$$) {
  my ($ufile,$vfile,$ofile,$uvar,$vvar) = @_;
  $uvar = "tauu" if !$uvar;
  $vvar = "tauv" if !$vvar;
  my @commands;
  my @options;
  push @options, "\'ufile=\"$ufile\"\'";
  push @options, "\'vfile=\"$vfile\"\'";
  push @options, "\'ofile=\"$ofile\"\'";
  push @options, "\'uvar=\"$uvar\"\'";
  push @options, "\'vvar=\"$vvar\"\'";
  push @options, "\'maskvar=\"sftlf\"\'";
  push @commands, "ncl @options $vortScript";  # vortScript is global
  push @commands, "ncrename -h -v vort,curl $ofile";
  return @commands;
}

#---------------------------------------------
# determine model descriptor name
# (as used in the figure label)

sub determine_model_descriptor {
  my $file = shift;
  my $idir = shift;
  my $descriptor;
  if (`ncdump -h $file` =~ /:title = "(.*)" ;/) {
     $descriptor = $1; 
  } else {
    # second, try to read it from the path name
    if ($idir =~ /.*\/(.*?)\/(.*?)\/pp.*\//) {
       my $desc1 = $1; my $desc2 = $2; 
       if ($2 =~ /^gfdl.ncrc/) {
          $descriptor = $desc1;
       } else {
          $descriptor = $desc2;
       }   
    }   
  }
  return $descriptor;
}


