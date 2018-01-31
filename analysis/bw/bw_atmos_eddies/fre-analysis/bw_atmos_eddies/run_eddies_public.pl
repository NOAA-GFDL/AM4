#!/usr/bin/perl
#-----------------------------------------------------------------------
# Computes and plots eddy variances and covariances for both model and
# observation (MERRA or NCEP Reanalysis) for standard atmospheric fields
# on several pressure levels.
#
# Input data must be daily time series with CMIP-type names.
# The standard CMIP6 input files with GFDL-style naming will be looked for.
#-----------------------------------------------------------------------
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
Getopt::Long::Configure("bundling");

# internal flags
my $FAST = 1;
my $TEST = 0;

# default plots and observed data sets
my @default_plots = qw/ uu10 tt10 tt250 uu250 uv250 vt250 vv250 zz250 tt500 wt500 zz500 ww500 uu850 vv850 tt850 qq850 uv850 vt850 vq850 wt850 wq850 ww850 zz850 /;
my @datasets = qw/ merra ncep2 /;

# required modules
my @required_modules = qw/fre ncarg/;
push @required_modules, "ifort" if $FAST;

# make sure the required modules are loaded from parent shell
my $platform = "desktop"; chomp $platform;
my $message = "";
my $public = "yes";
if (($platform eq "desktop") && ($public eq "no")) {
  my @error;
  foreach my $module (@required_modules) {
    push @error, $module if !grep /^$module\/.*/, split/:/, $ENV{"LOADEDMODULES"};
  }
  if (@error) {
    print STDERR "ERROR: the following modules are not loaded: @error\n"; die;
  }
} else {
if($public eq "no") {  die "Error: Invalid platform";}
}

#-----------------------
#  argument parsing
#-----------------------

my %Opt = (HELP=>0, SETUP_ONLY=>0, SKIP_SETUP=>0, VERBOSE=>0, WARN=>1, filter_opt=>0);
my $status = GetOptions ('h|help'           => \$Opt{HELP},
                         'S|setup_only'     => \$Opt{SETUP_ONLY},
                         'R|run_only'       => \$Opt{SKIP_SETUP},
                         'V|verbose+'       => \$Opt{VERBOSE},
                         'f|filter_opt=i'   => \$Opt{filter_opt},
                         'i|in_data_dir=s'  => \$Opt{in_data_dir},
                         'o|out_dir=s'      => \$Opt{out_dir},
                         'd|descriptor=s'   => \$Opt{descriptor},
                         's|staticfile=s'   => \$Opt{staticfile},
                         'y|year_range=s'   => \$Opt{year_range},
                         'c|data_chunks=s'  => \$Opt{data_chunks},
                         'D|datasets=s'     => \$Opt{datasets},
                         'p|plots=s'        => \$Opt{plots} );

#-------------------------
# argument error checking
#-------------------------

if (!$Opt{in_data_dir} && !$Opt{HELP}) {
   print STDERR "ERROR: no argument given for input directory\n";
   $Opt{HELP} = 1;
}

if (!$Opt{year_range} && !$Opt{HELP}) {
   print STDERR "ERROR: no argument given for yr1,yr2\n";
   $Opt{HELP} = 1;
} else {
   # user supplied year range 
   my @yrs = split /,/, $Opt{year_range};
   if (@yrs == 2) {
      ($Opt{yr1},$Opt{yr2}) = @yrs;
   } else {
      print STDERR "ERROR: invalid entry for yr1,yr2\n";
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

#------------------------------------------------------------
# determine data chunks (the years in each time series file)
#------------------------------------------------------------

if ($Opt{data_chunks}) {
  ($Opt{databegyr},$Opt{dataendyr},$Opt{datachunk}) = determine_data_chunks ($Opt{data_chunks},$Opt{yr1},$Opt{yr2});
} else {
  ($Opt{databegyr},$Opt{dataendyr},$Opt{datachunk}) = ($Opt{yr1},$Opt{yr2},$Opt{yr2}-$Opt{yr1}+1)
}
print "databegyr,dataendyr,datachunk: ".$Opt{databegyr}.", ".$Opt{dataendyr}.", ".$Opt{datachunk}."\n" if $Opt{VERBOSE} > 0;

#---------------------------------------------
# directory locations for this package
# package root is the location of this script
#----------------------------------------------

my $package_location = substr(abs_path($0),0,rindex(abs_path($0),"/"));
$ENV{BW_PACKAGE_ROOT} = $package_location;

# put color map directory in search path
$ENV{NCARG_COLORMAPS} = "$package_location/ncl/colormaps:".$ENV{NCARG_ROOT}."/lib/ncarg/colormaps";

# observed data directory
 my $obs_data_location = "/home/a1r/AM4/bw/fms/fre-analysis/test/DET/data/daily";
 print "$obs_data_location ";

#my $obs_data_location = "/archive/bw/fre-analysis/test/bw/bw_atmos_eddies/daily";
#my $obs_data_location = "/net2/bw/data/merra/daily/levels_3d/v";

# external scripts
my $date_script = "$package_location/shared/bw/file_dates.pl";
my $fvalue_script = "$package_location/shared/bw/ncfillvalue.pl";
my $label_script = "$package_location/shared/bw/pagelabels.pl";
my $ncl_script = "$package_location/ncl/eddies.ncl";
my $extract_script = "$package_location/ncl/extract.ncl";

# local directory locations
my $workdir = $ENV{TMPDIR} . "/bw_atmos_eddies";
mkdir $workdir if (!-e $workdir);
chdir $workdir;
mkdir "model" if (!-e "model");

#-------------------------------
# determine the component name
#-------------------------------

my $component_name = get_component_name($Opt{in_data_dir});
print "model component: $component_name\n" if $Opt{VERBOSE} > 0;

#----------------------------------------------
# determine model name (used as label in plot)
#----------------------------------------------

if (!$Opt{descriptor}) {
   $Opt{descriptor} = get_exper_descriptor($Opt{in_data_dir});
   print "model descriptor: ".$Opt{descriptor}."\n" if $Opt{VERBOSE} > 0;
}

#----------------------------------------------
# static file name (if not supplied)
#----------------------------------------------

if (!$Opt{staticfile}) {
   $Opt{staticfile} = get_staticfile_name($Opt{in_data_dir},$component_name);
   print "static file: ".$Opt{staticfile}."\n" if $Opt{VERBOSE} > 0;
}

#----------------------------------------------
# check frequency of input data (must be daily)
#----------------------------------------------

my $frequency = "null";
$frequency = $1 if ($Opt{in_data_dir} =~ /\/ts\/(.*)\/\d+yr/);
print "frequency: $frequency\n" if $Opt{VERBOSE} > 0 || $frequency ne "daily";
die "ERROR: input data set not daily timeseries data" if ($frequency ne "daily");

#------------------------
# which data sets to use
#------------------------

@datasets = split /,/, $Opt{datasets} if $Opt{datasets};
print STDOUT "datasets: @datasets\n" if $Opt{VERBOSE} > 0;

#-------------------------
# list of requested plots
#-------------------------
my @plot_list = @default_plots;
@plot_list = split /,/, $Opt{plots} if $Opt{plots};
my @plot_list_out;

# translate 'z' to 'h' in all plot requests
foreach (@plot_list) { s/z/h/g }

#-----------------------------------------------------------
# mapping between plot variable letters and variable names
#-----------------------------------------------------------

# between plot name -> 3D model/cmip variable name
my %cmip_letter_map = ( q=>"hus", t=>"ta", u=>"ua", v=>"va", w=>"wap", h=>"zg" );
my $cmip_suffix = "_unmsk";

# between plot name -> 3D observed variable name
my %obs_letter_map = ( q=>"qv", t=>"t", u=>"u", v=>"v", w=>"omg", h=>"h" );

#-------------------------
# NCL version information
#-------------------------

my $nclversion = `ncl -V`; chomp $nclversion;
die "ERROR: NCL version must be 6.2.1 or greater;" if ($nclversion lt "6.2.1");
print "NCL Version $nclversion\n" if $Opt{VERBOSE} > 0;


##############################################
########     PROCESS MODEL DATA      #########
##############################################

my @dmget_files;
my @plot_variables;
my @plots_skipped;
my %extracted_levels;

#------------------------------------
# list of date ranges for each file
#------------------------------------
my @date_list = split /\n/, `$date_script daily $Opt{databegyr} $Opt{dataendyr} $Opt{datachunk}`;
my $full_date_range = $date_list[0] . $date_list[$#date_list]; $full_date_range =~ s/-\d{16}//;
if ($Opt{VERBOSE} > 0) {
  print "Model Dates:\n";
  foreach (@date_list) {
    print "  $_\n";
  }
}

#-----------------------------------------------------------
# loop through requested plot 
# plot name is of the form: ab####
# where a & b are variable names, #### is the pressure level
#-----------------------------------------------------------
foreach my $plot (@plot_list) {
  my $lev = substr($plot,2);
  foreach my $pvar ( substr($plot,0,1), substr($plot,1,1) ) {

    # check for valid mapped letter
    if ($cmip_letter_map{$pvar}) {
      my $cvar = $cmip_letter_map{$pvar}.$cmip_suffix;  # cmip variable name
      my $ovar = $obs_letter_map{$pvar};                # obs variable name
      my $pvarlev = $pvar.$lev;                         # plot variable (extracted on level)

     #if (! -e "model/$component_name.$full_date_range.$pvarlev.nc") {
      if (! -e "model/$pvarlev.$full_date_range.nc") {
        my @arfiles = check_for_archive_files($Opt{in_data_dir},$component_name,\@date_list,$cvar);
        if (@arfiles == @date_list) {
          foreach my $arfile (@arfiles) {
            if (! -e "model/".tailname($arfile)) {
              push @dmget_files, $arfile if (! grep /$arfile/, @dmget_files);
            }
          }
          if (!$extracted_levels{$pvarlev}) {
            $extracted_levels{$pvarlev}{var} = $cvar;
            $extracted_levels{$pvarlev}{files} = [tailname(@arfiles)];
          }
          push @plot_list_out,  $plot    if (! grep /$plot/,    @plot_list_out);
          push @plot_variables, $pvarlev if (! grep /$pvarlev/, @plot_variables);
        } else {
          print STDERR "WARNING: archive files found for \"$cvar\" ... skipping plot $plot\n";
          push @plots_skipped, $plot if (! grep /$plot/, @plots_skipped);
        }
      } else {
        push @plot_list_out,  $plot    if (! grep /$plot/,    @plot_list_out);
        push @plot_variables, $pvarlev if (! grep /$pvarlev/, @plot_variables);
      }
    } else {
      die "ERROR: Invalid variable letter \"$pvar\" in plot field, $plot";
    }
  }
}

#-------------------------------------
# add pmask2 and/or ps to dmget list
#-------------------------------------
my $pres_vars_found = 0;
my %psfiles;
foreach my $var (qw/pmaskv2 ps/) {
  if (! -e "model/$var.$full_date_range.nc") {
    my @arfiles = check_for_archive_files($Opt{in_data_dir},$component_name,\@date_list,$var);
    if (@arfiles == @date_list) {
      foreach my $arfile (@arfiles) {
        if (! -e "model/".tailname($arfile)) {
          push @dmget_files, $arfile if (! grep /$arfile/, @dmget_files);
        }
      }
      $psfiles{$var} = [tailname(@arfiles)];
      $pres_vars_found++;
    }
  } else {
    $pres_vars_found++;
  }
}

# create an appropriate ps mask using zsurf (add staticfile to dmget list)
my $FLAG_ps_from_zs;
if ($pres_vars_found <= 0) {
  if (!-e "model/ps_est.nc") {
    if ($Opt{staticfile}) {
      push @dmget_files, $Opt{staticfile};
      $FLAG_ps_from_zs = 1;
    } else {
      die "Missing model ps files and no static file defined";
    }
  }
  print STDERR "WARNING: estimating surface pressure for masking using topography\n";
}


#---------------------------------------------------------
# dmget and copy the model files to working directory
#---------------------------------------------------------

if (@dmget_files) {
  print "dmgetting ".scalar(@dmget_files)." model data files\n";
  if ($Opt{VERBOSE} > 0) {
    foreach (@dmget_files) {
      print "  $_\n";
    }
  }
  system ("dmget @dmget_files") if !$TEST;
  print "copying files\n";
  print("gcp @dmget_files $workdir/model");
  system ("gcp @dmget_files $workdir/model") if !$TEST;
}

#-------------------------------------------
# remove skipped plots from the plot list
#-------------------------------------------
if (@plots_skipped) {
   print "plots skipped: @plots_skipped\n";
   foreach my $plot (@plots_skipped) {
      my $index = 0;
      $index++ until $plot_list[$index] eq $plot; 
      splice(@plot_list, $index, 1);
   }
   print "plot list @plot_list\n";
}

#---------------------------------------------------------------
# create approx psurf file (when pmaskv2 and ps are not found)
# use millibars since ncl code uses variable name for level
#---------------------------------------------------------------
if ($FLAG_ps_from_zs) {
   if (!-e "model/ps_est.nc") {
     my $cmd = "ps=1000*((288-zsurf*5e-3)/288)^6.852";
     print "ncap2 -v -s \"$cmd\" ".tailname($Opt{staticfile})." model/ps_est.nc\n" if ($TEST || $Opt{VERBOSE} > 0);
     system ("ncap2 -v -s \"$cmd\" ".tailname($Opt{staticfile})." model/ps_est.nc")  if !$TEST;
     unlink tailname($Opt{staticfile}) if !$TEST;
   }
}

# concatenate ps files
if ($pres_vars_found > 0) {
  foreach my $var (keys %psfiles) {
    my @files = @{$psfiles{$var}};
    foreach (@files) { $_ = "model/".$_; }  # prepend model directory
    concatenate_files (@files,"model/$var.$full_date_range.nc");
  }
}

#--------------------------------------
# EXTRACT LEVELS FROM CMIP 3D FILES
#--------------------------------------

extract_level_from_3d_data( "model", \%extracted_levels );


########################################################################
###############  OBSERVED DATA  ########################################
########################################################################

foreach my $data (@datasets) {
  my @dmget_files;
  my %extracted_levels;
  mkdir $data if (!-e $data);

  my $ardir = "$obs_data_location/$data";
  my ($oyr1,$oyr2) = determine_observation_years($ardir,$Opt{yr1},$Opt{yr2});
  my @date_list_obs = determine_obs_file_dates($ardir,$oyr1,$oyr2);
  my $obs_date_range = $date_list_obs[0] . $date_list_obs[$#date_list_obs]; $obs_date_range =~ s/-\d{16}//;
  if ($Opt{VERBOSE} > 0) {
    print uc($data)." Dates:\n";
    foreach (@date_list_obs) {
      print "  $_\n";
    }
    print "Observed ($data) date range: $obs_date_range\n";
  }

  # try to get/extract all plot variables available from model
  foreach my $pvarlev (@plot_variables) {
    my $pvar = substr($pvarlev,0,1);
    my $lev = substr($pvarlev,1);
    my $ovar = $obs_letter_map{$pvar};

    if (! -e "$data/$pvarlev.$obs_date_range.nc") {
      my @arfiles = check_for_archive_files($ardir,$data,$ovar,\@date_list_obs);
      if (@arfiles == @date_list_obs) {
        foreach my $arfile (@arfiles) {
          if (! -e "$data/".tailname($arfile)) {
            push @dmget_files, $arfile if (! grep /$arfile/, @dmget_files);
          }
        }
        if (!$extracted_levels{$pvarlev}) {
          $extracted_levels{$pvarlev}{var} = $ovar;
          $extracted_levels{$pvarlev}{files} = [tailname(@arfiles)];
        }
      } else {
        print STDERR "WARNING: archive files not found for \"$data\" data, variable \"$ovar\"\n";
       #push @plots_skipped, $plot if (! grep /$plot/, @plots_skipped);
      }
    } else {
     #push @plot_list_out,  $plot    if (! grep /$plot/,    @plot_list_out);
     #push @plot_variables, $pvarlev if (! grep /$pvarlev/, @plot_variables);
    }
  }

  if (@dmget_files) {
    print "dmgetting ".scalar(@dmget_files)." $data (obs) data files\n";
    if ($Opt{VERBOSE} > 0) {
      foreach (@dmget_files) {
        print "  $_\n";
      }
    }
    system ("dmget @dmget_files") if !$TEST;
    print "copying files\n";
    print("gcp @dmget_files $workdir/$data");
    system ("gcp @dmget_files $workdir/$data") if !$TEST;
  }

  extract_level_from_3d_data( $data, \%extracted_levels );
}

exit if $Opt{SETUP_ONLY} || $TEST;

################################################
##########     PLOTTING SECTION      ###########
################################################

my $subdirectory = "PLOTS";
mkdir $subdirectory if (!-e $subdirectory);

# output directory
if ($Opt{out_dir}) {
   $Opt{out_dir} .= "/atmos_" . $Opt{yr1} . "_" . $Opt{yr2} . "/Wyman.eddies";
}

# compile fortran codes
if ($FAST) {
   if (!-e "eddies.so") {
      system("mkdir -p wrapit");
      print("gcp $package_location/wrapit/{average1d,average3d,filter}.f wrapit");
      system("gcp $package_location/wrapit/{average1d,average3d,filter}.f wrapit");
      system("WRAPIT -in -n eddies wrapit/{average1d,average3d,filter}.f");
   }
}

my $descriptor = $Opt{descriptor};
my $mlabel = "\'mlab=\"$descriptor\"\'";
my $years = "yrbeg=".$Opt{yr1}." yrend=".$Opt{yr2};

my $seasons;
my $labels = "ANN,DJF,MAM,JJA,SON";
# band-pass filtering options
my @bp_options;
@bp_options = qw/False/      if ($Opt{filter_opt} == 0);
@bp_options = qw/True/       if ($Opt{filter_opt} == 1);
@bp_options = qw/False True/ if ($Opt{filter_opt} == 2);

# if the timing file exist then remove it
unlink "timing.out" if (-e "timing.out");

foreach my $obs (@datasets) {
   mkdir "$subdirectory/$obs" if (!-e "$subdirectory/$obs");
   my $obsOpt = "\'odir=\"$obs\"\'";
   $obsOpt .= " \'olab=\"".uc($obs)."\"\'";
   foreach my $bp (@bp_options) {
      my @plotfiles;
      my $bandpass = "\'bandpass=$bp\'";
      foreach my $outvar (@plot_list_out) {
         my $plotvar = "\'var=\"$outvar\"\'";
         # check if all the files needed are present
         if (needed_files_present($outvar,"model") && needed_files_present($outvar,$obs)) {
            print  "ncl -Q $plotvar $years $mlabel $obsOpt $bandpass $seasons $ncl_script\n";
            system("ncl -Q $plotvar $years $mlabel $obsOpt $bandpass $seasons $ncl_script");

            system("$label_script -O -l $labels eddies.$outvar.ps");
            if ($bp eq "True") {
               rename "eddies.$outvar.ps","$subdirectory/$obs/$outvar.bandpass.ps";
               push @plotfiles, "$subdirectory/$obs/$outvar.bandpass.ps";
            } else {
               rename "eddies.$outvar.ps","$subdirectory/$obs/$outvar.nofilter.ps";
               push @plotfiles, "$subdirectory/$obs/$outvar.nofilter.ps";
            }
         } else {
            print STDERR "WARNING: skipping plots for $obs/$outvar ... missing files\n"
         }
      }
      if ($Opt{out_dir}) {
          my $odir = $Opt{out_dir};
          system("gzip -r @plotfiles");
          s/\.ps$/.ps.gz/ for @plotfiles;
          print  "gcp -cd @plotfiles gfdl:$odir/$obs/\n";
          system("gcp -cd @plotfiles gfdl:$odir/$obs/");
         #system("gzip -r $subdirectory/$obs/*.ps");
         #print  "gcp -cd -r $subdirectory/$obs gfdl:$odir/\n";
         #system("gcp -cd -r $subdirectory/$obs gfdl:$odir/");
      }
   }
}


########################################################################
###################   E N D   O F   S C R I P T   ######################
########################################################################

sub determine_data_chunks {
  my ($chunks,$yr1,$yr2) = @_;
  my ($databegyr,$dataendyr,$datachunk);
  my @yrs = split /,/, $chunks;
  if (@yrs == 1) {
    $databegyr = $yr1;
    $dataendyr = $yr2;
    $datachunk = $yrs[0];
  } elsif (@yrs == 3) {
    ($databegyr,$dataendyr,$datachunk) = @yrs;
  } else {
    print STDERR "ERROR: invalid entry for data_chunk\n";
    usage(); exit 1;
  }
  return ($databegyr,$dataendyr,$datachunk);
}

#-----------------------------------------------------------------------
# determine the component name (from experiment directory path)

sub get_component_name {
  my $expDir = shift;
  my $fileType = "ts"; # make this an arg?
  if ($expDir =~ /pp.*\/(.*?)\/$fileType\//) {
    return $1; 
  } else {
    die "ERROR: could not determine component name for directory: $expDir";
  }   
}

#-----------------------------------------------------------------------
# determine model name (used as label in plot)

sub get_exper_descriptor ($;$) {
  my ($expDir,$modelFile) = shift;
  my $desc;

  # first, try to read it from the model file (if present)
  if ($modelFile) {
    if (`ncdump -h $modelFile` =~ /:title = "(.*)" ;/) {
      $desc = $1;
      return $desc;
    }
  }

  # second, try to read it from the path name
  if ($expDir =~ /.*\/(.*?)\/(.*?)\/pp.*\//) {
    my $desc1 = $1; 
    my $desc2 = $2; 
    if ($desc2 !~ /^gfdl/) {
      $desc = $desc2;
    } else {  #if ($desc1 !~ /^gfdl/) {
      $desc = $desc1;
    }
    return $desc if $desc;
  }   

  die "ERROR: could not determine experiment descriptor for directory: $expDir";
} 

#-----------------------------------------------------------------------
# determine static file (usually for land mask or topog)

sub get_staticfile_name {
  my $expDir = shift;
  my $comp = shift;
  if ($expDir =~ /^(.*)\/.*\/.*\/.*$/) {
    my $staticfile = "$1/$comp.static.nc";
    die "ERROR: static file \'$staticfile\' could not be located" if (!-e $staticfile);
    return $staticfile;
  } else {
    die "ERROR: could not determine static file name for directory: $expDir";
  }
}

#-----------------------------------------------------------------------
# input: /aaa/bbb/ccc.xxx
# returns: ccc.xxx
# NOTE: input can also be an array reference (outputs array)

sub tailname {
   my @tail = @_;
   if (@tail == 1) {
     my $tail = $tail[0];
     while ($tail =~ s/^.*\///) {}
     return $tail;
   } else {
     my @tails = @tail;
     foreach (@tails) {
       while ($_ =~ s/^.*\///) {}
     }
     return @tails;
   }
}

#-----------------------------------------------------------------------
# archive file names are created
# the file names that exist will be returned

sub check_for_archive_files {
  my ($dir,$part1,$part2,$part3) = @_;
  my @dates;
  if (ref $part2 && !ref $part3) {
    @dates = @$part2;
  } elsif (ref $part3 && !ref $part2) {
    @dates = @$part3;
  } else {
    die "ERROR in check_for_archive_files: no (or mulitiple) array references found";
  }

  my @files = ();
  foreach my $date (@dates) {
    if (ref $part2) {
      push @files, "$dir/$part1.$date.$part3.nc" if (-e "$dir/$part1.$date.$part3.nc");
    } else {
      print "Looking for $dir/$part1.$part2.$date.nc\n" if ($TEST || $Opt{VERBOSE} > 1);
      push @files, "$dir/$part1.$part2.$date.nc" if (-e "$dir/$part1.$part2.$date.nc");
    }
  }
  return @files;
}

#-----------------------------------------------------------------------
# extract all single-level variables from 3d data

sub extract_level_from_3d_data {
  my ($dir,$href) = @_;

  foreach my $pvar (keys %{$href}) {   # pvar = plot var name = letterLevel, e.g., w500, u250
    print "extract: pvar = $pvar\n";
    my $lev = substr($pvar,1);
    my $var = $href->{$pvar}{var};
    my @files = @{$href->{$pvar}->{files}};
    my $odate = "$1-$2" if ( "$files[0]$files[$#files]" =~ /\.(\d{8})-\d{8}\..*\.\d{8}-(\d{8})\./ );
    die "ERROR in extract_level_from_3d_data: ooutput date not determined" if !$odate;
    my $tsfile = "$pvar.$odate.nc";

    # extract requested level from all files
    my @pfiles;
    foreach my $ifile (@files) {
      my $ofile = $ifile; $ofile =~ s/$var/$pvar/;
      my @options;
      push @options, "\'ifile=\"$dir/$ifile\"\'";
      push @options, "\'ofile=\"$dir/$ofile\"\'";
      push @options, "level=$lev";
      print  "ncl -Q @options $extract_script\n" if ($TEST || $Opt{VERBOSE} > 0);
      system("ncl -Q @options $extract_script")  if !$TEST;
      if ("$var$lev" ne $pvar) {
        print  "ncrename -h -v $var$lev,$pvar $dir/$ofile\n" if ($TEST || $Opt{VERBOSE} > 0);
        system("ncrename -h -v $var$lev,$pvar $dir/$ofile")  if !$TEST;
      }
      push @pfiles, "$dir/$ofile";
    }

    # concatenate files into single time series
    concatenate_files(@pfiles,"$dir/$tsfile");
  }
  return 1;
}

#-----------------------------------------------------------------------
# concatenates netcdf file into a single timeseries
# the input files are then removed

sub concatenate_files {
  my @files = @_;
  die "ERROR: at least 2 arguments required" if @files < 2;
  my $nn = $#files;
  my $n = $nn-1;

  # concatenate files into single time series
  if (@files > 2) {
    print  "ncrcat -h @files[0..$n] $files[$nn]\n" if ($TEST || $Opt{VERBOSE} > 0);
    system("ncrcat -h @files[0..$n] $files[$nn]")  if !$TEST;
    print  "rm -f @files[0..$n]\n" if ($TEST || $Opt{VERBOSE} > 0);
    system("rm -f @files[0..$n]")  if !$TEST;
  } elsif ($files[0] ne $files[1]) {
    print  "mv -f $files[0] $files[1]\n" if ($TEST || $Opt{VERBOSE} > 0);
    system("mv -f $files[0] $files[1]")  if !$TEST;
  }
  return 1;
}

#-----------------------------------------------------------------------

sub get_observed_data_chunk {
  my $dir = shift;

  # find first and last available years
  my @files = glob "$dir/*.*.nc"; 
  my ($yrbeg,$yrend,$chunk);

  if ($files[0] =~ /\.(\d\d\d\d)0101-(\d\d\d\d)1231\./) {
    $yrbeg = $1;
    $chunk = $2;
  }
  $yrend = $1 if ($files[$#files] =~ /\.\d\d\d\d0101-(\d\d\d\d)1231\./);
  die "ERROR: cound not determine yrbeg or yrend in directory \"$dir\"" if (!$yrbeg || !$yrend || !$chunk);
  $chunk = ($chunk - $yrbeg) + 1;

  return ($yrbeg,$yrend,$chunk);
}

#-----------------------------------------------------------------------
# yrbeg, yrend = starting/ending year of available data
# y1, y2       = starting/ending year requested
# yr1, yr2     = starting/ending year to be used

sub determine_observation_years {
  my($dir,$y1,$y2) = @_;
 #print "determine_observation_years: dir, y1, y2 = $dir, $y1, $y2\n";

  # find first and last available years
  my ($yrbeg,$yrend,$chunk) = get_observed_data_chunk($dir);

  # determine new years
  my $yr1 = int($y1);
  my $yr2 = int($y2);
  if ($yr1 < $yrbeg) {
    $yr2 = $yr2 + ($yrbeg-$yr1);
    $yr1 = $yrbeg;
  }
  if ($yr2 > $yrend) {
    my $dyr1 = $yr1-$yrbeg;
    my $dyr2 = $yr2-$yrend;
    $dyr1 = $dyr2 if ($dyr2 < $dyr1);
    $yr1 = $yr1 - $dyr1 if ($dyr1 > 0);
    $yr2 = $yrend;
  }

  printf "NOTE: Adjusted observed starting date from %4.4d to %4.4d\n", $y1, $yr1 if ($yr1 != $y1);
  printf "NOTE: Adjusted observed ending date from %4.4d to %4.4d\n",   $y2, $yr2 if ($yr2 != $y2);
  return ($yr1,$yr2);
}

#-----------------------------------------------------------------------

sub determine_obs_file_dates {
  my ($dir,$yr1,$yr2) = @_;
  my ($yrbeg,$yrend,$chunk) = get_observed_data_chunk($dir);
  die "ERROR in determine_obs_file_dates: years outside data years" if ($yr1 < $yrbeg || $yr2 > $yrend);
  my @dates;
  while ($yrbeg < $yrend) {
    my $ylast = $yrbeg + ($chunk-1);
    $ylast = $yrend if ($ylast > $yrend);
    if (($yrbeg < $yr1 && $ylast < $yr1) || ($yrbeg > $yr2 && $ylast > $yr2)) {
    # print "  yrbeg, ylast = $yrbeg, $ylast (OUT)\n";
    } else {
    # print "  yrbeg, ylast = $yrbeg, $ylast (IN)\n";
      push @dates, "$yrbeg"."0101-$ylast"."1231";
    }
    $yrbeg = $ylast+1;
  }
  return @dates;
}

#-----------------------------------------------------------------------
# given a plotting variable (e.g. vv250) and data directory
# this fuction checks that the needed model and observed files are present

sub needed_files_present {
   my ($plot,$dir) = @_;
   # extract the level
   my $lev = substr($plot,2);
   # extract first and second letters
   if ( substr($plot,0,1) && substr($plot,1,1)) {
      foreach my $var ( substr($plot,0,1).$lev, substr($plot,1,1).$lev ) {
         my @f = glob "$dir/$var.*.nc";
         print "glob \"$dir/$var.*.nc\" = @f\n" if $Opt{VERBOSE} > 0;
         return 0 if (!@f);
      }
      return 1;
   } else {
     die "ERROR in needed_files_present: problem extracting variable names from $plot";
   }
}

#################################################################################

sub usage {
   my $cmdname = substr($0,rindex($0,"/")+1);
   my $csv_list = join ",", @default_plots;
   print STDERR "

[1mOVERVIEW[0m
   Computes and plots near atmospheric eddy variance and covariances.
   The model climatology is compared to observation (NCEP Reanalysis).

[1mUSAGE[0m
   $cmdname -i [4mINDIR[0m -y [4mYR1[0m,[4mYR2[0m [OPTIONS...]

[1mREQUIRED OPTIONS[0m
   [1m-i[0m, [1m--in_data_dir[0m [4mPATH[0m
         Path for the input directory containing the time series of the daily pressure level fields.

   [1m-y[0m, [1m--year_range[0m [4mYR1[0m,[4mYR2[0m
         The year range for the analyzed data.
         This can be a subset of the years of avaliable data.

[1mOPTIONAL[0m
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
         If 4mDATABEGYR[0m,[4mDATAENDYR[0m,[4mDATACHUNK[0m are not given, then they are derived from
         [4mYR1[0m,[4mYR2[0m with the assumption that all the data will be analyzed and that there
         is one file. If the optional form of this is used [1m-c[0m [4mDATACHUNK[0m, then
         [4mDATABEGYR[0m and [4mDATAENDYR[0m will be set to [4mYR1[0m and [4mYR2[0m, respectively.

   [1m-f[0m, [1m--filter_opt[0m [4mLIST[0m
         Options for bandpass filtering: 0 = no filter; 1 = 2-7 day bandpass filter; 2 = both.
         Note, bandpass filtering can take a long time.

   [1m-p[0m, [1m--plots[0m [4mLIST[0m
         Comma-separated list of fields to be plotted.  The possible variables are:
         $csv_list.

";
  exit 1;
}

