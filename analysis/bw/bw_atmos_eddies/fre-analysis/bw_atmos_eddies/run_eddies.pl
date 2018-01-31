#!/usr/bin/perl
#------------------------------------
# Computes and plots eddy variances and covariances for
# both model and observation (NCEP Reanalysis) for standard
# atmospheric fields on several pressure levels.
#------------------------------------
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
Getopt::Long::Configure("bundling");

my $fast_version = 1;

my @required_modules = qw/fre ncarg/;
push @required_modules, "ifort" if $fast_version;

# make sure the required modules are loaded from parent shell
my $platform = `gfdl_platform`; chomp $platform;
my $message = "";
if ($platform eq "hpcs-csc") {
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
   die "Invalid platform";
}

#print "$0 $*\n";
my $TEST = 0;

my @all_plots = qw/ qq850 tt200 tt500 tt850 uu200 uu850 uv200 uv850 vq850 vt200 vt850 vv200 vv850 wq850 wt500 wt850 ww500 ww850 zz200 zz500 zz850 /;
my @datasets = qw/ ncep2 merra /;
my @all_variable_list = qw/ u200 u850 v200 v850 omg500 omg850 t200 t500 t850 q850 h200 h500 h850 /;

# argument parsing
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
   print "ERROR: no argument given for input directory\n";
   $Opt{HELP} = 1;
}

#if (!$Opt{descriptor} && !$Opt{HELP}) {
#   print "ERROR: no argument given for descriptor\n";
#   $Opt{HELP} = 1;
#}

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

#--------------------------
# usage message then exit
#--------------------------

if ($Opt{HELP}) {
   usage("run_eddies.pl");
   exit 1;
}

#------------------------------------------------------------
# determine data chunks (the years in each time series file)
#------------------------------------------------------------

if ($Opt{data_chunks}) {
   my @yrs = split /,/, $Opt{data_chunks};
   if (@yrs == 1) {
      $Opt{databegyr} = $Opt{yr1};
      $Opt{dataendyr} = $Opt{yr2};
      $Opt{datachunk} = $yrs[0];
   } elsif (@yrs == 3) {
      ($Opt{databegyr},$Opt{dataendyr},$Opt{datachunk}) = @yrs;
   } else {
      print "ERROR: invalid entry for data_chunk\n";
      $Opt{HELP} = 1;
   }
} else {
   $Opt{databegyr} = $Opt{yr1};
   $Opt{dataendyr} = $Opt{yr2};
   $Opt{datachunk} = ($Opt{dataendyr} - $Opt{databegyr}) + 1
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
my $obs_data_location = $ENV{FRE_ANALYSIS_ARCHIVE}."/bw/bw_atmos_eddies";

# external scripts
my $date_script = "$package_location/shared/bw/file_dates.pl";
my $fvalue_script = "$package_location/shared/bw/ncfillvalue.pl";
my $label_script = "$package_location/shared/bw/pagelabels.pl";
my $ncl_script = "$package_location/ncl/eddies.ncl";
my $extract_script = "$package_location/ncl/extract.ncl";

#-------------------------------------
# list of date ranges for each file
#-------------------------------------

my @date_list = split /\n/, `$date_script daily $Opt{databegyr} $Opt{dataendyr} $Opt{datachunk}`;
my $full_date_range = $date_list[0] . $date_list[$#date_list]; $full_date_range =~ s/-\d{16}//;
if ($Opt{VERBOSE} > 0) {
  print "Dates:\n";
  foreach (@date_list) {
    print "  $_\n";
  }
}

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

#------------------------
# which data sets to use
#------------------------

@datasets = split /,/, $Opt{datasets} if $Opt{datasets};
print "datasets: @datasets\n";

#-------------------------
# list of requested plots
#-------------------------
my @plot_list = @all_plots;
@plot_list = split /,/, $Opt{plots} if $Opt{plots};
my @plot_list_out;

# translate 'z' to 'h' in all plot requests
foreach (@plot_list) { s/z/h/g }

#-----------------------------------------------------------
# mapping between plot  variable letters and variable names
#-----------------------------------------------------------
# between plot name -> gfdl variable name
# also used for observation variable name
my %gfdl_letter_map = ( q=>"q", t=>"t", u=>"u", v=>"v", w=>"omg", h=>"h" );

# between plot name -> 3D cmip variable name
my %cmip_letter_map = ( q=>"hus", t=>"ta", u=>"ua", v=>"va", w=>"wap", h=>"zg" );
my $cmip_suffix = "_unmsk";

#---------------------------------------------
# create a list of GFDL variable names needed
#---------------------------------------------
my @dmget_variables;
my @obs_variables;
my @plots_skipped;
my %extract_lev;

foreach my $plot (@plot_list) {
  my $lev = substr($plot,2,3);
  foreach my $pvar ( substr($plot,0,1), substr($plot,1,1) ) {
    my $ovar = $gfdl_letter_map{$pvar}.$lev;
    my $cplot = $plot; $cplot =~ s/200/250/; # cmip 250hPa only

    # FIRST: check for cmip names
    if ($cmip_letter_map{$pvar}) {
      my $cvar = $cmip_letter_map{$pvar}.$cmip_suffix;
      my $xvar = $pvar.$lev;
      if (! -e "model/$component_name.$full_date_range.$cvar.nc") {
        print "model/$component_name.$full_date_range.$cvar.nc does not exist\n" if $Opt{VERBOSE} > 1;
        if (check_for_archive_files($Opt{in_data_dir},$component_name,\@date_list,$cvar) == @date_list) {
          push @dmget_variables, $cvar  if (! grep /$cvar/, @dmget_variables);
          push @plot_list_out, $cplot if (! grep /$cplot/, @plot_list_out);
          $extract_lev{$xvar} = $cvar if !$extract_lev{$xvar};
          $ovar =~ s/200/250/; # cmip 250hPa only
          push @obs_variables, $ovar if (! grep /$ovar/, @obs_variables);
          next;
        } else {
         #push @plots_skipped, $plot if (! grep /$plot/, @plots_skipped);
        }
      } else {
        push @dmget_variables, $cvar  if (! grep /$cvar/, @dmget_variables);
        push @plot_list_out, $cplot if (! grep /$cplot/, @plot_list_out);
        $extract_lev{$xvar} = $cvar if !$extract_lev{$xvar};
        $ovar =~ s/200/250/; # cmip 250hPa only
        push @obs_variables, $ovar if (! grep /$ovar/, @obs_variables);
        next;
      }
    }

    # SECOND: check for for gfdl names
    if ($gfdl_letter_map{$pvar}) {
      my $gvar = $gfdl_letter_map{$pvar}."$lev";
      if (! -e "model/$component_name.$full_date_range.$gvar.nc") {
        print "model/$component_name.$full_date_range.$gvar.nc does not exist\n" if $Opt{VERBOSE} > 1;
        if (check_for_archive_files($Opt{in_data_dir},$component_name,\@date_list,$gvar) == @date_list) {
          push @dmget_variables, $gvar  if (! grep /$gvar/, @dmget_variables);
          push @obs_variables, $ovar if (! grep /$ovar/, @obs_variables);
          push @plot_list_out, $plot   if (! grep /$plot/, @plot_list_out);
          next;
        } else {
          print "WARNING: no archive files found for \"$pvar\" ... skipping plot $plot\n";
          push @plots_skipped, $plot if (! grep /$plot/, @plots_skipped);
        }
      } else {
        push @dmget_variables, $gvar  if (! grep /$gvar/, @dmget_variables);
        push @obs_variables, $ovar if (! grep /$ovar/, @obs_variables);
        push @plot_list_out, $plot   if (! grep /$plot/, @plot_list_out);
      }
    }

    if (!$cmip_letter_map{$pvar} && !$gfdl_letter_map{$pvar}) {
      die "ERROR: Invalid variable letter \"$pvar\" in plot field, $plot";
    }
  }
}

# add pmask2 and/or ps to dmget list
my $pres_vars_found = 0;
foreach my $var (qw/pmaskv2 ps/) {
  if (! -e "model/$var.$full_date_range.nc") {
    if (check_for_archive_files($Opt{in_data_dir},$component_name,\@date_list,$var) == @date_list) {
      push @dmget_variables, $var if (! grep /$var/, @dmget_variables);
      $pres_vars_found++;
    }
  } else {
    $pres_vars_found++;
  }
}

print "variables needed: @dmget_variables\n";

#-------------------------------
# local directory locations
#-------------------------------

my $workdir = $ENV{TMPDIR} . "/bw_atmos_eddies";
mkdir $workdir if (!-e $workdir);
chdir $workdir;
mkdir "model" if (!-e "model");

#-------------------------------
# dmget and copy model files
#-------------------------------
if (!$Opt{SKIP_SETUP}) {

my $FLAG_ps_from_zs = 0;

my @dmget_file_list;
foreach my $var (@dmget_variables) {
   if (! -e "model/$var.$full_date_range.nc") {
      push @dmget_file_list, check_for_archive_files($Opt{in_data_dir},$component_name,\@date_list,$var);
   }
}

# create an appropriate ps mask using zsurf (add staticfile to dmget list)
if ($pres_vars_found <= 0) {
  if (!-e "model/ps_est.nc") {
    if ($Opt{staticfile}) {
      push @dmget_file_list, $Opt{staticfile};
      $FLAG_ps_from_zs = 1;
    } else {
      die "Missing model ps files and no static file defined";
    }
  }
  print "WARNING: estimating surface pressure for masking using topography\n";
}

print "dmget variables: @dmget_variables\n";

# remove skipped plots from the plot list
if (@plots_skipped) {
   print "plots skipped: @plots_skipped\n";
   foreach my $plot (@plots_skipped) {
      my $index = 0;
      $index++ until $plot_list[$index] eq $plot; 
      splice(@plot_list, $index, 1);
   }
   print "plot list @plot_list\n";
}

# dmget and copy the file to working directory
if (@dmget_file_list) {
   print "dmgetting ".scalar(@dmget_file_list)." model files\n";
   if ($Opt{VERBOSE} > 0) {
     foreach (@dmget_file_list) {
       print "$_\n";
     }
   }
   system ("dmget @dmget_file_list") if !$TEST;
   print "copying files\n";
   system ("gcp @dmget_file_list $workdir") if !$TEST;
}


#----------------------------------
# concatenate model files together
#----------------------------------

foreach my $var (@dmget_variables) {
   my @commands;
   if (! -e "model/$var.$full_date_range.nc") {
      my @files;
      foreach my $date (@date_list) {
         push @files, "$component_name.$date.$var.nc";
      }
      my $tsfile = "$component_name.$full_date_range.$var.nc";
      if (@files > 1 || $files[0] ne $tsfile) {
        print "ncrcat -h @files $tsfile\n" if ($TEST || $Opt{VERBOSE} > 0);
        system ("ncrcat -h @files $tsfile") if !$TEST;
        unlink @files if !$TEST;
      }
      # rename file and move to model directory
      print "rename $tsfile,\"model/$var.$full_date_range.nc\"\n" if ($TEST || $Opt{VERBOSE} > 0);
      rename $tsfile,"model/$var.$full_date_range.nc" if !$TEST;
   }
}

# create approx psurf file (when pmaskv2 and ps are not found)
# use millibars since ncl code uses variable name for level
if ($FLAG_ps_from_zs) {
   if (!-e "model/ps_est.nc") {
     my $cmd = "ps=1000*((288-zsurf*5e-3)/288)^6.852";
     print   "ncap2 -v -s \"$cmd\" ".tailname($Opt{staticfile})." model/ps_est.nc\n" if ($TEST || $Opt{VERBOSE} > 0);
     system ("ncap2 -v -s \"$cmd\" ".tailname($Opt{staticfile})." model/ps_est.nc")  if !$TEST;
     unlink tailname($Opt{staticfile}) if !$TEST;
   }
}

#<<<<< EXTRACT LEVELS FROM CMIP 3D FILES >>>>>>>>
foreach my $pvarlev (keys %extract_lev) {
  my $cvar = $extract_lev{$pvarlev};
  my $lev  = substr($pvarlev,1,3);
  $lev = "250" if ($lev eq "200"); # 200hPa level should be replaced with 250hPa
  my $gvar = substr($pvarlev,0,1).$lev;
  my $ifile = "model/$cvar.$full_date_range.nc";
  my $ofile = "model/$gvar.$full_date_range.nc";
  if (! -e $ofile) {
    my @options;
    push @options, "\'ifile=\"$ifile\"\'";
    push @options, "\'ofile=\"$ofile\"\'";
    push @options, "level=$lev";
    print  "ncl @options $extract_script\n" if ($TEST || $Opt{VERBOSE} > 0);
    system("ncl @options $extract_script")  if !$TEST;
    print  "ncrename -h -v $cvar$lev,$gvar $ofile\n" if ($TEST || $Opt{VERBOSE} > 0);
    system("ncrename -h -v $cvar$lev,$gvar $ofile")  if !$TEST;
  }
}

#-------------------------------
# copy observed files
#-------------------------------

foreach my $data (@datasets) {
   mkdir $data   if (!-e $data);
   my $observed_directory = "$obs_data_location/$data";
   my %observed_dates = ( ncep =>"19710101-20141231",
                          ncep2=>"19790101-20151231", 
                          merra=>"19790101-20131231" );
   my @dmfiles;
   my @commands;
   foreach my $ovar (@obs_variables) {
     #my $ovar = $var; $ovar =~ s/w/omg/;
      my $pvar = $ovar; $pvar =~ s/omg/w/;
      my $ofile = "$ovar.".$observed_dates{$data}.".nc";
      my $dmfile = "$observed_directory/$ofile";
      if (-e $dmfile) {
        my $pfile = "$data/$pvar.".$observed_dates{$data}.".nc";
        if (!-e $pfile) {
          push @dmfiles, $dmfile;
          if ($ovar ne $pvar) {
            push @commands, "mv $data/$ofile $pfile";
            push @commands, "ncrename -h -v $ovar,$pvar $pfile";
          }
        }
      } else {
        print "WARNING: $dmfile does not exist\n";
      }
   }
   # execute the commands
   if (@dmfiles) {
      print "dmgetting ".scalar(@dmfiles)." observed $data data files\n";
      system("dmget @dmfiles");
      print "copying $data data files\n";
      system("gcp @dmfiles $data/");
      foreach my $cmd (@commands) {
         print "$cmd\n" if $Opt{VERBOSE} > 0 || $TEST;
         system ($cmd) if !$TEST;
      }
   }
}

} # SKIP_SETUP

exit if $Opt{SETUP_ONLY} || $TEST;
#-------------------------------
# plotting section
#-------------------------------

my $subdirectory = "PLOTS";
mkdir $subdirectory if (!-e $subdirectory);

# output directory
if ($Opt{out_dir}) {
   $Opt{out_dir} .= "/atmos_" . $Opt{yr1} . "_" . $Opt{yr2} . "/Wyman.eddies";
}

# compile fortran codes
if ($fast_version) {
   if (!-e "eddies.so") {
      system("mkdir -p wrapit");
      system("gcp $package_location/wrapit/{average1d,average3d,filter}.f wrapit");
      system("WRAPIT -in -n eddies wrapit/{average1d,average3d,filter}.f");
   }
}

my $descriptor = $Opt{descriptor};
my $mlabel = "\'mlab=\"$descriptor\"\'";
my $years = "yrbeg=".$Opt{yr1}." yrend=".$Opt{yr2};

my $seasons;
my $labels = "ANN,DJF,MAM,JJA,SON";
#$seasons = "\'seasons=\"1\"\'";
#$labels = "DJF";


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
         if (needed_files_present($outvar,$obs)) {
            print  "ncl $plotvar $years $mlabel $obsOpt $bandpass $seasons $ncl_script\n";
            system("ncl $plotvar $years $mlabel $obsOpt $bandpass $seasons $ncl_script");

            system("$label_script -O -l $labels eddies.$outvar.ps");
            if ($bp eq "True") {
               rename "eddies.$outvar.ps","$subdirectory/$obs/$outvar.bandpass.ps";
               push @plotfiles, "$subdirectory/$obs/$outvar.bandpass.ps";
            } else {
               rename "eddies.$outvar.ps","$subdirectory/$obs/$outvar.nofilter.ps";
               push @plotfiles, "$subdirectory/$obs/$outvar.nofilter.ps";
            }
         } else {
            print "WARNING: skipping plots for $obs/$outvar ... missing files\n"
         }
      }
      if ($Opt{out_dir}) {
          my $odir = $Opt{out_dir};
          system("gzip -r @plotfiles");
          s/\.ps$/.ps.gz/ for @plotfiles;
          print "gcp -cd @plotfiles gfdl:$odir/$obs/\n";
          system("gcp -cd @plotfiles gfdl:$odir/$obs/");
         #system("gzip -r $subdirectory/$obs/*.ps");
         #print "gcp -cd -r $subdirectory/$obs gfdl:$odir/\n";
         #system("gcp -cd -r $subdirectory/$obs gfdl:$odir/");
      }
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

#----------------------------------------------
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

#----------------------------------------------
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

#----------------------------------------------
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

#----------------------------------------------

sub check_for_archive_files {
  my ($dir,$comp,$dateref,$var) = @_;
  my @files = ();
  foreach my $date (@$dateref) {
    push @files, "$dir/$comp.$date.$var.nc" if (-e "$dir/$comp.$date.$var.nc");
  }
  if (scalar(@files) == scalar(@$dateref)) {
    return @files;
  } else {
    return ();
  }
}

#----------------------------------------------
# given a variable name and the plot list
# this function returns the plots associated with the given variable

# sub plots_for_this_variable {
#    my $gfdl_var = shift;
#    my $plot_list_ref = shift;
#    my $current_list_ref = shift;
#    my ($var,$lev) = split_var_name($gfdl_var);
#    
#    my @skip;
#    foreach my $plot (@$plot_list_ref) {
#       # only return if not in the current list
#       next if (grep /$plot/, @$current_list_ref);
#       if ($var eq $letter_map{substr($plot,0,1)} || $var eq $letter_map{substr($plot,1,1)}) {
#          push @skip, $plot if ($lev == substr($plot,2,3));  # match
#       }
#    }
#    return @skip;
# }

#----------------------------------------------

sub split_var_name {
   my $var = shift;
   if ($var =~ /([a-z]+)(\d+)/) {
      return ($1,$2);
   } else {
      die "Can not parse variable name: $var";
   }
}

#----------------------------------------------
# given a plotting variable (e.g. vv200) and the observation directory
# this fuction checks that the needed model and observed files are present

sub needed_files_present {
   my ($plot,$obs) = @_;
   # extract the level
   my $lev = substr($plot,2,3);
   # extract first and second letters
   if ( substr($plot,0,1) && substr($plot,1,1)) {
      foreach my $var ( substr($plot,0,1).$lev, substr($plot,1,1).$lev ) {
         my @a = glob "model/$var.*.nc";
         my @b = glob "$obs/$var.*.nc";
         print "glob \"model/$var.*.nc\" =  @a\n" if $Opt{VERBOSE} > 0;
         print "glob \"$obs/$var.*.nc\" =  @b\n"  if $Opt{VERBOSE} > 0;
         return 0 if (!@a);
         return 0 if (!@b);
      }
      return 1;
   } else {
     die "ERROR in needed_files_present: problem extracting variable names from $plot";
   }
}

#################################################################################

sub usage {
   my $cmdname = shift;
   my $cvs_list = join ",", @all_plots;
   print "

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
         If [4mDATABEGYR[0m,[4mDATAENDYR[0m,[4mDATACHUNK[0m are not given, then they are derived from
         [4mYR1[0m,[4mYR2[0m with the assumption that all the data will be analyzed and that there
         is one file. If the optional form of this is used [1m-c[0m [4mDATACHUNK[0m, then
         [4mDATABEGYR[0m and [4mDATAENDYR[0m will be set to [4mYR1[0m and [4mYR2[0m, respectively.

   [1m-f[0m, [1m--filter_opt[0m [4mLIST[0m
         Options for bandpass filtering: 0 = no filter; 1 = 2-7 day bandpass filter; 2 = both.
         Note, bandpass filtering can take a long time.

   [1m-p[0m, [1m--plots[0m [4mLIST[0m
         Comma-separated list of fields to be plotted.  The possible variables are:
         $cvs_list.

   By default, model data is compared versus NCEP Reanalysis I data.

[1mEXAMPLE[0m
   $cmdname -i /archive/user/myExperName/pp/atmos_daily_plev/ts/daily/1yr \
            -o /nbhome/\$USER/figures/myExperName/analysis \
            -d myExperName -y 1981,2000 -c 1981,2000,10

";
   return 0;
}


