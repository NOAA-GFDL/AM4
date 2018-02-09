#!/usr/bin/perl
use strict;
use Getopt::Std;

my $help = 0;
our ($opt_b,$opt_C,$opt_D,$opt_f,$opt_t,$opt_c);

unless (getopts 'bfDCc:t:') {
  $help = 1;
}

if (scalar(@ARGV) == 0) {
  $help = 1;
}

my $fname = $ARGV[0];

# optional arguments

if ($help) {
   usage("ncdate.pl");
   exit 1;
}

# netcdf must be loaded
my $platform = `gfdl_platform`; chomp $platform;
if ($platform eq "hpcs-csc") {
   my $error = 0;
   my @loaded_modules = split /:/, $ENV{'LOADEDMODULES'};
   foreach my $mod(qw/netcdf nco/) {
      my @mods = grep /$mod\/.*/, @loaded_modules;
      if (scalar @mods == 0) {
         print "ERROR: $mod module not loaded\n";
         $error++;
      }
   }
   if ($error) {
     die;
   }
} elsif ($platform eq "desktop") {
   # should have everything loaded
} else {
   die "Platform not supported.";
}

# initialize dates
my ($dump,$time) = dates_init($fname);

# debug info
if ($opt_D) {
   print "timename = $time->{'name'} \n";
   print "numtime = $time->{'size'} \n";
   print "units = $time->{'units'} \n";
   print "calendar = $time->{'calendar'} \n";
   print "bounds = $time->{'bounds'} \n";
   print "bounds axis = $time->{'axisbnds'} \n";
   print "baseyear = $time->{'year'} \n";
}

# first and last time index
my @index = ( 1, $time->{'size'} );
# substitute user specified indices
if ($opt_t) {
   my @newIndex = split /,/,$opt_t;
   if (scalar @newIndex <= 2) {
      for (my $i=0; $i<@newIndex; $i++) {
         if ($newIndex[$i] >= $index[0] && $newIndex[$i] <= $index[1]) {
            $index[$i] = $newIndex[$i];
         }
      }
   }
}

# extract either time bounds or time values
my @T;
if ($opt_b) {
   $T[@T] = `ncks -H -F -C -s "%g\n" -v $time->{'bounds'} -d $time->{'name'},$index[0] -d $time->{'axisbnds'},1 $fname`;
   $T[@T] = `ncks -H -F -C -s "%g\n" -v $time->{'bounds'} -d $time->{'name'},$index[1] -d $time->{'axisbnds'},2 $fname`;
} else {
   $T[@T] = `ncks -H -F -C -s "%g\n" -v $time->{'name'} -d $time->{'name'},$index[0] $fname`;
   $T[@T] = `ncks -H -F -C -s "%g\n" -v $time->{'name'} -d $time->{'name'},$index[1] $fname`;
}
foreach my $t (@T) {
   my @date = get_current_date($t,$time);
   print  "@date \n";
}

# check for contiguous time values
# get all time axis values
if ($opt_C) {
   # check that time axis is increasing
   my @T = split /\n/, `ncks -H -F -C -s "%g\n" -v $time->{'name'} -d $time->{'name'},1,$time->{'size'} $fname`;
   for (my $i=1; $i<3; $i++) {
      if ($T[$i] <= $T[$i-1]) {
         die "Time axis not increasing.";
      }
   }
   # check time bounds
   if ($time->{'bounds'} && $time->{'axisbnds'}) {
      my @TB = split /\n/, `ncks -H -F -C -s "%g\n" -v $time->{'bounds'} -d $time->{'name'},1,$time->{'size'} -d $time->{'axisbnds'},1 $fname`;
      my @TE = split /\n/, `ncks -H -F -C -s "%g\n" -v $time->{'bounds'} -d $time->{'name'},1,$time->{'size'} -d $time->{'axisbnds'},2 $fname`;
      for (my $i=1; $i<@TB; $i++) {
         if ($TB[$i] != $TE[$i-1]) {
            die "Time bounds are not continuous.";
         }
      }
      print "Time bounds are continuous.\n";
   } else {
      print "Time bounds not checked.\n";
   }
}

if ($opt_D) {
   my @times = split /\n/, `ncks -H -F -s "%g\n" -v $time->{'name'} -d $time->{'name'},1,$time->{'size'} $fname`;
   foreach my $t (@times) {
      my @date = get_current_date($t,$time);
      printf "%10f %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d\n", $t, @date;
   }
}


###############################################################

sub dates_init {
   my $file = shift;
   my %time;
   my $dump = `ncdump -h $file`;
   # time axis info
   if ($dump =~ /\t(.+) = UNLIMITED ; \/\/ \((\d+) currently\)/) {
      my $name = $1;
      $time{'size'} = $2;
      $time{'name'} = $name;
      # time units
      if ($dump =~ /\t\t$name:units = "(.*)"/) {
         my $units = lc($1);
         # calendar
         if ($dump =~ /\t\t$name:calendar = "(.*)"/) {
            $time{'calendar'} = lc($1);
         } else {
            $time{'calendar'} = "";
            $time{'calendar'} = $opt_c if ($opt_c);
         }
         # axis bounds
         if ($dump =~ /\t\t$name:bounds = "(.*)"/) {
            $time{'bounds'} = lc($1);
            $dump =~ /\t\w+ $time{'bounds'}\(.+, (.+)\)/;  $time{'axisbnds'} = lc($1);
         } else {
            $time{'bounds'} = "";
            $time{'axisbnds'} = "";
         }
         ($time{'units'},$time{'year'},$time{'month'},$time{'day'},$time{'hour'},$time{'minute'},$time{'second'}) = decode_units($units);
         return ($dump,\%time);
      } else {
         die "Could not determine units for time axis.";
      }
   } else {
      die "Could not determine if there was a time axis";
   }
}

sub decode_units {
   my $units = shift;
   my @W = split /\s+/, $units; # units since yyyy-mm-dd hh:mm:ss
   if ($W[1] ne "since" || scalar(@W) < 3) {
      print STDERR "$W[0] $W[1] $W[2] ... \n";
      die "Could not decode units";
   }
   my @D = split /-/, $W[2];
   if (scalar @D != 3) {
      print STDERR "$W[2]\n";
      die "Could not decode base date";
   }
   my @T = (0,0,0);
   if (scalar @W  > 3) {
      @T = split /:/, $W[3];
      if (scalar @T != 3) {
         print STDERR "$W[3]\n";
         die "Could not decode base time";
      }
   }
   return ($W[0],$D[0],$D[1],$D[2],$T[0],$T[1],$T[2]);
}

sub get_current_date {
   my $deltaTime = shift;
   my $timeRef = shift;
   my @basedate = ($timeRef->{'year'}, $timeRef->{'month'},  $timeRef->{'day'},
                   $timeRef->{'hour'}, $timeRef->{'minute'}, $timeRef->{'second'});
   return get_date($deltaTime,\@basedate,$timeRef);
}

sub get_date {
   my $deltaTime = shift;
   my $dateRef = shift;
   my $timeRef = shift;
   my @date = @$dateRef;

   my $idate = -1; 
   if ($timeRef->{'units'} eq 'seconds') {
      $idate = 5;
   } elsif ($timeRef->{'units'} eq 'minutes') {
      $deltaTime = $deltaTime * 60.;
      $idate = 5;
   } elsif ($timeRef->{'units'} eq 'hours') {
      $deltaTime = $deltaTime * 3600.;
      $idate = 5;
   } elsif ($timeRef->{'units'} eq 'days') {
      $deltaTime = $deltaTime * 86400.;
      $idate = 5;
   } elsif ($timeRef->{'units'} eq 'months') {
      $idate = 1;
   } elsif ($timeRef->{'units'} eq 'years') {
      $idate = 0;
   } else {
      print "units = ".$timeRef->{'units'}."\n";
      die "Could not decode units";
   }   

   while ($idate >= 0) {
      $date[$idate] = $date[$idate] + $deltaTime;
      if ($idate == 0) { last; }
      # truncate seconds->minutes or minutes->hours
      if ($idate > 3) {
         $deltaTime = int($date[$idate]/60);
         $date[$idate] = $date[$idate] - ($deltaTime*60);
      }
      # truncate hours->days
      if ($idate == 3) {
         $deltaTime = int($date[$idate]/24);
         $date[$idate] = $date[$idate] - ($deltaTime*24);
      }
      # truncate days->months
      if ($idate == 2) {
         my ($moreMonths,$lessDays) = days2months(\@date,$timeRef->{'calendar'});
         $deltaTime = $moreMonths;
         $date[$idate] = $date[$idate] - $lessDays;
      }
      # truncate months-->years
      if ($idate == 1) {
         $deltaTime = int(($date[$idate]-1)/12);
         $date[$idate] = $date[$idate] - ($deltaTime*12);
      }
      $idate--;
   }
   return @date;
}

sub days2months {
   my $dateRef = shift;
   my $cal = shift;
   my $year = $dateRef->[0];
   my $month = $dateRef->[1];
   my $day = $dateRef->[2];
   my $numMonths = 0;
   # determine number of days in month for this calendar
   my $nday = number_of_days_in_month($year,$month,$cal);
   while ($day > $nday) {
#if ($DEBUG) {
#   print "year=$year ";
#   print "month=$month ";
#   print "day=$day ";
#   print "nday=$nday \n";
#}
      $day = $day - $nday;
      $numMonths++;
      $month++;
      if ($month == 13) {
         $month = 1;
         $year++;
      }
      $nday = number_of_days_in_month($year,$month,$cal);
   }
   return ($numMonths,$dateRef->[2]-$day);
}

# return the number of days for the given month, year and calendar
# args: year,month,calendar; returns: integer days in month
sub number_of_days_in_month {
   my $yr = shift;
   my $mo = shift;
   my $cal = shift;
   my @mdays = (31,28,31,30,31,30,31,31,30,31,30,31);
   if ($cal eq "constant") {
       return 30;
   }

   my $ndays = $mdays[$mo-1];
   if ($cal eq "noleap") {
      return $ndays;
   } elsif ($cal eq "julian" || $cal eq "gregorian") {
      if ($mo == 2) {
         $ndays = $ndays + yleap($yr,$cal);
      }
      return $ndays;
   } else {
      print "cal: $cal\n";
      die "ERROR: unable to determine calendar";
   }
}


sub yleap {
   my $yr = shift;
   my $cal = shift;
   if ($cal eq "noleap" || $cal eq "constant") {
      return 0;
   }
   my $yleap = 0;
   if ($yr%4 == 0) {
      $yleap = 1;
   }
   if ($cal eq "julian") {
      return $yleap;
   }
   if ($yr%100 == 0) {
      $yleap = 0;
   }
   if ($yr%400 == 0) {
      $yleap = 1;
   }
   return $yleap;
}

sub usage {
   my $command = shift;
   print "
USAGE
   $command [-bfCD] file
\nOPTIONS
   -b  print bounds
   -C  check time axis continuity (time bounds must be present)
   -f  format date output
   -D  debug
\nEXAMPLE
   $command -C myfile.nc \n\n";
   return 1;
}
