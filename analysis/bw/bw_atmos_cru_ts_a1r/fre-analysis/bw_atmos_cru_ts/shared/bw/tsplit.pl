#!/usr/bin/perl
use strict;
use Getopt::Std;

my $help = 0;
our ($opt_c);

unless (getopts 'c:') {
   $help = 1;
}

usage() if (scalar @ARGV == 0 || $help);
my $file = $ARGV[0];
#my $period = "month";
#$period = $ARGV[1] if (scalar @ARGV > 1);
#usage() if ($period ne "month" && $period ne "year");

# time/calendar initialization
my $dump = `ncdump -h $file`;
my @dump = split /\n/, $dump;
my %timeInfo = calendar();
my %dateBase = decode_units($timeInfo{'units'});

# time axis values
my $times = `ncks -H -F -s "%f\n" -v $timeInfo{'name'} -d $timeInfo{'name'},1,$timeInfo{'size'} $file`;
my @timeValues = split /\n/, $times;

# initial date
my @basedate = ($dateBase{'year'},  $dateBase{'month'}, $dateBase{'day'},
                $dateBase{'hour'},  $dateBase{'minute'},$dateBase{'second'});
my $DEBUG = 0;
if ($DEBUG) {
   my $time0 = 0;
   printf "%10f %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d\n", $time0, @basedate;
}

my $time0 = $timeValues[0];
my @date = get_current_date(\@basedate,$timeInfo{'calendar'},$dateBase{'units'},$timeValues[0]);
if ($DEBUG) {
   printf "%10f %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d\n", $timeValues[0], @date;
}
my $month = $date[1];
my @countList;
my $count = 1;
my $countListSum = 0;
for (my $i=1; $i<scalar(@timeValues); $i++) {
   @date = get_current_date(\@date,$timeInfo{'calendar'},$dateBase{'units'},$timeValues[$i]-$timeValues[$i-1]);
   if ($date[1] != $month) {
      if ($date[2] == 1 && $date[3] == 0 && $date[4] == 0 &&$date[5] == 0) {
         $count++;
         if ($DEBUG) {
            printf "%10f %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d\n", $timeValues[$i], @date;
            print "------END-OF-MONTH----- ($count days)\n";
         }
         $countList[@countList] = $count;
         $countListSum = $countListSum + $count;
         $count = 0;
         $month = $date[1];
      } else {
         $countList[@countList] = $count;
         $countListSum = $countListSum + $count;
         if ($DEBUG) {
            print "------END-OF-MONTH----- ($count days)\n";
            printf "%10f %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d\n", $timeValues[$i], @date;
         }
         $count = 1;
         $month = $date[1];
      }
   } else {
         $count++;
         if ($DEBUG) {
            printf "%10f %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d\n", $timeValues[$i], @date;
         }
   }
}

# cleanup
if ($count > 0) {
   if ($DEBUG) {
      print "------END-OF-MONTH----- ($count days)\n";
   }
   $countList[@countList] = $count;
   $countListSum = $countListSum + $count;
}
   
# output
if ($countListSum != $timeInfo{'size'}) {
    print "time length: $timeInfo{'size'} ne count: $countListSum\n";
    die "ERROR: unable to split into months";
} else {
    foreach (@countList) {
       print "$_\n";
    }
}



# base date

sub calendar {
   my ($timeName,$timeSize,$timeUnits,$timeCal);
   foreach (@dump) {
      # time axis name and size
      if (/UNLIMITED/) {
         s/\t//g;
         my @W = split /\s+/;
         $timeName = $W[0];
         $timeSize = substr($W[5],1);
      }
      # locate time variable attributes
      if (/^\t\t$timeName:units = /) {
         /"(.*)"/;
         $timeUnits = lc($1);
      }
      if (/^\t\t$timeName:calendar = /) {
         /"(.*)"/;
         $timeCal = lc($1);
      }
   }
   if (!$timeCal) {
      if ($opt_c) {
         $timeCal = $opt_c;
      } else {
         print STDERR "Warning: setting calendar to julian\n";
         $timeCal = "julian";
      }
   }
   return ('name'=> $timeName, 'units' => $timeUnits, 'size' => $timeSize, 'calendar' => $timeCal);
}

sub decode_units {
   my $units = shift;
   my @W = split /\s+/, $units;
   if ($W[1] ne "since") {
      print "$W[0] $W[1] $W[2] ... \n";
      die "Could not decode units";
   }
   my @D = split /-/, $W[2];
   my @T = split /:/, $W[3];
   return ( 'units'=>$W[0], 'year'=>$D[0], 'month'=>$D[1],  'day'=>$D[2],
                            'hour'=>$T[0], 'minute'=>$T[1], 'second'=>$T[2] );
}

sub get_current_date {
   my $refDate = shift;
   my $calendar = shift;
   my $units = shift;
   my $deltaTime = shift;
   my @date = @$refDate;

   my $idate = -1;
   if ($units eq 'seconds') {
      $idate = 5;
   } elsif ($units eq 'minutes') {
      $deltaTime = $deltaTime * 60.;
      $idate = 5;
   } elsif ($units eq 'hours') {
      $deltaTime = $deltaTime * 3600.;
      $idate = 5;
   } elsif ($units eq 'days') {
      $deltaTime = $deltaTime * 86400.;
      $idate = 5;
   } elsif ($units eq 'months') {
      $idate = 1;
   } elsif ($units eq 'years') {
      $idate = 0;
   } else {
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
         my ($moreMonths,$lessDays) = days2months(\@date,$calendar);
         $deltaTime = $moreMonths;
         $date[$idate] = $date[$idate] - $lessDays 
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

sub days2months {
   my $date = shift;
   my $cal = shift;
   my $year = $date->[0];
   my $month = $date->[1];
   my $day = $date->[2];
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
   return ($numMonths,$date->[2]-$day);
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
   print "\n[1mOVERVIEW[0m
  Returns the number of records per month.\n
[1mUSAGE[0m
  tsplit.pl FILE\n\n";
   exit 1
}
