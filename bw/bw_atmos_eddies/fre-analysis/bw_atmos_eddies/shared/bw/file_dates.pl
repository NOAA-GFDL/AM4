#!/usr/bin/perl
use strict;
use Getopt::Long;
Getopt::Long::Configure("bundling");

#-------------------
# argument parsing
#-------------------

my %Opt = (HELP=>0);
my $status = GetOptions ('h|help'   => \$Opt{HELP},
                         'f|freq=s' => \$Opt{freq},
                         'Y|yr1=i'  => \$Opt{yr1},
                         'Z|yr2=i'  => \$Opt{yr2});

usage() if $Opt{HELP};

# argument checking
my ($type,$databegyr,$dataendyr,$datachunk);
if (scalar(@ARGV) == 3 && $Opt{freq}) {
   $type = $Opt{freq};
   ($databegyr,$dataendyr,$datachunk) = @ARGV;
} elsif (scalar(@ARGV) == 4 && !$Opt{freq}) {
   ($type,$databegyr,$dataendyr,$datachunk) = @ARGV;  # backward compatibility
} else {
   usage();
}

my ($mdh1, $mdh2);
if ($type eq "monthly" || $type eq "month") {
   $mdh1 = "01";
   $mdh2 = "12";
} elsif ($type eq "daily" || $type eq "day") {
   $mdh1 = "0101";
   $mdh2 = "1231";
} elsif ($type eq "hour" || $type eq "hr") {
   $mdh1 = "010100";
   $mdh2 = "123123";
} elsif ($type eq "year" || $type eq "yr") {
   $mdh1 = "";
   $mdh2 = "";
} else {
   print STDERR "ERROR: Invalid time interval designation\n";
   print STDERR "       Must be hour, daily, month or year\n";
   exit 1;
}

my $y1 = $databegyr;
for (my $y1 = $databegyr; $y1 <= $dataendyr; $y1 = $y1+$datachunk) {
   my $y2 = ($y1 + $datachunk) - 1;
   if ($Opt{yr1}) {
      next if $Opt{yr1} > $y2;
   }
   if ($Opt{yr2}) {
      last if $Opt{yr2} < $y1;
   }
   print  sprintf "%4.4d%s-%4.4d%s\n", $y1, $mdh1, $y2, $mdh2;
}

  
sub usage {
   my $scriptname = "file_dates.pl";
   print "
Purpose:  Returns the dates for time series files.

Syntax :  $scriptname -f freq [-Y yr1] [-Z yr2] databegyr dataendyr datachunk

          freq        = hour, daily, month or year
          yr1,yr2     = first and last year of analysis
          databegyr = first year of archive timeseries data
          dataendyr = last year of archive timeseries data
          datachunk = number of years in each archive timeseries data file

Example:  $scriptname -f month -Y 1986 -Z 1995 1981 2000 10

 output:  198101-199012
          199101-200012

";
   exit 1;
}
