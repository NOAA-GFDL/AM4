#!/usr/bin/perl
use strict;
use Getopt::Std;

my $help = 0;
our ($opt_d,$opt_t,$opt_v);
unless (getopts 'dtv:') {$help = 1;}

if (scalar @ARGV != 1 && scalar @ARGV != 2) {
   print "ERROR: No file name\n";
   usage();
   exit 1;
}
my $file = $ARGV[0];
my $ofile = "";
if (scalar @ARGV == 2) {
   $ofile = $ARGV[1];
}

if (!-e $file) {
   die ("File does not exist");
}

# debug mode turns on test-mode (no files changed)
if ($opt_d) {
   $opt_t = 1;
}

# split comma separated variable list
my @variables;
if ($opt_v) {
   @variables = split /,/, $opt_v;
}

# dump netcdf file
my $dump = `ncdump -h $file`;
my @dump = split /\n/, $dump;

# create a list of all variable and their type
my @varList;
foreach (@dump) {
    if (/^\t[a-z]* \w*\(.*\) ;/) {
       # extract variable type and name
       /^\t([a-z]*) (\w*)\((.*)\)/;
       my %V = ();
       $V{'type'} = $1;
       $V{'name'} = $2;
       push @varList, \%V;
    }
    # done when we reach the global attributes
    if (/ global attributes:/) {
       last;
    }
}

# debug: print variables found
if ($opt_d) {
   print "# vars = ".scalar(@varList)."\n";
}

my @ncatted_commands;
foreach my $ref (@varList) {
   my $n = $ref->{'name'};
   if (checkName($n,\@variables)) {
      # debug: print missing value and fillvalue
      if ($opt_d) {
         print "var: $n\n";
         if ($dump =~ /\n\t\t$n:missing_value = /) {
            $dump =~ /\n\t\t$n:missing_value = (.*) ;/;
            my $mval = $1;
            print "$n: missing_value = $mval\n";
         }
         if ($dump =~ /\n\t\t$n:_FillValue = /) {
            print "$n: _FillValue found\n";
         }
      }
      # check if fillvalue attribute should be added
      if ($dump =~ /\n\t\t$n:missing_value = / && $dump !~ /\n\t\t$n:_FillValue = /) {
         $dump =~ /\n\t\t$n:missing_value = (.*) ;/;
         my $mval = $1;
         my $t;
         if ($ref->{'type'} eq "double" || $ref->{'type'} eq "float" ||
             $ref->{'type'} eq "long"   || $ref->{'type'} eq "short") {
            $t = substr $ref->{'type'},0,1;
         } else {
            die "Cannot handle data type for variable $n";
         }
         my $value = $1; $value =~ s/$t$//; # remove type-letter from end of value
         push @ncatted_commands, "-a _FillValue,$n,c,$t,$value";
      }
   }
}

if (scalar @ncatted_commands > 0) {
   my $ncatted_opt = join ' ', @ncatted_commands;
   if ($opt_t) {
      print "ncatted -h $ncatted_opt $file $ofile\n";
   } else {
      system ("ncatted -h $ncatted_opt $file $ofile");
   }
}

#################################################

sub checkName {
   my $var = shift;
   my $vref = shift;
   if (scalar @$vref == 0) {
      return 1;
   }
   foreach (@$vref) {
      if ($_ eq $var) {
         return 1;
      }
   }
   return 0;
}

sub usage {
   print "
ncfillvalue.pl [-d] [-t] [-v var1,var2,...] file [outfile]

        -v   Comma separated list of variables (default all variables)
        -d   debug (implies test mode)
        -t   test-mode (nothing changed)
      file   netcdf file name

";
}
