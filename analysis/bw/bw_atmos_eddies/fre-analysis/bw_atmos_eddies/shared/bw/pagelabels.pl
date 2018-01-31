#!/usr/bin/perl
use Getopt::Std;
use strict;

our($opt_h,$opt_O,$opt_l);
my $help = 0;
my ($ifile,$ofile);
my $numPages;

my (@labels);

unless (getopts 'hOl:') {
   $help = 1;
}

if (scalar(@ARGV) == 2) {
   ($ifile,$ofile) = @ARGV;
} elsif (scalar(@ARGV) == 1) {
   $ifile = $ARGV[0];
   $ofile = $ifile;
} else {
   $help = 1;
   print "Input and output files must be specified!\n";
}

if ($help) {
   my $cmd = substr($0,rindex($0,"/")+1);
   print "USAGE: $cmd [-O] -l LAB1,LAB2,.. IN.nc [OUT.nc]\n";
   exit 1;
}

die "Must specify labels" if (! $opt_l);
die "Input file does not exists!" if (! -e $ifile);
die "Can not overwrite file without the -O option!" if (-e $ofile && !$opt_O);

my $text = `cat $ifile`; chomp $text;
my @lines = split /\n/, $text;

# number of pages
if ($text =~ /%%Pages: (\d+)\n/) {
   $numPages = $1;
} else {
   die "Could not determine number of pages";
}

# user specified labels
my @labels = split /,/, $opt_l;
die "Number of labels not equal number of pages" if (@labels != $numPages);

# output file
open (OUT,"> $ofile") || "Cannot open $ofile for output";

my $i = 1;
foreach my $line (@lines) {
   if ($line =~ /%%Page: \d+ \d+/) {
      $line =~ s/%%Page: \d+ /%%Page: $labels[$i-1] /;
      $i++;
   }
   print OUT "$line\n";
}
      
