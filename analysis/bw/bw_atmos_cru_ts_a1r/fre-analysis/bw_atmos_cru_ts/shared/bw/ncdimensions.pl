#!/usr/bin/perl
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
Getopt::Long::Configure("bundling");

# default arguments
my %Opt = (HELP=>0, VERBOSE=>0, BOUNDS=>0);

# argument parsing
my $status = GetOptions ('h|help'       => \$Opt{HELP},
                         'V|verbose'    => \$Opt{VERBOSE},
                         'b|bounds'     => \$Opt{BOUNDS},
                         'v|variable=s' => \$Opt{variable});

# help message
if ($Opt{HELP}) {
   usage();
   exit 1;
}

# error checks
die "ERROR: no variable argument"          if !$Opt{variable};
die "ERROR: incorrect number of arguments" if (@ARGV != 1);

my $name = $Opt{variable};
my $file = $ARGV[0];

die "ERROR: file does not exist" if (!-e $file);

# dump the file header into a string variable
my $dump = `ncdump -h $file` || die "ERROR: could not ncdump file";
my $unlim = get_unlim_name($dump);
print STDERR "unlim = $unlim\n" if $Opt{VERBOSE};

my @axes;
my @bnds;

if ( $dump =~ /\t\w+ $name\((.+)\)/ ) {
     # parse the dimensions for this variable
     @axes = split /, /, $1;
     if ($Opt{BOUNDS}) {
       print STDERR "axes found = @axes\n" if $Opt{VERBOSE};
       # check for bounds attribute for each dimension
       foreach my $axis (@axes) {
         if ( $dump =~ /\t\t$axis:bounds = "(.*)"/ ) {
            my $bnd = $1;
            push @bnds, $bnd if ( $dump =~ /\t\w+ $bnd\(.+\)/ );
         } elsif ( $dump =~ /\t\t$axis:climatology = "(.*)"/ ) {
            print STDERR "WARNING: using climatology attribute as bounds for non-unlimited dimension\n" if ($axis ne $unlim);
            my $bnd = $1;
            push @bnds, $bnd if ( $dump =~ /\t\w+ $bnd\(.+\)/ );
         } else {
            print STDERR "No bounds found for dimension: $axis\n" if $Opt{VERBOSE};
         }
       }
     }
} else {
     print STDERR "WARNING: variable not found in file\n";
}

# output
if ($Opt{BOUNDS}) {
   print "@bnds\n" if (@bnds);
} else {
   print "@axes\n" if (@axes);
}

########################################################

sub get_unlim_name {
   my $dmp = shift;
   my $unlim;
   if ( $dump =~ /\t(.+) = UNLIMITED ; \/\/ \((\d+) currently\)/ ) {
      $unlim = $1;
   }
   return $unlim;
}

sub usage {
   my $cmdname = substr($0,rindex($0,"/")+1);
   print "
[1mOVERVIEW[0m
   Prints the dimension names for a variable in a netCDF file.

[1mUSAGE[0m
   $cmdname [-h] [-V] [-b] -v VAR in.nc

[1mREQUIRED[0m
   [1m-v[0m, [1m--variable[0m VAR    Variable name to process
   in.nc                 Input file name

[1mOPTIONAL FLAGS[0m
   [1m-h[0m, [1m--help[0m      Help message (this message)
   [1m-V[0m, [1m--verbose[0m   Additional messages output to STDERR
   [1m-b[0m, [1m--bounds[0m    Output the bounds name for each dimension (if available) instead of the dimension name
   
";
   return 0;
}


