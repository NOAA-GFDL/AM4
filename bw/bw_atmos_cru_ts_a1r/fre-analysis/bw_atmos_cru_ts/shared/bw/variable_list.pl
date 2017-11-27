#!/usr/bin/perl
use strict;

my $vars = $ARGV[0];
my $file = $ARGV[1];
my @outlist;
my $DEBUG = 0;
my $dump = `ncdump -h $file`;

foreach my $var ( split /,/,$vars ) {

   # add this variable to the list
   print "Adding var: $var\n" if ($DEBUG);
   push @outlist, $var;

   # find axes for this variable
   if ($dump =~ /\t\w+ $var\((.+)\)/) {
      my @axes = split /, /, $1;
      print "Checking axes: @axes\n" if ($DEBUG);

      # for each axes find edges, bounds, climatology names (var:att = "value" ;)
      foreach my $axis (@axes) {
         print "axis = $axis\n" if ($DEBUG);
         foreach my $att (qw/edges bounds climatology/) {
            if ($dump =~ /\t\t$axis:$att = "(.+)"/) {
               print "Adding att: $att   value: $1\n" if ($DEBUG);
               push @outlist, $1 if var_does_not_exist($1);
            }
         }
      }
   }

   # check variable for time average info
   if ($dump =~ /\t\t$var:time_avg_info = "(.+)"/) {
      foreach my $att ( split /, /, $1 ) {
         push @outlist, $att if var_does_not_exist($att);
      }
   }
}

my $outlist = join ",", @outlist;
print "$outlist\n";


sub var_does_not_exist {
   my $thisVar = shift;
   foreach my $var (@outlist) {
       return 0 if ($var eq $thisVar);
   }
   return 1;
}
       
