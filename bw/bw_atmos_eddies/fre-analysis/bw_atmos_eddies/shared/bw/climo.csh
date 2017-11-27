#!/bin/csh -f
#
# Very simple monthly climatology script
#

set infile = $1
set outfile = $2

set timename = `ncdump -h $infile | grep UNLIMITED | awk '{print $1}'`
@ nt = `ncdump -h $infile | grep UNLIMITED | awk '{print $6}' | cut -c2-` - 1 

# loop through months
foreach ns ( 0 1 2 3 4 5 6 7 8 9 10 11)
  ncks -h -d $timename,$ns,$nt,12 $infile .climo.$infile.nc
  set mm = `printf "%2.2d" $ns`
  timavg.csh -amb -o .climo$mm.$infile.nc .climo.$infile.nc
  rm -f .climo.$infile.nc
end
ncrcat -h .climo??.$infile.nc $outfile
rm -f .climo??.$infile.nc

