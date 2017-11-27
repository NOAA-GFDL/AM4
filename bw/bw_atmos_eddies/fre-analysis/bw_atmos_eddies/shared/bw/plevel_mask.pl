#!/usr/bin/perl
use strict;

# load module NCL if needed (must be version 6.2.1 or greater)
if (!grep /^ncarg\/.*/, split/:/, $ENV{"LOADEDMODULES"}) {
  load_module_init();
  print STDOUT "NOTE: loading module ncarg/6.2.1\n";
  module("load", "ncarg/6.2.1");
}
my $nclversion = `ncl -V`; chomp $nclversion;
die "ERROR: NCL version must be 6.2.1 or greater;" if ($nclversion lt "6.2.1");

# input and output file names
die "ERROR: missing input and out file names" if (@ARGV != 2); 
my ($ifile,$ofile) = @ARGV;

# run the script
write_script_file("pmask.ncl");
my $command = "ncl -Q \'ifile=\"$ifile\"\' \'ofile=\"$ofile\"\' pmask.ncl";
print STDOUT "NCL Version $nclversion\n";
print STDOUT "$command\n";
system ($command);
unlink "pmask.ncl";

#------------------------------------------------------
# writes the NCL file

sub write_script_file {
  my $file = shift;
  unlink $file if (-e $file);
  open (OUT,"> $file") || die "Cannot open $file for output";
  print OUT '
begin
  fi = addfile(ifile,"r")
  varnames = getfilevarnames(fi)
  ps = fi->ps
  ; copy input file format
  format = systemfunc("ncdump -k "+ifile)
  format = str_sub_str(str_sub_str(str_capital(str_sub_str(format,"-"," "))," ",""),"Model","")
  setfileoption("nc","Format",format)
  setfileoption("nc","HeaderReserveSpace",16384)
  ; open the output file
  fo = addfile(ofile,"c")
  filedimdef(fo,ps!0,-1,True)
  ; copy bounds
  do i = 0, 2
    if (isatt(ps&$ps!i$,"bounds")) then
      bnds = ps&$ps!i$@bounds
      fo->$bnds$ = fi->$bnds$
    end if
  end do
  ; loop thru all variables (process 4-D variables with pressure dimension)
  ndone = 0
  do iv = 0, dimsizes(varnames)-1
    dimnames = getfilevardims(fi,varnames(iv))
    if (dimsizes(dimnames) .eq. 4) then
      if (isfilevaratt(fi,dimnames(1),"long_name")) then
        if (fi->$dimnames(1)$@long_name .eq. "pressure") then
          data = fi->$varnames(iv)$
          plev = data&$dimnames(1)$
          do k = 0, dimsizes(plev)-1
            data(:,k,:,:) = where(plev(k) .gt. ps, data@_FillValue, data(:,k,:,:))
          end do
          ; output variable name (remove trailing "_substr")
          outname = varnames(iv)
          index = str_index_of_substr(outname,"_",-1)
          if (.not.ismissing(index)) then
            outname = str_get_cols(outname,0,index-1)
          end if
          if (isatt(data,"time_avg_info")) then
            delete(data@time_avg_info)
          end if 
          print(varnames(iv)+" ---> "+outname)
          ; output file
          fo->$outname$ = data
          ndone = ndone+1
          delete([/data,plev/])
        end if
      end if
    end if
    delete(dimnames)
  end do
  delete(fi)
  delete(fo)
  if (ndone .eq. 0) then
    print("No fields processed, output file not created")
    if (isfilepresent(ofile)) then
      system("/bin/rm -f "+ofile)
    end if
  end if
end
';
  return 0;
}

#------------------------------------------------------
# sets up environment for module load using perl

sub load_module_init {
  if (! defined $ENV{MODULE_VERSION}) { 
    $ENV{MODULE_VERSION_STACK}="3.1.6";
    $ENV{MODULE_VERSION}="3.1.6";
  } else {
    $ENV{MODULE_VERSION_STACK}=$ENV{MODULE_VERSION};
  }
  $ENV{MODULESHOME} = "/usr/local/Modules/".$ENV{MODULE_VERSION};
  if (! defined $ENV{MODULEPATH} ) { 
    $ENV{MODULEPATH} = `sed 's/#.*$//' $ENV{MODULESHOME}/init/.modulespath | awk 'NF==1{printf("%s:",$1)}'` 
  }
  if (! defined $ENV{LOADEDMODULES} ) { 
    $ENV{LOADEDMODULES} = ""; 
  }
  return 0;
}

#------------------------------------------------------

sub module {
  my $exec_prefix = "/usr/local/Modules/".$ENV{MODULE_VERSION};
  eval `$exec_prefix/bin/modulecmd perl @_`;
}

