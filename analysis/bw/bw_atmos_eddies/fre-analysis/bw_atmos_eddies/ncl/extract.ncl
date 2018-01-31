
function getpresdim(fi[1]:file,var[1]:string)
local dims, id, axis
begin
  dims = getfilevardims(fi,var)
  axis = -1
  do id = 0, dimsizes(dims)-1
    if (isfilevaratt(fi,dims(id),"axis")) then
      if (fi->$dims(id)$@axis .eq. "Z") then
        axis = id
        break
      end if
    else if (isfilevaratt(fi,dims(id),"cartesian_axis")) then
      if (fi->$dims(id)$@cartesian_axis .eq. "Z") then
        axis = id
        break
      end if
    end if
    end if

    if (isfilevaratt(fi,dims(id),"long_name")) then
      if (str_lower(fi->$dims(id)$@long_name) .ne. "pressure") then
        axis = -1
      end if
    end if
  end do
  axis@_FillValue = -1
  return axis
end

; INPUT: ifile[1]:string, ofile[1]:string, level[1]:integer
begin
 ;ifile = "va_unmsk.19810101-19831231.nc"
 ;ofile = "out."+level+".nc"
 ;level = 250
  verbose = 2
  ; open file
  fi = addfile(ifile,"r")
  dimensions = getvardims(fi)
  variables = getfilevarnames(fi)
  ; copy input file format
  format = systemfunc("ncdump -k "+ifile)
  format = str_sub_str(str_sub_str(str_capital(str_sub_str(format,"-"," "))," ",""),"Model","")
  setfileoption("nc","Format",format)
  setfileoption("nc","HeaderReserveSpace",16384)
  ; open the output file
  system("rm -f "+ofile)
  fo = addfile(ofile,"c")
  ; define unlimited dimension
  do id = 0, dimsizes(dimensions)-1
    if (isunlimited(fi,dimensions(id))) then
      filedimdef(fo,dimensions(id),-1,True)
    end if
  end do
  ; loop through all variables
  do iv = 0, dimsizes(variables)-1
    if (any(variables(iv) .eq. dimensions)) then
      continue
    end if 
    kpdim = getpresdim(fi,variables(iv))
    if (.not.ismissing(kpdim)) then
      newName = variables(iv)+level
      data = fi->$variables(iv)$
      plevels = data&$data!kpdim$
      print("Processing: "+variables(iv)+" --> "+newName+" (axis="+data!kpdim+")")
      dims = dimsizes(data)
      rank = dimsizes(dims)
      ; convert level to vertical dimemsion units
      xlevel = level
      if (isatt(data&$data!kpdim$,"units")) then
        units = str_lower(data&$data!kpdim$@units)
        if (any(units .eq. (/"pa","pascal","pascals"/))) then
          xlevel = xlevel*100
        end if
      else
        if (tofloat(max(plevels)) .gt. 50000.) then
          xlevel = xlevel*100
        end if
      end if
      ; find exact level
      kplev = ind(xlevel .eq. plevels)
      print("  kplev="+kplev)
      if (ismissing(kplev)) then
        print ("ERROR: extract_level: exact level not found: xlevel="+xlevel)
        system("rm -f "+ofile)
        status_exit(1)
      end if
      ; extract data
      if (rank .eq. 1) then
        y = new(1,typeof(x))
        y = data(kplev)
      else if (rank .eq. 2 .and. kpdim .eq. 0) then
        y = data(kplev,:)
      else if (rank .eq. 3 .and. kpdim .eq. 0) then
        y = data(kplev,:,:)
      else if (rank .eq. 4 .and. kpdim .eq. 0) then
        y = data(kplev,:,:,:)
      else if (rank .eq. 4 .and. kpdim .eq. 1) then
        y = data(:,kplev,:,:)
      else if (rank .eq. 5 .and. kpdim .eq. 0) then
        y = data(kplev,:,:,:,:)
      else if (rank .eq. 5 .and. kpdim .eq. 1) then
        y = data(:,kplev,:,:,:)
      else if (rank .eq. 5 .and. kpdim .eq. 2) then
        y = data(:,:,kplev,:,:)
      else
        print ("ERROR: extract_level: unsupported rank ("+rank+") and/or kpdim ("+kpdim+")")
        status_exit(1)
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      fo->$newName$ = y
      delete(y)
    else
      if (verbose .gt. 0) then
        print("Copying: "+variables(iv))
      end if
      fo->$variables(iv)$ = fi->$variables(iv)$
    end if
  end do
end

