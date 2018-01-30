
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  compute monthly-mean climatology from a monthly timeseries

undef ("monthly_climatology")
function monthly_climatology(x:numeric,verbose[1]:integer)
local dimx,rank,nyrs,yrmo,dimc,clim,wts,n,k,ks,ke
begin
  dimx = dimsizes(x)
  rank = dimsizes(dimx)
  nyrs = dimx(0)/12
  yrmo = cd_calendar(x&$x!0$,-1) ; yearmo (e.g., 197901)

  wts = new(dimx(0),float)

  if (verbose .ge. 1) then
    print("Time series from "+yrmo(0)+" to "+yrmo(dimx(0)-1))
  end if

  ; weight in days for each month
  do n = 0, nyrs-1
     wts(12*n:12*n+11) = (/31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./)
     if (mod(toint(yrmo(12*n+1)/100),4) .eq. 0) then  ; leap-year
        wts(12*n+1) = wts(12*n+1)+1.0
     end if
  end do

  ; define output array (12 climatological records)
  dimc = dimx
  dimc(0) = 12
  clim = new(dimc,typeof(x))

  ; monthly averages (temporary)
  do ks = 0, 11
     ke = ks+(nyrs-1)*12
     if (rank .eq. 2) then
        clim(ks,:) = dim_avg_wgt_n_Wrap(x(ks:ke:12,:),wts(ks:ke:12),0,0)
     else if (rank .eq. 3) then
        clim(ks,:,:) = dim_avg_wgt_n_Wrap(x(ks:ke:12,:,:),wts(ks:ke:12),0,0)
     else if (rank .eq. 4) then
        clim(ks,:,:,:) = dim_avg_wgt_n_Wrap(x(ks:ke:12,:,:,:),wts(ks:ke:12),0,0)
     end if
     end if
     end if
  end do

  clim!0 = x!0
  clim&$clim!0$ = fspan(0.5,11.5,12)
  clim&$clim!0$@units = "months since 0000-01-01 00:00:00"
  if (isatt(x&$x!0$,"long_name")) then
     clim&$clim!0$@long_name = x&$x!0$@long_name
  else
     clim&$clim!0$@long_name = "Climatological time"
  end if
  return clim
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  returns True if the variable has NO unlimited dimension

undef("static_variable")
function static_variable (fi:file,variable:string,timeDimName:string)
local dimensions, n
begin
  dimensions = getfilevardims(fi,variable)
  do n = 0, dimsizes(dimensions)-1
   ;if (isunlimited(fi,dimensions(n))) then
    if (dimensions(n) .eq. timeDimName) then
      return False
    end if
  end do
  return True
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

begin

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;  input arguments  ;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  usage = False
  if (.not.isvar("ifile")) then
    print("ERROR: ifile not specified")
    usage = True
  end if
  if (.not.isvar("ofile")) then
    print("ERROR: ofile not specified")
    usage = True
  end if
  if (.not.isvar("yr1")) then
    yr1 = 1979
  end if
  if (.not.isvar("yr2")) then
    yr2 = 2008
  end if
  if (.not.isvar("verbose")) then
    verbose = 1
  end if
  if (usage) then
    print("USAGE: ncl ifile=STRING ofile=STRING yr1=INT yr2=INT compute_monthly_climatology.ncl")
    exit
  end if

  if (yr2 .lt. yr1) then
     print("ERROR: negative year range specified: yr2 < yr1")
     status_exit(1)
  end if

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;  open the input file ;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  fi = addfile(ifile,"r")

  ; get time coord name
  if (.not.isvar("timeName")) then
    dimNames = getvardims(fi)
    do i = 0, dimsizes(dimNames)-1
      if (isunlimited(fi,dimNames(i))) then
        timeName = dimNames(i)
      end if
    end do
  end if
  if (.not.isvar("timeName")) then
    print ("ERROR: no time coordinate for input file")
    status_exit(1)
  end if
  ; read time coordinate data
  time = fi->$timeName$
  ; create list of variable not processed (usually time related)
  ; check for time bounds attribute
  if (isatt(time,"bounds")) then
     if (verbose .ge. 2) then
        print ("banning variables: "+timeName+", "+time@bounds)
     end if
     banned_variables = (/timeName,"average_T1","average_T2","average_DT",time@bounds/)
  else
     banned_variables = (/timeName,"average_T1","average_T2","average_DT"/)
  end if

  ; determine start and end indices (jan yr1 - dec yr2)
  yrmo = cd_calendar(time,-1)
  do i = 0, dimsizes(yrmo)-1
    if (yrmo(i) .eq. yr1*100+1) then
      kb = i
    end if
    if (yrmo(i) .eq. yr2*100+12) then
      ke = i
    end if
  end do
  if (.not.isvar("kb") .or. .not.isvar("ke")) then
    print ("ERROR: unable to determine start/end time indices for model "+modelNames(n))
    status_exit(1)
  end if
  if (ke-kb .ne. 12*(yr2-yr1+1)-1) then
    print ("ERROR in determine_start_end: Start/end times could not be determined")
    print ("Dates: "+yrmo(0)+" to "+yrmo(ny))
    status_exit(1)
  end if

  if (verbose .ge. 1) then
    print("Climatology computed from "+yr1+" to "+yr2)
    print("Time dimension ("+timeName+") indices used "+kb+" to "+ke+" (full indices are 0 to "+(dimsizes(time)-1)+")")
  end if

  ; open output file
  system("/bin/rm -f "+ofile)
  fo = addfile(ofile,"c")

  ; explicitly declare file definition mode (improve efficiency)
  setfileoption(fo,"DefineMode",True)

  ; predefine the coordinate variables and their dimensionality
  ; copy dimensions from input file
  dNames = getvardims(fi)
  dSizes = getfiledimsizes(fi)
  dUnlim = new(dimsizes(dNames),logical)
  dMask  = new(dimsizes(dNames),logical)
  dUnlim = False
  dMask  = False
  do n = 0, dimsizes(dNames)-1
   ;if (isunlimited(fi,dNames(n))) then
    if (dNames(n) .eq. timeName) then
      dSizes(n) = -1
      dUnlim(n) = True
    end if
  end do

  ; get variable names
  variables = getfilevarnames(fi)

  ; loop through variables and flag dimensions used
  do n = 0, dimsizes(variables)-1
    ; skip variables that are banned or also dimensions
    if (any(variables(n) .eq. banned_variables) .or. any(variables(n) .eq. dNames)) then
      continue
    end if
    ; flag dimensions used
    if (verbose .ge. 3) then
      print("variable: "+variables(n))
    end if
    varDims = getfilevardims(fi,variables(n))
    do m = 0, dimsizes(varDims)-1
      if (varDims(m) .ne. "ncl_scalar") then ; scalar dim can not be defined by user
        dMask(ind(dNames .eq. varDims(m))) = True
      end if
    end do
    delete(varDims)
  end do

  ; define array of indices of used dimensions
  indx = new(num(dMask),integer)
  m = 0
  do n = 0, dimsizes(dMask)-1
    if (dMask(n)) then
       indx(m) = n
       m = m+1
    else
       if (verbose .ge. 2) then
         print ("Skip dimension: "+dNames(n))
       end if
    end if
  end do

  ; print dimension info
  do n = 0, dimsizes(indx)-1
    if (verbose .ge. 2) then
      print("dimension: "+dNames(indx(n))+", size="+dSizes(indx(n))+", unlim="+dUnlim(indx(n)))
    end if
  end do

  ; define the dimensions
  filedimdef(fo,dNames(indx),dSizes(indx),dUnlim(indx))
  delete([/dSizes,dUnlim,dMask/])


  ; predefine the the dimensionality of the variables to be written out
  do n = 0, dimsizes(variables)-1
    ; skip banned variables
    if (any(variables(n) .eq. banned_variables)) then
      continue
    end if
    ; skip variable/dimensions not used
    if (any(variables(n) .eq. dNames) .and. .not.any(variables(n) .eq. dNames(indx))) then
      continue
    end if
    ; special treatment of time axis
    if (variables(n) .eq. timeName) then
       filevardef(fo, variables(n), "float", getfilevardims(fi,variables(n)))
       if (verbose .ge. 1) then
         print ("Registering variable: "+variables(n)+"  dims="+str_join(getfilevardims(fi,variables(n)),",")+"  type=float")
       end if
    else
       type = getfilevartypes(fi,variables(n))
       filevardef(fo, variables(n), type, getfilevardims(fi,variables(n)))
       if (verbose .ge. 1) then
         print ("Registering variable: "+variables(n)+"  dims="+str_join(getfilevardims(fi,variables(n)),",")+"  type="+type)
       end if
    end if
  end do

  ; copy attributes associated with each variable to the file
  ; all attributes associated with each variable will be copied
  atts_to_skip = (/"_FillValue","missing_value"/)
  do n = 0, dimsizes(variables)-1
    ; skip banned variables
    if (any(variables(n) .eq. banned_variables)) then
      continue
    end if
    ; skip variable/dimensions not used
    if (any(variables(n) .eq. dNames) .and. .not.any(variables(n) .eq. dNames(indx))) then
      continue
    end if

    varAtts = 0
    atts = getfilevaratts(fi,variables(n))
    if (variables(n) .ne. timeName) then
      noskip = True
      ; if the variable has no attributes there will actually be one missing attribute
      if (dimsizes(atts) .eq. 1) then
        if (ismissing(atts)) then
          noskip = False
        end if
      end if
      natts = 0
      if (noskip) then
        if (verbose .ge. 2) then
          print ("Registering attributes for variable "+variables(n))
        end if
        do ii = 1, dimsizes(atts)
          i = dimsizes(atts)-ii
          if (.not.any(atts_to_skip .eq. atts(i))) then
            thisAtt = fi->$variables(n)$@$atts(i)$
            varAtts@$atts(i)$ = thisAtt
            if (verbose .ge. 3) then
               print ("    "+atts(i)+":"+thisAtt)
            end if
            natts = natts+1
            delete(thisAtt)
          end if
        end do
      end if
    else ; variables(n) .eq. timeName
        if (isatt(fi->$variables(n)$,"long_name")) then
           varAtts@long_name = fi->$variables(n)$@long_name
        else
           varAtts@long_name = "Climatological time"
        end if
        varAtts@units = "months since 0000-01-01 00:00:00"
        varAtts@climatological_years = yr1+" to "+yr2
        natts = 3
    end if
    if (natts .gt. 0) then
      filevarattdef(fo, variables(n), varAtts)
    end if
    delete([/varAtts,atts/])
  end do
  if (verbose .ge. 2) then
    print ("Done registering variable attributess")
  end if

  ;;;;;;;;;;;;;;;;;;;;;;
  ; write data
  ; loop thru variables

  do n = 0, dimsizes(variables)-1
    ; skip banned variables
    if (any(variables(n) .eq. banned_variables)) then
      continue
    end if
    ; skip variable/dimensions not used
    if (any(variables(n) .eq. dNames) .and. .not.any(variables(n) .eq. dNames(indx))) then
      continue
    end if

    ; write data
    if (static_variable(fi,variables(n),timeName)) then
      if (verbose .ge. 1) then
        print ("writing static variable: "+variables(n))
      end if
      fo->$variables(n)$ = fi->$variables(n)$
    else if (variables(n) .eq. timeName) then
      if (verbose .ge. 1) then
        print ("writing time dimension variable: "+variables(n))
      end if
      fo->$variables(n)$ = fspan(0.5,11.5,12)
    else
      if (verbose .ge. 1) then
        print ("writing time-varying variable: "+variables(n))
      end if
      dims = getfilevardims(fi,variables(n))
      rank = dimsizes(dims)
      if (rank .eq. 2) then
        fo->$variables(n)$ = monthly_climatology(fi->$variables(n)$(kb:ke,:),verbose)
      else if (rank .eq. 3) then
        fo->$variables(n)$ = monthly_climatology(fi->$variables(n)$(kb:ke,:,:),verbose)
      else if (rank .eq. 4) then
        fo->$variables(n)$ = monthly_climatology(fi->$variables(n)$(kb:ke,:,:,:),verbose)
      else
        print("ERROR: invalid rank: only variables with rank = 2,3,4 allowed")
        status_exit(1)
      end if
      end if
      end if
      delete([/dims,rank/])
    end if
    end if
  end do
  delete([/dNames,indx/])

end


