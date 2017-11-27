
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  compute monthly-mean climatology from a monthly timeseries

undef ("monthly_timeseries")
function monthly_timeseries(x:numeric,opt:logical)
local dimx,dimo,rank,time,ntime,nyrs,yrmody,xo,n,k,ks,ke,yr,mo,verbose,bounds
begin
  verbose = 0
  bounds = False
  if (opt) then
     if (isatt(opt,"verbose")) then
        verbose = opt@verbose
     end if
     if (isatt(opt,"bounds")) then
        bounds = opt@bounds
     end if
  end if

  dimx = dimsizes(x)
  rank = dimsizes(dimx)
  time = x&$x!0$
  ntime = dimx(0)
  yrmody = cd_calendar(time,-2) ; year,month,day (e.g., 19790101)
  nyrs = (toint(yrmody(ntime-1)/10000) - toint(yrmody(0)/10000)) + 1

  if (bounds .and. rank .ne. 2) then
     print("ERROR: invalid rank for bounds computation")
     status_exit(1)
  end if

  if (verbose .ge. 1) then
    print("Time series from "+yrmody(0)+" to "+yrmody(ntime-1))
    print("nyrs = "+nyrs)
  end if

  dimo = dimx
  dimo(0) = nyrs*12
  xo = new(dimo,typeof(x))
  to = new(nyrs*12,typeof(time))

  ; monthly averages
  ks = 0
  k = 0
  do n = 1, dimx(0)-1
     yr = toint(yrmody(n)/10000)
     mo = toint((yrmody(n)-yr*10000)/100)
     dy = yrmody(n) - (yr*10000 + mo*100)
     if (dy .eq. 1 .or. n .eq. ntime-1) then
        ke = n
        if (dy .eq. 1) then
          ke = ke-1
        end if
        if (verbose .ge. 3) then
           if (dy .eq. 1) then
              mo = mo-1
           end if
           if (mo .eq. 0) then
              mo = 12
              yr = yr-1
           end if
           print("yr/mo: "+yr+"/"+mo+", ndy: "+(ke-ks+1)+", ks:ke = "+ks+":"+ke)
        end if

        ; monthly average
        if (.not.bounds) then
           if (rank .eq. 1) then
              xo(k) = dim_avg_n_Wrap(x(ks:ke),0)
           else if (rank .eq. 2) then
              xo(k,:) = dim_avg_n_Wrap(x(ks:ke,:),0)
           else if (rank .eq. 3) then
              xo(k,:,:) = dim_avg_n_Wrap(x(ks:ke,:,:),0)
           else if (rank .eq. 4) then
              xo(k,:,:,:) = dim_avg_n_Wrap(x(ks:ke,:,:,:),0)
           end if
           end if
           end if
           end if
        else
           xo(k,0) = x(ks,0)
           xo(k,1) = x(ke,1)
        end if
        to(k) = dim_avg_Wrap(time(ks:ke))
        ks = ke+1
        k = k+1
     end if
  end do

  ; time values
  if (verbose .ge. 3) then
    printVarSummary(to)
  end if

  xo!0 = x!0
  xo&$xo!0$ = to
  return xo
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
 ;if (.not.isvar("yr1")) then
 ;  yr1 = 1979
 ;end if
 ;if (.not.isvar("yr2")) then
 ;  yr2 = 2008
 ;end if
  if (.not.isvar("verbose")) then
    verbose = 1
  end if
  if (usage) then
    print("USAGE: ncl ifile=STRING ofile=STRING yr1=INT yr2=INT compute_monthly_timeseries.ncl")
    exit
  end if

  if (isvar("yr1") .and. isvar("yr2")) then
    if (yr2 .lt. yr1) then
      print("ERROR: negative year range specified: yr2 < yr1")
      status_exit(1)
    end if
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
  banned_variables = (/timeName,"average_T1","average_T2","average_DT"/)

  ; check for time bounds attribute
  timeBounds = "_null_"
  if (isatt(time,"bounds")) then
     timeBounds = time@bounds
  end if

  ; determine start and end indices (jan 01 yr1 - dec 31 yr2)
  yrmody = cd_calendar(time,-2)

  if (isvar("yr1")) then
    do i = 0, dimsizes(yrmody)-1
      if (yrmody(i) .eq. yr1*10000+101) then
        kb = i
      end if
    end do
  else
    yr1 = toint(yrmody(0)/10000)
    if (yrmody(0) .eq. yr1*10000+101) then
      kb = 0
    end if
  end if

  if (isvar("yr2")) then
    do i = 0, dimsizes(yrmody)-1
      if (yrmody(i) .eq. yr2*10000+1231) then
        ke = i
      end if
    end do
  else
    yr2 = toint(yrmody(dimsizes(yrmody)-1)/10000)
    if (yrmody(dimsizes(yrmody)-1) .eq. yr2*10000+1231) then
      ke = dimsizes(yrmody)-1
    end if
  end if

  if (.not.isvar("kb") .or. .not.isvar("ke")) then
    print ("ERROR: unable to determine start/end time indices for model "+modelNames(n))
    status_exit(1)
  end if
 ;if (ke-kb .ne. 12*(yr2-yr1+1)-1) then
   ;print ("ERROR in determine_start_end: Start/end times could not be determined")
   ;print ("Dates: "+yrmody(0)+" to "+yrmody(dimsizes(yrmody)-1))
   ;status_exit(1)
 ;end if

  if (verbose .ge. 1) then
    print("Time series computed from "+yr1+" to "+yr2)
    print("File dates: "+yrmody(0)+" to "+yrmody(dimsizes(yrmody)-1))
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
   ;if (variables(n) .ne. timeName) then
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
   ;else ; variables(n) .eq. timeName
   ;    if (isatt(fi->$variables(n)$,"long_name")) then
   ;       varAtts@long_name = fi->$variables(n)$@long_name
   ;    else
   ;       varAtts@long_name = "Climatological time"
   ;    end if
   ;    varAtts@units = "months since 0000-01-01 00:00:00"
   ;    varAtts@climatological_years = yr1+" to "+yr2
   ;    natts = 3
   ;end if
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
  ;;;;;;;;;;;;;;;;;;;;;;

  Options = True

  do n = 0, dimsizes(variables)-1
    ; skip banned variables
    if (any(variables(n) .eq. banned_variables)) then
      continue
    end if
    ; skip variable/dimensions not used
    if (any(variables(n) .eq. dNames) .and. .not.any(variables(n) .eq. dNames(indx))) then
      continue
    end if

    ;;;;;;;;;;;;;;
    ; write data
    ;;;;;;;;;;;;;;

    ; copy static variables
    if (static_variable(fi,variables(n),timeName)) then
      if (verbose .ge. 1) then
        print ("writing static variable: "+variables(n))
      end if
      fo->$variables(n)$ = fi->$variables(n)$

    ; average and write time dimension values
    else if (variables(n) .eq. timeName) then
      if (verbose .ge. 1) then
        print ("writing time dimension variable: "+variables(n))
      end if
      Options@verbose = 0
      Options@bounds = False
      fo->$variables(n)$ = monthly_timeseries(fi->$variables(n)$(kb:ke),Options)

    ; time bound values
    else if (variables(n) .eq. timeBounds) then
      dims = getfilevardims(fi,variables(n))
      rank = dimsizes(dims)
      if (rank .ne. 2) then
        print("ERROR: invalid rank for bounds variable: "+variables(n))
        status_exit(1)
      end if
      if (verbose .ge. 1) then
        print ("writing time bounds variable: "+variables(n))
      end if
      Options@verbose = 0
      Options@bounds = True ; flag for bounds output
      fo->$variables(n)$ = monthly_timeseries(fi->$variables(n)$(kb:ke,:),Options)

    ; average and write all other non-static variable
    else
      if (verbose .ge. 1) then
        print ("writing time-varying variable: "+variables(n))
      end if
      dims = getfilevardims(fi,variables(n))
      rank = dimsizes(dims)
      Options@verbose = verbose
      Options@bounds = False
      if (rank .eq. 2) then
        fo->$variables(n)$ = monthly_timeseries(fi->$variables(n)$(kb:ke,:),Options)
      else if (rank .eq. 3) then
        fo->$variables(n)$ = monthly_timeseries(fi->$variables(n)$(kb:ke,:,:),Options)
      else if (rank .eq. 4) then
        fo->$variables(n)$ = monthly_timeseries(fi->$variables(n)$(kb:ke,:,:,:),Options)
      else
        print("ERROR: invalid rank: only variables with rank = 2,3,4 allowed")
        status_exit(1)
      end if
      end if
      end if
      end if
      delete([/dims,rank/])
    end if
    end if
  end do
  delete([/dNames,indx/])

end


