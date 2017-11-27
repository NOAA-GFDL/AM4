
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  compute vorticity

undef ("compute_vort")
function compute_vort(u:numeric,v:numeric,lm:logical)
local uu,vv,dims,rank,vort
begin
  uu = u
  vv = v
  if (dimsizes(dimsizes(lm)) .eq. 2) then
    dims = dimsizes(uu)
    rank = dimsizes(dims)
    ; mask is True where you want data
    uu = where(conform(uu,lm,(/rank-2,rank-1/)), uu, uu@_FillValue)
    vv = where(conform(vv,lm,(/rank-2,rank-1/)), vv, vv@_FillValue)
  end if
 ;vort = u
  vort = uv2vr_cfd( uu, vv, uu&$uu!1$, uu&$uu!2$, 0 ) * 1.e6
  copy_VarCoords(uu,vort)
  if (isatt(u,"cell_methods")) then
    vort@cell_methods = uu@cell_methods
  end if
  return vort
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
;;  return the name of the time coordinate

undef("get_time_coordinate_name")
function get_time_coordinate_name (fi:file)
local dnames,i,timeName
begin
  dnames = getvardims(fi)
  do i = 0, dimsizes(dnames)-1
    if (isunlimited(fi,dnames(i))) then
      timeName = dnames(i)
    end if
  end do
  if (.not.isvar("timeName")) then
    print ("ERROR in get_time_coordinate_name: no time coordinate for file: "+getfilepath(fi))
    status_exit(1)
  end if
  return timeName
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

begin

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;  input arguments  ;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  usage = False
  if (.not.isvar("ufile")) then
    print("ERROR: ufile not specified")
    usage = True
  end if
  if (.not.isvar("vfile")) then
    print("ERROR: vfile not specified")
    usage = True
  end if
  if (.not.isvar("ofile")) then
    print("ERROR: ofile not specified")
    usage = True
  end if
  if (.not.isvar("uvar")) then
    print("ERROR: uvar not specified")
    usage = True
  end if
  if (.not.isvar("vvar")) then
    print("ERROR: vvar not specified")
    usage = True
  end if
  if (.not.isvar("verbose")) then
    verbose = 1
  end if
  if (usage) then
    print("USAGE: ncl ufile=STRING vfile=STRING ofile=STRING uvar=STRING vvar=STRING compute_vorticity.ncl")
    exit
  end if


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;  open the input file ;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  fu = addfile(ufile,"r")

  ; check that variable "uvar" exists
  if (.not.isfilevar(fu,uvar)) then
    print("ERROR: variable "+uvar+" does not exist in file "+ufile)
    status_exit(1)
  end if

  ; get the dimensions of the zonal wind variable
  udims = getfilevardimsizes(fu,uvar)

  ; get time coord name
  timeName = get_time_coordinate_name(fu)

  ; mask land areas?
  if (isvar("maskvar")) then
    datamask = where( fu->$maskvar$ .ge. 10., False, True )  ; True where you want data
  else
    datamask = False ; scalar = no masking
  end if

  ; create list of variable not processed (usually time related)
  if ( ufile .eq. vfile) then
    banned_variables = (/timeName,"average_T1","average_T2","average_DT",vvar/)
  else
    banned_variables = (/timeName,"average_T1","average_T2","average_DT"/)
  end if

  ; check for time bounds attribute
  time = fu->$timeName$
  timeBounds = "_null_"
  if (isatt(time,"bounds")) then
     timeBounds = time@bounds
  end if

  ; open meridional wind file
  fv = addfile(vfile,"r")

  ; check that variable "vvar" exists in file
  if (.not.isfilevar(fv,vvar)) then
    print("ERROR: variable "+vvar+" does not exist in file "+vfile)
    status_exit(1)
  end if

  ; and has the same size as variable "uvar"
  vdims = getfilevardimsizes(fv,vvar)
  if (sum(where(udims .ne. vdims, 1, 0)) .gt. 0) then
    print("ERROR: variables "+uvar+" and "+vvar+" do not have the same size")
    print("udims= "+udims)
    print("vdims= "+vdims)
    status_exit(1)
  end if

  if (verbose .ge. 1) then
    yrmo = cd_calendar(time,-1)
    print("File dates: "+yrmo(0)+" to "+yrmo(dimsizes(yrmo)-1))
  end if

  ; open output file
  system("/bin/rm -f "+ofile)
  fo = addfile(ofile,"c")

  ; explicitly declare file definition mode (improve efficiency)
  setfileoption(fo,"DefineMode",True)

  ; predefine the coordinate variables and their dimensionality
  ; copy dimensions from input file
  dNames = getvardims(fu)
  dSizes = getfiledimsizes(fu)
  dUnlim = new(dimsizes(dNames),logical)
  dMask  = new(dimsizes(dNames),logical)
  dUnlim = False
  dMask  = False
  do n = 0, dimsizes(dNames)-1
    if (dNames(n) .eq. timeName) then
      dSizes(n) = -1
      dUnlim(n) = True
    end if
  end do

  ; get variable names
  variables = getfilevarnames(fu)

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
    varDims = getfilevardims(fu,variables(n))
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
    ; special treatment of zonal wind variable
    if (variables(n) .eq. uvar) then
      filevardef(fo, "vort", "float", getfilevardims(fu,variables(n)))
      if (verbose .ge. 1) then
        print ("Registering variable: vort,  dims="+str_join(getfilevardims(fu,variables(n)),",")+",  type=float")
      end if
    else
    ; if ( ufile .eq. vfile .and. variables(n) .eq. vvar) then
    ;   ; skip if ufile=vfile and vvar
    ; else
        type = getfilevartypes(fu,variables(n))
        filevardef(fo, variables(n), type, getfilevardims(fu,variables(n)))
        if (verbose .ge. 1) then
          print ("Registering variable: "+variables(n)+",  dims="+str_join(getfilevardims(fu,variables(n)),",")+",  type="+type)
        end if
    ; end if
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
    atts = getfilevaratts(fu,variables(n))
    if (variables(n) .ne. uvar) then
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
            thisAtt = fu->$variables(n)$@$atts(i)$
            varAtts@$atts(i)$ = thisAtt
            if (verbose .ge. 3) then
               print ("    "+atts(i)+":"+thisAtt)
            end if
            natts = natts+1
            delete(thisAtt)
          end if
        end do
      end if
      if (natts .gt. 0) then
        filevarattdef(fo, variables(n), varAtts)
      end if

    else ; variables(n) .eq. uvar
        varAtts@comment = "relative vorticity using surface wind stress components"
        varAtts@units = "1/s"
        varAtts@long_name = "Vorticity"
        filevarattdef(fo, "vort", varAtts)
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
    if (static_variable(fu,variables(n),timeName)) then
      if (verbose .ge. 1) then
        print ("writing static variable: "+variables(n))
      end if
      fo->$variables(n)$ = fu->$variables(n)$

    ; compute vorticity
    else if (variables(n) .eq. uvar) then
      dims = getfilevardims(fu,variables(n))
      rank = dimsizes(dims)
      if (rank .lt. 3) then
         print("ERROR: invalid rank for velocity component variable: "+variables(n))
         status_exit(1)
      end if
      if (verbose .ge. 1) then
        print ("computing vorticity from: "+variables(n)+" and "+vvar)
      end if
      fo->vort = compute_vort (fu->$variables(n)$, fv->$vvar$, datamask)
      delete([/dims,rank/])

    ; write all other non-static variables
    else
      if (verbose .ge. 1) then
        print ("writing time-varying variables: "+variables(n))
      end if
      fo->$variables(n)$ = fu->$variables(n)$
    end if
    end if
  end do
  delete([/dNames,indx/])

end


