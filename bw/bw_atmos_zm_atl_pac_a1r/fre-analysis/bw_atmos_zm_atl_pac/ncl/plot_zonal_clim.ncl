
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$BW_PACKAGE_ROOT/ncl/compute_quartiles.ncl"
load "$BW_PACKAGE_ROOT/ncl/shared.ncl"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;; start of script ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

begin

if (.not.isvar("yr1")) then
   yr1 = 1979
end if
if (.not.isvar("yr2")) then
   yr2 = 2008
end if
if (.not.isvar("model_descriptor")) then
   model_descriptor = "c96L48_awgp2a (1981-2000)"
end if
if (.not.isvar("psfile")) then
   psfile = "regave"
else
   ; strip ending .ps
   if (str_get_cols(psfile,-3,-1) .eq. ".ps") then
      psfile = str_get_cols(psfile,0,-4)
   end if
end if

if (.not.isvar("ensemble")) then
  ensemble = True   ; ensemble mean of CMIP5 models?
end if

if (.not.isvar("variable")) then
  variable = "pr"   ; pr,tauu
end if

if (.not.isvar("ANOMALY")) then
  ANOMALY = False
end if

if (.not.isvar("ALL_OBS")) then
   ALL_OBS     = False ; when True also plot COADS & ERS data
end if                 ; by default MERRA & ECMWF are plotted (only valid for tauu)

USE_NEW_COLORS = True  ; when False the original Kushner colors are used (only valid for tauu)

if (ensemble) then
  ; ensemble mean of models
  files = systemfunc("/bin/ls cmip5/"+variable+"_Amon_*_amip_r0i0p0_"+yr1+"01-"+yr2+"12-clim.nc")
  numFiles = dimsizes(files)
  cmip5_descriptor = "CMIP5 ("+numFiles+" models, "+yr1+"-"+yr2+")"
else
  ; use all model realizations
  files = systemfunc("/bin/ls models/"+variable+"_Amon_*_amip_r*_"+yr1+"01-"+yr2+"12-clim.nc")
  numFiles = dimsizes(files)
  cmip5_descriptor = "CMIP5 ("+numFiles+" runs, "+yr1+"-"+yr2+")"
end if
if (numFiles .eq. 0) then
   print("ERROR: no cmip5 model files found")
   exit
end if
modelNames = new(numFiles,string)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; output grid latitudes - when plotting cmip5 models
; all cmip5 model will be interpolated onto this grid
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  grid_delta = 0.5

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; store information about regions in data structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  REGIONS = True
  if (variable .eq. "pr") then
    REGIONS@name     = (/"pacific", "atlantic", "tropics"/)
    REGIONS@latbeg   = (/ -70.    ,  -70.     ,   -20.   /)
    REGIONS@latend   = (/  80.    ,   80.     ,    20.   /)
    REGIONS@lonbeg   = (/ 115.    ,  -70.     ,     0.   /)
    REGIONS@lonend   = (/ 280.    ,   20.     ,   360.   /)
    REGIONS@dimavg   = (/   2     ,    2      ,     1    /) ; 1=meridional avg, 2=zonal avg
    REGIONS@lndsea   = (/   1     ,    1      ,     1    /) ; land=0,sea=1
    REGIONS@longname = (/"Pacific (115E-100W)", "Atlantic (70W-20E)", "Tropical Ocean (20S-20N)"/)
  else if (variable .eq. "curl") then
    REGIONS@name     = (/"pacific", "atlantic", "global"/)
    REGIONS@latbeg   = (/ -70.    ,  -70.     ,   -70.   /)
    REGIONS@latend   = (/  70.    ,   70.     ,    70.   /)
    REGIONS@lonbeg   = (/ 115.    ,  -70.     ,     0.   /)
    REGIONS@lonend   = (/ 280.    ,   20.     ,   360.   /)
    REGIONS@dimavg   = (/   2     ,    2      ,     2    /) ; 1=meridional avg, 2=zonal avg
    REGIONS@lndsea   = (/   1     ,    1      ,     1    /) ; land=0,sea=1
    REGIONS@longname = (/"Pacific (115E-100W)", "Atlantic (70W-20E)", "Global Ocean"/)
   ;REGIONS@longname = (/"Pacific (115E-100W)", "Atlantic (70W-20E)", "Tropical Ocean (20S-20N)"/)
  else if (variable .eq. "tauu") then
    REGIONS@name     = (/"pacific", "atlantic", "global"/)
    REGIONS@latbeg   = (/ -70.    ,  -70.     ,   -55.   /)
    REGIONS@latend   = (/  80.    ,   80.     ,    80.   /)
    REGIONS@lonbeg   = (/ 115.    ,  -70.     ,     0.   /)
    REGIONS@lonend   = (/ 280.    ,   20.     ,   360.   /)
    REGIONS@dimavg   = (/   2     ,    2      ,     2    /) ; 1=meridional avg, 2=zonal avg
    REGIONS@lndsea   = (/   1     ,    1      ,     0    /) ; land=0,sea=1
    REGIONS@longname = (/"Pacific (115E-100W)", "Atlantic (70W-20E)", "Land (Global)"/)
  else
    print("ERROR: invalid variable: "+variable)
    exit
  end if
  end if
  end if


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  set up data structures to hold CMIP5 data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  
  CMIP5GRID = True
  CMIP5DATA = True
  CMIP5DATA@descriptor = cmip5_descriptor
  do reg = 0, dimsizes(REGIONS@name)-1
    regionName = REGIONS@name(reg)
    ; define grid
    if (REGIONS@dimavg(reg) .eq. 1) then
      npts = toint(((REGIONS@lonend(reg)-REGIONS@lonbeg(reg))/grid_delta)+0.5)
      if (REGIONS@lonend(reg)-REGIONS@lonbeg(reg) .lt. 359.999) then
         npts = npts + 1
         CMIP5GRID@$regionName$ = fspan(REGIONS@lonbeg(reg),REGIONS@lonend(reg),npts)
         print("defining  CMIP5GRID@"+regionName)
      else
         delh = (REGIONS@lonend(reg)-REGIONS@lonbeg(reg))/tofloat(2*npts)  ; delta/2
         CMIP5GRID@$regionName$ = fspan(REGIONS@lonbeg(reg)+delh,REGIONS@lonend(reg)-delh,npts)
         print("defining  CMIP5GRID@"+regionName)
         cyclic_attrib = regionName+"_cyclic"
         CMIP5GRID@$cyclic_attrib$ = True
         print("defining  CMIP5GRID@"+cyclic_attrib)
      end if
    else if (REGIONS@dimavg(reg) .eq. 2) then
      npts = toint(((REGIONS@latend(reg)-REGIONS@latbeg(reg))/grid_delta)+0.5) + 1
      CMIP5GRID@$regionName$ = fspan(REGIONS@latbeg(reg),REGIONS@latend(reg),npts)
      print("defining  CMIP5GRID@"+regionName)
    end if
    end if
      
    rdat = new((/numFiles,5,npts/),float)
    rdat = rdat@_FillValue
    print("defining  CMIP5DATA@"+regionName)
    CMIP5DATA@$regionName$ = rdat
    delete(rdat)
  end do

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;  LOOP THROUGH CMIP5 MODELS  ;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  do n = 0, numFiles-1
    nameSplit = str_split(files(n),"_")
    modelNames(n) = nameSplit(2)
    print ("Model: "+modelNames(n))

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; read the data (yr1-yr2)
    ; reading from climatology file
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;
    fm = addfile(files(n),"r")
    data = fm->$variable$
    if (any(variable .eq. (/"tauu","curl"/)) .and. modelNames(n) .eq. "BNU-ESM") then
      data = data*(-1.0)
    end if
    if (ANOMALY) then
      data = anomaly_Wrap(fm,data,yr1,yr2)
    end if
    if (variable .eq. "curl") then
      data = data*(-1)
    end if

    opt = True
    opt@yr1 = yr1
    opt@yr2 = yr2
    sdata = monthly_clim_to_seasonal(data,opt)

    ; read/setup land/sea masks (already as fraction, but not for cmip models)
    land = fm->sftlf
    if (max(land) .gt. 50.) then
      land = land / 100.
    end if
    diml = dimsizes(land)
    landsea = new((/2,diml(0),diml(1)/),typeof(land))
    landsea(0,:,:) = land
    landsea(1,:,:) = land
    landsea(1,:,:) = 1.-landsea(1,:,:)

    ; regional/zonal averages (dimensioned [season, latitude])
    do reg = 0, dimsizes(REGIONS@name)-1
      regionName = REGIONS@name(reg)
      flag = False
      if (reg .eq. 2) then
         flag = False
      end if
      rdat = CMIP5DATA@$regionName$
      dat = regional_average ( fm, sdata, landsea(REGIONS@lndsea(reg),:,:), REGIONS@latbeg(reg), REGIONS@latend(reg), \
                                                                            REGIONS@lonbeg(reg), REGIONS@lonend(reg), REGIONS@dimavg(reg), flag )
      cyclic = False
      if (isatt(CMIP5GRID,regionName+"_cyclic")) then
         cyclic = True
      end if
      rdat(n,:,:) = linint1_n_Wrap (dat&$dat!1$, dat, cyclic, CMIP5GRID@$regionName$, 0, 1) ; season x latitude  (5 x npts)
      CMIP5DATA@$regionName$ = rdat
      delete([/dat,rdat/])
    end do

    delete([/fm,data,land,landsea,sdata/])

  end do

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; compute quartiles for cmip5 models
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  QUARTILES = compute_quartiles(REGIONS@name, CMIP5DATA)
  


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  set up observed data sets
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  OBSDATA = True

  ; precipitation
  if (variable .eq. "pr") then
    OBSDATA@files  = (/"cmap", "gpcp", "merra"/)
    OBSDATA@land   = (/  True, True, True /) ; obs file has land data
    OBSDATA@colors = (/"darkgreen","blue","green"/) ; red is reserved for model
    obs_file_var = variable

  ; eastward wind stress
  else if (variable .eq. "tauu" .or. variable .eq. "curl") then
    if (ALL_OBS) then
      OBSDATA@files  = (/"ecmwf", "merra",  "coads",    "ers"  /)
      OBSDATA@land   = (/  True,    True,    False ,    False  /) ; obs file has land data
      if (USE_NEW_COLORS) then
        OBSDATA@colors = (/ "blue", "green", "darkgreen", "brown4"/) ; red is reserved for model
      else
        OBSDATA@colors = (/"blue",  "green", "red", "darkgreen"/) ; Kushner plot colors
      end if
    else
      OBSDATA@files  = (/"ecmwf", "merra"/)
      OBSDATA@land   = (/  True,   True /) ; obs file has land data
      OBSDATA@colors = (/ "blue", "green"/) ; red is reserved for model
    end if
   ;obs_file_var = "wind_stress"
    obs_file_var = variable
  end if
  end if

  ; loop though observed data sets
  do i = 0, dimsizes(OBSDATA@files)-1
    print ("OBS: "+OBSDATA@files(i))

    ; read obs data
    fm = addfile("obs/"+obs_file_var+"."+OBSDATA@files(i)+".climatology.nc","r")
    data = fm->$variable$
    if (ANOMALY) then
      data = anomaly_Wrap(fm,data,0,0)
    end if
    if (variable .eq. "curl") then
      data = data*(-1)
    end if
    ; read 2d land fraction (NOT percent)
    land = fm->sftlf
    if (max(land) .gt. 50.) then
      land = land/100.
    end if
    diml = dimsizes(land)
    landsea = new((/2,diml(0),diml(1)/),typeof(land))
    landsea(0,:,:) = land
    landsea(1,:,:) = land
    landsea(1,:,:) = 1.-landsea(1,:,:)
    sdata = monthly_clim_to_seasonal(data,False)

    do reg = 0, dimsizes(REGIONS@name)-1
      dat = regional_average ( fm, sdata, landsea(REGIONS@lndsea(reg),:,:), REGIONS@latbeg(reg), REGIONS@latend(reg), \
                                                                            REGIONS@lonbeg(reg), REGIONS@lonend(reg), REGIONS@dimavg(reg), False )
      ; save obs data for this region
      attname = OBSDATA@files(i)+"_"+REGIONS@name(reg)
      OBSDATA@$attname$ = dat
      ; save obs grid for this region
      attname = OBSDATA@files(i)+"_"+REGIONS@name(reg)+"_grid"
      OBSDATA@$attname$ = dat&$dat!1$
      delete(dat)
    end do

    delete([/fm,data,land,landsea,sdata/])
  end do


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   set up model data set
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  MODELDATA = True
  MODELDATA@descriptor = model_descriptor

  if (USE_NEW_COLORS) then
    MODELDATA@color = "red"
  else
    MODELDATA@color = "black"
  end if

  fm = addfile("model_climatology.nc","r")
  print ("Model: "+MODELDATA@descriptor)
  data = fm->$variable$
  if (ANOMALY) then
    data = anomaly_Wrap(fm,data,0,0)
  end if
  if (variable .eq. "tauu") then
    data = data*(-1)   ; this is assuming GFDL output files
  end if               ; where only the variable has been renamed
  land = fm->sftlf
  if (max(land) .gt. 50.) then
    land = land/100.
  end if
  diml = dimsizes(land)
  landsea = new((/2,diml(0),diml(1)/),typeof(land))
  landsea(0,:,:) = land
  landsea(1,:,:) = land
  landsea(1,:,:) = 1.-landsea(1,:,:)
  opt = False
  sdata = monthly_clim_to_seasonal(data,opt)
  do reg = 0, dimsizes(REGIONS@name)-1
    dat = regional_average ( fm, sdata, landsea(REGIONS@lndsea(reg),:,:), REGIONS@latbeg(reg), REGIONS@latend(reg), \
                                                                          REGIONS@lonbeg(reg), REGIONS@lonend(reg), REGIONS@dimavg(reg), False )
    regionName = REGIONS@name(reg)
    MODELDATA@$regionName$ = dat
    attname = regionName+"_grid"
    MODELDATA@$attname$ = dat&$dat!1$
    delete(dat)
  end do
  delete([/fm,data,land,landsea/])


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;; PLOTTING SECTION ;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;     create 3 panel plots for each season
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  OPTIONS = True
  OPTIONS@ANOMALY = ANOMALY

  ; open postscript file for output
  wks = gsn_open_wks("ps",psfile)

  ; plot curves for all seasons
  plot_curves ( wks, variable, REGIONS, MODELDATA, OBSDATA, CMIP5DATA, CMIP5GRID, QUARTILES, OPTIONS )


end
