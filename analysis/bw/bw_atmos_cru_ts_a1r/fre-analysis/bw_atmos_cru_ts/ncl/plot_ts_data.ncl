
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$BW_PACKAGE_ROOT/shared/bw/ncl/define_regions.ncl"
load "$BW_PACKAGE_ROOT/shared/bw/ncl/horiz_interp.ncl"
load "$BW_PACKAGE_ROOT/shared/bw/ncl/interp_Wrap.ncl"
load "$BW_PACKAGE_ROOT/ncl/compute_climatology_statistics.ncl"
load "$BW_PACKAGE_ROOT/ncl/contour_levels.ncl"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("initialize_resources")
function initialize_resources (x[*][*]:numeric)
local res, dimx, numLat, numLon, minLat, maxLat, minLon, maxLon, aspect, do_polar
begin

  lineLabelsOn = False

  ; lat/lon limits
  dimx = dimsizes(x)
  numLat = dimx(0)
  numLon = dimx(1)
  minLat = x&$x!0$(0)
  maxLat = x&$x!0$(numLat-1)
  minLon = x&$x!1$(0)
  maxLon = x&$x!1$(numLon-1)
  aspect = (maxLon-minLon) / (maxLat-minLat)

  ; initialize plotting resources
  res = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@mpFillOn = False
  res@mpPerimDrawOrder = "PostDraw"

  res@mpFillOn = True
  res@mpOceanFillColor       = "lightblue"
  res@mpInlandWaterFillColor = "lightblue"
  res@mpLandFillColor        = 0

  res@mpGridAndLimbOn = True
  res@mpGridMaskMode  = "MaskNotOcean"
  res@mpGridLineDashPattern  = 2  ; dotted

  if (numLon*(maxLon - minLon)/(numLon-1) + 0.0001 .gt. 360.) then
    res@gsnAddCyclic = True
    res@mpOutlineBoundarySets = "Geophysical"
  else
    res@gsnAddCyclic = False
    res@mpOutlineBoundarySets = "National"
    res@mpOutlineSpecifiers = "AllUSStates"
  end if

  ; polar stereographic?
  do_polar = False
  if (res@gsnAddCyclic .and. maxLat-minLat .le. 90.) then
    if (maxLat + (maxLat - minLat)/(numLat-1) .ge. 90.) then
      do_polar = True
      res@gsnPolar = "NH"
      res@mpCenterLonF = 270
    else if (minLat - (maxLat - minLat)/(numLat-1) .le. -90.) then
      do_polar = True
      res@gsnPolar = "SH"
      res@mpCenterLonF = 90
    end if
    end if
  end if

  if (do_polar) then
    aspect = 1.2
    res@gsnRightStringParallelPosF = 1.300
    res@gsnLeftStringParallelPosF = -0.300
    res@gsnRightStringOrthogonalPosF = -0.030
    res@gsnLeftStringOrthogonalPosF = -0.030
  end if

  ; global lat/lon plot start at 20W
  if (res@gsnAddCyclic .and. .not.do_polar) then
    minLon = -20
    maxLon = 340
  end if

  ; contour fill mode for high resolution grids
  if (numLat*numLon .gt. 800000) then
    print ("Using raster fill mode (npts="+numLat*numLon+")")
    res@cnFillMode = "RasterFill"
    res@cnRasterMinCellSizeF = .0005
    res@cnRasterSmoothingOn = True
  end if

  ; set up plot shape/aspect
  if (aspect .gt. 1.1 .and. aspect .lt. 4.00) then
    res@vpWidthF      = 1.00
    res@vpHeightF     = 0.52
  else if (aspect .le. 1.1) then
    res@mpShapeMode   = "FreeAspect"
    res@vpWidthF      = 0.60
    res@vpHeightF     = 0.52
    res@gsnRightStringParallelPosF = 1.100
  else
    res@mpShapeMode   = "FreeAspect"
    res@vpWidthF      = 1.00
    res@vpHeightF     = 0.25
  end if
  end if

  res@mpMinLatF = minLat
  res@mpMaxLatF = maxLat
  if (.not. do_polar) then
    res@mpMinLonF = minLon   ; -20
    res@mpMaxLonF = maxLon   ; 340
    res@mpCenterLonF = (minLon + maxLon) * 0.5  ; 160
  end if

  if (lineLabelsOn) then
    res@cnLinesOn      = True
    res@cnLineLabelsOn = True                      ; turn on line labels
    res@cnLineLabelPlacementMode = "Randomized"    ; constant,randomized, or computed
    res@cnLineLabelDensityF = 1.00                 ; normal = 1.0
    res@cnLineLabelBackgroundColor = -1
    res@cnLineLabelFontHeightF = .008
    res@cnLabelMasking = True
    res@gsnContourNegLineDashPattern = 5
    res@cnLineLabelFont = "simplex_roman"
  else
    res@cnLinesOn      = False
    res@cnLineLabelsOn = False
  end if

  ; labelbar resources (only if no contour line labels)
  if (.not.lineLabelsOn) then
    res@lbLabelBarOn = True
    res@lbAutoManage = False
    res@lbOrientation  = "Vertical"
    res@lbRightMarginF = .00
    res@lbLabelFontHeightF = 0.013
    res@lbLabelAutoStride    = True ; optimal labels
  else
    res@lbLabelBarOn = False
  end if

  ; plot title resources
  res@gsnLeftStringFontHeightF = .020  ; .015
  res@gsnRightStringFontHeightF = .020  ; .015
  res@gsnLeftString     = ""
  res@gsnRightString    = ""
  res@tiMainString = " "
  res@tiMainOn = False

  return res
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; compute area weights using latitude

undef("area_weight")
function area_weight (fi:file,x:numeric)
local dimo,dimw,kind,dtr,bnd,latb,wt
begin
  
  dimx = dimsizes(x)
  kind = dimsizes(dimx)
  dtr = acos(0.0)/90.
  dimo = (/dimx(kind-2),dimx(kind-1)/)

  ; check for latitude bounds (will be used for grid cell area if it exists)
  if (isatt(x&$x!(kind-2)$,"bounds")) then
    bnds = x&$x!(kind-2)$@bounds
    latb = fi->$bnds$
    wt = conform_dims(dimo, tofloat(sin(latb(:,1)*dtr)-sin(latb(:,0)*dtr)), 0)
    print ("  using lat bounds variable ("+bnds+") for weight")
  else
    wt = conform_dims(dimo, tofloat(cos(x&$x!(kind-2)$*dtr)), 0)
   ;print ("  using coslat for weight")
  end if
 ;copy_VarCoords(x,wt)
  wt!0 = x!(kind-2)
  wt!1 = x!(kind-1)
  wt&$wt!0$ = x&$x!(kind-2)$
  wt&$wt!1$ = x&$x!(kind-1)$

  return wt
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("compute_statistics")
function compute_statistics (x1[*][*]:numeric, x2[*][*]:numeric, wt[*][*]:numeric)
local stats
begin

  stats = True
  stats@mean1 = wgt_areaave2 (x1,wt,0)
  stats@mean2 = wgt_areaave2 (x2,wt,0)
  stats@rmse  = wgt_arearmse2(x1,x2,wt,0)
  stats@corr  = pattern_cor  (x1,x2,wt,0)  ; centered correlation (matches Kushner plots)
  
  return stats
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;     start of main script     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

begin

  ; input variables and error checking
  error = 0

  if (.not.isvar("mfile")) then
    print("ERROR: no model file (mfile) specified")
    error = 1
  end if
  if (.not.isvar("ofile")) then
    print("ERROR: no obs file (ofile) specified")
    error = 1
  end if
  if (.not.isvar("var")) then
    print("ERROR: no variable name (var) specified")
    error = 1
  end if
  if (error .eq. 1) then
    exit
  end if

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;      optional arguments      ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;  mlab = model descriptor/label
  ;  olab = obs descriptor/label
  ;  mdate = model date label
  ;  odate = obs date label
  ;  rfile = region definition file
  ;  regrid = regridding type ('conserve', 'bilinear', or default 'best')
  ;  stats  = statistics type plotted, comma separated list ('mean', 'stdev', 'trend', or 'all')
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if (.not.isvar("rfile")) then
    regionFile = "$BW_PACKAGE_ROOT/ncl/regions.txt"
  else
    if (fileexists(rfile)) then  ; ncarg/6.2.1 req.
      regionFile = rfile
    else
      print("ERROR: region file does not exist or is not readable")
      status_exit(1)
    end if
  end if

  if (.not.isvar("stats")) then
    stats = "mean"
  end if

  ; generate statistics name list
  ; order MUST be mean, stdev, trend
  slist = str_lower(str_split(stats,","))
  if (any(slist .eq. "all")) then
    stats_names = (/"Mean","StDev","Trend"/)
  else
    stats_names = slist
    nstat = 0
    if (any(slist .eq. "mean")) then
      stats_names(nstat) = "Mean"
      nstat = nstat + 1
    end if
    if (any(slist .eq. "stdev")) then
      stats_names(nstat) = "StDev"
      nstat = nstat + 1
    end if
    if (any(slist .eq. "trend")) then
      stats_names(nstat) = "Trend"
      nstat = nstat + 1
    end if
  end if

; all titles save as attributes for easy recall
  Titles = True
  Titles@pr = "Precipitation (mm/dy)"
  Titles@tas = "Near Surface Temperature (C)"
  Titles@tasmean = "Near Surface Mean Temperature (C)"
  Titles@tasdiff = "Near Surface Diurnal Temperature Range (C)"
  Titles@tasmax  = "Near Surface Maximum Temperature (C)"
  Titles@tasmin  = "Near Surface Minimum Temperature (C)"
  Titles@wet = "Frequency of Wet Days (%)"
  Titles@clt = "Total Cloud Amount (%)"

  if (isatt(Titles,var)) then
     varTitle = Titles@$var$
  else
     varTitle = var
  end if

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; open monthly time series files
  fmod = addfile(mfile,"r")
  fobs = addfile(ofile,"r")

  ;;;;;;;; model ;;;;;;;;;
  print ("Read model data")
  moddata = fmod->$var$
  modmask = fmod->sftlf ; fraction
 ; change_units
  if (var .eq. "pr") then
    moddata = moddata*86400.
  else if (any(var .eq. (/"tas","tasmax","tasmin","tasmean"/))) then
    moddata = moddata-273.16
  else if (var .eq. "wet") then
    moddata = moddata*100.
  end if
  end if
  end if

; land/sea mask must be fraction (horiz_interp.ncl requirement)
  if (max(modmask) .gt. 75.0) then
    print("NOTE: converting model land/sea mask from percent to fraction")
    modmask = modmask*0.01
  end if

  ;;;;;;;; observation ;;;;;;;;;
  obsdata = fobs->$var$
  obsarea = area_weight(fobs,obsdata)

  obslon = obsdata&$obsdata!2$
  obslat = obsdata&$obsdata!1$
  ; fix for cru data
  if (obslon(0) .lt. -170.) then
    obsdata = lonFlip(obsdata)
    obsarea = lonFlip(obsarea)
  end if
  global = True
  ; fix for prism data
  if (obslon(0) .lt. 0. .and. obslon(dimsizes(obslon)-1).lt.0.) then
    obsdata&$obsdata!2$ = obsdata&$obsdata!2$ + 360.
    obsarea&$obsarea!1$ = obsarea&$obsarea!1$ + 360.
    global = False
  end if
 ;printVarSummary(obsdata)
 ;change_units(var,obsdata)
  if (var .eq. "pr") then
    obsdata = obsdata*86400.
  else if (any(var .eq. (/"tas","tasmax","tasmin","tasmean"/))) then
    obsdata = obsdata-273.16
  else if (var .eq. "wet") then
    obsdata = obsdata*100.
  end if
  end if
  end if

  ; compute climatology by season (annual,djf,mam,jja,son)
  ; for different statistics (mean, stdev, trend)
  ; do this on the native grid
  ; output dimensions: (season,stat,lat,lon)
  opt = True
  opt@stats = stats_names
  obsclim = compute_climatology_statistics(obsdata,opt)
  mclim   = compute_climatology_statistics(moddata,opt)
  ; set mostly ocean values to missing (avoids contamination of land-average from ocean)
  critmask = 0.67
  print("NOTE: masking model grid boxes with land fraction < "+toint(critmask*100.)+"%")
  modmask = where( modmask .lt. critmask, 0.0, modmask )
; mclim = where( conform(mclim,modmask,(/2,3/)) .lt. critmask, mclim@_FillValue, mclim )

  delete([/obsdata,moddata/])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;; interpolate model data to obs grid ;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  opt = True
  opt@critpc = 10.0
  if (isvar("regrid")) then
    if (regrid .eq. "bilinear") then
      opt@bilinearOn = True
    else if (regrid .eq. "conserve") then
      opt@bilinearOn = False
    end if
    end if
  end if

  modclim = obsclim
  modclim = interp_Wrap ( mclim, modmask, modclim, opt )

  if (interp_type_is_bilinear( mclim, modclim, opt )) then
     info_label = "Bilinear regrid to observation grid"
  else
     info_label = "Conservative regrid to observation grid"
  end if

  ; info about wet day cutoff (from comment att),  add to regrid label
  if (var .eq. "wet" .and. isfilevaratt(fmod,var,"comment")) then
     comments = str_split(fmod->$var$@comment," ")
     info_label = info_label + ", wet day cutoff = "+comments(str_match_ind(comments,"ge")+1)+" mm/dy"
    ;info_label = info_label + ", wet day cutoff = "+str_match_regex(str_split(fmod->$var$@comment," "),"[0-9]")
  end if

 ;printVarSummary(obsclim)
 ;printVarSummary(modclim)
  delete(mclim)

  ; interp model mask to obs grid
; mland = linint2( modmask&$modmask!1$, modmask&$modmask!0$, modmask, True, obsclim&$obsclim!2$, obsclim&$obsclim!1$, 0 )

  ; seasonal averages
  season_names = (/"Annual","DJF","MAM","JJA","SON"/)
 ;modvar = monthly_clim_to_seasonal(mvar,opt)
 ;obsvar = monthly_clim_to_seasonal(obsclim,opt)

  ; apply masking
; common_mask = where( mland .le. 0.33, True, False )
  common_mask = where(ismissing(modclim(0,0,:,:)), True, False )
  common_mask = where(ismissing(obsclim(0,0,:,:)), True, common_mask )
  do stat = 0, nstat-1
  do sea = 0, 4
    modclim(stat,sea,:,:) = where( common_mask, modclim@_FillValue, modclim(stat,sea,:,:) )
    obsclim(stat,sea,:,:) = where( common_mask, obsclim@_FillValue, obsclim(stat,sea,:,:) )
  end do
  end do

  ; set up geographic regions
  REGIONS = define_regions(regionFile)
  if (REGIONS@numReg .eq. 0) then
    print("ERROR: number of regions is zero")
    status_exit(1)
  end if

; grid spacing hack   GBL NAM SAM EUR EUA ASA SEA AUS AFR
  gridLatSpacing = (/  30, 20, 30, 10, 20, 30, 10, 10, 30/)
  gridLonSpacing = (/  30, 30, 20, 20, 30, 30, 20, 30, 30/)

; initialize statistics output file
  if (isvar("statsfile") .and. any(stats_names .eq. "Mean")) then
    quote = str_get_dq()  ; will be used later
    if (.not.fileexists(statsfile)) then
      if (isvar("mlab")) then
        text = "model="+mlab
        if (isvar("mdate")) then
          text = text + " ("+mdate+")"
        end if
        system("perl -e 'print "+quote+text+"\n"+quote+"'>"+statsfile)
      end if
      if (isvar("olab")) then
        text = "obs="+olab
        if (isvar("odate")) then
          text = text + " ("+odate+")"
        end if
        system("perl -e 'print "+quote+text+"\n"+quote+"'>>"+statsfile)
      end if
    end if
  end if

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; loop thru regions ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;

  do reg = 0, REGIONS@numReg-1

    ; define the geographic region (lon/lat)
    regName  = REGIONS@id(reg)
    regTitle = REGIONS@title(reg)
    minLon   = REGIONS@lonbeg(reg)
    maxLon   = REGIONS@lonend(reg)
    minLat   = REGIONS@latbeg(reg)
    maxLat   = REGIONS@latend(reg)
   ;if (any(regName .eq. (/"ASA","GBL","EUA"/))) then
   ;   continue
   ;end if
    print ("regName: "+regName+"; lon: "+minLon+" to "+maxLon+"; lat: "+minLat+" to "+maxLat)

    ; shift data longitude for this region 
    if (global) then
      modclim = my_lonPivot(modclim,minLon)
      obsclim = my_lonPivot(obsclim,minLon)
      obsarea = my_lonPivot(obsarea,minLon)
    end if

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;; loop thru statistics ;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

     do stat = 0, nstat-1

      ; open output file for postscript plots
      wks = gsn_open_wks ("ps", var+"."+str_lower(stats_names(stat))+"."+regName+".ps")

      res = initialize_resources(modclim(0,0,{minLat:maxLat},{minLon:maxLon}))
      res@cnFillOn        = True  ; turn on color fill
      res@cnInfoLabelOn   = False
      res@cnLevelSelectionMode = "ExplicitLevels"
      res@gsnRightStringFontHeightF = .017

      res@mpGridLatSpacingF = gridLatSpacing(reg)
      res@mpGridLonSpacingF = gridLonSpacing(reg)

      ;;;;;;;;;;;;;;;;;;;;;;;;;
      ;;; loop thru seasons ;;;
      ;;;;;;;;;;;;;;;;;;;;;;;;;

      do sea = 0, 4
        plot = new(3,graphic)

        ;;;;;;;;;;;;;;;;;;;;
        ; global statistics
        ;;;;;;;;;;;;;;;;;;;;
        statistics = compute_statistics(modclim(stat,sea,{minLat:maxLat},{minLon:maxLon}), \
                                        obsclim(stat,sea,{minLat:maxLat},{minLon:maxLon}), \
                                        obsarea({minLat:maxLat},{minLon:maxLon}))
        avgm = str_squeeze(sprintf("%7.3f",statistics@mean1))
        avgo = str_squeeze(sprintf("%7.3f",statistics@mean2))
        bias = str_squeeze(sprintf("%7.3f",statistics@mean1-statistics@mean2))
        rmse = str_squeeze(sprintf("%7.3f",statistics@rmse))
        corr = str_squeeze(sprintf("%6.4f",statistics@corr))
        print("reg = "+regName+"   stat = "+stats_names(stat)+"   ssn = "+sea+"   avgm = "+avgm)
        print("reg = "+regName+"   stat = "+stats_names(stat)+"   ssn = "+sea+"   avgo = "+avgo)
        print("reg = "+regName+"   stat = "+stats_names(stat)+"   ssn = "+sea+"   bias = "+bias)
        print("reg = "+regName+"   stat = "+stats_names(stat)+"   ssn = "+sea+"   rmse = "+rmse)
        print("reg = "+regName+"   stat = "+stats_names(stat)+"   ssn = "+sea+"   corr = "+corr)

        ; statistics file (can not do for range of levels)
        if (isvar("statsfile") .and. stats_names(stat) .eq. "Mean") then 
          ssn_name = str_lower(str_get_cols(season_names(sea),0,2))
          ascii_record = str_join((/var,regName,ssn_name,avgm,bias,rmse,corr/),"; ")
          system("perl -e 'print "+quote+ascii_record+"\n"+quote+"'>>"+statsfile)
        end if

        ; colors and palette for upper and middle plots
       ;res@cnFillPalette = "matlab_jet"
        res@cnLevels = contour_levels( var, stats_names(stat), "field" )
        if ( abs(res@cnLevels(0)+res@cnLevels(dimsizes(res@cnLevels)-1)) .lt. 1.e-6 ) then
          res@cnFillPalette = "BW_BlueWhiteOrangeRed"
        else
          res@cnFillPalette = "BW_precip3_16lev"
        end if

        ;;;;;;;;;;;;;;;
        ; upper frame
        ;;;;;;;;;;;;;;;
        if (isvar("mlab")) then
          res@gsnLeftString  = "MODEL: "+mlab
        else
          res@gsnLeftString  = "MODEL"
        end if
        if (isvar("mdate")) then
          res@gsnRightString = "a="+avgm+" ("+mdate+")"
        else
          res@gsnRightString = "a="+avgm
        end if
        if (isatt(res,"gsnPolar")) then
          plot(0) = gsn_csm_contour_map_polar(wks, modclim(stat,sea,{minLat:maxLat},{minLon:maxLon}), res )
        else
          ; hokie fix for long mlab string and narrow plots
         ;print("vpWidthF: "+res@vpWidthF)
          if (res@vpWidthF .lt. 0.75 .and. strlen(res@gsnLeftString) .gt. 30) then
            res@gsnLeftStringParallelPosF = (30-strlen(res@gsnLeftString))*.018
            res@gsnRightStringParallelPosF = 1.1 - (30-strlen(res@gsnLeftString))*.015
            print("NOTE: adjusted model/date labels x-postion by "+res@gsnLeftStringParallelPosF)
          end if
          plot(0) = gsn_csm_contour_map_ce(wks, modclim(stat,sea,{minLat:maxLat},{minLon:maxLon}), res )
          if (res@vpWidthF .lt. 0.75 .and. strlen(res@gsnLeftString) .gt. 30) then
            res@gsnLeftStringParallelPosF = 0.0 ; default
            res@gsnRightStringParallelPosF = 1.1 ; default
          end if
        end if

        ;;;;;;;;;;;;;;;;
        ; middle frame
        ;;;;;;;;;;;;;;;;
        if (isvar("olab")) then
          res@gsnLeftString  = "OBS: "+olab
        else
          res@gsnLeftString  = "OBS"
        end if
        if (isvar("odate")) then
          res@gsnRightString = "a="+avgo+" ("+odate+")"
        else
          res@gsnRightString = "a="+avgo
        end if
        if (isatt(res,"gsnPolar")) then
          plot(1) = gsn_csm_contour_map_polar(wks, obsclim(stat,sea,{minLat:maxLat},{minLon:maxLon}), res )
        else
          plot(1) = gsn_csm_contour_map_ce(wks, obsclim(stat,sea,{minLat:maxLat},{minLon:maxLon}), res)
        end if
        delete(res@cnLevels)

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ; difference (model-obs)
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;
        diff = modclim(stat,sea,:,:) - obsclim(stat,sea,:,:)
        diff!0 = modclim!2
        diff&$diff!0$ = modclim&$modclim!2$
        diff!1 = modclim!3
        diff&$diff!1$ = modclim&$modclim!3$
        ; labels
        res@gsnLeftString  = "MODEL - OBS"
        rightString = "~F33~D~F21~="+str_squeeze(bias)
        res@gsnRightString = str_concat((/rightString, \
                            "   r="+str_squeeze(corr)+ \
                            "   E="+str_squeeze(rmse)/))

        ; colors and palette for difference plot
        res@cnLevels = contour_levels( var, stats_names(stat), "diff" )
        res@cnFillPalette = "BW_BlueWhiteOrangeRed"

        if (isatt(res,"gsnPolar")) then
          res@gsnLeftStringParallelPosF = -0.200 + .016*(12 - strlen(res@gsnLeftString))
          res@gsnRightStringParallelPosF = 1.350
          plot(2) = gsn_csm_contour_map_polar(wks, diff({minLat:maxLat},{minLon:maxLon}), res)
        else
          plot(2) = gsn_csm_contour_map_ce(wks, diff({minLat:maxLat},{minLon:maxLon}), res)
        end if
        delete(res@cnLevels)

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ; plot the panel
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;
        resP                  = True
        resP@gsnPaperOrientation = "portrait"
        resP@gsnPanelLabelBar = False                   ; add common colorbar
        resP@gsnMaximize      = True
        resP@gsnFrame   = False
        resP@gsnDraw    = True
        if (isatt(res,"gsnPolar")) then
          resP@gsnPanelYWhiteSpacePercent = 1.0
        else
          resP@gsnPanelYWhiteSpacePercent = 2.0
        end if

       ;resP@gsnPanelDebug = True
        resP@gsnPanelBottom = .025
        resP@gsnPanelXWhiteSpacePercent = 0.0
        resP@txFontColor = 0
        resP@txString = "Null" ; trick to allow room for titles
        gsn_panel(wks,plot,(/3,1/),resP)

        ; title (and sub-title)
        txres = True
        txres@txFontHeightF = 0.022
        txYF = 0.97
        txXF = 0.49
        gsn_text_ndc (wks, season_names(sea)+" "+varTitle+", "+stats_names(stat), txXF, txYF, txres)

        ; informational title (interpolation type and "to" grid, also wet day cutoff)
        txres = True
        txres@txFontHeightF = 0.008
        txres@txFontColor = "grey75"
        txres@txJust = "BottomCenter"
        txres@txFontThicknessF = 0.5 ; probably does nothing 
        txYF = 0.001
        txXF = 0.49
        gsn_text_ndc (wks, info_label, txXF, txYF, txres)
        frame(wks)
     
        delete([/plot,resP,txres/])
   
      end do ; sea
      delete([/wks,res/])
    end do ; stat
  end do ; reg

  
 end


