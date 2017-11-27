
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$BW_PACKAGE_ROOT/ncl/compute_quartiles.ncl"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; shift longitude for regional average

undef("regAve_lonPivot")
function regAve_lonPivot (x:numeric, pLon:numeric)
local dimx,rank,lonName,temp,xlon,dlon,indP,mlon,indL,n,tlon,indt
begin

  dimx = dimsizes(x)
  rank = dimsizes(dimx)
  if (rank .lt. 2 .or. rank .gt. 3) then
     print ("ERROR: regAve_lonPivot: rank must be 2 or 3")
     exit
  end if
  lonName = x!(rank-1)
  temp = x
  xlon = x&$lonName$
  dlon = xlon-pLon
  dlon = where(dlon.lt.-360.,dlon+360.,dlon)
  dlon = where(dlon.gt. 360.,dlon-360.,dlon)
  indP = minind(where(dlon.ge.0,dlon,361.))

  mlon = dimx(rank-1)
  indL = mlon-1
  n    = indL-indP

  if (rank .eq. 2) then
     temp(:,0:n) = (/ x(:,indP:indL) /)
     if (indP.gt.0) then
        temp(:,n+1:) = (/ x(:,0:indP-1) /)
     end if
  end if
  if (rank .eq. 3) then
     temp(:,:,0:n) = (/ x(:,:,indP:indL) /)
     if (indP.gt.0) then
        temp(:,:,n+1:) = (/ x(:,:,0:indP-1) /)
     end if
  end if

  tlon = new(mlon,typeof(xlon))
  tlon(0:n) = (/ xlon(indP:indL) /)
  if (indP.gt.0) then
     tlon(n+1:) = (/ xlon(0:indP-1) /)
  end if
  delete (tlon@_FillValue)

  if (tlon(0).ge.0.) then  ; (say) 20,25,...,350,355,0,5,..
      indt = ind(tlon.lt.tlon(0))
      if (.not.all(ismissing(indt))) then
          tlon(indt) = (/ tlon(indt) + 360. /)
      end if
  end if
  if (tlon(0).ge.180. .or. tlon(0).eq.360.) then
      tlon = (/ tlon -360. /)
  end if
  tlon = where(tlon.lt.tlon(0),tlon+360.,tlon)
  copy_VarAtts (xlon,tlon)

  temp&$lonName$ = tlon
  return temp
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  compute grid box weights based on latitude bounds if available 
;  otherwise based on cosine latitude.
;  returns a 1d weight

undef("grid_box_weights")
function grid_box_weights(fi:file,x[*][*][*]:numeric)
local name,dtr,bnds,wt
begin
  dtr = atan(1.)/45.
  if (isatt(x&$x!1$,"bounds")) then
    name = x&$x!1$@bounds
    bnds = fi->$name$ * dtr
    wt = tofloat(sin(bnds(:,1)*dtr)-sin(bnds(:,0)*dtr))
  else
    wt = tofloat(cos(x&$x!1$*dtr))
  end if
  return wt
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;  compute the zonal or meridional average for a region  ;;;;;;;;
; numeric mask (land vs. sea) is applied to the data before the average
; the average is computed where >1% of the data is available
; if dim = 1 then averaging is along the latitude axis
; if dim = 2 then averaging is along the longitude axis (zonal average)

undef("regional_average")
function regional_average (fi:file,x[*][*][*]:numeric, xm[*][*]:numeric, y0[1]:numeric, y1[1]:numeric, x0[1]:numeric, x1[1]:numeric, dim[1]:integer, opt:logical )
local xt,mt,dtr,name,bnds,xs,ms,xa,mm,ml,m1,m1s
begin
   xt = x  ; temporary of x
   mt = xm
   if (dim .eq. 1) then
     wt = grid_box_weights(fi,x) ; apply grid box weights
     mt = mt*conform(mt,wt,0)
   end if
   xt = xt * conform(x,mt,(/1,2/))  ; apply mask
   xt = regAve_lonPivot(xt,x0)     ; pivot on longitude
   mt = regAve_lonPivot(mt,x0)

   xs = dim_sum_n_Wrap(xt(:,{y0:y1},{x0:x1}),dim) ; summation across latitude or longitude
   ms = dim_sum_n_Wrap(mt  ({y0:y1},{x0:x1}),dim-1)
   m1 = mt
   if (dim .eq. 1) then
     m1 = conform(mt,wt,0)
   else
     m1 = 1.0
   end if
   m1s = dim_sum_n_Wrap(m1({y0:y1},{x0:x1}),dim-1) * 0.01 ; 1% of the data must be present
   ml = where(ms .lt. m1s, False, True)
   ms = where(ms .lt. m1s, 1.0, ms)
   xa = xs
   xa = where(conform(xa,ml,1), xs/conform(xa,ms,1), xa@_FillValue)
   ; optional printout for debugging
   if (opt) then
     printVarSummary(xa)
     printVarSummary(ml)
     print("x  (mn/mx) = "+min(x)+", "+max(x))
     print("xs (mn/mx) = "+min(xs)+", "+max(xs))
     print("wt (mn/mx) = "+min(wt)+", "+max(wt))
     print("mt (mn/mx) = "+min(mt)+", "+max(mt))
     print("ms (mn/mx) = "+min(ms)+", "+max(ms))
     print("xa (mn/mx) = "+min(xa)+", "+max(xa))
     exit
   end if
   return xa
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Compute weights for a monthly climatology given the starting
; and ending years. A standard calendar is assumed.
; Note: -12 < shift < 12
; For weights from Dec-Nov, shift = -1
; For weights from Feb-Jan, shift =  1

undef("monthly_clim_weights") 
function monthly_clim_weights (yrbeg[1]:integer, yrend[1]:integer, Opt:logical)
local  fdays,shft,wt,mo,yr,nyrs,nleap
begin
  ; non-leap year weights (Jan-Dec)
  fdays = (/31.,28.25,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./)
  ; shifted weights ?
  shft = 0
  if (Opt) then
    if (isatt(Opt,"shift")) then
      shft = Opt@shift
    end if
  end if
  if (shft .lt. 0) then
    shft = shft + 12
  end if
  if (shft .lt. 0 .or. shft .gt. 11) then
    print("ERROR: in monthly_clim_weights: incorrect month shift")
    exit
  end if

  wt = new(12,float)
  do mo = 0, 11
    wt(mo) = fdays(mod(mo+shft,12))
  end do
  ; compute correct weight for Feb (assuming Gregorian/Julian calendar)
  if (yrbeg .gt. 0 .and. yrend .gt. 0) then
    nyrs = (yrend-yrbeg)+1
    if (nyrs .lt. 0) then
      print("ERROR: invalid start and end year in monthly_clim_to_seasonal")
      exit
    end if
    nleap = 0
    do yr = yrbeg, yrend
      if (mod(yr,4) .eq. 0) then
        nleap = nleap+1
      end if
    end do
    wt(mod(13-shft,12)) = 28. + (tofloat(nleap)/tofloat(nyrs))
  end if
  return wt
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; convert a monthly climatology (J,F,M,A,...), i.e., 12 time levels
; to a seasonal climatology (Ann,DJF,MAM,JJA,SON), i.e., 5 time levels

undef("monthly_clim_to_seasonal")
function monthly_clim_to_seasonal (x[12][*][*]:numeric,opt:logical)
local dimx,wm,xm,xs,k,ks,ke,nyrs,nleap,yr
begin

  ; note month wts are shifted (Dec-Nov)
  res = True
  res@shift = -1
  if (isatt(opt,"yr1") .and. isatt(opt,"yr2")) then
    wm = monthly_clim_weights( opt@yr1, opt@yr2, res )
  else
    wm = monthly_clim_weights( 0, 0, res )
  end if

  dimx = dimsizes(x)
  xs = new((/5,dimx(1),dimx(2)/),typeof(x))
  xs!0 = "season"
  xs&season = ispan(0,4,1)
  xs&season@long_name = "Season number (0=ann,1=djf,2=mam,3=jja,4=son)"
  xs&season@units = "None"
  ; shift monthly average (dec-nov)
  xm = x
  xm(0,:,:) = x(11,:,:)
  do k = 1, 11
     xm(k,:,:) = x(k-1,:,:)
  end do
  ; annual average
  xs(0,:,:) = dim_avg_wgt_n_Wrap(x,wm,0,0)
  ; season averages (djf=0:2,mam=3:5,etc.)
  do k = 0, 3
     ks = 3*k
     ke = ks+2
     xs(k+1,:,:) = dim_avg_wgt_n_Wrap(xm(ks:ke,:,:),wm(ks:ke),0,0)
  end do
  return xs
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
undef("simple_legend")
function simple_legend (wks:graphic,plot:graphic,labels[*]:string,x[1]:numeric,y[1]:numeric,res)
local dummy,nlab,str,lineLength,fontHeight,fontOffset,lineHeight,n,resl,rest,yn,xl,yl,xt,yt
begin
  dummy = True ; variable to store polyline overlays
  nlab = dimsizes(labels)

  ; resources
  lineLength = .050
  fontHeight = .020
  if (isatt(res,"lineLength")) then
     lineLength = res@lineLength
  end if
  if (isatt(res,"fontHeight")) then
     fontHeight = res@fontHeight
  end if
  if (isatt(res,"fontOffset")) then
     fontOffset = res@fontOffset
  else
     fontOffset = fontHeight*0.60
  end if
  if (isatt(res,"lineHeight")) then
     lineHeight = res@lineHeight * fontHeight
  else
     lineHeight = 1.60 * fontHeight
  end if

  xl = new(2,typeof(x))
  yl = new(2,typeof(y))
  xt = new(1,typeof(x))
  yt = new(1,typeof(y))
  xndc = new(1,typeof(x))
  yndc = new(1,typeof(y))
  datatondc(plot, x, y, xndc, yndc)

  resl = True
  rest = True

  do n = 0, nlab-1
     yn = yndc - (tofloat(n)+0.5)*lineHeight

     ;;;;;;; plot the line ;;;;;;;;
     ndctodata(plot, (/xndc,xndc+lineLength/), (/yn,yn/), xl, yl)
     ; line resources
     if (isatt(res,"lineColors")) then
        resl@gsLineColor = res@lineColors(n)
     end if
     if (isatt(res,"lineThicknesses")) then
        resl@gsLineThicknessF = res@lineThicknesses(n)
     end if
     if (isatt(res,"lineDashPatterns")) then
        resl@gsLineDashPattern = res@lineDashPatterns(n)
     end if
     if (isatt(res,"lineDashSegLen")) then
        resl@gsLineDashSegLenF = res@lineDashSegLen
     end if

     str = unique_string("polyline")
     dummy@$str$ = gsn_add_polyline(wks, plot, xl, yl, resl)

     ;;;;;;; plot the text ;;;;;;;;
     ndctodata(plot, xndc+lineLength+fontOffset, yn, xt, yt)
     ; text resources
     rest@txJust = "CenterLeft"
     rest@txFontHeightF = fontHeight
     if (isatt(res,"fontColor")) then
        rest@gsFontColor = res@fontColor
     end if

     str = unique_string("text")
     dummy@$str$ = gsn_add_text(wks, plot, labels(n), xt, yt, rest)
  end do

  return dummy
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Compute anomaly of monthly climatology (on dimension 1)
;   input:  month x lat x lon
;
undef("anomaly_Wrap")
function anomaly_Wrap ( fi:file, x[12][*][*]:numeric, yr1:integer, yr2:integer )
local wt,tav,area,ave,xa
begin

  ; time average
  wt = monthly_clim_weights (yr1, yr2, False)  ; => wt(mon)
  tav = dim_avg_wgt_n_Wrap (x, wt, 8, 0)       ; => tav(lat,lon)

  area = grid_box_weights(fi,x) ; area(lat)
  ave = wgt_areaave2(tav, conform(tav,area,0), 0) ; => ave[1]

  xa = x
  xa = x - ave

  return xa
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;    Create a three panel plot for each season
;

undef("plot_curves")
procedure plot_curves ( wks:graphic, variable:string, REGIONS:logical, \
                       MODELDATA:logical, OBSDATA:logical, CMIP5DATA:logical, CMIP5GRID:logical, \
                       QUARTILES:logical, Options:logical)
local shadeExtremes, dashPattern, resf, seasonNames, axis_names, axis_units, plot, dummy, res, \
      funcCode, greekTau, subX, units, plot_title, scale_factor, ssn, reg, regionName, \
      trXMinF, trXMaxF, trYMinF, trYMaxF, ANOMALY, legend_labels, legend_colors, legend_linethickness, legend_dashes, \
      qdata, qdims, npoly, xpoly, ypoly, gridvalues, str, i, j, resL, attname, odata, ogrid, mdata, mgrid, ngrid, gridbnds, \
      ybeg, yend, xbeg, xend, xleg, yleg, lgres, resP
begin

  ; optional arguments
  shadeExtremes = True
  ANOMALY = False

  if (Options) then
    if (isatt(Options,"ANOMALY")) then
      ANOMALY = Options@ANOMALY
    end if
    if (isatt(Options,"shadeExtremes")) then
      shadeExtremes = Options@shadeExtremes
    end if
  end if

  if (shadeExtremes) then
     dashPattern = 0
  else
     dashPattern = 2
  end if

  ; fill resources
  resf = True
  resf@tfPolyDrawOrder = "Predraw"

  ; labels for plot
  seasonNames = (/"Annual", "DJF", "MAM", "JJA", "SON"/)

  axis_names = (/ "lon", "lat" /)
  axis_units = (/ "degrees_east", "degrees_north" /)

  ; define colors for future use
  gsn_define_colormap( wks, array_append_record((/"white","black","grey65","grey75","grey85",MODELDATA@color/),OBSDATA@colors,0) )

  plot = new(3,graphic)
  dummy = True ; variable to store polyline overlays

  ; general resources for xy plots
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@vpWidthF = 0.90
  res@vpHeightF = 0.45

  res@xyLineColors = (/"grey75","grey75"/)
  res@xyLineThicknesses = (/ 1., 1. /)
  res@xyDashPatterns = (/dashPattern,dashPattern/)
  res@xyLineDashSegLenF = .12

  ; labeling (y-axis and plot title)
  funcCode = "%"
  res@tiYAxisFuncCode = funcCode
  if (variable .eq. "tauu") then
    greekTau = funcCode+"F33"+funcCode+"t"+funcCode+"F21"+funcCode
    subX = funcCode+"B"+funcCode+"x"+funcCode+"N"+funcCode
    units = "Nt/m"+funcCode+"S"+funcCode+"2"+funcCode+"N"+funcCode
    res@tiYAxisString = greekTau+subX+" ("+units+")"
    plot_title = "Zonal Average Eastward Wind Stress ("+units+")" 
  else if (variable .eq. "curl") then
    greekTau = funcCode+"F33"+funcCode+"t"+funcCode+"F21"+funcCode
   ;subX = funcCode+"B"+funcCode+"x"+funcCode+"N"+funcCode
    units = "N/m"+funcCode+"S"+funcCode+"3"+funcCode+"N"+funcCode
    res@tiYAxisString = "curl("+greekTau+")*10"+funcCode+"S"+funcCode+"6"+funcCode+"N"+funcCode+" ("+units+")"
    plot_title = "Zonal Average Curl of Wind Stress ("+units+")" 
  else if (variable .eq. "pr") then
    units = "mm/d"
    res@tiYAxisString = units
    plot_title = "Precipitation Rate ("+units+")" 
    scale_factor = 86400.
  end if
  end if
  end if

;;;;;;;;;;;;;;;;;;;;;
; loop over seasons
;;;;;;;;;;;;;;;;;;;;;
  do ssn = 0, 4

;;;;;;;;;;;;;;;;;;;;;;;;;;
;    loop over region
;;;;;;;;;;;;;;;;;;;;;;;;;;
     do reg = 0, 2
        regionName = REGIONS@name(reg)
        print("season: "+seasonNames(ssn)+", region: "+REGIONS@longname(reg))

        ; same limits as Kushner plot (for Pacific/Atlantic)
        if (REGIONS@dimavg(reg) .eq. 1) then
          trXMinF = REGIONS@lonbeg(reg)
          trXMaxF = REGIONS@lonend(reg)
        else if (REGIONS@dimavg(reg) .eq. 2) then
          trXMinF = REGIONS@latbeg(reg)
          trXMaxF = REGIONS@latend(reg)
        end if
        end if

        ; same y-axis limits as Kushner plot
        if (variable .eq. "tauu") then
          if (reg .eq. 0 .or. reg .eq. 1) then
            trYMinF = -.14
            trYMaxF = .20
          else if (reg .eq. 2) then
            trYMinF = -.18
            trYMaxF = .34
          end if
          end if
        else if (variable .eq. "curl") then
          if (reg .eq. 0) then
            trYMinF = -.18
            trYMaxF = .20
          else if (reg .eq. 1) then
            trYMinF = -.34
            trYMaxF = .16
          else if (reg .eq. 2) then
            trYMinF = -.22
            trYMaxF = .22
          end if
          end if
          end if
        else if (variable .eq. "pr") then
          if (ANOMALY) then
            trYMinF = -3.0
            trYMaxF = 9.0
          else
            trYMinF = 0.0
            trYMaxF = 12.0
          end if
        end if
        end if
        end if

        res@trXMinF = trXMinF
        res@trXMaxF = trXMaxF
        res@trYMinF = trYMinF
        res@trYMaxF = trYMaxF

        ; legend labels/colors
        if (reg .lt. 2) then
           legend_labels = array_append_record(str_upper(OBSDATA@files),(/CMIP5DATA@descriptor,MODELDATA@descriptor/),0)
           legend_colors = array_append_record(str_upper(OBSDATA@colors),(/"grey65",MODELDATA@color/),0)
        else
           legend_labels = array_append_record(str_upper(OBSDATA@files (ind(OBSDATA@land))),(/CMIP5DATA@descriptor,MODELDATA@descriptor/),0)
           legend_colors = array_append_record(str_upper(OBSDATA@colors(ind(OBSDATA@land))),(/"grey65",MODELDATA@color/),0)
        end if
        legend_linethickness = new(dimsizes(legend_labels),float)
        legend_dashes = new(dimsizes(legend_labels),integer)
        legend_linethickness = 2.
        legend_dashes = 0
        legend_linethickness(dimsizes(legend_labels)-2) = 1.    ; cmip5
        legend_dashes(dimsizes(legend_labels)-2) = dashPattern  ; cmip5

        ; load quartile data for this region (and season)
        qdata = QUARTILES@$regionName$
        qdims = dimsizes(qdata) ; quartile x ssn x npts

        ; if any quartile is missing set them all missing
        ; only need to check this season
        do i = 0, qdims(2)-1
          if (any(ismissing(qdata(:,ssn,i)))) then
            qdata(:,ssn,i) = qdata@_FillValue
          end if
        end do

        ; scale data if necessary
        if (isvar("scale_factor")) then
          qdata = qdata*scale_factor
        end if

        ; polygon for shading
        npoly = qdims(2)
        xpoly = new(npoly*2,typeof(qdata))
        ypoly = new(npoly*2,typeof(qdata))

        xpoly(0:npoly-1) = CMIP5GRID@$regionName$
        xpoly(npoly:2*npoly-1:-1) = CMIP5GRID@$regionName$

        print("  plotting cmip5: "+CMIP5DATA@descriptor)
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ; plot min/max of cmip5 models
        res@gsnLeftStringFontHeightF = .023
        res@gsnRightStringFontHeightF = .023
        res@gsnLeftString = REGIONS@longname(reg)
        res@gsnRightString = seasonNames(ssn)
        gridvalues = CMIP5GRID@$regionName$ ; x-axis values
        gridvalues@units = axis_units(REGIONS@dimavg(reg)-1) ; add units to label axis
        plot(reg) = gsn_csm_xy(wks, gridvalues, qdata(0:4:4,ssn,:), res)

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ; add shading for extremes
        if (shadeExtremes) then
           ypoly(0:npoly-1) = qdata(0,ssn,:)
           ypoly(npoly:2*npoly-1:-1) = qdata(4,ssn,:)
           resf@gsFillColor = "grey85"
           str = unique_string("polyline")
           dummy@$str$ = gsn_add_polygon(wks, plot(reg), xpoly, ypoly, resf)
           if (reg .eq. 3) then ; turned off
             do i = 0, npoly-1
               str = sprinti("%3i",i)
               do j = 0, 4
                 if (j .ne. 2) then
                   if (.not. ismissing(qdata(j,ssn,i))) then
                     str = str + sprintf("%10.4f",qdata(j,ssn,i))
                   else
                     str = str + "    ------"
                   end if
                 end if
               end do
               print("i="+str)
             end do
           end if
        end if

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ; add shading of middle two quartiles for cmip5 models
        ypoly(0:npoly-1) = qdata(1,ssn,:)
        ypoly(npoly:2*npoly-1:-1) = qdata(3,ssn,:)
        resf@gsFillColor = "grey65"
        str = unique_string("polyline")
        dummy@$str$ = gsn_add_polygon(wks, plot(reg), xpoly, ypoly, resf)
        delete([/xpoly,ypoly/])

        ;;;;;;;;;;;;;;;;;;;;;;;;
        ;  add lines for obs
        ;;;;;;;;;;;;;;;;;;;;;;;;

        resL = True
        ; loop over obs data sets
        do i = 0, dimsizes(OBSDATA@files)-1
           ; region #2 is land - skip ocean only data sets
           if (reg .eq. 2 .and. .not.OBSDATA@land(i)) then
              continue
           end if
           print("  plotting obs: "+OBSDATA@files(i))
           attname = OBSDATA@files(i)+"_"+regionName
           odata = OBSDATA@$attname$  ; season x latitude
           ; scale data if necessary
           if (isvar("scale_factor")) then
             odata = odata*scale_factor
           end if
           attname = OBSDATA@files(i)+"_"+regionName+"_grid"
           ogrid = OBSDATA@$attname$  ; 1d grid
           resL@gsLineColor = OBSDATA@colors(i)
           resL@gsLineThicknessF = 2.0
           resL@tfPolyDrawOrder = "PostDraw"
           ; add line (need to add legend)
           str = unique_string("polyline")
           dummy@$str$ = gsn_add_polyline(wks, plot(reg), ogrid, odata(ssn,:), resL)
           delete([/odata,ogrid/])
        end do

        ;;;;;;;;;;;;;;;;;;;;;;;;
        ;  add line for model
        ;;;;;;;;;;;;;;;;;;;;;;;;
          print("  plotting model: "+MODELDATA@descriptor)
          resL@gsLineColor = MODELDATA@color
          resL@gsLineThicknessF = 2.0
          resL@tfPolyDrawOrder = "PostDraw"
          mdata = MODELDATA@$regionName$   ; season x latitude
          ; scale data if necessary
          if (isvar("scale_factor")) then
            mdata = mdata*scale_factor
          end if
          print("    mn,mx="+min(mdata(ssn,:))+","+max(mdata(ssn,:)))
          attname = regionName+"_grid"
          mgrid = MODELDATA@$attname$
          str = unique_string("polyline")
          dummy@$str$ = gsn_add_polyline(wks, plot(reg), mgrid, mdata(ssn,:), resL)
          delete([/mdata,mgrid/])

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ; add polygons for missing data
        ; use the cmip5 quartile data
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        ; define grid edges
        ngrid = dimsizes(gridvalues)
        gridbnds = new(ngrid+1,float)
        gridbnds(0) = 1.5*gridvalues(0) - 0.5*gridvalues(1)
        gridbnds(ngrid) = 1.5*gridvalues(ngrid-1) - 0.5*gridvalues(ngrid-2)
        gridbnds(1:ngrid-1) = 0.5*(gridvalues(0:ngrid-2)+gridvalues(1:ngrid-1))
        delete(gridvalues)

        resf = True
        resf@gsFillColor = "white"

        ybeg = trYMinF
        yend = trYMaxF
        xbeg = new(1,float)
        xend = new(1,float)
        do i = 0, ngrid-1
          if (ismissing(qdata(0,ssn,i)) .and. ismissing(xbeg)) then
            xbeg = gridbnds(i)
          else if (.not.ismissing(qdata(0,ssn,i)) .and. .not.ismissing(xbeg)) then
            xend = gridbnds(i+1)
            str = unique_string("polygon")
            dummy@$str$ = gsn_add_polygon(wks, plot(reg), (/xbeg,xend,xend,xbeg/), (/ybeg,ybeg,yend,yend/), resf)
            xbeg = xbeg@_FillValue
          end if
          end if
        end do
        delete([/qdata,gridbnds/])

        ; finally add reference lines
        resL = True
        resL@gsLineColor = "black"
        resL@gsLineThicknessF = 1.0
        resL@gsLineDashPattern = 2
        resL@gsLineDashSegLenF = 0.07
        resL@tfPolyDrawOrder = "PostDraw"
        if (trYMinF .lt. -0.0001 .and. trYMaxF.gt.0.0001) then
         ;resL@gsnYRefLine = 0.0  ; add a reference line at y=0.0
          str = unique_string("polyline")
          dummy@$str$ = gsn_add_polyline(wks, plot(reg), (/trXMinF,trXMaxF/), (/0.,0./), resL)
        end if
        if (trXMinF .lt. -0.1 .and. trXMaxF.gt.0.1) then
          str = unique_string("polyline")
          dummy@$str$ = gsn_add_polyline(wks, plot(reg), (/0.,0./), (/trYMinF,trYMaxF/), resL)
        else if (trXMinF .lt. 179. .and. trXMaxF.gt.181.) then
          str = unique_string("polyline")
          dummy@$str$ = gsn_add_polyline(wks, plot(reg), (/180.,180./), (/trYMinF,trYMaxF/), resL)
        else
          print ("trXMinF, trXMaxF = "+trXMinF+", "+trXMaxF)
        end if
        end if
        delete(resL)

        ;;;;;;;;;;;;;;;;;;;;;;;;
        ;  add/attach a legend
        ;;;;;;;;;;;;;;;;;;;;;;;;
        ; vary the starting point depending on variable
        ; may need to adjust this
        if (variable .eq. "tauu") then
          xleg = trXMinF + 0.31*(trXMaxF-trXMinF)
        else if (variable .eq. "curl") then
          xleg = trXMinF + 0.20*(trXMaxF-trXMinF)
        else
          xleg = trXMinF + 0.05*(trXMaxF-trXMinF)
        end if  
        end if  
        yleg = trYMinF + 0.92*(trYMaxF-trYMinF)
        lgres = True
        lgres@lineColors = legend_colors(::-1)
        lgres@lineThicknesses = legend_linethickness(::-1)
        lgres@lineDashPatterns = legend_dashes(::-1)
        lgres@lineDashSegLen = .050 ; not same as plot seg length? (.12)
        lgres@fontHeight = .015
        legend = simple_legend (wks,plot(reg),legend_labels(::-1),xleg,yleg,lgres)
        delete(lgres@lineColors)
        delete(lgres@lineThicknesses)
        delete(lgres@lineDashPatterns)
        delete([/legend_labels,legend_colors,legend_linethickness,legend_dashes/])

     end do

     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ; assemble the 3 plots into one panel
     ; add title to top of panel
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     resP                  = True
     resP@gsnPaperOrientation = "portrait"
     resP@gsnPanelLabelBar = False
     resP@gsnMaximize      = True
     resP@gsnPanelBottom   = 0.05
     resP@txFuncCode = funcCode
     resP@txString   = plot_title
     resP@txFontHeightF = .017
     resP@txPosXF = .53
     resP@gsnPanelYF = (/ .90, .59, .28/)
     gsn_panel(wks,plot,(/3,1/),resP)

  end do

end
