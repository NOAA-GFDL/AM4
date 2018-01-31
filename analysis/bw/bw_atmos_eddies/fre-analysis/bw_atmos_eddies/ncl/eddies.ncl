load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"
load "$BW_PACKAGE_ROOT/ncl/contributed.ncl"
load "$BW_PACKAGE_ROOT/ncl/lookup2.ncl"

external EDDIES "./eddies.so"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;  L O C A L    F U N C T I O N S  ;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef ("determine_start_end")
function determine_start_end (t:numeric,y1:integer,y2:integer)
local n, ny, y, yb, ye, ym1, ym2
begin
  y = cd_calendar(t,-2)
  ym1 = y1*10000+101
  ym2 = y2*10000+1231
  ny = dimsizes(y)-1
  yb = -1
  ye = -1
  do n=0,ny
     if (y(n) .eq. ym1) then
        yb = n
     end if
     if (y(n) .eq. ym2) then
        ye = n
     end if
  end do
 ;if (ye-yb .ne. 12*(y2-y1+1)-1) then
  if (yb .eq. -1 .or. ye .eq. -1) then
     print ("ERROR in determine_start_end: Start/end times could not be determined")
     print ("Dates: "+y(0)+" to "+y(ny))
     print ("ym1,ym2: "+ym1+", "+ym2)
     exit
  end if
  return (/yb,ye/)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("average_daily_to_monthly")
function average_daily_to_monthly (x:numeric, t:numeric, opt:logical)
local dimx,ks,ke,n,y,tn
begin
  dimx = toint(dimsizes(x))
  dimy = dimx
  dimy(0) = opt@nmonth
  y = new(dimy,typeof(x))
  y@_FillValue = x@_FillValue
  tf = new(opt@nmonth,float)
  mdays = new(opt@nmonth,integer)
  mdays = (opt@ke - opt@ks) + 1

  EDDIES::average3d(dimx(2),dimx(1),dimx(0),x,x@_FillValue, mdays,opt@nmonth, y)
  EDDIES::average1d(dimx(0),tofloat(t),mdays,opt@nmonth,tf)

  tn = new(opt@nmonth,typeof(t))
  tn = tf

  copy_VarAtts(t,tn)
  y!0 = t!0
  y&$y!0$ = tn
  return y
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("time_indexing_init")
function time_indexing_init(t:numeric)
local d,ntime,y1,y2,m1,m2,nmonth,ks,ke,opt,n,tn
begin
  opt = True
; pre-compute time axis indexing
  d = cd_calendar(t,-2)
  ntime = dimsizes(d)
  y1 = d(0)/10000
  y2 = d(ntime-1)/10000
  m1 = (d(0) - y1*10000)/100
  m2 = (d(ntime-1) - y2*10000)/100
  opt@nmonth = 12*(y2-y1) + (m2-m1) + 1
  print("y1,y2:"+y1+","+y2)
  print("m1,m2:"+m1+","+m2)
  print("nmonth:"+opt@nmonth)
  ks = new(opt@nmonth,integer)
  ke = new(opt@nmonth,integer)
  ks = ind(d-(d/100)*100 .eq. 1)
  ke(:opt@nmonth-2) = ks(1:)-1
  ke(opt@nmonth-1) = ntime-1

  ; new time coordinates
  tn = new(opt@nmonth,typeof(t))
  do n = 0, opt@nmonth-1
     tn(n) = dim_avg(t(ks(n):ke(n)))
  end do

  opt@ks = ks
  opt@ke = ke
  opt@tn = tn
  return opt
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("read_observed_data")
function read_observed_data (f:file,v:string,y1:integer,y2:integer)
local x,y,t,d,indx,dat
begin
  ; first check if latitude (dim 1) should be flipped
  d = getfilevardims(f,v)
  t = f->$d(0)$
  if (dimsizes(d) .eq. 3) then
    y = f->$d(1)$
    x = f->$d(2)$
  end if
  if (isatt(t,"calendar")) then
    t@calendar = str_lower(t@calendar) ; cd_calendar only recognizes lower case
  end if
  indx = determine_start_end(t,y1,y2)
  print ("obs time index: "+indx(0)+" to "+indx(1)+" (npts="+(indx(1)-indx(0)+1)+")")

  if (dimsizes(d) .eq. 3) then
     if (y(1) .gt. y(0)) then
        dat = f->$v$(indx(0):indx(1),:,:)
     else
        dat = f->$v$(indx(0):indx(1),::-1,:)
     end if
     ; check if longitude should be shifted
     if (x(0) .lt. -174.) then
        print ("Shifting longitude on observed data")
        dat = lonFlip(dat)
     end if
  else
     dat = f->$v$(indx(0):indx(1))
  end if
  return dat
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("var_letter_lookup")
function var_letter_lookup (a:string)
local b
begin
  if (any(a .eq. (/"u","v","t","q","h","w"/))) then
     b = a
; else if (a .eq. "w") then
;    b = "omg"
  else if (a .eq. "z") then
     b = "h"
  else
     print ("ERROR: invalid variable letter")
     exit
  end if
  end if
; end if

  return b
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("setup_color_maps")
function setup_color_maps(wks:graphic,var:string)
local name, rgb
begin
  name = str_get_cols(var,0,1)
  rgb = new(2,string)
  ; defaults
  rgb(0) = "local_default_grey"
  rgb(1) = "blue_grey_red"
  if (name .eq. "vq") then
    ;rgb(0) = "light_cyan_yellow"
     rgb(0) = "brown_grey_green"
     rgb(1) = "brown_grey_green"
  end if
  if (name .eq. "vt") then
    ;rgb(0) = "light_cyan_yellow"
    ;rgb(1) = "light_cyan_yellow"
     rgb(0) = "blue_grey_red"
     rgb(1) = "blue_grey_red"
  end if
  if (name .eq. "wq") then
     rgb(0) = "brown_grey_green"
     rgb(1) = "brown_grey_green"
  end if
  if (name .eq. "wt") then
     rgb(0) = "light_blue_red"
     rgb(1) = "light_blue_red"
  end if
  if (name .eq. "uv") then
    ;rgb(0) = "light_cyan_yellow"
    ;rgb(1) = "light_cyan_yellow"
     rgb(0) = "blue_grey_red"
     rgb(1) = "blue_grey_red"
  end if

  return rgb
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; compute statistics for model vs obs
; model mean, bias, rmse, corr

undef("compute_statistics")
function compute_statistics(mdat[*][*]:numeric,odat[*][*]:numeric,opt:logical)
local dtr,area,stats
begin

  dtr = acos(0.0)/90.
  if (opt) then
    if (isatt(opt,"lat_bnds")) then
      bnds = onedtond(opt@lat_bnds,(/dimsizes(opt@lat_bnds)/2,2/))
      area = conform(mdat, tofloat(sin(bnds(:,1)*dtr)-sin(bnds(:,0)*dtr)), 0)
    end if
  end if
  if (.not.isvar("area")) then
     area = conform(mdat, tofloat(cos(mdat&$mdat!0$*dtr)), 0)
  end if

  stats = new(4,float)
  stats(0) = wgt_areaave2 (mdat, area, 0)
  stats(1) = wgt_areaave2 (odat, area, 0)
  stats(2) = wgt_arearmse2(mdat, odat, area, 0)
  stats(3) = pattern_cor(mdat, odat, area, 0)

  return stats
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; get unlimited dimension name

undef ("get_unlim_dimension_name")
function get_unlim_dimension_name (f:file)
local dimNames, timeName, i
begin

  dimNames = getvardims(f)

  do i = 0, dimsizes(dimNames)-1
    if (isunlimited(f,dimNames(i))) then 
      timeName = dimNames(i)
    end if
  end do

  if (.not.isvar("timeName")) then 
     if (any(dimNames .eq. "time")) then 
        timeName = "time" ; fix for hadisst (and maybe other data sets)
     else 
        print ("ERROR: no unlimited dimension")
        status_exit(1)
     end if
  end if

  return timeName
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  determine the starting/ending years for a data file
;  given the target start/end years
;
;  y1, y2 = target (model analysis) years
;  ofile  = observational file name 
;  NOTE: internal parameter "freq"

function determine_observation_years( y1:integer, y2:integer, ofile:string )
local fd,times,dates,timeName,n,nn,yrbeg,yrend,yr1,yr2,dyr1,dyr2,freq,mfreq,mbeg,mend
begin

  freq = -2                        ; daily = -2;    monthly = -1
  mfreq = toint(10^(abs(freq)*2))  ; daily = 10000; monthly = 100
  mbeg = (/  1,  101 /)
  mend = (/ 12, 1231 /)

  ; dates for data file
  fd = addfile(ofile,"r")
  timeName = get_unlim_dimension_name(fd)
  times = fd->$timeName$ ; read time values
  if (isatt(times,"calendar")) then 
    times@calendar = str_lower(times@calendar)
  end if
  dates = cd_calendar(times,freq) ; yr-mo-dy
  delete(fd)
 ;print ("dates: "+dates(0)+" to "+dates(dimsizes(dates)-1))

  yrbeg = -1 
  yrend = -1 
  ; find year of the first day of january (yrbeg)
  do n = 0,  dimsizes(dates)-1
     if (dates(n) % mfreq .eq. mbeg(abs(freq)-1)) then 
        yrbeg = dates(n)/mfreq
        break
     end if
  end do
  ; find year of last december (yrend)
  do nn = 1, dimsizes(dates)
     n = dimsizes(dates) - nn
     if (dates(n) % mfreq .eq. mend(abs(freq)-1)) then
        yrend = dates(n)/mfreq
        break
     end if
  end do
  delete(dates)

  ; determine new years
  yr1 = toint(y1)
  yr2 = toint(y2)
  if (yr1 .lt. yrbeg) then
    yr2 = yr2 + (yrbeg-yr1)
    yr1 = yrbeg
  end if
  if (yr2 .gt. yrend) then
    dyr1 = yr1-yrbeg
    dyr2 = yr2-yrend
    if (dyr2 .lt. dyr1) then
      dyr1 = dyr2
    end if
    if (dyr1 .gt. 0) then
      yr1 = yr1 - dyr1
    end if
    yr2 = yrend
  end if

  if (yr1 .ne. y1) then
    print("NOTE: Adjusted observed starting date from "+sprinti("%4.4i",y1)+" to "+yr1)
  end if
  if (yr2 .ne. y2) then
    print("NOTE: Adjusted observed ending date from "+sprinti("%4.4i",y2)+" to "+yr2)
  end if

  return (/yr1,yr2/)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  initalize plotting resources

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

  ; -- background fill for oceans --
  res@mpFillOn = False ; set to True for land-only
  res@mpOceanFillColor       = "lightblue"
  res@mpInlandWaterFillColor = "lightblue"
  res@mpLandFillColor        = 0 

  ; -- grid lines over oceans --
  res@mpGridAndLimbOn = False ; set to True for land-only
  res@mpGridMaskMode  = "MaskNotOcean"
  res@mpGridLineDashPattern  = 2  ; dotted

  ; -- geographic outlines --
  if (numLon*(maxLon - minLon)/(numLon-1) + 0.0001 .gt. 360.) then
    res@gsnAddCyclic = True
    res@mpOutlineBoundarySets = "Geophysical"
  else
    res@gsnAddCyclic = False
    res@mpOutlineBoundarySets = "National"
    res@mpOutlineSpecifiers = "AllUSStates"
  end if

  ; -- polar stereographic (when global hemisperic data) --
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

  ; -- when NOT ploar, global lat/lon plot start at 20W --
  if (res@gsnAddCyclic .and. .not.do_polar) then
    minLon = -20
    maxLon = 340
  end if

  ; -- contour fill mode for high resolution grids (saves time/space) --
  if (numLat*numLon .gt. 800000) then
    print ("Using raster fill mode (npts="+numLat*numLon+")")
    res@cnFillMode = "RasterFill"
    res@cnRasterMinCellSizeF = .0005
    res@cnRasterSmoothingOn = True
  end if

  ; -- set up plot shape/aspect --
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
    res@mpMinLonF = minLon
    res@mpMaxLonF = maxLon
    res@mpCenterLonF = (minLon + maxLon) * 0.5
  end if

  ; -- line labels (hard-coded above) --
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

  ; -- labelbar resources (only if no contour line labels) --
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

  ; -- plot title resources --
  res@gsnLeftStringFontHeightF = .020
  res@gsnRightStringFontHeightF = .020
  res@gsnLeftString     = ""
  res@gsnRightString    = ""
  res@tiMainString = " "
  res@tiMainOn = False

  return res
end




wcStart = True
wcTime = True
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;  S T A R T   O F  N C L   S C R I P T  ;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

begin
 ; timing
  wcStart@total = systemfunc("date +%s")

 ; required variables
  if (.not.isvar("var")) then
     print ("ERROR: variable not specifed")
     exit
  end if
 ; var = "uv250"

 ; optional variables
  if (.not.isvar("mdir")) then
     mdir = "model"
  end if
  if (.not.isvar("odir")) then
     odir = "obs"
  end if
  if (.not.isvar("bandpass")) then
     bandpass = False
  end if
  if (.not.isvar("stdev")) then
     stdev = False
  end if
  if (.not.isvar("yrbeg")) then
     yrbeg = 1979
  end if
  if (.not.isvar("yrend")) then
     yrend = 2008 
  end if
  if (.not.isvar("verbose")) then
     verbose = False
  end if
  if (.not.isvar("seasons")) then
     ssn_list = (/ 0, 1, 2, 3, 4 /)
  else
     ssn_list = toint(str_split(seasons,","))
  end if

;;;;;;;;;;;;;;;;;;;;;;;
; variable names

  level = str_get_cols(var,2,-1)
  var1 = var_letter_lookup(str_get_cols(var,0,0)) + level
  var2 = var_letter_lookup(str_get_cols(var,1,1)) + level
  
  variable_label = setup_variable_label(var)

 ;cin = fspan(20,160,8)
 ;ci = fspan(10,50,5)
 ;cind = array_append_record(-ci(::-1),ci,0)

  model_date_label = sprinti("%4.4i",yrbeg) + " to " + sprinti("%4.4i",yrend)

  if (bandpass) then
     filter_label = "2 to 7 Day Filter"
     filter_wghts = filwgts_lanczos (21, 2, 1./7., 1./2., 1.)
  else
     filter_label = "Unfiltered"
  end if

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read observation (entire file)

  wcStart@readObs = systemfunc("date +%s")
  ofile = systemfunc("/bin/ls "+odir+"/"+var1+".*.nc")
  if (dimsizes(ofile) .ne. 1) then
     print ("ERROR: did not find single obs file, instead found "+dimsizes(ofile)+" files")
     exit
  end if
  print ("ofile: "+ofile)
  obsyears = determine_observation_years(yrbeg,yrend,ofile)
  obs_date_label = obsyears(0) + " to " + obsyears(1)
  print("obs years: "+obs_date_label)
  fobs = addfile(ofile,"r")
  ; first check if latitude (dim 1) should be flipped
  ovdims = getfilevardims(fobs,var1)
 ;print (ovdims)
; otime = fobs->$ovdims(0)$
 ;olat = fobs->$ovdims(1)$
 ;olon = fobs->$ovdims(2)$
  otime   = read_observed_data(fobs,ovdims(0),obsyears(0),obsyears(1))
  oOpt    = time_indexing_init(otime)
  obsdat1 = read_observed_data(fobs,var1,obsyears(0),obsyears(1))
  olat = obsdat1&$obsdat1!1$
  olon = obsdat1&$obsdat1!2$
  delete(fobs)
; printVarSummary(otime)
; printVarSummary(obsdat1)

  if (var1 .ne. var2) then
     ofile = systemfunc("/bin/ls "+odir+"/"+var2+".*.nc")
     if (dimsizes(ofile) .ne. 1) then
        print ("ERROR: did not find single obs file(2), instead found "+dimsizes(ofile)+" files")
        exit
     end if
     print ("ofile: "+ofile)
     fobs = addfile(ofile,"r")
     obsdat2 = read_observed_data(fobs,var2,obsyears(0),obsyears(1))
     delete(fobs)
  end if
  wcTime@readObs = stringtoint(systemfunc("date +%s")) - stringtoint(wcStart@readObs)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; setup model file for reading

  wcStart@readModel = systemfunc("date +%s")
  mfile = systemfunc("/bin/ls "+mdir+"/"+var1+".*.nc")
  if (dimsizes(mfile) .ne. 1) then
     print ("ERROR: did not find single model file, instead found "+dimsizes(mfile)+" files")
     exit
  end if
  print ("mfile: "+mfile)
  fmod1 = addfile(mfile,"r")
  dnames = getfilevardims(fmod1,var1)
  mtime = fmod1->$dnames(0)$
  mlat = fmod1->$dnames(1)$
  mlon = fmod1->$dnames(2)$
  if (isatt(mtime,"calendar")) then
    mtime@calendar = str_lower(mtime@calendar) ; cd_calendar only recognizes lower case
  end if
  index = determine_start_end(mtime,yrbeg,yrend)
  print ("model time index: "+index(0)+" to "+index(1)+" (npts="+(index(1)-index(0)+1)+")")
  mOpt = time_indexing_init(mtime(index(0):index(1)))

  ; create new time dimension
  ntime = mOpt@tn
  copy_VarAtts(mtime,ntime)

  ; open second file if necessary
  if (var1 .ne. var2) then
     mfile = systemfunc("/bin/ls "+mdir+"/"+var2+".*.nc")
     if (dimsizes(mfile) .ne. 1) then
        print ("ERROR: did not find single model file(2), instead found "+dimsizes(mfile)+" files")
        exit
     end if
     print ("mfile: "+mfile)
     fmod2 = addfile(mfile,"r")
  end if

  ; check for latitude bounds
  ; only needed for computing area weighted statistics
  stRes = False
  if (isatt(mlat,"bounds")) then
    stRes = True
    stRes@lat_bnds = ndtooned(fmod1->$mlat@bounds$)
  end if
  wcTime@readModel = stringtoint(systemfunc("date +%s")) - stringtoint(wcStart@readModel)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; process observed data on observed grid (interpolate to model grid later)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if (bandpass) then
     if (verbose) then
        print ("Bandpass filter observed data")
     end if
     dimobs = toint(dimsizes(obsdat1))
     wcStart@bpObs = systemfunc("date +%s")
     EDDIES::filter(dimobs(2),dimobs(1),dimobs(0),obsdat1,obsdat1@_FillValue,filter_wghts,21)
    ;obsdat1 = band_pass_latlon_time(obsdat1,1,(/2,7,21/),False)
     if (var1 .ne. var2) then
        EDDIES::filter(dimobs(2),dimobs(1),dimobs(0),obsdat2,obsdat2@_FillValue,filter_wghts,21)
        ;obsdat2 = band_pass_latlon_time(obsdat2,1,(/2,7,21/),False)
     end if
     wcTime@bpObs = stringtoint(systemfunc("date +%s")) - stringtoint(wcStart@bpObs)
  end if
 ;print (var+": obsdat1(mn,mx) = "+min(obsdat1)+", "+max(obsdat1))
  if (verbose) then
     print ("Compute monthly mean variance of observed data")
  end if
  wcStart@varObs = systemfunc("date +%s")
  avg1 = average_daily_to_monthly(obsdat1,otime,oOpt)
 ;print (var+": avg1(mn,mx) = "+min(avg1)+", "+max(avg1))
  if (var1 .ne. var2) then
     avg2 = average_daily_to_monthly(obsdat2,otime,oOpt)
     vobsgrid = average_daily_to_monthly(obsdat1*obsdat2,otime,oOpt)-avg1*avg2
  else
     vobsgrid = average_daily_to_monthly(obsdat1^2,otime,oOpt)-avg1^2
  end if
  wcTime@varObs = stringtoint(systemfunc("date +%s")) - stringtoint(wcStart@varObs)
  vobsgrid!0 = dnames(0)
  vobsgrid&$vobsgrid!0$ = ntime
  vobsgrid!1 = ovdims(1)
  vobsgrid&$vobsgrid!1$ = olat
  vobsgrid!2 = ovdims(2)
  vobsgrid&$vobsgrid!2$ = olon
  delete([/obsdat1,avg1/])
  if (var1 .ne. var2) then
     delete([/obsdat2,avg2/])
  end if

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; setup output variance grids and pressure masking file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  vdsize = filevardimsizes(fmod1,var1)
  vdsize(0) = mOpt@nmonth
  vmod = new(vdsize,getfilevartypes(fmod1,var1))
  vdsize(0) = oOpt@nmonth
  vobs = new(vdsize,getfilevartypes(fmod1,var1))

  ; masking file
  quote = str_get_dq()
  USE_PS_EST = False
  pmvar = "pmaskv2"
  pmask = systemfunc("/usr/bin/perl -e '@f=glob "+quote+mdir+"/pmaskv2.*.nc"+quote+";print "+quote+"@f\n"+quote+"'")
  if (any(.not.isfilepresent(pmask))) then
   ; surface pressure
    pmvar = "ps"
    pmask = systemfunc("/usr/bin/perl -e '@f=glob "+quote+mdir+"/ps.*.nc"+quote+";print "+quote+"@f\n"+quote+"'")
    if (any(.not.isfilepresent(pmask))) then
      ; estimated surface pressure
      pmask = systemfunc("/usr/bin/perl -e '@f=glob "+quote+mdir+"/ps_est.nc"+quote+";print "+quote+"@f\n"+quote+"'")
      if (all(isfilepresent(pmask)) .and. dimsizes(pmask) .eq. 1) then
        USE_PS_EST = True
      end if
    end if
  end if

  if (all(isfilepresent(pmask)) .and. dimsizes(pmask) .eq. 1) then
     print ("Open pmask file: "+pmask)
     fmsk = addfile(pmask,"r")
  else
     print ("ERROR: did not find single pressure masking file, instead found "+dimsizes(pmask)+" files")
     exit
  end if

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; process model data looping through latitude to save memory
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  jdim = 19
  js = 0

  wcTime@readModel = 0
  wcTime@readMask = 0
  wcTime@varModel = 0
  if (bandpass) then
    wcTime@bpModel = 0
  end if

; loop over model latitudes
  do while (js .le. dimsizes(mlat)-1)
     je = min((/js+jdim,dimsizes(mlat)-1/))
     print("j = "+je)
     if (verbose) then
        print ("Reading model data")
     end if
     wcStart@readModel = systemfunc("date +%s")
     dailyTS1 = fmod1->$var1$(:,js:je,:)
     if (var1 .ne. var2) then
        dailyTS2 = fmod2->$var2$(:,js:je,:)
     end if
     wcTime@readModel = wcTime@readModel + (stringtoint(systemfunc("date +%s")) - stringtoint(wcStart@readModel))
     ; band-pass filter?
     if (bandpass) then
        if (verbose) then
           print ("Bandpass filter model data")
        end if
        dim = toint(dimsizes(dailyTS1))
        wcStart@bpModel = systemfunc("date +%s")
        EDDIES::filter(dim(2),dim(1),dim(0),dailyTS1,dailyTS1@_FillValue,filter_wghts,21)
       ;dailyTS1 = band_pass_latlon_time(dailyTS1,1,(/2,7,21/),False)
        if (var1 .ne. var2) then
           EDDIES::filter(dim(2),dim(1),dim(0),dailyTS2,dailyTS2@_FillValue,filter_wghts,21)
          ;dailyTS2 = band_pass_latlon_time(dailyTS2,1,(/2,7,21/),False)
        end if
        wcTime@bpModel = wcTime@bpModel + (stringtoint(systemfunc("date +%s")) - stringtoint(wcStart@bpModel))
     end if
     ; read surface pressure for masking (band-pass also?)
     if (verbose) then
        print ("Read pressure mask and compute monthly mean")
     end if
     wcStart@readMask = systemfunc("date +%s")
     if (.not.USE_PS_EST) then
       pmsk = average_daily_to_monthly(fmsk->$pmvar$(:,js:je,:),mtime(index(0):index(1)),mOpt)
       if (isfilevaratt(fmsk,pmvar,"units")) then
         if (any(str_lower(fmsk->$pmvar$@units) .eq. (/"pa","pascals"/))) then
           pmsk = pmsk*0.01
         end if
       end if
     else
       dim = toint(dimsizes(dailyTS1))
       pmsk = conform_dims( (/mOpt@nmonth,dim(1),dim(2)/), fmsk->ps(js:je,:), (/1,2/) )
     end if
     wcTime@readMask = wcTime@readMask + (stringtoint(systemfunc("date +%s")) - stringtoint(wcStart@readMask))

     if (verbose) then
        print ("Compute monthly mean variance of model data")
     end if
     wcStart@varModel = systemfunc("date +%s")
     avg1 = average_daily_to_monthly(dailyTS1,mtime(index(0):index(1)),mOpt)
     if (var1 .ne. var2) then
        avg2 = average_daily_to_monthly(dailyTS2,mtime(index(0):index(1)),mOpt)
        vmod(:,js:je,:) = where( pmsk .ge. tofloat(level), \
                                 average_daily_to_monthly(dailyTS1*dailyTS2,mtime(index(0):index(1)),mOpt)-avg1*avg2, \
                                 dailyTS1@_FillValue )
        delete([/dailyTS1,avg1/])
        delete([/dailyTS2,avg2/])
     else
        vmod(:,js:je,:) = where( pmsk .ge. tofloat(level), \
                                 average_daily_to_monthly(dailyTS1^2,mtime(index(0):index(1)),mOpt)-avg1^2, \
                                 dailyTS1@_FillValue )
        delete([/dailyTS1,avg1/])
     end if
     wcTime@varModel = wcTime@varModel + (stringtoint(systemfunc("date +%s")) - stringtoint(wcStart@varModel))

     js = je+1
  end do

  ; copy dimensions
  vmod!0 = dnames(0)
  vmod&$vmod!0$ = ntime
  vmod!1 = dnames(1)
  vmod&$vmod!1$ = mlat
  vmod!2 = dnames(2)
  vmod&$vmod!2$ = mlon
 ;printVarSummary(vmod)

  ; close the input files
  delete(fmod1)
  if (var1 .ne. var2) then
     delete(fmod2)
  end if

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; interpolate observation to model grid
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; observation
  if (verbose) then
     print ("Interpolation of observed data")
  end if
  wcStart@interp = systemfunc("date +%s")
  vobsgrid@_FillValue = vmod@_FillValue
  vobs = linint2_Wrap(vobsgrid&$vobsgrid!2$, vobsgrid&$vobsgrid!1$, vobsgrid, True, vmod&$vmod!2$, vmod&$vmod!1$, 0)
  vobs = where( ismissing(vmod), vobs@_FillValue, vobs )
  wcTime@interp = stringtoint(systemfunc("date +%s")) - stringtoint(wcStart@interp)

  ;;;;;;;;;;;;;;;;;;;;;
  ; standard deviation
  ;;;;;;;;;;;;;;;;;;;;;

  if (stdev) then
     vmod = vmod^(0.5)
     vobs = vobs^(0.5)
  end if

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; compute climatology
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  season_names = (/"Annual", "DJF", "MAM", "JJA", "SON"/)
  wcStart@climo = systemfunc("date +%s")
  mclim = month_to_climo5_Wrap(vmod,vmod&$vmod!0$)
  oclim = month_to_climo5_Wrap(vobs,vobs&$vobs!0$)
  wcTime@climo = stringtoint(systemfunc("date +%s")) - stringtoint(wcStart@climo)
 ;printVarSummary(mclim)
 ;print ("max:"+max(mclim(1,:,:)))
  delete([/vmod,vobs/])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; scale data before plotting
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  scale = variable_scale_factor(var)
  if (scale .ne. 1.0) then
     print ("scale factor for "+var+": "+scale)
     mclim = mclim * scale
     oclim = oclim * scale
  end if

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; difference (model-obs) for plotting
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  dclim = mclim - oclim
  copy_VarMeta(mclim,dclim)

  if (verbose) then
     print("mclim(0,0) = "+mclim(0,0,0)+"   (FillValue = "+mclim@_FillValue+")")
     print("oclim(0,0) = "+oclim(0,0,0)+"   (FillValue = "+oclim@_FillValue+")")
     print("dclim(0,0) = "+dclim(0,0,0)+"   (FillValue = "+dclim@_FillValue+")")
  end if

  ; statistics
  print (var+": model(mn,mx) = "+min(mclim)+", "+max(mclim))
  print (var+":   obs(mn,mx) = "+min(oclim)+", "+max(oclim))
  print (var+":  diff(mn,mx) = "+min(dclim)+", "+max(dclim))
     
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;; Graphics Section ;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wcStart@graphics = systemfunc("date +%s")
  print("Plotting data, file = eddies."+var+".ps")
  wks = gsn_open_wks("ps","eddies."+var)
  rgb = setup_color_maps(wks,var)
  plot = new(3,graphic)

  ; general resources for plots
  res = initialize_resources(mclim(0,:,:))
  res@cnFillOn        = True  ; turn on color fill
  res@cnInfoLabelOn   = False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@gsnRightStringFontHeightF = .016
 ;res@gsnRightString = filter_label + "~C~" + model_date_label

  ; resources for panel plots
  resP = True
  resP@gsnMaximize = True
  resP@gsnFrame = False
  resP@gsnDraw  = True
  resP@txFontColor = 0
  resP@txFontHeightF = 0.030
  resP@txString = "Null" ; trick to allow room for titles

  ; resources for main titles
  resT = True
  resT@txFontHeightF = 0.022
  resS = True
  resS@txFontHeightF = 0.014

  ; initialize ascii output file for statistics
  ascii_record = "field          season        mean       bias       rmse       corr"
 ;system("echo "+ascii_record+">statistics."+var+".out")
  if (bandpass) then
    ascii_filter = " (bandpass)"
    bp = "bp"
  else
    ascii_filter = "           "
    bp = ""
  end if
  system("perl -e 'print "+quote+ascii_record+"\n"+quote+"'>>statistics."+var+bp+".out")  ; 'quote' defined above

  ; setup contour levels
  ; the first time thru season loop
  setupContours = True

  ;;;;;;;; loop thru seasons ;;;;;;;
  do n = 0, dimsizes(ssn_list)-1
     sea = ssn_list(n)

     ; compute statistics for this season
     ; write to stdout and ascii file
     stats = compute_statistics(mclim(sea,:,:),oclim(sea,:,:),stRes)
     bias = stats(0)-stats(1)
     print ("Season: "+sea+"   mean="+stats(0)+", bias="+bias+", rmse="+stats(2)+", corr="+stats(3))
     ascii_record_part1 = var+ascii_filter+"   "+str_lower(str_get_cols(season_names(sea),0,2))
     ascii_record_part2 = sprintf("%16.8g",stats(0))+sprintf("%16.8g",bias)+sprintf("%16.8g",stats(2))+sprintf("%12.8f",stats(3))
    ;system("echo "+ascii_record_part1+ascii_record_part2+">>statistics."+var+".out")
     system("perl -e 'print "+quote+ascii_record_part1+ascii_record_part2+"\n"+quote+"'>>statistics."+var+bp+".out")  ; 'quote' defined above

     res@cnFillPalette = rgb(0)
     
     ; try to set contour levels (first thru season loop)
     if (setupContours) then
        cin = setup_contour_levels(var,True,bandpass)
        if (dimsizes(cin) .gt. 1) then
           res@cnLevelSelectionMode = "ExplicitLevels"
           res@cnLevels = cin
        else
           print(var+": using AutomaticLevels standard contour plots")
        end if
        delete(cin)
     else
        res@cnLevelSelectionMode = "ExplicitLevels"
        if (isatt(res,"cnLevels")) then
           delete(res@cnLevels)
        end if
        res@cnLevels = contourLevels
     end if

     ; no contour lines for omega
     if (var1 .eq. "w" .and. var2 .eq. "w") then
        res@cnLinesOn = False
     else
        res@cnLinesOn = True
     end if

     ;;;;;;;;; top plot (model) ;;;;;;;;;
     res@gsnLeftString = "Model"
     if (isvar("mlab")) then
        res@gsnLeftString = res@gsnLeftString+": "+mlab
     end if
     res@gsnRightString = "a="+str_squeeze(sprintf("%7.3f",stats(0)))+" ("+model_date_label+")"
     plot(0) = gsn_csm_contour_map_ce(wks,mclim(sea,:,:),res)

     ; first time through retrieve contour levels
     if (setupContours) then
        getvalues plot@contour
           "cnLevels": contourLevels
        end getvalues
       ;print ("cnLevels: "+contourLevels)
     end if

     ;;;;;;;;; middle plot (obs) ;;;;;;;;;
     res@cnLevelSelectionMode = "ExplicitLevels"
     res@cnLevels = contourLevels
     res@gsnLeftString = "OBS"
     if (isvar("olab")) then
        res@gsnLeftString = res@gsnLeftString+": "+olab
     end if
     res@gsnRightString = "a="+str_squeeze(sprintf("%7.3f",stats(1)))+" ("+obs_date_label+")"
     plot(1) = gsn_csm_contour_map_ce(wks,oclim(sea,:,:),res)

     ;;;;;;;;; bottom plot (model-obs) ;;;;;;;;;
     res@cnFillPalette = rgb(1)
     if (isatt(res,"cnLevels")) then
        delete(res@cnLevels) ; reverts to "AutomaticLevels"
     end if
     
     ; try to set contour levels (first-time thru season loop)
     if (setupContours) then
        cin = setup_contour_levels(var,False,bandpass)
        if (dimsizes(cin) .gt. 1) then
           res@cnLevelSelectionMode = "ExplicitLevels"
           res@cnLevels = cin
        else
           symMinMaxPlt(dclim(sea,:,:),12,False,res)
           print(var+": using symMinMaxPlt to calculate symmetric contours for difference plots")
        end if
        delete(cin)
     else
        res@cnLevelSelectionMode = "ExplicitLevels"
        res@cnLevels = contourDiffLevels
     end if

     ; no contour lines for anything with omega
     if (var1 .eq. "w" .or. var2 .eq. "w") then
        res@cnLinesOn = False
     else
        res@cnLinesOn = True
     end if

     res@gsnLeftString = "Model - OBS"
     rightString = "~F33~D~F21~="+str_squeeze(sprintf("%7.3f",bias))
     res@gsnRightString = str_concat((/rightString, \
                        "   r="+str_squeeze(sprintf("%6.4f",stats(3)))+ \
                        "   E="+str_squeeze(sprintf("%7.3f",stats(2)))/))

     plot(2) = gsn_csm_contour_map_ce(wks,dclim(sea,:,:),res)

     ; first time through retrieve contour levels
     if (setupContours) then
        getvalues plot@contour
           "cnLevels": contourDiffLevels
        end getvalues
       ;print ("diff cnLevels: "+contourDiffLevels)
     end if

     gsn_panel(wks,plot,(/3,1/),resP)
     gsn_text_ndc(wks, season_names(sea) + " " + variable_label, 0.50, 0.97, resT)
     if (bandpass) then
       gsn_text_ndc(wks, filter_label, 0.50, 0.94, resS)
     end if
     frame(wks)
     setupContours = False
  end do

  ;;;;;;;;;;;;;;;;;;
  ; timing section
  ;;;;;;;;;;;;;;;;;;

  wcTime@graphics = stringtoint(systemfunc("date +%s")) - stringtoint(wcStart@graphics)
  wcTime@total = stringtoint(systemfunc("date +%s")) - stringtoint(wcStart@total)

  wcLabels = True
  wcLabels@total = "Total Runtime"
  wcLabels@graphics = "Graphics"
  wcLabels@readObs = "Obs: Read Data"
  wcLabels@bpObs = "Obs: Bandpass"
  wcLabels@varObs = "Obs: Averaging"
  wcLabels@readModel = "Model: Read Data"
  wcLabels@bpModel = "Model: Bandpass"
  wcLabels@varModel = "Model: Averaging"

  atts = getvaratts(wcTime)
  if (dimsizes(atts) .gt. 0) then
    if (isatt(wcTime,"total")) then
      if (wcTime@total .gt. 0) then
        total_runtime = wcTime@total
      end if
    end if
    NLINE = chartostring(inttochar(10)) ; new line character
    output = ""
    do n = 0, dimsizes(atts)-1
       if (isatt(wcLabels,atts(n))) then
          label = wcLabels@$atts(n)$
       else
          label = atts(n)
       end if
       output = output + str_insert(str_get_cols(label,0,15),"   ",16) + sprinti("%4.1i",wcTime@$atts(n)$) + " sec"
       if (isvar("total_runtime")) then
         pct = toint(100.*(tofloat(wcTime@$atts(n)$)/tofloat(total_runtime))+0.5)
         output = output + sprinti("%7.1i",pct) + "%" + NLINE
       else
         output = output + NLINE
       end if
    end do
    print(NLINE+output)
    ; also output timing in ascii file
    output = str_sub_str(output,NLINE,"\n")
    header = "\ntiming for "+var
    if (bandpass) then
      header = header + " (with bandpass filter)"
    end if
    system("perl -e 'print "+quote+header+output+quote+"'>>timing.out")  ; 'quote' defined above
  end if

end
