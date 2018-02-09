
;************************************************************
; month_to_season_Wrap (xMon)
;    xMon = Time series of monthly data of rank 3 and any type.
;           Must start with January; end with December.
; output
;    xSea = Time series of seasonal data, dimensions (year,season, ...)
;           Rank will be input rank+1.
;           Seasons are: ANN, DJF, MAM, JJA, SON (first DJF is missing).
;************************************************************

undef("month_to_season_Wrap")
function month_to_season_Wrap (xMon:numeric)
local dimx,dimo,rank,nyrs,yrs,wts,yearmo,year1,xSea,yr,sea,ka,ks,ke,kws,kwe,miss,kdec
begin
  dimx = dimsizes(xMon)
  rank = dimsizes(dimx)
  if (rank.lt.3 .or. rank.gt.3) then
     print ("ERROR in month_to_season_Wrap: input array rank must be 3")
     exit
  end if
  if ((dimx(0)%12).ne.0) then
     print ("ERROR in month_to_season_Wrap: dimension must be a multiple of 12")
     exit
  end if
  nyrs = dimx(0)/12
  miss = 0

  ; define output array
  dimo = new(rank+1,integer)
  dimo(0) = nyrs
  dimo(1) = 5
  dimo(2:rank) = dimx(1:rank-1)
  xSea = new(dimo,typeof(xMon))

  yrs = new(nyrs,float)
  wts = new(13,float) ; dec-dec
  yearmo = cd_calendar(xMon&$xMon!0$,-1)
  year1 = toint(yearmo(0)/100)

  do yr = 0, nyrs-1
     ; save the year
     yrs(yr) = tofloat(year1+yr)
     ; month weights
     wts = (/31.,31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./)
     if (mod(year1+yr,4) .eq. 0) then
        wts(2) = 29.
     end if
     ;annual
     ka = yr*12 ;jan
     xSea(yr,0,:,:) = dim_avg_wgt_n_Wrap(xMon(ka:ka+11,:,:),wts(1:12), miss,0)
     ;seasons
     do sea = 1, 4
        ks = max((/ka+3*sea-4,0/))
        ke = ka+3*sea-2
        kwe = 3*sea-1
        kws = kwe - (ke-ks)
        xSea(yr,sea,:,:) = dim_avg_wgt_n_Wrap(xMon(ks:ke,:,:),wts(kws:kwe), miss,0)
     end do

     ; special case for first DJF
     ; use last December
     ; this will make the mean correct
     if (yr .eq. 0) then
        kdec = dimx(0)-1
        xSea(yr,1,:,:) = (xMon(0,:,:)*wts(1) + xMon(1,:,:)*wts(2) + xMon(kdec,:,:)*wts(0)) / sum(wts(0:2))
     end if
  
  end do

  ; dimensions
  xSea!0 = "year"
  xSea&year = yrs
  xSea!1 = "season"
  xSea&season = (/0,1,2,3,4/)
  xSea!2 = xMon!1
  xSea&$xSea!2$ = xMon&$xMon!1$
  xSea!3 = xMon!2
  xSea&$xSea!3$ = xMon&$xMon!2$

  return xSea
end

;************************************************************

undef("season_weights_in_days")
function season_weights_in_days (xSea[*][5][*][*]:numeric)
local dimx,year,nyrs,yr,wtx
begin
  dimx = dimsizes(xSea)
  wtx = new(dimx(0:1),float)
  year = xSea&$xSea!0$
  nyrs = dimsizes(year)
  do yr = 0, nyrs-1
    wtx(yr,:) = (/365.,90.,92.,92.,91./) ; non-leap year
    if (mod(year(yr),4) .eq. 0) then
        wtx(yr,0) = 366.
        wtx(yr,1) = 91.
    end if
  end do

  return wtx
end

;************************************************************

undef("compute_climatology_statistics")
function compute_climatology_statistics (xMon[*][*][*]:numeric,opt:logical)
local dimx,dimo,nyrs,wt,xSea,xStat,xmiss,xmask,trdx,sea,ngood,stats,nstats,istat
begin

  ; optional arguments
  if (opt) then
    if (isatt(opt,"stats")) then
      stats = opt@stats
    end if 
  end if
  if (.not.isvar("stats")) then
    stats = (/"Mean"/)
  end if
  nstats = dimsizes(stats)

  xSea = month_to_season_Wrap(xMon) ; xSea(year,season,lat,lon)
  wt = season_weights_in_days(xSea)
  dimx = dimsizes(xSea)  ; (yr,sea,lat,lon)
  nyrs = dimx(0)
  ngood = nyrs-1
  dimo = dimx
  dimo(0) = nstats ; number of statistics (mean, stdev, trend)
  xStat = new(dimo,float)
  istat = 0

  ; mean
  if (any(stats .eq. "Mean")) then
    do sea = 0, 4
      xStat(istat,sea,:,:) = dim_avg_wgt_n_Wrap(xSea(:,sea,:,:), wt(:,sea), ngood, 0)
    end do
    istat = istat + 1
  end if
  xSea(0,1,:,:) = xSea@_FillValue ; first DJF missing
  wt(0,1) = 0.0

  if (any(stats .eq. "StDev")) then
    xStat(istat,:,:,:) = dim_stddev_n_Wrap (xSea, 0)
    istat = istat + 1
  end if

  if (any(stats .eq. "Trend")) then
    trdx = dtrend_msg_n(fspan(0,nyrs-1,nyrs), xSea, True, True, 0)
    xStat(istat,:,:,:) = onedtond(trdx@slope, (/dimo(1),dimo(2),dimo(3)/)) * nyrs
  end if

  ; apply mask
  xStat@_FillValue = xSea@_FillValue
  xmiss = nyrs-1
  xmask = where(dim_sum_n(where(ismissing(xSea),0,1),0) .lt. xmiss, False, True)
  xStat = where(conform(xStat,xmask,(/1,2,3/)),xStat,xStat@_FillValue)

  ; dimensions (lat,lon are most important)
  xStat!0 = "stat"
  xStat&stat = ispan(0,nstats-1,1)
  xStat!1 = "season"
  xStat&season = ispan(0,4,1)
  xStat!2 = xMon!1
  xStat&$xStat!2$ = xMon&$xMon!1$
  xStat!3 = xMon!2
  xStat&$xStat!3$ = xMon&$xMon!2$

  return xStat
end

