
;************************************************************
; month_weights_in_days (xTime)
;    xTime = numeric time dimension
; Returns the number of days in the month.
; xTime is assumed to be consecutive months.
; First month is assumed to be January, last month December.
; The Julian/Gregorian calendar is assumed.
;************************************************************
undef("month_weights_in_days")
function month_weights_in_days (xTime:numeric)
local dimx,wts,nyrs,yearmo,n,yr1 
begin

  dimx = dimsizes(xTime)
  if ((dimx(0)%12).ne.0) then
     print ("ERROR in month_weights_in_days: numbr of time values must be a multiple of 12")
     exit
  end if
  wts = new(dimx,float)
  nyrs = dimx(0)/12

  yearmo = cd_calendar(xTime,-1)
  yr1 = toint(yearmo(0)/100)
  do n = 0, nyrs-1
     wts(12*n:12*n+11) = (/31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./)
     if (mod(yr1+n,4) .eq. 0) then
        wts(12*n+1) = wts(12*n+1)+1.0
     end if
  end do

  return wts
end

;************************************************************
undef("month_to_climo5_Wrap")
function month_to_climo5_Wrap (xMon:numeric, xTime:numeric)
local dimx,rank,nyrs,miss,ymiss,wts,dimOut,xSea,mAv,mWt,k,ks,ke,j,js,je
begin
  dimx = dimsizes(xMon)
  rank = dimsizes(dimx)
  if (rank.lt.2 .or. rank.gt.5) then
     print ("ERROR in month_to_climo5_Wrap: input array rank must be between 2 and 5")
     exit
  end if
  if ((dimx(0)%12).ne.0) then
     print ("ERROR in month_to_climo5_Wrap: dimension must be a multiple of 12")
     exit
  end if
  nyrs = toint(dimx(0))/12
  miss = 2*toint(dimx(0))/3  ; number of non-missing points for annual average (2/3 of all points)
  ymiss = 2*nyrs/3           ; number of non-missing points for monthly averages

; weights for each month (Julian/Gregorian assumed)
  wts = month_weights_in_days(xTime)

  dimOut = dimx
  dimOut(0) = 5 
  if (typeof(xMon) .eq. "short" .or. (typeof(xMon) .eq. "double" .and. rank .ge. 3)) then
     xSea = new(dimOut,float)
  else
     xSea = new(dimOut,typeof(xMon))
  end if

  xSea!0 = "season"
  xSea&season = ispan(0,4,1)
  xSea&season@long_name = "Season number (0=ann,1=djf,2=mam,3=jja,4=son)"
  xSea&season@units = "None"
 ;do k = 1, rank-1
 ;   xSea!k = xMon!k
 ;   xSea&$xSea!k$ = xMon&$xMon!k$
 ;end do

  ; annual average
  if (rank.eq.2) then
     if (typeof(xMon) .eq. "short") then
        xSea(0,:) = dim_avg_wgt_n_Wrap(short2flt(xMon),wts,miss,0)
     else
        xSea(0,:) = dim_avg_wgt_n_Wrap(xMon,wts,miss,0)
     end if
  else if (rank.eq.3) then
     if (typeof(xMon) .eq. "short") then
        xSea(0,:,:) = dim_avg_wgt_n_Wrap(short2flt(xMon),wts,miss,0)
     else if (typeof(xMon) .eq. "double" .and. typeof(xSea) .eq. "float") then
        xSea(0,:,:) = dim_avg_wgt_n_Wrap(doubletofloat_Wrap(xMon),wts,miss,0)
     else
        xSea(0,:,:) = dim_avg_wgt_n_Wrap(xMon,wts,miss,0)
     end if
     end if
  else if (rank.eq.4) then
     if (typeof(xMon) .eq. "short") then
        xSea(0,:,:,:) = dim_avg_wgt_n_Wrap(short2flt(xMon),wts,miss,0)
     else if (typeof(xMon) .eq. "double" .and. typeof(xSea) .eq. "float") then
        xSea(0,:,:,:) = dim_avg_wgt_n_Wrap(doubletofloat_Wrap(xMon),wts,miss,0)
     else
        xSea(0,:,:,:) = dim_avg_wgt_n_Wrap(xMon,wts,miss,0)
     end if
     end if
  else if (rank.eq.5) then
     if (typeof(xMon) .eq. "short") then
        xSea(0,:,:,:,:) = dim_avg_wgt_n_Wrap(short2flt(xMon),wts,miss,0)
     else if (typeof(xMon) .eq. "double" .and. typeof(xSea) .eq. "float") then
        xSea(0,:,:,:,:) = dim_avg_wgt_n_Wrap(doubletofloat_Wrap(xMon),wts,miss,0)
     else
        xSea(0,:,:,:,:) = dim_avg_wgt_n_Wrap(xMon,wts,miss,0)
     end if
     end if
  end if
  end if
  end if
  end if

  mWt = new(12,typeof(wts))
  dimOut(0) = 12
  mAv = new(dimOut,typeof(xSea))
  mAv!0 = "month"
 ;do k = 1, rank-1
 ;   mAv!k = xMon!k
 ;   mAv&$mAv!k$ = xMon&$xMon!k$
 ;end do

; month averages first (0=dec, 1=jan, ....)
  do ks = 0, 11
     ke = ks+(nyrs-1)*12
     k = mod(ks+1,12)
     if (rank.eq.2) then
        if (typeof(xMon) .eq. "short") then
           mAv(k,:) = dim_avg_wgt_n_Wrap(short2flt(xMon(ks:ke:12,:)),wts(ks:ke:12),ymiss,0)
        else
           mAv(k,:) = dim_avg_wgt_n_Wrap(xMon(ks:ke:12,:),wts(ks:ke:12),ymiss,0)
        end if
     else if (rank.eq.3) then
        if (typeof(xMon) .eq. "short") then
           mAv(k,:,:) = dim_avg_wgt_n_Wrap(short2flt(xMon(ks:ke:12,:,:)),wts(ks:ke:12),ymiss,0)
        else if (typeof(xMon) .eq. "double" .and. typeof(mAv) .eq. "float") then
           mAv(k,:,:) = dim_avg_wgt_n_Wrap(doubletofloat_Wrap(xMon(ks:ke:12,:,:)),wts(ks:ke:12),ymiss,0)
        else
           mAv(k,:,:) = dim_avg_wgt_n_Wrap(xMon(ks:ke:12,:,:),wts(ks:ke:12),ymiss,0)
        end if
        end if
     else if (rank.eq.4) then
        if (typeof(xMon) .eq. "short") then
           mAv(k,:,:,:) = dim_avg_wgt_n_Wrap(short2flt(xMon(ks:ke:12,:,:,:)),wts(ks:ke:12),ymiss,0)
        else if (typeof(xMon) .eq. "double" .and. typeof(mAv) .eq. "float") then
           mAv(k,:,:,:) = dim_avg_wgt_n_Wrap(doubletofloat_Wrap(xMon(ks:ke:12,:,:,:)),wts(ks:ke:12),ymiss,0)
        else
           mAv(k,:,:,:) = dim_avg_wgt_n_Wrap(xMon(ks:ke:12,:,:,:),wts(ks:ke:12),ymiss,0)
        end if
        end if
     else if (rank.eq.5) then
        if (typeof(xMon) .eq. "short") then
           mAv(k,:,:,:,:) = dim_avg_wgt_n_Wrap(short2flt(xMon(ks:ke:12,:,:,:,:)),wts(ks:ke:12),ymiss,0)
        else if (typeof(xMon) .eq. "double" .and. typeof(mAv) .eq. "float") then
           mAv(k,:,:,:,:) = dim_avg_wgt_n_Wrap(doubletofloat_Wrap(xMon(ks:ke:12,:,:,:,:)),wts(ks:ke:12),ymiss,0)
        else
           mAv(k,:,:,:,:) = dim_avg_wgt_n_Wrap(xMon(ks:ke:12,:,:,:,:),wts(ks:ke:12),ymiss,0)
        end if
        end if
     end if
     end if
     end if
     end if
     mWt(k) = sum(wts(ks:ke:12))
  end do

; season averages
  do j = 0, 3
     js = 3*j
     je = js+2
     if (rank.eq.2) then
        xSea(j+1,:) = dim_avg_wgt_n_Wrap(mAv(js:je,:),mWt(js:je),2,0)
     else if (rank.eq.3) then
        xSea(j+1,:,:) = dim_avg_wgt_n_Wrap(mAv(js:je,:,:),mWt(js:je),2,0)
     else if (rank.eq.4) then
        xSea(j+1,:,:,:) = dim_avg_wgt_n_Wrap(mAv(js:je,:,:,:),mWt(js:je),2,0)
     else if (rank.eq.5) then
        xSea(j+1,:,:,:,:) = dim_avg_wgt_n_Wrap(mAv(js:je,:,:,:,:),mWt(js:je),2,0)
     end if
     end if
     end if
     end if
  end do

  return xSea
end
;************************************************************
