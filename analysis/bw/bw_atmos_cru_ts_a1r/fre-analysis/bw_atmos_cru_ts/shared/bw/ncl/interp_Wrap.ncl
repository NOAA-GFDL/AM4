;************************************************************
; interpolation with weighting
;   dataIn = data and grid to be interpolated
;   dataOut = target grid
;   opt = optional arguments
;************************************************************
undef("interp_type_is_bilinear")
function interp_type_is_bilinear ( din:numeric, dout:numeric, opt:logical )
local dimo,dimi,ranko,ranki, dloni,dlati,dlati,dlato, bilinear
begin

  dimi = dimsizes(din)
  dimo = dimsizes(dout)
  ranki = dimsizes(dimi)
  ranko = dimsizes(dimo)

  dloni = din&$din!(ranki-1)$(1) - din&$din!(ranki-1)$(0)
  dlati = din&$din!(ranki-2)$(1) - din&$din!(ranki-2)$(0)
  dlono = dout&$dout!(ranko-1)$(1) - dout&$dout!(ranko-1)$(0)
  dlato = dout&$dout!(ranko-2)$(1) - dout&$dout!(ranko-2)$(0)

  ; resolution: high -> low  .OR.  low -> high
  if ( dloni .lt. dlono .and. dlati .lt. dlato ) then
     bilinear = False   ; high to low res
  else
     bilinear = True    ; low to high res
  end if

  if (isatt(opt,"bilinearOn")) then
     bilinear = opt@bilinearOn
  end if

  return bilinear
end


;************************************************************
; interpolation with weighting
;   dataIn = data and grid to be interpolated
;   dataOut = target grid
;   opt = optional arguments
;************************************************************
undef("interp_Wrap")
function interp_Wrap ( dataIn:numeric, maskIn:numeric, dataOut:numeric, opt:logical )
local dimIn,dimOut,rankIn,rankOut,latIn,lonIn,latOut,lonOut,dlatIn,dlonIn,dlonOut,dlatOut,cyclic,bilinear,critpc,quiet,x,w,xo,wo,i,opti,dw
begin
  ; dimensions and rank
  dimIn = dimsizes(dataIn)
  dimOut = dimsizes(dataOut)
  rankIn = dimsizes(dimIn)
  rankOut = dimsizes(dimOut)
  dimMask = dimsizes(maskIn)
  rankMask = dimsizes(dimMask)
  if (rankOut .lt. 2 .or. rankIn .lt. 2 .or. rankOut .ne. rankIn) then
     print ("ERROR: in interp_wgt_Wrap: incorrect number of dimensions for data")
     status_exit(1)
  end if
  if (rankMask .ne. 1 .and. rankMask .ne. 2) then
     print ("ERROR: in interp_wgt_Wrap: incorrect number of dimensions for mask")
     status_exit(1)
  end if

  latIn  = dataIn&$dataIn!(rankIn-2)$
  lonIn  = dataIn&$dataIn!(rankIn-1)$
  latOut = dataOut&$dataOut!(rankOut-2)$
  lonOut = dataOut&$dataOut!(rankOut-1)$

  dlatIn = latIn(1) - latIn(0)
  dlonIn = lonIn(1) - lonIn(0)
  dlatOut = latOut(1) - latOut(0)
  dlonOut = lonOut(1) - lonOut(0)

  ; is the input grid cyclic?
  if (360. - lonIn(dimsizes(lonIn)-1) + lonIn(0) - dlonIn  .lt. dlonIn/10.) then
     cyclic = True
  else
     cyclic = False
  end if

  ; determine default interpolation type
  ;   bilinear:     low res -> high res
  ;   conservative: high res -> low res

  bilinear = interp_type_is_bilinear( dataIn, dataOut, opt )

  if (isatt(opt,"critpc")) then
     critpc = tofloat(opt@critpc)
  else
     critpc = 50.0
  end if

  if (isatt(opt,"quietOn")) then
     quiet = opt@quietOn
  else
     quiet = False
  end if

  if (.not.quiet) then
     print ("In (dlon,dlat): "+dlonIn+", "+dlatIn)
     print ("Out(dlon,dlat): "+dlonOut+", "+dlatOut)
  end if

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;   bi-linear interpolation with mask ;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if (bilinear) then
     ; bilinear interpolation
     ; weighting of data for missing values
     if (any(ismissing(dataIn)) .or. rankMask .eq. 2) then
        if (.not.quiet) then
           print ("using bi-linear interpolation with missing values, critpc = "+critpc)
        end if
        if (rankMask .eq. 2) then
           w = conform(dataIn,maskIn,(/rankIn-2,rankIn-1/)) 
           w = where(ismissing(dataIn),0.,w)
        else
           w = where(ismissing(dataIn),0.,1.)
        end if
        x = where(ismissing(dataIn),0.,dataIn*w)
        xo = linint2_Wrap(lonIn, latIn, x, cyclic, lonOut, latOut, 0)
        wo = linint2_Wrap(lonIn, latIn, w, cyclic, lonOut, latOut, 0)
        wo = where(wo .gt. 0.01*critpc, wo, wo@_FillValue)
        xo = where(ismissing(wo), xo@_FillValue, xo/wo)
     else
        if (.not.quiet) then
           print ("using bi-linear interpolation (no missing values)")
        end if
        xo = linint2_Wrap(lonIn, latIn, dataIn, cyclic, lonOut, latOut, 0)
     end if
     
    ; copy the dimensions
     if (rankIn .gt. 2) then
        do i = 0, rankIn-3
           xo!i = dataIn!i
           xo&$xo!i$ = dataIn&$dataIn!i$
        end do
     end if
     xo!(rankIn-2) = dataOut!(rankOut-2)
     xo&$xo!(rankIn-2)$ = dataOut&$dataOut!(rankOut-2)$
     xo!(rankIn-1) = dataOut!(rankOut-1)
     xo&$xo!(rankIn-1)$ = dataOut&$dataOut!(rankOut-1)$

  else
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;   area conserving interpolation  ;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

     if (.not.quiet) then
        print ("using area conserving interpolation (critpc="+critpc+")")
     end if
     opti = True
     opti@crit = critpc
     opti@verbose = 1
    ;xo = dataOut
     if (rankMask .eq. 1) then
        xo = horiz_interp_Wrap(lonIn, latIn, dataIn, lonOut, latOut, opti)
     else
        dw = (/rankIn-2,rankIn-1/) ; dimensions that mask conform to
        x = dataIn
        x = x * conform(x,maskIn,dw)
        xo = horiz_interp_Wrap(lonIn, latIn, x, lonOut, latOut, opti)
        wo = horiz_interp_Wrap(lonIn, latIn, maskIn, lonOut, latOut, False)
        wo = where(wo .gt. 0.01*critpc, wo, wo@_FillValue)
        xo = where(ismissing(conform(xo,wo,dw)), xo@_FillValue, xo/conform(xo,wo,dw))
     end if
        

  end if

  return xo
end

