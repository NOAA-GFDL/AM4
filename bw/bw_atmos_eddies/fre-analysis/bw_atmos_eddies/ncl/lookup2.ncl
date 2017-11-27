
undef("dspan")
function dspan (x[*]:numeric)
local dimx,xx
begin
  dimx = dimsizes(x)
  xx = tofloat(array_append_record(-x(::-1),x,0))
  return xx
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("setup_contour_levels")
function setup_contour_levels (var:string,clevs:logical,bpass:logical)
local ci
begin
   cnLevs = True
  ;----------------------------------------------
  ; standard contour levels with bandpass filter
  ;----------------------------------------------
   if (clevs .and. bpass) then
      cnLevs@uu200 = fspan(10,90,9)
      cnLevs@uu250 = cnLevs@uu200
      cnLevs@uu850 = fspan(3,30,10)
      cnLevs@uv200 = dspan(fspan(3,39,6))
      cnLevs@uv250 = cnLevs@uv200
      cnLevs@uv850 = dspan(fspan(2,14,7))
      cnLevs@vv200 = fspan(20,180,9)
      cnLevs@vv250 = cnLevs@vv200
      cnLevs@vv850 = fspan(3,36,12)
      cnLevs@vq850 = dspan(fspan(1,5,5))
      cnLevs@vt200 = dspan((/1,2,3,4,6,8,10/))
      cnLevs@vt250 = cnLevs@vt200
      cnLevs@vt500 = dspan((/1,2,3,4,6,8,10/))
      cnLevs@vt850 = dspan(fspan(2,12,6))
      cnLevs@ww500 = fspan(5,50,10)
      cnLevs@ww850 = fspan(2,20,10)
      cnLevs@wq850 = dspan(fspan(20,120,6))
      cnLevs@wt500 = dspan(fspan(40,200,5))
      cnLevs@wt850 = dspan((/10,30,50,100,200/))
      cnLevs@qq850 = fspan(0.50,4,8)
      cnLevs@qt850 = fspan(-1,9,10)
      cnLevs@tt200 = fspan(1,8,8)
      cnLevs@tt250 = cnLevs@tt200
      cnLevs@tt500 = fspan(1,12,12)
      cnLevs@tt850 = fspan(1,12,12)
      cnLevs@hh200 = fspan(10,120,12)
      cnLevs@hh250 = cnLevs@hh200
      cnLevs@hh500 = fspan(5,60,12)
      cnLevs@hh850 = tofloat((/2,4,6,10,15,20,25,30/))
      cnLevs@hh10  = fspan(20,160,8)
      cnLevs@tt10  = fspan(1,14,14)
      cnLevs@uu10  = fspan(2,20,10)

  ;------------------------------------------------
  ; difference contour levels with bandpass filter
  ;------------------------------------------------
   else if (.not.clevs .and. bpass) then
      cnLevs@uu200 = dspan(fspan(5,30,6))
      cnLevs@uu250 = cnLevs@uu200
      cnLevs@uu850 = dspan(fspan(2,12,6))
      cnLevs@uv200 = dspan(fspan(3,21,4))
      cnLevs@uv250 = cnLevs@uv200
      cnLevs@uv850 = dspan(fspan(1,9,5))
      cnLevs@vv200 = dspan((/10,20,40,60/))
      cnLevs@vv250 = cnLevs@vv200
      cnLevs@vv850 = dspan((/2,4,6,10,15,20/))
      cnLevs@vq850 = dspan(fspan(0.5,2,4))
      cnLevs@vt200 = dspan(fspan(1,5,5))
      cnLevs@vt250 = cnLevs@vt200
      cnLevs@vt500 = dspan((/0.2,0.6,1.,1.5,2.,3./))
      cnLevs@vt850 = dspan(fspan(1,5,5))
      cnLevs@ww500 = dspan(fspan(4,20,5))
      cnLevs@ww850 = dspan(fspan(2,10,5))
      cnLevs@wq850 = dspan(fspan(10,60,6))
      cnLevs@wt500 = dspan(fspan(10,70,4))
      cnLevs@wt850 = dspan(fspan(10,70,4))
      cnLevs@qq850 = dspan((/0.25,0.50,1.0,1.5,2.0/))
      cnLevs@qt850 = dspan((/0.5,1,2,3,4,5/))
      cnLevs@tt200 = dspan(fspan(0.5,2.5,5))
      cnLevs@tt250 = cnLevs@tt200
      cnLevs@tt500 = dspan((/0.5,1,2,3,4,5/))
      cnLevs@tt850 = dspan(fspan(1,6,6))
      cnLevs@hh200 = dspan((/5,10,20,30,40/))
      cnLevs@hh250 = cnLevs@hh200
      cnLevs@hh500 = dspan((/2,4,8,12,16/))
      cnLevs@hh850 = dspan(fspan(1,9,5))
      cnLevs@hh10  = dspan(fspan(10,70,7))
      cnLevs@tt10  = dspan(fspan(1,8,8))
      cnLevs@uu10  = dspan(fspan(2,10,5))

  ;------------------------------------------
  ; standard contour levels with NO filter
  ;------------------------------------------
   else if (clevs .and. .not.bpass) then
      ; contour plots - unfiltered
      cnLevs@uu200 = fspan(30,300,10)
      cnLevs@uu250 = cnLevs@uu200
      cnLevs@uu850 = fspan(10,100,10)
      cnLevs@uv200 = dspan(fspan(10,90,5))
      cnLevs@uv250 = cnLevs@uv200
      cnLevs@uv850 = dspan(fspan(5,30,6))
      cnLevs@vv200 = fspan(30,360,12)
      cnLevs@vv250 = cnLevs@vv200
      cnLevs@vv850 = fspan(10,100,10)
      cnLevs@vq850 = dspan(fspan(2,10,5))
      cnLevs@vt850 = dspan(fspan(5,25,5))
      cnLevs@wq850 = dspan((/25,50,100,150,200/))
     ;cnLevs@wq850 = tofloat((/-200,-150,-100,-50,-25,0,25,50/))
      cnLevs@wt500 = dspan(fspan(50,300,6))
      cnLevs@wt850 = dspan(fspan(50,300,6))
      cnLevs@ww500 = fspan(5,50,10)
      cnLevs@ww850 = fspan(5,50,10)
      cnLevs@qq850 = tofloat((/1,2,3,5,7,9,12,15/))
      cnLevs@qt850 = fspan(-4,7,12)
      cnLevs@tt200 = fspan(2,20,10)
      cnLevs@tt250 = cnLevs@tt200
      cnLevs@tt500 = fspan(3,30,10)
      cnLevs@tt850 = fspan(5,50,10)
      cnLevs@hh200 = tofloat((/5,15,30,60,100,150,200,250,300,400/))
      cnLevs@hh250 = cnLevs@hh200
      cnLevs@hh500 = tofloat((/5,10,25,50,75,100,125,150,200,250,300/))
      cnLevs@hh850 = tofloat((/2,5,10,20,40,60,90,120,160/))
      cnLevs@hh10  = tofloat((/200,400,600,800,1000,1200,1600,2000,2400/))
      cnLevs@tt10  = fspan(10,100,10)
      cnLevs@uu10  = fspan(30,300,10)

  ;------------------------------------------
  ; difference contour levels with NO filter
  ;     (difference plots - unfiltered)
  ;------------------------------------------
   else if (.not.clevs .and. .not.bpass) then
      cnLevs@uu200 = dspan(fspan(10,90,5))
      cnLevs@uu250 = cnLevs@uu200
      cnLevs@uu850 = dspan(fspan(4,20,5))
      cnLevs@uv200 = dspan(fspan(10,50,5))
      cnLevs@uv250 = cnLevs@uv200
      cnLevs@uv850 = dspan(fspan(2,18,5))
      cnLevs@vv200 = dspan((/10,20,40,60,80,100/))
      cnLevs@vv250 = cnLevs@vv200
      cnLevs@vv850 = dspan(fspan(5,25,5))
      cnLevs@vq850 = dspan(fspan(1,5,5))
      cnLevs@vt850 = dspan(fspan(3,15,5))
      cnLevs@wq850 = dspan((/10,20,50,100/))
     ;cnLevs@wq850 = dspan((/5,10,25,50/))
      cnLevs@wt500 = dspan((/10,20,50,100/))
      cnLevs@wt850 = dspan((/10,20,50,100/))
      cnLevs@ww500 = dspan((/2,5,10,20/))
      cnLevs@ww850 = dspan((/2,5,10,20/))
      cnLevs@qq850 = dspan((/1,2,3,5,7,10/))
      cnLevs@qt850 = dspan((/0.25,0.5,1,2,3/))
      cnLevs@tt200 = dspan(fspan(1,6,6))
      cnLevs@tt250 = cnLevs@tt200
      cnLevs@tt500 = dspan(fspan(1,6,6))
      cnLevs@tt850 = dspan((/2,5,10,15,20,30/))
      cnLevs@hh200 = dspan(fspan(25,125,5))
      cnLevs@hh250 = cnLevs@hh200
      cnLevs@hh500 = dspan(fspan(20,100,5))
      cnLevs@hh850 = dspan(fspan(10,50,5))
      cnLevs@hh10  = dspan((/100,200,400,600,800,1000/))
      cnLevs@tt10  = dspan(fspan(5,25,5))
      cnLevs@uu10  = dspan(fspan(30,150,5))
   end if
   end if
   end if
   end if
   
   if (isatt(cnLevs,var)) then
      ci = cnLevs@$var$
   else
      ci = (/0.0/)
   end if
   return ci
end


undef("setup_variable_label")
function setup_variable_label (var:string)
local name,prime,lev,ch,i,lab,units,label
begin
  prime = "'"
  name = str_get_cols(var,0,1)
  lev = str_get_cols(var,2,-1)
  ch = new(2,string)
  do i = 0, 1
     ch(i) = str_get_cols(var,i,i)
     if (ch(i) .eq. "t") then
        ch(i) = "T"
     end if
  end do
  lab = lev + "hPa " + ch(0) + prime + ch(1) + prime

  units = ""
  if (name .eq. "uu" .or. name .eq. "vv" .or. name .eq. "uv") then
     units = "(m~S~2~N~ s~S~-2~N~)"
  end if
  if (name .eq. "zz" .or. name .eq. "hh") then
     units = "(dm~S~2~N~)"
  end if
  if (name .eq. "tt") then
     units = "(K~S~2~N~)"
  end if
  if (name .eq. "ww") then
     units = "(10~S~3~N~ mb~S~2~N~ dy~S~-2~N~)"
  end if
  if (name .eq. "qq") then
     units = "(g~S~2~N~ kg~S~-2~N~)"
  end if
  if (name .eq. "qt") then
     units = "(g kg~S~-1~N~ K)"
  end if

  label = lab + " " + units
  return label
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("variable_scale_factor")
function variable_scale_factor (var:string)
local ch,i
begin
  scale = 1.0
  do i = 0, 1
     ch = str_get_cols(var,i,i)
     if (ch .eq. "q") then
        scale = scale * 1000.   ; kg/kg to g/kg
     end if
     if (ch .eq. "w") then
        scale = scale * 864.    ; Pa/s to hPa/dy
     end if
     if (ch .eq. "z" .or. ch .eq. "h") then
        scale = scale * 0.1     ; m to dm (decimeters)
     end if
  end do
  ; additional scaling for ww
  if (str_get_cols(var,0,1) .eq. "ww") then
     scale = scale * 0.001
  end if
  return scale
end

