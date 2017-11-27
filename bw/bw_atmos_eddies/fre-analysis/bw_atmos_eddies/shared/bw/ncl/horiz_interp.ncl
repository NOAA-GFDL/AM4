
external HORIZ_INTERP "./horiz_interp.so"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("horiz_interp_lon_axis")
function horiz_interp_lon_axis (nx[1]:numeric, x:numeric)
local dimx,xb
begin
  ; check the size of the input x dimension
  dimx = dimsizes(x)
  if (dimsizes(dimx) .eq. 1 .and. dimx .eq. nx) then
    ; average from mid-point value to grid box boundary
    ; then apply cyclic continuity for end points
    xb = new((/nx,2/),double)
    if (typeof(x) .eq. "float") then
      xb(1:nx-1,0) = flt2dble((x(0:nx-2)+x(1:nx-1))*0.5)
      xb(0,0) = flt2dble(1.5*x(0)-0.5*x(1))
    else
      xb(1:nx-1,0) = (x(0:nx-2)+x(1:nx-1))*0.5
      xb(0,0) = 1.5*x(0)-0.5*x(1)
    end if
    xb(0:nx-2,1) = xb(1:nx-1,0)
    xb(nx-1,1) = xb(0,0) + 360.d
  else if (dimsizes(dimx) .eq. 1 .and. dimx .eq. nx+1) then
    ; reshape boundary point array
    xb = new((/nx,2/),double)
    if (typeof(x) .eq. "float") then
      xb(0:nx-1,0) = flt2dble(x(0:nx-1))
      xb(0:nx-1,1) = flt2dble(x(1:nx))
    else
      xb(0:nx-1,0) = x(0:nx-1)
      xb(0:nx-1,1) = x(1:nx)
    end if
  else if (dimsizes(dimx) .eq. 2 .and. dimx(0) .eq. nx .and. dimx(1) .eq. 2) then
    if (typeof(x) .eq. "float") then
      xb = flt2dble(x)
    else
      xb = x
    end if
  else
    print("ERROR in horiz_interp: invalid array size for input longitude")
    status_exit(1)
  end if
  end if
  end if

  ; define dimensions if input axis was 1d
  if (dimsizes(dimx) .eq. 1) then
    xb!0 = x!0
    xb!1 = "bnds"
   ;xb&$xb!0$ = x&$x!0$
    if (isatt(x,"units")) then
      xb@units = x@units
    end if
  end if

  return xb
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("horiz_interp_lat_axis")
function horiz_interp_lat_axis (ny[1]:numeric, y:numeric)
local dimy,yb
begin
  ; check the size of the input y dimension
  dimy = dimsizes(y)
  if (dimsizes(dimy) .eq. 1 .and. dimy .eq. ny) then
    ; average from mid-point value to grid box boundary
    ; then set end points to pole values
    yb = new((/ny,2/),double)
    if (typeof(y) .eq. "float") then
      yb(1:ny-1,0) = flt2dble((y(0:ny-2)+y(1:ny-1))*0.5)
    else
      yb(1:ny-1,0) = (y(0:ny-2)+y(1:ny-1))*0.5
    end if
    yb(0:ny-2,1) = yb(1:ny-1,0)
    yb(0,0) = -90.
    yb(ny-1,1) = 90.
  else if (dimsizes(dimy) .eq. 1 .and. dimy .eq. ny+1) then
    ; reshape boundary point array
    yb = new((/ny,2/),double)
    if (typeof(y) .eq. "float") then
      yb(0:ny-1,0) = flt2dble(y(0:ny-1))
      yb(0:ny-1,1) = flt2dble(y(1:ny))
    else
      yb(0:ny-1,0) = y(0:ny-1)
      yb(0:ny-1,1) = y(1:ny)
    end if
  else if (dimsizes(dimy) .eq. 2 .and. dimy(0) .eq. ny .and. dimy(1) .eq. 2) then
    if (typeof(y) .eq. "float") then
      yb = flt2dble(y)
    else
      yb = y
    end if
  else
    print("ERROR in horiz_interp: invalid array size for input longitude")
    status_exit(1)
  end if
  end if
  end if

  ; define dimensions if input axis was 1d
  if (dimsizes(dimy) .eq. 1) then
    yb!0 = y!0
    yb!1 = "bnds"
   ;yb&$yb!0$ = y&$y!0$
    if (isatt(y,"units")) then
      yb@units = y@units
    end if
  end if

  return yb
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("horiz_interp_axis_size")
function horiz_interp_axis_size(dat:numeric,axis[1]:string)
local dims, size,eps,test,del
begin
  ; check the size of the dimensions
  dims = dimsizes(dat)
  if (dimsizes(dims) .eq. 2) then
    size = dims(0)
  else if (dimsizes(dims) .eq. 1) then
    eps =1.e-4
    if (str_lower(axis) .eq. "x") then
      ; cyclic continuity test
      test = dat(dims-1) - dat(0) - 360.
      del = dat(dims-1) - dat(dims-2)
      if (abs(test) .le. eps) then
        size = dims-1
      else if (abs(test+del) .le. eps) then
        size = dims
      else
        print("ERROR in horiz_interp: could not determine length of 1-dimensional x-axis")
        status_exit(1)
      end if
      end if
    else if (str_lower(axis) .eq. "y") then
      ; check first and last values for poles
      if (abs(dat(0)+90.) .le. eps .and. abs(dat(dims-1)-90.) .le. eps) then
        size = dims-1
      else
        size = dims
      end if
    else
      print("ERROR in horiz_interp: invalid axis identifier: must be 'x' or 'y'")
      status_exit(1)
    end if
    end if
  else
    print("ERROR in horiz_interp: invalid axis shape (must 1d or 2d)")
    status_exit(1)
  end if
  end if
  return size
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("horiz_interp_Wrap")
function horiz_interp_Wrap(xi:numeric, yi:numeric, fi:numeric, xo:numeric, yo:numeric, opt:logical)
local dimi,rank,nxi,nyi,nk,nl,nxo,nyo,din,dout,dimo,xbi,ybi,xbo,ybo,crit,miss,verb,fo
begin

  ; check the rank of the input data
  dimi = dimsizes(fi)
  rank = dimsizes(dimi)
  nxi = dimi(rank-1)
  nyi = dimi(rank-2)
  if (rank .eq. 2) then
     nl = 1
     nk = 1
     din = new((/nl,nk,nyi,nxi/),typeof(fi))
     din(0,0,:,:) = fi
  else if (rank .eq. 3) then
     nl = 1
     nk = dimi(0)
     din = new((/nl,nk,nyi,nxi/),typeof(fi))
     din(0,:,:,:) = fi
  else if (rank .eq. 4) then
     nl = dimi(0)
     nk = dimi(1)
     din = new((/nl,nk,nyi,nxi/),typeof(fi))
     din = fi
  else
    print("ERROR in horiz_interp: invalid rank ("+rank+") for input data")
    status_exit(1)
  end if
  end if
  end if
  din@_FillValue = fi@_FillValue

  ; check the size of the input dimensions
  ; returned axes will be double precision
  xbi = horiz_interp_lon_axis (nxi, xi)
  ybi = horiz_interp_lat_axis (nyi, yi)

  ; check the size of the output dimensions
  nxo = horiz_interp_axis_size(xo,"x")
  nyo = horiz_interp_axis_size(yo,"y")

  ; reshaped output axes
  ; returned axes will be double precision
  xbo = horiz_interp_lon_axis (nxo, xo)
  ybo = horiz_interp_lat_axis (nyo, yo)

  ; setup output data
  dout = new((/nl,nk,nyo,nxo/),typeof(din))
  
  ; optional arguments
  ; defaults
  crit = new(1,typeof(din))
  crit = 0.50
  verb = 0

  if (opt) then
    if (isatt(opt,"crit")) then
      crit = opt@crit
    end if
    if (isatt(opt,"verbose")) then
      verb = opt@verbose
    end if
  end if
  ; more values
  miss = din@_FillValue
  error = 0

  ; debug
  if (verb .gt. 0) then
    print("Arg types: xbi="+typeof(xbi)+", ybi="+typeof(ybi)+", din="+typeof(din))
    print("Arg types: xbo="+typeof(xbo)+", ybo="+typeof(ybo)+", dout="+typeof(dout))
    print("Arg types: crit="+typeof(crit)+",miss="+typeof(miss))
  end if

  ; call the appropriate interpolation scheme
  ; axes are always double precision
  if (typeof(din) .eq. "float") then
    HORIZ_INTERP::horiz_interp_float(nxi,nyi,nk,nl,xbi,ybi,din,nxo,nyo,xbo,ybo,dout,crit,miss,verb,error)
  else if (typeof(din) .eq. "double") then
    HORIZ_INTERP::horiz_interp_double(nxi,nyi,nk,nl,xbi,ybi,din,nxo,nyo,xbo,ybo,dout,crit,miss,verb,error)
  else
    print("ERROR: invalid input data type")
    status_exit(1)
  end if
  end if

  if (error .gt. 0) then
    print("ERROR in horiz_interp: bad exit status")
    status_exit(1)
  end if

  ; lat/lon for output (same type as input)
  lon = new(nxo,typeof(xo))
  lat = new(nyo,typeof(yo))
  if (typeof(xo) .eq. "float") then
    lon = tofloat((xbo(:,0)+xbo(:,1))*0.5)
  else
    lon = (xbo(:,0)+xbo(:,1))*0.5 ; probably double
  end if
  if (typeof(yo) .eq. "float") then
    lat = tofloat((ybo(:,0)+ybo(:,1))*0.5)
  else
    lat = (ybo(:,0)+ybo(:,1))*0.5
  end if

  ; prepare the output data
  if (rank .eq. 2) then
    fo = new((/nyo,nxo/),typeof(fi))
    fo = dout(0,0,:,:)
    fo!0 = "lat"
    fo!1 = "lon"
    fo&lat = lat
    fo&lon = lon
    fo@_FillValue = miss
    return fo
  else if (rank .eq. 3) then
    fo = new((/nk,nyo,nxo/),typeof(fi))
    fo = dout(0,:,:,:)
    fo!0 = din!0
    fo&$fo!0$ = din&$din!0$
    fo!1 = "lat"
    fo!2 = "lon"
    fo&lat = lat
    fo&lon = lon
    fo@_FillValue = miss
    return fo
  else ; rank eq 4
    dout!0 = din!0
    dout&$dout!0$ = din&$din!0$
    dout!1 = din!1
    dout&$dout!1$ = din&$din!1$
    dout!2 = "lat"
    dout!3 = "lon"
    dout&lat = lat
    dout&lon = lon
    dout@_FillValue = miss
    return dout
  end if
  end if

end

