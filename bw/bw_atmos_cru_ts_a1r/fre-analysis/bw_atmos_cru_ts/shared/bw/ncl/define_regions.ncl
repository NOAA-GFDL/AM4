 
undef("define_regions")
function define_regions(ffile:string)
local fRegions,REGIONS,i,n,fields,x,y,z,id,name,xmin,xmax,ymin,ymax
begin
; read regions
  fRegions = systemfunc("cat "+ffile)

  id   = new(dimsizes(fRegions),string)
  name = new(dimsizes(fRegions),string)
  xmin = new(dimsizes(fRegions),float)
  xmax = new(dimsizes(fRegions),float)
  ymin = new(dimsizes(fRegions),float)
  ymax = new(dimsizes(fRegions),float)

  n = 0
  do i = 0, dimsizes(fRegions)-1
    fields = str_split(fRegions(i),"|")
    if (dimsizes(fields).ge.3) then
      z = str_split(str_squeeze(fields(0)),":")
      x = str_split(str_squeeze(fields(1)),":")
      y = str_split(str_squeeze(fields(2)),":")
      id  (n) = z(0)
      name(n) = z(1)
      xmin(n) = tofloat(x(0))
      xmax(n) = tofloat(x(1))
      ymin(n) = tofloat(y(0))
      ymax(n) = tofloat(y(1))
      n = n+1
    end if
  end do

  REGIONS = True
  REGIONS@numReg = n
  REGIONS@id = id(0:n-1)
  REGIONS@title = name(0:n-1)
  REGIONS@lonbeg = xmin(0:n-1)
  REGIONS@lonend = xmax(0:n-1)
  REGIONS@latbeg = ymin(0:n-1)
  REGIONS@latend = ymax(0:n-1)
  return REGIONS
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("my_lonPivot")
function my_lonPivot (x:numeric, pivotLon:numeric)
local dimx,rank,lonName,temp,xlon,dlon,indP,mlon,indL,n,tlon,indt
begin
  dimx = dimsizes(x)
  rank = dimsizes(dimx)
  if (rank .lt. 2 .or. rank .gt. 5) then
    print ("ERROR: my_lonPivot: dims .ne. 2,3,4,5")
    exit
  end if
 ;print ("lonPivot= "+pivotLon)

  lonName = x!(rank-1)
  temp    = x 
  xlon    = x&$lonName$
  dlon    = xlon-pivotLon
  dlon    = where(dlon.lt.-360.,dlon+360.,dlon)
  dlon    = where(dlon.gt. 360.,dlon-360.,dlon)
  indP    = minind(where(dlon.ge.0,dlon,361.))

  mlon = dimx(rank-1)
  indL = mlon-1
  n    = indL-indP

 ;print("mlon= "+mlon)
 ;print("indP= "+indP)
 ;print("indL= "+indL)
 ;print("n= "+n)

  if (rank .eq. 2) then
     temp(:,0:n)  = (/ x(:,indP:indL)/)
     if (indP.gt.0) then
        temp(:,n+1:) = (/ x(:,0:indP-1) /)
     end if
  end if
  if (rank .eq. 3) then
     temp(:,:,0:n)  = (/ x(:,:,indP:indL)/)
     if (indP.gt.0) then
        temp(:,:,n+1:) = (/ x(:,:,0:indP-1) /)
     end if
  end if
  if (rank .eq. 4) then
     temp(:,:,:,0:n)  = (/ x(:,:,:,indP:indL)/)
     if (indP.gt.0) then
        temp(:,:,:,n+1:) = (/ x(:,:,:,0:indP-1) /)
     end if
  end if
  if (rank .eq. 5) then
     temp(:,:,:,:,0:n)  = (/ x(:,:,:,:,indP:indL)/)
     if (indP.gt.0) then
        temp(:,:,:,:,n+1:) = (/ x(:,:,:,:,0:indP-1) /)
     end if
  end if

  tlon       = new ( mlon, typeof(xlon) )
  tlon(0:n)  = (/ xlon(indP:indL) /)
  if (indP.gt.0) then
     tlon(n+1:) = (/ xlon(0:indP-1)/)
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
  copy_VarAtts   (xlon,tlon)

  temp&$lonName$ = tlon           ; new coord variable

  temp@lonPivot  = "reordered via lonPivot [NCL]: pivotLon="+pivotLon
  return (temp)
end

