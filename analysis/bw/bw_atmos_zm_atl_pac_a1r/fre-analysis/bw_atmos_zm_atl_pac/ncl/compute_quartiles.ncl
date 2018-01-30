;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; returns quartile values from sorted index list
; input dimensions:  runs x season x points

undef("quartile_values")
function quartile_values (is[*][*][*]:integer,x[*][*][*]:numeric)
local dimx,dimz,j,q,w,y,z
begin
   dimx = dimsizes(x) ; # of runs is first dimension
   runs = dimx(0)
   dimz = dimx
   dimz(0) = 5
   z = new(dimz,typeof(x))
   L = new(3,integer)
   y = new(3,float)
   w = new(3,float)

   q = (/.25,.50,.75/) ; min,max are first,last quartiles
   y = (runs-1)*q
   L = toint(y)
   w = y-L 

   do j = 0, dimx(2)-1
   do k = 0, dimx(1)-1
      LS0 = is(0,k,j)
      LS  = is(L,k,j)
      LS1 = is(L+1,k,j)
      LSX = is(runs-1,k,j)
      z(0,k,j) = x(LS0,k,j) ; minimum
      z(1:3,k,j) = x(LS,k,j)*(1.-w) + x(LS1,k,j)*w
      z(4,k,j) = x(LSX,k,j) ; maximum
   end do
   end do
   return z
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; compute quartile values across runs
; input:  x(runs,:,:,:)
; output: q(5,:)
;         for minimum,1/4,median,3/4,maximum

undef("compute_quartiles")
function compute_quartiles (regions[*]:string, xRegions:logical )
local nRegion,qOut,n,x,isort,xsort
begin
  ; determine number of regions
  nRegion = dimsizes(regions)

  qOut = True

  ; loop thru regions
  do n = 0, nRegion-1
    print("region: "+regions(n))
    x = xRegions@$regions(n)$
    isort = dim_pqsort_n(x,1,0)
    xsort = quartile_values(isort,x)  ; 5 x ssn x npts
    qOut@$regions(n)$ = xsort
    delete([/isort,x,xsort/])
  end do

  return qOut
end

