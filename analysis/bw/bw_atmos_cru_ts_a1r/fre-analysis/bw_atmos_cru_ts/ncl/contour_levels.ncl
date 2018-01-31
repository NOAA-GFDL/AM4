
undef("symetric_levels")
function symetric_levels ( values[*]:numeric )
local allvalues
begin
   if (abs(values(0)) .lt. 1.e-5) then
      allvalues = array_append_record( -values(1::-1), values, 0 )
   else
      allvalues = array_append_record( -values(::-1), values, 0 )
   end if
   return allvalues
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("contour_levels")
function contour_levels ( var:string, stat:string, plotType:string )
local cnLevels, statType, levels
begin

  statType = str_lower(stat)
  cnLevels = True

  if ( statType .eq. "mean" .and. plotType .eq. "field" ) then
    cnLevels@pr = (/0.5,1,2,3,4,5,6,7,9,11,13,15/)  ; from Klein plots
    cnLevels@tas = symetric_levels((/2,6,10,14,18,24,30/))
    cnLevels@tasmean = cnLevels@tas
    cnLevels@tasmax = cnLevels@tas
    cnLevels@tasmin = cnLevels@tas
    cnLevels@tasdiff = (/2,4,6,8,10,12,14,16/)
    cnLevels@clt = (/10,20,30,40,50,60,70,80,90/)
    cnLevels@wet = cnLevels@clt

  else if ( statType .eq. "stdev" .and. plotType .eq. "field" ) then
    cnLevels@pr = (/0.25,0.5,0.75,1.0,1.5,2,3,4/)
    cnLevels@tas = (/0.25,0.5,0.75,1.0,1.25,1.5,2.0,2.5,3.0/)
    cnLevels@tasmean = cnLevels@tas
    cnLevels@tasmax = cnLevels@tas
    cnLevels@tasmin = cnLevels@tas
    cnLevels@tasdiff = cnLevels@tas
    cnLevels@clt = (/1,2,3,4,5,7,10,15/)
    cnLevels@wet = (/1,2,3,4,5,7,10,15,20/)

  else if ( statType .eq. "trend" .and. plotType .eq. "field" ) then
    cnLevels@pr = symetric_levels((/0.5,1,2,3,4,5/))
    cnLevels@tas = symetric_levels((/0.5, 1.0, 1.5, 2.0, 3.0, 4.0/))
    cnLevels@tasmean = cnLevels@tas
    cnLevels@tasmax = cnLevels@tas
    cnLevels@tasmin = cnLevels@tas
    cnLevels@tasdiff = cnLevels@tas
    cnLevels@clt = symetric_levels((/2,5,15,20,25/))
    cnLevels@wet = symetric_levels((/2,5,15,20,25/))

  else if ( statType .eq. "mean" .and. plotType .eq. "diff" ) then
    cnLevels@pr = symetric_levels((/0.5,1,2,3,5,7/))
    cnLevels@tas = symetric_levels((/1,2,4,6,8,10/))
    cnLevels@tasmean = cnLevels@tas
    cnLevels@tasmax = cnLevels@tas
    cnLevels@tasmin = cnLevels@tas
    cnLevels@tasdiff = cnLevels@tas
    cnLevels@clt = symetric_levels((/5,10,15,20,30,40/))
    cnLevels@wet = symetric_levels((/5,10,15,20,25,30/))

  else if ( statType .eq. "stdev" .and. plotType .eq. "diff" ) then
    cnLevels@pr = symetric_levels((/0.10,0.25,0.5,1,2,4/))
    cnLevels@tas = symetric_levels((/0.25,0.5,1,2,4/))
    cnLevels@tasmean = cnLevels@tas
    cnLevels@tasmax = cnLevels@tas
    cnLevels@tasmin = cnLevels@tas
    cnLevels@tasdiff = symetric_levels((/0.25,0.5,1.0,1.5,2,3/))
    cnLevels@clt = symetric_levels((/1,2,4,6,8,10/))
    cnLevels@wet = symetric_levels((/1,2,4,6,10,15/))

  else if ( statType .eq. "trend" .and. plotType .eq. "diff" ) then
    cnLevels@pr = symetric_levels((/0.25,0.5,1,2,4/))
    cnLevels@tas = symetric_levels((/0.5,1.0,1.5,2,3,4,5/))
    cnLevels@tasmean = cnLevels@tas
    cnLevels@tasmax = cnLevels@tas
    cnLevels@tasmin = cnLevels@tas
    cnLevels@tasdiff = cnLevels@tas
    cnLevels@clt = symetric_levels((/5,10,15,20,25/))
    cnLevels@wet = symetric_levels((/5,10,15,20,25/))

  end if
  end if
  end if
  end if
  end if
  end if

  if (isatt(cnLevels,var)) then
    levels = cnLevels@$var$
    return levels
  else
    print("ERROR: no contour levels found, invalid field: "+var+", statType: "+statType+", plotType: "+plotType)
    status_exit(1)
  end if
end

