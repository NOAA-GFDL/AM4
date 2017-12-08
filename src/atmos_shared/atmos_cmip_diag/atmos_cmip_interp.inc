
        if (logp <= pm(1)) then 
          if (extrap) then
            ap(i,j,kp) = a(i,j,1) + (a(i,j,3)-a(i,j,1)) * (logp-pm(1))/(pm(3)-pm(1))
          else
            ap(i,j,kp) = a(i,j,1)
          endif
        else if (logp >= pm(km)) then 
          if (extrap) then
            ap(i,j,kp) = a(i,j,km) + (a(i,j,km)-a(i,j,km-2)) * (logp-pm(km))/(pm(km)-pm(km-2))
          else
            ap(i,j,kp) = a(i,j,km)
          endif
        else 
          do k = 1, km-1
            if (logp <= pm(k+1) .and. logp >= pm(k)) then 
              ap(i,j,kp) = a(i,j,k) + (a(i,j,k+1)-a(i,j,k)) * (logp-pm(k))/(pm(k+1)-pm(k))
              go to 1000
            endif
          enddo
        endif
1000  continue

