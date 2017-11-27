C NCLFORTSTART
      subroutine average1d (nt,a,md,nm,b)
      integer nt, nm, md(nm)
      real a(nt), b(nm)
C NCLEND
      integer k,ks,ke,ka,num

      ke = 0
      do k = 1, nm
         b(k) = 0.0
         ks = ke+1
         ke = ke+md(k)
         do ka = ks, ke
            b(k) = b(k) + a(ka)
         enddo
         b(k) = b(k)/float(md(k))
      enddo

      return
      end subroutine average1d
