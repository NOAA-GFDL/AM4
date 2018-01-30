C NCLFORTSTART
      subroutine average3d (nx,ny,nt,a,am,md,nm,b)
      integer nx,ny,nt,nm,md(nm)
      real a(nx,ny,nt), b(nx,ny,nm), am
C NCLEND
      integer i,j,k,ks,ke,ka,num(nx,ny)

      ke = 0
      do k = 1, nm

         do j = 1, ny
         do i = 1, nx
            num(i,j) = 0
            b(i,j,k) = 0.0
         enddo
         enddo

         ks = ke+1
         ke = ke+md(k)
         do ka = ks, ke
            do j = 1, ny
            do i = 1, nx
               if (a(i,j,ka) .ne. am) then
                  b(i,j,k) = b(i,j,k) + a(i,j,ka)
                  num(i,j) = num(i,j) + 1
               endif
            enddo
            enddo
         enddo
         do j = 1, ny
         do i = 1, nx
            if (num(i,j) .ge. (md(k)+1)/2) then
               b(i,j,k) = b(i,j,k)/float(num(i,j))
            else
               b(i,j,k) = am
            endif
         enddo
         enddo
      enddo

      return
      end subroutine average3d
