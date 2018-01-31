C NCLFORTSTART
      subroutine filter (nx,ny,nz,a,am,f,nf)
      integer nx,ny,nz,nf
      real a(nx,ny,nz), am, f(nf)
C NCLEND
      real b(nx,ny,nz)
      integer i,j,k,k1,k2,nfh

      nfh = nf/2
      do k = 1, nz
         k1 = k-nfh
         k2 = k+nfh
         if (k1 .ge. 1 .and. k2 .le. nz) then
            do j = 1, ny
            do i = 1, nx
               b(i,j,k) = 0.0
            enddo
            enddo
            do j = 1, ny
            do i = 1, nx
               do n = 0, nf-1
                  if (a(i,j,k1+n) .ne. am) then
                     b(i,j,k) = b(i,j,k) + f(n+1)*a(i,j,k1+n)
                  else
                    b(i,j,k) = am
                    exit
                  endif
               enddo
            enddo
            enddo
         else
            do j = 1, ny
            do i = 1, nx
               b(i,j,k) = am
            enddo
            enddo
         endif
      enddo
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         a(i,j,k) = b(i,j,k)
      enddo
      enddo
      enddo
      return
      end subroutine filter
