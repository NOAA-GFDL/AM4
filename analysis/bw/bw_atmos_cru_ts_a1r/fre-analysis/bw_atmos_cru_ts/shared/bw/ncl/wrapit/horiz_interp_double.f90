
subroutine horiz_interp_double ( nxi, nyi, nk, nl, blon_in, blat_in, data_in, &
                          nxo, nyo, blon_out, blat_out, data_out, &
                          critical_frac, missing_value, verbose, error_status )
implicit none

integer, intent(in)  :: nxi, nyi, nk, nl, nxo, nyo
real*8, intent(in)   :: blon_in(2,nxi) , blat_in(2,nyi)
real*8, intent(in)   :: blon_out(2,nxo), blat_out(2,nyo)
real*8, intent(in)   :: data_in(nxi,nyi,nk,nl)
real*8, intent(out)  :: data_out(nxo,nyo,nk,nl)
real*8, intent(in)   :: critical_frac, missing_value
integer, intent(in)  :: verbose
integer, intent(out) :: error_status

!-----------------------------------------------------------------------
 
include 'horiz_interp.inc'

!-----------------------------------------------------------------------

 end subroutine horiz_interp_double

