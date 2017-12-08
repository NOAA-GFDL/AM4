#ifdef CLUBB

module CLUBB_3D_var

! This module contains definitions of 3D arrays for 
! higher order terms
!-----------------------------------------------------------------------

use parameters_model, only: & 
       sclr_dim

    use constants_clubb, only:  & 
        rt_tol, &
        thl_tol, &
        w_tol_sqd

implicit none
public :: setup_CLUBB_3D_var, &
              cleanup_CLUBB_3D_var
	   
real, allocatable, dimension(:,:,:,:), public :: &	  ! save high-res results 
           RH_crit_clubb_3D     ! (spatial 3D) critical relative humidity for droplet and ice nucleation

real, allocatable, dimension(:,:,:), public :: &   !higher order terms
           upwp_3D,            &
           vpwp_3D,            &
           up2_3D,              &
           vp2_3D,              &
           wp2_3D,             &
           wprtp_3D,           &
           wpthlp_3D,          &
           rtp2_3D,              &
           thlp2_3D,            &
           rtpthlp_3D,          &
           wp3_3D

contains
!-----------------------------------------------------------------------
!  Allocates higher order terms
!  for the HOC model code
!-----------------------------------------------------------------------
subroutine setup_CLUBB_3D_var(id, jd, nzmax )

implicit none
! Input Variables
integer, intent(in) :: id, jd, nzmax
           allocate(   upwp_3D(1:id, 1:jd, 1:nzmax)    )
           allocate(   vpwp_3D(1:id, 1:jd, 1:nzmax)    )
           allocate(   up2_3D(1:id, 1:jd, 1:nzmax)   )
           allocate(   vp2_3D(1:id, 1:jd, 1:nzmax)   )
           allocate(   wp2_3D(1:id, 1:jd, 1:nzmax)    )
           allocate(   wprtp_3D(1:id, 1:jd, 1:nzmax)    )
           allocate(   wpthlp_3D(1:id, 1:jd, 1:nzmax)   )
           allocate(   rtp2_3D(1:id, 1:jd, 1:nzmax)   )
           allocate(   thlp2_3D(1:id, 1:jd, 1:nzmax)   )
           allocate(   rtpthlp_3D(1:id, 1:jd, 1:nzmax)   )
           allocate(   wp3_3D(1:id, 1:jd, 1:nzmax)    )

          if( sclr_dim > 0) then
              allocate(   RH_crit_clubb_3D(1:id, 1:jd, 1:nzmax, 2 ) )
              RH_crit_clubb_3D = 1.0
           endif

           wp2_3D       = w_tol_sqd * 800
           wp3_3D       = 0.0
           upwp_3D      = 0.0
           vpwp_3D      = 0.0
           wprtp_3D     = 0.0
           wpthlp_3D    = 0.0
           rtp2_3D        = rt_tol**2
           thlp2_3D       = thl_tol**2
           rtpthlp_3D     = 0.0
           up2_3D         = w_tol_sqd * 800
           vp2_3D         = w_tol_sqd * 800

return
end subroutine setup_CLUBB_3D_var


!------------------------------------------------------------------------
subroutine cleanup_CLUBB_3D_var( )
!       Description:
!       Subroutine to deallocate variables defined in module global
!------------------------------------------------------------------------
implicit none

! --- Deallocate --- 
           deallocate(   upwp_3D    )
           deallocate(   vpwp_3D    )
           deallocate(   up2_3D   )
           deallocate(   vp2_3D   )
           deallocate(   wp2_3D    )
           deallocate(   wprtp_3D    )
           deallocate(   wpthlp_3D   )
           deallocate(   rtp2_3D   )
           deallocate(   thlp2_3D   )
           deallocate(   rtpthlp_3D   )
           deallocate(   wp3_3D    )
 
            if( sclr_dim > 0) then
               deallocate(  RH_crit_clubb_3D )
           endif

return
end subroutine cleanup_CLUBB_3D_var

end module CLUBB_3D_var
#endif
