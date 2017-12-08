#ifdef CLUBB
!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module advance_sclrm_Nd_module

  ! Description:
  ! Contains the advance_sclrm_Nd_module scheme (droplet number concentration prediction Nd).

  ! References:
  ! None
  !-----------------------------------------------------------------------

    use clubb_precision, only:  & 
      time_precision ! Variable(s)

    use grid_class, only: & 
        gr,       & ! Variable(s)
        zm2zt, & ! Procedure(s)
        zt2zm      ! Procedure(s)

   use parameters_model, only: &
      sclr_dim, cloud_frac_min

  implicit none
   ! Parameter Constants
    real, parameter, private :: & 
    epsilon_Nd = 1.0e-12

  public  :: advance_sclrm_Nd_semi_implicit, &
             advance_sclrm_Nd_upwind, &
             advance_sclrm_Nd_diffusion_OG

  private :: sclrm_Nd_lhs, &
             sclrm_Nd_rhs, &
             sclrm_Nd_solve


  contains


subroutine advance_sclrm_Nd_diffusion_OG( dt, sclrm,  sclrm_trsport_only,  Kh_zm, cf_host, err_code) 
   ! Description: droplet number diffusion scheme
   ! Advance the droplet number concentration term by one timestep using diffusion scheme.

    ! References:
    ! Ovtchinnikov and Ghan (JGR, 2005)
     real(kind=time_precision), intent(in) ::  & 
      dt            ! Current timestep size    [s]

     real, intent(inout), dimension( gr%nz , sclr_dim ) :: &
      sclrm    ! Passive scalar forcing.        [{units vary}/s]    
  
     real, intent(INout), dimension( gr%nz , sclr_dim ) :: & 
      sclrm_trsport_only    ! Passive scalar concentration due to pure transport [{units vary}/s]
      
     real, intent(in), dimension( gr%nz ) :: &
       Kh_zm ! diffusion coefficient at momentum levels 
 
      real, intent(in), dimension( gr%nz ) :: & 
       cf_host   ! cloud fraction at thermo-dynamic levels    



      integer, intent(out) :: &
       err_code       ! clubb_singular_matrix when matrix is singular  
  
      integer k, km1, kp1, nmax

      integer  i_sclr_dim

      real, dimension( gr%nz-1 )   :: a1, b1, c1, a2, c2
      real, dimension( gr%nz-1 )   :: aa, bb, cc, qq
      real, dimension( gr%nz )      :: Kdiff
      real, dimension( gr%nz-1 )   :: cf, nm

 !--------------------------- Begin Code ------------------------------------
 ! initialization
         a1 = 0.0
         b1 = 0.0
         c1 = 0.0

         a2 = 0.0
         c2 = 0.0

         aa = 0.0
         bb = 0.0
         cc  = 0.0

         qq = 0.0
      
      nmax = gr%nz-1
      do k=1,   nmax           
         cf( k )     = cf_host( k+1 )
         Kdiff( k ) = Kh_zm( k)
      enddo

      do   i_sclr_dim = 1 ,  sclr_dim
        do k=1,   nmax           
           nm( k )   = sclrm( k+1, i_sclr_dim)
        enddo

        do k=1,nmax
           km1 = max(k-1,1)
           kp1 = min(k+1,nmax)
           a1(k) = - Kdiff(k) * gr%invrs_dzt(k+1) * gr%invrs_dzm(k) * dt
           b1(k) =   Kdiff(k) * gr%invrs_dzt(k+1) * gr%invrs_dzm(k) * dt &
                       +Kdiff(k+1) * gr%invrs_dzt(k+1) *gr%invrs_dzm(k+1)  * dt
           c1(k) = - Kdiff(k+1) * gr%invrs_dzt(k+1) *gr%invrs_dzm(k+1)  * dt

           if ( cf(km1) > cloud_frac_min ) then
             a2(k) = max(cf(km1)-cf(k),0.0)/cf(km1) &
                       * Kdiff(k) * gr%invrs_dzt(k+1) * gr%invrs_dzm(k) * dt
           else
             a2(k) = 0.0
           end if

           if ( cf(kp1) > cloud_frac_min ) then
             c2(k) = max(cf(kp1)-cf(k),0.0)/cf(kp1) &
                       * Kdiff(k+1) * gr%invrs_dzt(k+1) *gr%invrs_dzm(k+1)  * dt
           else
             c2(k) = 0.0 
           end if
           
        end do

!calculate the scalar concentration due to the pure transport/diffusion (supposed to conserved)           
        aa = a1
        bb = 1.0 + b1
        cc = c1
        call tridiag(nmax,aa,bb,cc,nm,qq)
        do k=1,   nmax
            sclrm_trsport_only( k+1, i_sclr_dim)   =   sclrm_trsport_only( k+1, i_sclr_dim)  &
                                                                         + nm( k ) -  sclrm( k+1 , i_sclr_dim)
        enddo
        sclrm_trsport_only( 1 , i_sclr_dim)  = sclrm_trsport_only( 2 , i_sclr_dim)  
     

!calculate the scalar concentration due to both pure transport/diffusion and evaporation           
        aa = a1 + a2
        bb = 1.0 + b1
        cc = c1 + c2
        do k=1,   nmax
            nm( k )   = sclrm( k+1, i_sclr_dim)
        enddo
        call tridiag(nmax,aa,bb,cc,nm,qq)
        do k=1,   nmax         
          sclrm( k+1, i_sclr_dim)   = nm( k )         
        enddo
        sclrm( 1, i_sclr_dim)       = sclrm( 2, i_sclr_dim)

      enddo !  i_sclr_dim
      err_code = 0

      return
endsubroutine advance_sclrm_Nd_diffusion_OG
!============================================================

!------------------------------------------------------------------------------------------
! Triadiagonal solver from Durran (1999) pages 440-441

       subroutine tridiag(jmx,a,b,c,f,q)
       implicit none

       real a(*),b(*),c(*),f(*),q(*),p
       integer j,jmx

       c(jmx)=0.

!      Forward elimination step

       q(1)=-c(1)/b(1)
       f(1)= f(1)/b(1)

       do j=2,jmx
         p= 1.0/( b(j)+a(j)*q(j-1) )
         q(j)= -c(j)*p
         f(j)=( f(j)-a(j)*f(j-1) )*p
       end do

!      Backward pass

       do j=jmx-1,1,-1
         f(j)=f(j)+q(j)*f(j+1)
       end do

       return
 endsubroutine tridiag
!============================================================

  
  !==========================================================
  subroutine advance_sclrm_Nd_upwind( dt, sclrm, rcm, wprcp_zm, cf, err_code ) 
   ! Description: upwind scheme
   ! Advance the droplet number concentration term by one timestep using upwind scheme.

    ! References:
    ! wpncp = wprcp *Nd /(rc + epsilon_Nd)
     real(kind=time_precision), intent(in) ::  & 
      dt            ! Current timestep size    [s]

     real, intent(inout), dimension( gr%nz , sclr_dim ) :: & 
      sclrm    ! Passive scalar forcing.        [{units vary}/s]

     real, intent(in), dimension( gr%nz ) :: & 
      rcm    ! cloud water content.        [ kg/kg ]

     real, intent(in), dimension( gr%nz ) :: & 
      wprcp_zm   ! w'rc'. at momentum levels       [ (kg m) /(kg s) ]

     real, intent(in), dimension( gr%nz ) :: & 
      cf   ! cloud fraction at thermo-dynamic levels    

     integer, intent(out) :: &
      err_code       ! clubb_singular_matrix when matrix is singular


     ! Local variables
     real,  dimension( gr%nz ) :: & 
       flux_nd, Nc_cld
 
     integer :: k, km1, i_sub,  n_sub  ! Array indices

     real :: dt_sub
     real :: sum_sclrm
 !--------------------------- Begin Code ------------------------------------
         dt_sub = dt
         do k = 2, gr%nz-1, 1
            if( rcm(k) .gt. 1.0e-6 ) then
                dt_sub = min( dt_sub, &
                   abs(1.0/gr%invrs_dzm(k)/( wprcp_zm(k) + 0.0  )*rcm(k)) )
            endif
         enddo
!         print*, ' dt_sub = ', dt_sub, dt

         n_sub = floor(dt/dt_sub)
         nc_cld(:) = sclrm( : , 1 )/(cf(:)+epsilon_Nd)

          do i_sub =1, n_sub
! initialization of flux_nd
           flux_nd = 0.0
           do k = 2, gr%nz-1
! use upwind scheme
             if( rcm(k) .gt. 1.0e-6 ) then
                if( wprcp_zm(k) >= 0.0)    then
                      flux_nd(k) = wprcp_zm(k) * sclrm( k, 1 )/(rcm(k)+epsilon_Nd)
                      flux_nd(k) =  flux_nd(k) * cf(k)
                endif ! wprcp_zm(k) > 0.0

                if( wprcp_zm(k) < 0.0) then
                      flux_nd(k) = wprcp_zm(k) * sclrm( k+1, 1 )/(rcm(k+1)+epsilon_Nd)
                      flux_nd(k) =  flux_nd(k) * cf(k)
                endif ! wprcp_zm(k) < 0.0
             endif ! rcm(k) .gt. 1.0e-6

           enddo ! k = 2, gr%nz-1, 1

           do k = 2, gr%nz-1
                sclrm( k, 1 ) = sclrm( k, 1 ) &
                                       - dt_sub * 1.0 &
                                        * ( flux_nd(k) - flux_nd(k-1) ) &
					  *gr%invrs_dzt(k)
           enddo
           sclrm( 1, 1 )              = sclrm( 2, 1 )
           sclrm( gr%nz, 1 )  = sclrm( gr%nz-1, 1 )

         enddo ! do i_sub =1, n_sub


         write(*,100) ' k ', 'h (km)', 'cf', 'wprcp',  'rcm', ' w_hat ', 'Nd_cld', ' flux ',  'd_flux' 

         flux_nd = 0.0
         do k = 2, gr%nz-1, 1
! use upwind scheme
             if( rcm(k) .gt. 1.0e-6 ) then
                if( wprcp_zm(k) >= 0.0)    then
                      flux_nd(k) = wprcp_zm(k) * sclrm( k, 1 )/(rcm(k)+epsilon_Nd)
                      flux_nd(k) =  flux_nd(k) * cf(k)
                endif ! wprcp_zm(k) > 0.0

                if( wprcp_zm(k) < 0.0) then
                      flux_nd(k) = wprcp_zm(k) * sclrm( k+1, 1 )/(rcm(k+1)+epsilon_Nd)
                      flux_nd(k) =  flux_nd(k) * cf(k)
                endif ! wprcp_zm(k) < 0.0
             endif ! rcm(k) .gt. 1.0e-6

         enddo ! k = 2, gr%nz-1, 1

         do k = 2, gr%nz-1

            if(k>12 .and. k<43)   then
               write(*, 101)  k , gr%zt(k)*0.001,  cf(k),  wprcp_zm(k),  (rcm(k)+epsilon_Nd)*1.0e6, &
                 wprcp_zm(k)/(rcm(k)+epsilon_Nd),  sclrm( k, 1 )/(cf(k)+epsilon_Nd)*1.0e-6, &
                 flux_nd(k)*1.0e-6, -( flux_nd(k) - flux_nd(k-1) )*1.0e-6
            endif
100      format (A4, 2x, A10,  2x,  A12,  8(5x,  A18) )
101      format (I4, 2x,  f10.2, 2x, f12.6, 8(2x, 1pe15.5) )

              sclrm( k, 1 ) = sclrm( k, 1 ) &
                                    - (dt - n_sub*dt_sub) * 1.0 &
                                        * ( flux_nd(k) - flux_nd(k-1) ) &
					  *gr%invrs_dzt(k)

         enddo
         sclrm( 1, 1 )              = sclrm( 2, 1 )
         sclrm( gr%nz, 1 )  = sclrm( gr%nz-1, 1 )

         err_code = 0

         return
endsubroutine advance_sclrm_Nd_upwind
 !============================================================


  
 !=============================================================
  subroutine advance_sclrm_Nd_semi_implicit( dt, sclrm, rcm, wprcp_zm, err_code )
  
   ! Description: semi-implicit scheme
   ! Advance the droplet number concentration term by one timestep using semi-implicit scheme.

    ! References:
    ! wpncp = wprcp *Nd /(rc + epsilon_Nd)
     real(kind=time_precision), intent(in) ::  & 
      dt            ! Current timestep size    [s]

     real, intent(inout), dimension( gr%nz , sclr_dim ) :: & 
      sclrm    ! Passive scalar forcing.        [{units vary}/s]

     real, intent(in), dimension( gr%nz ) :: & 
      rcm    ! cloud water content.        [ kg/kg ]

     real, intent(in), dimension( gr%nz ) :: & 
      wprcp_zm    ! w'rc' at momentum levels.        [ (kg m) /(kg s) ]

    real, dimension( 3 , gr%nz ) :: &
      lhs ! The implicit part of the tridiagonal matrix              [units vary]

    real, dimension( gr%nz ) :: &
      rhs,     &! The explicit part of the tridiagonal matrix        [units vary]
      solution  ! The solution to the tridiagonal matrix             [units vary]

    integer, intent(out) :: &
      err_code       ! clubb_singular_matrix when matrix is singular

    ! Local variables
    real :: rcond ! Estimate of the reciprocal of the condition number on the LHS matrix

   !--------------------------- Begin Code ------------------------------------

      call sclrm_Nd_lhs(dt,  rcm, wprcp_zm,  &   !  Intent (In)
                                 lhs)                            !  Intent (out)

      call sclrm_Nd_rhs(dt,  sclrm, rcm, wprcp_zm,  &   !  Intent (In)
                        rhs)                            !  Intent (out)

       call sclrm_Nd_solve(  lhs, rhs, &    !  Intent (In)
                   solution, err_code, rcond )  !  Intent (out)

        sclrm( 1:gr%nz, 1 ) = solution(1:gr%nz)

 return
endsubroutine advance_sclrm_Nd_semi_implicit
 !=============================================================================




  subroutine sclrm_Nd_lhs(dt,  rcm, wprcp_zm,  lhs)
   ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      dt               ! Timestep                                 [s]

    real, intent(in), dimension(gr%nz) :: & 
      rcm,      &  !cloud water content at thermo-dynamic levels [kg/kg]
      wprcp_zm       ! w'r_c' at momentum levels [(kg/kg) m/s]

! intermediate arrays
    real, dimension(gr%nz) :: & 
      rcm_zm     !cloud water content at momentum levels [kg/kg]

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2,    & ! Thermodynamic main diagonal index.
      km1_tdiag = 3       ! Thermodynamic subdiagonal index.

    ! Output Variable
    real, dimension(3,gr%nz), intent(out) :: &
      lhs           ! Implicit contributions to xm (tridiagonal matrix)

    ! Local Variables
    integer :: k, km1  ! Array indices

    ! --- Begin Code ---

    rcm_zm = zt2zm( rcm )

    ! Initialize the LHS array.
    lhs = 0.0

    do k = 2, gr%nz-1, 1
! Define index
      km1 = max( k-1, 1 )

     ! LHS time tendency.
      lhs( k_tdiag, k )  &
      = real( lhs(k_tdiag,k) + ( 1.0 / dt ) )

      lhs( k_tdiag, k )  &
      = real( lhs( k_tdiag, k ) + &
      0.5 * wprcp_zm(k) * ( gr%zt(k+1)- gr%zm(k) ) / &
      ( gr%zm(k) - gr%zm(km1) )/( gr%zt(k+1)- gr%zt(k) ) / &
      (  rcm_zm(k) + epsilon_Nd ) )

      lhs( k_tdiag, k )  &
      = real( lhs( k_tdiag, k ) - &
      0.5 * wprcp_zm(km1) * ( gr%zm(km1)- gr%zt(km1) ) / &
      ( gr%zm(k) - gr%zm(km1) )/( gr%zt(k)- gr%zt(km1) ) / &
      (  rcm_zm(km1) + epsilon_Nd ) )

      lhs( kp1_tdiag, k )  &
      = real( lhs( kp1_tdiag, k ) - &
      0.5 * wprcp_zm(km1) * ( gr%zt(k)- gr%zm(km1) ) / &
      ( gr%zm(k) - gr%zm(km1) )/( gr%zt(k)- gr%zt(km1) ) / &
      (  rcm_zm(km1) + epsilon_Nd ) )

      lhs( km1_tdiag, k )  &
      = real( lhs( km1_tdiag, k ) + &
      0.5 * wprcp_zm(k) * ( gr%zm(k) -  gr%zt(k) ) / &
      ( gr%zm(k) - gr%zm(km1) )/( gr%zt(k+1)- gr%zt(k) ) / &
      (  rcm_zm(k) + epsilon_Nd ) )

    enddo ! k = 2 .. gr%nz-1

    ! Lower Boundary
     k = 1
     lhs( k_tdiag, k ) = 1.0

    ! Upper Boundary
     k = gr%nz
     lhs( k_tdiag, k ) = 1.0

  return
  end subroutine sclrm_Nd_lhs
!=============================================================================



  subroutine sclrm_Nd_rhs(dt,  sclrm, rcm, wprcp_zm,  rhs)
   ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      dt               ! Timestep                                 [s]

   ! Input Variables
    real, intent(in), dimension(gr%nz,sclr_dim) ::  & 
      sclrm

    real, intent(in), dimension(gr%nz) :: & 
      rcm,      &  !cloud water content at thermo-dynamic levels [kg/kg]
      wprcp_zm       ! w'r_c'  at thermo-dynamic levels  [(kg/kg) m/s]

! intermediate arrays
    real, dimension(gr%nz) :: & 
      rcm_zm     !cloud water content at momentum levels [kg/kg]

    ! Output Variable
    real, dimension(gr%nz), intent(out) :: &
      rhs           ! Right-hand side of band diag. matrix. (LAPACK)

    ! Local Variables
    integer :: k, km1  ! Array indices

    ! --- Begin Code ---

    rcm_zm = zt2zm( rcm )

    ! Initialize the LHS array.
    rhs = 0.0

    do k = 2, gr%nz-1, 1
! Define index
      km1 = max( k-1, 1 )

       rhs( k )   =   rhs( k ) &
        + sclrm( k , 1 )/dt

       rhs( k )  &
      =  rhs( k )  &
        - 0.5 * wprcp_zm(k) * ( gr%zt(k+1)- gr%zm(k) ) * sclrm( k , 1 ) / &
      ( gr%zm(k) - gr%zm(km1) )/( gr%zt(k+1)- gr%zt(k) ) / &
      (  rcm_zm(k) + epsilon_Nd )
      
       rhs( k )  &
      =  rhs( k )  &
        + 0.5 * wprcp_zm(km1) * ( gr%zm(km1)- gr%zt(km1) ) * sclrm( k , 1 )/ &
      ( gr%zm(k) - gr%zm(km1) )/( gr%zt(k)- gr%zt(km1) ) / &
      (  rcm_zm(km1) + epsilon_Nd )
  

       rhs( k )  &
      =  rhs( k )  &
         - 0.5 * wprcp_zm(k) * ( gr%zm(k) -  gr%zt(k) ) * sclrm( k+1 , 1 )/ &
      ( gr%zm(k) - gr%zm(km1) )/( gr%zt(k+1)- gr%zt(k) ) / &
      (  rcm_zm(k) + epsilon_Nd )


      rhs( k )  &
      =  rhs( k )  &
         +0.5 * wprcp_zm(km1) * ( gr%zt(k)- gr%zm(km1) ) * sclrm( k-1 , 1 )/ &
      ( gr%zm(k) - gr%zm(km1) )/( gr%zt(k)- gr%zt(km1) ) / &
      (  rcm_zm(km1) + epsilon_Nd )
    enddo ! k = 2 .. gr%nz-1

    ! Boundary Conditions

    ! Lower Boundary
     k = 1
     rhs( k ) = sclrm( k , 1 )

    ! Upper Boundary
     k = gr%nz
     rhs( k ) = sclrm( k , 1 )

  return
  end subroutine sclrm_Nd_rhs
!===============================================================================



 subroutine sclrm_Nd_solve(  lhs, rhs, solution, err_code, rcond )
    ! Description:
    !   Solve for droplet number concentration (Nd) using the band diagonal solver.
    !   A Crank-Nicholson time-stepping algorithm
    ! is used in solving the turbulent advection term.

    ! The rate of change of an eddy-scalar quantity, xm, is:
    !
    ! d(Nd)/dt = - d(Nd'w')/dz 
    !  Nd'w'= rc'w'  *Nd /(rc + epsilon_Nd)

    ! Boundary Conditions:
    !
    ! An eddy-scalar quantity is not allowed to flux out the upper boundary.
    ! Thus, a zero-flux boundary condition is used for the upper boundary in the
    ! eddy-diffusion equation. 
    ! The lower boundary condition is a 
    ! zero-flux  boundary condition.
    !------------------------------------------------------------------------

    use grid_class, only: & 
      gr ! Variable(s)

     use lapack_wrap, only:  & 
      tridag_solve, & ! Procedure(s)
      tridag_solvex

    implicit none

    ! Constant parameters

    integer, parameter :: &
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2,    & ! Thermodynamic main diagonal index.
      km1_tdiag = 3       ! Thermodynamic subdiagonal index.

   ! Input/Output Variables
    real, dimension( 3 , gr%nz), intent(inout) :: & 
      lhs  ! Implicit contributions to Nd (band diag. matrix in LAPACK storage)

    real, dimension( gr%nz ), intent(inout) ::  & 
      rhs      ! Right-hand side of band diag. matrix. (LAPACK storage)

    real, dimension( gr%nz ), intent(out) ::  & 
      solution ! Solution to band diagonal system (LAPACK storage)

    integer, intent(out) :: & 
      err_code ! clubb_singular_matrix when matrix is singular

    ! Local variables
    real :: rcond ! Estimate of the reciprocal of the condition number on the LHS matrix
    !-----------------------------------------------------------------------

      call tridag_solvex( "Solve_Nd_HOC", gr%nz, 1, &                        !  Intent (In)
                         lhs(kp1_tdiag,:),  lhs(k_tdiag,:), lhs(km1_tdiag,:), rhs, & !  Intent(Inout)
                         solution, rcond, err_code )                                            !  Intent(Out)
 

 return
 end subroutine sclrm_Nd_solve
!===============================================================================

end module advance_sclrm_Nd_module
#endif
