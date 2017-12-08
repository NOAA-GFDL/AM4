                     MODULE aerosol_cloud_mod

!------------------------------------------------------------------------
!
!           aerosol_cloud_mod retrieves the aerosol information
!           needed by the model cloud scheme
!
!-----------------------------------------------------------------------

use fms_mod,                   only :  error_mesg, FATAL, mpp_pe,   &
                                       mpp_root_pe, open_namelist_file, &
                                       check_nml_error, close_file,  &
                                       write_version_number, file_exist, &
                                       stdlog
use constants_mod,             ONLY :  grav, cp_air, rdgas, rvgas, tfreeze
use time_manager_mod,          only :  time_type
use aerosol_types_mod,         ONLY :  aerosol_type
use mpp_mod,                   only :  mpp_clock_id, mpp_clock_begin,  &
                                       mpp_clock_end, CLOCK_LOOP,  &
                                       input_nml_file
use polysvp_mod,               ONLY :  polysvp_l, polysvp_i, polysvp_init,&
                                       polysvp_end
use aerosol_params_mod,        only :  aerosol_params_init, Nfact_du1,  &
                                       Nfact_du2, Nfact_du3, Nfact_du4, &
                                       Nfact_du5, rbar_du1, rbar_du2,   &
                                       rbar_du3, rbar_du4, rbar_du5
USE  aer_ccn_act_mod,          ONLY :  aer_ccn_act_wpdf_m, &
                                       aer_ccn_act_init,  aer_ccn_act_end
USE  ice_nucl_mod,             ONLY :  ice_nucl_wpdf, ice_nucl_wpdf_init, &
                                       ice_nucl_wpdf_end
USE lscloud_types_mod,         ONLY :  lscloud_types_init, &
                                       diag_id_type, diag_pt_type, &
                                       atmos_state_type, lscloud_nml_type,&
                                       lsc_constants_type, particles_type
use lscloud_debug_mod,         only :  aerosol_cloud_debug1,   &
                                       aerosol_cloud_debug2, &
                                       lscloud_debug_init
use moist_proc_utils_mod,      only :  mp_input_type, mp_conv2ls_type, &
                                       mp_nml_type, mp_lsdiag_type, &
                                       mp_lsdiag_control_type
use physics_radiation_exch_mod, only : exchange_control_type

implicit none 
private

!--------------------------------------------------------------------------
!---interfaces-------------------------------------------------------------

public   aerosol_cloud_init, aerosol_cloud_end
public   determine_available_aerosol, determine_activated_aerosol
public   aerosol_effects

!--------------------------------------------------------------------------
!---version number---------------------------------------------------------

Character(len=128) :: Version = '$Id$'
Character(len=128) :: Tagname = '$Name$'

!--------------------------------------------------------------------------
!---namelist---------------------------------------------------------------

integer :: rh_act_opt = 2           ! option as to which rh to use for ice
                                    ! nucleation. 1 => use gridbox mean;
                                    ! 2 => use rh in non-cloudy portion of
                                    ! box when the cloudy part of grid box
                                    ! is less than cf_thresh_nucl; use grid
                                    ! box mean otherwise.
real    :: sea_salt_scale_onl = 1.  ! scaling factor to convert seasalt 
                                    ! aerosol tracer to seasalt available
                                    ! for use as condensation nuclei
logical :: reproduce_rk = .true.    ! option to force use of same entries in
                                    ! aerosol activation tables as were used
                                    ! in legacy code; new code fixes an 
                                    ! index offset error.
real    :: var_limit_ice = -999.    ! lower limit to the vertical velocity
                                    ! variance of the pdf used for ice
                                    ! nucleation; default value implies 
                                    ! that no limit is imposed.
real    :: cf_thresh_nucl = 0.98    ! threshold cloud fraction below which
                                    ! the rh of the non-cloudy portion of 
                                    ! the gribox is used in determining
                                    ! ice nucleation (when rh_act_opt is >1)
!  <DATA NAME="var_limit" UNITS=" (m**2)/(s**2)" TYPE="real" DEFAULT="0.0">
!   minimum value of the variance in the vertical velocity pdf
!  </DATA>

real    :: var_limit = 0.0
integer :: var_limit_opt = 1
integer :: up_strat_opt = 1

namelist / aerosol_cloud_nml / rh_act_opt, sea_salt_scale_onl, &
                               reproduce_rk, var_limit_ice, &
                               var_limit,  var_limit_opt, up_strat_opt, &
                               cf_thresh_nucl



!-------------------------------------------------------------------------

real, parameter    :: d622 = rdgas / rvgas
real, parameter    :: d378 = 1. - d622
integer, parameter :: n_totmass  = 4
integer, parameter :: n_imass = 12

LOGICAL            :: do_var_limit_ice
INTEGER            :: wpdf_offs

integer            :: aero_init, aero_effects, aero_dust, aero_loop1, &
                      aero_loop2, aero_loop3

logical            :: module_is_initialized = .false.
real               :: missing_value = -1.e30

real    :: qmin
logical :: do_liq_num
logical :: do_ice_num
logical :: do_ice_nucl_wpdf
logical :: do_dust_berg      
logical :: use_online_aerosol
logical :: use_sub_seasalt       
real    :: sea_salt_scale
real    :: om_to_oc
logical :: do_pdf_clouds
logical :: do_mg_microphys, do_mg_ncar_microphys, do_ncar_microphys
logical :: total_activation

logical :: debug


                             CONTAINS


!#########################################################################

subroutine aerosol_cloud_init (Constants_lsc, Nml_lsc, Nml_mp, Exch_ctrl)

!------------------------------------------------------------------------
type(lsc_constants_type),    intent(in) :: Constants_lsc
type(lscloud_nml_type),      intent(in) :: NmL_lsc
type(mp_nml_type),           intent(in) :: NmL_mp
type(exchange_control_type), intent(in) :: Exch_ctrl


!-----local variables

      integer :: unit, io, ierr, logunit

!-----------------------------------------------------------------------
      if (module_is_initialized) return

!----------------------------------------------------------------------
!    define module variables that are obtained from other modules.
!----------------------------------------------------------------------
      qmin = Exch_ctrl%qmin
      do_liq_num = Exch_ctrl%do_liq_num
      do_ice_num = Exch_ctrl%do_ice_num
      do_ice_nucl_wpdf = Nml_lsc%do_ice_nucl_wpdf
      do_dust_berg = Nml_lsc%do_dust_berg
      use_online_aerosol = Nml_mp%use_online_aerosol
      use_sub_seasalt = Nml_mp%use_sub_seasalt
      sea_salt_scale = Nml_mp%sea_salt_scale
      om_to_oc = Nml_mp%om_to_oc
      do_pdf_clouds = Nml_lsc%do_pdf_clouds
      do_mg_microphys = Constants_lsc%do_mg_microphys
      do_mg_ncar_microphys = Constants_lsc%do_mg_ncar_microphys
      do_ncar_microphys = Constants_lsc%do_ncar_microphys
      total_activation = Constants_lsc%total_activation

!-------------------------------------------------------------------------
!    process namelist.
!-------------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=aerosol_cloud_nml, iostat=io)
      ierr = check_nml_error(io,'aerosol_cloud_nml')
#else
      if ( file_exist('input.nml')) then
        unit = open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=aerosol_cloud_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'aerosol_cloud_nml')
        enddo
10      call close_file (unit)
      endif
#endif
 
!-------------------------------------------------------------------------
!    write version and namelist to standard log.
!-------------------------------------------------------------------------
      call write_version_number ( version, tagname )
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe()) &
                                  write (logunit, nml=aerosol_cloud_nml)

      if (var_limit_ice .EQ. -999. ) THEN
        do_var_limit_ice = .false.
      ELSE
        do_var_limit_ice = .true.
      END IF

!-------------------------------------------------------------------------
!    define offset into aerosol activation tables.
!-------------------------------------------------------------------------
      IF (.NOT. (do_mg_microphys .OR.   & 
                 do_ncar_microphys .OR.   & 
                 do_mg_ncar_microphys)  & 
                                       .AND. reproduce_rk) THEN
        wpdf_offs =0
      ELSE 
        wpdf_offs = 1
      END IF

!-------------------------------------------------------------------------
!    be sure needed modules have been initialized.
!-------------------------------------------------------------------------
      call lscloud_types_init
      call polysvp_init (Nml_lsc)
      call aer_ccn_act_init
      call ice_nucl_wpdf_init
      call lscloud_debug_init (debug)

!-------------------------------------------------------------------------
!    initialize clocks.
!-------------------------------------------------------------------------
      aero_init = mpp_clock_id ('aerosol_cloud:  initial',   &
                                                      grain = CLOCK_LOOP)
      aero_effects = mpp_clock_id ('aerosol_cloud:  effects',  &
                                                      grain = CLOCK_LOOP)
      aero_dust = mpp_clock_id ('aerosol_cloud:  dust   ',   &
                                                      grain = CLOCK_LOOP)
      aero_loop1 = mpp_clock_id ('aerosol_cloud:  loop1  ',   &
                                                      grain = CLOCK_LOOP)
      aero_loop2 = mpp_clock_id ('aerosol_cloud:  loop2  ',  &
                                                      grain = CLOCK_LOOP)
      aero_loop3 = mpp_clock_id ('aerosol_cloud:  loop3  ',  &
                                                      grain = CLOCK_LOOP)

!-------------------------------------------------------------------------
!    declare module initialized.
!-------------------------------------------------------------------------
      module_is_initialized = .true.


end subroutine aerosol_cloud_init



!#####################################################################

subroutine determine_available_aerosol (idim, jdim, kdim, Lsdiag_mp,      &
                                        Lsdiag_mp_control, &
                                        Atmos_state, Particles, Aerosol)

!-------------------------------------------------------------------------
integer,                    intent(in)     :: idim, jdim, kdim
type(mp_lsdiag_type),       intent(inout)  :: Lsdiag_mp
type(mp_lsdiag_control_type), intent(inout):: Lsdiag_mp_control
type(atmos_state_type),     intent(inout)  :: Atmos_state
type(particles_type),       intent(inout)  :: Particles  
TYPE(aerosol_type),         INTENT (in)    :: Aerosol  


!-------------------------------------------------------------------------
!----local variables
      INTEGER :: i, j, k

!---------------------------------------------------------------------
!    call aerosol_effects to determine aerosol fields that will impact the 
!    cloud / ice droplet numbers and the bergeron process, if these effects
!    have been activated.
!---------------------------------------------------------------------
      call mpp_clock_begin (aero_effects)
      if (do_liq_num .or. do_dust_berg) then
        call aerosol_effects (idim, jdim, kdim,    &
                              Lsdiag_mp_control%n_diag_4d,  &
                              Atmos_state%pthickness, &
                              Particles%concen_dust_sub,&
                              Particles%totalmass1, Particles%imass1, &
                              Aerosol, Lsdiag_mp%diag_4d,  & 
                              Lsdiag_mp_control%diag_id,   &
                              Lsdiag_mp_control%diag_pt       )
      endif
      call mpp_clock_end   (aero_effects)

!-------------------------------------------------------------------------
!    define the number and mean radius of dust particles present.
!-------------------------------------------------------------------------
      call mpp_clock_begin (aero_dust)
      do k = 1,kdim
        do j=1,jdim
          do i = 1,idim
            Particles%ndust(i,j,k) =    &
                       Nfact_du1*Particles%imass1(i,j,k,8) + &
                       Nfact_du2*Particles%imass1(i,j,k,9) + &
                       Nfact_du3*Particles%imass1(i,j,k,10) +&
                       Nfact_du4*Particles%imass1(i,j,k,11) +&
                       Nfact_du5*Particles%imass1(i,j,k,12)

            Particles%rbar_dust(i,j,k) = (  &
                       Nfact_du1*Particles%imass1(i,j,k,8)*rbar_du1  +  &
                       Nfact_du2*Particles%imass1(i,j,k,9)*rbar_du2  +  &
                       Nfact_du3*Particles%imass1(i,j,k,10)*rbar_du3  + &
                       Nfact_du4*Particles%imass1(i,j,k,11)*rbar_du4  + &
                       Nfact_du5*Particles%imass1(i,j,k,12)*rbar_du5 ) / &
                                       MAX (Particles%ndust(i,j,k),1.e-10)
          end do
        end do        
      end do        
      call mpp_clock_end (aero_dust)


!-------------------------------------------------------------------------


END SUBROUTINE  determine_available_aerosol

!########################################################################

subroutine determine_activated_aerosol (   &
                     idim, jdim, kdim, n_diag_4d, C2ls_mp, Input_mp,  &
                     Atmos_state, Particles, qa_upd,  &
                                               diag_4d, diag_id, diag_pt )

!-------------------------------------------------------------------------
integer,                    intent(in)     :: idim, jdim, kdim, n_diag_4d
type(atmos_state_type),     intent(inout)  :: Atmos_state
type(mp_input_type),        intent(in)     :: Input_mp
type(mp_conv2ls_type),      intent(in)     :: c2ls_mp
type(particles_type),       intent(inout)  :: Particles  
REAL, DIMENSION(:,:,:),     INTENT(IN )    :: qa_upd
REAL, DIMENSION( :,:,:,0:), INTENT(INOUT ) :: diag_4d
TYPE(diag_id_type),         intent(in)     :: diag_id
TYPE(diag_pt_type),         intent(in)     :: diag_pt


!-------------------------------------------------------------------------
!----local variables
      real, dimension(idim, jdim, kdim)             ::                 &
                                 up_strat, wp2, wp2i, wp2t, thickness, &
                                 ni_sulf, ni_dust, ni_bc, cf, u_i, u_l, &
                                 eslt, esit, qvsl, qvsi, qs_d, qvt
      INTEGER :: i, j, k


!-------------------------------------------------------------------------
!    do the following calculations if droplet number is being predicted:
!-------------------------------------------------------------------------
      if (do_liq_num) then 

!-------------------------------------------------------------------------
!    compute the relevant upward velocity for droplet / ice activation.
!-------------------------------------------------------------------------
        call mpp_clock_begin (aero_loop1)
        if (up_strat_opt == 1) then   ! cjg
          do k = 1,kdim
            do j=1,jdim
              do i = 1,idim
                up_strat(i,j,k) = -1.*(((Input_mp%omega(i,j,k) + &
                                      grav*C2ls_mp%mc_full(i,j,k))/ &
                                      Atmos_state%airdens(i,j,k)/grav) + &
                                      Input_mp%radturbten(i,j,k)* &
                                                              cp_air/grav)
              end do
            end do
          end do

        elseif (up_strat_opt == 2) then   ! cjg
          do k = 1,kdim
            do j=1,jdim
              do i = 1,idim
                up_strat(i,j,k) = -1.*( Input_mp%omega(i,j,k) &
                                     / (Atmos_state%airdens(i,j,k)*grav) )
              end do
            end do
          end do
        endif   ! cjg
        call mpp_clock_end (aero_loop1)

!-------------------------------------------------------------------------
!    define the layer thickness and variance of the vertical velocity. 
!    call aer_ccn_act_wpdf_m to determine number of activated droplets.
!    this call need not be made if the activation vertical velocity is
!    downward, the rotstayn-klein microphysics is active, and pdf_clouds
!    are not activated; in such a case, no particles are activated. 
!-------------------------------------------------------------------------
        call mpp_clock_begin (aero_loop2)
        if (var_limit_opt == 1) then   ! cjg
          do k = 1,kdim
            do j=1,jdim
              do i = 1,idim
                if ( (do_pdf_clouds) .or.  &
                   (do_mg_microphys)  .or.  &
                   (do_ncar_microphys)  .or.  &
                   (do_mg_ncar_microphys)  .or.  &
! cjg: total activation for RK
                   (total_activation) .or.  &  
                   (up_strat(i,j,k) >= 0.0) )  then
                  thickness(i,j,k) = Input_mp%pmass(i,j,k)/  &
                                                Atmos_state%airdens(i,j,k)
                  wp2t(i,j,k) = 2.0/(3.0*0.548**2)* &
                      (0.5*(Input_mp%diff_t(i,j,k) +  &
                                 Input_mp%diff_t(i,j,min(k+1,KDIM)))/&
                                                     thickness(i,j,k) )**2
                  wp2(i,j,k) = MAX (wp2t(i,j,k), var_limit**2)
                  if(diag_id%subgrid_w_variance > 0)   &
                        diag_4d(i,j,k,diag_pt%subgrid_w_variance) =   &
                                                           wp2(i,j,k)**0.5
               
                  call aer_ccn_act_wpdf_m   &
                     (Input_mp%tin(i,j,k), Input_mp%pfull(i,j,k), &
                      up_strat(i,j,k), wp2(i,j,k), wpdf_offs,   &
                      Particles%totalmass1(i,j,k,:),    &
                                                   Particles%drop1(i,j,k))

                else
                  Particles%drop1(i,j,k) = 0.                
                endif
              end do
            end do
          end do

        elseif (var_limit_opt == 2) then   ! cjg
          do k = 1,kdim
            do j=1,jdim
              do i = 1,idim

                wp2(i,j,k) = var_limit**2
                if (diag_id%subgrid_w_variance > 0)   &
                  diag_4d(i,j,k,diag_pt%subgrid_w_variance) =   &
                                                         wp2(i,j,k)**0.5
                call aer_ccn_act_wpdf_m   &
                    (Input_mp%tin(i,j,k), Input_mp%pfull(i,j,k), &
                     up_strat(i,j,k), wp2(i,j,k), wpdf_offs,   &
                     Particles%totalmass1(i,j,k,:), Particles%drop1(i,j,k))

              end do
            end do
          end do
        endif   ! cjg
        call mpp_clock_end   (aero_loop2)

!-------------------------------------------------------------------------
!    if ice nuclei activation is being done via a vertical velocity pdf,
!    execute the following code.
!-------------------------------------------------------------------------
        call mpp_clock_begin (aero_loop3)
        if (do_ice_num) then
        if (do_ice_nucl_wpdf) THEN
          do k = 1,kdim
            do j = 1,jdim
              do i = 1,idim

!------------------------------------------------------------------------
!    ice nuclei activation may only occur at temps below -5 C.
!------------------------------------------------------------------------
                if (Input_mp%tin(i,j,k) .LT. tfreeze - 5.) then

!-----------------------------------------------------------------------
!    place a lower limit on the velocity pdf variance, if desired.
!-----------------------------------------------------------------------
                  IF (do_var_limit_ice) THEN
                    wp2i(i,j,k) = MAX (wp2t(i,j,k), var_limit_ice**2)
                  ELSE
                    wp2i(i,j,k) = wp2(i,j,k)
                  END IF

!------------------------------------------------------------------------
!    define the saturation specific humidities over liquid and ice.
!------------------------------------------------------------------------
                  eslt(i,j,k) = polysvp_l(Input_mp%tin(i,j,k))
                  qs_d(i,j,k) = Input_mp%pfull(i,j,k) - d378*eslt(i,j,k)
                  qs_d(i,j,k) = max(qs_d(i,j,k),eslt(i,j,k))
                  qvsl(i,j,k) = 0.622 *eslt(i,j,k)/qs_d(i,j,k)

                  esit(i,j,k) = polysvp_i(Input_mp%tin(i,j,k))
                  qs_d(i,j,k) = Input_mp%pfull(i,j,k) - d378*esit(i,j,k)
                  qs_d(i,j,k) = max(qs_d(i,j,k),esit(i,j,k))
                  qvsi(i,j,k) = 0.622 *esit(i,j,k)/qs_d(i,j,k)

!------------------------------------------------------------------------
!    define the relative humidities wrt ice and liquid to be used in the 
!    calculation of ice nuclei activation. it may vary based on the nml 
!    variable rh_act_opt.
!------------------------------------------------------------------------
                  IF (rh_act_opt .EQ. 1) THEN
                    qvt(i,j,k) = Input_mp%qin(i,j,k)
                    cf(i,j,k) = 0.
                  ELSE
                    cf(i,j,k) = qa_upd(i,j,k) +     &
                                    C2ls_mp%convective_humidity_area(i,j,k)
                    IF (cf(i,j,k) .LT. cf_thresh_nucl) THEN
                      qvt(i,j,k) =  (Input_mp%qin(i,j,k) -   &
                                      cf(i,j,k)*Atmos_State%qs(i,j,k))/ &
                                                           (1. - cf(i,j,k))
                    ELSE
                      qvt(i,j,k) =  Input_mp%qin(i,j,k)
                    ENDIF
                  END IF

                  if (qvt(i,j,k) .LE. 0.0) then
                     qvt(i,j,k) =  MAX(Input_mp%qin(i,j,k), qmin)
                  endif
                  u_i(i,j,k) =  qvt(i,j,k)/qvsi(i,j,k)
                  u_l(i,j,k) =  qvt(i,j,k)/qvsl(i,j,k)

!--------------------------------------------------------------------------
!    if debugging is active and the relative humidity exceeds 200%, output
!    relevant variables.
!--------------------------------------------------------------------------
                  if (debug .and. Input_mp%tin(i,j,k) .lt. 260. .and.  &
                                            u_i(i,j,k) .gt. 200. ) then
                     call aerosol_cloud_debug1 (    &
                           i, j, k, u_i(i,j,k), Atmos_state%qs(i,j,k),   &
                           qvsi(i,j,k), qvt(i,j,k), Input_mp%qin(i,j,k), &
                           cf(i,j,k),                                   &
                           C2ls_mp%convective_humidity_area(i,j,k),   &
                           qa_upd(i,j,k),   &
                           C2ls_mp%convective_humidity_ratio(i,j,k) )
                  endif

!-------------------------------------------------------------------------
!    call ice_nucl_wpdf to obtain number of activated ice crystals.
!-------------------------------------------------------------------------
                  call ice_nucl_wpdf (    &
                          Input_mp%tin(i,j,k), u_i(i,j,k), u_l(i,j,k),&
                          up_strat(i,j,k), wp2i(i,j,k),  & 
                          Input_mp%zfull(i,j,k),    &
                          Particles%totalmass1(i,j,k,:), &
                          Particles%imass1(i,j,k,:), n_totmass, n_imass, &
                          Particles%crystal1(i,j,k), &
                          Particles%drop1(i,j,k), Particles%hom(i,j,k), &
                          Atmos_state%rh_crit(i,j,k),  &
                          Atmos_state%rh_crit_min(i,j,k), ni_sulf(i,j,k), &
                          ni_dust(i,j,k), ni_bc(i,j,k))
                else
                  ni_sulf(i,j,k) = 0.
                  ni_dust(i,j,k) = 0.
                  ni_bc  (i,j,k) = 0.
                  cf(i,j,k) = missing_value
                  u_i(i,j,k) =  missing_value            
                  u_l(i,j,k) =  missing_value
                endif

!-------------------------------------------------------------------------
!    define the critical relative humidity that was used for ice nuclei 
!    activation, averaged over all spectral intervals (rh_crit). also define
!    the minimum value from any of the spectral regions (rh_crit_min).
!-------------------------------------------------------------------------
                IF (Input_mp%tin(i,j,k) .LT. 250.) THEN 
                  Atmos_state%rh_crit(i,j,k) =     &
                                      MAX(Atmos_state%rh_crit(i,j,k), 1.) 
                  Atmos_state%rh_crit_min(i,j,k) =  &
                                  MAX(Atmos_state%rh_crit_min(i,j,k), 1.) 

!-------------------------------------------------------------------------
!    if debugging is active, and the minimum is over 250%, output relevant
!    information.
!-------------------------------------------------------------------------
                  IF ( debug ) THEN
                    IF (Atmos_state%rh_crit_min(i,j,k) .GT. 2.5  ) THEN
                      call aerosol_cloud_debug2 (i,j,k,   &
                                           Atmos_state%rh_crit_min(i,j,k))
                    END IF
                  END IF
                ELSE

!-------------------------------------------------------------------------
!    if temp over 250 K, set values to 1.0.
!-------------------------------------------------------------------------
                  Atmos_state%rh_crit(i,j,k) = 1.
                  Atmos_state%rh_crit_min(i,j,k) = 1.
                END IF
              end do
            end do
          end do

!-------------------------------------------------------------------------
!    define various desired diagnostics.
!-------------------------------------------------------------------------

          if ( diag_id%imass7 > 0 )    &   
                  diag_4d(:,:,:,diag_pt%imass7) = Particles%imass1(:,:,:,7)

          if(diag_id%potential_crystals > 0)   &
                  diag_4d(:,:,:,diag_pt%potential_crystals) =  &
                                                       Particles%crystal1

          if ( diag_id%rhcrit > 0 )     &  
                  diag_4d(:,:,:,diag_pt%rhcrit) = 100.*Atmos_state%rh_crit

          if ( diag_id%rhcrit_min > 0 )     &  
                  diag_4d(:,:,:,diag_pt%rhcrit_min) =   &
                                            100.*Atmos_state%rh_crit_min

          if ( diag_id%ndust1 > 0 )     &  
                  diag_4d(:,:,:,diag_pt%ndust1) = Nfact_du1 *   &
                                                  Particles%imass1(:,:,:,8)

          if ( diag_id%ndust2 > 0 )     &  
                  diag_4d(:,:,:,diag_pt%ndust2) = Nfact_du2 *    &
                                                  Particles%imass1(:,:,:,9)

          if ( diag_id%ndust3 > 0 )     &  
                  diag_4d(:,:,:,diag_pt%ndust3) = Nfact_du3 *   &
                                                Particles%imass1(:,:,:,10)
         
          if ( diag_id%ndust4 > 0 )     &  
                  diag_4d(:,:,:,diag_pt%ndust4) = Nfact_du4 *   &
                                                Particles%imass1(:,:,:,11)
    
          if ( diag_id%ndust5 > 0 )    &
                  diag_4d(:,:,:,diag_pt%ndust5) = Nfact_du5 *  &
                                                Particles%imass1(:,:,:,12)

          if ( diag_id%ni_dust > 0 )     &
                  diag_4d(:,:,:,diag_pt%ni_dust) = ni_dust

          if ( diag_id%ni_bc > 0 ) diag_4d(:,:,:,diag_pt%ni_bc) = ni_bc

          if ( diag_id%ni_sulf > 0 )      & 
                  diag_4d(:,:,:,diag_pt%ni_sulf) = ni_sulf

          if (rh_act_opt .ne. 1) then
            if ( diag_id%cfin > 0 ) diag_4d(:,:,:,diag_pt%cfin) = cf
          end if

          if ( diag_id%rhiin > 0 ) diag_4d(:,:,:,diag_pt%rhiin) = u_i

          if ( diag_id%rhlin > 0 ) diag_4d(:,:,:,diag_pt%rhlin) = u_l

        END IF ! do_ice_nucl_wpdf
        END IF ! do_ice_num

        if(diag_id%potential_droplets > 0)   &
                      diag_4d(:,:,:,diag_pt%potential_droplets) =   &
                                                        Particles%drop1

        call mpp_clock_end (aero_loop3)
      end if  ! (do_liq_num)

!-------------------------------------------------------------------------


END SUBROUTINE  determine_activated_aerosol


!------------------------------------------------------------------------

subroutine  aerosol_cloud_end

      module_is_initialized = .false.

      call polysvp_end
      call aer_ccn_act_end
      call ice_nucl_wpdf_end

end subroutine  aerosol_cloud_end

!#########################################################################

 subroutine aerosol_effects (idim,jdim,kdim,n_diag_4d, pthickness,  &
                            concen_dust_sub, totalmass1, imass1, Aerosol, &
                            diag_4d, diag_id, diag_pt )

INTEGER,                                   INTENT (in)   :: idim, jdim,  &
                                                            kdim, n_diag_4d
REAL, DIMENSION(idim,jdim,kdim),           INTENT(in )   :: pthickness
REAL, DIMENSION(idim,jdim,kdim),           INTENT(inout) :: concen_dust_sub
REAL, DIMENSION(idim,jdim,kdim,n_totmass), INTENT(inout) :: totalmass1
REAL, DIMENSION(idim,jdim,kdim,n_imass),   INTENT(inout) :: imass1
TYPE(aerosol_type),                        INTENT (in)   :: Aerosol  
REAL, DIMENSION(idim,jdim,kdim,0:n_diag_4d),    &
                                           INTENT(INOUT) :: diag_4d
TYPE(diag_id_type),                        intent(in)    :: diag_id
TYPE(diag_pt_type),                        intent(in)    :: diag_pt

!-----------------------------------------------------------------------
!   local variables:

      REAL, DIMENSION(idim,jdim,kdim) :: concen_all_sub, concen_ss_sub, &
                                         concen_ss_sup
      INTEGER :: i, j, k, na, s
      LOGICAL :: used



!-------------------------------------------------------------------------
!    initialize local accumulation arrays. 
!-------------------------------------------------------------------------
      if (use_online_aerosol) then
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              concen_ss_sub(i,j,k) = 0.
              concen_ss_sup(i,j,k) = 0.
              concen_all_sub(i,j,k) = 0.
            end do
          end do
        end do
      endif

!-------------------------------------------------------------------------
!    define the aerosol content in various categories that may be used for
!    droplet or ice crystal nucleation. save various diagnostics as they
!    are calculated.
!-------------------------------------------------------------------------
      if (do_liq_num) then
        if (use_online_aerosol) then
          do na = 1,size(Aerosol%aerosol,4)               
            if (trim(Aerosol%aerosol_names(na)) == 'so4' .or. &
                trim(Aerosol%aerosol_names(na)) == 'so4_anthro' .or.&
                trim(Aerosol%aerosol_names(na)) == 'so4_natural')  then
              do k=1,kdim
                do j=1,jdim
                  do i=1,idim
                    totalmass1(i,j,k,1) = totalmass1(i,j,k,1) + &
                                          Aerosol%aerosol(i,j,k,na)
                  end do
                end do
              end do

            else if(trim(Aerosol%aerosol_names(na)) == 'omphilic' .or.&
                   trim(Aerosol%aerosol_names(na)) == 'omphobic')  then
              do k=1,kdim
                do j=1,jdim
                  do i=1,idim
                    totalmass1(i,j,k,4) = totalmass1(i,j,k,4) +  &
                                          Aerosol%aerosol(i,j,k,na)

                  end do
                end do
              end do

!-----------------------------------------------------------------------
! h1g, 2015-09-18
!    for fast aerosol, seasalt and dust names are changed from seasalt1,
!    seasalt2, ... dust1, dust2, ... to seasalt_aitken, ..., 
!    dust_mode1_of_2
!-----------------------------------------------------------------------
            else if(trim(Aerosol%aerosol_names(na)) == 'seasalt1' &
                 .or.trim(Aerosol%aerosol_names(na)) ==   &
                                             'seasalt_aitken' ) then 
              do k=1,kdim
                do j=1,jdim
                  do i=1,idim
                    concen_ss_sub(i,j,k) = concen_ss_sub(i,j,k) +  &
                                           Aerosol%aerosol(i,j,k,na)
                    imass1(i,j,k,1) = Aerosol%aerosol(i,j,k,na)/  &
                                                         pthickness(i,j,k)
                  end do
                end do
              end do

            else if(trim(Aerosol%aerosol_names(na)) == 'seasalt2' &
                 .or.trim(Aerosol%aerosol_names(na)) ==    &
                                                   'seasalt_fine' ) then
              do k=1,kdim
                do j=1,jdim
                  do i=1,idim
                    concen_ss_sub(i,j,k) = concen_ss_sub(i,j,k) +  &
                                           Aerosol%aerosol(i,j,k,na)
                    imass1(i,j,k,2) = Aerosol%aerosol(i,j,k,na)/   &
                                                         pthickness(i,j,k)
                  end do
                end do
              end do

            else if(trim(Aerosol%aerosol_names(na)) == 'seasalt3' &
                .or.trim(Aerosol%aerosol_names(na)) ==    &
                                                 'seasalt_coarse') then
              do k=1,kdim
                do j=1,jdim
                  do i=1,idim
                    concen_ss_sup(i,j,k) = concen_ss_sup(i,j,k) +  &
                                           Aerosol%aerosol(i,j,k,na)
                    imass1(i,j,k,3) = Aerosol%aerosol(i,j,k,na)/  &
                                                         pthickness(i,j,k)
                  end do
                end do
              end do

            else if(trim(Aerosol%aerosol_names(na)) == 'seasalt4' ) then 
              do k=1,kdim
                do j=1,jdim
                  do i=1,idim
                    concen_ss_sup(i,j,k) = concen_ss_sup(i,j,k) +  &
                                           Aerosol%aerosol(i,j,k,na)
                    imass1(i,j,k,4) = Aerosol%aerosol(i,j,k,na)/  &
                                                         pthickness(i,j,k)
                  end do
                end do
              end do

            else if(trim(Aerosol%aerosol_names(na)) == 'seasalt5' ) then
              do k=1,kdim
                do j=1,jdim
                  do i=1,idim
                    concen_ss_sup(i,j,k) = concen_ss_sup(i,j,k) +  &
                                           Aerosol%aerosol(i,j,k,na)
                    imass1(i,j,k,5) = Aerosol%aerosol(i,j,k,na)/   &
                                                         pthickness(i,j,k)
                  end do
                end do
              end do
             
            else if(trim(Aerosol%aerosol_names(na)) == 'bcphilic') then 
              do k=1,kdim
                do j=1,jdim
                  do i=1,idim
                    concen_all_sub(i,j,k) = concen_all_sub(i,j,k) +  &
                                            Aerosol%aerosol(i,j,k,na)
                    imass1(i,j,k,6) = Aerosol%aerosol(i,j,k,na)/ &
                                                         pthickness(i,j,k)
                  end do
                end do
              end do

            else if(trim(Aerosol%aerosol_names(na)) == 'bcphobic' ) then 
              do k=1,kdim
                do j=1,jdim
                  do i=1,idim
                    concen_all_sub(i,j,k) = concen_all_sub(i,j,k) +  &
                                            Aerosol%aerosol(i,j,k,na)
                  end do
                end do
              end do

            else if(trim(Aerosol%aerosol_names(na)) == 'dust1' &
                 .or. trim(Aerosol%aerosol_names(na)) ==     &
                                                  'dust_mode1_of_2' ) then
              do k=1,kdim
                do j=1,jdim
                  do i=1,idim
                    concen_all_sub(i,j,k) = concen_all_sub(i,j,k) +  &
                                            Aerosol%aerosol(i,j,k,na)
                    imass1(i,j,k,8) = Aerosol%aerosol(i,j,k,na)/  &
                                                         pthickness(i,j,k)
! according to Paul dust 1  and dust2 are sub- 2.5 micron
                    imass1(i,j,k,7) =    imass1(i,j,k,7) +   &
                                              Aerosol%aerosol(i,j,k,na)
                    if (do_dust_berg) then
                      concen_dust_sub(i,j,k) =   concen_dust_sub(i,j,k) + &
                                                 Aerosol%aerosol(i,j,k,na)
                    endif
                  end do
                end do
              end do

            else if(trim(Aerosol%aerosol_names(na)) == 'dust2' ) then   
              do k=1,kdim
                do j=1,jdim
                  do i=1,idim
                    concen_all_sub(i,j,k) = concen_all_sub(i,j,k) +  &
                                            Aerosol%aerosol(i,j,k,na)
                    imass1(i,j,k,9) = Aerosol%aerosol(i,j,k,na)/  &
                                                         pthickness(i,j,k)
! according to Paul dust 1  and dust2 are sub- 2.5 micron
                    imass1(i,j,k,7) =    imass1(i,j,k,7) +   &
                                              Aerosol%aerosol(i,j,k,na)
                    if (do_dust_berg) then
                      concen_dust_sub(i,j,k) = concen_dust_sub(i,j,k) + &
                                               Aerosol%aerosol(i,j,k,na)
                    endif
                  end do
                end do
              end do

            else if(trim(Aerosol%aerosol_names(na)) == 'dust3' ) then 
              do k=1,kdim
                do j=1,jdim
                  do i=1,idim
                    concen_all_sub(i,j,k) = concen_all_sub(i,j,k) +  &
                                            Aerosol%aerosol(i,j,k,na)
                    imass1(i,j,k,10) = Aerosol%aerosol(i,j,k,na)/   &
                                                         pthickness(i,j,k)
                    if (do_dust_berg) then
                      concen_dust_sub(i,j,k) = concen_dust_sub(i,j,k) + &
                                               Aerosol%aerosol(i,j,k,na)
                    endif
                  end do
                end do
              end do

            else if (trim(Aerosol%aerosol_names(na)) == 'dust4' &
                .or. trim(Aerosol%aerosol_names(na)) ==    &
                                                 'dust_mode2_of_2') then
              do k=1,kdim
                do j=1,jdim
                  do i=1,idim
                    imass1(i,j,k,11) = Aerosol%aerosol(i,j,k,na)/   &   
                                                         pthickness(i,j,k)
                  end do
                end do
              end do

            else if (trim(Aerosol%aerosol_names(na)) == 'dust5') then
              do k=1,kdim
                do j=1,jdim
                  do i=1,idim
                    imass1(i,j,k,12) = Aerosol%aerosol(i,j,k,na)/ &      
                                                         pthickness(i,j,k)
                  end do
                end do
              end do

            endif
          end do
          
          do k=1,kdim
            do j=1,jdim
              do i=1,idim
                totalmass1(i,j,k,3) = concen_ss_sub(i,j,k)
                totalmass1(i,j,k,2) = concen_all_sub(i,j,k) + &
                                      totalmass1(i,j,k,4) + &
                                      concen_ss_sub(i,j,k)
                if (use_sub_seasalt) then
                else
                  totalmass1(i,j,k,3) = concen_ss_sub(i,j,k) +  &
                                                  concen_ss_sup(i,j,k)
                endif
                if (sea_salt_scale_onl .NE. 1. ) then
                  totalmass1(i,j,k,3) =     sea_salt_scale_onl*  &
                                                       totalmass1(i,j,k,3)
                end if
                if (diag_id%sulfate > 0) then
                  diag_4d(i,j,k, diag_pt%sulfate) =   &
                                             0.7273*totalmass1(i,j,k,1)/  &
                                                    pthickness(i,j,k)*1.0e9
                endif
         
                if (diag_id%seasalt_sub > 0) then
                  diag_4d(i,j,k,diag_pt%seasalt_sub) =    &
                                                  concen_ss_sub(i,j,k)/  &
                                                    pthickness(i,j,k)*1.0e9
                endif 
                if (diag_id%seasalt_sup > 0) then
                  diag_4d(i,j,k, diag_pt%seasalt_sup) =    &
                                                 concen_ss_sup(i,j,k)/  &
                                                    pthickness(i,j,k)*1.0e9
                endif 
              end do
            end do
           end do

         else  ! (use_online_aerosol)

           if (do_dust_berg) then

!     YMice submicron dust (NO. 14 to NO. 18)
             do s = 14,18
               do k=1,kdim
                 do j=1,jdim
                   do i=1,idim
                     concen_dust_sub(i,j,k) = concen_dust_sub(i,j,k)+ &
                                              Aerosol%aerosol(i,j,k,s)
                   end do
                 end do
               end do
             end do
           endif

           if (diag_id%sulfate > 0) then
!     anthro. and natural sulfate concentration (ug so4/m3)
             do k=1,kdim
               do j=1,jdim
                 do i=1,idim
                   diag_4d(i,j,k, diag_pt%sulfate) = 0.7273*   &
                            (Aerosol%aerosol(i,j,k,1) +      &
                                            Aerosol%aerosol(i,j,k,2))/&
                                                pthickness(i,j,k)*1.0e9
                 end do
               end do
             end do
           endif

!offline
! NO. 1 natural Sulfate; NO. 2 anthro. sulfate; NO. 3 Sea Salt; NO. 4 Or        ganics
           do k=1,kdim
             do j=1,jdim
               do i=1,idim
                 totalmass1(i,j,k,1) = Aerosol%aerosol(i,j,k,2)
                 totalmass1(i,j,k,2) = Aerosol%aerosol(i,j,k,1)
                 totalmass1(i,j,k,3) = sea_salt_scale*  &
                                       Aerosol%aerosol(i,j,k,5)
                 totalmass1(i,j,k,4) = om_to_oc*  &
                                       Aerosol%aerosol(i,j,k,3)
               end do
             end do
           end do
         endif ! (use_online_aerosol)

        do na = 1, 4
          do k=1,kdim
            do j=1,jdim
              do i=1,idim
                totalmass1(i,j,k,na) = totalmass1(i,j,k,na)/  &
                                       pthickness(i,j,k)*1.0e9*1.0e-12
              end do
            end do
          end do
        end do

        do k=1,kdim
          do j=1,jdim
            do i=1,idim
! submicron dust concentration (ug/m3) (NO. 2 to NO. 4)
              imass1(i,j,k,7) = imass1(i,j,k,7)/pthickness(i,j,k)*1.0e9 
            end do
          end do
        end do


        if (do_dust_berg) then
! submicron dust concentration (ug/m3) (NO. 2 to NO. 4)
          do k=1,kdim
            do j=1,jdim
              do i=1,idim
                concen_dust_sub(i,j,k) = concen_dust_sub(i,j,k)/ &
                                                   pthickness(i,j,k)*1.0e9 
              end do
            end do
          end do
        endif

        if (diag_id%om > 0) then
          do k=1,kdim
            do j=1,jdim
              do i=1,idim
!RSH the 1.0e12 factor here is to counter the  1.0e-12 factor applied 
! above to totalmass1. The 1.0e9 factor (kg -> ug) is already in totalmass1.
                diag_4d(i,j,k,diag_pt%om) = totalmass1(i,j,k,2)*1.0e12
              end do
            end do
          end do
        endif

      endif  ! (do_liq_num)

!----------------------------------------------------------------------


end subroutine aerosol_effects

!#########################################################################



                    END MODULE aerosol_cloud_mod
