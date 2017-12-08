                      module tiedtke_macro_mod

use mpp_mod,                    only : input_nml_file
use fms_mod,                    only : FATAL, error_mesg, mpp_pe, &
                                       mpp_root_pe, open_namelist_file, &
                                       check_nml_error, close_file,  &
                                       write_version_number, file_exist, &
                                       stdlog
use constants_mod,              only : GRAV, TFREEZE, CP_AIR, HLV, HLS, &
                                       RVGAS, pi
use beta_dist_mod,              only : incomplete_beta, beta_dist_init, &
                                       beta_dist_end
use lscloud_types_mod,          only : lscloud_types_init, &
                                       diag_id_type, diag_pt_type, &
                                       atmos_state_type,  &
                                       lscloud_nml_type, &
                                       cloud_state_type, particles_type,&
                                       lsc_constants_type, &
                                       cloud_processes_type
use lscloud_debug_mod,          only : lscloud_debug_init,  &
                                       macrophysics_debug
use polysvp_mod,                only : polysvp_init, polysvp_end,  &
                                       compute_qs_x1
use moist_proc_utils_mod,       only : mp_input_type, mp_conv2ls_type, &
                                       mp_nml_type, mp_lsdiag_type,   &
                                       mp_lsdiag_control_type
use physics_radiation_exch_mod, only : exchange_control_type

implicit none
private 

!-----------------------------------------------------------------------
!---interfaces----------------------------------------------------------
public  tiedtke_macro_init, tiedtke_macro, tiedtke_macro_diagnostics,  &
        tiedtke_macro_end, tiedtke_macro_time_vary

private tiedtke_macro_nopdf_nosuper, tiedtke_macro_nopdf_super,   &
        tiedtke_macro_Single_Gaussian_pdf, erf, &
        tiedtke_macro_pdf, ppm2m_sak, kmppm_sak, steepz_sak

!-----------------------------------------------------------------------
!---version number------------------------------------------------------

Character(len=128) :: Version = '$Id$'
Character(len=128) :: Tagname = '$Name$'

!-----------------------------------------------------------------------
!---namelist------------------------------------------------------------

!  <DATA NAME="U00" UNITS="fraction" TYPE="real" DEFAULT="0.80">
!   threshold rel hum for cloud formation by large-scale condensation.
!  </DATA
!  <DATA NAME="u00_profile" TYPE="logical"  DEFAULT=".false.">
!   the low-level u00 ECMWF profile should be applied?
!  </DATA>
!  <DATA NAME="eros_scale" UNITS="1/sec" TYPE="real" DEFAULT="1.E-06">
!   normal erosion rate constant for cloud destruction (default = 1.E-06)
!  </DATA>
!  <DATA NAME="eros_choice" TYPE="logical"  DEFAULT=".false.">
!   enhanced erosion should occur in turbulent or convective conditions ?
!  </DATA>
!  <DATA NAME="eros_scale_c" UNITS="1/sec" TYPE="real"  DEFAULT="8.E-06">
!   erosion rate constant for cloud destruction for convective conditions.
!  </DATA>
!  <DATA NAME="eros_scale_t" UNITS="1/sec" TYPE="real" DEFAULT="5.E-05">
!   erosion rate constant for cloud destruction for turbulent conditions.
!  </DATA>
!  <DATA NAME="mc_thresh" UNITS="kg/m2/sec" TYPE="real" DEFAULT="0.001">
!   convective mass-flux threshold for enhanced erosion to turn on.
!  </DATA
!  <DATA NAME="diff_thresh" UNITS="m2/s" TYPE="real" DEFAULT="1.0">
!   diffusion coefficient threshold for enhanced erosion to turn on.
!  </DATA>
!  <DATA NAME="efact" UNITS="dimensionless" TYPE="real" DEFAULT="0.0">
!   factor used in formula to enhance cloud erosion as a function of height
!  </DATA>

real    :: U00 =0.80
logical :: u00_profile = .true.
real    :: eros_scale     =  1.E-06
logical :: eros_choice    =  .false.
real    :: eros_scale_c   =  8.E-06
real    :: eros_scale_t   =  5.E-05
real    :: mc_thresh      =  0.001   
real    :: diff_thresh   =  0.1      
real    :: efact = 0.

logical :: add_ahuco = .true.              ! include convective cloud area
                                           ! in grid box when determining 
                                           ! grid box mean relative 
                                           ! humidity required to allow
                                           ! large-scale condensation
logical :: use_qabar = .true.              ! use time-mean cloud area 
                                           ! during the timestep in
                                           ! calculating cloud processes;
                                           ! otherwise use cloud area at
                                           ! end of step
logical :: rk_repartition_first = .false.  ! repartition the condensation 
                                           ! between liquid and ice
                                           ! based on bergeron relation
                                           ! before computing microphysics?
logical :: do_aero_eros = .false.          ! enhance cloud erosion based
                                           ! on number of droplets
real    :: ae_lb = 0.8                     ! used in aero_eros calculation
real    :: ae_ub = 1.2                     ! used in aero_eros calculation
real    :: ae_N_lb = 0.                    ! used in aero_eros calculation
real    :: ae_N_ub = 150.                  ! used in aero_eros calculation
logical :: include_neg_mc = .false.
logical :: Single_Gaussion_pdf = .false.

namelist / tiedtke_macro_nml /  &
                       U00, U00_profile, eros_scale, eros_choice,      &
                       eros_scale_c, eros_scale_t, mc_thresh, diff_thresh,&
                       efact, add_ahuco, use_qabar, rk_repartition_first, &
                       do_aero_eros, ae_lb, ae_ub, ae_N_lb, ae_N_ub, &
                       include_neg_mc, Single_Gaussion_pdf
 
!------------------------------------------------------------------------
!    module variables obtained from other modules
!------------------------------------------------------------------------
real    :: qmin
real    :: qcvar
logical :: limit_conv_cloud_frac
real    :: Dmin
real    :: cfact
integer :: super_ice_opt
logical :: do_pdf_clouds
integer :: betaP
real    :: qthalfwidth
integer :: nsublevels
integer :: kmap
integer :: kord
logical :: pdf_org
real    :: dtcloud, inv_dtcloud
logical :: do_rk_microphys
logical :: total_activation
logical :: debug

!------------------------------------------------------------------------

logical :: module_is_initialized = .false.





                               CONTAINS



!#######################################################################

SUBROUTINE tiedtke_macro_init (Nml_mp, Nml_lsc, Constants_lsc, Exch_ctrl)

!-----------------------------------------------------------------------
type(lscloud_nml_type),      intent(in) :: Nml_lsc
type(lsc_constants_type),    intent(in) :: Constants_lsc
type(mp_nml_type),           intent(in) :: Nml_mp
type(exchange_control_type), intent(in) :: Exch_ctrl
 
!-----------------------------------------------------------------------
!---local variables----

      integer :: unit, io, ierr, logunit

!-----------------------------------------------------------------------
      if (module_is_initialized) return

!-------------------------------------------------------------------------
!    process namelist.
!-------------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=tiedtke_macro_nml, iostat=io)
      ierr = check_nml_error(io,'tiedtke_macro_nml')
#else
      if ( file_exist('input.nml')) then
        unit = open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=tiedtke_macro_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'tiedtke_macro_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!-------------------------------------------------------------------------
!    write version and namelist to standard log.
!-------------------------------------------------------------------------
      call write_version_number ( version, tagname )
      logunit = stdlog()
      if ( mpp_pe() == mpp_root_pe() ) &
                        write ( logunit, nml=tiedtke_macro_nml )

!----------------------------------------------------------------------
!    save variables that will be needed later in this module.
!----------------------------------------------------------------------
      qmin = Exch_ctrl%qmin
      qcvar = Exch_ctrl%qcvar
      limit_conv_cloud_frac = Nml_mp%limit_conv_cloud_frac
      Dmin = Nml_lsc%Dmin
      cfact = Nml_lsc%cfact
      super_ice_opt = Nml_lsc%super_ice_opt
      do_pdf_clouds = Nml_lsc%do_pdf_clouds
      betaP = Nml_lsc%betaP
      qthalfwidth = Nml_lsc%qthalfwidth
      nsublevels = Nml_lsc%nsublevels
      kmap = Nml_lsc%kmap
      kord = Nml_lsc%kord
      pdf_org = Nml_lsc%pdf_org
      do_rk_microphys = Constants_lsc%do_rk_microphys
      total_activation = Constants_lsc%total_activation

!------------------------------------------------------------------------
!    be sure needed modules are initialized.
!------------------------------------------------------------------------
      call lscloud_types_init
      call lscloud_debug_init (debug)
      call polysvp_init (Nml_lsc)
      if (do_pdf_clouds) then
        call beta_dist_init
      endif

      module_is_initialized = .true.

!-----------------------------------------------------------------------


END SUBROUTINE tiedtke_macro_init


!########################################################################

subroutine tiedtke_macro_time_vary (dtcloud_in)

real, intent(in) :: dtcloud_in

!---------------------------------------------------------------------
!    save variables needed during prognostic loop.
!---------------------------------------------------------------------
      dtcloud = dtcloud_in
      inv_dtcloud = 1./dtcloud

end subroutine tiedtke_macro_time_vary 



!########################################################################

SUBROUTINE tiedtke_macro (     &
             idim, jdim, kdim, C2ls_mp, Input_mp, Atmos_state,    &
             Cloud_state, ST, SQ, Cloud_processes, Particles, Lsdiag_mp, &
                                                       Lsdiag_mp_control)

!------------------------------------------------------------------------
INTEGER,                         INTENT(IN )   :: idim, jdim, kdim 
type(mp_conv2ls_type),           intent(inout) :: C2ls_mp    
type(mp_input_type),             intent(inout) :: Input_mp    
type(atmos_state_type),          intent(inout) :: Atmos_state
type(cloud_state_type),          intent(inout) :: Cloud_state
REAL, dimension(idim,jdim,kdim), INTENT(IN )   :: ST, SQ
type(cloud_processes_type),      intent(inout) :: Cloud_processes
type(particles_type),            intent(inout) :: Particles
type(mp_lsdiag_type),            intent(inout) :: Lsdiag_mp    
type(mp_lsdiag_control_type),    intent(inout) :: Lsdiag_mp_control    

!------------------------------------------------------------------------
!----local variables----
!
!       SA             cloud area tendency
!       U00p           critical relative humidity fraction which may be 
!                      a function of pressure 
!       U00pr          gridbox mean rh needed for condensation, after
!                      accounting for convective cloud area
!       liq_frac       fraction of condensation which is liquid
!       ice_frac       fraction of condensation which is frozen
!       erosion_scale  cloud erosion scale
!       edum           factor used in aero_eros calculation
!       dum            ratio of bergeron to condensation tendencies
!       D2_dts         tendency due to bergeron process
!       A_plus_B_s     sum of vapor diffusion factor and thermal conduct-
!                      ivity factor
!       mdum           factor used in aero_eros calculation
!       bdum           factor used in aero_eros calculation
!       i,j,k          do-loop indices
!-------------------------------------------------------------------------

      real, dimension(idim, jdim,kdim) :: SA, U00p, U00pr, liq_frac, &
                                          ice_frac, erosion_scale, edum
      real                             :: dum, D2_dts, A_plus_B_s, mdum, &
                                          bdum
      INTEGER                          :: i, j ,k


!-----------------------------------------------------------------------
!    define a local copy of the current cloud area tendency. this will
!    be modified in this routine and the modified version returned to the
!    Cloud_state%SA_out variable upon completion.
!-----------------------------------------------------------------------
      SA = Cloud_state%SA_out

!------------------------------------------------------------------------
!    process the non-convective condensation for pdf clouds. 
!                                                                      
!                NON-CONVECTIVE CONDENSATION                          
!                                                                      
!                STATISTICAL CLOUD FRACTION                            
!                                                                      
!------------------------------------------------------------------------
      if (do_pdf_clouds) then
        if (Single_Gaussion_pdf) then
!RSH remove constants, nml, otun as args.
!         call nc_cond_Single_Gaussian_pdf(idim, jdim, kdim, Nml, 
!         call nc_cond_Single_Gaussian_pdf(idim, jdim, kdim,      
          call tiedtke_macro_Single_Gaussian_pdf(idim, jdim, kdim,    & 
!                         Constants, Atmos_state, &
                       Atmos_state, Input_mp, &
                       Cloud_state, ST, SQ, Cloud_processes, &
                       Particles, Lsdiag_mp_control%n_diag_4d,   &
                       Lsdiag_mp%diag_4d, &
!                      Lsdiag_mp%diag_id, Lsdiag_mp%diag_pt, otun, SA)
                       Lsdiag_mp_control%diag_id,   &
                       Lsdiag_mp_control%diag_pt,       SA)
        else
          call tiedtke_macro_pdf (     &
              idim, jdim, kdim, Atmos_state, Input_mp, Cloud_state, ST,  &
              SQ, Cloud_processes, Particles, Lsdiag_mp_control%n_diag_4d,&
              Lsdiag_mp%diag_4d, Lsdiag_mp_control%diag_id,   &
              Lsdiag_mp_control%diag_pt, SA)
        endif
      else  !(do_pdf_clouds)

!------------------------------------------------------------------------
!    process the non-convective condensation for non-pdf clouds. some
!    additional calculations are needed.
!
!                 NON-CONVECTIVE CONDENSATION                          
!                                                                     
!                 TIEDTKE (1993) CLOUD FRACTION                    
!                                                                    
!       ANALYTIC INTEGRATION OF SATURATED VOLUME FRACTION EQUATION     
!
!       Do non-convective condensation following Tiedtke, pages 3044-5.
!       In this formulation stratiform clouds are only formed/destroyed 
!       when there is upward or downward motion to support/destroy it. 
!-------------------------------------------------------------------------

!-----------------------------------------------------------------------
!    calculate aerosol erosion term which optionally may be used.
!-----------------------------------------------------------------------
        if (do_aero_eros ) then 
          mdum = (ae_ub - ae_lb)/(ae_N_ub - ae_N_lb)
          bdum = ae_lb - mdum*ae_N_lb
          edum = mdum*Particles%drop1 + bdum 
        else
          edum = 1.
        endif

!------------------------------------------------------------------------
!    compute pressure dependent critical relative humidity for onset of
!    condensation (U00p) following ECMWF formula if desired. otherwise,
!    it is given by the nml variable u00.
!------------------------------------------------------------------------
        U00p = U00
        if (u00_profile) then
          DO k=1,kdim
            where (Input_mp%pfull(:,:,k).gt.   &
                          0.8*Input_mp%phalf(:,:,KDIM+1)) 
              U00p(:,:,k) = U00 + (1.- U00)* &
                         (((Input_mp%pfull(:,:,k) -   &
                            (0.8*Input_mp%phalf(:,:,KDIM+1)))/  &
                               (0.2*Input_mp%phalf(:,:,KDIM+1)) )**2.)
            end where
          END DO
        endif       

!------------------------------------------------------------------------
!    modify u00p to account for humidity in convective system.  see 
!    "Tiedtke u00 adjustment" notes, 10/22/02 -- ljd
!    u00p = critical rh for condensation in the grid box, accounting for 
!    the fact that ahuco of the box is saturated (has convective cloud) 
!    better written : ahuco*100% + (1-ahuco)*u00p = gridboxmean rh
!------------------------------------------------------------------------
        IF (add_ahuco) THEN
          u00pr = u00p + (1. - u00p)*C2ls_mp%convective_humidity_area
          u00p = u00pr
        END IF

!-----------------------------------------------------------------------
!    Theory for eros_scale
!
!    If eros_choice equals false, then a single erosion time scale
!    is used in all conditions (eros_scale).  If eros_choice equals
!    true then it is assumed that the timescale for turbulent 
!    evaporation is a function of the conditions in the grid box.  
!    Specifically, if the flow is highly turbulent then the scale is 
!    short, and eros_scale is large.  Likewise if convection is 
!    occurring, then it is assumed that the erosion term is larger 
!    than backround conditions. 
!
!    Here are the typical values for the timescales and the 
!    switches used (subject to changes via namelist):
!
!         Mixing type      eros_scale (sec-1)          Indicator
!       ----------------   ------------------     --------------------
!
!       Background            1.e-06              always present
!       Convective layers     5.e-06              Mc > Mc_thresh
!       Turbulent  layers     5.e-05              diff_t > diff_thresh
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!    define array with background erosion scale.
!-----------------------------------------------------------------------
        erosion_scale =     eros_scale
       
!-----------------------------------------------------------------------
!    Do enhanced erosion in convective or turbulent layers?
!
!                   IMPORTANT NOTE
!                
!    note that convection is considered first, so that if turbulence and 
!    convection occur in the same layer, the erosion rate for turbulence 
!    is selected.                
!-----------------------------------------------------------------------
        if (eros_choice) then
          do k=1,kdim
            do j=1,jdim
              do i=1,idim

!-----------------------------------------------------------------------
!    Enhanced erosion in convective layers
!-----------------------------------------------------------------------
                if (include_neg_mc) then
                  if (abs(C2ls_mp%mc_full(i,j,k)) .gt. mc_thresh) then 
                    erosion_scale(i,j,k) = eros_scale_c
                  endif
                else
                  if (C2ls_mp%mc_full(i,j,k) .gt. mc_thresh) then
                    erosion_scale(i,j,k) = eros_scale_c
                  endif
                endif
                if ((Input_mp%diff_t(i,j,K) .gt. diff_thresh) .or.&
                    (Input_mp%diff_t(i,j,min(k+1,KDIM)) .gt.  &
                                                       diff_thresh) ) then
                  erosion_scale(i,j,K) = eros_scale_t
                endif 
              END DO
            END DO
          END DO
        end if   !(eros_choice)        

!------------------------------------------------------------------------
!    enhance erosion scale if that option has been selected.
!------------------------------------------------------------------------
        IF (do_aero_eros) erosion_scale = edum*erosion_scale

!------------------------------------------------------------------------
!     determine condensation for different supersaturation assumptions.
!------------------------------------------------------------------------
        IF (super_ice_opt .LT. 1 ) THEN

!------------------------------------------------------------------------
!    if supersaturation is not to be allowed in the clear-sky portion
!    of the grid box, call subroutine tiedtke_macro_nopdf_nosuper to 
!    calculate the condensation.
!------------------------------------------------------------------------
          call tiedtke_macro_nopdf_nosuper (    &
               idim, jdim, kdim, C2ls_mp, Input_mp, Atmos_state,   &
               Cloud_state, ST, SQ, Cloud_processes, Particles,   &
               Lsdiag_mp_control%n_diag_4d, Lsdiag_mp%diag_4d,   &
               Lsdiag_mp_control%diag_id, Lsdiag_mp_control%diag_pt,  &
                                       edum, SA, U00p, erosion_scale)
        else

!------------------------------------------------------------------------
!    if supersaturation is to be allowed in the clear-sky portion of the
!    gridbox, call subroutine tiedtke_macro_nopdf_super to calculate the 
!    condensation.
!------------------------------------------------------------------------
          call tiedtke_macro_nopdf_super  (     &
               idim, jdim, kdim, C2ls_mp, Input_mp, Atmos_state,   &
               Cloud_state, ST, SQ, Cloud_processes, Particles,   &
               Lsdiag_mp_control%n_diag_4d, Lsdiag_mp%diag_4d,   &
               Lsdiag_mp_control%diag_id, Lsdiag_mp_control%diag_pt,  &
                                            edum, SA, U00p, erosion_scale)
        endif
      endif  ! (do_pdf_clouds)

!-----------------------------------------------------------------------
!    the next step is the apportionment of the non-convective condensation !    between liquid and ice phases. Following the suggestion of Rostayn
!    (2000), all condensation at temperatures greater than -40C is in 
!    liquid form as ice nuclei are generally limited in the atmosphere.
!    The droplets may subsequently be converted to ice by the 
!    Bergeron-Findeison mechanism.  
!
!    one problem with this formulation is that the proper saturation vapor 
!    pressure is not used for cold clouds as it should be liquid saturation
!    in the case of first forming liquid, but change to ice saturation as 
!    the cloud glaciates.  The current use of ice saturation beneath -20C 
!    thus crudely mimics the result that nearly all stratiform clouds are
!    glaciated for temperatures less than -15C.
!
!    in the case of large-scale evaporation (dcond_ls<0.), it is assumed 
!    that cloud liquid will evaporate faster than cloud ice because if 
!    both are present in the same volume the saturation vapor pressure over
!    the droplet is higher than that over the ice crystal.
!
!    the fraction of large-scale condensation that is liquid is stored in 
!    the temporary variable liq_frac.
!-----------------------------------------------------------------------
      liq_frac = 1.
      ice_frac = 0.

!-----------------------------------------------------------------------
!    for cases of cloud condensation where temperatures are less than -40C 
!    create only ice,
!-----------------------------------------------------------------------
      where (Cloud_processes%dcond_ls .ge. 0. .and.   &
                                Input_mp%tin .lt. tfreeze - 40.)
        liq_frac = 0.
        ice_frac = 1.
      endwhere

!-----------------------------------------------------------------------
!    for cases of cloud evaporation of mixed phase clouds, set liquid 
!    evaporation to preferentially occur first.
!-----------------------------------------------------------------------
      where ((Cloud_processes%dcond_ls .lt. 0.) .and.   &
                     (Cloud_state%ql_upd .gt. qmin)  .and.        &
                            (Cloud_state%qi_upd .gt. qmin))   
        liq_frac = min(-1.*Cloud_processes%dcond_ls, Cloud_state%ql_upd)/ &
                              max(-1.*Cloud_processes%dcond_ls, qmin)
        ice_frac = 1. - liq_frac
      end where

!-----------------------------------------------------------------------
!    do evaporation of pure ice cloud.
!-----------------------------------------------------------------------
      where ( (Cloud_processes%dcond_ls .lt. 0.) .and.   &
                 (Cloud_state%ql_upd .le. qmin)  .and.            &
                          (Cloud_state%qi_upd .gt. qmin) )
        liq_frac = 0.
        ice_frac = 1.
      end where
        
!-----------------------------------------------------------------------
!    calculate partitioning among liquid and ice to dcond_ls.    
!-----------------------------------------------------------------------
      Cloud_processes%dcond_ls_ice = ice_frac*Cloud_processes%dcond_ls
      Cloud_processes%dcond_ls     = liq_frac*Cloud_processes%dcond_ls   

!------------------------------------------------------------------------
!    output debug diagnostics, if desired.
!------------------------------------------------------------------------
      if (debug) then
        call macrophysics_debug (Cloud_processes, Cloud_state)
      end if

!-------------------------------------------------------------------------
!    if  desired, re-partition dcond_ls based on Bergeron process 
!    calculation. ?????
!-------------------------------------------------------------------------
      IF (rk_repartition_first) THEN
        DO k = 1,kdim
          DO j = 1,jdim
            DO i=1,idim
              Cloud_processes%dcond_ls_tot(i,j,k) =  &
                    Cloud_processes%dcond_ls(i,j,k) +   &
                               Cloud_processes%dcond_ls_ice(i,j,k)
              A_plus_B_s =   &
                ( (HLV/0.024/Input_mp%tin(i,j,k))*   &
                             ((HLV/RVGAS/Input_mp%tin(i,j,k))-1.) ) + &
                (rvgas*Input_mp%tin(i,j,k)*Input_mp%pfull(i,j,k)/ &
                                         2.21/Atmos_state%esat0(i,j,k))
              if ( (Input_mp%tin(i,j,k) < TFREEZE)    &
                  .and. (CLoud_state%ql_upd(i,j,k) + max  &
                      (Cloud_processes%dcond_ls(i,j,k), 0.) > qmin) &
                  .and. (CLoud_state%qa_upd(i,j,k) .gt. qmin))  then 
                D2_dts = dtcloud*Cloud_state%qa_upd(i,j,k)* &
                              ((cfact*1000.*exp((12.96*0.0125*   &
                          (TFREEZE - Input_mp%tin(i,j,k))) - 0.639)/  &
                                Atmos_state%airdens(i,j,k))**(2./3.))* &
                         7.8* ((max(CLoud_state%qi_upd(i,j,k)/  &
                         Cloud_state%qa_upd(i,j,k),1.E-12*cfact*1000.*   &
                          exp((12.96*0.0125*(TFREEZE -   &
                                      Input_mp%tin(i,j,k)))-0.639)   &
                      /Atmos_state%airdens(i,j,k)))**(1./3.))*0.0125*  &
                     (TFREEZE - Input_mp%tin(i,j,k))/((700.**(1./3.))*&
                       A_plus_B_s)
              else
                D2_dts = 0.0        
              end if
              if (Cloud_processes%dcond_ls_tot(i,j,k).ne.0.) then
                dum = D2_dts/Cloud_processes%dcond_ls_tot(i,j,k)
              endif
              dum = max(dum, 0.)
              dum = min(dum, 1.)
              Cloud_processes%dcond_ls_ice(i,j,k) =   &
                                  dum*Cloud_processes%dcond_ls_tot(i,j,k)
              Cloud_processes%dcond_ls(i,j,k) =   &
                           (1. - dum)*Cloud_processes%dcond_ls_tot(i,j,k)
            END DO
          END DO
        END DO
      END IF
 

!----------------------------------------------------------------------
!    return new value of cloud area tendency to permanent location.
!----------------------------------------------------------------------
      Cloud_state%SA_out = SA

!----------------------------------------------------------------------
!    define the relative variance of the subgrid scale cloud water 
!    distribution. For now it is nml-supplied (through Exch_ctrl), may 
!    later be defined based on atmospheric conditions.
!----------------------------------------------------------------------
      Cloud_state%relvarn = qcvar

!-----------------------------------------------------------------------
!    the next step is to compute semi-implicit qa,ql,qi which are used in 
!    the r-k microphysics.  this gives a somewhat implicitness to the
!    scheme. in this calculation an estimate is made of what the cloud 
!    fields would be in the absence of cloud microphysics and cloud 
!    erosion.
!
!    in the case of the Tiedtke cloud scheme, the mean cloud condensate is 
!    incremented if large-scale condensation is occurring. For cloud 
!    fraction, the value from the analytic integration above is used.
!
!    for the statistical cloud scheme these are set equal to the values 
!    diagnosed from the beta-distribution apart from the corrections for 
!    mixed phase clouds.
!-----------------------------------------------------------------------
      if (do_rk_microphys) then
        Cloud_state%qa_mean(:,:,:) = Cloud_state%qa_upd(:,:,:)

        if (.not. do_pdf_clouds) then
          Cloud_state%ql_mean(:,:,:) = Cloud_state%ql_upd(:,:,:) +  &
                             max(Cloud_Processes%dcond_ls(:,:,:)    ,0.)
          Cloud_state%qi_mean(:,:,:) = Cloud_state%qi_upd(:,:,:) +  &
                             max(Cloud_Processes%dcond_ls_ice(:,:,:),0.)
        else
          Cloud_state%ql_mean(:,:,:) = max(Cloud_state%ql_upd(:,:,:) + &
                          Cloud_processes%dcond_ls(:,:,:), qmin)
          Cloud_state%qi_mean(:,:,:) = max(Cloud_state%qi_upd(:,:,:) + &
                          Cloud_processes%dcond_ls_ice(:,:,:), qmin) 
        end if

!------------------------------------------------------------------------
!    define the mean droplet number after this timestep's activation. 
!    distinguish two cases: total activation, in which the percent of
!    available droplets activated is equal to cloud fraction, and partial 
!    activation, in which droplet number increases only by the amount that
!    the cloud area has increased during the timestep. this process only 
!    applies to r-k microphysics. (ncar flavors handle this process within
!    the microphysics.)
!    delta_cf:  A_dt * (1.-qabar)   where A_dt = A*dt , A source rate
!    Eq. 7 of Yi's 2007 paper
!------------------------------------------------------------------------
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              if (total_activation) then
                Cloud_state%qn_mean(i,j,k) =  MAX( 0.,   &
                   Cloud_state%qa_upd(i,j,k)*Particles%drop1(i,j,k)*1.e6/ &
                                               Atmos_state%airdens(i,j,k))
              else if (Cloud_processes%da_ls(i,j,k) > 0.0 ) then
                Cloud_state%qn_mean(i,j,k) =   &
                    Cloud_state%qn_upd(i,j,k) +   &
                        max(Cloud_processes%delta_cf(i,j,k),0.)*  &
                           Particles%drop1(i,j,k)*1.e6/  &
                                               Atmos_state%airdens(i,j,k)
              else
                Particles%drop1(i,j,k) = 0.                
                Cloud_state%qn_mean(i,j,k) = Cloud_state%qn_upd(i,j,k)
              end if
            end do
          end do
        end do
      else  ! do_rk
        Cloud_state%qa_mean(:,:,:) = Cloud_state%qa_upd(:,:,:)

        if (.not. do_pdf_clouds) then
          Cloud_state%ql_mean(:,:,:) = Cloud_state%ql_upd(:,:,:) +  &
                             max(Cloud_Processes%dcond_ls(:,:,:)    ,0.)
          Cloud_state%qi_mean(:,:,:) = Cloud_state%qi_upd(:,:,:) +  &
                             max(Cloud_Processes%dcond_ls_ice(:,:,:),0.)
        else
          Cloud_state%ql_mean(:,:,:) = max(Cloud_state%ql_upd(:,:,:) + &
                          Cloud_processes%dcond_ls(:,:,:), qmin)
          Cloud_state%qi_mean(:,:,:) = max(Cloud_state%qi_upd(:,:,:) + &
                          Cloud_processes%dcond_ls_ice(:,:,:), qmin) 
        end if

!------------------------------------------------------------------------
!    define the mean droplet number after this timestep's activation. 
!    distinguish two cases: total activation, in which the percent of
!    available droplets activated is equal to cloud fraction, and partial 
!    activation, in which droplet number increases only by the amount that
!    the cloud area has increased during the timestep. this process only 
!    applies to r-k microphysics. (ncar flavors handle this process within
!    the microphysics.)
!    delta_cf:  A_dt * (1.-qabar)   where A_dt = A*dt , A source rate
!    Eq. 7 of Yi's 2007 paper
!------------------------------------------------------------------------
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              if (total_activation) then
                Cloud_state%qn_mean(i,j,k) =  MAX( 0.,   &
                   Cloud_state%qa_upd(i,j,k)*Particles%drop1(i,j,k)*1.e6/ &
                                               Atmos_state%airdens(i,j,k))
              else if (Cloud_processes%da_ls(i,j,k) > 0.0 ) then
                Cloud_state%qn_mean(i,j,k) =   &
                    Cloud_state%qn_upd(i,j,k) +   &
                        max(Cloud_processes%delta_cf(i,j,k),0.)*  &
                           Particles%drop1(i,j,k)*1.e6/  &
                                               Atmos_state%airdens(i,j,k)
              else
                Cloud_state%qn_mean(i,j,k) = Cloud_state%qn_upd(i,j,k)
              end if
            end do
          end do
        end do
      endif  ! do_rk

!----------------------------------------------------------------------!


END SUBROUTINE tiedtke_macro


!#########################################################################

subroutine tiedtke_macro_diagnostics (     &
                       idim, jdim, kdim, n_diag_4d, diag_id, diag_4d,    &
                       diag_pt, Cloud_processes, Cloud_state)

!------------------------------------------------------------------------
!    subroutine tiedtke_macro_diagnostics calculates requested diagnostics
!    related to the tiedtke macrophysics.
!------------------------------------------------------------------------

integer,                                 intent(in)    :: idim,jdim,kdim
INTEGER,                                 INTENT(IN)    :: n_diag_4d
type(diag_id_type),                      intent(in)    :: diag_id
type(diag_pt_type),                      intent(in)    :: diag_pt
REAL,     &
  dimension(idim,jdim,kdim,0:n_diag_4d), INTENT(INOUT) :: diag_4d
type(cloud_processes_type),              intent(in)    :: Cloud_processes
type(cloud_state_type),                  intent(in)    :: Cloud_state  

!-----------------------------------------------------------------------
!   local variables:

      integer :: i,j,k

!-----------------------------------------------------------------------
!    save a diagnostic for change in cloud fraction, if desired.
!-----------------------------------------------------------------------
      if ( diag_id%delta_cf > 0 ) then
        diag_4d(:,:,:,diag_pt%delta_cf) = Cloud_processes%delta_cf
      end if

!-----------------------------------------------------------------------
!    compute other diagnostics for cloud fraction.
!-----------------------------------------------------------------------
      if (diag_id%aall > 0)   &
                diag_4d(:,:,:,diag_pt%aall) = Cloud_state%qa_mean(:,:,:)
      if (diag_id%aliq > 0 .or. diag_id%rvolume > 0) then
        where (Cloud_State%ql_mean(:,:,:) .gt. qmin) &
                diag_4d(:,:,:,diag_pt%aliq) = Cloud_state%qa_mean(:,:,:)
      end if
      if (diag_id%aice > 0 .or. diag_id%vfall > 0) then
        where (Cloud_state%qi_mean(:,:,:) .gt. qmin)  &
                 diag_4d(:,:,:,diag_pt%aice) = Cloud_state%qa_mean(:,:,:)
      end if              

      if (diag_id%dcond > 0 )   &
           diag_4d(:,:,:,diag_pt%dcond) = Cloud_processes%dcond_ls(:,:,:)

      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            if (diag_id%potential_droplets > 0 .and.   &
                      Cloud_processes%da_ls(i,j,k) <= 0.0)   &
                 diag_4d(i,j,k,diag_pt%potential_droplets) = 0.0
            if (diag_id%subgrid_w_variance > 0 .and.   &
                      Cloud_processes%da_ls(i,j,k) <= 0.0)   &
                 diag_4d(i,j,k,diag_pt%subgrid_w_variance) = 0.0
          end do
        end do
      end do

!-----------------------------------------------------------------------

end subroutine tiedtke_macro_diagnostics



!########################################################################

SUBROUTINE tiedtke_macro_end

!----------------------------------------------------------------------

      module_is_initialized = .false.

!------------------------------------------------------------------------
!    terminate the modules used by tiedtke_macro_mod.
!------------------------------------------------------------------------
      call polysvp_end
      if ( do_pdf_clouds ) then
        call beta_dist_end
      end if

!-----------------------------------------------------------------------

END SUBROUTINE tiedtke_macro_end



!########################################################################

SUBROUTINE tiedtke_macro_nopdf_nosuper (    &
            idim, jdim, kdim, C2ls_mp, Input_mp, Atmos_state,   &
            Cloud_state, ST, SQ, Cloud_processes, Particles, n_diag_4d, &
            diag_4d, diag_id, diag_pt, edum, SA, U00p, erosion_scale)

!------------------------------------------------------------------------
!   subroutine tiedtke_macro_nopdf_nosuper calculates non-convective
!   condensation while not allowing supersaturation in the clear-sky 
!   portion of grid boxes (as in original Tiedtke scheme).
!------------------------------------------------------------------------

INTEGER,                         INTENT(IN )   :: idim, jdim, kdim      
type(mp_input_type),             intent(inout) :: Input_mp   
type(mp_conv2ls_type),           intent(inout) :: C2ls_mp   
type(atmos_state_type),          intent(inout) :: Atmos_state
type(cloud_state_type),          intent(inout) :: Cloud_state
type(cloud_processes_type),      intent(inout) :: Cloud_processes
type(particles_type),            intent(inout) :: Particles
REAL, dimension(idim,jdim,kdim), INTENT(IN)    :: ST, SQ, edum, U00p,  &
                                                  erosion_scale
REAL, dimension(idim,jdim,kdim), INTENT(INOUT) :: SA
INTEGER,                         INTENT(IN )   :: n_diag_4d
REAL,  &
 dimension(idim, jdim, kdim, 0:n_diag_4d),   &
                                 INTENT(INOUT) :: diag_4d
TYPE(diag_id_type),              intent(in)    :: diag_id
TYPE(diag_pt_type),              intent(in)    :: diag_pt

!-------------------------------------------------------------------------
!----local variables----

!-------------------------------------------------------------------------
!       dqs_ls         change in qs due to large 
!                      scale processes
!       A_dt           product of A and dtcloud in     dimensionless in 
!                      in the analytic integration     qa integration
!                      of the qa equation, or C and
!                      dtcloud in the analytic         kg condensate/
!                      integration of the ql and qi    kg air in ql or 
!                      equations.                      qi integration
!
!       B_dt           product of B and dtcloud in     dimensionless in
!                      in the analytic integration     qa, ql, and qi
!                      of the qa equation, or D and    integration
!                      dtcloud in the analytic         
!                      integration of the ql and qi    
!                      equations.                      
!
!       qa0            value of cloud fraction or      dimensionless or
!                      cloud condensate at the         kg condensate /
!                      initial time                    kg air
!        
!       qa1            value of cloud fraction or      dimensionless or
!                      cloud condensate at the final   kg condensate /
!                      time                            kg air
!
!       qaeq           equilibrium value of cloud      dimensionless or
!                      fraction or cloud condensate    kg condensate /
!                      that the analytic integration   kg air
!                      approaches                   
!       qabar          mean value of cloud fraction    dimensionless or
!                      or cloud condensate over the    kg condensate /
!                      t0 to t0 + dtcloud interval     kg air
!
!------------------------------------------------------------------------
     
      real, dimension(idim, jdim,kdim)   :: dqs_ls
      real, dimension(idim, jdim,kdim)   :: A_dt, B_dt
      real, dimension(idim, jdim,kdim)   :: qa1, qa0, qabar
      real, dimension(idim, jdim,kdim)   :: qaeq
      real, dimension(idim, jdim,kdim)   :: tmp1
      REAL                               :: eslt, qvs, qs_d, qvi, esit, ul
      REAL                               :: ttmp, qtmp, qs_t, dqsdT1,  &
                                            qvmax, esat0, gamma1,  &
                                            tmp1s, qs_l, qs_i
      INTEGER                            :: ns, id
      INTEGER                            :: i,j,k

!-----------------------------------------------------------------------
!
!       The first step is to compute the change in qs due to large-
!       scale processes, dqs_ls.   In Tiedtke, it has contributions from 
!       large-scale uplift, convection induced compensating subsidence,
!       turbulence cooling and radiative cooling.  dqs_ls has the form:
!
!               (((omega+ grav*Mc)/airdens/cp)+radturbten)*dqsdT*dtcloud
!   (6) dqs_ls= --------------------------------------------------------
!                  1.  +   ( qa +  (da_ls/2.) ) * gamma
!
!       Here da_ls is the increase in cloud fraction due to non-
!       convective processes.  Because this increase is also a function
!       of dqs_ls, a quadratic equation must be solved for dqs_ls in
!       the case that da_ls is not equal to zero.
!
!------------------------------------------------------------------------
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            dqs_ls(i,j,k) = ((                    &
              (Input_mp%omega(i,j,k) + grav*C2ls_mp%mc_full(i,j,k))/  &
                                Atmos_state%airdens(i,j,k)/cp_air)  +   &
                         Input_mp%radturbten(i,j,k))* &
                                         dtcloud*Atmos_state%dqsdT(i,j,k)
            if (dqs_ls(i,j,k) .le. 0. .and.    &
                       Atmos_state%U_ca(i,j,k) .ge. U00p(i,j,k) .and.   &
                                Cloud_state%qa_upd(i,j,k) .lt. 1.)  then
              tmp1(i,j,k) = sqrt( (1. + Cloud_state%qa_upd(i,j,k)*  &
                                       Atmos_state%gamma(i,j,k))**2. -  &
                                  (1. - Cloud_state%qa_upd(i,j,k))*  &
                                     (1. - Cloud_state%qa_upd(i,j,k))*&
                                 Atmos_state%gamma(i,j,k)*dqs_ls(i,j,k)/&
                                          Atmos_state%qs(i,j,k)/      &
                    max(1.-Atmos_state%U_ca(i,j,k), qmin) ) -   &
                        (1. + Cloud_state%qa_upd(i,j,k)*  &
                                               Atmos_state%gamma(i,j,k))
              tmp1(i,j,k) = -1.*tmp1(i,j,k)/((1. -   &
                                       Cloud_state%qa_upd(i,j,k))*   &
                               (1. - Cloud_state%qa_upd(i,j,k))*  &
                                            Atmos_state%gamma(i,j,k)/&
                                Atmos_state%qs(i,j,k)/   &
                            max(1. - Atmos_state%U_ca(i,j,k), qmin)/2.)
              dqs_ls(i,j,k) = min(tmp1(i,j,k),dqs_ls(i,j,k)/(1. +   &
                              0.5*(1. + Cloud_State%qa_upd(i,j,k))*  &
                                             Atmos_state%gamma(i,j,k)))
            else
              dqs_ls(i,j,k) = dqs_ls(i,j,k)/(1. +   &
                   Cloud_state%qa_upd(i,j,k)*Atmos_state%gamma(i,j,k))
            endif
          end do
        end do
      end do

!------------------------------------------------------------------------
!       The next step is to compute the change in saturated volume
!       fraction due to non-convective condensation, da_ls.   This 
!       occurs in two conditions:
!
!       (a)  dqs_ls < 0. and U00 < U < 1., where U00 is the threshold
!            relative humidity for non-convective condensation. Note 
!            that if U is greater than or equal to 1., ideally qa = 1,
!            and da_ls = 0.  However this may not be the case for 
!            numerical reasons so this must be assured after analytic 
!            integration of the qa equation.
!
!            For these cases the change in saturated volume fraction is:
!
!   (7)      da_ls = - (1.-qa)*(1.-qa)*dqs_ls/2./qs/(1.-U)
!
!            This formula arises from the assumption that vapor is uni-
!            formly distributed in the range [qv_clr - (qs - qv_clr),qs]
!            where qv_clr is the amount of vapor in the unsaturated 
!            volume and is given from the following equation:
!
!   (8)      qv  =   qa * qs      +   (1.-qa) * qv_clr
!          
!            Implicit in equation (7) is the following assumption:
!            As qsat changes, the distribution of qv+ql+qi 
!            remains constant.  That is as qsat rises, portions where
!            qv+ql+qi > qsat+dqsat remain saturated.  This can only
!            occur if it is assumed that ql+qi evaporate-sublimate or
!            condense-deposit to keep qv = qsat. 
!
!       (b)  dqs_ls > 0.  Ideally some portion of the cloud should
!            evaporate however this is not accounted for at present.
!            
!------------------------------------------------------------------------
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            if ((dqs_ls(i,j,k) .le. 0. .and.    &
               Atmos_state%U_ca(i,j,k) .ge. U00p(i,j,k)) .and. &
                (Cloud_state%qa_upd(i,j,k) +  &
                  C2ls_mp%convective_humidity_area(i,j,k) .le. 1.))   then
              Cloud_processes%da_ls(i,j,k) =    &
                      -0.5*(1. - Cloud_state%qa_upd(i,j,k) -  &
                               C2ls_mp%convective_humidity_area(i,j,k))*  &
                            (1. - Cloud_state%qa_upd(i,j,k) -    &
                              C2ls_mp%convective_humidity_area(i,j,k) )* &
                                   dqs_ls(i,j,k)/Atmos_state%qs(i,j,k) /  &
                                  max(1.-Atmos_state%U_ca(i,j,k), qmin)
            else
              Cloud_processes%da_ls(i,j,k) = 0.
            endif
          end do
        end do
      end do

!------------------------------------------------------------------------
!       Turbulent erosion of clouds
!
!       As in Tiedtke (1993) this is calculated using the eros_scale
!       parameter as:
!
!   (9) dql/dt    =  - qa * eros_scale * (qs - qv) * (ql/ ql+qi )
!
!  (10) dqi/dt    =  - qa * eros_scale * (qs - qv) * (qi/ ql+qi )
!
!  (11) dqa/dt    =  - qa * eros_scale * (qs - qv) * (qa/ ql+qi )
!
!       for which the erosion sink term (B in equation 13) is
!
!  (12) B = qa * eros_scale * (qs - qv) / (ql + qi)  
!
!------------------------------------------------------------------------
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            if (Cloud_state%ql_upd(i,j,k) .gt. qmin .or.   &
                        Cloud_state%qi_upd (i,j,k).gt. qmin) then
              Cloud_processes%D_eros(i,j,k) =   &
                    Cloud_state%qa_upd(i,j,k) * erosion_scale(i,j,k) *   &
                               dtcloud * Atmos_state%qs(i,j,k) *      &
                      (1.-Atmos_state%U_ca(i,j,k)) /   &
                   (Cloud_state%qi_upd(i,j,k) + CLoud_state%ql_upd(i,j,k))

              if (Input_mp%pfull(i,j,k) .gt. 400.e02) then
                Cloud_processes%D_eros(i,j,k) =  &
                     Cloud_processes%D_eros(i,j,k) + efact*  &
                                      Cloud_processes%D_eros(i,j,k)* &
                          ((Input_mp%pfull(i,j,kdim) -   &
                                               Input_mp%pfull(i,j,k))/  &
                                (Input_mp%pfull(i,j,kdim) - 400.e02))
              else
                Cloud_processes%D_eros(i,j,k)=  &
                      Cloud_processes%D_eros(i,j,k) +  &
                                     efact*CLoud_processes%D_eros(i,j,k)

              endif
            else
              Cloud_processes%D_eros(i,j,k) = 0.
            endif
          END DO    
        END DO    
      END DO    

!------------------------------------------------------------------------
!       The next step is to analytically integrate the saturated volume
!       fraction equation.  This follows the Tiedtke approach
!
!       The qa equation is written in the form:
!
!  (13) dqa/dt    =   (1.-qa) * A   -  qa * B 
!
!       Note that over the physics time step, A, B are assumed to be 
!       constants.
!
!       Defining qa(t) = qa0 and qa(t+dtcloud) = qa1, the analytic
!       solution of the above equation is:
!
!  (14) qa1 = qaeq -  (qaeq - qa0) * exp (-(A+B)*dtcloud)
! 
!       where qaeq is the equilibrium cloud fraction that is approached
!       with an time scale of 1/(A+B),
!
!  (15) qaeq  =  A/(A+B)
!
!
!       To diagnose the magnitude of each of the right hand terms of
!       (13) integrated over the time step, define the average cloud
!       fraction in the interval t to t + dtcloud qabar as:
!
!  (16) qabar  = qaeq - [ (qa1-qa0) / ( dtcloud * (A+B) ) ]
! 
!       from which the magnitudes of the A and B terms integrated
!       over the time step are:
!
!       A * (1-qabar)    and    -B * (qabar)
!
!       Additional notes on this analytic integration:
!
!       1.   For large-scale cloud formation or destruction from 
!            the dqs_ls term the contributions to A or B are defined
!            from:
!
!  (19)      A_ls * (1. - qa) = da_ls / dtcloud      if da_ls >= 0.
! 
!  (20)      B_ls * qa        = da_ls / dtcloud      if da_ls < 0.
!
!
!       3.   Qa goes to zero only in the case of ql and qi less than or
!            equal to     qmin; see 'cloud destruction code' near the end of 
!            this loop over levels.
!------------------------------------------------------------------------
        
!------------------------------------------------------------------------
!    compute A_dt; This is assigned to the large-scale source term
!    following (18). Reset B_dt.
!------------------------------------------------------------------------
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            A_dt(i,j,k) = Cloud_processes%da_ls(i,j,k)/   &
                              max((1.-Cloud_state%qa_upd(i,j,k)), qmin)
            B_dt(i,j,k) = Cloud_processes%D_eros(i,j,k)
  
!------------------------------------------------------------------------
!    do analytic integration.      
!------------------------------------------------------------------------
            if ( (A_dt(i,j,k) .gt. Dmin) .or.   &
                 (B_dt(i,j,k) .gt. Dmin) )  then
              qa0(i,j,k) = Cloud_state%qa_upd(i,j,k)
              qaeq(i,j,k) = A_dt(i,j,k)/(A_dt(i,j,k)  + B_dt(i,j,k))
              qa1(i,j,k) = qaeq(i,j,k) - (qaeq(i,j,k) - qa0(i,j,k)) * &
                                      exp ( -1.*(A_dt(i,j,k)+B_dt(i,j,k)) )
              qabar(i,j,k) = qaeq(i,j,k) - ((qa1(i,j,k) - qa0(i,j,k))/  &
                                               (A_dt(i,j,k) + B_dt(i,j,k)))
            else
              qa0(i,j,k)   = Cloud_state%qa_upd(i,j,k)
              qaeq(i,j,k)  = qa0(i,j,k)
              qa1(i,j,k)   = qa0(i,j,k)   
              qabar(i,j,k) = qa0(i,j,k)  
            endif 
          END DO    
        END DO    
      END DO    

      do k=1,kdim
        do j=1,jdim
          do i=1,idim
!------------------------------------------------------------------------
!    save some diagnostics.
!------------------------------------------------------------------------
            if ( (A_dt(i,j,k) .gt.     Dmin) .and.   &
                 (B_dt(i,j,k) .gt.     Dmin) )  then
              if (diag_id%qadt_lsform + diag_id%qa_lsform_col > 0) then
                diag_4d(i,j,k,diag_pt%qadt_lsform) =  A_dt(i,j,k)*  &
                                            (1.-qabar(i,j,k))*inv_dtcloud
              end if
              if ( diag_id%qadt_eros + diag_id%qa_eros_col > 0 ) then
                diag_4d(i,j,k,diag_pt%qadt_eros)  =    &
                     ((qa1(i,j,k) - qa0(i,j,k))*inv_dtcloud )- &
                                     diag_4d(i,j,k,diag_pt%qadt_lsform)   
              end if
            else if (A_dt(i,j,k) .gt. Dmin) then
              if (diag_id%qadt_lsform + diag_id%qa_lsform_col > 0) then
                diag_4d(i,j,k,diag_pt%qadt_lsform) =   &
                          (qa1(i,j,k) - qa0(i,j,k))*inv_dtcloud
              end if 
            else if (B_dt(i,j,k) .gt. Dmin)  then
              if ( diag_id%qadt_eros + diag_id%qa_eros_col > 0 ) then
                diag_4d(i,j,k,diag_pt%qadt_eros)  =   &
                          (qa1(i,j,k) - qa0(i,j,k))*inv_dtcloud
              end if 
            endif
          END DO    
        END DO    
      END DO    

!-----------------------------------------------------------------------
!   define value needed for diagnostic calculation.
!-----------------------------------------------------------------------
      if (diag_id%qadt_ahuco + diag_id%qa_ahuco_col > 0) then
        diag_4d(:,:,:,diag_pt%qadt_ahuco) = qa1(:,:,:)
      end if

!------------------------------------------------------------------------
!    limit cloud area to be no more than that which is not being
!    taken by convective clouds.
!------------------------------------------------------------------------
      if (limit_conv_cloud_frac) then
        qa1 = MIN(qa1, 1.0 -C2ls_mp%convective_humidity_area)
      endif
                 
!------------------------------------------------------------------------
!    complete calculation of diagnostic.
!------------------------------------------------------------------------
      if (diag_id%qadt_ahuco + diag_id%qa_ahuco_col > 0) then
        diag_4d(:,:,:,diag_pt%qadt_ahuco) =    &
                  (qa1(:,:,:) - diag_4d(:,:,:,diag_pt%qadt_ahuco))* & 
                                                              inv_dtcloud 
      end if

!------------------------------------------------------------------------
!    set total tendency term and update cloud fraction    
!------------------------------------------------------------------------
      SA = (SA + qa1) - qa0
      Cloud_state%qa_upd = qa1
        
!------------------------------------------------------------------------
!       The next step is to calculate the change in condensate
!       due to non-convective condensation, dcond_ls. Note that this is
!       not the final change but is used only to apportion condensate
!       change between phases. According to Tiedtke 1993 this takes the
!       form:
!
!  (21) dcond_ls = -1. * (qa +  0.5*da_ls) * dqs_ls
!
!       Here the 0.5*da_ls represents using a midpoint cloud fraction.
!       This is accomplished by using the variable qabar.
!
!       Also, save term needed when predicted droplets is active (tmp5).
!------------------------------------------------------------------------
      IF (use_qabar) THEN
        Cloud_processes%dcond_ls = -1.*qabar*dqs_ls
        Cloud_processes%delta_cf = A_dt*(1.-qabar)
      ELSE
        Cloud_processes%dcond_ls = -1.*Cloud_state%qa_upd*dqs_ls   
        Cloud_processes%delta_cf = A_dt*(1.-Cloud_state%qa_upd)
      END IF      

!-----------------------------------------------------------------------


end subroutine tiedtke_macro_nopdf_nosuper 


!########################################################################

SUBROUTINE tiedtke_macro_nopdf_super (   &
                   idim, jdim, kdim, C2ls_mp,Input_mp, Atmos_state,    &
                   Cloud_state, ST, SQ, Cloud_processes, Particles,  &
                   n_diag_4d, diag_4d, diag_id, diag_pt, edum, SA,   &
                                                      U00p, erosion_scale) 

!------------------------------------------------------------------------
!   subroutine tiedtke_macro_nopdf_super calculates non-convective
!   condensation while allowing ice supersaturation in the clear-sky 
!   portion of grid boxes, as is needed in order to activate ice crystals 
!   (see Tompkins et al 2007, Salzmann et al , 2010).
!------------------------------------------------------------------------

INTEGER,                            INTENT(IN )    :: idim, jdim, kdim  
type(atmos_state_type),             intent(inout)  :: Atmos_state
type(mp_input_type),                intent(in)     :: Input_mp   
type(mp_conv2ls_type),              intent(in)     :: C2ls_mp   
type(cloud_state_type),             intent(inout)  :: Cloud_state
type(cloud_processes_type),         intent(inout)  :: Cloud_processes
type(particles_type),               intent(inout)  :: Particles
REAL, dimension(idim,jdim,kdim),    INTENT(IN )    :: ST, SQ, edum, U00p,&
                                                      erosion_scale
REAL, dimension(idim,jdim,kdim),    INTENT(INout ) :: SA
INTEGER,                            INTENT(IN )    :: n_diag_4d
REAL, dimension(idim, jdim, kdim, 0:n_diag_4d),        &
                                    INTENT(INOUT ) :: diag_4d
TYPE(diag_id_type),                 intent(in)     :: diag_id
TYPE(diag_pt_type),                 intent(in)     :: diag_pt

!-------------------------------------------------------------------------
!-----local variables-----
! 
!       A_dt           product of A and dtcloud in     dimensionless in 
!                      in the analytic integration     qa integration
!                      of the qa equation, or C and
!                      dtcloud in the analytic         kg condensate/
!                      integration of the ql and qi    kg air in ql or 
!                      equations.                      qi integration
!
!       B_dt           product of B and dtcloud in     dimensionless in
!                      in the analytic integration     qa, ql, and qi
!                      of the qa equation, or D and    integration
!                      dtcloud in the analytic         
!                      integration of the ql and qi    
!                      equations.                      
!
!       qa0            value of cloud fraction or      dimensionless or
!                      cloud condensate at the         kg condensate /
!                      initial time                    kg air
!        
!       qa1            value of cloud fraction or      dimensionless or
!                      cloud condensate at the final   kg condensate /
!                      time                            kg air
!
!       qaeq           equilibrium value of cloud      dimensionless or
!                      fraction or cloud condensate    kg condensate /
!                      that the analytic integration   kg air
!                      approaches                   
!       qabar          mean value of cloud fraction    dimensionless or
!                      or cloud condensate over the    kg condensate /
!                      t0 to t0 + dtcloud interval     kg air
!       qve            Tiedtke environmenmtal 
!                      humidity
!
!------------------------------------------------------------------------
      real, dimension(idim, jdim,kdim) :: deltpg, dqs_ls, drhcqs_ls, dum,&
                                          qve 
      real, dimension(idim, jdim,kdim) :: A_dt, B_dt, qa1, qa0,  &
                                          qabar, qaeq, qa_t, tmp0, tmp1, &
                                          drhcqsdT, beta, ttmp, qtmp,  &
                                          qs_t, qs_l, qs_i, dqsdT1, gamma1
      real, dimension(idim,jdim)       :: qagtmp,qcgtmp,qvgtmp
      REAL                             :: eslt, qvs, qs_d, qvi, esit, ul
      REAL                             :: qvmax, esat0, tmp1s            
      INTEGER                          :: ns, id
      INTEGER                          :: i,j,k

!------------------------------------------------------------------------
!    calculate the derivative of the threshold relative humidity over water!    (for temps > 233) and over ice (temps < 233) with respect to 
!    temperature, and the L/Cp drhqsdT term (beta) in the
!    temperature equation to be used in this calculation. for
!    the water case (T >= 233), no supersaturation is allowed. For the ice
!    case, superstaurations up to rh_crit are permitted.
!------------------------------------------------------------------------
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            if (Input_mp%tin(i,j,k) .GE. TFREEZE - 40. ) then 
              drhcqsdT(i,j,k) =  HLV*Atmos_state%qsl(i,j,k)/  &
                                        (RVGAS*Input_mp%tin(i,j,k)**2)
              beta(i,j,k) = drhcqsdT(i,j,k)*HLV/CP_AIR
            else
              drhcqsdT(i,j,k) = Atmos_state%rh_crit(i,j,k)*   &
                          HLS*Atmos_state%qsi(i,j,k)/  &
                                       (RVGAS*Input_mp%tin(i,j,k)**2) +  &
                           Particles%hom(i,j,k)*Atmos_state%qsi(i,j,k)*  &
                       (2.*0.0073*(Input_mp%tin(i,j,k) - TFREEZE) + 1.466)
              beta(i,j,k) = drhcqsdT(i,j,k)*HLS/CP_AIR
            endif 
          end do
        end do
      end do

      do k=1,kdim
        do j=1,jdim
          do i=1,idim

!------------------------------------------------------------------------
!    calculate environmental specific humidity (qve) based on input cloud 
!    fraction. take into account convective cloud presence.
!------------------------------------------------------------------------
            qa_t(i,j,k) =    &
              MAX(MIN(Cloud_State%qa_upd(i,j,k) +   &
                          C2ls_mp%convective_humidity_area(i,j,k), 1.),0.)
            qve(i,j,k) =  (Input_mp%qin(i,j,k) - qa_t(i,j,k)*  &
                            Atmos_state%qs(i,j,k) ) /   &
                                       MAX((1. - qa_t(i,j,k) ), qmin)

!------------------------------------------------------------------------
!       The first step is to compute the change in qs due to large-
!       scale processes, dqs_ls.   In Tiedtke, it has contributions from 
!       large-scale uplift, convection induced compensating subsidence,
!       turbulence cooling and radiative cooling.  dqs_ls has the form:
!
!               (((omega+ grav*Mc)/airdens/cp)+radturbten)*dqsdT*dtcloud
!   (6) dqs_ls= --------------------------------------------------------
!                  1.  +   ( qa +  (da_ls/2.) ) * gamma
!
!       Here da_ls is the increase in cloud fraction due to non-
!       convective processes.  Because this increase is also a function
!       of dqs_ls, a quadratic equation must be solved for dqs_ls in
!       the case that da_ls is not equal to zero.
!------------------------------------------------------------------------
            dum(i,j,k) = (((Input_mp%omega(i,j,k) + grav*  &
                             C2ls_mp%mc_full(i,j,k))/ &
                                 Atmos_state%airdens(i,j,k)/cp_air) + &
                                   Input_mp%radturbten(i,j,k))*dtcloud
            dqs_ls(i,j,k) = dum(i,j,k) * Atmos_state%dqsdT(i,j,k)
            drhcqs_ls(i,j,k) = dum(i,j,k) *drhcqsdT(i,j,k)
          end do
        end do
      end do

 !------------------------------------------------------------------------
 !       The next step is to compute the change in saturated volume
 !       fraction due to non-convective condensation, da_ls.   This 
 !       occurs in two conditions:
 !
 !       (a)  dqs_ls < 0. and U00 < U < 1., where U00 is the threshold
 !            relative humidity for non-convective condensation. Note 
 !            that if U is greater than or equal to 1., ideally qa = 1,
 !            and da_ls = 0.  However this may not be the case for 
 !            numerical reasons so this must be assured after analytic 
 !            integration of the qa equation.
 !
 !            For these cases the change in saturated volume fraction is:
 !
 !   (7)      da_ls = - (1.-qa)*(1.-qa)*dqs_ls/2./qs/(1.-U)
 !
 !            This formula arises from the assumption that vapor is uni-
 !            formly distributed in the range [qv_clr - (qs - qv_clr),qs]
 !            where qv_clr is the amount of vapor in the unsaturated 
 !            volume and is given from the following equation:
 !
 !   (8)      qv  =   qa * qs      +   (1.-qa) * qv_clr
 !          
 !            Implicit in equation (7) is the following assumption:
 !            As qsat changes, the distribution of qv+ql+qi 
 !            remains constant.  That is as qsat rises, portions where
 !            qv+ql+qi > qsat+dqsat remain saturated.  This can only
 !            occur if it is assumed that ql+qi evaporate-sublimate or
 !            condense-deposit to keep qv = qsat. 
 !
 !       (b)  dqs_ls > 0.  Ideally some portion of the cloud should
 !            evaporate however this is not accounted for at present.
 !----------------------------------------------------------------------- 
      where (dqs_ls .le. 0. .and.    &
                  qve .ge. U00p*Atmos_state%rh_crit_min*Atmos_state%qs &
                                                       .and. qa_t .lt. 1.)
        tmp0 =  (1 + Atmos_state%gamma*qa_t)*drhcqs_ls - beta*qa_t* &
                                                                    dqs_ls
        tmp1 = SQRT( (1.+qa_t*Atmos_state%gamma)**2. - (1.-qa_t)* &
                         beta * tmp0 /MAX( (  Atmos_state%rh_crit *  &
                          Atmos_state%qs - qve ), qmin ) ) -  &
                                            (1.+qa_t*Atmos_State%gamma) 
        tmp1 = -1.*tmp1/((1. - qa_t)*beta/(2.*MAX(   &
                   (Atmos_state%rh_crit*Atmos_State%qs - qve), qmin )) )
        drhcqs_ls = min(tmp1, tmp0 / MAX(1.+ Atmos_state%gamma*qa_t +&
                                         beta*0.5*(1. - qa_t), qmin))
        Cloud_processes%da_ls = -0.5*(1. - qa_t)*drhcqs_ls/   &
             MAX( Atmos_state%rh_crit*Atmos_state%qs - qve, qmin)
        dqs_ls = (dqs_ls - Atmos_state%gamma*0.5*  &
                  Cloud_processes%da_ls*drhcqs_ls)/   &
                                          (1. + qa_t*Atmos_state%gamma)
      elsewhere
        drhcqs_ls = 0.
        Cloud_processes%da_ls = 0.
        dqs_ls = dqs_ls/(1. + qa_t*Atmos_state%gamma)
      endwhere

!------------------------------------------------------------------------
!       Turbulent erosion of clouds
!
!       As in Tiedtke (1993) this is calculated using the eros_scale
!       parameter as:
!
!   (9) dql/dt    =  - qa * eros_scale * (qs - qv) * (ql/ ql+qi )
!
!  (10) dqi/dt    =  - qa * eros_scale * (qs - qv) * (qi/ ql+qi )
!
!  (11) dqa/dt    =  - qa * eros_scale * (qs - qv) * (qa/ ql+qi )
!
!       for which the erosion sink term (B in equation 13) is
!
!  (12) B = qa * eros_scale * (qs - qv) / (ql + qi)  
!
!------------------------------------------------------------------------
      DO k=1,kdim
        DO j=1,jdim
          DO i=1,idim
            if (Cloud_state%ql_upd(i,j,k) .gt. qmin .or. &
                           Cloud_state%qi_upd (i,j,k).gt. qmin) then
              Cloud_processes%D_eros(i,j,k) =    &
                        Cloud_state%qa_upd(i,j,k)*erosion_scale(i,j,k)*   &
                                                           dtcloud*    &
                   (MAX(Atmos_state%qs(i,j,k) - Input_mp%qin(i,j,k),&
                         0.)/(1. + Atmos_state%gamma(i,j,k)) ) /      &
                    (Cloud_state%qi_upd(i,j,k) + Cloud_state%ql_upd(i,j,k))

              if (Input_mp%pfull(i,j,k) .gt. 400.e02) then
                Cloud_processes%D_eros(i,j,k) =   &
                     Cloud_processes%D_eros(i,j,k) + efact*  &
                                Cloud_processes%D_eros(i,j,k)*  &
                    ((Input_mp%pfull(i,j,kdim) - Input_mp%pfull(i,j,k))/  &
                                (Input_mp%pfull(i,j,kdim) - 400.e02))

              else
                Cloud_processes%D_eros(i,j,k) =  &
                     Cloud_processes%D_eros(i,j,k) + efact*  &
                                          Cloud_processes%D_eros(i,j,k)
              endif
            else
              Cloud_processes%D_eros(i,j,k) = 0.
            endif
          END DO    
        END DO    
      END DO    

!------------------------------------------------------------------------ 
!       The next step is to analytically integrate the saturated volume
!       fraction equation.  This follows the Tiedtke approach
!
!       The qa equation is written in the form:
!
!  (13) dqa/dt    =   (1.-qa) * A   -  qa * B 
!
!       Note that over the physics time step, A, B are assumed to be 
!       constants.
!
!       Defining qa(t) = qa0 and qa(t+dtcloud) = qa1, the analytic
!       solution of the above equation is:
!
!  (14) qa1 = qaeq -  (qaeq - qa0) * exp (-(A+B)*dtcloud)
! 
!       where qaeq is the equilibrium cloud fraction that is approached
!       with an time scale of 1/(A+B),
!
!  (15) qaeq  =  A/(A+B)
!
!
!       To diagnose the magnitude of each of the right hand terms of
!       (13) integrated over the time step, define the average cloud
!       fraction in the interval t to t + dtcloud qabar as:
!
!  (16) qabar  = qaeq - [ (qa1-qa0) / ( dtcloud * (A+B) ) ]
! 
!       from which the magnitudes of the A and B terms integrated
!       over the time step are:
!
!       A * (1-qabar)    and    -B * (qabar)
!
!       Additional notes on this analytic integration:
!
!       1.   For large-scale cloud formation or destruction from 
!            the dqs_ls term the contributions to A or B are defined
!            from:
!
!  (19)      A_ls * (1. - qa) = da_ls / dtcloud      if da_ls >= 0.
! 
!  (20)      B_ls * qa        = da_ls / dtcloud      if da_ls < 0.
!
!
!       3.   Qa goes to zero only in the case of ql and qi less than or
!            equal to     qmin; see 'cloud destruction code' near the end of 
!            this loop over levels.
!------------------------------------------------------------------------
        
!-------------------------------------------------------------------------
!    compute A_dt; This is assigned to the large-scale source term
!    following (18). Reset B_dt.
!-------------------------------------------------------------------------
   
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            A_dt(i,j,k) = Cloud_processes%da_ls(i,j,k)/   &
                                            max((1.-qa_t(i,j,k)), qmin)
            B_dt(i,j,k) = Cloud_processes%D_eros(i,j,k)
            if ( (A_dt(i,j,k) .gt. Dmin) .or.   &
                                      (B_dt(i,j,k) .gt. Dmin) ) then 
              qa0(i,j,k)   = Cloud_state%qa_upd(i,j,k)
              qaeq(i,j,k)  =                                     &
                                 A_dt(i,j,k)/(A_dt(i,j,k) + B_dt(i,j,k))
              qa1(i,j,k)  = qaeq(i,j,k) - (qaeq(i,j,k) - qa0(i,j,k))* &
                                     exp(-1.*(A_dt(i,j,k) + B_dt(i,j,k)) )
              qabar(i,j,k) = qaeq(i,j,k) - ((qa1(i,j,k) - qa0(i,j,k))/  &
                                         (A_dt(i,j,k) + B_dt(i,j,k)))
            else
              qa0(i,j,k)   = Cloud_state%qa_upd(i,j,k)
              qaeq(i,j,k)  = qa0(i,j,k)   
              qa1(i,j,k)   = qa0(i,j,k)   
              qabar(i,j,k) = qa0(i,j,k)  
            endif
          end do
        end do
      end do

!-------------------------------------------------------------------------
!    output some diagnostics.
!-------------------------------------------------------------------------
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            if ( (A_dt(i,j,k) .gt. Dmin) .and.   &
                 (B_dt(i,j,k) .gt. Dmin) )  then
              if (diag_id%qadt_lsform + diag_id%qa_lsform_col > 0) then
                diag_4d(i,j,k,diag_pt%qadt_lsform) =  A_dt(i,j,k)*  &
                              (1. - qabar(i,j,k))*inv_dtcloud
              end if
              if ( diag_id%qadt_eros + diag_id%qa_eros_col > 0 ) then
                diag_4d(i,j,k,diag_pt%qadt_eros)  =   &
                     ((qa1(i,j,k) - qa0(i,j,k))*inv_dtcloud ) - &
                                       diag_4d(i,j,k,diag_pt%qadt_lsform)  
              end if
            else if (A_dt(i,j,k) .gt. Dmin) then
              if (diag_id%qadt_lsform + diag_id%qa_lsform_col > 0) then
                diag_4d(i,j,k,diag_pt%qadt_lsform) =    &
                          (qa1(i,j,k) - qa0(i,j,k))*inv_dtcloud
              end if 
            else if (B_dt(i,j,k) .gt. Dmin)  then
              if ( diag_id%qadt_eros + diag_id%qa_eros_col > 0 ) then
                diag_4d(i,j,k,diag_pt%qadt_eros)  =    &
                           (qa1(i,j,k) - qa0(i,j,k))*inv_dtcloud
              end if 
            endif
          end do
        end do
      end do

!-----------------------------------------------------------------------
!   define value needed for diagnostic calculation.
!-----------------------------------------------------------------------
      if ( diag_id%qadt_ahuco + diag_id%qa_ahuco_col > 0 ) then
        diag_4d(:,:,:,diag_pt%qadt_ahuco)  = qa1
      end if

!-------------------------------------------------------------------------
!    limit cloud area to be no more than that which is not being
!    taken by convective clouds.
!-------------------------------------------------------------------------
      if (limit_conv_cloud_frac) then
        qa1 = MIN(qa1, 1.0 - C2ls_mp%convective_humidity_area)
      endif

!------------------------------------------------------------------------
!    complete calculation of diagnostic.
!------------------------------------------------------------------------
      if ( diag_id%qadt_ahuco + diag_id%qa_ahuco_col > 0 ) then
        diag_4d(:,:,:,diag_pt%qadt_ahuco)  =   &
             (qa1(:,:,:) - diag_4d(:,:,:,diag_pt%qadt_ahuco))*inv_dtcloud 
      end if
                 
!-------------------------------------------------------------------------
!    set total tendency term and update cloud fraction.    
!-------------------------------------------------------------------------
      SA = (SA + qa1) - qa0
      Cloud_state%qa_upd = qa1
      Cloud_processes%delta_cf = MAX(qa1 - qa0 , 0.)

!-------------------------------------------------------------------------
!       The next step is to calculate the change in condensate
!       due to non-convective condensation, dcond_ls. Note that this is
!       not the final change but is used only to apportion condensate
!       change between phases. According to Tiedtke 1993 this takes the
!       form:
!
!  (21) dcond_ls = -1. * (qa +  0.5*da_ls) * dqs_ls
!
!       Here the 0.5*da_ls represents using a midpoint cloud fraction.
!       This is accomplished by using the variable qabar.

!cms but qabar is not limited to area outside  conv. clouds ...
!!      dcond_ls = -1. * MIN(qabar, 1.- ahuco)  * dqs_ls
!-----------------------------------------------------------------------
      Cloud_processes%dcond_ls = -1.*qa_t*dqs_ls -    &
           0.5*MAX(MIN(Cloud_processes%da_ls, 1. - qa_t),0.)*drhcqs_ls 

!-----------------------------------------------------------------------
!     compute qs and condensation using updated temps and vapor. 
!     see Tompkins et al., 2007
!-----------------------------------------------------------------------
      ttmp=Input_mp%tin + ST        
      qtmp= Input_mp%qin + SQ         
      CALL compute_qs_x1 (idim, jdim, kdim, ttmp, Input_mp%pfull, &
                                       qs_t, qs_l, qs_i, dqsdT1, gamma1 )
      DO k=1,kdim
        do j=1,jdim
          DO i=1, idim
            qvmax = qs_t(i,j,k)*(qa_t(i,j,k) + (1. -  qa_t(i,j,k))*  &
                                             Atmos_state%rh_crit(i,j,k) )
            tmp1s =  max(0., (qtmp(i,j,k) - qvmax) )/(1. + gamma1(i,j,k))
            ! limit 
            IF (Cloud_processes%dcond_ls(i,j,k) .GT. 0. .OR.    &
                                                    tmp1s .GT. 0. )  THEN
              Cloud_processes%dcond_ls(i,j,k) =   &
                          MAX(Cloud_processes%dcond_ls(i,j,k), tmp1s) 
            ELSE
            !don't evaporate into saturated grid cells
            !CHECK
              if ( qve(i,j,k) .lt. qs_t(i,j,k)) then
                tmp1s =  max(0.,( qs_t(i,j,k) - qtmp(i,j,k) ))/   &
                               (1. + gamma1(i,j,k)) ! evap .le. qs - qtmp
                Cloud_processes%dcond_ls(i,j,k) = MAX( -1.*tmp1s,  &
                                         Cloud_processes%dcond_ls(i,j,k))
              end if
            END IF
          END DO
        END DO
      END DO

!-------------------------------------------------------------------  


end SUBROUTINE tiedtke_macro_nopdf_super 


!#########################################################################

SUBROUTINE tiedtke_macro_pdf (  &
             idim, jdim, kdim, Atmos_state, Input_mp, Cloud_state, ST,   &
             SQ, Cloud_processes, Particles, n_diag_4d, diag_4d,    &
                                                     diag_id, diag_pt, SA)

!----------------------------------------------------------------------!
!                                                                      !
!                     STATISTICAL CLOUD SCHEME                         !
!                                                                      !
!-----------------------------------------------------------------------

integer,                           intent(in )   :: idim, jdim, kdim     
type(atmos_state_type),            intent(inout) :: Atmos_state
type(mp_input_type),               intent(inout) :: Input_mp    
type(cloud_state_type),            intent(inout) :: Cloud_state
type(cloud_processes_type),        intent(inout) :: Cloud_processes
type(particles_type),              intent(inout) :: Particles
real, dimension(idim,jdim,kdim),   intent(in )   :: ST, SQ
real, dimension(idim,jdim,kdim),   intent(inout) :: SA
integer,                           intent(in )   :: n_diag_4d
real, dimension( idim, jdim, kdim, 0:n_diag_4d ),         &
                                   intent(inout) :: diag_4d
type(diag_id_type),                intent(in)    :: diag_id
type(diag_pt_type),                intent(in)    :: diag_pt

!-----------------------------------------------------------------------
!-----local variables
 
!-----------------------------------------------------------------------
!                STATISTICAL CLOUD SCHEME VARIABLES
!
!
!       qcg            equilibrium value of cloud      kg condensate /
!                      condensate that PDF clouds      kg air
!                      wants           
!
!       qag            equilibrium value of cloud      dimensionless 
!                      fraction for statistical 
!                      cloud scheme           
!
!       qs_norm        the difference between the      dimensionless
!                      saturation specific humidity    
!                      and qtmin normalized by deltaQ
!
!       icbp           the value of the incomplete     dimensionless
!                      beta function evaluated with
!                      x=qs_norm, p=betaP, and q=betaP
!
!       icbp1          the value of the incomplete     dimensionless
!                      beta function evaluated with                  
!                      x=qs_norm, p=betaP+1, and 
!                      q=betaP
!
!       qtbar          total water specific humidity   kg water /
!                      which is equal to the sum of    kg air
!                      liquid water, ice water, and
!                      water vapor
!
!       deltaQ         the width of the total water    kg water /
!                      subgrid distribution (= qtmax   kg air
!                      minus qtmin)
!
!       qtmin          the minimum value to the total  kg water /
!                      sub-grid scale distribution     kg air
!
      real, dimension(idim, jdim,kdim)   :: qa1, qa0, qcg, qag
      real, dimension(idim, jdim,kdim)   :: qta, qtqsa
      real, dimension(idim,jdim)         :: qagtmp, qcgtmp, qvgtmp, qtbar,&
                                            deltaQ, qtmin, qs_norm     
      real                               :: icbp, icbp1
      integer                            :: i,j,k

 
!-----------------------------------------------------------------------
!    Compute total water and excess saturation.
!-----------------------------------------------------------------------
      if (pdf_org) THEN
        qta(:,:,:) = max(qmin, Input_mp%qin +  &
                                  Cloud_state%ql_in + Cloud_state%qi_in)
      else
        qta(:,:,:) = max(qmin, Input_mp%qin + SQ +  &
                             Cloud_state%ql_upd + Cloud_state%qi_upd)
      endif 

      qtqsa(:,:,:) = qta(:,:,:) - Atmos_state%qs
        

#ifdef SKIP
      qcg(:,:,:) = 0.
      Cloud_processes%qvg(:,:,:) = 0.
      qag(:,:,:) = 0.
#endif
       
!-----------------------------------------------------------------------
!    set Tiedtke erosion term to zero.
!-----------------------------------------------------------------------
      Cloud_processes%D_eros = 0.
        
!-----------------------------------------------------------------------
!    compute pdf cloud fraction and condensate
! 
!    Note that the SYMMETRIC beta distribution is used here.
!
!
!    Initialize grid-box mean values of cloud fraction (qag),
!    cloud condensate(qcg), and clear sky water vapor (qvg)
!-----------------------------------------------------------------------
     ks_loop: DO k=1,kdim
       !! qcg(:,k) = 0.
        ! qvg(:,k) = 0.
        !!qag (:,k)= 0.
        
#ifdef SKIP
        !Create loop over sub-levels within a grid box
   
        sublevel_loop: do ns = 1, nsublevels
        
             !calculate normalized vertical level
             ! 0. = top of gridbox
             ! 1. = bottom of gridbox
        
             pnorm =  (real(ns) - 0.5 )/real(nsublevels)
        
             !First step is to calculating the minimum (qtmin)
             !of the total water distribution and 
             !the width of the qt distribution (deltaQ)
             !
             !For diagnostic variance this is set to (1.-qthalfwidth)*qtbar
             !and 2*qthalfwidth*qtbar, respectively, where qtbar is the
             !mean total water in the grid box.        
             !
             !

             qtbar = qta4(2,:,:,k) + pnorm*((qta4(3,:,:,k) -   &
                                qta4(2,:,:,k)) + qta4(4,:,:,k)*(1-pnorm) )
             
             qtbar  = max(qmin ,qtbar )
             deltaQ  = 2.*qthalfwidth *qtbar
             qtmin = (1.- qthalfwidth )*qtbar 
#endif
             qtbar = qta(:,:,k)
             qtbar  = max(qmin ,qtbar )
             deltaQ  = 2.*qthalfwidth *qtbar
             qtmin = (1.-qthalfwidth )*qtbar
        
             !From this the variable normalized saturation specific
             !humidity qs_norm is calculated.
             !
             !  qs_norm = (qs(Tl) - qtmin)/(qtmax-qtmin)
             !
             !          = 0.5  - (qtbar - qs(Tl))/deltaQ
             !
             !Note that if qs_norm > 1., the grid box is fully clear.
             !If qs_norm < 0., the grid box is fully cloudy.
        
!            qs_norm = qtqsa4(2,:,:,k)+  &
!                      pnorm*( (qtqsa4(3,:,:,k)-qtqsa4(2,:,:,k)) + &
!                      qtqsa4(4,:,:,k)*(1-pnorm) )

             qs_norm = qtqsa(:,:,k)
             qs_norm = 0.5 - ( qs_norm/deltaQ )
             
             !Calculation of cloud fraction (qagtmp), cloud condensate 
             !(qcgtmp), and water vapor in clear air part of the grid 
             !box (qvgtmp)
             !
             !Formulas (from Tompkins, and personal derivations):
             !
             !  Define icbp  = incomplete_beta(qs_norm,p,q)
             !         icbp1 = incomplete_beta(qs_norm,p+1,q)
             !
             !  qagtmp = 1. - icbp
             !
             !  qcgtmp = aThermo * {  (qtbar-qtmin)*(1.-icbp1) - 
             !                       qs_norm*deltaQ*(1.-icbp ) }
             !
             !
             !  qvgtmp = qtmin + (p/(p+q))*(icbp1/icbp)*deltaQ
             !
             !  
             ! where aThermo = 1./(1.+(L/cp)*dqsdT)
             !
             ! note that in the qvg formula below the factor of 0.5
             ! is equal to (p/(p+q)).
             !

             do j=1,jdim
             
             do i = 1,idim
        
             if (qs_norm(i,j).le.1.) then
                 
                 icbp = incomplete_beta(max(0.,qs_norm(i,j)), &
                                      p = betaP    , q =     betaP)
                 icbp1= incomplete_beta(max(0.,qs_norm(i,j)), &
                                      p = betaP + 1, q =     betaP)
                 qagtmp(i,j) = 1.-icbp
                 qcgtmp(i,j) = (qtbar(i,j)-qtmin(i,j))*(1.-icbp1)&
                               - qs_norm(i,j)*deltaQ(i,j)*(1.-icbp)    
                 qcgtmp(i,j) = qcgtmp(i,j)/   &
                                          (1.+Atmos_state%gamma(i,j,k))
                 qvgtmp(i,j) = qtmin(i,j) + &
                               0.5*(icbp1/max(icbp, qmin))*deltaQ(i,j)
             
                 !bound very very small cloud fractions which may
                 !cause negative cloud condensates due to roundoff 
                 !errors or similar errors in the beta table lookup.
                 if((qagtmp(i,j).lt.0.).or.(qcgtmp(i,j).le.0.))then
                      qagtmp(i,j) = 0.
                      qcgtmp(i,j) = 0.
                      qvgtmp(i,j) = qtbar(i,j)
                 end if
                 
             else             
                 qagtmp(i,j) = 0.
                 qcgtmp(i,j) = 0.
                 qvgtmp(i,j) = qtbar(i,j)             
             end if
             
             enddo
             enddo
              
#ifdef SKIP
             !sum vertically
             !
             !note special averaging of clear-sky water vapor
             !this is weighting clear-sky relative humidity by the 
             !clear-sky fraction
         
             qag(:,:,k) = qag(:,:,k) + qagtmp
             qcg(:,:,k) = qcg(:,:,k) + qcgtmp
             Cloud_processes%qvg(:,:,k) =    &
                     Cloud_processes%qvg(:,:,k)+ &
                               (1.-qagtmp)*min(max(qvgtmp/max(qmin, &
                                   (qtbar+((qs_norm-0.5)*deltaQ))),0.),1.)
             
        enddo sublevel_loop!for number of sublevels loop

        


        !compute grid-box average cloud fraction, cloud condensate
        !and water vapor
        
        if (nsublevels.gt.1) then
             qag (:,:,k)= qag(:,:,k) / real(nsublevels)
             qcg(:,:,k) = qcg(:,:,k) / real(nsublevels)
             
             !note special averaging of clear-sky water vapor
            
             do j = 1,jdim
             do id = 1,idim
                  if ((1.-qag(id,j,k)).gt. qmin) then
                    Cloud_processes%qvg(id,j,k) =  &
                       Cloud_processes%qvg(id,j,k)/real(nsublevels)/&
                                                        (1.-qag(id,j,k))
                    Cloud_processes%qvg(id,j,k) =  &
                              Cloud_processes%qvg(id,j,k)*   &
                                                 Atmos_state%qs(id,j,k)
                  else
                    Cloud_processes%qvg(id,j,k) =   &
                                                Atmos_state%qs(id,j,k)
                  end if
             enddo
             enddo
          
             
!       else
#endif
             qag(:,:,k) = qagtmp
             qcg(:,:,k) = qcgtmp
             Cloud_processes%qvg(:,:,k) = qvgtmp
!        end if
     
       
   END DO ks_loop
        
        !do adjustment of cloud fraction
        qa0 = Cloud_state%qa_in
        qa1 = qag

        !set total tendency term and update cloud fraction    
        SA   = (SA   + qa1) - qa0
        Cloud_state%qa_upd     = qa1


        !define da_ls and tmp5 needed when do_liq_num = .true. (cjg)
        Cloud_processes%da_ls = max(qa1-qa0,0.)
        Cloud_processes%delta_cf = max(qa1-qa0,0.)

        !compute large-scale condensation / evaporation
        Cloud_processes%dcond_ls = qcg -    &
                              (Cloud_state%ql_upd + Cloud_state%qi_upd)



   IF ( .NOT. pdf_org ) THEN
   !!!! INVESTIGATE!!!
   !make sure super/subsat is not created here
   ! this is different from the original PDF assumption
   !  (as is saturation adjustment) 

     Cloud_processes%dcond_ls = MAX( ((Input_mp%qin + SQ -  &
                Atmos_State%qs)/(1.+Atmos_state%gamma)),    &
                                               Cloud_processes%dcond_ls)
       
      
   endif 

!-----------------------------------------------------------------------
!     Diagnostics
!-----------------------------------------------------------------------

      if (diag_id%qadt_lsform + diag_id%qa_lsform_col > 0) then
          diag_4d(:,:,:,diag_pt%qadt_lsform ) =  max(qa1 - qa0, 0.)*  &
                                                               inv_dtcloud
      end if
      if (diag_id%qadt_lsdiss + diag_id%qa_lsdiss_col > 0) then
          diag_4d(:,:,:,diag_pt%qadt_lsdiss ) =  max(qa0 - qa1, 0.)* &
                                                               inv_dtcloud
      end if

!------------------------------------------------------------------------
 


end SUBROUTINE tiedtke_macro_pdf 


!#########################################################################

!SUBROUTINE nc_cond_Single_Gaussian_pdf (idim, jdim, kdim, Nml, Constants, Atmos_state,   &
SUBROUTINE tiedtke_macro_Single_Gaussian_pdf (idim, jdim, kdim,  Atmos_state,   &
                        Input_mp, &
                        Cloud_state, ST, SQ, Cloud_processes, Particles, &
                        n_diag_4d, diag_4d, diag_id, diag_pt,       SA)  

!----------------------------------------------------------------------!
!                                                                      !
!                     STATISTICAL CLOUD SCHEME                         !
!                                                                      !
!-----------------------------------------------------------------------

integer,                           intent(in)    :: idim, jdim, kdim     
!type(strat_nml_type),              intent(in)    :: Nml
type(atmos_state_type),            intent(inout) :: Atmos_state
type(mp_input_type),               intent(inout) :: Input_mp
type(cloud_state_type),            intent(inout) :: Cloud_state
type(cloud_processes_type),        intent(inout) :: Cloud_processes
type(particles_type),              intent(inout) :: Particles
!type(strat_constants_type),        intent(inout) :: Constants  
real, dimension(idim,jdim,kdim),   intent(in)    :: ST, SQ
real, dimension(idim,jdim,kdim),   intent(inout) :: SA
integer,                           intent(in)    :: n_diag_4d
real, dimension(idim,jdim,kdim,0:n_diag_4d),         &
                                   intent(inout) :: diag_4d
type(diag_id_type),                intent(in)    :: diag_id
type(diag_pt_type),                intent(in)    :: diag_pt
!integer,                           intent(in)    :: otun

!-----------------------------------------------------------------------
!-----local variables---------------------------------------------------
 
!-----------------------------------------------------------------------
!                STATISTICAL CLOUD SCHEME VARIABLES
!
!
!       qcg            total condensate diagnosed        kg condensate /
!                      from Single-Gaussian PDF          kg air        
!
!       qag            cloud fraction diagnosed          dimensionless 
!                      from Single-Gaussian PDF                                  
!
!
!       qtbar          total water specific humidity     kg water /
!                      which is equal to the sum of      kg air
!                      liquid water, ice water, and
!                      water vapor
!
!       s              s = (qt -qs)/(1 + L/Cp * dqsdT)   kg condensate/kg air
!                     
!       stddev_s       standard deviation of s           kg condensate/kg air
!
!       if ignoring liquid potential temperature variance,  
!       and correlation between liquid potential temperature
!       and total water, 
!
!       stddev_s = stddev_qt /(1 + L/Cp * dqsdT)
!                = stddev_qt /(1 + Atmos_state%gamma ) 


      real, dimension(idim,jdim,kdim)    :: qa1, qa0, qcg, qag
      real, dimension(idim,jdim,kdim)    :: qta, qtqsa
 
!---> h1g, 2015-07-22
     real, dimension(idim,jdim,kdim)    ::  s,   stdev_s
     real                               ::  s_mellor_tol, zeta
     integer                            ::  i, j, k
!-----------------------------------------------------------------------
!     Compute total water and excess saturation
!-----------------------------------------------------------------------
      if (pdf_org) then
!RSH  change to Input_mp%qin ??
!       qta(:,:,:) = max(qmin, Atmos_state%qv_in +  &
        qta(:,:,:) = max(qmin, Input_mp%qin +  &
                                     Cloud_state%ql_in + Cloud_state%qi_in)
      else
!       qta(:,:,:) = max(qmin, Atmos_State%qv_in + SQ +  &
        qta(:,:,:) = max(qmin, Input_mp%qin + SQ +  &
                                     Cloud_state%ql_upd + Cloud_state%qi_upd)
      endif
      qtqsa(:,:,:) = qta(:,:,:) - Atmos_state%qs      

      s =  qtqsa(:,:,:) / (1.0 +  Atmos_state%gamma )

!---> h1g, 2015-07-23
! 99.7% data are within mean +- 3 stdev 
      stdev_s =     qthalfwidth * qta * 0.333
      stdev_s = stdev_s / ( 1.0 +  Atmos_state%gamma ) 

      s_mellor_tol = 1.e-8
      do k= 1,kdim
        do  j = 1,jdim
          do i = 1,idim
            if( stdev_s(i,j,k)  > s_mellor_tol ) then
              zeta = s(i,j,k)/stdev_s(i,j,k)
              qag(i,j,k)  = 0.5 * ( 1.0 + erf( zeta /sqrt(2.0) ) )
              qcg(i,j,k)  = s(i,j,k) * qag(i,j,k)  &
                          + stdev_s(i,j,k)/sqrt(2.0*pi) * exp( -0.5 *zeta*zeta )
            else
              if ( s(i,j,k) < 0.0 ) then
                qag(i,j,k)  = 0.0
                qcg(i,j,k)  = 0.0
              else
                qag(i,j,k)  = 1.0
                qcg(i,j,k)  = s(i,j,k)
              endif  
            endif

!---> h1g, 2015-07-23
!     qtbar:     grid-average total specific humidity
!     qv,clr:    clear-sky  water vapor specific humidity
!     qc,cld:    in-cloud  cloud condensate
!     qs:        satuation vapor specific humidity
!
!     qtbar = (1.-CF) * qv,clr + CF * ( qc,cld + qs )
!     qv,clr = (qtbar - CF * ( qc,cld + qs ) )/(1.-CF)      
!<--- h1g, 2015-07-23
            Cloud_processes%qvg(i,j,k) = (qta(i,j,k) - qcg(i,j,k) - qag(i,j,k)*Atmos_state%qs(i,j,k))  &
                                        /(max( (1.0 - qag(i,j,k)), qmin) )
          enddo
        enddo
      enddo        

      !do adjustment of cloud fraction
      qa0 = Cloud_state%qa_in
      qa1 = qag

      !set total tendency term and update cloud fraction    
      SA   = (SA   + qa1) - qa0
      Cloud_state%qa_upd     = qa1

      !define da_ls and tmp5 needed when do_liq_num = .true. (cjg)
      Cloud_processes%da_ls = max(qa1-qa0,0.)
      Cloud_processes%delta_cf = max(qa1-qa0,0.)

      !compute large-scale condensation / evaporation
      Cloud_processes%dcond_ls = qcg -    &
                              (Cloud_state%ql_upd + Cloud_state%qi_upd)


!--->h1g, the following is inherited from nc_cond_pdf, 2015-07-22
      if ( .not. pdf_org ) then
      !!!! INVESTIGATE!!!
      ! make sure super/subsat is not created here
      ! this is different from the original PDF assumption
      ! (as is saturation adjustment) 

!RSH    Cloud_processes%dcond_ls = MAX( ((Atmos_state%qv_in + SQ -  &
        Cloud_processes%dcond_ls = MAX( ((Input_mp%qin      + SQ -  &
                Atmos_State%qs)/(1.+Atmos_state%gamma)),    &
                                               Cloud_processes%dcond_ls)
  
      end if
!<--- h1g, 2015-07-22

!-----------------------------------------------------------------------
!     Diagnostics
!-----------------------------------------------------------------------

      if (diag_id%qadt_lsform + diag_id%qa_lsform_col > 0) then
          diag_4d(:,:,:,diag_pt%qadt_lsform ) =  max(qa1 - qa0, 0.)*  &
                                                     inv_dtcloud 
      end if
      if (diag_id%qadt_lsdiss + diag_id%qa_lsdiss_col > 0) then
          diag_4d(:,:,:,diag_pt%qadt_lsdiss ) =  max(qa0 - qa1, 0.)* &
                                                     inv_dtcloud
      end if
!------------------------------------------------------------------------

end SUBROUTINE tiedtke_macro_Single_Gaussian_pdf 

!#####################################################################
 
  function erf(x)
!--------------------------------------------------------------
! This numerical recipes routine calculates the error function.
!--------------------------------------------------------------
 
   real :: erf
   real, intent(in) :: x
   real :: t,z
 
   z=abs(x)
   t=1./(1.+0.5*z)
 
    erf= 1.0 - t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t*      &
         (.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t*    &
         (1.48851587+t*(-.82215223+t*.17087277)))))))))
 
    if (x.lt.0.) erf= -erf
 
   end function erf


!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  ppm2m_sak --- Piecewise parabolic method for fields
!
! !INTERFACE:
 subroutine ppm2m_sak(a4, delp, km, kmap, i1, i2, iv, kord)

!implicit none

! !INPUT PARAMETERS:
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: i1      ! Starting longitude
 integer, intent(in):: i2      ! Finishing longitude
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: kmap    ! partial remap to start
 integer, intent(in):: kord    ! Order (or more accurately method no.):
                               ! 
 real, intent(in):: delp(i1:i2,km)     ! layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real, intent(inout):: a4(4,i1:i2,km)  ! Interpolated values

! !DESCRIPTION:
!
!   Perform the piecewise parabolic method 
! 
! !REVISION HISTORY: 
!   ??.??.??    Lin        Creation
!   02.04.04    Sawyer     Newest release from FVGCM
! 
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! local arrays:
      real   dc(i1:i2,km)
      real   h2(i1:i2,km)
      real delq(i1:i2,km)
      real  df2(i1:i2,km)
      real   d4(i1:i2,km)

! local scalars:
      integer i, k, km1, lmt
      integer it
      real fac
      real a1, a2, c1, c2, c3, d1, d2
      real qmax, qmin, cmax, cmin
      real qm, dq, tmp
      real qmp, pmp
      real lac

      km1 = km - 1
       it = i2 - i1 + 1

      do k=max(2,kmap-2),km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo
 
      do k=max(2,kmap-2),km1
         do i=i1,i2
            c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1)
            c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k)
            tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            qmax = max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1)) - a4(1,i,k)
            qmin = a4(1,i,k) - min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))
             dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp)
            df2(i,k) = tmp
         enddo
      enddo

!-----------------------------------------------------------
! 4th order interpolation of the provisional cell edge value
!-----------------------------------------------------------

      do k=max(3,kmap), km1
      do i=i1,i2
        c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
        a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
        a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
        a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *    &
                  ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                                delp(i,k-1)*a1*dc(i,k  ) )
      enddo
      enddo

      if(km>8 .and. kord>3) call steepz_sak(i1, i2, km, kmap, a4, df2, dc, delq, delp, d4)

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      if ( kmap <= 2 ) then
      do i=i1,i2
         d1 = delp(i,1)
         d2 = delp(i,2)
         qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
         dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
         c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
         c3 = dq - 0.5*c1*(d2* (5.*d1+d2)-3.*d1**2)
         a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
         a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
         dc(i,1) =  a4(1,i,1) - a4(2,i,1)
! No over- and undershoot condition
         cmax = max(a4(1,i,1), a4(1,i,2))
         cmin = min(a4(1,i,1), a4(1,i,2))
         a4(2,i,2) = max(cmin,a4(2,i,2))
         a4(2,i,2) = min(cmax,a4(2,i,2))
      enddo
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1**2)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
         dc(i,km) = a4(3,i,km) -  a4(1,i,km)
! No over- and under-shoot condition
         cmax = max(a4(1,i,km), a4(1,i,km1))
         cmin = min(a4(1,i,km), a4(1,i,km1))
         a4(2,i,km) = max(cmin,a4(2,i,km))
         a4(2,i,km) = min(cmax,a4(2,i,km))
      enddo

      do k=max(1,kmap),km1
         do i=i1,i2
            a4(3,i,k) = a4(2,i,k+1)
         enddo
      enddo

! Enforce monotonicity of the "slope" within the top layer
      if ( kmap <= 2 ) then
      do i=i1,i2
         if ( a4(2,i,1) * a4(1,i,1) <= 0. ) then 
              a4(2,i,1) = 0.
                dc(i,1) = a4(1,i,1)
         endif
         if ( dc(i,1) * (a4(2,i,2) - a4(1,i,1)) <= 0. ) then
! Setting DC==0 will force piecewise constant distribution after
! calling kmppm_sak
              dc(i,1) = 0.
         endif
      enddo
      endif

! Enforce constraint on the "slope" at the surface

      do i=i1,i2
         if( a4(3,i,km) * a4(1,i,km) <= 0. ) then
!            a4(3,i,km) = 0.
!              dc(i,km) =  -a4(1,i,km)
               dc(i,km) = 0.
         endif
         if( dc(i,km) * (a4(1,i,km) - a4(2,i,km)) <= 0. ) then
             dc(i,km) = 0.
         endif
      enddo
 
!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
      if ( kmap <= 2 ) then
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
            call kmppm_sak(dc(i1,k), a4(1,i1,k), it, 0)
      enddo
      endif

      if(kord >= 7) then
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      do k=max(2,kmap-1), km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2
!           h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))
!    &               / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )
!    &               * delp(i,k)**2
! Method#3
            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      if( kord == 7 ) then
         fac = 1.5           ! original quasi-monotone
      else
         fac = 0.125         ! full monotone
      endif

      do k=max(3,kmap), km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(3,i,k) = min(max(a4(3,i,k), qmin), qmax)
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(2,i,k) = min(max(a4(2,i,k), qmin), qmax)
! Recompute A6
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to ensure positivity when kord=7
         if (iv == 0 .and. kord == 7) then
             call kmppm_sak(dc(i1,k), a4(1,i1,k), it, 2)
         endif
      enddo

      else
 
         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv == 0) lmt = min(2, lmt)

      do k=max(3,kmap), km-2
      if( kord /= 4) then
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
      endif
         call kmppm_sak(dc(i1,k), a4(1,i1,k), it, lmt)
      enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm_sak(dc(i1,k), a4(1,i1,k), it, 0)
      enddo
!EOC
 end subroutine ppm2m_sak
!-----------------------------------------------------------------------


!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  kmppm_sak --- Perform piecewise parabolic method in vertical
!
! !INTERFACE:
 subroutine kmppm_sak(dm, a4, itot, lmt)

!implicit none

! !INPUT PARAMETERS:
      real, intent(in):: dm(*)     ! the linear slope
      integer, intent(in) :: itot      ! Total Longitudes
      integer, intent(in) :: lmt       ! 0: Standard PPM constraint
                                       ! 1: Improved full monotonicity constraint (Lin)
                                       ! 2: Positive definite constraint
                                       ! 3: do nothing (return immediately)
! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout) :: a4(4,*)   ! PPM array
                                           ! AA <-- a4(1,i)
                                           ! AL <-- a4(2,i)
                                           ! AR <-- a4(3,i)
                                           ! A6 <-- a4(4,i)

! !DESCRIPTION:
!
! !REVISION HISTORY: 
!    00.04.24   Lin       Last modification
!    01.03.26   Sawyer    Added ProTeX documentation
!    02.04.04   Sawyer    Incorporated newest FVGCM version
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      real, parameter:: r12 = 1./12.
      real qmp
      real da1, da2, a6da
      real fmin
      integer i

! Developer: S.-J. Lin, NASA-GSFC
! Last modified: Apr 24, 2000

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      do i=1,itot
      if(dm(i) == 0.) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt == 1) then

! Improved full monotonicity constraint (Lin 2003)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, itot
           qmp = 2.*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = 3.*( 2.*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt == 2) then

! Positive definite constraint
      do i=1,itot
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
         if( fmin < 0. ) then
         if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
         elseif(a4(3,i) > a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

!EOC
 end subroutine kmppm_sak
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  steepz_sak --- Calculate attributes for PPM
!
! !INTERFACE:
 subroutine steepz_sak(i1, i2, km, kmap, a4, df2, dm, dq, dp, d4)

!  implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: km                   ! Total levels
      integer, intent(in) :: kmap                 ! 
      integer, intent(in) :: i1                   ! Starting longitude
      integer, intent(in) :: i2                   ! Finishing longitude
      real, intent(in) ::  dp(i1:i2,km)       ! grid size
      real, intent(in) ::  dq(i1:i2,km)       ! backward diff of q
      real, intent(in) ::  d4(i1:i2,km)       ! backward sum:  dp(k)+ dp(k-1) 
      real, intent(in) :: df2(i1:i2,km)       ! first guess mismatch
      real, intent(in) ::  dm(i1:i2,km)       ! monotonic mismatch

! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout) ::  a4(4,i1:i2,km)  ! first guess/steepened

!
! !DESCRIPTION:
!   This is complicated stuff related to the Piecewise Parabolic Method
!   and I need to read the Collela/Woodward paper before documenting
!   thoroughly.
!
! !REVISION HISTORY: 
!   ??.??.??    Lin?       Creation
!   01.03.26    Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, k
      real alfa(i1:i2,km)
      real    f(i1:i2,km)
      real  rat(i1:i2,km)
      real  dg2

! Compute ratio of dq/dp
      do k=max(2,kmap-1),km
         do i=i1,i2
            rat(i,k) = dq(i,k-1) / d4(i,k)
         enddo
      enddo

! Compute F
      do k=max(2,kmap-1),km-1
         do i=i1,i2
            f(i,k) =   (rat(i,k+1) - rat(i,k))                          &
                     / ( dp(i,k-1)+dp(i,k)+dp(i,k+1) )
         enddo
      enddo

      do k=max(3,kmap),km-2
         do i=i1,i2
         if(f(i,k+1)*f(i,k-1)<0. .and. df2(i,k)/=0.) then
            dg2 = (f(i,k+1)-f(i,k-1))*((dp(i,k+1)-dp(i,k-1))**2          &
                   + d4(i,k)*d4(i,k+1) )
            alfa(i,k) = max(0., min(0.5, -0.1875*dg2/df2(i,k))) 
         else
            alfa(i,k) = 0.
         endif
         enddo
      enddo

      do k=max(4,kmap+1),km-2
         do i=i1,i2
            a4(2,i,k) = (1.-alfa(i,k-1)-alfa(i,k)) * a4(2,i,k) +         &
                        alfa(i,k-1)*(a4(1,i,k)-dm(i,k))    +             &
                        alfa(i,k)*(a4(1,i,k-1)+dm(i,k-1))
         enddo
      enddo

!EOC
 end subroutine steepz_sak
!----------------------------------------------------------------------- 




!#######################################################################


                end module tiedtke_macro_mod
