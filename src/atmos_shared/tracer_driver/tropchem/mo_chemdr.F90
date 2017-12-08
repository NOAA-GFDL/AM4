
      module mo_chemdr_mod

!<f1p
      use fms_mod,            only : FATAL, error_mesg, mpp_pe, mpp_root_pe
      use tropchem_types_mod, only : tropchem_opt,tropchem_diag, small_value
      use mo_chem_utls_mod,   only : get_spc_ndx      
      use tracer_manager_mod, only : get_tracer_index        
      use field_manager_mod,  only : MODEL_ATMOS
      use mpp_mod,            only : mpp_clock_id,         &
                                     mpp_clock_begin,      &
                                     mpp_clock_end
!>
      implicit none

 !<f1p
      integer :: nh4no3_ndx, hno3_ndx
      integer :: nh4_ndx, nh3_ndx
      integer :: so2_ndx, so4_ndx
      integer :: ox_ndx, o3s_ndx
      integer :: o3s_e90_ndx, e90_ndx
      integer :: implicit_clock_id, photo_clock_id, cloud_clock_id
!>
      private
      public :: chemdr
      public :: chemdr_init !<f1p>

character(len=128), parameter :: version     = '$Id$'
character(len=128), parameter :: tagname     = '$Name$'
logical                       :: module_is_initialized = .false.

      contains

! <SUBROUTINE NAME="chemdr">
!   <OVERVIEW>
!     Tropospheric chemistry driver
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine calculates chemical production and loss of trace species,
!     using the MOZART chemical mechanism and solver.
!     Species and reactions can be modified using MOZART chemical pre-processor.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call chemdr( vmr, Time, lat, lon, delt, ps, ptop, pmid, pdel, &
!                  zma, zi, cldfr, cwat, tfld, inv_data, sh, &
!                  albedo, coszen, esfact, &
!                  prod_out, loss_out, sulfate, psc, &
!                  do_interactive_h2o, solar_phase, imp_slv_nonconv, plonl,&
!                  prod_ox, loss_ox )
!   </TEMPLATE>
!   <IN NAME="Time" TYPE="time_type">
!     Model time
!   </INOUT>
!   <IN NAME="lat" TYPE="real" DIM="(:)">
!     The latitudes for the local domain.
!   </IN>
!   <IN NAME="lon" TYPE="real" DIM="(:)">
!     The longitudes for the local domain.
!   </IN>
!   <IN NAME="delt" TYPE="real">
!     Model timestep (s)
!   </IN>
!   <IN NAME="pmid" TYPE="real" DIM="(:,:)">
!     Pressure on the model full levels (Pa)
!   </IN>
!   <IN NAME="pdel" TYPE="real" DIM="(:,:)">
!     Pressure thickness of model layers (Pa)
!   </IN>
!   <IN NAME="ps" TYPE="real" DIM="(:)">
!     Model surface pressure (Pa)
!   </IN>
!   <IN NAME="ptop" TYPE="real" DIM="(:)">
!     Pressure at model top (Pa)
!   </IN>
!   <IN NAME="zma" TYPE="real" DIM="(:,:)">
!     Model full level absolute geopotential heights (m)
!   </IN>
!   <IN NAME="zi" TYPE="real" DIM="(:,:)">
!     Model half level absolute geopotential heights (m)
!   </IN>
!   <IN NAME="coszen" TYPE="real" DIM="(:)">
!     Cosine of solar zenith angle
!   </IN>
!   <IN NAME="albedo" TYPE="real" DIM="(:)">
!     Surface albedo
!   </IN>
!   <IN NAME="cldfr" TYPE="real" DIM="(:,:)">
!     Cloud fraction
!   </IN>
!   <IN NAME="cwat" TYPE="real" DIM="(:,:)">
!     Cloud liquid+ice water (kg/kg)
!   </IN>
!   <IN NAME="cwat" TYPE="real" DIM="(:,:)">
!     Cloud liquid+ice water (kg/kg)
!   </IN>
!   <IN NAME="tfld" TYPE="real" DIM="(:,:)">
!     Model temperature (K)
!   </IN>
!   <IN NAME="sh" TYPE="real" DIM="(:,:)">
!     Model specific humidity (kg/kg)
!   </IN>
!   <IN NAME="sulfate" TYPE="real" DIM="(:,:)">
!     Off-line sulfate aerosol volume mixing ratio (mol/mol)
!   </IN>
!   <IN NAME="psc" TYPE="psc_type">
!     Polar stratospheric cloud amounts
!   </IN>
!   <IN NAME="solar_phase" TYPE="real">
!     Solar cycle phase (1=max, 0=min)
!   </IN>
!   <IN NAME="esfact" TYPE="real">
!     Earth-sun distance factor (r_avg/r)^2
!   </IN>
!   <IN NAME="inv_data" TYPE="real" DIM="(:,:,:)">
!     Volume mixing ratios of "invariant" species (mol/mol)
!   </IN>
!   <IN NAME="do_interactive_h2o" TYPE="logical">
!     Include water vapor sources/sinks?
!   </IN>
!   <INOUT NAME="vmr" TYPE="real" DIM="(:,:,:)">
!     Trace species volume mixing ratio (mol/mol)
!   </INOUT>
!   <OUT NAME="prod_out" TYPE="real" DIM="(:,:,:)">
!     Trace species photochemical production rates (mol/mol/s)
!   </OUT>
!   <OUT NAME="loss_out" TYPE="real" DIM="(:,:,:)">
!     Trace species photochemical loss rates (mol/mol/s)
!   </OUT>
!   <OUT NAME="imp_slv_nonconv" TYPE="real" DIM="(:,:)">
!     Flag for implicit solver non-convergence (fraction)
!   </OUT>
!   <OUT NAME="prod_ox" TYPE="real" DIM="(:,:)">
!     output for ox production
!   </OUT>
!   <OUT NAME="loss_ox" TYPE="real" DIM="(:,:)">
!     output for ox loss
!   </OUT>
! prod_ox and loss_ox are added to chemdr. (jmao, 1/1/2011)
! r,phalf, pwt and j_ndx are added for fastjx.(jmao,1/1/2011)
      subroutine chemdr( vmr, &
                         r, &
                         phalf,&
                         pwt , &
                         j_ndx, &
                         Time, &
                         lat, lon, &
                         delt, &
                         ps, ptop, pmid, pdel, &
                         zma, zi, &
                         cldfr, cwat, fliq, tfld, inv_data, sh, &
                         albedo, coszen, esfact, &
                         prod_out, loss_out, jvals_out, rate_const_out, sulfate, psc, &
                         do_interactive_h2o, solar_phase, imp_slv_nonconv, &
                         plonl, prod_ox, loss_ox, e90_vmr, &
                         e90_tropopause_vmr, &
                         co2, &
                         trop_diag_array, trop_option,trop_diag) !>
!-----------------------------------------------------------------------
!     ... Chem_solver advances the volumetric mixing ratio
!         forward one time step via a combination of explicit,
!         ebi, hov, fully implicit, and/or rodas algorithms.
!-----------------------------------------------------------------------

#ifndef AM3_CHEM
      use chem_mods_mod,    only : indexm, nadv_mass, phtcnt, gascnt, rxntot, clscnt1, clscnt4, clscnt5, &
                                   ncol_abs, grpcnt, nfs, extcnt, hetcnt
#else
      use AM3_chem_mods_mod,only : indexm, nadv_mass, phtcnt, gascnt, rxntot, clscnt1, clscnt4, clscnt5, &
                                   ncol_abs, grpcnt, nfs, extcnt, hetcnt
#endif
      use mo_photo_mod,     only : set_ub_col, setcol, photo, &
                                   sundis
      use mo_exp_sol_mod,   only : exp_sol
      use mo_imp_sol_mod,   only : imp_sol
      use mo_rodas_sol_mod, only : rodas_sol
      use mo_usrrxt_mod,    only : usrrxt
      use mo_setinv_mod,    only : setinv
#ifndef AM3_CHEM
      use mo_setrxt_mod,    only : setrxt
      use mo_adjrxt_mod,    only : adjrxt
      use mo_phtadj_mod,    only : phtadj
#else
      use AM3_mo_setrxt_mod,    only : setrxt
      use AM3_mo_adjrxt_mod,    only : adjrxt
      use AM3_mo_phtadj_mod,    only : phtadj
#endif
      use mo_setsox_mod,    only : setsox
#ifndef AM3_CHEM
      use mo_fphoto_mod,    only : fphoto
#else
      use AM3_fphoto_mod,    only : fphoto
#endif
      use mo_chem_utls_mod, only : inti_mr_xform, adjh2o, negtrc, mmr2vmr, vmr2mmr !<f1p get_spc_ndx, get_grp_mem_ndx >
      use time_manager_mod, only : time_type
      use strat_chem_utilities_mod, only : psc_type


      implicit none

!-----------------------------------------------------------------------
!        ... Dummy arguments
!-----------------------------------------------------------------------
      type(time_type), intent(in) :: Time             ! time
      real,    intent(in) ::  lat(:), lon(:)          ! latitude, longitude
      integer, intent(in) ::  plonl
      real,    intent(in) ::  delt                    ! timestep in seconds
      real, intent(inout) ::  vmr(:,:,:)              ! transported species ( vmr )
      integer, intent(in) ::  j_ndx             
      real, intent(in)    ::  r(:,:,:)                ! the original r species
      real, intent(in)    ::  co2(:,:)
      real, intent(inout) ::  trop_diag_array(:,:,:)
!>
      real, dimension(:), intent(in) :: &
                              ps, &                   ! surface press ( pascals )
                              ptop, &                 ! model top pressure (pascals)
!                             oro, &                  ! surface orography flag
                              albedo, &               ! surface albedo
                              coszen                  ! cosine of solar zenith angle
!                             tsurf, &                ! surface temperature
!                             phis, &                 ! surf geopot
!                             cldtop                  ! cloud top level ( 1 ... plev )
      real, dimension(:,:), intent(in) :: &
                              pmid, &                 ! midpoint press ( pascals )
                              pdel, &                 ! delta press across midpoints
                              zma, &                  ! abs geopot height at midpoints ( m )
                              cldfr, &                ! cloud fraction
!                             cmfdqr, &               ! dq/dt for convective rainout
!                             nrain, &                ! release of strt precip ( 1/s )
!                             nevapr, &               ! evap precip ( 1/s )
                              cwat, &                 ! total cloud water (kg/kg)
                              fliq, &                 !liquid fraction
	                          tfld, &                 ! midpoint temperature
                              sh, &                   ! specific humidity ( kg/kg )
                              sulfate, &              ! sulfate aerosol
                              phalf,   &              ! pressure at boundaries (Pa)
                              pwt
      type(psc_type), intent(in) :: &
                              psc                     ! polar stratospheric clouds (PSCs)
      real, intent(in) ::     solar_phase, &          ! solar cycle phase (1=max, 0=min)
                              esfact                  ! earth-sun distance factor (r_avg/r)^2
      real, dimension(:,:), intent(in) :: &
                              zi                      ! abs geopot height at interfaces ( m )
      real, dimension(:,:,:), intent(out) :: &
                              prod_out, &             ! chemical production rate
                              loss_out                ! chemical loss rate
      real, dimension(:,:),   intent(out) :: &        
                              prod_ox, &              ! chemical production rate, defined as (nlon,nlev)(jmao,1/1/2011)
                              loss_ox                 ! chemical loss rate(jmao,1/1/2011)      
      real, dimension(:,:,:), intent(out) :: &
                              jvals_out               ! photolysis rates (J-values, s^-1)
      real, dimension(:,:,:), intent(out) :: &
                              rate_const_out          ! kinetic rxn rate constants (cm^3 molec^-1 s^-1 for 2nd order)
      real, dimension(:,:,:), intent(in) :: &
                              inv_data                ! invariant species
      real, dimension(:,:), intent(out) :: &
                              imp_slv_nonconv         ! flag for implicit solver non-convergence (fraction)
      logical, intent(in) ::  do_interactive_h2o      ! include h2o sources/sinks
      real, dimension(:,:), intent(in) :: &
                              e90_vmr                ! e90 conc 
      real,    intent(in) ::  e90_tropopause_vmr     ! e90 tropopause threshold 

!<f1p
      type(tropchem_opt),  intent(in) :: trop_option
      type(tropchem_diag), intent(in) :: trop_diag
!>

!-----------------------------------------------------------------------
!             ... Local variables
!-----------------------------------------------------------------------
      integer, parameter :: inst = 1, avrg = 2
      integer  :: k,i,n
      integer  :: troplev(plonl)
      real, parameter    :: ztrop_low = 5.  
                                     ! lowest tropopause level allowed (km)
      real, parameter    :: ztrop_high = 20.  
                                    ! highest tropopause level allowed (km)
      real, parameter    :: max_dtdz   = 2.  
                             ! max dt/dz for tropopause level (degrees k/km)
      real     :: dt
      real, dimension(plonl,SIZE(vmr,2))    :: k_loss_ox   
                                       ! Loss rate coefficient of Ox (s-1)
      integer,dimension(plonl,SIZE(vmr,2)) :: strate90 
      real     ::  invariants(plonl,SIZE(vmr,2),max(1,nfs))
      real     ::  col_dens(plonl,SIZE(vmr,2),max(1,ncol_abs))                  ! column densities (molecules/cm^2)
      real     ::  col_delta(plonl,0:SIZE(vmr,2),max(1,ncol_abs))               ! layer column densities (molecules/cm^2)
      real     ::  het_rates(plonl,SIZE(vmr,2),max(1,hetcnt))
      real     ::  extfrc(plonl,SIZE(vmr,2),max(1,extcnt))
      real     ::  reaction_rates(plonl,SIZE(vmr,2),rxntot)


      integer   :: e90_tropk
      real, dimension(plonl,SIZE(vmr,2)) :: &
                   h2ovmr, &             ! water vapor volume mixing ratio
                   mbar, &               ! mean wet atmospheric mass ( amu )
                   zmid                  ! midpoint geopotential in km
      real, dimension(plonl,SIZE(zi,2)) :: &
                   zint                  ! interface geopotential in km
      integer :: plev, plevp, plnplv, num_invar
      integer :: nstep

      if ( .not. module_is_initialized ) then
         call error_mesg ('chemdr','chemdr_init must be called first.', FATAL)
      end if
      

      plev = SIZE(vmr,2)
      plevp = SIZE(zi,2)
      plnplv = plonl*plev
      num_invar = SIZE(invariants,3)
      nstep = 0
      
      
!-----------------------------------------------------------------------      
!        ... Initialize xform between mass and volume mixing ratios
!-----------------------------------------------------------------------      
      call inti_mr_xform( sh, mbar, plonl )
!-----------------------------------------------------------------------      
!        ... Xform from mmr to vmr
!-----------------------------------------------------------------------      
!     call mmr2vmr( vmr, mmr, mbar, plonl )
!-----------------------------------------------------------------------      
!        ... Xform water vapor from mmr to vmr and adjust in stratosphere
!-----------------------------------------------------------------------      
      call adjh2o( h2ovmr, sh, mbar, vmr, do_interactive_h2o, plonl )
!-----------------------------------------------------------------------      
!        ... Xform geopotential height from m to km 
!            and pressure from hPa to mb
!-----------------------------------------------------------------------      
      do k = 1,plev
         zmid(:,k) = 1.e-3 * zma(:,k)
         zint(:,k) = 1.e-3 * zi(:,k)
      end do
      zint(:,plevp) = 1.e-3 * zi(:,plevp)

      if( nfs > 0 ) then
!-----------------------------------------------------------------------      
!        ... Set the "invariants"
!-----------------------------------------------------------------------      
         call setinv( invariants, tfld, h2ovmr, pmid, inv_data, &
                      do_interactive_h2o, plonl )
      end if
      if( ncol_abs > 0 .and. phtcnt > 0 ) then
!-----------------------------------------------------------------------      
!        ... Xform family ox assuming that all ox is o3
!-----------------------------------------------------------------------      
!        ox_ndx = get_spc_ndx( 'OX' )
!        if( ox_ndx > 0 ) then
!           o3_ndx = get_grp_mem_ndx( 'O3' )
!           if( o3_ndx > 0 ) then
!              vmr(:,:,ox_ndx) = mbar(:,:) * mmr(:,:,ox_ndx) / nadv_mass(o3_ndx)
!           end if
!        end if
!-----------------------------------------------------------------------      
!        ... Set the column densities at the upper boundary
!-----------------------------------------------------------------------      
         call set_ub_col( col_delta, vmr, invariants, pdel, ptop, plonl )
      end if
      if( gascnt > 0 ) then
!-----------------------------------------------------------------------      
!       ...  Set rates for "tabular" and user specified reactions
!-----------------------------------------------------------------------      
         call setrxt( reaction_rates, tfld, invariants(:,:,indexm), plonl, plev, plnplv )
!        call sulf_interp( lat, ip, pmid, caldayn, sulfate, plonl )
         call usrrxt( reaction_rates, tfld, invariants, h2ovmr, pmid, &
                      invariants(:,:,indexm), sulfate, psc, vmr, sh, delt, plonl , r, &
                      trop_diag_array, trop_option, trop_diag)
!-----------------------------------------------------------------------      
!       ...  Save reaction rate constants for diagnostic output
!-----------------------------------------------------------------------      
         rate_const_out(:,:,:) = reaction_rates(:,:,phtcnt+1:rxntot)
!-----------------------------------------------------------------------      
!       ...  History output for instantaneous reaction rates
!-----------------------------------------------------------------------      
!        do file = 1,moz_file_cnt
!           if( hfile(file)%wrhstts .and. hfile(file)%histout_cnt(9,inst) > 0 ) then
!              do m = 1,hfile(file)%histout_cnt(9,inst)
!                  fldname = hfile(file)%hist_inst(hfile(file)%histout_ind(9,inst)+m-1)
!                  hndx    = hfile(file)%inst_map(hfile(file)%histout_ind(9,inst)+m-1)
!                 call outfld( fldname, reaction_rates(1,1,hndx+phtcnt), plonl, ip, lat, file )
!              end do
!           end if
!-----------------------------------------------------------------------      
!       ...  History output for time averaged reaction rates
!-----------------------------------------------------------------------      
!           do m = 1,hfile(file)%histout_cnt(9,avrg)
!               fldname = hfile(file)%hist_timav(hfile(file)%histout_ind(9,avrg)+m-1)
!               hndx = hfile(file)%timav_map(hfile(file)%histout_ind(9,avrg)+m-1)
!              call outfld( fldname, reaction_rates(1,1,hndx+phtcnt), plonl, ip, lat, file )
!           end do
!        end do
         call adjrxt( reaction_rates, invariants, invariants(:,:,indexm), &
                      plnplv )
      end if
      
      if( phtcnt > 0 ) then
!-----------------------------------------------------------------------
!        ... Compute the photolysis rates
!-----------------------------------------------------------------------      
         if( ncol_abs > 0 ) then
!-----------------------------------------------------------------------      
!             ... Set the column densities
!-----------------------------------------------------------------------      
            call setcol( col_delta, col_dens, pdel, plonl )
         end if
!-----------------------------------------------------------------------      
!             ... Calculate the surface albedo
!-----------------------------------------------------------------------      
!        call srfalb( lat, ip, albs, caldayn, tsurf, plonl )
!-----------------------------------------------------------------------      
!             ... Calculate the photodissociation rates
!-----------------------------------------------------------------------      
!<f1p: replace cwat by cliq+cice
   call mpp_clock_begin(photo_clock_id)
         if (.not. trop_option%do_fastjx_photo) then
           call photo( reaction_rates(:,:,:phtcnt), pmid, pdel, tfld, zmid, &
                     col_dens, &
                     coszen, albedo, &
                     cwat, cldfr, &
                     esfact, solar_phase, plonl )
         else    
#ifndef AM3_CHEM
            call fphoto( reaction_rates(:,:,:phtcnt), &
                         pmid, pdel, &
                         tfld, zmid, &
                         col_dens, &
                         coszen, albedo, &
                         j_ndx, &
                         cwat, &
                         cldfr,  &
                         esfact, solar_phase, plonl,&           
                         trop_option%use_lsc_in_fastjx, &
                         phalf,&
                         zi,&
                         pwt , &
                         sh, &
                         r, &            
                         Time, & 
                         trop_option%time_varying_solarflux)
#else
            call fphoto( reaction_rates(:,:,:phtcnt), &
                         pmid, pdel, &
                         tfld, zmid, &
                         col_dens, &
                         coszen, albedo, &
                         j_ndx, &
                         cwat, &
                         cldfr,  &
                         esfact, solar_phase, plonl,&
                         trop_option%use_lsc_in_fastjx, &
                         phalf,&
                         zi,&
                         pwt , &
                         sh, &
                         r)
#endif
         end if
   call mpp_clock_end(photo_clock_id)

!-----------------------------------------------------------------------      
!       ...  History output for instantaneous photo rates
!-----------------------------------------------------------------------      
!        do file = 1,moz_file_cnt
!           if( hfile(file)%wrhstts .and. hfile(file)%histout_cnt(8,inst) > 0 ) then
!              do m = 1,hfile(file)%histout_cnt(8,inst)
!                 fldname = hfile(file)%hist_inst(hfile(file)%histout_ind(8,inst)+m-1)
!                  hndx    = hfile(file)%inst_map(hfile(file)%histout_ind(8,inst)+m-1)
!                 call outfld( fldname, reaction_rates(1,1,hndx), plonl, ip, lat, file )
!              end do
!           end if
!-----------------------------------------------------------------------      
!       ...  History output for time averaged photo rates
!-----------------------------------------------------------------------      
!           do m = 1,hfile(file)%histout_cnt(8,avrg)
!               fldname = hfile(file)%hist_timav(hfile(file)%histout_ind(8,avrg)+m-1)
!               hndx    = hfile(file)%timav_map(hfile(file)%histout_ind(8,avrg)+m-1)
!              call outfld( fldname, reaction_rates(1,1,hndx), plonl, ip, lat, file )
!           end do
!        end do
!-----------------------------------------------------------------------      
!       ...  Save photolysis rates for diagnostic output
!-----------------------------------------------------------------------      
         jvals_out(:,:,:) = reaction_rates(:,:,:phtcnt)
!-----------------------------------------------------------------------      
!             ... Adjust the photodissociation rates
!-----------------------------------------------------------------------      
         call phtadj( reaction_rates, invariants, invariants(:,:,indexm), &
                      plnplv )
      end if
      if( hetcnt > 0 ) then
!-----------------------------------------------------------------------
!        ... Compute the heterogeneous rates at time = t(n+1)
!-----------------------------------------------------------------------      
!        call sethet( het_rates, pmid, lat, zmid, phis, &
!                     tfld, cmfdqr, nrain, nevapr, delt, &
!                     invariants(1,1,indexm), vmr, plonl )
!-----------------------------------------------------------------------      
!       ...  History output for instantaneous wet removal rates
!-----------------------------------------------------------------------      
!        do file = 1,moz_file_cnt
!           if( hfile(file)%wrhstts .and. hfile(file)%histout_cnt(10,inst) > 0 ) then
!              do m = 1,hfile(file)%histout_cnt(10,inst)
!                 fldname = hfile(file)%hist_inst(hfile(file)%histout_ind(10,inst)+m-1)
!                     hndx    = hfile(file)%inst_map(hfile(file)%histout_ind(10,inst)+m-1)
!                    call outfld( fldname, het_rates(1,1,hndx), plonl, ip, lat, file )
!              end do
!           end if
!-----------------------------------------------------------------------      
!       ...  History output for time averaged wet removal rates
!-----------------------------------------------------------------------      
!           do m = 1,hfile(file)%histout_cnt(10,avrg)
!               fldname = hfile(file)%hist_timav(hfile(file)%histout_ind(10,avrg)+m-1)
!               hndx    = hfile(file)%timav_map(hfile(file)%histout_ind(10,avrg)+m-1)
!              call outfld( fldname, het_rates(1,1,hndx), plonl, ip, lat, file )
!           end do
!        end do
      end if
      if( extcnt > 0 ) then
!-----------------------------------------------------------------------
!        ... Compute the extraneous frcing at time = t(n+1)
!-----------------------------------------------------------------------      
!        call setext( extfrc, lat, ip, zint, cldtop, plonl )
!-----------------------------------------------------------------------      
!       ...  History output for instantaneous external forcing rates
!-----------------------------------------------------------------------      
!        do file = 1,moz_file_cnt
!           if( hfile(file)%wrhstts .and. hfile(file)%histout_cnt(11,inst) > 0 ) then
!              do m = 1,hfile(file)%histout_cnt(11,inst)
!                  fldname = hfile(file)%hist_inst(hfile(file)%histout_ind(11,inst)+m-1)
!                  hndx = hfile(file)%inst_map(hfile(file)%histout_ind(11,inst)+m-1)
!                 call outfld( fldname, extfrc(1,1,hndx), plonl, ip, lat, file )
!              end do
!           end if
!-----------------------------------------------------------------------      
!       ...  History output for time averaged external forcing rates
!-----------------------------------------------------------------------      
!           do m = 1,hfile(file)%histout_cnt(11,avrg)
!               fldname = hfile(file)%hist_timav(hfile(file)%histout_ind(11,avrg)+m-1)
!               hndx    = hfile(file)%timav_map(hfile(file)%histout_ind(11,avrg)+m-1)
!              call outfld( fldname, extfrc(1,1,hndx), plonl, ip, lat, file )
!           end do
!        end do
!        do m = 1,max(1,extcnt)
!           do k = 1,SIZE(vmr,2)
!               extfrc(:,k,m) = extfrc(:,k,m) / invariants(:,k,indexm)
!           end do
!        end do
      end if
      if( grpcnt > 0 ) then
!-----------------------------------------------------------------------
!        ... Set the group ratios
!-----------------------------------------------------------------------      
!        call set_grp_ratios( group_ratios, reaction_rates, vmr, mmr, nas, &
!                              mbar, invariants, plonl )
!-----------------------------------------------------------------------
!             ... Modify the reaction rate of any reaction
!           with group member or proportional reactant(s)
!-----------------------------------------------------------------------
!        call rxt_mod( reaction_rates, het_rates, group_ratios, plnplv )
      end if

!=======================================================================
!        ... Call the class solution algorithms
!=======================================================================
      if( clscnt1 > 0 .and. rxntot > 0 ) then
!-----------------------------------------------------------------------
!        ... Solve for "explicit" species
!-----------------------------------------------------------------------
         call exp_sol( vmr, reaction_rates, &
                       het_rates, extfrc, &
                       nstep, delt, &
!                      invariants(1,1,indexm), &
                       prod_out, loss_out, &
                       plonl, plnplv )
      end if
      if( clscnt4 > 0 .and. rxntot > 0 ) then
!-----------------------------------------------------------------------
!        ... Solve for "Implicit" species
!-----------------------------------------------------------------------
!prod_ox and loss_ox are added to implicit solver. (jmao,1/1/2011)
   call mpp_clock_begin(implicit_clock_id)
         call imp_sol( vmr, reaction_rates, &
                       het_rates, extfrc, &
                       nstep, delt, &
!                      invariants(1,1,indexm), &
                       lat, lon, &
                       prod_out, loss_out, trop_option%check_convergence, &
                       imp_slv_nonconv, &
                       plonl, plnplv, &
                       prod_ox, loss_ox)
   call mpp_clock_end(implicit_clock_id)
      end if
      if( clscnt5 > 0 .and. rxntot > 0 ) then
!-----------------------------------------------------------------------
!        ... Solve for "Rodas" species
!-----------------------------------------------------------------------
         call rodas_sol( vmr, reaction_rates, &
                         het_rates, extfrc, &
                         nstep, delt, &
!                        invariants(1,1,indexm), &
                         plonl, plnplv )
      end if

!-----------------------------------------------------------------------
!       ... Assign O3strat to O3 at and above the tropopause and 
!           let it undergo loss processes below the tropopause
!-----------------------------------------------------------------------
      do i = 1,plonl
        do k = plev-1,2,-1
           if (zmid(i,k) < ztrop_low ) then
               cycle
           else if( zmid(i,k) > ztrop_high ) then
               troplev(i)    = k
               exit
           end if
           dt = tfld(i,k) - tfld(i,k-1)
           if( dt < max_dtdz*(zmid(i,k-1) - zmid(i,k)) ) then
              troplev(i)    = k
              exit
           end if
       end do
     end do
     do i = 1,plonl
       do k = 1,plev
         if (k <= troplev(i)) then
            vmr(i,k,o3s_ndx) = vmr(i,k,ox_ndx)
         else
            k_loss_ox(i,k) = loss_ox(i,k) / (max( vmr(i,k,ox_ndx), 1.0e-20     ))
            vmr(i,k,o3s_ndx) = vmr(i,k,o3s_ndx)*exp( -delt*k_loss_ox(i,k) )
         end if
       end do
     end do

   !O3-Strat with tropopause diagnosed by e90 
   !Criteria defining "stratospheric":
   ! (1) e90 < e90_tropopause_vmr, excluding points when
   ! (2) k >= e90_tropk and O3<100 ppbv, so to avoid labeling
   !     aged trop air in Hadley cell[Prather et al,2011]
   e90_ndx = get_tracer_index( MODEL_ATMOS,'e90' )
   o3s_e90_ndx = get_spc_ndx( 'O3S_E90' )
   ox_ndx = get_spc_ndx( 'O3' )

   if( e90_ndx > 0 .and. o3s_e90_ndx > 0 ) then

       strate90(:,:) = 0  ! default to trop air

       do i = 1,plonl
          e90_tropk = 0
          do k = 1, plev
           if ( e90_vmr(i,k) < e90_tropopause_vmr .and. pmid(i,k)*0.01 < 600.0 ) then
              strate90(i,k) = 1  ! strat air
              if (e90_tropk>0 .and. vmr(i,k,ox_ndx)*1.0e9<100.0) then
                 strate90(i,k) = 0
              endif
           else 
              strate90(i,k) = 0  ! trop air
              if (e90_tropk==0) then
                  e90_tropk=k
              endif
           endif
          end do
       end do
 
       do i = 1,plonl
          do k = 1,plev 
            if (strate90(i,k)==1 ) then
               vmr(i,k,o3s_e90_ndx) = vmr(i,k,ox_ndx)
            else
               k_loss_ox(i,k) = loss_ox(i,k) / (max( vmr(i,k,ox_ndx), 1.0e-20 ))
               vmr(i,k,o3s_e90_ndx) = vmr(i,k,o3s_e90_ndx)*exp( -delt*k_loss_ox(i,k) )
            end if
          end do
       end do 
 
   end if

!-----------------------------------------------------------------------
!       ... Heterogeneous+cloud chemistry
!-----------------------------------------------------------------------
   call mpp_clock_begin(cloud_clock_id)
      if( so2_ndx > 0 .and. so4_ndx > 0 ) then
         call setsox( pmid, plonl, delt, tfld, sh, &
                      cwat, fliq, cldfr, invariants(:,:,indexm), &
                      vmr, co2, &
                      trop_diag_array,trop_option,trop_diag)
      end if
   call mpp_clock_end(cloud_clock_id)

!-----------------------------------------------------------------------      
!         ... Check for negative values and reset to zero
!-----------------------------------------------------------------------      
!     call negtrc( lat, 'After chemistry ', vmr, plonl )
!-----------------------------------------------------------------------      
!         ... Output instantaneous "wet" advected volume mixing
!-----------------------------------------------------------------------      
!     do file = 1,moz_file_cnt
!        if( hfile(file)%wrhstts .and. hfile(file)%histout_cnt(1,inst) > 0 ) then
!           do m = 1,hfile(file)%histout_cnt(1,inst)
!               fldname = hfile(file)%hist_inst(hfile(file)%histout_ind(1,inst)+m-1)
!               hndx    = hfile(file)%inst_map(hfile(file)%histout_ind(1,inst)+m-1)
!              call outfld( fldname, vmr(1,1,hndx), plonl, ip, lat, file )
!           end do
!        end if
!-----------------------------------------------------------------------      
!         ... Output time averaged "wet" advected volume mixing ratios
!-----------------------------------------------------------------------      
!        do m = 1,hfile(file)%histout_cnt(1,avrg)
!            fldname = hfile(file)%hist_timav(hfile(file)%histout_ind(1,avrg)+m-1)
!            hndx    = hfile(file)%timav_map(hfile(file)%histout_ind(1,avrg)+m-1)
!           call outfld( fldname, vmr(1,1,hndx), plonl, ip, lat, file )
!        end do
!     end do
!-----------------------------------------------------------------------      
!         ... Output instantaneous "wet" non-advected volume mixing
!-----------------------------------------------------------------------      
!     group_write(:moz_file_cnt) = hfile(:moz_file_cnt)%wrhstts .and. &
!                                  hfile(:moz_file_cnt)%histout_cnt(2,inst) > 0
!     if( ANY( group_write(:moz_file_cnt) ) .or. &
!         ANY( hfile(:moz_file_cnt)%histout_cnt(2,avrg) > 0 ) ) then
!        call mak_grp_vmr( vmr, group_ratios(1,1,1), group_vmr(1,1,1), plonl )
!     end if
!     do file = 1,moz_file_cnt
!        if( hfile(file)%wrhstts .and. hfile(file)%histout_cnt(2,inst) > 0 ) then
!           do m = 1,hfile(file)%histout_cnt(2,inst)
!               fldname = hfile(file)%hist_inst(hfile(file)%histout_ind(2,inst)+m-1)
!               hndx    = hfile(file)%inst_map(hfile(file)%histout_ind(2,inst)+m-1)
!              call outfld( fldname, group_vmr(1,1,hndx), plonl, ip, lat, file )
!           end do
!        end if
!-----------------------------------------------------------------------      
!         ... Output time averaged "wet" non-advected volume mixing ratios
!-----------------------------------------------------------------------      
!        do m = 1,hfile(file)%histout_cnt(2,avrg)
!            fldname = hfile(file)%hist_timav(hfile(file)%histout_ind(2,avrg)+m-1)
!            hndx    = hfile(file)%timav_map(hfile(file)%histout_ind(2,avrg)+m-1)
!           call outfld( fldname, group_vmr(1,1,hndx), plonl, ip, lat, file )
!        end do
!     end do
!-----------------------------------------------------------------------      
!         ... Xform from vmr to mmr
!-----------------------------------------------------------------------      
!     call vmr2mmr( vmr, mmr, nas, group_ratios, mbar, plonl )

      end subroutine chemdr
!</SUBROUTINE>

!</SUBROUTINE>
!<f1p
      subroutine chemdr_init(trop_option)

        type(tropchem_opt), intent(in) :: trop_option        

        integer :: n
        real    :: nit_scale_factor

        nh4no3_ndx   = get_spc_ndx('NH4NO3')
        nh4_ndx      = get_spc_ndx('NH4')
        nh3_ndx      = get_spc_ndx('NH3')
        hno3_ndx     = get_spc_ndx('HNO3')
        so2_ndx      = get_spc_ndx( 'SO2' )
        so4_ndx      = get_spc_ndx( 'SO4' )
        ox_ndx       = get_spc_ndx( 'O3' )
        o3s_ndx      = get_spc_ndx( 'O3S' )

        e90_ndx = get_tracer_index( MODEL_ATMOS,'e90' )
        o3s_e90_ndx = get_spc_ndx( 'O3S_E90' )

        implicit_clock_id = mpp_clock_id('Chemistry: Implicit solver')
        photo_clock_id = mpp_clock_id('Chemistry: Photolysis')
        cloud_clock_id = mpp_clock_id('Chemistry: Cloud Chemistry')

        module_is_initialized = .true.
     
      end subroutine chemdr_init
!>
      end module mo_chemdr_mod
