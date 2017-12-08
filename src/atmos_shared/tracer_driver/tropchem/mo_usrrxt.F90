      module mo_usrrxt_mod

      use sat_vapor_pres_mod, only : compute_qs
      use constants_mod, only : rdgas, rvgas
      use strat_chem_utilities_mod, only : psc_type, strat_chem_get_gamma, &
                                           strat_chem_get_hetrates
      use field_manager_mod,  only : MODEL_ATMOS   
      use tracer_manager_mod, only : get_tracer_index,  query_method
      use field_manager_mod,  only: parse             
      use tropchem_types_mod, only : tropchem_opt, tropchem_diag
      use fms_mod,    only : open_file, close_file

implicit none
      public :: usrrxt_init, usrrxt
      public :: HET_CHEM_LEGACY, HET_CHEM_J1M

      private
      integer, parameter     :: HET_CHEM_LEGACY    = 1
      integer, parameter     :: HET_CHEM_J1M       = 2

      integer :: uo_o2_ndx, uno2_no3_ndx, un2o5_ndx, uoh_hno3_ndx, uho2_no2_ndx, uhno4_ndx,&
                 uco_oha_ndx, uho2_ho2_ndx, upan_f_ndx, upan_b_ndx, umpan_f_ndx, umpan_b_ndx, &
                 n2o5h_ndx, no3h_ndx, ho2h_ndx, no2h_ndx, nh3h_ndx, uoh_xooh_ndx, uoh_acet_ndx, &
                 uoh_dms_ndx,&
                 so4_ndx, bc1_ndx, bc2_ndx, oc1_ndx, oc2_ndx, soa_ndx,nh4_ndx,nh4no3_ndx,&
                 ssa_ndx(5), dust_ndx(5),&
                 h2o_ndx, hcl_ndx, clono2_ndx, hbr_ndx, &
                 strat37_ndx, strat38_ndx, strat72_ndx, strat73_ndx, strat74_ndx, &
                 strat75_ndx, strat76_ndx, strat77_ndx, strat78_ndx, strat79_ndx, &
                 strat80_ndx

      real, parameter :: d622 = rdgas/rvgas
      real, parameter :: d378 = 1. - d622
! setup for heteorogenous chemistry (jmao, 02/28/2012)      
      integer, parameter           ::       A_HET     =1000     !dim = no. of Aerosol/cloud Mie sets FOR AM3(input data)     
      character(len=78)            ::      TITLE0_HET           !Title of the first line in am3_scat.dat
      character(len=20)            ::      TITL_HET(A_HET)      ! 
      integer                      ::      NAA_HET              !Number of categories for scattering phase functions for am3, is 880
      real(kind=8)                 ::      MEE_HET(5,A_HET)     !Aerosol mass extinction efficiency, MEE*colume mass =od
      real(kind=8)                 ::      SAA_HET(5,A_HET)     !Single scattering albedo
      real(kind=8)                 ::      PAA_HET(8,5,A_HET)   !Phase function: first 8 terms of expansion      
      real(kind=8)                 ::      WAA_HET(5,A_HET)     !Wavelengths for the NK supplied phase functions(nm)
      real(kind=8)                 ::      RAA_HET(A_HET)

character(len=128), parameter :: version     = '$Id$'
character(len=128), parameter :: tagname     = '$Name$'
logical                       :: module_is_initialized = .false.

      contains

      subroutine usrrxt_init( verbose, trop_option )
!-----------------------------------------------------------------
!        ... Intialize the user reaction constants module
!-----------------------------------------------------------------

      use mo_chem_utls_mod, only : get_rxt_ndx, get_spc_ndx

      implicit none

!-----------------------------------------------------------------------
!       ... Dummy arguments
!-----------------------------------------------------------------------
      integer,          intent(in) :: verbose
      type(tropchem_opt), intent(in) :: trop_option

!-----------------------------------------------------------------
!        ... local variables
!-----------------------------------------------------------------
      integer :: funit, I, J, K

      uo_o2_ndx = get_rxt_ndx( 'uo_o2' )
      uno2_no3_ndx = get_rxt_ndx( 'uno2_no3' )
      un2o5_ndx = get_rxt_ndx( 'un2o5' )
      uoh_hno3_ndx = get_rxt_ndx( 'uoh_hno3' )
      uho2_no2_ndx = get_rxt_ndx( 'uho2_no2' )
      uhno4_ndx = get_rxt_ndx( 'uhno4' )
      uco_oha_ndx = get_rxt_ndx( 'uco_oha' )
      uho2_ho2_ndx = get_rxt_ndx( 'uho2_ho2' )
      upan_f_ndx = get_rxt_ndx( 'upan_f' )
      upan_b_ndx = get_rxt_ndx( 'upan_b' )
      umpan_f_ndx = get_rxt_ndx( 'umpan_f' )
      umpan_b_ndx = get_rxt_ndx( 'umpan_b' )
      n2o5h_ndx = get_rxt_ndx( 'n2o5h' )
      no3h_ndx = get_rxt_ndx( 'no3h' )
      ho2h_ndx = get_rxt_ndx( 'ho2h' )
      no2h_ndx = get_rxt_ndx( 'no2h' )
      nh3h_ndx = get_rxt_ndx( 'nh3h' )
      uoh_xooh_ndx = get_rxt_ndx( 'uoh_xooh' )
      uoh_acet_ndx = get_rxt_ndx( 'uoh_acet' )
      uoh_dms_ndx = get_rxt_ndx( 'uoh_dms' )
      strat37_ndx = get_rxt_ndx( 'strat37' )
      strat38_ndx = get_rxt_ndx( 'strat38' )
      strat72_ndx = get_rxt_ndx( 'strat72' )
      strat73_ndx = get_rxt_ndx( 'strat73' )
      strat74_ndx = get_rxt_ndx( 'strat74' )
      strat75_ndx = get_rxt_ndx( 'strat75' )
      strat76_ndx = get_rxt_ndx( 'strat76' )
      strat77_ndx = get_rxt_ndx( 'strat77' )
      strat78_ndx = get_rxt_ndx( 'strat78' )
      strat79_ndx = get_rxt_ndx( 'strat79' )
      strat80_ndx = get_rxt_ndx( 'strat80' )

      !Note here we cannot use get_spc_ndx, because aerosols are not
      !tracnam, which is for get_spc_ndx in mo_chem_utls.F90. (jmao, 03/16/2012)
!      so4_ndx     = get_tracer_index(MODEL_ATMOS,'so4')
if (trop_option%het_chem .eq. HET_CHEM_LEGACY) then
      so4_ndx = get_spc_ndx( 'SO4' )
elseif ( trop_option%het_chem .eq. HET_CHEM_J1M) then
! check wiki for this, search get_tracer_index.
      so4_ndx     = get_tracer_index(MODEL_ATMOS,'so4')
      nh4_ndx     = get_tracer_index(MODEL_ATMOS,'nh4')
      nh4no3_ndx  = get_tracer_index(MODEL_ATMOS,'nh4no3')
      bc1_ndx     = get_tracer_index(MODEL_ATMOS,'bcphob')
      bc2_ndx     = get_tracer_index(MODEL_ATMOS,'bcphil')
      oc1_ndx     = get_tracer_index(MODEL_ATMOS,'omphob')
      oc2_ndx     = get_tracer_index(MODEL_ATMOS,'omphil')
      soa_ndx     = get_tracer_index(MODEL_ATMOS,'SOA')
      ssa_ndx(1)    = get_tracer_index(MODEL_ATMOS,'ssalt1')
      ssa_ndx(2)    = get_tracer_index(MODEL_ATMOS,'ssalt2')
      ssa_ndx(3)    = get_tracer_index(MODEL_ATMOS,'ssalt3')
      ssa_ndx(4)    = get_tracer_index(MODEL_ATMOS,'ssalt4')
      ssa_ndx(5)    = get_tracer_index(MODEL_ATMOS,'ssalt5')
      dust_ndx(1)   = get_tracer_index(MODEL_ATMOS,'dust1')
      dust_ndx(2)   = get_tracer_index(MODEL_ATMOS,'dust2') 
      dust_ndx(3)   = get_tracer_index(MODEL_ATMOS,'dust3')
      dust_ndx(4)   = get_tracer_index(MODEL_ATMOS,'dust4')
      dust_ndx(5)   = get_tracer_index(MODEL_ATMOS,'dust5')
end if      
      h2o_ndx = get_spc_ndx( 'H2O' )
      hcl_ndx = get_spc_ndx( 'HCl' )
      clono2_ndx = get_spc_ndx( 'ClONO2' )
      hbr_ndx = get_spc_ndx( 'HBr' )

      if (verbose >= 2) then
      write(*,*) 'usrrxt_init: diagnostics '
      write(*,'(10i5)') uo_o2_ndx, uno2_no3_ndx, un2o5_ndx, uoh_hno3_ndx, uho2_no2_ndx, uhno4_ndx, &
                 uco_oha_ndx, uho2_ho2_ndx, upan_f_ndx, upan_b_ndx, umpan_f_ndx, umpan_b_ndx, &
                 n2o5h_ndx, no3h_ndx, uoh_xooh_ndx, uoh_acet_ndx, &
                 uoh_dms_ndx, nh3h_ndx, &
                 strat37_ndx, strat38_ndx, strat72_ndx, strat73_ndx, strat74_ndx, &
                 strat75_ndx, strat76_ndx, strat77_ndx, strat78_ndx, strat79_ndx, &
                 strat80_ndx
      end if
!----------------------------------------------------------------------
!Read a new set of aerosol density and effective radius from FASTJX (jmao, 03/14/2012)
!------------------------------------------------------------------------
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)        
!     funit  Channel number for reading data file                     
!     NAA      Number of categories for scattering phase functions      
!     QAA      Aerosol scattering phase functions                       
!     NK       Number of wavelengths at which functions supplied (set as
!     WAA      Wavelengths for the NK supplied phase functions          
!     PAA      Phase function: first 8 terms of expansion               
!     RAA      Effective radius associated with aerosol type            
!     SAA      Single scattering albedo                                 
!-----------------------------------------------------------------------
      funit = open_file(FILE='INPUT/am3_uptake.dat',form='formatted',action='read',threading='multi', &
                        dist=.false.)

      read (funit,'(i4,a78)') NAA_HET, TITLE0_HET 
      if (NAA_HET .gt. A_HET) then 
         write(*,*) 'ATMOS:fastjx_init: too many scat-data sets for AM3:', NAA_HET, A_HET 
         stop 
      endif
 
      do J = 1,NAA_HET
         read (funit,'(3x,a20,32x,f6.3)')                       &
              &      TITL_HET(J),RAA_HET(J)                                
         ! ver 6.0 extend to 5 ref wavelengths for mie-sca
         do K = 1,5 
            read (funit,'(f4.0,2e14.6,7f6.3,1x,f7.3,f8.4)')              &
                 WAA_HET(K,J),MEE_HET(K,J),SAA_HET(K,J),(PAA_HET(I,K,J),I=2,8)                   
            PAA_HET(1,K,J) = 1.d0 
         enddo
       enddo  
       
       call close_file(funit,dist=.false.)

    end subroutine usrrxt_init

      subroutine usrrxt( rxt, temp, invariants, h2ovmr, pmid, m, &
                         sulfate, psc, qin, sh, delt, plonl ,r, trop_diag_array, trop_option, trop_diag)
!-----------------------------------------------------------------
!        ... set the user specified reaction rates
!-----------------------------------------------------------------

#ifndef AM3_CHEM
      use chem_mods_mod, only : nfs, rxntot, indexh2o
#else
      use AM3_chem_mods_mod, only : nfs, rxntot, indexh2o
#endif
      use constants_mod, only : PI

      implicit none

!-----------------------------------------------------------------
!        ... dummy arguments
!-----------------------------------------------------------------
      integer, intent(in) :: plonl
      type(tropchem_opt), intent(in)  :: trop_option
      type(tropchem_diag), intent(in) :: trop_diag
      real, intent(inout)             :: trop_diag_array(:,:,:)
      real,    intent(in) :: qin(:,:,:)            ! transported species ( vmr )
      real,    intent(in) :: temp(:,:), &          ! temperature
                             m(:,:), &             ! total atm density
                             sulfate(:,:), &       ! sulfate aerosol vmr
                             h2ovmr(:,:), &        ! water vapor vmr
                             pmid(:,:), &          ! midpoint pressure in pa
                             sh(:,:), &            ! specific humidity
                             invariants(:,:,:)     ! invariants density
      type(psc_type), intent(in) :: &
                             psc                   ! polar stratospheric clouds (PSCs)
      real,    intent(in) :: delt                  ! timestep (sec)
      real, intent(inout) :: rxt(:,:,:)            ! gas phase rates
      real,    intent(in) :: r(:,:,:)              ! r is for aerosol concentrations(jmao, 03/18/2012)
      
!-----------------------------------------------------------------
!        ... local variables
!-----------------------------------------------------------------
      real, parameter :: boltz = 1.38044e-16            ! erg / k
      real, parameter :: avo   = 6.023e23               ! molecules/mole
!-----------------------------------------------------------------
!        ... density of sulfate aerosol
!-----------------------------------------------------------------
!     real, parameter :: gam1 = 0.04                    ! n2o5+sul ->2hno3
!     real, parameter :: gam1 = 0.10                    ! n2o5+sul ->2hno3
!     real, parameter :: gam4 = 0.05                    ! NH3 +SUL ->NH4SO4 (Dentener 1994)
      real, parameter :: wso4 = 98.
      real, parameter :: den  = 1.15                    ! each molecule of so4(aer) density g/cm3
!-------------------------------------------------
!         ... volume of sulfate particles
!           assuming mean rm 
!           continient 0.05um  0.07um  0.09um
!           ocean      0.09um  0.25um  0.37um
!                      0.16um                  blake jgr,7195, 1995
!-------------------------------------------------
      real, parameter :: rm1  = 0.16*1.e-4                   ! mean radii in cm
      real, parameter :: fare = 4.*3.14*rm1*rm1              ! each mean particle(r=0.1u) area   cm2/cm3
      real, parameter :: dg   = 0.1                          ! mole diffusion =0.1 cm2 (Dentener, 1993)

      integer  ::  i, k, ilev, n
      real     ::  amas
      real, dimension( SIZE(temp,1) ) :: &
                   tp, &                    ! 300/t
                   tinv, &                  ! 1/t
                   ko, &
                   kinf, &
                   fc, &
                   xr, &                    ! factor to increase particle radii depending on rel hum
                   sur, &                   ! sulfate particle surface area (cm^2/cm^3)
                   exp_fac, &               ! vector exponential
                   tmp_hcl, &               ! temporary array for HCl VMR
                   tmp_clono2, &            ! temporary array for ClONO2 VMR
                   tmp_hbr                  ! temporary array for HBr VMR
      integer, parameter :: nhet_strat = 9
      real, dimension(SIZE(temp,1), nhet_strat) :: &
                   gamma_strat, &           ! reaction probabilities for strat. het. rxns
                   hetrates_strat           ! effective second-order reaction rates for strat. het. rxns
      integer ::   plev
      real, dimension(SIZE(temp,1),SIZE(temp,2)) :: &
                   relhum                   ! relative humidity
      INTEGER :: tmp_indexh2o

      integer, parameter :: naero_het = 18
      real, parameter :: mw_HO2 = 33.
      real, parameter :: mw_n2o5= 108.
      real, parameter :: mw_no3 = 62.
      real, parameter :: mw_no2 = 46.
      real, parameter :: mw_nh3 = 17.
      real, dimension(size(qin,1), naero_het)::drymass_het,&
                                        rd_het,re_het,sfca_het
      real :: uptk_het
      real :: gam_n2o5, gam_no3, gam_nh3

      plev = SIZE(temp,2)
      ilev = SIZE(temp,1)
      
      amas = 4.*PI*rm1**3*den/3.            ! each mean particle(r=0.1u) mass (g)
!-----------------------------------------------------------------
!        ... o + o2 + m --> o3 + m
!-----------------------------------------------------------------
      do k = 1,plev
         tinv(:)           = 1. / temp(:,k)
         tp(:)             = 300. * tinv(:)
         if( uo_o2_ndx > 0 ) then
            rxt(:,k,uo_o2_ndx) = 6.e-34 * tp(:)**2.4
         end if
!-----------------------------------------------------------------
!        ... n2o5 + m --> no2 + no3 + m
!-----------------------------------------------------------------
         if( un2o5_ndx > 0 ) then
            if( uno2_no3_ndx > 0 ) then
               rxt(:,k,un2o5_ndx) = rxt(:,k,uno2_no3_ndx) * 3.704e26 * exp( -11000.*tinv(:) )
            else
               rxt(:,k,un2o5_ndx) = 0.
            end if
         end if

!-----------------------------------------------------------------
!        set rates for:
!         ... hno3 + oh --> no3 + h2o
!           ho2no2 + m --> ho2 + no2 + m
!           co + oh --> co2 + ho2
!-----------------------------------------------------------------
         if( uoh_hno3_ndx > 0 ) then
            ko(:) = m(:,k) * 6.5e-34 * exp( 1335.*tinv(:) )
            ko(:) = ko(:) / (1. + ko(:)/(2.7e-17*exp( 2199.*tinv(:) )))
            rxt(:,k,uoh_hno3_ndx) = ko(:) + 2.4e-14*exp( 460.*tinv(:) )
         end if
         if( uhno4_ndx > 0 ) then
            if( uho2_no2_ndx > 0 ) then
               rxt(:,k,uhno4_ndx) = rxt(:,k,uho2_no2_ndx) * exp( -10900.*tinv(:) )/ 2.1e-27
            else
               rxt(:,k,uhno4_ndx) = 0.
            end if
         end if
!        if( uco_oha_ndx > 0 ) then
!           rxt(:,k,uco_oha_ndx) = 1.5e-13 * (1. + 6.e-7*boltz*m(:,k)*temp(:,k))
!        end if

!-----------------------------------------------------------------
!            ... mco3 + no2 -> mpan
!-----------------------------------------------------------------
         if( umpan_f_ndx > 0 ) then
            rxt(:,k,umpan_f_ndx) = 9.3e-12 * tp(:) / m(:,k)
         end if

!-----------------------------------------------------------------
!        ... pan + m --> ch3co3 + no2 + m
!-----------------------------------------------------------------
         exp_fac(:) = exp( -14000.*tinv(:) )
         if( upan_b_ndx > 0 ) then
            if( upan_f_ndx > 0 ) then
               rxt(:,k,upan_b_ndx) = rxt(:,k,upan_f_ndx) * 1.111e28 * exp_fac(:)
            else
               rxt(:,k,upan_b_ndx) = 0.
            end if
         end if

!-----------------------------------------------------------------
!        ... mpan + m --> mco3 + no2 + m
!-----------------------------------------------------------------
         if( umpan_b_ndx > 0 ) then
            if( umpan_f_ndx > 0 ) then
               rxt(:,k,umpan_b_ndx) = rxt(:,k,umpan_f_ndx) * 1.111e28 * exp_fac(:)
            else
               rxt(:,k,umpan_b_ndx) = 0.
            end if
         end if

!-----------------------------------------------------------------
!       ... xooh + oh -> h2o + oh
!   This reaction has been removed in AM4. We put it here only for AM3.
!-----------------------------------------------------------------
         if( uoh_xooh_ndx > 0 ) then
            rxt(:,k,uoh_xooh_ndx) = temp(:,k)**2 * 7.69e-17 * exp( 253.*tinv(:) )
         end if

!-----------------------------------------------------------------
!       ... ch3coch3 + oh -> ro2 + h2o
!-----------------------------------------------------------------
         if( uoh_acet_ndx > 0 ) then
            rxt(:,k,uoh_acet_ndx) = 3.82e-11 * exp( -2000.*tinv(:) ) &
                                 + 1.33e-13
         end if
!-----------------------------------------------------------------
!        ... Cl2O2 + M -> 2*ClO + M
!-----------------------------------------------------------------
         if( strat38_ndx > 0 ) then
            if( strat37_ndx > 0 ) then
               rxt(:,k,strat38_ndx) = rxt(:,k,strat37_ndx) * 1.075e27 * exp( -8835.*tinv(:) )
            else
               rxt(:,k,strat38_ndx) = 0.
            end if
         end if

if (trop_option%het_chem .eq. HET_CHEM_LEGACY) then
!-----------------------------------------------------------------
!        ... ho2 + ho2 --> h2o2
!        note: this rate involves the water vapor number density
!-----------------------------------------------------------------
         if( uho2_ho2_ndx > 0 ) then
            if( indexh2o > 0 ) then
               tmp_indexh2o = indexh2o
               fc(:)   = 1. + 1.4e-21 * invariants(:,k,tmp_indexh2o) * exp( 2200.*tinv(:) )
            else if( h2o_ndx > 0 ) then
               fc(:)   = 1. + 1.4e-21 * qin(:,k,h2o_ndx) * m(:,k) * exp( 2200.*tinv(:) )
            else
               fc(:) = 1.
            end if
            ko(:)   = 3.5e-13 * exp( 430.*tinv(:) )
            kinf(:) = 1.7e-33 * m(:,k) * exp( 1000.*tinv(:) )
            rxt(:,k,uho2_ho2_ndx) = (ko(:) + kinf(:)) * fc(:)
         end if
!-----------------------------------------------------------------
!       ... DMS + OH -> .75 * SO2
!-----------------------------------------------------------------
         if( uoh_dms_ndx > 0 ) then
            ko(:) = 1. + 5.0e-30 * exp( 6280.*tinv(:) ) * m(:,k) * 0.21
            rxt(:,k,uoh_dms_ndx) = 1.0e-39 * exp( 5820.*tinv(:) ) &
                                 * m(:,k) * 0.21 / ko(:)
         end if


         if( n2o5h_ndx > 0 .or. no3h_ndx > 0 .or. nh3h_ndx > 0 ) then
!-----------------------------------------------------------------
!         ... n2o5 --> 2*hno3
!             no3 --> hno3
!This reaction is updated from JPL11. (jmao,04/30/2013)         
!-----------------------------------------------------------------
!        ... first compute the relative humidity
!-----------------------------------------------------------------
!           call aqsat( temp(1,k), pmid(1,k), satv, satq, plonl, &
!                       plonl, 1, 1, 1 )
!           relhum(:) = .622 * h2ovmr(:,k) / satq(:)
!           relhum(:) = max( 0.,min( 1.,relhum(:) ) )
            call rh_calc( pmid(:,k), temp(:,k), sh(:,k), relhum(:,k) )
!-------------------------------------------------------------------------
!         ... estimate humidity effect on aerosols (from shettle and fenn, 1979)
!           xr is a factor of the increase aerosol radii with hum (hum=0., factor=1)
!-------------------------------------------------------------------------
            xr(:)     = .999151 + relhum(:,k)*(1.90445 + relhum(:,k)*(-6.35204 + relhum(:,k)*5.32061))
!-------------------------------------------------------------------------
!         ... estimate sulfate particles surface area (cm2/cm3) in each grid
!-------------------------------------------------------------------------
            gam_n2o5 =  trop_option%gN2O5
            gam_no3  =  trop_option%gNO3
            gam_nh3  =  trop_option%gNH3

            if( so4_ndx > 0 ) then
               sur(:)    = qin(:,k,so4_ndx)
            else
               sur(:)    = sulfate(:,k)
            end if
            sur(:)    = sur(:)*m(:,k)/avo*wso4 &              ! xform mixing ratio to g/cm3
                        / amas &                                    ! xform g/cm3 to num particels/cm3
                        * fare &                                    ! xform num particels/cm3 to cm2/cm3
                        * xr(:)*xr(:)                               ! humidity factor
!-----------------------------------------------------------------
!        ... compute the "aerosol" reaction rates
!-----------------------------------------------------------------
!             k = gam * a * velo/4
!
!       where velo = sqrt[ 8*bk*t/pi/(w/av) ]
!             bk = 1.381e-16
!             av = 6.02e23
!             w  = 108 (n2o5)  ho2(33)  ch2o (30)  nh3(15)  
!
!       so that velo = 1.40e3*sqrt(t)  (n2o5)   gama=0.1
!       so that velo = 2.53e3*sqrt(t)  (ho2)    gama>0.2
!       so that velo = 2.65e3*sqrt(t)  (ch2o)   gama>0.022
!       so that velo = 3.75e3*sqrt(t)  (nh3)    gama=0.4
!--------------------------------------------------------
!           xr(:) = .25 * gam1 * sur(:) * 1.40e3 * sqrt( temp(:,k) )
            xr(:) = 1./(rm1/dg + 4./(gam_n2o5+1.e-30)/(1.40e3 * sqrt( temp(:,k))))*sur(:)
            if( n2o5h_ndx > 0 ) then
               rxt(:,k,n2o5h_ndx) = xr(:)
            end if
            if( no3h_ndx > 0 ) then
               rxt(:,k,no3h_ndx) = xr(:)
            end if
            if( nh3h_ndx > 0 ) then
               rxt(:,k,nh3h_ndx) = &
                  1./(rm1/dg + 4./(gam_nh3+1.e-30)/(3.75e3 * sqrt( temp(:,k))))*sur(:)
            end if
         end if
elseif ( trop_option%het_chem .eq. HET_CHEM_J1M) then
!-----------------------------------------------------------------
!        ... ho2 + ho2 --> h2o2
!        note: this rate involves the water vapor number density
!        This reaction is updated from JPL11. (jmao,04/30/2013)         
!-----------------------------------------------------------------
         if( uho2_ho2_ndx > 0 ) then
            if( indexh2o > 0 ) then
               tmp_indexh2o = indexh2o
               fc(:)   = 1. + 1.4e-21 * invariants(:,k,tmp_indexh2o) * exp( 2200.*tinv(:) )
            else if( h2o_ndx > 0 ) then
               fc(:)   = 1. + 1.4e-21 * qin(:,k,h2o_ndx) * m(:,k) * exp( 2200.*tinv(:) )
            else
               fc(:) = 1.
            end if
            ko(:)   = 3.0e-13 * exp( 460.*tinv(:) )
            kinf(:) = 2.1e-33 * m(:,k) * exp( 920.*tinv(:) )
            rxt(:,k,uho2_ho2_ndx) = (ko(:) + kinf(:)) * fc(:)
         end if
!-----------------------------------------------------------------
!       ... DMS + OH -> .75 * SO2 
!       This reaction is updated from JPL11. (jmao,04/30/2013)         
!-----------------------------------------------------------------
         if( uoh_dms_ndx > 0 ) then
            ko(:) = 1. + 1.05e-5 * exp( 3644.*tinv(:) ) * 0.21
            rxt(:,k,uoh_dms_ndx) = 8.2e-39 * exp( 5376.*tinv(:) ) & 
                                 * m(:,k) * 0.21 / ko(:)
         end if

!-----------------------------------------------------------------
!        ... first compute the relative humidity
!-----------------------------------------------------------------
            call rh_calc( pmid(:,k), temp(:,k), sh(:,k), relhum(:,k) )
!-------------------------------------------------------------------------
!         ... estimate humidity effect on aerosols (from shettle and fenn, 1979)
!           xr is a factor of the increase aerosol radii with hum (hum=0., factor=1)
!-------------------------------------------------------------------------
            xr(:)     = .999151 + relhum(:,k)*(1.90445 + relhum(:,k)*(-6.35204 + relhum(:,k)*5.32061))
!My new calculation for HO2,NO3,N2O5 uptake(jmao, 03/23/2012)
!usr16: N2O5 -> 2HNO3
!usr17: NO3 -> HNO3
!usr18: HO2 -> water
!usr19: NO2 -> 0.5HNO3 + 0.5HONO
!----------------------------------------------------------------------------------------
            ! calculate surface area for each kind of aerosol (total 18)
            call set_aerosol(r(:,k,:),relhum(:,k),m(:,k),drymass_het(:,:),&
                 rd_het(:,:),re_het(:,:),sfca_het(:,:))            


            do i=1,ilev
               if( n2o5h_ndx > 0 ) then
                  rxt(i,k,n2o5h_ndx)=0.
                 if ( trop_option%gN2O5 .gt. 0. ) then
                    do n=1, naero_het
                       uptk_het = 0.
                        call calc_hetrate(sfca_het(i,n),re_het(i,n)*1.D-4,m(i,k),trop_option%gN2O5, &
                          sqrt( temp(i,k)),sqrt(mw_n2o5),uptk_het)
                        rxt(i,k,n2o5h_ndx) = rxt(i,k,n2o5h_ndx) + uptk_het
!                     if (trop_diag%ind_gn2o5 .gt. 0) then
!                        trop_diag_array(i,k,trop_diag%ind_gn2o5) = gam_n2o5
!                     end if
                    end do
                 end if
               end if

               if( no3h_ndx > 0 ) then
                  rxt(i,k,no3h_ndx)=0.
                  if ( trop_option%gNO3 .gt. 0. ) then
                     do n=1, naero_het
                        uptk_het = 0.
                        ! we need to make sure the effective radius unit is cm.
                        call calc_hetrate(sfca_het(i,n),re_het(i,n)*1.D-4,m(i,k),trop_option%gNO3, &
                             sqrt( temp(i,k)),sqrt(mw_no3),uptk_het)
                        rxt(i,k,no3h_ndx) = rxt(i,k,no3h_ndx) + uptk_het
                     end do
                  end if
               end if

               if ( ho2h_ndx > 0) then
                  rxt(i,k,ho2h_ndx)=0.
                  if ( trop_option%gHO2 .gt. 0. ) then
                     do n=1, naero_het
                        uptk_het = 0.
                        ! we need to make sure the effective radius unit is cm.
                        call calc_hetrate(sfca_het(i,n),re_het(i,n)*1.D-4,m(i,k),trop_option%gHO2, &
                             sqrt( temp(i,k)),sqrt(mw_ho2),uptk_het)
                        rxt(i,k,ho2h_ndx) = rxt(i,k,ho2h_ndx) + uptk_het
                     end do
                  end if
               end if

            if ( no2h_ndx > 0) then
               rxt(i,k,no2h_ndx)=0.             
               if ( trop_option%gNO2 .gt. 0. ) then
                  do n=1, naero_het
                     uptk_het = 0.
                  ! we need to make sure the effective radius unit is cm.
                     call calc_hetrate(sfca_het(i,n),re_het(i,n)*1.D-4,m(i,k),trop_option%gNO2, &
                          sqrt( temp(i,k)),sqrt(mw_no2),uptk_het)
                     rxt(i,k,no2h_ndx) = rxt(i,k,no2h_ndx) + uptk_het
                  end do
               end if
            end if

            !NH3+SO4->NH4SO4,(NH4)2SO4
            if( nh3h_ndx > 0 ) then
               rxt(i,k,nh3h_ndx)=0.             
               if ( trop_option%gNH3 .gt. 0. ) then
                  do n=1, naero_het
                     uptk_het = 0.
                  ! we need to make sure the effective radius unit is cm.
                     call calc_hetrate(sfca_het(i,n),re_het(i,n)*1.D-4,m(i,k),trop_option%gNH3, &
                          sqrt( temp(i,k)),sqrt(mw_nh3),uptk_het   )
                     rxt(i,k,nh3h_ndx) = rxt(i,k,nh3h_ndx) + uptk_het
                  end do
               end if
            end if
           end do ! (ilev)
end if ! (trop_option%het_chem)
         if( strat72_ndx > 0 .or. strat73_ndx > 0 .or. strat74_ndx > 0 .or. &
             strat75_ndx > 0 .or. strat76_ndx > 0 .or. strat77_ndx > 0 .or. &
             strat78_ndx > 0 .or. strat79_ndx > 0 .or. strat80_ndx > 0 ) then

            if (hcl_ndx>0) then
               tmp_hcl(:) = qin(:,k,hcl_ndx)
            else
               tmp_hcl(:) = 0.
            end if
            if (clono2_ndx>0) then
               tmp_clono2(:) = qin(:,k,clono2_ndx)
            else
               tmp_clono2(:) = 0.
            end if
            if (hbr_ndx>0) then
               tmp_hbr(:) = qin(:,k,hbr_ndx)
            else
               tmp_hbr(:) = 0.
            end if
            call strat_chem_get_gamma( temp(:,k), pmid(:,k), m(:,k), &
                                       tmp_hcl, tmp_clono2, &
                                       psc, k, gamma_strat )
            call strat_chem_get_hetrates( temp(:,k), tmp_hcl, tmp_hbr, h2ovmr(:,k), &
                                          m(:,k), psc, gamma_strat, k, delt, hetrates_strat )
            
            if( strat72_ndx > 0 ) then
               rxt(:,k,strat72_ndx) = hetrates_strat(:,1)
            end if
            if( strat73_ndx > 0 ) then
               rxt(:,k,strat73_ndx) = hetrates_strat(:,2)
            end if
            if( strat74_ndx > 0 ) then
               rxt(:,k,strat74_ndx) = hetrates_strat(:,3)
            end if
            if( strat75_ndx > 0 ) then
               rxt(:,k,strat75_ndx) = hetrates_strat(:,4)
            end if
            if( strat76_ndx > 0 ) then
               rxt(:,k,strat76_ndx) = hetrates_strat(:,5)
            end if
            if( strat77_ndx > 0 ) then
               rxt(:,k,strat77_ndx) = hetrates_strat(:,6)
            end if
            if( strat78_ndx > 0 ) then
               rxt(:,k,strat78_ndx) = hetrates_strat(:,7)
            end if
            if( strat79_ndx > 0 ) then
               rxt(:,k,strat79_ndx) = hetrates_strat(:,8)
            end if
            if( strat80_ndx > 0 ) then
               rxt(:,k,strat80_ndx) = hetrates_strat(:,9)
            end if

         end if

      end do

      end subroutine usrrxt

      subroutine rh_calc(pmid, temp, sh, rh)
              
        implicit none
        
        real, intent(in), dimension(:) :: pmid, temp, sh
        real, intent(out), dimension(:) :: rh
        
!-----------------------------------------------------------------------
!       Calculate RELATIVE humidity.
!       This is calculated according to the formula:
!
!       RH   = qv / (epsilon*esat/ [pfull  -  (1.-epsilon)*esat])
!
!       Where epsilon = Rdgas/RVgas = d622
!
!       and where 1- epsilon = d378
!
!       Note that rh does not have its proper value
!       until all of the following code has been executed.  That
!       is, rh is used to store intermediary results
!       in forming the full solution.
!-----------------------------------------------------------------------
        
!-----------------------------------------------------------------------
!calculate water saturated specific humidity
!-----------------------------------------------------------------------
        call compute_qs (temp, pmid, rh, q = sh)
        
!-----------------------------------------------------------------------
!calculate rh
!-----------------------------------------------------------------------
        rh(:)= sh(:) / rh(:)
        
      end subroutine rh_calc

      subroutine set_aerosol(r_,rh,airdensity,drymass, rd, re, sfc_area)
!----------------------------------------------------------------
!     set aerosol information for heteorogenous calculation,
!      including dry radius(rd), effective radius(re, with rh correction), and surface area  
!----------------------------------------------------------------
        implicit none
        real, intent(in)        :: r_(:,:)              !species, second dimension is species index
        real, intent(in)        :: rh(:),airdensity(:)  !relative humidity
        real, intent(out)       :: drymass(:,:)         !aerosol dry mass (g/cm3)
        real, intent(out)       :: rd(:,:), re(:,:), sfc_area(:,:)  !second dimension is aerosol index
 
!-----------------------------------------------------------------
!     local parameter variables
!-----------------------------------------------------------------

        integer,  parameter :: nso4=26, nbcphob = 1, nbcphil=18, nocphob =1, nocphil=21, nssalt=130, ndust=8
        integer,  parameter :: nso4_rh = 26, nssalt_rh = 26
        integer,  parameter :: so4_rh1(nso4_rh)=(/30,35,40,45,50,55,60,65,70,75,80,82,84,86,88,90,91,92,93,94,95,96,97,98,99,100/)
        integer,  parameter :: bcphil_rh1(nbcphil)=(/0,70,75,80,82,84,86,88,90,91,92,93,94,95,96,97,98,99/)
        integer,  parameter :: ocphil_rh1(nocphil)=(/0,55,60,65,70,75,80,82,84,86,88,90,91,92,93,94,95,96,97,98,99/)
        integer,  parameter :: ssalt_rh1(nssalt_rh)=(/0,47,50,55,60,65,70,72,74,76,78,80,82,84,86,88,90,91,92,93,94,95,96,97,98,99/)
!these densities have unit of g/cm3.
        real,  parameter    :: denso4 = 1.74, denbc = 1.0, denocphob = 1.0, denocphil = 1.5 
        real,  parameter    :: denssalt(5) = (/2.2774594,2.2774594,2.2774594, 2.2774594,    2.2774594/)
        real,  parameter    :: dendust(8) = (/2.5,     2.5,     2.5,     2.5,     2.6,    2.6,    2.6,     2.6  /)
        real,  parameter   ::  fractiondust(8) = (/ 0.01053,0.08421,0.25263,0.65263,1.,1.,1.,1. /)
!-----------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------
!        real, dimension(size(r,1),size(r,2))  ::      so4, &           ! ammonium sulfate (VMR)
!                                                      bc1, &           ! Hydrophobic black carbon (MMR)
!                                                      bc2, &           ! Hydrophilic black carbon (MMR)
!                                                      oc1, &           ! Hydrophobic organic carbon (MMR)
!                                                      oc2, &           ! Hydrophilic organic carbon (MMR)
!                                                      soa, &           ! Secondary organic aerosol (MMR)
!                                                      ssa1, &          ! sea salt 1 (MMR)
!                                                      ssa2, &          ! sea salt 2 (MMR)
!                                                      ssa3, &          ! sea salt 3 (MMR)
!                                                      ssa4, &          ! sea salt 4 (MMR)
!                                                      ssa5, &          ! sea salt 5 (MMR)
!                                                      dust11, &        ! dust1 (MMR)
!                                                      dust12, &        ! dust1 (MMR)
!                                                      dust13, &        ! dust1 (MMR)
!                                                      dust14, &        ! dust1 (MMR)
!                                                      dust1, &         ! dust1 (MMR)
!                                                      dust2, &         ! dust2 (MMR) 
!                                                      dust3, &         ! dust3 (MMR)
!                                                      dust4, &         ! dust4 (MMR)
!                                                      dust5, &         ! dust5 (MMR)
        real, dimension(size(r_,1))        ::     rh_het
        integer, dimension(size(r_,1))     ::     irh
        integer, dimension(size(r_,1),size(r_,2))     ::     aeroindx!to save index
        real, parameter :: avo   = 6.023e23               ! molecules/mole
        integer   :: i, st1

        rh_het(:) = rh(:) *100.
        drymass(:,:)=0.
        sfc_area(:,:)=0.
        rd(:,:)=0.
        re(:,:)=0.
!----------------------------------------------------------------
!     SO4 
!----------------------------------------------------------------
!----------------------------------------------------------------
!     convert mmr to g/cm3 (drymass=aerop*m(:)/avo*28.97)
!=g/g*molecules/cm3 / (molecules/mol) *g/mol=g/cm3
!----------------------------------------------------------------   
!        drymass(:,1) = r_(:,so4_ndx)*132.*airdensity(:)/avo     !VMR => g/cm3
! here to take into account nitrate aerosols.
         drymass(:,1) = (r_(:,so4_ndx)*132.+r_(:,nh4no3_ndx)*80.) * &
                airdensity(:)/avo     !VMR => g/cm3
        irh(:)  = 0
        call find_indx(so4_rh1(:), rh_het(:), irh(:))
        aeroindx(:,1) = irh(:)
        rd(:,1)=RAA_HET(1)
        call give_indx(RAA_HET(:), aeroindx(:,1),re(:,1))
        call calc_sfc(re(:,1)*1.0D-4,rd(:,1)*1.0D-4,denso4,drymass(:,1),sfc_area(:,1))
!----------------------------------------------------------------
!     BC 
!----------------------------------------------------------------
!     --phobic,bc1
        st1=nso4+1
        aeroindx(:,2)=st1
        drymass(:,2) = r_(:,bc1_ndx)*airdensity(:)/avo*28.97 !MMR=> g/cm3
        rd(:,2)=RAA_HET(st1) !rd and re have unit of um.
        re(:,2)=RAA_HET(st1)
        call calc_sfc(re(:,2)*1.0D-4,rd(:,2)*1.0D-4,denbc,drymass(:,2),sfc_area(:,2))
!     --philic,bc2
        st1=st1+nbcphob
        irh(:)  = 0
        call find_indx(bcphil_rh1(:), rh_het(:), irh(:))
        aeroindx(:,3)= st1 + irh(:)
        drymass(:,3)=r_(:,bc2_ndx)*airdensity(:)/avo*28.97 !MMR=> g/cm3
        rd(:,3)=RAA_HET(st1)
        call give_indx(RAA_HET(:), aeroindx(:,3),re(:,3))
        call calc_sfc(re(:,3)*1.0D-4,rd(:,3)*1.0D-4,denbc,drymass(:,3),sfc_area(:,3))              
!----------------------------------------------------------------
!     OC(SOA is also taken into account)
!----------------------------------------------------------------
!     --phobic
        st1 = st1 + nbcphil
        aeroindx(:,4) = st1 
        drymass(:,4) = r_(:,oc1_ndx)*airdensity(:)/avo*28.97 !MMR=> g/cm3
        rd(:,4)=RAA_HET(st1)
        re(:,4)=RAA_HET(st1)
        call calc_sfc(re(:,4)*1.0D-4,rd(:,4)*1.0D-4,denocphob,drymass(:,4),sfc_area(:,4))
!     --philic      
        st1 = st1 + nocphob
        irh(:)  = 0
        call find_indx(so4_rh1(:), rh_het(:), irh(:))
        aeroindx(:,5)= st1 + irh(:)
        drymass(:,5)=(r_(:,oc2_ndx) + r_(:,soa_ndx))*airdensity(:)/avo*28.97 !MMR=> g/cm3
        rd(:,5)=RAA_HET(st1)
        call give_indx(RAA_HET(:), aeroindx(:,5),re(:,5))
        call calc_sfc(re(:,5)*1.0D-4,rd(:,5)*1.0D-4,denocphil,drymass(:,5),sfc_area(:,5))
!----------------------------------------------------------------
!     sea salt
!----------------------------------------------------------------
        st1 = st1 + nocphil
        irh(:)  = 0
        call find_indx(ssalt_rh1(:), rh_het(:), irh(:))      
        do i= 1, 5
           aeroindx(:,5+i) = st1 + (i-1)*nssalt_rh + irh(:)
           drymass(:,5+i)=r_(:,ssa_ndx(i))*airdensity(:)/avo*28.97 !MMR=> g/cm3 
           rd(:,5+i) = RAA_HET(st1 + (i-1)*nssalt_rh)
           call give_indx(RAA_HET(:), aeroindx(:,5+i),re(:,5+i))
           call calc_sfc(re(:,5+i)*1.0D-4,rd(:,5+i)*1.0D-4,denssalt(i),drymass(:,5+i),sfc_area(:,5+i))        
        end do
!----------------------------------------------------------------
!     dust, no hygroscopicity correction
!----------------------------------------------------------------
        st1 = st1 + nssalt -1
        do i= 1, 4
           aeroindx(:,10+i) = st1 + i
           drymass(:,10+i)=r_(:,dust_ndx(1)) * fractiondust(i)*airdensity(:)/avo*28.97
           rd(:,10+i)=RAA_HET(st1+i)
           re(:,10+i)=RAA_HET(st1+i)
           call calc_sfc(re(:,10+i)*1.0D-4,rd(:,10+i)*1.0D-4,dendust(i),drymass(:,10+i),sfc_area(:,10+i)) 
        end do
        do i= 5, 8
           aeroindx(:,10+i) = st1 + i
           drymass(:,10+i)=r_(:,dust_ndx(i-3))*airdensity(:)/avo*28.97 
           rd(:,10+i)=RAA_HET(st1+i)
           re(:,10+i)=RAA_HET(st1+i)
           call calc_sfc(re(:,10+i)*1.0D-4,rd(:,10+i)*1.0D-4,dendust(i),drymass(:,10+i),sfc_area(:,10+i))        
        end do

      end subroutine set_aerosol

      subroutine calc_sfc(aero_re, aero_rd, aero_density,aero_mass,aero_sfca)
        implicit none
        real*8, intent(in):: aero_re(:), aero_rd(:), aero_mass(:)
        real*8, intent(in):: aero_density !density for each aerosol
        !here re is effective radius and rd is the dry radius
        real*8, intent(out) :: aero_sfca(:)

        real*8, dimension(size(aero_re)) :: scaler, scalevol
! aero_sfca (cm2/cm3), aero_radius(cm),aero_density(g/cm3),aero_mass(g/cm3)     
        scaler(:)=aero_re(:)/aero_rd(:)
        scalevol(:)=scaler(:)**3
        aero_sfca=3*aero_mass(:)*scalevol(:)/aero_re(:)/aero_density
        
      end subroutine calc_sfc
      
      subroutine calc_hetrate(aero_area,aero_radius,denair,stkcf,stk,sqm,hetrate)
!******************************************************************************
!  This subroutine calculates the 1st-order loss rate of species on 
!  wet aerosol surface. This subroutine is from arsl1k.f in GEOS-Chem.
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) AREA   (REAL*8) : sfc area of wet aerosols/volume of air [cm2/cm3]
!  (2 ) RADIUS (REAL*8) : radius of wet aerosol [cm], order of 0.01-10 um;
!                           note that radius here is Rd, not Ro
!  (3 ) DENAIR (REAL*8) : Density of air [#/cm3]
!  (4 ) STKCF  (REAL*8) : Sticking coefficient [unitless], order of 0.1
!  (5 ) STK    (REAL*8) : Square root of temperature [K]
!  (6 ) SQM    (REAL*8) : Square root of molecular weight [g/mole]
!
!  References:
!  ============================================================================
!  The 1st-order loss rate on wet aerosol (Dentener's Thesis, p. 14)
!  is computed as:
!
!      ARSL1K [1/s] = area / [ radius/dfkg + 4./(stkcf * xmms) ]        
! 
!  where XMMS = Mean molecular speed [cm/s] = sqrt(8R*TK/pi/M) for Maxwell 
!        DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
!        R=8.314 kg m2/s2 mol K=83140000 g cm2/s2 mol K
!******************************************************************************
        IMPLICIT NONE
        ! Arguments
        REAL*8, INTENT(IN) :: STKCF, aero_AREA, aero_RADIUS, STK, SQM, DENAIR
        REAL*8, INTENT(OUT) :: hetrate
        ! Local variables
        REAL*8             :: DFKG

        !=================================================================
        ! ARSL1K begins here!
        !=================================================================
        if ( aero_AREA .le. 0d0 .or. aero_RADIUS .le. 1d-30 ) then

           ! Use default value if AREA is negative
           ! or if RADIUS is either zero or negative
           hetrate = 1.D-30

        else 

           ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
           DFKG  = 9.45D17/DENAIR * STK * SQRT(3.472D-2 + 1.D0/(SQM*SQM))

           ! Compute ARSL1K according to the formula listed above
           hetrate = aero_AREA / ( aero_RADIUS/DFKG + 2.749064E-4*SQM/(STKCF*STK) )

        end if

      end subroutine calc_hetrate

        
      subroutine find_indx(t,tt, it)
        !this is routine is to find the index of rh for the whole array.
        implicit none

        integer, intent(in)   :: t(:)         ! defined categories
        real, intent(in)      :: tt(:)        ! model simulated values
        integer, intent(out)  :: it(:)        ! index


        integer      :: N_   ! number of parameters
        integer      :: i
        real, dimension(size(t,1)-1) :: dt 


        N_ = size(t,1)

        do i = 1, N_-1
           dt(i) = (t(i+1)-t(i))/2.      
        end do

        where (tt(:) .le. t(1) + dt(1)  )
           it(:) = 1
        endwhere

        where (tt(:) .gt. t(N_-1) + dt(N_-1)  )
           it(:) = N_
        endwhere

        do i=1,N_-2
           where(tt(:) .gt. t(i) + dt(i) .and. tt(:) .le. t(i+1) + dt(i+1)) 
              it(:) = i+1
           endwhere
        end do

      end subroutine find_indx


      subroutine give_indx(raa,rindx,ree)
        implicit none

        real, intent(in) :: raa(:) !data read from am3_scat.dat
        integer, intent(in) :: rindx(:) !index
        real, intent(out):: ree(:)

        integer :: N_
        integer    :: i,tmp
        
        N_=size(rindx,1)
        do i= 1, N_
           tmp=rindx(i)
           ree(i)=raa(tmp)
        end do
      end subroutine give_indx
      
      end module mo_usrrxt_mod
