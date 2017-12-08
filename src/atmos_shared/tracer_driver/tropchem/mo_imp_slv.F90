
      module mo_imp_sol_mod

#ifndef AM3_CHEM
      use chem_mods_mod,        only : clscnt4
#else
      use AM3_chem_mods_mod,    only : clscnt4
#endif
      use constants_mod,    only : PI

      implicit none

      private
      public :: imp_slv_init, imp_sol

!     save

      integer, parameter ::  inst = 1, avrg = 2
      real, parameter    ::  rel_err      = 1.e-3
      real, parameter    ::  high_rel_err = 1.e-4
!-----------------------------------------------------------------------
!               newton-raphson iteration limits
!-----------------------------------------------------------------------
!     integer, parameter :: cut_limit    = 5
      integer, parameter :: cut_limit    = 8

      integer :: ox_ndx,o1d_ndx,h2o_ndx
      integer :: oh_ndx, ho2_ndx, ch3o2_ndx, po2_ndx, ch3co3_ndx, &
                 c2h5o2_ndx, isopo2_ndx, macro2_ndx, mco3_ndx, c3h7o2_ndx, &
                 ro2_ndx, xo2_ndx, no_ndx, no2_ndx, no3_ndx, n2o5_ndx, &
                 c2h4_ndx, c3h6_ndx, isop_ndx, mvk_ndx, c10h16_ndx, &
                 eo2_ndx, isopnbo2_ndx, iepoxoo_ndx, &
                 mobaoo_ndx, ino2_ndx, mvko2_ndx,&
                 macrno2_ndx, mao3_ndx, maopo2_ndx, ato2_ndx,&
                 isopnb_ndx, macr_ndx, moba_ndx
      integer :: op_ho2_ndx, op_mo2_ndx, op_ch3co3_ndx, op_po2_ndx, op_c2h5o2_ndx, &
                 op_isopo2_ndx, op_macro2_ndx, op_mco3_ndx, op_c3h7o2_ndx, op_ro2_ndx, &
                 op_xo2_ndx, op_isopnbo2_ndx, op_iepoxo2_ndx, op_ino2_ndx, op_mvko2_ndx,&
                 op_macrno2_ndx, op_mao3_ndx, op_maopo2_ndx, op_ato2_ndx, op_eo2_ndx
      integer :: ol_o1d_ndx, ol_oh_ndx, ol_ho2_ndx, ol_c3h6_ndx, ol_isop_ndx, &
                 ol_c2h4_ndx, ol_mvk_ndx, ol_macr_ndx, ol_c10h16_ndx, uoh_no2_ndx, &
                 n2o5h_ndx, no3h_ndx, no2h_ndx, ol_isopnb_ndx
      logical :: do_ox_pl = .true.
      integer :: verbose
      real    :: r2d

      type hst_pl
         integer  ::  cnt(2)
         logical  ::  do_hst(2)
      end type hst_pl

      real, private                      ::   small
      real, private                      ::   epsilon(clscnt4)
      type(hst_pl), private, allocatable ::   imp_hst_prod(:)
      type(hst_pl), private, allocatable ::   imp_hst_loss(:)
      logical, private, allocatable      ::   factor(:)

character(len=128), parameter :: version     = '$Id$'
character(len=128), parameter :: tagname     = '$Name$'
logical                       :: module_is_initialized = .false.

      contains

      subroutine imp_slv_init( verbose_in, retain_cm3_bugs )
!-----------------------------------------------------------------------
!        ... initialize the implict solver
!-----------------------------------------------------------------------

#ifndef AM3_CHEM
      use chem_mods_mod,  only : clscnt4, endrun, implicit
      use mo_grid_mod,    only : pcnstm1
#else
      use AM3_chem_mods_mod,  only : clscnt4, endrun, implicit
      use AM3_mo_grid_mod,    only : pcnstm1
#endif
      use mo_chem_utls_mod, only : get_spc_ndx, get_rxt_ndx

      implicit none

!-----------------------------------------------------------------------
!       ... Dummy arguments
!-----------------------------------------------------------------------
      integer,          intent(in) :: verbose_in
      logical,          intent(in) :: retain_cm3_bugs

!-----------------------------------------------------------------------
!        ... local variables
!-----------------------------------------------------------------------
      integer :: m, astat
      integer :: wrk(32)
      real    :: eps(pcnstm1)
      character(len=128) ::  msg

      allocate( factor(implicit%iter_max),stat=astat )
      if( astat /= 0 ) then
         write(msg,*) 'imp_slv_init: failed to allocate factor array; error = ',astat
         call endrun(msg)
      end if
      factor(:) = .true.
      eps(:)    = rel_err

      if (retain_cm3_bugs) then
        ox_ndx = get_spc_ndx( 'OX' )
      else
        ox_ndx = get_spc_ndx ('O3')
      endif
      o1d_ndx = get_spc_ndx('O1D')
      h2o_ndx = get_spc_ndx('H2O')
      if( ox_ndx > 0 ) then
         eps(ox_ndx) = high_rel_err
      else
         m = get_spc_ndx( 'O3' )
         if( m > 0 ) then
            eps(m) = high_rel_err
         end if
      end if
      m = get_spc_ndx( 'N' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'NO' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'NO2' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'NO3' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'HNO3' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'HO2NO2' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'N2O5' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'ClONO2' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'BrONO2' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'OH' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'HO2' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'Cl' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'ClO' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'Br' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'BrO' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if

      do m = 1,max(1,clscnt4)
         epsilon(m) = eps(implicit%clsmap(m))
      end do
#ifndef AM3_CHEM
     if( ox_ndx > 0 ) then
         op_ho2_ndx = get_rxt_ndx( 'op_ho2' )
         op_mo2_ndx = get_rxt_ndx( 'op_mo2' )
         op_eo2_ndx = get_rxt_ndx( 'op_eo2' )
         op_po2_ndx = get_rxt_ndx( 'op_po2' )
         op_ch3co3_ndx = get_rxt_ndx( 'op_ch3co3' )
         op_c2h5o2_ndx = get_rxt_ndx( 'op_c2h5o2' )
         op_c3h7o2_ndx = get_rxt_ndx( 'op_c3h7o2' )
         op_isopo2_ndx = get_rxt_ndx( 'op_isopo2' )
         op_isopnbo2_ndx = get_rxt_ndx( 'op_isopnbo2' )
         op_iepoxo2_ndx = get_rxt_ndx( 'op_iepoxo2' )
         op_ino2_ndx = get_rxt_ndx( 'op_ino2' )
         op_mvko2_ndx = get_rxt_ndx( 'op_mvko2' )
         op_macro2_ndx = get_rxt_ndx( 'op_macro2' )
         op_macrno2_ndx = get_rxt_ndx( 'op_macrno2' )
         op_mao3_ndx = get_rxt_ndx( 'op_mao3' )
         op_maopo2_ndx = get_rxt_ndx( 'op_maopo2' )
         op_ato2_ndx = get_rxt_ndx( 'op_ato2' )
         wrk(1:17) = (/ op_ho2_ndx, op_mo2_ndx, op_eo2_ndx, op_po2_ndx, op_ch3co3_ndx, &
                        op_c2h5o2_ndx, op_c3h7o2_ndx, op_isopo2_ndx, op_isopnbo2_ndx, op_iepoxo2_ndx, &
                        op_ino2_ndx, op_mvko2_ndx, op_macro2_ndx, &
                        op_macrno2_ndx, op_mao3_ndx, op_maopo2_ndx, op_ato2_ndx /)
         if( any( wrk(1:17) < 1 ) ) then
            do_ox_pl = .false.
         end if
         if( do_ox_pl ) then
            ol_o1d_ndx = get_rxt_ndx( 'ol_o1d' )
            ol_oh_ndx = get_rxt_ndx( 'ol_oh' )
            ol_ho2_ndx = get_rxt_ndx( 'ol_ho2' )
            ol_c2h4_ndx = get_rxt_ndx( 'ol_c2h4' )
            ol_c3h6_ndx = get_rxt_ndx( 'ol_c3h6' )
            ol_isop_ndx = get_rxt_ndx( 'ol_isop' )
            ol_isopnb_ndx = get_rxt_ndx( 'ol_isopnb' )
            ol_mvk_ndx = get_rxt_ndx( 'ol_mvk' )
            ol_macr_ndx = get_rxt_ndx( 'ol_macr' )
            ol_c10h16_ndx = get_rxt_ndx( 'ol_c10h16' )
            uoh_no2_ndx = get_rxt_ndx( 'uoh_no2' )
            n2o5h_ndx = get_rxt_ndx( 'n2o5h' )
            no3h_ndx = get_rxt_ndx( 'no3h' )
            no2h_ndx = get_rxt_ndx( 'no2h' )
            wrk(1:13) = (/ ol_o1d_ndx, ol_oh_ndx, ol_ho2_ndx, ol_c2h4_ndx, ol_c3h6_ndx, &
                           ol_isop_ndx, ol_isopnb_ndx, ol_mvk_ndx, ol_macr_ndx, &
                           uoh_no2_ndx, n2o5h_ndx, no3h_ndx, &
                           no2h_ndx /)
            if( any( wrk(1:13) < 1 ) ) then
               do_ox_pl = .false.
            end if
         end if

         if( do_ox_pl ) then
            oh_ndx = get_spc_ndx( 'OH' )
            ho2_ndx = get_spc_ndx( 'HO2' )
            ch3o2_ndx = get_spc_ndx( 'CH3O2' )
            eo2_ndx = get_spc_ndx( 'EO2' )
            po2_ndx = get_spc_ndx( 'PO2' )
            ch3co3_ndx = get_spc_ndx( 'CH3CO3' )
            c2h5o2_ndx = get_spc_ndx( 'C2H5O2' )
            c3h7o2_ndx = get_spc_ndx( 'C3H7O2' )
            isopo2_ndx = get_spc_ndx( 'ISOPO2' )
            isopnbo2_ndx = get_spc_ndx( 'ISOPNBO2' )
            iepoxoo_ndx = get_spc_ndx( 'IEPOXOO' )
            ino2_ndx = get_spc_ndx( 'INO2' )
            mvko2_ndx = get_spc_ndx( 'MVKO2' )
            macro2_ndx = get_spc_ndx( 'MACRO2' )
            macrno2_ndx = get_spc_ndx( 'MACRNO2' )
            mao3_ndx = get_spc_ndx( 'MAO3' )
            maopo2_ndx = get_spc_ndx( 'MAOPO2' )
            ato2_ndx = get_spc_ndx( 'ATO2' )
            no_ndx = get_spc_ndx( 'NO' )
            no2_ndx = get_spc_ndx( 'NO2' )
            no3_ndx = get_spc_ndx( 'NO3' )
            n2o5_ndx = get_spc_ndx( 'N2O5' )
            c2h4_ndx = get_spc_ndx( 'C2H4' )
            c3h6_ndx = get_spc_ndx( 'C3H6' )
            isop_ndx = get_spc_ndx( 'ISOP' )
            isopnb_ndx = get_spc_ndx( 'ISOPNB' )
            mvk_ndx = get_spc_ndx( 'MVK' )
            macr_ndx = get_spc_ndx( 'MACR' )
! Here I removed terpenes to test the code (jmao, 05/29/2013)
            c10h16_ndx = get_spc_ndx( 'C10H16' )
            wrk(1:29) = (/ oh_ndx, ho2_ndx, ch3o2_ndx, eo2_ndx, po2_ndx, ch3co3_ndx, &
                           c2h5o2_ndx, c3h7o2_ndx, isopo2_ndx, isopnbo2_ndx,&
                           iepoxoo_ndx, &
                           ino2_ndx, mvko2_ndx, macro2_ndx, macrno2_ndx, mao3_ndx, maopo2_ndx, ato2_ndx,& 
                           no_ndx, no2_ndx, no3_ndx, n2o5_ndx, &
                           c2h4_ndx, c3h6_ndx, isop_ndx, isopnb_ndx, mvk_ndx, macr_ndx, &
                           c10h16_ndx /)
            if( any( wrk(1:29) < 1 ) ) then
               do_ox_pl = .false.
            end if
         end if
      else
         do_ox_pl = .false.
      end if



#else
       if( ox_ndx > 0 ) then
         op_ho2_ndx = get_rxt_ndx( 'op_ho2' )
         op_mo2_ndx = get_rxt_ndx( 'op_mo2' )
         op_po2_ndx = get_rxt_ndx( 'op_po2' )
         op_ch3co3_ndx = get_rxt_ndx( 'op_ch3co3' )
         op_c2h5o2_ndx = get_rxt_ndx( 'op_c2h5o2' )
         op_isopo2_ndx = get_rxt_ndx( 'op_isopo2' )
         op_macro2_ndx = get_rxt_ndx( 'op_macro2' )
         op_mco3_ndx = get_rxt_ndx( 'op_mco3' )
         op_c3h7o2_ndx = get_rxt_ndx( 'op_c3h7o2' )
         op_ro2_ndx = get_rxt_ndx( 'op_ro2' )
         op_xo2_ndx = get_rxt_ndx( 'op_xo2' )
         wrk(1:11) = (/ op_ho2_ndx, op_mo2_ndx, op_po2_ndx, op_ch3co3_ndx, op_c2h5o2_ndx, &
                        op_isopo2_ndx, op_macro2_ndx, op_mco3_ndx, op_c3h7o2_ndx, op_ro2_ndx, op_xo2_ndx /)
         if( any( wrk(1:11) < 1 ) ) then
            do_ox_pl = .false.
         end if
         if( do_ox_pl ) then
            ol_o1d_ndx = get_rxt_ndx( 'ol_o1d' )
            ol_oh_ndx = get_rxt_ndx( 'ol_oh' )
            ol_ho2_ndx = get_rxt_ndx( 'ol_ho2' )
            ol_c3h6_ndx = get_rxt_ndx( 'ol_c3h6' )
            ol_isop_ndx = get_rxt_ndx( 'ol_isop' )
            ol_c2h4_ndx = get_rxt_ndx( 'ol_c2h4' )
            ol_mvk_ndx = get_rxt_ndx( 'ol_mvk' )
            ol_macr_ndx = get_rxt_ndx( 'ol_macr' )
            ol_c10h16_ndx = get_rxt_ndx( 'ol_terp' )
            uoh_no2_ndx = get_rxt_ndx( 'uoh_no2' )
            n2o5h_ndx = get_rxt_ndx( 'n2o5h' )
            no3h_ndx = get_rxt_ndx( 'no3h' )
            wrk(1:12) = (/ ol_o1d_ndx, ol_oh_ndx, ol_ho2_ndx, ol_c3h6_ndx, ol_isop_ndx, &
                           ol_c2h4_ndx, ol_mvk_ndx, ol_macr_ndx, ol_c10h16_ndx, uoh_no2_ndx, &
                           n2o5h_ndx, no3h_ndx /)
            if( any( wrk(1:12) < 1 ) ) then
               do_ox_pl = .false.
            end if
         end if

         if( do_ox_pl ) then
            oh_ndx = get_spc_ndx( 'OH' )
            ho2_ndx = get_spc_ndx( 'HO2' )
            ch3o2_ndx = get_spc_ndx( 'CH3O2' )
            po2_ndx = get_spc_ndx( 'PO2' )
            ch3co3_ndx = get_spc_ndx( 'CH3CO3' )
            c2h5o2_ndx = get_spc_ndx( 'C2H5O2' )
            macro2_ndx = get_spc_ndx( 'MACRO2' )
            mco3_ndx = get_spc_ndx( 'MCO3' )
            c3h7o2_ndx = get_spc_ndx( 'C3H7O2' )
            ro2_ndx = get_spc_ndx( 'RO2' )
            xo2_ndx = get_spc_ndx( 'XO2' )
            no_ndx = get_spc_ndx( 'NO' )
            no2_ndx = get_spc_ndx( 'NO2' )
            no3_ndx = get_spc_ndx( 'NO3' )
            n2o5_ndx = get_spc_ndx( 'N2O5' )
            c2h4_ndx = get_spc_ndx( 'C2H4' )
            c3h6_ndx = get_spc_ndx( 'C3H6' )
            isop_ndx = get_spc_ndx( 'ISOP' )
            isopo2_ndx = get_spc_ndx( 'ISOPO2' )
            mvk_ndx = get_spc_ndx( 'MVK' )
            c10h16_ndx = get_spc_ndx( 'C10H16' )
            wrk(1:21) = (/ oh_ndx, ho2_ndx, ch3o2_ndx, po2_ndx, ch3co3_ndx, &
                           c2h5o2_ndx, macro2_ndx, mco3_ndx, c3h7o2_ndx, ro2_ndx, &
                           xo2_ndx, no_ndx, no2_ndx, no3_ndx, n2o5_ndx, &
                           c2h4_ndx, c3h6_ndx, isop_ndx, isopo2_ndx, mvk_ndx, c10h16_ndx /)
            if( any( wrk(1:21) < 1 ) ) then
               do_ox_pl = .false.
            end if
         end if
      else
         do_ox_pl = .false.
      end if
#endif

      small = 1.e6 * tiny( small )
      r2d = 180./PI

      verbose = verbose_in

      end subroutine imp_slv_init

      subroutine imp_sol( base_sol, reaction_rates, &
                          het_rates, extfrc, &
                          nstep, delt, &
                          lat, lon, &
                          prod_out, loss_out, check_convergence, &
                          non_convergence, &
                          plonl, plnplv, &
                          prod_ox, loss_ox)
!-----------------------------------------------------------------------
!              ... imp_sol advances the volumetric mixing ratio
!           forward one time step via the fully implicit euler scheme.
!           this source is meant for small l1 cache machines such as
!           the intel pentium and itanium cpus
!-----------------------------------------------------------------------

      use m_tracname_mod,        only : tracnam
#ifndef AM3_CHEM
      use chem_mods_mod,         only : clscnt4, imp_nzcnt, clsze, rxntot, hetcnt, extcnt, implicit
      use mo_grid_mod,           only : pcnstm1
      use mo_indprd_mod,         only : indprd
      use mo_imp_prod_loss_mod,  only : imp_prod_loss
      use mo_imp_lin_matrix_mod, only : imp_linmat
      use mo_imp_nln_matrix_mod, only : imp_nlnmat
      use mo_imp_factor_mod,     only : imp_lu_fac
      use mo_imp_solve_mod,      only : imp_lu_slv
#else
      use AM3_chem_mods_mod,         only : clscnt4, imp_nzcnt, clsze, rxntot, hetcnt, extcnt, implicit
      use AM3_mo_grid_mod,           only : pcnstm1
      use AM3_mo_indprd_mod,         only : indprd
      use AM3_mo_imp_prod_loss_mod,  only : imp_prod_loss
      use AM3_mo_imp_lin_matrix_mod, only : imp_linmat
      use AM3_mo_imp_nln_matrix_mod, only : imp_nlnmat
      use AM3_mo_imp_factor_mod,     only : imp_lu_fac
      use AM3_mo_imp_solve_mod,      only : imp_lu_slv
#endif

      implicit none

!-----------------------------------------------------------------------
!             ... dummy args
!-----------------------------------------------------------------------
      integer, intent(in)    :: nstep                     ! time step index (zero based)
      real,    intent(in)    :: lat(:)                    ! latitude
      real,    intent(in)    :: lon(:)                    ! longitude
      integer, intent(in)    :: plonl                     ! longitude tile dimension
      integer, intent(in)    :: plnplv                    ! plonl*plev
      real,    intent(in)    :: delt                      ! time step (seconds)
      real,    intent(in)    :: reaction_rates(plnplv,rxntot)
      real,    intent(in)    :: het_rates(plnplv,max(1,hetcnt)), &
                                extfrc(plnplv,max(1,extcnt))
      logical, intent(in)    :: check_convergence         ! check for convergence ?
      real,    intent(out)   :: non_convergence(plnplv)   ! flag for implicit solver non-convergence (fraction)
      real,    intent(inout) :: base_sol(plnplv,pcnstm1)
      real,    intent(inout) :: prod_out(plnplv,pcnstm1),loss_out(plnplv,pcnstm1)
      real,    intent(inout) :: prod_ox(plnplv),loss_ox(plnplv)
!-----------------------------------------------------------------------
!             ... local variables
!-----------------------------------------------------------------------
      type hst_buff
         real, pointer, dimension(:,:) :: buff
      end type hst_buff

      integer ::   nr_iter, &
                   lev, &
                   indx, &
                   i, &
                   j, &
                   k, &
                   m, &
                   cut_cnt, stp_con_cnt
      real :: interval_done, dt, dti
      real :: max_delta(max(1,clscnt4))
      real, dimension(max(1,imp_nzcnt)) :: &
                   sys_jac, &
                   lin_jac
      real, dimension(max(1,clscnt4)) :: &
                   solution, &
                   forcing, &
                   iter_invariant, &
                   prod, &
                   loss
      real :: lrxt(max(1,rxntot))
      real :: lsol(max(1,pcnstm1))
      real :: lhet(max(1,hetcnt))
      real, dimension(plnplv,max(1,clscnt4)) :: &
                   ind_prd
      logical ::   convergence
      logical ::   frc_mask
      logical ::   converged(max(1,clscnt4))
      integer :: indx_old

      if( implicit%indprd_cnt > 0 ) then
!-----------------------------------------------------------------------
!        ... class independent forcing
!-----------------------------------------------------------------------
         call indprd( 4, ind_prd, base_sol, extfrc, reaction_rates )
      else
         do m = 1,max(1,clscnt4)
            ind_prd(:,m) = 0.
         end do
      end if

!-----------------------------------------------------------------------
!        ... Initialize production/loss diagnostics
!-----------------------------------------------------------------------
      do k = 1,clscnt4
         j = implicit%clsmap(k)
         prod_out(:,j) = 0.
         loss_out(:,j) = 0.
      end do
!for ox budget (jmao, 1/1/2011)
      prod_ox(:) = 0.
      loss_ox(:) = 0.
      non_convergence(:) = 0.


level_loop : &
!++lwh
!     do lev = 1,plev
      do lev = 1,plnplv/plonl
!--lwh
lon_tile_loop : &
         do i = 1,plonl
            indx = (lev - 1)*plonl + i
            indx_old = 0
!-----------------------------------------------------------------------
!        ... transfer from base to local work arrays
!-----------------------------------------------------------------------
            do m = 1,rxntot
               lrxt(m) = reaction_rates(indx,m)
            end do
            if( hetcnt > 0 ) then
               do m = 1,hetcnt
                  lhet(m) = het_rates(indx,m)
                end do
            end if
!-----------------------------------------------------------------------
!        ... time step loop
!-----------------------------------------------------------------------
            dt            = delt
            cut_cnt       = 0
            stp_con_cnt   = 0
            interval_done = 0.
time_step_loop : &
            do
               dti = 1. / dt
!-----------------------------------------------------------------------
!        ... transfer from base to local work arrays
!-----------------------------------------------------------------------
               do m = 1,pcnstm1
                  lsol(m) = base_sol(indx,m)
               end do
!-----------------------------------------------------------------------
!        ... transfer from base to class array
!-----------------------------------------------------------------------
               do k = 1,clscnt4
                  j = implicit%clsmap(k)
                  m = implicit%permute(k)
                  solution(m) = lsol(j)
               end do
!-----------------------------------------------------------------------
!        ... set the iteration invariant part of the function f(y)
!        ... if there is "independent" production put it in the forcing
!-----------------------------------------------------------------------
               if( implicit%indprd_cnt > 0 .or. extcnt > 0 ) then
                  do m = 1,clscnt4
                     iter_invariant(m) = dti * solution(m) + ind_prd(indx,m)
                  end do
               else
                  do m = 1,clscnt4
                     iter_invariant(m) = dti * solution(m)
                  end do
               end if
!-----------------------------------------------------------------------
!        ... the linear component
!-----------------------------------------------------------------------
               if( implicit%lin_rxt_cnt > 0 ) then
                  call imp_linmat( lin_jac, lsol, lrxt, lhet )
               else
                  do j = 1,clscnt4
                     m = implicit%diag_map(j)
                     lin_jac(m) = -dti
                  end do
               end if

!=======================================================================
!        the newton-raphson iteration for f(y) = 0
!=======================================================================
iter_loop : &
               do nr_iter = 1,implicit%iter_max
!-----------------------------------------------------------------------
!        ... the non-linear component
!-----------------------------------------------------------------------
                  if( factor(nr_iter) ) then
                     if( implicit%nln_rxt_cnt > 0 ) then
                        call imp_nlnmat( sys_jac, lsol, lrxt, lin_jac, dti )
                     else
                        do m = 1,imp_nzcnt
                           sys_jac(m) = lin_jac(m)
                        end do
                     end if
!-----------------------------------------------------------------------
!         ... factor the "system" matrix
!-----------------------------------------------------------------------
                     call imp_lu_fac( sys_jac )
                  end if
!-----------------------------------------------------------------------
!           ... form f(y)
!-----------------------------------------------------------------------
                  call imp_prod_loss( prod, loss, lsol, lrxt, lhet )
                  do m = 1,clscnt4
                     forcing(m) = solution(m)*dti - (iter_invariant(m) + prod(m) - loss(m))
                  end do
!-----------------------------------------------------------------------
!         ... solve for the mixing ratio at t(n+1)
!-----------------------------------------------------------------------
                  call imp_lu_slv( sys_jac, forcing )
                  do m = 1,clscnt4
                     solution(m) = solution(m) + forcing(m)
                  end do
!-----------------------------------------------------------------------
!            ... convergence measures
!-----------------------------------------------------------------------
                  if( nr_iter > 1 ) then
                     do k = 1,clscnt4
                        m = implicit%permute(k)
                        if( abs(solution(m)) > 1.e-40 ) then
                           max_delta(k) = abs( forcing(m)/solution(m) )
                        else
                           max_delta(k) = 0.
                        end if
                     end do
                  end if
!-----------------------------------------------------------------------
!           ... limit iterate
!-----------------------------------------------------------------------
                  where( solution(:) < 0. )
                     solution(:) = 0.
                     endwhere
!-----------------------------------------------------------------------
!           ... transfer latest solution back to work array
!-----------------------------------------------------------------------
                  do k = 1,clscnt4
                     j = implicit%clsmap(k)
                     m = implicit%permute(k)
                     lsol(j) = solution(m)
                  end do
!-----------------------------------------------------------------------
!            ... check for convergence
!-----------------------------------------------------------------------
                  if( nr_iter > 1 ) then
                     do k = 1,clscnt4
                        m = implicit%permute(k)
                        frc_mask = abs( forcing(m) ) > small
                        if( frc_mask ) then
                           converged(k) =  abs(forcing(m)) <= epsilon(k)*abs(solution(m))
                        else
                           converged(k) =  .true.
                        end if
                     end do
                     convergence = all( converged(:) )
                     if( convergence ) then
                        exit
                     end if
                  end if
               end do iter_loop

!-----------------------------------------------------------------------
!            ... check for newton-raphson convergence
!-----------------------------------------------------------------------
               if( .not. convergence ) then
!-----------------------------------------------------------------------
!           ... non-convergence
!-----------------------------------------------------------------------
!                  if( pdiags%imp_slv ) then
!                     write(*,'('' imp_sol: Time step '',1p,e21.13,'' Failed to converge'')') dt
!                  end if
                  stp_con_cnt = 0
                  if( cut_cnt < cut_limit ) then
                     cut_cnt = cut_cnt + 1
                     dt = .5 * dt
                     cycle
                  else
!                    write(*,'('' imp_sol: failed to converge @ (lon,lat,lev,dt) = '',3i5,1p,e21.13)') indx,lat,lev,dt
                     non_convergence(indx) = non_convergence(indx) + dt/delt
                     if (indx_old /= indx) then
                        if (verbose >= 3) then
                        write(*,105) lon(i)*r2d,lat(i)*r2d,lev,dt
 105                    format('imp_sol: failed to converge @ (lon,lat,lev,dt) = ', 2f8.2,i5,1p,f12.4)
                        end if
                        if (verbose >= 4) then
                           do m = 1,clscnt4
                              if( .not. converged(m)) then
                                 write(*,'(1x,a8,1x,1pe10.3)') &
                                    tracnam(implicit%clsmap(m)), max_delta(m)
                              end if
                           end do
                        end if
                     else
                        if (verbose >= 4) then
                        write(*,105) lon(i)*r2d,lat(i)*r2d,lev,dt
                        end if
                     end if
                     indx_old = indx
                  end if
               end if
!-----------------------------------------------------------------------
!           ... check for interval done
!-----------------------------------------------------------------------
               interval_done = interval_done + dt
               if( abs( delt - interval_done ) <= .0001 ) then
                  exit
               else
!-----------------------------------------------------------------------
!           ... transfer latest solution back to base array
!-----------------------------------------------------------------------
                  if( convergence ) then
                     stp_con_cnt = stp_con_cnt + 1
                  end if
                  if ( convergence .or. .not. check_convergence ) then
                     do m = 1,pcnstm1
                        base_sol(indx,m) = lsol(m)
                     end do
!++lwh
!-----------------------------------------------------------------------
!        ... Production/loss diagnostics
!-----------------------------------------------------------------------
                     do k = 1,clscnt4
                        j = implicit%clsmap(k)
                        m = implicit%permute(k)
                        prod_out(indx,j) = prod_out(indx,j) + prod(m) * dt/delt
                        loss_out(indx,j) = loss_out(indx,j) + loss(m) * dt/delt
                     end do
!--lwh
                     if( stp_con_cnt >= 2 ) then
                        dt = 2.*dt
                        stp_con_cnt = 0
                        cut_cnt = max( 0,cut_cnt-1 )
                     end if
                  endif
                  dt = min( dt,delt-interval_done )
               end if
            end do time_step_loop
!-----------------------------------------------------------------------
!           ... transfer latest solution back to base array
!-----------------------------------------------------------------------
            if ( convergence .or. .not. check_convergence ) then
               do k = 1,clscnt4
                  j = implicit%clsmap(k)
                  m = implicit%permute(k)
                  base_sol(indx,j) = solution(m)
!++lwh
!-----------------------------------------------------------------------
!        ... Production/loss diagnostics
!-----------------------------------------------------------------------
                  prod_out(indx,j) = prod_out(indx,j) + prod(m) * dt/delt &
                                                      + ind_prd(indx,m)
                  loss_out(indx,j) = loss_out(indx,j) + loss(m) * dt/delt
!--lwh
               end do
            end if
#ifndef AM3_CHEM
     if( do_ox_pl ) then
          prod_ox(indx) = &
           (reaction_rates(indx,op_ho2_ndx)*base_sol(indx,ho2_ndx) &
           + reaction_rates(indx,op_mo2_ndx) *base_sol(indx,ch3o2_ndx) &
           + reaction_rates(indx,op_eo2_ndx) *base_sol(indx,eo2_ndx) &
           + reaction_rates(indx,op_po2_ndx) *base_sol(indx,po2_ndx) &
           + reaction_rates(indx,op_ch3co3_ndx) *base_sol(indx,ch3co3_ndx) &
           + reaction_rates(indx,op_c2h5o2_ndx) *base_sol(indx,c2h5o2_ndx) &
           + reaction_rates(indx,op_c3h7o2_ndx) *base_sol(indx,c3h7o2_ndx) &
           + .9*reaction_rates(indx,op_isopo2_ndx) *base_sol(indx,isopo2_ndx) &
           + 1.6*reaction_rates(indx,op_isopnbo2_ndx) *base_sol(indx,isopnbo2_ndx) &
           + reaction_rates(indx,op_iepoxo2_ndx) *base_sol(indx,iepoxoo_ndx) &
           + 1.3*reaction_rates(indx,op_ino2_ndx) *base_sol(indx,ino2_ndx) &
           + .965*reaction_rates(indx,op_mvko2_ndx) *base_sol(indx,mvko2_ndx) &
           + .97*reaction_rates(indx,op_macro2_ndx) *base_sol(indx,macro2_ndx) &
           + 1.85*reaction_rates(indx,op_macrno2_ndx) *base_sol(indx,macrno2_ndx) &
           + reaction_rates(indx,op_mao3_ndx) *base_sol(indx,mao3_ndx) &
           + reaction_rates(indx,op_maopo2_ndx) *base_sol(indx,maopo2_ndx) &
           + .96*reaction_rates(indx,op_ato2_ndx) *base_sol(indx,ato2_ndx) ) * base_sol(indx,no_ndx)

          !-----------------------------------------------------------------------
          !         ... ozone destruction (only valid for the troposphere!)
          !             also include ox loss from no2+oh, n2o5+aerosol, no3+aerosol
          !-----------------------------------------------------------------------
          loss_ox(indx) = reaction_rates(indx,ol_oh_ndx) *base_sol(indx,oh_ndx) &
          + reaction_rates(indx,ol_ho2_ndx) *base_sol(indx,ho2_ndx) &
          + reaction_rates(indx,ol_c2h4_ndx) *base_sol(indx,c2h4_ndx) &
          + reaction_rates(indx,ol_c3h6_ndx) *base_sol(indx,c3h6_ndx) &
          + reaction_rates(indx,ol_isop_ndx) *base_sol(indx,isop_ndx) &
          + reaction_rates(indx,ol_isopnb_ndx) *base_sol(indx,isopnb_ndx) &
          + .235*reaction_rates(indx,ol_c10h16_ndx)*base_sol(indx,c10h16_ndx) &
          + reaction_rates(indx,ol_mvk_ndx)*base_sol(indx,mvk_ndx) &
          + reaction_rates(indx,ol_macr_ndx)*base_sol(indx,macr_ndx) &
          + (reaction_rates(indx,uoh_no2_ndx) * base_sol(indx,no2_ndx) * base_sol(indx,oh_ndx) &
             + 3. * reaction_rates(indx,n2o5h_ndx) * base_sol(indx,n2o5_ndx) &
             + 2. * reaction_rates(indx,no3h_ndx) * base_sol(indx,no3_ndx)  &
             + .5* reaction_rates(indx,no2h_ndx) * base_sol(indx,no2_ndx) ) &
            /max( base_sol(indx,ox_ndx),1.e-20 )
         ! here I sperate the o1d term because we need ozone implicitly from o1d.
            loss_ox(indx)=loss_ox(indx)*base_sol(indx,ox_ndx)&
          +reaction_rates(indx,ol_o1d_ndx)*base_sol(indx,o1d_ndx)*base_sol(indx,h2o_ndx)
    end if
#else

     if( do_ox_pl ) then
             prod_ox(indx) = &
              (reaction_rates(indx,op_ho2_ndx)*base_sol(indx,ho2_ndx) &
              + reaction_rates(indx,op_mo2_ndx) *base_sol(indx,ch3o2_ndx) &
              + reaction_rates(indx,op_po2_ndx) *base_sol(indx,po2_ndx) &
              + reaction_rates(indx,op_ch3co3_ndx) *base_sol(indx,ch3co3_ndx) &
              + reaction_rates(indx,op_c2h5o2_ndx) *base_sol(indx,c2h5o2_ndx) &
              + .88*reaction_rates(indx,op_isopo2_ndx)*base_sol(indx,isopo2_ndx) &
              + .985*reaction_rates(indx,op_macro2_ndx)*base_sol(indx,macro2_ndx) &
              + reaction_rates(indx,op_mco3_ndx)*base_sol(indx,mco3_ndx) &
              + reaction_rates(indx,op_c3h7o2_ndx)*base_sol(indx,c3h7o2_ndx) &
              + reaction_rates(indx,op_ro2_ndx)*base_sol(indx,ro2_ndx) &
              + reaction_rates(indx,op_xo2_ndx)*base_sol(indx,xo2_ndx)) * base_sol(indx,no_ndx)
             !-----------------------------------------------------------------------
             !         ... ozone destruction (only valid for the troposphere!)
             !             also include ox loss from no2+oh, n2o5+aerosol, no3+aerosol
             !-----------------------------------------------------------------------
             loss_ox(indx) = reaction_rates(indx,ol_oh_ndx) *base_sol(indx,oh_ndx) &
             + reaction_rates(indx,ol_ho2_ndx) *base_sol(indx,ho2_ndx) &
             + reaction_rates(indx,ol_c2h4_ndx) *base_sol(indx,c2h4_ndx) &
             + reaction_rates(indx,ol_c3h6_ndx) *base_sol(indx,c3h6_ndx) &
             + .9*reaction_rates(indx,ol_isop_ndx) *base_sol(indx,isop_ndx) &
             + .8*(reaction_rates(indx,ol_mvk_ndx)*base_sol(indx,mvk_ndx) &
                   + reaction_rates(indx,ol_macr_ndx)*base_sol(indx,macro2_ndx)) &
             + .235*reaction_rates(indx,ol_c10h16_ndx)*base_sol(indx,c10h16_ndx) &
             + (reaction_rates(indx,uoh_no2_ndx) * base_sol(indx,no2_ndx) * base_sol(indx,oh_ndx) &
                + 3. * reaction_rates(indx,n2o5h_ndx) * base_sol(indx,n2o5_ndx) &
                + 2. * reaction_rates(indx,no3h_ndx) * base_sol(indx,no3_ndx)) &
               /max( base_sol(indx,ox_ndx),1.e-20 )
             ! here I sperate the o1d term because we need ozone implicitly from o1d.
               loss_ox(indx)=loss_ox(indx)*base_sol(indx,ox_ndx)&
             +reaction_rates(indx,ol_o1d_ndx)*base_sol(indx,o1d_ndx)*base_sol(indx,h2o_ndx)
             prod_ox(indx)=max(prod_ox(indx),0.0)
             loss_ox(indx)=max(loss_ox(indx),0.0)
       end if

#endif
         end do lon_tile_loop
      end do level_loop


      end subroutine imp_sol

      end module mo_imp_sol_mod
