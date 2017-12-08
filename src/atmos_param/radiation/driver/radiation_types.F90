module radiation_types_mod
#include <fms_platform.h>
!---------------------------------------------------------------------

use mpp_mod, only: mpp_pe

!---- module data ----
use atmos_co2_mod,      only: atmos_co2_rad, co2_radiation_override
use block_control_mod,  only: block_control_type
use constants_mod,      only: WTMAIR, WTMCO2
use field_manager_mod,  only: MODEL_ATMOS
use mpp_domains_mod,    only: mpp_global_sum, domain2D, &
                              NON_BITWISE_EXACT_SUM, BITWISE_EXACT_SUM
use time_manager_mod,   only: time_type
use tracer_manager_mod, only: get_number_tracers, get_tracer_index, NO_TRACER

!---- public data ----
!---
!--- RADIATION STRUCTURE
!--- radiation control structure
 public radiation_input_control_type
 type radiation_input_control_type
     integer :: sphum                    ! index for specific humidity tracer
     logical :: phys_hydrostatic         ! hydrostratic flag for physics
     logical :: do_uni_zfull             ! miz
     type (domain2D) :: domain
 end type radiation_input_control_type

!--- atmosphere inputs block structure
 public radiation_input_block_type
 type radiation_input_block_type
     real, dimension(:,:),     _ALLOCATABLE :: phis   _NULL
     real, dimension(:,:,:),   pointer :: t      => null(), &
                                          pe     => null(), &
                                          peln   => null(), &
                                          delp   => null(), &
                                          delz   => null(), &
                                          p_full => null(), &
                                          p_half => null(), &
                                          z_full => null(), &
                                          z_half => null()
     real, dimension(:,:,:,:), pointer :: q      => null()
 end type radiation_input_block_type

!--- radiation global quantities (cannot be blocked)
 public radiation_input_glblqty_type
 type radiation_input_glblqty_type
     real                                :: atm_mass = -999. ! global mean atmospheric mass for CF/CMOR diagnostics
     real, dimension(:),         pointer :: gavg_q => null()  !_NULL
     real, dimension(:,:),  _ALLOCATABLE :: area   _NULL  ! atmospheric grid areas (for compute_g_avg)
     real, dimension(:,:),  _ALLOCATABLE :: pref   _NULL  ! reference pressure levels (for diagnostics)
 end type radiation_input_glblqty_type

!--- Radiation type definition
 public radiation_type
 type radiation_type
     type (radiation_input_control_type)                           :: control
     type (radiation_input_block_type), dimension(:), _ALLOCATABLE :: block _NULL
     type (radiation_input_glblqty_type)                           :: glbl_qty
 end type radiation_type

public :: compute_g_avg, alloc_radiation_type, dealloc_radiation_type

contains

 subroutine compute_g_avg(Time, tracer_name, Radiation, Atm_block, esm2_bugs)
    type(time_type),      intent(in)    :: Time
    character(len=*),     intent(in)    :: tracer_name
    type(radiation_type), intent(inout) :: Radiation
    type(block_control_type), intent(in) :: Atm_block
    logical, optional, intent(in) :: esm2_bugs
!------------------------------------------------------------
    real, dimension(Atm_block%iec-Atm_block%isc+1, &
                    Atm_block%jec-Atm_block%jsc+1) :: psfc_sum, qp_sum
    real :: qp, s1, bug_mult
    integer :: nb, ibs, ibe, jbs, jbe
    integer :: j, i, k, jb, ib, idx, npz
    logical :: esm2_bugs_local

!---the following is needed in order to allow ESM2 to reproduce a bug that
!---existed in the fv-latlon core (atmos_fv_dynamics::fv_physics.F90)
!---in which the dry mass mixing ratio was not used resulting in the following:
!--- qp = qp + q(co2)*(pe(k+1)-pe(k))
      bug_mult=1.0
      esm2_bugs_local=.false.
      if (PRESENT(esm2_bugs)) then
        if (esm2_bugs) then
          bug_mult=0.0
          esm2_bugs_local=.true.
        endif
      endif

    npz = Atm_block%npz

    Radiation%glbl_qty%gavg_q = 0.

!------------------------------------------------------------------------
!---compute global mean atmospheric mass to be used later in diagnostics
!------------------------------------------------------------------------
    psfc_sum = 0.
!$OMP parallel do default(shared) private(nb, ibs, ibe, jbs, jbe, j, i, jb, ib)
    do nb = 1, Atm_block%nblks
      ibs = Atm_block%ibs(nb) - Atm_block%isc + 1
      ibe = Atm_block%ibe(nb) - Atm_block%isc + 1
      jbs = Atm_block%jbs(nb) - Atm_block%jsc + 1
      jbe = Atm_block%jbe(nb) - Atm_block%jsc + 1

      do j = jbs, jbe
        jb = j-jbs+1
        do i = ibs, ibe
          ib = i-ibs+1
          psfc_sum(i,j) = Radiation%block(nb)%pe(ib,npz+1,jb)*Radiation%glbl_qty%area(i,j)
        enddo
      enddo
    enddo
    if (esm2_bugs_local) then
      Radiation%glbl_qty%atm_mass = mpp_global_sum(Radiation%control%domain, psfc_sum, &
                                      flags=BITWISE_EXACT_SUM)
    else
      Radiation%glbl_qty%atm_mass = REAL(mpp_global_sum(Radiation%control%domain, psfc_sum, &
                                      flags=NON_BITWISE_EXACT_SUM),KIND=4)
    endif

!------------------------------------------------------------------------
!---check to override predicted global pressure-weighted rad co2
!------------------------------------------------------------------------
    idx = get_tracer_index(MODEL_ATMOS, trim(tracer_name))
    if(idx /= NO_TRACER .and. co2_radiation_override) then
      call atmos_co2_rad(Time, Radiation%glbl_qty%gavg_q(idx))
    elseif (idx /= NO_TRACER) then

      npz = Atm_block%npz
      qp_sum = 0.
!$OMP parallel do default(shared) private(ibs, ibe, jbs, jbe, i, j, qp)
      do nb = 1, Atm_block%nblks
        ibs = Atm_block%ibs(nb) - Atm_block%isc + 1
        ibe = Atm_block%ibe(nb) - Atm_block%isc + 1
        jbs = Atm_block%jbs(nb) - Atm_block%jsc + 1
        jbe = Atm_block%jbe(nb) - Atm_block%jsc + 1

        do j = jbs, jbe
          jb = j-jbs+1
          do i = ibs, ibe
            ib = i-ibs+1
!---------------------------------------------------------------------
!  define pressure-weighted column mean value of dry mass mixing
!  ratio for tracer idx. assumption is that the tracer field q
!  is a moist mass mixing ratio. convert to dry mass mixing ratio by
!  dividing by (1 - qh2o).
!---------------------------------------------------------------------
            qp = 0.0
            do k=1,npz
              qp = qp + Radiation%block(nb)%q(ib,jb,k,idx)* &
                        (Radiation%block(nb)%pe(ib,k+1,jb)-Radiation%block(nb)%pe(ib,k,jb)) &
                        /(1.-bug_mult*Radiation%block(nb)%q(ib,jb,k,Radiation%control%sphum))
            enddo
            qp_sum(i,j) = qp * Radiation%glbl_qty%area(i,j)
           enddo
         enddo
       enddo
       if (esm2_bugs_local) then
         s1 = mpp_global_sum(Radiation%control%domain, qp_sum, flags=BITWISE_EXACT_SUM)
       else
         s1 = REAL(mpp_global_sum(Radiation%control%domain, qp_sum, flags=NON_BITWISE_EXACT_SUM),KIND=4)
       endif
       Radiation%glbl_qty%gavg_q(idx) = s1 / Radiation%glbl_qty%atm_mass
!---------------------------------------------------------------------
!    convert the tracer dry mass mixing ratio to the dry volume
!    mixing ratio.
!---------------------------------------------------------------------
       if (trim(tracer_name).eq.'co2') then
          Radiation%glbl_qty%gavg_q(idx) = Radiation%glbl_qty%gavg_q(idx)*WTMAIR/WTMCO2
       end if
    endif

 end subroutine compute_g_avg


 subroutine alloc_radiation_type (Radiation, Atm_block, p_hydro, do_uni_zfull) !miz
  type (radiation_type), intent(inout) :: Radiation
  type (block_control_type), intent(in) :: Atm_block
  logical,               intent(in)    :: p_hydro, do_uni_zfull !miz
!--- local varialbes
  integer :: n, ix, jx, npz, nt_prog

!---allocate/set control data
   npz = Atm_block%npz
   Radiation%control%phys_hydrostatic = p_hydro
   Radiation%control%do_uni_zfull = do_uni_zfull !miz

   call get_number_tracers(MODEL_ATMOS, num_prog=nt_prog)
   Radiation%control%sphum = get_tracer_index(MODEL_ATMOS, 'sphum')

!---allocate global_quantities
!--- allocate and set pref
   allocate (Radiation%glbl_qty%pref(npz+1,2))
   Radiation%glbl_qty%pref =0.
!--- allocate and set gavg_q
   allocate (Radiation%glbl_qty%gavg_q(nt_prog))
   Radiation%glbl_qty%gavg_q =0.
!--- allocate and set area
   ix = Atm_block%iec-Atm_block%isc+1
   jx = Atm_block%jec-Atm_block%jsc+1
   allocate (Radiation%glbl_qty%area(ix,jx))
   Radiation%glbl_qty%area = -999.

!---allocate input block data
   allocate (Radiation%block(Atm_block%nblks))
   do n = 1,Atm_block%nblks
     ix = Atm_block%ibe(n)-Atm_block%ibs(n)+1
     jx = Atm_block%jbe(n)-Atm_block%jbs(n)+1
     allocate (Radiation%block(n)%phis   (ix,jx),           &
               Radiation%block(n)%t      (ix,jx,npz),       &
               Radiation%block(n)%q      (ix,jx,npz,nt_prog), &
               Radiation%block(n)%pe     (ix,npz+1,jx),     &
               Radiation%block(n)%peln   (ix,npz+1,jx),     &
               Radiation%block(n)%delp   (ix,jx,npz),       &
               Radiation%block(n)%delz   (ix,jx,npz),       &
               Radiation%block(n)%p_full (ix,jx,npz),       &
               Radiation%block(n)%p_half (ix,jx,npz+1),     &
               Radiation%block(n)%z_full (ix,jx,npz),       &
               Radiation%block(n)%z_half (ix,jx,npz+1)      )
     Radiation%block(n)%phis   =0.
     Radiation%block(n)%t      =0.
     Radiation%block(n)%q      =0.
     Radiation%block(n)%pe     =0.
     Radiation%block(n)%peln   =0.
     Radiation%block(n)%delp   =0.
     Radiation%block(n)%delz   =0.
     Radiation%block(n)%p_full =0.
     Radiation%block(n)%p_half =0.
     Radiation%block(n)%z_full =0.
     Radiation%block(n)%z_half =0.
   enddo

 end subroutine alloc_radiation_type

 subroutine dealloc_radiation_type (Radiation)
   type (radiation_type), intent(inout) :: Radiation
!--- local variables
   integer :: n

!---deallocate global quantity
   deallocate (Radiation%glbl_qty%pref, &
               Radiation%glbl_qty%area, &
               Radiation%glbl_qty%gavg_q)
!---deallocate input block data
   do n = 1, size(Radiation%block,1)
     deallocate (Radiation%block(n)%phis,   &
                 Radiation%block(n)%t,      &
                 Radiation%block(n)%q,      &
                 Radiation%block(n)%pe,     &
                 Radiation%block(n)%peln,   &
                 Radiation%block(n)%delp,   &
                 Radiation%block(n)%delz,   &
                 Radiation%block(n)%p_full, &
                 Radiation%block(n)%p_half, &
                 Radiation%block(n)%z_full, &
                 Radiation%block(n)%z_half  )
   enddo
   deallocate (Radiation%block)

 end subroutine dealloc_radiation_type


end module radiation_types_mod
