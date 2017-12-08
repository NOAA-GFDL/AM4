module physics_types_mod
#include <fms_platform.h>
!---------------------------------------------------------------------
!---- module data ----
use block_control_mod,  only: block_control_type
use field_manager_mod,  only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_index, get_number_tracers

!---- public data ----
 public physics_control_type
 type physics_control_type
     integer :: sphum
     integer :: nsphum, nql, nqi, nqa, nqn, nqni, nqr, nqs, nqg
     integer :: num_prog_tracers
     logical :: hydrostatic, phys_hydrostatic, do_uni_zfull !miz
     logical :: l_host_applies_sfc_fluxes
     logical :: use_tau
     logical, dimension(:), _ALLOCATABLE :: cloud_tracer _NULL
 end type physics_control_type

!--- atmosphere inputs block structure
 public physics_input_block_type
 type physics_input_block_type
     real, dimension(:,:),     _ALLOCATABLE :: phis _NULL
     real, dimension(:,:,:,:), _ALLOCATABLE :: tmp_4d  _NULL
     real, dimension(:,:,:),   pointer :: u      => null(), &
                                          v      => null(), &
                                          w      => null(), &
                                          t      => null(), &
                                          omega  => null(), &
                                          pe     => null(), &
                                          peln   => null(), &
                                          delp   => null(), &
                                          delz   => null(), &
                                          p_full => null(), &
                                          p_half => null(), &
                                          z_full => null(), &
                                          z_half => null(), &
                                          ! used for leapfrog scheme
                                          um     => null(), &
                                          vm     => null(), &
                                          tm     => null()
     real, dimension(:,:,:,:), pointer :: q      => null(), &
                                          qm     => null()
 end type physics_input_block_type

!--- atmosphere global quantities (cannot be blocked)
 public physics_global_quantity_type
 type physics_global_quantity_type
     real,    dimension(:,:),  _ALLOCATABLE :: pref _NULL
 end type physics_global_quantity_type

!--- tendency block structure
 public physics_tendency_block_type
 type physics_tendency_block_type
     real, dimension(:,:,:),   pointer :: u_dt => null(), &
                                          v_dt => null(), &
                                          t_dt => null()
     real, dimension(:,:,:,:), pointer :: q_dt => null(), &
                                          qdiag => null()
 end type physics_tendency_block_type

!--- Physics type definition
 public physics_type
 type physics_type
     type (physics_control_type)                                 :: control
     type (physics_input_block_type), dimension(:), _ALLOCATABLE :: block _NULL
     type (physics_global_quantity_type)                         :: glbl_qty
 end type physics_type

!--- Physics tendency type definition
 public physics_tendency_type
 type physics_tendency_type
     type (physics_tendency_block_type), dimension(:), _ALLOCATABLE :: block _NULL
 end type physics_tendency_type

 public phys_mp_exch_type
 type phys_mp_exch_type
      real, dimension(:,:,:), pointer :: diff_t => null()
      real, dimension(:,:,:), pointer :: radturbten => null()
      real, dimension(:,:,:), pointer :: diff_cu_mo => null()
      real, dimension(:,:,:), pointer :: diff_t_clubb => null()
      real, dimension(:,:  ), pointer :: cush        => null()
      real, dimension(:,:  ), pointer :: cbmf        => null()
      real, dimension(:,:  ), pointer :: pbltop      => null()
   logical, dimension(:,:  ), pointer :: convect     => null()
      real, dimension(:,:  ), pointer :: tdt_shf     => null()
      real, dimension(:,:  ), pointer :: qdt_lhf     => null() 
      real, dimension(:,:  ), pointer :: hmint       => null() 
      real, dimension(:,:  ), pointer :: cgust       => null() 
      real, dimension(:,:  ), pointer :: tke         => null() 
      real, dimension(:,:  ), pointer :: pblhto      => null() 
      real, dimension(:,:  ), pointer :: rkmo        => null() 
      real, dimension(:,:  ), pointer :: taudpo     => null() 
   integer, dimension(:,:,:), pointer :: exist_shconv=> null() 
   integer, dimension(:,:,:), pointer :: exist_dpconv=> null() 
      real, dimension(:,:,:), pointer :: pblht_prev => null() 
      real, dimension(:,:,:), pointer :: hlsrc_prev => null() 
      real, dimension(:,:,:), pointer :: qtsrc_prev => null() 
      real, dimension(:,:,:), pointer ::  cape_prev => null() 
      real, dimension(:,:,:), pointer ::   cin_prev => null() 
      real, dimension(:,:,:), pointer ::   tke_prev => null() 
 end type phys_mp_exch_type

public Phys2cosp_type
type Phys2cosp_type 
   real, dimension(:,:,:), pointer :: temp_last=>NULL()
   real, dimension(:,:,:), pointer :: q_last=>NULL()
   real, dimension(:,:,:), pointer :: p_full=>NULL()
   real, dimension(:,:,:), pointer :: p_half=>NULL()
   real, dimension(:,:,:), pointer :: z_full=>NULL()
   real, dimension(:,:,:), pointer :: z_half=>NULL()
   real, dimension(:,:,:), pointer ::      u=>NULL()
   real, dimension(:,:,:), pointer ::      v=>NULL()
   real, dimension(:,:  ), pointer :: frac_land=>NULL()
   real, dimension(:,:  ), pointer :: lat=>NULL()
   real, dimension(:,:  ), pointer :: lon=>NULL()
end type Phys2cosp_type 

public precip_flux_type
type precip_flux_type
   real, dimension(:,:,:),        pointer :: fl_lsrain=>NULL()
   real, dimension(:,:,:),        pointer :: fl_lssnow=>NULL()
   real, dimension(:,:,:),        pointer :: fl_lsgrpl=>NULL()
   real, dimension(:,:,:),        pointer :: fl_ccrain=>NULL()
   real, dimension(:,:,:),        pointer :: fl_ccsnow=>NULL()
   real, dimension(:,:,:),        pointer :: fl_donmca_rain=>NULL()
   real, dimension(:,:,:),        pointer :: fl_donmca_snow=>NULL()
end type precip_flux_type



public :: alloc_physics_type, dealloc_physics_type, &
          alloc_physics_tendency_type, dealloc_physics_tendency_type

contains

 subroutine alloc_physics_type (Physics, Atm_block, p_hydro, hydro, do_uni_zfull) !miz
  type (physics_type), intent(inout) :: Physics
  type (block_control_type), intent(in) :: Atm_block
  logical,               intent(in) :: p_hydro, hydro, do_uni_zfull !miz
!--- local varialbes
  integer :: n, ix, jx, npz, nt_tot, nt_prog

!---control data
   npz = Atm_block%npz
   Physics%control%phys_hydrostatic = p_hydro
   Physics%control%hydrostatic = hydro
   Physics%control%do_uni_zfull = do_uni_zfull !miz

   call get_number_tracers(MODEL_ATMOS, num_tracers=nt_tot, num_prog=nt_prog)
   Physics%control%sphum = get_tracer_index(MODEL_ATMOS, 'sphum')

!---allocate global quantities
!--- set pref
   allocate (Physics%glbl_qty%pref(npz+1,2))
   Physics%glbl_qty%pref = 0.

!---allocate input block data
   allocate (Physics%block(Atm_block%nblks))
   do n = 1, Atm_block%nblks
     ix = Atm_block%ibe(n)-Atm_block%ibs(n)+1
     jx = Atm_block%jbe(n)-Atm_block%jbs(n)+1
     allocate (Physics%block(n)%phis   (ix,jx),             &
               Physics%block(n)%tmp_4d (ix,jx,npz,nt_prog+1:nt_tot), &
               Physics%block(n)%u      (ix,jx,npz),         &
               Physics%block(n)%v      (ix,jx,npz),         &
               Physics%block(n)%w      (ix,jx,npz),         &
               Physics%block(n)%t      (ix,jx,npz),         &
               Physics%block(n)%q      (ix,jx,npz,nt_prog), &
               Physics%block(n)%omega  (ix,jx,npz),         &
               Physics%block(n)%pe     (ix,npz+1,jx),       &
               Physics%block(n)%peln   (ix,npz+1,jx),       &
               Physics%block(n)%delp   (ix,jx,npz),         &
               Physics%block(n)%delz   (ix,jx,npz),         &
               Physics%block(n)%p_full (ix,jx,npz),         &
               Physics%block(n)%p_half (ix,jx,npz+1),       &
               Physics%block(n)%z_full (ix,jx,npz),         &
               Physics%block(n)%z_half (ix,jx,npz+1)        )

     Physics%block(n)%phis   = 0.
     Physics%block(n)%u      = 0.
     Physics%block(n)%v      = 0.
     Physics%block(n)%w      = 0.
     Physics%block(n)%t      = 0.
     Physics%block(n)%q      = 0.
     Physics%block(n)%omega  = 0.
     Physics%block(n)%pe     = 0.
     Physics%block(n)%peln   = 0.
     Physics%block(n)%delp   = 0.
     Physics%block(n)%delz   = 0.
     Physics%block(n)%p_full = 0.
     Physics%block(n)%p_half = 0.
     Physics%block(n)%z_full = 0.
     Physics%block(n)%z_half = 0.
   enddo

 end subroutine alloc_physics_type

 subroutine dealloc_physics_type (Physics)
   type (physics_type), intent(inout) :: Physics
!--- local variables
   integer :: n

!---deallocate global quantities
   deallocate (Physics%glbl_qty%pref)
!---deallocate input block data
   do n = 1, size(Physics%block,1)
     deallocate (Physics%block(n)%phis,   &
                 Physics%block(n)%u,      &
                 Physics%block(n)%v,      &
                 Physics%block(n)%w,      &
                 Physics%block(n)%t,      &
                 Physics%block(n)%q,      &
                 Physics%block(n)%omega,  &
                 Physics%block(n)%pe,     &
                 Physics%block(n)%peln,   &
                 Physics%block(n)%delp,   &
                 Physics%block(n)%delz,   &
                 Physics%block(n)%p_full, &
                 Physics%block(n)%p_half, &
                 Physics%block(n)%z_full, &
                 Physics%block(n)%z_half  )
     if (_ALLOCATED(Physics%block(n)%tmp_4d)) &
        deallocate(Physics%block(n)%tmp_4d)
   enddo
   deallocate (Physics%block)
 end subroutine dealloc_physics_type


 subroutine alloc_physics_tendency_type (Physics_tendency, Atm_block)
   type (physics_tendency_type), intent(inout) :: Physics_tendency
   type (block_control_type),    intent(in)    :: Atm_block
!--- local variables
   integer :: n, ix, jx, npz, nt_tot, nt_prog

    call get_number_tracers(MODEL_ATMOS, num_tracers=nt_tot, num_prog=nt_prog)
    allocate ( Physics_tendency%block(Atm_block%nblks) )
    npz = Atm_block%npz
    do n = 1,Atm_block%nblks
      ix = Atm_block%ibe(n) - Atm_block%ibs(n) + 1
      jx = Atm_block%jbe(n) - Atm_block%jbs(n) + 1
      allocate ( Physics_tendency%block(n)%u_dt(ix,jx,npz), &
                 Physics_tendency%block(n)%v_dt(ix,jx,npz), &
                 Physics_tendency%block(n)%t_dt(ix,jx,npz), &
                 Physics_tendency%block(n)%q_dt(ix,jx,npz,nt_prog), &
                 Physics_tendency%block(n)%qdiag(ix,jx,npz,nt_prog+1:nt_tot) )

      Physics_tendency%block(n)%u_dt = 0.
      Physics_tendency%block(n)%v_dt = 0.
      Physics_tendency%block(n)%t_dt = 0.
      Physics_tendency%block(n)%q_dt = 0.
      Physics_tendency%block(n)%qdiag = 0.
    enddo

 end subroutine alloc_physics_tendency_type


 subroutine dealloc_physics_tendency_type (Physics_tendency)
   type (physics_tendency_type), intent(inout) :: Physics_tendency
!--- local variables
   integer :: n

    do n = 1, size(Physics_tendency%block)
      deallocate (Physics_tendency%block(n)%u_dt, &
                  Physics_tendency%block(n)%v_dt, &
                  Physics_tendency%block(n)%t_dt, &
                  Physics_tendency%block(n)%q_dt, &
                  Physics_tendency%block(n)%qdiag)
    enddo
    deallocate (Physics_tendency%block)

 end subroutine dealloc_physics_tendency_type

end module physics_types_mod
