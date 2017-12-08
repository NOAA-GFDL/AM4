 
module aerosolrad_types_mod

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!------ public data structures ------
!--------------------------------------------------------------------

public aerosolrad_control_type

type aerosolrad_control_type
    logical :: do_aerosol
    logical :: do_swaerosol
    logical :: do_lwaerosol
    logical :: volcanic_sw_aerosols
    logical :: volcanic_lw_aerosols
    logical :: do_swaerosol_forcing
    logical :: do_lwaerosol_forcing
    integer :: indx_swaf
    integer :: indx_lwaf
    integer :: num_sw_aerosol_bands
    integer :: num_lw_aerosol_bands
end type aerosolrad_control_type

!--------------------------------------------------------------------

public   aerosolrad_diag_type

type aerosolrad_diag_type
     real, dimension(:,:,:,:,:), pointer  :: extopdep=>NULL(), &
                                             absopdep=>NULL(), &
                                             asymdep=>NULL()

     real, dimension(:,:,:,:), pointer  :: extopdep_vlcno=>NULL(), &
                                           absopdep_vlcno=>NULL(), &
                                           lw_extopdep_vlcno=>NULL(), &
                                           lw_absopdep_vlcno=>NULL()

     real, dimension(:,:,:), pointer  :: sw_heating_vlcno=>NULL()

     real, dimension(:,:,:,:),   pointer  :: sw_ext=>NULL(), &
                                             sw_ssa=>NULL(), &
                                             sw_asy=>NULL(), &
                                             lw_ext=>NULL(), &
                                             lw_ssa=>NULL(), &
                                             lw_asy=>NULL()

     contains
        procedure :: alloc=>aerosolrad_diag_alloc
        procedure :: dealloc=>aerosolrad_diag_dealloc
end type aerosolrad_diag_type

!---------------------------------------------------------------------

CONTAINS

!#####################################################################
!#####################################################################
!              constructor/destructor routines
!#####################################################################
!#####################################################################

subroutine aerosolrad_diag_alloc ( Aerosolrad_diags, &
                                   ix, jx, kx, nfld, &
                                   nfld_sw_ext, nfld_sw_ssa, nfld_sw_asy, &
                                   nfld_lw_ext, nfld_lw_ssa, nfld_lw_asy, &
                                   Aerosolrad_control)

class(aerosolrad_diag_type),    intent(out) :: Aerosolrad_diags
integer,                        intent(in)  :: ix, jx, kx, nfld 
integer,                        intent(in)  :: nfld_sw_ext, nfld_sw_ssa, nfld_sw_asy, &
                                               nfld_lw_ext, nfld_lw_ssa, nfld_lw_asy
type(aerosolrad_control_type),  intent(in)  :: Aerosolrad_control

!---------------------------------------------------------------------------
!
!  Number of volcanic aerosol fields
!     nfld_sw_ext   number of fields contained in supplemental sw_ext file
!     nfld_sw_ssa   number of fields contained in supplemental sw_ssa file
!     nfld_sw_asy   number of fields contained in supplemental sw_asy file
!     nfld_lw_ext   number of fields contained in supplemental lw_ext file
!     nfld_lw_ssa   number of fields contained in supplemental lw_ssa file
!     nfld_lw_asy   number of fields contained in supplemental lw_asy file
!
!---------------------------------------------------------------------------

      allocate (Aerosolrad_diags%extopdep (ix, jx, kx, nfld, 10)) 
      allocate (Aerosolrad_diags%absopdep (ix, jx, kx, nfld, 10)) 
      allocate (Aerosolrad_diags%asymdep  (ix, jx, kx, nfld, 10)) 
      Aerosolrad_diags%extopdep = 0.0
      Aerosolrad_diags%absopdep = 0.0
      Aerosolrad_diags%asymdep = 0.0

      if (Aerosolrad_control%volcanic_sw_aerosols) then 
        allocate (Aerosolrad_diags%extopdep_vlcno   (ix, jx, kx, 3))
        allocate (Aerosolrad_diags%absopdep_vlcno   (ix, jx, kx, 3))
        allocate (Aerosolrad_diags%sw_heating_vlcno (ix, jx, kx)) 
        Aerosolrad_diags%extopdep_vlcno = 0.0
        Aerosolrad_diags%absopdep_vlcno = 0.0
        Aerosolrad_diags%sw_heating_vlcno = 0.0
        allocate (Aerosolrad_diags%sw_ext (ix, jx, kx, nfld_sw_ext))
        allocate (Aerosolrad_diags%sw_ssa (ix, jx, kx, nfld_sw_ssa))
        allocate (Aerosolrad_diags%sw_asy (ix, jx, kx, nfld_sw_asy))
      endif

      if (Aerosolrad_control%volcanic_lw_aerosols) then 
        allocate (Aerosolrad_diags%lw_extopdep_vlcno (ix, jx, kx+1, 2))
        allocate (Aerosolrad_diags%lw_absopdep_vlcno (ix, jx, kx+1, 2))
        Aerosolrad_diags%lw_extopdep_vlcno = 0.0
        Aerosolrad_diags%lw_absopdep_vlcno = 0.0
        allocate (Aerosolrad_diags%lw_ext (ix, jx, kx, nfld_lw_ext))
        allocate (Aerosolrad_diags%lw_ssa (ix, jx, kx, nfld_lw_ssa))
        allocate (Aerosolrad_diags%lw_asy (ix, jx, kx, nfld_lw_asy))
      endif

!--------------------------------------------------------------------

end subroutine aerosolrad_diag_alloc

!#####################################################################

subroutine aerosolrad_diag_dealloc (Aerosolrad_diags)

class(aerosolrad_diag_type), intent(inout) :: Aerosolrad_diags
!--------------------------------------------------------------------

      deallocate (Aerosolrad_diags%extopdep)
      deallocate (Aerosolrad_diags%absopdep)
      deallocate (Aerosolrad_diags%asymdep)

  !BW if (Aerosolrad_control%volcanic_sw_aerosols) then 
      if (associated(Aerosolrad_diags%extopdep_vlcno)) then
        deallocate (Aerosolrad_diags%extopdep_vlcno)
        deallocate (Aerosolrad_diags%absopdep_vlcno)
        deallocate (Aerosolrad_diags%sw_heating_vlcno)
        deallocate (Aerosolrad_diags%sw_ext)
        deallocate (Aerosolrad_diags%sw_ssa)
        deallocate (Aerosolrad_diags%sw_asy)
      endif

  !BW if (Aerosolrad_control%volcanic_lw_aerosols) then 
      if (associated(Aerosolrad_diags%lw_extopdep_vlcno)) then
        deallocate (Aerosolrad_diags%lw_extopdep_vlcno)
        deallocate (Aerosolrad_diags%lw_absopdep_vlcno)
        deallocate (Aerosolrad_diags%lw_ext)
        deallocate (Aerosolrad_diags%lw_ssa)
        deallocate (Aerosolrad_diags%lw_asy)
      endif

!--------------------------------------------------------------------

end subroutine aerosolrad_diag_dealloc

!#####################################################################

end module aerosolrad_types_mod

