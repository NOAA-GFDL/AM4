 
module shortwave_types_mod

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!-------  version number --------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'

!---------------------------------------------------------------------
!-------  interfaces --------

public :: assignment(=)

interface assignment(=)
   module procedure sw_output_type_eq
end interface

!---------------------------------------------------------------------
!------- public derived-types ------
public :: sw_output_type

type sw_output_type
     real, dimension(:,:,:), pointer :: dfsw=>NULL(),   &
                                        ufsw=>NULL(),  &
                                        fsw=>NULL(),   &
                                        hsw=>NULL()   
     real, dimension(:,:,:), pointer :: dfswcf=>NULL(),   &
                                        ufswcf=>NULL(),&
                                        fswcf=>NULL(),  &
                                        hswcf=>NULL()
      real, dimension(:,:), pointer :: dfsw_vis_sfc=>NULL(),   &
                                       ufsw_vis_sfc=>NULL()
      real, dimension(:,:), pointer :: dfsw_dir_sfc=>NULL()
      real, dimension(:,:), pointer :: ufsw_dir_sfc=>NULL()
      real, dimension(:,:), pointer :: dfsw_dir_sfc_clr=>NULL()
      real, dimension(:,:), pointer :: dfsw_dif_sfc=>NULL(),   &
                                       ufsw_dif_sfc=>NULL()
      real, dimension(:,:), pointer :: dfsw_dif_sfc_clr=>NULL()
      real, dimension(:,:), pointer :: dfsw_vis_sfc_dir=>NULL()
      real, dimension(:,:), pointer :: ufsw_vis_sfc_dir=>NULL()
      real, dimension(:,:), pointer :: dfsw_vis_sfc_clr=>NULL()
      real, dimension(:,:), pointer :: dfsw_vis_sfc_dif=>NULL(),   &
                                       ufsw_vis_sfc_dif=>NULL()
      real, dimension(:,:,:), pointer :: bdy_flx=>NULL()
      real, dimension(:,:,:), pointer :: bdy_flx_clr=>NULL()

      contains
         procedure :: alloc => shortwave_output_alloc
         procedure :: dealloc => shortwave_output_dealloc
end type sw_output_type

!---------------------------------------------------------------------
!------- private data ------

logical :: module_is_initialized=.false.   ! module is initialized ?

!---------------------------------------------------------------------

CONTAINS

!#####################################################################
!
!                     PUBLIC SUBROUTINES
!
!#####################################################################

subroutine sw_output_type_eq(sw_output_out,sw_output_in)

   type(sw_output_type), intent(inout) :: sw_output_out
   type(sw_output_type), intent(in)    :: sw_output_in

!  Need to add error trap to catch unallocated sw_output_in
   sw_output_out%fsw              = sw_output_in%fsw
   sw_output_out%dfsw             = sw_output_in%dfsw
   sw_output_out%ufsw             = sw_output_in%ufsw
   sw_output_out%hsw              = sw_output_in%hsw
   sw_output_out%dfsw_dir_sfc     = sw_output_in%dfsw_dir_sfc
   sw_output_out%ufsw_dir_sfc     = sw_output_in%ufsw_dir_sfc
   sw_output_out%dfsw_dif_sfc     = sw_output_in%dfsw_dif_sfc
   sw_output_out%ufsw_dif_sfc     = sw_output_in%ufsw_dif_sfc
   sw_output_out%dfsw_vis_sfc     = sw_output_in%dfsw_vis_sfc
   sw_output_out%ufsw_vis_sfc     = sw_output_in%ufsw_vis_sfc
   sw_output_out%ufsw_vis_sfc_dir = sw_output_in%ufsw_vis_sfc_dir
   sw_output_out%dfsw_vis_sfc_dir = sw_output_in%dfsw_vis_sfc_dir
   sw_output_out%dfsw_vis_sfc_dif = sw_output_in%dfsw_vis_sfc_dif
   sw_output_out%ufsw_vis_sfc_dif = sw_output_in%ufsw_vis_sfc_dif
   sw_output_out%bdy_flx          = sw_output_in%bdy_flx
   if (ASSOCIATED(sw_output_in%fswcf))then
       sw_output_out%fswcf            = sw_output_in%fswcf
       sw_output_out%dfswcf           = sw_output_in%dfswcf
       sw_output_out%ufswcf           = sw_output_in%ufswcf
       sw_output_out%hswcf            = sw_output_in%hswcf
       sw_output_out%dfsw_dir_sfc_clr = sw_output_in%dfsw_dir_sfc_clr
       sw_output_out%dfsw_dif_sfc_clr = sw_output_in%dfsw_dif_sfc_clr
       sw_output_out%dfsw_vis_sfc_clr = sw_output_in%dfsw_vis_sfc_clr
       sw_output_out%bdy_flx_clr      = sw_output_in%bdy_flx_clr
   endif

end subroutine sw_output_type_eq

!#####################################################################
!
!                     PRIVATE SUBROUTINES ?????
!
!#####################################################################
!######################################################################
! <SUBROUTINE NAME="shortwave_output_alloc">
!  <OVERVIEW>
!   Code that allocates and initializes shortwave output variables
!  </OVERVIEW>
!  <DESCRIPTION>
!   Shortwave_output_alloc allocates and initializes the components
!   of the sw_output_type variable Sw_output, which is used to hold
!   output data from shortwave_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call shortwave_output_alloc (ix, jx, kx, Sw_output)
!  </TEMPLATE>
!  <IN NAME="ix" TYPE="integer">
!   x dimention of the radiation grid where shortwave output is desired
!  </IN>
!  <IN NAME="jx" TYPE="integer">
!   y dimention of the radiation grid where shortwave output is desired
!  </IN>
!  <IN NAME="kx" TYPE="integer">
!   z dimention of the radiation grid where shortwave output is desired
!  </IN>
!  <INOUT NAME="Sw_output" TYPE="sw_output_type">
!   shortwave radiation output variable
!  </INOUT>
! </SUBROUTINE>
!
subroutine shortwave_output_alloc (Sw_output, ix, jx, kx, &
                                   do_totcld_forcing)

!--------------------------------------------------------------------
!    shortwave_output_alloc allocates and initializes the components
!    of the sw_output_type variable Sw_output, which is used to hold
!    output data from shortwave_driver_mod.
!--------------------------------------------------------------------

class(sw_output_type), intent(inout) ::  Sw_output 
integer,              intent(in)     ::  ix, jx, kx
logical,              intent(in)     ::  do_totcld_forcing

!-------------------------------------------------------------------
!  intent(in) variables:
!
!    ix, jx, kx   dimensions of the radiation grid on which output 
!                 will be produced
!
!  intent(inout) variables:
!
!      Sw_output  sw_output_type variable containing shortwave 
!                 radiation output data 
!
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!    allocate and initialize fields to contain net(up-down) sw flux 
!    (fsw), upward sw flux (ufsw), downward sw flux(dfsw) at flux 
!    levels and sw heating in model layers (hsw).
!--------------------------------------------------------------------
      allocate (Sw_output%fsw  (ix, jx, kx+1) )
      allocate (Sw_output%ufsw (ix, jx, kx+1) )
      allocate (Sw_output%dfsw (ix, jx, kx+1) )
      allocate (Sw_output%hsw  (ix, jx, kx  ) )
      allocate (Sw_output%dfsw_dir_sfc (ix, jx) )
      allocate (Sw_output%ufsw_dir_sfc (ix, jx) )
      allocate (Sw_output%ufsw_dif_sfc (ix, jx) )
      allocate (Sw_output%dfsw_dif_sfc (ix, jx) )
      allocate (Sw_output%dfsw_vis_sfc (ix, jx  ) )
      allocate (Sw_output%ufsw_vis_sfc (ix, jx  ) )
      allocate (Sw_output%dfsw_vis_sfc_dir (ix, jx  ) )
      allocate (Sw_output%ufsw_vis_sfc_dir (ix, jx  ) )
      allocate (Sw_output%dfsw_vis_sfc_dif (ix, jx  ) )
      allocate (Sw_output%ufsw_vis_sfc_dif (ix, jx  ) )
      allocate (Sw_output%bdy_flx          (ix, jx, 4) )

      Sw_output%fsw   (:,:,:) = 0.0
      Sw_output%dfsw  (:,:,:) = 0.0
      Sw_output%ufsw  (:,:,:) = 0.0
      Sw_output%hsw   (:,:,:) = 0.0
      Sw_output%dfsw_dir_sfc = 0.0
      Sw_output%ufsw_dir_sfc = 0.0
      Sw_output%dfsw_dif_sfc  = 0.0
      Sw_output%ufsw_dif_sfc = 0.0
      Sw_output%dfsw_vis_sfc = 0.
      Sw_output%ufsw_vis_sfc = 0.
      Sw_output%dfsw_vis_sfc_dir = 0.
      Sw_output%ufsw_vis_sfc_dir = 0.
      Sw_output%dfsw_vis_sfc_dif = 0.
      Sw_output%ufsw_vis_sfc_dif = 0.
      Sw_output%bdy_flx(:,:,:) = 0.0       

!---------------------------------------------------------------------
!    if the cloud-free values are desired, allocate and initialize 
!    arrays for the fluxes and heating rate in the absence of clouds.
!----------------------------------------------------------------------
      if (do_totcld_forcing) then
        allocate (Sw_output%fswcf  (ix, jx, kx+1) )
        allocate (Sw_output%dfswcf (ix, jx, kx+1) )
        allocate (Sw_output%ufswcf (ix, jx, kx+1) )
        allocate (Sw_output%hswcf  (ix, jx, kx  ) )
        allocate (Sw_output%dfsw_dir_sfc_clr (ix, jx) )
        allocate (Sw_output%dfsw_dif_sfc_clr (ix, jx) )
        allocate (Sw_output%dfsw_vis_sfc_clr (ix, jx  ) )
        allocate (Sw_output%bdy_flx_clr      (ix, jx, 4) )

        Sw_output%fswcf (:,:,:) = 0.0
        Sw_output%dfswcf(:,:,:) = 0.0
        Sw_output%ufswcf(:,:,:) = 0.0
        Sw_output%hswcf (:,:,:) = 0.0
        Sw_output%dfsw_dir_sfc_clr = 0.0
        Sw_output%dfsw_dif_sfc_clr  = 0.0
        Sw_output%dfsw_vis_sfc_clr = 0.
        Sw_output%bdy_flx_clr (:,:,:) = 0.0
      endif

!--------------------------------------------------------------------

end subroutine shortwave_output_alloc

!###################################################################
! <SUBROUTINE NAME="shortwave_output_dealloc">
!  <OVERVIEW>
!   Code that allocates and initializes shortwave output variables
!  </OVERVIEW>
!  <DESCRIPTION>
!   Shortwave_output_dealloc allocates and initializes the components
!   of the sw_output_type variable Sw_output, which is used to hold
!   output data from shortwave_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call shortwave_output_dealloc (Sw_output)
!  </TEMPLATE>
!  <INOUT NAME="Sw_output" TYPE="sw_output_type">
!   shortwave radiation output variable
!  </INOUT>
! </SUBROUTINE>
!
subroutine shortwave_output_dealloc (Sw_output)

!--------------------------------------------------------------------
!    shortwave_output_dealloc deallocates the components
!    of the sw_output_type variable Sw_output, which is used to hold
!    output data from shortwave_driver_mod.
!--------------------------------------------------------------------

class(sw_output_type), intent(inout)  ::  Sw_output

!-------------------------------------------------------------------
!  intent(inout) variables:
!
!      Sw_output  sw_output_type variable containing shortwave 
!                 radiation output data 
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    deallocate fields to contain net(up-down) sw flux 
!    (fsw), upward sw flux (ufsw), downward sw flux(dfsw) at flux 
!    levels and sw heating in model layers (hsw).
!--------------------------------------------------------------------
      deallocate (Sw_output%fsw)
      deallocate (Sw_output%ufsw)
      deallocate (Sw_output%dfsw)
      deallocate (Sw_output%hsw)
      deallocate (Sw_output%dfsw_dir_sfc)
      deallocate (Sw_output%ufsw_dir_sfc)
      deallocate (Sw_output%ufsw_dif_sfc)
      deallocate (Sw_output%dfsw_dif_sfc)
      deallocate (Sw_output%dfsw_vis_sfc)
      deallocate (Sw_output%ufsw_vis_sfc)
      deallocate (Sw_output%dfsw_vis_sfc_dir)
      deallocate (Sw_output%ufsw_vis_sfc_dir)
      deallocate (Sw_output%dfsw_vis_sfc_dif)
      deallocate (Sw_output%ufsw_vis_sfc_dif)
      deallocate (Sw_output%bdy_flx)
!---------------------------------------------------------------------
!    if the cloud-free values are desired, allocate and initialize 
!    arrays for the fluxes and heating rate in the absence of clouds.
!----------------------------------------------------------------------
      if (ASSOCIATED(Sw_output%hswcf)) then
        deallocate (Sw_output%fswcf)
        deallocate (Sw_output%dfswcf)
        deallocate (Sw_output%ufswcf)
        deallocate (Sw_output%hswcf)
        deallocate (Sw_output%dfsw_dir_sfc_clr)
        deallocate (Sw_output%dfsw_dif_sfc_clr)
        deallocate (Sw_output%dfsw_vis_sfc_clr)
        deallocate (Sw_output%bdy_flx_clr)
      endif

!--------------------------------------------------------------------

end subroutine shortwave_output_dealloc

!####################################################################

end module shortwave_types_mod

