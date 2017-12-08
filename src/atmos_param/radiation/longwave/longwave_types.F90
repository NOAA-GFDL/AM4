 
module longwave_types_mod

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
   module procedure lw_output_type_eq
   module procedure lw_diagnostics_type_eq
end interface

!---------------------------------------------------------------------

public lw_diagnostics_type

type lw_diagnostics_type
     real, dimension(:,:),   pointer   :: flx1e1=>NULL(),  &
                                          gxcts=>NULL()
     real, dimension(:,:,:), pointer   :: flx1e1f=>NULL(),  &
                                          excts=>NULL(),&
                                          fctsg=>NULL()
     real, dimension(:,:,:,:), pointer :: fluxn=>NULL(),   &
                                          fluxncf=>NULL(),   &
                                          exctsn=>NULL(),  &
                                          cts_out=>NULL(), &
                                          cts_outcf=>NULL()

     contains
        procedure :: alloc => longwave_diag_alloc
        procedure :: dealloc => longwave_diag_dealloc
end type lw_diagnostics_type

!-------------------------------------------------------------------

public lw_output_type

type lw_output_type
     real, dimension(:,:,:), pointer :: heatra=>NULL(), &
                                        flxnet=>NULL(),  &
                                        heatracf=>NULL(), &
                                        flxnetcf=>NULL()
     real, dimension(:,:,:), pointer :: bdy_flx=>NULL(), &
                                        bdy_flx_clr=>NULL()

     contains
        procedure :: alloc => longwave_output_alloc
        procedure :: dealloc => longwave_output_dealloc
end type lw_output_type

!------------------------------------------------------------------

CONTAINS

!#####################################################################
!
!                     PUBLIC SUBROUTINES
!
!#####################################################################

! <SUBROUTINE NAME="longwave_output_alloc">
!  <OVERVIEW>
!   Subroutine to allocate output variables from longwave calculation
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine allocates and initializes the components
!    of the lw_output_type variable Lw_output which holds the longwave
!    output needed by radiation_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_output_alloc (Lw_output, ix, jx, kx)
!  </TEMPLATE>
!  <IN NAME="ix" TYPE="integer">
!   Dimension 1 length of radiation arrays to be allocated
!  </IN>
!  <IN NAME="jx" TYPE="integer">
!   Dimension 2 length of radiation arrays to be allocated
!  </IN>
!  <IN NAME="kx" TYPE="integer">
!   Dimension 3 length of radiation arrays to be allocated
!  </IN>
!  <OUT NAME="Lw_output" TYPE="lw_output_type">
!   lw_output_type variable containing longwave 
!                   radiation output data
!  </OUT>
! </SUBROUTINE>
!
subroutine longwave_output_alloc (Lw_output, ix, jx, kx,  &
                                  do_totcld_forcing)

!--------------------------------------------------------------------
!    longwave_driver_alloc allocates and initializes the components
!    of the lw_output_type variable Lw_output which holds the longwave
!    output needed by radiation_driver_mod.
!--------------------------------------------------------------------

class(lw_output_type),     intent(inout) :: Lw_output
integer,                   intent(in)    :: ix, jx, kx
logical,                   intent(in)    :: do_totcld_forcing

!--------------------------------------------------------------------
!   intent(in) variables:
!
!     ix,jx,kx      (i,j,k) dimensions of current physics window 
!     do_totcld_forcing   True if arrays for clear sky output
!                         are allocated
!
!   intent(inout) variables:
!
!     Lw_output     lw_output_type variable containing longwave 
!                   radiation output data 
!  
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    allocate and initialize arrays to hold net longwave fluxes and 
!    the longwave heating rate at each gridpoint. if the
!    cloud-forcing calculation is to be done, also allocate and init-
!    ialize arrays for fluxes and heating rates without clouds.
!-------------------------------------------------------------------
      allocate (Lw_output%flxnet ( ix, jx, kx+1) )
      allocate (Lw_output%heatra ( ix, jx, kx  ) )
      allocate (Lw_output%bdy_flx( ix, jx, 7) )    ! change by Xianglei Huang
      Lw_output%flxnet(:,:,:) = 0.0
      Lw_output%heatra(:,:,:) = 0.0
      Lw_output%bdy_flx (:,:,:) = 0.0      
      if (do_totcld_forcing)  then
        allocate (Lw_output%flxnetcf   ( ix, jx, kx+1) )
        allocate (Lw_output%heatracf   ( ix, jx, kx  ) )
        allocate (Lw_output%bdy_flx_clr( ix, jx, 7) )   ! change by Xianglei Huang
        Lw_output%flxnetcf(:,:,:) = 0.0
        Lw_output%heatracf(:,:,:) = 0.0
        Lw_output%bdy_flx_clr (:,:,:) = 0.0      
      endif
    
!--------------------------------------------------------------------

end subroutine longwave_output_alloc


!#####################################################################
! <SUBROUTINE NAME="longwave_output_dealloc">
!  <OVERVIEW>
!   Subroutine to deallocate output variables from longwave calculation
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine allocates and initializes the components
!    of the lw_output_type variable Lw_output which holds the longwave
!    output needed by radiation_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_output_alloc (Lw_output)
!  </TEMPLATE>
!  <OUT NAME="Lw_output" TYPE="lw_output_type">
!   lw_output_type variable containing longwave 
!                   radiation output data
!  </OUT>
! </SUBROUTINE>
!
subroutine longwave_output_dealloc (Lw_output)

!--------------------------------------------------------------------
!    longwave_output_dealloc deallocates the components
!    of the lw_output_type variable Lw_output.
!--------------------------------------------------------------------

class(lw_output_type),      intent(inout) :: Lw_output

!--------------------------------------------------------------------
!
!   intent(inout) variables:
!
!      Lw_output    lw_output_type variable containing longwave 
!                   radiation output data 
!  
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    deallocate arrays to hold net longwave fluxes and 
!    the longwave heating rate at each gridpoint
!-------------------------------------------------------------------
      deallocate (Lw_output%flxnet)
      deallocate (Lw_output%heatra)
      deallocate (Lw_output%bdy_flx)
      if (ASSOCIATED(Lw_output%flxnetcf))  then
        deallocate (Lw_output%flxnetcf)
        deallocate (Lw_output%heatracf)
        deallocate (Lw_output%bdy_flx_clr)
      endif
    
!--------------------------------------------------------------------

end subroutine longwave_output_dealloc

!#####################################################################
! <SUBROUTINE NAME="longwave_diag_dealloc">
!  
!   <OVERVIEW>
!     A routine to deallocate arrays allocated temporarily for
!     longwave diagnostics during radiation calculation.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine deallocates arrays used in longwave 
!     diagnostics and cloud space parameters used in the
!     lacis-hansen formulation.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call longwave_diag_dealloc (Lw_diagnostics)
!   </TEMPLATE>
!
!   <IN NAME="Lw_diagnostics" TYPE="lw_diagnostics_type">
!     Desired diagnostics from longwave_driver
!     so they may be passed to radiation_diag_mod
!   </IN>
! </SUBROUTINE>

subroutine longwave_diag_dealloc (Lw_diagnostics)

!---------------------------------------------------------------------
!    deallocate_arrays deallocates the array cpomponents of local
!    derived-type variables.
!---------------------------------------------------------------------

class(lw_diagnostics_type),  intent(inout)  :: Lw_diagnostics

!---------------------------------------------------------------------
!  intent(inout) variables:
!
!         Lw_diagnostics      lw_diagnostics_type variable to hold
!                             desired diagnostics from longwave_driver
!                             so they may be passed to 
!                             radiation_diag_mod
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!    deallocate the components of Lw_diagnostics.
!--------------------------------------------------------------------
      deallocate (Lw_diagnostics%flx1e1)
      deallocate (Lw_diagnostics%fluxn )
      deallocate (Lw_diagnostics%cts_out)
      deallocate (Lw_diagnostics%cts_outcf)
      deallocate (Lw_diagnostics%gxcts )
      deallocate (Lw_diagnostics%excts )
      deallocate (Lw_diagnostics%exctsn)
      deallocate (Lw_diagnostics%fctsg )
        deallocate (Lw_diagnostics%flx1e1f)
      if (ASSOCIATED(Lw_diagnostics%fluxncf)) then
        deallocate (Lw_diagnostics%fluxncf)
      endif

!--------------------------------------------------------------------

end subroutine longwave_diag_dealloc

!#####################################################################

! <SUBROUTINE NAME="longwave_diag_alloc">
!  <OVERVIEW>
!   Subroutine to allocate variables needed for longwave diagnostics
!  </OVERVIEW>
!  <TEMPLATE>
!   call longwave_diag_alloc (Lw_diagnostics, ix, jx, kx)
!  </TEMPLATE>
!  <IN NAME="ix" TYPE="integer">
!   Dimension 1 length of radiation arrays to be allocated
!  </IN>
!  <IN NAME="jx" TYPE="integer">
!   Dimension 2 length of radiation arrays to be allocated
!  </IN>
!  <IN NAME="kx" TYPE="integer">
!   Dimension 3 length of radiation arrays to be allocated
!  </IN>
!  <INOUT NAME="Lw_diagnostics" TYPE="lw_diagnostics_type">
!   lw_diagnostics_type variable containing longwave 
!                   radiation output data
!  </INOUT>
! </SUBROUTINE>
!
subroutine longwave_diag_alloc (Lw_diagnostics, ix, jx, kx, nbtrge, nbly, &
                                do_totcld_forcing)

!--------------------------------------------------------------------
!    longwave_diag_dealloc allocates and initializes the components of the 
!    lw_diagnostics_type variable Lw_diagnostics which holds diagnostic
!    output generated by sealw99_mod.  
!--------------------------------------------------------------------

class(lw_diagnostics_type), intent(inout) :: Lw_diagnostics
integer,                    intent(in)    :: ix, jx, kx, nbtrge, nbly
logical,                    intent(in)    :: do_totcld_forcing

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      ix,jx,kx     (i,j,k) lengths of radiation arrays to be allocated
!
!
!   intent(inout) variables:
!
!      Lw_diagnostics
!                   lw_diagnostics_type variable containing diagnostic
!                   longwave output used by the radiation diagnostics
!                   module
!  
!---------------------------------------------------------------------
!  module variable used (now arguments)
!
!   NBTRGE
!   NBLY
!
!--------------------------------------------------------------------
!    allocate (and initialize where necessary) lw_diagnostics_type 
!    component arrays.
!---------------------------------------------------------------------

      allocate ( Lw_diagnostics%flx1e1   (ix, jx                ) )
      allocate ( Lw_diagnostics%fluxn    (ix, jx, kx+1, 6+NBTRGE) )
      allocate (Lw_diagnostics%cts_out   (ix, jx, kx,   6       ) )
      allocate (Lw_diagnostics%cts_outcf (ix, jx, kx,   6       ) )
      allocate (Lw_diagnostics%gxcts     (ix, jx                ) )
      allocate (Lw_diagnostics%excts     (ix, jx, kx            ) )
      allocate (Lw_diagnostics%exctsn    (ix, jx, kx,   NBLY    ) )
      allocate (Lw_diagnostics%fctsg     (ix, jx,       NBLY    ) )

      Lw_diagnostics%flx1e1   = 0.
      Lw_diagnostics%cts_out    = 0.
      Lw_diagnostics%cts_outcf = 0.
      Lw_diagnostics%gxcts    = 0.
      Lw_diagnostics%excts  = 0.
      Lw_diagnostics%exctsn   = 0.
      Lw_diagnostics%fctsg   = 0.

      Lw_diagnostics%fluxn  (:,:,:,:) = 0.0

      if (do_totcld_forcing) then
        allocate ( Lw_diagnostics%fluxncf (ix, jx, kx+1, 6+NBTRGE) )
        Lw_diagnostics%fluxncf(:,:,:,:) = 0.0
      endif

      allocate( Lw_diagnostics%flx1e1f  (ix, jx,       NBTRGE  ) )
      Lw_diagnostics%flx1e1f  = 0.

!--------------------------------------------------------------------

end subroutine longwave_diag_alloc

!#####################################################################

subroutine lw_output_type_eq (Lw_output_out, Lw_output_in)

   type(lw_output_type), intent(inout) :: Lw_output_out
   type(lw_output_type), intent(in)    :: Lw_output_in

!  Need to add error trap to catch unallocated Lw_output_in
   Lw_output_out%heatra        = Lw_output_in%heatra
   Lw_output_out%flxnet        = Lw_output_in%flxnet
   Lw_output_out%bdy_flx       = Lw_output_in%bdy_flx
   if (ASSOCIATED(Lw_output_in%heatracf))then
       Lw_output_out%heatracf          = Lw_output_in%heatracf
       Lw_output_out%flxnetcf          = Lw_output_in%flxnetcf
       Lw_output_out%bdy_flx_clr       = Lw_output_in%bdy_flx_clr
   endif

end subroutine lw_output_type_eq

!##################################################################

subroutine lw_diagnostics_type_eq (Lw_diagnostics_out, Lw_diagnostics_in)

   type(lw_diagnostics_type), intent(inout) :: Lw_diagnostics_out
   type(lw_diagnostics_type), intent(in)    :: Lw_diagnostics_in

!  Need to add error trap to catch unallocated Lw_diagnostics_in
   Lw_diagnostics_out%flx1e1    = Lw_diagnostics_in%flx1e1
   Lw_diagnostics_out%cts_out   = Lw_diagnostics_in%cts_out
   Lw_diagnostics_out%gxcts     = Lw_diagnostics_in%gxcts
   Lw_diagnostics_out%excts     = Lw_diagnostics_in%excts
   Lw_diagnostics_out%exctsn    = Lw_diagnostics_in%exctsn
   Lw_diagnostics_out%fctsg     = Lw_diagnostics_in%fctsg
   Lw_diagnostics_out%flx1e1f   = Lw_diagnostics_in%flx1e1f
   Lw_diagnostics_out%fluxn     = Lw_diagnostics_in%fluxn
   Lw_diagnostics_out%cts_outcf = Lw_diagnostics_in%cts_outcf
   
   if (ASSOCIATED(Lw_diagnostics_in%fluxncf)) then
      Lw_diagnostics_out%fluxncf = Lw_diagnostics_in%fluxncf
   end if

end subroutine lw_diagnostics_type_eq

!####################################################################

end module longwave_types_mod

