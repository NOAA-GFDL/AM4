
module radiative_gases_types_mod

!--------------------------------------------------------------------

use time_manager_mod, only: time_type

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------

!---------------------------------------------------------------------
!-------  interfaces --------

public :: assignment(=)

interface assignment(=)
  module procedure radiative_gases_type_eq
end interface

!---------------------------------------------------------------------
!------- public derived-types ------

public   radiative_gases_type

type radiative_gases_type
     real, dimension(:,:,:), pointer :: qo3=>NULL()
     real                            :: rrvch4, rrvn2o, rrvco2,    &   
                                        rrvf11, rrvf12, rrvf113,  &
                                        rrvf22, rf11air, rf12air,  &
                                        rf113air, rf22air, &
                                        co2_for_last_tf_calc,  &
                                        co2_tf_offset, &
                                        co2_for_next_tf_calc, &
                                        ch4_for_last_tf_calc,  &
                                        ch4_tf_offset, &
                                        ch4_for_next_tf_calc, &
                                        n2o_for_last_tf_calc,  &
                                        n2o_tf_offset, &
                                        n2o_for_next_tf_calc, &
                                        ch4_for_tf_calc, &
                                        n2o_for_tf_calc, &
                                        co2_for_tf_calc
     logical                         :: time_varying_co2,  &
                                        time_varying_f11, &
                                        time_varying_f12,  &
                                        time_varying_f113, &
                                        time_varying_f22,  &
                                        time_varying_ch4, &
                                        time_varying_n2o, &
                                        use_ch4_for_tf_calc, &
                                        use_n2o_for_tf_calc, &
                                        use_co2_for_tf_calc, &
                                        use_model_supplied_co2
     type(time_type)                 :: Co2_time, Ch4_time, N2o_time
     contains
        procedure :: alloc => radiative_gases_alloc
        procedure :: dealloc => radiative_gases_dealloc
end type radiative_gases_type

!---------------------------------------------------------------------

CONTAINS

!#####################################################################
!
!                     PUBLIC SUBROUTINES
!
!#####################################################################

subroutine radiative_gases_alloc (Rad_gases, ix, jx, kx)

class(radiative_gases_type), intent(inout) :: Rad_gases
integer,                     intent(in)    :: ix, jx, kx

!--------------------------------------------------------------------
!    allocate an array in a radiative_gases_type variable to hold the
!    model ozone field at the current time. call ozone_driver to define
!    this field for use in the radiation calculation.
!--------------------------------------------------------------------

      allocate (Rad_gases%qo3(ix,jx,kx))
      Rad_gases%qo3 = 0.

end subroutine radiative_gases_alloc

!####################################################################
! <SUBROUTINE NAME="radiative_gases_dealloc">
!  <OVERVIEW>
!    radiative_gases_end is the destructor for radiative_gases_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    radiative_gases_end is the destructor for radiative_gases_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiative_gases_dealloc (Rad_gases)
!  </TEMPLATE>
!  <INOUT NAME="Rad_gases" TYPE="radiative_gases_type">
!   radiative_gases_type variable containing the radi-
!              ative gas input fields needed by the radiation package
!  </INOUT>
! </SUBROUTINE>

subroutine radiative_gases_dealloc (Rad_gases)

!--------------------------------------------------------------------

class(radiative_gases_type), intent(inout)  :: Rad_gases

!---------------------------------------------------------------------
!  intent(inout) variable:
!
!     Rad_gases
!
!---------------------------------------------------------------------
!    deallocate the variables in Rad_gases.
!--------------------------------------------------------------------

      deallocate (Rad_gases%qo3)

!--------------------------------------------------------------------

end subroutine radiative_gases_dealloc 

!####################################################################

subroutine radiative_gases_type_eq (Rad_gases_out, Rad_gases_in)
type(radiative_gases_type), intent(inout) :: Rad_gases_out
type(radiative_gases_type), intent(in)    :: Rad_gases_in

!---------------------------------------------------------------------
!   copy all the data except ozone
!   that will allocated and assigned elsewhere
!---------------------------------------------------------------------

      Rad_gases_out%ch4_tf_offset = Rad_gases_in%ch4_tf_offset
      Rad_gases_out%n2o_tf_offset = Rad_gases_in%n2o_tf_offset
      Rad_gases_out%co2_tf_offset = Rad_gases_in%co2_tf_offset
      Rad_gases_out%ch4_for_next_tf_calc = Rad_gases_in%ch4_for_next_tf_calc
      Rad_gases_out%n2o_for_next_tf_calc = Rad_gases_in%n2o_for_next_tf_calc
      Rad_gases_out%co2_for_next_tf_calc = Rad_gases_in%co2_for_next_tf_calc
      Rad_gases_out%rrvch4  = Rad_gases_in%rrvch4
      Rad_gases_out%rrvn2o  = Rad_gases_in%rrvn2o
      Rad_gases_out%rrvf11  = Rad_gases_in%rrvf11
      Rad_gases_out%rrvf12  = Rad_gases_in%rrvf12
      Rad_gases_out%rrvf113 = Rad_gases_in%rrvf113
      Rad_gases_out%rrvf22  = Rad_gases_in%rrvf22
      Rad_gases_out%rrvco2  = Rad_gases_in%rrvco2
      Rad_gases_out%time_varying_co2  = Rad_gases_in%time_varying_co2
      Rad_gases_out%time_varying_ch4  = Rad_gases_in%time_varying_ch4
      Rad_gases_out%time_varying_n2o  = Rad_gases_in%time_varying_n2o
      Rad_gases_out%time_varying_f11  = Rad_gases_in%time_varying_f11
      Rad_gases_out%time_varying_f12  = Rad_gases_in%time_varying_f12
      Rad_gases_out%time_varying_f113 = Rad_gases_in%time_varying_f113
      Rad_gases_out%time_varying_f22  = Rad_gases_in%time_varying_f22
      Rad_gases_out%Co2_time = Rad_gases_in%Co2_time
      Rad_gases_out%Ch4_time = Rad_gases_in%Ch4_time
      Rad_gases_out%N2o_time = Rad_gases_in%N2o_time
      Rad_gases_out%use_model_supplied_co2 = Rad_gases_in%use_model_supplied_co2

      Rad_gases_out%co2_for_last_tf_calc = Rad_gases_in%co2_for_last_tf_calc
      Rad_gases_out%ch4_for_last_tf_calc = Rad_gases_in%ch4_for_last_tf_calc
      Rad_gases_out%n2o_for_last_tf_calc = Rad_gases_in%n2o_for_last_tf_calc

!---------------------------------------------------------------------

end subroutine radiative_gases_type_eq

!#####################################################################

end module radiative_gases_types_mod

