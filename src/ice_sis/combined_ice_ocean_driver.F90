module combined_ice_ocean_driver
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of the Sea Ice Simulator (version 1).           *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!----------------------------------------------------------------------------
! This module provides a common interface for jointly stepping SIS2 and MOM6,
! and will evolve as a platform for tightly integrating the ocean and sea ice
! models.  This particular version is a null version for interface
! compatibility with SIS2, and should not be called to drive SIS1.

use fms_mod, only : error_mesg, FATAL, WARNING
use time_manager_mod, only : time_type
use ice_model_mod,   only : ice_data_type
use ocean_model_mod, only : ocean_public_type, ocean_state_type, ice_ocean_boundary_type

implicit none ; private

public update_slow_ice_and_ocean, ice_ocean_driver_init, ice_ocean_driver_end

type, public :: ice_ocean_driver_type ; private
  logical :: CS_is_initialized = .false.
end type ice_ocean_driver_type

contains

!>   This subroutine initializes the combined ice ocean coupling control type.
subroutine ice_ocean_driver_init(CS, Time_init, Time_in)
  type(ice_ocean_driver_type), pointer       :: CS        !< The control structure for combined ice-ocean driver
  type(time_type),             intent(in)    :: Time_init !< The start time for the coupled model's calendar.
  type(time_type),             intent(in)    :: Time_in   !< The time at which to initialize the coupled model.

  if (associated(CS)) then
    call error_mesg("ice_ocean_driver_init", &
           "ice_ocean_driver_init called with an associated "// &
           "ice_ocean_driver_type. Model is already initialized.", WARNING)
    return
  endif

  call error_mesg("ice_ocean_driver_init", &
         "ice_ocean_driver_init called in SIS1.  The dummy version of "// &
         "ice_ocean_driver for SIS1 does not work and should not be used.", &
         FATAL)

  ! allocate(CS)

end subroutine ice_ocean_driver_init

!>   The subroutine update_slow_ice_and_ocean uses the forcing already stored in
!! the ice_data_type to advance both the sea-ice (and icebergs) and ocean states
!! for a time interval coupling_time_step.  The dummy version in SIS1 does nothing.
subroutine update_slow_ice_and_ocean(CS, Ice, Ocn, Ocean_sfc, Ice_ocean_boundary, &
                               time_start_update, coupling_time_step)
  type(ice_ocean_driver_type), &
                           pointer       :: CS   !< The control structure for this driver
  type(ice_data_type),     intent(inout) :: Ice  !< The publicly visible ice data type
  type(ocean_state_type),  pointer       :: Ocn  !< The internal ocean state and control structures
  type(ocean_public_type), intent(inout) :: Ocean_sfc !< The publicly visible ocean surface state type
  type(ice_ocean_boundary_type), &
                           intent(inout) :: Ice_ocean_boundary !< A structure containing the various forcing
                                                               !! fields going from the ice to the ocean
  type(time_type),         intent(in)    :: time_start_update  !< The time at the beginning of the update step
  type(time_type),         intent(in)    :: coupling_time_step !< The amount of time over which to advance
                                                               !! the ocean and ice

  call error_mesg("update_slow_ice_and_ocean", &
        "update_slow_ice_and_ocean called in SIS1.  The dummy version of "// &
        "update_slow_ice_and_ocean for SIS1 does not work and should not be used.", &
        FATAL)

end subroutine update_slow_ice_and_ocean


!>   The subroutine ice_ocean_driver_end terminates the model run, saving
!! the ocean and slow ice states in restart files and deallocating any data
!! associated with the ocean and slow ice.
subroutine ice_ocean_driver_end(CS, Ice, Ocean_sfc, Ocn, Time)
  type(ice_ocean_driver_type), pointer   :: CS   !< The control structure for combined ice-ocean driver
  type(ice_data_type),     intent(inout) :: Ice  !< The publicly visible ice data type
  type(ocean_state_type),  pointer       :: Ocn  !< The internal ocean state and control structures
  type(ocean_public_type), intent(inout) :: Ocean_sfc !< The publicly visible ocean surface state type
  type(time_type),         intent(in)    :: Time !< The model time, used for writing restarts

  if (associated(CS)) deallocate(CS)

end subroutine ice_ocean_driver_end

end module combined_ice_ocean_driver
