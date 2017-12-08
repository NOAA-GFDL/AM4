module land_tile_selectors_mod

use fms_mod, only : error_mesg, WARNING
use land_data_mod, only : log_version

implicit none
private

! ==== public interface ======================================================
public :: tile_selector_type      !

! selector tags
public :: SEL_SOIL, SEL_VEGN, SEL_LAKE, SEL_GLAC, SEL_SNOW, SEL_CANA, SEL_HLSP

public :: tile_selectors_init     ! initialize module
public :: tile_selectors_end      ! clean up after ourselves

public :: register_tile_selector  ! register selector for diag field
public :: selector_suffix         ! return suffix for the field name

public :: selectors               ! array of selectors
public :: n_selectors             ! number of available selectors
! ==== end of public interface ===============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'land_tile_selectors_mod'
#include "../shared/version_variable.inc"

integer, parameter :: SEL_LEN           = 16  ! max length of the selector name
integer, parameter :: SEL_LONG_NAME_LEN = 128 ! max name of the selector long name
integer, parameter :: INIT_SEL_SIZE     = 1   ! initial size of the array of selectors

! tags for tile-specific diagnostic selectors
integer, parameter :: SEL_SOIL = 1
integer, parameter :: SEL_VEGN = 2
integer, parameter :: SEL_LAKE = 3
integer, parameter :: SEL_GLAC = 4
integer, parameter :: SEL_SNOW = 5
integer, parameter :: SEL_CANA = 6
integer, parameter :: SEL_HLSP = 7

! ==== derived types =========================================================
type :: tile_selector_type
   character(len=SEL_LEN)           :: name =''          ! name of the selector
   character(len=SEL_LONG_NAME_LEN) :: long_name = ''    ! long name of the selector
   logical :: is_default = .FALSE.     ! if true, returned suffix is set to empty string --
                                       ! this is to handle default areas for soil/glac/vegn diagnostics
   logical :: area_is_static = .FALSE. ! if true, area does not change in time (e.g. area of land, soil, or lakes)

   integer :: tag = 0 ! tag of the model component
   integer :: idata1=0, idata2=0 ! integer data
   integer :: rdata1=0, rdata2=0 ! real data
   integer :: area_id = -1 ! diag manger ID of the area field associated with this filter
end type tile_selector_type

! ==== module private data ===================================================
logical :: module_is_initialized = .false.

! ==== module public data ====================================================
! array of registered selectors: not protected because land_tile_diag sets
! is_default and area_id
type(tile_selector_type), pointer :: selectors(:) => NULL()
integer, protected :: n_selectors = 0 ! number of registered selectors

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! ============================================================================
subroutine tile_selectors_init()

  if (module_is_initialized) return

  module_is_initialized = .true.
  call log_version(version, module_name, &
  __FILE__)

  allocate (selectors(INIT_SEL_SIZE))
  n_selectors = 0 ! initialize number of registered selectors
  ! register couple of default selectors (for all tiles and for each tile)
end subroutine tile_selectors_init


! ============================================================================
subroutine tile_selectors_end()

  module_is_initialized = .false.
  deallocate(selectors)
  n_selectors = 0
end subroutine tile_selectors_end


! ============================================================================
! registers a selector to be used for diagnostic output
subroutine register_tile_selector( name, long_name, tag, idata1, idata2, rdata1, rdata2, area_depends_on_time )
  character(len=*), intent(in) :: name
  character(len=*), intent(in), optional :: long_name
  integer, intent(in), optional :: tag
  integer, intent(in), optional :: idata1, idata2
  real,    intent(in), optional :: rdata1, rdata2
  logical, intent(in), optional :: area_depends_on_time

  ! ---- local vars
  type(tile_selector_type), pointer :: new_selectors(:)
  integer :: i

  ! check for conflict of names -- presumably, if the selector was already
  ! registered, then it is an error to register it again
  do i = 1, n_selectors
     if (trim(name)==trim(selectors(i)%name)) then
        call error_mesg(module_name,'attempt to register selector "'&
             //trim(name)//'" which has already been registered',WARNING)
        return ! just skip it
     endif
  enddo

  ! allocate additional space for selectors if necessary
  if(n_selectors >= size(selectors)) then
     allocate(new_selectors(max(n_selectors*2,1)))
     new_selectors(1:n_selectors) = selectors(1:n_selectors)
     deallocate(selectors)
     selectors => new_selectors
  endif

  ! set up the selector values
  n_selectors = n_selectors + 1
  selectors(n_selectors)%name = name
  if (present(long_name)) &
       selectors(n_selectors)%long_name = long_name
  if (present(tag)) selectors(n_selectors)%tag = tag
  if (present(idata1)) selectors(n_selectors)%idata1 = idata1
  if (present(idata2)) selectors(n_selectors)%idata2 = idata2
  if (present(rdata1)) selectors(n_selectors)%rdata1 = rdata1
  if (present(rdata2)) selectors(n_selectors)%rdata2 = rdata2
  if (present(area_depends_on_time)) &
          selectors(n_selectors)%area_is_static = .not.area_depends_on_time
end subroutine register_tile_selector


! ============================================================================
! returns variable suffix for given selector
function selector_suffix(selector)
  character(len=SEL_LEN+1) :: selector_suffix
  type(tile_selector_type), intent(in) :: selector

  if(selector%is_default) then
     selector_suffix = ''
  else
     selector_suffix = '_'//trim(selector%name)
  endif
end function selector_suffix

end module land_tile_selectors_mod
