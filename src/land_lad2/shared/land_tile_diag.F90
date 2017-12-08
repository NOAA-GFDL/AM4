module land_tile_diag_mod

use mpp_mod,            only : mpp_sum
use mpp_efp_mod,        only : mpp_reproducing_sum
use mpp_domains_mod,    only : mpp_pass_ug_to_sg
use time_manager_mod,   only : time_type
use diag_axis_mod,      only : get_axis_length
use diag_manager_mod,   only : register_diag_field, register_static_field, &
     diag_field_add_attribute, diag_field_add_cell_measures, send_data
use diag_util_mod,      only : log_diag_field_info
use fms_mod,            only : error_mesg, string, FATAL

use land_debug_mod, only : set_current_point, check_var_range
use land_tile_selectors_mod, only : tile_selectors_init, tile_selectors_end, &
     tile_selector_type, register_tile_selector, selector_suffix, &
     n_selectors, selectors
use land_tile_mod,      only : land_tile_type, diag_buff_type, land_tile_list_type, &
     land_tile_enum_type, first_elmt, loop_over_tiles, &
     land_tile_map, tile_is_selected, fptr_i0, fptr_r0, fptr_r0i
use land_data_mod,      only : lnd, lnd_sg, log_version, land_data_type
use land_debug_mod,     only : check_var_range, set_current_point
use tile_diag_buff_mod, only : diag_buff_type, realloc_diag_buff
use land_debug_mod,     only : check_var_range, set_current_point

implicit none
private


! ==== public interface ======================================================
public :: tile_diag_init
public :: tile_diag_end

public :: set_default_diag_filter ! set the default filter for the consequent
          ! register_tiled_diag_field calls. The default filter will handle the
          ! fields that appear in diag_table without subsampling suffix.

public :: diag_buff_type

public :: register_tiled_area_fields
public :: register_tiled_diag_field
public :: register_tiled_static_field
public :: add_tiled_diag_field_alias
public :: add_tiled_static_field_alias

public :: send_tile_data
public :: send_tile_data_r0d_fptr, send_tile_data_r1d_fptr
public :: send_tile_data_i0d_fptr

public :: dump_tile_diag_fields
public :: send_global_land_diag

! codes of tile aggregation operations
public :: OP_AVERAGE, OP_SUM, OP_VAR, OP_STD

interface send_tile_data
   module procedure send_tile_data_0d
   module procedure send_tile_data_1d
   module procedure send_tile_data_0d_array
end interface

public :: get_area_id

! name of the table used for CMOR-compatible variables
character(*), public, parameter :: cmor_name='cmor_land'
real, public, parameter      :: cmor_mrsos_depth=0.1 ! depth of mrsos soil moisture averaging, m
! ==== end of public interface ===============================================


! ==== module constants ======================================================
character(len=*), parameter :: mod_name = 'land_tile_diag_mod'
#include "../shared/version_variable.inc"

integer, parameter :: INIT_FIELDS_SIZE     = 1     ! initial size of the fields array
integer, parameter :: BASE_TILED_FIELD_ID  = 65536 ! base value for tiled field
! ids, to distinguish them from regular diagnostic fields. All IDs of tiled
! (that is, registered by register_*tiled_field functions are larger than
! BASE_TILED_FIELD_ID)
integer, parameter :: MIN_DIAG_BUFFER_SIZE = 1     ! min size of the per-tile diagnostic buffer
! operations used for tile data aggregation
integer, parameter :: OP_AVERAGE = 0 ! weighted average of tile values
integer, parameter :: OP_SUM     = 1 ! sum of all tile values
integer, parameter :: OP_VAR     = 2 ! variance of tile values
integer, parameter :: OP_STD     = 3 ! standard deviation of tile values

integer, parameter :: FLD_STATIC    = 0
integer, parameter :: FLD_DYNAMIC   = 1
integer, parameter :: FLD_LIKE_AREA = 2

! ==== derived types =========================================================
type :: tiled_diag_field_type
   integer, pointer :: ids(:) => NULL()
   integer :: offset ! offset of the field data in the buffer
   integer :: size   ! size of the field data in the per-tile buffers
   integer :: op     ! aggregation operation
   integer :: static ! static/dynamic indicator, one of FLD_STATIC/FLD_DYNAMIC/FLD_LIKE_AREA
   logical :: fill_missing ! if TRUE, missing values (e.g. ocean points) are
                     ! filled with zeros, as per CMIP requirements
   integer :: n_sends! number of data points sent to the field since last dump
   integer :: alias = 0 ! ID of the first alias in the chain
   character(32) :: module,name ! for debugging purposes only
end type tiled_diag_field_type


! ==== module data ===========================================================
logical :: module_is_initialized = .false.

! list of registered fields
type(tiled_diag_field_type), pointer :: fields(:) => NULL()
integer :: n_fields       = 0 ! current number of diag fields
integer :: current_offset = 1 ! current total size of the diag fields per tile



contains



! ============================================================================
subroutine tile_diag_init()

  if (module_is_initialized) return

  module_is_initialized = .true.
  call log_version(version, mod_name, &
  __FILE__)

  ! initialize diag selectors
  call tile_selectors_init()
  call register_tile_selector('land',area_depends_on_time=.FALSE.)

  ! initialize global data
  allocate(fields(INIT_FIELDS_SIZE))
  n_fields       = 0
  current_offset = 1

end subroutine tile_diag_init



! ============================================================================
subroutine tile_diag_end()

  integer :: i

  ! deallocate global data
  do i = 1, n_fields
     deallocate(fields(i)%ids)
  end do
  deallocate(fields)

  ! destroy selectors
  call tile_selectors_end()

  module_is_initialized = .false.

end subroutine tile_diag_end

! ============================================================================
! give a name of the diagnostic selector, returns id of area variable associated
! with this selector
function get_area_id(name); integer get_area_id
   character(*), intent(in) :: name ! name of the selector

   integer :: i

   do i = 1, n_selectors
      if (trim(name)==trim(selectors(i)%name)) then
         get_area_id = selectors(i)%area_id
         return
      endif
   enddo
   call error_mesg(mod_name,&
      'diag filter "'//trim(name)//'" was not found among registered diag filters', FATAL)
end function get_area_id

! ============================================================================
! sets the default tile diagnostic selector: all tiled diag fields registered
! after calling this function (until the next call) will be associated with the
! specified selector
subroutine set_default_diag_filter(name)
   character(*), intent(in) :: name ! name of the selector

   integer :: i

   selectors(:)%is_default = .FALSE.
   do i = 1, n_selectors
      if (trim(name)==trim(selectors(i)%name)) then
         selectors(i)%is_default = .TRUE.
         return
      endif
   enddo
   call error_mesg(mod_name,&
      'diag filter "'//trim(name)//'" was not found among registered diag filters', FATAL)
end subroutine set_default_diag_filter


! ============================================================================
! area and frac are special, because:
! * static/dynamic is controlled by the selectors, not by the registration
! * frac is always associated with the land area
! * frac is average over land area, despite being calculated with OP_SUM
! TODO: somehow make sure that frac_land is the fraction of the area in the
!       grid cell, not 1
subroutine register_tiled_area_fields(module_name, axes, init_time, &
     id_area, id_frac)

  character(len=*), intent(in)  :: module_name
  integer,          intent(in)  :: axes(:)
  type(time_type),  intent(in)  :: init_time
  integer,          intent(out) :: id_area, id_frac

  integer :: i_area, k

  ! register areas for all tiles
  id_area = reg_field(FLD_LIKE_AREA, module_name, 'area', init_time, axes, &
         'area in the grid cell', 'm2', missing_value=-1.0, op=OP_SUM)
  if (id_area>0) then
     call add_cell_methods(id_area,'area: sum')
  endif
  ! store the ids of area for each of the selectors
  if (id_area > 0) then
     i_area = id_area - BASE_TILED_FIELD_ID
     do k = 1, n_selectors
        selectors(k)%area_id = fields(i_area)%ids(k)
     enddo
  endif
  id_frac = reg_field(FLD_LIKE_AREA, module_name, 'frac', init_time, axes, &
         'fraction of land area', 'unitless', missing_value=-1.0, op=OP_SUM, &
         standard_name='area_fraction')
  if (id_frac > 0) then
     call add_cell_measures(id_frac, get_area_id('land'))
     call add_cell_methods(id_frac,'area: mean')
  endif
end subroutine register_tiled_area_fields


! ============================================================================
! sets cell_measures for all selectors of the tiled diag field
subroutine add_cell_measures(id, area)
   integer, intent(in) :: id ! id of the tiled diag field
   integer, intent(in), optional :: area ! id of the area diag field.
     ! If "area" is present, cell_measures for all selectors are set to indicated
     ! output field; otherwise diag field for each selector is associated with the
     ! area for the selector.

   integer :: i,k, area_id

   if (id<=0) return ! do nothing if the field is not registered

   i = id - BASE_TILED_FIELD_ID
   do k = 1, n_selectors
      if (present(area)) then
         area_id = area
      else
         area_id = selectors(k)%area_id
      endif
      if (fields(i)%ids(k)>0) then
         call diag_field_add_cell_measures(fields(i)%ids(k),area_id)
      endif
   enddo
end subroutine add_cell_measures


! ============================================================================
! adds cell_methods attribute to all selectors of specified tiled field
! if cell_method is present, its value is used for all selectors; otherwise
! the value is deduced based on aggreagtion opertaion code
subroutine add_cell_methods(id,cell_methods)
  integer, intent(in) :: id ! id of the tiled diagnostic field
  character(*), intent(in), optional :: cell_methods ! value cell_method attribute

  integer :: i,k
  character(64) :: cell_methods_

  if (id<=0) return ! do nothing if the field is not registered

  i = id - BASE_TILED_FIELD_ID
  if (present(cell_methods)) then
     cell_methods_ = cell_methods
  else
     select case (fields(i)%op)
     case (OP_AVERAGE)
        cell_methods_ = 'area: mean'
     case (OP_SUM)
        cell_methods_ = 'area: sum'
     case (OP_VAR)
        cell_methods_ = 'area: variance (of sub-grid tile values)'
     case (OP_STD)
        cell_methods_ = 'area: standard_deviation (of sub-grid tile values)'
     case default
        call error_mesg('add_cell_methods', 'unknown aggregation operation for field "'//&
               trim(fields(i)%module)//'/'//trim(fields(i)%name)//'"',FATAL)
     end select
  endif

  do k = 1, n_selectors
     if (fields(i)%ids(k)>0) then
        call diag_field_add_attribute(fields(i)%ids(k),'cell_methods',trim(cell_methods_))
     endif
  enddo
end subroutine add_cell_methods


! ============================================================================
function register_tiled_diag_field(module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op, standard_name, fill_missing) result (id)

  integer :: id

  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  integer,          intent(in), optional :: op ! aggregation operation code
  character(len=*), intent(in), optional :: standard_name
  logical,          intent(in), optional :: fill_missing

  id = reg_field(FLD_DYNAMIC, module_name, field_name, init_time, axes, long_name, &
         units, missing_value, range, op=op, standard_name=standard_name, fill_missing=fill_missing)
  call add_cell_measures(id)
  call add_cell_methods(id)
end function register_tiled_diag_field

! ============================================================================
function register_tiled_static_field(module_name, field_name, axes, &
     long_name, units, missing_value, range, require, op, standard_name, fill_missing) result (id)

  integer :: id

  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  logical,          intent(in), optional :: require
  integer,          intent(in), optional :: op ! aggregation operation code
  character(len=*), intent(in), optional :: standard_name
  logical,          intent(in), optional :: fill_missing

  ! --- local vars
  type(time_type) :: init_time

  id = reg_field(FLD_STATIC, module_name, field_name, init_time, axes, long_name, &
         units, missing_value, range, require, op, standard_name=standard_name, &
         fill_missing=fill_missing)
  call add_cell_measures(id)
  call add_cell_methods(id)
end function register_tiled_static_field


! ============================================================================
subroutine add_tiled_static_field_alias(id0, module_name, field_name, axes, &
     long_name, units, missing_value, range, op, standard_name)
  integer,          intent(inout) :: id0 ! id of the original diag field on input;
   ! if negative then it may be replaced with the alias id on output
  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  integer,          intent(in), optional :: op ! aggregation operation code
  character(len=*), intent(in), optional :: standard_name

  ! --- local vars
  type(time_type) :: init_time

  call reg_field_alias(id0, FLD_STATIC, module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op, standard_name=standard_name)
end subroutine add_tiled_static_field_alias


! ============================================================================
subroutine add_tiled_diag_field_alias(id0, module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op, standard_name, fill_missing)
  integer,          intent(inout) :: id0 ! id of the original diag field on input;
   ! if negative then it may be replaced with the alias id on output
  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  integer,          intent(in), optional :: op ! aggregation operation code
  character(len=*), intent(in), optional :: standard_name
  logical,          intent(in), optional :: fill_missing

  call reg_field_alias(id0, FLD_DYNAMIC, module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op, standard_name=standard_name, &
     fill_missing=fill_missing)
end subroutine add_tiled_diag_field_alias

! ============================================================================
subroutine reg_field_alias(id0, static, module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, op, standard_name, fill_missing)


  integer,          intent(inout) :: id0 ! id of the original diag field on input;
  integer,          intent(in) :: static
   ! if negative then it may be replaced with the alias id on output
  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  integer,          intent(in), optional :: op ! aggregation operation code
  character(len=*), intent(in), optional :: standard_name
  logical,          intent(in), optional :: fill_missing

  ! local vars
  integer :: id1
  integer :: ifld0, ifld1

  if (id0>0) then
    ifld0 = id0-BASE_TILED_FIELD_ID
    if (ifld0<1.or.ifld0>n_fields) &
       call error_mesg(mod_name, 'incorrect index ifld0 '//string(ifld0)//&
                    ' in definition of tiled diag field alias "'//&
                    trim(module_name)//'/'//trim(field_name)//'"', FATAL)
    id1 = reg_field(static, module_name, field_name, init_time, axes, long_name, &
          units, missing_value, range, op=op, offset=fields(ifld0)%offset, &
          fill_missing=fill_missing)
    call add_cell_measures(id1)
    call add_cell_methods(id1)
    if (id1>0) then
      ifld1 = id1-BASE_TILED_FIELD_ID
      ! check that sizes of the fields are identical
      if (fields(ifld0)%size/=fields(ifld1)%size) &
         call error_mesg(mod_name, 'sizes of diag field "'//              &
           trim(fields(ifld0)%module)//'/'//trim(fields(ifld0)%name)//    &
           '" and its alias "'//trim(module_name)//'/'//trim(field_name)//&
           '" are not the same', FATAL)
      ! check that "static" status of the fields is the same
      if(fields(ifld0)%static.ne.fields(ifld1)%static) &
         call error_mesg(mod_name,                                        &
           'attempt to register alias "'//trim(module_name)//'/'//trim(field_name)//   &
           '" of field "'//trim(fields(ifld0)%module)//'/'//trim(fields(ifld0)%name)// &
           '" with different static/dynamic status',&
           FATAL)
      ! copy alias field from the original into the alias, to preserve the chain
      fields(ifld1)%alias = fields(ifld0)%alias
      ! update alias field in the head of alias chain
      fields(ifld0)%alias = ifld1
    endif
  else
    ! the "main" field has not been registered, so simply register the alias
    ! as a diag field
    id0 = reg_field(static, module_name, field_name, init_time, axes, long_name, &
          units, missing_value, range, op=op, standard_name=standard_name,&
          fill_missing=fill_missing)
    call add_cell_measures(id0)
    call add_cell_methods(id0)
  endif
end subroutine reg_field_alias

! ============================================================================
! provides unified interface for registering a diagnostic field with full set
! of selectors
function reg_field(static, module_name, field_name, init_time, axes, &
     long_name, units, missing_value, range, require, op, offset, &
     area, cell_methods, standard_name, fill_missing) result(id)

  integer :: id

  integer,          intent(in) :: static
  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  logical,          intent(in), optional :: require
  integer,          intent(in), optional :: op
  integer,          intent(in), optional :: offset
  character(len=*), intent(in), optional :: area ! name of the area associated with this field, if not default
  character(len=*), intent(in), optional :: cell_methods ! cell_methods associated with this field, if not default
  character(len=*), intent(in), optional :: standard_name
  logical,          intent(in), optional :: fill_missing ! if true, missing values (e.g. ocean points)
                                                         ! are filled with zeros, as per CMIP requirements

  ! ---- local vars
  integer, pointer :: diag_ids(:) ! ids returned by FMS diag manager for each selector
  integer :: i, j
  type(tiled_diag_field_type), pointer :: new_fields(:)
  ! ---- global vars: n_fields, fields, current_offset -- all used and updated

  ! log diagnostic field information
  call log_diag_field_info ( module_name, trim(field_name), axes, long_name, units,&
                             missing_value, range, dynamic=(static.ne.FLD_STATIC) )

  ! go through all possible selectors and try to register a diagnostic field
  ! with the name derived from field name and selector; if any of the
  ! registrations succeeds, return a tiled field id, otherwise return 0.
  ! Note that by design one of the selectors have empty name and selects all
  ! the tiles.
  id = 0
  allocate(diag_ids(n_selectors))

  do i = 1, n_selectors
     diag_ids(i) = reg_field_set(static, selectors(i), &
          module_name, field_name, axes, &
          init_time, long_name, units, missing_value, range, require, &
          standard_name=standard_name)
  enddo

  if(any(diag_ids>0)) then
     ! if any of the field+selector pairs was found for this field, an entry
     ! must be added to the table of tile diagnostic fields

     ! if there is not enough slots in the field table to add another one,
     ! allocate more space
     if(n_fields>=size(fields)) then
        allocate(new_fields(max(2*n_fields,1)))
        new_fields(1:n_fields) = fields(1:n_fields)
        deallocate(fields)
        fields => new_fields
     endif
     ! add the current field to the field table
     n_fields = n_fields+1
     id       = n_fields
     ! set the array of FMS diagnostic field IDs for each selector
     fields(id)%ids => diag_ids
     ! set the field offset in the diagnostic buffers
     if (present(offset)) then
        fields(id)%offset = offset
     else
        fields(id)%offset = current_offset
     endif
     ! calculate field size per tile and increment current offset to
     ! reserve space in per-tile buffers. We assume that the first two axes
     ! are horizontal coordinates, so their size is not taken into account
     fields(id)%size = 1
     do i = 2, size(axes(:))
        fields(id)%size = fields(id)%size * get_axis_length(axes(i))
     enddo
     ! if offset is present in the list of the arguments, it means that we do not
     ! want to increase the current_offset -- this is an alias field
     if (.not.present(offset)) &
        current_offset = current_offset + fields(id)%size
     ! store the code of the requested tile aggregation operation
     if(present(op)) then
        fields(id)%op = op
     else
        fields(id)%op = OP_AVERAGE
     endif
     ! store the static field flag
     fields(id)%static = static
     ! zero out the number of data points sent to the field
     fields(id)%n_sends = 0
     ! store the name of the field -- for now, only to be able to see what it is
     ! in the debugger
     fields(id)%module=module_name
     fields(id)%name=field_name
     ! store the filler flag
     fields(id)%fill_missing = .FALSE.
     if(present(fill_missing))fields(id)%fill_missing = fill_missing
     ! increment the field id by some (large) number to distinguish it from the
     ! IDs of regular FMS diagnostic fields
     id = id + BASE_TILED_FIELD_ID
  else
     deallocate(diag_ids)
  endif

end function reg_field


! ============================================================================
! provides unified interface for registering a diagnostic field with a given
! selector, whether static or time-dependent
function reg_field_set(static, sel, module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, require, area, standard_name) result (id)

  integer :: id

  integer,          intent(in) :: static
  type(tile_selector_type), intent(in) :: sel
  character(len=*), intent(in) :: module_name
  character(len=*), intent(in) :: field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), intent(in), optional :: long_name
  character(len=*), intent(in), optional :: units
  real,             intent(in), optional :: missing_value
  real,             intent(in), optional :: range(2)
  logical,          intent(in), optional :: require
  integer,          intent(in), optional :: area
  character(len=*), intent(in), optional :: standard_name

  character(len=128) :: fname
  logical :: static_

  ! form field name as concatenation of name of the field and selector suffix
  fname = trim(field_name)//trim(selector_suffix(sel))

  ! try registering diagnostic field with FMS diagnostic manager.
  select case(static)
  case (FLD_STATIC)
     static_=.TRUE.
  case (FLD_DYNAMIC)
     static_=.FALSE.
  case (FLD_LIKE_AREA)
     static_ = sel%area_is_static
  end select
  if (static_) then
     id = register_static_field ( module_name, fname,   &
          axes, long_name, units, missing_value, range, require, &
          do_not_log=.TRUE., area=area, standard_name=standard_name )
  else
     id = register_diag_field ( module_name,  fname,   &
          axes, init_time, long_name, units, missing_value, range, &
          mask_variant=.true., do_not_log=.TRUE., area=area, &
          standard_name=standard_name )
  endif

end function reg_field_set


! ============================================================================
subroutine send_tile_data_0d(id, x, buffer)
  integer, intent(in) :: id
  real   , intent(in) :: x
  type(diag_buff_type), intent(inout) :: buffer

  integer :: idx, i
#ifdef DEBUG_LAND_TILE_DIAG
  real*4 r4
#endif

  if (id <= 0) return

  ! reallocate diagnostic buffer according to the current number and size of
  ! tiled diag fields
  call realloc_diag_buff(buffer,current_offset)

  ! calculate offset for the current diagnostic field in the buffer
  i = id - BASE_TILED_FIELD_ID
  idx = fields(i)%offset

#ifdef DEBUG_LAND_TILE_DIAG
  ! DEBUG : check validity of input data
  call check_var_range(x,REAL(-HUGE(r4)),REAL(HUGE(r4)),'send_tile_data_0d',fields(i)%name, FATAL)
#endif

  ! store the diagnostic data
  buffer%data(idx) = x
  buffer%mask(idx) = .TRUE.

  ! increment sent data counter
  fields(i)%n_sends = fields(i)%n_sends + 1
  ! increment sent data counter in all aliases
  do while(fields(i)%alias>0)
    i=fields(i)%alias
    fields(i)%n_sends = fields(i)%n_sends + 1
  enddo
end subroutine send_tile_data_0d

! ============================================================================
subroutine send_tile_data_1d(id, x, buffer)
  integer, intent(in) :: id
  real   , intent(in) :: x(:)
  type(diag_buff_type), intent(inout) :: buffer

  integer :: is, ie, i
#ifdef DEBUG_LAND_TILE_DIAG
  real*4 r4
#endif

  if (id <= 0) return

  ! reallocate diagnostic buffer according to the current number and size of
  ! tiled diag fields
  call realloc_diag_buff(buffer, current_offset)

  ! calculate offset for the current diagnostic field in the buffer
  i = id - BASE_TILED_FIELD_ID ! index in the array of fields
  is = fields(i)%offset ; ie = is+fields(i)%size-1

#ifdef DEBUG_LAND_TILE_DIAG
  ! DEBUG : check validity of input data
  call check_var_range(x,REAL(-HUGE(r4)),REAL(HUGE(r4)),'send_tile_data_1d',fields(i)%name, FATAL)
#endif

  ! store the data
  buffer%data(is:ie) = x(:)
  buffer%mask(is:ie) = .TRUE.

  ! increment sent data counter
  fields(i)%n_sends = fields(i)%n_sends + 1
  ! increment sent data counter in all aliases
  do while(fields(i)%alias>0)
    i=fields(i)%alias
    fields(i)%n_sends = fields(i)%n_sends + 1
  enddo
end subroutine send_tile_data_1d

! NOTE: 2-d fields can be handled similarly to 1-d with reshape

! ============================================================================
subroutine send_tile_data_0d_array(id, x, send_immediately)
  integer, intent(in) :: id
  real   , intent(in) :: x(:,:)
  logical, intent(in), optional :: send_immediately ! if true, send data to diag_manager
    ! right avay, instead of waiting for the next tiled diag fields dump

  integer :: l,k
  type(land_tile_enum_type)     :: ce
  type(land_tile_type), pointer :: tileptr

  ce = first_elmt( land_tile_map )
  do while(loop_over_tiles(ce,tileptr,l,k))
     call send_tile_data(id,x(l,k),tileptr%diag)
  enddo
  if(present(send_immediately)) then
     ! TODO: perhaps need to add time to the arguments, instead of using lnd%time?
     ! not clear if this will have any effect
     if(send_immediately) call dump_tile_diag_field(land_tile_map, id, lnd%time)
  endif
end subroutine send_tile_data_0d_array

! ============================================================================
subroutine send_tile_data_r0d_fptr(id, fptr)
  integer, intent(in) :: id
  procedure(fptr_r0)  :: fptr

  type(land_tile_enum_type)     :: ce      ! tile list enumerator
  type(land_tile_type), pointer :: tileptr ! pointer to tile
  real                , pointer :: ptr     ! pointer to the data element within a tile

  if(id <= 0) return
  ce = first_elmt( land_tile_map )
  do while(loop_over_tiles(ce,tileptr))
     call fptr(tileptr,ptr)
     if(associated(ptr)) call send_tile_data(id,ptr,tileptr%diag)
  enddo
end subroutine send_tile_data_r0d_fptr


! ============================================================================
subroutine send_tile_data_r1d_fptr(id, fptr)
  integer, intent(in) :: id
  procedure(fptr_r0i) :: fptr

  type(land_tile_enum_type)     :: ce      ! tile list enumerator
  type(land_tile_type), pointer :: tileptr ! pointer to tile
  real                , pointer :: ptr     ! pointer to the data element within a tile
  real,             allocatable :: buffer(:) ! buffer for accumulating data
  integer :: i, i2,i3, k
  logical :: have_data

  if(id <= 0) return
  i = id - BASE_TILED_FIELD_ID ! index in the array of fields

  allocate(buffer(fields(i)%size))
  ce = first_elmt( land_tile_map, ls=lnd%ls )
  do while(loop_over_tiles(ce, tileptr,i2,i3))
#ifdef DEBUG_LAND_TILE_DIAG
     call set_current_point(i2,i3)
#endif
     have_data = .FALSE.
     do k = 1,fields(i)%size
        call fptr(tileptr,k,ptr)
        if(associated(ptr)) then
            buffer(k) = ptr
            have_data = .TRUE.
        else
            buffer(k) = 0.0
        endif
     enddo
     if (have_data) call send_tile_data(id,buffer,tileptr%diag)
  enddo
  deallocate(buffer)
end subroutine send_tile_data_r1d_fptr


! ============================================================================
subroutine send_tile_data_i0d_fptr(id, fptr)
  integer, intent(in) :: id
  procedure(fptr_i0)  :: fptr

  type(land_tile_enum_type)     :: ce      ! tile list enumerator
  type(land_tile_type), pointer :: tileptr ! pointer to tile
  integer             , pointer :: ptr     ! pointer to the data element within a tile

  if(id <= 0) return
  ce = first_elmt( land_tile_map )
  do while(loop_over_tiles(ce,tileptr))
     call fptr(tileptr,ptr)
     if(associated(ptr)) call send_tile_data(id,real(ptr),tileptr%diag)
  enddo
end subroutine send_tile_data_i0d_fptr


! ============================================================================
!pass in land_tile_map into this routine to temporarily solve the crash issue
!with Intel compiler when running with multiple openmp threads.
subroutine dump_tile_diag_fields(land_tile_map,time)
  type(land_tile_list_type), intent(in) :: land_tile_map(:) ! map of tiles
  type(time_type)          , intent(in) :: time       ! current time

  ! ---- local vars
  integer :: ifld ! field number
  integer :: isel ! selector number
  type(land_tile_enum_type)     :: ce
  type(land_tile_type), pointer :: tile
  integer :: total_n_sends(n_fields)

  total_n_sends(:) = fields(1:n_fields)%n_sends
  call mpp_sum(total_n_sends, n_fields, pelist=lnd%pelist)

!$OMP parallel do default(none) shared(land_tile_map,n_fields,total_n_sends,n_selectors,fields,selectors,time)
  do ifld = 1, n_fields
     if (total_n_sends(ifld) == 0) cycle ! no data to send
     do isel = 1, n_selectors
        if (fields(ifld)%ids(isel) <= 0) cycle
        call dump_diag_field_with_sel (land_tile_map, fields(ifld)%ids(isel), &
             fields(ifld), selectors(isel), time )
     enddo
  enddo
  ! zero out the number of data points sent to the field
  fields(1:n_fields)%n_sends=0

  ! all the data are sent to the output, so set the data presence tag to FALSE
  ! in all diag buffers in preparation for the next time step
  ce = first_elmt(land_tile_map)
  do while(loop_over_tiles(ce,tile))
    tile%diag%mask(:) = .FALSE.
  enddo
end subroutine dump_tile_diag_fields

! ============================================================================
! dumps a single field
! TODO: perhaps need dump aliases as well
! TODO: perhaps total_n_sends check can be removed to avoid communication
! pass in land_tile_map into this routine to temporarily solve the crash issue
! with Intel compiler when running with multiple openmp threads.
subroutine dump_tile_diag_field(land_tile_map, id, time)
  type(land_tile_list_type), intent(in) :: land_tile_map(:)   ! map of tiles
  integer, intent(in) :: id ! diag id of the field
  type(time_type), intent(in) :: time       ! current time

  ! ---- local vars
  integer :: ifld ! field number
  integer :: isel ! selector number
  type(land_tile_enum_type)     :: ce
  type(land_tile_type), pointer :: tile
  integer :: total_n_sends

  if (id<=0) return ! do nothing if field not registered

  ifld = id-BASE_TILED_FIELD_ID
  if (ifld<1.or.ifld>n_fields) &
     call error_mesg(mod_name, 'incorrect field id '//string(id)//' in dump_tile_diag_field ', FATAL)

  total_n_sends = fields(ifld)%n_sends
  call mpp_sum(total_n_sends, pelist=lnd%pelist)

  if (total_n_sends == 0) return ! no data to send
!$OMP parallel do default(none) shared(land_tile_map,n_selectors,fields,ifld,selectors,time) private(isel)
  do isel = 1, n_selectors
     if (fields(ifld)%ids(isel) <= 0) cycle
     call dump_diag_field_with_sel (land_tile_map, fields(ifld)%ids(isel), &
          fields(ifld), selectors(isel), time )
  enddo
  ! zero out the number of data points sent to the field
  fields(ifld)%n_sends=0

  ! all the data are sent to the output, so set the data presence tag to FALSE
  ! in all diag buffers in preparation for the next time step
  ce = first_elmt(land_tile_map)
  do while(loop_over_tiles(ce,tile))
    tile%diag%mask(fields(ifld)%offset:fields(ifld)%offset+fields(ifld)%size-1) = .FALSE.
  enddo

end subroutine dump_tile_diag_field

! ============================================================================
subroutine dump_diag_field_with_sel(land_tile_map, id, field, sel, time)
  type(land_tile_list_type)  , intent(in) :: land_tile_map(:)
  integer                    , intent(in) :: id
  type(tiled_diag_field_type), intent(in) :: field
  type(tile_selector_type)   , intent(in) :: sel
  type(time_type)            , intent(in) :: time ! current time

! NOTE that passing in land_tile_map (despite the fact that this array is also globally
! available) is a work around (apparent) compiler issue, when with multiple openmp threads
! *and* debug flags Intel compilers (15 and 16) report index errors, as if global
! land_tile_map array started from 1,1 instead of is,js. Passing it in as argument solves
! this issue.

  ! ---- local vars
  integer :: l ! iterators
  integer :: ks,ke ! array boundaries
  integer :: ls, le
  logical :: used ! value returned from send_data (ignored)
  real, allocatable :: buffer(:,:), weight(:,:), var(:,:)
  logical, allocatable :: mask(:,:)
  type(land_tile_enum_type)     :: ce
  type(land_tile_type), pointer :: tile

  ! calculate array boundaries
  ls = lbound(land_tile_map,1); le = ubound(land_tile_map,1)
  ks = field%offset   ; ke = field%offset + field%size - 1

  ! allocate and initialize temporary buffers
  allocate(buffer(ls:le,ks:ke), weight(ls:le,ks:ke), mask(ls:le,ks:ke))
  weight(:,:) = 0.0
  buffer(:,:) = 0.0

  ! accumulate data
  ce = first_elmt(land_tile_map, ls=ls)
  do while(loop_over_tiles(ce, tile, l))
    if ( size(tile%diag%data) < ke )       cycle ! do nothing if there is no data in the buffer
    if ( .not.tile_is_selected(tile,sel) ) cycle ! do nothing if tile is not selected
    select case (field%op)
    case (OP_AVERAGE,OP_VAR,OP_STD)
       where(tile%diag%mask(ks:ke))
          buffer(l,:) = buffer(l,:) + tile%diag%data(ks:ke)*tile%frac
       end where
       weight(l,:) = weight(l,:) + tile%frac
    case (OP_SUM)
       where(tile%diag%mask(ks:ke))
          buffer(l,:) = buffer(l,:) + tile%diag%data(ks:ke)
       end where
       weight(l,:) = 1
    end select
  enddo

  ! normalize accumulated data
  mask = (weight>0)
  where (mask) buffer=buffer/weight

  if (field%op == OP_VAR.or.field%op==OP_STD) then
     ! second loop to process the variance and standard deviation diagnostics.
     ! it may be possible to calc. var and std in one pass with weighted incremental
     ! algorithm from http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
     ! code her is more straightforward. buffer(:,:,:) already contains the mean,
     ! and weight(:,:,:) -- sum of  tile fractions
     allocate(var(ls:le,ks:ke))
     var(:,:) = 0.0
     ! the loop is somewhat different from the first, for no particular reason:
     ! perhaps this way is better for performance?
     do l = ls, le
        ce = first_elmt(land_tile_map(l))
        do while(loop_over_tiles(ce,tile))
           if ( size(tile%diag%data) < ke )       cycle ! do nothing if there is no data in the buffer
           if ( .not.tile_is_selected(tile,sel) ) cycle ! do nothing if tile is not selected
           where(tile%diag%mask(ks:ke))
              var(l,:) = var(l,:) + tile%frac*(tile%diag%data(ks:ke)-buffer(l,:))**2
           end where
        enddo
     enddo
     ! renormalize the variance or standard deviation. note that weight is
     ! calculated in the first loop
     select case (field%OP)
     case (OP_VAR)
         where (mask) buffer = var/weight
     case (OP_STD)
         where (mask) buffer = sqrt(var/weight)
     end select
     deallocate(var)
  endif

  if (field%fill_missing) then
      where (.not. mask) buffer = 0.0
      mask = .true.
  endif
  ! send diag field
  used = send_data(id,buffer,time,mask=mask)

  ! clean up temporary data
  deallocate(buffer,weight,mask)

end subroutine dump_diag_field_with_sel

  !#######################################################################
  !> \brief Send out the land model field on unstructured grid for global integral
  logical function send_global_land_diag( id, diag, Time, tile, mask, Land )
  integer,                 intent(in) :: id
  real,    dimension(:,:), intent(in) :: diag, tile
  type(time_type),         intent(in) :: Time
  logical, dimension(:,:), intent(in) :: mask
  type(land_data_type),    intent(in) :: Land

  real,    dimension(size(diag,1),1)    :: diag_ug, tile_ug, area_ug
  logical, dimension(size(mask,1))    :: mask_ug
  integer :: k
  real    :: area_sum, diag_sum

    ! sum over tiles on unstructured grid
    diag_ug = 0.0
    tile_ug = 0.0
    do k = 1, size(diag,2)
      where (mask(:,k))
        diag_ug(:,1) = diag_ug(:,1) + diag(:,k)*tile(:,k)
        tile_ug(:,1) = tile_ug(:,1) + tile(:,k)
      endwhere
    enddo
    ! average on unstructured grid
    where (tile_ug > 0.0)
      diag_ug = diag_ug/tile_ug
    endwhere
    mask_ug(:) = ANY(mask,dim=2)

    where(mask_ug)
       diag_ug(:,1) = diag_ug(:,1) * lnd%area
       area_ug(:,1) = lnd%area
    elsewhere
       diag_ug(:,1) = 0.0
       area_ug(:,1) = 0.0
    endwhere

    area_sum = mpp_reproducing_sum(diag_ug) !, overflow_check=.true.)
    diag_sum = mpp_reproducing_sum(area_ug) !, overflow_check=.true.)

    send_global_land_diag = send_data( id, diag_sum/area_sum, Time)

  end function send_global_land_diag


end module land_tile_diag_mod
