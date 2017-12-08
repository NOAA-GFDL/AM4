module atmos_global_diag_mod

!----------------------------------------------------------------------
!  Module for computed globally averaged atmospheric quantities
!  and then registering and sending this data to the diag_manager.
!----------------------------------------------------------------------

use mpp_mod,          only: input_nml_file
use mpp_domains_mod,  only: domain2d, mpp_global_sum, BITWISE_EFP_SUM, &
                            null_domain2d, operator(.eq.)
use fms_mod,          only: open_namelist_file, check_nml_error, &
                            close_file, stdlog, mpp_pe, mpp_root_pe, &
                            write_version_number, file_exist, &
                            error_mesg, FATAL, WARNING
use time_manager_mod, only: time_type
use diag_manager_mod, only: register_diag_field, send_data, &
                            DIAG_FIELD_NOT_FOUND
!use diag_data_mod,    only: null_axis_id, &
use diag_data_mod,    only: CMOR_MISSING_VALUE
use diag_axis_mod,    only: get_domain2d

!-----------------------------------------------------------------------

implicit none
private

!-----------------------------------------------------------------------
! public interfaces

public :: atmos_global_diag_init, &
          register_global_diag_field, &
          get_global_diag_field_id, &
          buffer_global_diag, &
          send_global_diag, &
          atmos_global_diag_end

!-----------------------------------------------------------------------

interface send_global_diag
   module procedure send_global_diag_data
   module procedure send_global_diag_buffer
end interface

!-----------------------------------------------------------------------
! private data structure

type atmos_global_diag_type
  character(len=128) :: field_name
  integer            :: field_id
  real, pointer      :: buffer(:,:)
  type(time_type)    :: Time
  logical            :: use_buffer
  logical            :: use_masking
end type atmos_global_diag_type

!-----------------------------------------------------------------------
! private module data

real, allocatable, dimension(:,:) :: area_g
real :: area_g_sum
integer :: id, jd
!integer :: null_axis(1)
real :: missing_value = CMOR_MISSING_VALUE

character(len=12) :: mod_name = 'atmos_global'
type(domain2d), save :: Domain2

integer :: num_fields = 0
type(atmos_global_diag_type), allocatable :: fields(:)

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

logical :: module_is_initialized=.false.

!-----------------------------------------------------------------------
! namelist

integer :: max_fields = 100

namelist /atmos_global_diag_nml/ max_fields


CONTAINS

!#######################################################################

subroutine atmos_global_diag_init (axes, area)
integer, intent(in), dimension(:)   :: axes
real,    intent(in), dimension(:,:) :: area

!-----------------------------------------------------------------------
! local data

integer :: iunit, ierr, io

!-----------------------------------------------------------------------

  if (module_is_initialized) then
    call error_mesg ('atmos_global_diag_mod', &
                     'module has already been initialized', WARNING)
    return
  endif

!-----------------------------------------------------------------------
!----- read namelist -----
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=atmos_global_diag_nml, iostat=io)
  ierr = check_nml_error (io, 'atmos_global_diag_nml')
#else
  if (file_exist('input.nml') ) then
    iunit = open_namelist_file()
    ierr=1
    do while (ierr /= 0)
      read (iunit, nml=atmos_global_diag_nml, iostat=io, end=10)
      ierr = check_nml_error (io, 'atmos_global_diag_nml')
    enddo
10  call close_file (iunit)
  endif
#endif

!----- write version and namelist to log file -----

  iunit = stdlog()
  call write_version_number ( version, tagname )
  if (mpp_pe() == mpp_root_pe()) write (iunit, nml=atmos_global_diag_nml)

!-----------------------------------------------------------------------
!----- allocate storage for fields -----
  allocate(fields(max_fields))
  num_fields = 0

!----- determine the atmospheric 2D domain from the axes using diag_manager -----
  Domain2 = get_domain2d(axes)
  if (Domain2 .eq. null_domain2d) then
    call error_mesg ('atmos_global_diag_mod', &
                     'could not determine 2d domain', FATAL)
  endif

!----- save cell areas and the computed global area ----
  id = size(area,1); jd = size(area,2)
  allocate(area_g(id,jd))
  area_g = area
  area_g_sum = mpp_global_sum (Domain2, area_g, flags=BITWISE_EFP_SUM)

! null axis for global fields in diag_manager
! null_axis(1) = null_axis_id

  module_is_initialized = .true.

!-----------------------------------------------------------------------

end subroutine atmos_global_diag_init

!#######################################################################

function register_global_diag_field (field_name, Time_init, long_name, &
                                     units, standard_name, realm, buffer, use_masking )
character(len=*), intent(in) :: field_name
type(time_type),  intent(in) :: Time_init
character(len=*), intent(in), optional :: long_name, units, standard_name, realm
logical,          intent(in), optional :: buffer, use_masking

!-----------------------------------------------------------------------

integer :: register_global_diag_field
integer :: field_index, field_id

!-----------------------------------------------------------------------

  if (.not.module_is_initialized) call error_mesg &
      ('atmos_global_diag_mod', 'module has not been initialized', FATAL)

  register_global_diag_field = 0

  field_index = find_field_index(field_name)

  if (field_index > 0) call error_mesg ('atmos_global_diag_mod', &
                 'trying to register field already registered', FATAL)

  ! register field with diag_manager
  field_id = register_diag_field (mod_name, field_name, Time_init, &
                                  long_name, units, standard_name=standard_name, &
                                  realm=realm)
  if (field_id == DIAG_FIELD_NOT_FOUND) return

  ! register field with this module
  num_fields = num_fields+1
  if (num_fields > max_fields) call error_mesg ('atmos_global_diag_mod', &
              'exceeded max_fields, increase value in namelist', FATAL)
  register_global_diag_field = num_fields

  fields(num_fields)%field_name = trim(field_name)
  fields(num_fields)%field_id = field_id
  fields(num_fields)%use_masking = .false.
  fields(num_fields)%use_buffer = .false.

  if (present(buffer)) then
    if (buffer) then
      allocate(fields(num_fields)%buffer(id,jd))
      fields(num_fields)%buffer = missing_value
      fields(num_fields)%use_buffer = .true.
      if (present(use_masking)) fields(num_fields)%use_masking = use_masking
    endif
  endif

  return

!-----------------------------------------------------------------------

end function register_global_diag_field

!#######################################################################

function get_global_diag_field_id (field_num)
integer, intent(in) :: field_num

integer :: get_global_diag_field_id
!-----------------------------------------------------------------------
  get_global_diag_field_id = DIAG_FIELD_NOT_FOUND

  if (.not.module_is_initialized) call error_mesg &
      ('atmos_global_diag_mod', 'module has not been initialized', FATAL)

  if (field_num == 0) return
  if (field_num < 0 .or. field_num > num_fields) call error_mesg &
      ('atmos_global_diag_mod', 'invalid field number in get_global_diag_field_id', FATAL)

  get_global_diag_field_id = fields(field_num)%field_id

!-----------------------------------------------------------------------

end function get_global_diag_field_id

!#######################################################################

subroutine buffer_global_diag (field_num, data, Time, is_in, js_in, mask)
integer,           intent(in) :: field_num
real,              intent(in) :: data(:,:)
type(time_type),   intent(in) :: Time
integer, optional, intent(in) :: is_in, js_in
logical, optional, intent(in) :: mask(:,:)
!-----------------------------------------------------------------------

integer :: is, ie, js, je

!-----------------------------------------------------------------------

  if (.not.module_is_initialized) call error_mesg &
      ('atmos_global_diag_mod', 'module has not been initialized', FATAL)

  if (field_num == 0) return
  if (field_num < 0 .or. field_num > num_fields) call error_mesg &
      ('atmos_global_diag_mod', 'invalid field number in buffer_global_data', FATAL)

  if (.not.fields(field_num)%use_buffer) call error_mesg ('atmos_global_diag_mod', &
      'buffer not allocated in buffer_global_data for field='//trim(fields(field_num)%field_name), FATAL)

  is = 1; if (present(is_in)) is = is_in
  js = 1; if (present(js_in)) js = js_in
  ie = is + size(data,1) -1
  je = js + size(data,2) -1

  if (present(mask)) then
    if (.not.fields(field_num)%use_masking) call error_mesg &
        ('atmos_global_diag_mod', 'attempt to pass mask field without '// &
         'setting use_masking=True', FATAL)
    where (mask) fields(field_num)%buffer(is:ie,js:je) = data
  else
    fields(field_num)%buffer(is:ie,js:je) = data
  endif
  fields(field_num)%Time = Time

!-----------------------------------------------------------------------
  
end subroutine buffer_global_diag

!#######################################################################

function send_global_diag_data (field_num, data, Time, area, mask)
integer,           intent(in) :: field_num
real,              intent(in) :: data(:,:)
type(time_type),   intent(in) :: Time
real,    optional, intent(in) :: area(:,:)
logical, optional, intent(in) :: mask(:,:)

logical :: send_global_diag_data
!-----------------------------------------------------------------------

real, dimension(id,jd) :: data_gbl, area_gbl
real :: gbl_sum, area_sum

!-----------------------------------------------------------------------

  if (.not.module_is_initialized) call error_mesg &
      ('atmos_global_diag_mod', 'module has not been initialized', FATAL)

  if (field_num == 0) return

  if (field_num < 0 .or. field_num > num_fields) call error_mesg &
      ('atmos_global_diag_mod', 'invalid field number in send_global_diag_data', FATAL)

  if (fields(field_num)%use_buffer) call error_mesg ('atmos_global_diag_mod', &
       'buffer allocated in send_global_diag_data for field='//trim(fields(field_num)%field_name), FATAL)

  if (.not.present(area) .and. .not.present(mask)) then
    data_gbl = data*area_g
    area_sum = area_g_sum
    
  else if (present(area) .and. present(mask)) then
    where (mask)
      data_gbl = data*area
      area_gbl = area
    elsewhere
      data_gbl = 0.0
      area_gbl = 0.0
    endwhere
    area_sum = mpp_global_sum (Domain2, area_gbl, flags=BITWISE_EFP_SUM)

  else if (present(area) .and. .not.present(mask)) then
    data_gbl = data*area
    area_sum = mpp_global_sum (Domain2, area, flags=BITWISE_EFP_SUM)
 
  else if (.not.present(area) .and. present(mask)) then
    where (mask)
      data_gbl = data
      area_gbl = area_g
    elsewhere
      data_gbl = 0.0
      area_gbl = 0.0
    endwhere
    area_sum = mpp_global_sum (Domain2, area_gbl, flags=BITWISE_EFP_SUM)
  endif

  gbl_sum = mpp_global_sum (Domain2, data_gbl, flags=BITWISE_EFP_SUM)
  send_global_diag_data = send_data (fields(field_num)%field_id, gbl_sum/area_sum, Time)

!-----------------------------------------------------------------------

end function send_global_diag_data

!#######################################################################

function send_global_diag_buffer (field_num)
integer,         intent(in) :: field_num
logical :: send_global_diag_buffer
!-----------------------------------------------------------------------

real, dimension(id,jd) :: data_g, area_m
real :: gbl_sum, area_sum

!-----------------------------------------------------------------------

  if (.not.module_is_initialized) call error_mesg &
      ('atmos_global_diag_mod', 'module has not been initialized', FATAL)

  if (field_num == 0) return

  if (field_num < 0 .or. field_num > num_fields) call error_mesg &
      ('atmos_global_diag_mod', 'invalid field number in send_global_diag_buffer', FATAL)

  if (.not.fields(field_num)%use_buffer) call error_mesg ('atmos_global_diag_mod', &
      'buffer not allocated in send_global_diag_buffer for field='//trim(fields(field_num)%field_name), FATAL)

  ! weight data with area
  data_g = fields(field_num)%buffer*area_g
  area_sum = area_g_sum

  !--- when masking cells, need to compute new global area ---
  if (fields(field_num)%use_masking) then
    where(fields(field_num)%buffer .eq. missing_value)
      data_g = 0.0
      area_m = 0.0
    elsewhere
      area_m = area_g
    endwhere
    ! global sum of masked area
    area_sum = mpp_global_sum(Domain2, area_m, flags=BITWISE_EFP_SUM)
  endif

  ! sum of data*area
  gbl_sum = mpp_global_sum (Domain2, data_g, flags=BITWISE_EFP_SUM)

  send_global_diag_buffer = send_data (fields(field_num)%field_id, gbl_sum/area_sum, &
                                      fields(field_num)%Time)

!-----------------------------------------------------------------------

end function send_global_diag_buffer

!#######################################################################

subroutine atmos_global_diag_end

integer :: ind

!-----------------------------------------------------------------------

  if (.not.module_is_initialized) then
    call error_mesg ('atmos_global_diag_mod', &
          'module has not been initialized, '// &
          'nothing to do in atmos_global_diag_end', WARNING)
    return
  endif

  do ind = 1, num_fields
    if (fields(ind)%use_buffer) then
       deallocate(fields(ind)%buffer)
       fields(ind)%use_buffer = .false.
    endif
  enddo
  deallocate(fields)

  module_is_initialized = .false.

!-----------------------------------------------------------------------

end subroutine atmos_global_diag_end

!#######################################################################
! private routine

function find_field_index (field_name)
character(len=*), intent(in) :: field_name
integer :: find_field_index
integer :: ind

  find_field_index = 0

  do ind = 1, num_fields
    if (trim(field_name) .eq. trim(fields(ind)%field_name)) then
      find_field_index = ind
      return
    endif
  enddo

end function find_field_index

!#######################################################################

end module atmos_global_diag_mod

