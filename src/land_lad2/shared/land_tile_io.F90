module land_tile_io_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use mpp_mod, only : mpp_send, mpp_recv, mpp_sync
use mpp_mod, only : COMM_TAG_1,  COMM_TAG_2,  COMM_TAG_3,  COMM_TAG_4
use mpp_mod, only : COMM_TAG_5,  COMM_TAG_6,  COMM_TAG_7,  COMM_TAG_8
use mpp_mod, only : COMM_TAG_9,  COMM_TAG_10
use fms_mod, only : error_mesg, FATAL, NOTE, mpp_pe, get_mosaic_tile_file
use fms_io_mod, only : restart_file_type, free_restart_type, &
     get_instance_filename, read_data
use time_manager_mod, only : time_type
use data_override_mod, only : data_override_ug
use mpp_domains_mod,   only : mpp_pass_SG_to_UG
use nf_utils_mod, only : nfu_inq_dim, nfu_inq_var, nfu_def_dim, nfu_def_var, &
     nfu_get_var, nfu_put_var, nfu_put_att
use land_io_mod, only : print_netcdf_error, read_field, input_buf_size, new_land_io
use land_tile_mod, only : land_tile_type, land_tile_list_type, land_tile_enum_type, &
     first_elmt, loop_over_tiles, &
     tile_exists_func, fptr_i0, fptr_i0i, fptr_r0, fptr_r0i, fptr_r0ij, fptr_r0ijk, &
     land_tile_map

use land_data_mod, only  : lnd_sg, lnd
use land_utils_mod, only : put_to_tiles_r0d_fptr

use fms_io_mod, only: fms_io_unstructured_save_restart
use fms_io_mod, only: fms_io_unstructured_register_restart_axis
use fms_io_mod, only: fms_io_unstructured_register_restart_field
use fms_io_mod, only: CIDX,ZIDX,CCIDX
use fms_io_mod, only: fms_io_unstructured_get_field_size
use fms_io_mod, only: fms_io_unstructured_read
use fms_io_mod, only: fms_io_unstructured_field_exist

implicit none
private

! ==== public interfaces =====================================================
! restart i/o subroutines: those use CF "compression by gathering" technique
! to pack tile data.
public :: land_restart_type
public :: init_land_restart, open_land_restart, save_land_restart, free_land_restart
public :: add_restart_axis
public :: add_tile_data, add_int_tile_data, add_scalar_data
public :: get_tile_data, get_int_tile_data, get_scalar_data
public :: field_exists
public :: gather_tile_index

public :: read_field

public :: create_tile_out_file

! auxiliary subroutines
public :: get_tile_by_idx
public :: print_netcdf_error

! ==== end of public interfaces ==============================================
interface create_tile_out_file
   module procedure create_tile_out_file_idx_old
   module procedure create_tile_out_file_idx_new
end interface

interface add_tile_data
   module procedure add_tile_data_r0d_fptr_r0
   module procedure add_tile_data_r0d_fptr_r0i
   module procedure add_tile_data_r1d_fptr_r0i
   module procedure add_tile_data_r1d_fptr_r0ij
   module procedure add_tile_data_r2d_fptr_r0ij
   module procedure add_tile_data_r2d_fptr_r0ijk
end interface

interface add_int_tile_data
   module procedure add_tile_data_i0d_fptr_i0
   module procedure add_tile_data_i1d_fptr_i0i
end interface

interface get_tile_data
   module procedure get_tile_data_r0d_fptr_r0
   module procedure get_tile_data_r0d_fptr_r0i
   module procedure get_tile_data_r1d_fptr_r0i
   module procedure get_tile_data_r1d_fptr_r0ij
   module procedure get_tile_data_r2d_fptr_r0ij
   module procedure get_tile_data_r2d_fptr_r0ijk
end interface

interface get_int_tile_data
   module procedure get_tile_data_i0d_fptr_i0
   module procedure get_tile_data_i1d_fptr_i0i
end interface

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'land_tile_io_mod'

! name of the "compressed" dimension (and dimension variable) in the output
! netcdf files -- that is, the dimensions written out using compression by
! gathering, as described in CF conventions. See subroutines write_tile_data,
! read_tile_data, read_unpack_tile_data, write_cohort_data
character(len=*), parameter :: tile_index_name = 'tile_index'

! ==== module types ==========================================================
type axis
   character(128) :: name = ''  ! name of the axis
   integer        :: len = -1   ! length of the axis
end type axis
! land restart type encapsulates the data needed for the land restarts
type land_restart_type
   type(restart_file_type) :: rhandle ! fms_io restart file data type
   logical :: should_free_rhandle = .FALSE.

   character(267) :: basename ='' ! name of the restart file
   character(267) :: filename ='' ! name of the restart file after adding PE number and such
   integer :: ncid = -1 ! netcdf id, used only for old land io
   integer, allocatable :: tidx(:) ! tile index
   integer, allocatable :: cidx(:) ! vegetation cohort index
   integer :: tile_dim_length = -1! length of tile dimension

   ! axis information
   integer    :: nax=0
   type(axis) :: ax(5)
end type land_restart_type

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ==============================================================================
! given a file name and tile indicator funcrion, creates land restart structure
subroutine init_land_restart(restart,filename,tile_exists,tile_dim_length)
  type(land_restart_type), intent(out) :: restart
  character(*),            intent(in)  :: filename ! name of the file to create
  procedure(tile_exists_func) :: tile_exists ! existence detector function:
      ! returns true if specific tile exists (hence should be written to restart)
  integer,                 intent(in)  :: tile_dim_length ! length of the tile dimension
      ! in the output file

  restart%basename=filename
  ! TODO: determine tile_dim_length inside this subroutine. It is equal to the
  ! max number of tiles per grid cell
  restart%tile_dim_length = tile_dim_length
  ! allocate and fill tile compression index
  call gather_tile_index(tile_exists,restart%tidx)

  if (new_land_io) then
     call create_tile_out_file_idx_new(restart%rhandle,restart%basename,restart%tidx, &
          restart%tile_dim_length)
     restart%should_free_rhandle = .TRUE.
  else
     call create_tile_out_file_idx_old(restart%ncid,'RESTART/'//trim(restart%basename), &
          restart%tidx, restart%tile_dim_length, lnd%coord_glon, lnd%coord_glat)
  endif
end subroutine init_land_restart

! ==============================================================================
subroutine open_land_restart(restart,filename,restart_exists)
  type(land_restart_type), intent(out) :: restart
  character(*),            intent(in)  :: filename
  logical,                 intent(out) :: restart_exists

  ! ---- local vars
  integer :: flen(4) ! length of the index
  logical :: found   ! true if field exists
  integer :: ierr

  restart%basename = filename
  call get_input_restart_name(restart%basename,restart_exists,restart%filename)
  if (.not.restart_exists) return

  if (new_land_io) then
    !Get the size of the tile dimension from the file.
     call fms_io_unstructured_get_field_size(filename, "tile", flen, lnd%domain, &
                                             field_found=found)
     if (.not. found) then
         call error_mesg("open_land_restart", "dimension 'tile' not found in file '" &
                         //trim(filename)//"'.", FATAL)
     endif
     restart%tile_dim_length = flen(1)

    !Get the size of the tile index dimension from the file.
     call fms_io_unstructured_get_field_size(filename, "tile_index", flen, lnd%domain, &
                                             field_found=found)
     if (.not. found) then
         call error_mesg("open_land_restart", "'tile_index' not found in file '" &
                         //trim(filename)//"'.", FATAL)
     endif
     allocate(restart%tidx(flen(1)))

     ! Read in the tile_index field from the file.
     call fms_io_unstructured_read(filename, "tile_index", restart%tidx, lnd%domain, &
                        timelevel=1)

     ! Get the size of the cohort_index dimension from the file.
     call fms_io_unstructured_get_field_size(restart%basename, "cohort_index", flen, &
                        lnd%domain, field_found=found)
     if (found) then
        ! Read in the cohort_index field from the file.
        allocate(restart%cidx(flen(1)))
        call fms_io_unstructured_read(restart%basename, "cohort_index", restart%cidx, &
                                      lnd%domain, timelevel=1)
     endif
     ! TODO: possibly make tile index and cohort index names parameters in this module
     !       just constants, no sense to make them namelists vars
  else ! old i/o
      __NF_ASRT__(nf_open(restart%filename,NF_NOWRITE,restart%ncid))
      ierr = nfu_inq_dim(restart%ncid,'tile',len=restart%tile_dim_length)
      if (ierr/=NF_NOERR) call error_mesg('open_land_restart', &
           'dimension "tile" not found in file "'//trim(filename)//'"', FATAL)
  endif
end subroutine open_land_restart

! ==============================================================================
subroutine save_land_restart(restart)
  type(land_restart_type), intent(inout) :: restart

  if (restart%should_free_rhandle) then
       call fms_io_unstructured_save_restart(restart%rhandle)
  endif
end subroutine save_land_restart

! ==============================================================================
! restet land restart data type to its initial state
subroutine free_land_restart(restart)
  type(land_restart_type), intent(inout) :: restart

  if (restart%should_free_rhandle) call free_restart_type(restart%rhandle)
  restart%should_free_rhandle = .FALSE.
  if (restart%ncid>0) then
     __NF_ASRT__(nf_close(restart%ncid))
     restart%ncid = -1
  endif
  restart%basename = ''
  restart%filename = ''
  if (allocated(restart%tidx)) deallocate(restart%tidx)
  if (allocated(restart%cidx)) deallocate(restart%cidx)
  restart%tile_dim_length = -1
end subroutine free_land_restart

! ==============================================================================
subroutine add_restart_axis(restart,name,data,cartesian,units,longname,sense)
  type(land_restart_type),    intent(inout) :: restart
  character(len=*),           intent(in)    :: name
  real,                       intent(in)    :: data(:)
  character(len=*),           intent(in)    :: cartesian
  character(len=*), optional, intent(in)    :: units, longname
  integer,          optional, intent(in)    :: sense
  ! how to get rid of "cartesian" attribute? In some cases (carbon cohort) it is far from obvious what it should be.

  integer :: n
  real, pointer :: data_(:)

  if (new_land_io) then
     allocate(data_(size(data)))
     data_(:) = data(:)
     call fms_io_unstructured_register_restart_axis(restart%rhandle, restart%basename, &
               name, data_, cartesian, lnd%domain, units=units, longname=longname, sense=sense)
  else
     if (mpp_pe()==lnd%io_pelist(1)) then
        __NF_ASRT__(nfu_def_dim(restart%ncid,name,data(:),longname,units))
        if (present(sense)) then
           if (sense<0) then
              __NF_ASRT__(nfu_put_att(restart%ncid,name,'positive','down'))
           endif
        endif
     endif
  endif
  ! record dimension information for future use
  n = restart%nax+1; restart%nax = n
  restart%ax(n)%name = name
  restart%ax(n)%len  = size(data)
end subroutine add_restart_axis

! ==============================================================================
logical function field_exists(restart,name)
  type(land_restart_type), intent(in) :: restart
  character(len=*),        intent(in) :: name

  if (new_land_io) then
     field_exists = fms_io_unstructured_field_exist(restart%basename, name, domain=lnd%domain)
  else
     field_exists = (nfu_inq_var(restart%ncid,trim(name))==NF_NOERR)
  endif
end function field_exists

! ==============================================================================
subroutine add_scalar_data(restart,varname,datum,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  integer,          intent(in) :: datum
  character(len=*), intent(in), optional :: units, longname

  integer :: id_restart, ierr

  if (new_land_io) then
     id_restart = fms_io_unstructured_register_restart_field(restart%rhandle, &
            restart%basename, varname, datum, lnd%domain, longname=longname, units=units)
  else
     if(mpp_pe()==lnd%io_pelist(1)) then
        ierr = nf_redef(restart%ncid)
        __NF_ASRT__(nfu_def_var(restart%ncid,varname,NF_INT,long_name=longname,units=units))
        ierr = nf_enddef(restart%ncid)
        __NF_ASRT__(nfu_put_var(restart%ncid,varname,datum))
     end if
  endif
end subroutine add_scalar_data

subroutine add_tile_data_i0d_fptr_i0(restart,varname,fptr,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(fptr_i0)           :: fptr ! subroutine returning pointer to the data
  character(len=*), intent(in), optional :: units, longname

  integer :: id_restart
  integer, pointer :: data(:)

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_r0d_fptr_r0', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)
  allocate(data(size(restart%tidx)))
  call gather_tile_data_i0d(fptr,restart%tidx,data)
  if (new_land_io) then
     id_restart = fms_io_unstructured_register_restart_field(restart%rhandle, restart%basename, &
            varname, data, (/CIDX/), lnd%domain, longname=longname, units=units, &
            restart_owns_data=.true.)
  else ! old land io
     call write_tile_data_i1d(restart%ncid,varname,data,longname,units)
     deallocate(data)
  endif
end subroutine add_tile_data_i0d_fptr_i0

subroutine add_tile_data_r0d_fptr_r0(restart,varname,fptr,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(fptr_r0)           :: fptr ! subroutine returning pointer to the data
  character(len=*), intent(in), optional :: units, longname

  integer :: id_restart
  real, pointer :: data(:)

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_r0d_fptr_r0', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)
  allocate(data(size(restart%tidx)))
  call gather_tile_data_r0d(fptr,restart%tidx,data)
  if (new_land_io) then
     id_restart = fms_io_unstructured_register_restart_field(restart%rhandle, restart%basename, &
            varname, data, (/CIDX/), lnd%domain, longname=longname, units=units, &
            restart_owns_data=.true.)
  else ! old land io
     call write_tile_data_r1d(restart%ncid,varname,data,longname,units)
     deallocate(data)
  endif
end subroutine add_tile_data_r0d_fptr_r0

subroutine add_tile_data_r0d_fptr_r0i(restart,varname,fptr,index,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(fptr_r0i)          :: fptr ! subroutine returning pointer to the data
  integer ,         intent(in) :: index ! index of the fptr array element to write
  character(len=*), intent(in), optional :: units, longname

  integer :: id_restart
  real, pointer :: data(:)

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_r0d_fptr_r0', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)
  allocate(data(size(restart%tidx)))
  call gather_tile_data_r0d_idx(fptr,index,restart%tidx,data)
  if (new_land_io) then
     id_restart = fms_io_unstructured_register_restart_field(restart%rhandle, restart%basename, &
            varname, data, (/CIDX/), lnd%domain, longname=longname, units=units, &
            restart_owns_data=.true.)
  else ! old land io
     call write_tile_data_r1d(restart%ncid,varname,data,longname,units)
     deallocate(data)
  endif
end subroutine add_tile_data_r0d_fptr_r0i


! given restart and name of the dimension, returns the size of this dimension
integer function dimlen(restart,dimname)
  type(land_restart_type), intent(in) :: restart ! restart data structure
  character(*),            intent(in) :: dimname ! name of the dimansion

  integer :: i
  dimlen = -1
  do i = 1,restart%nax
     if (trim(restart%ax(i)%name) == trim(dimname)) then
        dimlen=restart%ax(i)%len
        return
     endif
  enddo
  if (dimlen<1) call error_mesg('dimlen', 'axis "'//trim(dimname)//'" not found', FATAL)
end


subroutine add_tile_data_i1d_fptr_i0i(restart,varname,zdim,fptr,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: zdim      ! name of the z-dimension
  procedure(fptr_i0i)          :: fptr    ! subroutine returning pointer to the data
  character(len=*), intent(in), optional :: units, longname

  integer :: id_restart
  integer, pointer :: data(:,:) ! needs to be pointer; we are passing ownership to restart object
  integer :: i,nlev

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_i1d_fptr_i0i', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)

  nlev = dimlen(restart,zdim)
  allocate(data(size(restart%tidx),nlev))
  call gather_tile_data_i1d(fptr,restart%tidx,data)

  if (new_land_io) then
     id_restart = fms_io_unstructured_register_restart_field(restart%rhandle, restart%basename, &
            varname, data, (/CIDX,ZIDX/), lnd%domain, longname=longname, units=units, &
            restart_owns_data=.true.)
  else ! old land io
     call write_tile_data_i2d(restart%ncid,varname,data,zdim,longname,units)
     deallocate(data)
  endif
end subroutine add_tile_data_i1d_fptr_i0i

subroutine add_tile_data_r1d_fptr_r0i(restart,varname,zdim,fptr,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: zdim      ! name of the z-dimension
  procedure(fptr_r0i)          :: fptr    ! subroutine returning pointer to the data
  character(len=*), intent(in), optional :: units, longname

  integer :: id_restart
  real, pointer :: data(:,:) ! needs to be pointer; we are passing ownership to restart object
  integer :: i,nlev

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_r0d_fptr_r0i', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)

  nlev = dimlen(restart,zdim)
  allocate(data(size(restart%tidx),nlev))
  call gather_tile_data_r1d(fptr,restart%tidx,data)

  if (new_land_io) then
     ! checking name of the dimension here is a dirty trick, which is sure to
     ! bite us in the future, but it is necessary because fms_io has no way to
     ! figure out what the additional dimension of the variable is. A better way
     ! to fix that is to rewrite fms_io so that it allows to specify the dimensions
     ! of the variable in a sane way.
     if (trim(zdim)=='soilCCohort') then
        ! write (*,*) 'writing "',trim(varname),'" with C_CC'
        id_restart = fms_io_unstructured_register_restart_field(restart%rhandle, restart%basename, &
              varname, data, (/CIDX,CCIDX/), lnd%domain, longname=longname, units=units, &
              restart_owns_data=.true.)
     else
        id_restart = fms_io_unstructured_register_restart_field(restart%rhandle, restart%basename, &
              varname, data, (/CIDX,ZIDX/), lnd%domain, longname=longname, units=units, &
              restart_owns_data=.true.)
     endif
  else ! old land io
     call write_tile_data_r2d(restart%ncid,varname,data,zdim,longname,units)
     deallocate(data)
  endif
end subroutine add_tile_data_r1d_fptr_r0i

subroutine add_tile_data_r1d_fptr_r0ij(restart,varname,zdim,fptr,index,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: zdim      ! name of the z-dimension
  procedure(fptr_r0ij)         :: fptr    ! subroutine returning pointer to the data
  integer         , intent(in) :: index   ! index of the array element to write
  character(len=*), intent(in), optional :: units, longname

  integer :: id_restart
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real, pointer :: data(:,:) ! needs to be pointer; we are passing ownership to restart object
  real, pointer :: ptr ! pointer to the tile data
  integer :: i,n,nlev

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_r0d_fptr_r0i', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)

  nlev = dimlen(restart,zdim)
  allocate(data(size(restart%tidx),nlev))
  data = NF_FILL_DOUBLE

  ! gather data into an array along the tile dimension. It is assumed that
  ! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(restart%tidx)
     call get_tile_by_idx(restart%tidx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                          lnd%ls, lnd%gs, lnd%ge, tileptr)
     do n = 1,nlev
        call fptr(tileptr,n,index, ptr)
        if(associated(ptr)) then
           data(i,n) = ptr
        endif
     enddo
  enddo
  if (new_land_io) then
     ! checking name of the dimension here is a dirty trick, which is sure to
     ! bite us in the future, but it is necessary because fms_io has no way to
     ! figure out what the additional dimension of the variable is. A better way
     ! to fix that is to rewrite fms_io so that it allows to specify the dimensions
     ! of the variable in a sane way.
     if (trim(zdim)=='soilCCohort') then
        ! write (*,*) 'writing "',trim(varname),'" with C_CC'
        id_restart = fms_io_unstructured_register_restart_field(restart%rhandle, restart%basename, &
                varname, data, (/CIDX,CCIDX/), lnd%domain, longname=longname, units=units, &
                restart_owns_data=.true.)
     else
        id_restart = fms_io_unstructured_register_restart_field(restart%rhandle, restart%basename, &
                varname, data, (/CIDX,ZIDX/), lnd%domain, longname=longname, units=units, &
                restart_owns_data=.true.)
     endif
  else ! old land io
     call write_tile_data_r2d(restart%ncid,varname,data,zdim,longname,units)
     deallocate(data)
  endif
end subroutine add_tile_data_r1d_fptr_r0ij

subroutine add_tile_data_r2d_fptr_r0ij(restart,varname,dim1,dim2,fptr,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: dim1,dim2 ! names of extra dimensions
  procedure(fptr_r0ij)         :: fptr    ! subroutine returning pointer to the data
  character(len=*), intent(in), optional :: units, longname

  integer :: id_restart
  real, pointer :: data(:,:,:) ! needs to be pointer; we are passing ownership to restart object
  integer :: i,dim1len,dim2len

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_r2d_fptr_r0ijk', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)

  dim1len = dimlen(restart,dim1)
  dim2len = dimlen(restart,dim2)
  allocate(data(size(restart%tidx),dim1len,dim2len))
  call gather_tile_data_r2d(fptr,restart%tidx,data)

  if (new_land_io) then
     id_restart = fms_io_unstructured_register_restart_field(restart%rhandle, restart%basename, &
            varname, data, (/CIDX,ZIDX,CCIDX/), lnd%domain, longname=longname, units=units, &
            restart_owns_data=.true.)
  else ! old land io
     call write_tile_data_r3d(restart%ncid,varname,data,dim1,dim2,longname,units)
     deallocate(data)
  endif
end subroutine add_tile_data_r2d_fptr_r0ij

subroutine add_tile_data_r2d_fptr_r0ijk(restart,varname,dim1,dim2,fptr,index,longname,units)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: dim1,dim2 ! names of extra dimensions
  procedure(fptr_r0ijk)        :: fptr    ! subroutine returning pointer to the data
  integer         , intent(in) :: index   ! index of the array element to write
  character(len=*), intent(in), optional :: units, longname

  integer :: id_restart
  real, pointer :: data(:,:,:) ! needs to be pointer; we are passing ownership to restart object
  integer :: i,dim1len,dim2len

  if (.not.allocated(restart%tidx)) call error_mesg('add_tile_data_r2d_fptr_r0ijk', &
        'tidx not allocated: looks like land restart was not initialized',FATAL)

  dim1len = dimlen(restart,dim1)
  dim2len = dimlen(restart,dim2)
  allocate(data(size(restart%tidx),dim1len,dim2len))
  call gather_tile_data_r2d_idx(fptr,index,restart%tidx,data)

  if (new_land_io) then
     id_restart = fms_io_unstructured_register_restart_field(restart%rhandle, restart%basename, &
            varname, data, (/CIDX,ZIDX,CCIDX/), lnd%domain, longname=longname, units=units, &
            restart_owns_data=.true.)
  else ! old land io
     call write_tile_data_r3d(restart%ncid,varname,data,dim1,dim2,longname,units)
     deallocate(data)
  endif
end subroutine add_tile_data_r2d_fptr_r0ijk

! =============================================================================
subroutine get_scalar_data(restart,varname,datum)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  integer,          intent(out) :: datum

  if (new_land_io) then
     call read_data(restart%basename,varname,datum,domain=lnd_sg%domain)
  else
     __NF_ASRT__(nfu_get_var(restart%ncid,varname,datum))
  endif
end subroutine get_scalar_data

subroutine get_tile_data_i0d_fptr_i0(restart,varname,fptr)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(fptr_i0)           :: fptr    ! subroutine returning pointer to the data

  ! ---- local vars
  integer, allocatable :: r(:) ! input data buffer
  logical :: found

  if (new_land_io) then
     allocate(r(size(restart%tidx)))
     call fms_io_unstructured_read(restart%basename, varname, r, lnd%domain, timelevel=1)
     call assemble_tiles_i0d(fptr,restart%tidx,r)
     deallocate(r)
  else ! old land io
     call read_tile_data_i0d_fptr(restart%ncid,varname,fptr)
  endif
end subroutine get_tile_data_i0d_fptr_i0

subroutine get_tile_data_r0d_fptr_r0(restart,varname,fptr)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(fptr_r0)           :: fptr    ! subroutine returning pointer to the data

  ! ---- local vars
  real, allocatable :: r(:) ! input data buffer
  logical :: found

  if (new_land_io) then
     allocate(r(size(restart%tidx)))
     call fms_io_unstructured_read(restart%basename, varname, r, lnd%domain, timelevel=1)
     call assemble_tiles_r0d(fptr,restart%tidx,r)
     deallocate(r)
  else ! old land io
     call read_tile_data_r0d_fptr_r0(restart%ncid,varname,fptr)
  endif
end subroutine get_tile_data_r0d_fptr_r0

subroutine get_tile_data_r0d_fptr_r0i(restart,varname,fptr,index)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  procedure(fptr_r0i)          :: fptr    ! subroutine returning pointer to the data
  integer,          intent(in) :: index   ! index where to read the data

  ! ---- local vars
  real, allocatable :: r(:) ! input data buffer
  logical :: found

  if (new_land_io) then
     allocate(r(size(restart%tidx)))
     call fms_io_unstructured_read(restart%basename, varname, r, lnd%domain, timelevel=1)
     call assemble_tiles_r0d_idx(fptr,index,restart%tidx,r)
     deallocate(r)
  else ! old land io
     call read_tile_data_r0d_fptr_r0i(restart%ncid,varname,fptr,index)
  endif
end subroutine get_tile_data_r0d_fptr_r0i

subroutine get_tile_data_r1d_fptr_r0i(restart,varname,zdim,fptr)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: zdim ! name of the z-dimension
  procedure(fptr_r0i)          :: fptr    ! subroutine returning pointer to the data

  ! ---- local vars
  integer :: flen(4) ! size of the input field
  real, allocatable :: r(:,:) ! input data buffer
  logical :: found

  if (new_land_io) then
    !Get the size of z-dimension from the file.
     call fms_io_unstructured_get_field_size(restart%basename, zdim, flen, lnd%domain, &
                                             field_found=found)
     if (.not. found) then
         call error_mesg("get_tile_data_r0d_fptr_r0i", &
               "axis '"//trim(zdim)//"' was not found in file '"//trim(restart%basename)//"'.", &
               FATAL)
     endif

    !Read in the field from the file.
     allocate(r(size(restart%tidx),flen(1)))
     call fms_io_unstructured_read(restart%basename, varname, r, lnd%domain, timelevel=1)
     call assemble_tiles_r1d(fptr,restart%tidx,r)
     deallocate(r)
  else ! old land io
     call read_tile_data_r1d_fptr_r0i(restart%ncid,varname,fptr)
  endif
end subroutine get_tile_data_r1d_fptr_r0i

subroutine get_tile_data_i1d_fptr_i0i(restart,varname,zdim,fptr)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: zdim ! name of the z-dimension
  procedure(fptr_i0i)          :: fptr    ! subroutine returning pointer to the data

  ! ---- local vars
  integer :: flen(4) ! size of the input field
  integer, allocatable :: r(:,:) ! input data buffer
  logical :: found

  if (new_land_io) then
    !Get the size of z-dimension from the file.
     call fms_io_unstructured_get_field_size(restart%basename, zdim, flen, lnd%domain, &
                                             field_found=found)
     if (.not. found) then
         call error_mesg("get_tile_data_i1d_fptr_i0i", &
                         "axis '"//trim(zdim)//"' was not found in file '"//trim(restart%basename)//"'.", &
                         FATAL)
     endif

    !Read in the field data from the file.
     allocate(r(size(restart%tidx),flen(1)))
     call fms_io_unstructured_read(restart%basename, varname, r, lnd%domain, timelevel=1)
     call assemble_tiles_i1d(fptr,restart%tidx,r)
     deallocate(r)
  else ! old land io
     call read_tile_data_i1d_fptr_i0i(restart%ncid,varname,fptr)
  endif
end subroutine get_tile_data_i1d_fptr_i0i

subroutine get_tile_data_r1d_fptr_r0ij(restart,varname,zdim,fptr,index)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: zdim ! name of the z-dimension
  procedure(fptr_r0ij)         :: fptr    ! subroutine returning pointer to the data
  integer ,         intent(in) :: index

  ! ---- local vars
  integer :: flen(4) ! size of the input field
  real, allocatable :: r(:,:) ! input data buffer
  logical :: found

  if (new_land_io) then
    !Get the size of z-dimension from the file.
     call fms_io_unstructured_get_field_size(restart%basename, zdim, flen, lnd%domain, &
                                             field_found=found)
     if (.not. found) then
         call error_mesg("get_tile_data_r1d_fptr_r0ij", &
                         "axis '"//trim(zdim)//"' was not found in file '"//trim(restart%basename)//"'.", &
                         FATAL)
     endif

    !Read in the field from the file.
     allocate(r(size(restart%tidx),flen(1)))
     call fms_io_unstructured_read(restart%basename, varname, r, lnd%domain, timelevel=1)
     call assemble_tiles_r1d_idx(fptr,index,restart%tidx,r)
     deallocate(r)
  else ! old land io
     call read_tile_data_r1d_fptr_r0ij(restart%ncid,varname,fptr,index)
  endif
end subroutine get_tile_data_r1d_fptr_r0ij

subroutine get_tile_data_r2d_fptr_r0ij(restart,varname,dim1,dim2,fptr)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: dim1,dim2 ! names of the dimensions
  procedure(fptr_r0ij)         :: fptr    ! subroutine returning pointer to the data

  ! ---- local vars
  integer :: flen(4),n,m ! size of the input field
  real, allocatable :: r(:,:,:) ! input data buffer
  logical :: found

  if (new_land_io) then
    !Get the size of the first dimension of the field.
     call fms_io_unstructured_get_field_size(restart%basename, dim1, flen, lnd%domain, &
                                             field_found=found)
     if (.not. found) then
         call error_mesg("get_tile_data_r0d_fptr_r0i", &
              "axis '"//trim(dim1)//"' was not found in file '"//trim(restart%basename)//"'.", &
              FATAL)
     endif
     n = flen(1)

    !Get the size of the second dimension of the field.
     call fms_io_unstructured_get_field_size(restart%basename, dim2, flen, lnd%domain, &
                                             field_found=found)
     if (.not. found) then
         call error_mesg("get_tile_data_r0d_fptr_r0i", &
              "axis '"//trim(dim2)//"' was not found in file '"//trim(restart%basename)//"'.", &
              FATAL)
     endif
     m = flen(1)

    !Read in the field data from the file.
     allocate(r(size(restart%tidx),n,m))
     call fms_io_unstructured_read(restart%basename, varname, r, lnd%domain, timelevel=1)
     call assemble_tiles_r2d(fptr,restart%tidx,r)
     deallocate(r)
  else ! old land io
     call read_tile_data_r2d_fptr_r0ij(restart%ncid,varname,fptr)
  endif
end subroutine get_tile_data_r2d_fptr_r0ij

subroutine get_tile_data_r2d_fptr_r0ijk(restart,varname,dim1,dim2,fptr,index)
  type(land_restart_type), intent(inout) :: restart
  character(len=*), intent(in) :: varname ! name of the variable to write
  character(len=*), intent(in) :: dim1,dim2 ! names of the extra dimensions
  procedure(fptr_r0ijk)        :: fptr    ! subroutine returning pointer to the data
  integer,          intent(in) :: index   ! index where to read the data

  ! ---- local vars
  integer :: flen(4), m, n ! size of the input field
  real, allocatable :: r(:,:,:) ! input data buffer
  logical :: found

  if (new_land_io) then
    !Get the size of the z-dimension from the file.
     call fms_io_unstructured_get_field_size(restart%basename, dim1, flen, lnd%domain, &
                                             field_found=found)
     if (.not. found) then
         call error_mesg("get_tile_data_r0d_fptr_r0i", &
              "axis '"//trim(dim1)//"' was not found in file '"//trim(restart%basename)//"'.", &
              FATAL)
     endif
     n = flen(1)

    !Get the size of the 3rd dimension of the field.
     call fms_io_unstructured_get_field_size(restart%basename, dim2, flen, lnd%domain, &
                                             field_found=found)
     if (.not. found) then
         call error_mesg("get_tile_data_r0d_fptr_r0i", &
              "axis '"//trim(dim1)//"' was not found in file '"//trim(restart%basename)//"'.", &
              FATAL)
     endif
     m = flen(1)

    !Read in the field from the file.
     allocate(r(size(restart%tidx),n,m))
     call fms_io_unstructured_read(restart%basename, varname, r, lnd%domain, timelevel=1)
     call assemble_tiles_r2d_idx(fptr,index,restart%tidx,r)
     deallocate(r)
  else ! old land io
     call read_tile_data_r2d_fptr_r0ijk(restart%ncid,varname,fptr,index)
  endif
end subroutine get_tile_data_r2d_fptr_r0ijk

! =============================================================================
! given a generic name of the restart file, checks if a file with one of the
! possible restarts file names exists, and if it does returns the tile-qualified
! (or tile- and processor-qualified) name of the restart.
subroutine get_input_restart_name(name, restart_exists, actual_name, new_land_io)
  character(*), intent(in)  :: name        ! "generic" name of the restart
  logical     , intent(out) :: restart_exists ! TRUE if any file found
  character(*), intent(out) :: actual_name ! name of the found file, if any
  logical, intent(in), optional :: new_land_io

  ! ---- local vars
  character(6) :: PE_suffix ! PE number
  character(len=256) :: distributed_name

  ! Build the restart file name.
  call get_instance_filename(trim(name), actual_name)
  call get_mosaic_tile_file(trim(actual_name),actual_name, lnd%domain)
  ! we can't use fms file_exist function here, because it lies: it checks not
  ! just the original name, but the name with PE suffix, and returns true if
  ! either of those exist
  inquire (file=trim(actual_name), exist=restart_exists)
  if (.not.restart_exists) then
     ! try the name with current PE number attached
     write(PE_suffix,'(".",I4.4)') lnd%io_id
     distributed_name = trim(actual_name)//trim(PE_suffix)
     inquire (file=trim(distributed_name), exist=restart_exists)
     if(present(new_land_io)) then
       if(.not.new_land_io) actual_name = trim(distributed_name)
     else
       ! if new_land_io is not present then revert to behavior of previous revision. That is, as if new_land_io=.false.
       actual_name = trim(distributed_name)
     endif
  endif

end subroutine get_input_restart_name

! =============================================================================
! this subroutine creates netcdf file for output of tiled data using "compression
! by gathering," as described in CF conventions, and creates coordinate system
! necessary for write_tile_data subroutines.
! In particular:
!  "compressed" dimension and integer variable with appropriate attributes, with
! variable filled with packing indices
!   horizontal dimensions "lat" and "lon" with associated variable, and boundaries,
! describing global grid
!   dimension "tile," without any associated variable. length of this dimension is
! equal to current global max of number of tiles per grid cell
!
! The file is actually created only by root processor of our io_domain; the rest
! of the processors just open the created file in NOWRITE mode.
subroutine create_tile_out_file_idx_old(ncid, name, tidx, tile_dim_length, glon, glat, reserve)
  integer          , intent(out) :: ncid      ! resulting NetCDF id
  character(len=*) , intent(in)  :: name      ! name of the file to create
  real             , intent(in)  :: glon(:)   ! longitudes of the grid centers
  real             , intent(in)  :: glat(:)   ! latitudes of the grid centers
  integer          , intent(in)  :: tidx(:)   ! integer compressed index of tiles
  integer          , intent(in)  :: tile_dim_length ! length of tile axis
  integer, optional, intent(in)  :: reserve   ! amount of space to reserve for
  ! header expansion. This subroutine and following calls to write_tile_data
  ! will work even if this is set to 0, but for efficiency it is useful to
  ! specify some non-zero value, so that netcdf library do not have to rewrite
  ! entire file each time a new variable is added. Default value (8K) should do
  ! fine in most cases.

  ! ---- local vars
  integer        :: reserve_  ! local value of space to reserve at the end of NetCDF header
  character(256) :: full_name ! full name of the file, including the processor number
  character(6)   :: PE_suffix ! PE number
  integer, allocatable :: ntiles(:) ! list of land tile numbers for each of PEs in io_domain
  integer, allocatable :: tidx2(:)  ! array of tile indices from all PEs in io_domain
  integer :: p ! io_domain PE iterator
  integer :: k ! current index in tidx2 array for receive operation
  integer :: i
  integer :: iret

  ! form the full name of the file
  call get_instance_filename(trim(name), full_name)
  call get_mosaic_tile_file(trim(full_name),full_name,lnd%domain)
  if (lnd%append_io_id) then
      write(PE_suffix,'(".",I4.4)') lnd%io_id
  else
      PE_suffix = ''
  endif

  full_name = trim(full_name)//trim(PE_suffix)
  if(tile_dim_length<=0) &
       call error_mesg('create_tile_out_file','tile axis length must be positive', FATAL)

  if (mpp_pe()/=lnd%io_pelist(1)) then
     ! if current PE doesn't do io, we just send the data to the processor that
     ! does
     call mpp_send(size(tidx), plen=1,          to_pe=lnd%io_pelist(1), tag=COMM_TAG_1)
     call mpp_send(tidx(1),    plen=size(tidx), to_pe=lnd%io_pelist(1), tag=COMM_TAG_2)
  else
     ! gather an array of tile sizes from all processors in our io_domain
     allocate(ntiles(size(lnd%io_pelist)))
     ntiles(1) = size(tidx)
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(ntiles(p), from_pe=lnd%io_pelist(p), glen=1, tag=COMM_TAG_1)
     enddo
     ! gather tile indices from all processors in our io_domain
     allocate(tidx2(sum(ntiles(:))))
     tidx2(1:ntiles(1))=tidx(:)
     k=ntiles(1)+1
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(tidx2(k), from_pe=lnd%io_pelist(p), glen=ntiles(p), tag=COMM_TAG_2)
        k = k+ntiles(p)
     enddo
     ! create netcdf file
#ifdef use_netCDF3
     __NF_ASRT__(nf_create(full_name,NF_CLOBBER,ncid))
#elif use_LARGEFILE
     __NF_ASRT__(nf_create(full_name,ior(NF_64BIT_OFFSET,NF_CLOBBER),ncid))
#else
     __NF_ASRT__(nf_create(full_name,ior(NF_NETCDF4,NF_CLASSIC_MODEL),ncid))
#endif

     ! create lon, lat, dimensions and variables
     __NF_ASRT__(nfu_def_dim(ncid,'lon' ,glon(:) ,'longitude','degrees_east'))
     __NF_ASRT__(nfu_def_dim(ncid,'lat' ,glat(:) ,'latitude','degrees_north'))

     iret=nfu_def_dim(ncid,'tile',(/(p,p=1,tile_dim_length)/),'tile number within grid cell')
     __NF_ASRT__(iret)
     ! the size of tile dimension really does not matter for the output, but it does
     ! matter for uncompressing utility, since it uses it as a size of the array to
     ! unpack to
     ! create tile index dimension and variable
     __NF_ASRT__(nfu_def_dim(ncid,tile_index_name,tidx2,'compressed land point index'))
     __NF_ASRT__(nfu_put_att(ncid,tile_index_name,'compress','tile lat lon'))
     __NF_ASRT__(nfu_put_att(ncid,tile_index_name,'valid_min',0))

     ! determine the local value of space reserved in the header; by default 16K
     reserve_ = 1024*16
     if(present(reserve)) reserve_ = reserve

     ! end definition mode, reserving some space for future additions
     ! this call also commits the changes to the disk
     __NF_ASRT__(nf__enddef(ncid,reserve_,4,0,4))
     ! arguments are ncid,h_minfree,v_align,v_minfree,r_align; default is (ncid,0,4,0,4).
     ! The above call reserves some space at the end of the netcdf header for
     ! future expansion without library's having to rewrite the entire file. See
     ! manual pages netcdf(3f) or netcdf(3) for more information.
  endif

  call mpp_sync()
end subroutine create_tile_out_file_idx_old

subroutine create_tile_out_file_idx_new(rhandle,name,tidx,tile_dim_length,zaxis_data,soilCCohort_data)
  type(restart_file_type), intent(inout) :: rhandle     ! restart file handle
  character(len=*),      intent(in)  :: name                ! name of the file to create
  integer              , intent(in)  :: tidx(:)             ! integer compressed index of tiles (local)
  integer              , intent(in)  :: tile_dim_length     ! length of tile axis
  real,        optional, intent(in)  :: zaxis_data(:)       ! data for the Z-axis
  real,        optional, intent(in)  :: soilCCohort_data(:)

  ! ---- local vars
  character(256) :: file_name ! full name of the file, including the processor number

  ! form the full name of the file
  call get_instance_filename(trim(name), file_name)
  call get_mosaic_tile_file(trim(file_name),file_name,lnd%domain)

  call fms_io_unstructured_register_restart_axis(rhandle, file_name, "lon", lnd%coord_glon, "X", &
          lnd%domain, units="degrees_east", longname="longitude")
  call fms_io_unstructured_register_restart_axis(rhandle, file_name, "lat", lnd%coord_glat, "Y", &
          lnd%domain, units="degrees_north", longname="latitude")
  ! the size of tile dimension really does not matter for the output, but it does
  ! matter for uncompressing utility, since it uses it as a size of the array to
  ! unpack to create tile index dimension and variable.
  call fms_io_unstructured_register_restart_axis(rhandle, file_name, trim(tile_index_name), &
          tidx, "tile lat lon", "C", tile_dim_length, lnd%domain, dimlen_name="tile", &
          dimlen_lname="tile number within grid cell", longname="compressed land point index", imin=0)
  if (present(zaxis_data)) then
      call fms_io_unstructured_register_restart_axis(rhandle, file_name, "zfull", zaxis_data, "Z", &
          lnd%domain, units="m", longname="full level", sense=-1)
  endif

  if (present(soilCCohort_data)) then
      call fms_io_unstructured_register_restart_axis(rhandle, file_name, "soilCCohort", soilCCohort_data, "CC", &
          lnd%domain, longname="Soil carbon cohort")
  endif
end subroutine create_tile_out_file_idx_new

! ============================================================================
! given a tile existence detection function, allocates and fills the tile index vector
subroutine gather_tile_index(tile_exists,idx)
  procedure(tile_exists_func) :: tile_exists   ! existence detector function:
         ! returns true if specific tile exists (hence should be written to restart)
  integer, allocatable,  intent(out) :: idx(:) ! rank local tile index vector

  ! ---- local vars
  type(land_tile_enum_type) :: ce ! tile list enumerator
  type(land_tile_type), pointer :: tile
  integer :: i,j,k,n

  ! count total number of tiles in this domain
  ce = first_elmt(land_tile_map)
  n  = 0
  do while (loop_over_tiles(ce,tile))
     if (tile_exists(tile)) n = n+1
  end do

  ! calculate compressed tile index to be written to the restart file;
  allocate(idx(max(n,1))); idx(:) = -1 ! set init value to a known invalid index
  ce = first_elmt(land_tile_map, lnd%ls)
  n = 1
  do while (loop_over_tiles(ce,tile,i=i,j=j,k=k))
     if(tile_exists(tile)) then
        idx(n) = (k-1)*lnd%nlon*lnd%nlat + (j-1)*lnd%nlon + (i-1)
        n = n+1
     endif
  end do
end subroutine gather_tile_index

! ============================================================================
! given compressed index, sizes of the global grid, 2D array of tile lists
! and the lower boundaries of this array returns a pointer to the tile
! corresponding to the compressed index, or NULL is the index is outside
! current domain, or such tile does not exist.
subroutine get_tile_by_idx(idx,nlon,nlat,tiles,ls,gs,ge,ptr)
   integer, intent(in) :: idx ! index
   integer, intent(in) :: nlon, nlat
   integer, intent(in) :: ls,gs,ge
   type(land_tile_list_type), intent(in) :: tiles(ls:)
   type(land_tile_type), pointer :: ptr

   ! ---- local vars
   integer :: g,k,npts,l
   type(land_tile_enum_type) :: ce

   ptr=>null()

   if(idx<0) return ! negative indices do not correspond to any tile

   ! given tile idx, calculate global lon, lat, and tile indices
   k = idx
   npts = nlon*nlat
   g = modulo(k,npts)+1 ; k = k/npts
   ! do nothing if the indices is outside of our domain
   if (g<gs.or.g>ge) return ! skip points outside of domain
   ! loop through the list of tiles at the given point to find k+1st tile
   l = lnd%l_index(g)
   if(l < lnd%ls .OR. l > lnd%le) then
      print*, " l= ", l, g, lnd%ls, lnd%le, mpp_pe()
      call error_mesg("land_tile_io", "l < lnd%ls .OR. l > lnd%le", FATAL)
   endif
   ce = first_elmt (tiles(l))
   do while(loop_over_tiles(ce, ptr).and.k>0)
      k = k-1
   enddo
   ! NOTE that at the end of the loop (that is, if there are less tiles in the list
   ! then requested by the idx), loop_over_tiles(ce,ptr) returns NULL

end subroutine get_tile_by_idx


! ============================================================================
! given the netcdf file id, name of the variable, and accessor subroutine, this
! subroutine reads integer 2D data (a scalar value per grid cell, that's why
! there is 0d in the name of this subroutine) and assigns the input values to
! each tile in respective grid cell.
subroutine read_tile_data_i0d_fptr(ncid,name,fptr)
   integer     , intent(in) :: ncid ! netcdf file id
   character(*), intent(in) :: name ! name of the variable to read
   procedure(fptr_i0)       :: fptr ! subroutine returning the pointer to the
                                    ! data to be written

   ! ---- local constants
   character(*), parameter :: module_name='read_tile_data_i0d_fptr'
   ! ---- local vars
   integer :: ndims     ! number of the variable dimensions
   integer :: dimids(1) ! IDs of the variable dimensions
   integer :: dimlen(1) ! size of the variable dimensions
   character(NF_MAX_NAME) :: idxname ! name of the index variable
   integer, allocatable :: idx(:) ! storage for compressed index
   integer, allocatable :: x1d(:) ! storage for the data
   integer :: i,j,bufsize
   integer :: varid, idxid
   type(land_tile_type), pointer :: tileptr ! pointer to tile
   integer, pointer :: ptr

   ! get the number of variable dimensions, and their lengths
   __NF_ASRT__(nfu_inq_var(ncid,name,id=varid,ndims=ndims,dimids=dimids,dimlens=dimlen))
   if(ndims/=1) then
      call error_mesg(module_name,'variable "'//trim(name)//'" has incorrect number of dimensions -- must be 1-dimensional', FATAL)
   endif
   ! get the name of compressed dimension and ID of corresponding variable
   __NF_ASRT__(nfu_inq_dim(ncid,dimids(1),idxname))
   __NF_ASRT__(nfu_inq_var(ncid,idxname,id=idxid))
   ! allocate input buffers for compression index and the variable
   bufsize = min(input_buf_size,dimlen(1))
   allocate(idx(bufsize),x1d(bufsize))
   ! read the input buffer-by-buffer
   do j = 1,dimlen(1),bufsize
      ! read the index variable
      __NF_ASRT__(nf_get_vara_int(ncid,idxid,(/j/),(/min(bufsize,dimlen(1)-j+1)/),idx))
      ! read the data
      __NF_ASRT__(nf_get_vara_int(ncid,varid,(/j/),(/min(bufsize,dimlen(1)-j+1)/),x1d))
      ! distribute the data over the tiles
      do i = 1, min(input_buf_size,dimlen(1)-j+1)
         call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                              lnd%ls,lnd%gs,lnd%ge, tileptr)
         call fptr(tileptr, ptr)
         if(associated(ptr)) ptr = x1d(i)
      enddo
   enddo
   ! release allocated memory
   deallocate(idx,x1d)
end subroutine read_tile_data_i0d_fptr

subroutine read_tile_data_r0d_fptr_r0(ncid,name,fptr)
   integer     , intent(in) :: ncid ! netcdf file id
   character(*), intent(in) :: name ! name of the variable to read
   procedure(fptr_r0)       :: fptr ! subroutine returning the pointer to the
                                    ! data to be written
   ! ---- local constants
   character(*), parameter :: module_name='read_tile_data_r0d_fptr_r0'
   ! ---- local vars
   integer :: ndims     ! number of the variable dimensions
   integer :: dimids(1) ! IDs of the variable dimensions
   integer :: dimlen(1) ! size of the variable dimensions
   character(NF_MAX_NAME) :: idxname ! name of the index variable
   integer, allocatable :: idx(:) ! storage for compressed index
   real   , allocatable :: x1d(:) ! storage for the data
   integer :: i, j, bufsize
   integer :: varid, idxid
   type(land_tile_type), pointer :: tileptr ! pointer to tile
   real, pointer :: ptr

   ! get the number of variable dimensions, and their lengths
   __NF_ASRT__(nfu_inq_var(ncid,name,id=varid,ndims=ndims,dimids=dimids,dimlens=dimlen))
   if(ndims/=1) then
      call error_mesg(module_name,'variable "'//trim(name)//'" has incorrect number of dimensions -- must be 1-dimensional', FATAL)
   endif
   ! get the name of compressed dimension and ID of corresponding variable
   __NF_ASRT__(nfu_inq_dim(ncid,dimids(1),idxname))
   __NF_ASRT__(nfu_inq_var(ncid,idxname,id=idxid))
   ! allocate input buffers for compression index and the variable
   bufsize=min(input_buf_size,dimlen(1))
   allocate(idx(bufsize),x1d(bufsize))
   ! read the input buffer-by-buffer
   do j = 1,dimlen(1),bufsize
      ! read the index variable
      __NF_ASRT__(nf_get_vara_int(ncid,idxid,(/j/),(/min(bufsize,dimlen(1)-j+1)/),idx))
      ! read the data
      __NF_ASRT__(nf_get_vara_double(ncid,varid,(/j/),(/min(bufsize,dimlen(1)-j+1)/),x1d))
      ! distribute the data over the tiles
      do i = 1, min(bufsize,dimlen(1)-j+1)
         call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                              lnd%ls,lnd%gs,lnd%ge, tileptr)
         call fptr(tileptr, ptr)
         if(associated(ptr)) ptr = x1d(i)
      enddo
   enddo
   ! release allocated memory
   deallocate(idx,x1d)
end subroutine read_tile_data_r0d_fptr_r0

subroutine read_tile_data_r0d_fptr_r0i (ncid,name,fptr,index)
   integer     , intent(in) :: ncid ! netcdf file id
   character(*), intent(in) :: name ! name of the variable to read
   procedure(fptr_r0i)      :: fptr ! subroutine returning the pointer to the
                                    ! data to be written
   integer     , intent(in) :: index ! index where to read the data

   ! ---- local constants
   character(*), parameter :: module_name='read_tile_data_r0d_fptr_r0i'
   ! ---- local vars
   integer :: ndims     ! number of the variable dimensions
   integer :: dimids(1) ! IDs of the variable dimensions
   integer :: dimlen(1) ! size of the variable dimensions
   character(NF_MAX_NAME) :: idxname ! name of the index variable
   integer, allocatable :: idx(:) ! storage for compressed index
   real   , allocatable :: x1d(:) ! storage for the data
   integer :: i,j,bufsize
   integer :: varid, idxid
   type(land_tile_type), pointer :: tileptr ! pointer to tile
   real, pointer :: ptr

   ! get the number of variable dimensions, and their lengths
   __NF_ASRT__(nfu_inq_var(ncid,name,id=varid,ndims=ndims,dimids=dimids,dimlens=dimlen))
   if(ndims/=1) then
      call error_mesg(module_name,'variable "'//trim(name)//'" has incorrect number of dimensions -- must be 1-dimensional', FATAL)
   endif
   ! get the name of compressed dimension and ID of corresponding variable
   __NF_ASRT__(nfu_inq_dim(ncid,dimids(1),idxname))
   __NF_ASRT__(nfu_inq_var(ncid,idxname,id=idxid))
   ! allocate input buffers for compression index and the variable
   bufsize=min(input_buf_size,dimlen(1))
   allocate(idx(bufsize),x1d(bufsize))
   ! read the input buffer-by-buffer
   do j = 1,dimlen(1),bufsize
      ! read the index variable
      __NF_ASRT__(nf_get_vara_int(ncid,idxid,(/j/),(/min(bufsize,dimlen(1)-j+1)/),idx))
      ! read the data
      __NF_ASRT__(nf_get_vara_double(ncid,varid,(/j/),(/min(bufsize,dimlen(1)-j+1)/),x1d))
      ! distribute the data over the tiles
      do i = 1, min(bufsize,dimlen(1)-j+1)
         call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                              lnd%ls,lnd%gs,lnd%ge, tileptr)
         call fptr(tileptr, index, ptr)
         if(associated(ptr)) ptr = x1d(i)
      enddo
   enddo
   ! release allocated memory
   deallocate(idx,x1d)
end subroutine read_tile_data_r0d_fptr_r0i

subroutine read_tile_data_i1d_fptr_i0i(ncid,name,fptr)
   integer     , intent(in) :: ncid ! netcdf file id
   character(*), intent(in) :: name ! name of the variable to read
   procedure(fptr_i0i)      :: fptr ! subroutine returning the pointer to the data

   ! ---- local constants
   character(*), parameter :: module_name='read_tile_data_i1d_fptr_i0i'
   ! ---- local vars
   integer :: ndims     ! number of the variable dimensions
   integer :: dimids(2) ! IDs of the variable dimensions
   integer :: dimlen(2) ! size of the variable dimensions
   character(NF_MAX_NAME) :: idxname ! name of the index variable
   integer, allocatable :: idx(:)   ! storage for compressed index
   integer, allocatable :: x1d(:)   ! storage for the data
   integer :: i, j, n, bufsize
   integer :: varid,idxid
   integer :: start(2), count(2) ! input slab parameters
   type(land_tile_type), pointer :: tileptr ! pointer to tile
   integer, pointer :: ptr

   ! get the number of variable dimensions, and their lengths
   __NF_ASRT__(nfu_inq_var(ncid,name,id=varid,ndims=ndims,dimids=dimids,dimlens=dimlen))
   if(ndims/=2) then
      call error_mesg(module_name,'variable "'//trim(name)//'" has incorrect number of dimensions -- must be 2-dimensional', FATAL)
   endif
   ! get the name of compressed dimension and ID of corresponding variable
   __NF_ASRT__(nfu_inq_dim(ncid,dimids(1),idxname))
   __NF_ASRT__(nfu_inq_var(ncid,idxname,id=idxid))
   ! allocate input buffers for compression index and the variable
   bufsize=min(input_buf_size,dimlen(1))
   allocate(idx(bufsize),x1d(bufsize*dimlen(2)))
   ! read the input buffer-by-buffer
   do j = 1,dimlen(1),bufsize
      ! set up slab parameters
      start(1) = j ; count(1) = min(bufsize,dimlen(1)-j+1)
      start(2) = 1 ; count(2) = dimlen(2)
      ! read the index variable
      __NF_ASRT__(nf_get_vara_int(ncid,idxid,start(1),count(1),idx))
      ! read the data
      __NF_ASRT__(nf_get_vara_int(ncid,varid,start,count,x1d))
      ! distribute the data over the tiles
      do i = 1, min(bufsize,dimlen(1)-j+1)
         call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                              lnd%ls,lnd%gs, lnd%ge, tileptr)
         do n = 1,count(2)
            call fptr(tileptr, n, ptr)
            if(associated(ptr)) ptr = x1d(i+count(1)*(n-1))
         enddo
      enddo
   enddo
   ! release allocated memory
   deallocate(idx,x1d)
end subroutine read_tile_data_i1d_fptr_i0i

subroutine read_tile_data_r1d_fptr_r0i(ncid,name,fptr)
   integer     , intent(in) :: ncid ! netcdf file id
   character(*), intent(in) :: name ! name of the variable to read
   procedure(fptr_r0i)      :: fptr ! subroutine returning the pointer to the data

   ! ---- local constants
   character(*), parameter :: module_name='read_tile_data_r1d_fptr_r0i'
   ! ---- local vars
   integer :: ndims     ! number of the variable dimensions
   integer :: dimids(2) ! IDs of the variable dimensions
   integer :: dimlen(2) ! size of the variable dimensions
   character(NF_MAX_NAME) :: idxname ! name of the index variable
   integer, allocatable :: idx(:)   ! storage for compressed index
   real   , allocatable :: x1d(:)   ! storage for the data
   integer :: i, j, n, bufsize
   integer :: varid,idxid
   integer :: start(2), count(2) ! input slab parameters
   type(land_tile_type), pointer :: tileptr ! pointer to tile
   real, pointer :: ptr

   ! get the number of variable dimensions, and their lengths
   __NF_ASRT__(nfu_inq_var(ncid,name,id=varid,ndims=ndims,dimids=dimids,dimlens=dimlen))
   if(ndims/=2) then
      call error_mesg(module_name,'variable "'//trim(name)//'" has incorrect number of dimensions -- must be 2-dimensional', FATAL)
   endif
   ! get the name of compressed dimension and ID of corresponding variable
   __NF_ASRT__(nfu_inq_dim(ncid,dimids(1),idxname))
   __NF_ASRT__(nfu_inq_var(ncid,idxname,id=idxid))
   ! allocate input buffers for compression index and the variable
   bufsize=min(input_buf_size,dimlen(1))
   allocate(idx(bufsize),x1d(bufsize*dimlen(2)))
   ! read the input buffer-by-buffer
   do j = 1,dimlen(1),bufsize
      ! set up slab parameters
      start(1) = j ; count(1) = min(bufsize,dimlen(1)-j+1)
      start(2) = 1 ; count(2) = dimlen(2)
      ! read the index variable
      __NF_ASRT__(nf_get_vara_int(ncid,idxid,start(1),count(1),idx))
      ! read the data
      __NF_ASRT__(nf_get_vara_double(ncid,varid,start,count,x1d))
      ! distribute the data over the tiles
      do i = 1, min(bufsize,dimlen(1)-j+1)
         call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                              lnd%ls,lnd%gs,lnd%ge, tileptr)
         do n = 1,count(2)
            call fptr(tileptr, n, ptr)
            if(associated(ptr)) ptr = x1d(i+count(1)*(n-1))
         enddo
      enddo
   enddo
   ! release allocated memory
   deallocate(idx,x1d)
end subroutine read_tile_data_r1d_fptr_r0i

subroutine read_tile_data_r1d_fptr_r0ij(ncid,name,fptr,index)
   integer     , intent(in) :: ncid ! netcdf file id
   character(*), intent(in) :: name ! name of the variable to read
   procedure(fptr_r0ij)     :: fptr ! subroutine returning the pointer to the data
   integer     , intent(in) :: index

   ! ---- local constants
   character(*), parameter :: module_name='read_tile_data_r1d_fptr_r0ij'
   ! ---- local vars
   integer :: ndims     ! number of the variable dimensions
   integer :: dimids(2) ! IDs of the variable dimensions
   integer :: dimlen(2) ! size of the variable dimensions
   character(NF_MAX_NAME) :: idxname ! name of the index variable
   integer, allocatable :: idx(:)   ! storage for compressed index
   real   , allocatable :: x1d(:)   ! storage for the data
   integer :: i, j, n, bufsize
   integer :: varid,idxid
   integer :: start(2), count(2) ! input slab parameters
   type(land_tile_type), pointer :: tileptr ! pointer to tile
   real, pointer :: ptr

   ! get the number of variable dimensions, and their lengths
   __NF_ASRT__(nfu_inq_var(ncid,name,id=varid,ndims=ndims,dimids=dimids,dimlens=dimlen))
   if(ndims/=2) then
      call error_mesg(module_name,'variable "'//trim(name)//'" has incorrect number of dimensions -- must be 2-dimensional', FATAL)
   endif
   ! get the name of compressed dimension and ID of corresponding variable
   __NF_ASRT__(nfu_inq_dim(ncid,dimids(1),idxname))
   __NF_ASRT__(nfu_inq_var(ncid,idxname,id=idxid))
   ! allocate input buffers for compression index and the variable
   bufsize=min(input_buf_size,dimlen(1))
   allocate(idx(bufsize),x1d(bufsize*dimlen(2)))
   ! read the input buffer-by-buffer
   do j = 1,dimlen(1),bufsize
      ! set up slab parameters
      start(1) = j ; count(1) = min(bufsize,dimlen(1)-j+1)
      start(2) = 1 ; count(2) = dimlen(2)
      ! read the index variable
      __NF_ASRT__(nf_get_vara_int(ncid,idxid,start(1),count(1),idx))
      ! read the data
      __NF_ASRT__(nf_get_vara_double(ncid,varid,start,count,x1d))
      ! distribute the data over the tiles
      do i = 1, min(bufsize,dimlen(1)-j+1)
         call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                              lnd%ls,lnd%gs,lnd%ge, tileptr)
         do n = 1,count(2)
            call fptr(tileptr, n, index, ptr)
            if(associated(ptr)) ptr = x1d(i+count(1)*(n-1))
         enddo
      enddo
   enddo
   ! release allocated memory
   deallocate(idx,x1d)
end subroutine read_tile_data_r1d_fptr_r0ij

subroutine read_tile_data_r2d_fptr_r0ij (ncid,name,fptr)
   integer     , intent(in) :: ncid ! netcdf file id
   character(*), intent(in) :: name ! name of the variable to read
   procedure(fptr_r0ij)     :: fptr ! subroutine returning the pointer to the data

   ! ---- local constants
   character(*), parameter :: module_name='read_tile_data_r2d_fptr_r0ij'
   ! ---- local vars
   integer :: ndims     ! number of the variable dimensions
   integer :: dimids(3) ! IDs of the variable dimensions
   integer :: dimlen(3) ! size of the variable dimensions
   character(NF_MAX_NAME) :: idxname ! name of the index variable
   integer, allocatable :: idx(:) ! storage for compressed index
   real   , allocatable :: x1d(:) ! storage for the data
   integer :: i,j,m,n,bufsize
   integer :: varid, idxid
   integer :: start(3), count(3) ! input slab parameters
   type(land_tile_type), pointer :: tileptr ! pointer to tile
   real, pointer :: ptr

   ! get the number of variable dimensions, and their lengths
   __NF_ASRT__(nfu_inq_var(ncid,name,id=varid,ndims=ndims,dimids=dimids,dimlens=dimlen))
   if(ndims/=3) then
      call error_mesg(module_name,'variable "'//trim(name)//'" has incorrect number of dimensions -- must be 3-dimensional', FATAL)
   endif
   ! get the name of compressed dimension and ID of corresponding variable
   __NF_ASRT__(nfu_inq_dim(ncid,dimids(1),idxname))
   __NF_ASRT__(nfu_inq_var(ncid,idxname,id=idxid))
   ! allocate input buffers for compression index and the variable
   bufsize=min(input_buf_size,dimlen(1))
   allocate(idx(bufsize),x1d(bufsize*dimlen(2)*dimlen(3)))
   ! read the input buffer-by-buffer
   do j = 1,dimlen(1),bufsize
      ! set up slab parameters
      start(1) = j ; count(1) = min(bufsize,dimlen(1)-j+1)
      start(2) = 1 ; count(2) = dimlen(2)
      start(3) = 1 ; count(3) = dimlen(3)
      ! read the index variable
      __NF_ASRT__(nf_get_vara_int(ncid,idxid,start(1),count(1),idx))
      ! read the data
      __NF_ASRT__(nf_get_vara_double(ncid,varid,start,count,x1d))
      ! distribute the data over the tiles
      do i = 1, min(bufsize,dimlen(1)-j+1)
         call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                              lnd%ls,lnd%gs,lnd%ge, tileptr)
         do m = 1,dimlen(2)
         do n = 1,dimlen(3)
            call fptr(tileptr, m, n, ptr)
            if(associated(ptr)) ptr = x1d(i+count(1)*(m-1)+count(1)*count(2)*(n-1))
         enddo
         enddo
      enddo
   enddo
   ! release allocated memory
   deallocate(idx,x1d)
end subroutine read_tile_data_r2d_fptr_r0ij

subroutine read_tile_data_r2d_fptr_r0ijk (ncid,name,fptr,index)
   integer     , intent(in) :: ncid ! netcdf file id
   character(*), intent(in) :: name ! name of the variable to read
   procedure(fptr_r0ijk)    :: fptr ! subroutine returning the pointer to the data
   integer, intent(in)      :: index

   ! ---- local constants
   character(*), parameter :: module_name='read_tile_data_r2d_fptr_r0ijk'
   ! ---- local vars
   integer :: ndims     ! number of the variable dimensions
   integer :: dimids(3) ! IDs of the variable dimensions
   integer :: dimlen(3) ! size of the variable dimensions
   character(NF_MAX_NAME) :: idxname ! name of the index variable
   integer, allocatable :: idx(:) ! storage for compressed index
   real   , allocatable :: x1d(:) ! storage for the data
   integer :: i,j,m,n,bufsize
   integer :: varid, idxid
   integer :: start(3), count(3) ! input slab parameters
   type(land_tile_type), pointer :: tileptr ! pointer to tile
   real, pointer :: ptr

   ! get the number of variable dimensions, and their lengths
   __NF_ASRT__(nfu_inq_var(ncid,name,id=varid,ndims=ndims,dimids=dimids,dimlens=dimlen))
   if(ndims/=3) then
      call error_mesg(module_name,'variable "'//trim(name)//'" has incorrect number of dimensions -- must be 3-dimensional', FATAL)
   endif
   ! get the name of compressed dimension and ID of corresponding variable
   __NF_ASRT__(nfu_inq_dim(ncid,dimids(1),idxname))
   __NF_ASRT__(nfu_inq_var(ncid,idxname,id=idxid))
   ! allocate input buffers for compression index and the variable
   bufsize=min(input_buf_size,dimlen(1))
   allocate(idx(bufsize),x1d(bufsize*dimlen(2)*dimlen(3)))
   ! read the input buffer-by-buffer
   do j = 1,dimlen(1),bufsize
      ! set up slab parameters
      start(1) = j ; count(1) = min(bufsize,dimlen(1)-j+1)
      start(2) = 1 ; count(2) = dimlen(2)
      start(3) = 1 ; count(3) = dimlen(3)
      ! read the index variable
      __NF_ASRT__(nf_get_vara_int(ncid,idxid,start(1),count(1),idx))
      ! read the data
      __NF_ASRT__(nf_get_vara_double(ncid,varid,start,count,x1d))
      ! distribute the data over the tiles
      do i = 1, min(bufsize,dimlen(1)-j+1)
         call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                              lnd%ls,lnd%gs,lnd%ge, tileptr)
         do m = 1,dimlen(2)
         do n = 1,dimlen(3)
            call fptr(tileptr, m, n, index, ptr)
            if(associated(ptr)) ptr = x1d(i+count(1)*(m-1)+count(1)*count(2)*(n-1))
         enddo
         enddo
      enddo
   enddo
   ! release allocated memory
   deallocate(idx,x1d)
end subroutine read_tile_data_r2d_fptr_r0ijk


! ============================================================================
! The subroutines write_tile_data_* below write tiled data (that is, data provided
! in arrays congruous with current tiling) to NetCDF files using "compression
! by gathering" (see CF conventions). They assume that the compressed dimension
! is already created, has certain name (see parameter "tile_index_name" at the beginning
! of this file), is written in the same order.

! ============================================================================
! writes out 1-d integer tiled data using "compression by gathering"
subroutine write_tile_data_i1d(ncid,name,data,long_name,units)
  integer         , intent(in) :: ncid
  character(len=*), intent(in) :: name
  integer         , intent(in) :: data(:) ! data to write
  character(len=*), intent(in), optional :: units, long_name

  ! local vars
  integer :: varid,iret,p,k
  integer, allocatable :: buffer(:) ! data buffers
  integer, allocatable :: ntiles(:) ! list of land tile numbers for each of PEs in io_domain

  if (mpp_pe()/=lnd%io_pelist(1)) then
     call mpp_send(size(data), plen=1,       to_pe=lnd%io_pelist(1), tag=COMM_TAG_3)
     call mpp_send(data(1), plen=size(data), to_pe=lnd%io_pelist(1), tag=COMM_TAG_4)
  else
     allocate(ntiles(size(lnd%io_pelist)))
     ntiles(1) = size(data)
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(ntiles(p), from_pe=lnd%io_pelist(p), glen=1, tag=COMM_TAG_3)
     enddo
     ! gather data from all processors in io_domain
     allocate(buffer(sum(ntiles(:))))
     buffer(1:ntiles(1)) = data(:)
     k=ntiles(1)+1
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(buffer(k), glen=ntiles(p), from_pe=lnd%io_pelist(p), tag=COMM_TAG_4)
        k = k+ntiles(p)
     enddo
     ! create variable, if it does not exist
     if(nf_inq_varid(ncid,name,varid)/=NF_NOERR) then
        __NF_ASRT__(nfu_def_var(ncid,name,NF_INT,(/tile_index_name/),long_name,units,varid))
     endif
     ! write data
     iret = nf_enddef(ncid) ! ignore errors (file may be in data mode already)
     __NF_ASRT__(nf_put_var_int(ncid,varid,buffer))
  endif
  ! wait for all PEs to finish: necessary because mpp_send doesn't seem to
  ! copy the data, and therefore on non-root io_domain PE there would be a chance
  ! that the data and mask are destroyed before they are actually sent.
  call mpp_sync()
end subroutine write_tile_data_i1d

! ============================================================================
! writes out 1-d real tiled data using "compression by gathering"
subroutine write_tile_data_r1d(ncid,name,data,long_name,units)
  integer         , intent(in) :: ncid    ! netcdf ID
  character(len=*), intent(in) :: name    ! name of the variable
  real            , intent(in) :: data(:) ! data to write
  character(len=*), intent(in), optional :: units, long_name ! attributes

  ! ---- local vars
  integer :: varid,iret,p,k
  real,    allocatable :: buffer(:) ! data buffer
  integer, allocatable :: ntiles(:) ! list of land tile numbers for each of PEs in io_domain

  ! if our PE doesn't do IO (that is, it isn't the root io_domain processor),
  ! simply send data size and data the root IO processor
  if (mpp_pe()/=lnd%io_pelist(1)) then
     call mpp_send(size(data), plen=1,          to_pe=lnd%io_pelist(1), tag=COMM_TAG_5)
     call mpp_send(data(1),    plen=size(data), to_pe=lnd%io_pelist(1), tag=COMM_TAG_6)
  else
     allocate(ntiles(size(lnd%io_pelist)))
     ntiles(1) = size(data)
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(ntiles(p), from_pe=lnd%io_pelist(p), glen=1, tag=COMM_TAG_5)
     enddo
     ! gather data from the processors in io_domain
     allocate(buffer(sum(ntiles(:))))
     buffer(1:ntiles(1)) = data(:)
     k=ntiles(1)+1
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(buffer(k), glen=ntiles(p), from_pe=lnd%io_pelist(p), tag=COMM_TAG_6)
        k = k+ntiles(p)
     enddo
     ! create variable, if it does not exist
     if(nf_inq_varid(ncid,name,varid)/=NF_NOERR) then
        __NF_ASRT__(nfu_def_var(ncid,name,NF_DOUBLE,(/tile_index_name/),long_name,units,varid))
     endif
     ! write data
     iret = nf_enddef(ncid) ! ignore errors (file may be in data mode already)
     __NF_ASRT__(nf_put_var_double(ncid,varid,buffer))
     deallocate(buffer,ntiles)
  endif
  ! wait for all PEs to finish: necessary because mpp_send doesn't seem to
  ! copy the data, and therefore on non-root io_domain PE there would be a chance
  ! that the data and mask are destroyed before they are actually sent.
  call mpp_sync()
end subroutine write_tile_data_r1d

! ============================================================================
! writes out 2-d integer tiled data using "compression by gathering". The dimension
! of the data is (tile,z), and both tile and z dimensions are assumed to be
! already created
subroutine write_tile_data_i2d(ncid,name,data,zdim,long_name,units)
  integer         , intent(in) :: ncid ! netcdf id
  character(len=*), intent(in) :: name ! name of the variable to write
  character(len=*), intent(in) :: zdim ! name of the z-dimension
  integer         , intent(inout) :: data(:,:) ! (tile,z)
  character(len=*), intent(in), optional :: units, long_name
  ! data and mask are "inout" to save the memory on send-receive buffers. On the
  ! root io_domain PE mask is destroyed and data is filled with the information
  ! from other PEs in our io_domain. On other PEs these arrays reman intact.

  ! local vars
  integer :: varid,iret,p,i,k
  character(NF_MAX_NAME)::dimnames(2)
  integer, allocatable :: buff2(:,:), buff1(:) ! send/receive buffers
  integer, allocatable :: ntiles(:)   ! list of land tile numbers for each of PEs in io_domain

  ! if our PE doesn't do io (that is, it isn't the root io_domain processor),
  ! simply send the data and mask of valid data to the root IO processor
  if (mpp_pe()/=lnd%io_pelist(1)) then
     call mpp_send(size(data,1), plen=1,          to_pe=lnd%io_pelist(1), tag=COMM_TAG_7)
     call mpp_send(data(1,1),    plen=size(data), to_pe=lnd%io_pelist(1), tag=COMM_TAG_8)
  else
     allocate(ntiles(size(lnd%io_pelist)))
     ntiles(1) = size(data,1)
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(ntiles(p), from_pe=lnd%io_pelist(p), glen=1, tag=COMM_TAG_7)
     enddo
     allocate(buff2(sum(ntiles),size(data,2)),buff1(maxval(ntiles)*size(data,2)))
     ! gather data from the processors in our io_domain
     buff2(1:ntiles(1),:) = data(:,:)
     k=ntiles(1)
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(buff1(1), glen=ntiles(p)*size(data,2), from_pe=lnd%io_pelist(p), tag=COMM_TAG_8)
        do i = 1,size(data,2)
           buff2(k+1:k+ntiles(p),i) = buff1((i-1)*ntiles(p)+1:i*ntiles(p))
        enddo
     enddo

     ! create variable, if it does not exist
     if(nf_inq_varid(ncid,name,varid)/=NF_NOERR) then
        dimnames(1) = tile_index_name
        dimnames(2) = zdim
        __NF_ASRT__(nfu_def_var(ncid,name,NF_DOUBLE,dimnames,long_name,units,varid))
     endif
     ! write data
     iret = nf_enddef(ncid) ! ignore errors: its OK if file is in data mode already
     __NF_ASRT__(nf_put_var_double(ncid,varid,buff2))
     deallocate(buff2,buff1,ntiles)
  endif
  ! wait for all PEs to finish: necessary because mpp_send doesn't seem to
  ! copy the data, and therefore on non-root io_domain PE there would be a chance
  ! that the data and mask are destroyed before they are actually sent.
  call mpp_sync()
end subroutine write_tile_data_i2d

! ============================================================================
! writes out 2-d real tiled data using "compression by gathering". The dimension
! of the data is (tile,z), and both tile and z dimensions are assumed to be
! already created
subroutine write_tile_data_r2d(ncid,name,data,zdim,long_name,units)
  integer         , intent(in) :: ncid ! netcdf id
  character(len=*), intent(in) :: name ! name of the variable to write
  character(len=*), intent(in) :: zdim ! name of the z-dimension
  real            , intent(in) :: data(:,:) ! (tile,z)
  character(len=*), intent(in), optional :: units, long_name

  ! local vars
  integer :: varid,iret,p,i,k
  character(NF_MAX_NAME)::dimnames(2)
  real,    allocatable :: buff2(:,:), buff1(:) ! send/receive buffers
  integer, allocatable :: ntiles(:)   ! list of land tile numbers for each of PEs in io_domain

  ! if our PE doesn't do io (that is, it isn't the root io_domain processor),
  ! simply send the data and mask of valid data to the root IO processor
  if (mpp_pe()/=lnd%io_pelist(1)) then
     call mpp_send(size(data,1), plen=1,          to_pe=lnd%io_pelist(1), tag=COMM_TAG_7)
     call mpp_send(data(1,1),    plen=size(data), to_pe=lnd%io_pelist(1), tag=COMM_TAG_8)
  else
     allocate(ntiles(size(lnd%io_pelist)))
     ntiles(1) = size(data,1)
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(ntiles(p), from_pe=lnd%io_pelist(p), glen=1, tag=COMM_TAG_7)
     enddo
     allocate(buff2(sum(ntiles),size(data,2)),buff1(maxval(ntiles)*size(data,2)))
     ! gather data from the processors in our io_domain
     buff2(1:ntiles(1),:) = data(:,:)
     k=ntiles(1)
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(buff1(1), glen=ntiles(p)*size(data,2), from_pe=lnd%io_pelist(p), tag=COMM_TAG_8)
        do i = 1,size(data,2)
           buff2(k+1:k+ntiles(p),i) = buff1((i-1)*ntiles(p)+1:i*ntiles(p))
        enddo
        k = k+ntiles(p)
     enddo

     ! create variable, if it does not exist
     if(nf_inq_varid(ncid,name,varid)/=NF_NOERR) then
        dimnames(1) = tile_index_name
        dimnames(2) = zdim
        __NF_ASRT__(nfu_def_var(ncid,name,NF_DOUBLE,dimnames,long_name,units,varid))
     endif
     ! write data
     iret = nf_enddef(ncid) ! ignore errors: its OK if file is in data mode already
     __NF_ASRT__(nf_put_var_double(ncid,varid,buff2))
     deallocate(buff2,buff1,ntiles)
  endif
  ! wait for all PEs to finish: necessary because mpp_send doesn't seem to
  ! copy the data, and therefore on non-root io_domain PE there would be a chance
  ! that the data and mask are destroyed before they are actually sent.
  call mpp_sync()
end subroutine write_tile_data_r2d

! ============================================================================
! writes out 3-d real tiled data using "compression by gathering". The dimension
! of the data is (tile,z,cohort), and both tile and z dimensions are assumed to be
! already created
subroutine write_tile_data_r3d(ncid,name,data,dim1,dim2,long_name,units)
  integer         , intent(in) :: ncid ! netcdf id
  character(len=*), intent(in) :: name ! name of the variable to write
  character(len=*), intent(in) :: dim1,dim2 ! names of dimensions
  real            , intent(in) :: data(:,:,:) ! (tile,dim1,dim2)
  character(len=*), intent(in), optional :: units, long_name

  ! local vars
  integer :: varid,iret,p,i,j,k,n
  character(NF_MAX_NAME)::dimnames(3)
  real, allocatable :: buff3(:,:,:),buff1(:) ! send/receive buffers
  integer, allocatable :: ntiles(:)   ! list of land tile numbers for each of PEs in io_domain

  ! if our PE does not do io (that is, it is not the root io_domain processor),
  ! simply send the data and mask of valid data to the root IO processor
  if (mpp_pe()/=lnd%io_pelist(1)) then
     call mpp_send(size(data,1), plen=1,          to_pe=lnd%io_pelist(1), tag=COMM_TAG_9)
     call mpp_send(data(1,1,1),  plen=size(data), to_pe=lnd%io_pelist(1), tag=COMM_TAG_10)
  else
     allocate(ntiles(size(lnd%io_pelist)))
     ntiles(1) = size(data,1)
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(ntiles(p), from_pe=lnd%io_pelist(p), glen=1, tag=COMM_TAG_9)
     enddo
     allocate(buff3(sum(ntiles),size(data,2),size(data,3)),&
              buff1(maxval(ntiles)*size(data,2)*size(data,3)))
     ! gather data from the processors in our io_domain
     buff3(1:ntiles(1),:,:) = data(:,:,:)
     k=ntiles(1)
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(buff1(1), glen=ntiles(p)*size(data,2)*size(data,3), from_pe=lnd%io_pelist(p),&
                      tag=COMM_TAG_10)
        n = 0
        do i = 1,size(data,2)
        do j = 1,size(data,3)
           buff3(k+1:k+ntiles(p),i,j) = buff1(n*ntiles(p)+1:n*ntiles(p))
           n = n+ntiles(p)
        enddo
        enddo
        k = k+ntiles(p)
     enddo
     ! create variable, if it does not exist
     if(nf_inq_varid(ncid,name,varid)/=NF_NOERR) then
        dimnames(1) = tile_index_name
        dimnames(2) = dim1
        dimnames(3) = dim2
        __NF_ASRT__(nfu_def_var(ncid,name,NF_DOUBLE,dimnames,long_name,units,varid))
     endif
     ! write data
     iret = nf_enddef(ncid) ! ignore errors: its OK if file is in data mode already
     __NF_ASRT__(nf_put_var_double(ncid,varid,buff3))
  endif
  ! wait for all PEs to finish: necessary because mpp_send does not seem to
  ! copy the data, and therefore on non-root io_domain PE there would be a chance
  ! that the data and mask are destroyed before they are actually sent.
  call mpp_sync()
end subroutine write_tile_data_r3d

! ============================================================================
subroutine gather_tile_data_i0d(fptr,idx,data)
  procedure(fptr_i0)  :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  integer, intent(out) :: data(:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  integer, pointer :: ptr ! pointer to the tile data
  integer :: i

  data = NF_FILL_INT

! gather data into an array along the tile dimension. It is assumed that
! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                          lnd%ls,lnd%gs,lnd%ge, tileptr)
     call fptr(tileptr, ptr)
     if(associated(ptr)) data(i)=ptr
  enddo
end subroutine gather_tile_data_i0d

subroutine gather_tile_data_r0d(fptr,idx,data)
  procedure(fptr_r0) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real, intent(out) :: data(:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i

  data = NF_FILL_DOUBLE

! gather data into an array along the tile dimension. It is assumed that
! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                          lnd%ls,lnd%gs,lnd%ge, tileptr)
     call fptr(tileptr, ptr)
     if(associated(ptr)) data(i)=ptr
  enddo
end subroutine gather_tile_data_r0d

subroutine gather_tile_data_r0d_idx(fptr,n,idx,data)
  procedure(fptr_r0i) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: n   ! additional index argument for fptr
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real, intent(out) :: data(:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i

  data = NF_FILL_DOUBLE

! gather data into an array along the tile dimension. It is assumed that
! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                          lnd%ls,lnd%gs,lnd%ge, tileptr)
     call fptr(tileptr, n, ptr)
     if(associated(ptr)) data(i)=ptr
  enddo
end subroutine gather_tile_data_r0d_idx

subroutine gather_tile_data_r1d(fptr,idx,data)
  procedure(fptr_r0i) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real, intent(out) :: data(:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i,j

  data = NF_FILL_DOUBLE

! gather data into an array along the tile dimension. It is assumed that
! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                          lnd%ls,lnd%gs,lnd%ge, tileptr)
     do j = 1,size(data,2)
        call fptr(tileptr, j, ptr)
        if(associated(ptr)) data(i,j)=ptr
     enddo
  enddo
end subroutine gather_tile_data_r1d

subroutine gather_tile_data_i1d(fptr,idx,data)
  procedure(fptr_i0i) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  integer, intent(out) :: data(:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  integer, pointer :: ptr ! pointer to the tile data
  integer :: i,j

  data = NF_FILL_INT

! gather data into an array along the tile dimension. It is assumed that
! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                          lnd%ls,lnd%gs,lnd%ge, tileptr)
     do j = 1,size(data,2)
        call fptr(tileptr, j, ptr)
        if(associated(ptr)) data(i,j)=ptr
     enddo
  enddo
end subroutine gather_tile_data_i1d

subroutine gather_tile_data_r2d(fptr,idx,data)
  procedure(fptr_r0ij):: fptr ! subroutine returning the pointer to the data to be written
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real, intent(out) :: data(:,:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i,k,m

  data = NF_FILL_DOUBLE

! gather data into an array along the tile dimension. It is assumed that
! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                          lnd%ls,lnd%gs,lnd%ge, tileptr)
     do k=1,size(data,2)
     do m=1,size(data,3)
        call fptr(tileptr, k, m, ptr)
        if(associated(ptr)) data(i,k,m)=ptr
     enddo
     enddo
  enddo
end subroutine gather_tile_data_r2d

subroutine gather_tile_data_r2d_idx(fptr,n,idx,data)
  procedure(fptr_r0ijk) :: fptr ! subroutine returning the pointer to the data
  integer, intent(in) :: n ! additional index argument for fptr
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real, intent(out) :: data(:,:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i,k,m

  data = NF_FILL_DOUBLE

! gather data into an array along the tile dimension. It is assumed that
! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                          lnd%ls,lnd%gs,lnd%ge, tileptr)
     do k=1,size(data,2)
     do m=1,size(data,3)
        call fptr(tileptr, k, m, n, ptr)
        if(associated(ptr)) data(i,k,m)=ptr
     enddo
     enddo
  enddo
end subroutine gather_tile_data_r2d_idx

subroutine assemble_tiles_i0d(fptr,idx,data)
  procedure(fptr_i0) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  integer, intent(in) :: data(:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  integer, pointer :: ptr ! pointer to the tile data
  integer :: i

! distribute the data over the tiles
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                          lnd%ls,lnd%gs,lnd%ge, tileptr)
     call fptr(tileptr, ptr)
     if(associated(ptr)) ptr=data(i)
  enddo
end subroutine assemble_tiles_i0d

subroutine assemble_tiles_r0d(fptr,idx,data)
  procedure(fptr_r0) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real,    intent(in) :: data(:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i

! distribute the data over the tiles
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                          lnd%ls,lnd%gs,lnd%ge, tileptr)
     call fptr(tileptr, ptr)
     if(associated(ptr)) ptr=data(i)
  enddo
end subroutine assemble_tiles_r0d

subroutine assemble_tiles_r0d_idx(fptr,n,idx,data)
  procedure(fptr_r0i) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: n ! additional index argument for fptr
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real,    intent(in) :: data(:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i

! distribute the data over the tiles
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                          lnd%ls,lnd%gs,lnd%ge, tileptr)
     call fptr(tileptr, n, ptr)
     if(associated(ptr)) ptr=data(i)
  enddo
end subroutine assemble_tiles_r0d_idx

subroutine assemble_tiles_i1d(fptr,idx,data)
  procedure(fptr_i0i) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  integer, intent(in) :: data(:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  integer, pointer :: ptr ! pointer to the tile data
  integer :: i,j

! distribute the data over the tiles
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                          lnd%ls,lnd%gs,lnd%ge, tileptr)
     do j = 1,size(data,2)
        call fptr(tileptr, j, ptr)
        if(associated(ptr)) ptr=data(i,j)
     enddo
  enddo
end subroutine assemble_tiles_i1d

subroutine assemble_tiles_r1d(fptr,idx,data)
  procedure(fptr_r0i) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real,    intent(in) :: data(:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i,j

! distribute the data over the tiles
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                          lnd%ls,lnd%gs,lnd%ge, tileptr)
     do j = 1,size(data,2)
        call fptr(tileptr, j, ptr)
        if(associated(ptr)) ptr=data(i,j)
     enddo
  enddo
end subroutine assemble_tiles_r1d

subroutine assemble_tiles_r1d_idx(fptr,n,idx,data)
  procedure(fptr_r0ij) :: fptr ! subroutine returning pointer to the data
  integer, intent(in) :: n ! additional index argument for fptr
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real,    intent(in) :: data(:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: i,j

! distribute the data over the tiles
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                          lnd%ls,lnd%gs,lnd%ge, tileptr)
     do j = 1,size(data,2)
        call fptr(tileptr, j, n, ptr)
        if(associated(ptr)) ptr=data(i,j)
     enddo
  enddo
end subroutine assemble_tiles_r1d_idx

subroutine assemble_tiles_r2d(fptr,idx,data)
  procedure(fptr_r0ij):: fptr ! subroutine returning the pointer to the data to be written
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real,    intent(in) :: data(:,:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real, pointer :: ptr ! pointer to the tile data
  integer :: i,k,m

! distribute the data over the tiles
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                          lnd%ls,lnd%gs,lnd%ge, tileptr)
     do k=1,size(data,2)
     do m=1,size(data,3)
        call fptr(tileptr, k, m, ptr)
        if(associated(ptr)) ptr=data(i,k,m)
     enddo
     enddo
  enddo
end subroutine assemble_tiles_r2d

subroutine assemble_tiles_r2d_idx(fptr,n,idx,data)
  procedure(fptr_r0ijk) :: fptr ! subroutine returning the pointer to the data
  integer, intent(in) :: n ! additional index argument for fptr
  integer, intent(in) :: idx(:)  ! local vector of tile indices
  real,    intent(in) :: data(:,:,:) ! local tile data

  ! ---- local vars
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real, pointer :: ptr ! pointer to the tile data
  integer :: i,k,m

! distribute the data over the tiles
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,land_tile_map,&
                          lnd%ls,lnd%gs,lnd%ge, tileptr)
     do k=1,size(data,2)
     do m=1,size(data,3)
        call fptr(tileptr, k, m, n, ptr)
        if(associated(ptr)) ptr=data(i,k,m)
     enddo
     enddo
  enddo
end subroutine assemble_tiles_r2d_idx

! ============================================================================
subroutine override_tile_data_r0d_fptr(fieldname,fptr,time,override)
  character(len=*), intent(in)   :: fieldname ! field to override
  procedure(fptr_r0) :: fptr ! subroutine returning pointer to the data
  type(time_type),  intent(in)   :: time      ! model time
  logical, optional, intent(out) :: override  ! true if the field has been
                                              ! overridden successfully

  ! ---- local vars
  real    :: data1D(lnd%ls:lnd%le)  ! storage for the input data
  logical :: override_

  call data_override_ug('LND',fieldname,data1D, time, override_ )
  if(present(override)) override=override_
  if(.not.override_) return ! do nothing if the field was not overridden

  ! distribute the data over the tiles
  call put_to_tiles_r0d_fptr(data1d,land_tile_map,fptr)

end subroutine override_tile_data_r0d_fptr

end module
