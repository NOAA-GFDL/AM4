module static_vegn_mod

use constants_mod,      only : pi
use mpp_mod,            only : mpp_max, mpp_sum
use mpp_io_mod, only : fieldtype, axistype, mpp_get_atts, mpp_open, MPP_RDONLY, &
      MPP_NETCDF, MPP_MULTI, MPP_SINGLE, mpp_get_axis_by_name, default_axis, &
      mpp_get_info, mpp_get_times, mpp_get_fields, mpp_get_axis_data, mpp_get_axis_data, &
      validtype, mpp_is_valid, mpp_get_time_axis
use fms_io_mod, only : restart_file_type, set_domain, nullify_domain, &
     get_file_name
use time_manager_mod,   only : time_type, set_date, time_type_to_real, &
     get_calendar_type, valid_calendar_types, operator(-), get_date
use get_cal_time_mod,   only : get_cal_time

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod,            only : error_mesg, FATAL, NOTE, &
     mpp_pe, file_exist, close_file, check_nml_error, stdlog, lowercase, &
     mpp_root_pe, get_mosaic_tile_file, fms_error_handler
use time_interp_mod,    only : time_interp
use diag_manager_mod,   only : get_base_date

use nf_utils_mod,       only : nfu_inq_dim, nfu_get_dim, nfu_def_dim, &
     nfu_inq_compressed_var, nfu_get_compressed_rec, nfu_validtype, &
     nfu_get_valid_range, nfu_is_valid, nfu_put_rec, nfu_put_att
use land_data_mod,      only : lnd_sg, log_version, lnd
use land_io_mod,        only : print_netcdf_error, new_land_io
use land_numerics_mod,  only : nearest
use land_tile_io_mod,   only : create_tile_out_file, gather_tile_index
use land_tile_mod,      only : land_tile_map, land_tile_type, land_tile_enum_type, first_elmt, &
     tail_elmt, next_elmt, current_tile, operator(/=), nitems
use vegn_cohort_mod,    only : vegn_cohort_type
use cohort_io_mod,      only : create_cohort_dimension_new, create_cohort_dimension_orig, gather_cohort_data, &
     write_cohort_data_i0d, write_cohort_data_r0d, gather_cohort_index

use fms_io_mod, only: fms_io_unstructured_register_restart_axis
use fms_io_mod, only: fms_io_unstructured_register_restart_field
use fms_io_mod, only: fms_io_unstructured_save_restart
use fms_io_mod, only: HIDX
use fms_io_mod, only: fms_io_unstructured_get_field_size
use fms_io_mod, only: fms_io_unstructured_read

implicit none
private

! ==== public interface =====================================================
public :: read_static_vegn_namelist
public :: static_vegn_init
public :: static_vegn_end

public :: read_static_vegn
public :: write_static_vegn
! ==== end of public interface ==============================================

! ==== module constants =====================================================
character(len=*), parameter :: module_name = 'static_vegn_mod'
#include "../shared/version_variable.inc"

! ==== module data ==========================================================
logical :: module_is_initialized = .FALSE.
integer :: ncid  ! netcdf id of the input file
integer, allocatable, dimension(:) :: species, status
real,    allocatable, dimension(:) :: bl, blv, br, bsw, bwood, bliving
type(time_type),allocatable :: time_line(:) ! time line of input data
type(time_type)             :: ts,te        ! beginning and end of time interval
integer, allocatable :: map_i(:), map_j(:)! remapping arrays: for each of the
     ! land grid cells in current domain they hold indices of corresponding points
     ! in the input grid.
type(time_type) :: base_time ! model base time for static vegetation output
type(fieldtype), allocatable :: Fields(:)
integer :: input_unit, ispecies, ibl, iblv, ibr, ibsw, ibwood, ibliving, istatus

type(restart_file_type) :: static_veg_file ! handle of output file, for new IO
integer :: ncid2 ! netcdf id of the output file, for old IO
integer :: tile_dim_length ! length of tile dimension in output files. global max of number of tiles per gridcell
integer, allocatable :: cidx(:) ! cohort compression index, local for current PE

! ---- namelist variables ---------------------------------------------------
logical :: use_static_veg = .FALSE.
character(len=512) :: input_file = & ! name of input file for static vegetation
     "INPUT/static_veg_data.nc"
character(len=10)  :: timeline   = 'normal' ! type of timeline ('normal' or 'loop')
integer, dimension(6) :: &
     start_loop = (/1,1,1,0,0,0/), & ! beginning of the time loop
     end_loop   = (/1,1,1,0,0,0/)    ! end of the time loop
logical :: fill_land_mask = .FALSE. ! if true, all the vegetation points on the
     ! map are filled with the information from static vegetation data, using
     ! nearest point remap; otherwise only the points that overlap with valid
     ! static vegetation data are overridden.
logical :: write_static_veg = .FALSE. ! if true, the state of vegetation is saved
     ! periodically for future use as static vegetation input
character(16) :: static_veg_freq = 'daily' ! or 'monthly', or 'annual'
     ! specifies the frequency for writing the static vegetation data file

namelist/static_veg_nml/use_static_veg,input_file,timeline,start_loop,end_loop,&
     fill_land_mask, write_static_veg, static_veg_freq

! ==== NetCDF declarations ==================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),__FILE__,__LINE__)

contains

! ===========================================================================
subroutine read_static_vegn_namelist(static_veg_used)
  logical, intent(out) :: static_veg_used

  ! ---- local vars
  integer :: unit, ierr, io

  call log_version(version, module_name, &
  __FILE__)

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=static_veg_nml, iostat=io)
  ierr = check_nml_error(io, 'static_veg_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=static_veg_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'static_veg_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif

  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=static_veg_nml)
  endif

  if (    (trim(static_veg_freq)=='daily') &
      .or.(trim(static_veg_freq)=='monthly') &
      .or.(trim(static_veg_freq)=='annual') ) then
     ! static_veg_freq is OK -- do nothing
  else
     call error_mesg('static_vegn_init','option static_veg_freq="'&
          //trim(static_veg_freq)&
          //'" is invalid, use "daily", "monthly", or "annual"', FATAL)
  endif

  static_veg_used = use_static_veg
end subroutine read_static_vegn_namelist


! ===========================================================================
subroutine static_vegn_init( )

  ! ---- local vars
  integer :: unlimdim, timelen, timeid
  integer :: i,j,k,iret
  character(len=NF_MAX_NAME) :: dimname  ! name of the dimension variable : time, lon, and lat
  integer                    :: ndims    ! rank of input vars
  integer                    :: dimids (NF_MAX_VAR_DIMS) ! netcdf IDs of input var dimensions
  integer                    :: dimlens(NF_MAX_VAR_DIMS) ! sizes of respective dimensions
  real, allocatable          :: t(:)     ! temporary real timeline
  character(len=256)         :: units    ! units of time in the file
  character(len=256)         :: calendar ! calendar of the data
  real, allocatable          :: in_lon(:)! longitude coordinates in input file
  real, allocatable          :: in_lat(:)! latitude coordinates in input file
  logical, allocatable       :: mask(:,:)! mask of valid points in input data
  integer, allocatable       :: data(:,:,:,:) ! temporary array used to calculate the mask of
                                         ! valid input data
  logical                    :: has_records ! true if input variable has records
  integer :: m, n, siz(4), ndim, nvar, natt, l
  character(len=1024) :: actual_input_file, actual_input_file2
  logical :: input_is_multiface ! TRUE if the input files are face-specific
  logical :: found_file, read_dist, io_domain_exist
  type(axistype) :: Lon_axis, Lat_axis, Tile_axis, Cohort_axis
  type(axistype) :: Time_axis
  character(len=256) :: name
  integer, allocatable :: cidx(:), idata(:)

  if(module_is_initialized) return

  if(use_static_veg) then

     ! SET UP LOOP BOUNDARIES
     ts = set_date(start_loop(1),start_loop(2),start_loop(3), start_loop(4),start_loop(5),start_loop(6))
     te = set_date(end_loop(1)  ,end_loop(2)  ,end_loop(3)  , end_loop(4)  ,end_loop(5)  ,end_loop(6)  )

     if(new_land_io) then
     ! OPEN INPUT FILE
       if(file_exist(trim(input_file),no_domain=.true.)) then
          call mpp_open(input_unit, trim(input_file), action=MPP_RDONLY, form=MPP_NETCDF, &
                     threading=MPP_MULTI, fileset=MPP_SINGLE)
          call error_mesg('static_vegn_init','Reading global static vegetation file "'&
               //trim(input_file)//'"', NOTE)
          input_is_multiface = .FALSE.
          actual_input_file = input_file
       else
          if(lnd%nfaces==1) then
             ! for 1-face grid we cannot use multi-face input, even if it exists
             call error_mesg('static_vegn_init','input file "'//trim(input_file)&
                     //'" does not exist', FATAL)
          else
             ! if there is more then one face, try opening face-specific input with consideration of io_layout
             call get_mosaic_tile_file(trim(input_file),actual_input_file2,.FALSE.,lnd_sg%domain)
             found_file = get_file_name(input_file, actual_input_file, read_dist, io_domain_exist, &
                             domain=lnd_sg%domain)
             if(.not.found_file) call error_mesg('static_vegn_init','"'//trim(actual_input_file2)// &
                '" and corresponding distributed file are not found', FATAL)
             if(read_dist) then
                call mpp_open(input_unit, trim(actual_input_file), action=MPP_RDONLY, form=MPP_NETCDF, &
                                 threading=MPP_MULTI, fileset=MPP_MULTI, domain=lnd_sg%domain)
             else
                call mpp_open(input_unit, trim(actual_input_file), action=MPP_RDONLY, form=MPP_NETCDF, &
                              threading=MPP_MULTI, fileset=MPP_SINGLE)
             endif
             call error_mesg('static_vegn_init','Reading face-specific vegetation file "'&
                  //trim(actual_input_file)//'"', NOTE)
             input_is_multiface = .TRUE.
          endif
       endif

       call mpp_get_info(input_unit,ndim,nvar,natt,timelen)
       call mpp_get_time_axis(input_unit, Time_axis)

       ! READ TIME AXIS DATA
       allocate(t(timelen))
       call mpp_get_times(input_unit, t)

       allocate(Fields(nvar))
       call mpp_get_fields(input_unit, Fields)
       do i=1,nvar
         call mpp_get_atts(Fields(i), name=name)
         select case (name)
         case ('species')
           ispecies = i
         case ('bl')
           ibl = i
         case ('blv')
           iblv = i
         case ('br')
           ibr = i
         case ('bsw')
           ibsw = i
         case ('bwood')
           ibwood = i
         case ('bliving')
           ibliving = i
         case ('status')
           istatus = i
         end select
       enddo

       ! GET UNITS OF THE TIME AND CALENDAR OF THE DATA
       units = ' '
       calendar = 'JULIAN'
       call mpp_get_atts(Time_axis, units=units, calendar=calendar)

       ! CONVERT TIME TO THE FMS TIME_TYPE AND STORE IT IN THE TIMELINE FOR THE DATA SET
       allocate(time_line(timelen))
       do i = 1, size(t)
          ! set the respective value in the timeline
          time_line(i) = get_cal_time(t(i),units,calendar)
       enddo

       ! READ HORIZONTAL COORDINATES
       Lon_axis = mpp_get_axis_by_name(input_unit,'lon')
       call mpp_get_atts(Lon_axis, len=dimlens(1))
       allocate(in_lon(dimlens(1)))
       call mpp_get_axis_data(Lon_axis, in_lon)

       Lat_axis = mpp_get_axis_by_name(input_unit,'lat')
       call mpp_get_atts(Lat_axis, len=dimlens(2))
       allocate(in_lat(dimlens(2)))
       call mpp_get_axis_data(Lat_axis, in_lat)

       in_lon = in_lon*PI/180.0 ; in_lat = in_lat*PI/180.0

       ! COMPUTE INDEX REMAPPING ARRAY
       allocate(map_i(lnd%ls:lnd%le))
       allocate(map_j(lnd%ls:lnd%le))
       map_i = -1
       map_j = -1
       if( .not. input_is_multiface ) then
          allocate(mask(size(in_lon),size(in_lat)))
          mask = .false.

          if(fill_land_mask) then
             ! READ THE FIRST RECORD AND CALCULATE THE MASK OF THE VALID INPUT DATA
             Tile_axis   = mpp_get_axis_by_name(input_unit,'tile')
             call mpp_get_atts(Tile_axis, len=dimlens(3))
             Cohort_axis = mpp_get_axis_by_name(input_unit,'cohort')
             call mpp_get_atts(Cohort_axis, len=dimlens(4))
             ! Note: The input file used for initial testing had
             ! lon = 144, lat = 90, tile = 2, cohort = 1

             call fms_io_unstructured_get_field_size(trim(input_file), "cohort_index", siz, &
                                                     lnd%domain)
             allocate(cidx(siz(1)), idata(siz(1)))

             call fms_io_unstructured_read(trim(input_file), "cohort_index", cidx, lnd%domain, timelevel=1)
             call fms_io_unstructured_read(trim(input_file), "species", idata, lnd%domain, timelevel=1)
             do n = 1,size(cidx)
                m = cidx(n)
                i = modulo(m,dimlens(1))+1
                m = m/dimlens(1)
                j = modulo(m,dimlens(2))+1
                m = m/dimlens(2)
                ! k = modulo(m,dimlens(3))+1 ! This is how to get tile number, if it were needed.
                m = m/dimlens(3)
                   ! L = m+1  ! This is how to get cohort number, if it were needed. No need to do
                              ! modulo with dimlens(4) because at this point m is always < dimlens(4)
                if(idata(n)>=0 .or. mask(i,j)) then
                   mask(i,j) = .TRUE. ! If species exists in any cohort of this grid cell then mask is .TRUE.
                endif
             enddo
             deallocate(idata)
          else
             mask(:,:) = .TRUE.
          endif
       endif
     else ! original code below here. i.e. if(.not.new_land_io)
       ! OPEN INPUT FILE
       if (nf_open(input_file,NF_NOWRITE,ncid)/=NF_NOERR) then
          if(lnd%nfaces==1) then
             ! for 1-face grid we cannot use multi-face input, even if it exists
             call error_mesg('static_vegn_init','input file "'//trim(input_file)&
                     //'" does not exist', FATAL)
          else
             ! if there is more then one face, try opening face-specific input
             call get_mosaic_tile_file(trim(input_file),actual_input_file,lnd%domain)
             if (nf_open(actual_input_file,NF_NOWRITE,ncid)/=NF_NOERR) then
                call error_mesg('static_vegn_init','Neither "'//trim(input_file)&
                     //'" nor "'//trim(actual_input_file)//'" files exist', FATAL)
             else
                call error_mesg('static_vegn_init','Reading face-specific vegetation file "'&
                     //trim(actual_input_file)//'"', NOTE)
                input_is_multiface = .TRUE.
             endif
          endif
       else
          call error_mesg('static_vegn_init','Reading global static vegetation file "'&
               //trim(input_file)//'"', NOTE)
          input_is_multiface = .FALSE.
          actual_input_file = input_file
       endif

       ! READ TIME AXIS DATA
       if(nf_inq_unlimdim( ncid, unlimdim )/=NF_NOERR) then
          call error_mesg('static_vegn_init',&
              'Input file "'//trim(actual_input_file)//'" lacks record dimension.', FATAL)
       endif
       __NF_ASRT__(nf_inq_dimname ( ncid, unlimdim, dimname ))
       __NF_ASRT__(nf_inq_varid   ( ncid, dimname, timeid ))
       __NF_ASRT__(nf_inq_dimlen( ncid, unlimdim, timelen ))
       allocate (time_line(timelen), t(timelen))
       __NF_ASRT__(nf_get_var_double (ncid, timeid, t ))

       ! GET UNITS OF THE TIME
       units = ' '
       if (nf_get_att_text(ncid, timeid,'units',units)/=NF_NOERR) then
          call error_mesg('static_vegn_init',&
              'Cannot read time  units from file "'//trim(actual_input_file)//'"', FATAL)
       endif

       ! GET CALENDAR OF THE DATA
       calendar = ' '
       iret = nf_get_att_text(ncid, timeid, 'calendar',calendar)
       if(iret/=NF_NOERR) &
            iret = nf_get_att_text(ncid, timeid,'calendar_type',calendar)
       if(iret/=NF_NOERR) &
            calendar='JULIAN' ! use model calendar? how to get the name of the model calendar?

       ! CONVERT TIME TO THE FMS TIME_TYPE AND STORE IT IN THE TIMELINE FOR THE
       ! DATA SET
       do i = 1, size(t)
          ! set the respective value in the timeline
          time_line(i) = get_cal_time(t(i),units,calendar)
       enddo

       ! READ HORIZONTAL COORDINATES
       iret = nfu_inq_compressed_var(ncid,'species',ndims=ndims,dimids=dimids,dimlens=dimlens,&
            has_records=has_records)
       if (iret/=NF_NOERR) then
          call error_mesg('static_vegn_init',&
             'Cannot read compression information from file "'//trim(actual_input_file)//&
             '": check that all dimensions listed in "compress" attributes are present in the file.', FATAL)
       endif
       __NF_ASRT__(iret)
       allocate(in_lon(dimlens(1)),in_lat(dimlens(2)))
       __NF_ASRT__(nfu_get_dim(ncid,dimids(1),in_lon)) ! get longitude
       __NF_ASRT__(nfu_get_dim(ncid,dimids(2),in_lat)) ! get latitude
       in_lon = in_lon*PI/180.0 ; in_lat = in_lat*PI/180.0

       ! COMPUTE INDEX REMAPPING ARRAY
       allocate(map_i(lnd%ls:lnd%le))
       allocate(map_j(lnd%ls:lnd%le))
       allocate(mask(size(in_lon),size(in_lat)))

       map_i = -1
       map_j = -1
       mask = .false.

       if(fill_land_mask) then
          ! CALCULATE THE DIMENSIONS OF THE BUFFER FOR THE INPUT DATA
          if (has_records) ndims=ndims-1
          do i = ndims+1,4
             dimlens(i) = 1
          enddo
          ! READ THE FIRST RECORD AND CALCULATE THE MASK OF THE VALID INPUT DATA
          allocate(data(dimlens(1),dimlens(2),dimlens(3),dimlens(4)))
          !             lon        lat        tile       cohort
          data(:,:,:,:) = -1
          __NF_ASRT__(nfu_get_compressed_rec(ncid,'species',1,data))
          do j = 1,size(data,2)
          do i = 1,size(data,1)
             mask(i,j) = any(data(i,j,:,:)>=0)
          enddo
          enddo
          deallocate(data)
       else
          mask(:,:) = .TRUE.
       endif
    endif

    if(input_is_multiface) then
       ! check that the sizes of input data and the model data are the same
       if(dimlens(1)/=lnd%nlon.or.dimlens(2)/=lnd%nlat) then
          call error_mesg('static_vegn_init','size of face-specific static vegetation '&
               //'data isn''t the same as the size of the mosaic face', FATAL)
       endif
       ! in case of multi-face input, we do not do any remapping
       do l = lnd%ls,lnd%le
          map_i(l) = lnd%i_index(l); map_j(l) = lnd%j_index(l)
       enddo
    else
       ! do the nearest-point remapping
       do l = lnd%ls,lnd%le
          call nearest(mask,in_lon,in_lat,lnd%lon(l),lnd%lat(l),map_i(l),map_j(l))
       enddo
    endif
    deallocate (in_lon,in_lat)
    if(allocated(mask)) deallocate(mask)
    deallocate(t)
  endif

  if(write_static_veg) &
      call init_writing_static_veg()

  module_is_initialized = .true.
end subroutine static_vegn_init


! create output file for static vegetation
subroutine init_writing_static_veg()
  integer, allocatable :: tidx(:)
  integer :: k, l, csize, iret
  integer :: year, month, day, hour, minute, sec ! components of base date
  character(len=256) :: units ! units of time in the output file

  ! count all land tiles and determine the length of tile dimension
  ! sufficient for the current domain
  tile_dim_length = 0
  do l = lnd%ls, lnd%le
     k = nitems(land_tile_map(l))
     tile_dim_length = max(tile_dim_length,k)
  enddo

  ! [1.1] calculate the tile dimension length by taking the max across all domains
  call mpp_max(tile_dim_length)
  call gather_tile_index(vegn_tile_exists,tidx)
  call gather_cohort_index(tile_dim_length, cidx)

  csize = size(cidx)
  allocate(species(csize),status(csize),bl(csize), blv(csize), br(csize), bsw(csize), &
           bwood(csize), bliving(csize))

  call gather_tile_index(vegn_tile_exists, tidx)

  call get_base_date(year,month,day,hour,minute,sec)
  base_time = set_date(year, month, day, hour, minute, sec)
  write(units, 11) year, month, day, hour, minute, sec

  if(new_land_io) then
     call create_tile_out_file(static_veg_file, 'static_veg_out.nc', tidx, tile_dim_length)
     call create_cohort_dimension_new(static_veg_file, cidx, 'static_veg_out.nc', tile_dim_length)
     call fms_io_unstructured_register_restart_axis(static_veg_file, "static_veg_out.nc", &
            "time",(/0.0/), "T", lnd%domain, units=units, calendar=valid_calendar_types(get_calendar_type()))
     k = fms_io_unstructured_register_restart_field(static_veg_file, "static_veg_out.nc", &
            "species", species, (/HIDX/), lnd%domain, longname="vegetation species")
     k = fms_io_unstructured_register_restart_field(static_veg_file, "static_veg_out.nc", &
            "bl", bl, (/HIDX/), lnd%domain, longname="biomass of leaves per individual", units="kg C/m2")
     k = fms_io_unstructured_register_restart_field(static_veg_file, "static_veg_out.nc", &
            "blv", blv, (/HIDX/), lnd%domain, longname="biomass of virtual leaves (labile store) per individual", units="kg C/m2")
     k = fms_io_unstructured_register_restart_field(static_veg_file, "static_veg_out.nc", &
            "br", br, (/HIDX/), lnd%domain, longname="biomass of fine roots per individual", units="kg C/m2")
     k = fms_io_unstructured_register_restart_field(static_veg_file, "static_veg_out.nc", &
            "bsw", bsw, (/HIDX/), lnd%domain, longname="biomass of sapwood per individual", units="kg C/m2")
     k = fms_io_unstructured_register_restart_field(static_veg_file, "static_veg_out.nc", &
            "bwood", bwood, (/HIDX/), lnd%domain, longname="biomass of heartwood per individual", units="kg C/m2")
     k = fms_io_unstructured_register_restart_field(static_veg_file, "static_veg_out.nc", &
            "bliving", bliving, (/HIDX/), lnd%domain, longname="total living biomass per individual", units="")
     k = fms_io_unstructured_register_restart_field(static_veg_file, "static_veg_out.nc", &
            "status", status, (/HIDX/), lnd%domain, longname="leaf status", units="")
     call fms_io_unstructured_save_restart(static_veg_file, directory="", time_level=-1.0)
  else
     call create_tile_out_file(ncid2,'static_veg_out.nc', tidx, tile_dim_length, &
                               lnd%coord_glon, lnd%coord_glat)
     ! create compressed dimension for vegetation cohorts
     call create_cohort_dimension_orig(ncid2, cidx, tile_dim_length)
     ! get the base date of the simulation
     if(mpp_pe()==lnd%io_pelist(1)) then
        ! create time axis, on root IO processors only
        __NF_ASRT__(nfu_def_dim(ncid2,'time',NF_UNLIMITED,NF_DOUBLE,units=trim(units)))
        ! add calendar attribute to the time axis
        iret=nfu_put_att(ncid2,'time','calendar',trim(valid_calendar_types(get_calendar_type())))
        __NF_ASRT__(iret)
     endif
  endif
11 format('days since ', i4.4, '-', i2.2, '-', i2.2, ' ', i2.2, ':', i2.2, ':', i2.2)
end subroutine init_writing_static_veg

! ===========================================================================
subroutine static_vegn_end()

  if(new_land_io) return

  if(use_static_veg) then
     __NF_ASRT__(nf_close(ncid))
     deallocate(time_line,map_i,map_j)
  endif
  if(write_static_veg .and.  mpp_pe()==lnd%io_pelist(1) ) then
     __NF_ASRT__(nf_close(ncid2))
  endif
  module_is_initialized = .false.
end subroutine static_vegn_end

! ===========================================================================
subroutine read_static_vegn (time, err_msg)
  type(time_type), intent(in)    :: time
  character(len=*), intent(out), optional :: err_msg

  ! ---- local vars
  integer :: index1, index2 ! result of time interpolation (only index1 is used)
  real    :: weight         ! another result of time interp, not used
  character(len=256) :: msg
  integer :: siz(4)
  integer, allocatable :: cidx(:), idata(:)
  real,    allocatable :: rdata(:)

  if(.not.use_static_veg)return;

  msg = ''
  !   time_interp to find out the index of the current time interval
  if (timeline == 'loop') then
     call time_interp(time, ts, te, time_line, weight, index1, index2, &
                      correct_leap_year_inconsistency=.true., err_msg=msg)
  else if (timeline == 'normal') then
     call time_interp(time, time_line, weight, index1, index2, err_msg=msg)
  else
     call error_mesg(module_name,'timeline option "'//trim(timeline)// &
          '" is incorrect, use "normal" or "loop"', FATAL)
  endif
  if(msg /= '') then
    if(fms_error_handler('read_static_vegn','Message from time_interp: '//trim(msg),err_msg)) return
  endif

  ! read the data into cohort variables
  if(new_land_io) then
     call fms_io_unstructured_get_field_size(trim(input_file), "cohort_index", siz, lnd%domain)
     allocate(cidx(siz(1)), idata(siz(1)), rdata(siz(1)))
     call fms_io_unstructured_read(trim(input_file), "cohort_index", cidx, lnd%domain, timelevel=index1)
     call fms_io_unstructured_read(trim(input_file), "species", idata, lnd%domain, timelevel=index1)
     call read_remap_cohort_data_i0d_new(Fields(ispecies), cohort_species_ptr, map_i, map_j, cidx, idata)
     call fms_io_unstructured_read(trim(input_file), "bl", rdata, lnd%domain, timelevel=index1)
     call read_remap_cohort_data_r0d_new(Fields(ibl), cohort_bl_ptr, map_i, map_j, cidx, rdata)
     call fms_io_unstructured_read(trim(input_file), "blv", rdata, lnd%domain, timelevel=index1)
     call read_remap_cohort_data_r0d_new(Fields(iblv), cohort_blv_ptr, map_i, map_j, cidx, rdata)
     call fms_io_unstructured_read(trim(input_file), "br", rdata, lnd%domain, timelevel=index1)
     call read_remap_cohort_data_r0d_new(Fields(ibr), cohort_br_ptr, map_i, map_j, cidx, rdata)
     call fms_io_unstructured_read(trim(input_file), "bsw", rdata, lnd%domain, timelevel=index1)
     call read_remap_cohort_data_r0d_new(Fields(ibsw), cohort_bsw_ptr, map_i, map_j, cidx, rdata)
     call fms_io_unstructured_read(trim(input_file), "bwood", rdata, lnd%domain, timelevel=index1)
     call read_remap_cohort_data_r0d_new(Fields(ibwood), cohort_bwood_ptr, map_i, map_j, cidx, rdata)
     call fms_io_unstructured_read(trim(input_file), "bliving", rdata, lnd%domain, timelevel=index1)
     call read_remap_cohort_data_r0d_new(Fields(ibliving), cohort_bliving_ptr, map_i, map_j, cidx, rdata)
     call fms_io_unstructured_read(trim(input_file), "status", idata, lnd%domain, timelevel=index1)
     call read_remap_cohort_data_i0d_new(Fields(istatus), cohort_status_ptr, map_i, map_j, cidx, idata)
     deallocate(cidx, idata, rdata)
  else
     call read_remap_cohort_data_i0d_fptr(ncid, 'species' , cohort_species_ptr , map_i, map_j, index1)
     call read_remap_cohort_data_r0d_fptr(ncid, 'bl'      , cohort_bl_ptr      , map_i, map_j, index1)
     call read_remap_cohort_data_r0d_fptr(ncid, 'blv'     , cohort_blv_ptr     , map_i, map_j, index1)
     call read_remap_cohort_data_r0d_fptr(ncid, 'br'      , cohort_br_ptr      , map_i, map_j, index1)
     call read_remap_cohort_data_r0d_fptr(ncid, 'bsw'     , cohort_bsw_ptr     , map_i, map_j, index1)
     call read_remap_cohort_data_r0d_fptr(ncid, 'bwood'   , cohort_bwood_ptr   , map_i, map_j, index1)
     call read_remap_cohort_data_r0d_fptr(ncid, 'bliving' , cohort_bliving_ptr , map_i, map_j, index1)
     call read_remap_cohort_data_i0d_fptr(ncid, 'status'  , cohort_status_ptr  , map_i, map_j, index1)
  endif

  ! derived variables will be updated in update_land_bc_fast
end subroutine read_static_vegn


! ===========================================================================
subroutine write_static_vegn()

  real :: t ! time in output units
  integer :: rec ! number of record to write
  ! components of the date
  integer :: second, minute, hour, day0, day1, month0, month1, year0, year1

  if(.not.write_static_veg) return;

  ! get components of calendar dates for this and previous time step
  call get_date(lnd%time,             year0,month0,day0,hour,minute,second)
  call get_date(lnd%time-lnd%dt_fast, year1,month1,day1,hour,minute,second)

  if (.not.((trim(static_veg_freq)=='daily'  .and.  day1/=day0)   &
        .or.(trim(static_veg_freq)=='monthly'.and.month1/=month0) &
        .or.(trim(static_veg_freq)=='annual' .and. year1/=year0))) return

  t = (time_type_to_real(lnd%time)-time_type_to_real(base_time))/86400
  call gather_cohort_data(cohort_species_ptr,cidx,tile_dim_length,species)
  call gather_cohort_data(cohort_bl_ptr,cidx,tile_dim_length,bl)
  call gather_cohort_data(cohort_blv_ptr,cidx,tile_dim_length,blv)
  call gather_cohort_data(cohort_br_ptr,cidx,tile_dim_length,br)
  call gather_cohort_data(cohort_bsw_ptr,cidx,tile_dim_length,bsw)
  call gather_cohort_data(cohort_bwood_ptr,cidx,tile_dim_length,bwood)
  call gather_cohort_data(cohort_bliving_ptr,cidx,tile_dim_length,bliving)
  call gather_cohort_data(cohort_status_ptr,cidx,tile_dim_length,status)
  if(new_land_io) then
     call fms_io_unstructured_save_restart(static_veg_file, directory="", append=.true., time_level=t)
  else
     ! get the current number of records in the output file, rec is only needed by the io_pelist root pe.
     rec = 0
     ! create new record in the output file and store current value of time
     if(mpp_pe()==lnd%io_pelist(1)) then
          __NF_ASRT__(nfu_inq_dim(ncid2,'time',rec))
          rec = rec+1
        __NF_ASRT__(nfu_put_rec(ncid2,'time',rec,t))
     endif
     ! write static vegetation data
     call write_cohort_data_i0d(ncid2,'species', species, 'vegetation species',record=rec)
     call write_cohort_data_r0d(ncid2,'bl', bl, &
          'biomass of leaves per individual','kg C/m2', record=rec)
     call write_cohort_data_r0d(ncid2,'blv', blv, &
          'biomass of virtual leaves (labile store) per individual','kg C/m2',record=rec)
     call write_cohort_data_r0d(ncid2,'br', br, &
          'biomass of fine roots per individual','kg C/m2', record=rec)
     call write_cohort_data_r0d(ncid2,'bsw', bsw, &
          'biomass of sapwood per individual','kg C/m2', record=rec)
     call write_cohort_data_r0d(ncid2,'bwood', bwood, &
          'biomass of heartwood per individual','kg C/m2', record=rec)
     call write_cohort_data_r0d(ncid2,'bliving', bliving, &
          'total living biomass per individual','kg C/m2', record=rec)
     call write_cohort_data_i0d(ncid2,'status',  status, &
          'leaf status', record=rec)
  endif
end subroutine write_static_vegn


! ============================================================================
#define F90_TYPE       integer
#define READ_REMAP_SUB read_remap_cohort_data_i0d_fptr
#include "read_remap_cohort_data.inc"

#define F90_TYPE       real
#define READ_REMAP_SUB read_remap_cohort_data_r0d_fptr
#include "read_remap_cohort_data.inc"

! ============================================================================
#define F90_TYPE       integer
#define READ_REMAP_SUB read_remap_cohort_data_i0d_new
#include "read_remap_cohort_data_new.inc"

#define F90_TYPE       real
#define READ_REMAP_SUB read_remap_cohort_data_r0d_new
#include "read_remap_cohort_data_new.inc"
! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function vegn_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   vegn_tile_exists = associated(tile%vegn)
end function vegn_tile_exists

! ============================================================================
! cohort accessor functions: given a pointer to cohort, return a pointer to a
! specific member of the cohort structure
#define DEFINE_COHORT_ACCESSOR(xtype,x) subroutine cohort_ ## x ## _ptr(c,p);\
type(vegn_cohort_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%x;\
end subroutine

DEFINE_COHORT_ACCESSOR(integer,species)
DEFINE_COHORT_ACCESSOR(real,bl)
DEFINE_COHORT_ACCESSOR(real,br)
DEFINE_COHORT_ACCESSOR(real,blv)
DEFINE_COHORT_ACCESSOR(real,bsw)
DEFINE_COHORT_ACCESSOR(real,bwood)
DEFINE_COHORT_ACCESSOR(real,bliving)
DEFINE_COHORT_ACCESSOR(integer,status)

end module static_vegn_mod
